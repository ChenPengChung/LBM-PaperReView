!======================================================================
! Periodic Hill Structured Grid Generator
! Using TTM (Thompson-Thames-Mastin) Poisson Equation Solver
!
! Based on the ERCOFTAC UFR 3-30 Periodic Hill geometry
! Hill profile from Almeida et al. (1993), Mellen et al. (2000)
!
! Polynomial variable: physical x coordinate (NOT x/h)
! Breakpoints: x = 0, 9, 14, 20, 30, 40, 54
! Symmetry about x = Lx/2 = 126
!======================================================================
program periodic_hill_grid
  implicit none
  integer, parameter :: dp = selected_real_kind(15)

  ! ---- Geometry parameters ----
  real(dp), parameter :: h_hill = 28.0_dp
  real(dp), parameter :: Lx = 9.0_dp * h_hill        ! 252
  real(dp), parameter :: Ly = 3.036_dp * h_hill       ! 85.008
  real(dp), parameter :: Lx_half = 0.5_dp * Lx        ! 126

  ! ---- Grid dimensions ----
  integer, parameter :: Nx = 160
  integer, parameter :: Ny = 80

  ! ---- Grid arrays ----
  real(dp) :: xg(Nx, Ny), yg(Nx, Ny)

  ! ---- Bottom/top wall x-coordinates ----
  real(dp) :: x_bot(Nx), y_bot(Nx)
  real(dp) :: x_top(Nx)

  ! ---- Source term arrays ----
  real(dp) :: Qsrc(Nx, Ny)

  ! ---- Solver parameters ----
  real(dp), parameter :: omega_sor = 1.4_dp
  integer,  parameter :: max_iter  = 20000
  real(dp), parameter :: tol_conv  = 1.0e-6_dp
  real(dp), parameter :: beta_wall = 2.5_dp   ! tanh stretching param

  write(*,*) '================================================='
  write(*,*) ' Periodic Hill Structured Grid Generator'
  write(*,*) ' TTM Poisson Equation Solver'
  write(*,*) '================================================='
  write(*,'(A,I4,A,I4)') '  Grid: Nx=', Nx, '  Ny=', Ny
  write(*,'(A,F8.2,A,F8.2)') '  Domain: [0,', Lx, '] x [0,', Ly, ']'
  write(*,'(A,F8.2)') '  Hill height h = ', h_hill
  write(*,*) ''

  ! Step 1: Compute bottom wall (arc-length parameterization)
  call compute_bottom_wall()

  ! Step 2: Set up all four boundaries
  call setup_boundaries()

  ! Step 3: TFI (Transfinite Interpolation) initialization
  call tfi_init()

  ! Step 4: Compute source terms from TFI grid
  call compute_source_terms()

  ! Step 5: Poisson (TTM) iterative solver
  call poisson_solve()

  ! Step 6: Check grid quality
  call check_jacobian()

  ! Step 7: Output files
  call write_plot3d()
  call write_grid_dat()

  write(*,*) ''
  write(*,*) 'Grid generation complete!'
  write(*,*) '  Plot3D file: periodic_hill.xyz'
  write(*,*) '  Data file:   grid_data.dat'

contains

  !--------------------------------------------------------------------
  ! Hill profile function: returns y given physical x in [0, Lx]
  ! Uses correct ERCOFTAC polynomial coefficients (physical x variable)
  !--------------------------------------------------------------------
  function hill_profile(xp) result(yh)
    real(dp), intent(in) :: xp
    real(dp) :: yh
    real(dp) :: xx

    ! Apply symmetry about Lx/2
    xx = xp
    if (xx > Lx_half) xx = Lx - xx

    ! Piecewise polynomial (6 pieces + flat)
    if (xx < 0.0_dp) then
      yh = h_hill
    else if (xx < 9.0_dp) then
      yh = min(28.0_dp, 2.800000000000e+01_dp &
           + 0.000000000000e+00_dp * xx &
           + 6.775070969851e-03_dp * xx**2 &
           - 2.124527775800e-03_dp * xx**3)
    else if (xx < 14.0_dp) then
      yh = 2.507355893131e+01_dp &
           + 9.754803562315e-01_dp * xx &
           - 1.016116352781e-01_dp * xx**2 &
           + 1.889794677828e-03_dp * xx**3
    else if (xx < 20.0_dp) then
      yh = 2.579601052357e+01_dp &
           + 8.206693007457e-01_dp * xx &
           - 9.055370274339e-02_dp * xx**2 &
           + 1.626510569859e-03_dp * xx**3
    else if (xx < 30.0_dp) then
      yh = 4.046435022819e+01_dp &
           - 1.379581654948e+00_dp * xx &
           + 1.945884504128e-02_dp * xx**2 &
           - 2.070318932190e-04_dp * xx**3
    else if (xx < 40.0_dp) then
      yh = 1.792461334664e+01_dp &
           + 8.743920332081e-01_dp * xx &
           - 5.567361123058e-02_dp * xx**2 &
           + 6.277731764683e-04_dp * xx**3
    else if (xx < 54.0_dp) then
      yh = max(0.0_dp, 5.639011190988e+01_dp &
           - 2.010520359035e+00_dp * xx &
           + 1.644919857549e-02_dp * xx**2 &
           + 2.674976141766e-05_dp * xx**3)
    else
      yh = 0.0_dp
    end if

  end function hill_profile

  !--------------------------------------------------------------------
  ! Compute bottom wall points using arc-length parameterization
  ! This naturally clusters points where the hill slope is steep
  !--------------------------------------------------------------------
  subroutine compute_bottom_wall()
    integer, parameter :: Nfine = 20000
    real(dp) :: xf(Nfine), sf(Nfine)
    real(dp) :: dx_fine, dy, s_total, s_target
    integer :: i, k

    ! Compute arc length along bottom wall at fine resolution
    dx_fine = Lx / real(Nfine - 1, dp)
    xf(1) = 0.0_dp
    sf(1) = 0.0_dp

    do k = 2, Nfine
      xf(k) = real(k - 1, dp) * dx_fine
      dy = hill_profile(xf(k)) - hill_profile(xf(k-1))
      sf(k) = sf(k-1) + sqrt(dx_fine**2 + dy**2)
    end do

    s_total = sf(Nfine)
    write(*,'(A,F10.3)') '  Bottom wall arc length = ', s_total

    ! Distribute Nx points uniformly in arc length
    x_bot(1)  = 0.0_dp
    x_bot(Nx) = Lx

    do i = 2, Nx - 1
      s_target = real(i - 1, dp) / real(Nx - 1, dp) * s_total

      ! Binary search for efficiency
      do k = 2, Nfine
        if (sf(k) >= s_target) then
          x_bot(i) = xf(k-1) + (s_target - sf(k-1)) &
                     / (sf(k) - sf(k-1)) * dx_fine
          exit
        end if
      end do
    end do

    ! Compute y coordinates on bottom wall
    do i = 1, Nx
      y_bot(i) = hill_profile(x_bot(i))
    end do

    ! Top wall: uniform x distribution
    do i = 1, Nx
      x_top(i) = real(i - 1, dp) / real(Nx - 1, dp) * Lx
    end do

    write(*,'(A,F10.4,A,F10.4)') '  x_bot range: [', x_bot(1), ', ', x_bot(Nx), ']'

  end subroutine compute_bottom_wall

  !--------------------------------------------------------------------
  ! Set up all four boundary conditions
  !--------------------------------------------------------------------
  subroutine setup_boundaries()
    integer :: i, j
    real(dp) :: eta, s
    real(dp) :: ybot_0

    ybot_0 = hill_profile(0.0_dp)  ! = h_hill = 28

    ! Bottom boundary (j = 1): hill profile
    do i = 1, Nx
      xg(i, 1) = x_bot(i)
      yg(i, 1) = y_bot(i)
    end do

    ! Top boundary (j = Ny): flat at y = Ly
    do i = 1, Nx
      xg(i, Ny) = x_top(i)
      yg(i, Ny) = Ly
    end do

    ! Left boundary (i = 1): vertical from (0, h) to (0, Ly)
    do j = 1, Ny
      eta = real(j - 1, dp) / real(Ny - 1, dp)
      s = tanh_stretch_wall(eta)
      xg(1, j) = 0.0_dp
      yg(1, j) = ybot_0 + s * (Ly - ybot_0)
    end do
    yg(1, 1)  = ybot_0  ! exact bottom
    yg(1, Ny) = Ly       ! exact top

    ! Right boundary (i = Nx): periodic (same y as left, x = Lx)
    do j = 1, Ny
      xg(Nx, j) = Lx
      yg(Nx, j) = yg(1, j)
    end do

    write(*,*) '  Boundaries initialized.'

  end subroutine setup_boundaries

  !--------------------------------------------------------------------
  ! Tanh stretching: clusters near eta=0 (bottom wall)
  ! Maps [0,1] -> [0,1] with fine spacing near 0
  !--------------------------------------------------------------------
  function tanh_stretch_wall(eta) result(s)
    real(dp), intent(in) :: eta
    real(dp) :: s
    s = 1.0_dp - tanh(beta_wall * (1.0_dp - eta)) / tanh(beta_wall)
  end function tanh_stretch_wall

  !--------------------------------------------------------------------
  ! Transfinite Interpolation (TFI) for initial grid
  !--------------------------------------------------------------------
  subroutine tfi_init()
    integer :: i, j
    real(dp) :: xi, eta

    do j = 2, Ny - 1
      do i = 2, Nx - 1
        xi  = real(i - 1, dp) / real(Nx - 1, dp)
        eta = real(j - 1, dp) / real(Ny - 1, dp)

        xg(i,j) = (1.0_dp - eta) * xg(i,1) + eta * xg(i,Ny) &
                 + (1.0_dp - xi)  * xg(1,j) + xi  * xg(Nx,j) &
                 - (1.0_dp - xi) * (1.0_dp - eta) * xg(1,1)   &
                 - xi * (1.0_dp - eta) * xg(Nx,1)              &
                 - (1.0_dp - xi) * eta * xg(1,Ny)              &
                 - xi * eta * xg(Nx,Ny)

        yg(i,j) = (1.0_dp - eta) * yg(i,1) + eta * yg(i,Ny) &
                 + (1.0_dp - xi)  * yg(1,j) + xi  * yg(Nx,j) &
                 - (1.0_dp - xi) * (1.0_dp - eta) * yg(1,1)   &
                 - xi * (1.0_dp - eta) * yg(Nx,1)              &
                 - (1.0_dp - xi) * eta * yg(1,Ny)              &
                 - xi * eta * yg(Nx,Ny)
      end do
    end do

    write(*,*) '  TFI initialization done.'

  end subroutine tfi_init

  !--------------------------------------------------------------------
  ! Compute Q source terms (Thomas-Middlecoff method)
  ! Preserves wall-normal clustering during Poisson smoothing
  !--------------------------------------------------------------------
  subroutine compute_source_terms()
    integer :: i, j
    real(dp) :: x_eta, y_eta, x_ee, y_ee
    real(dp) :: Qbot_i, Qtop_i, eta, phi_b, phi_t
    real(dp), parameter :: decay_rate = 5.0_dp

    Qsrc = 0.0_dp

    do i = 2, Nx - 1
      ! Bottom boundary Q (j=1): one-sided 2nd-order differences
      x_eta = (-3.0_dp*xg(i,1) + 4.0_dp*xg(i,2) - xg(i,3)) * 0.5_dp
      y_eta = (-3.0_dp*yg(i,1) + 4.0_dp*yg(i,2) - yg(i,3)) * 0.5_dp
      x_ee  = xg(i,1) - 2.0_dp*xg(i,2) + xg(i,3)
      y_ee  = yg(i,1) - 2.0_dp*yg(i,2) + yg(i,3)

      if (x_eta**2 + y_eta**2 > 1.0e-20_dp) then
        Qbot_i = -(x_ee*x_eta + y_ee*y_eta) / (x_eta**2 + y_eta**2)
      else
        Qbot_i = 0.0_dp
      end if

      ! Top boundary Q (j=Ny)
      x_eta = (3.0_dp*xg(i,Ny) - 4.0_dp*xg(i,Ny-1) + xg(i,Ny-2)) * 0.5_dp
      y_eta = (3.0_dp*yg(i,Ny) - 4.0_dp*yg(i,Ny-1) + yg(i,Ny-2)) * 0.5_dp
      x_ee  = xg(i,Ny) - 2.0_dp*xg(i,Ny-1) + xg(i,Ny-2)
      y_ee  = yg(i,Ny) - 2.0_dp*yg(i,Ny-1) + yg(i,Ny-2)

      if (x_eta**2 + y_eta**2 > 1.0e-20_dp) then
        Qtop_i = -(x_ee*x_eta + y_ee*y_eta) / (x_eta**2 + y_eta**2)
      else
        Qtop_i = 0.0_dp
      end if

      ! Exponential decay interpolation to interior
      do j = 1, Ny
        eta  = real(j - 1, dp) / real(Ny - 1, dp)
        phi_b = exp(-decay_rate * eta)
        phi_t = exp(-decay_rate * (1.0_dp - eta))
        Qsrc(i,j) = Qbot_i * phi_b + Qtop_i * phi_t
      end do
    end do

    write(*,*) '  Source terms computed.'

  end subroutine compute_source_terms

  !--------------------------------------------------------------------
  ! TTM Poisson solver with SOR (Successive Over-Relaxation)
  !
  ! Solves: alpha*x_xixi - 2*beta*x_xieta + gamma*x_etaeta
  !         + J^2*(P*x_xi + Q*x_eta) = 0
  ! and similarly for y.
  !--------------------------------------------------------------------
  subroutine poisson_solve()
    integer :: i, j, iter
    real(dp) :: x_xi, y_xi, x_eta, y_eta
    real(dp) :: ac, bc, gc, Jac2
    real(dp) :: rhs_x, rhs_y, diag
    real(dp) :: xold, yold, max_res, res

    write(*,*) ''
    write(*,*) '  Starting Poisson (TTM) iteration...'
    write(*,'(A,F6.3,A,ES10.3)') '    omega=', omega_sor, '  tol=', tol_conv

    do iter = 1, max_iter
      max_res = 0.0_dp

      do j = 2, Ny - 1
        do i = 2, Nx - 1
          ! Central differences for first derivatives
          x_xi  = 0.5_dp * (xg(i+1,j) - xg(i-1,j))
          y_xi  = 0.5_dp * (yg(i+1,j) - yg(i-1,j))
          x_eta = 0.5_dp * (xg(i,j+1) - xg(i,j-1))
          y_eta = 0.5_dp * (yg(i,j+1) - yg(i,j-1))

          ! Metric coefficients
          ac = x_eta**2 + y_eta**2       ! alpha
          bc = x_xi*x_eta + y_xi*y_eta   ! beta
          gc = x_xi**2 + y_xi**2         ! gamma

          ! Jacobian squared
          Jac2 = (x_xi*y_eta - x_eta*y_xi)**2

          ! Diagonal coefficient
          diag = 2.0_dp * (ac + gc)
          if (abs(diag) < 1.0e-20_dp) cycle

          ! RHS for x-equation
          rhs_x = ac * (xg(i+1,j) + xg(i-1,j)) &
                + gc * (xg(i,j+1) + xg(i,j-1)) &
                - 0.5_dp * bc * (xg(i+1,j+1) - xg(i-1,j+1) &
                                - xg(i+1,j-1) + xg(i-1,j-1))
          ! Source term: + J^2 * Q * x_eta
          rhs_x = rhs_x + Jac2 * Qsrc(i,j) * x_eta

          ! RHS for y-equation
          rhs_y = ac * (yg(i+1,j) + yg(i-1,j)) &
                + gc * (yg(i,j+1) + yg(i,j-1)) &
                - 0.5_dp * bc * (yg(i+1,j+1) - yg(i-1,j+1) &
                                - yg(i+1,j-1) + yg(i-1,j-1))
          ! Source term: + J^2 * Q * y_eta
          rhs_y = rhs_y + Jac2 * Qsrc(i,j) * y_eta

          ! SOR update
          xold = xg(i,j)
          yold = yg(i,j)

          xg(i,j) = (1.0_dp - omega_sor) * xold &
                   + omega_sor * rhs_x / diag
          yg(i,j) = (1.0_dp - omega_sor) * yold &
                   + omega_sor * rhs_y / diag

          res = max(abs(xg(i,j) - xold), abs(yg(i,j) - yold))
          max_res = max(max_res, res)
        end do
      end do

      if (mod(iter, 1000) == 0 .or. iter == 1) then
        write(*,'(A,I6,A,ES12.5)') '    Iter ', iter, &
              '  max_res = ', max_res
      end if

      if (max_res < tol_conv) then
        write(*,'(A,I6,A)') '    Converged after ', iter, ' iterations.'
        return
      end if
    end do

    write(*,*) '    WARNING: Did NOT converge!'
    write(*,'(A,ES12.5)') '    Final residual = ', max_res

  end subroutine poisson_solve

  !--------------------------------------------------------------------
  ! Check grid quality (Jacobian sign)
  !--------------------------------------------------------------------
  subroutine check_jacobian()
    integer :: i, j, neg_count
    real(dp) :: Jac, Jmin, Jmax_v
    real(dp) :: dx_xi, dy_xi, dx_eta, dy_eta

    neg_count = 0
    Jmin = huge(1.0_dp)
    Jmax_v = -huge(1.0_dp)

    do j = 1, Ny - 1
      do i = 1, Nx - 1
        dx_xi  = xg(i+1,j) - xg(i,j)
        dy_xi  = yg(i+1,j) - yg(i,j)
        dx_eta = xg(i,j+1) - xg(i,j)
        dy_eta = yg(i,j+1) - yg(i,j)
        Jac = dx_xi * dy_eta - dx_eta * dy_xi

        Jmin = min(Jmin, Jac)
        Jmax_v = max(Jmax_v, Jac)
        if (Jac <= 0.0_dp) neg_count = neg_count + 1
      end do
    end do

    write(*,*) ''
    write(*,*) '  Grid quality:'
    write(*,'(A,ES12.5)') '    Min Jacobian = ', Jmin
    write(*,'(A,ES12.5)') '    Max Jacobian = ', Jmax_v
    write(*,'(A,I6)')     '    Negative Jacobian cells = ', neg_count

    if (neg_count > 0) then
      write(*,*) '    WARNING: Grid has crossings!'
    else
      write(*,*) '    OK: All Jacobians positive.'
    end if

  end subroutine check_jacobian

  !--------------------------------------------------------------------
  ! Write Plot3D formatted grid file
  !--------------------------------------------------------------------
  subroutine write_plot3d()
    integer :: i, j

    open(unit=10, file='periodic_hill.xyz', status='replace', &
         form='formatted')
    write(10, *) 1         ! one block
    write(10, *) Nx, Ny    ! dimensions
    ! Write x coordinates
    write(10, '(5ES20.12)') ((xg(i,j), i=1,Nx), j=1,Ny)
    ! Write y coordinates
    write(10, '(5ES20.12)') ((yg(i,j), i=1,Nx), j=1,Ny)
    close(10)

  end subroutine write_plot3d

  !--------------------------------------------------------------------
  ! Write simple data file for Python visualization
  !--------------------------------------------------------------------
  subroutine write_grid_dat()
    integer :: i, j

    open(unit=11, file='grid_data.dat', status='replace')
    write(11, '(2I6)') Nx, Ny
    do j = 1, Ny
      do i = 1, Nx
        write(11, '(2ES22.14)') xg(i,j), yg(i,j)
      end do
    end do
    close(11)

  end subroutine write_grid_dat

end program periodic_hill_grid
