#!/usr/bin/env python3
"""
Periodic Hill Structured Grid Generator
Using TTM (Thompson-Thames-Mastin) Poisson Equation Solver

Based on the ERCOFTAC UFR 3-30 Periodic Hill geometry
Hill profile from Almeida et al. (1993), Mellen et al. (2000)

Polynomial variable: physical x coordinate (NOT x/h)
Breakpoints: x = 0, 9, 14, 20, 30, 40, 54
Symmetry about x = Lx/2 = 126
"""

import numpy as np
import struct
import time

# ======== Geometry parameters ========
h_hill = 28.0
Lx = 9.0 * h_hill        # 252
Ly = 3.036 * h_hill       # 85.008
Lx_half = 0.5 * Lx        # 126

# ======== Grid dimensions ========
Nx = 160
Ny = 80

# ======== Solver parameters ========
omega_sor = 1.4
max_iter = 20000
tol_conv = 1.0e-6
beta_wall = 2.5   # tanh stretching param


def hill_profile(xp):
    """Hill profile function: returns y given physical x in [0, Lx].
    Uses correct ERCOFTAC polynomial coefficients (physical x variable).
    Vectorized version."""
    xp = np.atleast_1d(np.asarray(xp, dtype=np.float64))
    yh = np.zeros_like(xp)

    for idx in range(len(xp)):
        xx = xp[idx]
        # Apply symmetry about Lx/2
        if xx > Lx_half:
            xx = Lx - xx

        if xx < 0.0:
            yh[idx] = h_hill
        elif xx < 9.0:
            yh[idx] = min(28.0, 2.800000000000e+01
                          + 0.000000000000e+00 * xx
                          + 6.775070969851e-03 * xx**2
                          - 2.124527775800e-03 * xx**3)
        elif xx < 14.0:
            yh[idx] = (2.507355893131e+01
                       + 9.754803562315e-01 * xx
                       - 1.016116352781e-01 * xx**2
                       + 1.889794677828e-03 * xx**3)
        elif xx < 20.0:
            yh[idx] = (2.579601052357e+01
                       + 8.206693007457e-01 * xx
                       - 9.055370274339e-02 * xx**2
                       + 1.626510569859e-03 * xx**3)
        elif xx < 30.0:
            yh[idx] = (4.046435022819e+01
                       - 1.379581654948e+00 * xx
                       + 1.945884504128e-02 * xx**2
                       - 2.070318932190e-04 * xx**3)
        elif xx < 40.0:
            yh[idx] = (1.792461334664e+01
                       + 8.743920332081e-01 * xx
                       - 5.567361123058e-02 * xx**2
                       + 6.277731764683e-04 * xx**3)
        elif xx < 54.0:
            yh[idx] = max(0.0, 5.639011190988e+01
                          - 2.010520359035e+00 * xx
                          + 1.644919857549e-02 * xx**2
                          + 2.674976141766e-05 * xx**3)
        else:
            yh[idx] = 0.0

    return yh if len(yh) > 1 else yh[0]


def hill_profile_scalar(xp):
    """Scalar version for inner loops."""
    xx = xp
    if xx > Lx_half:
        xx = Lx - xx

    if xx < 0.0:
        return h_hill
    elif xx < 9.0:
        return min(28.0, 2.800000000000e+01
                   + 0.000000000000e+00 * xx
                   + 6.775070969851e-03 * xx**2
                   - 2.124527775800e-03 * xx**3)
    elif xx < 14.0:
        return (2.507355893131e+01
                + 9.754803562315e-01 * xx
                - 1.016116352781e-01 * xx**2
                + 1.889794677828e-03 * xx**3)
    elif xx < 20.0:
        return (2.579601052357e+01
                + 8.206693007457e-01 * xx
                - 9.055370274339e-02 * xx**2
                + 1.626510569859e-03 * xx**3)
    elif xx < 30.0:
        return (4.046435022819e+01
                - 1.379581654948e+00 * xx
                + 1.945884504128e-02 * xx**2
                - 2.070318932190e-04 * xx**3)
    elif xx < 40.0:
        return (1.792461334664e+01
                + 8.743920332081e-01 * xx
                - 5.567361123058e-02 * xx**2
                + 6.277731764683e-04 * xx**3)
    elif xx < 54.0:
        return max(0.0, 5.639011190988e+01
                   - 2.010520359035e+00 * xx
                   + 1.644919857549e-02 * xx**2
                   + 2.674976141766e-05 * xx**3)
    else:
        return 0.0


def tanh_stretch_wall(eta):
    """Tanh stretching: clusters near eta=0 (bottom wall)."""
    return 1.0 - np.tanh(beta_wall * (1.0 - eta)) / np.tanh(beta_wall)


def compute_bottom_wall():
    """Compute bottom wall points using arc-length parameterization."""
    Nfine = 20000
    dx_fine = Lx / (Nfine - 1)

    xf = np.linspace(0, Lx, Nfine)
    yf = hill_profile(xf)

    # Arc length
    sf = np.zeros(Nfine)
    for k in range(1, Nfine):
        dy = yf[k] - yf[k - 1]
        sf[k] = sf[k - 1] + np.sqrt(dx_fine**2 + dy**2)

    s_total = sf[-1]
    print(f"  Bottom wall arc length = {s_total:.3f}")

    # Distribute Nx points uniformly in arc length
    x_bot = np.zeros(Nx)
    y_bot = np.zeros(Nx)
    x_bot[0] = 0.0
    x_bot[-1] = Lx

    for i in range(1, Nx - 1):
        s_target = i / (Nx - 1) * s_total
        k = np.searchsorted(sf, s_target)
        if k >= Nfine:
            k = Nfine - 1
        if k < 1:
            k = 1
        x_bot[i] = xf[k - 1] + (s_target - sf[k - 1]) / (sf[k] - sf[k - 1]) * dx_fine

    for i in range(Nx):
        y_bot[i] = hill_profile_scalar(x_bot[i])

    # Top wall: uniform x distribution
    x_top = np.linspace(0, Lx, Nx)

    print(f"  x_bot range: [{x_bot[0]:.4f}, {x_bot[-1]:.4f}]")
    return x_bot, y_bot, x_top


def setup_boundaries(xg, yg, x_bot, y_bot, x_top):
    """Set up all four boundary conditions."""
    ybot_0 = hill_profile_scalar(0.0)  # = h_hill = 28

    # Bottom boundary (j=0): hill profile
    xg[:, 0] = x_bot
    yg[:, 0] = y_bot

    # Top boundary (j=Ny-1): flat at y = Ly
    xg[:, -1] = x_top
    yg[:, -1] = Ly

    # Left boundary (i=0): vertical from (0, h) to (0, Ly)
    for j in range(Ny):
        eta = j / (Ny - 1)
        s = tanh_stretch_wall(eta)
        xg[0, j] = 0.0
        yg[0, j] = ybot_0 + s * (Ly - ybot_0)
    yg[0, 0] = ybot_0
    yg[0, -1] = Ly

    # Right boundary (i=Nx-1): periodic (same y as left, x = Lx)
    xg[-1, :] = Lx
    yg[-1, :] = yg[0, :]

    print("  Boundaries initialized.")


def tfi_init(xg, yg):
    """Transfinite Interpolation (TFI) for initial grid."""
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            xi = i / (Nx - 1)
            eta = j / (Ny - 1)

            xg[i, j] = ((1 - eta) * xg[i, 0] + eta * xg[i, -1]
                         + (1 - xi) * xg[0, j] + xi * xg[-1, j]
                         - (1 - xi) * (1 - eta) * xg[0, 0]
                         - xi * (1 - eta) * xg[-1, 0]
                         - (1 - xi) * eta * xg[0, -1]
                         - xi * eta * xg[-1, -1])

            yg[i, j] = ((1 - eta) * yg[i, 0] + eta * yg[i, -1]
                         + (1 - xi) * yg[0, j] + xi * yg[-1, j]
                         - (1 - xi) * (1 - eta) * yg[0, 0]
                         - xi * (1 - eta) * yg[-1, 0]
                         - (1 - xi) * eta * yg[0, -1]
                         - xi * eta * yg[-1, -1])

    print("  TFI initialization done.")


def compute_source_terms(xg, yg):
    """Thomas-Middlecoff Q source terms with exponential decay."""
    Qsrc = np.zeros((Nx, Ny))
    decay_rate = 5.0

    for i in range(1, Nx - 1):
        # Bottom boundary Q (j=0): one-sided 2nd-order differences
        x_eta = (-3.0 * xg[i, 0] + 4.0 * xg[i, 1] - xg[i, 2]) * 0.5
        y_eta = (-3.0 * yg[i, 0] + 4.0 * yg[i, 1] - yg[i, 2]) * 0.5
        x_ee = xg[i, 0] - 2.0 * xg[i, 1] + xg[i, 2]
        y_ee = yg[i, 0] - 2.0 * yg[i, 1] + yg[i, 2]

        denom = x_eta**2 + y_eta**2
        if denom > 1e-20:
            Qbot_i = -(x_ee * x_eta + y_ee * y_eta) / denom
        else:
            Qbot_i = 0.0

        # Top boundary Q (j=Ny-1)
        x_eta = (3.0 * xg[i, -1] - 4.0 * xg[i, -2] + xg[i, -3]) * 0.5
        y_eta = (3.0 * yg[i, -1] - 4.0 * yg[i, -2] + yg[i, -3]) * 0.5
        x_ee = xg[i, -1] - 2.0 * xg[i, -2] + xg[i, -3]
        y_ee = yg[i, -1] - 2.0 * yg[i, -2] + yg[i, -3]

        denom = x_eta**2 + y_eta**2
        if denom > 1e-20:
            Qtop_i = -(x_ee * x_eta + y_ee * y_eta) / denom
        else:
            Qtop_i = 0.0

        # Exponential decay interpolation
        for j in range(Ny):
            eta = j / (Ny - 1)
            phi_b = np.exp(-decay_rate * eta)
            phi_t = np.exp(-decay_rate * (1.0 - eta))
            Qsrc[i, j] = Qbot_i * phi_b + Qtop_i * phi_t

    print("  Source terms computed.")
    return Qsrc


def poisson_solve(xg, yg, Qsrc):
    """TTM Poisson solver with SOR."""
    print()
    print("  Starting Poisson (TTM) iteration...")
    print(f"    omega={omega_sor:.3f}  tol={tol_conv:.1e}")

    t0 = time.time()

    for iteration in range(1, max_iter + 1):
        max_res = 0.0

        for j in range(1, Ny - 1):
            for i in range(1, Nx - 1):
                # Central differences for first derivatives
                x_xi = 0.5 * (xg[i + 1, j] - xg[i - 1, j])
                y_xi = 0.5 * (yg[i + 1, j] - yg[i - 1, j])
                x_eta = 0.5 * (xg[i, j + 1] - xg[i, j - 1])
                y_eta = 0.5 * (yg[i, j + 1] - yg[i, j - 1])

                # Metric coefficients
                ac = x_eta**2 + y_eta**2      # alpha
                bc = x_xi * x_eta + y_xi * y_eta  # beta
                gc = x_xi**2 + y_xi**2        # gamma

                # Jacobian squared
                Jac2 = (x_xi * y_eta - x_eta * y_xi)**2

                # Diagonal coefficient
                diag = 2.0 * (ac + gc)
                if abs(diag) < 1e-20:
                    continue

                # RHS for x-equation
                rhs_x = (ac * (xg[i + 1, j] + xg[i - 1, j])
                          + gc * (xg[i, j + 1] + xg[i, j - 1])
                          - 0.5 * bc * (xg[i + 1, j + 1] - xg[i - 1, j + 1]
                                        - xg[i + 1, j - 1] + xg[i - 1, j - 1]))
                rhs_x += Jac2 * Qsrc[i, j] * x_eta

                # RHS for y-equation
                rhs_y = (ac * (yg[i + 1, j] + yg[i - 1, j])
                          + gc * (yg[i, j + 1] + yg[i, j - 1])
                          - 0.5 * bc * (yg[i + 1, j + 1] - yg[i - 1, j + 1]
                                        - yg[i + 1, j - 1] + yg[i - 1, j - 1]))
                rhs_y += Jac2 * Qsrc[i, j] * y_eta

                # SOR update
                xold = xg[i, j]
                yold = yg[i, j]

                xg[i, j] = (1.0 - omega_sor) * xold + omega_sor * rhs_x / diag
                yg[i, j] = (1.0 - omega_sor) * yold + omega_sor * rhs_y / diag

                res = max(abs(xg[i, j] - xold), abs(yg[i, j] - yold))
                max_res = max(max_res, res)

        if iteration % 500 == 0 or iteration == 1:
            elapsed = time.time() - t0
            print(f"    Iter {iteration:6d}  max_res = {max_res:.5e}  ({elapsed:.1f}s)")

        if max_res < tol_conv:
            elapsed = time.time() - t0
            print(f"    Converged after {iteration} iterations. ({elapsed:.1f}s)")
            return

    elapsed = time.time() - t0
    print(f"    WARNING: Did NOT converge! Final residual = {max_res:.5e} ({elapsed:.1f}s)")


def check_jacobian(xg, yg):
    """Check grid quality (Jacobian sign)."""
    neg_count = 0
    Jmin = np.inf
    Jmax_v = -np.inf

    for j in range(Ny - 1):
        for i in range(Nx - 1):
            dx_xi = xg[i + 1, j] - xg[i, j]
            dy_xi = yg[i + 1, j] - yg[i, j]
            dx_eta = xg[i, j + 1] - xg[i, j]
            dy_eta = yg[i, j + 1] - yg[i, j]
            Jac = dx_xi * dy_eta - dx_eta * dy_xi

            Jmin = min(Jmin, Jac)
            Jmax_v = max(Jmax_v, Jac)
            if Jac <= 0.0:
                neg_count += 1

    print()
    print("  Grid quality:")
    print(f"    Min Jacobian = {Jmin:.5e}")
    print(f"    Max Jacobian = {Jmax_v:.5e}")
    print(f"    Negative Jacobian cells = {neg_count}")

    if neg_count > 0:
        print("    WARNING: Grid has crossings!")
    else:
        print("    OK: All Jacobians positive.")

    return neg_count


def write_plot3d(xg, yg, filename):
    """Write Plot3D formatted grid file (ASCII)."""
    with open(filename, 'w') as f:
        f.write("1\n")             # one block
        f.write(f"{Nx} {Ny}\n")    # dimensions
        # Write x coordinates
        vals = []
        for j in range(Ny):
            for i in range(Nx):
                vals.append(f"{xg[i, j]:.12e}")
                if len(vals) == 5:
                    f.write(" ".join(vals) + "\n")
                    vals = []
        if vals:
            f.write(" ".join(vals) + "\n")
            vals = []
        # Write y coordinates
        for j in range(Ny):
            for i in range(Nx):
                vals.append(f"{yg[i, j]:.12e}")
                if len(vals) == 5:
                    f.write(" ".join(vals) + "\n")
                    vals = []
        if vals:
            f.write(" ".join(vals) + "\n")

    print(f"  Plot3D file written: {filename}")


def write_grid_dat(xg, yg, filename):
    """Write simple data file for visualization."""
    with open(filename, 'w') as f:
        f.write(f"{Nx:6d}{Ny:6d}\n")
        for j in range(Ny):
            for i in range(Nx):
                f.write(f"{xg[i, j]:22.14e}{yg[i, j]:22.14e}\n")
    print(f"  Data file written: {filename}")


# ======== Main Execution ========
if __name__ == "__main__":
    print("=================================================")
    print(" Periodic Hill Structured Grid Generator")
    print(" TTM Poisson Equation Solver")
    print("=================================================")
    print(f"  Grid: Nx={Nx}  Ny={Ny}")
    print(f"  Domain: [0, {Lx:.2f}] x [0, {Ly:.2f}]")
    print(f"  Hill height h = {h_hill:.2f}")
    print()

    # Step 1: Bottom wall (arc-length parameterization)
    x_bot, y_bot, x_top = compute_bottom_wall()

    # Step 2: Set up grid arrays and boundaries
    xg = np.zeros((Nx, Ny))
    yg = np.zeros((Nx, Ny))
    setup_boundaries(xg, yg, x_bot, y_bot, x_top)

    # Step 3: TFI initialization
    tfi_init(xg, yg)

    # Step 4: Source terms
    Qsrc = compute_source_terms(xg, yg)

    # Step 5: Poisson solver
    poisson_solve(xg, yg, Qsrc)

    # Step 6: Check quality
    neg_count = check_jacobian(xg, yg)

    # Step 7: Output
    write_plot3d(xg, yg, "periodic_hill.xyz")
    write_grid_dat(xg, yg, "grid_data.dat")

    # Save numpy arrays for visualization
    np.savez("periodic_hill_grid.npz", xg=xg, yg=yg, Nx=Nx, Ny=Ny)

    print()
    print("Grid generation complete!")
    print("  Plot3D file: periodic_hill.xyz")
    print("  Data file:   grid_data.dat")
    print("  NumPy file:  periodic_hill_grid.npz")
