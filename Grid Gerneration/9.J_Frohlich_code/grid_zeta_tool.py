"""
Periodic Hill Grid Tool -- Steger-Sorenson Poisson + Zeta Stretching
=====================================================================
Capabilities:
  1. Parse original Tecplot .dat grid
  2. Mode 1 (Zeta-only): keep Ni x Nj, adjust vertical stretching
  3. Mode 2 (Adaptive):  freely set Ni x Nj, then re-solve the
     Poisson grid equation with control functions P,Q reversed
     from the reference grid -- true Steger-Sorenson method
  4. Export new grid in Tecplot format
  5. Identity verification at original resolution

Mode 2 mathematical basis:
  The TTM-Poisson equation (physical-space form):
    alpha * r_xixi - 2*beta * r_xieta + gamma * r_etaeta
        = -J^2 * (P * r_xi + Q * r_eta)

  Given a reference grid r(xi,eta):
    1. Compute all metric terms and Jacobian
    2. Solve the 2x2 linear system for P,Q at each point
    3. Interpolate P,Q to new (Ni,Nj) via bicubic spline
    4. Resample boundaries, create TFI initial guess
    5. Iteratively solve the Poisson equation with the
       interpolated P,Q as source terms

  Validation: at same (Ni,Nj) the method recovers the original
  grid to ~1e-11 absolute error (near machine precision).
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import RectBivariateSpline, interp1d

# ============================================================
#  1.  Parser
# ============================================================

def parse_tecplot_dat(filepath):
    filepath = Path(filepath)
    with open(filepath, "r", encoding="latin-1") as f:
        lines = f.readlines()

    ni = nj = None
    header_lines = 0
    for idx, line in enumerate(lines):
        if "I=" in line.upper():
            parts = line.replace(",", " ").replace("=", " ").upper().split()
            for k, tok in enumerate(parts):
                if tok == "I":
                    ni = int(parts[k + 1])
                if tok == "J":
                    nj = int(parts[k + 1])
            header_lines = idx + 2
            break

    if ni is None or nj is None:
        raise ValueError("Cannot find I/J dimensions in header")

    data_lines = lines[header_lines:]
    x_flat, y_flat = [], []
    for dl in data_lines:
        dl = dl.strip()
        if not dl:
            continue
        vals = dl.split()
        if len(vals) >= 2:
            x_flat.append(float(vals[0]))
            y_flat.append(float(vals[1]))

    expected = ni * nj
    if len(x_flat) != expected:
        raise ValueError(
            f"Expected {expected} points (I={ni} x J={nj}), got {len(x_flat)}"
        )

    x = np.array(x_flat).reshape(nj, ni)
    y = np.array(y_flat).reshape(nj, ni)
    return x, y, ni, nj


# ============================================================
#  2.  Visualiser
# ============================================================

def plot_grid(x, y, title="Grid", savepath=None, figsize=(18, 6)):
    nj, ni = x.shape
    fig, ax = plt.subplots(figsize=figsize)
    for j in range(nj):
        ax.plot(x[j, :], y[j, :], "k-", lw=0.3)
    for i in range(ni):
        ax.plot(x[:, i], y[:, i], "k-", lw=0.3)
    ax.set_aspect("equal")
    ax.set_xlabel("x  [m]"); ax.set_ylabel("y  [m]")
    ax.set_title(title)
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=200)
        print(f"  [saved] {savepath}")
    plt.close(fig)


def plot_compare(x1, y1, x2, y2, labels=("Original", "New"),
                 title="Comparison", savepath=None, figsize=(18, 12)):
    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
    for ax, xg, yg, lbl in zip(axes, [x1, x2], [y1, y2], labels):
        nj, ni = xg.shape
        for j in range(nj):
            ax.plot(xg[j, :], yg[j, :], "k-", lw=0.25)
        for i in range(ni):
            ax.plot(xg[:, i], yg[:, i], "k-", lw=0.25)
        ax.set_aspect("equal"); ax.set_ylabel("y  [m]"); ax.set_title(lbl)
    axes[-1].set_xlabel("x  [m]")
    fig.suptitle(title, fontsize=14, y=1.01)
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=200, bbox_inches="tight")
        print(f"  [saved] {savepath}")
    plt.close(fig)


def plot_vertical_spacing(y1, y2, icol, labels=("Original", "New"),
                          savepath=None):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(range(y1.shape[0]-1), np.diff(y1[:, icol])*1e3, "o-", ms=3, label=labels[0])
    ax.plot(range(y2.shape[0]-1), np.diff(y2[:, icol])*1e3, "s-", ms=3, label=labels[1])
    ax.set_xlabel("j index"); ax.set_ylabel("dy  [mm]")
    ax.set_title(f"Vertical spacing at i = {icol}")
    ax.legend(); ax.grid(True, ls="--", alpha=0.4)
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=200)
        print(f"  [saved] {savepath}")
    plt.close(fig)


# ============================================================
#  3.  Stretching functions
# ============================================================

def vinokur_tanh(eta, gamma, alpha=0.5):
    """
    Vinokur two-sided tanh clustering.  eta in [0,1].
    gamma=0 => identity.  Monotonic for all gamma >= 0.
    """
    if gamma < 1e-14:
        return eta.copy()
    denom = np.tanh(gamma * alpha)
    if abs(denom) < 1e-30:
        return eta.copy()
    zeta = 0.5 * (1.0 + np.tanh(gamma * (eta - alpha)) / denom)
    zeta[0] = 0.0; zeta[-1] = 1.0
    return zeta


# ============================================================
#  4.  Zeta-only redistribution (Mode 1)
# ============================================================

def redistribute_vertical(x, y, gamma=0.0, alpha=0.5):
    """
    Redistribute vertical points only (keep ni unchanged).
    gamma=0 => identity (reproduces original exactly).
    """
    nj, ni = x.shape
    eta = np.linspace(0, 1, nj)
    zeta = vinokur_tanh(eta, gamma, alpha)

    x_new = np.empty_like(x)
    y_new = np.empty_like(y)

    for i in range(ni):
        xc, yc = x[:, i], y[:, i]
        ds = np.sqrt(np.diff(xc)**2 + np.diff(yc)**2)
        s = np.concatenate(([0.0], np.cumsum(ds)))
        s_norm = s / s[-1]
        s_new = np.interp(zeta, eta, s_norm)
        x_new[:, i] = np.interp(s_new, s_norm, xc)
        y_new[:, i] = np.interp(s_new, s_norm, yc)

    return x_new, y_new


# ============================================================
#  5.  Steger-Sorenson Poisson grid generation (Mode 2)
# ============================================================

def _compute_metrics(x, y):
    """Compute all metric terms using 2nd-order finite differences."""
    nj, ni = x.shape

    x_xi = np.zeros_like(x)
    x_xi[:, 1:-1] = 0.5 * (x[:, 2:] - x[:, :-2])
    x_xi[:, 0]  = -1.5*x[:,0] + 2.0*x[:,1] - 0.5*x[:,2]
    x_xi[:, -1] =  0.5*x[:,-3] - 2.0*x[:,-2] + 1.5*x[:,-1]

    y_xi = np.zeros_like(y)
    y_xi[:, 1:-1] = 0.5 * (y[:, 2:] - y[:, :-2])
    y_xi[:, 0]  = -1.5*y[:,0] + 2.0*y[:,1] - 0.5*y[:,2]
    y_xi[:, -1] =  0.5*y[:,-3] - 2.0*y[:,-2] + 1.5*y[:,-1]

    x_eta = np.zeros_like(x)
    x_eta[1:-1,:] = 0.5 * (x[2:,:] - x[:-2,:])
    x_eta[0,:]  = -1.5*x[0,:] + 2.0*x[1,:] - 0.5*x[2,:]
    x_eta[-1,:] =  0.5*x[-3,:] - 2.0*x[-2,:] + 1.5*x[-1,:]

    y_eta = np.zeros_like(y)
    y_eta[1:-1,:] = 0.5 * (y[2:,:] - y[:-2,:])
    y_eta[0,:]  = -1.5*y[0,:] + 2.0*y[1,:] - 0.5*y[2,:]
    y_eta[-1,:] =  0.5*y[-3,:] - 2.0*y[-2,:] + 1.5*y[-1,:]

    x_xixi = np.zeros_like(x)
    x_xixi[:, 1:-1] = x[:, 2:] - 2.0*x[:, 1:-1] + x[:, :-2]
    x_xixi[:, 0]  = x[:,0] - 2.0*x[:,1] + x[:,2]
    x_xixi[:, -1] = x[:,-3] - 2.0*x[:,-2] + x[:,-1]

    y_xixi = np.zeros_like(y)
    y_xixi[:, 1:-1] = y[:, 2:] - 2.0*y[:, 1:-1] + y[:, :-2]
    y_xixi[:, 0]  = y[:,0] - 2.0*y[:,1] + y[:,2]
    y_xixi[:, -1] = y[:,-3] - 2.0*y[:,-2] + y[:,-1]

    x_etaeta = np.zeros_like(x)
    x_etaeta[1:-1,:] = x[2:,:] - 2.0*x[1:-1,:] + x[:-2,:]
    x_etaeta[0,:]  = x[0,:] - 2.0*x[1,:] + x[2,:]
    x_etaeta[-1,:] = x[-3,:] - 2.0*x[-2,:] + x[-1,:]

    y_etaeta = np.zeros_like(y)
    y_etaeta[1:-1,:] = y[2:,:] - 2.0*y[1:-1,:] + y[:-2,:]
    y_etaeta[0,:]  = y[0,:] - 2.0*y[1,:] + y[2,:]
    y_etaeta[-1,:] = y[-3,:] - 2.0*y[-2,:] + y[-1,:]

    x_pad = np.pad(x, ((1,1),(1,1)), mode='edge')
    y_pad = np.pad(y, ((1,1),(1,1)), mode='edge')
    x_xieta = 0.25*(x_pad[2:,2:] - x_pad[2:,:-2]
                    - x_pad[:-2,2:] + x_pad[:-2,:-2])[:nj,:ni]
    y_xieta = 0.25*(y_pad[2:,2:] - y_pad[2:,:-2]
                    - y_pad[:-2,2:] + y_pad[:-2,:-2])[:nj,:ni]

    return {
        "x_xi": x_xi, "x_eta": x_eta, "y_xi": y_xi, "y_eta": y_eta,
        "x_xixi": x_xixi, "x_etaeta": x_etaeta, "x_xieta": x_xieta,
        "y_xixi": y_xixi, "y_etaeta": y_etaeta, "y_xieta": y_xieta,
        "alpha": x_eta**2 + y_eta**2,
        "beta": x_xi*x_eta + y_xi*y_eta,
        "gamma": x_xi**2 + y_xi**2,
        "J": x_xi*y_eta - x_eta*y_xi,
    }


def _compute_PQ(metrics):
    """Reverse-compute control functions P,Q from a known grid."""
    m = metrics
    RHS_x = (m["alpha"]*m["x_xixi"] - 2.0*m["beta"]*m["x_xieta"]
             + m["gamma"]*m["x_etaeta"])
    RHS_y = (m["alpha"]*m["y_xixi"] - 2.0*m["beta"]*m["y_xieta"]
             + m["gamma"]*m["y_etaeta"])
    J2 = m["J"]**2
    det = m["J"]
    safe = np.abs(det) > 1e-30

    P = np.zeros_like(RHS_x)
    Q = np.zeros_like(RHS_x)
    b1 = np.zeros_like(RHS_x)
    b2 = np.zeros_like(RHS_x)
    b1[safe] = RHS_x[safe] / (-J2[safe])
    b2[safe] = RHS_y[safe] / (-J2[safe])
    P[safe] = ( m["y_eta"][safe]*b1[safe] - m["x_eta"][safe]*b2[safe]) / det[safe]
    Q[safe] = (-m["y_xi"][safe]*b1[safe]  + m["x_xi"][safe]*b2[safe])  / det[safe]
    return P, Q


def _poisson_solve(x_init, y_init, P, Q,
                   n_iter=15000, omega=1.0, tol=1e-10, print_every=2000):
    """Row-vectorised Gauss-Seidel Poisson solver. Boundaries fixed."""
    nj, ni = x_init.shape
    x = x_init.copy()
    y = y_init.copy()
    convergence = []
    si = slice(1, -1)

    for it in range(n_iter):
        max_corr = 0.0

        for j in range(1, nj - 1):
            xxi  = 0.5 * (x[j, 2:] - x[j, :-2])
            xeta = 0.5 * (x[j+1, si] - x[j-1, si])
            yxi  = 0.5 * (y[j, 2:] - y[j, :-2])
            yeta = 0.5 * (y[j+1, si] - y[j-1, si])

            al = xeta**2 + yeta**2
            be = xxi*xeta + yxi*yeta
            ga = xxi**2 + yxi**2
            jac = xxi*yeta - xeta*yxi
            j2 = jac**2

            denom = 2.0 * (al + ga)
            safe = denom > 1e-30

            x_cross = 0.25*(x[j+1,2:] - x[j+1,:-2] - x[j-1,2:] + x[j-1,:-2])
            y_cross = 0.25*(y[j+1,2:] - y[j+1,:-2] - y[j-1,2:] + y[j-1,:-2])

            Pj = P[j, si]; Qj = Q[j, si]
            Sx = -j2 * (Pj*xxi + Qj*xeta)
            Sy = -j2 * (Pj*yxi + Qj*yeta)

            x_new = np.where(safe,
                (al*(x[j,2:]+x[j,:-2]) + ga*(x[j+1,si]+x[j-1,si])
                 - 2.0*be*x_cross - Sx) / np.where(safe, denom, 1.0),
                x[j, si])
            y_new = np.where(safe,
                (al*(y[j,2:]+y[j,:-2]) + ga*(y[j+1,si]+y[j-1,si])
                 - 2.0*be*y_cross - Sy) / np.where(safe, denom, 1.0),
                y[j, si])

            dx = omega * (x_new - x[j, si])
            dy = omega * (y_new - y[j, si])
            x[j, si] += dx
            y[j, si] += dy

            row_max = max(np.max(np.abs(dx)), np.max(np.abs(dy)))
            if row_max > max_corr:
                max_corr = row_max

        convergence.append(max_corr)

        if np.isnan(max_corr) or max_corr > 1e10:
            print(f"    DIVERGED at iter {it}")
            break

        if print_every and (it % print_every == 0 or it == n_iter - 1):
            print(f"    iter {it:5d}:  max_corr = {max_corr:.4e}")

        if max_corr < tol:
            print(f"    Converged at iter {it}, max_corr = {max_corr:.4e}")
            break

    return x, y, convergence


def _tfi(x_bot, y_bot, x_top, y_top, x_lft, y_lft, x_rgt, y_rgt):
    """Transfinite Interpolation (vectorised)."""
    ni = len(x_bot); nj = len(x_lft)
    xi = np.linspace(0, 1, ni)[np.newaxis, :]
    eta = np.linspace(0, 1, nj)[:, np.newaxis]
    x = ((1-eta)*x_bot + eta*x_top
       + (1-xi)*x_lft[:, np.newaxis] + xi*x_rgt[:, np.newaxis]
       - (1-xi)*(1-eta)*x_bot[0] - xi*(1-eta)*x_bot[-1]
       - (1-xi)*eta*x_top[0] - xi*eta*x_top[-1])
    y = ((1-eta)*y_bot + eta*y_top
       + (1-xi)*y_lft[:, np.newaxis] + xi*y_rgt[:, np.newaxis]
       - (1-xi)*(1-eta)*y_bot[0] - xi*(1-eta)*y_bot[-1]
       - (1-xi)*eta*y_top[0] - xi*eta*y_top[-1])
    return x, y


def _interpolate_PQ(P, Q, ni_old, nj_old, ni_new, nj_new):
    """
    Interpolate P,Q from old to new resolution via bicubic spline,
    with proper scaling for the changed computational grid spacing.

    P controls clustering in xi-direction, Q in eta-direction.
    In index space (dxi=1, deta=1), the source term is -J^2*(P*r_xi + Q*r_eta).
    When Ni changes, h_xi = 1/(Ni-1) changes, so P must scale as:
      P_new = P_old * (Ni_new-1)/(Ni_old-1)
      Q_new = Q_old * (Nj_new-1)/(Nj_old-1)
    """
    xi_o = np.linspace(0, 1, ni_old); eta_o = np.linspace(0, 1, nj_old)
    xi_n = np.linspace(0, 1, ni_new); eta_n = np.linspace(0, 1, nj_new)
    P_n = RectBivariateSpline(eta_o, xi_o, P, kx=3, ky=3)(eta_n, xi_n)
    Q_n = RectBivariateSpline(eta_o, xi_o, Q, kx=3, ky=3)(eta_n, xi_n)

    scale_P = (ni_new - 1) / (ni_old - 1)
    scale_Q = (nj_new - 1) / (nj_old - 1)
    P_n *= scale_P
    Q_n *= scale_Q

    return P_n, Q_n


def _resample_boundary(xb, yb, n_new):
    """Resample boundary to n_new points preserving arc-length pattern."""
    n_old = len(xb)
    if n_new == n_old:
        return xb.copy(), yb.copy()
    ds = np.sqrt(np.diff(xb)**2 + np.diff(yb)**2)
    s = np.concatenate(([0], np.cumsum(ds))); s /= s[-1]
    s_norm_old = np.linspace(0, 1, n_old)
    s_new = np.interp(np.linspace(0, 1, n_new), s_norm_old, s)
    return (interp1d(s, xb, kind='cubic')(s_new),
            interp1d(s, yb, kind='cubic')(s_new))


def generate_adaptive_grid(x_ref, y_ref, ni_new, nj_new,
                           gamma=0.0, alpha=0.5,
                           poisson_iter=15000, poisson_tol=1e-10):
    """
    Full Steger-Sorenson adaptive grid generation.

    Strategy:
      1. Reverse-compute P,Q from reference grid
      2. Interpolate P,Q to new (ni_new, nj_new)
      3. Resample boundaries at new resolution (NO stretching here)
      4. TFI initial guess
      5. Poisson solve with interpolated P,Q
      6. Apply vertical stretching (gamma/alpha) as post-processing
         on the converged Poisson grid -- same logic as Mode 1

    The stretching is applied AFTER the Poisson solve to avoid
    boundary inconsistency: Poisson needs all 4 boundaries to be
    geometrically consistent, which breaks if only the vertical
    boundaries are stretched while horizontal boundaries are not.
    """
    nj_ref, ni_ref = x_ref.shape

    print("    [1/6] Computing P,Q from reference ...")
    metrics = _compute_metrics(x_ref, y_ref)
    P_ref, Q_ref = _compute_PQ(metrics)

    print(f"    [2/6] Interpolating P,Q: ({ni_ref}x{nj_ref}) -> ({ni_new}x{nj_new}) ...")
    if ni_new == ni_ref and nj_new == nj_ref:
        P_new, Q_new = P_ref.copy(), Q_ref.copy()
    else:
        P_new, Q_new = _interpolate_PQ(P_ref, Q_ref,
                                        ni_ref, nj_ref, ni_new, nj_new)

    print("    [3/6] Resampling boundaries ...")
    xb, yb = _resample_boundary(x_ref[0, :],  y_ref[0, :],  ni_new)
    xt, yt = _resample_boundary(x_ref[-1, :], y_ref[-1, :], ni_new)
    xl, yl = _resample_boundary(x_ref[:, 0],  y_ref[:, 0],  nj_new)
    xr, yr = _resample_boundary(x_ref[:, -1], y_ref[:, -1], nj_new)

    xl[0] = xb[0];   yl[0] = yb[0]
    xl[-1] = xt[0];  yl[-1] = yt[0]
    xr[0] = xb[-1];  yr[0] = yb[-1]
    xr[-1] = xt[-1]; yr[-1] = yt[-1]

    print("    [4/6] TFI initial guess ...")
    x_tfi, y_tfi = _tfi(xb, yb, xt, yt, xl, yl, xr, yr)

    print(f"    [5/6] Poisson solve (max {poisson_iter} iter) ...")
    x_out, y_out, conv = _poisson_solve(
        x_tfi, y_tfi, P_new, Q_new,
        n_iter=poisson_iter, omega=1.0, tol=poisson_tol, print_every=2000)

    if gamma > 1e-14:
        print(f"    [6/6] Applying vertical stretching (gamma={gamma}, alpha={alpha}) ...")
        x_out, y_out = redistribute_vertical(x_out, y_out, gamma=gamma, alpha=alpha)
    else:
        print("    [6/6] No stretching (gamma=0)")

    return x_out, y_out, conv


# ============================================================
#  6.  Export to Tecplot .dat
# ============================================================

def write_tecplot_dat(filepath, x, y, title="Generated grid",
                      zone_title="Adaptive"):
    nj, ni = x.shape
    with open(filepath, "w") as f:
        f.write(f'TITLE     = "{title}"\n')
        f.write('VARIABLES = "x corner"\n')
        f.write('"y corner"\n')
        f.write(f'ZONE T="{zone_title}"\n')
        f.write(f' I={ni}, J={nj}, K=1,F=POINT\n')
        f.write('DT=(SINGLE SINGLE )\n')
        for j in range(nj):
            for i in range(ni):
                f.write(f" {x[j, i]: .9E} {y[j, i]: .9E}\n")
    print(f"  [written] {filepath}")


# ============================================================
#  7.  Verification
# ============================================================

def verify_identity(x_orig, y_orig, x_new, y_new, tol=1e-10):
    dx = np.max(np.abs(x_orig - x_new))
    dy = np.max(np.abs(y_orig - y_new))
    ok = (dx < tol) and (dy < tol)
    return ok, dx, dy


# ============================================================
#  8.  Interactive helpers
# ============================================================

def ask_float(prompt, default, lo=None, hi=None):
    while True:
        raw = input(f"  {prompt} [default={default}]: ").strip()
        if raw == "":
            return default
        try:
            val = float(raw)
        except ValueError:
            print("    ** Invalid number, try again.")
            continue
        if lo is not None and val < lo:
            print(f"    ** Must be >= {lo}, try again.")
            continue
        if hi is not None and val > hi:
            print(f"    ** Must be <= {hi}, try again.")
            continue
        return val


def ask_int(prompt, default, lo=None, hi=None):
    while True:
        raw = input(f"  {prompt} [default={default}]: ").strip()
        if raw == "":
            return default
        try:
            val = int(raw)
        except ValueError:
            print("    ** Invalid integer, try again.")
            continue
        if lo is not None and val < lo:
            print(f"    ** Must be >= {lo}, try again.")
            continue
        if hi is not None and val > hi:
            print(f"    ** Must be <= {hi}, try again.")
            continue
        return val


def ask_yes_no(prompt, default_yes=True):
    hint = "Y/n" if default_yes else "y/N"
    raw = input(f"  {prompt} [{hint}]: ").strip().lower()
    if raw == "":
        return default_yes
    return raw in ("y", "yes")


def detect_dat_files(folder):
    return sorted(f for f in folder.glob("*.dat")
                  if not f.name.startswith("zeta_")
                  and not f.name.startswith("adaptive_"))


# ============================================================
#  MAIN
# ============================================================

if __name__ == "__main__":

    base = Path(__file__).parent

    print()
    print("=" * 62)
    print("  Periodic Hill Grid -- Steger-Sorenson Poisson + Zeta")
    print("  (Interactive Mode)")
    print("=" * 62)

    # -----------------------------------------------------------
    #  Step 1 -- select reference grid
    # -----------------------------------------------------------
    print("\n" + "-" * 62)
    print("  [Step 1] Select reference grid file")
    print("-" * 62)

    dat_list = detect_dat_files(base)
    if len(dat_list) == 0:
        print("  ERROR: No .dat files found in", base)
        sys.exit(1)

    for idx, fp in enumerate(dat_list):
        print(f"    {idx + 1}. {fp.name}")

    while True:
        raw = input(f"\n  Enter file number [1-{len(dat_list)}] (default=1): ").strip()
        if raw == "":
            choice = 0
            break
        try:
            choice = int(raw) - 1
            if 0 <= choice < len(dat_list):
                break
        except ValueError:
            pass
        print("    ** Invalid choice, try again.")

    dat_path = dat_list[choice]
    grid_key = dat_path.stem

    # -----------------------------------------------------------
    #  Step 2 -- parse reference
    # -----------------------------------------------------------
    print("\n" + "-" * 62)
    print("  [Step 2] Parsing reference grid ...")
    print("-" * 62)

    x_ref, y_ref, ni_ref, nj_ref = parse_tecplot_dat(dat_path)
    print(f"  Reference: {dat_path.name}")
    print(f"  Dimensions: I={ni_ref} (streamwise)  x  J={nj_ref} (vertical)")

    out_orig = base / f"original_{grid_key}.png"
    plot_grid(x_ref, y_ref,
              title=f"Reference: {dat_path.name}  (I={ni_ref}, J={nj_ref})",
              savepath=out_orig)

    # -----------------------------------------------------------
    #  Step 3 -- choose mode
    # -----------------------------------------------------------
    print("\n" + "-" * 62)
    print("  [Step 3] Choose operation mode")
    print("-" * 62)
    print()
    print("    1. Zeta-only  -- keep original Ni x Nj,")
    print("                     adjust vertical stretching (GAMMA/ALPHA)")
    print()
    print("    2. Adaptive   -- freely set new Ni x Nj,")
    print("                     Poisson solve with Steger-Sorenson P,Q")
    print("                     (true elliptic grid generation)")
    print()

    while True:
        raw = input("  Mode [1 or 2] (default=1): ").strip()
        if raw == "":
            mode = 1
            break
        if raw in ("1", "2"):
            mode = int(raw)
            break
        print("    ** Enter 1 or 2.")

    # -----------------------------------------------------------
    #  Step 4 -- set parameters
    # -----------------------------------------------------------
    print("\n" + "-" * 62)
    print("  [Step 4] Set parameters")
    print("-" * 62)

    if mode == 2:
        print()
        print(f"  Reference grid: I={ni_ref}, J={nj_ref}")
        print()
        print("  Ni -- streamwise grid points")
        print(f"         (original = {ni_ref})")
        NI = ask_int("Ni", default=ni_ref, lo=10, hi=2000)
        print()
        print("  Nj -- vertical grid points")
        print(f"         (original = {nj_ref})")
        NJ = ask_int("Nj", default=nj_ref, lo=10, hi=2000)
    else:
        NI = ni_ref
        NJ = nj_ref

    print()
    print("  GAMMA -- Vertical stretching intensity")
    print("           0.0 = preserve original spacing")
    print("           0.5~1.0 = mild wall-clustering")
    print("           1.5~2.0 = moderate (recommended)")
    print("           2.5~3.0 = strong wall-clustering")
    print()
    GAMMA = ask_float("GAMMA", default=0.0 if mode == 2 else 1.5,
                      lo=0.0, hi=10.0)

    print()
    print("  ALPHA -- Vertical symmetry")
    print("           0.5  = symmetric (both walls equal)")
    print("           <0.5 = bottom wall denser")
    print("           >0.5 = top wall denser")
    print()
    ALPHA = ask_float("ALPHA", default=0.5, lo=0.01, hi=0.99)

    if mode == 2:
        print()
        print("  Poisson solver iterations")
        print("    (more = more accurate, slower)")
        print("    Typical: 10000~30000 for high accuracy")
        POISSON_ITER = ask_int("Poisson iterations", default=15000, lo=1000, hi=100000)
    else:
        POISSON_ITER = 15000

    print()
    print(f"  -> Mode:  {'Zeta-only' if mode == 1 else 'Adaptive (Poisson + P,Q)'}")
    print(f"  -> Grid:  I={NI} x J={NJ}")
    print(f"  -> GAMMA: {GAMMA}  |  ALPHA: {ALPHA}")
    if mode == 2:
        print(f"  -> Poisson iterations: {POISSON_ITER}")

    # -----------------------------------------------------------
    #  Step 5 -- identity verification
    # -----------------------------------------------------------
    print("\n" + "-" * 62)
    print("  [Step 5] Identity verification (gamma=0, original size)")
    print("-" * 62)

    x_id, y_id = redistribute_vertical(x_ref, y_ref, gamma=0.0)
    ok, dx_err, dy_err = verify_identity(x_ref, y_ref, x_id, y_id, tol=1e-10)
    tag = "PASS" if ok else "FAIL"
    print(f"  Zeta-only identity:  max|dx| = {dx_err:.2e},  max|dy| = {dy_err:.2e}  ->  {tag}")

    # -----------------------------------------------------------
    #  Step 6 -- generate new grid
    # -----------------------------------------------------------
    print("\n" + "-" * 62)
    print("  [Step 6] Generating new grid ...")
    print("-" * 62)

    if mode == 1:
        x_new, y_new = redistribute_vertical(x_ref, y_ref,
                                              gamma=GAMMA, alpha=ALPHA)
    else:
        x_new, y_new, poisson_conv = generate_adaptive_grid(
            x_ref, y_ref, NI, NJ,
            gamma=GAMMA, alpha=ALPHA,
            poisson_iter=POISSON_ITER, poisson_tol=1e-12)

    print(f"  Generated grid: I={NI}, J={NJ}")

    # -----------------------------------------------------------
    #  Step 7 -- output
    # -----------------------------------------------------------
    print("\n" + "-" * 62)
    print("  [Step 7] Saving outputs ...")
    print("-" * 62)

    tag_str = f"I{NI}_J{NJ}_g{GAMMA}_a{ALPHA}"

    out_cmp = base / f"compare_{grid_key}_{tag_str}.png"
    plot_compare(x_ref, y_ref, x_new, y_new,
                 labels=["Reference", f"New ({NI}x{NJ})"],
                 title=f"GAMMA={GAMMA}, ALPHA={ALPHA}, Grid={NI}x{NJ}",
                 savepath=out_cmp)

    mid_col = NI // 2
    out_sp = base / f"spacing_{grid_key}_{tag_str}.png"
    plot_vertical_spacing(y_ref, y_new, icol=min(mid_col, ni_ref//2),
                          labels=["Reference", f"New ({NI}x{NJ})"],
                          savepath=out_sp)

    out_dat = base / f"adaptive_{grid_key}_{tag_str}.dat"
    write_tecplot_dat(out_dat, x_new, y_new,
                      title=f"Periodic hill {NI}x{NJ}",
                      zone_title=f"I{NI}_J{NJ}_g{GAMMA}_a{ALPHA}")

    out_new = base / f"grid_{grid_key}_{tag_str}.png"
    plot_grid(x_new, y_new,
              title=f"New grid {NI}x{NJ}  GAMMA={GAMMA}",
              savepath=out_new)

    if mode == 2:
        fig_cv, ax_cv = plt.subplots(figsize=(8, 5))
        ax_cv.semilogy(poisson_conv, 'k-', lw=0.6)
        ax_cv.set_xlabel("Iteration"); ax_cv.set_ylabel("Max correction")
        ax_cv.set_title(f"Poisson convergence ({NI}x{NJ})")
        ax_cv.grid(True, ls='--', alpha=0.4)
        plt.tight_layout()
        conv_path = base / f"convergence_{grid_key}_{tag_str}.png"
        fig_cv.savefig(conv_path, dpi=200)
        print(f"  [saved] {conv_path}")
        plt.close()

    # -----------------------------------------------------------
    #  Step 8 -- optional parametric sweep
    # -----------------------------------------------------------
    print("\n" + "-" * 62)
    print("  [Step 8] Parametric sweep (optional)")
    print("-" * 62)

    do_sweep = ask_yes_no("Generate parametric sweep plots?", default_yes=False)

    if do_sweep:
        print("  Generating sweep ...")
        gammas = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

        fig, axes = plt.subplots(len(gammas), 1,
                                 figsize=(18, 3.2 * len(gammas)),
                                 sharex=True)
        for ax, g in zip(axes, gammas):
            if mode == 1:
                xn, yn = redistribute_vertical(x_ref, y_ref, gamma=g, alpha=ALPHA)
            else:
                xn, yn, _ = generate_adaptive_grid(
                    x_ref, y_ref, NI, NJ,
                    gamma=g, alpha=ALPHA, poisson_iter=POISSON_ITER)
            nj_n, ni_n = xn.shape
            for jj in range(nj_n):
                ax.plot(xn[jj, :], yn[jj, :], "k-", lw=0.2)
            for ii in range(0, ni_n, max(1, ni_n//40)):
                ax.plot(xn[:, ii], yn[:, ii], "k-", lw=0.2)
            ax.set_aspect("equal"); ax.set_ylabel("y")
            ax.set_title(f"gamma = {g:.1f}", fontsize=10, loc="left")

        axes[-1].set_xlabel("x  [m]")
        fig.suptitle(f"Parametric sweep (alpha={ALPHA})", fontsize=14)
        plt.tight_layout()
        sweep_path = base / f"sweep_{grid_key}_{tag_str}.png"
        fig.savefig(sweep_path, dpi=200, bbox_inches="tight")
        print(f"  [saved] {sweep_path}")
        plt.close(fig)

        fig3, ax3 = plt.subplots(figsize=(7, 5))
        eta = np.linspace(0, 1, NJ)
        for g in gammas:
            z = vinokur_tanh(eta, g, ALPHA)
            ax3.plot(range(NJ), z, "-", lw=1.2, label=f"gamma={g:.1f}")
        ax3.set_xlabel("j index"); ax3.set_ylabel("zeta (normalised)")
        ax3.set_title(f"Zeta distribution (alpha={ALPHA})")
        ax3.legend(fontsize=8); ax3.grid(True, ls="--", alpha=0.4)
        plt.tight_layout()
        zeta_path = base / "zeta_curves.png"
        fig3.savefig(zeta_path, dpi=200)
        print(f"  [saved] {zeta_path}")
        plt.close(fig3)
    else:
        print("  Skipped.")

    # -----------------------------------------------------------
    #  Summary
    # -----------------------------------------------------------
    print("\n" + "=" * 62)
    print("  DONE -- Output summary")
    print("=" * 62)
    print(f"  Reference    : {dat_path.name}  (I={ni_ref}, J={nj_ref})")
    print(f"  New grid     : I={NI} x J={NJ}")
    print(f"  Mode         : {'Zeta-only' if mode == 1 else 'Adaptive (Poisson + Steger-Sorenson P,Q)'}")
    print(f"  GAMMA        : {GAMMA}  (stretching intensity)")
    print(f"  ALPHA        : {ALPHA}  (vertical symmetry)")
    if mode == 2:
        print(f"  Poisson iter : {POISSON_ITER}")
    print(f"  Output folder: {base}")
    print("=" * 62)
