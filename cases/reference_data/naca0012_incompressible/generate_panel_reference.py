"""
Reference-data generator: incompressible inviscid NACA0012 by a 2D
Hess-Smith panel method (constant-strength source per panel + single
constant vortex strength, Kutta closure).

Provenance / method (see README.md):
  - Geometry: the SAME closed-TE NACA0012 coordinate set the solver
    meshes use (pyfp3d.meshgen.planar.naca0012_coordinates, last
    coefficient -0.1036), so the FEM-vs-panel comparison of gate G2.3 is
    method-vs-method on identical geometry, not geometry-vs-geometry.
  - Discretization: N cosine-clustered nodes, panels traversed clockwise
    (TE -> lower -> LE -> upper -> TE), control points at panel midpoints,
    outward normal n = (-t_y, t_x).
  - Influence: analytic constant-source / constant-vortex panel velocity
    in panel-local coordinates; self-influence by the y -> 0+ limit
    (source: normal sigma/2; vortex: tangential -tau/2 on the outward side).
  - System: N tangency equations + Kutta (V_t,first + V_t,last = 0).
  - cl primary from Cp integration (cl = -oint Cp n dS . lift_dir);
    cross-checked against Kutta-Joukowski 2*Gamma/(U c) and against panel
    count convergence (N = 100..800, Richardson-style tail reported).

Writes: cl_reference.csv, cp_alpha4.csv, convergence.csv.
Run:    python generate_panel_reference.py
"""

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from pyfp3d.meshgen.planar import naca0012_coordinates  # noqa: E402


def _panel_geometry(n_half: int):
    """Clockwise closed polygon: TE -> lower -> LE -> upper -> TE."""
    coords = naca0012_coordinates(n_half=n_half)  # TE->upper->LE->lower->TE
    pts = coords[::-1].copy()                     # TE->lower->LE->upper->TE
    x0, y0 = pts[:-1, 0], pts[:-1, 1]
    x1, y1 = pts[1:, 0], pts[1:, 1]
    dx, dy = x1 - x0, y1 - y0
    length = np.hypot(dx, dy)
    t_hat = np.stack([dx / length, dy / length], axis=1)
    n_hat = np.stack([-t_hat[:, 1], t_hat[:, 0]], axis=1)  # outward (cw poly)
    xm, ym = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    return x0, y0, t_hat, n_hat, xm, ym, length


def _influence_matrices(x0, y0, t_hat, n_hat, xm, ym, length):
    """(u, v) at every control point per unit source / vortex on each panel.

    Panel-local fields (density 1, panel along local x in [0, L], local +y
    along the outward normal):
        source: u' = ln(r1^2/r2^2)/(4 pi),  v' = (th2 - th1)/(2 pi)
        vortex (counterclockwise-positive): (u', v') = (-v'_src, u'_src)
    with r1, th1 = atan2(y', x') from the panel start and r2,
    th2 = atan2(y', x' - L) from the panel end. Self-influence is the
    outward-side (y' -> 0+) limit: source v' = +1/2 (outflow), vortex
    u' = -1/2.
    """
    n = len(xm)
    dxg = xm[:, None] - x0[None, :]
    dyg = ym[:, None] - y0[None, :]
    tx, ty = t_hat[:, 0], t_hat[:, 1]
    xl = dxg * tx[None, :] + dyg * ty[None, :]
    yl = -dxg * ty[None, :] + dyg * tx[None, :]

    r1sq = xl**2 + yl**2
    r2sq = (xl - length[None, :]) ** 2 + yl**2
    th1 = np.arctan2(yl, xl)
    th2 = np.arctan2(yl, xl - length[None, :])

    us_l = np.log(r1sq / r2sq) / (4.0 * np.pi)
    vs_l = (th2 - th1) / (2.0 * np.pi)

    diag = np.arange(n)
    us_l[diag, diag] = 0.0
    # Control point ON its own panel: at exactly y' = 0 the atan2 pair
    # already evaluates to the outward-side (+y') limit th1 = 0, th2 = pi,
    # i.e. v' = +1/2; set it explicitly anyway to be robust against -0.0.
    vs_l[diag, diag] = 0.5

    uv_l = -vs_l
    vv_l = us_l

    tx, ty = t_hat[:, 0], t_hat[:, 1]
    us = us_l * tx[None, :] - vs_l * ty[None, :]
    vs = us_l * ty[None, :] + vs_l * tx[None, :]
    uv = uv_l * tx[None, :] - vv_l * ty[None, :]
    vv = uv_l * ty[None, :] + vv_l * tx[None, :]
    return us, vs, uv, vv


def solve_hess_smith(alpha_deg: float, n_half: int):
    """Returns dict with cl (Cp-integrated), cl_kj, cp at control points."""
    x0, y0, t_hat, n_hat, xm, ym, length = _panel_geometry(n_half)
    n = len(xm)
    us, vs, uv, vv = _influence_matrices(x0, y0, t_hat, n_hat, xm, ym, length)

    a = np.deg2rad(alpha_deg)
    vinf = np.array([np.cos(a), np.sin(a)])

    A = np.zeros((n + 1, n + 1))
    rhs = np.zeros(n + 1)
    # Tangency: n_i . V_i = 0.
    A[:n, :n] = n_hat[:, 0:1] * us + n_hat[:, 1:2] * vs
    A[:n, n] = (n_hat[:, 0:1] * uv + n_hat[:, 1:2] * vv).sum(axis=1)
    rhs[:n] = -(n_hat @ vinf)
    # Kutta: tangential velocities on the two TE panels cancel.
    for i in (0, n - 1):
        A[n, :n] += t_hat[i, 0] * us[i] + t_hat[i, 1] * vs[i]
        A[n, n] += (t_hat[i, 0] * uv[i] + t_hat[i, 1] * vv[i]).sum()
    rhs[n] = -(t_hat[0] @ vinf) - (t_hat[n - 1] @ vinf)

    sol = np.linalg.solve(A, rhs)
    sigma, tau = sol[:n], sol[n]

    vt = (
        t_hat[:, 0] * (us @ sigma + tau * uv.sum(axis=1) + vinf[0])
        + t_hat[:, 1] * (vs @ sigma + tau * vv.sum(axis=1) + vinf[1])
    )
    cp = 1.0 - vt**2

    # Cp-integrated force (chord = 1): dC = -Cp n dS.
    cf = -(cp * length) @ n_hat
    lift_dir = np.array([-np.sin(a), np.cos(a)])
    cl = float(cf @ lift_dir)
    # Kutta-Joukowski cross-check: the sheet's total circulation is
    # Gamma_ccw = tau * perimeter (ccw-positive kernel), and
    # L' = -rho U Gamma_ccw, so cl = -2 tau * perimeter (U = c = 1).
    cl_kj = float(-2.0 * tau * length.sum())
    return {
        "cl": cl, "cl_kj": cl_kj, "cp": cp, "xm": xm, "ym": ym,
        "n_panels": n,
    }


def main():
    out = Path(__file__).parent
    alphas = [0.0, 2.0, 4.0, 6.0]
    n_half_levels = [51, 101, 201, 401]   # -> ~100/200/400/800 panels

    conv_rows = []
    for a in alphas:
        for nh in n_half_levels:
            r = solve_hess_smith(a, nh)
            conv_rows.append([a, r["n_panels"], r["cl"], r["cl_kj"]])
            print(f"alpha={a:4.1f}  N={r['n_panels']:4d}  "
                  f"cl={r['cl']:.6f}  cl_KJ={r['cl_kj']:.6f}")

    with open(out / "convergence.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["alpha_deg", "n_panels", "cl_cp_integrated", "cl_kutta_joukowski"])
        w.writerows(conv_rows)

    with open(out / "cl_reference.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["alpha_deg", "cl"])
        for a in alphas:
            r = solve_hess_smith(a, n_half_levels[-1])
            w.writerow([a, f"{r['cl']:.6f}"])

    r4 = solve_hess_smith(4.0, n_half_levels[-1])
    upper = r4["ym"] > 0
    with open(out / "cp_alpha4.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_c", "cp", "surface"])
        for xm, cp, up in zip(r4["xm"], r4["cp"], upper):
            w.writerow([f"{xm:.6f}", f"{cp:.6f}", "upper" if up else "lower"])
    print(f"wrote {out}/cl_reference.csv, cp_alpha4.csv, convergence.csv")


if __name__ == "__main__":
    main()
