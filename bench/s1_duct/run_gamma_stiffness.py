"""GS1.3b step 1: the Gamma-row stiffness as a MATRIX quantity, both closure
renderings, three mesh levels.

GS1.3 measured the closure amplification by perturbing the Kutta target
(tip_taper = 1+eps) and reading dGamma: A = 15 / 22 / 33 on coarse / medium /
fine at SUBSONIC M0.5, i.e. A ~ h^(-0.56). That is a response measurement on
one rendering. This script computes the same thing directly from the Jacobian
blocks, which (a) cross-validates it, (b) works for the `pressure` rendering
where no clean perturbation vehicle exists, and (c) exposes the MECHANISM.

Algebra (newton.py::kutta_blocks): the coupled system is

    J_ff dphi + B dGamma = -R
    K   dphi + D dGamma  = -F

The driver eliminates the second row (K~ = -D^-1 K, F~ = -D^-1 F) and solves
(J_ff + B K~) dphi = -R - B F~ via Woodbury with

    S_code = I + K~ J_ff^-1 B          and      S_Gamma = D . S_code

so a closure-equation perturbation dF moves Gamma by dGamma = S_Gamma^-1 dF.

Rendering-independent yardstick -- "the closure condition is violated by a
RELATIVE amount eps":

  * probe:    F = jump - Gamma, so dF = eps*Gamma and D = -1
              =>  A = 1 / |S_code|                     (dimensionless already)
  * pressure: F_raw = |q_u|^2 - |q_l|^2, so dF_raw = eps*u_inf^2
              =>  A = u_inf^2 / (|D_raw| * |S_code| * Gamma)

Mechanism read-out: for the probe S_code = 1 + K~ J^-1 B, and A -> inf means
K~ J_ff^-1 B -> -1, i.e. the field feedback cancels the explicit -Gamma term.
The printed `feedback` column is exactly that product.

Outputs: results/gs1_3b_stiffness.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                      # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients        # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting           # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA = 1.25
LEVELS = (("coarse", 0.02), ("medium", 0.01), ("fine", 0.005))
H_REF = {"coarse": 0.02, "medium": 0.01, "fine": 0.005}
UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR, U_INF = 1.5, 0.95, 3.0, 0.05, 1.0


def stiffness(ws, phi_cut, gamma):
    """(S_code, D_raw, feedback) at the given state."""
    phi_free = np.asarray(phi_cut, dtype=np.float64)[:ws.n_red][ws.free].copy()
    R, F, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT, M_CAP,
                                   RHO_FLOOR)
    J_ff, B = ws.assemble_coupled(state, UPWIND_C, M_CRIT, RHO_FLOOR)
    K_t, _ = ws.kutta_blocks(state, F)
    JB = spla.spsolve(J_ff.tocsc(), B.toarray())
    if JB.ndim == 1:
        JB = JB[:, None]
    feedback = np.asarray(K_t @ JB, dtype=np.float64).reshape(ws.n_st, ws.n_st)
    S_code = np.eye(ws.n_st) + feedback
    D_raw = None
    if ws.kutta_estimator == "pressure":
        _, D_raw = ws.cvs.newton_rows(state["phi_cut"])
        D_raw = np.asarray(D_raw, dtype=np.float64).reshape(ws.n_st, ws.n_st)
    return S_code, D_raw, feedback, float(np.max(np.abs(R)))


def main():
    rows = []
    for m_inf, tag in ((0.50, "subsonic"), (0.7875, "transonic")):
        for level, _ in LEVELS:
            path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
            if not path.exists():
                continue
            mc, wc = cut_wake(read_mesh(path))
            dz = float(np.ptp(mc.nodes[:, 2]))
            for est in ("probe", "pressure"):
                t0 = time.perf_counter()
                try:
                    r = solve_newton_lifting(
                        mc, wc, m_inf=m_inf, alpha_deg=ALPHA,
                        upwind_c=UPWIND_C, m_crit=M_CRIT, freeze_tol=1e-6,
                        freeze_refresh_max=8, precond="direct",
                        direct_refactor_every=4, n_newton_max=80,
                        kutta_estimator=est)
                except Exception as exc:                       # noqa: BLE001
                    print(f"  {m_inf} {level} {est}: FAILED "
                          f"{type(exc).__name__}", flush=True)
                    continue
                ws = r["workspace"]
                g = np.atleast_1d(np.asarray(r["gamma"], dtype=np.float64))
                S, D_raw, fb, rmax = stiffness(ws, r["phi"], g)
                s_abs = float(np.abs(np.linalg.det(S))) if ws.n_st > 1 \
                    else float(abs(S[0, 0]))
                if est == "probe":
                    amp = 1.0 / s_abs
                    d_abs = 1.0
                else:
                    d_abs = float(abs(D_raw[0, 0])) if ws.n_st == 1 \
                        else float(np.abs(np.linalg.det(D_raw)))
                    amp = U_INF ** 2 / (d_abs * s_abs * float(g[0]))
                f = wall_force_coefficients(
                    mc.nodes, mc.elements, mc.boundary_faces["wall"],
                    r["phi"], alpha_deg=ALPHA, s_ref=dz, m_inf=m_inf)
                row = dict(regime=tag, m_inf=m_inf, level=level,
                           h_wall=H_REF[level], estimator=est,
                           n_dof=len(mc.nodes), converged=r["converged"],
                           n_newton=r["n_newton"],
                           res_at_state=rmax,
                           gamma=round(float(g[0]), 8),
                           cl_p=round(f["cl"], 6),
                           S_code=round(s_abs, 8),
                           D_raw=round(d_abs, 8),
                           feedback=round(float(fb[0, 0]), 8),
                           amplification=round(amp, 4),
                           wall_s=round(time.perf_counter() - t0, 1))
                print(f"  {tag:9s} {level:7s} {est:9s} "
                      f"|S|={s_abs:.4e} feedback={fb[0, 0]:+.6f} "
                      f"D_raw={d_abs:.4e} gamma={g[0]:.6f} "
                      f"-> A={amp:9.3f}  (conv={r['converged']}, "
                      f"cl={f['cl']:.5f})", flush=True)
                rows.append(row)

    with open(OUT / "gs1_3b_stiffness.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1_3b_stiffness.csv")

    print("\n=== headline: amplification A vs h, per rendering ===")
    for tag in ("subsonic", "transonic"):
        for est in ("probe", "pressure"):
            sub = [r for r in rows if r["regime"] == tag
                   and r["estimator"] == est]
            if len(sub) < 2:
                continue
            sub.sort(key=lambda r: -r["h_wall"])
            aa = [r["amplification"] for r in sub]
            hh = [r["h_wall"] for r in sub]
            ratios = [aa[i + 1] / aa[i] for i in range(len(aa) - 1)]
            p = (np.polyfit(np.log(hh), np.log(aa), 1)[0]
                 if len(aa) > 2 else float("nan"))
            print(f"  {tag:9s} {est:9s} A = "
                  + " -> ".join(f"{a:.2f}" for a in aa)
                  + "   ratios " + ", ".join(f"{x:.2f}" for x in ratios)
                  + f"   fit A ~ h^{p:.2f}")
    print("\n  cross-check vs GS1.3 eps-perturbation (probe, subsonic): "
          "14.94 / 22.05 / 32.44")


if __name__ == "__main__":
    main()
