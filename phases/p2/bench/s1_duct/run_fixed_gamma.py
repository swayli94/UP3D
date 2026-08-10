"""GS1.6 diagnostic: the fixed-Gamma discriminator applied to h-refinement.

Criterion 0 of GS1.6 refuted the entropy-fix premise: the exact-solution nozzle
stays O(h) convergent at every shock strength up to M_shock 1.60 (|err| ratios
0.46 / 0.49 / 0.53 / 0.53), so a strong ISENTROPIC shock is not what breaks the
airfoil. The remaining evidence says the divergence needs BOTH lift and shock
strength -- each alone converges:

    M0.60 alpha1.25 (lift, no shock)   cl +0.23 % over the last refinement
    M0.80 alpha0    (shock, no lift)   x_shock 0.4640/0.5016/0.5130, converging
    M0.7875 alpha1.25 (both)           cl 0.3725/0.5234/0.7150, diverging

This script splits the two halves of the coupled problem. It solves the FIELD
only, with the circulation held FIXED at a prescribed value (no Kutta row at
all), and refines. Phase one used the same fixed-Gamma trick (A2) to separate a
measurement artefact from flow content; here it separates:

  * if x_shock and cl CONVERGE at fixed Gamma, the field solve is healthy and
    the divergence lives entirely in how Gamma is DETERMINED (closure + wake
    model) -- and that is where S1 must work;
  * if they still misbehave, the field solve itself is at fault at strong shock
    with a lifting flow, which the nozzle cannot reproduce.

Outputs: results/gs1_6_fixed_gamma.csv
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
from pyfp3d.physics.isentropic import mach_squared_field       # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve              # noqa: E402
from pyfp3d.post.shock import shock_report                     # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients        # noqa: E402
from pyfp3d.solve.newton import NewtonWorkspace                # noqa: E402
from pyfp3d.solve.picard import solve_subsonic_lifting         # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA = 1.25
UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR = 1.5, 0.95, 3.0, 0.05
#: the coarse converged circulation at each condition (GS1.3b / GS1.2b tables);
#: held FIXED across all three mesh levels so only the field is refined.
GAMMA_FIX = {0.7875: 0.189159, 0.80: 0.229487, 0.72: 0.121833}
LEVELS = ("coarse", "medium", "fine")


def field_newton(ws, phi_free, gamma, n_max=60, tol=1e-9, verbose=False):
    """Newton on the FIELD block only (Gamma frozen): J_ff dphi = -R_free."""
    R, F, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT, M_CAP,
                                   RHO_FLOOR)
    hist = [float(np.max(np.abs(R)))]
    for it in range(n_max):
        J_ff, _ = ws.assemble_coupled(state, UPWIND_C, M_CRIT, RHO_FLOOR)
        try:
            d = spla.spsolve(J_ff.tocsc(), -R)
        except Exception:                                      # noqa: BLE001
            return phi_free, state, hist, "linear solve failed"
        if not np.all(np.isfinite(d)):
            return phi_free, state, hist, "non-finite step"
        lam, best, r0 = 1.0, None, hist[-1]
        for _ in range(10):
            trial = phi_free + lam * d
            Rt, Ft, st = ws.eval_residual(trial, gamma, UPWIND_C, M_CRIT,
                                          M_CAP, RHO_FLOOR)
            rt = float(np.max(np.abs(Rt)))
            if np.isfinite(rt) and rt < r0:
                break
            if best is None or (np.isfinite(rt) and rt < best[0]):
                best = (rt, trial, Rt, Ft, st)
            lam *= 0.5
        else:
            if best is None:
                return phi_free, state, hist, "line search"
            rt, trial, Rt, Ft, st = best
        phi_free, R, F, state = trial, Rt, Ft, st
        hist.append(rt)
        if verbose:
            print(f"      {it:2d}: |R|={rt:.3e} lam={lam:g}", flush=True)
        if rt < tol:
            return phi_free, state, hist, "tol"
    return phi_free, state, hist, "cap"


def main():
    rows = []
    for m_inf in (0.7875, 0.72):
        gamma = np.array([GAMMA_FIX[m_inf]], dtype=np.float64)
        print(f"\n=== M {m_inf}, Gamma FIXED at {gamma[0]:.6f} ===", flush=True)
        for level in LEVELS:
            path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
            if not path.exists():
                continue
            mc, wc = cut_wake(read_mesh(path))
            dz = float(np.ptp(mc.nodes[:, 2]))
            ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA)
            ws.set_mach(m_inf)
            t0 = time.perf_counter()
            seed = solve_subsonic_lifting(
                mc, wc, m_inf=m_inf, alpha_deg=ALPHA, upwind_c=UPWIND_C,
                m_crit=M_CRIT, m_cap=M_CAP, n_picard_max=5, tol_rho=1e-3)
            phi_free = np.asarray(seed["phi"],
                                  dtype=np.float64)[:ws.n_red][ws.free].copy()
            phi_free, state, hist, why = field_newton(ws, phi_free, gamma)
            wall = time.perf_counter() - t0
            phi_cut = state["phi_cut"]
            rep = shock_report(wall_cp_curve(mc, phi_cut, z=0.5 * dz,
                                            m_inf=m_inf), m_inf)
            f = wall_force_coefficients(mc.nodes, mc.elements,
                                        mc.boundary_faces["wall"], phi_cut,
                                        alpha_deg=ALPHA, s_ref=dz,
                                        m_inf=m_inf)
            row = dict(m_inf=m_inf, gamma_fixed=round(float(gamma[0]), 6),
                       level=level, n_dof=len(mc.nodes),
                       reason=why, n_newton=len(hist) - 1,
                       res_final=hist[-1],
                       cl_p=round(f["cl"], 6),
                       x_shock=rep["upper"].get("x_shock"),
                       n_cells=rep["upper"].get("n_cells"),
                       n_limited=int(state["n_limited"]),
                       n_floored=int(state["n_floored"]),
                       wall_s=round(wall, 1))
            row["m_max"] = round(float(np.sqrt(np.max(
                mach_squared_field(state["q2l"], m_inf)))), 6)
            print(f"   {level:7s} {why:20s} |R|={hist[-1]:.2e} "
                  f"cl={row['cl_p']:.6f} x_shock={row['x_shock']} "
                  f"M_max={row['m_max']} nfl={row['n_floored']} "
                  f"wall={row['wall_s']}s", flush=True)
            rows.append(row)

    with open(OUT / "gs1_6_fixed_gamma.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1_6_fixed_gamma.csv")

    print("\n=== headline: at FIXED Gamma, does the field converge? ===")
    for m_inf in (0.7875, 0.72):
        sub = [r for r in rows if r["m_inf"] == m_inf]
        if len(sub) < 3:
            continue
        xs = [r["x_shock"] for r in sub]
        cl = [r["cl_p"] for r in sub]
        if all(x is not None for x in xs):
            d1, d2 = xs[1] - xs[0], xs[2] - xs[1]
            print(f"  M {m_inf}: x_shock " + " -> ".join(f"{x:.4f}" for x in xs)
                  + f"   deltas {d1:+.4f} / {d2:+.4f}"
                  + (f"   ratio {d2 / d1:.3f}" if d1 else ""))
        d1, d2 = cl[1] - cl[0], cl[2] - cl[1]
        print(f"          cl_p    " + " -> ".join(f"{c:.4f}" for c in cl)
              + f"   deltas {100 * d1 / cl[0]:+.2f} % / "
              f"{100 * d2 / cl[1]:+.2f} %"
              + (f"   ratio {d2 / d1:.3f}" if d1 else ""))


if __name__ == "__main__":
    main()
