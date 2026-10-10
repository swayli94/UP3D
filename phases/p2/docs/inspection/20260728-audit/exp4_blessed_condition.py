"""AUDIT 2026-07-28 / experiment 4 -- the same dissipation / cross-path check
at the project's OWN blessed gate condition.

tests/test_p8_newton.py re-specced the medium-mesh transonic gate away from the
canonical M0.80 to M0.7875, on the stated grounds that M0.80 medium "sits at
the edge of the FP non-uniqueness fold". Fair enough -- so this experiment asks
the same two questions AT THE BLESSED CONDITION:

  (1) how much does the answer move with the artificial-dissipation constant C?
  (2) do the Newton and Picard paths land on the same solution?

If the spread is small here, the M0.80 behaviour is a fold artifact. If it is
comparable to M0.80, then the transonic answer is set by the tuning generally.

Outputs: results/exp4_blessed.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                      # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve              # noqa: E402
from pyfp3d.post.shock import shock_report                     # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients        # noqa: E402
from pyfp3d.solve.continuation import solve_transonic_lifting  # noqa: E402
from pyfp3d.solve.newton import solve_newton_transonic         # noqa: E402

M_INF, ALPHA = 0.7875, 1.25          # the blessed G8.1 medium condition
MESH = "cases/meshes/naca0012_2.5d/medium.msh"
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)


def post(mc, phi, gamma, dz):
    curve = wall_cp_curve(mc, phi, z=0.5 * dz, m_inf=M_INF)
    rep = shock_report(curve, M_INF)
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], phi,
                                alpha_deg=ALPHA, s_ref=dz, m_inf=M_INF)
    return dict(cl_p=round(f["cl"], 5),
                cl_kj=round(2.0 * float(np.atleast_1d(gamma)[0]), 5),
                x_shock_up=rep["upper"].get("x_shock"),
                n_cells_up=rep["upper"].get("n_cells"),
                x_shock_lo=rep["lower"].get("x_shock"),
                cp_min_up=round(rep["upper"].get("cp_min", float("nan")), 4))


def main():
    mesh = read_mesh(REPO / MESH)
    mc, wc = cut_wake(mesh)
    dz = float(np.ptp(mc.nodes[:, 2]))
    print(f"medium mesh: {len(mc.nodes)} nodes, M={M_INF}, alpha={ALPHA}",
          flush=True)
    rows = []

    for C in (1.0, 1.5, 2.0, 3.0):
        t0 = time.perf_counter()
        try:
            r = solve_newton_transonic(
                mc, wc, m_inf=M_INF, alpha_deg=ALPHA, m_start=0.70,
                dm=0.025, dm_min=0.003, upwind_c=C, freeze_tol=1e-6,
                newton_kw=dict(freeze_refresh_max=8, precond="direct",
                               n_newton_max=60))
        except Exception as exc:                              # noqa: BLE001
            print(f"  newton C={C}: FAILED {type(exc).__name__}", flush=True)
            rows.append(dict(case=f"newton_C{C}", C=C, converged=False,
                             error=type(exc).__name__))
            continue
        row = dict(case=f"newton_C{C}", C=C, converged=r["converged"],
                   wall_s=round(time.perf_counter() - t0, 1),
                   m_max=round(float(np.sqrt(r["mach2_max"])), 4),
                   froze=r["froze"], error="")
        row.update(post(mc, r["phi"], r["gamma"], dz))
        print("  ", row, flush=True)
        rows.append(row)

    t0 = time.perf_counter()
    r = solve_transonic_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                max_gamma_evals=12, n_picard_eval=800)
    row = dict(case="picard", C=1.5, converged=r["converged"],
               wall_s=round(time.perf_counter() - t0, 1),
               m_max=round(float(np.sqrt(r["mach2_max"])), 4),
               froze=False, error="",
               kutta_mismatch=r["kutta_mismatch"])
    row.update(post(mc, r["phi"], r["gamma"], dz))
    print("  ", row, flush=True)
    rows.append(row)

    keys = sorted({k for row in rows for k in row})
    with open(OUT / "exp4_blessed.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["case"] + [k for k in keys
                                                      if k != "case"])
        w.writeheader()
        w.writerows(rows)
    xs = [r["x_shock_up"] for r in rows if r.get("x_shock_up")]
    cls = [r["cl_p"] for r in rows if r.get("cl_p")]
    print(f"\nx_shock spread {min(xs):.4f}..{max(xs):.4f} "
          f"({max(xs)-min(xs):.4f} chord); cl_p spread {min(cls):.4f}.."
          f"{max(cls):.4f} ({100*(max(cls)-min(cls))/min(cls):.1f} %)")


if __name__ == "__main__":
    main()
