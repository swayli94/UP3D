"""AUDIT 2026-07-28 / experiment 3 -- is the transonic answer mesh-converging
to a dissipation-INDEPENDENT solution?

Experiment 1 showed that on the coarse 2.5-D NACA0012 mesh the shock position
spans ~0.09 chord and cl spans ~50 % as the artificial-density constant C is
varied over 1.0..3.0, with every run reported "converged". That is expected of
ANY artificial-viscosity scheme at fixed h; the question that decides whether
it is acceptable is whether the spread SHRINKS under mesh refinement (a
consistent O(h) dissipation) or stays put (an inconsistent operator, i.e. the
answer is set by the tuning, forever).

Runs the same Newton recipe at C in {1.0, 1.5, 3.0} on coarse/medium/fine.

Outputs: results/exp3_dissipation.csv
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
from pyfp3d.solve.newton import solve_newton_transonic         # noqa: E402

M_INF, ALPHA = 0.80, 1.25
LEVELS = ["coarse", "medium", "fine"]
CS = [1.0, 1.5, 3.0]
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    rows = []
    for level in LEVELS:
        path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
        if not path.exists():
            print(f"skip {level}: no mesh")
            continue
        mesh = read_mesh(path)
        mc, wc = cut_wake(mesh)
        dz = float(np.ptp(mc.nodes[:, 2]))
        print(f"\n=== {level}: {len(mc.nodes)} nodes ===", flush=True)
        for C in CS:
            t0 = time.perf_counter()
            try:
                r = solve_newton_transonic(
                    mc, wc, m_inf=M_INF, alpha_deg=ALPHA, m_start=0.70,
                    dm=0.025, dm_min=0.003, upwind_c=C, freeze_tol=1e-6,
                    newton_kw=dict(freeze_refresh_max=8, precond="direct",
                                   n_newton_max=60))
            except Exception as exc:                          # noqa: BLE001
                print(f"  C={C}: FAILED {type(exc).__name__}: {exc}",
                      flush=True)
                rows.append(dict(level=level, n_dof=len(mc.nodes), C=C,
                                 converged=False,
                                 error=f"{type(exc).__name__}"))
                continue
            dt = time.perf_counter() - t0
            curve = wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=M_INF)
            rep = shock_report(curve, M_INF)
            f = wall_force_coefficients(mc.nodes, mc.elements,
                                        mc.boundary_faces["wall"], r["phi"],
                                        alpha_deg=ALPHA, s_ref=dz,
                                        m_inf=M_INF)
            row = dict(level=level, n_dof=len(mc.nodes), C=C,
                       converged=r["converged"], wall_s=round(dt, 2),
                       cl_p=round(f["cl"], 5),
                       cl_kj=round(2.0 * float(r["gamma"][0]), 5),
                       x_shock_up=rep["upper"].get("x_shock"),
                       n_cells_up=rep["upper"].get("n_cells"),
                       x_shock_lo=rep["lower"].get("x_shock"),
                       m_max=round(float(np.sqrt(r["mach2_max"])), 4),
                       froze=r["froze"], error="")
            print("  ", row, flush=True)
            rows.append(row)

    keys = sorted({k for row in rows for k in row})
    with open(OUT / "exp3_dissipation.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["level", "C"]
                           + [k for k in keys if k not in ("level", "C")])
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT / "exp3_dissipation.csv")

    # spread per level -- the headline of this experiment
    print("\nspread of x_shock / cl_p across C at each level:")
    for level in LEVELS:
        sub = [r for r in rows if r.get("level") == level
               and r.get("x_shock_up") is not None]
        if len(sub) < 2:
            continue
        xs = [r["x_shock_up"] for r in sub]
        cls = [r["cl_p"] for r in sub]
        print(f"  {level:6s} n_dof={sub[0]['n_dof']:6d}  "
              f"x_shock {min(xs):.4f}..{max(xs):.4f} (spread {max(xs)-min(xs):.4f})  "
              f"cl_p {min(cls):.4f}..{max(cls):.4f} "
              f"(spread {100*(max(cls)-min(cls))/min(cls):.1f} %)")


if __name__ == "__main__":
    main()
