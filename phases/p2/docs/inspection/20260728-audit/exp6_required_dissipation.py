"""AUDIT 2026-07-28 / experiment 6 -- how much artificial dissipation would be
needed to hit the canonical Euler anchor, and is that amount mesh-independent?

Experiments 1/3/4 showed the converged transonic answer moves monotonically
with the artificial-density constant C, and that the shipped default C = 1.5
puts NACA0012 M0.80/alpha1.25 at cl 0.459 / shock 0.658 -- versus the canonical
Euler anchor cl ~ 0.363 / shock 0.60-0.63 (the same anchor the project's own
G4.1 reference CSV is built on).

`kernels/upwind.py`'s own docstring records that a face-neighbour hop on these
prism-split sliver tets reaches only 0.25-0.39 of an element's streamwise
extent, i.e. the upwind distance in rho~ = rho - nu*Dl*drho/dl is delivered
SHORT. If that is the mechanism, the C needed to reproduce the anchor should be
a few times the literature value C in [1,2] -- and should be roughly the same
number on both meshes. This experiment measures it.

Outputs: results/exp6_required_C.csv
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

ALPHA = 1.25
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)

CASES = [
    ("coarse", 0.80, 4.0), ("coarse", 0.80, 5.0), ("coarse", 0.80, 6.0),
    ("coarse", 0.7875, 1.5), ("coarse", 0.7875, 3.0),
    ("medium", 0.7875, 4.0), ("medium", 0.7875, 5.0),
    ("medium", 0.80, 2.0),
]


def main():
    cache = {}
    rows = []
    for level, m_inf, C in CASES:
        if level not in cache:
            mesh = read_mesh(REPO / f"cases/meshes/naca0012_2.5d/{level}.msh")
            cache[level] = cut_wake(mesh)
        mc, wc = cache[level]
        dz = float(np.ptp(mc.nodes[:, 2]))
        t0 = time.perf_counter()
        try:
            r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=ALPHA,
                                       m_start=0.70, dm=0.025, dm_min=0.003,
                                       upwind_c=C, freeze_tol=1e-6,
                                       newton_kw=dict(freeze_refresh_max=8,
                                                      precond="direct",
                                                      n_newton_max=60))
        except Exception as exc:                              # noqa: BLE001
            print(f"  {level} M{m_inf} C={C}: FAILED {type(exc).__name__}",
                  flush=True)
            rows.append(dict(level=level, m_inf=m_inf, C=C, converged=False,
                             error=type(exc).__name__))
            continue
        dt = time.perf_counter() - t0
        done = [lr["m"] for lr in r["level_results"] if lr["converged"]]
        curve = wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=m_inf)
        rep = shock_report(curve, m_inf)
        f = wall_force_coefficients(mc.nodes, mc.elements,
                                    mc.boundary_faces["wall"], r["phi"],
                                    alpha_deg=ALPHA, s_ref=dz, m_inf=m_inf)
        row = dict(level=level, n_dof=len(mc.nodes), m_inf=m_inf, C=C,
                   converged=r["converged"], wall_s=round(dt, 1),
                   m_reached=max(done) if done else None,
                   cl_p=round(f["cl"], 5),
                   cl_kj=round(2 * float(r["gamma"][0]), 5),
                   x_shock_up=rep["upper"].get("x_shock"),
                   n_cells_up=rep["upper"].get("n_cells"),
                   m_max=round(float(np.sqrt(r["mach2_max"])), 4), error="")
        print("  ", row, flush=True)
        rows.append(row)
    keys = sorted({k for r_ in rows for k in r_})
    with open(OUT / "exp6_required_C.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["level", "m_inf", "C"]
                           + [k for k in keys
                              if k not in ("level", "m_inf", "C")])
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT / "exp6_required_C.csv")


if __name__ == "__main__":
    main()
