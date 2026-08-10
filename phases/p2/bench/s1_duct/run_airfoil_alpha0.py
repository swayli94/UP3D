"""GS1.1 escalation: the SAME dissipation/mesh sweep on the production airfoil
mesh at ZERO incidence -- the discriminator between "flux operator on sliver
tets" and "lift / Kutta coupling".

Why: the nozzle bench (run_nozzle.py) measured the shipped walk flux to be
CONSISTENT -- shock error O(h), dissipation sensitivity O(h), unaffected by mesh
jitter. Yet on the lifting NACA0012 at M0.80/alpha1.25 the audit measured a
26-47 % cl spread over C and a +40 % cl jump under one mesh refinement, with no
clamps active. Something the duct does not contain is responsible. The two
candidates, in order:

  (a) the mesh + oblique flow: production 2.5-D airfoil meshes are prism-split
      slivers with edge lengths spanning ~300:1 and the flow crosses cells
      obliquely, so the walk's upstream distance varies far more than in the
      duct (where the flow is aligned with the grid lines);
  (b) the lift / Kutta coupling: phase one measured that the Kutta map
      amplifies a closure shift into Gamma by 1/(1-b) ~ 14x (P14 / G14.7), so a
      sub-cell shock movement can produce a large lift change.

At alpha = 0 the NACA0012 is symmetric, Gamma == 0, and the Kutta row is inert:
(b) is switched OFF while (a) is unchanged. So:

  * if the shock is mesh-convergent and C-insensitive at alpha = 0, the flux
    operator is fine on production meshes too and the airfoil pathology lives
    in the lift/Kutta coupling -> S1's problem statement must be rewritten;
  * if it is pathological at alpha = 0, the mesh/obliqueness interaction with
    the walk is the root cause -> GS1.2 proceeds as planned.

Conditions: M 0.80 and 0.82, alpha 0. Both are BELOW the documented
full-potential non-uniqueness band for the NACA0012 at low lift
(M ~ 0.82-0.85, design.md Sec 1), so the test is not run inside the
literature's ambiguous zone.

Outputs: results/gs1_1_alpha0.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                      # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve              # noqa: E402
from pyfp3d.post.shock import shock_report                     # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients        # noqa: E402
from pyfp3d.solve.newton import solve_newton_transonic         # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA = 0.0
MACHS = (0.80, 0.82)
CS = (1.0, 1.5, 2.0, 3.0)
LEVELS = ("coarse", "medium")
RECIPE = dict(m_start=0.70, dm=0.025, dm_min=0.003, freeze_tol=1e-6,
              newton_kw=dict(freeze_refresh_max=8, precond="direct",
                             n_newton_max=60))


def main():
    rows = []
    for level in LEVELS:
        path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
        if not path.exists():
            print(f"skip {level}: mesh missing")
            continue
        mc, wc = cut_wake(read_mesh(path))
        dz = float(np.ptp(mc.nodes[:, 2]))
        print(f"\n=== {level}: {len(mc.nodes)} nodes ===", flush=True)
        for m_inf in MACHS:
            for C in CS:
                t0 = time.perf_counter()
                try:
                    r = solve_newton_transonic(mc, wc, m_inf=m_inf,
                                               alpha_deg=ALPHA, upwind_c=C,
                                               **RECIPE)
                    err = ""
                except Exception as exc:                       # noqa: BLE001
                    rows.append(dict(level=level, n_dof=len(mc.nodes),
                                     m_inf=m_inf, C=C, converged=False,
                                     error=type(exc).__name__))
                    print("  ", rows[-1], flush=True)
                    continue
                wall = time.perf_counter() - t0
                rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                                m_inf=m_inf), m_inf)
                f = wall_force_coefficients(mc.nodes, mc.elements,
                                            mc.boundary_faces["wall"],
                                            r["phi"], alpha_deg=ALPHA,
                                            s_ref=dz, m_inf=m_inf)
                up, lo = rep["upper"], rep["lower"]
                row = dict(
                    level=level, n_dof=len(mc.nodes), m_inf=m_inf, C=C,
                    converged=r["converged"], wall_s=round(wall, 2),
                    x_shock_up=up.get("x_shock"), x_shock_lo=lo.get("x_shock"),
                    n_cells_up=up.get("n_cells"),
                    cp_min_up=round(up.get("cp_min", float("nan")), 5),
                    cl_p=round(f["cl"], 6),
                    gamma=round(float(r["gamma"][0]), 8),
                    m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                    n_limited=r["n_limited"], n_floored=r["n_floored"],
                    res_final=r["residual_history"][-1], error=err)
                print("  ", {k: row[k] for k in
                             ("m_inf", "C", "converged", "x_shock_up",
                              "x_shock_lo", "m_max", "gamma", "n_floored",
                              "n_limited")}, flush=True)
                rows.append(row)

    with open(OUT / "gs1_1_alpha0.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1_1_alpha0.csv")

    print("\n=== headline: upper-surface shock x/c at alpha = 0 ===")
    for m_inf in MACHS:
        print(f"  M {m_inf}")
        print(f"    {'level':8s} " + " ".join(f"{'C=' + str(c):>9s}"
                                             for c in CS) + "    spread")
        for level in LEVELS:
            vals = []
            for C in CS:
                s = [r for r in rows if r.get("level") == level
                     and r.get("m_inf") == m_inf and r.get("C") == C]
                v = s[0].get("x_shock_up") if s else None
                vals.append(v)
            cells = [f"{v:9.4f}" if v is not None else f"{'--':>9s}"
                     for v in vals]
            good = [v for v in vals if v is not None]
            spread = (max(good) - min(good)) if len(good) > 1 else float("nan")
            print(f"    {level:8s} " + " ".join(cells) + f"    {spread:.4f}")


if __name__ == "__main__":
    main()
