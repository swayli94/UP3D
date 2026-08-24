"""AUDIT 2026-07-28 / experiment 2 -- where the wall clock goes, and what the
preconditioner choice costs.

The product target is "wing-body in ~2 CPU-minutes including mesh generation"
(BLWF class). The project's own only speed gate is G8.2: ONERA M6 medium
end-to-end < 300 s. This experiment runs the project's OWN G8.2 recipe
(tests/test_p8_newton.py::NEWTON_M6_RECIPE) verbatim and changes exactly ONE
thing -- `precond`: "direct" (the shipped choice: a lagged sparse LU of the
Newton Jacobian) vs "amg" -- so the cost difference is attributable to the
linear-algebra strategy alone, not to the language, the mesh or the physics.

Outputs: results/exp2_cost.csv
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
from pyfp3d.post.surface import (planform_area,                # noqa: E402
                                 wall_force_coefficients)
from pyfp3d.solve.newton import solve_newton_transonic         # noqa: E402

OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)

M_INF, ALPHA = 0.84, 3.06          # the G8.2 condition


def recipe(precond):
    """tests/test_p8_newton.py::NEWTON_M6_RECIPE with `precond` swapped."""
    return dict(dm=0.05, dm_min=0.01, freeze_tol=1e-6, intermediate_tol=1e-5,
                newton_kw=dict(freeze_refresh_max=8, precond=precond,
                               direct_refactor_every=1000, n_newton_max=60,
                               farfield_spanwise_gamma=True))


RUNS = [
    ("m6_coarse_direct", "cases/meshes/onera_m6/coarse.msh", "direct"),
    ("m6_coarse_amg", "cases/meshes/onera_m6/coarse.msh", "amg"),
    ("m6_medium_direct", "cases/meshes/onera_m6/medium.msh", "direct"),
    ("m6_medium_amg", "cases/meshes/onera_m6/medium.msh", "amg"),
]


def run(tag, mesh_rel, precond):
    t_all = time.perf_counter()
    mesh = read_mesh(REPO / mesh_rel)
    t0 = time.perf_counter()
    mc, wc = cut_wake(mesh)
    t_cut = time.perf_counter() - t0
    n_dof = len(mc.nodes)
    print(f"\n=== {tag}: {n_dof} nodes, {len(mc.elements)} tets, "
          f"cut_wake {t_cut:.1f} s ===", flush=True)
    t0 = time.perf_counter()
    try:
        r = solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                   verbose=True, **recipe(precond))
    except Exception as exc:                                # noqa: BLE001
        print("FAILED:", type(exc).__name__, exc, flush=True)
        return dict(case=tag, precond=precond, n_dof=n_dof,
                    n_tets=len(mc.elements), t_cut_wake=round(t_cut, 2),
                    wall_solve_s=round(time.perf_counter() - t0, 2),
                    converged=False,
                    error=f"{type(exc).__name__}: {exc}"[:160])
    wall = time.perf_counter() - t0
    tm = r["timings_total"]
    n_steps = sum(len(lr["residual_history"]) - 1 for lr in r["level_results"])
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], r["phi"],
                                alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
    row = dict(
        case=tag, precond=precond, n_dof=n_dof, n_tets=len(mc.elements),
        t_cut_wake=round(t_cut, 2), wall_solve_s=round(wall, 2),
        wall_total_s=round(time.perf_counter() - t_all, 2),
        converged=r["converged"], error="",
        n_levels=len(r["level_results"]), n_newton_total=n_steps,
        s_per_newton=round(wall / max(n_steps, 1), 3),
        t_seed=round(tm.get("seed", 0.0), 2),
        t_assembly=round(tm.get("assembly", 0.0), 2),
        t_precond=round(tm.get("precond", 0.0), 2),
        t_linsolve=round(tm.get("linsolve", 0.0), 2),
        t_residual=round(tm.get("residual", 0.0), 2),
        n_gmres=r["n_gmres_total"], n_refactor=r["n_refactor"],
        cl_p=round(f["cl"], 6), gamma_mean=round(float(np.mean(r["gamma"])), 6),
        m_max=round(float(np.sqrt(r["mach2_max"])), 4),
        n_limited=r["n_limited"], n_floored=r["n_floored"],
        res_final=r["residual_history"][-1])
    print("  ", row, flush=True)
    return row


def main():
    rows = [run(*r) for r in RUNS]
    keys = sorted({k for row in rows for k in row})
    with open(OUT / "exp2_cost.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["case"] + [k for k in keys
                                                     if k != "case"])
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT / "exp2_cost.csv")


if __name__ == "__main__":
    main()
