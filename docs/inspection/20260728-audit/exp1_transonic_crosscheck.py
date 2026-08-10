"""AUDIT 2026-07-28 / experiment 1 -- 2.5-D transonic cross-path consistency
and artificial-dissipation sensitivity.

Question (fresh eyes, no prior conclusions assumed): on ONE mesh at ONE
condition (NACA0012 coarse, M 0.80, alpha 1.25 -- the G4.1 gate condition),
do the two nonlinear drivers agree, and how much does the answer depend on
the tuned artificial-density constant C (upwind_c) and the switch M_crit?

Outputs: results/exp1_crosscheck.csv
Run:  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
          python docs/inspection/20260728-audit/exp1_transonic_crosscheck.py
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

MESH = REPO / "cases/meshes/naca0012_2.5d/coarse.msh"
M_INF, ALPHA = 0.80, 1.25
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)


def diagnose(mc, phi, gamma, tag, wall_s, extra):
    dz = float(np.ptp(mc.nodes[:, 2]))
    curve = wall_cp_curve(mc, phi, z=0.5 * dz, m_inf=M_INF)
    rep = shock_report(curve, M_INF)
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], phi,
                                alpha_deg=ALPHA, s_ref=dz, m_inf=M_INF)
    row = dict(case=tag, wall_s=round(wall_s, 2),
               cl_p=f["cl"], cd_p=f.get("cd"),
               cl_kj=2.0 * float(np.atleast_1d(gamma)[0]),
               gamma=float(np.atleast_1d(gamma)[0]),
               x_shock_up=rep["upper"].get("x_shock"),
               n_cells_up=rep["upper"].get("n_cells"),
               monotone_up=rep["upper"].get("monotone"),
               x_shock_lo=rep["lower"].get("x_shock"),
               cp_min_up=rep["upper"].get("cp_min"))
    row.update(extra)
    print("  ", {k: (round(v, 5) if isinstance(v, float) else v)
                 for k, v in row.items()})
    return row


def main():
    mesh = read_mesh(MESH)
    mc, wc = cut_wake(mesh)
    print(f"mesh: {len(mc.nodes)} nodes, {len(mc.elements)} tets")
    rows = []

    # --- Newton path, sweeping the artificial-dissipation constant C -------
    for C in (1.0, 1.5, 2.0, 3.0):
        t0 = time.perf_counter()
        r = solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                   m_start=0.70, dm=0.025, dm_min=0.003,
                                   upwind_c=C, freeze_tol=1e-6,
                                   newton_kw=dict(freeze_refresh_max=8,
                                                  precond="direct",
                                                  n_newton_max=60))
        dt = time.perf_counter() - t0
        rows.append(diagnose(
            mc, r["phi"], r["gamma"], f"newton_C{C}", dt,
            dict(converged=r["converged"],
                 res_final=r["residual_history"][-1],
                 res_unfrozen=r.get("residual_unfrozen"),
                 froze=r["froze"], m_max=float(np.sqrt(r["mach2_max"])),
                 n_limited=r["n_limited"], n_floored=r["n_floored"],
                 upwind_c=C, m_crit=0.95)))

    # --- Newton path, sweeping the switch threshold M_crit ----------------
    for mc_switch in (0.90, 1.00):
        t0 = time.perf_counter()
        r = solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                   m_start=0.70, dm=0.025, dm_min=0.003,
                                   upwind_c=1.5, m_crit=mc_switch,
                                   freeze_tol=1e-6,
                                   newton_kw=dict(freeze_refresh_max=8,
                                                  precond="direct",
                                                  n_newton_max=60))
        dt = time.perf_counter() - t0
        rows.append(diagnose(
            mc, r["phi"], r["gamma"], f"newton_Mc{mc_switch}", dt,
            dict(converged=r["converged"],
                 res_final=r["residual_history"][-1],
                 res_unfrozen=r.get("residual_unfrozen"),
                 froze=r["froze"], m_max=float(np.sqrt(r["mach2_max"])),
                 n_limited=r["n_limited"], n_floored=r["n_floored"],
                 upwind_c=1.5, m_crit=mc_switch)))

    # --- Picard path (the recipe that produced the committed G4.1 gate) ---
    t0 = time.perf_counter()
    r = solve_transonic_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                max_gamma_evals=12, n_picard_eval=800)
    dt = time.perf_counter() - t0
    rows.append(diagnose(
        mc, r["phi"], r["gamma"], "picard_G4.1_recipe", dt,
        dict(converged=r["converged"], res_final=None,
             res_unfrozen=None, froze=False,
             m_max=float(np.sqrt(r["mach2_max"])),
             n_limited=r["n_limited"], n_floored=r["n_floored"],
             upwind_c=1.5, m_crit=0.95,
             kutta_mismatch=r["kutta_mismatch"],
             n_picard_total=r["n_picard_total"])))

    keys = sorted({k for row in rows for k in row})
    with open(OUT / "exp1_crosscheck.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["case"] + [k for k in keys
                                                      if k != "case"])
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT / "exp1_crosscheck.csv")


if __name__ == "__main__":
    main()
