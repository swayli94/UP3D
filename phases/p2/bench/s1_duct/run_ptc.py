"""GS1.2b: can pseudo-transient continuation replace the Mach ramp?

Pre-registered in phases/p2/docs/dev_phase_two/20260728-2030-s1-ptc.md.

The ramp was measured to fail under refinement (GS1.2 Q2): at h_wall = 0.005c it
cannot reach M0.7875/alpha1.25 (dies at M0.741) nor M0.80/alpha0 (dies at
M0.784), while 0.02c and 0.01c reach both. PTC (SER schedule, new
`ptc_lambda0`) solves AT the target Mach from a cold start instead, so there is
no continuation path to fall off.

Legs, in increasing order of what they prove:

  1. calibration (coarse): sweep ptc_lambda0 to find the working range, and
     check the PTC answer equals the ramp answer at the same condition.
  2. medium: same check (the ramp works here, so this is an equivalence test).
  3. fine: the acceptance leg -- the ramp CANNOT get here; PTC must.

Outputs: results/gs1_2b_ptc.csv
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
from pyfp3d.solve.newton import (solve_newton_lifting,         # noqa: E402
                                 solve_newton_transonic)

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

LAMBDAS = (0.1, 1.0, 10.0)
CASES = [("M0.7875_a1.25", 0.7875, 1.25), ("M0.80_a0", 0.80, 0.0)]


def post(mc, r, m_inf, alpha, dz):
    rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=m_inf),
                       m_inf)
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], r["phi"],
                                alpha_deg=alpha, s_ref=dz, m_inf=m_inf)
    return rep, f


def row_from(tag, level, mc, wc, m_inf, alpha, dz, r, wall, kind, lam=None):
    rep, f = post(mc, r, m_inf, alpha, dz)
    return dict(case=tag, level=level, kind=kind, ptc_lambda0=lam,
                n_dof=len(mc.nodes), converged=r["converged"],
                accept_reason=r.get("accept_reason"),
                n_newton=r.get("n_newton"),
                wall_s=round(wall, 1),
                cl_p=round(f["cl"], 6),
                gamma=round(float(r["gamma"][0]), 8),
                x_shock=rep["upper"].get("x_shock"),
                m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                n_limited=r["n_limited"], n_floored=r["n_floored"],
                res_final=r["residual_history"][-1],
                ptc_lambda_last=(r.get("ptc_history") or [None])[-1])


def main():
    rows = []
    levels = [("coarse", LAMBDAS), ("medium", (1.0,)), ("fine", (1.0,))]
    for level, lams in levels:
        path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
        if not path.exists():
            print(f"skip {level}")
            continue
        mc, wc = cut_wake(read_mesh(path))
        dz = float(np.ptp(mc.nodes[:, 2]))
        print(f"\n=== {level}: {len(mc.nodes)} nodes ===", flush=True)
        for tag, m_inf, alpha in CASES:
            # --- reference: the shipped Mach ramp -------------------------
            t0 = time.perf_counter()
            try:
                rr = solve_newton_transonic(
                    mc, wc, m_inf=m_inf, alpha_deg=alpha, m_start=0.70,
                    dm=0.025, dm_min=0.003, upwind_c=1.5, freeze_tol=1e-6,
                    newton_kw=dict(freeze_refresh_max=8, precond="direct",
                                   direct_refactor_every=4, n_newton_max=60))
                done = [lr["m"] for lr in rr["level_results"]
                        if lr["converged"]]
                row = row_from(tag, level, mc, wc, m_inf, alpha, dz, rr,
                               time.perf_counter() - t0, "ramp")
                row["m_reached"] = max(done) if done else None
            except Exception as exc:                           # noqa: BLE001
                row = dict(case=tag, level=level, kind="ramp",
                           converged=False, error=type(exc).__name__)
            print("   ramp ", {k: row.get(k) for k in
                               ("case", "converged", "m_reached", "cl_p",
                                "x_shock", "wall_s")}, flush=True)
            rows.append(row)

            # --- PTC at the target Mach, cold start ----------------------
            for lam in lams:
                t0 = time.perf_counter()
                try:
                    rp = solve_newton_lifting(
                        mc, wc, m_inf=m_inf, alpha_deg=alpha, upwind_c=1.5,
                        m_crit=0.95, freeze_tol=1e-6, ptc_lambda0=lam,
                        freeze_refresh_max=8, precond="direct",
                        direct_refactor_every=4, n_newton_max=200)
                    row = row_from(tag, level, mc, wc, m_inf, alpha, dz, rp,
                                   time.perf_counter() - t0, "ptc", lam)
                    row["m_reached"] = m_inf if rp["converged"] else None
                except Exception as exc:                       # noqa: BLE001
                    row = dict(case=tag, level=level, kind="ptc",
                               ptc_lambda0=lam, converged=False,
                               error=f"{type(exc).__name__}: {exc}"[:120])
                print(f"   ptc l={lam:<5g}", {k: row.get(k) for k in
                                              ("converged", "n_newton",
                                               "cl_p", "x_shock", "m_max",
                                               "wall_s", "res_final")},
                      flush=True)
                rows.append(row)

    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(OUT / "gs1_2b_ptc.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1_2b_ptc.csv")

    print("\n=== headline: ramp vs PTC ===")
    for level in ("coarse", "medium", "fine"):
        for tag, _, _ in CASES:
            sub = [r for r in rows if r.get("level") == level
                   and r.get("case") == tag]
            if not sub:
                continue
            print(f"  {level:7s} {tag:14s}")
            for r in sub:
                lam = "" if r.get("ptc_lambda0") is None \
                    else f" l={r['ptc_lambda0']:g}"
                print(f"     {r['kind']:4s}{lam:8s} conv="
                      f"{str(r.get('converged')):5s} "
                      f"m_reached={r.get('m_reached')} "
                      f"cl_p={r.get('cl_p')} x_shock={r.get('x_shock')} "
                      f"wall={r.get('wall_s')}s")


if __name__ == "__main__":
    main()
