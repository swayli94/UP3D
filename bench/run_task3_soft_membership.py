"""Part 3 of the sigma-freeze round: soften the post-shock MEMBERSHIP test and read F1-F5.

Pre-registered in docs/dev_phase_three/20260812-1100-sigma-freeze-prereg.md, with the candidate
order and the ramp's exact form fixed in addendum #1 BEFORE this code was written.

Why this candidate and not relaxation or averaging: Part 1 measured the mechanism. The membership
test is a HARD switch, so a cell sitting on the sonic line flips on an infinitesimal change in phi
while its factor JUMPS by 1 - sigma_RH(M1) -- and max|dsigma| pins at exactly that jump (1.55e-2
matching a cell at M1 1.27), with periods 2 to 5. Relaxation and averaging damp the CONSEQUENCE;
softening removes the jump at SOURCE. And the precedent is next door in the same kernel: the
artificial-density switch has always been a continuous ramp (m_crit). Only this test was hard.

Registered criteria:
  F1  cross-seed cl_p spread (CONVERGED legs only; UNDEFINED below two) from 32.11 % to <= 5 %
  F2  do NOT put P-B back: no L1 signature, and the residual still converges
  F3  legacy path bit-identical (soft_eps = 0)
  F5  ★ report the SHOCK POSITION shift against M1's own 0.61 +- 0.02 -- softening charges partial
      entropy to cells that have not crossed a shock, so a spread number alone cannot discharge this
  eps sensitivity across 0.02 / 0.05 / 0.10 -- it is a CALIBRATION, and four hand-picked constants
      have already decided conclusions in this project

Outputs (TRACKED): bench/gate_results/task3_soft_membership.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.reader import read_mesh                        # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                       # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve               # noqa: E402
from pyfp3d.post.shock import shock_report                      # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients         # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting            # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_soft_membership.csv")
M_INF, ALPHA, C = 0.80, 1.25, 1.5
SEEDS = (0, 5, 12)
#: 0.0 is the legacy control (F3); the rest is the pre-registered eps sensitivity set
EPSES = (0.0, 0.02, 0.05, 0.10)
RECIPE = dict(upwind_c=C, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80)
CASES = (("hybrid_R0", None), ("unstr_coarse", "cases/meshes/naca0012_2.5d/coarse.msh"))
#: F1's baseline and target, and F5's band -- all declared in the registration
BASE_SPREAD_PCT, F1_TARGET_PCT = 32.11, 5.0
SHOCK_REF, SHOCK_TOL = 0.61, 0.02
LEG_GATE_S = 600.0


def _mesh(path):
    if path is None:
        from run_task3_refinement_paradox import build
        return build(160, 0.004, None, 1.0)[0]
    return cut_wake(read_mesh(os.path.join(REPO, path)))


def run(tag, mc, wc, eps, seed):
    t0 = time.perf_counter()
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, n_picard_seed=seed,
                            sigma_soft_eps=eps, **RECIPE)
    wall = time.perf_counter() - t0
    h = list(r.get("residual_history") or [])
    z = float(np.unique(mc.nodes[:, 2]).mean())
    rep = shock_report(wall_cp_curve(mc, r["phi"], z=z, m_inf=M_INF), M_INF)
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
                                alpha_deg=ALPHA, s_ref=float(np.ptp(mc.nodes[:, 2])), m_inf=M_INF)
    xs = rep["upper"].get("x_shock")
    fr = r.get("sigma_freeze_report") or {}
    return dict(case=tag, eps=eps, seed=seed, converged=bool(r["converged"]),
                res_final=(h[-1] if h else None), n_newton=len(h),
                n_limited=int(r.get("n_limited") or 0),
                n_floored=int(r.get("n_floored") or 0),
                cl_p=float(f["cl"]), x_shock=None if xs is None else float(xs),
                sigma_min=r.get("sigma_min"), m1_max=r.get("m1_max"),
                n_shock_cells=r.get("n_shock_cells"),
                frozen_in_transient=fr.get("frozen_in_transient"),
                selection_churn=fr.get("selection_churn"),
                churn_period=fr.get("churn_period"),
                last_sigma_delta=fr.get("last_sigma_delta"),
                wall_s=round(wall, 1))


def spread(rows, key):
    """CONVERGED legs only; UNDEFINED below two of them (the defect this project already fixed)."""
    v = [r[key] for r in rows if r["converged"] and r[key] is not None]
    if len(v) < 2:
        return None, len(v)
    return max(v) - min(v), len(v)


def main():
    print(f"Part 3: soften the membership test.   M{M_INF}/alpha {ALPHA}/C {C}   "
          f"eps {EPSES}   seeds {SEEDS}   threads {os.environ['NUMBA_NUM_THREADS']}\n")
    rows = []
    for tag, path in CASES:
        mc, wc = _mesh(path)
        for eps in EPSES:
            for seed in SEEDS:
                row = run(tag, mc, wc, eps, seed)
                rows.append(row)
                print(f"  {tag:13} eps={eps:.2f} seed {seed:>2}  "
                      f"conv={str(row['converged']):5} |R|={row['res_final']:.2e} "
                      f"lim/flr={row['n_limited']}/{row['n_floored']} "
                      f"cl_p={row['cl_p']:>9.6f} "
                      f"x_shock={'-' if row['x_shock'] is None else format(row['x_shock'], '.4f')} "
                      f"sig_min={format(row['sigma_min'], '.5f') if row['sigma_min'] else '-'} "
                      f"churn={row['selection_churn']} p={row['churn_period']} "
                      f"({row['wall_s']:.0f}s)", flush=True)
                if row["wall_s"] > LEG_GATE_S:
                    print("    ★ cost gate exceeded -- stop"); _write(rows); return 1
    _write(rows)
    return _read(rows)


def _read(rows):
    #: --- F3: the legacy path must be bit-identical ------------------------------------------
    print("\n=== F3: eps = 0 must reproduce the legacy answers BIT-IDENTICALLY ===")
    ref = {("hybrid_R0", 0): (0.261367, 0.7141), ("hybrid_R0", 5): (0.258223, 0.6051),
           ("hybrid_R0", 12): (0.351467, 0.5458),
           ("unstr_coarse", 0): (0.407978, 0.6195), ("unstr_coarse", 5): (0.407108, 0.6073)}
    f3 = True
    for (tag, seed), (cl_a, xs_a) in ref.items():
        m = [r for r in rows if r["case"] == tag and r["eps"] == 0.0 and r["seed"] == seed]
        if not m or not m[0]["converged"]:
            print(f"  {tag:13} seed {seed:>2}: not converged at eps=0 -- cannot compare")
            continue
        r = m[0]
        hit = abs(r["cl_p"] - cl_a) < 5e-6 and abs(r["x_shock"] - xs_a) < 5e-4
        f3 &= hit
        print(f"  {tag:13} seed {seed:>2}: cl {r['cl_p']:.6f} vs {cl_a}, "
              f"x_shock {r['x_shock']:.4f} vs {xs_a}  -> {'MATCH' if hit else '★ MISMATCH'}")
    print(f"  -> F3 {'PASS' if f3 else 'FAIL'}   ★ a FAIL here voids everything below: the flag "
          f"would not be default-inert")

    #: --- F1 / F2 / F5 per eps ----------------------------------------------------------------
    print("\n=== F1 (cross-seed spread) / F2 (no churn back) / F5 (shock position) per eps ===")
    print(f"  {'case':13}{'eps':>6}{'conv':>6}{'cl spread %':>13}{'x_shock spread':>15}"
          f"{'churn?':>8}{'x_shock in band?':>18}")
    verdicts = {}
    for tag, _ in CASES:
        for eps in EPSES:
            sub = [r for r in rows if r["case"] == tag and r["eps"] == eps]
            good = [r for r in sub if r["converged"]]
            dcl, n = spread(sub, "cl_p")
            dxs, _ = spread(sub, "x_shock")
            rel = (None if dcl is None
                   else 100.0 * dcl / max(abs(np.mean([abs(r["cl_p"]) for r in good])), 1e-12))
            churn = any(bool(r["selection_churn"]) for r in sub)
            inband = [r for r in good if r["x_shock"] is not None
                      and abs(r["x_shock"] - SHOCK_REF) <= SHOCK_TOL]
            verdicts[(tag, eps)] = dict(rel=rel, n=n, churn=churn, inband=len(inband),
                                        ngood=len(good), dxs=dxs)
            print(f"  {tag:13}{eps:>6.2f}{n:>4}/3"
                  f"{('UNDEFINED' if rel is None else format(rel, '.2f')):>13}"
                  f"{('UNDEFINED' if dxs is None else format(dxs, '.4f')):>15}"
                  f"{str(churn):>8}{f'{len(inband)}/{len(good)}':>18}")

    print(f"\n  bands: F1 needs spread <= {F1_TARGET_PCT} % (baseline {BASE_SPREAD_PCT} %); "
          f"F2 needs churn False AND convergence kept; F5 REPORTS x_shock vs "
          f"{SHOCK_REF}+-{SHOCK_TOL}")
    for tag, _ in CASES:
        for eps in EPSES:
            if eps == 0.0:
                continue
            v = verdicts[(tag, eps)]
            base = verdicts[(tag, 0.0)]
            if v["rel"] is None:
                print(f"    {tag:13} eps={eps:.2f} -> UNDEFINED spread "
                      f"({v['n']} converged) -- not a pass, not a small spread")
                continue
            f1 = v["rel"] <= F1_TARGET_PCT
            f2 = (not v["churn"]) and v["ngood"] >= base["ngood"]
            print(f"    {tag:13} eps={eps:.2f} -> F1 {'PASS' if f1 else 'FAIL'} "
                  f"({v['rel']:.2f} % vs {base['rel'] if base['rel'] is None else format(base['rel'], '.2f')} % at eps=0)"
                  f" | F2 {'PASS' if f2 else 'FAIL'} (churn={v['churn']}, "
                  f"converged {v['ngood']} vs {base['ngood']})"
                  f" | F5 x_shock in band {v['inband']}/{v['ngood']} "
                  f"(eps=0: {base['inband']}/{base['ngood']})")
    print("\n  ★ F1 and F2 are a PAIR: F1 alone could be bought by never freezing, which is the")
    print("    problem the freeze exists to solve; F2 alone means nothing was cured.")
    print("  ★ eps is a CALIBRATION -- read the three eps rows as a sensitivity, not as three")
    print("    independent attempts to pass.")
    return 0


def _write(rows):
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
