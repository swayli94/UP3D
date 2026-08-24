"""Part 1 of the sigma-freeze round: REPRODUCE the limit cycle before analysing it.

Pre-registered in phases/p3/docs/dev_phase_three/20260812-1100-sigma-freeze-prereg.md, committed before this
file was written.

Why reproduction has to come first. The user asked why the limit cycle exists. The entry check found
that there are TWO diseases with OPPOSITE treatments:

  P-A "frozen too early": sigma pinned while the field is still in a large transient (sigma_min
      swinging 0.02-0.99). Measured last round. Wants refreshing to CONTINUE.
  P-B "tail selection churn": n_shock 73 <-> 74 with sigma_min moving in the SEVENTH digit and the
      residual stalled at ~5e-6. Recorded ONLY in a code comment -- never reproduced here. Wants
      refreshing to STOP.

`_SIGMA_REFRESH_MAX = 8` is a compromise that guarantees neither. And analysing a phenomenon one has
only read about in prose is analysing the prose -- so this script tries to make P-B happen, at the
exact condition the comment records, by removing the refresh cap.

★ No library edit for that: the cap is a module constant read at RUNTIME in the driver loop, so it is
monkeypatched in this process only.

★ L1 is defined so that "it looks like it oscillates" cannot pass -- four conditions together, and
the PERIOD must be reported, because the project has already published a conclusion off a period-3
cycle whose descent10 came out 2.0021 and fell on the wrong side of a hand-picked threshold.

Outputs (TRACKED): bench/gate_results/task3_sigma_freeze_part1.csv
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
from pyfp3d.solve import newton as _newton                      # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting            # noqa: E402
from run_task3_refinement_paradox import build                  # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_sigma_freeze_part1.csv")
RECIPE = dict(m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80)
#: (tag, mesh, m_inf, alpha, upwind_c) -- the FIRST leg is the condition the comment records for
#: P-B, reproduced verbatim; the others are the conditions already in hand this phase.
LEGS = (("comment_M07875", "cases/meshes/naca0012_2.5d/coarse.msh", 0.7875, 1.25, 1.5),
        ("unstr_coarse_M080", "cases/meshes/naca0012_2.5d/coarse.msh", 0.80, 1.25, 1.5),
        ("hybrid_R0_M080", None, 0.80, 1.25, 1.5))
SEEDS = (0, 5, 12)
TAIL = 10          # the "tail" window, fixed in advance
MAX_DISTINCT = 3   # L1: at most this many distinct n_shock values recur in the tail
LEG_GATE_S = 600.0


def classify_tail(hist, res):
    """L1 / L2 / L3 on ONE leg, by the pre-registered four-part test.

    hist: the per-refresh (sigma_min, n_shock, m1_max, sigma_delta, converged) records.
    res:  the residual history.
    Returns (label, detail dict). "It looks like it oscillates" cannot pass this.
    """
    if len(hist) < TAIL + 2:
        return "L3", dict(why=f"only {len(hist)} refreshes -- tail window ({TAIL}) not available")
    tail = hist[-TAIL:]
    nsh = [int(h[1]) for h in tail]
    dsig = [float(h[3]) for h in tail]
    distinct = sorted(set(nsh))
    #: (i) small recurring value set for the post-shock COUNT
    c_set = len(distinct) <= MAX_DISTINCT and len(nsh) > len(distinct)
    #: (ii) the set actually keeps changing (a constant count is not churn)
    c_move = len(distinct) > 1
    #: (iii) max|dsigma| does NOT go to zero
    c_dsig = min(dsig) > 1e-8 and (dsig[-1] >= 0.1 * max(dsig))
    #: (iv) the residual stopped descending over the same window
    r_tail = np.asarray(res[-TAIL:], dtype=float) if len(res) >= TAIL else np.asarray(res)
    c_res = bool(len(r_tail) >= 3 and r_tail[-1] >= 0.5 * r_tail[0])
    #: the PERIOD, reported not assumed: smallest p in 1..5 with nsh[i] == nsh[i-p] throughout
    period = None
    for p in range(1, 6):
        if len(nsh) > p and all(nsh[i] == nsh[i - p] for i in range(p, len(nsh))):
            period = p
            break
    detail = dict(tail_nshock=";".join(map(str, nsh)), distinct_nshock=len(distinct),
                  tail_dsigma_min=min(dsig), tail_dsigma_max=max(dsig),
                  period=period, c_set=c_set, c_move=c_move, c_dsig=c_dsig, c_res=c_res)
    if c_set and c_move and c_dsig and c_res:
        return "L1", detail
    if not c_move and not c_dsig and c_res is False:
        return "L2", detail
    return "L3", detail


def main():
    #: ★ the cap is removed IN THIS PROCESS ONLY -- no library edit
    orig = _newton._SIGMA_REFRESH_MAX
    _newton._SIGMA_REFRESH_MAX = 10 ** 9
    print(f"Part 1: reproduce P-B.   _SIGMA_REFRESH_MAX {orig} -> {_newton._SIGMA_REFRESH_MAX} "
          f"(monkeypatched, this process only)   threads {os.environ['NUMBA_NUM_THREADS']}\n")
    rows = []
    try:
        for tag, path, m_inf, alpha, cc in LEGS:
            if path is None:
                mc, wc = build(160, 0.004, None, 1.0)[0]
            else:
                mc, wc = cut_wake(read_mesh(os.path.join(REPO, path)))
            for seed in SEEDS:
                t0 = time.perf_counter()
                r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=alpha, upwind_c=cc,
                                         n_picard_seed=seed, **RECIPE)
                wall = time.perf_counter() - t0
                hist = r.get("sigma_history") or []
                res = list(r.get("residual_history") or [])
                label, det = classify_tail(hist, res)
                ws = r["workspace"]
                sets = ws.shock_set_history
                #: how many cells actually flip between the last two recorded sets
                flips = (len(np.setxor1d(sets[-1], sets[-2])) if len(sets) >= 2 else None)
                row = dict(leg=tag, m_inf=m_inf, seed=seed, label=label,
                           converged=bool(r["converged"]), n_refresh=len(hist),
                           n_newton=len(res), res_final=(res[-1] if res else None),
                           sigma_min=r.get("sigma_min"), m1_max=r.get("m1_max"),
                           n_shock_cells=r.get("n_shock_cells"),
                           flips_last_two=flips, wall_s=round(wall, 1), **det)
                rows.append(row)
                print(f"  {tag:18} seed {seed:>2}  {label}  conv={str(r['converged']):5} "
                      f"refresh={len(hist):>3}/{len(res):>3}steps  |R|={row['res_final']:.2e}  "
                      f"n_shock tail=[{det['tail_nshock']}]  period={det['period']}  "
                      f"dsig {det['tail_dsigma_min']:.2e}..{det['tail_dsigma_max']:.2e}  "
                      f"flips={flips}  ({wall:.0f}s)", flush=True)
                if wall > LEG_GATE_S:
                    print("    ★ cost gate exceeded -- stop"); _write(rows); return 1
    finally:
        _newton._SIGMA_REFRESH_MAX = orig
    _write(rows)

    print("\n=== Part 1 reading (bands fixed in the registration section 1) ===")
    l1 = [r for r in rows if r["label"] == "L1"]
    print(f"  L1 (P-B reproduced): {len(l1)} of {len(rows)} legs")
    for r in l1:
        print(f"    {r['leg']:18} seed {r['seed']:>2}  period {r['period']}  "
              f"n_shock tail [{r['tail_nshock']}]  dsigma {r['tail_dsigma_min']:.2e}.."
              f"{r['tail_dsigma_max']:.2e}  |R| {r['res_final']:.2e}")
    if l1:
        print("  -> P-B REPRODUCED. The mechanism question now has a subject; Part 3 candidate 1")
        print("     (freeze the SET, keep the magnitude live) keeps its target.")
    else:
        print("  -> ★ L2/L3: P-B did NOT reproduce at these conditions with the cap removed.")
        print("     Per the registration that is NOT a kill, and it does NOT license hunting for")
        print("     a case elsewhere that proves it. What it does mean is recorded honestly:")
        print("     _SIGMA_REFRESH_MAX's own justification needs re-measuring, and Part 3 loses")
        print("     candidate 1's target (its premise was 'the set limit-cycles').")
        for r in rows:
            print(f"    {r['leg']:18} seed {r['seed']:>2}  {r['label']}  "
                  f"set={r['c_set']} move={r['c_move']} dsig={r['c_dsig']} res={r['c_res']}  "
                  f"distinct={r['distinct_nshock']}  refresh={r['n_refresh']}")
    return 0


def _write(rows):
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
