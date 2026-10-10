"""How many times does one streamline get charged for a shock? A literal hop-by-hop count.

Pre-registered in phases/p3/docs/dev_phase_three/20260812-0700-sigma-charge-count-prereg.md, committed before
the record-only store and this file were written.

The setting: sigma_e is the PRODUCT of per-element factors s = sigma_RH(M1) along the donor chain
starting at e. Since every factor obeys M1 <= the flow's peak Mach, sigma_accumulated >=
sigma_RH(m_max)^N -- which forces N >= 7.0 / 8.9 / 12.2 on legs whose peak Mach is only 1.33-1.40.
That part is arithmetic, not a hypothesis. What is hypothetical is WHY, and the answer left standing
(a shocked donor cycle is excluded, because such a state cannot be reported converged) is: several
DISTINCT sonic crossings lie on one acyclic chain. This script counts them.

★ The walker SHARES NO CODE with transport_sigma. That is not stylistic: the project has already
paid for it. Of three candidate termination criteria for that kernel, TWO WRONG ONES passed every
hand-built case, including five deliberately constructed cycle variants; what killed them was a
test that walked each donor chain hop by hop and compared. So the oracle here is the literal walk,
and G1 requires it to reproduce ws.sigma_frozen at rtol 1e-12 -- ★ NOT bit-equality, because the
fast algorithm multiplies the same factors in a different ORDER and multiplication is not
associative. Widening that rtol to make it pass is forbidden by the kill criterion.

Outputs (TRACKED): bench/gate_results/task3_sigma_charge_count.csv
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

from pyfp3d.kernels.entropy import total_pressure_ratio            # noqa: E402
from pyfp3d.mesh.reader import read_mesh                           # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                          # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients            # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting               # noqa: E402
from run_task3_refinement_paradox import build                     # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_sigma_charge_count.csv")
M_INF, ALPHA, C = 0.80, 1.25, 1.5
SEEDS = (0, 5, 12)
RECIPE = dict(upwind_c=C, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80)
CASES = (("hybrid_R0", None), ("unstr_coarse", "cases/meshes/naca0012_2.5d/coarse.msh"))
RTOL_ORACLE = 1e-12          # fixed by the registration -- may NOT be widened
OUTLIER_SIGMA, CONTROL_SIGMA = 0.85, 0.95      # band D1's two sides, declared in advance
LEG_GATE_S = 600.0


def walk(e, s, up, cap):
    """LITERAL hop-by-hop walk from `e` to its root: no memo, no doubling, no shared code.

    Returns (product, n_factors_below_one, n_hops, cycle_hit). The convention is the one the
    kernel's own test pins (`test_sigma_is_the_product_along_the_donor_chain`): the chain STARTS AT
    e inclusive and ends at the root, where up[root] == root.
    """
    prod, cnt, hops, cur = 1.0, 0, 0, int(e)
    seen = set()
    while True:
        prod *= s[cur]
        if s[cur] < 1.0:
            cnt += 1
        nxt = int(up[cur])
        if nxt == cur:
            return prod, cnt, hops, False
        #: ★ cycle protection AND cycle REPORTING. The registration argues no shocked cycle can
        #: survive here (such a state cannot be reported converged), but an argument is not a
        #: measurement, so a cycle is detected and reported rather than assumed away.
        if cur in seen or hops >= cap:
            return prod, cnt, hops, True
        seen.add(cur)
        cur = nxt
        hops += 1


def walk_all(s, up, cap):
    """Same literal walk for every element, with a memo of COMPLETED results only.

    The memo caches the outcome of a walk already performed; it does not borrow the kernel's
    doubling idea. Cycle-hit chains are never memoised.
    """
    n = len(s)
    memo_p = np.full(n, np.nan)
    memo_c = np.full(n, -1, dtype=np.int64)
    prods = np.empty(n)
    cnts = np.zeros(n, dtype=np.int64)
    cyc = np.zeros(n, dtype=bool)
    for e in range(n):
        path = []
        cur, prod, cnt, hops, hit = int(e), 1.0, 0, 0, False
        while True:
            if memo_c[cur] >= 0:
                prod *= memo_p[cur]; cnt += int(memo_c[cur]); break
            path.append(cur)
            prod *= s[cur]
            if s[cur] < 1.0:
                cnt += 1
            nxt = int(up[cur])
            if nxt == cur:
                break
            if hops >= cap or nxt in path:
                hit = True; break
            cur = nxt; hops += 1
        prods[e], cnts[e], cyc[e] = prod, cnt, hit
        if not hit:
            #: memoise only the head (safe and enough -- the tail is memoised on its own turn)
            memo_p[e], memo_c[e] = prods[e], cnts[e]
    return prods, cnts, cyc


def main():
    print(f"charge-count diagnostic   M{M_INF}/alpha {ALPHA}/C {C}   seeds {SEEDS}   "
          f"threads {os.environ['NUMBA_NUM_THREADS']}   (production default -- no instrument)\n")
    rows = []
    for tag, path in CASES:
        if path is None:
            mc, wc = build(160, 0.004, None, 1.0)[0]
        else:
            mc, wc = cut_wake(read_mesh(os.path.join(REPO, path)))
        for seed in SEEDS:
            t0 = time.perf_counter()
            r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                     n_picard_seed=seed, **RECIPE)
            wall = time.perf_counter() - t0
            ws = r["workspace"]
            f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"],
                                        r["phi"], alpha_deg=ALPHA,
                                        s_ref=float(np.ptp(mc.nodes[:, 2])), m_inf=M_INF)
            row = dict(case=tag, seed=seed, converged=bool(r["converged"]),
                       cl_p=float(f["cl"]), sigma_min=r.get("sigma_min"),
                       m1_max=r.get("m1_max"), n_shock_cells=r.get("n_shock_cells"),
                       sigma_transport_converged=r.get("sigma_transport_converged"),
                       n_limited=int(r.get("n_limited") or 0),
                       n_floored=int(r.get("n_floored") or 0),
                       #: ★ the decisive extra field. _SIGMA_REFRESH_MAX = 8: sigma is refreshed
                       #: for the first 8 Newton steps and then HELD, deliberately (the post-shock
                       #: SET is a discrete selection that limit-cycles, so freezing it makes the
                       #: tail a fixed smooth system Newton can finish). The side effect this round
                       #: measures: the held value is whatever the flow looked like AT THE FREEZE
                       #: MOMENT, and where the iteration is at step 8 depends on the seed.
                       n_sigma_refresh=r.get("n_sigma_refresh"),
                       n_newton=len(r.get("residual_history") or []),
                       #: the m1_max trajectory over the refreshes, so "frozen at a transient" is
                       #: read off the history rather than inferred from one endpoint
                       m1_max_traj=";".join(
                           f"{h[2]:.4f}" for h in (r.get("sigma_history") or [])),
                       sigma_min_traj=";".join(
                           f"{h[0]:.5f}" for h in (r.get("sigma_history") or [])),
                       wall_s=round(wall, 1))

            s, up, sig = ws.ent._s, ws.upstream_sigma, ws.sigma_frozen
            if s is None or up is None or sig is None:
                print(f"  {tag:13} seed {seed:>2}  ★ PREMISE FAIL: sigma arrays unreachable "
                      f"(s={s is not None}, up={up is not None}, sigma={sig is not None})")
                row.update(oracle_max_rel=None, n_at_min=None, n_max=None)
                rows.append(row); continue

            n = len(s)
            prods, cnts, cyc = walk_all(s, up, cap=n)
            #: G1 oracle -- the walk's product against the kernel's accumulated sigma
            rel = np.abs(prods - sig) / np.maximum(np.abs(sig), 1e-300)
            imin = int(np.argmin(sig))
            p_i, c_i, h_i, cyc_i = walk(imin, s, up, cap=n)     # uncached, for the binding element
            row.update(n_elements=n, n_roots=int(np.count_nonzero(up == np.arange(n))),
                       oracle_max_rel=float(rel.max()),
                       n_at_min=int(c_i), n_at_min_hops=int(h_i),
                       n_max=int(cnts.max()), n_cycles=int(np.count_nonzero(cyc)),
                       cycle_at_min=bool(cyc_i),
                       n_charged=int(np.count_nonzero(s < 1.0)),
                       chain_len_p95=float(np.percentile(
                           [walk(e, s, up, n)[2] for e in range(0, n, max(1, n // 400))], 95)))
            rows.append(row)
            print(f"  {tag:13} seed {seed:>2}  conv={str(row['converged']):5} "
                  f"sig_min={row['sigma_min']:.5f} m1_max={row['m1_max']:.4f}  "
                  f"charged={row['n_charged']:>4}  ★ N_at_min={row['n_at_min']:>3} "
                  f"N_max={row['n_max']:>3}  oracle_rel={row['oracle_max_rel']:.2e} "
                  f"cycles={row['n_cycles']}  refresh={row['n_sigma_refresh']}"
                  f"/{row['n_newton']}steps  ({wall:.0f}s)", flush=True)
            if wall > LEG_GATE_S:
                print("    ★ cost gate exceeded -- stop"); _write(rows); return 1
    _write(rows)
    return _read(rows)


def _read(rows):
    good = [r for r in rows if r["converged"] and r.get("oracle_max_rel") is not None]

    print("\n=== G1 oracle: the literal walk vs the kernel's accumulated sigma ===")
    worst = max((r["oracle_max_rel"] for r in good), default=float("inf"))
    print(f"  worst relative difference over all elements, all legs: {worst:.3e} "
          f"(threshold {RTOL_ORACLE:.0e}, FIXED by the registration)")
    if not (worst <= RTOL_ORACLE):
        print("  -> ★ G1 FAIL. Kill clause 1: every count below is VOID, and the rtol may NOT be")
        print("     widened to make this pass. Fix the tool, not the threshold.")
        return 1
    print("  -> PASS")

    print("\n=== N2 consistency assert (my algebra vs my walker) ===")
    ok2 = True
    for r in good:
        sm, m1 = r["sigma_min"], r["m1_max"]
        s1 = total_pressure_ratio(m1) if m1 and m1 > 1.0 else 1.0
        bound = (np.log(sm) / np.log(s1)) if (0 < sm < 1 and s1 < 1) else 0.0
        hit = r["n_at_min"] >= np.floor(bound)
        ok2 &= hit
        print(f"  {r['case']:13} seed {r['seed']:>2}  sigma_min {sm:.5f}  m1_max {m1:.4f}  "
              f"bound N >= {bound:5.1f}   measured N_at_min = {r['n_at_min']:>3}  "
              f"-> {'ok' if hit else '★ BELOW BOUND -- tool or algebra is wrong'}")
    if not ok2:
        print("  -> ★ a measured count below the arithmetic bound is a TOOL failure, not a "
              "physics finding. Stopping.")
        return 1

    print("\n=== RECORDED: the sigma FREEZE moment (_SIGMA_REFRESH_MAX = 8) ===")
    for r in good:
        print(f"  {r['case']:13} seed {r['seed']:>2}  refreshes {r['n_sigma_refresh']}"
              f" of {r['n_newton']} Newton steps   m1_max over refreshes: "
              f"{r['m1_max_traj']}")
    print("  ★ sigma is HELD after the 8th refresh by design. If m1_max at the LAST refresh")
    print("    exceeds the converged field's peak Mach, the correction riding the final answer")
    print("    was computed from a TRANSIENT -- and where the iteration sits at step 8 depends")
    print("    on the seed. That is a path dependence, read here, not inferred.")

    print("\n=== D1 / D2 / D3 (registration section 4; N_at_min is binding) ===")
    print(f"  {'case':13}{'seed':>5}{'sigma_min':>11}{'N_at_min':>10}{'N_max':>7}"
          f"{'charged':>9}{'cl_p':>10}")
    for r in good:
        print(f"  {r['case']:13}{r['seed']:>5}{r['sigma_min']:>11.5f}{r['n_at_min']:>10}"
              f"{r['n_max']:>7}{r['n_charged']:>9}{r['cl_p']:>10.6f}")
    out = [r for r in good if r["sigma_min"] <= OUTLIER_SIGMA]
    ctl = [r for r in good if r["sigma_min"] >= CONTROL_SIGMA]
    print(f"\n  outliers (sigma_min <= {OUTLIER_SIGMA}): "
          f"{[(r['case'], r['seed'], r['n_at_min']) for r in out]}")
    print(f"  controls (sigma_min >= {CONTROL_SIGMA}): "
          f"{[(r['case'], r['seed'], r['n_at_min']) for r in ctl]}")
    if out and ctl and min(r["n_at_min"] for r in out) >= 5 \
            and max(r["n_at_min"] for r in ctl) <= 2:
        print("\n  -> D1: MULTIPLICITY is the discriminator. The target is the CHARGING RULE")
        print("     (candidates: charge a shock once per connected region; or bound the")
        print("     accumulated loss along a chain physically). Register separately.")
    elif good and (max(r["n_at_min"] for r in good) - min(r["n_at_min"] for r in good)) <= 1:
        print("\n  -> D2: counts agree within 1 while sigma_min does not. It is PER-CHARGE")
        print("     MAGNITUDE, not count -- the target becomes the knee walk's choice of M1.")
    else:
        print("\n  -> D3: neither band cleanly. RECORDED, no direction claimed.")
    return 0


def _write(rows):
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
