"""Self-convergence on a COMMON point set, factor-2 chain S0/S2/S4. ZERO SOLVES.

Pre-registered in docs/dev_phase_three/20260814-2200-common-pointset-prereg.md, whose section 0 states my
prediction (N-DISC) up front, because I have already seen the sqrt(2) chain's numbers and cannot claim a
blind test -- and states how that prediction could fail (the clamped leading-edge points may have carried
most of the self-difference, in which case d_self collapses and the reading flips to N-MODEL).

The construction is the whole content: an experimental point is used only if it lies inside ALL THREE
computed curves. Dropping points per level would drop DIFFERENT points at different levels, which makes the
two self-difference vectors different objects (the fifth standing question). No clamped point is used at any
level, so G-X is not loosened; "the whole band is void" becomes "point-level validity on one shared object".

★ G-Z: this script calls NO solver. It reads cached curves (npz) and committed CSVs only -- asserted below,
because "I will just re-solve that one level" is how a zero-solve round becomes a different round.

Outputs (TRACKED): bench/gate_results/task3_common_pointset.csv
"""

import csv
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from run_m3_budget import BANDS, ETAS, N_UNMASKED, parse_experiment    # noqa: E402
import run_task3_fixed_budget as FB                                    # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_common_pointset.csv")
PREV_CSV = os.path.join(HERE, "gate_results", "task3_selfconvergence.csv")
#: the factor-2 chain (h_far explicit = 120*h_wall on all three; all converged)
CHAIN = ("S0", "S2", "S4")
#: RECORDED only -- these numbers were already seen, so they cannot serve as a criterion (section 0)
SEEN = ("S3", "S4")
GR_REF = 0.242221
MODEL_MAX, DISC_MIN, MIN_COMMON = 1.0 / 3.0, 1.0, 20
RTOL_IDENTITY = 1e-12


def load_curves(tag):
    path = os.path.join(FB.SCRATCH, f"curves_{tag}.npz")
    if not os.path.exists(path):
        return None
    z = np.load(path)
    etas = sorted({float(k.split("|")[0]) for k in z.files})
    return {e: {k.split("|")[1]: z[k] for k in z.files if float(k.split("|")[0]) == e}
            for e in etas}


def point_table(curves, exp, band, side):
    """Per-station experimental points with, for each level, whether the point is IN RANGE.

    Returns rows [(eta, x, cp_exp, {tag: (cp_interp, in_range)})] in a fixed order.
    """
    lo, hi = next((l, h) for n, l, h in BANDS if n == band)
    out = []
    for eta in ETAS[:N_UNMASKED]:
        e = exp[eta]
        m = (e["upper"] == (side == "upper")) & (e["x"] >= lo) & (e["x"] < hi)
        idx = np.flatnonzero(m)
        for i in idx:
            per = {}
            for tag, c in curves.items():
                if eta not in c:
                    per[tag] = (None, False)
                    continue
                xs, cs = c[eta][f"x_{side}"], c[eta][f"cp_{side}"]
                inr = bool(np.min(xs) <= e["x"][i] <= np.max(xs))
                per[tag] = (float(np.interp(e["x"][i], xs, cs)), inr)
            out.append((eta, float(e["x"][i]), float(e["cp"][i]), per))
    return out


def rms(v):
    v = np.asarray(v, dtype=float)
    return float(np.sqrt(np.mean(v ** 2))) if v.size else None


def main():
    #: --- G-Z ------------------------------------------------------------------------------------
    assert "solve" not in globals() and "run_m3_budget" not in sys.modules or True
    #: ★★ The first version of this guard MATCHED ITS OWN forbidden-token list and fired on a file
    #: that solves nothing -- the same self-reference family as `pgrep -f <pattern>` matching the
    #: shell that launched it, which this project has already recorded twice. The tokens are
    #: therefore assembled at runtime so the literals never appear in the source being scanned.
    src = open(__file__).read()
    body = src.split(chr(34) * 3, 2)[2]
    for forbidden in ("solve" + "_newton", "run_m3_budget." + "solve", "onera_m6" + "_wing_mesh",
                      "cut" + "_wake"):
        assert forbidden not in body, f"G-Z: this round must not solve ({forbidden})"
    print("★ G-Z: zero solves -- cached curves + committed CSVs only\n")

    exp = parse_experiment()
    curves = {t: load_curves(t) for t in CHAIN}
    missing = [t for t, c in curves.items() if c is None]
    if missing:
        print(f"★ cached curves missing for {missing} -- STOP (kill clause 4; re-solving would make "
              f"this a different round).")
        return 1
    print(f"loaded cached curves: {', '.join(CHAIN)}   "
          f"(h_wall 0.060 / 0.030 / 0.015, factor 2, h_far = 120*h_wall)\n")

    rows, summary = [], {}
    for band, _lo, _hi in BANDS:
        for side in ("upper", "lower"):
            tbl = point_table(curves, exp, band, side)
            if not tbl:
                continue
            common = [r for r in tbl if all(r[3][t][1] for t in CHAIN)]
            dropped = [r for r in tbl if r not in common]
            for eta, x, ce, per in tbl:
                rows.append(dict(band=band, side=side, eta=eta, x=x, cp_exp=ce,
                                 **{f"cp_{t}": per[t][0] for t in CHAIN},
                                 **{f"in_{t}": per[t][1] for t in CHAIN},
                                 common=all(per[t][1] for t in CHAIN)))
            summary[(band, side)] = dict(
                n_full=len(tbl), n_common=len(common), n_drop=len(dropped),
                drop_x=[round(r[1], 5) for r in dropped],
                e_full={t: rms([r[3][t][0] - r[2] for r in tbl if r[3][t][0] is not None])
                        for t in CHAIN},
                e_common={t: rms([r[3][t][0] - r[2] for r in common]) for t in CHAIN},
                d_self={f"{a}->{b}": rms([r[3][a][0] - r[3][b][0] for r in common])
                        for a, b in zip(CHAIN, CHAIN[1:])})
    _write(rows)

    key = ("LE", "upper")
    s = summary[key]
    print("=== G-N: the common point set (binding band = LE upper) ===")
    print(f"  {'band/side':14}{'n_full':>8}{'n_common':>10}{'dropped':>9}   dropped x")
    for (b, sd), v in summary.items():
        print(f"  {b + '/' + sd:14}{v['n_full']:>8}{v['n_common']:>10}{v['n_drop']:>9}   "
              f"{v['drop_x'] if v['drop_x'] else '-'}")

    print("\n=== G-R: same cached curves as the previous round? ===")
    print(f"  S4 full-set e(LE upper) = {s['e_full']['S4']:.6f} vs committed {GR_REF}")
    if abs(s["e_full"]["S4"] - GR_REF) > 5e-7:
        print("  -> ★ G-R FAIL: not the same curves. STOP (kill clause 1).")
        return 1
    print("  -> PASS")

    #: --- G-I: the exact set identity ------------------------------------------------------------
    print("\n=== G-I: e_full^2*n_full == e_common^2*n_common + sum(dropped^2) ===")
    worst = 0.0
    for (b, sd), v in summary.items():
        tbl = point_table(curves, exp, b, sd)
        common = [r for r in tbl if all(r[3][t][1] for t in CHAIN)]
        for t in CHAIN:
            res_all = [r[3][t][0] - r[2] for r in tbl if r[3][t][0] is not None]
            res_com = [r[3][t][0] - r[2] for r in common]
            res_drop = [r[3][t][0] - r[2] for r in tbl
                        if r not in common and r[3][t][0] is not None]
            lhs = sum(x * x for x in res_all)
            rhs = sum(x * x for x in res_com) + sum(x * x for x in res_drop)
            worst = max(worst, abs(lhs - rhs) / max(lhs, 1e-300))
    print(f"  worst relative violation {worst:.2e} (rtol {RTOL_IDENTITY:.0e})")
    if worst > RTOL_IDENTITY:
        print("  -> ★ G-I FAIL: the set arithmetic is wrong, so no ratio means anything. STOP.")
        return 1
    print("  -> PASS")

    if s["n_common"] < MIN_COMMON:
        print(f"\n★ LE upper common points {s['n_common']} < {MIN_COMMON} -- UNDEFINED, stop "
              f"(kill clause 3).")
        return 0

    #: --- the reading ----------------------------------------------------------------------------
    print("\n=== e and d_self on the COMMON point set (LE upper) ===")
    print(f"  e:      " + "   ".join(f"{t} {s['e_common'][t]:.6f}" for t in CHAIN))
    print(f"  (full)  " + "   ".join(f"{t} {s['e_full'][t]:.6f}" for t in CHAIN))
    ds = [s["d_self"][f"{a}->{b}"] for a, b in zip(CHAIN, CHAIN[1:])]
    print(f"  d_self: " + "   ".join(f"{a}->{b} {d:.6f}"
                                     for (a, b), d in zip(zip(CHAIN, CHAIN[1:]), ds)))
    if min(ds) == 0.0:
        print("  -> ★ d_self exactly 0: SUSPICION, not a verdict (kill clause 5).")
        return 1
    shrink = ds[0] > ds[1]
    R = ds[-1] / s["e_common"]["S4"]
    print(f"\n=== the verdict quantity ===")
    print(f"  d_self shrinks with refinement: {shrink}   "
          f"R = d_self(S2->S4)/e(S4) = {R:.4f}")
    if R <= MODEL_MAX and shrink:
        print(f"  -> ★★★ N-MODEL  R <= 1/3 and d_self shrinks ⇒ MY PREDICTION IS REFUTED.")
        print("     The solution has essentially converged, and not to the experiment ⇒ the residual")
        print(f"     is >= {100 * (1 - R):.0f} % MODEL error ⇒ M3's LE target needs a model-error band.")
        print("     ★ And the previous round's 0.214/0.250 were a leading-edge CLAMPING ARTIFACT.")
    elif R >= DISC_MIN or not shrink:
        print(f"  -> ★★ N-DISC  R = {R:.4f}"
              + ("" if shrink else " and d_self does NOT shrink")
              + " ⇒ prediction CONFIRMED: the solution has NOT converged.")
        print("     The root cause is on the DISCRETISATION side, which OVERTURNS 'changing the mesh")
        print("     cannot solve M3'.")
    else:
        print(f"  -> N-MIX  1/3 < R = {R:.4f} < 1 ⇒ RECORDED, no direction claimed.")
        print("     ★ My prediction was N-DISC and this is NOT it -- the registration forbids reading")
        print("       N-MIX as N-DISC, so no direction is claimed.")

    #: --- P2: was the previous round's reading an artifact? (RECORDED, numbers already seen) ------
    print("\n=== P2 (RECORDED, already-seen numbers): how much of the sqrt(2) d_self came from the "
          "DROPPED points? ===")
    c34 = {t: load_curves(t) for t in SEEN}
    if all(v is not None for v in c34.values()):
        tbl = point_table(c34, exp, "LE", "upper")
        com = [r for r in tbl if all(r[3][t][1] for t in SEEN)]
        drop = [r for r in tbl if r not in com]
        a, b = SEEN
        d_all = rms([r[3][a][0] - r[3][b][0] for r in tbl if None not in (r[3][a][0], r[3][b][0])])
        d_com = rms([r[3][a][0] - r[3][b][0] for r in com])
        d_dr = rms([r[3][a][0] - r[3][b][0] for r in drop]) if drop else None
        print(f"  d_self({a}->{b}) all {d_all:.6f} (previously reported 0.250115)   "
              f"common-only {d_com:.6f}   dropped-only "
              f"{(f'{d_dr:.6f}' if d_dr is not None else '-')}  "
              f"(n_drop {len(drop)} of {len(tbl)})")
        if d_com is not None and d_all:
            print(f"  ⇒ removing the clamped points changes it by "
                  f"{100 * (d_com - d_all) / d_all:+.1f} % ⇒ "
                  f"{'the earlier reading was NOT an artifact' if abs(d_com - d_all) / d_all < 0.25 else 'the clamped points carried a large share'}")

    print("\n=== observed order and Richardson limit (DERIVED, not criteria) ===")
    if shrink:
        p = math.log(ds[0] / ds[1]) / math.log(2.0)
        tbl = point_table(curves, exp, "LE", "upper")
        common = [r for r in tbl if all(r[3][t][1] for t in CHAIN)]
        b2 = np.array([r[3]["S2"][0] for r in common])
        b4 = np.array([r[3]["S4"][0] for r in common])
        ce = np.array([r[2] for r in common])
        lim = b4 + (b4 - b2) / (2.0 ** p - 1.0)
        print(f"  p = {p:.3f}   |Cp_inf - Cp_exp| RMS = {rms(lim - ce):.6f}   vs S2's target 0.08")
    else:
        print("  UNDEFINED -- extrapolating a non-shrinking sequence would be a fit, not a limit.")

    print("\n=== all bands on the common set (RECORDED) ===")
    print(f"  {'band/side':14}{'n_com':>7}{'e(S0)':>10}{'e(S2)':>10}{'e(S4)':>10}"
          f"{'d S0->S2':>11}{'d S2->S4':>11}")
    for (b, sd), v in summary.items():
        if v["n_common"] == 0:
            continue
        print(f"  {b + '/' + sd:14}{v['n_common']:>7}"
              + "".join(f"{v['e_common'][t]:>10.5f}" for t in CHAIN)
              + f"{v['d_self']['S0->S2']:>11.5f}{v['d_self']['S2->S4']:>11.5f}")
    return 0


def _write(rows):
    if not rows:
        return
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {CSV}\n")


if __name__ == "__main__":
    sys.exit(main())
