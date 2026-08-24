"""Is the non-convergence confined to the leading-edge NEAR-SINGULAR region? ZERO SOLVE.

Pre-registered in phases/p3/docs/dev_phase_three/20260815-1800-nearsingular-partition-prereg.md.

The explanation under test, from already-committed numbers: at the finest level's worst station one wall cell
is about 0.0102 chord, the sonic front foot sits at 0.00405 chord (inside the FIRST cell), and the most
upstream scored experimental point is at 0.00057 chord -- a small fraction of one cell from the leading edge.
A P1 field cannot represent a stagnation-to-supersonic transition inside one cell, and refinement RELOCATES
that failure instead of removing it.

★★★ Its one falsifiable prediction: split the band and the INNER (near-singular) region should not converge
while the OUTER region should. If the outer region also fails, the explanation is wrong.

★★ No new threshold: R = d_self/e and its bounds (<= 1/3 converged, >= 1 not converged) were registered in
the self-convergence round; this round only computes them PER REGION.

★★★ The split is x_foot(S2) -- the COARSER level's sonic foot -- so both levels share ONE split point and
therefore ONE point set. Same logic as the common-point-set fix; the fifth standing question again.

★ Controls: LE lower and MID upper, already measured to converge normally, are partitioned identically. If
their outer regions come back non-convergent, the partition implementation is suspect, not the hypothesis.

Outputs (TRACKED): bench/gate_results/task3_nearsingular_partition.csv
"""

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.meshgen.wing3d import B_SEMI, chord_at                    # noqa: E402
from run_m3_budget import BANDS, ETAS, M_INF, parse_experiment        # noqa: E402
import run_task3_nonconv_discriminator as ND                         # noqa: E402
import run_task3_shape_and_taper as ST                               # noqa: E402
import run_task3_tipstations_mechanism as TM                         # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_nearsingular_partition.csv")
CHAIN, PAIR = ("S0", "S2", "S4"), ("S2", "S4")
UNMASKED = (0.20, 0.44, 0.65, 0.80, 0.90)
GR7 = ST.GR7
#: reused verbatim from the self-convergence registration -- no new numbers
R_CONV, R_NOCONV = 1.0 / 3.0, 1.0
N_MIN, RTOL_MASS = 10, 1e-12
#: measured LE-band wall tangential spacing of each level (metres), from the committed G-M tables
LE_HT = {"S2": 0.010508, "S4": 0.004979}


def regionals(xs, ce, per, split):
    """(e_S2, e_S4, d_self, R, n) on the sub-set selected by `split`."""
    if not np.any(split):
        return None
    d = per[PAIR[0]][split] - per[PAIR[1]][split]
    e2 = float(np.sqrt(np.mean((per[PAIR[0]][split] - ce[split]) ** 2)))
    e4 = float(np.sqrt(np.mean((per[PAIR[1]][split] - ce[split]) ** 2)))
    ds = float(np.sqrt(np.mean(d ** 2)))
    return dict(n=int(split.sum()), e_S2=e2, e_S4=e4, d_self=ds,
                R=(ds / e4 if e4 else None), conv_ratio=(e4 / e2 if e2 else None),
                sumsq=float(np.sum(d ** 2)))


def pool(parts):
    """Pool sub-results over stations: RMS pooling by n, and the summed d^2 for G-I."""
    parts = [p for p in parts if p]
    if not parts:
        return None
    n = sum(p["n"] for p in parts)
    f = lambda k: float(np.sqrt(sum(p[k] ** 2 * p["n"] for p in parts) / n))   # noqa: E731
    e2, e4, ds = f("e_S2"), f("e_S4"), f("d_self")
    return dict(n=n, e_S2=e2, e_S4=e4, d_self=ds, R=(ds / e4 if e4 else None),
                conv_ratio=(e4 / e2 if e2 else None),
                sumsq=sum(p["sumsq"] for p in parts))


def analyse(curves, exp, band, side, split_mode, cp_star):
    """Partition one (band, side) and pool over the unmasked stations."""
    per_station, inner, outer, full = {}, [], [], []
    for eta in UNMASKED:
        got = TM.band_points(curves, exp, eta, band=band, side=side)
        if got[0] is None:
            continue
        xs, ce, per = got[0], got[1], got[2]
        if split_mode == "foot":
            xs_split = TM.geometry(curves[PAIR[0]][eta], cp_star)[0]
        else:                                            # sensitivity: 2 local cells
            xs_split = 2.0 * LE_HT[PAIR[0]] / float(chord_at(eta * B_SEMI))
        if xs_split is None:
            continue
        m_in = xs < xs_split
        ri = regionals(xs, ce, per, m_in)
        ro = regionals(xs, ce, per, ~m_in)
        rf = regionals(xs, ce, per, np.ones_like(m_in))
        per_station[eta] = dict(x_split=float(xs_split), inner=ri, outer=ro, full=rf)
        if ri:
            inner.append(ri)
        if ro:
            outer.append(ro)
        if rf:
            full.append(rf)
    return dict(per_station=per_station, inner=pool(inner), outer=pool(outer), full=pool(full))


def main():
    src = open(__file__).read().split(chr(34) * 3, 2)[2]
    for forbidden in ("solve" + "_newton", "run_m3_budget." + "solve"):
        assert forbidden not in src, f"G-Z: no solver ({forbidden})"
    print("★ G-Z: zero solves\n")
    exp = parse_experiment()
    cp_star = ND.cp_critical(M_INF)
    curves = {}
    for t in CHAIN:
        c, _r = TM.solve_all_stations(t, exp)
        if c is None:
            print(f"★ {t} curves unavailable -- STOP")
            return 1
        curves[t] = c

    #: --- G-R ------------------------------------------------------------------------------------
    print("=== G-R: full-band per-station d_self must reproduce bit-for-bit ===")
    for eta in ETAS:
        got = TM.band_points(curves, exp, eta)
        if got[0] is None:
            continue
        r = float(np.sqrt(np.mean((got[2][PAIR[0]] - got[2][PAIR[1]]) ** 2)))
        if abs(r - GR7[eta]) > 5e-7:
            print(f"  eta {eta:.2f}: {r:.6f} vs {GR7[eta]:.6f} ★ FAIL -> STOP (kill clause 1)")
            return 1
    print("  -> PASS (all seven)")

    rows = []
    A = analyse(curves, exp, "LE", "upper", "foot", cp_star)

    #: --- G-I ------------------------------------------------------------------------------------
    lhs = A["inner"]["sumsq"] + A["outer"]["sumsq"]
    rhs = A["full"]["sumsq"]
    rel = abs(lhs - rhs) / max(rhs, 1e-300)
    print(f"\n=== G-I: sum d^2(inner) + sum d^2(outer) == sum d^2(full) ===")
    print(f"  {lhs:.12e} vs {rhs:.12e}   rel {rel:.2e} (rtol {RTOL_MASS:.0e}) -> "
          f"{'PASS' if rel <= RTOL_MASS else '★ FAIL'}")
    if rel > RTOL_MASS:
        print("  -> ★ G-I FAIL: the partition is not exclusive-and-exhaustive. STOP (kill clause 2).")
        return 1

    print("\n=== the split point (binding: x_foot(S2), the COARSER level's sonic foot) ===")
    print(f"  {'eta':>6}{'x_split':>10}{'n_inner':>9}{'n_outer':>9}")
    for eta, v in A["per_station"].items():
        print(f"  {eta:>6.2f}{v['x_split']:>10.5f}"
              f"{(v['inner']['n'] if v['inner'] else 0):>9}"
              f"{(v['outer']['n'] if v['outer'] else 0):>9}")
    ni = A["inner"]["n"] if A["inner"] else 0
    no = A["outer"]["n"] if A["outer"] else 0
    print(f"  pooled: n_inner {ni}, n_outer {no}   (need >= {N_MIN} each)")
    if ni < N_MIN or no < N_MIN:
        print(f"  -> ★ UNDEFINED: pooled sample too thin (kill clause 3). The sampling cannot")
        print("     resolve this boundary; no direction claimed.")
        _write(rows, A, None)
        return 0

    print("\n=== binding: LE upper, pooled over the 5 unmasked stations ===")
    print(f"  {'region':8}{'n':>5}{'e(S2)':>10}{'e(S4)':>10}{'e4/e2':>9}{'d_self':>10}{'R':>9}"
          f"   reading")
    for name in ("inner", "outer", "full"):
        v = A[name]
        rr = v["R"]
        note = ("★ NOT converged (R >= 1)" if rr >= R_NOCONV
                else "★ converged (R <= 1/3)" if rr <= R_CONV else "mixed")
        print(f"  {name:8}{v['n']:>5}{v['e_S2']:>10.6f}{v['e_S4']:>10.6f}"
              f"{v['conv_ratio']:>9.3f}{v['d_self']:>10.6f}{rr:>9.4f}   {note}")
        rows.append(dict(kind="binding", region=name, **v))

    Ri, Ro = A["inner"]["R"], A["outer"]["R"]
    print(f"\n  R_inner = {Ri:.4f}   R_outer = {Ro:.4f}")
    if Ro <= R_CONV and Ri >= R_NOCONV:
        print("  -> ★★★ P-CONF  the non-convergence is CONFINED to the near-singular region and the")
        print("     rest of the band converges normally ⇒ 'the LE band is unreachable' becomes an")
        print("     EXPLAINED fact, and M3 should be reported per region (targets unchanged).")
    elif Ro >= R_NOCONV:
        print("  -> ★★ P-REFUTE  the OUTER region does not converge either ⇒ MY EXPLANATION IS WRONG:")
        print("     the non-convergence is NOT confined to the near-singular region, and the")
        print("     'close it out' route does not hold. Ten rounds of synthesis need reorganising.")
    else:
        print(f"  -> P-MIX  R_outer = {Ro:.4f} sits between 1/3 and 1 ⇒ RECORDED, no direction claimed.")

    #: --- controls (tool self-check, pre-registered) ----------------------------------------------
    print("\n=== ★ tool self-check: CONTROLS already measured to converge normally ===")
    print(f"  {'band/side':12}{'region':8}{'n':>5}{'e(S4)':>10}{'d_self':>10}{'R':>9}   reading")
    for band, side in (("LE", "lower"), ("MID", "upper")):
        C = analyse(curves, exp, band, side, "foot", cp_star)
        for name in ("inner", "outer"):
            v = C[name]
            if not v:
                print(f"  {band + '/' + side:12}{name:8}   (empty)")
                continue
            rr = v["R"]
            note = ("★ NOT converged" if rr >= R_NOCONV
                    else "converged" if rr <= R_CONV else "mixed")
            print(f"  {band + '/' + side:12}{name:8}{v['n']:>5}{v['e_S4']:>10.6f}"
                  f"{v['d_self']:>10.6f}{rr:>9.4f}   {note}")
            rows.append(dict(kind="control", band=band, side=side, region=name, **v))
        if C["outer"] and C["outer"]["R"] >= R_NOCONV:
            print(f"  ★ {band}/{side} OUTER also non-convergent -> the PARTITION is suspect, not the")
            print("    hypothesis (kill clause 4): check the implementation before reading the verdict.")

    #: --- sensitivity (reported, may NOT change the verdict) --------------------------------------
    print("\n=== sensitivity: x_split = 2 local cells (REPORTED ONLY, cannot change the verdict) ===")
    S = analyse(curves, exp, "LE", "upper", "cells", cp_star)
    if S["inner"] and S["outer"]:
        print(f"  {'region':8}{'n':>5}{'e(S4)':>10}{'d_self':>10}{'R':>9}")
        for name in ("inner", "outer"):
            v = S[name]
            print(f"  {name:8}{v['n']:>5}{v['e_S4']:>10.6f}{v['d_self']:>10.6f}{v['R']:>9.4f}")
            rows.append(dict(kind="sensitivity", region=name, **v))
        print(f"  x_split per station: "
              + ", ".join(f"{e:.2f}:{v['x_split']:.4f}" for e, v in S["per_station"].items()))
    else:
        print("  (one region empty at this split -- reported as unusable)")

    print("\n=== per-station (REPORTED, not judged -- 2-4 points per region) ===")
    print(f"  {'eta':>6}{'x_split':>10}{'in n':>6}{'in R':>9}{'out n':>7}{'out R':>9}")
    for eta, v in A["per_station"].items():
        ri, ro = v["inner"], v["outer"]
        print(f"  {eta:>6.2f}{v['x_split']:>10.5f}"
              f"{(ri['n'] if ri else 0):>6}{(ri['R'] if ri and ri['R'] else float('nan')):>9.3f}"
              f"{(ro['n'] if ro else 0):>7}{(ro['R'] if ro and ro['R'] else float('nan')):>9.3f}")
        rows.append(dict(kind="station", eta=eta, x_split=v["x_split"],
                         inner_n=(ri["n"] if ri else 0), inner_R=(ri["R"] if ri else None),
                         outer_n=(ro["n"] if ro else 0), outer_R=(ro["R"] if ro else None)))
    _write(rows, A, None)
    return 0


def _write(rows, A, _unused):
    if not rows:
        rows = [dict(kind="undefined", note="pooled sample too thin")]
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
