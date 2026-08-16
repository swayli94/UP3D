"""Is the LE error a CONSTANT OFFSET or SCATTER? One solve, one identity, no fit.

Pre-registered in docs/dev_phase_three/20260814-0100-le-offset-scatter-prereg.md.

Established before this round: the LE deficit is a LEVEL error, not a shape or position one -- eta 0.90
is the worst station while sitting on the FLATTEST part of the Cp curve, so position cannot explain it
there even in principle, and a chordwise shift scan put its optimum at delta = 0. The candidate list is
down to the G1.6 family.

But a level error still splits into two different diseases: a CONSTANT OFFSET points at local
loading/circulation, while SCATTER points at wall recovery and faceted geometry. This decides where the
next effort goes.

★ The quantity is an IDENTITY, not a fit: with d = Cp_comp - Cp_exp, RMS^2 = bias^2 + var, so
f_bias = bias^2 / RMS^2 is dimensionless and bounded in [0, 1] -- no hand-picked physical threshold is
needed, which is the absolute-threshold defect this session has already tripped over. G-I asserts the
identity to rtol 1e-12, because if it does not hold no f_bias means anything.

★ The bias SIGN is half the diagnosis: bias < 0 means the computed Cp is more negative, i.e. the
leading edge is over-accelerating there.

Outputs (TRACKED): bench/gate_results/task3_le_offset_scatter.csv
"""

import csv
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from run_m3_budget import (ALPHA, B_SEMI, BANDS, ETAS, M_INF,       # noqa: E402
                           N_UNMASKED, parse_experiment, solve)

CSV = os.path.join(HERE, "gate_results", "task3_le_offset_scatter.csv")
#: G-R: this leg's committed pooled LE upper RMS, bit-for-bit
POOLED_LE_REF = 0.242221
O1_MIN, O2_MAX, RTOL_IDENTITY = 0.7, 0.3, 1e-12
LEG_GATE_S = 900.0


def decompose(curves, exp, eta, band, side):
    """bias / var / rms for one (station, band, side). The identity is the point, not a fit."""
    lo, hi = next((l, h) for n, l, h in BANDS if n == band)
    e = exp[eta]
    m = (e["upper"] == (side == "upper")) & (e["x"] >= lo) & (e["x"] < hi)
    if not np.any(m):
        return None
    cp_i = np.interp(e["x"][m], curves[eta][f"x_{side}"], curves[eta][f"cp_{side}"])
    d = cp_i - e["cp"][m]
    bias = float(np.mean(d))
    var = float(np.mean((d - bias) ** 2))
    rms = float(np.sqrt(np.mean(d ** 2)))
    return dict(eta=eta, band=band, side=side, n=int(m.sum()), bias=bias,
                sd=float(np.sqrt(var)), rms=rms,
                f_bias=(bias * bias / (rms * rms) if rms > 0 else None),
                identity_rel=(abs(bias * bias + var - rms * rms) / max(rms * rms, 1e-300)))


def main():
    print("resolved threads: " + ", ".join(
        f"{k}={os.environ.get(k)}" for k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS",
                                             "OPENBLAS_NUM_THREADS")))
    print(f"load average: {os.getloadavg()}\n")
    exp = parse_experiment()
    mc, wc = cut_wake(read_mesh(os.path.join(REPO, "cases", "meshes", "onera_m6", "medium.msh")))

    t0 = time.perf_counter()
    r = solve(mc, wc, entropy=True, kutta="pressure", taper=True, probe_seed=0, taper_rc=0.05)
    wall = time.perf_counter() - t0
    print(f"solve: converged={r.get('converged')} |R|={r['residual_history'][-1]:.3e} ({wall:.0f}s)")
    if not r.get("converged"):
        print("★ not converged -- nothing to decompose. STOP (kill clause 1).")
        return 1
    if wall > LEG_GATE_S:
        print(f"★ cost gate {wall:.0f}s -- stop"); return 1

    phi = np.asarray(r["phi"])
    etas = ETAS[:N_UNMASKED]
    curves = {eta: section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF) for eta in ETAS}

    rows = [d for eta in etas for band, _, _ in BANDS for side in ("upper", "lower")
            if (d := decompose(curves, exp, eta, band, side)) is not None]
    _write(rows)

    #: --- G-I: the identity, or nothing below means anything ---------------------------------
    worst = max(x["identity_rel"] for x in rows)
    print(f"\n=== G-I: bias^2 + var == RMS^2   worst relative violation {worst:.2e} "
          f"(rtol {RTOL_IDENTITY:.0e}) ===")
    if worst > RTOL_IDENTITY:
        print("  -> ★ G-I FAIL: the decomposition is wrong, so no f_bias is meaningful. STOP.")
        return 1
    print("  -> PASS")

    #: --- G-R: same leg as the previous rounds ------------------------------------------------
    ss = sum(x["rms"] ** 2 * x["n"] for x in rows if x["band"] == "LE" and x["side"] == "upper")
    nn = sum(x["n"] for x in rows if x["band"] == "LE" and x["side"] == "upper")
    pooled = float(np.sqrt(ss / nn))
    print(f"\n=== G-R: pooled LE upper RMS {pooled:.6f} vs committed {POOLED_LE_REF} ===")
    if abs(pooled - POOLED_LE_REF) > 5e-7:
        print("  -> ★ G-R FAIL: this is not the same leg. STOP and align.")
        return 1
    print("  -> PASS (same leg)")

    #: --- the reading -------------------------------------------------------------------------
    #: ★★ 1/n is printed next to f_bias BY CONSTRUCTION, because comparing f_bias against ZERO is
    #: the wrong comparison: for residuals with NO bias, bias has variance sigma^2/n while rms^2 is
    #: about sigma^2, so E[f_bias] ~ 1/n. Without this column the TE/lower band -- which has n = 1
    #: at every station, making f_bias identically 1.0 -- reads as a beautiful "pure offset"
    #: finding. It is arithmetic. Same family as the project's recorded vacuous PASS where
    #: `d2h == 0` counted as success because the baseline was already 0.
    print("\n=== decomposition, LE upper (binding station = eta 0.90) ===")
    print(f"  {'eta':>6}{'n':>4}{'bias':>11}{'sd':>11}{'rms':>11}{'f_bias':>9}{'1/n':>8}"
          f"   note")
    le_up = [x for x in rows if x["band"] == "LE" and x["side"] == "upper"]
    for x in le_up:
        inv = 1.0 / x["n"]
        note = ("VACUOUS (n=1: f_bias is identically 1)" if x["n"] == 1
                else "consistent with PURE SCATTER" if x["f_bias"] <= inv
                else "slightly above 1/n" if x["f_bias"] <= 2 * inv
                else "★ a REAL bias")
        print(f"  {x['eta']:>6.2f}{x['n']:>4}{x['bias']:>11.6f}{x['sd']:>11.6f}"
              f"{x['rms']:>11.6f}{x['f_bias']:>9.3f}{inv:>8.3f}   {note}")
    signs = {int(np.sign(x["bias"])) for x in le_up}
    print(f"\n  G-S (RECORDED): bias signs across stations = {sorted(signs)}"
          + ("   ★ CONSISTENT -- a global level error is possible" if len(signs) == 1
             else "   ★ SIGNS FLIP -- this is NOT one global level error"))

    print("\n  other bands (RECORDED):")
    for band in ("MID", "TE"):
        for side in ("upper", "lower"):
            sub = [x for x in rows if x["band"] == band and x["side"] == side]
            if not sub:
                continue
            fb = np.mean([x["f_bias"] for x in sub])
            #: ★ per-station 1/n comparison, not the mean f_bias alone -- see the note above
            above = sum(1 for x in sub if x["n"] > 1 and x["f_bias"] > 2.0 / x["n"])
            vac = sum(1 for x in sub if x["n"] == 1)
            print(f"    {band:4}/{side:6} mean f_bias {fb:.3f}   "
                  f"bias range [{min(x['bias'] for x in sub):+.4f}, "
                  f"{max(x['bias'] for x in sub):+.4f}]   "
                  f"{above}/{len(sub)} stations show a REAL bias"
                  + (f"   ★ {vac} VACUOUS (n=1)" if vac else ""))

    b = next(x for x in le_up if abs(x["eta"] - 0.90) < 1e-9)
    print(f"\n=== O1 / O2 / O3 (binding: eta 0.90 LE upper) ===")
    print(f"  f_bias = {b['f_bias']:.3f}   bias = {b['bias']:+.6f}   sd = {b['sd']:.6f}")
    if b["f_bias"] >= O1_MIN:
        print(f"  -> ★ O1  f_bias >= {O1_MIN} ⇒ the error is mostly a LEVEL OFFSET")
        print(f"     bias is {'NEGATIVE (computed Cp more negative -> the LE OVER-accelerates)' if b['bias'] < 0 else 'POSITIVE (computed Cp less negative -> the LE UNDER-accelerates)'}")
        print("     ⇒ points at LOCAL LOADING / CIRCULATION, not at wall recovery.")
        print("     ★ and this OVERTURNS the current 'only the G1.6 family remains' state.")
    elif b["f_bias"] <= O2_MAX:
        print(f"  -> ★ O2  f_bias <= {O2_MAX} ⇒ the error is mostly SCATTER")
        print("     ⇒ points at WALL RECOVERY / faceted geometry, i.e. the G1.6 family --")
        print("       consistent with it being the only candidate left, which TIGHTENS that")
        print("       conclusion rather than merely repeating it.")
    else:
        print(f"  -> O3  {O2_MAX} < f_bias < {O1_MIN} ⇒ MIXED. RECORDED, no direction claimed;")
        print(f"     quantitatively: |bias| {abs(b['bias']):.6f} against sd {b['sd']:.6f}.")
    return 0


def _write(rows):
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    keys = list(rows[0])
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
