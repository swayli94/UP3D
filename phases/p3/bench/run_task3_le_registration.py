"""How much of the LE 70 % is just CHORDWISE MISREGISTRATION? One parameter, one solve.

Pre-registered in phases/p3/docs/dev_phase_three/20260813-1700-le-registration-prereg.md, committed before this
file.

Entry check (in the registration): three hypotheses for the LE deficit are already REFUTED -- the
Kutta form, LE resolution (a major negative: refining LE alone gives at best -8.3 % and is
NON-MONOTONE, and a resolution error cannot be non-monotone), and geometry via the h_te split field
(bit-identical). The tip route was tested with a gate its own author recorded as mis-designed. None of
those is re-measured here.

What has never been asked: the LE band is x/c in [0, 0.15] on the UPPER surface, where the experimental
Cp is very steep, so a small chordwise misregistration inflates the RMS enormously. Every earlier round
asked why the AMPLITUDE is wrong. This asks how much is merely POSITION.

★★★ The interpretation ban, fixed before any number: "a shift explains the error" does NOT mean the
error is absent. A displaced suction peak is still a solver error. This RE-ATTRIBUTES it from the
amplitude family (G1.6's faceted geometry, P1 wall recovery) to the position family (leading-edge
geometric discretisation, flow direction, circulation) -- which imply completely different work.

★ Sign convention, stated so a sign error cannot invert the reading: delta > 0 samples the COMPUTED
curve at x_experiment + delta, i.e. it is equivalent to shifting the computed curve UPSTREAM by delta.

★ np.interp CLAMPS at the ends rather than extrapolating, so an out-of-range point would silently take
an endpoint value. Validity is therefore checked EXPLICITLY and an out-of-range delta is declared
INVALID (registration, question 2).

Outputs (TRACKED): bench/gate_results/task3_le_registration.csv
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

CSV = os.path.join(HERE, "gate_results", "task3_le_registration.csv")
LEVEL = "medium"
#: declared in advance: scan range and step
DELTAS = np.round(np.arange(-0.020, 0.02001, 0.001), 6)
#: S2's OWN registered intermediate target -- not a number picked here
TARGET_RMS, H1_MAX_SHIFT, H2_FLOOR = 0.08, 0.010, 0.15
LEG_GATE_S = 900.0


def le_band_rms(curves, exp, etas, delta, band="LE", side="upper"):
    """Pooled RMS over `etas` for one band/side, sampling the computed curve at x_exp + delta.

    Returns (rms, n_points, valid). `valid` is False if ANY sampled abscissa falls outside the
    computed curve's range -- np.interp would clamp there, silently substituting an endpoint value.
    """
    lo, hi = next((l, h) for n, l, h in BANDS if n == band)
    ss, nn, valid = 0.0, 0, True
    for eta in etas:
        e = exp[eta]
        m = (e["upper"] == (side == "upper")) & (e["x"] >= lo) & (e["x"] < hi)
        if not np.any(m):
            continue
        cx = curves[eta][f"x_{side}"]
        xs = e["x"][m] + delta
        if xs.min() < float(np.min(cx)) or xs.max() > float(np.max(cx)):
            valid = False
        cp_i = np.interp(xs, cx, curves[eta][f"cp_{side}"])
        ss += float(np.sum((cp_i - e["cp"][m]) ** 2))
        nn += int(m.sum())
    return (float(np.sqrt(ss / nn)) if nn else None), nn, valid


def main():
    print("resolved threads: " + ", ".join(
        f"{k}={os.environ.get(k)}" for k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS",
                                             "OPENBLAS_NUM_THREADS")))
    print(f"load average: {os.getloadavg()}\n")
    exp = parse_experiment()
    path = os.path.join(REPO, "cases", "meshes", "onera_m6", f"{LEVEL}.msh")
    mc, wc = cut_wake(read_mesh(path))

    t0 = time.perf_counter()
    r = solve(mc, wc, entropy=True, kutta="pressure", taper=True, probe_seed=0)
    wall = time.perf_counter() - t0
    print(f"solve: converged={r.get('converged')} |R|={r['residual_history'][-1]:.3e} "
          f"({wall:.0f}s)")
    if not r.get("converged"):
        print("★ not converged -- no curves to analyse. STOP (kill clause 1).")
        return 1
    if wall > LEG_GATE_S:
        print(f"★ cost gate {wall:.0f}s > {LEG_GATE_S:.0f}s -- stop")
        return 1

    phi = np.asarray(r["phi"])
    etas = ETAS[:N_UNMASKED]
    curves = {eta: section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
              for eta in ETAS}

    #: the steepness that motivates the whole question, measured from the EXPERIMENT
    slopes = []
    for eta in etas:
        e = exp[eta]
        m = (e["upper"]) & (e["x"] >= 0.0) & (e["x"] < 0.15)
        x, cp = e["x"][m], e["cp"][m]
        o = np.argsort(x)
        slopes.extend(np.abs(np.diff(cp[o]) / np.maximum(np.diff(x[o]), 1e-12)))
    slope_med = float(np.median(slopes))
    print(f"experiment LE-band |dCp/dx| median = {slope_med:.2f} per chord   "
          f"=> 0.001 c of shift ~ {0.001 * slope_med:.4f} in Cp")

    rows = []
    for d in DELTAS:
        rec = dict(delta=float(d))
        for band in ("LE", "MID", "TE"):
            rms, n, ok = le_band_rms(curves, exp, etas, float(d), band=band, side="upper")
            rec[f"rms_{band}_upper"] = rms
            rec[f"n_{band}_upper"] = n
            rec[f"valid_{band}"] = ok
        rec["valid"] = rec["valid_LE"]
        rows.append(rec)
    #: per-station optimum -- ★ if the stations disagree on the SIGN it is not a global shift
    per_station = {}
    for eta in etas:
        best = None
        for d in DELTAS:
            rms, n, ok = le_band_rms(curves, exp, [eta], float(d))
            if ok and rms is not None and (best is None or rms < best[1]):
                best = (float(d), rms)
        per_station[eta] = best
    _write(rows, per_station, slope_med)

    valid = [x for x in rows if x["valid"] and x["rms_LE_upper"] is not None]
    base = next(x for x in rows if x["delta"] == 0.0)
    print(f"\nvalid deltas: {len(valid)}/{len(rows)}   baseline (delta=0) LE upper RMS "
          f"{base['rms_LE_upper']:.6f}")
    if len(valid) < len(DELTAS) // 2:
        print("★ fewer than half the deltas are valid -- the scan is meaningless. STOP (clause 2).")
        return 1

    best = min(valid, key=lambda x: x["rms_LE_upper"])
    print(f"best: delta = {best['delta']:+.3f} c   LE upper RMS {best['rms_LE_upper']:.6f}   "
          f"(MID {best['rms_MID_upper']:.6f}, TE {best['rms_TE_upper']:.6f})")
    print(f"      at delta=0: MID {base['rms_MID_upper']:.6f}, TE {base['rms_TE_upper']:.6f}")

    print("\nper-station optimum (★ sign disagreement means it is NOT a global misregistration):")
    signs = set()
    for eta, b in per_station.items():
        if b is None:
            print(f"  eta {eta}: no valid delta"); continue
        signs.add(np.sign(b[0]) if b[0] != 0 else 0)
        print(f"  eta {eta:.2f}  best delta {b[0]:+.3f} c   LE RMS {b[1]:.6f}")
    print(f"  distinct signs among station optima: {sorted(signs)}")

    print("\n=== H1 / H2 / H3 (registration section 5) ===")
    on_edge = abs(best["delta"]) >= float(DELTAS.max()) - 1e-12
    if on_edge:
        print(f"  -> H3   the minimum sits ON the scan boundary ({best['delta']:+.3f} c) -- "
              f"boundary-limited, and clause 3 forbids extrapolating. RECORDED.")
    elif abs(best["delta"]) <= H1_MAX_SHIFT and best["rms_LE_upper"] <= TARGET_RMS:
        print(f"  -> ★ H1   |delta| {abs(best['delta']):.3f} <= {H1_MAX_SHIFT} and RMS "
              f"{best['rms_LE_upper']:.4f} <= {TARGET_RMS} ⇒ POSITION dominates the LE deficit")
        print("     ★ and this does NOT mean the error is absent: a displaced peak is still a")
        print("       solver error. It re-attributes the work from amplitude to position.")
    elif min(x["rms_LE_upper"] for x in valid) > H2_FLOOR:
        print(f"  -> ★ H2   no shift in the scan brings LE below {H2_FLOOR} (best "
              f"{best['rms_LE_upper']:.4f}) ⇒ position CANNOT explain it; the deficit is genuine")
        print("     amplitude error, and the next suspect is the G1.6 family (faceted geometry +")
        print("     P1 wall recovery) -- for which P11 already measured the curved-element route")
        print("     NEGATIVE, so this would be a hard conclusion.")
    else:
        print(f"  -> H3   best RMS {best['rms_LE_upper']:.4f} at delta {best['delta']:+.3f} c "
              f"falls between the bands. RECORDED, no direction claimed.")
    return 0


def _write(rows, per_station, slope_med):
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    p2 = CSV.replace(".csv", "_stations.csv")
    with open(p2, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["eta", "best_delta", "best_rms", "exp_LE_slope_median"])
        for eta, b in per_station.items():
            w.writerow([eta, (b[0] if b else None), (b[1] if b else None), slope_med])
    print(f"\nwrote {CSV}\nwrote {p2}")


if __name__ == "__main__":
    sys.exit(main())
