"""Is the LE-band error DISCRETISATION or MODEL error? A self-convergence (Cauchy) read.

Pre-registered in docs/dev_phase_three/20260814-1800-selfconvergence-prereg.md.

Every earlier round measured e(h) = RMS[Cp_h - Cp_exp]. Whether the solution has CONVERGED needs the
solution against itself: d_self(h) = RMS[Cp_h - Cp_finer] at the same experimental abscissae, so both
quantities are Cp RMS values in the same units and the verdict is their RATIO -- no physical threshold.

  d_self shrinks and is far below e  -> the solution converged to something that is NOT the experiment
                                       => the residual is MODEL error, and M3's target needs a band
  d_self does not shrink / ~ e       -> the solution has not converged => the root cause is discretisation

★★ The family must be repaired first: the default h_far = min(2.5, 120 h_wall) CLAMP breaks self-similarity
above h_wall 0.0208 (generate_onera_m6._level_params, recorded 2026-07-13 / P13-G13.3: "coarse is NOT on the
same refinement ray ... any THREE-POINT RICHARDSON is INVALID"). So h_far is passed EXPLICITLY as
120 h_wall on a factor-sqrt(2) ray. Two anchors follow for free: S4's h_far equals the unclamped default so
S4 must bit-reproduce C_medium (G-R), and S2's h_far differs from the clamped default so S2-vs-C_coarse
measures what the clamp does to the LE band -- never quantified before.

★ G-X: np.interp CLAMPS outside its range instead of extrapolating. That already invalidated a leg in the
LE registration round, so every (level, band, side) records whether the experimental abscissae actually lie
inside the computed curve, and an invalid one is dropped rather than silently flattened.

Outputs (TRACKED): bench/gate_results/task3_selfconvergence.csv
"""

import csv
import math
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.wake_cut import cut_wake                             # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                              # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                  # noqa: E402
from pyfp3d.post.surface import planform_area                         # noqa: E402
from run_m3_budget import (ALPHA, BANDS, ETAS, M_INF, N_UNMASKED,     # noqa: E402
                           parse_experiment, solve)
import run_task3_fixed_budget as FB                                   # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_selfconvergence.csv")
PREV = os.path.join(HERE, "gate_results", "task3_fixed_budget.csv")
#: G-R: S4 uses h_far = 120*0.015 = 1.80, which IS the unclamped default -> same mesh as C_medium
GR_REF = 0.242221
#: the factor-sqrt(2) ray, h_far EXPLICIT (section 2 of the registration)
RAY = tuple((f"S{i}", h) for i, h in enumerate((0.0600, 0.0424, 0.0300, 0.0212, 0.0150)))
MODEL_MAX, DISC_MIN = 1.0 / 3.0, 1.0
LEG_GATE_S, TOTAL_GATE_S = 900.0, 3600.0


def knobs(h_wall):
    return dict(h_wall=h_wall, h_edge=0.5 * h_wall, h_te=0.5 * h_wall, h_far=120.0 * h_wall)


def curves_of(mc, phi, tag):
    """★ CACHE BEFORE YOU REPORT (a project discipline bought with a 40-minute solve): the curves go
    to disk as they are produced, so a crash in the reporting layer costs nothing instead of the
    whole ray. Reloaded on a re-run."""
    path = os.path.join(FB.SCRATCH, f"curves_{tag}.npz")
    if os.path.exists(path):
        z = np.load(path)
        etas = sorted({float(k.split("|")[0]) for k in z.files})
        print(f"  [cached curves] {tag}", flush=True)
        return {e: {k.split("|")[1]: z[k] for k in z.files if float(k.split("|")[0]) == e}
                for e in etas}
    out = {}
    for eta in ETAS[:N_UNMASKED]:
        try:
            out[eta] = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
        except Exception:                                             # noqa: BLE001
            pass
    try:
        os.makedirs(FB.SCRATCH, exist_ok=True)
        np.savez(path, **{f"{e}|{k}": np.asarray(v) for e, c in out.items()
                          for k, v in c.items() if np.ndim(v) >= 1})
    except Exception as exc:                                          # noqa: BLE001
        print(f"  (curves not cached: {exc})")
    return out


def sampled(curves, exp, band, side):
    """Cp of this level at the experimental abscissae of one (band, side), pooled over stations.

    Returns (vector, exp vector, valid). ★ G-X: np.interp clamps outside its range, so `valid` is
    False if any abscissa falls outside the computed curve -- a clamped value is not a reading.
    """
    lo, hi = next((l, h) for n, l, h in BANDS if n == band)
    cp, ce, valid = [], [], True
    for eta in ETAS[:N_UNMASKED]:
        if eta not in curves or eta not in exp:
            valid = False
            continue
        e = exp[eta]
        m = (e["upper"] == (side == "upper")) & (e["x"] >= lo) & (e["x"] < hi)
        if not np.any(m):
            continue
        xs, cs = curves[eta][f"x_{side}"], curves[eta][f"cp_{side}"]
        if e["x"][m].min() < np.min(xs) or e["x"][m].max() > np.max(xs):
            valid = False
        cp.append(np.interp(e["x"][m], xs, cs))
        ce.append(e["cp"][m])
    if not cp:
        return None, None, False
    return np.concatenate(cp), np.concatenate(ce), valid


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def main():
    print("resolved threads: " + ", ".join(
        f"{k}={os.environ.get(k)}" for k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS",
                                             "OPENBLAS_NUM_THREADS")))
    print(f"load average: {os.getloadavg()}")
    print(f"★ G-C: tip_cap = {FB.TIP_CAP};  h_far passed EXPLICITLY as 120*h_wall "
          f"(the clamp breaks the ray)\n")
    exp = parse_experiment()
    rows, store, elapsed = [], {}, [0.0]

    for tag, h in RAY:
        if elapsed[0] > TOTAL_GATE_S:
            print(f"  ★ total gate exceeded -- {tag} NOT run (kill clause 4)")
            break
        k = knobs(h)
        mc, wc = cut_wake(FB.cached_mesh(tag, k["h_wall"], k["h_edge"], k["h_te"], k["h_far"]))
        mrow = FB.mesh_row(mc)
        s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        t0 = time.perf_counter()
        r = solve(mc, wc, entropy=True, kutta="pressure", taper=True, probe_seed=0, taper_rc=0.05)
        wall = time.perf_counter() - t0
        elapsed[0] += wall
        conv = bool(r.get("converged"))
        res = float(r.get("residual_history", [float("nan")])[-1])
        row = dict(leg=tag, tip_cap=FB.TIP_CAP, **k, converged=conv, res_final=res,
                   n_limited=r.get("n_limited"), n_floored=r.get("n_floored"),
                   solve_s=round(wall, 1), near_pct=100.0 * float(
                       (FB.cell_sizes(mc)[1] < 1.0).sum()) / len(mc.elements), **mrow)
        if conv:
            row.update(FB.flow_row(mc, wc, r, s_ref, exp))
            store[tag] = curves_of(mc, np.asarray(r["phi"]), tag)
        print(f"  {tag} h_wall {h:.4f}  tets {mrow['n_tet']:>8}  conv={conv} |R|={res:.2e} "
              f"lim={r.get('n_limited')} flr={r.get('n_floored')}  "
              f"LE_up {row.get('rms_LE_upper', float('nan')):.6f}  ({wall:.0f}s)", flush=True)
        rows.append(row)
        _write(rows)
        if wall > LEG_GATE_S:
            print(f"  ★ leg gate {wall:.0f}s exceeded -- stopping the ray here (kill clause 4)")
            break

    return report(rows, store, exp)


def report(rows, store, exp):
    by = {x["leg"]: x for x in rows}
    print("\n=== G-R: S4 must bit-reproduce C_medium (its h_far 1.80 IS the unclamped default) ===")
    s4 = by.get("S4")
    if s4 is None or not s4["converged"]:
        print("  -> ★ S4 missing or not converged. STOP (kill clause 1).")
        return 1
    print(f"  S4 LE upper {s4['rms_LE_upper']:.6f} vs committed C_medium {GR_REF}  "
          f"tets {s4['n_tet']}")
    if abs(s4["rms_LE_upper"] - GR_REF) > 5e-7:
        print("  -> ★ G-R FAIL: the instrument moved. STOP and align (kill clause 1).")
        return 1
    print("  -> PASS (same leg)")

    #: --- e(h) and d_self on adjacent pairs ------------------------------------------------------
    print("\n=== e(h) per level, and d_self on adjacent pairs (LE upper = binding) ===")
    print(f"  {'level':6}{'h_wall':>9}{'tets':>9}{'e(LE_up)':>11}{'d_self':>11}{'valid':>8}"
          f"{'AR':>8}{'near %':>8}")
    seq = [t for t, _h in RAY if t in store]
    pairs, elines = [], {}
    for i, tag in enumerate(seq):
        v, ce, ok = sampled(store[tag], exp, "LE", "upper")
        e = rms(v - ce) if ok and v is not None else None
        elines[tag] = (e, ok)
        d = dv = None
        #: ★ DEFECT FIXED: pair only levels ADJACENT ON THE RAY. Using consecutive entries of `seq`
        #: silently promoted a factor-2 gap to an "adjacent" pair whenever a level failed to
        #: converge (S1 did), which would have mixed sqrt(2) and 2 steps in one d_self chain and
        #: corrupted both the monotonicity read and the observed order. Not a criterion change.
        ray = [t for t, _h in RAY]
        j = ray.index(tag)
        nxt = ray[j + 1] if j + 1 < len(ray) else None
        if nxt is not None and nxt in store:
            a, _c1, ok1 = sampled(store[tag], exp, "LE", "upper")
            b, _c2, ok2 = sampled(store[nxt], exp, "LE", "upper")
            if a is not None and b is not None and len(a) == len(b):
                d, dv = rms(a - b), bool(ok1 and ok2)
                pairs.append((tag, nxt, d, dv))
        x = by[tag]
        print(f"  {tag:6}{x['h_wall']:>9.4f}{x['n_tet']:>9}"
              f"{(e if e is not None else float('nan')):>11.6f}"
              f"{(d if d is not None else float('nan')):>11.6f}"
              f"{str(dv) if d is not None else '-':>8}{x['ar_max']:>8.2f}{x['near_pct']:>8.2f}")

    good = [p for p in pairs if p[3]]
    if len(good) < 2:
        print(f"\n  -> ★ only {len(good)} valid adjacent pair(s): UNDEFINED, no direction claimed "
              f"(kill clause 2/3).")
        return 0
    ds = [p[2] for p in good]
    mono = all(ds[i] > ds[i + 1] for i in range(len(ds) - 1))
    e_fine, ok_fine = elines[seq[-1]]
    R = ds[-1] / e_fine if e_fine else float("nan")
    print(f"\n=== the verdict quantity: R = d_self(finest pair) / e(finest) ===")
    print(f"  d_self chain: " + " -> ".join(f"{d:.6f}" for d in ds)
          + f"   monotone shrinking: {mono}")
    print(f"  e(finest) = {e_fine:.6f}   R = {R:.4f}")
    if len({round(d, 12) for d in ds}) == 1:
        print("  -> ★ d_self identical across pairs: SUSPICION, not a verdict (kill clause 5).")
        return 1
    if R <= MODEL_MAX and mono:
        print(f"  -> ★★★ N-MODEL  R = {R:.4f} <= 1/3 and d_self shrinks monotonically ⇒")
        print("     the solution has essentially CONVERGED, and not to the experiment.")
        print(f"     At least {100 * (1 - R):.0f} % of the residual is MODEL error ⇒ M3's LE target")
        print("     needs a model-error band (the way every Track V gate carries the A4 band).")
    elif R >= DISC_MIN or not mono:
        print(f"  -> ★★ N-DISC  R = {R:.4f}"
              + ("" if mono else " and d_self does NOT shrink monotonically")
              + " ⇒ the solution has NOT converged;")
        print("     the root cause is on the DISCRETISATION side. This OVERTURNS 'changing the mesh")
        print("     cannot solve M3' -- the stronger result, and it must be reported as such.")
    else:
        print(f"  -> N-MIX  1/3 < R = {R:.4f} < 1 ⇒ RECORDED, no direction claimed.")

    #: --- Richardson limit, only if d_self is monotone --------------------------------------------
    print("\n=== Richardson limit (reported ONLY if d_self is monotone) ===")
    if mono and len(good) >= 2:
        p = math.log(ds[-2] / ds[-1]) / math.log(math.sqrt(2.0))
        a, _c, _ok = sampled(store[seq[-2]], exp, "LE", "upper")
        b, ce, _ok2 = sampled(store[seq[-1]], exp, "LE", "upper")
        f = 2.0 ** (p / 2.0)
        lim = b + (b - a) / (f - 1.0) if abs(f - 1.0) > 1e-9 else b
        print(f"  observed order p = {p:.3f} (DERIVED, not a criterion)")
        print(f"  |Cp_inf - Cp_exp| RMS = {rms(lim - ce):.6f}   vs S2's target 0.08")
        print(f"  ⇒ {'ABOVE the target -- unreachable by refinement on this ray' if rms(lim - ce) > 0.08 else 'below the target'}")
    else:
        print("  UNDEFINED (d_self not monotone) -- extrapolating a non-monotone sequence would be")
        print("  a fit, not a limit.")

    #: --- the free by-product: what the h_far clamp does to the LE band ---------------------------
    print("\n=== RECORDED (free): the h_far CLAMP's effect on the LE band, first quantification ===")
    try:
        prev = {r["leg"]: r for r in csv.DictReader(open(PREV))}
        cc = prev["C_coarse"]
        s2 = by.get("S2")
        if s2 and s2["converged"]:
            print(f"  same h_wall 0.030: h_far 2.50 (CLAMPED, C_coarse) LE_up "
                  f"{float(cc['rms_LE_upper']):.6f}, {cc['n_tet']} tets")
            print(f"                     h_far 3.60 (on-ray, S2)      LE_up "
                  f"{s2['rms_LE_upper']:.6f}, {s2['n_tet']} tets")
            print(f"  ⇒ the clamp changes the LE-band RMS by "
                  f"{100 * (float(cc['rms_LE_upper']) - s2['rms_LE_upper']) / s2['rms_LE_upper']:+.1f} % "
                  f"-- and it makes coarse LOOK better/worse accordingly")
    except Exception as exc:                                          # noqa: BLE001
        print(f"  (previous CSV unavailable: {exc})")

    print("\n=== other bands (RECORDED) ===")
    for band in ("LE", "MID", "TE"):
        for side in ("upper", "lower"):
            line = []
            for i, tag in enumerate(seq):
                if i + 1 >= len(seq):
                    continue
                a, ca, o1 = sampled(store[tag], exp, band, side)
                b, _cb, o2 = sampled(store[seq[i + 1]], exp, band, side)
                if a is None or b is None or len(a) != len(b):
                    continue
                line.append(f"{tag}->{seq[i+1]} d {rms(a - b):.4f}"
                            + ("" if o1 and o2 else "(INVALID)"))
            ef, okf = None, True
            v, ce, okf = sampled(store[seq[-1]], exp, band, side)
            if v is not None:
                ef = rms(v - ce)
            print(f"  {band:4}/{side:6} e(finest) {(ef if ef else float('nan')):.4f}"
                  + ("" if okf else " (INVALID)") + "   " + "  ".join(line))

    print("\n=== G-S: the solution moves along the ray (RECORDED) ===")
    for tag in seq:
        x = by[tag]
        print(f"  {tag} M_max {x['m_max']:.4f} (inboard {x['m_max_inboard']:.4f})  "
              f"m1_max {x.get('m1_max')}  sigma_min {x.get('sigma_min')}  cl_p {x['cl_p']:.6f}")
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


if __name__ == "__main__":
    sys.exit(main())
