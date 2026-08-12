"""N1 confirm the near-body axis IN ADVANCE / N2 the flat-cap hypothesis / N3 price the MID-for-LE trade.

Pre-registered in docs/dev_phase_three/20260814-1400-nearbody-axis-prereg.md.

Last round found LE-band error monotone in the near-body cell share -- but on an axis identified AFTER
the data, because the arm registered as "bulk-heavy" had actually moved cells outward (the per-shell
count guard caught it). So the axis is registered here first and swept with ONE independent variable,
h_edge, with h_wall searched to hold the cell budget and h_te PINNED (CLAUDE.md: h_edge sizes the LE and
the TE together, so a free h_te reads the P13 tip trailing-edge singularity instead of the LE band).

★★★ The near-body share s and the cell count are both OUTPUTS of the mesh generator. The legs are defined
by knobs; the criterion is a relation between two MEASURED quantities. Nothing here sets s.

Machinery is imported from run_task3_fixed_budget rather than re-typed, so the mesh builder, the guards
and the recipe cannot drift between the two rounds -- and D2/D4 reuse cached meshes by parameter key,
which makes D2 a bit-identical reproduction of last round's control leg (guard G-R).

Outputs (TRACKED): bench/gate_results/task3_nearbody_axis.csv
"""

import csv
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
from pyfp3d.post.surface import planform_area                         # noqa: E402
from run_m3_budget import (ALPHA, BANDS, ETAS, N_UNMASKED,            # noqa: E402
                           parse_experiment, solve)
#: imported, never re-typed -- the builder, the guards and the search are the previous round's
import run_task3_fixed_budget as FB                                   # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_nearbody_axis.csv")
PREV_CSV = os.path.join(HERE, "gate_results", "task3_fixed_budget.csv")
#: G-R: D2 is last round's control leg by construction (same knobs, same mesh key, same recipe)
CTRL_LE_REF, CTRL_TETS = 0.372531, 68624
#: N1: one independent variable. (tag, h_edge, h_wall seed or None = pinned 0.030 with no search)
DOSES = (("D1_edge0200", 0.0200, 0.028),
         ("D2_edge0150", 0.0150, None),
         ("D3_edge0100", 0.0100, 0.032),
         ("D4_edge0075", 0.0075, 0.034),
         ("D5_edge0050", 0.0050, 0.040))
H_TE = 0.015
#: N2: GS2.1's own LE-only arm, cap changed and nothing else. Control = C_medium (previous CSV).
N2_LEG = ("N2_le_only_round", 0.015, 0.00375, 0.0075)
GS21_FLAT_REDUCTION_PCT = 8.3
N2A_MIN, N2B_MAX = 15.0, 8.3
SPREAD_MIN_FRAC = 0.10
MIN_LEGS = 4
TOTAL_GATE_S = 3600.0


def band_counts(exp):
    """n per (eta, band, side) -- a property of the EXPERIMENT file and the band edges only, not of
    any solution, so the pooled all-band RMS can be assembled exactly from the per-band RMS values."""
    out = {}
    for eta in ETAS[:N_UNMASKED]:
        e = exp[eta]
        for name, lo, hi in BANDS:
            for side in ("upper", "lower"):
                m = (e["upper"] == (side == "upper")) & (e["x"] >= lo) & (e["x"] < hi)
                out[f"{name}_{side}"] = out.get(f"{name}_{side}", 0) + int(m.sum())
    return out


def pooled_all(row, counts):
    """N3's key quantity: if pooled does NOT improve while LE does, the LE gain is cells being
    shuffled between bands rather than accuracy being bought."""
    ss = nn = 0.0
    for k, n in counts.items():
        v = row.get(f"rms_{k}")
        if v is None or n == 0:
            continue
        ss += float(v) ** 2 * n
        nn += n
    return float(np.sqrt(ss / nn)) if nn else None


def near_share(mc):
    _h, r = FB.cell_sizes(mc)
    return 100.0 * float((r < 1.0).sum()) / len(mc.elements)


def run(tag, h_wall, h_edge, h_te, exp, counts, elapsed, budget_ok=None):
    if elapsed[0] > TOTAL_GATE_S:
        print(f"  ★ total gate exceeded -- {tag} NOT run (kill clause 4)")
        return None
    mc, wc = cut_wake(FB.cached_mesh(tag, h_wall, h_edge, h_te, None))
    mrow = FB.mesh_row(mc)
    s = near_share(mc)
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    t0 = time.perf_counter()
    r = solve(mc, wc, entropy=True, kutta="pressure", taper=True, probe_seed=0, taper_rc=0.05)
    wall = time.perf_counter() - t0
    elapsed[0] += wall
    conv = bool(r.get("converged"))
    res = float(r.get("residual_history", [float("nan")])[-1])
    row = dict(leg=tag, tip_cap=FB.TIP_CAP, h_wall=h_wall, h_edge=h_edge, h_te=h_te,
               near_pct=s, converged=conv, res_final=res, n_limited=r.get("n_limited"),
               n_floored=r.get("n_floored"), solve_s=round(wall, 1), budget_met=budget_ok, **mrow)
    if conv:
        row.update(FB.flow_row(mc, wc, r, s_ref, exp))
        row["rms_pooled"] = pooled_all(row, counts)
    print(f"  {tag:16} tets {mrow['n_tet']:>7} near {s:5.2f} %  conv={conv} |R|={res:.2e} "
          f"lim={r.get('n_limited')} flr={r.get('n_floored')}  "
          f"LE_up {row.get('rms_LE_upper', float('nan')):.6f}  ({wall:.0f}s)", flush=True)
    return row


def main():
    print("resolved threads: " + ", ".join(
        f"{k}={os.environ.get(k)}" for k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS",
                                             "OPENBLAS_NUM_THREADS")))
    print(f"load average: {os.getloadavg()}")
    print(f"★ G-C: tip_cap = {FB.TIP_CAP} on every leg\n")
    exp = parse_experiment()
    counts = band_counts(exp)
    print(f"  band counts (from the experiment file, solution-independent): {counts}\n")
    rows, elapsed = [], [0.0]

    print("=== N1: one independent variable (h_edge) at the fixed budget ===")
    for tag, h_edge, seed in DOSES:
        if seed is None:
            h_wall, ok = 0.030, True                     # D2 = last round's control, no search
        else:
            h_wall, _he, n, ok, _hist = FB.search_budget(
                tag, None, f"abs:{h_edge}", seed, CTRL_TETS, H_TE)
            if not ok:
                print(f"  ★ {tag}: budget NOT met after {FB.MAX_GENS} generations "
                      f"({n} vs {CTRL_TETS}) -- recorded (kill clause 3)")
        row = run(tag, h_wall, h_edge, H_TE, exp, counts, elapsed, ok)
        if row is not None:
            rows.append(row)
            _write(rows)

    #: --- G-R before any axis is read ------------------------------------------------------------
    d2 = next((x for x in rows if x["leg"].startswith("D2")), None)
    print("\n=== G-R: D2 must bit-reproduce last round's control leg ===")
    if d2 is None or not d2["converged"]:
        print("  -> ★ D2 missing or not converged. STOP (kill clause 1).")
        return 1
    print(f"  D2 LE upper {d2['rms_LE_upper']:.6f} vs committed {CTRL_LE_REF}   "
          f"tets {d2['n_tet']} vs {CTRL_TETS}")
    if abs(d2["rms_LE_upper"] - CTRL_LE_REF) > 5e-7 or d2["n_tet"] != CTRL_TETS:
        print("  -> ★ G-R FAIL: the instrument moved. STOP and align (kill clause 1).")
        _write(rows)
        return 1
    print("  -> PASS (same leg)")

    print("\n=== N2: GS2.1's LE-only arm, ROUND cap, no budget constraint ===")
    tag, hw, he, hte = N2_LEG
    row = run(tag, hw, he, hte, exp, counts, elapsed, None)
    if row is not None:
        rows.append(row)
    _write(rows)
    return report(rows, d2, counts)


def report(rows, d2, counts):
    doses = [x for x in rows if x["leg"].startswith("D")]
    good = [x for x in doses if x["converged"]
            and x["ar_max"] <= FB.AR_FACTOR * d2["ar_max"]
            and x["dih_min"] >= FB.DIH_FACTOR * d2["dih_min"]]
    excl = [x for x in doses if x not in good]

    print("\n=== N1: the axis, read on the pre-registered relation ===")
    print(f"  {'leg':16}{'h_edge':>9}{'h_wall':>9}{'tets':>8}{'budget':>9}{'near %':>9}"
          f"{'LE_upper':>11}{'AR':>8}{'dih':>7}{'quality':>10}")
    for x in sorted(doses, key=lambda y: y["near_pct"]):
        q = x in good
        print(f"  {x['leg']:16}{x['h_edge']:>9.4f}{x['h_wall']:>9.5f}{x['n_tet']:>8}"
              f"{100 * (x['n_tet'] / CTRL_TETS - 1):>+8.1f}%{x['near_pct']:>9.2f}"
              f"{x.get('rms_LE_upper', float('nan')):>11.6f}{x['ar_max']:>8.2f}"
              f"{x['dih_min']:>7.3f}{'OK' if q else '★EXCLUDED':>10}")
    for x in excl:
        print(f"  ★ excluded: {x['leg']} (conv={x['converged']}, AR {x['ar_max']:.2f} vs "
              f"{FB.AR_FACTOR}x{d2['ar_max']:.2f}, dih {x['dih_min']:.3f} vs "
              f"{FB.DIH_FACTOR}x{d2['dih_min']:.3f}) -- RECORDED as evidence of an ALLOCATION "
              f"CEILING, kept out of the axis reading")

    if len(good) < MIN_LEGS:
        print(f"\n  -> ★ only {len(good)} quality+converged legs (< {MIN_LEGS}): "
              f"N1 UNDEFINED, no direction claimed (kill clause 2).")
    else:
        seq = sorted(good, key=lambda y: y["near_pct"])
        le = [x["rms_LE_upper"] for x in seq]
        dec = all(le[i] > le[i + 1] for i in range(len(le) - 1))
        inc = all(le[i] < le[i + 1] for i in range(len(le) - 1))
        spread = (max(le) - min(le)) / d2["rms_LE_upper"]
        if len(set(round(v, 12) for v in le)) == 1:
            print("\n  -> ★ every leg bit-identical: NOT N1c, SUSPICION (kill clause 5).")
            return 1
        print(f"\n  monotone decreasing in near %: {dec}   increasing: {inc}   "
              f"spread {100 * spread:.1f} % of control (needs >= {100 * SPREAD_MIN_FRAC:.0f} %)")
        if dec and spread >= SPREAD_MIN_FRAC:
            print("  -> ★★ N1a  the axis is CONFIRMED with the criterion fixed in advance ⇒")
            print("     next spend is a fixed-budget graded recipe.")
        elif inc:
            print("  -> ★★ N1b  MONOTONE INCREASING: last round's direction is REVERSED. Report it")
            print("     and stop further allocation spending.")
        else:
            print("  -> N1c  RECORDED, no direction claimed "
                  f"({'non-monotone' if not dec else 'spread too small'}).")

    print("\n=== N2: is 'LE resolution is refuted' a flat-cap artifact? ===")
    n2 = next((x for x in rows if x["leg"].startswith("N2")), None)
    prev = {r["leg"]: r for r in csv.DictReader(open(PREV_CSV))}
    cm = prev.get("C_medium")
    if n2 is None or not n2["converged"] or cm is None:
        print("  -> ★ N2 leg missing or not converged -- UNDEFINED, no claim.")
    else:
        base = float(cm["rms_LE_upper"])
        red = 100.0 * (base - n2["rms_LE_upper"]) / base
        print(f"  control C_medium {base:.6f} (previous round, same script/code/threads)")
        print(f"  N2 (h_edge 0.0075 -> 0.00375, round cap) {n2['rms_LE_upper']:.6f}  "
              f"tets {n2['n_tet']} (vs {cm['n_tet']})")
        print(f"  reduction {red:+.1f} %   vs GS2.1's flat-cap {-GS21_FLAT_REDUCTION_PCT:+.1f} %")
        if red >= N2A_MIN:
            print("  -> ★ N2a supports 'the refutation was a flat-cap artifact'")
        elif red <= N2B_MAX:
            print("  -> ★★ N2b does NOT support it ⇒ A3's gain came from the BUDGET COMPENSATION,")
            print("     not from the cap. The attribution has to be rewritten that way.")
        else:
            print("  -> N2c RECORDED, no claim")
        print("  ★ one arm measures MAGNITUDE only -- GS2.1 had several legs and read monotonicity;")
        print("    this cannot and does not test that (pre-registered).")

    print("\n=== N3 (RECORDED): is the LE gain accuracy, or cells shuffled between bands? ===")
    print(f"  {'leg':16}{'near %':>9}{'LE_upper':>11}{'LE_lower':>11}{'MID_upper':>11}"
          f"{'TE_upper':>11}{'POOLED':>11}")
    for x in sorted([y for y in rows if y["converged"]], key=lambda y: y["near_pct"]):
        print(f"  {x['leg']:16}{x['near_pct']:>9.2f}{x['rms_LE_upper']:>11.6f}"
              f"{x['rms_LE_lower']:>11.6f}{x['rms_MID_upper']:>11.6f}"
              f"{x['rms_TE_upper']:>11.6f}{x['rms_pooled']:>11.6f}")
    if len(good) >= 2:
        seq = sorted(good, key=lambda y: y["near_pct"])
        po = [x["rms_pooled"] for x in seq]
        print(f"\n  pooled monotone decreasing with near %: "
              f"{all(po[i] > po[i + 1] for i in range(len(po) - 1))}   "
              f"({po[0]:.6f} -> {po[-1]:.6f}, {100 * (po[-1] - po[0]) / po[0]:+.1f} %)")
        print("  ★ if pooled does NOT improve while LE does, the LE gain is a RESHUFFLE between")
        print("    bands, not accuracy bought -- this sentence belongs next to N1's verdict.")

    print("\n=== G-S: the solution moves too (RECORDED) ===")
    for x in sorted([y for y in rows if y["converged"]], key=lambda y: y["near_pct"]):
        print(f"  {x['leg']:16} M_max {x['m_max']:.4f} (inboard {x['m_max_inboard']:.4f})  "
              f"m1_max {x.get('m1_max')}  sigma_min {x.get('sigma_min')}  cl_p {x['cl_p']:.6f}")
    print("\n=== G-M: LE spacing per leg ===")
    for x in sorted(rows, key=lambda y: y["near_pct"]):
        print(f"  {x['leg']:16} LE h_t {x['le_ht']:.6f}  h_n {x['le_hn']:.6f}  "
              f"aniso {x['le_aniso']:.3f}")
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
