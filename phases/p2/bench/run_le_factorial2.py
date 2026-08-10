"""LE-mechanism factorial, second version: DISPLACEMENT at M0.75, with `h_far` as the
bulk arm.

Pre-registered in docs/dev_phase_two/20260806-0000-le-mechanism-prereg2.md and its two
addenda, all committed before this file was written. What each addendum changed:

  #1 (user)  M0.50 downgrades the QUESTION, not the condition -- with no supersonic cell
             the upwinding is inert and the entropy correction is bit-identical ON/OFF,
             so a subsonic run measures the pure P1 wall error P11 already
             characterised. The condition must carry a shock, and it must be MEASURED
             rather than assumed (NEWTON_M6_RECIPE's claim that M0.70 is subcritical
             turned out to be false: M_max 1.5358, 214 shock cells).
  #2         step 0 measured the window as {0.70, 0.75}; M0.75 is chosen (557 shock
             cells, shock at x/c 0.15, refined leg still 0/0 clamps). And h_far is NOT a
             clean bulk-only knob either (G1 fails at 1e-12: LE spacing moves 0.37 %),
             so instead of relaxing that threshold the contamination gets a measurable
             BOUND from the LE-only arm's own 26.2 % change -- 70.8x larger.

The measure is DISPLACEMENT against the baseline leg, not error against the experiment:
this round asks which refinement direction MOVES the LE band, which needs no experiment
and sidesteps scoring an inviscid solver against a viscous one. The interpolation,
stations and band split are run_m3_budget's, so displacement and error are computed at
the same points by the same code.

Outputs (TRACKED): bench/gate_results/le_factorial2.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.shock import shock_report                          # noqa: E402
from run_le_response import le_geometry                             # noqa: E402
from run_le_factorial import le_face_count                          # noqa: E402
from run_le_window import solve_at                                  # noqa: E402
from run_m3_budget import BANDS, ETAS, N_UNMASKED, parse_experiment  # noqa: E402
from run_seed_exposure import clamp_map                             # noqa: E402

#: ★ addendum #3: `PYFP3D_LE_HFAR` selects the far-field step of the refined arms, so
#: the third round reuses this script rather than forking it. 1.2 was round two (whose
#: all-scales leg hit m_cap); 1.8 is round three's milder step. The CRITERION is
#: untouched -- only the leg moves -- and BOTH refined arms move together, because
#: comparing a 1.2 displacement against a 1.8 one would be mixing step sizes.
H_FAR_FINE = float(os.environ.get("PYFP3D_LE_HFAR", "1.8"))
CSV = os.path.join(_GATE,
                   f"le_factorial2_hfar{H_FAR_FINE:g}.csv".replace(".", "p", 1)
                   if H_FAR_FINE != 1.2 else "le_factorial2.csv")
M_INF = 0.75                       # addendum #2: measured as the highest in-window Mach
#: (tag, h_le, h_far). h_wall is 0.020 on every leg -- `build` fixes it -- so the only
#: things that move are the LE spacing and the far field.
LEGS = (("G00_base",     0.010,  2.4),
        ("G10_faronly",  0.010,  H_FAR_FINE),
        ("G01_leonly",   0.0075, 2.4),
        ("G11_both",     0.0075, H_FAR_FINE))
LEG_GATE_S = 500.0                 # addendum #2 (c), from step 0's measured 320-443 s
TOTAL_GATE_S = 2000.0


def cp_at_exp_points(mc, phi, exp):
    """Each station's computed Cp sampled at the EXPERIMENT's x locations, per side.

    Same interpolation and same side mapping as run_m3_budget.band_rms -- reused so that
    a displacement and an error are computed at identical points, which is the only way
    the two can be compared later.
    """
    out = {}
    for eta in ETAS[:N_UNMASKED]:
        c = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
        e = exp[eta]
        for want_upper in (True, False):
            side = "upper" if want_upper else "lower"
            m = e["upper"] == want_upper
            if not np.any(m):
                continue
            out[(eta, side)] = (e["x"][m],
                                np.interp(e["x"][m], c[f"x_{side}"],
                                          c[f"cp_{side}"]))
    return out


def band_err(cp_leg, exp):
    """RMS of (Cp_leg - Cp_EXPERIMENT) per band -- addendum #4's RECORDED quantity.

    Displacement (band_disp) says how much a refinement MOVES the solution; this says
    whether it moves toward the experiment. The distinction is not academic here: GS2.1
    measured LE refinement changing the ERROR by only -8.3 % while this round measures
    the same direction DISPLACING the solution by 0.14332. Solution movement and error
    reduction are demonstrably different quantities on this band.
    """
    acc = {}
    for (eta, side), (x, cp) in cp_leg.items():
        e = exp[eta]
        m = e["upper"] == (side == "upper")
        d = cp - e["cp"][m]
        for name, lo, hi in BANDS:
            b = (x >= lo) & (x < hi)
            ss, nn = acc.get(f"{name}_{side}", (0.0, 0))
            acc[f"{name}_{side}"] = (ss + float(np.sum(d[b] ** 2)), nn + int(b.sum()))
    out = {}
    for name, _, _ in BANDS:
        for side in ("upper", "lower"):
            ss, nn = acc.get(f"{name}_{side}", (0.0, 0))
            out[f"{name}_{side}"] = (ss / nn) ** 0.5 if nn else float("nan")
    return out


def band_disp(cp_leg, cp_base):
    """RMS of (Cp_leg - Cp_base) per band and side, pooled over the unmasked stations."""
    acc = {}
    for (eta, side), (x, cp) in cp_leg.items():
        xb, cpb = cp_base[(eta, side)]
        assert np.array_equal(x, xb), "the two legs are not sampled at the same x"
        d = cp - cpb
        for name, lo, hi in BANDS:
            m = (x >= lo) & (x < hi)
            ss, nn = acc.get(f"{name}_{side}", (0.0, 0))
            acc[f"{name}_{side}"] = (ss + float(np.sum(d[m] ** 2)),
                                     nn + int(m.sum()))
    out = {}
    for name, _, _ in BANDS:
        for side in ("upper", "lower"):
            ss, nn = acc.get(f"{name}_{side}", (0.0, 0))
            out[f"{name}_{side}"] = (ss / nn) ** 0.5 if nn else float("nan")
    tot_ss = sum(v[0] for v in acc.values())
    tot_n = sum(v[1] for v in acc.values())
    out["pooled"] = (tot_ss / max(tot_n, 1)) ** 0.5
    return out


def main():
    exp = parse_experiment()
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    print(f"LE factorial v2 at M{M_INF} (prereg 20260806-0000 + addenda #1/#2): "
          f"which refinement direction MOVES the LE band?\n")
    rows, cps, spent = [], {}, 0.0
    for tag, h_le, h_far in LEGS:
        t0 = time.perf_counter()
        mc, wc = cut_wake(_build(h_le, h_far))
        ht, hn, aniso = le_geometry(mc)
        n_le, n_wall = le_face_count(mc), len(mc.boundary_faces["wall"])
        t_gen = time.perf_counter() - t0
        t1 = time.perf_counter()
        r = solve_at(mc, wc, M_INF)
        wall = time.perf_counter() - t1
        spent += wall
        hist = np.asarray(r.get("residual_history", []), dtype=float)
        res = float(hist[-1]) if len(hist) else float("nan")
        nlim, nflr = int(r.get("n_limited") or 0), int(r.get("n_floored") or 0)
        usable = bool(r.get("converged")) and res < 1e-9 and nlim == 0 and nflr == 0
        phi = np.asarray(r["phi"])
        cps[tag] = cp_at_exp_points(mc, phi, exp) if usable else None
        loc = clamp_map(mc, phi, M_INF)
        sh = {}
        for eta in (0.44, 0.65, 0.90):
            c = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
            sh[eta] = shock_report(c, M_INF)["upper"].get("x_shock")
        row = dict(tag=tag, h_le=h_le, h_far=h_far, n_tet=len(mc.elements),
                   n_wall=n_wall, n_le_faces=n_le, le_ht=round(ht, 10),
                   converged=bool(r.get("converged")), res_final=res,
                   n_limited=nlim, n_floored=nflr, usable=usable,
                   m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                   n_shock_cells=r.get("n_shock_cells"),
                   fallback_fired=bool(r.get("seed_fallback", {}).get("fired")),
                   **{f"x_shock_{e}": sh[e] for e in sh},
                   t_gen_s=round(t_gen, 1), wall_s=round(wall, 1),
                   clamp_frac_LE=loc["clamp_frac_LE_0_15"],
                   clamp_frac_MID=loc["clamp_frac_MID_15_85"],
                   x_of_peak_q2=loc["x_of_peak_q2"])
        rows.append(row)
        print(f"  {tag:14} h_le {h_le:<7} h_far {h_far:<5} tets {row['n_tet']:>7} "
              f"LE faces {n_le:>5} h_t {ht:.6f}  conv={str(row['converged']):5} "
              f"|R|={res:.1e} lim/flr={nlim}/{nflr} M_max={row['m_max']:.4f} "
              f"shock={row['n_shock_cells']} ({t_gen:.0f}+{wall:.0f}s)", flush=True)
        if not usable:
            print("    ★ NOT USABLE (G2): its displacement is not a reading", flush=True)
        if wall > LEG_GATE_S:
            print(f"  ★ leg gate {wall:.0f} s > {LEG_GATE_S:.0f} s -- stopping",
                  flush=True)
            break
        if spent > TOTAL_GATE_S:
            print(f"  ★ total gate {spent:.0f} s > {TOTAL_GATE_S:.0f} s -- stopping",
                  flush=True)
            break

    base = cps.get("G00_base")
    if base is not None:
        #: addendum #4: the baseline's own error, so every leg's ΔE has a reference
        eb = band_err(base, exp)
        for r in rows:
            if r["tag"] == "G00_base":
                r.update({f"err_{k}": round(v, 6) for k, v in eb.items()})
        print(f"\n  [addendum #4, RECORDED] base LE-upper error vs EXPERIMENT "
              f"{eb['LE_upper']:.5f}  MID_up {eb['MID_upper']:.5f}  "
              f"TE_up {eb['TE_upper']:.5f}")
        # G5: the zero test. Differencing the baseline against itself must be exactly 0.
        z = band_disp(base, base)
        assert all(v == 0.0 for k, v in z.items()), f"G5 zero test failed: {z}"
        print("\n  G5 zero test (base minus itself): exactly 0 -> PASS")
        for row in rows:
            c = cps.get(row["tag"])
            if c is None or row["tag"] == "G00_base":
                continue
            d = band_disp(c, base)
            row.update({f"disp_{k}": round(v, 6) for k, v in d.items()})
            e = band_err(c, exp)
            row.update({f"err_{k}": round(v, 6) for k, v in e.items()})
            print(f"  disp {row['tag']:14} LE_up {d['LE_upper']:.5f} "
                  f"LE_lo {d['LE_lower']:.5f} MID_up {d['MID_upper']:.5f} "
                  f"TE_up {d['TE_upper']:.5f} pooled {d['pooled']:.5f}")
            print(f"  err  {row['tag']:14} LE_up {e['LE_upper']:.5f} "
                  f"(base {eb['LE_upper']:.5f}, "
                  f"{100*(e['LE_upper']-eb['LE_upper'])/eb['LE_upper']:+.1f} %)  "
                  f"MID_up {e['MID_upper']:.5f} TE_up {e['TE_upper']:.5f}")

    keys = sorted({k for r in rows for k in r})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")
    _verdict(rows)
    return 0


def _build(h_le, h_far):
    """run_le_window.build with h_le as a parameter (it fixes H_LE = 0.010)."""
    from pyfp3d.meshgen.wing3d import MAC, onera_m6_wing_mesh
    from run_le_window import H_WALL
    return onera_m6_wing_mesh(h_wall=H_WALL, h_edge=h_le, h_te=0.5 * H_WALL,
                              h_wake=3.0 * H_WALL, h_far=h_far,
                              r_far=15.0 * MAC, tip_cap="round", embed_wake=True)


def _verdict(rows):
    by = {r["tag"]: r for r in rows}
    need = ("G00_base", "G10_faronly", "G01_leonly", "G11_both")
    print("\n=== verdict (criteria fixed in the registration sec 4 + addendum #2b) ===")
    if any(t not in by for t in need):
        print(f"  incomplete: {[t for t in need if t not in by]} missing -- "
              f"nothing read"); return
    bad = [t for t in need if not by[t]["usable"]]
    if bad:
        print(f"  G2 FAIL on {bad} -- their displacements are not readings")
        print("  -> RECORDED, no conclusion"); return
    b = by["G00_base"]
    d10 = by["G10_faronly"].get("disp_LE_upper")
    d01 = by["G01_leonly"].get("disp_LE_upper")
    d11 = by["G11_both"].get("disp_LE_upper")
    if None in (d10, d01, d11):
        print("  displacements missing -- nothing read"); return
    # the far arm's own LE contamination, and the LE arm's, as MEASURED spacings
    c10 = abs(by["G10_faronly"]["le_ht"] - b["le_ht"]) / b["le_ht"]
    c01 = abs(by["G01_leonly"]["le_ht"] - b["le_ht"]) / b["le_ht"]
    print(f"  LE-upper displacement: far-only {d10:.5f} | le-only {d01:.5f} | "
          f"both {d11:.5f}")
    print(f"  LE spacing changed by: far-only {100*c10:.3f} % | "
          f"le-only {100*c01:.3f} %  (ratio {c01/max(c10,1e-12):.1f}x)")
    #: addendum #2b: the bound the far arm's own contamination could explain, by LINEAR
    #: extrapolation from the le-only arm -- an ASSUMPTION, labelled as one, with a 5x
    #: margin. We may never claim the far arm is clean; G1 already failed.
    bound = 5.0 * d01 * (c10 / max(c01, 1e-12))
    print(f"  contamination bound (ASSUMES local linearity, 5x margin): {bound:.5f}")
    ctrl = (by["G10_faronly"].get("disp_MID_upper", 0.0) >= 0.8 * d10
            and by["G10_faronly"].get("disp_TE_upper", 0.0) >= 0.8 * d10)
    if ctrl:
        print("  ★ CONTROL CLAUSE FIRED: far refinement moves MID and TE as much as LE")
        print("    -> a GLOBAL displacement, not LE-specific pollution. P/N NOT read.")
        return
    if abs(d10) >= 0.60 * abs(d11) and abs(d10) > bound:
        print("  => P: the LE deficit is POLLUTION-DOMINATED -- the far field moves the")
        print("     LE band by more than its own LE contamination can explain, and by")
        print("     >= 60 % of what all-scales does. S2's target moves to the bulk/far")
        print("     mesh, and this joins P11's sphere mechanism as one finding.")
    elif abs(d10) < 0.10 * abs(d11):
        print("  => N: far-field pollution REFUTED too. Neither LE-local (v1 measured")
        print("     -1.55 %) nor far-field => a third mechanism, own round.")
    else:
        print("  => I: between the thresholds -> RECORDED, no conclusion")


if __name__ == "__main__":
    sys.exit(main())
