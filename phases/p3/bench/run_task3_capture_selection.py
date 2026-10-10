"""Change the CAPTURE without changing the MAGNITUDE: SELECT a computed sigma field, never blend one.

★★★ NOT RUNNABLE AT HEAD (2026-08-16) -- and that is the intended end state. K3 fired: the selection
rule is HARMFUL (balanced panel [0, 5], cl spread 1.21 % -> 7.47 %, one seed lost), so the route was
killed, and at the phase-3 close-out the `capture_select` / `capture_select_abs` knobs were REMOVED
from pyfp3d/solve/newton.py rather than left as default-OFF options -- this project's own rule is that
a knob kept "in case we need it again" is how temporary knobs become permanent (the same disposal the
temporary sigma_scale instrument got at its registered expiry). This file therefore raises TypeError
on the solve call and is kept as the PROVENANCE of its committed CSV, not as a live harness.

  to reproduce: git checkout c38f9a6 -- pyfp3d/solve/newton.py   (then restore with
                git checkout HEAD -- pyfp3d/solve/newton.py, NEVER the bare form -- see CLAUDE.md)

Pre-registered in phases/p3/docs/dev_phase_three/20260812-2100-capture-selection-prereg.md, with addendum #1
(the relative-delta rule) committed after the smoke test and BEFORE any band was read.

Why selection and not smoothing: averaging sigma fields is NOT magnitude-neutral (the mean of 0.743
and 0.987 is 0.865, weaker than either input) and relaxation is exponentially-weighted averaging, so
both re-import the confound the previous round was spent removing. Adopting one refresh's OWN output
unchanged is neutral by construction, which makes it the only available experiment that varies
capture-uniformity alone.

★ Addendum #1's correction, kept visible here because it is the round's main methodological point:
scoring candidates by the RAW |dsigma| is biased toward weak corrections -- |dsigma| is minimised by
refreshes where sigma barely acts, i.e. the early ones before the shock forms, and it duly picked
refresh 2 (sigma_min 0.98975) landing on x_shock 0.5601 against the isentropic 0.5605. The primary
rule therefore normalises by (1 - sigma_min); the biased rule is kept as a RECORDED control arm so
that "the bias is real" is evidence rather than an assertion.

Bands (unchanged from the registration): K1 = m1_max spread down AND cl spread meets F1' (relative to
each family's own baseline) -> capture-uniformity alone suffices. K2 = m1_max down but cl not ->
magnitude is also involved. K3 = spread >= baseline or fewer seeds converge -> the route dies. K4 ->
RECORDED. G-M = the adopted field must be BIT-IDENTICAL to a candidate. Spreads use a BALANCED PANEL.

Outputs (TRACKED): bench/gate_results/task3_capture_selection.csv
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

CSV = os.path.join(HERE, "gate_results", "task3_capture_selection.csv")
M_INF, ALPHA, C = 0.80, 1.25, 1.5
SEEDS = (0, 5, 12)
#: (label, capture_select, capture_select_abs) -- "abs" is the recorded control arm
ARMS = (("off", False, False), ("rel", True, False), ("abs", True, True))
RECIPE = dict(upwind_c=C, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80)
CASES = (("hybrid_R0", None), ("unstr_coarse", "cases/meshes/naca0012_2.5d/coarse.msh"))
F1_ABS_PCT, F1_REL = 5.0, 1.0 / 3.0
SHOCK_REF, SHOCK_TOL = 0.61, 0.02
LEG_GATE_S = 600.0


def _mesh(path):
    if path is None:
        from run_task3_refinement_paradox import build
        return build(160, 0.004, None, 1.0)[0]
    return cut_wake(read_mesh(os.path.join(REPO, path)))


def run(tag, mc, wc, arm, sel, absr, seed):
    t0 = time.perf_counter()
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, n_picard_seed=seed,
                            capture_select=sel, capture_select_abs=absr, **RECIPE)
    wall = time.perf_counter() - t0
    h = list(r.get("residual_history") or [])
    z = float(np.unique(mc.nodes[:, 2]).mean())
    rep = shock_report(wall_cp_curve(mc, r["phi"], z=z, m_inf=M_INF), M_INF)
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
                                alpha_deg=ALPHA, s_ref=float(np.ptp(mc.nodes[:, 2])), m_inf=M_INF)
    xs = rep["upper"].get("x_shock")
    fr = r.get("sigma_freeze_report") or {}
    hist = r.get("sigma_history") or []
    i = r.get("capture_adopted_index")
    #: ★ the ADOPTED refresh's own diagnostics. The result dict's sigma_min / m1_max are the LAST
    #: refresh's -- reporting those would describe a state the solve did not freeze, which the smoke
    #: test caught before this experiment ran.
    ad = hist[i] if (i is not None and 0 <= i < len(hist)) else None
    ws = r["workspace"]
    gm = (None if not sel or ws.sigma_best is None
          else bool(np.array_equal(ws.sigma_frozen, ws.sigma_best)))
    return dict(case=tag, arm=arm, seed=seed, converged=bool(r["converged"]),
                res_final=(h[-1] if h else None), n_newton=len(h),
                cl_p=float(f["cl"]), x_shock=None if xs is None else float(xs),
                adopted_index=i, adopted_score=r.get("capture_adopted_delta"),
                adopted_sigma_min=(None if ad is None else float(ad[0])),
                adopted_m1_max=(None if ad is None else float(ad[2])),
                last_sigma_min=r.get("sigma_min"), last_m1_max=r.get("m1_max"),
                gm_bit_identical=gm, selection_churn=fr.get("selection_churn"),
                n_limited=int(r.get("n_limited") or 0),
                n_floored=int(r.get("n_floored") or 0), wall_s=round(wall, 1))


def main():
    print(f"capture selection   M{M_INF}/alpha {ALPHA}/C {C}   arms {[a[0] for a in ARMS]}   "
          f"seeds {SEEDS}   threads {os.environ['NUMBA_NUM_THREADS']}\n")
    rows = []
    for tag, path in CASES:
        mc, wc = _mesh(path)
        for arm, sel, absr in ARMS:
            for seed in SEEDS:
                row = run(tag, mc, wc, arm, sel, absr, seed)
                rows.append(row)
                ai = row["adopted_index"]
                print(f"  {tag:13} {arm:4} seed {seed:>2}  conv={str(row['converged']):5} "
                      f"|R|={row['res_final']:.2e} cl_p={row['cl_p']:>9.6f} "
                      f"x_shock={'-' if row['x_shock'] is None else format(row['x_shock'], '.4f')} "
                      f"adopted={str(ai):>4} "
                      f"(sig {'-' if row['adopted_sigma_min'] is None else format(row['adopted_sigma_min'], '.4f')}"
                      f" m1 {'-' if row['adopted_m1_max'] is None else format(row['adopted_m1_max'], '.4f')})"
                      f" G-M={row['gm_bit_identical']} ({row['wall_s']:.0f}s)", flush=True)
                if row["wall_s"] > LEG_GATE_S:
                    print("    ★ cost gate exceeded -- stop"); _write(rows); return 1
    _write(rows)
    return _read(rows)


def _read(rows):
    print("\n=== G-M (magnitude neutrality, INPUT-side): adopted field == a candidate, bitwise ===")
    gms = [r for r in rows if r["gm_bit_identical"] is not None]
    bad = [r for r in gms if not r["gm_bit_identical"]]
    print(f"  {len(gms) - len(bad)}/{len(gms)} adopted fields are bit-identical to their candidate")
    if bad:
        print("  -> ★ G-M FAIL: this was implemented as SMOOTHING, not selection. Kill clause 1.")
        return 1
    print("  -> PASS (nothing was blended, so the correction's magnitude formula is untouched)")

    print("\n=== F3: the off arm must bit-reproduce the legacy answers ===")
    ref = {("hybrid_R0", 0): 0.261367, ("hybrid_R0", 5): 0.258223, ("hybrid_R0", 12): 0.351467,
           ("unstr_coarse", 0): 0.407978, ("unstr_coarse", 5): 0.407108}
    f3 = True
    for (tag, seed), cl_a in ref.items():
        m = [r for r in rows if r["case"] == tag and r["arm"] == "off" and r["seed"] == seed]
        if not m or not m[0]["converged"]:
            print(f"  {tag:13} seed {seed:>2}: not converged"); continue
        hit = abs(m[0]["cl_p"] - cl_a) < 5e-6
        f3 &= hit
        print(f"  {tag:13} seed {seed:>2}: {m[0]['cl_p']:.6f} vs {cl_a}  "
              f"-> {'MATCH' if hit else '★ MISMATCH'}")
    print(f"  -> F3 {'PASS' if f3 else 'FAIL'}")
    if not f3:
        return 1

    print("\n=== K1-K4 on a BALANCED PANEL (seeds converging in BOTH arms) ===")
    for tag, _ in CASES:
        off = {r["seed"]: r for r in rows if r["case"] == tag and r["arm"] == "off"}
        print(f"\n  {tag}")
        for arm in ("rel", "abs"):
            trt = {r["seed"]: r for r in rows if r["case"] == tag and r["arm"] == arm}
            panel = [s for s in SEEDS if off.get(s, {}).get("converged")
                     and trt.get(s, {}).get("converged")]
            if len(panel) < 2:
                print(f"    {arm:4} -> K4  balanced panel has {len(panel)} seed(s) -> UNDEFINED "
                      f"(NOT 'a small spread')")
                continue

            def sp(d, key):
                v = [d[s][key] for s in panel if d[s][key] is not None]
                return (max(v) - min(v)) if len(v) >= 2 else None

            def rel(d):
                v = [abs(d[s]["cl_p"]) for s in panel]
                return 100.0 * (max(v) - min(v)) / max(float(np.mean(v)), 1e-12)

            r_off, r_trt = rel(off), rel(trt)
            #: the mechanism variable: how much the ADOPTED shock strength varies across seeds
            m_off = sp(off, "last_m1_max")
            m_trt = sp(trt, "adopted_m1_max") or sp(trt, "last_m1_max")
            f1p = r_trt <= r_off * F1_REL and r_trt <= F1_ABS_PCT
            churn = any(bool(trt[s]["selection_churn"]) for s in panel)
            n_off = sum(1 for s in SEEDS if off.get(s, {}).get("converged"))
            n_trt = sum(1 for s in SEEDS if trt.get(s, {}).get("converged"))
            inb_off = sum(1 for s in panel if abs(off[s]["x_shock"] - SHOCK_REF) <= SHOCK_TOL)
            inb_trt = sum(1 for s in panel if abs(trt[s]["x_shock"] - SHOCK_REF) <= SHOCK_TOL)
            if r_trt >= r_off or n_trt < n_off:
                band = "K3"
                why = (f"spread {r_trt:.2f} % >= baseline {r_off:.2f} %" if r_trt >= r_off
                       else f"converged {n_trt} < {n_off}")
            elif m_trt is not None and m_off is not None and m_trt < m_off and f1p:
                band, why = "★ K1", f"m1_max spread {m_off:.3f} -> {m_trt:.3f} and F1' met"
            elif m_trt is not None and m_off is not None and m_trt < m_off:
                band, why = "K2", (f"m1_max spread down ({m_off:.3f} -> {m_trt:.3f}) but F1' NOT "
                                   f"met ({r_trt:.2f} % > {r_off * F1_REL:.2f} %)")
            else:
                band, why = "K4", "neither band cleanly"
            print(f"    {arm:4} panel {panel}  cl spread {r_off:.2f} % -> {r_trt:.2f} %  "
                  f"(F1' needs <= {r_off * F1_REL:.2f} % and <= {F1_ABS_PCT} %)")
            print(f"         m1_max spread {'-' if m_off is None else format(m_off, '.3f')} -> "
                  f"{'-' if m_trt is None else format(m_trt, '.3f')}   churn={churn}   "
                  f"F5 in band {inb_off} -> {inb_trt} of {len(panel)}   "
                  f"converged {n_off} -> {n_trt}")
            print(f"         -> {band}   {why}")
    print("\n  ★ the 'abs' arm is the RECORDED control: it is registered as BIASED toward weak")
    print("    corrections, so a better-looking number there is not a better result.")
    return 0


def _write(rows):
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
