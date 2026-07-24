"""V5.1d -- the near-band seed: does a quadratic basin exist ADJACENT to the floor?

Binding text: cases/analysis/v5_1d_near_band_window/PRE_REGISTRATION.md
(committed BEFORE the first execution, per discipline). Design inputs: the
committed GV5.1c VERDICT (../v5_1c_above_band_window/VERDICT.md -- the
above-band window read: NO quadratic regime anywhere above the floor; the
trajectory STALLS mid-range at F_BL ~ 1e-2, 3-4 decades above the floor,
never reaching the band) and the IBL-floor diagnosis
(../v5_ibl_floor/results/findings.md).

Protocol (the GV5.1c protocol verbatim, NEW near-band windows):
  base seed = the GV5.1 amended protocol verbatim (the GV5.1b runner's
  build_loose_state IMPORTED via the GV5.1c runner, not mirrored --
  HEAD-regen loose-converged state, wiring guard |dcl_k0| <= 1e-8); the
  tight pack is built on the UNPERTURBED state (every frozen operator
  identical to GV5.1b's / GV5.1c's); perturbed seed = U[:,0] * (1+eps)
  at the FREE (non-inflow) BL nodes, eps calibrated by the GV5.1c
  deterministic log10-bisection (IMPORTED) into the pre-registered
  NEAR-BAND windows T1=[1e-4,1e-3] (primary; 5.8-58x the medium floor
  band / 3.2-31.6x the coarse band -- the decades directly above the
  band, BELOW the GV5.1c mid-range stall ~1e-2) / T2=[1e-3,1e-2]
  (escalation, fires ONLY if the T1 trajectory yields <3 above-band
  contraction triples; one re-seed max per level);
  solve = the GV5.1b machinery verbatim: newton_tight(scaling="rowcol",
  lm_damping=True, floor_stop=True; max_iter=10, tol=1e-8,
  tol_abs=1e-10, line_search=True). Assembly untouched; no solver-side
  edits of any kind.

Bands (the pre-registered ones):
  (a) implementation exactness: the committed suite green (the tight
      fleet 28 + tests/test_v5_above_band_seed.py 9 + the new synthetic
      tests test_v5_near_band_seed.py) + live identities on the
      PERTURBED-seed J (T1): e1 <= 1e-12; e2 <= max(1e-10,
      10*kappa_1(J)*eps) -- the cond-aware tolerance carried in (the
      GV5.1b 2026-07-24 adjudication); e3 <= 1e-6.
  (b) the near-band quadratic-basin window (medium binding, coarse
      recorded): p on F_BL max-norm contraction triples with both
      predecessors above the floor band (coarse 3.16e-5 / medium
      1.71e-5, the GV5.1b bands verbatim); PASS = >= 3 triples (T1;
      pooled with T2 if the escalation fires) with median p in
      [1.5, 2.5]; converged outright -> trivially PASS (kept verbatim
      for protocol identity; with the floor in place the expected
      terminations are floor_reached / cap); < 3 triples after the
      escalation -> the pre-registered RECORDED fallback (available
      triples + contraction factors + the log-log regression slope
      quoted; NOT a gate crash). BAND-ENTRY iteration elevated to a
      key datum: GV5.1c never entered the band from above.
  (c) counts (recorded; aspirational N_polish <= 2x loose; the seeds
      are deliberately off-point -- N_polish measures the near-band
      traversal, not production convergence).

Environment: executed under a TEMPORARY 8-thread constraint (this
session only, user-directed 2026-07-24): the runner keeps the 16
defaults below, the constraint is applied via the environment at
execution time and recorded in the artifacts; wall times are NOT
comparable to the 16-thread ledger entries.

Run:  python cases/analysis/v5_1d_near_band_window/run.py
      python cases/analysis/v5_1d_near_band_window/run.py --levels coarse
(partial runs merge summary.csv by (band, level); a full run
regenerates every artifact from scratch)
"""

import argparse
import csv
import importlib.util
import os
import sys
import time

# resolve pyfp3d from THIS worktree (the site-packages editable install may
# point at a sibling worktree)
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..")))

# runner defaults (the ledger-standard 16); the temporary session
# constraint (8) is applied via the environment, never here
os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np

from pyfp3d.viscous import tight_driver as td
from tests.v5_state import M_CAP, M_CRIT, RHO_FLOOR, UPWIND_C

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

GV51C_RUN = os.path.join(HERE, "..", "v5_1c_above_band_window", "run.py")
GV51C_RESULTS = os.path.join(HERE, "..", "v5_1c_above_band_window",
                             "results")


def _load_gv51c():
    """The COMMITTED GV5.1c runner module -- imported, not mirrored, so
    the seed/calibration helpers (and through it the GV5.1b protocol
    functions) cannot drift."""
    spec = importlib.util.spec_from_file_location("gv51c_run", GV51C_RUN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gc = _load_gv51c()
gb = gc.gb                            # the GV5.1b runner module

FLOOR_BAND = gb.FLOOR_BAND            # coarse 3.16e-5 / medium 1.71e-5
TARGETS = {"T1": (1.0e-4, 1.0e-3), "T2": (1.0e-3, 1.0e-2)}  # NEAR-BAND
BAND_A_TOL = gb.BAND_A_TOL            # (1e-12, 1e-10, 1e-6)
THREADS = os.environ.get("NUMBA_NUM_THREADS", "16")

SUMMARY = []                          # (band, level, key_result, verdict)


# ---------------------------------------------------------------------------
# summary/compare writers (the v5_tight_coupling merge discipline, runner-local)
# ---------------------------------------------------------------------------


def _write_lines(name, header, lines):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as f:
        f.write(header + "\n")
        for ln in lines:
            f.write(ln + "\n")
    print(f"  wrote {path}", flush=True)


def _write_summary():
    rows = {}
    path = os.path.join(RESULTS, "summary.csv")
    if os.path.exists(path):
        with open(path) as f:
            lines = [ln.rstrip("\n") for ln in f if ln.strip()]
        for ln in lines[1:]:
            flds = ln.split(",", 2)
            if len(flds) >= 2:
                rows[(flds[0], flds[1])] = ln
    for t in SUMMARY:
        rows[(t[0], t[1])] = ",".join(str(x) for x in t)
    order = {"(a)": 0, "(b)": 1, "(c)": 2, "calib": 3, "wall": 4}
    out = sorted(rows.values(),
                 key=lambda ln: (order.get(ln.split(",", 1)[0], 99),
                                 ln.split(",", 2)[1]))
    _write_lines("summary.csv", "band,level,key_result,verdict", out)


def _record(band, level, key_result, verdict):
    SUMMARY.append((band, level, key_result, verdict))
    print(f"  [{verdict}] {band} {level}: {key_result}", flush=True)
    _write_summary()


def _merge_compare(level, rows):
    """compare.csv: one row per (level, metric); fresh rows replace the
    same pair, other levels' rows kept."""
    header = "level,metric,gv5_1d,committed_reference,note"
    kept = {}
    path = os.path.join(RESULTS, "compare.csv")
    if os.path.exists(path):
        with open(path) as f:
            lines = [ln.rstrip("\n") for ln in f if ln.strip()]
        for ln in lines[1:]:
            flds = ln.split(",", 2)
            if len(flds) >= 2:
                kept[(flds[0], flds[1])] = ln
    for r in rows:
        kept[(level, str(r[0]))] = ",".join(str(x) for x in (level,) + r)
    order = {"coarse": 0, "medium": 1}
    out = sorted(kept.values(),
                 key=lambda ln: (order.get(ln.split(",", 1)[0], 9),
                                 ln.split(",", 2)[1]))
    _write_lines("compare.csv", header, out)


# ---------------------------------------------------------------------------
# the near-band specific read: band entry
# ---------------------------------------------------------------------------


def band_entry_iter(e_seq, band):
    """First iteration index whose F_BL max-norm sits AT/BELOW the floor
    band, or None -- the key new datum of this gate (the GV5.1c
    above-band trajectories never entered the band)."""
    return next((i for i, v in enumerate(e_seq) if v <= band), None)


# ---------------------------------------------------------------------------
# live reads on the pack (CSV history writer, runner-local; the physics
# helpers -- fbl_of / band_a_live_at / calibrate_eps / perturb_delta /
# regression_slope / contraction_factors / window_verdict -- are the
# GV5.1c runner's, IMPORTED)
# ---------------------------------------------------------------------------


def _history_final(results_dir, level):
    """A committed final row (f_bl_max, merit) -- read, never recomputed."""
    with open(os.path.join(
            results_dir, f"newton_history_{level}.csv")) as f:
        last = list(csv.DictReader(f))[-1]
    return float(last["f_bl_max"]), float(last["merit"])


def _write_history(tag, hist, band, term):
    """newton_history_{tag}.csv -- the GV5.1b columns verbatim."""
    e_bl = [float(h["block_max"][2]) for h in hist]
    m_seq = [float(h["merit"]) for h in hist]
    p_bl = dict(gb._p_series(e_bl, band=band))
    p_m = dict(gb._p_series(m_seq))
    lines = []
    for i, h in enumerate(hist):
        bm = h["block_max"]
        lines.append(",".join([
            str(h["iter"]),
            f"{bm[0]:.6e}", f"{bm[1]:.6e}", f"{bm[2]:.6e}",
            f"{h['merit']:.6e}",
            "" if h.get("mu") is None else f"{h['mu']:.6e}",
            "" if h["lam"] is None else f"{h['lam']:.6f}",
            "" if i not in p_bl else f"{p_bl[i]:.4f}",
            "" if i not in p_m else f"{p_m[i]:.4f}",
            f"{h['ds_change']:.6e}",
            term if i == len(hist) - 1 else "",
            f"{h['wall_s']:.2f}"]))
    _write_lines(
        f"newton_history_{tag}.csv",
        "iter,f_phi_max,f_gamma_max,f_bl_max,merit,mu,lam,p_fbl,"
        "p_merit,ds_change,termination,wall_s", lines)


def _polish(pack, x0, tag):
    """One GV5.1b-machinery polish from x0 (the pre-registered flags
    verbatim); returns (res_n, wall_s)."""
    print(f"    polish [{tag}]: newton_tight(scaling=rowcol, "
          "lm_damping, floor_stop; max_iter=10) ---", flush=True)
    t0 = time.perf_counter()
    res_n = td.newton_tight(pack, x0=x0, max_iter=10, tol=1e-8,
                            tol_abs=1e-10, line_search=True, verbose=True,
                            scaling="rowcol", lm_damping=True,
                            floor_stop=True)
    wall = time.perf_counter() - t0
    print(f"    polish [{tag}] in {wall:.0f}s: converged="
          f"{res_n['converged']} termination={res_n['termination']} "
          f"n_iter={res_n['n_iter']}", flush=True)
    _write_history(tag, res_n["history"], FLOOR_BAND[tag.split("_")[0]],
                   res_n["termination"])
    return res_n, wall


# ---------------------------------------------------------------------------
# the legs
# ---------------------------------------------------------------------------


def run_leg(level):
    t_leg = time.perf_counter()
    band = FLOOR_BAND[level]
    print(f"=== GV5.1d leg [{level}]: near-band seed polish (floor band "
          f"{band:.2e}; T1 {TARGETS['T1']} / T2 {TARGETS['T2']}) ===",
          flush=True)

    st = gb.build_loose_state(level)   # the amended seed + wiring guard
    res = st["loose_res"]
    loose = gb._loose_committed(level)

    print(f"--- [{level}] tight pack ---", flush=True)
    t0 = time.perf_counter()
    pack = td.build_tight_pack(st, UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
    print(f"    pack in {time.perf_counter() - t0:.0f}s", flush=True)

    x_base = pack.x_base()
    fbl_base = gc.fbl_of(pack, x_base)
    floor = gb._gv51_seed_ibl_floor(level)
    print(f"    unperturbed seed F_BL={fbl_base:.3e} = "
          f"{fbl_base / floor:.2f}x the diagnosed floor {floor:.3e} "
          f"(the GV5.1b constructional read: inside the band {band:.2e})",
          flush=True)

    inflow_mask = np.asarray(st["inflow_mask"], dtype=bool)

    def fbl_eps(eps):
        return gc.fbl_of(pack, gc.perturb_delta(
            x_base, pack.n_free, pack.n_st, inflow_mask, eps))

    # -- T1 calibration ------------------------------------------------
    print(f"--- [{level}] T1 seed calibration (NEAR-BAND target window "
          f"{TARGETS['T1']}) ---", flush=True)
    t0 = time.perf_counter()
    eps1, f1, trace1, ok1, why1 = gc.calibrate_eps(fbl_eps, TARGETS["T1"])
    calib_wall = time.perf_counter() - t0
    _write_lines(
        f"seed_calibration_{level}.csv",
        "seed,eval,eps,f_bl_max,window_lo,window_hi",
        [f"T1,{k},{eps:.6e},{f:.6e},{TARGETS['T1'][0]:.1e},"
         f"{TARGETS['T1'][1]:.1e}" for k, eps, f in trace1])
    print(f"    T1 calibration in {calib_wall:.0f}s "
          f"({len(trace1)} evals): {why1}", flush=True)
    if not ok1:
        _record("calib", level,
                f"T1 calibration FAILED: {why1} -- the pre-registered "
                f"abort: leg stopped, user adjudication | unperturbed "
                f"seed F_BL={fbl_base:.3e} ({fbl_base / floor:.2f}x "
                f"floor)", "RECORDED")
        for b in ("(a)", "(b)", "(c)"):
            _record(b, level, "not reached (calibration abort)",
                    "RECORDED")
        return
    _record("calib", level,
            f"T1 eps={eps1:.3e} -> seed F_BL={f1:.3e} in "
            f"[{TARGETS['T1'][0]:.0e},{TARGETS['T1'][1]:.0e}] "
            f"({f1 / band:.2e}x the floor band; {len(trace1)} evals) | "
            f"unperturbed seed F_BL={fbl_base:.3e} "
            f"({fbl_base / floor:.2f}x floor)", "RECORDED")

    x1 = gc.perturb_delta(x_base, pack.n_free, pack.n_st, inflow_mask,
                          eps1)

    # -- band (a) ------------------------------------------------------
    e1, e2, e3, tol_e2, a_ok = gc.band_a_live_at(pack, x1)
    _record("(a)", level,
            f"live identities on the perturbed T1-seed J: e1={e1:.3e} "
            f"e2={e2:.3e} (<= {tol_e2:.1e} cond-aware, carried in) "
            f"e3={e3:.3e} (the suite gates the machine-precision "
            "algebra + the mu schedule)",
            "PASS" if a_ok else "FAIL")

    # -- T1 polish -----------------------------------------------------
    res1, wall1 = _polish(pack, x1, level)
    hist1 = res1["history"]
    e_seq1 = [float(h["block_max"][2]) for h in hist1]
    m_seq1 = [float(h["merit"]) for h in hist1]
    triples1 = gb._p_series(e_seq1, band=band)
    converged = bool(res1["converged"])
    n_band1 = band_entry_iter(e_seq1, band)

    # -- T2 escalation (fires only on < 3 triples, not converged) ------
    triples2, fired_t2 = [], False
    eps2 = f2 = None
    res2 = None
    wall2 = 0.0
    n_band2 = None
    if not converged and len(triples1) < 3:
        fired_t2 = True
        print(f"--- [{level}] T2 escalation ({len(triples1)} < 3 "
              f"triples on T1; target window {TARGETS['T2']}) ---",
              flush=True)
        eps2, f2, trace2, ok2, why2 = gc.calibrate_eps(
            fbl_eps, TARGETS["T2"])
        with open(os.path.join(
                RESULTS, f"seed_calibration_{level}.csv"), "a") as f:
            for k, eps, f in trace2:
                f.write(f"T2,{k},{eps:.6e},{f:.6e},"
                        f"{TARGETS['T2'][0]:.1e},{TARGETS['T2'][1]:.1e}\n")
        print(f"    T2 calibration: {why2}", flush=True)
        if ok2:
            x2 = gc.perturb_delta(x_base, pack.n_free, pack.n_st,
                                  inflow_mask, eps2)
            res2, wall2 = _polish(pack, x2, f"{level}_t2")
            e_seq2 = [float(h["block_max"][2]) for h in res2["history"]]
            triples2 = gb._p_series(e_seq2, band=band)
            n_band2 = band_entry_iter(e_seq2, band)
            converged = converged or bool(res2["converged"])
        else:
            _record("calib", f"{level}-T2",
                    f"T2 calibration FAILED: {why2} (T1 read stands)",
                    "RECORDED")
    pooled = list(triples1) + list(triples2)

    # -- band (b) ------------------------------------------------------
    verdict_b, key_b = gc.window_verdict(pooled, converged, level)
    p_str = ("; ".join(f"i{i}:{p:.2f}" for i, p in triples1)
             if triples1 else "n/a (no above-band triple on T1)")
    if triples2:
        p_str += (" | T2: " + "; ".join(f"i{i}:{p:.2f}"
                                        for i, p in triples2))
    pm_str = ("; ".join(f"i{i}:{p:.2f}"
                        for i, p in gb._p_series(m_seq1))
              or "n/a")
    slope, n_pairs = gc.regression_slope(e_seq1, band)
    cf = gc.contraction_factors(e_seq1, band)
    cf_str = "; ".join(f"i{i}:{d:.2f}dex" for i, d in cf) or "n/a"
    slope_str = (f"{slope:.2f} over {n_pairs} pairs"
                 if slope is not None else f"n/a ({n_pairs} pairs)")
    fin = hist1[-1]
    fin_fbl, fin_merit = float(fin["block_max"][2]), float(fin["merit"])
    g51b_fbl, g51b_merit = gc._gv51b_final(level)
    key_b_full = (
        f"{key_b} | T1 F_BL p (above-band): {p_str} | merit p (T1): "
        f"{pm_str} | contraction factors (T1): {cf_str} | regression "
        f"slope (T1): {slope_str} | BAND-ENTRY iter: T1="
        f"{'none' if n_band1 is None else n_band1}"
        + (f" T2={'none' if n_band2 is None else n_band2}"
           if fired_t2 else "")
        + f" (GV5.1c never entered the band from above) | T1 final F_BL "
        f"{fin_fbl:.3e} vs floor {floor:.3e} ({fin_fbl / floor:.2f}x), "
        f"final merit {fin_merit:.3e} vs GV5.1b committed "
        f"{g51b_merit:.3e}, termination={res1['termination']}")
    _record("(b)", level, key_b_full, verdict_b)

    # -- band (c) ------------------------------------------------------
    n_polish = res1["n_iter"] + (res2["n_iter"] if res2 else 0)
    n_total = res.n_outer + n_polish
    asp = 2 * loose["n_outer"]
    key_c = (f"N_polish={n_polish} (T1 {res1['n_iter']} term="
             f"{res1['termination']}"
             + (f" + T2 {res2['n_iter']} term={res2['termination']}"
                if res2 else "")
             + f"; band-entry iter="
             f"{'none' if n_band1 is None else n_band1}) | N_total="
             f"{n_total} (regen loose {res.n_outer} + {n_polish}) vs "
             f"committed loose {loose['n_outer']} | aspirational "
             f"N_polish <= {asp}: "
             f"{'met' if n_polish <= asp else 'NOT met'}"
             + (" (3+: recorded; user adjudication)" if n_polish >= 3
                else "")
             + " | deliberately off-point seeds: N_polish measures the "
             "near-band traversal, not production convergence")
    _record("(c)", level, key_c, "RECORDED")

    # -- compare.csv ---------------------------------------------------
    retries = sum(h.get("mu_retries", 0) for h in hist1)
    mus = [h["mu"] for h in hist1 if h.get("mu") is not None]
    g51c_fbl, g51c_merit = _history_final(GV51C_RESULTS, level)
    cmp_rows = [
        ("unperturbed seed |F_BL|_max (x the diagnosed floor)",
         f"{fbl_base:.6e} ({fbl_base / floor:.2f}x)", f"floor {floor:.3e}",
         "the GV5.1b constructional read (inside the band)"),
        ("T1 eps / seed |F_BL|_max (NEAR-BAND)", f"{eps1:.6e} / {f1:.6e}",
         f"window [{TARGETS['T1'][0]:.0e},{TARGETS['T1'][1]:.0e}]",
         f"{f1 / band:.2e}x the floor band; {len(trace1)} calibration "
         "evals; GV5.1c's T1 was [5e-2,5e-1]"),
        ("band_a_e1 diag(rn)@Jsc@diag(cn) vs J (relerr)",
         f"{e1:.3e}", f"<= {BAND_A_TOL[0]:.0e}",
         "on the perturbed T1-seed J"),
        ("band_a_e2 mu=0 damped vs undamped splu (relerr)",
         f"{e2:.3e}", f"<= {tol_e2:.1e}",
         "cond-aware tolerance carried in (the GV5.1b adjudication)"),
        ("band_a_e3 dx = C dy round-trip (relerr)",
         f"{e3:.3e}", f"<= {BAND_A_TOL[2]:.0e}", ""),
        ("T1 F_BL p above the floor band", p_str, "",
         "the band-(b) read"),
        ("T1 merit p (all contractions)", pm_str, "",
         "recorded alongside per the prereg"),
        ("T1 contraction factors log10(e_prev/e_i)", cf_str, "",
         "fallback context"),
        ("T1 regression slope log e_{i+1} vs log e_i", slope_str, "",
         "fallback context (>= 3 pairs)"),
        ("above-band triples T1 / pooled", f"{len(triples1)} / "
         f"{len(pooled)}", ">= 3 -> median-p band applies",
         f"T2 escalation fired: {fired_t2}"
         + (f" (eps={eps2:.3e}, seed F_BL={f2:.3e})" if res2 else "")),
        ("BAND-ENTRY iter T1 (T2 if fired)",
         f"{'none' if n_band1 is None else n_band1}"
         + (f" ({'none' if n_band2 is None else n_band2})"
            if fired_t2 else ""),
         "GV5.1c: never (both levels)",
         "the key new datum -- a near-band trajectory reaching the band "
         "has traversed the near-band region"),
        ("T1 final |F_BL|_max vs floor", f"{fin_fbl:.6e}",
         f"GV5.1b {g51b_fbl:.6e} / GV5.1c {g51c_fbl:.6e}",
         f"{fin_fbl / floor:.2f}x the diagnosed floor"),
        ("T1 final merit", f"{fin_merit:.6e}",
         f"GV5.1b {g51b_merit:.6e} / GV5.1c {g51c_merit:.6e}",
         "vs the committed finals (read; never recomputed)"),
        ("T1 termination / converged",
         f"{res1['termination']} / {res1['converged']}",
         "GV5.1b: floor_reached (medium) / cap (coarse); GV5.1c: cap "
         "(both)", ""),
        ("N_polish (T1[+T2]) / N_total vs committed loose outer",
         f"{n_polish} / {n_total}", loose["n_outer"],
         f"regen loose n_outer={res.n_outer}"),
        ("mu schedule: final mu / total reject-retries (T1)",
         f"{mus[-1]:.3e} / {retries}" if mus else "n/a", "", ""),
        ("loose committed: cl / ds_max / ibl_floor",
         "", f"{loose['cl']:.8f} / {loose['ds_max']:.6e} / "
         f"{loose['ibl_floor']:.3e}", "read; never recomputed"),
        ("polish wall time (s) T1[+T2]",
         f"{wall1:.0f}+{wall2:.0f}" if fired_t2 else f"{wall1:.0f}", "",
         f"threads={THREADS} (temporary 8-thread session constraint; "
         "runner default 16)"),
    ]
    _merge_compare(level, cmp_rows)
    _record("wall", level,
            f"leg total {time.perf_counter() - t_leg:.0f}s (loose regen "
            f"+ pack + calibration + band (a) + polish"
            f"{'+T2' if fired_t2 else ''}) | threads={THREADS} "
            "(temporary 8-thread session constraint; runner default 16)",
            "RECORDED")
    print(f"GV5.1d [{level}] done in "
          f"{time.perf_counter() - t_leg:.0f}s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=["coarse", "medium"],
                    choices=["coarse", "medium"],
                    help="mesh level(s) to run (default: both)")
    args = ap.parse_args()

    print(f"GV5.1d runner: threads={THREADS} (temporary 8-thread "
          "session constraint; runner default 16)", flush=True)
    for level in args.levels:
        run_leg(level)

    n_fail = sum(1 for *_, v in SUMMARY if v == "FAIL")
    n_pass = sum(1 for *_, v in SUMMARY if v == "PASS")
    n_rec = sum(1 for *_, v in SUMMARY if v == "RECORDED")
    print(f"\nGV5.1d: {n_pass} PASS / {n_fail} FAIL / {n_rec} RECORDED",
          flush=True)
    if n_fail:
        print("HONEST FAIL -- see summary.csv", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
