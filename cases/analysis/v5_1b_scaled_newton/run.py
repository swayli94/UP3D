"""V5.1b -- the scaled + damped augmented Newton, the pre-floor quadratic
window.

Binding text: cases/analysis/v5_1b_scaled_newton/PRE_REGISTRATION.md
(committed at 8b7793f BEFORE the first execution, per discipline). Design
inputs: the committed IBL-floor diagnosis
(cases/analysis/v5_ibl_floor/results/findings.md). Solver-internal only:
the assembled F/J are bit-identical to GV5.1 (the committed FD verdicts
stand) -- this runner exercises newton_tight(scaling="rowcol",
lm_damping=True, floor_stop=True), the pre-registered R/C equilibration +
Levenberg mu schedule + floor-reached stop.

Protocol = the GV5.1 amended protocol verbatim:
  seed = HEAD-regenerated loose-converged state (the committed GV3.1
  recipe; build_loose_converged_state IMPORTED from
  ../v5_tight_coupling/run.py so the seed construction cannot drift),
  wiring guard |dcl_k0| <= 1e-8 per leg (abort + record on failure);
  augmented Newton as a POLISH (max_iter=10, tol=1e-8, tol_abs=1e-10);
  loose numbers read -- never recomputed -- from
  ../v3_loose_coupling/results/, GV5.1 finals and seed IBL floors read
  from ../v5_tight_coupling/results/.

Bands (the pre-registered ones):
  (a) implementation exactness: the algebraic identities live on the
      seed J (diag(rn) @ Jsc @ diag(cn) == J; mu=0 damped step ==
      undamped splu step; dx = C dy round-trip), PASS/FAIL; the test
      suite gates the machine-precision algebra + the mu schedule
      (tests/test_v5_tight_scaled.py).
  (b) pre-floor quadratic window: p on the F_BL max-norm sequence,
      measured on iterates with F_BL above the pre-registered floor
      band (coarse 3.16e-5, medium 1.71e-5 = 10x the diagnosed
      standalone floors); PASS = >= 3 consecutive contractions with
      median p in [1.5, 2.5] (medium binding, coarse recorded);
      converged outright -> trivially PASS; no slope-2 window -> the
      pre-registered fallback read (full history + partial success:
      descent below the GV5.1 committed final merit, floor-reached vs
      lambda-collapse), NOT a gate crash.
  (c) counts (recorded; aspirational): N_polish to enter the floor
      band (or converge), N_total = N_loose + N_polish vs the
      committed loose outer counts (4 coarse / 5 medium).

The standalone k=1-seed retry (the S3 fixture pack) is recorded only,
no band (the PRE_REGISTRATION).

Run:  python cases/analysis/v5_1b_scaled_newton/run.py
      python cases/analysis/v5_1b_scaled_newton/run.py --levels coarse
      python cases/analysis/v5_1b_scaled_newton/run.py --no-k1
(partial runs merge summary.csv by (band, level); a full run
regenerates every artifact from scratch)
"""

import argparse
import csv
import importlib.util
import os
import re
import sys
import time

# resolve pyfp3d from THIS worktree (the site-packages editable install may
# point at a sibling worktree)
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..")))

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.viscous import tight_driver as td
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    build_airfoil_case,
    make_picard_lifting_driver,
    run_loose_coupling,
)
from tests.v5_state import (
    ALPHA,
    M_CAP,
    M_CRIT,
    M_INF,
    NACA_DIR,
    RE,
    RHO_FLOOR,
    UPWIND_C,
    build_k1_state,
    build_naca_case,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

V3_RESULTS = os.path.join(HERE, "..", "v3_loose_coupling", "results")
GV51_RESULTS = os.path.join(HERE, "..", "v5_tight_coupling", "results")
GV51_RUN = os.path.join(HERE, "..", "v5_tight_coupling", "run.py")

# the pre-registered floor bands (10x the diagnosed standalone floors)
FLOOR_BAND = {"coarse": 3.16e-5, "medium": 1.71e-5}
RECIPE_TOL_K0 = 1e-8       # the GV5.1 wiring guard |dcl_k0|

# band (a) live-identity thresholds: e1/e2 are solve-free / same-solve
# algebra (machine precision); e3 is the round-trip through a cond ~
# 1e10 system (cond-amplified)
BAND_A_TOL = (1.0e-12, 1.0e-10, 1.0e-6)

SUMMARY = []               # (band, level, key_result, verdict)


def _write_summary():
    """summary.csv: one row per (band, level). Merge-written after every
    record: fresh rows replace the same (band, level) pair, rows of legs
    not re-run are kept (the v5_tight_coupling _merge_tagged discipline)
    -- a late crash loses only the legs that did not complete."""
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
    order = {"(a)": 0, "(b)": 1, "(c)": 2, "k1": 3, "wall": 4}
    out = sorted(rows.values(),
                 key=lambda ln: (order.get(ln.split(",", 1)[0], 99),
                                 ln.split(",", 2)[1]))
    _write_lines("summary.csv", "band,level,key_result,verdict", out)


def _record(band, level, key_result, verdict):
    SUMMARY.append((band, level, key_result, verdict))
    print(f"  [{verdict}] {band} {level}: {key_result}", flush=True)
    _write_summary()


def _write_lines(name, header, lines):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as f:
        f.write(header + "\n")
        for ln in lines:
            f.write(ln + "\n")
    print(f"  wrote {path}", flush=True)


def _load_gv51_recipe():
    """build_loose_converged_state from the COMMITTED GV5.1 runner
    (../v5_tight_coupling/run.py) -- imported, not mirrored, so the seed
    construction cannot drift from the amended-protocol recipe."""
    spec = importlib.util.spec_from_file_location("gv51_run", GV51_RUN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.build_loose_converged_state


def _gv51_final(level):
    """The GV5.1 committed amended final row (f_bl_max, merit) -- read,
    never recomputed."""
    with open(os.path.join(
            GV51_RESULTS, f"gv5_1_newton_history_{level}.csv")) as f:
        last = list(csv.DictReader(f))[-1]
    return float(last["f_bl_max"]), float(last["merit"])


def _gv51_seed_ibl_floor(level):
    """The seed IBL floor from the COMMITTED GV5.1 summary.csv (the
    'seed IBL floor (loose-final...)' row) -- read, never recomputed."""
    with open(os.path.join(GV51_RESULTS, "summary.csv")) as f:
        for ln in f:
            if (f",{level} seed IBL floor" in ln
                    and "final_residual=" in ln):
                m = re.search(r"final_residual=([0-9.eE+-]+)", ln)
                return float(m.group(1))
    raise RuntimeError(f"no committed GV5.1 seed-IBL-floor row for {level}")


def _loose_committed(level):
    """The committed GV3.1 loose history (read, never recomputed):
    n_outer, final cl, final ds_max, final ibl_final_residual."""
    with open(os.path.join(V3_RESULTS, f"gv3_1_history_{level}.csv")) as f:
        hist = list(csv.DictReader(f))
    last = hist[-1]
    return {
        "n_outer": max(int(r["k"]) for r in hist),
        "cl": float(last["cl"]),
        "ds_max": float(last["ds_max"]),
        "ibl_floor": float(last["ibl_final_residual"]),
        "k0_cl": float(hist[0]["cl"]),
    }


# ---------------------------------------------------------------------------
# state construction (the GV5.1 amended seed, verbatim)
# ---------------------------------------------------------------------------


def build_loose_state(level):
    """The loose loop's converged state, regenerated with the committed
    GV3.1 recipe (../v5_tight_coupling/run.py phase 1, verbatim: picard
    driver, omega=1.0, n_outer_max=10, tol_ds=1e-3), wiring guard
    included (converged + |dcl_k0| <= 1e-8 vs the committed GV3.1
    history, read -- never recomputed). Same body as the committed
    v5_ibl_floor runner's build_loose_state."""
    print(f"--- [{level}] loose-loop regeneration (committed GV3.1 "
          "recipe, the GV5.1 amended seed) ---", flush=True)
    mc, wc = cut_wake(read_mesh(os.path.join(NACA_DIR, f"{level}.msh")))
    cfg = CouplingConfig(re_chord=RE, m_inf=M_INF, alpha_deg=ALPHA)
    case = build_airfoil_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg)
    dz = float(np.ptp(mc.nodes[:, 2]))
    s_ref = 1.0 * dz

    from pyfp3d.post.surface import wall_force_coefficients

    def probe(phi, gamma, k):
        f = wall_force_coefficients(
            mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
            alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
        return {"cl": f["cl"], "cd_p": f["cd_pressure"]}

    t0 = time.perf_counter()
    driver = make_picard_lifting_driver(mc, wc, M_INF, ALPHA)
    res = run_loose_coupling(driver, case, cfg, probe=probe)
    wall = time.perf_counter() - t0
    loose = _loose_committed(level)
    d_k0 = abs(float(res.history[0]["cl"]) - loose["k0_cl"])
    guard_ok = bool(res.converged) and d_k0 <= RECIPE_TOL_K0
    print(f"    loose regen in {wall:.0f}s: converged={res.converged} "
          f"n_outer={res.n_outer} cl={float(res.history[-1]['cl']):.8f} "
          f"| wiring guard: dcl_k0={d_k0:.3e} (<= {RECIPE_TOL_K0:.0e}) "
          f"-> {'OK' if guard_ok else 'RECIPE ERROR'}", flush=True)
    if not guard_ok:
        raise RuntimeError(
            f"[{level}] loose-regen wiring guard failed (converged="
            f"{res.converged}, dcl_k0={d_k0:.3e}) -- the pre-registered "
            "abort: leg stopped, user adjudication")
    build_lc = _load_gv51_recipe()
    st = build_lc(mc, wc, cfg, case, res)
    st["loose_res"] = res
    st["level"] = level
    return st


# ---------------------------------------------------------------------------
# band (a): the live identities on the seed J
# ---------------------------------------------------------------------------


def band_a_live(pack):
    """The band-(a) algebraic identities evaluated on the assembled seed
    operators (assembly bit-identical to GV5.1): e1 = max-norm rel error
    of diag(rn) @ Jsc @ diag(cn) vs J; e2 = the mu=0 damped step vs the
    undamped splu step; e3 = the dx = C dy round-trip. Returns
    (e1, e2, e3) (e2/e3 nan if the seed J is splu-singular -- recorded,
    the polish itself starts at mu > 0 and is unaffected)."""
    x0 = pack.x_base()
    F0 = td.augmented_residual(pack, x0)
    J0 = td.augmented_jacobian(pack, x0)
    rn, cn, Jsc = td.equilibrate_rc(J0)
    E = (sp.diags(rn) @ Jsc @ sp.diags(cn) - J0).tocsr()
    e1 = float(abs(E.data).max() / abs(J0.data).max())
    try:
        d_ref = -sla.splu(J0.tocsc()).solve(F0)
        scale = max(float(np.max(np.abs(d_ref))), 1e-300)
        d0 = td.scaled_damped_step(J0, F0, 0.0)
        e2 = float(np.max(np.abs(d0 - d_ref)) / scale)
        d_rc = td.scaled_damped_step(J0, F0, 0.0, scaling="rowcol")
        e3 = float(np.max(np.abs(d_rc - d_ref)) / scale)
    except RuntimeError as exc:
        print(f"    band (a) live: seed J splu-singular ({exc}) -- "
              "e2/e3 recorded nan", flush=True)
        e2 = e3 = float("nan")
    ok = (e1 <= BAND_A_TOL[0]
          and (np.isnan(e2) or e2 <= BAND_A_TOL[1])
          and (np.isnan(e3) or e3 <= BAND_A_TOL[2]))
    print(f"    band (a) live on the seed J: e1={e1:.3e} (<="
          f"{BAND_A_TOL[0]:.0e}) e2={e2:.3e} (<={BAND_A_TOL[1]:.0e}) "
          f"e3={e3:.3e} (<={BAND_A_TOL[2]:.0e}) -> "
          f"{'PASS' if ok else 'FAIL'}", flush=True)
    return e1, e2, e3, ok


# ---------------------------------------------------------------------------
# band (b): the pre-floor quadratic window read
# ---------------------------------------------------------------------------


def _p_series(e, band=None):
    """p_i = log(e_i/e_{i-1}) / log(e_{i-1}/e_{i-2}) on the contraction
    steps (e_i < e_{i-1} < e_{i-2}); when band is not None, only
    contractions whose two predecessor iterates sit ABOVE the band are
    kept (the pre-registered floor-band filter). Returns
    list[(iter, p)]."""
    out = []
    for i in range(2, len(e)):
        if not (0.0 < e[i] < e[i - 1] < e[i - 2]):
            continue
        if band is not None and min(e[i - 2], e[i - 1]) <= band:
            continue
        out.append((i, float(np.log(e[i] / e[i - 1])
                             / np.log(e[i - 1] / e[i - 2]))))
    return out


def _longest_run(ps):
    """The longest iteration-consecutive sub-run of [(iter, p)]."""
    best, cur = [], []
    for item in ps:
        if cur and item[0] != cur[-1][0] + 1:
            if len(cur) > len(best):
                best = cur
            cur = []
        cur.append(item)
    if len(cur) > len(best):
        best = cur
    return best


def band_b_read(hist, level, converged):
    """The pre-registered band-(b) read on the F_BL max-norm sequence.
    Returns (verdict, key_result, run) with verdict in
    PASS/FAIL/RECORDED (coarse always RECORDED; the no-window fallback
    RECORDED per the pre-registered 'not a gate crash')."""
    band = FLOOR_BAND[level]
    e = [float(h["block_max"][2]) for h in hist]
    m = [float(h["merit"]) for h in hist]
    ps = _p_series(e, band=band)
    ps_merit = _p_series(m)
    run = _longest_run(ps)
    p_str = ("; ".join(f"i{i}:{p:.2f}" for i, p in ps)
             if ps else "n/a (no above-band contraction)")
    pm_str = ("; ".join(f"i{i}:{p:.2f}" for i, p in ps_merit)
              if ps_merit else "n/a")
    print(f"    band (b) [{level}]: floor band={band:.2e} | F_BL p "
          f"(above-band): {p_str} | merit p: {pm_str}", flush=True)
    if converged:
        return "PASS", ("converged outright -- the band is trivially "
                        f"PASS | F_BL p: {p_str} | merit p: {pm_str}"), run
    if len(run) >= 3:
        med = float(np.median([p for _, p in run]))
        hit = 1.5 <= med <= 2.5
        verdict = ("RECORDED" if level == "coarse"
                   else ("PASS" if hit else "FAIL"))
        return verdict, (f">= 3 consecutive above-band contractions "
                         f"(iters {run[0][0]}..{run[-1][0]}); median "
                         f"p={med:.2f} vs [1.5..2.5] | F_BL p: {p_str} "
                         f"| merit p: {pm_str}"), run
    return "RECORDED", ("no slope-2 window above the floor band -- the "
                        f"pre-registered fallback read | F_BL p: {p_str} "
                        f"| merit p: {pm_str}"), run


# ---------------------------------------------------------------------------
# the legs
# ---------------------------------------------------------------------------


def run_leg(level):
    t_leg = time.perf_counter()
    print(f"=== GV5.1b leg [{level}]: scaled + damped augmented Newton "
          f"polish from the GV5.1 amended seed (floor band "
          f"{FLOOR_BAND[level]:.2e}) ===", flush=True)

    st = build_loose_state(level)
    res = st["loose_res"]
    loose = _loose_committed(level)

    print(f"--- [{level}] tight pack ---", flush=True)
    t0 = time.perf_counter()
    pack = td.build_tight_pack(st, UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
    print(f"    pack in {time.perf_counter() - t0:.0f}s", flush=True)

    e1, e2, e3, a_ok = band_a_live(pack)
    _record("(a)", level,
            f"live identities on the seed J: e1={e1:.3e} e2={e2:.3e} "
            f"e3={e3:.3e} (the suite gates the machine-precision "
            "algebra + the mu schedule)", "PASS" if a_ok else "FAIL")

    print(f"--- [{level}] polish: newton_tight(scaling=rowcol, "
          "lm_damping, floor_stop; max_iter=10) ---", flush=True)
    t0 = time.perf_counter()
    res_n = td.newton_tight(pack, max_iter=10, tol=1e-8, tol_abs=1e-10,
                            line_search=True, verbose=True,
                            scaling="rowcol", lm_damping=True,
                            floor_stop=True)
    wall_newton = time.perf_counter() - t0
    hist = res_n["history"]
    term = res_n["termination"]
    print(f"    polish in {wall_newton:.0f}s: converged="
          f"{res_n['converged']} termination={term} "
          f"n_iter={res_n['n_iter']}", flush=True)

    # -- newton_history_{level}.csv (the pre-registered columns) ------
    band = FLOOR_BAND[level]
    e_bl = [float(h["block_max"][2]) for h in hist]
    m_seq = [float(h["merit"]) for h in hist]
    p_bl = dict(_p_series(e_bl, band=band))
    p_m = dict(_p_series(m_seq))
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
        f"newton_history_{level}.csv",
        "iter,f_phi_max,f_gamma_max,f_bl_max,merit,mu,lam,p_fbl,"
        "p_merit,ds_change,termination,wall_s", lines)

    # -- band (b) ------------------------------------------------------
    verdict_b, key_b, run = band_b_read(hist, level, res_n["converged"])
    gv51_fbl, gv51_merit = _gv51_final(level)
    fin = hist[-1]
    fin_fbl = float(fin["block_max"][2])
    fin_merit = float(fin["merit"])
    floor = _gv51_seed_ibl_floor(level)
    below_gv51 = fin_merit < gv51_merit
    lam_last = fin["lam"] if fin["lam"] is not None else 0.0
    fallback_note = (
        f" | fallback read: final merit {fin_merit:.3e} "
        f"{'<' if below_gv51 else '>='} GV5.1 committed "
        f"{gv51_merit:.3e}; termination={term} (floor-reached rather "
        f"than lambda-collapse: lam_last={lam_last:.4f})"
        if "fallback" in key_b else "")
    _record("(b)", level, key_b + fallback_note, verdict_b)

    # -- band (c) ------------------------------------------------------
    n_to_band = next(
        (i for i, v in enumerate(e_bl) if v <= band), None)
    n_polish = res_n["n_iter"]
    n_total = res.n_outer + n_polish
    asp = 2 * loose["n_outer"]
    key_c = (f"N_polish={n_polish} (termination={term}; literal "
             f"band-entry iter={'none' if n_to_band is None else n_to_band}"
             f" -- the seed F_BL {e_bl[0]:.3e} already sits inside the "
             f"band) | N_total={n_total} (regen loose {res.n_outer} + "
             f"{n_polish}) vs committed loose {loose['n_outer']} | "
             f"aspirational N_polish <= {asp}: "
             f"{'met' if n_polish <= asp else 'NOT met'}"
             + (" (3+: recorded; user adjudication)"
                if n_polish >= 3 else ""))
    _record("(c)", level, key_c, "RECORDED")

    # -- compare.csv ---------------------------------------------------
    retries = sum(h.get("mu_retries", 0) for h in hist)
    mus = [h["mu"] for h in hist if h.get("mu") is not None]
    cmp_rows = [
        ("band_a_e1 diag(rn)@Jsc@diag(cn) vs J (relerr)",
         f"{e1:.3e}", f"<= {BAND_A_TOL[0]:.0e}", ""),
        ("band_a_e2 mu=0 damped vs undamped splu (relerr)",
         f"{e2:.3e}", f"<= {BAND_A_TOL[1]:.0e}", ""),
        ("band_a_e3 dx = C dy round-trip (relerr)",
         f"{e3:.3e}", f"<= {BAND_A_TOL[2]:.0e}",
         "cond-amplified through cond(J) ~ 1e10"),
        ("seed |F_BL|_max", f"{e_bl[0]:.6e}", "",
         "regen seed (HEAD trajectory)"),
        ("seed merit", f"{m_seq[0]:.6e}", "", ""),
        ("final |F_BL|_max", f"{fin_fbl:.6e}", f"{gv51_fbl:.6e}",
         "vs GV5.1 committed amended final"),
        ("final merit", f"{fin_merit:.6e}", f"{gv51_merit:.6e}",
         f"{'below' if below_gv51 else 'at/above'} GV5.1 committed"),
        ("final |F_BL|_max vs diagnosed floor", f"{fin_fbl:.6e}",
         f"{floor:.3e}",
         f"ratio {fin_fbl / floor:.3f}x the seed IBL floor; floor band "
         f"{band:.2e}"),
        ("termination", term, "GV5.1: cap (lambda-collapse crawl)",
         f"lam_last={lam_last:.4f}"),
        ("converged", res_n["converged"], "GV5.1: False", ""),
        ("F_BL p above the floor band",
         "; ".join(f"i{i}:{p:.2f}" for i, p in _p_series(e_bl, band=band))
         or "n/a", "", "the band-(b) read"),
        ("merit p (all contractions)",
         "; ".join(f"i{i}:{p:.2f}" for i, p in _p_series(m_seq))
         or "n/a", "", "recorded alongside per the prereg"),
        ("longest above-band contraction run", len(run), "",
         ">= 3 -> median-p band applies"),
        ("N_polish (polish iterations to termination)", n_polish, "",
         f"literal band-entry iter "
         f"{'none' if n_to_band is None else n_to_band}"),
        ("N_total = N_loose + N_polish vs committed loose outer",
         n_total, loose["n_outer"], f"regen loose n_outer={res.n_outer}"),
        ("ds_change_last (polish)", f"{res_n['ds_change_last']:.6e}",
         "", "tol_ds = 1e-3 cross-check"),
        ("mu schedule: final mu / total reject-retries",
         f"{mus[-1]:.3e} / {retries}" if mus else "n/a", "", ""),
        ("loose committed: cl / ds_max / ibl_floor",
         "", f"{loose['cl']:.8f} / {loose['ds_max']:.6e} / "
         f"{loose['ibl_floor']:.3e}", "read; never recomputed"),
        ("polish wall time (s)", f"{wall_newton:.0f}", "", ""),
    ]
    _merge_compare(level, cmp_rows)
    _record("wall", level,
            f"leg total {time.perf_counter() - t_leg:.0f}s (loose regen "
            "+ pack + band (a) + polish)", "RECORDED")
    print(f"GV5.1b [{level}] done in {time.perf_counter() - t_leg:.0f}s",
          flush=True)


def _merge_compare(level, rows):
    """compare.csv: one row per (level, metric); fresh rows replace the
    same pair, other levels' rows kept."""
    header = "level,metric,gv5_1b,gv5_1_or_loose_committed,note"
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
    order = {"coarse": 0, "medium": 1, "k1standalone": 2}
    out = sorted(kept.values(),
                 key=lambda ln: (order.get(ln.split(",", 1)[0], 9),
                                 ln.split(",", 2)[1]))
    _write_lines("compare.csv", header, out)


def run_k1standalone():
    """The standalone k=1-seed retry on the S3 fixture pack (the
    Stage-3 smoke state) -- recorded only, no band (the prereg)."""
    t0_leg = time.perf_counter()
    print("=== GV5.1b k=1 standalone: scaled + damped polish from the "
          "S3 k=1 seed (recorded only, no band) ===", flush=True)
    st = build_k1_state(build_naca_case())
    pack = td.build_tight_pack(st, UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
    res_n = td.newton_tight(pack, max_iter=10, tol=1e-8, tol_abs=1e-10,
                            line_search=True, verbose=True,
                            scaling="rowcol", lm_damping=True,
                            floor_stop=True)
    wall = time.perf_counter() - t0_leg
    hist = res_n["history"]
    term = res_n["termination"]
    m_seq = [float(h["merit"]) for h in hist]
    p_m = dict(_p_series(m_seq))
    lines = []
    for i, h in enumerate(hist):
        bm = h["block_max"]
        lines.append(",".join([
            str(h["iter"]),
            f"{bm[0]:.6e}", f"{bm[1]:.6e}", f"{bm[2]:.6e}",
            f"{h['merit']:.6e}",
            "" if h.get("mu") is None else f"{h['mu']:.6e}",
            "" if h["lam"] is None else f"{h['lam']:.6f}",
            "",  # no floor band on this leg -> no p_fbl read
            "" if i not in p_m else f"{p_m[i]:.4f}",
            f"{h['ds_change']:.6e}",
            term if i == len(hist) - 1 else "",
            f"{h['wall_s']:.2f}"]))
    _write_lines(
        "newton_history_coarse_k1standalone.csv",
        "iter,f_phi_max,f_gamma_max,f_bl_max,merit,mu,lam,p_fbl,"
        "p_merit,ds_change,termination,wall_s", lines)

    with open(os.path.join(
            GV51_RESULTS,
            "gv5_1_newton_history_coarse_k1seed.csv")) as f:
        k1_last = list(csv.DictReader(f))[-1]
    k1_fbl, k1_merit = float(k1_last["f_bl_max"]), float(k1_last["merit"])
    fin = hist[-1]
    fin_fbl, fin_merit = float(fin["block_max"][2]), float(fin["merit"])
    lam_last = fin["lam"] if fin["lam"] is not None else 0.0
    key = (f"termination={term} n_iter={res_n['n_iter']} "
           f"converged={res_n['converged']} | final |F_BL|_max "
           f"{fin_fbl:.3e} vs k1seed committed {k1_fbl:.3e} | final "
           f"merit {fin_merit:.3e} vs k1seed committed {k1_merit:.3e} "
           f"({'below' if fin_merit < k1_merit else 'at/above'}) | "
           f"lam_last={lam_last:.4f} | {wall:.0f}s")
    _record("k1", "coarse", key, "RECORDED")
    _merge_compare("k1standalone", [
        ("final |F_BL|_max", f"{fin_fbl:.6e}", f"{k1_fbl:.6e}",
         "vs GV5.1 k1seed committed iter-10"),
        ("final merit", f"{fin_merit:.6e}", f"{k1_merit:.6e}",
         f"{'below' if fin_merit < k1_merit else 'at/above'} k1seed"),
        ("termination", term, "GV5.1 k1seed: cap", ""),
        ("n_iter", res_n["n_iter"], 10, ""),
        ("wall time (s)", f"{wall:.0f}", "", ""),
    ])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=["coarse", "medium"],
                    choices=["coarse", "medium"],
                    help="mesh level(s) to run (default: both)")
    ap.add_argument("--no-k1", action="store_true",
                    help="skip the standalone k=1-seed leg")
    args = ap.parse_args()

    for level in args.levels:
        run_leg(level)
    if not args.no_k1:
        run_k1standalone()

    n_fail = sum(1 for *_, v in SUMMARY if v == "FAIL")
    n_pass = sum(1 for *_, v in SUMMARY if v == "PASS")
    n_rec = sum(1 for *_, v in SUMMARY if v == "RECORDED")
    print(f"\nGV5.1b: {n_pass} PASS / {n_fail} FAIL / {n_rec} RECORDED",
          flush=True)
    if n_fail:
        print("HONEST FAIL -- see summary.csv", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
