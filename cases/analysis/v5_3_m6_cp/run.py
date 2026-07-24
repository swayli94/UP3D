"""GV5.3 M6 wing direction+magnitude check vs committed Cp (Track V5).

Binding text: docs/roadmap/track_v.md GV5.3 (re-anchored 2026-07-22,
user-directed: the committed experiment holds Cp only -- no experimental
CL); pre-registered: cases/analysis/v5_3_m6_cp/PRE_REGISTRATION.md
(committed ba636d9 BEFORE the first code change).

  ONERA M6 (coarse recorded + medium binding), conforming path, TEST 2308
  verbatim M 0.8395 / alpha 3.06, Re_MAC = 11.72e6, forced transition
  x_tr/c = 0.05 both sides. The loose GV3.1 recipe verbatim (omega 1.0,
  <= 10 outer, tol_ds 1e-3) drives the GV5.0 wing case (tip band
  z > 0.95*b_semi pinned + mdot-masked). The FP driver is the P14
  transonic recipe verbatim (M0.70 probe seed -> the NEWTON_M6_RECIPE
  ramp imported from tests/test_p8_newton.py, pressure Kutta; warm outer
  solves; NO tip_taper so the k=0 solve anchors the committed P14
  numbers) + the pre-registered rescue chain (strict -> stall-accept ->
  warm continuation m_start=0.80 strict -> stall-accept; cold: ramp
  strict -> ramp stall-accept -> raise).

  Gates: (a) Delta-CL DOWN beyond the A4 floor (2.5%, medium binding,
  cl_KJ read; cl_p consistency quoted); (b) per-station Cp RMS to the
  committed 7-station experiment decreases vs the same-mesh k=0 inviscid
  baseline (5 unmasked stations binding: >= 4/5 + pooled decrease =
  PASS; eta 0.96/0.99 tip-masked RECORDED-only). (c) convergence/guards
  RECORDED. Wiring guards W1 (k=0 cl vs P14 anchors 1%), W2 (experiment
  side mapping), W3 (section extraction) raise = recipe error, not
  verdicts. Exits 1 iff any gate reads FAIL.

Run:  python cases/analysis/v5_3_m6_cp/run.py [--levels coarse medium]
"""

import argparse
import os
import sys
import time

# resolve pyfp3d from THIS worktree (the site-packages editable install may
# point at a sibling worktree)
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..")))

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at, x_le
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.surface import (
    cl_kj_3d,
    planform_area,
    wall_force_coefficients,
)
from pyfp3d.solve.newton import solve_newton_lifting, solve_newton_transonic
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    build_wing_case,
    run_loose_coupling,
)

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..")))
from tests.test_p8_newton import NEWTON_M6_RECIPE  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
M6_DIR = os.path.join(REPO, "cases", "meshes", "onera_m6")
EXP_FILE = os.path.join(REPO, "cases", "reference_data",
                        "onera_m6_experiment", "experiment-Cp.dat")

M_INF, ALPHA = 0.8395, 3.06            # TEST 2308 dataset label verbatim
MAC = 0.64607
RE_MAC = 11.72e6
RE_CHORD = RE_MAC / MAC                # per meter (meshes in NASA meters)
X_TR = 0.05                            # forced transition, both sides
TIP_FRAC = 0.05                        # tip mask = production tip_taper radius

ETAS = (0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99)
N_UNMASKED = 5                         # eta 0.96/0.99 inside the tip mask
A4_CL_FLOOR = 0.025                    # A4 medium peak-rel u_e input floor
A4_CP_SCALE = 0.05                     # ~2*|Cp|*2.5% Cp-scale annotation

# committed inviscid anchors (P14 pressure-Kutta, M0.84,
# cases/demo/p14_pressure_kutta/results/m084_pressure.csv); the 1% guard
# absorbs the dM = 0.0005 dataset-vs-anchor label difference
P14_ANCHOR = {"coarse": (0.262778, 0.268813),
              "medium": (0.277628, 0.282263)}

# production M6 Newton recipe (P14 demo, tier 1; the bridge verbatim)
M6_NEWTON_KW = dict(farfield_spanwise_gamma=True, precond="direct",
                    direct_refactor_every=1000, n_newton_max=60)

SUMMARY = []  # (gate, metric, band, measured, verdict)


def _record(gate, metric, band, measured, ok=None):
    verdict = "RECORDED" if ok is None else ("PASS" if ok else "FAIL")
    SUMMARY.append((gate, metric, band, measured, verdict))
    print(f"  [{verdict:8s}] {gate} {metric}: measured={measured} "
          f"(band: {band})", flush=True)


def _write_csv(name, header, rows):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write(",".join(str(x) for x in row) + "\n")
    print(f"  wrote {path}", flush=True)


# ---------------------------------------------------------------------------
# committed experiment: 7 Tecplot zones, (NP, X/L, Y/b, Z/L, Cp) per point;
# X/L = local x/c, Y/b = eta station, side via the Z/L sign (positive =
# upper/suction side -- W2-validated against the k=0 mapping at execution)
# ---------------------------------------------------------------------------

def parse_experiment(path=EXP_FILE):
    zones = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("ZONE"):
                zones.append([])
            elif s and not s.startswith(("TITLE", "VARIABLES")):
                p = s.split()
                if len(p) == 5:
                    zones[-1].append((float(p[1]), float(p[2]),
                                      float(p[3]), float(p[4])))
    exp = {}
    for z in zones:
        a = np.asarray(z)
        eta = round(float(a[0, 1]), 2)
        exp[eta] = {"x": a[:, 0], "zl": a[:, 2], "cp": a[:, 3],
                    "upper": a[:, 2] > 0.0}
    return exp


# ---------------------------------------------------------------------------
# FP driver: P14 transonic recipe verbatim + the pre-registered rescue chain.
# Per FP call, ordered cheap -> deep, FIRST success wins:
#   warm-started (k >= 1): single strict -> single stall-accept -> warm
#       continuation (m_start=0.80) strict -> continuation stall-accept ->
#       the GV3.3 loud raise.
#   cold start (k = 0):    the P14 chain (M0.70 probe seed ->
#       NEWTON_M6_RECIPE ramp, pressure Kutta) strict -> ramp stall-accept
#       -> raise.
# Every attempt's (path, accept_reason, converged) is logged for (c).
# ---------------------------------------------------------------------------

def make_m6_driver(mc, wc):
    log = []

    def _ramp_kw(pressure_seed=None, rhs=None, stall=False):
        kw = dict(NEWTON_M6_RECIPE["newton_kw"],
                  kutta_estimator="pressure", n_picard_seed=0)
        if pressure_seed is not None:
            kw.update(phi_init=pressure_seed[0], gamma_init=pressure_seed[1])
        if rhs is not None:
            kw["external_rhs"] = rhs
        if stall:
            kw["accept_on_stall"] = True
        return kw

    def _ramp_args():
        return {k: v for k, v in NEWTON_M6_RECIPE.items() if k != "newton_kw"}

    def _cold(rhs, stall):
        # P14 verbatim: M0.70 PROBE-Kutta seed, then the pressure ramp.
        # NO early return on a non-converged seed (addendum 2026-07-25 #1:
        # the ramp's level 0 re-converges the seed itself; the bridge's
        # M0.5 short-circuit poisoned the medium k=0 at the first
        # execution -- W1 caught it).
        r0 = solve_newton_lifting(mc, wc, m_inf=0.70, alpha_deg=ALPHA,
                                  **M6_NEWTON_KW)
        log.append(("cold_seed_m070", str(r0.get("accept_reason")),
                    bool(r0["converged"])))
        r = solve_newton_transonic(
            mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
            newton_kw=_ramp_kw((r0["phi"], r0["gamma"]), rhs, stall),
            **_ramp_args())
        log.append(("cold_ramp_stall_accept" if stall else "cold_ramp_strict",
                    str(r.get("accept_reason")), bool(r["converged"])))
        return r

    def _single(rhs, seed, stall):
        kw = dict(M6_NEWTON_KW, kutta_estimator="pressure", n_picard_seed=0,
                  phi_init=seed.phi, gamma_init=seed.gamma)
        if rhs is not None:
            kw["external_rhs"] = rhs
        if stall:
            kw["accept_on_stall"] = True
        r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw)
        log.append(("single_stall_accept" if stall else "single_strict",
                    str(r.get("accept_reason")), bool(r["converged"])))
        return r

    def _cont(rhs, seed, stall):
        r = solve_newton_transonic(
            mc, wc, m_inf=M_INF, alpha_deg=ALPHA, m_start=0.80,
            newton_kw=_ramp_kw((seed.phi, seed.gamma), rhs, stall),
            **_ramp_args())
        log.append(("continuation_stall_accept" if stall
                    else "continuation_strict",
                    str(r.get("accept_reason")), bool(r["converged"])))
        return r

    def solve(rhs, seed):
        warm = seed is not None and seed.phi is not None
        attempts = (
            [lambda: _single(rhs, seed, False),
             lambda: _single(rhs, seed, True),
             lambda: _cont(rhs, seed, False),
             lambda: _cont(rhs, seed, True)] if warm else
            [lambda: _cold(rhs, False), lambda: _cold(rhs, True)])
        r = attempts[0]()
        for attempt in attempts[1:]:
            if r["converged"]:
                break
            r = attempt()
        return r["phi"], r["gamma"], r

    return solve, log


# ---------------------------------------------------------------------------
# one level
# ---------------------------------------------------------------------------

def run_level(level, exp):
    print(f"--- GV5.3 {level}: ONERA M6 M={M_INF} alpha={ALPHA} "
          f"Re_MAC={RE_MAC:.3e} ---", flush=True)
    mc, wc = cut_wake(read_mesh(os.path.join(M6_DIR, f"{level}.msh")))
    wall = mc.boundary_faces["wall"]
    cfg = CouplingConfig(re_chord=RE_CHORD, m_inf=M_INF, alpha_deg=ALPHA,
                         x_tr_upper=X_TR, x_tr_lower=X_TR)
    case = build_wing_case(mc.nodes, mc.elements, wall, cfg,
                           x_le=x_le, chord_at=chord_at,
                           tip_mask_frac=TIP_FRAC)
    print(f"    IBL surface: {case.sm.n_node} nodes / {case.sm.n_tri} tris; "
          f"tip-masked {int(case.outflow_pin_surf.sum())}", flush=True)

    s_ref = planform_area(mc.nodes, wall)
    o = np.argsort(wc.station_z)
    k0 = {}

    def probe(phi, gamma, k):
        f = wall_force_coefficients(mc.nodes, mc.elements, wall, phi,
                                    alpha_deg=ALPHA, s_ref=s_ref,
                                    m_inf=M_INF)
        cl_kj = cl_kj_3d(np.asarray(gamma)[o], wc.station_z[o], s_ref,
                         B_SEMI)
        # progress logging only (the GV5.0 blind-run lesson: a multi-hour
        # loop with no per-iteration output is unmonitorable; not a
        # numerics change)
        print(f"    [k={k}] cl_p={f['cl']:.4f} cl_kj={float(cl_kj):.4f} "
              f"(+{time.perf_counter() - t0:.0f}s)", flush=True)
        if k == 0 and "phi" not in k0:
            k0.update(phi=phi.copy(), gamma=np.asarray(gamma).copy(),
                      cl_p=f["cl"], cl_kj=float(cl_kj))
        return {"cl_p": f["cl"], "cl_kj": float(cl_kj),
                "cd_p": f["cd_pressure"]}

    driver, path_log = make_m6_driver(mc, wc)
    t0 = time.perf_counter()
    res = run_loose_coupling(driver, case, cfg, probe=probe)
    wall_s = time.perf_counter() - t0
    print(f"    converged={res.converged} n_outer={res.n_outer} "
          f"wall={wall_s:.0f}s", flush=True)

    # -- W1 wiring guard: k=0 inviscid cl vs the committed P14 anchors ------
    a_p, a_kj = P14_ANCHOR[level]
    for name, got, ref in (("cl_p", k0["cl_p"], a_p),
                           ("cl_kj", k0["cl_kj"], a_kj)):
        rel = abs(got / ref - 1.0)
        if rel > 0.01:
            raise RuntimeError(
                f"W1 wiring guard: {level} k=0 {name} {got:.6f} vs the "
                f"committed P14 anchor {ref:.6f} (rel {rel:.3%} > 1%) -- "
                "the FP driver recipe drifted; recipe error, not a verdict")
    print(f"    [W1 ok] k=0 cl_p {k0['cl_p']:.4f} (anchor {a_p:.4f}), "
          f"cl_kj {k0['cl_kj']:.4f} (anchor {a_kj:.4f})", flush=True)

    # -- history CSV ---------------------------------------------------------
    head = ("k,ds_max,ds_change_rel,ds_neg_floored,mdot_max,ibl_n_iter,"
            "ibl_converged,ibl_final_residual,cl_p,cl_kj,cd_p")
    _write_csv(f"history_{level}.csv", head,
               [tuple(h.get(k, "") for k in head.split(","))
                for h in res.history])

    # -- section Cp at the 7 committed stations, k=0 vs terminal -------------
    states = (("k0", k0["phi"]), ("final", res.phi))
    curves = {}
    for tag, phi in states:
        curves[tag] = {}
        rows = []
        for eta in ETAS:
            sec = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI,
                                   u_inf=1.0, m_inf=M_INF)
            curves[tag][eta] = sec
            for side in ("upper", "lower"):
                for xv, cv in zip(sec[f"x_{side}"], sec[f"cp_{side}"]):
                    rows.append((f"{eta:.2f}", side, f"{xv:.5f}",
                                 f"{cv:.6f}"))
        _write_csv(f"cp_stations_{tag}_{level}.csv", "eta,side,x_over_c,cp",
                   rows)

    # -- per-station RMS vs experiment (+ W2 side-mapping guard) -------------
    def station_rms(phi_curves, eta, flip=False):
        e = exp[eta]
        tot, n = 0.0, 0
        for want_upper in (True, False):
            side = "upper" if want_upper != flip else "lower"
            m = e["upper"] == want_upper
            if not np.any(m):
                continue
            cp_i = np.interp(e["x"][m], phi_curves[eta][f"x_{side}"],
                             phi_curves[eta][f"cp_{side}"])
            tot += float(np.sum((cp_i - e["cp"][m]) ** 2))
            n += int(m.sum())
        return (tot / max(n, 1)) ** 0.5, n

    for eta, e in exp.items():
        i_le = int(np.argmax(e["cp"]))
        if e["x"][i_le] >= 0.05:
            raise RuntimeError(
                f"W2 wiring guard: eta={eta} max-Cp point at x/c="
                f"{e['x'][i_le]:.4f} (not the LE stagnation point) -- "
                "the experiment side mapping is broken")
    rms_flip = sum(station_rms(curves["k0"], eta, flip=True)[0]
                   for eta in ETAS[:N_UNMASKED])
    rms_keep = sum(station_rms(curves["k0"], eta, flip=False)[0]
                   for eta in ETAS[:N_UNMASKED])
    if rms_keep >= rms_flip:
        raise RuntimeError(
            f"W2 wiring guard: k=0 pooled RMS with the chosen side mapping "
            f"({rms_keep:.4f}) >= the flipped mapping ({rms_flip:.4f})")
    print(f"    [W2 ok] side mapping: pooled k0 RMS {rms_keep:.4f} < "
          f"flipped {rms_flip:.4f}", flush=True)

    rms_rows, dec, tot_sq, n_pts = [], 0, {"k0": 0.0, "final": 0.0}, 0
    for j, eta in enumerate(ETAS):
        masked = j >= N_UNMASKED
        r0_, n0 = station_rms(curves["k0"], eta)
        r1_, _ = station_rms(curves["final"], eta)
        d = r1_ - r0_
        better = d < 0.0
        inp = abs(d) < A4_CP_SCALE
        if not masked:
            dec += int(better)
            tot_sq["k0"] += r0_ ** 2 * n0
            tot_sq["final"] += r1_ ** 2 * n0
            n_pts += n0
        rms_rows.append((f"{eta:.2f}", int(masked), n0, f"{r0_:.5f}",
                         f"{r1_:.5f}", f"{d:+.5f}", int(better),
                         int(inp)))
    _write_csv(f"cp_rms_{level}.csv",
               "eta,tip_masked,n_points,rms_inviscid,rms_viscous,delta_rms,"
               "better,input_limited", rms_rows)
    pooled0 = (tot_sq["k0"] / max(n_pts, 1)) ** 0.5
    pooled1 = (tot_sq["final"] / max(n_pts, 1)) ** 0.5

    # -- (a) Delta-CL direction gate -----------------------------------------
    cl1_kj = float(res.history[-1]["cl_kj"])
    cl1_p = float(res.history[-1]["cl_p"])
    dcl_kj, dcl_p = cl1_kj - k0["cl_kj"], cl1_p - k0["cl_p"]
    rel_kj = abs(dcl_kj) / max(abs(k0["cl_kj"]), 1e-30)
    rel_p = abs(dcl_p) / max(abs(k0["cl_p"]), 1e-30)
    if dcl_kj < 0 and rel_kj > A4_CL_FLOOR:
        ok_a = True
        flag_a = "DOWN beyond the A4 floor"
    elif rel_kj <= A4_CL_FLOOR:
        ok_a = None
        flag_a = "input-limited (|dcl| <= A4 2.5% floor)"
    else:
        ok_a = False
        flag_a = "UP beyond the A4 floor"
    gate_a = "GV5.3(a)" if level == "medium" else "GV5.3(a)-coarse"
    _record(gate_a, f"{level} Delta-cl_KJ ({k0['cl_kj']:.4f} -> "
            f"{cl1_kj:.4f})",
            "DOWN > 2.5% A4 floor = PASS; <= 2.5% = RECORDED input-limited"
            if level == "medium" else "recorded (binding = medium)",
            f"{dcl_kj:+.4f} ({100 * rel_kj:.2f}%) {flag_a}",
            ok_a if level == "medium" else None)
    _record(gate_a, f"{level} Delta-cl_p consistency ({k0['cl_p']:.4f} -> "
            f"{cl1_p:.4f})", "recorded consistency read",
            f"{dcl_p:+.4f} ({100 * rel_p:.2f}%)"
            + (" -- estimator disagreement > A4 floor"
               if abs(rel_p - rel_kj) > A4_CL_FLOOR else ""))

    # -- (b) per-station Cp RMS direction gate --------------------------------
    gate_b = "GV5.3(b)" if level == "medium" else "GV5.3(b)-coarse"
    ok_b = (dec >= 4) and (pooled1 < pooled0)
    _record(gate_b, f"{level} unmasked-station RMS decreases "
            f"({dec}/5) + pooled ({pooled0:.4f} -> {pooled1:.4f})",
            ">= 4/5 stations + pooled decrease = PASS"
            if level == "medium" else "recorded (binding = medium)",
            f"{dec}/5, pooled {pooled1 - pooled0:+.4f}",
            ok_b if level == "medium" else None)
    for j, eta in enumerate(ETAS):
        if j < N_UNMASKED:
            continue
        r0_, _ = station_rms(curves["k0"], eta)
        r1_, _ = station_rms(curves["final"], eta)
        _record(gate_b, f"{level} masked eta={eta:.2f} RMS "
                f"({r0_:.4f} -> {r1_:.4f})",
                "recorded-only (tip mask: expect ~no change)",
                f"{r1_ - r0_:+.5f}"
                + (" -- LARGE SHIFT at a masked station, investigate"
                   if abs(r1_ - r0_) > 0.02 else ""))

    # -- (c) convergence and guards RECORDED ----------------------------------
    hist = res.history
    _record("GV5.3(c)" if level == "medium" else "GV5.3(c)-coarse",
            f"{level} loose loop n_outer / converged", "recorded",
            f"{res.n_outer} / {res.converged}")
    mdot = [float(h.get("mdot_max", 0.0)) for h in hist]
    floor = max((h.get("ibl_final_residual", 0.0) for h in hist[1:]),
                default=0.0)
    _record("GV5.3(c)" if level == "medium" else "GV5.3(c)-coarse",
            f"{level} IBL residual floor / mdot_max first->last",
            "recorded", f"{floor:.3e} / {mdot[0]:.3e} -> {mdot[-1]:.3e}")
    live = ~case.outflow_pin_surf
    U = res.U
    a_mag = float(np.max(np.abs(U[live, 1]))) if np.any(live) else 0.0
    b_mag = float(np.max(np.abs(U[live, 2]))) if np.any(live) else 0.0
    _record("GV5.3(c)" if level == "medium" else "GV5.3(c)-coarse",
            f"{level} crossflow max|B|/max|A|", "recorded",
            f"{b_mag:.4e} / {a_mag:.4e} = "
            f"{b_mag / max(a_mag, 1e-30):.4f}")
    n_stall = sum(1 for _, a, c in path_log if a == "stall" and c)
    n_cont = sum(1 for p, _, _ in path_log if "continuation" in p
                 or "ramp" in p)
    _record("GV5.3(c)" if level == "medium" else "GV5.3(c)-coarse",
            f"{level} fp calls (continuation / stall-accepted) / wall",
            "recorded",
            f"{len(path_log)} ({n_cont} / {n_stall}) / {wall_s:.0f}s "
            "@8-thread session constraint, non-comparable")

    _panels(level, exp, curves)
    _panel_cl(level, res)
    return {"level": level, "res": res, "wall_s": wall_s}


def _failure_row(level, err):
    """The pre-registered recipe-limit clause (mdot runaway / FP chain
    exhausted -> run_loose_coupling raised): the level reads RECORDED."""
    _record("GV5.3" if level == "medium" else "GV5.3-coarse",
            f"{level} recipe-limit (loop raised)", "recorded, gates "
            "unreadable at this level" if level == "medium" else "recorded",
            str(err)[:160])
    return {"level": level, "res": None, "wall_s": 0.0}


# ---------------------------------------------------------------------------
# panels
# ---------------------------------------------------------------------------

def _panels(level, exp, curves):
    fig, axes = plt.subplots(2, 4, figsize=(16, 6.4), sharex=True)
    for j, eta in enumerate(ETAS):
        ax = axes[j // 4, j % 4]
        e = exp[eta]
        up = e["upper"]
        ax.plot(e["x"][up], e["cp"][up], "ko", ms=2.5, alpha=0.8,
                label="exp" if j == 0 else None)
        ax.plot(e["x"][~up], e["cp"][~up], "ks", ms=2.5, alpha=0.8)
        for tag, ls, lb in (("k0", "-", "inviscid (k=0)"),
                            ("final", "--", "viscous (terminal)")):
            c = curves[tag][eta]
            ax.plot(c["x_upper"], c["cp_upper"], "C0", ls=ls, lw=1.1,
                    label=lb if j == 0 else None)
            ax.plot(c["x_lower"], c["cp_lower"], "C1", ls=ls, lw=1.1)
        ax.set_title(f"eta = {eta:.2f}"
                     + (" (tip-masked)" if j >= N_UNMASKED else ""),
                     fontsize=9)
        ax.invert_yaxis()
        ax.grid(alpha=0.3)
        if j % 4 == 0:
            ax.set_ylabel("Cp")
        if j // 4 == 1:
            ax.set_xlabel("x / c")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"GV5.3 M6 {level} M{M_INF} a{ALPHA}: Cp vs the committed "
                 "TEST 2308 experiment (circles/squares = exp upper/lower; "
                 "solid = k=0 inviscid, dashed = viscous terminal)")
    fig.tight_layout()
    path = os.path.join(RESULTS, f"cp_overlay_{level}.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}", flush=True)


def _panel_cl(level, res):
    h = res.history
    ks = [r["k"] for r in h]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharex=True)
    axes[0].plot(ks, [r.get("cl_p", np.nan) for r in h], "o-")
    axes[1].plot(ks, [r.get("cl_kj", np.nan) for r in h], "o-")
    axes[0].set_ylabel("cl_p (pressure integral)")
    axes[1].set_ylabel("cl_KJ (Gamma integration)")
    for ax in axes:
        ax.set_xlabel("outer iteration k (0 = inviscid baseline)")
        ax.grid(alpha=0.3)
    fig.suptitle(f"GV5.3 cl history of the loose loop (M6 {level} "
                 f"M{M_INF} a{ALPHA})")
    fig.tight_layout()
    path = os.path.join(RESULTS, f"cl_history_{level}.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}", flush=True)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=None,
                    choices=["coarse", "medium"])
    args = ap.parse_args()
    levels = args.levels or [
        lv for lv in ("coarse", "medium")
        if os.path.exists(os.path.join(M6_DIR, f"{lv}.msh"))]
    print(f"GV5.3 levels: {levels}", flush=True)
    exp = parse_experiment()
    missing = [e for e in ETAS if e not in exp]
    if missing:
        raise RuntimeError(f"experiment stations missing: {missing}")

    runs = []
    for lv in levels:
        try:
            runs.append(run_level(lv, exp))
        except RuntimeError as e:
            runs.append(_failure_row(lv, e))
            print(f"    LEVEL RECORDED (recipe-limit clause): {e}",
                  flush=True)

    _write_csv("summary.csv", "gate,metric,band,measured,verdict", SUMMARY)
    n_fail = sum(1 for *_, v in SUMMARY if v == "FAIL")
    n_pass = sum(1 for *_, v in SUMMARY if v == "PASS")
    n_rec = sum(1 for *_, v in SUMMARY if v == "RECORDED")
    print(f"\nGV5.3: {n_pass} PASS / {n_fail} FAIL / {n_rec} RECORDED",
          flush=True)
    if n_fail:
        print("HONEST FAIL -- see summary.csv", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
