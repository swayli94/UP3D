"""R15 -- the C=3.0 limit cycle and the rho_floor leg, with control arms.

Binding text: phases/p5/docs/dev_phase_five/20260823-0900-r15-prereg.md (committed first).

Feasibility resolved BEFORE registering: the floored cells cannot be located from what
the library returns (rho_tilde is not returned; the workspace exposes no rho/nu/floor
attribute). So the floor half is answered by WHEN, not WHERE, and the gap is recorded.

Run:  PYTHONNOUSERSITE=1 python bench/studies/r15_modes/run.py
"""
import csv
import os
import signal
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "bench"))
RES = os.path.join(HERE, "results")

from run_gs40h_strength_and_2p5d import NACA_KW                         # noqa: E402
from run_le14_common_root import classify_failure                       # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2                 # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry             # noqa: E402
from pyfp3d.mesh.reader import read_mesh                                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                               # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared, q2_at_mach   # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting                    # noqa: E402

LEVEL = "medium"
CANON = ("xcoarse", "coarse", "medium")
MESH = f"{ROOT}/cases/meshes/naca0012_2.5d/{LEVEL}.msh"
M_INF, ALPHA, GAMMA = 0.80, 1.25, 1.4
#: cheapest first (committed step counts 19 / 80 / 70)
LEGS = (("C3.0_s5_control", 3.0, 5), ("C3.0_s0_target", 3.0, 0),
        ("C1.0_s5_target", 1.0, 5))
ANCHOR = {"C3.0_s0_target": (0, 0, 80), "C1.0_s5_target": (2587, 391, 70)}
LEG_S, TOTAL_S = 25 * 60, 75 * 60          # G-KILL, now IN THE CODE
#: G-BANDS: must cover the whole x range, wake included (R14's defect)
SUMMARY = []


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}", flush=True)


def bands_for(x):
    """Bands covering the FULL x range, with an explicit wake band. G-BANDS asserts
    the fractions sum to 1 -- R14's four bands stopped at 1.01 and lost 15.4 %."""
    lo, hi = float(x.min()) - 1e-9, float(x.max()) + 1e-9
    edges = [lo, 0.0, 0.1, 0.5, 0.7, 1.0, hi] if hi > 1.0 else [lo, 0.0, 0.1, 0.5, 0.7, hi]
    names = ["UPSTREAM", "LE", "MID", "SHOCK", "TE", "WAKE"][:len(edges) - 1]
    out = []
    for nm, a, b in zip(names, edges[:-1], edges[1:]):
        if b <= a:
            continue
        out.append({"band": nm, "x_lo": round(a, 4), "x_hi": round(b, 4),
                    "n": int(((x >= a) & (x < b)).sum())})
    tot = sum(d["n"] for d in out)
    # the topmost band must be closed on the right
    out[-1]["n"] += int(len(x) - tot)
    for d in out:
        d["frac"] = d["n"] / len(x)
    s = sum(d["frac"] for d in out)
    assert abs(s - 1.0) < 1e-12, f"G-BANDS: fractions sum to {s}, not 1"
    return out


def diag(r):
    """Frozen-selection / conditioning diagnostics, Optional-safe (G-ARITY)."""
    k = ("assignment_cycle", "n_freeze_refresh", "n_freeze_reverts",
         "n_assignment_stale", "n_gmres_stalled", "n_gmres_total", "n_refactor",
         "accept_reason", "nu_max", "n_shock_cells", "n_sigma_refresh",
         "capture_n_refresh", "froze", "seed_fallback", "n_newton")
    return {n: (r[n] if n in r else "<absent>") for n in k}


def main():
    os.makedirs(RES, exist_ok=True)
    t0 = time.perf_counter()
    assert LEVEL in CANON, "G-SCOPE"
    for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        print(f"  G-THREADS  {v} = {os.environ.get(v, '<unset>')}")
    print(f"  G-THREADS  load average {os.getloadavg()}")
    print(f"  G-KILL     per leg {LEG_S}s, total {TOTAL_S}s -- ENFORCED IN CODE")
    print(f"  G-SCOPE    level={LEVEL}; fine is NOT touched", flush=True)

    mc, wc = cut_wake(read_mesh(MESH))
    nodes, el = mc.nodes, mc.elements
    B, _ = precompute_element_geometry(nodes, el)
    cent = nodes[el].mean(axis=1)
    cap_q2 = float(q2_at_mach(3.0, M_INF, GAMMA))
    signal.signal(signal.SIGALRM, _alarm)

    legs, brows = {}, []
    for name, C, seed in LEGS:
        if time.perf_counter() - t0 > TOTAL_S:
            print(f"  ★ G-KILL total budget spent; skipping {name}"); break
        kw = dict(NACA_KW); kw["upwind_c"] = C; kw["n_picard_seed"] = seed
        print(f"\n  === {name}  (C={C}, seed={seed}) ===", flush=True)
        t1 = time.perf_counter()
        signal.alarm(LEG_S)
        try:
            r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw)
        except Timeout:
            signal.alarm(0)
            print(f"  ★ G-KILL: {name} exceeded {LEG_S}s -- killed, recorded as a "
                  "reading, not retried", flush=True)
            legs[name] = {"leg": name, "C": C, "seed": seed, "status": "KILLED_TIMEOUT",
                          "wall_s": round(time.perf_counter() - t1, 1)}
            continue
        signal.alarm(0)
        w = time.perf_counter() - t1
        phi = np.asarray(r["phi"], float)
        hist = np.asarray(r["residual_history"], float)
        ch = np.asarray(r["clamp_history"], float)
        nl, nf = int(r["n_limited"]), int(r["n_floored"])
        np.savez_compressed(os.path.join(RES, f"{name}.npz"), phi=phi,
                            conv=bool(r["converged"]), nlim=nl, nflr=nf,
                            res_hist=hist, clamp_hist=ch,
                            gamma=np.asarray(r.get("gamma", []), float))
        g = np.empty((len(el), 3)); q2 = np.empty(len(el))
        element_velocity_q2(el, B, phi, g, q2)
        M = np.sqrt(np.maximum(mach_number_squared(q2, M_INF, GAMMA), 0.0))
        i = int(np.argmax(q2))
        te = (cent[:, 0] >= 0.80) & (cent[:, 0] <= 1.20)
        mode, ev = ("converged", "")
        if not bool(r["converged"]):
            mode, ev = classify_failure(hist, ch, np.asarray([], float), 0, "", nl, nf)[:2]
        rec = {"leg": name, "C": C, "seed": seed, "status": "OK", "wall_s": round(w, 1),
               "conv": bool(r["converged"]), "n_limited": nl, "n_floored": nf,
               "n_steps": len(hist), "res_last": hist[-1], "mode": mode, "evidence": ev,
               "maxM": round(float(M.max()), 4), "maxM_x": round(float(cent[i, 0]), 4),
               "te_maxM": round(float(M[te].max()), 4),
               "te_sup": int((M[te] > 1).sum()),
               "cap_ratio": round(float(q2.max() / cap_q2), 4), **diag(r)}
        legs[name] = rec
        print(f"  {w:7.1f}s conv={rec['conv']} clamps={nl}/{nf} steps={len(hist)} "
              f"mode={mode}\n           maxM {rec['maxM']} @x/c {rec['maxM_x']}  "
              f"TE maxM {rec['te_maxM']} (sup {rec['te_sup']})  q2/cap {rec['cap_ratio']}",
              flush=True)
        if name in ANCHOR:
            a = ANCHOR[name]
            ok = (nl, nf, len(hist)) == a
            _record(f"G-REPRO/{name}", "clamps and steps reproduce the committed run",
                    f"{a}", f"{(nl, nf, len(hist))}",
                    "PASS" if ok else "★ FAIL -- session difference, check before reading")
        # capped-cell bands, only where there ARE capped cells
        if nl:
            lim = q2 >= cap_q2 * (1.0 - 1e-12)
            bb = bands_for(cent[lim, 0])
            for d in bb:
                d.update(leg=name, n_capped=int(lim.sum()))
            brows += bb
            print("           capped bands: " + "  ".join(
                f"{d['band']} {100*d['frac']:.1f}%" for d in bb), flush=True)
        # clamp_history timeline (F-WHEN)
        if ch.size:
            steps_cap = np.where(ch[:, 0] > 0)[0]
            steps_flr = np.where(ch[:, 1] > 0)[0]
            rec["cap_steps"] = f"{steps_cap.min()}-{steps_cap.max()}" if steps_cap.size else "none"
            rec["flr_steps"] = f"{steps_flr.min()}-{steps_flr.max()}" if steps_flr.size else "none"
            rec["n_steps_capped"] = int(steps_cap.size)
            rec["n_steps_floored"] = int(steps_flr.size)
            print(f"           clamp timeline: capped on {rec['n_steps_capped']} steps "
                  f"({rec['cap_steps']}), floored on {rec['n_steps_floored']} steps "
                  f"({rec['flr_steps']})", flush=True)

    with open(os.path.join(RES, "legs.csv"), "w", newline="") as f:
        ks = sorted({k for d in legs.values() for k in d})
        w_ = csv.DictWriter(f, fieldnames=ks); w_.writeheader(); w_.writerows(legs.values())
    if brows:
        with open(os.path.join(RES, "capped_bands.csv"), "w", newline="") as f:
            ks = sorted({k for d in brows for k in d})
            w_ = csv.DictWriter(f, fieldnames=ks); w_.writeheader(); w_.writerows(brows)

    # ---- criteria ---------------------------------------------------------
    tgt, ctl = legs.get("C3.0_s0_target"), legs.get("C3.0_s5_control")
    if tgt and ctl and tgt.get("status") == "OK" and ctl.get("status") == "OK":
        fields = ("assignment_cycle", "n_freeze_refresh", "n_freeze_reverts",
                  "n_assignment_stale", "n_gmres_stalled")
        cmp = "; ".join(f"{k}: target {tgt[k]} vs control {ctl[k]}" for k in fields)
        churn = (tgt["assignment_cycle"] is True or
                 (isinstance(tgt["n_freeze_reverts"], (int, float)) and
                  isinstance(ctl["n_freeze_reverts"], (int, float)) and
                  tgt["n_freeze_reverts"] > ctl["n_freeze_reverts"]) or
                 (isinstance(tgt["n_freeze_refresh"], (int, float)) and
                  isinstance(ctl["n_freeze_refresh"], (int, float)) and
                  tgt["n_freeze_refresh"] > ctl["n_freeze_refresh"]))
        stall = isinstance(tgt["n_gmres_stalled"], (int, float)) and tgt["n_gmres_stalled"] > 0
        _record("L-MECH", "frozen-selection churn vs ill-conditioning, against the "
                "same-C converged control",
                "assignment_cycle True or freeze counts above control => churn;  "
                "n_gmres_stalled > 0 with freeze counts at control => ill-conditioning",
                cmp,
                "L-MECH: frozen-selection churn" if churn and not stall else
                "L-MECH: ill-conditioning" if stall and not churn else
                "L-MECH: UNDEFINED -- neither signature isolates (recorded, not attributed)")
        h = np.asarray(np.load(os.path.join(RES, "C3.0_s0_target.npz"))["res_hist"], float)
        r4 = np.round(h[-12:], 12)
        _, ev, d10, revis = classify_failure(h, np.zeros((len(h), 2)),
                                             np.asarray([], float), 0, "", 0, 0)
        _record("L-PERIOD", "does the residual tail revisit values", "RECORDED -- no "
                "threshold; a fixed one is a calibration",
                f"tail {' '.join(f'{v:.3e}' for v in h[-8:])};  {ev}", "RECORDED")
        _record("L-WHERE", "where max M sits with ZERO clamps", "RECORDED",
                f"target maxM {tgt['maxM']} @x/c {tgt['maxM_x']} (TE maxM {tgt['te_maxM']}, "
                f"sup {tgt['te_sup']}) vs control maxM {ctl['maxM']} @x/c {ctl['maxM_x']} "
                f"(TE {ctl['te_maxM']})", "RECORDED")

    f5 = legs.get("C1.0_s5_target")
    if f5 and f5.get("status") == "OK":
        _record("F-WHEN", "when floored appears, and whether it co-occurs with capped",
                "floored only on steps that are also capped => two ends of one "
                "divergence;  floored persisting without capped => two diseases",
                f"capped on {f5.get('n_steps_capped')} steps ({f5.get('cap_steps')}), "
                f"floored on {f5.get('n_steps_floored')} steps ({f5.get('flr_steps')})",
                "RECORDED (read in the verdict against the two branches)")
        cached = os.path.join(ROOT, "bench/studies/r14_medium_coverage/results/"
                                    "medium_c10_s0.npz")
        if os.path.exists(cached):
            d0 = np.load(cached)
            _record("F-SEED", "same C, only the seed differs",
                    "RECORDED -- what changes between seed 0 and seed 5 at C=1.0",
                    f"seed 0 (cached): clamps {int(d0['nlim'])}/{int(d0['nflr'])}, "
                    f"{len(d0['res_hist'])} steps;  seed 5: clamps "
                    f"{f5['n_limited']}/{f5['n_floored']}, {f5['n_steps']} steps, "
                    f"maxM {f5['maxM']} @x/c {f5['maxM_x']}", "RECORDED")
    _record("F-GAP", "can the floored cells be located from what the library returns",
            "resolved BEFORE registering: no",
            "rho_tilde is not returned and workspace exposes no rho/nu/floor/mask "
            "attribute; only clamp_history counts per step. CLAUDE.md's clamping "
            "signature asks for 'where those cells are' -- for n_floored that cannot "
            "be executed without a library change",
            "F-GAP: a named library instrumentation gap (NOT measured here)")

    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w_ = csv.writer(f); w_.writerow(["tag", "metric", "band", "measured", "verdict"])
        w_.writerows(SUMMARY)
    print(f"\n  {time.perf_counter() - t0:.1f} s total", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
