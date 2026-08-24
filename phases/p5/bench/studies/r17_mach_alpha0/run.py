"""R17 -- alpha = 0, rising M-infinity: a strong shock with ZERO wake jump.

Binding text: phases/p5/docs/dev_phase_five/20260823-1300-r17-prereg.md (committed first).

A symmetric section at alpha = 0 has Gamma identically zero, so M-infinity strengthens
the shock without introducing a wake jump. That is the single-variable isolation R16's
confound could not give.

Also repairs R16's own scope defect: its 0.50-0.70 "shock band" was set from alpha=1.25's
shock position, while at alpha=0 the pocket sits near x/c 0.39. The confound measure here
is max M UPSTREAM of the TE band, which finds the pocket wherever it is; both numbers are
reported so R16 stays comparable.

Run:  PYTHONNOUSERSITE=1 python bench/studies/r17_mach_alpha0/run.py
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
from failure_modes import classify_failure                       # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2                 # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry             # noqa: E402
from pyfp3d.mesh.reader import read_mesh                                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                               # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared, q2_at_mach   # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting                    # noqa: E402

LEVEL, C_ARM, SEED, ALPHA = "medium", 1.0, 0, 0.0
CANON = ("xcoarse", "coarse", "medium")
MESH = f"{ROOT}/cases/meshes/naca0012_2.5d/{LEVEL}.msh"
ENVELOPE_MAX = 0.87                       # G-ENVELOPE, the documented envelope
MACHS = (0.82, 0.84, 0.86)
CACHED_A0_M080 = os.path.join(ROOT, "bench/studies/r16_alpha_dose/results/alpha0.0.npz")
CACHED_CTRL = os.path.join(ROOT, "bench/studies/r14_medium_coverage/results/"
                                 "medium_c10_s0.npz")   # alpha 1.25, M 0.80
TE_LO, TE_HI = 0.80, 1.20
OLD_SH = (0.50, 0.70)                     # R16's band, reported for comparability
LEG_S, TOTAL_S = 25 * 60, 60 * 60
SUMMARY = []


class Timeout(Exception):
    pass


def _alarm(s, f):
    raise Timeout()


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}", flush=True)


def read_field(phi, m_inf, nodes, el, B, cent, cap_q2):
    """upM = max M upstream of the TE band -- finds the pocket wherever it sits.
    old_sh = R16's fixed 0.50-0.70 band, kept only for comparability."""
    g = np.empty((len(el), 3)); q2 = np.empty(len(el))
    element_velocity_q2(el, B, np.asarray(phi, float), g, q2)
    M = np.sqrt(np.maximum(mach_number_squared(q2, m_inf, 1.4), 0.0))
    up = cent[:, 0] < TE_LO
    te = (cent[:, 0] >= TE_LO) & (cent[:, 0] <= TE_HI)
    sh = (cent[:, 0] >= OLD_SH[0]) & (cent[:, 0] <= OLD_SH[1])
    i = int(np.argmax(q2[up]))
    xs_up = cent[up, 0]
    out = dict(upM=float(M[up].max()), upM_x=float(xs_up[i]),
               old_shM=float(M[sh].max()),
               te_maxM=float(M[te].max()), te_sup=int((M[te] > 1).sum()),
               maxM=float(M.max()), cap_ratio=float(q2.max() / cap_q2))
    #: G-SYM -- if a TE structure appears at alpha=0 it should be up/down symmetric
    if out["te_sup"]:
        y = cent[te & (M > 1), 1]
        out["te_sym"] = round(float(min((y > 0).mean(), (y < 0).mean())), 4)
    else:
        out["te_sym"] = None
    return out


def main():
    os.makedirs(RES, exist_ok=True)
    t0 = time.perf_counter()
    assert LEVEL in CANON, "G-SCOPE"
    for m in MACHS:
        assert m <= ENVELOPE_MAX, f"G-ENVELOPE: {m} > {ENVELOPE_MAX}"
    for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        print(f"  G-THREADS  {v} = {os.environ.get(v, '<unset>')}")
    print(f"  G-THREADS  load average {os.getloadavg()}")
    print(f"  G-ENVELOPE all M_inf <= {ENVELOPE_MAX}: {MACHS}")
    print(f"  G-KILL     per leg {LEG_S}s / total {TOTAL_S}s -- ENFORCED IN CODE")
    print(f"  G-RECIPE   NACA_KW + upwind_c={C_ARM} seed={SEED} alpha={ALPHA}; "
          "ONLY m_inf varies", flush=True)

    mc, wc = cut_wake(read_mesh(MESH))
    nodes, el = mc.nodes, mc.elements
    B, _ = precompute_element_geometry(nodes, el)
    cent = nodes[el].mean(axis=1)
    signal.signal(signal.SIGALRM, _alarm)

    # ---- G-UPM: recompute the control's upM from cache, do NOT reuse 2.5483 ----
    dc = np.load(CACHED_CTRL)
    ctrl = read_field(dc["phi"], 0.80, nodes, el, B, cent,
                      float(q2_at_mach(3.0, 0.80, 1.4)))
    upM_ref = ctrl["upM"]
    _record("G-UPM", "the control's upM, recomputed from cache (NOT the old band's 2.5483)",
            "must be recomputed -- R16's band was scoped to alpha=1.25's shock position",
            f"alpha=1.25 / M0.80: upM {upM_ref:.4f} @x/c {ctrl['upM_x']:.4f}  "
            f"(old 0.50-0.70 band gave {ctrl['old_shM']:.4f});  TE maxM "
            f"{ctrl['te_maxM']:.4f}, sup {ctrl['te_sup']}", "G-UPM: recomputed")

    rows = []
    # cached alpha=0 / M0.80 from R16
    if os.path.exists(CACHED_A0_M080):
        d0 = np.load(CACHED_A0_M080)
        f0 = read_field(d0["phi"], 0.80, nodes, el, B, cent,
                        float(q2_at_mach(3.0, 0.80, 1.4)))
        ch0 = np.asarray(d0["clamp_hist"], float)
        rows.append(dict(m_inf=0.80, source="cached R16", wall_s=66.8,
                         conv=bool(d0["conv"]), n_limited=int(d0["nlim"]),
                         n_floored=int(d0["nflr"]), n_steps=len(d0["res_hist"]),
                         steps_capped=int((ch0[:, 0] > 0).sum()),
                         steps_floored=int((ch0[:, 1] > 0).sum()),
                         **{k: (round(v, 4) if isinstance(v, float) else v)
                            for k, v in f0.items()}))

    for m in MACHS:
        if time.perf_counter() - t0 > TOTAL_S:
            print(f"  ★ G-KILL total budget spent; skipping M={m}"); break
        kw = dict(NACA_KW); kw["upwind_c"] = C_ARM; kw["n_picard_seed"] = SEED
        print(f"\n  === M_inf = {m}  (alpha = {ALPHA}) ===", flush=True)
        t1 = time.perf_counter(); signal.alarm(LEG_S)
        try:
            r = solve_newton_lifting(mc, wc, m_inf=m, alpha_deg=ALPHA, **kw)
        except Timeout:
            signal.alarm(0)
            print(f"  ★ G-KILL: M={m} exceeded {LEG_S}s -- killed, recorded", flush=True)
            rows.append(dict(m_inf=m, source="KILLED_TIMEOUT",
                             wall_s=round(time.perf_counter() - t1, 1)))
            continue
        signal.alarm(0)
        w = time.perf_counter() - t1
        phi = np.asarray(r["phi"], float)
        hist = np.asarray(r["residual_history"], float)
        ch = np.asarray(r["clamp_history"], float)
        nl, nf = int(r["n_limited"]), int(r["n_floored"])
        np.savez_compressed(os.path.join(RES, f"M{m}.npz"), phi=phi,
                            conv=bool(r["converged"]), nlim=nl, nflr=nf,
                            res_hist=hist, clamp_hist=ch,
                            gamma=np.asarray(r.get("gamma", []), float))
        fr = read_field(phi, m, nodes, el, B, cent, float(q2_at_mach(3.0, m, 1.4)))
        mode = "converged" if bool(r["converged"]) else classify_failure(
            hist, ch, np.asarray([], float), 0, "", nl, nf)[0]
        rows.append(dict(m_inf=m, source="new", wall_s=round(w, 1),
                         conv=bool(r["converged"]), n_limited=nl, n_floored=nf,
                         n_steps=len(hist), mode=mode,
                         steps_capped=int((ch[:, 0] > 0).sum()),
                         steps_floored=int((ch[:, 1] > 0).sum()),
                         **{k: (round(v, 4) if isinstance(v, float) else v)
                            for k, v in fr.items()}))
        print(f"  {w:7.1f}s conv={r['converged']} clamps={nl}/{nf} steps={len(hist)} "
              f"mode={mode}\n           upM {fr['upM']:.4f} @x/c {fr['upM_x']:.4f}  "
              f"TE maxM {fr['te_maxM']:.4f} (sup {fr['te_sup']})  "
              f"old-band {fr['old_shM']:.4f}  q2/cap {fr['cap_ratio']:.4f}\n"
              f"           G-CLAMPHIST: capped {int((ch[:, 0] > 0).sum())} steps, "
              f"floored {int((ch[:, 1] > 0).sum())} steps (final {nl}/{nf})", flush=True)

    with open(os.path.join(RES, "mach_sweep.csv"), "w", newline="") as f:
        ks = sorted({k for d in rows for k in d})
        w_ = csv.DictWriter(f, fieldnames=ks); w_.writeheader(); w_.writerows(rows)

    got = [r_ for r_ in rows if "upM" in r_]
    print(f"\n  {'M_inf':>7}{'conv':>7}{'clamps':>10}{'upM':>9}{'upM x/c':>10}"
          f"{'TE maxM':>10}{'TE sup':>8}{'wall_s':>9}")
    for r_ in got:
        print(f"  {r_['m_inf']:7.2f}{str(r_['conv']):>7}"
              f"{str(r_['n_limited'])+'/'+str(r_['n_floored']):>10}"
              f"{r_['upM']:9.4f}{r_['upM_x']:10.4f}{r_['te_maxM']:10.4f}"
              f"{r_['te_sup']:8d}{r_['wall_s']:9.1f}")

    # ---- the three branches ----------------------------------------------
    diverged = [r_ for r_ in got if r_["te_sup"] > 0 or r_["n_limited"] or r_["n_floored"]]
    matched = [r_ for r_ in got if r_["upM"] >= upM_ref]
    best = max((r_["upM"] for r_ in got), default=float("nan"))
    if diverged:
        v = ("★★★ X-EXCLUDE: an alpha=0 leg DIVERGES at the TE with Gamma == 0 => "
             "wake coupling is EXCLUDED (decisive, needs no strength match)")
    elif matched:
        v = ("★★ X-NECESSARY: upM reached the control's value with the TE still clean "
             "=> the wake jump is a NECESSARY condition")
    else:
        v = (f"X-INCONCLUSIVE: best upM {best:.4f} < control {upM_ref:.4f} inside the "
             "envelope, TE clean throughout -- NO direction may be read")
    _record("X-BRANCH", "does alpha=0 (Gamma == 0) reproduce the TE divergence",
            "TE divergence at any leg => EXCLUDE (decisive);  upM >= control with TE "
            "clean => wake jump NECESSARY;  else INCONCLUSIVE",
            "  ".join(f"M{r_['m_inf']}: upM {r_['upM']:.4f}, TE sup {r_['te_sup']}, "
                      f"clamps {r_['n_limited']}/{r_['n_floored']}" for r_ in got)
            + f";  control upM {upM_ref:.4f}", v)
    _record("X-TREND", "upM and TE maxM against M_inf", "RECORDED",
            "  ".join(f"M{r_['m_inf']}: upM {r_['upM']:.4f} TE {r_['te_maxM']:.4f}"
                      for r_ in got), "RECORDED")
    _record("G-CLAMPHIST", "final-step scalars vs steps that actually clamped",
            "RECORDED -- R15: 0/0 does not mean never clamped",
            "  ".join(f"M{r_['m_inf']}: final {r_['n_limited']}/{r_['n_floored']}, "
                      f"clamped on {r_['steps_capped']}/{r_['steps_floored']} steps"
                      for r_ in got), "RECORDED")

    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w_ = csv.writer(f); w_.writerow(["tag", "metric", "band", "measured", "verdict"])
        w_.writerows(SUMMARY)
    print(f"\n  {time.perf_counter() - t0:.1f} s total", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
