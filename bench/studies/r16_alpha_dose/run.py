"""R16 -- is the trailing-edge divergence lift/wake coupled? An alpha dose-response.

Binding text: phases/p5/docs/dev_phase_five/20260823-1100-r16-prereg.md (committed first).

Method is B23's: alpha zero clean and growing superlinearly means lift/wake coupled. The
confound is declared in advance -- alpha moves the wake jump AND the shock strength -- so
the shock-band max M is reported and W-DOSE is recorded as not isolable if it moves more
than 20%. One branch survives the confound: if alpha = 0 already diverges at the trailing
edge, wake coupling is excluded regardless.

Run:  PYTHONNOUSERSITE=1 python bench/studies/r16_alpha_dose/run.py
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

LEVEL, C_ARM, SEED = "medium", 1.0, 0
CANON = ("xcoarse", "coarse", "medium")
MESH = f"{ROOT}/cases/meshes/naca0012_2.5d/{LEVEL}.msh"
M_INF, GAMMA = 0.80, 1.4
ALPHAS = (0.0, 0.5)                       # 1.25 is cached from R14
CACHED_125 = os.path.join(ROOT, "bench/studies/r14_medium_coverage/results/"
                                "medium_c10_s0.npz")
TE_LO, TE_HI = 0.80, 1.20
SH_LO, SH_HI = 0.50, 0.70
CONFOUND_TOL = 0.20
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


def field_read(phi, nodes, el, B, cent, cap_q2):
    g = np.empty((len(el), 3)); q2 = np.empty(len(el))
    element_velocity_q2(el, B, np.asarray(phi, float), g, q2)
    M = np.sqrt(np.maximum(mach_number_squared(q2, M_INF, GAMMA), 0.0))
    te = (cent[:, 0] >= TE_LO) & (cent[:, 0] <= TE_HI)
    sh = (cent[:, 0] >= SH_LO) & (cent[:, 0] <= SH_HI)
    i = int(np.argmax(q2))
    return dict(maxM=float(M.max()), maxM_x=float(cent[i, 0]),
                te_maxM=float(M[te].max()), te_sup=int((M[te] > 1).sum()),
                sh_maxM=float(M[sh].max()), cap_ratio=float(q2.max() / cap_q2))


def main():
    os.makedirs(RES, exist_ok=True)
    t0 = time.perf_counter()
    assert LEVEL in CANON, "G-SCOPE"
    for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        print(f"  G-THREADS  {v} = {os.environ.get(v, '<unset>')}")
    print(f"  G-THREADS  load average {os.getloadavg()}")
    print(f"  G-KILL     per leg {LEG_S}s / total {TOTAL_S}s -- ENFORCED IN CODE")
    print(f"  G-RECIPE   NACA_KW + upwind_c={C_ARM} seed={SEED}; ONLY alpha varies")
    print(f"  G-SCOPE    level={LEVEL}; fine is NOT touched", flush=True)

    mc, wc = cut_wake(read_mesh(MESH))
    nodes, el = mc.nodes, mc.elements
    B, _ = precompute_element_geometry(nodes, el)
    cent = nodes[el].mean(axis=1)
    cap_q2 = float(q2_at_mach(3.0, M_INF, GAMMA))
    signal.signal(signal.SIGALRM, _alarm)

    rows = []
    # ---- the cached alpha = 1.25 control --------------------------------
    d = np.load(CACHED_125)
    ch = np.asarray(d["clamp_hist"], float)
    fr = field_read(d["phi"], nodes, el, B, cent, cap_q2)
    cs = np.where(ch[:, 0] > 0)[0]; fs = np.where(ch[:, 1] > 0)[0]
    rows.append(dict(alpha=1.25, source="cached R14", wall_s=1285.0,
                     conv=bool(d["conv"]), n_limited=int(d["nlim"]),
                     n_floored=int(d["nflr"]), n_steps=len(d["res_hist"]),
                     res_last=float(np.asarray(d["res_hist"], float)[-1]),
                     steps_capped=int(cs.size), steps_floored=int(fs.size),
                     **{k: round(v, 4) for k, v in fr.items()}))

    for a in ALPHAS:
        if time.perf_counter() - t0 > TOTAL_S:
            print(f"  ★ G-KILL total budget spent; skipping alpha={a}"); break
        kw = dict(NACA_KW); kw["upwind_c"] = C_ARM; kw["n_picard_seed"] = SEED
        print(f"\n  === alpha = {a} ===", flush=True)
        t1 = time.perf_counter(); signal.alarm(LEG_S)
        try:
            r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=a, **kw)
        except Timeout:
            signal.alarm(0)
            print(f"  ★ G-KILL: alpha={a} exceeded {LEG_S}s -- killed, recorded",
                  flush=True)
            rows.append(dict(alpha=a, source="KILLED_TIMEOUT",
                             wall_s=round(time.perf_counter() - t1, 1)))
            continue
        signal.alarm(0)
        w = time.perf_counter() - t1
        phi = np.asarray(r["phi"], float)
        hist = np.asarray(r["residual_history"], float)
        ch = np.asarray(r["clamp_history"], float)
        nl, nf = int(r["n_limited"]), int(r["n_floored"])
        np.savez_compressed(os.path.join(RES, f"alpha{a}.npz"), phi=phi,
                            conv=bool(r["converged"]), nlim=nl, nflr=nf,
                            res_hist=hist, clamp_hist=ch,
                            gamma=np.asarray(r.get("gamma", []), float))
        fr = field_read(phi, nodes, el, B, cent, cap_q2)
        #: G-CLAMPHIST -- the scalars are FINAL-STEP; the history is the process (R15)
        cs = np.where(ch[:, 0] > 0)[0]; fs = np.where(ch[:, 1] > 0)[0]
        mode = "converged" if bool(r["converged"]) else classify_failure(
            hist, ch, np.asarray([], float), 0, "", nl, nf)[0]
        rows.append(dict(alpha=a, source="new", wall_s=round(w, 1),
                         conv=bool(r["converged"]), n_limited=nl, n_floored=nf,
                         n_steps=len(hist), res_last=hist[-1], mode=mode,
                         steps_capped=int(cs.size), steps_floored=int(fs.size),
                         **{k: round(v, 4) for k, v in fr.items()}))
        print(f"  {w:7.1f}s conv={r['converged']} clamps={nl}/{nf} steps={len(hist)} "
              f"mode={mode}\n           maxM {fr['maxM']:.4f} @x/c {fr['maxM_x']:.4f}  "
              f"TE maxM {fr['te_maxM']:.4f} (sup {fr['te_sup']})  "
              f"SHOCK maxM {fr['sh_maxM']:.4f}  q2/cap {fr['cap_ratio']:.4f}\n"
              f"           G-CLAMPHIST: capped on {cs.size} steps, floored on "
              f"{fs.size} steps  (final scalars {nl}/{nf})", flush=True)

    with open(os.path.join(RES, "alpha_dose.csv"), "w", newline="") as f:
        ks = sorted({k for d_ in rows for k in d_})
        w_ = csv.DictWriter(f, fieldnames=ks); w_.writeheader(); w_.writerows(rows)

    got = {r_["alpha"]: r_ for r_ in rows if "te_maxM" in r_}
    print(f"\n  {'alpha':>6}{'conv':>7}{'clamps':>12}{'TE maxM':>10}{'TE sup':>8}"
          f"{'SHOCK maxM':>12}{'q2/cap':>12}")
    for a in sorted(got):
        r_ = got[a]
        print(f"  {a:6.2f}{str(r_['conv']):>7}"
              f"{str(r_['n_limited']) + '/' + str(r_['n_floored']):>12}"
              f"{r_['te_maxM']:10.4f}{r_['te_sup']:8d}{r_['sh_maxM']:12.4f}"
              f"{r_['cap_ratio']:12.4f}")

    # ---- W-CONFOUND (declared in advance) --------------------------------
    sh = [got[a]["sh_maxM"] for a in sorted(got)]
    rel = (max(sh) - min(sh)) / min(sh)
    isolable = rel <= CONFOUND_TOL
    _record("W-CONFOUND", "shock-band max M across the alpha legs (the confound)",
            f"<= {CONFOUND_TOL:.0%} => the dose-response is readable;  > that => "
            "W-DOSE is NOT isolable",
            "  ".join(f"a={a}: {got[a]['sh_maxM']:.4f}" for a in sorted(got))
            + f";  spread {100*rel:.1f} % (denominator = min)",
            "W-CONFOUND: isolable" if isolable else
            "★ W-CONFOUND FIRED: shock strength moves too much -- W-DOSE not isolable")

    # ---- W-DOSE ----------------------------------------------------------
    a0 = got.get(0.0)
    clean0 = bool(a0) and a0["te_sup"] == 0 and a0["n_limited"] == 0 and a0["n_floored"] == 0
    grows = all(got[a]["te_sup"] <= got[b]["te_sup"]
                for a, b in zip(sorted(got)[:-1], sorted(got)[1:]))
    if a0 and not clean0:
        verdict = ("★★★ W-DOSE: alpha=0 ALREADY diverges at the TE => wake/lift coupling "
                   "is EXCLUDED. ★ This branch is immune to W-CONFOUND.")
    elif clean0 and grows and isolable:
        verdict = "W-DOSE: lift/wake coupled (B23 class)"
    elif clean0 and grows and not isolable:
        verdict = ("W-DOSE: consistent with lift/wake coupling BUT NOT ISOLABLE "
                   "(W-CONFOUND fired) -- no attribution")
    else:
        verdict = "W-DOSE: UNDEFINED -- non-monotone or incomplete (recorded, not attributed)"
    _record("W-DOSE", "does the TE divergence scale with alpha",
            "alpha=0 clean and growing => lift/wake coupled;  alpha=0 already diverging "
            "=> wake coupling EXCLUDED (immune to the confound);  else UNDEFINED",
            "  ".join(f"a={a}: TE maxM {got[a]['te_maxM']:.4f}, sup {got[a]['te_sup']}, "
                      f"clamps {got[a]['n_limited']}/{got[a]['n_floored']}"
                      for a in sorted(got)),
            verdict)
    _record("W-COST", "wall time per leg", "RECORDED",
            "  ".join(f"a={a}: {got[a]['wall_s']}s" for a in sorted(got)), "RECORDED")

    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w_ = csv.writer(f); w_.writerow(["tag", "metric", "band", "measured", "verdict"])
        w_.writerows(SUMMARY)
    print(f"\n  {time.perf_counter() - t0:.1f} s total", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
