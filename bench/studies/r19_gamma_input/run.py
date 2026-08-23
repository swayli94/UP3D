"""R19 experiment -- separate upM from |Gamma| by making Gamma an INPUT.

Binding text: docs/dev_phase_five/20260823-1700-r19-prereg.md.

R18 found upM and |Gamma| BOTH separate TE-clean from TE-supersonic legs, are not
independent, and that no cached leg breaks their correlation. With gamma_target the
correlation is breakable by construction: prescribe Gamma and let upM follow M_inf.

★ With tip_taper = 0 the Kutta condition is NOT enforced -- these are fixed-Gamma
DIAGNOSTIC probes (A2's discriminator), not physical solutions.

Run:  PYTHONNOUSERSITE=1 python bench/studies/r19_gamma_input/run.py
"""
import csv
import itertools
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
from pyfp3d.kernels.gradient import element_velocity_q2                 # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry             # noqa: E402
from pyfp3d.mesh.reader import read_mesh                                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                               # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared, q2_at_mach   # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting                    # noqa: E402

LEVEL, C_ARM, SEED, ALPHA = "coarse", 1.0, 0, 1.25
CANON = ("xcoarse", "coarse", "medium")
GAMMAS = (0.0, 0.1, 0.2, 0.3, 0.4)
MACHS = (0.80, 0.84, 0.86)
ENVELOPE_MAX = 0.87
TE_LO, TE_HI = 0.80, 1.20
MATCH_TOL = 0.05
LEG_S, TOTAL_S = 25 * 60, 90 * 60
SUMMARY = []


class Timeout(Exception):
    pass


def _alarm(s, f):
    raise Timeout()


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}", flush=True)


def main():
    os.makedirs(RES, exist_ok=True)
    t0 = time.perf_counter()
    assert LEVEL in CANON, "G-SCOPE"
    for m in MACHS:
        assert m <= ENVELOPE_MAX, "G-ENVELOPE"
    print(f"  G-THREADS  load average {os.getloadavg()}")
    print(f"  G-SCOPE    level={LEVEL}; fine NOT touched;  G-ENVELOPE M<= {ENVELOPE_MAX}")
    print(f"  G-RECIPE   NACA_KW + upwind_c={C_ARM} seed={SEED} alpha={ALPHA}, "
          f"tip_taper=0 (pin live); axes = gamma_target x m_inf", flush=True)
    mc, wc = cut_wake(read_mesh(f"{ROOT}/cases/meshes/naca0012_2.5d/{LEVEL}.msh"))
    nodes, el = mc.nodes, mc.elements
    B, _ = precompute_element_geometry(nodes, el)
    cent = nodes[el].mean(axis=1)
    z = np.zeros(wc.n_stations)
    signal.signal(signal.SIGALRM, _alarm)

    rows = []
    for g, m in itertools.product(GAMMAS, MACHS):
        if time.perf_counter() - t0 > TOTAL_S:
            print("  ★ G-KILL total budget spent"); break
        kw = dict(NACA_KW); kw["upwind_c"] = C_ARM; kw["n_picard_seed"] = SEED
        t1 = time.perf_counter(); signal.alarm(LEG_S)
        try:
            r = solve_newton_lifting(mc, wc, m_inf=m, alpha_deg=ALPHA, tip_taper=z,
                                     gamma_target=np.full(wc.n_stations, g), **kw)
        except Timeout:
            signal.alarm(0)
            rows.append(dict(gamma_target=g, m_inf=m, status="KILLED")); continue
        signal.alarm(0)
        w = time.perf_counter() - t1
        phi = np.asarray(r["phi"], float)
        gg = np.empty((len(el), 3)); q2 = np.empty(len(el))
        element_velocity_q2(el, B, phi, gg, q2)
        M = np.sqrt(np.maximum(mach_number_squared(q2, m, 1.4), 0.0))
        up = cent[:, 0] < TE_LO
        te = (cent[:, 0] >= TE_LO) & (cent[:, 0] <= TE_HI)
        ch = np.asarray(r["clamp_history"], float)
        rows.append(dict(gamma_target=g, m_inf=m, status="OK", wall_s=round(w, 1),
                         conv=int(bool(r["converged"])),
                         n_limited=int(r["n_limited"]), n_floored=int(r["n_floored"]),
                         steps_capped=int((ch[:, 0] > 0).sum()),
                         gamma_out=round(float(np.asarray(r["gamma"])[0]), 9),
                         upM=round(float(M[up].max()), 4),
                         te_maxM=round(float(M[te].max()), 4),
                         te_sup=int((M[te] > 1).sum())))
        print(f"  G={g:.1f} M={m:.2f}  {w:6.1f}s conv={rows[-1]['conv']} "
              f"clamps={rows[-1]['n_limited']}/{rows[-1]['n_floored']}  "
              f"upM {rows[-1]['upM']:8.4f}  TE {rows[-1]['te_maxM']:8.4f} "
              f"(sup {rows[-1]['te_sup']:5d})", flush=True)

    ok = [r for r in rows if r.get("status") == "OK"]
    with open(os.path.join(RES, "grid.csv"), "w", newline="") as f:
        ks = sorted({k for d in rows for k in d})
        w_ = csv.DictWriter(f, fieldnames=ks); w_.writeheader(); w_.writerows(rows)

    # ---- P-MATCH: matched upM, different Gamma ---------------------------
    pairs = []
    for a, b in itertools.combinations(ok, 2):
        if a["gamma_target"] == b["gamma_target"]:
            continue
        lo = min(a["upM"], b["upM"])
        if lo > 0 and abs(a["upM"] - b["upM"]) / lo <= MATCH_TOL:
            pairs.append((a, b))
    dis = [(a, b) for a, b in pairs if (a["te_sup"] > 0) != (b["te_sup"] > 0)]
    _record("P-MATCH", "pairs matched in upM (<=5%) with DIFFERENT gamma_target",
            "TE status differs within a matched pair => Gamma drives it;  TE status "
            "agrees throughout => upM drives it",
            f"{len(pairs)} matched pairs, {len(dis)} of them disagree on TE status"
            + ("; " + "; ".join(f"G{a['gamma_target']}/M{a['m_inf']}(upM {a['upM']}, "
                                f"sup {a['te_sup']}) vs G{b['gamma_target']}/M{b['m_inf']}"
                                f"(upM {b['upM']}, sup {b['te_sup']})" for a, b in dis[:4])
               if dis else ""),
            "★ P-MATCH: Gamma drives it" if dis else
            "★★ P-MATCH: upM drives it (TE status constant at matched upM)"
            if pairs else "P-MATCH: no matched pairs -- UNDEFINED")

    # ---- P-MATCHG: same Gamma, different upM ----------------------------
    flips = []
    for g in GAMMAS:
        col = sorted([r for r in ok if r["gamma_target"] == g], key=lambda r: r["upM"])
        if len({r["te_sup"] > 0 for r in col}) > 1:
            flips.append((g, [(r["m_inf"], r["upM"], r["te_sup"]) for r in col]))
    _record("P-MATCHG", "at FIXED gamma_target, does raising upM flip TE status",
            "flips at fixed Gamma => upM drives it;  never flips => it does not",
            f"{len(flips)}/{len(GAMMAS)} gamma_target columns flip; "
            + "; ".join(f"G{g}: {v}" for g, v in flips[:3]),
            "★★ P-MATCHG: upM drives it" if flips else
            "P-MATCHG: upM alone never flips TE status")

    _record("P-GRID", "the full (gamma_target, m_inf) -> (upM, TE sup) map", "RECORDED",
            "; ".join(f"G{r['gamma_target']}/M{r['m_inf']}: upM {r['upM']}, "
                      f"sup {r['te_sup']}" for r in ok), "RECORDED")
    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w_ = csv.writer(f); w_.writerow(["tag", "metric", "band", "measured", "verdict"])
        w_.writerows(SUMMARY)
    print(f"\n  {time.perf_counter() - t0:.1f} s total", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
