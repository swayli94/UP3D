"""R22 -- the usable upwind_c interval at medium, and the cl spread across it.

Binding text: docs/dev_phase_five/20260824-0200-r22-prereg.md (committed first).

(c) cannot be populated at medium because C=1.0 clamps and C=3.0 limit-cycles. Turning
"cannot evaluate" into "here is the usable interval and the spread across it" is strictly
more information and does not re-spec the gate.

Run:  PYTHONNOUSERSITE=1 python bench/studies/r22_c_interval/run.py
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
from pyfp3d.mesh.reader import read_mesh                                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                               # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve                       # noqa: E402
from pyfp3d.post.shock import shock_report                              # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients                 # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting                    # noqa: E402

LEVEL, SEED, ALPHA, M_INF = "medium", 0, 1.25, 0.80
CANON = ("xcoarse", "coarse", "medium")
NEW_CS = (1.1, 1.25, 2.0, 2.5)
CACHED = {1.0: f"{ROOT}/bench/studies/r14_medium_coverage/results/medium_c10_s0.npz",
          1.5: f"{ROOT}/bench/studies/r12_h_pricing/results/medium.npz"}
SHOCK_REF, SHOCK_TOL = 0.62, 0.03      # ★ the COMMITTED reference, not 0.61 +- 0.02
LEG_S, TOTAL_S = 25 * 60, 90 * 60
IMPL, SUMMARY = {}, []


class Timeout(Exception):
    pass


def _alarm(s, f):
    raise Timeout()


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}", flush=True)


def spread(vals):
    """★ R13's lesson: a spread must be quoted WITH its denominator."""
    IMPL["E-VAR"] = True
    v = np.asarray(vals, float); d = v.max() - v.min()
    return dict(rel_min=d / abs(v.min()), rel_max=d / abs(v.max()),
                rel_mean=d / abs(v.mean()), lo=float(v.min()), hi=float(v.max()))


def classify_leg(hist, ch, nl, nf, conv):
    """★ never report conv=False; and report clamp_history too (R15)."""
    IMPL["E-MODE"] = True
    if conv and not nl and not nf:
        return "converged", ""
    m, ev = classify_failure(hist, ch, np.asarray([], float), 0, "", nl, nf)[:2]
    return m, ev


def main():
    os.makedirs(RES, exist_ok=True)
    t0 = time.perf_counter()
    assert LEVEL in CANON, "G-SCOPE"
    for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        print(f"  G-THREADS  {v} = {os.environ.get(v, '<unset>')}")
    print(f"  G-THREADS  load average {os.getloadavg()}")
    print(f"  G-KILL     per leg {LEG_S}s / total {TOTAL_S}s -- ENFORCED IN CODE")
    print(f"  G-RECIPE   NACA_KW + seed={SEED} alpha={ALPHA}; ONLY upwind_c varies",
          flush=True)
    mc, wc = cut_wake(read_mesh(f"{ROOT}/cases/meshes/naca0012_2.5d/{LEVEL}.msh"))
    dz = float(np.ptp(mc.nodes[:, 2]))
    signal.signal(signal.SIGALRM, _alarm)

    def post(phi):
        f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"],
                                    np.asarray(phi, float), alpha_deg=ALPHA,
                                    m_inf=M_INF, s_ref=dz)
        cur = wall_cp_curve(mc, np.asarray(phi, float), z=0.5 * dz, m_inf=M_INF)
        xs = shock_report(cur, M_INF)["upper"].get("x_shock")
        return float(f["cl"]), (float(xs) if xs is not None and np.isfinite(xs) else np.nan)

    rows = []
    for C, p in sorted(CACHED.items()):
        d = np.load(p)
        cl, xs = post(d["phi"])
        hist = np.asarray(d["res_hist"], float); ch = np.asarray(d["clamp_hist"], float)
        nl, nf = int(d["nlim"]), int(d["nflr"])
        m, ev = classify_leg(hist, ch, nl, nf, bool(d["conv"]))
        rows.append(dict(C=C, source="cached", conv=int(bool(d["conv"])), n_limited=nl,
                         n_floored=nf, usable=int(bool(d["conv"]) and not nl and not nf),
                         mode=m, evidence=ev, n_steps=len(hist),
                         steps_capped=int((ch[:, 0] > 0).sum()),
                         steps_floored=int((ch[:, 1] > 0).sum()),
                         cl_p=round(cl, 6), x_shock=round(xs, 6), wall_s=None))
    for C in NEW_CS:
        if time.perf_counter() - t0 > TOTAL_S:
            print("  ★ G-KILL total budget spent"); break
        kw = dict(NACA_KW); kw["upwind_c"] = C; kw["n_picard_seed"] = SEED
        print(f"\n  === upwind_c = {C} ===", flush=True)
        t1 = time.perf_counter(); signal.alarm(LEG_S)
        try:
            r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw)
        except Timeout:
            signal.alarm(0)
            print(f"  ★ G-KILL: C={C} exceeded {LEG_S}s -- killed, recorded", flush=True)
            rows.append(dict(C=C, source="KILLED", usable=0,
                             wall_s=round(time.perf_counter() - t1, 1))); continue
        signal.alarm(0)
        w = time.perf_counter() - t1
        hist = np.asarray(r["residual_history"], float)
        ch = np.asarray(r["clamp_history"], float)
        nl, nf = int(r["n_limited"]), int(r["n_floored"])
        np.savez_compressed(os.path.join(RES, f"C{C}.npz"),
                            phi=np.asarray(r["phi"], float), conv=bool(r["converged"]),
                            nlim=nl, nflr=nf, res_hist=hist, clamp_hist=ch,
                            gamma=np.asarray(r.get("gamma", []), float))
        cl, xs = post(r["phi"])
        m, ev = classify_leg(hist, ch, nl, nf, bool(r["converged"]))
        rows.append(dict(C=C, source="new", conv=int(bool(r["converged"])),
                         n_limited=nl, n_floored=nf,
                         usable=int(bool(r["converged"]) and not nl and not nf),
                         mode=m, evidence=ev, n_steps=len(hist),
                         steps_capped=int((ch[:, 0] > 0).sum()),
                         steps_floored=int((ch[:, 1] > 0).sum()),
                         cl_p=round(cl, 6), x_shock=round(xs, 6), wall_s=round(w, 1)))
        r_ = rows[-1]
        print(f"  {w:7.1f}s usable={r_['usable']} clamps={nl}/{nf} steps={len(hist)} "
              f"mode={m}\n           cl_p {cl:.6f}  x_shock {xs:.6f}  "
              f"G-CLAMPHIST capped/floored on {r_['steps_capped']}/{r_['steps_floored']} steps",
              flush=True)

    rows.sort(key=lambda r: r["C"])
    with open(os.path.join(RES, "c_interval.csv"), "w", newline="") as f:
        ks = sorted({k for r in rows for k in r})
        w_ = csv.DictWriter(f, fieldnames=ks); w_.writeheader(); w_.writerows(rows)
    print(f"\n  {'C':>6}{'usable':>8}{'clamps':>12}{'mode':>14}{'cl_p':>11}{'x_shock':>10}")
    for r in rows:
        print(f"  {r['C']:6.2f}{r.get('usable', 0):8d}"
              f"{str(r.get('n_limited', '-')) + '/' + str(r.get('n_floored', '-')):>12}"
              f"{str(r.get('mode', '-')):>14}{r.get('cl_p', float('nan')):11.6f}"
              f"{r.get('x_shock', float('nan')):10.6f}")

    ok = [r for r in rows if r.get("usable")]
    bad = [r for r in rows if not r.get("usable")]
    IMPL["E-EDGE"] = True
    lo_gap = max([r["C"] for r in bad if r["C"] < min(x["C"] for x in ok)], default=None)
    hi_gap = min([r["C"] for r in bad if r["C"] > max(x["C"] for x in ok)], default=None)
    _record("E-EDGE", "the usable upwind_c set at medium/seed 0",
            "report which two tested points the edges fall between -- ★ an edge is a GAP "
            "between test points, NOT a threshold (kill criterion 6)",
            f"usable C = {[r['C'] for r in ok]};  lower edge in ({lo_gap}, "
            f"{min(x['C'] for x in ok)}];  upper edge in [{max(x['C'] for x in ok)}, "
            f"{hi_gap});  gate asks for the whole of [1, 3]",
            "E-EDGE: usable interval is NARROWER than the gate's [1,3]"
            if (lo_gap is not None or hi_gap is not None) else
            "E-EDGE: the whole of [1,3] is usable")
    if len(ok) >= 2:
        s = spread([r["cl_p"] for r in ok])
        _record("E-VAR", "cl_p spread across the USABLE C set",
                "RECORDED, not a gate verdict (spec is the user's; and M0.80 medium is in "
                "the FOLD ZONE, discipline #4). (c)'s target is 3%",
                f"cl_p {s['lo']:.6f}..{s['hi']:.6f};  spread /min {100*s['rel_min']:.2f}% "
                f"/max {100*s['rel_max']:.2f}% /mean {100*s['rel_mean']:.2f}%  "
                f"over {len(ok)} usable C values",
                f"RECORDED: {100*s['rel_min']:.2f}% (/min) vs (c)'s 3% -- "
                + ("EXCEEDS" if s["rel_min"] > 0.03 else "INSIDE"))
        xs_ = spread([r["x_shock"] for r in ok])
        _record("E-XS", "x_shock spread across the usable set, vs the COMMITTED reference",
                f"reference {SHOCK_REF} +- {SHOCK_TOL} (shock_reference.csv; NOT the "
                "0.61 +- 0.02 that circulates in nine documents)",
                f"x_shock {xs_['lo']:.6f}..{xs_['hi']:.6f} (spread {100*xs_['rel_min']:.2f}% "
                f"/min);  all inside band: "
                f"{all(abs(r['x_shock'] - SHOCK_REF) <= SHOCK_TOL for r in ok)}", "RECORDED")
    _record("E-MODE", "failure mode of each unusable C (never conv=False) + clamp_history",
            "RECORDED",
            "; ".join(f"C={r['C']}: {r.get('mode')} ({r.get('n_limited')}/"
                      f"{r.get('n_floored')}, capped/floored on "
                      f"{r.get('steps_capped')}/{r.get('steps_floored')} steps)"
                      for r in bad), "RECORDED")
    reg = ("E-EDGE", "E-VAR", "E-MODE", "E-XS")
    IMPL["E-XS"] = IMPL.get("E-VAR", False)
    print("\n  G-CHECKOFF:")
    for c in reg:
        print(f"    {c:8} {'implemented' if IMPL.get(c) else '★ NOT IMPLEMENTED'}")
    _record("G-CHECKOFF", "every registered criterion has code", "all four",
            ", ".join(f"{c}={'yes' if IMPL.get(c) else 'NO'}" for c in reg),
            "PASS" if all(IMPL.get(c) for c in reg) else "★ FAIL")
    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w_ = csv.writer(f); w_.writerow(["tag", "metric", "band", "measured", "verdict"])
        w_.writerows(SUMMARY)
    print(f"\n  {time.perf_counter() - t0:.1f} s total", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
