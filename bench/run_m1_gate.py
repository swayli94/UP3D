"""Product metric M1, evaluated honestly (phase two GS1.5).

M1 (docs/dev_phase_two/roadmap.md 2): NACA0012 M0.80 / alpha 1.25 --
  (a) upper-surface shock at 0.61 +- 0.02 (the project's own Euler-anchored
      reference band, cases/reference_data/naca0012_m080/),
  (b) cl agreeing within 3 % between two mesh levels,
  (c) cl varying by less than 3 % over upwind_c in [1, 3].

Solved at the TARGET Mach directly (no Mach ramp -- GS1.2b measured the ramp to
reproduce the same answers where it works, to be 2.8x slower, and to fail
outright at the fine level).

This script exists so the FAIL is an artifact rather than a sentence: exit code
0 only if all three criteria hold. Decision D3 keeps M1 in the metric table with
its measured failure while M1a locks the in-envelope capability.

    python bench/run_m1_gate.py        # ~6 min
"""
import csv
import os
import sys
import time
from pathlib import Path

#: ★ PIN THE THREADS before numpy/numba import (2026-08-10, measured): the same M6 solve
#: took 566.6 s uncapped on 24 cores against 113.7 s capped at 8, same machine load --
#: a 5.0x oversubscription penalty. This script reports per-leg wall times, so leaving the
#: caps to the caller's shell makes those numbers unreproducible. `setdefault` so an
#: explicit export still wins; the resolved values are printed below.
THREAD_VARS = ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")
for _v in THREAD_VARS:
    os.environ.setdefault(_v, "8")

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                      # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve              # noqa: E402
from pyfp3d.post.shock import shock_report                     # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients        # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting           # noqa: E402
sys.path.insert(0, str(HERE))
from run_le14_common_root import classify_failure               # noqa: E402

M_INF, ALPHA = 0.80, 1.25
SHOCK_REF, SHOCK_TOL = 0.61, 0.02
CS = (1.0, 1.5, 3.0)
LEVELS = ("coarse", "medium")
#: ★ THE SEED IS AN AXIS, added 2026-08-05, because it moved a published number.
#: `solve_newton_lifting`'s n_picard_seed default went 5 -> 0 on 2026-08-02 ("noseed
#: global"), and this script does not pass it, so it silently inherited the change.
#: Measured consequence at M0.80 / alpha 1.25: coarse converges either way, but
#: MEDIUM dies at the m_cap limiter with seed 0 (|R| 3.29e-02, M_max exactly 3.0,
#: 7265 limited / 758 floored) and converges cleanly with seed 5 (|R| 2.85e-13,
#: 0/0, x_shock 0.6006) -- reproducing GS1b.11's published 0.6073 / 0.6006 exactly.
#: So the two seeds are BOTH measured here rather than one being chosen: 0 is what
#: production ships, 5 is what the published M1 numbers were produced on, and the
#: gap between them is itself the finding. The regression was invisible because the
#: matching test (test_p4_transonic::test_g41_transonic_medium_gate) is gated AND
#: xfail(strict=False) -- a non-strict xfail cannot detect a regression, and that
#: one was made non-strict for an unrelated reason (thread dependence).
SEEDS = (0, 5)


def main():
    print("threads (per-leg wall times are reported, so this is part of the result): "
          + ", ".join(f"{v}={os.environ[v]}" for v in THREAD_VARS))
    try:
        print(f"machine load average: {os.getloadavg()[0]:.1f} over {os.cpu_count()} cpus")
    except OSError:                                   # pragma: no cover - non-Linux
        pass
    # ★ CALL-SITE DEBT FIXED 2026-08-05. This script used to take --entropy,
    # defaulting to FALSE, so the plain invocation measured the ISENTROPIC density
    # long after the entropy correction became the library default (GS1b.11,
    # 2026-07-31). That is the same class of debt GS3.1 found for `precond`: the
    # library default had moved and the call sites had not. So the default leg now
    # passes NO density kwarg at all -- whatever the library ships is what M1 is
    # measured against -- and the isentropic control has to be asked for by name.
    #
    # The old artifact is instructive about why this matters: gate_results/
    # m1_gate_entropy.csv is dated 2026-07-29 15:15, which predates the sigma
    # freeze binding (GS1b.5), the m_cap escape fix (GS1b.10/11) and the
    # 2026-08-05 sigma-transport fix, and its medium legs show it -- C = 1.0
    # clamped at M_max 3.0 with cl 0.055, C = 3.0 "converged" at cl 0.0708 against
    # coarse's 0.370. Those numbers were never usable and the file kept them.
    #
    # gate_results/m1_gate.csv keeps its meaning (the isentropic-era record, which
    # roadmap 2 deliberately preserves); the library-default leg writes its own
    # file rather than overwriting it.
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--density", choices=("default", "isentropic"),
                    default="default",
                    help="'default' passes no density kwarg (= production); "
                         "'isentropic' is the control leg")
    a = ap.parse_args()
    kw_density = {} if a.density == "default" else {"entropy_correction": False}
    print(f"density = {a.density}"
          f"{'' if kw_density else '  (library default, no kwarg passed)'}")
    rows = []
    for level in LEVELS:
        path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
        if not path.exists():
            print(f"skip {level}: mesh missing")
            continue
        mc, wc = cut_wake(read_mesh(path))
        dz = float(np.ptp(mc.nodes[:, 2]))
        for seed in SEEDS:
          for C in CS:
            t0 = time.perf_counter()
            try:
                r = solve_newton_lifting(
                    mc, wc, m_inf=M_INF, alpha_deg=ALPHA, upwind_c=C,
                    m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
                    precond="direct", direct_refactor_every=4,
                    n_newton_max=80, n_picard_seed=seed, **kw_density)
                err = ""
            except Exception as exc:                           # noqa: BLE001
                rows.append(dict(level=level, C=C, n_picard_seed=seed,
                                 converged=False, error=type(exc).__name__))
                print(f"  {level} seed={seed} C={C}: FAILED "
                      f"{type(exc).__name__}", flush=True)
                continue
            #: ★ never report a bare conv=False (CLAUDE.md): the five modes have
            #: different signatures and completely different fixes, and correlating
            #: one label against a knob can be separating a MIXTURE of diseases.
            #: A converged leg gets no mode -- the classifier is a FAILURE
            #: classifier and says something arbitrary about a resolved trajectory.
            hist = np.asarray(r.get("residual_history", []), dtype=float)
            if r["converged"]:
                mode, ev = "", f"accept_reason={r.get('accept_reason')!r}"
            else:
                mode, ev, _d10, _rv = classify_failure(
                    hist, np.asarray(r.get("clamp_history", []), dtype=float),
                    np.asarray(r.get("F_history", []), dtype=float),
                    int(r.get("n_gmres_stalled") or 0),
                    str(r.get("accept_reason")),
                    int(r.get("n_limited") or 0), int(r.get("n_floored") or 0))
            rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                            m_inf=M_INF), M_INF)
            f = wall_force_coefficients(mc.nodes, mc.elements,
                                        mc.boundary_faces["wall"], r["phi"],
                                        alpha_deg=ALPHA, s_ref=dz,
                                        m_inf=M_INF)
            rows.append(dict(level=level, C=C, n_picard_seed=seed,
                             converged=r["converged"],
                             clamped=r.get("clamped"),
                             cl_p=round(f["cl"], 6),
                             x_shock=rep["upper"].get("x_shock"),
                             m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                             res_final=r["residual_history"][-1],
                             failure_mode=mode, mode_evidence=ev,
                             n_limited=r.get("n_limited"),
                             n_floored=r.get("n_floored"),
                             accept_reason=r.get("accept_reason"),
                             sigma_min=r.get("sigma_min"),
                             n_shock_cells=r.get("n_shock_cells"),
                             n_newton=len(hist),
                             #: ★ pre-registered instrumentation (2026-08-10, rule R1):
                             #: a seed-0 leg that converges may have been rescued by the
                             #: clamped-triggered fallback, and without this flag that
                             #: gets silently attributed to the SEED. Record which one
                             #: did the work.
                             fallback_fired=bool((r.get("seed_fallback") or {})
                                                 .get("fired")),
                             fallback_seed=(r.get("seed_fallback") or {}).get("seed"),
                             wall_s=round(time.perf_counter() - t0, 1),
                             error=err))
            print(f"  {level:7s} seed={seed} C={C:<4} "
                  f"conv={str(r['converged']):5s} "
                  f"cl={rows[-1]['cl_p']} x_shock={rows[-1]['x_shock']} "
                  f"M_max={rows[-1]['m_max']} |R|={rows[-1]['res_final']:.2e}"
                  f"{'  MODE=' + mode if mode else ''}", flush=True)

    # ★ `bench/results/` is gitignored (it holds the big bitcheck npz dumps), so
    # gate evidence goes to `bench/gate_results/`, which is TRACKED. Found
    # 2026-07-29: the GS1.5 close-out round file claimed the M1 FAIL artifact was
    # committed and it never was -- discipline #3 says a number living only in a
    # .md is not evidence, and that was exactly the situation.
    out_dir = HERE / "gate_results"
    out_dir.mkdir(exist_ok=True)
    out_csv = out_dir / ("m1_gate.csv" if a.density == "isentropic"
                         else "m1_gate_default.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)

    overall = True
    for seed in SEEDS:
        overall &= _criteria([r for r in rows if r.get("n_picard_seed") == seed],
                             seed)
    print(f"\nM1: {'PASS' if overall else 'FAIL'}")
    return 0 if overall else 1


def _criteria(rows, seed):
    print(f"\n=== M1 criteria, n_picard_seed = {seed}"
          f"{'  (library default)' if seed == 0 else '  (the published numbers)'}"
          " ===")
    ok = True
    # (a) shock band, default C, both levels
    for level in LEVELS:
        m = [r for r in rows if r.get("level") == level and r.get("C") == 1.5]
        if not m or m[0].get("x_shock") is None:
            print(f"  (a) {level}: no result -> FAIL")
            ok = False
            continue
        x, conv = m[0]["x_shock"], m[0]["converged"]
        good = conv and abs(x - SHOCK_REF) <= SHOCK_TOL
        print(f"  (a) {level:7s} x_shock {x:.4f} vs {SHOCK_REF}+-{SHOCK_TOL}"
              f"  conv={conv}  -> {'PASS' if good else 'FAIL'}")
        ok &= bool(good)
    # (b) two-level cl agreement at the default C
    cl = {}
    for level in LEVELS:
        m = [r for r in rows if r.get("level") == level and r.get("C") == 1.5
             and r.get("cl_p") is not None]
        if m:
            cl[level] = (m[0]["cl_p"], m[0]["converged"])
    if len(cl) == 2:
        d = (cl["medium"][0] - cl["coarse"][0]) / abs(cl["coarse"][0])
        good = abs(d) < 0.03 and all(c[1] for c in cl.values())
        print(f"  (b) cl coarse {cl['coarse'][0]:.4f} -> medium "
              f"{cl['medium'][0]:.4f} = {100 * d:+.1f} % (< 3 %)"
              f"  converged={[c[1] for c in cl.values()]}"
              f"  -> {'PASS' if good else 'FAIL'}")
        ok &= bool(good)
    else:
        print("  (b) FAIL: fewer than two levels produced a result")
        ok = False
    # (c) dissipation sensitivity per level
    for level in LEVELS:
        sub = [r for r in rows if r.get("level") == level
               and r.get("cl_p") is not None and r.get("converged")]
        if len(sub) < 2:
            print(f"  (c) {level}: {len(sub)} converged leg(s) -> FAIL")
            ok = False
            continue
        v = [r["cl_p"] for r in sub]
        spread = (max(v) - min(v)) / min(v)
        good = spread < 0.03
        print(f"  (c) {level:7s} cl over C in [1,3]: {min(v):.4f}..{max(v):.4f}"
              f" = {100 * spread:.1f} % (< 3 %) -> {'PASS' if good else 'FAIL'}")
        ok &= bool(good)

    return ok


if __name__ == "__main__":
    sys.exit(main())
