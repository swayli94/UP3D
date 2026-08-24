"""Product metric M1, evaluated honestly (phase two GS1.5).

M1 (phases/p2/docs/dev_phase_two/roadmap.md 2): NACA0012 M0.80 / alpha 1.25 --
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
#: ★★★ M1 RE-SPEC, 2026-08-24, BY USER RULING -- the FIRST relaxation of a target
#: number in this project (the phase-two plan recorded "目标数字一个都没有放松过").
#:   (a) shock band          -> DELETED as a criterion (kept below as a RECORDED reading)
#:   (b) two-level cl        -> 3 % relaxed to 20 %, and split out as its own metric M1b
#:   (c) C-sweep cl          -> 3 % relaxed to 20 %, and split out as its own metric M1c
#: ★★ WHY the old 3 % went: its only documented rationale was decidability, and its VALUE
#: was borrowed from M2's "cl within 3 % of experiment" target -- which the user retired as
#: meaningless, on grounds the project's own numbers support: the A4 input band is 2.5 %, so
#: a 3 % target sits at the measurement floor, and an inviscid full-potential model compared
#: against viscous experiment carries a model-form error larger than 3 % by construction
#: (GV5.2: every RAE2822 shock 0.06-0.10 c downstream). With that gone, the 3 % had no basis.
#: ★ The two criteria themselves are UNCHANGED in form and remain meaningful: (b) asks
#: whether the answer has stopped moving under refinement, (c) asks how much of the lift is
#: set by a stabilisation constant that does not exist in the true solution.
M1B_TOL = 0.20
M1C_TOL = 0.20
#: ★★ (a) is no longer a criterion, but the reading is still printed -- and printed against
#: the COMMITTED reference `cases/reference_data/naca0012_m080/shock_reference.csv`
#: (0.62 +- 0.03), NOT the 0.61 +- 0.02 that this constant used to carry and that circulates
#: in nine documents while appearing in no reference file (measured 2026-08-23).
#: ★★★ Consequence of deleting it, recorded because it is a real loss: (a) was M1's ONLY
#: externally anchored criterion. (b) and (c) are code-against-itself, so both can pass on a
#: solution that is uniformly wrong -- upwind_c = 1.10 at medium is a measured instance
#: (converged, 0 clamps, |R| 2.3e-13, Gamma off 7.6x, x_shock 0.657 OUTSIDE the band, and
#: (c) counts it as a legal leg). See bench/usability.py.
SHOCK_REF, SHOCK_TOL = 0.62, 0.03
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
#: ★★ A THIRD seed, opt-in via --third-seed (phase 3, 2026-08-11). Two seeds cannot see solution
#: NON-UNIQUENESS: measured at M0.80/alpha1.25 the discrete problem has SEVERAL roots and the seed
#: picks one, with cl_p spreading 23-32 % across {0, 5, 12} at a BIT-IDENTICAL mesh -- 8-11x the
#: 3 % of criteria (b)/(c). Seeds 0 and 5 happen to agree on coarse (0.02-2.76 %), so the gate's
#: own sampling is blind to it. Evidence:
#: phases/p3/docs/dev_phase_three/20260811-2100-task3-nonuniqueness-verdict.md.
#:
#: It is OPT-IN, and that is a recorded deviation from the registration's letter (addendum #2): the
#: default invocation must stay byte-compatible with the committed m1_gate_default.csv, which the
#: capability boundary cites. A --third-seed run writes its own file.
#:
#: ★ Seed 12 is DIAGNOSTIC ONLY -- it never enters the PASS/FAIL. Letting it vote would re-spec
#: M1, and the previous round handed M1's specification to the user.
THIRD_SEED = 12


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
    ap.add_argument("--third-seed", action="store_true",
                    help=f"also run n_picard_seed={THIRD_SEED} as a DIAGNOSTIC (never votes on "
                         f"PASS/FAIL) and write to a separate CSV; two seeds cannot see solution "
                         f"non-uniqueness")
    ap.add_argument("--density", choices=("default", "isentropic"),
                    default="default",
                    help="'default' passes no density kwarg (= production); "
                         "'isentropic' is the control leg")
    a = ap.parse_args()
    kw_density = {} if a.density == "default" else {"entropy_correction": False}
    seeds = SEEDS + ((THIRD_SEED,) if a.third_seed else ())
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
        for seed in seeds:
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
    #: ★ guard G3: never overwrite a committed record. m1_gate.csv's numbers are cited BY VALUE
    #: in the capability boundary (isentropic +101.4 % / 36.1 %), so a --third-seed run -- which
    #: has a different row count -- goes to its own file.
    _stem = "m1_gate" if a.density == "isentropic" else "m1_gate_default"
    out_csv = out_dir / f"{_stem}{'_3seed' if a.third_seed else ''}.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)

    _cross_seed_diagnostic(rows, seeds)

    overall = True
    #: ★ SEEDS, not `seeds`: the third seed is diagnostic and must not be able to move the verdict
    for seed in SEEDS:
        v = _criteria([r for r in rows if r.get("n_picard_seed") == seed], seed)
        #: ★ two INDEPENDENT metrics since the 2026-08-24 ruling; `overall` is kept only
        #: as a convenience exit code -- the per-metric verdicts are the result.
        overall &= all(v.values())
    print(f"\nM1: {'PASS' if overall else 'FAIL'}"
          f"{'   (seed %d ran as a DIAGNOSTIC and did not vote)' % THIRD_SEED if a.third_seed else ''}")
    return 0 if overall else 1


def _cross_seed_diagnostic(rows, seeds):
    """RECORDED, never a criterion: how far apart do the seeds land on the SAME mesh?

    ★ This exists because the data was already here and nothing looked at it. `_criteria` evaluates
    each seed's (a)/(b)/(c) separately and never asks whether the seeds AGREE -- so on the committed
    coarse rows, seed 0 and seed 5 differ by 2.76 % in cl_p at C=1.0, which is 92 % of the entire
    3 % budget of criteria (b) and (c), and that number was never printed. Same family as the
    project's own "a non-strict xfail cannot detect a regression": a criterion that does not look
    at something cannot fail on it.

    A spread needs TWO converged seeds or it is UNDEFINED -- reporting an under-populated spread as
    a small one is the vacuous-PASS mistake this project has already made once.
    """
    print(f"\n=== cross-seed consistency at a FIXED mesh (RECORDED -- not a criterion) ===")
    print(f"  seeds run: {list(seeds)}   (M1's verdict uses {list(SEEDS)} only)")
    print(f"  {'level':8}{'C':>5}{'conv':>7}{'cl_p spread':>14}{'rel %':>9}"
          f"{'x_shock spread':>16}")
    worst = 0.0
    for level in LEVELS:
        for C in CS:
            m = [r for r in rows if r.get("level") == level and r.get("C") == C]
            good = [r for r in m if r.get("converged") and r.get("cl_p") is not None]
            if len(good) < 2:
                print(f"  {level:8}{C:>5}{len(good):>4}/{len(m):<2}{'UNDEFINED':>14}"
                      f"{'-':>9}{'UNDEFINED':>16}")
                continue
            cl = [float(r["cl_p"]) for r in good]
            xs = [float(r["x_shock"]) for r in good if r.get("x_shock") is not None]
            dcl = max(cl) - min(cl)
            rel = 100.0 * dcl / max(abs(sum(cl) / len(cl)), 1e-12)
            worst = max(worst, rel)
            dxs = (max(xs) - min(xs)) if len(xs) >= 2 else None
            print(f"  {level:8}{C:>5}{len(good):>4}/{len(m):<2}{dcl:>14.6f}{rel:>9.2f}"
                  f"{('-' if dxs is None else f'{dxs:.4f}'):>16}")
    print(f"  worst cross-seed cl_p spread: {worst:.2f} %   "
          f"(criteria (b) and (c) both work at 3 %)")
    if worst > 3.0:
        print("  ★ the seed alone exceeds the (b)/(c) budget -- at a FIXED mesh, so it is not a")
        print("    resolution question. RECORDED: this does NOT change M1's verdict above.")


def _criteria(rows, seed):
    print(f"\n=== M1 criteria, n_picard_seed = {seed}"
          f"{'  (library default)' if seed == 0 else '  (the published numbers)'}"
          " ===")
    verdicts = {}
    # (a) RETIRED as a criterion by the 2026-08-24 ruling -- RECORDED only, and against
    # the COMMITTED band. It does not enter any PASS/FAIL below.
    for level in LEVELS:
        m = [r for r in rows if r.get("level") == level and r.get("C") == 1.5]
        if not m or m[0].get("x_shock") is None:
            print(f"  (a) {level}: no result -> RECORDED (no reading)")
            continue
        x, conv = m[0]["x_shock"], m[0]["converged"]
        inband = conv and abs(x - SHOCK_REF) <= SHOCK_TOL
        print(f"  (a) {level:7s} x_shock {x:.4f} vs committed {SHOCK_REF}+-{SHOCK_TOL}"
              f"  conv={conv}  in-band={inband}  -> RECORDED (criterion RETIRED)")
    # (b) two-level cl agreement at the default C
    cl = {}
    for level in LEVELS:
        m = [r for r in rows if r.get("level") == level and r.get("C") == 1.5
             and r.get("cl_p") is not None]
        if m:
            cl[level] = (m[0]["cl_p"], m[0]["converged"])
    if len(cl) == 2:
        d = (cl["medium"][0] - cl["coarse"][0]) / abs(cl["coarse"][0])
        good = abs(d) < M1B_TOL and all(c[1] for c in cl.values())
        print(f"  M1b  cl coarse {cl['coarse'][0]:.4f} -> medium "
              f"{cl['medium'][0]:.4f} = {100 * d:+.2f} % (< {100 * M1B_TOL:.0f} %)"
              f"  converged={[c[1] for c in cl.values()]}"
              f"  -> {'PASS' if good else 'FAIL'}")
        verdicts["M1b"] = bool(good)
    else:
        print("  M1b  FAIL: fewer than two levels produced a result")
        verdicts["M1b"] = False
    # (c) dissipation sensitivity per level
    m1c = True
    for level in LEVELS:
        sub = [r for r in rows if r.get("level") == level
               and r.get("cl_p") is not None and r.get("converged")]
        if len(sub) < 2:
            #: ★ this is a COVERAGE failure, not a precision one -- relaxing the tolerance
            #: cannot cure it, and saying so keeps the two apart in the output
            print(f"  M1c  {level:7s} {len(sub)} converged leg(s) -> FAIL "
                  "(COVERAGE, not precision -- a tolerance cannot fix it)")
            m1c = False
            continue
        v = [r["cl_p"] for r in sub]
        spread = (max(v) - min(v)) / min(v)
        good = spread < M1C_TOL
        print(f"  M1c  {level:7s} cl over C in [1,3]: {min(v):.4f}..{max(v):.4f}"
              f" = {100 * spread:.2f} % (< {100 * M1C_TOL:.0f} %) over {len(sub)} legs"
              f" -> {'PASS' if good else 'FAIL'}")
        m1c &= bool(good)
    verdicts["M1c"] = m1c
    print(f"  ⇒ M1b {'PASS' if verdicts['M1b'] else 'FAIL'};  "
          f"M1c {'PASS' if verdicts['M1c'] else 'FAIL'}   "
          "(two INDEPENDENT metrics since the 2026-08-24 ruling -- they ask different "
          "questions: convergence vs how much of the lift a stabilisation constant sets)")
    return verdicts


if __name__ == "__main__":
    sys.exit(main())
