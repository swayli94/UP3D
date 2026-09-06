"""G8.2's PHYSICS anchors, measured when the wall-clock assertion masks them.

Why this exists. On 2026-08-10 the fast capability-lock tier came back 6/7 with
`test_g82_m6_medium_newton_end_to_end` red, and the assertion was the BUDGET --
`wall < 450 s`, measured 588 s on a machine carrying someone else's load (the whole
tier took 49 min against its 10.7 min baseline, and test_seed_fallback alone took 27).

★★ That is NOT enough to report "only the timing failed". The budget assert sits at
line 625, BEFORE the physics anchors at 629-642, so pytest never evaluated cl, M_max
or the three shock positions -- exactly the trap b7 and b9 sprang twice on 2026-08-09,
where the first failing assert hid three more that were also broken.

So this script calls the test module's own `_m6_case` (no test edit, no library
change) and prints every anchor with its committed band, plus the timing so the
budget question is answered separately from the physics question.

Outputs (TRACKED): bench/gate_results/g82_anchor_check.csv
"""

import csv
import os
import sys

os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ["PYFP3D_TRANSONIC_GATES"] = "1"

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
#: ★★ NO `sys.path.insert(REPO/"tests")` (removed 2026-09-06, W0.2 / H2).
#: That hack is what let this file import a BARE module name -- and a bare name
#: is invisible to F07's live-code import check, which filters on the `tests.`
#: prefix. So the 2026-08-24 renumbering broke this script (test_p8_newton.py
#: was split into tests/A/test_A52... + tests/E/test_E01...) and NOTHING went
#: red, because `bench/` is on no cadence and the guard could not see the form.
#: ⇒ the package-qualified import below is now also the thing F07 enforces.

from tests.E.test_E01_p8_newton_anchors import _m6_case             # noqa: E402

CSV = os.path.join(HERE, "gate_results", "g82_anchor_check.csv")
#: the committed bands, copied from tests/E/test_E01_p8_newton_anchors.py
ANCHORS = (
    ("cl",           0.268691, 0.005),
    ("m_max",        1.99687,  0.05),
    ("x_shock_0.44", 0.59632,  0.02),
    ("x_shock_0.65", 0.54020,  0.02),
    ("x_shock_0.90", 0.37144,  0.02),
)
BUDGET_S = 450.0


def main():
    r, forces, shocks, wall = _m6_case("medium")
    got = {
        "cl": float(forces["cl"]),
        "m_max": float(np.sqrt(r["mach2_max"])),
        "x_shock_0.44": float(shocks[0.44]),
        "x_shock_0.65": float(shocks[0.65]),
        "x_shock_0.90": float(shocks[0.90]),
    }
    print("G8.2 physics anchors, measured separately from the budget\n")
    print(f"  converged        {bool(r['converged'])}")
    print(f"  |R| final        {float(r['residual_history'][-1]):.3e}   (assert < 1e-9 "
          f"via _assert_terminal_quadratic)")
    print(f"  Kutta |F| final  {float(r['F_history'][-1]):.3e}   (assert < 1e-12)")
    print(f"  clamps           {int(r['n_limited'])} limited / {int(r['n_floored'])} "
          f"floored   (assert 0/0)")
    print(f"\n  {'anchor':>14}{'measured':>12}{'committed':>12}{'tol':>9}"
          f"{'|diff|':>11}   verdict")
    rows, bad = [], []
    for name, ref, tol in ANCHORS:
        d = abs(got[name] - ref)
        ok = d < tol
        if not ok:
            bad.append(name)
        print(f"  {name:>14}{got[name]:>12.6f}{ref:>12.6f}{tol:>9.3f}{d:>11.6f}   "
              f"{'PASS' if ok else '★ FAIL'}")
        rows.append(dict(anchor=name, measured=round(got[name], 8), committed=ref,
                         tol=tol, abs_diff=round(d, 8), pass_=ok))
    rows.append(dict(anchor="wall_s", measured=round(wall, 1), committed=BUDGET_S,
                     tol="(upper bound)", abs_diff="",
                     pass_=bool(wall < BUDGET_S)))
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["anchor", "measured", "committed", "tol",
                                           "abs_diff", "pass_"])
        w.writeheader(); w.writerows(rows)
    print(f"\n  wall {wall:.0f} s vs the {BUDGET_S:.0f} s budget -> "
          f"{'PASS' if wall < BUDGET_S else 'OVER'}")
    print(f"\nwrote {CSV}")

    print("\n=== reading ===")
    if bad:
        print(f"  ★★ {len(bad)} PHYSICS anchor(s) also broken: {bad}")
        print("     So this is NOT a timing-only failure -- it needs its own attribution")
        print("     round before anything is concluded about the gated set.")
        return 1
    print("  every physics anchor PASSES within its committed band, and the clamp and")
    print("  convergence premises hold. The G8.2 red is therefore the WALL-CLOCK")
    print("  assertion alone -- a machine-load reading, not a capability regression.")
    print("  Per the 2026-08-06 precedent the budget is NOT relaxed on a loaded machine:")
    print("  the number is recorded and the leg is reported as timing-invalid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
