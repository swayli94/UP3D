"""The FAST capability-lock tier: the gated anchors that must be run every close-out.

Why this exists (2026-08-09). The 17 gated test files hold the project's absolute
capability anchors and take over four hours to run, so in practice they were not run on
any cadence -- and on 2026-08-06 the first full gated run since the round-tip switch of
2026-08-04 came back 7 failed. One of those (b9) had been red for SEVENTEEN DAYS because
B28 retired a premise in the demo and never carried it to the test. Debt in the gated set
is invisible by construction: the ungated suite stays green while capability locks rot.

So this is a subset chosen to be CHEAP enough to actually run, and it is deliberately
paired with an explicit statement of what it does NOT cover, because a fast tier that is
mistaken for full coverage would recreate the exact failure it is meant to prevent.

    PYFP3D_TRANSONIC_GATES=1 python bench/run_capability_locks.py

Outputs (TRACKED): bench/gate_results/capability_locks.csv
"""

import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
CSV = os.path.join(HERE, "gate_results", "capability_locks.csv")

#: ★★ PIN THE THREADS (added 2026-08-10, measured -- this tier shipped without it and
#: the omission produced a false red on its first real use). One of the locks here,
#: G8.2, asserts a WALL-CLOCK budget, and `test_p8_newton`'s own docstring records the
#: protocol that number was measured under: NUMBA_NUM_THREADS *and* OMP_NUM_THREADS
#: *and* OPENBLAS_NUM_THREADS all capped (CLAUDE.md discipline 1; uncapped measured
#: ~333 s against 252 s capped on the reference box). This script left the caps to
#: whatever shell invoked it, so a timing GATE was reading the caller's environment.
#: Measured cost of the omission, same test, same machine load (average 14.3):
#:     uncapped (24 cores, oversubscribed)   566.6 s   FAIL
#:     capped at 8 threads                   113.7 s   PASS   <- 5.0x
#: The physics anchors were bit-identical in both, so the red was purely environmental
#: (bench/gate_results/g82_anchor_check.csv). `setdefault` so an explicit export still
#: wins -- and the resolved values are PRINTED, because a timing gate whose environment
#: is not on the record is not reproducible.
THREAD_VARS = ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")
DEFAULT_THREADS = "8"
for _v in THREAD_VARS:
    os.environ.setdefault(_v, DEFAULT_THREADS)

#: (node id, what it locks). Chosen for anchors that are cheap; the heavy RAMPS are in
#: the excluded list below with their measured costs, so the trade-off is visible.
#: ★★ TWO LOCKS REMOVED 2026-08-10 by the phases/ reorganisation, and this is a
#: disposition, not an accident. `test_b9_wingbody_ls.py` and
#: `test_b22_ls_3d_anchors.py::test_m6_coarse_ramp_anchor` were level-set locks, and ruling
#: D5 abandoned that route -- both files are archived under phases/p1/tests/ and phase 3's
#: first task deletes them. Keeping them here would have meant pointing a live close-out
#: gate into an archive, which is worse than not having the gate: it would run, pass, and
#: assert capability for a route nobody maintains. What they locked is recorded in
#: phases/p2/docs/dev_phase_two/20260806-1200-b22-respec.md and the b9 re-spec round.
#: ★ THIRD lock removed 2026-08-10 by phase 3 task 1: `test_b7_onera_m6.py` (the M6
#: level-set 3-D machinery + transonic gate, whose three legs were the strict xfails)
#: is archived with the route. ⇒ this tier is now 4 groups, all CONFORMING:
#: M1a, the seed fallback, the wing-body conforming locks, and P8's G8.1/G8.2.
LOCKS = (
    ("tests/test_s1_m1a_envelope.py",
     "M1a: the in-envelope three-level convergence (re-spec'd 2026-08-05)"),
    ("tests/test_seed_fallback.py",
     "the cold-start seed fallback, incl. the real NACA M0.80 medium recovery"),
    ("tests/test_b9_wingbody_conforming.py",
     "wing-body conforming: junction loading, the B8 lift-loss detector"),
    ("tests/test_p8_newton.py",
     "G8.1/G8.2 conforming Newton anchors (G8.2 re-anchored 2026-08-06)"),
)

#: ★ WHAT THIS TIER DOES NOT COVER, with the measured reason. Read this before treating a
#: green run as "the gated set passes" -- that conflation is what this file exists to
#: prevent, and the numbers are from the 2026-08-06 full gated run and the re-spec rounds.
NOT_COVERED = (
    ("phases/p1/tests/test_b22_ls_3d_anchors.py::test_m6_medium_ramp_anchor  [ARCHIVED]",
     "~35 min. strict-xfail since 2026-08-09, and the xfail still RUNS the solve"),
    ("tests/test_p4_transonic.py", "~32 min, dominated by the G4.1 medium gate"),
    ("tests/test_p5_onera_m6.py", "45-75 min from scratch (M6 medium continuation+polish)"),
    ("phases/p1/tests/test_b14_schur_ls.py  [ARCHIVED]", "M6 medium Schur/lagged-LU arms"),
    ("phases/p1/tests/test_b15_ls_newton_freeze.py  [ARCHIVED]", "the M6 medium M0.84 freeze ramp (~515 s idle)"),
    ("tests/test_b18_wingbody_transonic.py", "wing-body transonic ramps"),
    ("[ARCHIVED, except test_p8_jacobian] test_b6_transonic, test_b16/b17, test_b11, test_b19",
     "passed on 2026-08-06; cost not individually measured"),
)


def main():
    if os.environ.get("PYFP3D_TRANSONIC_GATES") != "1":
        print("★ PYFP3D_TRANSONIC_GATES=1 is required, else every gated lock SKIPS "
              "and this reports a vacuous green.")
        return 2
    print("capability locks, FAST tier -- run this at every close-out")
    print("  threads (G8.2 asserts a wall-clock budget, so this is part of the result): "
          + ", ".join(f"{v}={os.environ[v]}" for v in THREAD_VARS))
    try:
        print(f"  machine load average: {os.getloadavg()[0]:.1f} over "
              f"{os.cpu_count()} cpus")
    except OSError:                                   # pragma: no cover - non-Linux
        pass
    print()
    rows, t_all = [], time.perf_counter()
    for node, what in LOCKS:
        t0 = time.perf_counter()
        p = subprocess.run([sys.executable, "-m", "pytest", node, "-q", "--tb=line"],
                           cwd=REPO, capture_output=True, text=True)
        dt = time.perf_counter() - t0
        tail = [l for l in p.stdout.splitlines() if " passed" in l or " failed" in l
                or " error" in l]
        summary = tail[-1].strip() if tail else "(no summary)"
        rows.append(dict(node=node, locks=what, returncode=p.returncode,
                         summary=summary, wall_s=round(dt, 1)))
        flag = "OK  " if p.returncode == 0 else "FAIL"
        print(f"  {flag} {dt:7.1f}s  {node}")
        print(f"        {summary}")
        if p.returncode != 0:
            for l in p.stdout.splitlines():
                if l.startswith("FAILED") or l.startswith("ERROR"):
                    print(f"        {l}")
    total = time.perf_counter() - t_all
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["node", "locks", "returncode", "summary",
                                           "wall_s"], extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    bad = [r for r in rows if r["returncode"] != 0]
    print(f"\n  total {total:.0f} s ({total/60:.1f} min);  "
          f"{len(rows) - len(bad)}/{len(rows)} groups green")
    print(f"  wrote {CSV}")
    print("\n  [ARCHIVED] = under phases/p1/ per ruling D5; deleted in phase 3, "
          "not runnable from tests/")
    print("\n★ NOT covered by this tier -- a green run here is NOT 'the gated set passes':")
    for node, why in NOT_COVERED:
        print(f"    {node}\n        {why}")
    print("  The FULL gated set is `PYFP3D_TRANSONIC_GATES=1 pytest tests/` and cost "
          "4 h 09 on 2026-08-06 at 8 threads under load. Run it at phase boundaries.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
