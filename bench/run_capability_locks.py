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
from pathlib import Path

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
#: M1a, the seed fallback, the wing-body conforming locks, P8's G8.1/G8.2, and
#: (added 2026-08-11) the wing-body transonic ceiling -- 5 groups.
#: ★★★ 2026-08-24 重编号（裁决 ④）：节点从 `tests/test_*.py` 改为 `tests/<CLASS>/test_<C><nn>_*.py`。
#: 分层现在由**目录**承载（裁决 ⑤）:
#:     pytest tests/A tests/B     未门控套件
#:     pytest tests/C tests/E     快层（本文件）
#:     pytest tests/D             门控全集
#:     tests/F                    可选开启（F 允许有默认开启子集）
#: ★★ 这一层的五组按新分类是 **1 个 A + 1 个 C-相关 + 3 个 E** —— 使用者的反问「快层不算 A 类吗」
#: 是对的：`seed_fallback` 与 `b9` 实测**全是机制/拓扑断言，零物理值** ⇒ 纯 A。
#: ★ 而 E01/E02/E03 是「等参照到位就搬空」的中转站：CFL3D Euler 一到，E01/E03 升格 D。
LOCKS = (
    ("tests/E/test_E03_m1a_level_lock.py",
     "M1a: the in-envelope three-level convergence (re-spec'd 2026-08-05)"),
    ("tests/A/test_A30_seed_fallback_contract.py",
     "the cold-start seed fallback, incl. the real NACA M0.80 medium recovery"),
    ("tests/A/test_A28_wingbody_topology.py",
     "wing-body conforming: junction loading, the B8 lift-loss detector"),
    ("tests/E/test_E01_p8_newton_anchors.py",
     "G8.1/G8.2 conforming Newton anchors (G8.2 re-anchored 2026-08-06)"),
    #: ★ ADDED 2026-08-11, and the reason is worth stating: this lock was written in
    #: phase 3 precisely so the conforming WING-BODY transonic capability (M0.84, cl_p
    #: 0.2738, 0 clamps) has an alarm -- it had none, while the capability boundary
    #: asserted it. Leaving it out of this tier would have meant a lock that runs only
    #: in the 2 h gated set, i.e. 6x less often than the tier it belongs to. Measured
    #: 652 s standalone; with the level-set locks gone this tier was down to 252 s, so
    #: adding it lands the tier back at the ~15 min scale it was designed around.
    ("tests/E/test_E02_wingbody_transonic_ceiling.py",
     "the conforming WING-BODY transonic ceiling: M0.84 reached, 0 clamps, cl_p"),
)

#: ★ WHAT THIS TIER DOES NOT COVER, with the measured reason. Read this before treating a
#: green run as "the gated set passes" -- that conflation is what this file exists to
#: prevent, and the numbers are from the 2026-08-06 full gated run and the re-spec rounds.
NOT_COVERED = (
    ("phases/p1/tests/test_b22_ls_3d_anchors.py::test_m6_medium_ramp_anchor  [ARCHIVED]",
     "~35 min. strict-xfail since 2026-08-09, and the xfail still RUNS the solve"),
    ("tests/D/test_D03_naca0012_m080_shock.py", "~32 min, dominated by the G4.1 medium gate"),
    ("tests/D/test_D04_onera_m6_experiment.py", "45-75 min from scratch (M6 medium continuation+polish)"),
    ("phases/p1/tests/test_b14_schur_ls.py  [ARCHIVED]", "M6 medium Schur/lagged-LU arms"),
    ("phases/p1/tests/test_b15_ls_newton_freeze.py  [ARCHIVED]", "the M6 medium M0.84 freeze ramp (~515 s idle)"),
    #: ★ GS4.0 2026-08-16: was listed as `tests/test_b18_wingbody_transonic.py`, a path
    #: that has not existed since phase 3 archived the file (its 4 legs were 100 %
    #: level-set). It sat here UNMARKED while every sibling carried [ARCHIVED], so this
    #: list over-stated what it had checked. Note what actually changed underneath: the
    #: CONFORMING half of b18's capability is no longer uncovered at all -- it is the
    #: `test_b32_wingbody_conforming_transonic.py` lock added to LOCKS above.
    ("phases/p1/tests/test_b18_wingbody_transonic.py  [ARCHIVED]",
     "level-set wing-body transonic ramps; the conforming half is now IN this tier (b32)"),
    ("[ARCHIVED, except test_p8_jacobian] test_b6_transonic, test_b16/b17, test_b11, test_b19",
     "passed on 2026-08-06; cost not individually measured"),
)


def _check_node_paths():
    """Every node named in LOCKS and NOT_COVERED must actually exist.

    ★ GS4.0 (2026-08-16). The exclusion list is deliberately treated as part of
    this script's OUTPUT, not as a comment -- the file's own header says a fast
    tier mistaken for full coverage would recreate the failure it exists to
    prevent. But it was only ever PRINTED, never validated, and the 2026-08-16
    audit found `tests/test_b18_wingbody_transonic.py` still listed as a live
    exclusion after phase 3 archived it -- unmarked, while every other archived
    entry carries `[ARCHIVED]`. A rotted list is worse than no list, because it
    reads as though someone checked.

    Entries tagged `[ARCHIVED]` are exempt from existence (that IS the tag's
    meaning), and the last NOT_COVERED entry is a prose roll-up rather than a
    path, so only entries that look like a path are checked.

    ★ Two constructive-pass holes closed on purpose (the pre-registration's
    8th question -- what does this give on an extreme sample?): an EMPTY list
    would make the loop vacuously true, and a list of nothing but `[ARCHIVED]`
    entries would too. Both are asserted, so this can never pass by being empty.
    ★ Self-reference (the 8th question's other half): this reads the LOCKS /
    NOT_COVERED tuples and the filesystem. It does not scan its own source, so
    it cannot match itself -- the failure mode that hit `pgrep -f`, the G-Z
    forbidden-word guard and the G-L load guard three times this season.
    """
    root = Path(__file__).resolve().parent.parent
    checked, missing = 0, []
    for node, _why in LOCKS + NOT_COVERED:
        if "[ARCHIVED]" in node:
            continue
        path = node.split("::")[0].strip()
        if not path.endswith(".py"):          # prose roll-up entry, not a path
            continue
        checked += 1
        if not (root / path).exists():
            missing.append(path)
    #: ★ `checked >= len(LOCKS)` ALONE passes vacuously on an empty list (0 >= 0).
    #: Caught by this round's own extreme-sample test before the guard shipped --
    #: which is the argument for asking the 8th question rather than trusting that
    #: a count comparison is a floor. A tier with no locks is a bug, not a pass.
    assert LOCKS, "the fast tier has no locks -- that is a bug, not a green run"
    assert checked >= len(LOCKS), (
        f"only {checked} path-like entries checked but there are {len(LOCKS)} "
        "locks -- the check is not covering what it claims to cover")
    if missing:
        raise SystemExit(
            "★ this tier's node list has rotted -- these paths do not exist:\n  "
            + "\n  ".join(missing)
            + "\n(if a file was archived, point the entry at its archived path "
              "AND tag it [ARCHIVED], the way its siblings are.)")
    return checked


def main():
    if os.environ.get("PYFP3D_TRANSONIC_GATES") != "1":
        print("★ PYFP3D_TRANSONIC_GATES=1 is required, else every gated lock SKIPS "
              "and this reports a vacuous green.")
        return 2
    #: ★ GS4.0: validate the node list BEFORE spending ~9 min, so a rotted entry costs
    #: a second rather than a whole tier run that then reports a coverage claim it
    #: cannot back. Raises SystemExit with the offending paths.
    n_checked = _check_node_paths()
    print("capability locks, FAST tier -- run this at every close-out")
    print(f"  node list validated: {n_checked} live paths exist "
          "(archived entries are tagged and exempt)")
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
