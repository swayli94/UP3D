# pyFP3D Roadmap — Track A (verification & analysis A1–…)

> New track, created 2026-07-15 (user-directed). Global working rules, gate-ID
> conventions and the track index live in [roadmap.md](../roadmap.md); the
> human-readable status snapshot is [overview.md](../overview.md). Evidence
> sections: [demo_report/track_a.md](../demo_report/track_a.md).

## Track A — Verification & analysis (created 2026-07-15; IN PROGRESS — A1)

**Purpose (user-arbitrated 2026-07-15):** a home for cross-cutting *verification
and analysis* work that is not a solver/mesh feature but measures the ones that
exist — profiling, method A/B comparisons, convergence studies, cost accounting.
Track A does not add physics; it adds instrumentation and evidence. Gate IDs are
`GA<phase>.<n>`.

Unlike Tracks P/M/B/V, a Track-A phase's runnable artifact lives under
**`cases/analysis/<phase>/`** (not `cases/demo/`), because it is an analysis
study rather than a per-phase capability demo — but it obeys the same evidence
discipline: a self-checking script, committed `results/*.png` + `*.csv`,
gitignored heavy `.npz`, and **a claim without a committed artifact is not
evidence** (the 2026-07-13 audit rule).

### A1 — Solver bottleneck study: conforming vs level-set × Picard vs Newton ✓ CLOSED 2026-07-16 (opened 2026-07-15)

**Deliverable:** a common timing-instrumentation layer on all four nonlinear
drivers, plus a controlled 2×2 benchmark that measures where the wall clock goes.

- **Instrumentation** (`pyfp3d/solve/timing.py` + edits to `picard.py`,
  `newton.py`, `picard_ls.py`, `newton_ls.py` and the four Mach-ramp wrappers):
  every driver now returns a `timings` dict on ONE canonical schema
  (`seed / assembly / precond / linsolve / residual / kutta / other / wall`),
  a `step_records` list (one dict per outer/Newton iteration: residual, Γ,
  clamp counts, linear solves+iterations, per-phase seconds, cumulative wall),
  and per-solve linear-algebra counts that were previously summed away. The
  ramp wrappers add per-level `wall_s` + `timings` + `step_records` and a
  `timings_total` (killing the "timings = final level only" footgun). Before
  A1 only conforming Newton had ANY `timings`; the conforming Picard ramp kept
  no per-level record at all. **Additive and inert** — no returned field that
  existed before A1 changed (conforming Newton keeps its legacy
  `jacobian`/`amg_setup`/`gmres` keys, aliased over the canonical buckets, so
  `cases/demo/p8_newton` needs no re-run).
- **Benchmark** (`cases/analysis/a1_solver_bottleneck/`): `run_a1.py` (ungated
  2.5-D NACA coarse, minutes) + `run_a1_m6.py` (gated 3-D ONERA M6 medium
  M0.84, ~1 h, `.npz`-cached). Both wake models on the SAME wake-embedded mesh
  for the headline 2×2 (isolates method from mesh), plus a dual-mesh leg on the
  wake-free mesh to price the mesh itself. Committed recipes reused verbatim —
  A1 measures cost, it does not tune.

**Gates:** ALL CLOSED 2026-07-16.
- [x] **GA1.1 — instrumentation faithful.** `other / wall < 5%` on every run in
  `a1_runs.csv` (the accounted phases explain ≥ 95% of wall; otherwise "the
  bottleneck is phase X" is not earned). Also asserted at test scale in
  `tests/test_a1_instrumentation.py`. **PASS** — worst `pct_other` 0.7% (2.5-D)
  and 0.7% (3-D).
- [x] **GA1.2 — instrumentation inert.** Full suite still **396 passed + 18
  skipped + 2 xfailed** plus the new A1 tests; no existing return-dict key
  changed. Primary regression `pytest tests/test_v0_freestream.py` green.
  **PASS** — **403 passed + 18 skipped + 2 xfailed** in 984.83 s @16 threads
  (= 396 + the 7 new A1 tests; baseline wall 988.73 s, unchanged: A1 adds tests
  and costs nothing).
  *Baseline reconciliation:* measured against the then-current **396**. A
  concurrent Track-M session re-spec'd the M2 body later the same day (+2
  tests, 396 → **398**) and recorded the post-A1 number as **405** = 398 + 7,
  which agrees with this run. If the suite reports 405, that is both changes
  landed, not a regression.
- [x] **GA1.3 — the four methods agree, so the comparison is meaningful.** On
  the shared wake-embedded mesh at M0.5 (matched vortex far field) the cl_KJ /
  Γ spread across the four methods is < 2%, and all runs converge. Without this
  the timing table would be comparing four different answers. **PASS** — spread
  **0.322%**, 6/6 subsonic runs converged.
- [x] **GA1.4 — the bottleneck is quantified.** `a1_runs.csv` names the
  dominant timing phase and its % of wall for all four methods × both regimes;
  `a1_levels.csv` gives the Mach-ramp anatomy. Descriptive gate — the
  deliverable is the measured number, not a threshold. **PASS** — 8/8 runs
  named; see the corrected headline in the ledger.
- [x] **GA1.5 — the harness reproduces committed history** (gated 3-D leg): M6
  medium M0.84 wall times land within ±25% of the committed anchors
  (conforming Newton ≈ **145 s** per G8.2-as-amended — see the correction below;
  LS Newton ≈ 657.4 s and LS Picard ≈ 2304.7 s per B15). A harness that cannot
  re-derive the known numbers cannot be trusted on the new ones. **PASS** —
  149 s (+2.9%) / 668 s (+1.6%) / 2337 s (+1.4%), and the physics reproduces
  *exactly*: cl_p 0.26464 & cl_KJ 0.26918 = the G8.2 locks 0.2646/0.2692;
  γ 0.088338 & 98 iters = GB15.4 digit for digit; cl_KJ 0.27648 = B15
  `PICARD_M6`, still parking on the plateau (conv ✗) as recorded.

**Corrections found BY GA1.5 (2026-07-16).** The gate paid for itself: all four
defects sat in paths the green 2.5-D `checks.csv` never executes. Full account in
[demo_report/track_a.md](../demo_report/track_a.md#ga15--3-d-reproduction-gated--pass-and-it-earned-its-keep).
1. **The conf_newton anchor 250 s was STALE** — it is the pre-P10 number. P10
   promoted `intermediate_tol=1e-5` into `NEWTON_M6_RECIPE` (239.5 → 140.3 s
   solve, "gated G8.2 now ~145 s", final level identical to 4 digits), and the
   script inlines the post-promotion recipe. A faithful 149 s run scored −40%
   against the stale anchor. Corrected to 145 s, provenance now in the source.
   *Anyone quoting "G8.2 = 250 s" is quoting a superseded number.*
2. **The LS legs were not running the committed M6 recipes** — `_bench.run_ls`
   carries the 2.5-D NACA defaults (`freeze_tol=1e-6` below the churn floor so
   the freeze never arms = the P9/G9.1 wall, vs `B_NEWTON_M6_DEFAULTS` 1e-3; no
   lagged LU = the B12/B13 splu trap; Picard budget 400 vs the committed 200,
   dm 0.05 vs 0.04). Against anchors that is not a reproduction. Fixed:
   `run_a1_m6.py` calls the solvers directly with `LS_*_M6_KW`.
3. **The 3-D leg wrote the 2.5-D sectional `cl_KJ = 2Γ`** (imported `cl_kj_3d`,
   never called it) — 0.17 where the committed answer is 0.2692, plus a
   `v6_consistency` derived from it, both headed for a committed CSV.
4. **The 3-D figures overwrote the 2.5-D ones** (shared `results/`, CSVs
   namespaced, PNGs not). Now `a1_*` vs `a1_m6_*` via a `prefix`.

**Scope guards (user-arbitrated 2026-07-15):**
- Full-timings-refactor instrumentation (native in the drivers), NOT an external
  monkeypatch probe.
- 2.5-D ungated + 3-D gated; the 3-D transonic story (the B15 Picard plateau) is
  the point, but it sits behind `PYFP3D_TRANSONIC_GATES=1`.
- Same-mesh headline + dual-mesh leg (isolate the method; price the mesh).
- A1 measures cost only. It makes NO grid-convergence or physics claim; GA1.3
  checks only that the four methods agree with each other. The M6 fine mesh (the
  `precond="direct"` 4h39m/26GB trap and the LS fine escape = B14) is out of
  scope. Conforming Picard in 3-D transonic (the P5 45-75 min beast, no M0.84
  anchor) is opt-in behind `PYFP3D_A1_CONF_PICARD_3D=1`.

---

## Progress ledger

### Track A — verification & analysis

Track status: **◐ IN PROGRESS** — created 2026-07-15. A1 ✓ CLOSED 2026-07-16
(GA1.1–GA1.5 all PASS). Next A1-successor phase not yet specified — A2 is
open for the user to scope.

| Phase | Status | Closed on | Notes |
|-------|--------|-----------|-------|
| A1 | ✓ | 2026-07-16 | Solver bottleneck study. Instrumentation `pyfp3d/solve/timing.py` + additive edits to the four drivers and four ramp wrappers (canonical `timings` schema, `step_records`, per-solve linear-algebra counts, per-level `wall_s`/`timings`/`timings_total`). Benchmark `cases/analysis/a1_solver_bottleneck/` (`run_a1.py` ungated 2.5-D ~5 min, `run_a1_m6.py` gated 3-D ~52 min cold). Tests `tests/test_a1_instrumentation.py` (7); suite 403+18+2. **Headline finding (mesh-scale dependent — quote WITH the mesh): in 3-D (M6 medium M0.84) BOTH Newton methods are PRECONDITIONER dominated (39.5% conforming / 42.6% LS, even with B12/B13 lagged LU; conforming Newton is 76.3% linear algebra), while the Picard warm-start seed — the 2.5-D headline at 71%/85% — falls to 13.7%/20.2%. The 2.5-D "the seed is the cost" result is a small-problem artifact and is NOT the lever in 3-D; it is measured support for B14 (Schur+AMG) aiming at the right target.** BOTH Picard methods are linear-solve dominated in both dimensions (61.2% in 3-D), and conforming Picard runs >1 linear solve per outer (2.33 — the inner Kutta secant on the frozen matrix), an overhead the LS implicit Kutta does not carry. GA1.5 reproduces G8.2/B15 digit for digit and **found 4 harness defects in the process** (stale 250 s anchor; LS legs not on the committed M6 recipes; 2.5-D `cl_KJ` written for a 3-D wing; 3-D figures clobbering the 2.5-D ones) — see the corrections block above. |
