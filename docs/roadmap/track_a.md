# pyFP3D Roadmap — Track A (verification & analysis A1–…)

> New track, created 2026-07-15 (user-directed). Global working rules, gate-ID
> conventions and the track index live in [roadmap.md](../roadmap.md); the
> human-readable status snapshot is [overview.md](../overview.md). Evidence
> sections: [demo_report/track_a.md](../demo_report/track_a.md).

## Track A — Verification & analysis (created 2026-07-15; A1 ✓, A2 ✓; no A3 scoped)

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

### A2 — TE/Kutta fidelity: conforming Γ(z) jitter + TE Cp jump, attribution study ✓ CLOSED 2026-07-17 (opened + scaffolded 2026-07-16; GA2.1–GA2.5 all met)

**Motivation (user-flagged 2026-07-16, from committed figures).** Two symptoms:
- **S1** — the conforming spanwise circulation Γ(z) carries station-to-station
  jitter (`a1_m6_spanwise.png`, P5 `g52_spanwise_*.png`) while the level-set
  path is smooth on BOTH mesh families (B7 `gamma_of_z.png`). 3-D only: the
  quasi-2D extrusion has a single Kutta station, so there is no spanwise
  quantity to be noisy in.
- **S2** — section Cp jumps unphysically at the trailing edge
  (`g51_sections_*.png`, `m6_cp_sections.png`); present in 2.5-D too
  (`a1_cp.png`), conforming much worse than level-set (B7 `section_cp.png`).

**Positioning against prior art (cited, not re-derived).** B7 committed the
S1 *contrast* (roughness = RMS 2nd difference / range: 0.0970 conforming-P5
vs 0.0079/0.0091 LS, 11–12×) and P13 committed the *prose attributions*
(Γ jitter = "conforming Kutta-probe placement"; TE Cp gap 0.14–0.22 =
potential-jump Kutta approximation + sharp-TE P1 floor) — **neither
attribution has intervention-grade evidence**, and
`INVESTIGATION_kutta_closure.md` left the probe degeneracy (35 upper + 41
lower of 166 stations share probes with a neighbour; off-plane by O(h)) as a
known-robustness item. A2's mandate is exactly the Track-A one: turn those
prose attributions into intervention-tested, decomposed, committed evidence.
A2 implements **no fixes** (Track A adds no physics).

**Mechanism map (code facts the hypotheses are grounded in).** Conforming
Kutta = per-station Γ whose target is a φ difference at *probe nodes one
edge off the TE* (`constraints/wake.py::kutta_targets`,
`mesh/wake_cut.py::_kutta_probe_nodes` — nearest-neighbour pick, off-plane
O(h) on the swept unstructured TE), stations independent. Level-set Kutta =
implicit pressure-equality row on wall-adjacent TE control volumes
(`wake/multivalued.py`, B4); no per-station machinery exists. Section Cp =
per-wall-triangle constant tangential-gradient Cp
(`post/section_cut.py::_wall_section_points`), raw (P5) or P6-smoothed (P8).

**Pre-registered hypotheses** (thresholds fixed before the intervention runs;
do not tune after seeing numbers):
- **H1 (S1 dominant cause = target-sampling noise).** Decided by the gated
  intervention discriminator **D = roughness(kutta_targets(φ(Γ_smooth))) /
  roughness(Γ_smooth)** after a fixed-Γ warm-start re-solve with a smooth
  diagnostic Γ (the T3/E methodology — explicitly NOT the refuted
  smoothing-as-fix route): **D > 3 ⇒ H1 confirmed** (the probe estimator
  regenerates jitter from a smooth field); **D < 1.5 ⇒ H1 refuted**; else
  mixed, report the split.
- **H2 (S1 secondary = incomplete per-station closure).** Measured by the
  |F_j|/max|Γ| census on committed final states; expected minor post-P5-polish
  (small |F| + persisting jitter ⇒ the noise lives in the target definition).
- **H3 (structural contrast + h-scaling).** One shared metric implementation
  must reproduce B7's committed numbers, then extend to all cached states;
  conforming jitter expected ~O(h).
- **H4 (S2 dominant = Kutta-form error).** Conforming TE pressure gap
  |Cp_u − Cp_l| = O(0.1–0.2), weakly variant-sensitive, slow in h; level-set
  ≈ 0 under the SAME estimator (its Kutta *is* pressure equality).
- **H5 (S2 secondary = last-point recovery artifact).** The last-point spike
  (deviation from the x/c∈[0.85,0.97] trend fit) is a P1 recovery artifact:
  shrinks under P6 passes / quadratic recovery / h, and exists on BOTH paths.

**Deliverable** — `cases/analysis/a2_te_kutta_fidelity/`: `_metrics.py`
(shared metric/estimator implementations; correctness anchored by GA2.1's
reproduction, no `pyfp3d/` edits at all), `run_a2.py` (zero-solve legs:
harvests the P5/B7/A1 local `.npz` caches — L0 unified jitter table, L1 probe
census + Kutta-closure census, L2z wall-Δφ(z; x/c) decay profile, L3 TE
gap/spike decomposition; ~90 s cold, level-set operator products cached to
gitignored `results/a2_cache_*.npz`), `run_a2_interventions.py` (GATED, the
GA2.2 discriminator: control solve + smooth-Γ fixed-Γ warm-start per the
committed T-series operation; **written, not yet run** — user-arbitrated).

**Gates (ALL MET — closed 2026-07-17):**
- [x] **GA2.1 — metric unification + reproduction.** The shared roughness
  implementation reproduces B7's committed 0.0970 / 0.0079 / 0.0091 within
  ±10% and the unified table covers all cached states. *Reproduced
  digit-for-digit (0.0970/0.0079/0.0091); table adds conf Picard medium
  0.0390, conf Newton medium 0.0365, LS Newton/Picard medium 0.0033/0.0033
  (`a2_jitter.csv`) — the 11–12× contrast holds at medium, and conforming
  jitter roughly halves coarse→medium (~O(h), H3).*
- [x] **GA2.2 — S1 attribution (H1 verdict). H1 CONFIRMED at both mesh scales:
  D = 7.33 (coarse) / 25.70 (medium).** Fed a smooth diagnostic Γ (roughness
  0.0126 coarse / 0.0014 medium — the finer mesh's fit is nearly flat) into a
  fixed-Γ warm-start re-solve, then read the probe Kutta targets back off the
  new field: they came back at **0.0921 / 0.0361 — i.e. essentially the
  original cached jitter (0.0970 / 0.0390)**, D = 7.33 / 25.70 (both far above
  the pre-registered >3 confirm threshold; the discriminator *grows* with
  refinement because the smooth input flattens ~O(h) while the estimator
  reproduces the cached jitter regardless). The control (fixed-Γ = cached)
  reproduced the cached targets to max dev 0.00 at both levels (the cached
  state is an exact fixed point — T2 path check clean).
  `checks_interventions.csv`, `a2_intervention.csv/png`. **The jitter is
  manufactured by the per-station probe-difference target estimator; it is not
  flow-field content and not incomplete closure.** *Zero-solve support that
  pointed here first: wall Δφ jitter at x/c=0.92 is 0.02×–0.07× the station-Γ
  jitter (confined to the TE-adjacent enforcement layer, `a2_decay.csv`);
  single-feature probe-geometry correlations are weak (top |r| 0.11–0.14,
  `a2_probe_census.csv`) — the noise is the combined sampling geometry, not one
  feature, which is why the intervention (not the correlations) earned the
  verdict.*
- [x] **GA2.3 — S1 closure completeness (H2).** |F_j| census on all cached
  conforming states. *max|F|/max|Γ| = 0.51% (P5 coarse), 0.59% (P5 medium),
  0.00% (Newton medium) — enforcement is complete; the jitter is in the
  targets, consistent with H1 (and now proven by GA2.2).*
- [x] **GA2.4 — S2 decomposition (H4/H5).** Gap + spike tables across
  {estimator × variant × level × method}. *SAME
  estimator (section last points, raw) conforming 0.318/0.228 vs level-set
  0.009/0.002 (34×/133×, coarse/medium); conforming all-station sweep median
  0.221/0.175 — squarely the P13 prose range 0.14–0.22; LS own-CV pressure
  gap 1.4e-9/4.9e-7 (its constraint, zeroed). Last-point spike exists on
  BOTH paths (conf coarse 0.147 → 0.081 under 0→2 P6 passes; LS coarse
  0.086–0.107) — the H4/H5 split (form error vs recovery artifact) is
  visible and quantified in `a2_spike.csv`.*
- [x] **GA2.5 — evidence dossier + routing.** demo_report/track_a.md A2
  section written (2026-07-17); every prose number backed by a committed CSV.
  **Routing:** S1's fix — replace the per-station probe-difference potential-
  jump estimator with a wall-adjacent-control-volume recovered-velocity /
  pressure-equality estimator (the B4 objects, ported to the conforming cut
  mesh; kills the single-node sampling noise AND closes the S2 TE Cp gap, one
  swap) — is a change to the conforming Kutta MECHANISM ⇒ **new Track P phase
  P14 (designed-not-started, trigger recorded)**. A2 implements nothing.

**Scope guards (user-arbitrated 2026-07-16):**
- Scaffold session (2026-07-16) = design + zero-solve legs; the gated
  intervention + closure followed 2026-07-16/17 on user direction. **B9
  remains NEXT** — A2 is an analysis interlude, not a queue change.
- No new from-scratch solves anywhere in A2; interventions are warm-start
  fixed-Γ only (coarse mandatory, medium behind `PYFP3D_A2_MEDIUM=1`).
- Dead routes stay closed: Γ_smooth is a *diagnostic input* under the T3/E
  precedent, NOT the refuted spanwise-smoothing fix; no B8 constraint-side
  tip cures; no G1.6 fix routes.
- States = the M6 wing family (P5/B7/A1 caches). Wing-body is B9 scope;
  P13 round-tip states can join at closure if needed.

## Progress ledger

### Track A — verification & analysis

Track status: **◐ IN PROGRESS** — created 2026-07-15. A1 ✓ CLOSED 2026-07-16
(GA1.1–GA1.5 all PASS). **A2 ✓ CLOSED 2026-07-17** (GA2.1–GA2.5; S1 = a
Kutta-probe measurement-operator artifact, D=7.33/25.70 coarse/medium; S2 =
potential-jump Kutta form error + P1 recovery artifact; fix routed to new
**P14**). B9 stays NEXT.

| Phase | Status | Closed on | Notes |
|-------|--------|-----------|-------|
| A2 | ✓ | 2026-07-17 | TE/Kutta fidelity attribution: conforming Γ(z) jitter (S1) + TE Cp jump (S2), conforming vs level-set. `cases/analysis/a2_te_kutta_fidelity/` (`_metrics.py` + `run_a2.py` zero-solve ~90 s + `run_a2_interventions.py` gated), figures `a2_jitter/decay/te_gap/spike/intervention.png` + CSVs + checks.csv (22 passed) + checks_interventions.csv (4 passed). No `pyfp3d/` edits; suite untouched. **Findings (all in committed CSVs):** **S1 SETTLED — the Γ(z) jitter is a measurement-operator artifact of the per-station probe-difference Kutta target estimator, NOT flow content and NOT unclosed stations.** Proof chain: closure |F|/max|Γ| ≤ 0.6% (not H2); wall Δφ jitter 0.02–0.07× the station-Γ jitter at x/c=0.92 (confined to the TE-adjacent enforcement layer); **decisive fixed-Γ discriminator D = roughness(probe targets read back on φ(Γ_smooth)) / roughness(Γ_smooth) = 7.33 (coarse) / 25.70 (medium)** — feeding a smooth Γ back-produces essentially the cached jitter (0.0921≈0.0970, 0.0361≈0.0390), controls exact (max dev 0.00). B7 roughness reproduced digit-for-digit; medium contrast conforming 0.0390/0.0365 vs LS 0.0033 (11–12×). **S2 DECOMPOSED:** dominant = potential-jump Kutta form error (conforming-only), SAME-estimator TE gap conforming 0.318/0.228 vs LS 0.009/0.002 (34×/133×), sweep median 0.221/0.175 = P13's prose 0.14–0.22 now measured; secondary = P1 last-point recovery artifact present on BOTH paths (conf coarse spike 0.147→0.081 over 0→2 P6 passes; LS 0.086–0.107). The 2.5-D `a1_cp.png` (S2 present, S1 absent — one Kutta station) confirms the two are distinct mechanisms. **Routing:** S1/S2 fix (probe-free wall-adjacent-CV pressure-equality Kutta estimator) → new **P14** (designed-not-started); A2 implements nothing. |
| A1 | ✓ | 2026-07-16 | Solver bottleneck study. Instrumentation `pyfp3d/solve/timing.py` + additive edits to the four drivers and four ramp wrappers (canonical `timings` schema, `step_records`, per-solve linear-algebra counts, per-level `wall_s`/`timings`/`timings_total`). Benchmark `cases/analysis/a1_solver_bottleneck/` (`run_a1.py` ungated 2.5-D ~5 min, `run_a1_m6.py` gated 3-D ~52 min cold). Tests `tests/test_a1_instrumentation.py` (7); suite 403+18+2. **Headline finding (mesh-scale dependent — quote WITH the mesh): in 3-D (M6 medium M0.84) BOTH Newton methods are PRECONDITIONER dominated (39.5% conforming / 42.6% LS, even with B12/B13 lagged LU; conforming Newton is 76.3% linear algebra), while the Picard warm-start seed — the 2.5-D headline at 71%/85% — falls to 13.7%/20.2%. The 2.5-D "the seed is the cost" result is a small-problem artifact and is NOT the lever in 3-D; it is measured support for B14 (Schur+AMG) aiming at the right target.** BOTH Picard methods are linear-solve dominated in both dimensions (61.2% in 3-D), and conforming Picard runs >1 linear solve per outer (2.33 — the inner Kutta secant on the frozen matrix), an overhead the LS implicit Kutta does not carry. GA1.5 reproduces G8.2/B15 digit for digit and **found 4 harness defects in the process** (stale 250 s anchor; LS legs not on the committed M6 recipes; 2.5-D `cl_KJ` written for a 3-D wing; 3-D figures clobbering the 2.5-D ones) — see the corrections block above. |
