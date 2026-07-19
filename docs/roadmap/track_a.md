# pyFP3D Roadmap — Track A (verification & analysis A1–…)

> New track, created 2026-07-15 (user-directed). Global working rules, gate-ID
> conventions and the track index live in [roadmap.md](../roadmap.md); the
> human-readable status snapshot is [overview.md](../overview.md). Evidence
> sections: [demo_report/track_a.md](../demo_report/track_a.md).

## Track A — Verification & analysis (created 2026-07-15; A1 ✓, A2 ✓, A3 ✓ CLOSED 2026-07-18 — this header had stale-read "no A3 scoped" until 2026-07-19)

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
  landed, not a regression. *(A3 erratum 2026-07-18: this reconciliation was
  correct WHEN WRITTEN — at the 15-MAC M2 body, +2 tests. The same-day
  far-field enlargement added one more lock, so the final same-day account is
  **+3 → 399 → 406**, which is what overview.md's lineage and every later
  baseline build on. Trust the lineage, not this line.)*
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
  *(Historical arbitration record: B9 duly closed 2026-07-17.)*
- No new from-scratch solves anywhere in A2; interventions are warm-start
  fixed-Γ only (coarse mandatory, medium behind `PYFP3D_A2_MEDIUM=1`).
- Dead routes stay closed: Γ_smooth is a *diagnostic input* under the T3/E
  precedent, NOT the refuted spanwise-smoothing fix; no B8 constraint-side
  tip cures; no G1.6 fix routes.
- States = the M6 wing family (P5/B7/A1 caches). Wing-body is B9 scope;
  P13 round-tip states can join at closure if needed.

### A3 — Response to the 2026-07-17 independent inspection: docs consistency + cross-path hardening ✓ CLOSED 2026-07-18 (opened same day, user-directed)

**Trigger.** Kimi Code CLI filed three independent audits on 2026-07-17
(`docs/inspection/`): a documentation-consistency review (17 findings), a code
review (C1–C7, T1/T2, P1/P2 + minors), and a roadmap/plan assessment. Every
finding was re-verified against the current tree on 2026-07-18 before any edit;
**all of them still stood** (only overview.md's baseline had been fixed in
passing by B16–B18). User arbitration: do all three layers, as a Track A phase.

**Scope guard.** A3 is a HARDENING phase: no converged result on any committed
mesh may move. Every code change is either dormant-path (fail-fast / disarm),
byte-inert by default (a new kwarg, a gamma default), or a loud guard that is a
no-op on committed data — each verified as such, not assumed.

- [x] **GA3.1 — documentation consistency: all 17 findings dispositioned.**
  Fixed: PROJECT_STRUCTURE footer (stale "B9 = NEXT", missing B14/B16–B18,
  baseline 421→463) + its four directory trees (meshgen/`fuselage.py`,
  `wingbody.py`, solve/`timing.py`, five missing mesh families, the dead
  `cases/test_*.py` line dropped, demo/tests trees re-pointed at their
  authoritative listings); **M2 closed** in track_m (title, solver-leg
  checkbox, ledger) matching agent-rules/track_b; overview snapshot date;
  GB15.3 timings → the committed CSV truth **41.9 / 7.5 (5.6×) / 6.5 s**
  (demo_report; the dated analysis snapshot got an erratum instead of a
  rewrite); a real **"Known gaps"** heading in PROJECT_STRUCTURE so CLAUDE.md,
  design.md §5.1.2 and four internal pointers resolve; four stale
  P9-means-curved-elements references annotated; renumber counts 2→3 per
  track; `track_a` added to three doc maps; four stale headers/ledger lines;
  both `cases/*/README.md` tables completed (13 + 3 missing directories);
  GA1.2's superseded 398/405 reconciliation erratum'd to the final 399/406;
  "42.6% → 2.6%" → **43.6% → 2.6%** with A1's independent 42.6% named
  separately; `docs/references/` marked gitignored; M2 regen cost unified;
  hard-coded `/home/lrz/code/UP3D` de-personalized.
  **Two deliberate deviations, both cost-discipline:** (a) the V14.6 dual pair
  (0.17/0.36 vs 0.15/0.34) is documented at its source rather than
  "fixed" — switching `p14/run_demo.py`'s rounded `LS_REF` to a CSV read would
  move the committed `checks.csv` and force a re-run of a heavy demo for a
  0.02pp rounding difference; (b) the dated `docs/analysis/` snapshot keeps its
  original prose plus an erratum, per the "dated snapshots are not maintained
  documents" rule.
- [x] **GA3.2 — C2/C3 backported to the conforming Newton** (`solve/newton.py`).
  C2: `r_level_best` is now reset at all three SELECTION-EPOCH boundaries
  (freeze-arm, revert, refresh), mirroring `newton_ls.py`'s B15 fix. Frozen and
  live residuals are not comparable — carrying the frozen best (~1e-11) across
  a refresh makes the fail-fast read a spurious 1e8× blow-up and kill a healthy
  freeze cycle. C3: new `freeze_max_reverts=3` + a `freeze_armed` flag —
  permanent disarm after repeated diverging freezes ("a freeze may only ever
  help; it can never cost convergence"). `n_freeze_reverts` was already counted
  and reported but never read. Both paths are DORMANT on every committed run
  (no conforming level has ever reverted more than twice, and the 2.5D floors
  sit below the 1e-3 absolute fail-fast floor) ⇒ suite unchanged. Flows through
  `solve_newton_transonic` via the existing `newton_kw` dict — no plumbing.
- [x] **GA3.3 — mesh reader hardened (C4/C5) + 3 tests.** C4 (**the dangerous
  one, silent**): a `.msh` naming only SOME physical surface groups dropped the
  unnamed groups' triangles entirely; for an imported mesh with wall/wake named
  but symmetry not, the chain is: symmetry triangles vanish →
  `wake_cut._sheet_free_edge_nodes` reads the sheet's symmetry edge as an
  interior free edge → root station not duplicated → **Γ(root) silently pinned
  to 0**, wrong aerodynamics with no error. Now every tag present in the file
  survives, unnamed ones under `surface_<id>`. C5: physical volume tags with no
  `$PhysicalNames` crashed `Mesh.validate()` ("Tag out of range") on legal
  input; `tag_names` is now padded with `volume_<i>` placeholders. Both
  reproduced on synthetic `.msh` and **verified to FAIL before the fix**
  (`git stash` A/B); a third test locks inertness on the committed
  sphere-shell. All in-repo meshes name everything ⇒ latent traps, not live
  bugs, which is exactly why they needed tests.
- [x] **GA3.4 — guards and small fixes.** C7a: `CutElementMap` now asserts
  every TE node owns an aux DOF (silently writing `reroute[-1]` redirected the
  LAST aux DOF's mass row onto a TE main row and left that DOF equation-less);
  ★ **scoped to `n_ext_dofs > 0` after the B1 suite caught it firing on a
  legitimate degenerate synthetic case** (a lone TE-fan tet with no cut
  elements at all — no multivalued system exists there). C7b: TE nodes are
  forced into the eps side-shift, closing the 1e-3·h (TE detection) vs 1e-6·h
  (side shift) window in which a TE node could be classed "−" and have its
  UPPER value used as the LOWER reference; **measured 0 nodes in that window on
  all six committed families ⇒ bit-identical**. C6: the far-field wake-master
  branch test drops its exact-float `dy == 0.0` (mask membership already
  implies cut membership; the generators only guarantee |y| < 1e-7·r_far, so a
  master at y=1e-9 took the wrong branch and was off by one full Γ). P1: `gamma`
  threaded through the conforming `section_cp_curve`/`wall_cp_curve`/
  `_wall_section_points` so `unified.section_cp` no longer silently applies
  γ=1.4 on one branch and the caller's γ on the other — the comparability
  contract `unified` exists for. T2 (`validate_coloring` return dropped, a
  no-op assert), F0 (the load-flaky 10% wall-clock bound → 20%, and the suite's
  only CWD-relative mesh path anchored to REPO_ROOT), plus the `surface_ls`
  "x/c" docstring that returns raw x.
- [x] **GA3.5 — T1: the B3 gate is now actually locked.** `test_b3_lifting`
  asserted 2% where design_track_b's B3 gate says **0.3%** (7× looser ⇒ the
  gate's semantics were not locked by any test). Measured 2026-07-18:
  **0.1441%** (embedded 0.141376 vs wake-free 0.141580). Tightened to the
  gate's 0.3% (~2× margin) and gate text ↔ test now cite each other.
- [x] **GA3.6 — C1 VERIFIED AND CONFIRMED ★ (the phase's one real finding).**
  `cases/analysis/c1_ls_jacobian_fd/run_check.py` — B6's FD gate's 3-D twin,
  reusing its harness verbatim so it measures the SHIPPED Newton system.
  **On ONERA M6 coarse at M0.70 the LS Newton Jacobian is NOT the derivative
  of its residual on mixed-side plain elements:**

  | probe | n DOFs | ‖Jv − FD‖/‖FD‖ |
  |---|---|---|
  | targeted (aux of cut nodes touching a mixed-side plain element) | 102 | **1.146e-01** |
  | control (all other aux) | 1509 | **6.33e-10** |
  | global free (= B6's own direction) | 12566 | 2.47e-03 |

  Eight orders of magnitude apart ⇒ the defect is exactly where C1 predicted
  and nowhere else. ★ **Adversarial discriminator (added after the first run,
  because "FD noise / a perturbation flipping upwind branches" is the obvious
  alternative explanation): the targeted error is eps-INDEPENDENT — 1.532e-01
  at eps = 1e-6, 1e-7 AND 1e-8, max/min = 1.00 across three decades.** FD
  non-smoothness cannot do that; a missing Jacobian term does. ★ The affected
  class is **not** just the beyond-tip strip: 3378 mixed-side plain elements vs
  428 beyond-tip. **Consequence (bounded, honest): the LS Newton degrades to a
  QUASI-Newton on 3-D meshes — R is untouched, so every converged LS state,
  every committed γ/cl and every gate number stands; what is affected is the
  convergence RATE and step count.** Why no gate saw it: B6's FD gate runs on
  the quasi-2D mesh, which structurally has no such elements, and the B15/B7 M6
  gates are convergence gates, not FD gates.
  **RECORDED, NOT FIXED** — the fix (side-aware column mapping) is a change to
  a shipped kernel that would move committed step-count trajectories, so it is
  the user's call as its own phase. Caveat recorded: measured at a seeded
  (|R| 4.8) state, which is legitimate for an FD identity that must hold at
  EVERY state, but a converged-state repeat is the natural follow-up.

**Suite:** 460 + 3 (reader) = **463 passed + 21 skipped + 2 xfailed**, **measured 1165.41 s @16 threads** (not inferred).
**No `pyfp3d/` numerics change**: every edit is dormant-path, default-inert, or
a measured no-op on committed data.

## Progress ledger

### Track A — verification & analysis

Track status: **◐ IN PROGRESS** — created 2026-07-15. **A3 ✓ CLOSED
2026-07-18** (GA3.1–GA3.6; response to the 2026-07-17 independent Kimi
inspection — 17 docs findings dispositioned, C2/C3 backported to the conforming
Newton, reader C4/C5 hardened, C6/C7/P1/T1/T2/F0 fixed, and **C1 VERIFIED AND
CONFIRMED**: the LS Newton Jacobian mismatches its residual on mixed-side plain
elements, RECORDED not fixed). A1 ✓ CLOSED 2026-07-16
(GA1.1–GA1.5 all PASS). **A2 ✓ CLOSED 2026-07-17** (GA2.1–GA2.5; S1 = a
Kutta-probe measurement-operator artifact, D=7.33/25.70 coarse/medium; S2 =
potential-jump Kutta form error + P1 recovery artifact; fix routed to new
**P14**). **A3 ◐ OPENED 2026-07-18** — response to the 2026-07-17 independent
inspection (docs/inspection/): documentation consistency + cross-path
hardening + the C1 Jacobian verification experiment. *(The "B9 stays NEXT"
that stood here is superseded: B9 closed 2026-07-17.)*

| Phase | Status | Closed on | Notes |
|-------|--------|-----------|-------|
| A3 | ✓ | 2026-07-18 | **Response to the 2026-07-17 independent Kimi inspection** (`docs/inspection/`, three audits: 17 docs-consistency findings, a code review C1–C7/T1/T2/P1/P2, a plan assessment). All findings re-verified against the current tree first — **all still stood**. **★ HEADLINE — C1 CONFIRMED by the verification the audit itself asked for.** `cases/analysis/c1_ls_jacobian_fd/run_check.py` (B6's FD gate's 3-D twin, reusing its harness verbatim), ONERA M6 coarse M0.70: targeted probe (aux DOFs of cut nodes touching a mixed-side plain element) ‖Jv−FD‖/‖FD‖ = **1.146e-01** vs control aux **6.33e-10** — eight orders apart, exactly where C1 predicted. ★ **Adversarial discriminator against the obvious alternative (FD noise / branch flipping): the error is eps-INDEPENDENT — 1.532e-01 at eps 1e-6, 1e-7 AND 1e-8, max/min = 1.00 across three decades** ⇒ a missing term, not non-smoothness. Class is **3378** mixed-side plain elements (vs 428 beyond-tip) — bigger than the audit's framing. **Bounded consequence: R is untouched ⇒ every converged LS state, γ, cl, M_max and gate number STANDS; what degrades is the convergence RATE — the LS Newton is a quasi-Newton in 3-D.** No gate could see it (B6's FD gate is quasi-2D, which structurally has no such elements; B7/B15's M6 gates are convergence gates). **RECORDED, NOT FIXED** — side-aware column mapping is a shipped-kernel change that would move committed step-count trajectories ⇒ its own phase, user's call. Caveat: measured at a seeded \|R\|=4.8 state (legitimate — the FD identity must hold at every state — but a converged-state repeat is the follow-up). **★ C2/C3 backported** to `solve/newton.py` (selection-epoch `r_level_best` reset at freeze-arm/revert/refresh; `freeze_max_reverts=3` permanent disarm) — both B15-era LS fixes that sat un-backported for three phases; dormant on every committed run. **★ Reader C4/C5** — C4 was the SILENT one: naming only some physical surface groups dropped the rest, and for an imported mesh the chain ends at **Γ(root) pinned to 0 with no error**; +3 tests, the two repro tests verified FAILING before the fix. **C6** (far-field master branch no longer depends on bit-exact `dy == 0.0`), **C7a** (TE-node aux-DOF assert; ★ scoped to `n_ext_dofs > 0` after the B1 suite caught the first version firing on a legitimate degenerate synthetic case), **C7b** (TE nodes forced into the eps shift — **measured 0 affected nodes on all six committed families**), **P1** (`gamma` threaded through the conforming section-Cp chain), **T2**, **F0** (10%→20% on the load-flaky timing bound + REPO_ROOT anchoring). **T1: the B3 gate is now actually locked** — test asserted 2% against a 0.3% gate text; measured **0.1441%**, tightened to 0.3%, gate ↔ test cross-cited. **Docs: 17/17 dispositioned** (15 fixed; 2 fixed-by-documenting — the V14.6 dual pair 0.17/0.36 vs 0.15/0.34 is rounding provenance, and switching p14's `LS_REF` to a CSV read would move committed evidence for 0.02pp; the dated analysis snapshot got an erratum not a rewrite). Includes the one authority-level contradiction (M2 solver leg ◐ vs CLOSED). **Process: close-out ritual extended to FIVE surfaces + a backport check** (CLAUDE.md step 5; agent-rules disciplines #9/#10) — the audit's findings were mostly close-out debt, concentrated in exactly the two surfaces the old two-item ritual omitted. **NO `pyfp3d/` numerics change.** Suite **463 + 21 + 2** (+3 reader). Response report: `docs/inspection/20260718-response-to-kimi-inspection.md`. ★ **Errata (2026-07-19, from the Kimi second-round audit):** (1) the **P2** physics-bounds contradiction is **BACKLOG, not fixed** — previously recorded only in the response doc; now also here and as a comment at `pyfp3d/physics/isentropic.py::validate_physics_bounds`; (2) the claim "M2 authority contradiction resolved in title, checkbox AND ledger" was overstated — the `track_m.md` ledger row (and the roadmap.md/overview.md M lines) still read ◐ until 2026-07-19. |
| A2 | ✓ | 2026-07-17 | TE/Kutta fidelity attribution: conforming Γ(z) jitter (S1) + TE Cp jump (S2), conforming vs level-set. `cases/analysis/a2_te_kutta_fidelity/` (`_metrics.py` + `run_a2.py` zero-solve ~90 s + `run_a2_interventions.py` gated), figures `a2_jitter/decay/te_gap/spike/intervention.png` + CSVs + checks.csv (22 passed) + checks_interventions.csv (4 passed). No `pyfp3d/` edits; suite untouched. **Findings (all in committed CSVs):** **S1 SETTLED — the Γ(z) jitter is a measurement-operator artifact of the per-station probe-difference Kutta target estimator, NOT flow content and NOT unclosed stations.** Proof chain: closure |F|/max|Γ| ≤ 0.6% (not H2); wall Δφ jitter 0.02–0.07× the station-Γ jitter at x/c=0.92 (confined to the TE-adjacent enforcement layer); **decisive fixed-Γ discriminator D = roughness(probe targets read back on φ(Γ_smooth)) / roughness(Γ_smooth) = 7.33 (coarse) / 25.70 (medium)** — feeding a smooth Γ back-produces essentially the cached jitter (0.0921≈0.0970, 0.0361≈0.0390), controls exact (max dev 0.00). B7 roughness reproduced digit-for-digit; medium contrast conforming 0.0390/0.0365 vs LS 0.0033 (11–12×). **S2 DECOMPOSED:** dominant = potential-jump Kutta form error (conforming-only), SAME-estimator TE gap conforming 0.318/0.228 vs LS 0.009/0.002 (34×/133×), sweep median 0.221/0.175 = P13's prose 0.14–0.22 now measured; secondary = P1 last-point recovery artifact present on BOTH paths (conf coarse spike 0.147→0.081 over 0→2 P6 passes; LS 0.086–0.107). The 2.5-D `a1_cp.png` (S2 present, S1 absent — one Kutta station) confirms the two are distinct mechanisms. **Routing:** S1/S2 fix (probe-free wall-adjacent-CV pressure-equality Kutta estimator) → new **P14** (designed-not-started); A2 implements nothing. |
| A1 | ✓ | 2026-07-16 | Solver bottleneck study. Instrumentation `pyfp3d/solve/timing.py` + additive edits to the four drivers and four ramp wrappers (canonical `timings` schema, `step_records`, per-solve linear-algebra counts, per-level `wall_s`/`timings`/`timings_total`). Benchmark `cases/analysis/a1_solver_bottleneck/` (`run_a1.py` ungated 2.5-D ~5 min, `run_a1_m6.py` gated 3-D ~52 min cold). Tests `tests/test_a1_instrumentation.py` (7); suite 403+18+2. **Headline finding (mesh-scale dependent — quote WITH the mesh): in 3-D (M6 medium M0.84) BOTH Newton methods are PRECONDITIONER dominated (39.5% conforming / 42.6% LS, even with B12/B13 lagged LU; conforming Newton is 76.3% linear algebra), while the Picard warm-start seed — the 2.5-D headline at 71%/85% — falls to 13.7%/20.2%. The 2.5-D "the seed is the cost" result is a small-problem artifact and is NOT the lever in 3-D; it is measured support for B14 (Schur+AMG) aiming at the right target.** BOTH Picard methods are linear-solve dominated in both dimensions (61.2% in 3-D), and conforming Picard runs >1 linear solve per outer (2.33 — the inner Kutta secant on the frozen matrix), an overhead the LS implicit Kutta does not carry. GA1.5 reproduces G8.2/B15 digit for digit and **found 4 harness defects in the process** (stale 250 s anchor; LS legs not on the committed M6 recipes; 2.5-D `cl_KJ` written for a 3-D wing; 3-D figures clobbering the 2.5-D ones) — see the corrections block above. |
