# Phase demo report — Track A (verification & analysis)

> New track file, created 2026-07-15. Scope, reproduce instructions and the
> honesty/evidence rule: see the [demo_report.md](../demo_report.md) index.
> Roadmap gates: [roadmap/track_a.md](../roadmap/track_a.md).

## Track A — verification & analysis (A1 ◐, opened 2026-07-15)

**What the track is.** Cross-cutting *measurement* of the machinery the other
tracks built — profiling, method A/B, cost accounting — not new physics. A1 is
the first phase: a shared timing-instrumentation layer on the four nonlinear
drivers, then a controlled 2×2 benchmark (conforming vs level-set wake) ×
(Picard vs Newton) that answers **where the wall clock goes**.

**Reproduce.** `python cases/analysis/a1_solver_bottleneck/run_a1.py` (ungated
2.5-D, ~5 min, matplotlib Agg). Exit 0 = all checks pass; `results/checks.csv`
holds each GA1 verdict. The gated 3-D leg is
`PYFP3D_TRANSONIC_GATES=1 python cases/analysis/a1_solver_bottleneck/run_a1_m6.py`
(~52 min from cold: 149 s + 668 s + 2337 s; `.npz`-cached per method, so a
re-run redraws figures in seconds — `PYFP3D_A1_RESOLVE=1` forces a re-solve).
Verdicts land in `results/checks_m6.csv`. Both legs must run at the **16-thread
cap** (`NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16`): GA1.5
compares wall clock against anchors measured at 16 threads, so a different thread
count measures SMT, not harness fidelity. Figures are namespaced `a1_*` (2.5-D)
vs `a1_m6_*` (3-D) — both legs share `results/`.

## Track A / A1 — solver bottleneck study (`cases/analysis/a1_solver_bottleneck/`, 2026-07-15)

### The instrumentation, and why it had to come first

Before A1 the four drivers were instrumented wildly unevenly: **only conforming
Newton reported any `timings` at all**, the conforming Picard Mach ramp kept *no
per-level record whatsoever*, and every driver computed its per-solve linear
iteration count and then summed it away. You could not say where a solve spent
its time. A1 adds `pyfp3d/solve/timing.py` — one canonical schema
(`seed / assembly / precond / linsolve / residual / kutta / other / wall`) — and
threads it through all four drivers and all four ramp wrappers, additively:
`step_records` (one dict per iteration), per-solve linear-algebra counts, and
per-level `wall_s` + `timings` + a ramp `timings_total`. The `other` term is the
honesty column — the wall time the named phases do **not** explain — and
**GA1.1** asserts it stays under 5%.

Inertness (**GA1.2**) is load-bearing: every change is additive, no pre-existing
return-dict key changed, conforming Newton keeps its legacy
`jacobian`/`amg_setup`/`gmres` keys (aliased over the canonical buckets) so
`cases/demo/p8_newton` needs no re-run, and the full suite is unchanged at
**396 + 18 + 2** (plus the 7 new `tests/test_a1_instrumentation.py`). The
level-set Newton residual evaluation *contains* an assembly (`_System.residual`
calls `assemble_matrix`), which is split out in place so the level-set path's
assembly is not silently charged to `residual` — without that the two wake
models would not be comparable, which is the whole point.

### GA1.3 — the four methods agree, so the comparison is meaningful

`results/a1_cp.png`, `a1_runs.csv`. On the SAME wake-embedded NACA coarse mesh
at **M∞0.50, α1.25°** (matched vortex far field) the four drivers converge to
the same circulation:

| method | Γ | cl_p | dominant phase | wall (s) |
|---|---|---|---|---|
| conforming Picard | 0.086894 | 0.1734 | linsolve (95.8%) | 6.2 |
| conforming Newton | 0.086894 | 0.1734 | **seed** (70.6%) | 2.7 |
| level-set Picard | 0.087173 | 0.17412 | linsolve (54.7%) | 1.7 |
| level-set Newton | 0.087174 | 0.17412 | **seed** (84.8%) | 2.4 |

Γ spread **0.32%** across all four (< 2%) — conforming vs level-set agree to the
known < 1%, and instrumentation perturbed none of them. This licenses reading
the timing table as a method comparison rather than four different answers.

### GA1.4 — where the wall clock actually goes (the bottleneck)

Two findings, both visible in `results/a1_time_breakdown.png` and
`a1_linear_solver.png`:

1. **Both Newton methods are dominated by their Picard warm-start seed**, not by
   Newton itself — 70.6% (conforming) / 84.8% (level-set) of subsonic wall is the
   seed. The Newton steps are cheap; the seed is the cost.
   **This does NOT survive the move to 3-D** — see GA1.5 below, where both Newton
   methods are *preconditioner* dominated and the seed falls to 13.7% / 20.2%.
   At 2.5-D coarse scale the seed dominates because everything else is trivially
   cheap; it is a small-problem artifact, not the lever. Quote this number only
   with its mesh attached.
2. **Both Picard methods are linear-solve dominated**, and conforming Picard runs
   **> 1 linear solve per outer** (the inner Kutta secant re-solving the frozen
   matrix — `lin_solves_per_outer` = 2.33 in `a1_runs.csv`), an overhead the
   level-set path's implicit Kutta does not carry (it has no Γ DOF and no secant).
   This one *does* transfer: Picard is linsolve-dominated in 3-D too (61.2%).

### Transonic (M∞0.80 ramp) — cost anatomy of a hard case

`results/a1_time_breakdown.png` (transonic panel), `a1_ramp.png`, `a1_runs.csv`.
NACA coarse, **M∞0.80, α1.25°**, each path's committed recipe (conforming Newton
= the G8.1 `precond="direct"` recipe; the Picard legs = the P4/B6 defaults):

| method | wall (s) | iterations | converged? | Γ | cl_p | dominant phase |
|---|---|---|---|---|---|---|
| conforming Picard | 160.9 | 10464 | parks (P4 stall state) | 0.1819 | 0.3572 | linsolve (38.3%) |
| **conforming Newton** | **4.7** | 57 | **✓** | 0.2338 | **0.459** | seed (57.1%) |
| level-set Picard | 30.3 | 464 | parks (B15 plateau) | 0.3008 | 0.5783 | linsolve (52.8%) |
| **level-set Newton** | **9.4** | 60 | **✓** | 0.2418 | 0.4745 | assembly (46.6%) |

The headline: **Newton is 34× faster than Picard on the conforming path AND
reaches the true discrete solution** (cl_p 0.459 = the committed G8.1 answer),
while conforming Picard burns 161 s / 10464 iterations and parks at the P4 stall
state (0.357) that P8 recorded is *not* a discrete solution. The level-set Picard
parks on the B15 shock-position plateau (γ 0.301, not converged) after its whole
budget; the level-set Newton reaches a strict solution in 9.8 s. Even the
converged Newton methods spend their time differently — the conforming one on its
Picard seed (56%), the level-set one on per-step assembly (47%, the re-linearised
cut-element operator), which is the lever for each.

At **NACA coarse M0.80** the four methods are near the transonic fold and do NOT
all reach a strict, agreeing solution — which is itself the point, and matches
the committed record (P8: "Picard states are NOT discrete solutions"; B15: the
LS Picard plateau). So the transonic leg is a **cost study of a hard case**, not
an agreement claim (GA1.3 agreement is asserted at subsonic only; transonic
convergence is *reported* in `a1_runs.csv`, not gated). `a1_convergence.png` puts
residual and Γ against **wall-clock** — the only view in which a 60-step Newton
solve and a multi-thousand-iteration Picard sweep are commensurable.

### GA1.5 — 3-D reproduction (gated) — PASS, and it earned its keep

`run_a1_m6.py`, ONERA M6 medium M0.84/α3.06, 16 threads on an idle box.
`results/a1_m6_runs.csv`, `checks_m6.csv`, `a1_m6_*.png`.

| method | wall | anchor | dev | committed cross-check | matches |
|---|---|---|---|---|---|
| conforming Newton | 149 s | 145 s | **+2.9%** | cl_p 0.26464, cl_KJ 0.26918 | G8.2 locks 0.2646 / 0.2692 |
| level-set Newton | 668 s | 657.4 s | **+1.6%** | γ 0.088338, 98 iters, conv ✓ | GB15.4 0.088338, 98 |
| level-set Picard | 2337 s | 2304.7 s | **+1.4%** | cl_KJ 0.27648, parks (conv ✗) | B15 `PICARD_M6` 0.27648, plateau |

The reproduction is **exact in the physics**, not merely inside the ±25% band:
every committed lift/circulation/iteration number comes back digit for digit, and
the two methods that should not converge (LS Picard on the B15 plateau) still
don't. GA1.1 holds in 3-D too (`pct_other` ≤ 0.7%).

**GA1.5 caught four real defects — it is the reason to have the gate.** All four
lived in paths the green 2.5-D `checks.csv` never executes:

1. **The conf_newton anchor was stale.** G8.2 closed at 249.2 s (2026-07-11), but
   P10 then PROMOTED `intermediate_tol=1e-5` into `NEWTON_M6_RECIPE` — 239.5 s →
   140.3 s solve, final level identical to 4 digits, "gated G8.2 now ~145 s"
   ([track_p.md](track_p.md) P10 A/B). The script inlines the *post*-promotion
   recipe, so 250 s was the wrong partner for it: a faithful 149 s run scored
   −40% and would have failed. Anchor corrected to 145 s **with its provenance in
   the source**, since a bare number is what went stale.
2. **The LS legs were not running the committed M6 recipes at all.**
   `_bench.run_ls` carries the 2.5-D NACA defaults, materially different on this
   mesh: `freeze_tol=1e-6` sits *below* the Mach-rising churn floor (2.7e-4 by
   M0.70) so the freeze never arms (the P9/G9.1 wall) where `B_NEWTON_M6_DEFAULTS`
   uses 1e-3; no `direct_refactor_every` means no lagged LU (the B12/B13 true-3D
   splu 18.6 s/step trap); and Picard would have used `n_outer_level=400` vs the
   committed 200 with `dm=0.05` vs 0.04 — and for a budget-bound run that parks on
   the plateau, the budget *is* the wall clock. The M6 leg now calls the solvers
   directly with the committed recipes (`LS_*_M6_KW`); `_bench`'s docstring records
   that `run_ls` is 2.5-D-only.
3. **The 3-D leg wrote a 2.5-D lift coefficient.** It imported `cl_kj_3d` and
   never called it, taking `_bench.add_forces`' sectional `cl_KJ = 2Γ` instead —
   on the M6 half wing that is meaningless (0.17 vs the true 0.2692), and
   `v6_consistency` is derived from it. Both would have been committed as
   evidence. Now integrated over the span per `post/surface.cl_kj_3d`.
4. **The 3-D figures silently overwrote the 2.5-D ones.** Both legs write into one
   `results/`, the CSVs were namespaced but the PNGs were not, so `a1_ramp.png` et
   al. — cited above as 2.5-D evidence — were replaced with M6 content. Figures are
   now `a1_*` (2.5-D) vs `a1_m6_*` (3-D) via a `prefix`.

Committed evidence is the PNG/CSV; the solutions are gitignored `.npz`.

### The 3-D bottleneck — and why the 2.5-D headline does not transfer

`results/a1_m6_time_breakdown.png`, `a1_m6_runs.csv`:

| method | seed | assembly | **precond** | linsolve | wall |
|---|---|---|---|---|---|
| conforming Newton | 13.7% | 1.4% | **39.5%** | 36.8% | 149 s |
| level-set Newton | 20.2% | 29.0% | **42.6%** | 6.5% | 668 s |
| level-set Picard | 0.0% | 26.1% | 11.1% | **61.2%** | 2337 s |

**In 3-D both Newton methods are preconditioner-dominated** — the sparse LU
factorization, 39.5% / 42.6% — *even with* B12/B13 lagged LU already lagging it
(4 and 17 refactorizations respectively). The Picard warm-start seed, the 2.5-D
headline at 71%/85%, collapses to **13.7% / 20.2%**. Add `linsolve` and the
conforming Newton run is **76.3% linear algebra** against 1.4% assembly.

So the actionable lever differs by dimension, and the 2.5-D answer is the wrong
one to optimise: at coarse 2.5-D scale the seed dominates only because the
factorization is trivially cheap there. This is measured support for **B14**
(Schur + AMG, designed-not-scheduled) being aimed at the right target — the
factorization is already the largest single phase at *medium*, before the fine
mesh where `precond="direct"` becomes the 4h39m/26GB trap outright. Picard,
by contrast, is linsolve-dominated in both dimensions (61.2% here), consistent
with the 2.5-D finding.

A1 stops at measuring this. It does not propose or cost a fix.

### Honest boundaries

- A1 measures **cost only** — no grid-convergence or physics claim. GA1.3 checks
  only that the four methods agree *with each other* at subsonic. The GA1.5
  digit-for-digit agreement with G8.2/B15 is a *reproduction* check on the
  harness, not independent confirmation that those numbers are physically right.
- **The dominant phase is mesh-scale dependent, and the 2.5-D and 3-D answers
  disagree** (seed vs precond for Newton). Neither is "the" bottleneck — quote
  each with its mesh. Only two runs per method were taken, one per regime; the
  phase splits carry ordinary run-to-run scatter (a few tenths of a percent on
  re-run) and no error bars are claimed.
- Wall times are **this box** (16C/32T Xeon 8369B, 16 threads, idle). The ±25%
  GA1.5 band is a harness-fidelity tolerance, not a portable performance spec.
- The transonic 2.5-D leg runs on the **coarse** mesh (committed, repeatable);
  the coarse fold means transonic Γ does not agree across methods and some
  methods park short of strict convergence. That is recorded, not hidden.
- The far-field BC matters: the level-set drivers' *defaults* differ (Picard
  vortex, Newton neumann) and neumann shifts Γ ≈ 4% — the harness forces a matched
  far field for the agreement check, and the M6 3-D leg uses neumann throughout
  (the B7 convention).

## Track A / A2 — TE/Kutta fidelity: Γ(z) jitter + TE Cp jump (`cases/analysis/a2_te_kutta_fidelity/`, 2026-07-16/17)

### What A2 answers

Two symptoms the user flagged on the committed figures: **S1** — the conforming
spanwise circulation Γ(z) carries station-to-station jitter (`a1_m6_spanwise.png`,
P5 `g52_spanwise_*.png`) while the level-set path is smooth on both mesh families
(B7 `gamma_of_z.png`); **S2** — section Cp jumps at the trailing edge
(`g51_sections_*.png`, `m6_cp_sections.png`), present in 2.5-D too (`a1_cp.png`),
conforming ≫ level-set (B7 `section_cp.png`). B7 had committed the S1 *contrast*
(11–12×) and P13 the *prose attributions* (probe placement; potential-jump-Kutta
gap 0.14–0.22) — neither with intervention-grade evidence. A2 turns those into
measured, decomposed, committed findings. **A2 adds no physics and edits no
`pyfp3d/` file** (the metric correctness anchor is GA2.1, below).

### Reproduce

`python cases/analysis/a2_te_kutta_fidelity/run_a2.py` (zero-solve, ~90 s cold —
harvests the P5/B7/A1 local `.npz` caches; level-set operator products cached to
gitignored `results/a2_cache_*.npz`; `results/checks.csv` = 22 PASS). The S1
verdict comes from the gated sibling
`PYFP3D_TRANSONIC_GATES=1 [PYFP3D_A2_MEDIUM=1] python …/run_a2_interventions.py`
(2 fixed-Γ warm-start solves per level, ~2 s coarse / ~14 s medium at the
16-thread cap; `results/checks_interventions.csv` = 4 PASS). Figures
`a2_jitter/decay/te_gap/spike/intervention.png`.

### GA2.1 — one metric, reproducing B7

`_metrics.roughness_d2` is B7's committed roughness (RMS 2nd difference / range)
lifted verbatim. It reproduces the committed numbers digit-for-digit
(`a2_jitter.csv`): P5 baseline **0.0970**, level-set M1 **0.00793**, M4
**0.00912**. Extending the one metric to every cached state:

| state | mesh | wake model | roughness r |
|---|---|---|---|
| conforming Picard (P5) | coarse | conforming | **0.0970** |
| conforming Picard (P5) | medium | conforming | **0.0390** |
| conforming Newton (A1/G8.2) | medium | conforming | **0.0365** |
| level-set Picard M1 (B7) | coarse | level-set | 0.00793 |
| level-set Picard M4 (B7) | coarse | level-set | 0.00912 |
| level-set Newton (A1/B15) | medium | level-set | 0.00326 |
| level-set Picard (A1/B13) | medium | level-set | 0.00328 |

The 11–12× contrast holds at medium (0.039 vs 0.0033), it is not a Picard-vs-Newton
thing (conforming Newton 0.0365 ≈ conforming Picard 0.0390), and conforming
roughness roughly halves coarse→medium — an ~O(h) artifact, not a fixed physical
feature.

### GA2.3 — the jitter is not unclosed stations (rules out H2)

Per-station Kutta residual census on the committed final states: max|F|/max|Γ|
= **0.51%** (P5 coarse), **0.59%** (P5 medium), **0.00%** (Newton medium). Every
station's secant closed. So the jitter is not stations that failed to converge —
it is in the *target* the closed secant converged to.

### GA2.2 — S1 SETTLED: the probe target estimator manufactures the jitter

The zero-solve legs point the finger: reading the wall potential jump
Δφ = φ_upper − φ_lower at x/c < 1 (a circulation proxy the wake constraint does
**not** pin — at the wake/TE node it equals the prescribed Γ identically, which is
why `kutta_targets` steps "one node off"), the jitter is a small fraction of the
station-Γ jitter and *grows toward the TE* (`a2_decay.csv`, local-fit residual /
range):

| x/c read on the wall | conf Picard coarse | conf Picard medium | conf Newton medium |
|---|---|---|---|
| 0.92 | 0.00084 | 0.00120 | 0.00120 |
| 0.98 | 0.00825 | 0.00283 | 0.00271 |
| 1.0 (= the station Γ) | 0.03493 | 0.01710 | 0.01596 |

At x/c=0.92 the field-side jitter is **0.02–0.07×** the station-Γ jitter: a few
cells up from the TE the circulation is smooth; the roughness is concentrated in
the last sliver where the probe estimator samples. Single-feature probe-geometry
correlations (probe distance, off-plane offset, shared-probe flag vs the
station's local-fit residual, `a2_probe_census.csv`) are all weak (top |r|
0.11–0.14) — no single feature predicts which station is noisiest, i.e. the noise
is the *combined* sampling geometry, which is why a correlation table cannot earn
the verdict and the intervention must.

**The decisive test** (`run_a2_interventions.py`, the P5 T3/E fixed-Γ
methodology): prescribe a *smooth* diagnostic Γ (a local-quadratic fit), warm-start
a fixed-Γ density re-solve, then read the probe Kutta targets back off the new
field. If the field carried the jitter, smooth Γ → smooth targets; if the probe
*estimator* manufactures it, smooth Γ → jittery targets. The discriminator

  D = roughness(kutta_targets(φ(Γ_smooth))) / roughness(Γ_smooth)

was pre-registered (>3 confirm / <1.5 refute) before the run:

| level | Γ_smooth in (r) | targets back (r) | cached (r) | **D** | control max-dev |
|---|---|---|---|---|---|
| coarse | 0.0126 | 0.0921 | 0.0970 | **7.33** | 0.00 |
| medium | 0.0014 | 0.0361 | 0.0390 | **25.70** | 0.00 |

The estimator hands back **essentially the original cached jitter** regardless of
how smooth the Γ that produced the field was — D = 7.33 / 25.70, both far above
the confirm threshold, and *growing* with refinement (the smooth input flattens
~O(h) while the estimator's output tracks the cached jitter). The control
(fixed-Γ = cached) reproduces the cached targets to machine zero, so the cached
state is an exact fixed point and the operation is clean (T2 path check).
`a2_intervention.png` is the one-figure proof: the read-back targets (red) trace
the cached jittery Γ (blue), not the smooth input (green).

**⇒ S1 is a measurement-operator artifact of the per-station probe-difference
potential-jump Kutta target estimator** (`constraints/wake.py::kutta_targets` +
`mesh/wake_cut.py::_kutta_probe_nodes`): each station independently samples φ at a
nearest-neighbour probe node one edge off the swept, unstructured TE — off-plane
by O(h), 35+41 of 166 stations sharing a probe with a neighbour — of a field that
varies sharply right at the TE. The level-set path has no per-station estimator
(Γ is one coupled solution mode; the Kutta is an implicit pressure-equality on
wall-adjacent control volumes), which is why it is smooth *by construction*, not
by tuning.

### GA2.4 — S2 DECOMPOSED: a Kutta-form error plus a recovery artifact

Two separable contributors, both measured (`a2_te_gap.csv`, `a2_spike.csv`):

1. **Dominant, conforming-only — the potential-jump Kutta is not pressure
   equality.** design.md 4.4 enforces equal potential jump; the physical Kutta is
   equal pressure (|q_u|²=|q_l|²). So the upper/lower recovered Cp genuinely differ
   at the TE. Under the **same** section-last-point estimator on both paths:
   conforming TE gap **0.318 / 0.228** (coarse/medium) vs level-set **0.009 /
   0.002** — 34× / 133×. The all-station conforming sweep median is **0.221 /
   0.175**, squarely P13's committed prose range 0.14–0.22 (now a measured number,
   mildly decreasing in h = a form error, not a blow-up). Level-set on its own
   control volumes closes its pressure-equality residual to 1e-9 / 5e-7.
2. **Secondary, both paths — a P1 last-point recovery artifact.** The very last
   1–2 points deviate from the x/c∈[0.85,0.97] trend fit by an amount that shrinks
   under P6 normal-gated smoothing (conforming coarse mean spike **0.147 → 0.081**
   over 0→2 passes) and is present on the level-set path too (0.086–0.107) — the
   same sliver-triangle gradient-recovery family as the G6.1 sawtooth, independent
   of the wake model.

**The two symptoms are distinct mechanisms**, and the 2.5-D `a1_cp.png` proves it:
it shows S2 (the TE Cp jump) but *no* S1 (a quasi-2D slab has one Kutta station,
so there is no spanwise quantity to be noisy in). S1 is the spanwise-sampling face
of the conforming Kutta; S2 is its potential-vs-pressure face.

### GA2.5 — routing (A2 implements nothing)

The indicated fix for S1 — and, because it is the physically-correct condition,
for S2's dominant term at the same time — is to replace the per-station
probe-difference potential-jump estimator with a **wall-adjacent control-volume
recovered-velocity / pressure-equality** estimator: the B4 objects
(`wake/multivalued.py::_build_te_control_volumes` + `te_velocities`, which recover
each side's velocity over a *consistent* TE-adjacent element set, not a single
off-TE node), ported to the conforming cut mesh. That removes the single-node
sampling degeneracy (S1) and, being equal-pressure rather than equal-potential,
closes the TE Cp gap (S2) — one estimator swap, both symptoms. This changes the
conforming solver's **Kutta enforcement mechanism**, so it is Track P work, not
Track A: filed as **P14 (designed-not-started)**. Dead routes stay closed —
spanwise-Γ smoothing (moves Γ off the self-consistent value), full-element-fan
recovery at the TE (B4: +11–15% wrong; wall-adjacent is <1%), and mere
better-probe-picking (a band-aid, not the mechanism). A2 recommends; P14 (if
opened) builds and gates.

### Honest boundaries

- A2 measures and attributes; it proves **no fix works** — the P14 route is the
  *indicated* one, to be gated on its own (Γ(z) smooth, lift unshifted, V6
  consistent), not a demonstrated cure. Porting B4's control volumes onto the
  conforming cut mesh (duplicated TE nodes, a Γ DOF, the tip region) has its own
  unexamined risks.
- The intervention is warm-start fixed-Γ, not from-scratch; it makes no recipe
  claim (the "warm-start ≠ from-scratch recipe" P5 lesson does not bite — there is
  no recipe here, only a diagnostic operation). Coarse + medium agree, which is
  the extent of the mesh-scale evidence.
- Γ_smooth is a **diagnostic input** under the P5 T3/E precedent, explicitly not
  the refuted smoothing-as-fix route.
