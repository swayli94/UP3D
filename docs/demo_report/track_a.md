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
