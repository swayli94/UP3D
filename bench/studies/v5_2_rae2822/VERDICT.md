# VERDICT — GV5.2: RAE2822 transonic VII vs committed experiment

**Outcome: band (b) FAIL** (medium binding: P1 terminal x_shock 0.6288
outside [0.495, 0.580], leg non-converged RECORDED; P2 loop runaway at
k = 4, RECORDED) — with the loose-recipe transonic-limit anatomy
recorded: 1/4 legs converged, every computed shock 0.06–0.10 c
DOWNSTREAM of the experimental bracket.

Executed 2026-07-24→25 under the temporary 8-thread session constraint
(wall times flagged non-comparable). Protocol, bands, failure clauses
and the three execution addenda: [PRE_REGISTRATION.md](PRE_REGISTRATION.md)
(committed d4302ce BEFORE code; addenda #1–#3 committed 62110bc /
8a50869 BEFORE the respective (re-)executions).

## 1. Question and protocol (as pre-registered)

Does the loose VII loop — the committed GV3.1 recipe verbatim
(ω = 1.0, ≤ 10 outer, tol_ds = 1e-3) with the GV3.2 transonic-point
Newton-driver protocol (NEWTON_ARGS imported, not re-invented) —
reproduce the committed RAE2822 experimental Cp at the two
dataset-labeled conditions P1 (M 0.725, α 2.55, Re 6.5e6) and P2
(M 0.73, α 3.19)? Forced transition x_tr/c = 0.03 both sides.
Coarse recorded, medium binding.

Execution-mechanics deviations from the first attempt, all
pre-registered in addenda BEFORE (re-)execution: the FP driver runs a
cheap→deep rescue chain (strict 1e-10 first; the library's Mach
continuation `solve_newton_transonic` from m_start = 0.70; the
honesty-guarded stall acceptance `accept_on_stall`) because the
single-shot Newton stalls at M ≥ 0.725 on this geometry (addenda
#2/#3 with the measured plateaus). The loose recipe, conditions,
metrics and bands are unchanged.

## 2. Meshes and band (a) — TE wedge pre-check (RECORDED, no fallback)

New family `cases/meshes/rae2822_2.5d/` (the M0-style embedded-wake
quasi-2D recipe mirrored from the NACA family): coarse 5560 nodes /
16236 tets, medium 20790 / 61494; stats + layer PNGs committed.
Geometry = the Cook/AGARD-AR-138 Table 6.1 ordinates (positive-DOWN
lower column; signature-locked: 12.1 % @ 37.9 %c thickness, 1.3 % @
75.7 %c camber; 0 contour self-intersections).

`results/te_wedge.csv`: mesh-crease wedge (A4 method, UNCUT mesh)
**9.46° coarse / 9.92° medium** vs the ordinate-fit measure
**12.91°** (x ≥ 0.95 secant sum; the gap is the reflex-camber
curvature over the fit window). The ≈6° quadratic-recovery guard
clears on both levels (`quadratic_available=True`) → the
pre-registered linear+smoothed fallback does NOT fire.

Two camber-specific library/runner defects were found and fixed
BEFORE execution (addendum #1): the `cut_wake` Kutta probe gained a
TE-wedge bisector-normal fallback (RAE2822's aft lower surface sits
ABOVE the chord line — both TE flank neighbours on the +y side;
fallback fires only where the old code raised; regression
`test_kutta_probes_cambered_te`); the runner's Cp side split switched
to the outward-normal idiom (a centroid-y split mislabels 9 aft lower
triangles as upper on coarse).

## 3. Band (b) — shock location (PASS/FAIL, medium binding): **FAIL**

`results/shock_{coarse,medium}.csv`, x_shock = windowed (x/c ∈
[0.2, 0.9]) compression-branch max dCp/dx on the loose-final upper
wall Cp; acceptance = experimental bracket ± 0.03 c:

| level  | point | x_shock | band           | in band |
|--------|-------|---------|----------------|---------|
| coarse | P1    | 0.6122  | [0.495, 0.580] | no      |
| coarse | P2    | 0.6520  | [0.520, 0.605] | no      |
| medium | P1    | 0.6288† | [0.495, 0.580] | no      |
| medium | P2    | — (§6 recipe-limit RECORDED) | [0.520, 0.605] | — |

† at the k = 10 capped (non-converged) loose-final state — the metric
definition reads the loose-final phi; the leg's non-convergence is
itself RECORDED (§5).

**FAIL**: at the binding level neither point lands in band (P1
terminal out of band; P2 unreadable). Every computed shock sits
DOWNSTREAM of its experimental bracket — by 0.06–0.08 c at P1 and
0.05–0.10 c at P2 — and the displacement grows with Mach/α. Coarse →
medium at P1 moves the shock slightly FURTHER aft (0.6122 → 0.6288),
with a deepening LE suction spike (§4): the miss is not a
coarse-mesh artifact. The direction (inviscid-family shock aft of the
viscous experiment) is the expected one for a full-potential solve
whose displacement-thickness feedback cannot pull the shock forward
at these conditions; the magnitude is the gate's honest negative
finding.

## 4. (c) Cp RMS (RECORDED; A4 medium input band ~2.5 % peak-rel u_e
annotated)

`results/cp_compare_{point}_{level}.csv` (per experimental station,
per side; computed Cp interpolated):

| level  | point | rms_upper | rms_lower |
|--------|-------|-----------|-----------|
| coarse | P1    | 0.1852    | 0.1460    |
| coarse | P2    | 0.1756    | 0.1177    |
| medium | P1    | 0.2654    | 0.1285    |

(RECORDED at the terminal state of non-converged legs where
applicable.) The residuals are an order of magnitude above anything
attributable to the A4 input band and are dominated by (i) the shock
displacement — the RMS integrates the Cp mismatch through the
misplaced jump — and (ii) the over-deep LE suction peak (medium P1
upper: diff −0.27 at x/c ≈ 0.009; the M_peak note in §5). Lower-side
agreement is markedly better than upper-side everywhere.

## 5. (d) Convergence and guards (RECORDED)

`results/convergence_{point}_{level}.csv`, `results/summary.csv`:

| level  | point | converged | outer | wall (8 thr) | fp calls (cont / stall-acc) | cl_final | M_peak (x/c) |
|--------|-------|-----------|-------|--------------|------------------------------|----------|--------------|
| coarse | P1    | **yes**   | 7     | 150 s        | 10 (1 / 2)                   | 0.9036   | 1.271 (0.155) |
| coarse | P2    | no (cap)  | 10    | 541 s        | 20 (2 / 8)                   | 0.9792   | **1.365** (0.294) — outside-envelope RECORDED |
| medium | P1    | no (cap)  | 10    | 1421 s       | 22 (2 / 10)                  | 0.9433   | **1.306** (0.034) — outside-envelope RECORDED |
| medium | P2    | **runaway k = 4** — §6 RECORDED | — | — | — | — | — |

- The loose fixed point does not contract at these transonic points:
  at the capped legs ds_change_rel oscillates (medium P1:
  0.06 → 0.49 → 0.11 → 0.32) instead of decaying, mdot_max grows
  (0.013 → 0.13), ds_neg_floored climbs (36 nodes at k = 10), and the
  IBL itself reaches the 100-iteration cap at every outer iteration.
  medium P2 diverges outright: mdot_max = 1.59 at k = 4 (the GV3.3
  transpiration-runaway class; the loud-fail guard fired; the §6
  clause reads the point RECORDED).
- coarse P2 hit an IBL `MatrixRankWarning` (exactly singular J)
  mid-loop and continued; its k = 10 IBL residual is NaN (recorded).
- The FP rescue chain (addenda #2/#3) was load-bearing everywhere:
  2 continuation cold starts per level, and 2/10 (P1 coarse) → 10/22
  (P1 medium) FP calls stall-accepted at |R| ≲ 1e-9 plateaus (Kutta
  constraint converged, zero limiter/floor activity at acceptance,
  accept_reason recorded in run.log). Strict 1e-10 remained the first
  choice at every call.
- No CL metric is registered for this gate; cl_final is quoted as
  convergence data only.

## 6. Reading

1. **Band (b): FAIL** at the binding level — the loose VII loop with
   the GV3.2 Newton protocol does not reproduce the experimental shock
   position at either RAE2822 point; shocks land 0.06–0.10 c too far
   aft, worsening with Mach/α and not improving with mesh refinement.
2. The failure is a RECIPE limit, recorded per the pre-registered
   clauses, not a code crash: 3 of 4 legs fail to converge in ≤ 10
   outer (one capped, one capped-with-singular-IBL, one runaway); the
   single converged leg (coarse P1) still shocks 0.06 c aft.
3. Anatomy for the follow-ups: the miss direction (inviscid-family
   shock aft of experiment) plus the loop's inability to contract says
   the displacement-thickness feedback is too weak/slow at these
   conditions through the loose update — consistent with the
   IBL-floor findings (GV5.1b/5.5) and motivating the tight/augmented
   path; the shock-cell Newton plateaus that forced the rescue chain
   mark the inner-solve robustness boundary at M ≥ 0.725 (2.5-D).
4. Registered-not-opened follow-ups (user adjudication): none opened
   by this gate; the data argues the next transonic-VII reads should
   come from the tight/augmented coupling, not from further loose-loop
   tuning.

## 7. Artifacts

- `cases/meshes/rae2822_2.5d/`: `rae2822.dat` (+sha256 in the
  generator header), `generate_rae2822.py`, `{coarse,medium}.msh`,
  `{coarse,medium}_stats.csv`, `{coarse,medium}_layer.png`.
- `bench/studies/v5_2_rae2822/`: PRE_REGISTRATION.md (+3 addenda),
  run.py, VERDICT.md, `results/te_wedge.csv`,
  `results/shock_{coarse,medium}.csv`,
  `results/cp_compare_{P1,P2}_{coarse,medium}.csv` (P2 medium absent —
  §6 clause), `results/convergence_{P1,P2}_{coarse,medium}.csv`
  (P2 medium absent), `results/summary.csv`. (run.log stays in the
  working tree — `*.log` is gitignored; the per-attempt fp-chain
  narrative it carries is summarized in §5 and the counts are in
  `results/summary.csv`.)
- Tests added: `tests/test_meshgen_rae2822.py` (5),
  `tests/test_p2_wake_cut.py::test_kutta_probes_cambered_te`.
- Library change: `pyfp3d/mesh/wake_cut.py` Kutta-probe
  bisector-normal fallback (fires only where the old code raised).
- Full-suite baseline: see the three ledger lines (filled with the
  baseline commit).
