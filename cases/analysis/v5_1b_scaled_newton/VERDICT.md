# GV5.1b VERDICT — scaled + damped augmented Newton: machinery delivered, the window question reframed

Gate: `docs/roadmap/track_v.md` V5 follow-up **GV5.1b** (2026-07-24).
Pre-registration: `PRE_REGISTRATION.md` (committed 8b7793f before the
first execution). Runner: `run.py` (regenerates every artifact in
`results/`; seeds by the committed recipe, wiring guard PASS both
legs — coarse 1.56e-12, medium 1.31e-9). Executed 2026-07-24 on
`kimi/track-v5-gv5-1b`. Row-level record: `results/summary.csv` —
**1 PASS / 1 FAIL / 7 RECORDED** (the FAIL is on a non-pre-registered
threshold; §3, user adjudication requested).

**Answer to the gate question, up front:** the scaled + damped
Newton machinery is delivered and exact (band (a) suite green: 28
tests; the medium live-check FAIL is a threshold calibration issue on
a cond ~ 1e10 solve, not an algebra error — §3). It does NOT unlock a
quadratic window — and the window question itself is **reframed, not
answered**: the amended-protocol seeds sit at ~1× the IBL floor
(coarse F_BL = 3.154e-6 vs floor 3.154e-6; medium 1.710e-6 vs
1.712e-6), i.e. INSIDE the pre-registered 10× floor band from
iteration 0, so no above-band contraction segment exists to read a
slope from (pre-registered fallback taken, RECORDED). What the gate
does establish: (i) the pre-registered floor-reached stop works —
medium terminates cleanly at iteration 5 instead of GV5.1's 10-step
λ-collapse crawl, at the same merit (9.074e-11 vs 9.025e-11); (ii) on
the harsh k=1 seed the scaled Newton descends markedly deeper than
GV5.1 (final F_BL 3.268e-6 vs 4.726e-6, −31 %; merit 2.25e-10 vs
5.28e-10, 2.3× below); (iii) coarse polish ends slightly deeper than
GV5.1 and still descending (merit 2.044e-10 < 2.068e-10, λ_last =
0.031 not collapsed); (iv) the damping arm is essentially inert — 0
rejection-retries across all three runs, μ monotone decay — the
active ingredient is the row/column scaling, consistent with the
diagnosis that globalization alone cannot cross the IBL floor.

## 1. Configuration as executed

Per the pre-registration, no deviations: assembly bit-identical to
GV5.1 (FD verdicts stand; the committed tight regression + the
legacy-path bit-regression against the committed k1seed history both
green). Solver-internal path: per-iteration row/column 2-norm
equilibration R, C of the full augmented J (zero-safe); Levenberg
damping (R·J·C + μI) δy = −R·F via splu, δx = C·δy; deterministic μ
schedule (μ₀ = 1e-6, ×10 on rejection bounded 1e2, ÷3 on accept
bounded 1e-12); P8/P14 backtracking + probe guard unchanged;
floor-reached stop = merit relative decrease < 1e-4 over 3
consecutive accepted steps. Flags `scaling / lm_damping /
floor_stop` default OFF (legacy path bit-reproduces the GV5.1
committed histories); the gate runner passes them explicitly.
States and protocol: the GV5.1 amended protocol verbatim (HEAD-regen
loose-converged seeds, polish read; coarse recorded, medium binding;
floor bands coarse 3.16e-5 / medium 1.71e-5).

## 2. Gate bands

| band | criterion | coarse | medium (binding) | verdict |
|---|---|---|---|---|
| (a) suite identities + regression | machine-precision algebra, μ schedule, tight 28 green | green | green | **PASS** (both) |
| (a) live-seed identities (extra runner check, thresholds NOT pre-registered) | e1 ≤ 1e-12, e2 ≤ 1e-10, e3 ≤ 1e-6 | e1 2.6e-16 / e2 9.2e-11 / e3 1.5e-11 | e1 2.3e-16 / e2 **1.96e-10** / e3 3.7e-10 | coarse PASS; **medium FAIL on e2 — §3, adjudication** |
| (b) pre-floor quadratic window | ≥ 3 contractions, median p ∈ [1.5, 2.5] above the floor band | no above-band segment (seed inside band) — fallback read | same | RECORDED (fallback) |
| (b) fallback: final merit vs GV5.1 committed | recorded | 2.044e-10 < 2.068e-10; termination cap, still descending | 9.074e-11 ≈ 9.025e-11; termination **floor_reached** (iter 5) | RECORDED |
| (c) counts | N_polish ≤ 2× loose (aspirational) | 10 vs 8 — NOT met (recorded) | 5 vs 10 — met (degenerate: band-entry iter 0) | RECORDED |
| k=1 standalone | recorded, no band | n_iter 10, F_BL 3.268e-6 (−31 % vs GV5.1 k1seed), merit 2.3× below | — | RECORDED |

Wall: coarse leg 317 s, medium leg 1629 s, k=1 113 s (machine under
external load; quote flagged).

## 3. band (a) medium e2 — threshold calibration, not algebra

The e2 live check compares the μ = 0 damped splu step against the
undamped splu step on the seed J. The two matrices are mathematically
identical; adding 0.0·I introduces explicit zero entries that change
SuperLU's column pivot order, so the two solves differ by roundoff
propagated through cond(J) ~ 1e10. The measured forward difference
1.96e-10 sits far inside the condition-amplified machine floor
(cond·eps·‖x‖ ~ 1e-6·‖x‖): the identity holds in the backward-error
sense. The ≤ 1e-10 forward-error threshold was chosen at
implementation time, NOT pre-registered, and was not adjusted after
the fact (discipline). The pre-registered band (a) gate — the
machine-precision identity suite on well-conditioned systems plus the
regression fleet — is green on both levels. **Adjudication
requested:** read band (a) medium as PASS under a cond-aware
threshold (recommended; the suite remains the binding gate and the
live check is re-issued as RECORDED with the backward-error note), or
keep the FAIL as the gate's honest answer on the live check.

## 4. band (b): why there is no window at these seeds

The pre-registration sized the floor band at 10× the diagnosed
floors, expecting polish room above it. In the event, the amended
seeds sit at 1.00× the floor (they ARE the loose-final states whose
IBL residual defines the floor), so every iterate lives inside the
band: F_BL never trades above 3.16e-5 / 1.71e-5, no ≥ 3-contraction
segment exists, and the median-p read is vacuous on both legs
(recorded merit-p sequences in `compare.csv` show no slope-2 either:
medium 0.52/0.03/1.64/0.54; coarse wanders 0.03–51). The window
question is therefore **open and reframed**: measuring a pre-floor
quadratic window requires a seed whose F_BL sits ABOVE the floor band
(e.g. an early loose iterate, or the k=1 state at 4.7e-6 — still only
1.5× floor on coarse; a genuinely above-band seed wants a deliberately
early loose state or a perturbed δ*). The k=1 standalone is the
current best evidence that the scaled Newton descends when given room
(−31 % F_BL over 10 iterations, merit 2.3× below GV5.1's crawl), but
10 iterations of slow descent is not a quadratic tail.

## 5. Consequences (recorded, not decided here)

1. **The scaled solve is kept** (flag-gated): it terminates
   gracefully (floor_reached replaces the λ-collapse crawl), descends
   deeper on harsh seeds, and its equilibration is a GV5.4
   preconditioner ingredient regardless of the window question.
2. **The window question moves to a seed-above-the-floor protocol**
   (candidate GV5.1c): pre-register an above-band seed (early loose
   iterate or perturbed-δ* state), read the slope-2 window there.
   Until then, "the augmented Newton converges quadratically before
   the floor" is UNTESTED, not false.
3. **Breaking the floor itself remains formulation work**: the
   diagnosis localized the floor residual in the TE-band (B, δ)
   equations; options (TE natural-outflow discretization, closure
   regularization) are out of GV5.1x scope and queue behind user
   adjudication.
4. Damping as implemented (Levenberg diagonal) is inert at these
   states (μ never engages) — keep the code path (it is exercised by
   the suite) but do not expect it to matter until an above-band
   seed exists. The V4-reopen trigger stays parked.
5. The band (a) §3 adjudication, once given, is recorded here and in
   the summary; the runner's live-check thresholds belong in the
   pre-registration of any GV5.1c.

## 6. Artifact index

- `results/summary.csv` — 1 PASS / 1 FAIL / 7 RECORDED + wall rows.
- `results/newton_history_{coarse,medium}.csv` — polish histories
  (F blocks, merit, μ, λ, p_merit, ds_change, termination).
- `results/newton_history_coarse_k1standalone.csv` — the k=1 record.
- `results/compare.csv` — 44 rows vs GV5.1 / loose committed (read,
  never recomputed): band-(a) identities with thresholds, merits,
  F_BL vs diagnosed floors, p sequences, μ schedule, counts.
- Implementation: `pyfp3d/viscous/tight_driver.py`
  (`equilibrate_rc`, `scaled_damped_step`, μ schedule helpers,
  `FloorStop`; flags default OFF = legacy path bit-identical);
  `tests/test_v5_tight_scaled.py` (8). Regression: tight fleet
  28 passed, twice (pre/post execution).
