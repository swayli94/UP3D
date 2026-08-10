# VERDICT — GV5.1d: the near-band seed — does a quadratic basin exist ADJACENT to the floor?

**Executed 2026-07-24 · 2 PASS / 1 FAIL / 7 RECORDED (honest FAIL).**
Binding text: `PRE_REGISTRATION.md` (committed pre-execution, branch
`kimi/track-v5-gv5-1d`). Runner: `run.py` (regenerates every artifact).
Machinery: the GV5.1b scaled + damped augmented Newton verbatim
(scaling=rowcol + lm_damping + floor_stop, max_iter=10); assembly
untouched; no solver-side edits. Environment: temporary 8-thread
session constraint (runner default 16; wall times NOT comparable to
the 16-thread ledger entries).

## §1 Question and protocol as executed

Does a quadratic contraction basin exist in the decade(s) DIRECTLY
above the floor band — the region the GV5.1c above-band trajectories
never reached (they stalled mid-range at F_BL ~ 1e-2, 4262×/12867×
the floor at the cap)? Seeds: the GV5.1c protocol verbatim (amended
base + δ·(1+ε) at the FREE BL nodes, deterministic log10-bisection)
with the near-band windows T1 = [1e-4, 1e-3] primary / T2 =
[1e-3, 1e-2] escalation (fires only on < 3 above-band triples).

As executed: **T1 calibrated both levels; T2 never fired** (≥ 3
above-band triples on both T1 trajectories, so the pooled read is the
T1 read). Wiring guards PASS both legs (|dcl_k0| = 1.56e-12 coarse /
1.54e-9 medium ≤ 1e-8). The medium loose regen landed on the **same
4th fixed point as GV5.1c at 8 threads** (cl 0.28245999, unperturbed
F_BL 1.824e-6 = 1.07× the committed floor — the GV5.1 §4 scatter,
flagged); coarse bit-identical (cl 0.26791639, 1.00× floor).

## §2 Seeds (calibration, deterministic)

| level | ε | seed F_BL | × band | × floor | evals |
|---|---|---|---|---|---|
| coarse | 1.0e1 | 1.711e-4 | 5.42× | 54× | 4 |
| medium | 5.62e1 | 6.02e-4 | 35.2× | 351× | 6 |

Both inside T1 = [1e-4, 1e-3] — the seeds sit 1–2 decades directly
above the band, BELOW the GV5.1c stall region (~1e-2), exactly as
pre-registered. Traces: `results/seed_calibration_{coarse,medium}.csv`.

## §3 Band (a) — implementation exactness: PASS both levels

Suite pre/post execution: 49 passed (tight fleet 33 +
`test_v5_above_band_seed.py` 9 + `test_v5_near_band_seed.py` 7) at 8
threads. Live identities on the perturbed T1-seed J:

| level | e1 (≤1e-12) | e2 (cond-aware tol) | tol_e2 | e3 (≤1e-6) |
|---|---|---|---|---|
| coarse | 2.62e-16 | 1.54e-11 | 9.2e-2 | 6.73e-12 |
| medium | 2.28e-16 | 1.20e-11 | 9.7e-2 | 7.13e-12 |

The cond-aware e2 tolerance (the GV5.1b adjudication, carried in)
passes with ~12 decades of margin.

## §4 Band (b) — the near-band window read: NO basin

Floor bands verbatim: coarse 3.16e-5 / medium 1.71e-5 (10× the
diagnosed floors 3.154e-6 / 1.712e-6).

**Coarse (RECORDED):** seed 1.711e-4 → iter 1: 8.66e-5 (a λ = 0.5
capped halving — the line-search cap, not Newton asymptotics) → then
the trajectory CRAWLS: 8.7e-5 → 7.59e-5 over 9 iterations with λ
collapsing to 6e-5–8e-3 (contraction factors ≤ 0.03 dex/step). Final
F_BL 7.593e-5 = **24.07× the floor (2.4× the band)**; band entry:
**none**; termination cap. Median p = 0.35 (8 triples); regression
slope 0.15 over 9 pairs; final merit 2.065e-7 vs the GV5.1b committed
2.044e-10.

**Medium (binding, honest FAIL):** seed 6.02e-4 → iter 1: F_BL moves
AWAY from the band, 6.0e-4 → 9.84e-4 (the accepted λ = 0.5 step buys
merit 1.82e-5 → 1.09e-5 by rebalancing the blocks — F_φ drops
7.8e-4 → 4.9e-4 — while the BL block max-norm grows) → then the same
crawl: 9.8e-4 → 8.43e-4 over 9 iterations, λ 0.003–0.06, contraction
factors ≤ 0.03 dex/step. Final F_BL 8.433e-4 = **492.6× the floor**;
band entry: **none**; termination cap. **Median p = 1.17 ∉
[1.5, 2.5] → FAIL** (8 triples: 1.01, 8.45, 0.06, 0.90, 1.34, 2.86,
0.37, 5.38 — plateau-roundoff scatter, not regimes); regression slope
0.88 over 9 pairs; final merit 8.610e-6 vs the GV5.1b committed
9.074e-11.

**μ schedule:** rejection retries = 0 on both legs, μ decaying
1e-6 → 5e-11 — the damping arm inert for the THIRD time (GV5.1b,
GV5.1c, GV5.1d): the line search carries all the globalization.

The T2 escalation did not fire (≥ 3 triples on both T1 legs).

## §5 Band (c) — counts (recorded)

N_polish = 10 (both levels, term=cap); N_total = 14 coarse (regen
loose 4 + 10) / 13 medium (3 + 10) vs committed loose 4/5;
aspirational N_polish ≤ 2× loose: met medium (10 ≤ 10), NOT met
coarse (10 > 8) — recorded with user adjudication per the honesty
note; the seeds are deliberately off-point: N_polish measures the
near-band traversal, not production convergence. Wall: 32 s coarse /
76 s medium at 8 threads (idle).

## §6 Interpretation

1. **The basin hypothesis is dead.** Seeds at 5.4× (coarse) and 35×
   (medium) the floor band — directly adjacent to it, below the
   GV5.1c stall region — show NO quadratic contraction: the
   line search collapses λ to 1e-3–1e-4 within 2–3 iterations and
   the residual creeps at ≤ 0.03 dex/step. Together with GV5.1c
   (far region: capped halvings then mid-range stall at ~1e-2), the
   whole above-floor range is now measured: **no slope-2 regime
   exists anywhere between 1e4× the floor and 24× the floor.** The
   stall is not a mid-range barrier with a basin below it — the
   flat/ragged merit neighborhood EXTENDS DOWN to within ~1.5
   decades of the floor itself.
2. **The near-band stall is qualitatively sharper than the mid-range
   one.** Medium's first accepted step moves F_BL AWAY from the band
   (merit bought by block rebalance, not BL descent) — the Newton
   direction at the near-band seeds is not a BL-descent direction at
   all until the line search throttles it. This is consistent with
   the diagnosis (findings Q2): the scaled (A, Ψ) stiffness
   1e5–1e7 means the merit is extremely flat along the directions
   the step must move to reduce the TE-band (rows 0/2) residual —
   the line search, not the model, is what "converges".
3. **Consequence for the program.** Basin hunting is exhausted:
   globalization (GV5.1b), above-band seeds (GV5.1c), and near-band
   seeds (GV5.1d) all fail to descend the last ~1.5 decades. The
   floor and its flat neighborhood are one obstacle, and it is
   formulation-level: **GV5.5 (the registered TE-band (B, δ)
   formulation item — mechanically rows 0 = x-momentum and 2 =
   kinetic-energy per the diagnosis's index naming) is now the only
   registered open route for the floor itself.** The GV5.4 block
   preconditioner work addresses cost, not the floor.
4. **V4-reopen trigger: not invoked** (pre-registered; the floor is
   already localized inside the BL block).
5. **8-thread caveat:** the medium fixed point scattered AGAIN (the
   4th fixed point cl 0.28245999, same as GV5.1c at 8 threads;
   coarse bit-identical). The reads quote the committed floors; the
   medium fixed-point migration expected once GV5.5 moves the floor
   will re-open the medium trajectory comparability (the GV5.1 §4
   caveat carried forward).

## §7 Artifacts

`results/seed_calibration_{coarse,medium}.csv` (the bisection traces)
· `results/newton_history_{coarse,medium}.csv` (F blocks, merit, μ,
λ, p_fbl, p_merit, ds_change, termination, wall) ·
`results/compare.csv` (vs GV5.1b / GV5.1c / loose committed — read,
never recomputed) · `results/summary.csv` (one row per band) ·
`run.py` + `run.log` · this VERDICT. Pre-registration:
`PRE_REGISTRATION.md` (committed pre-execution).
