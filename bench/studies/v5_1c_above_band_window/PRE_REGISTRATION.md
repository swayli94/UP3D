# PRE-REGISTRATION — GV5.1c: the above-band seed — a true read of the pre-floor slope-2 window

Committed BEFORE the first execution, per discipline. Gate:
`docs/roadmap/track_v.md` V5 follow-up **GV5.1c** (2026-07-24,
user-directed). Design inputs: the committed GV5.1b VERDICT
(`../v5_1b_scaled_newton/VERDICT.md` §4–5 — the amended seeds sit at
1.00× the diagnosed floor, INSIDE the 10× floor band from iteration 0,
so no above-band contraction segment exists by construction; the window
question is reframed to an above-band-seed protocol) and the IBL-floor
diagnosis (`../v5_ibl_floor/results/findings.md` =
`docs/design_track_v.md` §13).

**Question:** does the scaled + damped augmented Newton (the GV5.1b
machinery, delivered and exact) exhibit a quadratic contraction regime —
median p ∈ [1.5, 2.5] on the F_BL max-norm sequence — on iterates ABOVE
the pre-registered floor band? "The augmented Newton converges
quadratically before the floor" is UNTESTED, not false (GV5.1b VERDICT
§5 item 2); this gate tests it, once, with seeds that genuinely sit
above the band.

## Seed: the amended protocol + a calibrated multiplicative δ perturbation

1. Base state = the GV5.1 amended protocol verbatim (the GV5.1b runner's
   `build_loose_state` IMPORTED, not mirrored: HEAD-regenerated
   loose-converged state by the committed GV3.1 recipe; wiring guard
   converged + |dcl_k0| ≤ 1e-8 vs the committed GV3.1 history, read —
   never recomputed; abort + record on failure = the GV5.1 pattern).
   The tight pack is built on the UNPERTURBED state, so every frozen
   operator / edge datum is identical to GV5.1b's.
2. Perturbation: U[:, 0] (δ) ← δ·(1+ε), a single deterministic scalar
   ε > 0, applied at the FREE BL nodes only (the inflow Dirichlet band
   untouched — perturbing a pinned row would be a boundary-data
   violation, instantly undone on the first step while contaminating
   the seed norm). Rationale recorded: at a fixed point whose scaled
   Jacobian has no exact null directions (findings Q3) the asymptotic
   contraction order is seed-direction independent; δ is chosen because
   the floor residual itself lives in the (B, δ) equations (findings
   Q5) — the perturbation directly excites the floor-carrying block.
3. Calibration: ε set by a deterministic bisection on log10(ε) ∈
   [1e-8, 1e4] (≤ 20 residual evaluations on the fixed pack) placing
   the seed F_BL max-norm inside the target window:
   - **T1 = [5e-2, 5e-1]** (≈ 3e3–3e4 × the floor band) — the primary
     seed;
   - **T2 = [5e-1, 5e0]** — the escalation seed, run ONLY if the T1
     trajectory yields < 3 above-band contraction triples (one re-seed
     max per level).
   If the bracket cannot place F_BL in the window (F_BL(ε = 1e4) below
   it, or a grossly non-monotone response — recorded): leg aborted,
   RECORDED, user adjudication. The unperturbed seed F_BL is recorded
   as the ~1× floor baseline (the GV5.1b constructional read).

## Solve (the GV5.1b machinery verbatim)

`newton_tight(scaling="rowcol", lm_damping=True, floor_stop=True;
max_iter=10, tol=1e-8, tol_abs=1e-10, line_search=True)` from the
perturbed x0. **Assembly untouched** (the GV5.1 FD verdicts stand; the
committed suite re-greens; no solver-side edit of any kind in this
gate). Protocol: medium binding, coarse recorded — the GV5.1b
convention. The k=1 standalone leg is DROPPED (its GV5.1b record
stands; the calibrated above-band seeds supersede it as window
probes).

## Gate bands

- **(a) implementation exactness (PASS/FAIL).** The committed suite
  green: `tests/test_v5_tight_scaled.py` (8) + the tight fleet (28) +
  the new seed-helper tests `tests/test_v5_above_band_seed.py`
  (synthetic maps only — no heavy compute). Live identities on the
  PERTURBED-seed J (T1): e1 ≤ 1e-12 (solve-free algebra); e2 ≤
  **max(1e-10, 10·κ₁(J)·eps)** — the cond-aware tolerance
  PRE-REGISTERED this time (the GV5.1b 2026-07-24 user adjudication,
  VERDICT §3 / §5 item 5, carried in); e3 ≤ 1e-6.
- **(b) the pre-floor slope-2 window (medium binding, coarse
  recorded).** Floor band verbatim from GV5.1b: coarse 3.16e-5 /
  medium 1.71e-5 (10× the diagnosed floors). p_i =
  log(e_i/e_{i-1}) / log(e_{i-1}/e_{i-2}) on F_BL max-norm contraction
  triples with both predecessors above the band. **PASS = ≥ 3
  above-band triples (the T1 trajectory; pooled with T2's if the
  escalation fires) with median p ∈ [1.5, 2.5]**; converged outright →
  trivially PASS. If < 3 triples after the escalation: the
  pre-registered fallback — RECORDED, quoting the available triples,
  the per-step contraction factors e_{i-1}/e_i, and the least-squares
  slope of log e_{i+1} vs log e_i over above-band pairs (quoted when
  ≥ 3 pairs exist; recorded context, NOT the binding statistic); NOT a
  gate crash. Honestly recorded either way: termination class
  (converged / floor_reached / cap), band-entry iteration, final F_BL
  vs the diagnosed floor, final merit vs the GV5.1b committed finals
  (read from `../v5_1b_scaled_newton/results/`, never recomputed), the
  μ schedule + rejection retries, the λ history.
- **(c) counts (recorded; aspirational).** N_polish to enter the floor
  band (or terminate) and N_total = N_loose + N_polish vs the
  committed loose outer counts (4 coarse / 5 medium); aspirational
  N_polish ≤ 2× loose; 3+ recorded with user adjudication (the GV5.1
  honesty note carried over). The seeds are deliberately far: N_polish
  here measures the window traversal, not production convergence.

## Fallbacks and aborts (pre-registered)

- Wiring guard failure on either leg: abort that leg, record, user
  adjudication (the GV5.1 pattern).
- Calibration failure (bracket exhausted / grossly non-monotone
  response): abort that leg, RECORDED.
- A non-finite Newton step / line-search failure from the far seed:
  the trajectory up to the failure is read for band (b); the leg is
  RECORDED with the failure class, not a crash (the probe guard +
  accept-or-least-bad machinery makes this unlikely).
- If the trajectory contracts so fast that < 3 above-band triples
  exist even after the T2 escalation: the window answer is "the
  above-band quadratic segment is shorter than the read resolution" —
  RECORDED with the contraction factors + regression slope quoted.
- V4-reopen trigger: NOT invoked by this gate's failure (the floor is
  already localized inside the BL block; this gate measures the
  approach, it does not touch the floor).
- Breaking the floor itself is OUT of scope here — it is the
  separately registered TE-band (B, δ) formulation item
  (`docs/roadmap/track_v.md` V5 **GV5.5**, registered 2026-07-24,
  user-directed).

## Out of scope

Closure/formulation edits; TE-band discretization changes (GV5.5);
GMRES and block preconditioners (GV5.4); loose-loop criterion changes;
3-D M6 tight solves (the 2.5-D strip is the GV5.1x testbed);
solver-side changes of any kind (this gate exercises the delivered
GV5.1b machinery only).

## Environment note (temporary, this session only)

Executed under a **temporary 8-thread constraint** (user-directed
2026-07-24, this session only): `NUMBA_NUM_THREADS=8 /
OMP_NUM_THREADS=8 / OPENBLAS_NUM_THREADS=8` vs the ledger-standard 16,
applied via the environment at execution time (the runner keeps the
16 defaults untouched) and recorded in summary.csv / VERDICT; wall
times are flagged NOT comparable to the 16-thread ledger entries.

## Artifacts

`run.py` (regenerates everything; the GV5.1b runner's seed/protocol
functions imported, not mirrored) ·
`results/seed_calibration_{coarse,medium}.csv` (the bisection trace) ·
`results/newton_history_{coarse,medium}.csv` (+ `_t2` variants if the
escalation fires; F blocks, merit, μ, λ, p_fbl, p_merit, ds_change,
termination) · `results/compare.csv` (vs GV5.1b / GV5.1 / loose
committed — read, never recomputed) · `results/summary.csv` (one row
per band) · VERDICT at wrap-up.
