# PRE-REGISTRATION — GV5.1d: the near-band seed — does a quadratic basin exist ADJACENT to the floor?

Committed BEFORE the first execution, per discipline. Gate:
`docs/roadmap/track_v.md` V5 follow-up **GV5.1d** (2026-07-24,
user-directed: "开启 GV5.1d 近带种子"). Design inputs: the committed
GV5.1c VERDICT (`../v5_1c_above_band_window/VERDICT.md` — the
above-band window read: NO quadratic regime anywhere above the floor;
the clean-descent steps are line-search-capped halvings, then the
trajectory STALLS mid-range at F_BL ~ 1e-2, 3–4 decades above the
floor, never reaching the band — 4262×/12867× the floor at the cap;
its §5 follow-up: "whether a quadratic basin exists ADJACENT to the
floor = the near-band-seed follow-up question, candidate GV5.1d") and
the IBL-floor diagnosis (`../v5_ibl_floor/results/findings.md` =
`docs/design_track_v.md` §13).

**Question:** does the scaled + damped augmented Newton (the GV5.1b
machinery, delivered and exact, unchanged) exhibit a quadratic
contraction basin in the decade(s) DIRECTLY above the floor band — the
region the GV5.1c above-band trajectories never reached (they stalled
3–4 decades higher)? GV5.1c answered "no slope-2 above the floor" for
seeds ≥ ~1e4× the floor; it did NOT probe the near-band region: its
trajectories never descended below ~1.3e-2. This gate probes that
region once, with seeds calibrated to sit inside it.

## Seed: the amended protocol + the calibrated multiplicative δ perturbation (verbatim, new windows)

1. Base state = the GV5.1 amended protocol verbatim (the GV5.1b
   runner's `build_loose_state` IMPORTED via the GV5.1c runner, not
   mirrored: HEAD-regenerated loose-converged state; wiring guard
   converged + |dcl_k0| ≤ 1e-8 vs the committed GV3.1 history, read —
   never recomputed; abort + record on failure = the GV5.1 pattern).
   The tight pack is built on the UNPERTURBED state, so every frozen
   operator / edge datum is identical to GV5.1b's and GV5.1c's.
2. Perturbation: U[:, 0] (δ) ← δ·(1+ε) at the FREE BL nodes only (the
   inflow Dirichlet band untouched) — the GV5.1c perturbation verbatim
   (the GV5.1c runner's `perturb_delta` IMPORTED). The rationale
   carries in: at a fixed point whose scaled Jacobian has no exact
   null directions (findings Q3) the asymptotic contraction order is
   seed-direction independent; δ directly excites the floor-carrying
   (B, δ) block (findings Q5).
3. Calibration: the GV5.1c deterministic log10-bisection verbatim
   (`calibrate_eps` IMPORTED; bracket [1e-8, 1e4], ≤ 20 inf-safe
   evals), placing the seed F_BL max-norm inside the target window:
   - **T1 = [1e-4, 1e-3]** — the primary seed: 5.8–58× the medium
     floor band (58–585× the medium floor), 3.2–31.6× the coarse band
     (32–317× the coarse floor). This is the decade(s) directly above
     the band and BELOW the GV5.1c mid-range stall (~1e-2): the
     near-band window.
   - **T2 = [1e-3, 1e-2]** — the escalation seed, run ONLY if the T1
     trajectory yields < 3 above-band contraction triples (one re-seed
     max per level). Its upper edge touches the lower edge of the
     GV5.1c stall region; it is the resolution extension, not a new
     question.
   If the bracket cannot place F_BL in the window: leg aborted,
   RECORDED, user adjudication (the GV5.1c clause verbatim). The
   unperturbed seed F_BL is recorded as the ~1× floor baseline.

**8-thread scatter caveat (carried from GV5.1c):** the medium loose
fixed point scatters across thread counts (the GV5.1 §4 mechanism; a
4th fixed point cl 0.28245999 appeared at 8 threads). The wiring guard
still governs the seed; the floor/band ratios are quoted against the
committed diagnosed floors with the scatter flagged if it recurs.
Coarse was bit-identical across every code/env to date.

## Solve (the GV5.1b machinery verbatim)

`newton_tight(scaling="rowcol", lm_damping=True, floor_stop=True;
max_iter=10, tol=1e-8, tol_abs=1e-10, line_search=True)` from the
perturbed x0. **Assembly untouched; no solver-side edit of any kind**
(the committed suite re-greens). Protocol: medium binding, coarse
recorded. The k=1 standalone leg stays dropped (the GV5.1c clause).
Note (expectation, not a band): with the formulation floor in place,
outright convergence to tol 1e-8 below the floor is impossible — the
expected terminations are floor_reached / cap; the "converged
outright → trivially PASS" rule is kept verbatim for protocol
identity with GV5.1c.

## Gate bands

- **(a) implementation exactness (PASS/FAIL).** The committed suite
  green: the tight fleet (28) + `tests/test_v5_above_band_seed.py` (9)
  + the new `tests/test_v5_near_band_seed.py` (synthetic maps only —
  no heavy compute). Live identities on the PERTURBED-seed J (T1):
  e1 ≤ 1e-12; e2 ≤ **max(1e-10, 10·κ₁(J)·eps)** (the cond-aware
  tolerance, the GV5.1b 2026-07-24 adjudication carried in); e3 ≤
  1e-6.
- **(b) the near-band quadratic-basin window (medium binding, coarse
  recorded).** Floor band verbatim from GV5.1b: coarse 3.16e-5 /
  medium 1.71e-5 (10× the diagnosed floors). p_i =
  log(e_i/e_{i-1}) / log(e_{i-1}/e_{i-2}) on F_BL max-norm contraction
  triples with both predecessors above the band. **PASS = ≥ 3
  above-band triples (the T1 trajectory; pooled with T2's if the
  escalation fires) with median p ∈ [1.5, 2.5]**; converged outright →
  trivially PASS. If < 3 triples after the escalation: the
  pre-registered fallback — RECORDED, quoting the available triples,
  the per-step contraction factors, and the least-squares slope of
  log e_{i+1} vs log e_i over above-band pairs (quoted when ≥ 3 pairs
  exist; recorded context, NOT the binding statistic); NOT a gate
  crash. Honestly recorded either way: termination class, **band-entry
  iteration** (the key new datum — GV5.1c never entered the band from
  above; a near-band trajectory that enters the band and terminates
  floor_reached has traversed the near-band region, whatever the
  triple count), final F_BL vs the diagnosed floor, final merit vs
  the GV5.1b committed finals (read, never recomputed), the μ
  schedule + rejection retries, the λ history.
- **(c) counts (recorded; aspirational).** N_polish to enter the floor
  band (or terminate) and N_total = N_loose + N_polish vs the
  committed loose outer counts (4 coarse / 5 medium); aspirational
  N_polish ≤ 2× loose; 3+ recorded with user adjudication. The seeds
  are deliberately off-point: N_polish here measures the near-band
  traversal, not production convergence. (A fast near-band traversal
  — band entry within 2–3 iterations — is the basin signature this
  gate looks for; a slow one reprises the GV5.1c stall verdict at
  shorter range.)

## Fallbacks and aborts (pre-registered)

- Wiring guard failure on either leg: abort that leg, record, user
  adjudication (the GV5.1 pattern).
- Calibration failure (bracket exhausted / grossly non-monotone
  response): abort that leg, RECORDED.
- A non-finite Newton step / line-search failure: the trajectory up to
  the failure is read for band (b); the leg is RECORDED with the
  failure class, not a crash.
- If the trajectory contracts so fast that < 3 above-band triples
  exist even after the T2 escalation (a basin traversed in 1–2 steps):
  the window answer is "the near-band quadratic segment is shorter
  than the read resolution" — RECORDED with the contraction factors +
  regression slope quoted (the GV5.1c clause verbatim; band entry, if
  reached, is quoted as the traversal evidence).
- If the near-band trajectory STALLS without entering the band (the
  GV5.1c mid-range barrier extending down to the band's doorstep):
  that is the negative basin answer — honestly FAIL if ≥ 3 triples
  exist with median p outside [1.5, 2.5], else RECORDED per the
  fallback.
- V4-reopen trigger: NOT invoked by this gate's failure (the GV5.1c
  clause verbatim).
- Breaking the floor itself is OUT of scope here — it is **GV5.5**
  (registered 2026-07-24, user-directed; this gate measures with the
  floor in place).

## Out of scope

Closure/formulation edits; TE-band discretization changes (GV5.5);
GMRES and block preconditioners (GV5.4); loose-loop criterion changes;
3-D M6 tight solves; solver-side changes of any kind (this gate
exercises the delivered GV5.1b machinery only, exactly as GV5.1c did).

## Environment note (temporary, this session only)

Executed under a **temporary 8-thread constraint** (user-directed
2026-07-24, this session only): `NUMBA_NUM_THREADS=8 /
OMP_NUM_THREADS=8 / OPENBLAS_NUM_THREADS=8` vs the ledger-standard 16,
applied via the environment at execution time (the runner keeps the
16 defaults untouched) and recorded in summary.csv / VERDICT; wall
times are flagged NOT comparable to the 16-thread ledger entries.

## Artifacts

`run.py` (regenerates everything; the GV5.1c runner's seed/calibration
helpers and the GV5.1b runner's protocol functions imported, not
mirrored) · `results/seed_calibration_{coarse,medium}.csv` (the
bisection trace) · `results/newton_history_{coarse,medium}.csv` (+
`_t2` variants if the escalation fires; F blocks, merit, μ, λ, p_fbl,
p_merit, ds_change, termination) · `results/compare.csv` (vs GV5.1b /
GV5.1c / loose committed — read, never recomputed) ·
`results/summary.csv` (one row per band) · VERDICT at wrap-up.
