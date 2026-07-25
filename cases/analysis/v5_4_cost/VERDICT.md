# VERDICT — GV5.4: augmented-step cost on M6 medium with the block preconditioner

**0 PASS / 1 FAIL / 17 RECORDED** (runner exit 1, as pre-registered for a
(b) FAIL). Pre-registration: `PRE_REGISTRATION.md` (committed a0c2a5b
BEFORE the first code change) + Addenda 2026-07-25 #1–#4 (each committed
BEFORE the re-execution it enabled). Execution: 2026-07-25, the
temporary 8-thread session constraint (user-directed; every wall time
below is non-comparable to the 16-thread ledger entries).

## Question (binding: docs/roadmap/track_v.md GV5.4)

> augmented step wall-time ≤ ~2× the inviscid Newton step on M6 medium
> with the block preconditioner working; measured number recorded either
> way.

## Answers (medium, binding)

| item | verdict | measured |
|---|---|---|
| (a) augmented step wall vs inviscid step wall | **RECORDED** | **22.93 s / 3.05 s = 7.53×** — ABOVE the ≤ ~2× reference (recorded either way) |
| (b) the block preconditioner working (D5) | **FAIL** | **NOT-WORKING**: rung 1 (block-Jacobi) DIVERGED (rel_res 5.75e4 at the cap, step not accepted); rung 2 (exact-BL Schur) converged 1/4 steps (rel_res 2.66e-8, then stagnated 2.07e-7 / 6.13e-5 / 2.52e-6 at the 300-iter cap) |
| (c) diagnostics | RECORDED | below |

**(a)** is the cost of the pre-registered ladder **as executed** — since
(b) failed, four of the five measured steps include non-converged GMRES
work capped at 300 iterations; a converged-GMRES step would cost more
(iterations↑) or need a stronger preconditioner (setup↑). This is
exactly the "measured number recorded either way" the gate registers.

## (c) diagnostics (medium unless noted)

- **System**: 124,216 DOFs = n_free 62,820 + n_st 166 + 6·n_s 61,230
  (n_s 10,205; the pre-registered ~124k estimate confirmed; W2 PASS).
- **Inviscid anchor (D6, in-session @8t)**: the A1-chain ramp's final
  level 0.8395 — 13 steps, mean **3.05 s/step** (assembly 0.12 + precond
  0.81 + linsolve 2.05; 1 refactor). Committed cross-check (NON-binding,
  @16t): A1 12 steps / 52.14 s ≈ 4.35 s/step.
- **Seed (D2)**: the A1 `conf_newton` chain verbatim; final level
  converged strict (12 steps), cl_p 0.26429 vs the P14-committed probe
  G8.2 lock 0.2646 (rel 0.116 % — W1 PASS; the stale A1-era reading
  0.26918 is superseded, addendum #4). ONE standalone IBL solve:
  converged=False at the floor 2.889e-6 (expected, the GV5.1 finding).
- **Seed blocks**: |F_φ| 4.92e-5 (the transpiration kick at the k=0+IBL
  state), |F_Γ| 1.80e-16 (the converged probe state is exactly
  consistent), |F_BL| 2.889e-6 (the floor). F eval 0.3 s; J assembly
  1.07 s/step.
- **Per-step phase split** (5 steps): assembly 1.07 + precond setup 2.1
  + GMRES 14.6 + residual/backtrack ~5.0 ≈ 22.9 s. GMRES ≈ 52 ms/iter
  (each rung-2 iteration = 2 lu.solve + the sparse matvecs + an AMG-φ
  cycle).
- **"Measure before Schur" answer**: splu(J_BL,BL) setup = **1.8 s**
  (61,230² surface operator) — the BL elimination is CHEAP; the Krylov
  convergence is the bottleneck, not the Schur setup. AMG-φ setup 0.2 s,
  ILU-BL 2.5 s (rung 1).
- **GMRES**: rung 1 step 1 — 241 iters, info=5, rel_res **5.75e4**
  (block-Jacobi diverges: the φ–BL off-diagonal coupling is too strong
  for a diagonal preconditioner). Rung 2 steps 2–5 — 277 iters
  **converged 2.66e-8**, then 300-cap at 2.07e-7 / 6.13e-5 / 2.52e-6
  (the reduced-system GMRES stagnates above rtol 1e-8: plain AMG on
  J_φφ does not see the Schur correction J_hB·J_BB⁻¹·J_Bh).
- **Newton trajectory**: the IBL-floor signature — merit pinned at
  9.35e-9, |F_BL| at 2.888e-6, λ ≈ 1e-4 with steps 2–5 strictly
  merit-decreasing (accepted), step 1 a least-bad fallback (the
  diverged-GMRES step).
- **W3 FD guard (the wing augmented J on the true gate state)**: per-block
  median rel error φ 8.66e-12 / Γ 7.34e-12 / BL 0.00e+00 (the BL sample
  is dominated by the Dirichlet-pinned identity rows — tip + LE-band —
  whose FD error is exactly 0; the φ/Γ medians carry the wiring verdict)
  — PASS (60-row sample: 24 φ incl. wake-adjacent, 6 Γ incl. the tip
  station, 30 BL).
- **Memory / wall @8t**: ru_maxrss 1.7 GB; seed 107 s + IBL 1032 s +
  W3 ~300 s + Newton 115 s (total leg ~26 min; the run 1475 s incl.
  coarse).
- **Coarse shakedown** (narrative only, not gated): 27,047 DOFs; W1
  recorded (cl_p 0.25579, the committed −5.35 % mesh effect −4.97 %);
  W3 PASS; rung 1 converged 2/3 steps then capped → rung 2 clean; the
  coarse adjudication "rung 2 working" (RECORDED); ratio 20.58× — the
  coarse inviscid step (0.17 s) is too small for a meaningful ratio.

## Execution narrative (four pre-registered addenda, each committed before re-execution)

1. **#1** — coarse W1 fired (cl_p −4.97 % vs the A1 anchor): the A1
   anchor is a medium number; the coarse deviation matches the committed
   pressure mesh effect (−5.35 %). W1's tolerance scoped to medium.
2. **#2** — the medium M0.70 probe seed stalled strict (the committed
   plateau signature): seed non-convergence no longer raises (GV5.3's
   measured ramp re-convergence evidence). *Moot under #3.*
3. **#3** — W1 fired again at 1.816 % under the pressure-chain-shaped
   seed: the committed probe anchor belongs to **A1's own chain** (one
   `solve_newton_transonic(NEWTON_M6_RECIPE)`); the seed chain switched
   to A1 verbatim. *Fired again at the identical cl_p — the chain was
   never the driver.*
4. **#4** — anchor archaeology: **the A1-era 0.26918 is stale**; the
   P14 cross-model table (`cross_model_medium_m084.csv`, the V14.6
   source) locks the conforming probe at **cl_p 0.2646** (G8.2 lock).
   The in-session 0.264293 = 0.116 % below the P14 lock — today's probe
   branch IS the P14-locked branch. W1 re-anchored (tolerance
   unchanged), then PASS.

## Reading

- At the real 124k size the pre-registered ladder does **not** deliver a
  working block preconditioner (D5): block-Jacobi **diverges** on the
  coupled system, and exact-BL Schur is *nearly* there — one converged
  step at 277 iters, then stagnation just above the rtol-1e-8 budget.
  The Schur DIRECTION is right (its setup is 1.8 s — cheap); its
  reduced-space preconditioner (plain AMG on J_φφ) is what falls short.
- The 7.53× is dominated by Krylov iterations (~52 ms × ~300), not by
  setup or assembly. For scale: a DIRECT augmented solve is off the
  table — the 3-D LU fill makes the inviscid 60k refactor alone 14.7 s;
  at 124k it is memory-prohibitive — the GMRES ladder, even
  non-converged, is the only viable step so far.
- This is the honest negative the gate was registered to record either
  way: the augmented Newton at wing scale needs a stronger preconditioner
  before its cost can even be read against the ≤ ~2× band.

## Follow-ups (registered, NOT opened — user adjudication)

- The design-doc's **(A,Ψ)-structured BL treatment** + a **Schur-aware
  reduced-space preconditioner** (the correction J_hB·J_BB⁻¹·J_Bh is
  what plain AMG-φ misses; the stagnation at 1e-7..1e-5 is its
  signature).
- An **Eisenstat-Walker forcing** variant (η_k instead of the fixed
  rtol 1e-8; the P8 inviscid driver itself runs EW) — a relaxed-tolerance
  inexact-Newton cost read.
- The pressure-Kutta row wiring (pre-registration §7, unchanged).

## Artifacts

`results/summary.csv`; per level: `inviscid_anchor_{level}.csv` (the
final-level step_records), `steps_{level}.csv` (per-step phases, GMRES,
λ, merit), `gmres_{level}.csv` (per-call solver records),
`fd_guard_{level}.csv`, `cost_breakdown_{level}.png`; `run.log`.
