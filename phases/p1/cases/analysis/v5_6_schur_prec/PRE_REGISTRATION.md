# PRE-REGISTRATION — GV5.6: Schur-aware reduced-space preconditioner for the augmented step

Committed BEFORE the first code change (Track V discipline). **Opening
adjudicated 2026-07-25 (user)**: this gate opens the GV5.4 registered
follow-up *"the design-doc's (A,Ψ)-structured BL treatment + a Schur-aware
reduced-space preconditioner (the correction J_hB·J_BB⁻¹·J_Bh is what plain
AMG-φ misses; the stagnation at 1e-7..1e-5 is its signature)"*
(`bench/studies/v5_4_cost/VERDICT.md` "Follow-ups"; `phases/p1/docs/roadmap/track_v.md`
GV5.4 reading). The sibling registered follow-ups (the EW-forcing variant,
the pressure-Kutta row wiring, stale-preconditioner productization) stay
registered-not-opened (§7).

## 1. Question

On the GV5.4 system verbatim (ONERA M6, TEST 2308, medium binding, the
124,216-DOF augmented (φ, Γ, BL) system at the same A1 `conf_newton` seed):

- **(a)** What does ONE augmented Newton step cost with a Schur-aware
  block-preconditioned GMRES linear solve, relative to ONE inviscid Newton
  step measured in the same session? RECORDED vs the ≤ ~2× reference band,
  **the measured number recorded either way**.
- **(b)** Does the Schur-aware preconditioner WORK (the §4 binary
  adjudication D5, verbatim from GV5.4)? PASS = a rung works; FAIL = none
  does.
- **(c)** Diagnostics RECORDED: the correction assembly (nnz/density of Ĉ,
  the D_BB fallback count, t_corr), the AMG-Ŝ setup, per-step GMRES
  statistics, the Schur splu setup re-measured, memory.

**The falsifiable hypothesis the gate tests** (the registered basis): GV5.4
rung 2 stagnated at rel_res 2.07e-7 / 6.13e-5 / 2.52e-6 vs rtol 1e-8
because its reduced-space preconditioner bdiag(AMG(J_φφ), M_Γ) cannot
represent the Schur correction J_hB·J_BB⁻¹·J_Bh. Rung 3 puts a sparsified
copy of that correction INTO the AMG matrix. If the stagnation signature
persists unchanged, the hypothesis is wrong and the verdict records that —
the honest negative is a pre-registered outcome, not a recipe failure.

## 2. Conditions, meshes, anchors

Verbatim from GV5.4 (`bench/studies/v5_4_cost/PRE_REGISTRATION.md` §2):

- **M6 medium binding**: `cases/meshes/onera_m6/medium.msh` (gitignored;
  regenerate with `cases/meshes/onera_m6/generate_onera_m6.py`). The
  augmented system ≈ 62.7k φ + 166 Γ + 61,230 BL ≈ **124k DOFs** (the
  exact counts measured and recorded, W2). **Coarse = wiring shakedown
  only** (recorded, not gated).
- Condition = TEST 2308: M 0.8395, α 3.06, Re_MAC 11.72e6, x_tr/c = 0.05
  both sides, tip band z > 0.95·b_semi.
- Probe-branch anchor (W1): the P14-committed probe G8.2 lock **cl_p
  0.2646** (`cross_model_medium_m084.csv`; the GV5.4 addendum-#4 anchor),
  1.5 % tolerance, medium-binding.
- Inviscid step anchor: measured **in-session, same driver, same branch,
  same thread count** from the seed ramp's final-level `step_records`
  (D6-protocol). The GV5.4 committed numbers (22.93 s / 3.05 s = 7.53×,
  rung stagnation signature) are quoted as the cross-check — annotated
  @8 threads, NON-binding (D7).
- Seed = the A1 `conf_newton` chain VERBATIM (GV5.4 addendum #3) + ONE
  standalone IBL solve (its ~1e-6 floor expected, recorded).

## 3. Design decisions

- **D1 — system/seed/protocol = GV5.4 verbatim.** Same runner mechanics;
  `scaling="rowcol"`, `lm_damping=False`, `floor_stop=False` (mu ≡ 0 — a
  pure equilibrated Newton step); N = 5 measured steps (`max_iter=5`; the
  IBL floor caps the loop, an early termination is recorded with the
  measured n); GMRES budget rtol 1e-8, restart 60, maxiter 5 (≤ 300
  iterations), `on_fail="return"`; per-step linsolve wall cap 1800 s; no
  stale preconditioner reuse (every measured step builds its
  preconditioner fresh — the honest per-step cost).
- **D2 — rung 3 (binding, runs first): Schur-aware AMG on the reduced
  operator.** The GV5.4 rung-2 exact-BL Schur operator verbatim —
  splu(J_BL,BL) once per step, the matrix-free reduced operator
  S = J_hh − J_hB·J_BB⁻¹·J_Bh on (φ,Γ), exact back-substitution — with
  the reduced-space preconditioner bdiag(**AMG(Ŝ_φφ)**, M_Γ), where

      Ŝ_φφ = J_φφ − Ĉ_φφ,   Ĉ = J_hB · D_BB⁻¹ · J_Bh   (assembled explicitly)

  - **D_BB** = the per-node 6×6 block-diagonal of J_BL,BL (the node-major
    6i+k layout; extracted via the (6,6)-BSR view of the CSR block),
    inverted per node. Guard: a node block with rcond < 1e-12 falls back
    to its (zero-safe) diagonal inverse; the fallback count is recorded
    (the Dirichlet-pinned identity rows invert exactly — the expected
    fallback count is 0).
  - All blocks are extracted from the row+col equilibrated system (the
    §13-diagnosis frame: the genuine scaled (A,Ψ) stiffness there is
    1e5–1e7), so Ĉ and the exact operator live in the same frame.
  - **Physics**: D_BB⁻¹ is the local BL response — the quasi-simultaneous
    interaction matrix of the SKIPPED V4 route returning as preconditioner
    algebra (not as a solver). The design-doc's "(A,Ψ)-structured"
    guidance (§13: the BL preconditioner should be organized by the (A,Ψ)
    structure) is honored by the per-node BLOCK treatment — the full local
    6-variable coupling incl. (A,Ψ) — vs a scalar diagonal. What D_BB
    drops is the inter-node FE coupling; whether that sparsification
    suffices is exactly what the gate measures.
  - **M_Γ** = the exact inverse of the (Γ,Γ) block (unchanged from rung 2).
- **D3 — escalation rung 4: block-triangular full-system preconditioner.**
  GMRES on the full assembled 124k operator (cheap sparse matvec) with the
  block upper-triangular preconditioner P⁻¹: y_B = lu.solve(J_BL,BL, r_B)
  (the same per-step splu); y_h = P_hh(r_h − J_hB·y_B) with P_hh =
  bdiag(AMG(Ŝ_φφ), M_Γ). One application = 1 lu.solve + 1 AMG cycle +
  sparse matvecs (≈ rung 2's per-iteration cost), now carrying BOTH the
  J_hB direction and the Schur-aware φ cycle. Escalation rule = GV5.4's:
  any measured step's GMRES fails to converge within budget, a setup
  failure, or the linsolve cap trips ⇒ the current rung is adjudicated
  NOT-working and the remaining steps move up the ladder. Both rungs'
  numbers are recorded.
- **D4 — adjudication = the GV5.4 D5 binary verbatim.** A rung WORKS iff
  GMRES converges (info = 0) within budget on EVERY measured Newton step
  of that rung AND every computed step is accepted by the safety
  backtracking AND the W3 FD guard passed. Item (b) = PASS iff any rung
  works (the working rung reported); FAIL = neither does.
- **D5 — measurement = the GV5.4 D7/D8 protocol verbatim.** Per step:
  preconditioner setup + GMRES timed INSIDE the callback; assembly timed
  separately (2 samples); the residual + backtracking share = step wall
  minus the timed pieces. ratio = mean(augmented step wall) /
  mean(inviscid final-level step wall), both in-session. RECORDED with
  the within/above annotation vs ≤ ~2×. No PASS/FAIL on the ratio.
- **D6 — NO library change.** The ladder lives in the case runner; the
  committed GV5.4 `step_solve` injection point
  (`viscous/tight_driver.py:scaled_damped_step`/`newton_tight`) is reused
  as-is. The `pyfp3d/` tree diff = empty; the suite baseline is unchanged
  (652 passed + 25 skipped + 2 xfailed); W4 guards this.
- **D7 — threads.** The runner default **16** (the AGENTS.md cap; the
  temporary 8-thread session constraint of the GV5.1b–GV5.4 executions is
  not re-imposed for this gate). Wall times are comparable to the
  16-thread ledger entries; the GV5.4 8-thread numbers are quoted as a
  NON-binding cross-check only. The ratio is internally consistent (both
  sides measured in-session).

## 4. Wiring guards (raise = recipe error, NOT verdicts)

- **W1** — the probe seed at M 0.8395: the final level converged strict
  AND (medium-binding) |cl_p − 0.2646| ≤ 1.5 % vs the P14 probe G8.2 lock;
  coarse recorded with the committed −5.35 % mesh-effect cross-check
  (verbatim from GV5.4 addenda #1/#4).
- **W2** — pack DOF counts (medium): n_st = 166, 6·n_s = 61,230 (n_s =
  10,205), n_free ∈ [60k, 64k].
- **W3** — the 60-row sampled FD guard on the assembled wing augmented J
  at the seed (the GV5.4 sample: 24 φ incl. wake-adjacent, 6 Γ incl. the
  tip station, 30 BL; the h = 1e-5/1e-6/1e-7 ladder + the fallback row
  mask), per-block median relative error < 1e-6.
- **W4** (CI-side) — the `pyfp3d/` tree is untouched (git diff outside
  `cases/analysis/v5_6_schur_prec/` and the ledger files is empty): no
  suite re-baseline.

## 5. Verdict items + exit code

- (a) the cost ratio — RECORDED (annotation: within / above the ≤ ~2×
  reference).
- (b) the preconditioner adjudication — WORKING (the rung reported) =
  PASS; NOT-WORKING = FAIL.
- (c) diagnostics — RECORDED.
- The runner exits 1 iff (b) reads FAIL; the guards raise (recipe error).
  The cost number and the GMRES statistics are recorded either way.

## 6. Execution mechanics

- Coarse shakedown first (the same protocol, cheaper), then medium
  (binding). Expected cost @16 threads: the seed ramp ≈ 2-4 min (A1
  scaling) + the standalone IBL solve (≈ 8-15 min) + the FD guard (≈
  5-10 min) + 5 augmented steps (est. 1-15 min each) per level.
- Artifacts: `VERDICT.md`, `results/summary.csv`,
  `results/steps_{coarse,medium}.csv`, `results/gmres_{coarse,medium}.csv`
  (per-call records incl. t_corr / t_amg / nnz(Ĉ) / the D_BB fallback
  count in `notes`), `results/inviscid_anchor_{level}.csv`,
  `results/fd_guard_{level}.csv`, `cost_breakdown_{level}.png`, `run.log`.
- `--levels` for partial re-runs.

## 7. Follow-ups (registered, NOT opened — user adjudication)

- Productizing the working rung into `tight_driver` (only if a rung
  works), incl. stale-preconditioner reuse across Newton steps (the GV5.4
  §7 registration, unchanged).
- The Eisenstat-Walker forcing variant (η_k instead of the fixed rtol
  1e-8; the GV5.4 registered follow-up, unchanged) — a relaxed-tolerance
  inexact-Newton cost read.
- An inter-node (A,Ψ)-structured M_BB (station-strip block inverse or an
  ILU-based approximate inverse inside Ĉ) if both rungs fall short — the
  recorded escalation route beyond rung 4.
- The pressure-Kutta row wiring into the augmented un-eliminated layout
  (the GV5.4 §7 registration, unchanged).
