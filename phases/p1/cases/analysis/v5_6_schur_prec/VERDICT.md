# VERDICT — GV5.6: Schur-aware reduced-space preconditioner for the augmented step

**0 PASS / 1 FAIL / 17 RECORDED** (runner exit 1, as pre-registered for a
(b) FAIL). Pre-registration: `PRE_REGISTRATION.md` (committed `091f9fe`
BEFORE the first code change). Execution: 2026-07-25, the runner default
**16 threads** (the temporary 8-thread session constraint of the
GV5.1b–GV5.4 executions not re-imposed; the GV5.4 8-thread numbers quoted
as a NON-binding cross-check only). One clean execution — no addenda (the
wiring guards passed first try at both levels).

## Question (binding: phases/p1/docs/roadmap/track_v.md GV5.6)

> The GV5.4 system/seed/protocol verbatim: does the Schur-aware
> block-preconditioned GMRES linear solve WORK (the GV5.4 D5 binary), and
> what does one augmented Newton step cost vs one inviscid Newton step
> measured in-session (RECORDED vs ≤ ~2×, either way)?

## Answers (medium, binding)

| item | verdict | measured |
|---|---|---|
| (a) augmented step wall vs inviscid step wall | **RECORDED** | **23.88 s / 3.03 s = 7.87×** — ABOVE the ≤ ~2× reference (recorded either way; every step carries capped-GMRES work) |
| (b) the Schur-aware preconditioner working (D5) | **FAIL** | **NOT-WORKING**: rung 3 (Schur-aware AMG) — GMRES info=5 at 242/300 iters, rel_res **6.64e-01** on the only rung-3 step; rung 4 (block-triangular) — info=5 at the 300 cap on all 4 remaining steps, rel_res 6.82e-01 / 7.63e-01 / 8.62e-01 / 1.06e+00 (monotonically worsening) |
| (c) diagnostics | RECORDED | below |

## (c) diagnostics (medium unless noted)

- **System**: 124,216 DOFs = n_free 62,820 + n_st 166 + 6·n_s 61,230
  (n_s 10,205 — W2 PASS; the GV5.4 counts reproduced exactly).
- **Inviscid anchor (in-session @16t)**: the A1-chain ramp's final level
  0.8395 — 13 steps, mean **3.03 s/step** (assembly 0.11 + precond 0.81 +
  linsolve 2.04; 1 refactor). Cross-checks (NON-binding): GV5.4 3.05 s
  @8t (essentially thread-insensitive in this regime), A1 4.35 s @16t
  committed.
- **Seed (D1)**: the A1 `conf_newton` chain verbatim; final level
  converged strict (12 steps), cl_p 0.26429 vs the P14 probe G8.2 lock
  0.2646 (rel 0.116 % — W1 PASS). ONE standalone IBL solve:
  converged=False at the floor 2.875e-6 (expected, the GV5.1 finding).
- **Seed blocks**: |F_φ| 4.92e-5, |F_Γ| 2.12e-16, |F_BL| 2.875e-6 (the
  floor). F eval 0.6 s; J assembly 1.3 s/step (2 samples).
- **The correction assembly (the new diagnostics)**: nnz(Ĉ) =
  **597,154**, nnz(Ŝ_φφ) = **1,560,192** (Ĉ_φφ roughly doubles the φ-block
  density); t_corr = **0.2 s** (D_BB inversion + the triple product —
  negligible), t_amg(Ŝ) = 0.2 s, t_lu = **1.8 s** (the GV5.4
  "measure before Schur" number reproduced); **n_fallback = 0** — every
  one of the 10,205 per-node 6×6 blocks inverts exactly (the
  pre-registered expectation; the pinned identity rows included).
- **GMRES**: rung 3 step 1 — 242 iters, info=5, rel_res **6.64e-01**.
  For reference, GV5.4 rung 2 (plain AMG(J_φφ) on the same reduced
  operator, same seed) CONVERGED its first step (2.66e-8 @277 iters) and
  stagnated at 2.07e-7 / 6.13e-5 / 2.52e-6 afterwards. **The Schur-aware
  AMG matrix is catastrophically worse than plain AMG(J_φφ) at medium** —
  the pre-registered falsifiable hypothesis is falsified in this form
  (see the reading). Rung 4 steps 2–5 — 300/300/300/265 iters, info=5,
  rel_res 0.68 → 1.06 (worsening).
- **Per-step phase split** (5 steps): assembly 1.3 + precond setup 2.2 +
  GMRES 16.8 + residual/backtrack ~3.5 ≈ 23.9 s. The setup cost of the
  Schur-aware machinery is NOT the problem (2.2 s vs GV5.4 rung 2's
  2.1 s) — the Krylov non-convergence is.
- **Newton trajectory**: the IBL-floor signature intact — merit pinned at
  9.385e-09 (GV5.4: 9.35e-09), |F_BL| at 2.874e-6, λ ~ 1e-4, every step
  accepted (the committed accept-or-least-bad idiom). One benign
  RuntimeWarning (sqrt of a transient negative q² during the step-4
  backtracking; the least-bad line search rejects those tries).
- **W3 FD guard**: per-block median rel error φ 9.56e-12 / Γ 7.34e-12 /
  BL 0.00e+00 — PASS (the same 60-row sample as GV5.4).
- **Memory / wall @16t**: ru_maxrss 1.7 GB; seed 106 s + IBL 1000 s +
  W3 + Newton 120 s (the medium leg 1261 s total).
- **Coarse shakedown** (narrative only, not gated): 27,047 DOFs; W1
  recorded (cl_p 0.25579, the committed −5.35 % mesh-effect family); W3
  PASS (medians 5.09e-12 / 9.91e-12 / 1.53e-14); **rung 3 WORKS at
  coarse** — 5/5 steps GMRES info=0, 120–209 iters, rel_res ≤ 3.5e-9,
  no escalation (RECORDED; GV5.4's coarse rung 2 also worked); ratio
  18.98× (the 0.16 s inviscid step too small to read, the GV5.4 20.58×
  class).

## Reading

- **The pre-registered hypothesis is falsified in its naive form**:
  putting the per-node (quasi-simultaneous) sparsification Ĉ =
  J_hB·D_BB⁻¹·J_Bh of the Schur correction INTO the AMG matrix does not
  repair the reduced-space preconditioner at medium — it destroys it
  (rel_res 0.66 vs GV5.4 rung-2's 2e-7..6e-5 stagnation band). The
  coarse/medium split is the anatomy: the D_BB-local correction is
  adequate at coarse (rung 3 converges every step) but at medium either
  the inter-node BL coupling it drops becomes essential, or Ŝ_φφ = J_φφ
  − Ĉ_φφ loses the algebraic character (the diagonal-dominance /
  M-matrix-likeness) that made plain AMG(J_φφ) nearly work — Ĉ_φφ roughly
  doubles the block's density and subtracts a feedback term. The two
  discriminating measurements are committed (`results/gmres_*.csv`).
- **Rung 4 does not rescue it**: carrying BOTH coupling directions in the
  block-triangular preconditioner (with the same Schur-aware φ cycle)
  stagnates identically (0.68–1.06) — the failure is in the φ-cycle's
  quality on Ŝ_φφ, not in the block direction(s) carried.
- **The cost number is unchanged in kind**: 7.87× ≈ GV5.4's 7.53× (both
  dominated by capped-GMRES work; the per-step setup overhead of the
  Schur-aware machinery is +0.1 s — irrelevant). The ≤ ~2× band stays
  unreached; the honest negative the gate was registered to record.
- **What this closes and what it leaves**: the "AMG on a sparsified
  Schur" direction is now measured dead at wing scale (both block
  directions, D_BB-local physics). The escalation route that remains
  consistent with the diagnosis is an inter-node (A,Ψ)-structured M_BB
  inside Ĉ (the station-strip / ILU-based approximate inverse —
  pre-registration §7) or abandoning the corrected-AMG form for a
  fundamentally different reduced-space preconditioner; both stay
  registered-not-opened (user adjudication).

## Follow-ups (registered, NOT opened — user adjudication)

- An **inter-node (A,Ψ)-structured M_BB** (station-strip block inverse or
  an ILU-based approximate inverse inside Ĉ) — the recorded escalation
  route beyond rung 4.
- The **Eisenstat-Walker forcing** variant (η_k instead of the fixed rtol
  1e-8) — unchanged from GV5.4 §7; note it does not address a
  preconditioner-quality failure of this size (0.66 rel_res), only the
  near-miss stagnation class.
- Productizing into `tight_driver` — MOOT (no rung worked).
- The pressure-Kutta row wiring — unchanged from GV5.4 §7.

## Artifacts

`results/summary.csv`; per level: `inviscid_anchor_{level}.csv`,
`steps_{level}.csv`, `gmres_{level}.csv` (incl. t_lu / t_corr / t_amg /
nnz(Ĉ) / nnz(Ŝ) / n_fallback in `notes`), `fd_guard_{level}.csv`,
`cost_breakdown_{level}.png`; `run.log`.
