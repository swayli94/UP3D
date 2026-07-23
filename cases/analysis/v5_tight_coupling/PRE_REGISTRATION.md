# GV5.1 PRE-REGISTRATION — augmented Newton: exactness + convergence

Binding text: `docs/roadmap/track_v.md` **GV5.1** (with the 2026-07-22
user-directed pre-registered FD note): both coupling blocks FD-verified
(project Jacobian discipline; the B19/B31 FD-gate pattern) + quadratic tail
on the GV3.1 case; outer iterations ≤ half the V3 loose loop. Committed
before the first execution (project discipline; the GV5.0 precedent).

## Architecture (pinned)

Monolithic Newton over the augmented state **x = (φ_free, Γ, U)**:

- φ_free: inviscid full-potential DOFs (reduced, `NewtonWorkspace.free`);
  Γ: wake-station circulations (P8/P14 machinery, `solve/newton.py`).
- U: IBL state on wall surface nodes, 6 unknowns per node
  (δ, A, B, Ψ, Cτ1, Cτ2), `viscous/ibl3.py` layout.

Residual F(x) = (F_φ, F_Γ, F_BL):

- F_φ/F_Γ: the P8/P14 residual with the transpiration RHS threaded through
  the V2 `external_rhs` channel: Tᵀ(R_inv(φ,Γ) − b(ṁ(φ,U))) with
  ṁ = P · div_Γ(ρ_e(φ) · u_e(φ) · δ\*(U)) (P = tip/tail pin mask,
  linear diagonal 0/1; `viscous/transpiration.py`).
- F_BL: the IBL3 residual R_IBL(U; e(φ)) with e the per-node edge data
  (q, ρ_e, μ, M_e, direction frames) recovered from φ by
  `edge_velocity_per_zone` + isentropic relations (the GV5.0 driver chain).

Jacobian blocks (all assembled sparse):

- J_φφ, J_φΓ, J_Γφ, J_ΓΓ: existing P8/P14 blocks
  (`assemble_coupled` + `kutta_blocks`), **plus the augmentation**
  −Tᵀ·W·∂ṁ/∂φ through ρ_e(q²(φ))·u_e(φ) — explicit algebra (no closure
  calculus), included so the φ-block stays exact under the transpiration
  coupling (quasi-Newton otherwise).
- **J_φ,BL** = −Tᵀ·W·P·L(ρ_e,u_e)·D_δ\*(U): W = wall-load matrix of
  `assemble_transpiration_rhs` (linear in ṁ), L = surface-divergence
  operator (linear in δ\*), D_δ\* = per-node ∂δ\*/∂U (1×6 row from the
  closure packet's `douts[:, OUT_DS1, :]` — already shipped).
  All pieces exist; the block is assembled, not re-derived.
- **J_BL,φ** = (∂R_IBL/∂e)·(∂e/∂φ): ∂e/∂φ is the FIXED sparse linear
  recovery map (per-node 3×3 local frames + q=|u_e| + isentropic scalars;
  zones are state-independent — see the FD-note decision below);
  ∂R_IBL/∂e requires extending the closure packet with derivatives
  w.r.t. the edge inputs (q, ρ, μ, M) — the V1 packet explicitly does not
  produce them (`closures.py` "held fixed; derivatives w.r.t. them are NOT
  produced in V1") — and chaining them through the `_assemble` wiring
  wherever edge data appears (closure fluxes, the q²/ψ/u/w/ρ source terms,
  the s1 streamwise-diffusion direction, the diffusion scales).
  **This closure edge-derivative extension is the heavy piece of V5.**
- J_BL,BL: existing `IBL3Solver.residual_jacobian` (FD-verified in V1).

## Design decisions recorded per the gate's "the choice is recorded in the gate"

1. **FD-note decision (zone switch)**: the LE-band/elsewhere zone partition
   is decided by FIXED geometry (x/c < le_band_x at the node), NOT by the
   solution — so the φ → u_e recovery is a fixed sparse linear map with no
   kink; perturbing φ never moves a node across zones. Decision: **keep the
   discrete zone switch** (production behavior; smoothing the switch would
   change the physics to ease verification) and verify with **central FD
   plus kink-row masking** (the P8/B15 pattern: mask only rows whose
   perturbation genuinely crosses a non-smooth point; the sole candidate is
   q = |u_e| → 0 at the stagnation node). One-sided differences at
   boundary rows are retained as the fallback, pre-registered by the note,
   and fire only if masked rows exceed 2 % of the block.
2. **Tight-residual housekeeping**: the loose loop's δ\* under-relaxation
   (ω) is DROPPED (identity at convergence — the Newton map is the
   unrelaxed one); the max(δ̂\*, 0) floor is DROPPED (the closure packet
   already floors δ at DELTA_MIN internally with the derivative masked,
   which is FD-consistent by construction; the loose floor only guarded
   the relaxed fixed-point field); the tip/tail ṁ pin mask is KEPT
   (linear diagonal; it is a boundary condition of the coupled system).
3. **Linear algebra scope**: assembled sparse Jacobian + `splu` per Newton
   step at the 2.5-D GV3.1 scale (~30k DOFs medium). GMRES + block
   preconditioning (AMG-φ / ILU-BL) is DEFERRED to GV5.4 (the cost gate)
   per the deliverable's "measure, don't assume" — GV5.1 gates exactness
   and convergence, not iteration cost.
4. **Regime locks**: same as the loose loop (forced transition
   x_tr/c = 0.05 both sides, inflow Dirichlet band frozen at the seed's
   laminar/turbulent per-node regime, laminar stress pins,
   eps_diff = 0.005 / eps_diff_s = 0.02 — the load-bearing V1
   calibration).

## FD verification protocol (the B19/B31 pattern)

- Full-state random-direction central FD of the ASSEMBLED augmented
  Jacobian: v over (φ_free, Γ, U), max-normalized; fd =
  (F(x+εv)−F(x−εv))/(2ε) vs J@v; max-norm relative error, ε ladder
  (1e-5 / 1e-6 / 1e-7) with the sweet-spot recorded. Tolerance **< 1e-5**
  (the V1 IBL FD-gate discipline); unmasked rows are expected at
  1e-6–1e-8 (the P8/P14 precedent).
- Targeted per-block checks: J_φ,BL columns from unit U-perturbations;
  J_BL,φ columns from unit φ-perturbations (including columns whose
  recovery stencil touches the LE-band boundary — the decision-1 masking
  applies here); the J_φφ augmentation term separately against FD of
  F_φ w.r.t. φ at fixed U.
- State point: the loose-loop GV3.1 k=1 state (inviscid-converged φ + one
  IBL solve), coarse mesh (test-speed), repeated at the coupled solution.

## Convergence protocol (GV3.1 case)

- Case: the committed GV3.1 configuration — NACA0012 2.5-D wake-cut strip,
  M∞ = 0.5, α = 2°, Re 3.0e6, x_tr 0.05/0.05; **medium binding, coarse
  recorded** (GV3.1 discipline). Same inviscid driver recipe as GV3.1's
  Newton runs (P14 pressure-Kutta; the 2.5-D strip has no tip_taper).
- Seed: inviscid-converged (φ, Γ) + U = one IBL Newton solve on the
  inviscid u_e (identical to the loose loop's k=1 state before the FP
  re-solve — the loose FP re-solve is a fixed-point artifact the
  augmented system does not need).
- Globalization: backtracking line search on the augmented merit
  |F_φ|²+|F_Γ|²+|F_BL|² (the P8/P14 line-search idiom).
- Convergence criterion: max-norm scaled residual < 1e-8, with the
  δ\*-change of the last step recorded against tol_ds = 1e-3 as a
  cross-check.
- **Pass bands**: (a) both blocks FD-verified per the protocol above;
  (b) quadratic tail = a slope-2 regime in the residual-norm history over
  the last decades before the floor; (c) Newton iterations N_aug ≤ 2
  (= half the committed loose-loop count 5 on medium, floored).
  **Honesty note (pre-registered):** if N_aug lands at 3 the run is
  RECORDED and the letter-vs-spirit call goes to the user (the V4-skip
  adjudication precedent) — 3 < 5 still beats the loose loop but misses
  the halving letter; the FD exactness and the quadratic tail are this
  gate's teeth.

## Risks (recorded up front)

- The GV3.1 record shows the IBL inner solve stalling at a ~1.2e-6
  residual floor on the harsh k=1/k=2 states (max_iter=100). The augmented
  seed is exactly the k=1 state; if the floor is a genuine basin effect
  the first Newton steps carry it — recorded, not re-diagnosed here.
- `splu` memory at medium scale (~30k DOFs): expected fine at 2.5-D
  density; if it exceeds the machine, coarse becomes binding and medium is
  recorded with GMRES(fallback) — flagged in the VERDICT.
- The closure edge-derivative extension (q, ρ, μ, M chains through the
  eta-quadrature packet) is hand-derived calculus; the FD gate arbitrates
  (that is its purpose).

## Outputs (committed under `cases/analysis/v5_tight_coupling/results/`)

1. `gv5_1_fd_report.csv` — per-block and full-state FD errors, ε ladder,
   masked-row counts.
2. `gv5_1_newton_history_{coarse,medium}.csv` — per-iteration residual
   norms (F_φ/F_Γ/F_BL and total), step length, δ\*-change.
3. `gv5_1_compare.csv` — augmented vs loose: iteration counts, final δ\*
   max-diff, cl_p/cl_kj of both endpoints (the loose endpoint is the
   committed GV3.1 result, NOT recomputed).
4. `summary.csv` — pass/fail rows for the three pass bands + recorded
   rows; `VERDICT.md` after execution.
