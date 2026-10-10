# FINDINGS — V5 IBL-floor diagnosis (the GV5.1 follow-up)

Binding text: `../PRE_REGISTRATION.md` (committed at 53bf904 BEFORE the
first execution). Verdict class: **RECORDED** (diagnostic; no pass/fail
bands). One execution of `run.py` (2026-07-24, wall 2928 s under a
machine load ~90) regenerated every artifact in `results/` from scratch;
the two loose-converged seeds reproduced bit-identically across three
independent HEAD regenerations (coarse cl=0.26791639 n_outer=4, medium
cl=0.28137437 n_outer=3; wiring guards dcl_k0=1.56e-12 / 1.31e-9,
converged=True). Operator: J_BL,BL = the U-U block of
`tight_driver.augmented_jacobian` at `pack.x_base()` (edge data frozen at
pack base; the GV5.1 FD-verified semantics); dense SVD, exact to machine
precision at 1236²/2460².

States: **S1** = coarse loose-converged, **S2** = medium loose-converged
(every S2 number carries the GV5.1 §4 trajectory-scatter caveat; the
diagnosis conclusions below rest only on features common to S1 and S2),
**S3** = coarse k=1 fixture (the Stage-3 smoke baseline).

## Q1 — spectrum: the near-null cluster PERSISTS at the converged state

| state | n | σmax | σmin | cond | <1e-6·σmax | <1e-8 | <1e-10 |
|---|---|---|---|---|---|---|---|
| S3 (k=1) | 1236 | 1.0 | 2.45e-11 | 4.08e10 | 501 | 121 | 1 |
| S1 (conv) | 1236 | 1.0 | 7.59e-12 | 1.32e11 | 500 | 136 | 1 |
| S2 (conv) | 2460 | 1.0 | 2.52e-14 | 3.97e13 | 1082 | 541 | 12 |

The cluster is **not k=1-specific**: S1 reproduces the S3 spectrum almost
curve-for-curve (`spectrum.png` — the s1/s3 lines overlap), and S2 carries
a proportionally *larger* cluster (44 % of all σ below 1e-6·σmax).
σmax = 1.0 exactly at all three states — the identity Dirichlet/laminar
pin rows (ibl3.py `_apply_rows`) set the spectral scale; the "cluster" is
measured against it (see Q3).

## Q2 — null-space anatomy: turbulent (A, Ψ), mid-chord → TE

Energy of the top-20 smallest-σ right singular vectors, by state variable
(δ, A, B, Ψ, Cτ1, Cτ2):

- S1: **A = 0.427, Ψ = 0.558** (B = 0.014; δ, Cτ1, Cτ2 ≤ 3e-9).
- S2: **A = 0.358, Ψ = 0.635** (B = 0.006; δ, Cτ1, Cτ2 ≤ 2e-9).

Common feature: the near-null space is carried by the turbulent
profile-shape/crossflow variables **(A, Ψ)** — not by the thickness δ,
not by the stress Cτ (Cτ is pinned at laminar nodes; at turbulent nodes
it stays well-conditioned). Nodal support (`nullspace_map_{s1,s2}.png`):
distributed over the **whole turbulent region** x/c ≳ 0.05 with a
mid-chord bias (region mass: mid 0.86/0.82, TE band 0.14/0.18, LE bands
0.000), machine-zero (≤1e-31) in the laminar/inflow-pinned LE band,
near-symmetric upper/lower. Top-10 nodes hold only ~0.49/0.44 of the
mass — a *distributed* flat direction set, not a point singularity at
the TE or the LE pin edge.

## Q3 — scaling: the 4e10–4e13 is MOSTLY a scaling artifact; a real (A,Ψ) stiffness remains

One pass of row-then-column 2-norm equilibration R·J·C:

| state | cond raw | cond scaled | row-norm range | col-norm range | <1e-6·σmax raw → scaled |
|---|---|---|---|---|---|
| S3 | 4.08e10 | 2.10e4 | 8.4e4 | 1.4e6 | 501 → 0 |
| S1 | 1.32e11 | 7.43e5 | 7.7e4 | 3.0e6 | 500 → 0 |
| S2 | 3.97e13 | 1.14e7 | 3.3e5 | 7.6e5 | 1082 → 2 |

Discrimination (the pre-registered solver-vs-formulation question): the
raw near-null cluster is a **scaling artifact** — after equilibration
essentially *no* singular value sits below 1e-6·σmax (0/0/2), i.e. there
are no exact null directions (consistent with Q4's empty floor-active
set). But the residual conditioning is still 2e4 → 7e5 → 1e7 and grows
from S3 to S1 to S2: a **genuine stiffness** lives in the physics (Q2's
(A,Ψ) block). Consequence: a solver-level fix (equilibration) removes
the spectral artifact; the remaining ~1e5–1e7 is the formulation-level
target for GV5.1b (regularization / projected solve) and GV5.4 (block
preconditioner).

## Q4 — closure-floor active set: EMPTY at S1 and S2

Mirroring `closure_node` (closures.py:490-507) at the pack-base edge
data: **0 of 206 (S1) / 0 of 410 (S2)** nodes have DELTA_MIN active
(min δ = 1.62e-6 / 1.78e-6 — 2.2 decades above the 1e-8 floor) and none
have RE_D_MIN active (min re_d = 2.43 / 2.31 — 3.4 decades above the
1e-3 floor). (The remaining min/max clamps in closures.py are 1e-30/1e-12
division guards plus the inner ctau-solve step clamp — never active at a
physical state and not state maps.) Two consequences:

- The "floor-active rows have zero interior derivatives → candidate
  exact null directions" hypothesis is **not operative** at the
  converged states — the closure floors explain neither the Q1/Q2 null
  space nor the Q5 residual floor.
- The DELTA_MIN sensitivity at S1/S2 is **identically zero** by
  construction — the substitute evidence for Q6(b).

## Q5 — residual anatomy: TE band, (B, δ) equations, fully inside J's range

|F_BL|inf = 3.154e-6 (S1) / 1.710e-6 (S2) — matching the committed
loose-final IBL floors (3.154e-6 / 1.712e-6, HEAD regen). Support:
**trailing-edge band** — S1 top nodes x/c 0.96–0.98 (upper) plus one
mid node x/c 0.124 (lower); S2 top nodes x/c 0.98–1.00 (both sides).
Equation norm shares (δ, A, B, Ψ, Cτ1, Cτ2): **B = 0.83, δ = 0.48**
dominant at both states (A, Ψ ≈ 0.19; Cτ ≈ 0.01) — the floor residual is
carried by the crossflow-shape (B) and mass/defect (δ) equations, *not*
by the (A,Ψ) equations that carry the null space (Q2). Alignment with
the left singular vectors (`residual_alignment.png`): |u_iᵀF|/|F| over
the 20 smallest-σ left vectors is max 7.7e-3 (S1) / 6.0e-3 (S2) and
mostly 1e-7–1e-15 — **10+ orders of magnitude below** the alignment with
the well-conditioned left vectors (1e-3–1e-1). The residual is therefore
essentially entirely **inside the well-conditioned range** of J_BL,BL:
the floor is not a range-deficiency. Combined with Q2: the Jacobian
"sees" the residual, but eliminating it requires motion along the flat
(A,Ψ) directions — where the nonlinear residual does not follow the
linear model (Q7).

## Q6 — formulation sensitivity probes (S1): NOT an artificial-viscosity floor; DELTA_MIN deferred

(a) eps_diff probes (fresh IBL3Solver at scaled constructor kwargs,
solved from the same seed U0 = S1's U; `sensitivity_probes.csv`):

| probe | floor | ×base | near-null (<1e-6) | cond |
|---|---|---|---|---|
| base (×1, the Q7 re-solve) | 3.154e-6 | 1.0 | 500 | 1.38e11 |
| eps_diff ×0.5 | 3.062e-6 | 0.971 | 501 | 1.93e11 |
| eps_diff ×2 | 3.330e-6 | 1.056 | 498 | 6.45e10 |
| both ×0.5 (extra) | 3.056e-6 | 0.969 | 503 | 5.10e10 |
| both ×2 (extra) | 3.894e-6 | 1.235 | 483 | 2.52e10 |

A 4× change of the artificial viscosity moves the floor by ≤ 6 %
(≤ 24 % even when both coefficients scale) and the near-null count by
< 4 % — the floor is **not** an artificial-viscosity truncation floor,
and the near-null cluster is not diffusion-tunable either.

(b) DELTA_MIN probe: **DEFERRED** per the pre-registration fallback
clause (no clean in-run probe exists; no source edited anywhere in
this study).
DELTA_MIN/RE_D_MIN are module constants baked into the
`@njit(cache=True)` closures (closures.py:77-78, consumed in
`closure_node` :490-507): exposed through neither the `closure_all`
signature nor the `IBL3Solver` constructor, and numba freezes globals at
compile time — the pre-registered "logged local edit reverted in the
same run" would silently no-op against the cached compilation, so no
clean in-run probe exists. Substitute evidence: Q4's empty active set
makes the DELTA_MIN sensitivity identically zero at S1/S2.

## Q7 — pseudo-time stall anatomy (S1): the controller bottoms out with the residual frozen

Standalone `IBL3Solver.solve` from S1's U (verbose controller trace,
`pseudotime_stall.csv`): the pure steady residual starts at 3.154e-6 and
**never moves** (3.154000e-6 at every one of the 21 iterations; the few
accepted steps change it below 1e-12 relatively — motion along the
near-null manifold); 14 step rejections; cfl collapses 1.0 → 1.25e-1 →
… → 1.0e-3 = cfl_min and pins there; exit on n_fail > 10 (n_iter = 21,
converged = False). Read on the pre-registered dichotomy: the **step
controller bottoms out** — and because even arbitrarily small
pseudo-time steps find no descending direction, this is the
*formulation floor expressed through the controller* (a local minimum of
the steady residual norm on the flat (A,Ψ) manifold), not a
step-size tuning issue. Feeds the GV5.1b design directly: globalization
changes alone cannot pass this floor.

## Synthesis (features common to S1 and S2 only)

1. The near-null cluster persists at the loose-converged states (it is
   not a k=1 artifact) and lives in the turbulent **(A, Ψ)** variables,
   distributed mid-chord → TE.
2. Its raw magnitude (cond 4e10–4e13, 500–1082 sub-1e-6 σ's) is
   **mostly a scaling artifact** (equilibration: 0–2 sub-1e-6 σ's); a
   genuine stiffness of cond ~1e5–1e7 remains after scaling — the real
   GV5.1b/GV5.4 target.
3. The IBL floor residual lives at the **TE band** in the **(B, δ)**
   equations and lies essentially entirely inside J's range.
4. Neither the closure floors (inactive by 2+ decades) nor the
   artificial viscosity (≤ 6 %/4× floor movement) explains the floor.
5. The pseudo-time controller bottoms out with the residual frozen at
   the floor — a formulation-level floor, not a solver limitation.

## Honesty notes

- Dense SVD at 1236²/2460² is exact to machine precision; no iterative
  eigensolver tolerance involved.
- All S2 numbers carry the GV5.1 §4 trajectory-scatter caveat (the
  floor itself appeared on every trajectory; state-dependent numbers are
  quoted per-state). The HEAD regeneration is bit-identical run-to-run
  (verified across three independent executions of this study).
- Q6(b) is deferred with the reason recorded in
  `sensitivity_probes.csv`; no source was edited anywhere in this study.
- Extra artifacts beyond the pre-registered list:
  `residual_alignment_{s1,s2}.csv` (the raw alignment series behind
  `residual_alignment.png`).
- Wall-clock is trivial science-wise (seed regen 317 s/1712 s under a
  machine load ~90; dense SVDs ≤ 150 s); no heavy committed artifact was
  recomputed. `run.py` regenerates every artifact from scratch in one
  pass (`--states`/`--phases` allow partial re-runs that merge
  `summary.csv` by (question, states)).
