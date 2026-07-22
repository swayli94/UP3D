# GV2.1 VERDICT — transpiration channel exactness

- Date: 2026-07-22 · Branch: `kimi/track-v2-transpiration`
- Binding text: `docs/roadmap/track_v.md` GV2.1(a)–(c) (2026-07-22 re-spec)
- Pre-registration: [PRE_REGISTRATION.md](PRE_REGISTRATION.md) (committed
  before any gate execution; **no band changed during/after execution**)
- Evidence: `results/summary.csv`, `results/gv2_1a_refinement.csv`,
  `results/gv2_1c_fd.csv`, `results/gv2_1_panels.png` — regenerate with
  `python cases/analysis/v2_transpiration_channel/run.py` (exit 1 on any
  honest FAIL)
- **VERDICT: GV2.1 PASS — 23 PASS / 0 FAIL / 16 RECORDED.** V2 closable.

## Result table

| gate | metric (band) | measured | verdict |
|---|---|---|---|
| (a) | relmax strict decrease c>m>f (binding) | 2.1572e-02 > 6.8738e-03 > 2.2062e-03 | PASS |
| (a) | order coarse→medium (≥ 1.0) | 1.650 | PASS |
| (a) | order medium→fine (≥ 1.0) | 1.640 | PASS |
| (a) | relmax / relmax_wall / relL2 per level, n_cg, residuals | see gv2_1a_refinement.csv | RECORDED |
| (b) | solve_laplace zeros: phi + residual bit-identical | True | PASS |
| (b) | solve_subsonic zeros: phi + residual_history bit-identical | True | PASS |
| (b) | solve_subsonic_lifting zeros: phi + gamma + residual_history bit-identical | True | PASS |
| (b) | newton_lifting external_rhs=zeros: phi + gamma + residual_history bit-identical | True | PASS |
| (b) | ls_newton wall_rhs=zeros: phi_ext + residual_history bit-identical | True | PASS |
| (c) | J_ff / B bit-invariant under lagged b_ext | True / True | PASS |
| (c) | residual identity R′ = R − (Tᵀb)[free] (< 1e-12) | 0.0 exact | PASS |
| (c) | F invariant under b_ext | True | PASS |
| (c) | FD J_ff ×3 dirs + B gamma col (< 1e-5) | 6.6e-09 … 9.8e-09 / 7.2e-08 | PASS |
| (c) | nonzero load consistent solve | converged, ‖R‖ 2.69e-12, ‖ΔΓ‖ 2.98e-04 | RECORDED |

## Per-gate analysis

### (a) manufactured blowing — sign + accuracy pinned

The Fourier mode n = 2, v0 = 0.1 manufactured blowing assembles through
`viscous/transpiration.py::assemble_transpiration_rhs` (b = −load(ṁ),
blowing-positive) into `solve_laplace(body_source_rhs=...)` and reproduces
the analytic exterior Laplace solution φ = −(v0a/n)(a/r)ⁿcos nθ with the
max-norm error at the wall (the curved-boundary P1 approximation
dominates), strictly decreasing with order ≈ 1.65 — above the
pre-registered O(h) floor, below the O(h²) asymptote, consistent with the
G1.6/M0-cylinder family's curved-wall behavior. **The sign convention is
pinned by this clause**: a flipped sign lands relmax at O(2), two orders
outside the passing band — no silent sign error is possible.

### (b) ṁ = 0 bit-identity — all five driver legs

Every driver path that carries (or newly threads) the channel reproduces
its channel-absent result bit-for-bit under an explicit zero RHS:
`solve_laplace` (pre-existing hook), `solve_subsonic` and
`solve_subsonic_lifting` (V2 threading; the lifting path rides the
reduced_rhs Tᵀ reduction), `solve_newton_lifting`
(`NewtonWorkspace.external_rhs`), `solve_multivalued_newton` (existing
`b_base` slot). This is the pre-2026-07-22 sketch's "δ* = 0 bit-identical"
clause, re-homed to V2 by the re-phase.

### (c) Newton Jacobian exact under lagged ṁ

With a nonzero lagged b_ext (uniform blowing 0.01 at the wall) the coupled
Jacobian at the same state is **bit-identical** with and without the load
(b_ext enters only eval_residual's RHS, never assemble_coupled), the
residual identity R′_free = R_free − (Tᵀb_ext)[free] holds to 0.0 exact,
and the FD oracle of the b_ext-including residual still matches the
analytic J_ff / B to 6.6e-09–7.2e-08 (pre-registered band 1e-5; the
test_p8/test_b31 FD protocol). Tight coupling (∂ṁ/∂BL terms) is V5's
augmentation, deliberately NOT this channel.

## Implementation fixes made during gate execution

None. First execution passed every pre-registered band; no code, band, or
case was touched between pre-registration and this verdict.

## Numerical settings (as run)

- Thread cap 16 (NUMBA/OMP/OPENBLAS). Meshes: committed
  `cases/meshes/cylinder_2.5d/{coarse,medium,fine}.msh` (h_wall
  0.10/0.05/0.025) and `naca0012_2.5d/coarse.msh`.
- (a): solve_laplace rtol 1e-11, maxiter 3000; far-field Dirichlet =
  analytic.
- (b) cases: cylinder coarse base flow / cylinder M 0.3 / NACA cut M 0.3
  α 2° / NACA Newton M 0.5 α 2° / NACA LS M 0.3 α 2° neumann n_seed 10.
- (c): b_ext = uniform blowing ṁ = 0.01 at wall nodes; FD central
  difference eps = 1e-5 at the converged M 0.5 baseline (subcritical ⇒
  walk selection inert, FD-safe).
