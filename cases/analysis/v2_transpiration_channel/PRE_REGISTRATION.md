# GV2.1 Pre-registration (written before any gate execution)

Binding gate text: `docs/roadmap/track_v.md` GV2.1(a)–(c) (2026-07-22 re-spec,
V2 — Transpiration channel through all three drivers). User directive
2026-07-22 (Track V opening): **no gate re-spec** — the bands below are the
execution-level operationalization the gate text itself requires, not a
redefinition of the gate.

The pytest suite (`tests/test_v2_transpiration.py`,
`tests/test_v2_newton_rhs_channel.py`) covers the same clauses as fast
regression locks; THIS gate is the committed-evidence execution
(CSV/PNG artifacts + VERDICT) per the project evidence discipline.

## Gate (a) — manufactured blowing on the M0 cylinder

- Case: `cases/meshes/cylinder_2.5d/{coarse,medium,fine}.msh` (committed;
  a = 1, r_far = 20, single-layer 2.5-D). Fourier blowing mode n = 2,
  v0 = 0.1 (linear problem — amplitude arbitrary): nodal manufactured
  `m_dot = v0 cos(2θ)` on the wall, assembled by
  `viscous/transpiration.py::assemble_transpiration_rhs` (b = −load(ṁ),
  blowing-positive), solved by `solve_laplace(body_source_rhs=...)`
  (rtol = 1e-11, maxiter = 3000), far-field Dirichlet = the exact
  analytic value.
- Analytic solution (exterior Laplace with ∂φ/∂r = v0 cos(nθ) at r = a,
  decaying): φ = −(v0·a/n)·(a/r)ⁿ·cos(nθ)
  (`tests/mesh_utils.py::cylinder_blowing_phi_exact`).
- **Error measure (binding)**: relmax = max|φ − φ_exact| / max|φ_exact|
  over ALL volume nodes (the max sits at the wall — the curved-boundary
  P1 approximation dominates there).
- **Bands (binding)**:
  1. relmax strictly decreases coarse → medium → fine (h_wall halves per
     level: 0.10 / 0.05 / 0.025);
  2. measured order p = log2(e_h/e_{h/2}) ≥ 1.0 on BOTH successive pairs
     (gate text: "φ error O(h) vs analytic").
- **RECORDED (no band)**: relmax per level, wall-only relmax, relL2 per
  level, CG iteration counts, final residual norms. The sign convention
  check is implicit: a flipped transpiration sign lands relmax at O(2),
  two orders outside any passing band — the gate fails loudly, by design.

## Gate (b) — ṁ = 0 bit-identity on ALL drivers

A/B each driver: default (channel absent) vs an explicit ZERO assembled
transpiration RHS. **Assertion (binding)**: `np.array_equal` on the
solution (phi / phi_ext, plus gamma where the driver carries it) AND on
the recorded residual history. Legs:

1. `solve_laplace` — cylinder coarse, base-flow far field
   (tests/mesh_utils.py::cylinder_phi_exact Dirichlet), rtol 1e-11.
2. `solve_subsonic` — cylinder coarse, M = 0.3, freestream far field.
3. `solve_subsonic_lifting` — NACA0012 2.5-D coarse (wake-cut), M = 0.3,
   α = 2°.
4. conforming Newton `solve_newton_lifting` — NACA0012 2.5-D coarse,
   M = 0.5, α = 2° (subcritical ⇒ walk selection inert),
   external_rhs = zeros.
5. LS Newton `solve_multivalued_newton` — NACA0012 2.5-D coarse, M = 0.3,
   α = 2°, farfield = "neumann", n_seed = 10, wall_rhs = zeros.

## Gate (c) — conforming Newton Jacobian EXACT under lagged ṁ

- Case: NACA0012 2.5-D coarse, M = 0.5, α = 2°, at the converged baseline
  state (subcritical ⇒ FD-safe, the test_p8 protocol). Nonzero lagged
  b_ext = the assembled transpiration load of uniform blowing
  ṁ = 0.01 at the wall nodes.
- **Assertions (binding)**:
  1. structural: the coupled Jacobian assembled at the same state is
     BIT-IDENTICAL with and without b_ext (`np.array_equal` on J_ff and B
     data — b_ext never enters assemble_coupled);
  2. FD: central-difference JVPs of eval_residual (which INCLUDES the
     b_ext subtraction) match the analytic blocks — J_ff on 3 random
     directions and B on every Gamma column — rel err < 1e-5 (the
     test_p8/test_b31 FD discipline; measured values quoted, expected
     ~1e-6–1e-9);
  3. algebra: R′_free = R_free − (Tᵀ b_ext)[free] at the baseline state
     to max abs err < 1e-12.

## Numerical settings (recorded with the artifacts)

- Thread cap 16: NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
  OPENBLAS_NUM_THREADS=16.
- Newton leg: upwind_c = 1.5, m_crit = 0.95, m_cap = 3.0, rho_floor = 0.05,
  default precond/tolerances (tol_residual = 1e-10).
- LS leg: n_seed = 10, n_newton_max = 15, default precond (spsolve).
- Nonzero-load consistency runs (RECORDED, not banded): cylinder coarse
  M = 0.3 with ṁ = 0.02·cos(2θ); NACA Newton with the (c) load; LS with
  uniform ṁ = 0.005.
