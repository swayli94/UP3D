# pyFP3D Agent Rules

Current phase: **P3 closed (2026-07-07)** (P2 + M0 closed 2026-07-06; details in the roadmap ledger) — subsonic compressible delivered: NESTED density Picard `solve/picard.py::solve_subsonic{,_lifting}` (outer ρ update, inner P2 secant Kutta at frozen ρ — interleaving one Γ step per density step was measured to spike the residual 10× and was rejected; a loose *relative* inner CG tolerance with warm starts false-converges — the opt-in inexact scheme is a forcing term `atol = η‖b−Ax₀‖`, default η = 0; see roadmap G3.2), PG-scaled vortex far field (`constraints/dirichlet.py`, β stretches only the atan2 argument, β = 1 reduces bit-exactly), compressible Cp in post (`m_inf` params), and the **P1 assembly tech debt retired**: precomputed B_e/V_e + symbolic-pattern `elem_to_csr` + colored-`prange` kernels (`kernels/{gradient,jacobian,residual}.py`, `PicardOperator` per-mesh workspace; numba-jitted greedy coloring, same assignment as the old Python loop); the public `assemble_stiffness_matrix` now IS the fast path (5.7e-16 vs the retained serial reference kernels, bit-deterministic across calls/threads, ~160× hot reassembly on medium NACA). Gates: G3.1 0.32% vs PG-corrected (same-mesh comparison cancels the G1.6 wall bias); G3.2 cl −0.33% from the PG/Kármán–Tsien midpoint (`cases/reference_data/naca0012_m05/`), 15 density iterations, strictly monotone residual; G3.3 bitwise M∞ → 0 ⇔ Laplace for matrix, φ and Γ — enabled by seeding pyamg's spectral-radius RNG in `solve/linear.py::build_amg_preconditioner` (repeated identical solves differed at 2e-11 before; ALL solves are now bit-reproducible run-to-run). Suite: 117 passed + 2 xfailed, ~96 s (G3.2's medium nested solve ≈ 45 s of it). **Next priority: P4 (transonic artificial density: `kernels/upwind.py` upstream-element search + ν switch + ρ̃ per design.md §3, Mach continuation in `solve/continuation.py`, shock monitors; G4.2's ν ≡ 0 subcritical no-op must stay bit-identical to P3) and M1 meshing (needed by P5).** Also: P1 still open on one item — G1.6 sphere-Cp (strict xfail) awaits its **Option C acceptance re-spec** (geometry-consistent reference); curved/isoparametric wall elements separately scoped (design.md §5.1.2); G1.3/G1.4 oracles closed 2026-07-06 with negative results (Option A ceiling ≈ 11.3%), DP1 took the "> 5%" branch, G1.5 void; P1 gates renumbered 2026-07-06 (mapping in roadmap.md). Do not propose h-refinement, recovery tweaks, Nitsche, or boundary-data corrections for G1.6 — all ruled out with evidence — per [docs/roadmap.md](roadmap.md).

Authoritative docs:
- [docs/design.md](design.md): theory, discretization, architecture, and validation formulas.
- [docs/roadmap.md](roadmap.md): phases, gates, and progress ledger.

Hard rules:
1. After any kernel or assembly change, run `pytest tests/test_v0_freestream.py` first.
2. A phase closes only when its medium-mesh gate and the full coarse regression suite are green.
3. Numba kernels use SoA arrays only, no Python objects inside `@njit`, colored assembly for `prange`, zero allocation in hot loops, `cache=True`, and the `PYFP3D_NOJIT` switch.
4. All sparse linear algebra stays in SciPy/PyAMG.
5. Physics scalars live only in `physics/isentropic.py`.
6. Never edit files under `cases/reference_data/`.
7. Wake-cut changes must rerun topology asserts on all meshes in `cases/meshes/`.
8. Mesh generation stays vanilla Gmsh; node duplication happens only in the solver preprocessor.
9. Validate against full-potential references before Euler.
10. Keep commits phase-scoped and update the v2 progress ledger when a gate closes.
