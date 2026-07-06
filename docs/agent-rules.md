# pyFP3D Agent Rules

Current phase: **P2 + M0 closed (2026-07-06)** — wake cut / circulation / Kutta delivered on the M0 mesh: `mesh/wake_cut.py` (per-node flood-fill side classification; **TE nodes ARE duplicated** — the original “TE not duplicated” topology assert was re-specced with evidence: a single-valued TE tapers [φ] to zero over the first wake cell ≡ a point vortex parked at the TE, spurious suction ~Γ²/h that *diverges* under refinement; see roadmap P2 assert block), `constraints/wake.py` (master–slave elimination, Γ RHS-only), `constraints/dirichlet.py` (incompressible vortex far-field, branch cut on the wake sheet), `solve/picard.py::solve_laplace_lifting` (secant-accelerated Kutta — 2 updates; plain ω-relaxation would need O(100) since the measured map slope is b≈0.93), `post/section_cut.py` general z=const marching-tets path + sectional wall-Cp curves, `post/surface.py` triangle-wise force integration, Hess–Smith panel reference `cases/reference_data/naca0012_incompressible/`. Gates G2.1–G2.5 all green (medium cl −0.82% vs panel; Γ-vs-pressure cl 0.01%; G2.5(b) re-specced to p99-of-|w| ≥ 1st-order decay, measured ratio 2.05, stripe-free heatmap); M0's acceptance link closed with it (topology asserts pass on the whole mesh family, hard rule 7 test). Suite: 100 passed + 2 xfailed, ~38 s. **Next priority: P3 (subsonic compressible: Picard density outer loop + retiring the P1 assembly tech debt — precomputed B_e/V_e, colored `prange` assembly per design.md §7).** Also: P1 still open on one item — **P1 gates renumbered in workflow order 2026-07-06** (old → new: G1.3 → G1.2, G1.2-a0 → G1.3, G1.2-a → G1.4, G1.2-c → DP1, G1.2-b → G1.5, G1.2 → G1.6; mapping table in roadmap.md P1 section): G1.1 MMS and G1.2 CG-iteration gates closed; **G1.3 + G1.4 oracle experiments completed 2026-07-06 with negative results** (`solve/wall_correction.py`, `tests/test_wall_correction_cylinder.py`, `cases/demo/g14_sphere_oracle_experiment.py`; on body-fitted meshes the boundary-data defect is (near-)zero — exact per-facet flux vanishes by the divergence theorem — Option A ceiling ≈ 11.3% vs < 2% target; the cylinder case was re-framed as recovery-dominated, not crime-dominated), so DP1 took the "> 5%" branch and G1.5 is void; G1.6 sphere-Cp gate stays open (strict xfail) — its **Option C acceptance re-spec (geometry-consistent reference) is the open P1 item**, curved/isoparametric wall elements are separately scoped (design.md §5.1.2); do not propose h-refinement, recovery tweaks, Nitsche, or further boundary-data corrections — all ruled out with evidence — per [docs/roadmap.md](roadmap.md).

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
