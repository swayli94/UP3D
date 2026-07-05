# pyFP3D Agent Rules

Current phase: P0 (solver) + M0 (meshing), per [docs/roadmap.md](roadmap.md).

Authoritative docs:
- [docs/design.md](design.md): theory, discretization, architecture, and validation formulas.
- [docs/roadmap.md](roadmap.md): phases, gates, and progress ledger.

Hard rules:
1. After any kernel or assembly change, run `pytest cases/test_v0_freestream.py` first.
2. A phase closes only when its medium-mesh gate and the full coarse regression suite are green.
3. Numba kernels use SoA arrays only, no Python objects inside `@njit`, colored assembly for `prange`, zero allocation in hot loops, `cache=True`, and the `PYFP3D_NOJIT` switch.
4. All sparse linear algebra stays in SciPy/PyAMG.
5. Physics scalars live only in `physics/isentropic.py`.
6. Never edit files under `cases/reference_data/`.
7. Wake-cut changes must rerun topology asserts on all meshes in `cases/meshes/`.
8. Mesh generation stays vanilla Gmsh; node duplication happens only in the solver preprocessor.
9. Validate against full-potential references before Euler.
10. Keep commits phase-scoped and update the v2 progress ledger when a gate closes.
