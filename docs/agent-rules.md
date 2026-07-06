# pyFP3D Agent Rules

Current phase: **M0 (meshing) is the active priority** — re-specified 2026-07-06 as a single-layer extruded NACA0012 quasi-2D testbed (one cell layer in z, globally consistent prism→3-tet split, embedded wake; see roadmap.md Track M and gate G2.5). Also: P1 in progress (G1.1 MMS and G1.3 CG-iteration gates closed; G1.2 sphere-Cp gate open — root-caused to a curved-wall/flat-facet geometric consistency error in the volume solve, not the surface recovery scheme; see PROJECT_STRUCTURE.md "Known gaps" before proposing more h-refinement or recovery-scheme tweaks), P0 in progress (mesh infra done, M0 mesh family still missing) — per [docs/roadmap.md](roadmap.md).

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
