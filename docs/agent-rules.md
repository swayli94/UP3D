# pyFP3D Agent Rules

Current phase: **M0 mesh-side items delivered (2026-07-06)** — `pyfp3d/meshgen/` (single-layer extrusion, globally consistent prism→3-tet split, quad-split consistency assert), NACA0012 quasi-2D family (`cases/meshes/naca0012_2.5d/`) plus a cylinder-flow validation case (`cases/meshes/cylinder_2.5d/`, analytic Cp), tests in `tests/test_m0_*.py`. M0 closure waits on P2 (`mesh/wake_cut.py` topology asserts + gate G2.5, whose criterion (b) needs re-spec: solved-field spanwise noise is O(h) by construction for 3-tet prisms — see roadmap G2.5 evidence note). **Next priority: P2 (wake cut, circulation, Kutta) on the M0 mesh.** Also: P1 in progress (G1.1 MMS and G1.3 CG-iteration gates closed; G1.2 sphere-Cp gate open — root-caused to a curved-wall/flat-facet geometric consistency error in the volume solve, not the surface recovery scheme; fix routes now researched and documented (design.md §5.1, roadmap G1.2-a/b/c) — the next G1.2 action is the G1.2-a0 cylinder oracle pre-study (see roadmap G1.2-a0; docs only, no code yet); see PROJECT_STRUCTURE.md "Known gaps" before proposing more h-refinement or recovery-scheme tweaks) — per [docs/roadmap.md](roadmap.md).

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
