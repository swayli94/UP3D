# Project Structure

This document describes the directory layout and initialization status of pyFP3D.

## Directory Tree

```
pyfp3d/                    # Main package
├── __init__.py           # Package entry point, NOJIT mode flag
├── mesh/                 # Mesh I/O, topology, metrics, coloring
│   ├── __init__.py
│   ├── reader.py         # [P0] meshio → SoA arrays + boundary tags
│   ├── metrics.py        # [P0] volumes, gradients, face adjacency
│   ├── coloring.py       # [P0] element graph coloring
│   └── wake_cut.py       # [P2] node duplication from wake surface
├── physics/              # Physics constants and constitutive relations
│   ├── __init__.py
│   └── isentropic.py     # ✓ [P0] ρ(q²), M(q²), a(q²), Cp, etc. (complete)
├── kernels/              # Element-wise assembly kernels (Numba-jitted)
│   ├── __init__.py
│   ├── residual.py       # [P1] Laplace residual + stiffness assembly (done) → [P3] isentropic → [P4] Newton
│   └── upwind.py         # [P3] Artificial compressibility
├── solve/                # Linear and nonlinear solvers
│   ├── __init__.py
│   ├── linear.py         # [P1] Dirichlet elimination + CG/PyAMG preconditioner (done)
│   ├── picard.py         # [P1] Laplace driver (single linear solve, no outer loop) (done)
│   └── newton.py         # [P4] Newton method
└── post/                 # Post-processing
    ├── __init__.py
    ├── vtk_out.py        # [P0] Write .vtu for ParaView; also the PNG/CSV gate-artifact helpers
    │                       #      (export_error_heatmap, export_matplotlib_plot) live here, not
    │                       #      in a separate artifacts.py
    └── surface.py        # [P1] Nodal gradient recovery (volume-weighted), for surface Cp (done,
                            #      but see "Known gaps" -- boundary accuracy not yet gate-passing)

cases/                     # Test cases and reference data
├── meshes/               # Mesh families (coarse/medium/fine)
│   ├── sphere_shell/     # [P1] Gmsh sphere-shell case for gate G1.2 (coarse/medium generated)
│   ├── naca0012_2.5d/    # [M0] Extruded NACA0012 + wake surface -- not started
│   └── onera_m6/         # [M1] Swept wing -- not started
├── reference_data/       # Ground truth (DO NOT EDIT)
│   ├── sphere_incomp_cp.csv
│   ├── naca0012_cl.csv
│   └── ...
└── test_*.py             # [Deprecated] Integration tests (use tests/ now)

tests/                     # Unit and gate tests
├── conftest.py           # ✓ Pytest fixtures: artifacts_dir, mesh_dir, etc.
├── mesh_utils.py         # ✓ [P1] Dependency-free structured-cube + sphere-shell mesh generators
├── __init__.py
├── test_v0_freestream.py # ✓ [P0/P1] Primary regression test (incl. cut-free residual check)
├── test_mesh_*.py        # [P0] Gates G0.1–G0.4
├── test_laplace_*.py     # [P1] Gates G1.1–G1.3 -- NOT YET WRITTEN (see "Known gaps")
├── test_wake_*.py        # [P2] Gates G2.1–G2.3
└── test_transonic_*.py   # [P3] Gates G3.1–G3.2

artifacts/                 # Gate outputs (auto-generated, gitignored)
├── G0.1/                 # Volume conservation heatmap
├── G0.2/                 # Gradient recovery plots
├── G0.3/                 # Element coloring 3D render
└── ...

pyproject.toml            # ✓ Project metadata and dependencies
setup.py                  # ✓ Legacy setup (pyproject.toml preferred)
.copilot-instructions.md  # ✓ This AI agent's domain-specific rules
```

## Implementation Status

### ✓ Complete (P0)
- **pyfp3d/physics/isentropic.py** — All physics scalars, numba-jitted
- **pyfp3d/mesh/reader.py** — Mesh I/O (meshio → SoA), mesh validation, tagged round-trip write
- **pyfp3d/mesh/metrics.py** — Geometry: volumes, gradients, face adjacency (Numba)
- **pyfp3d/mesh/coloring.py** — Element graph coloring for @prange
- **pyfp3d/post/vtk_out.py** — VTK writer (point fields) + PNG/CSV artifact helpers
- **tests/conftest.py** — Pytest fixtures
- **tests/test_v0_freestream.py** — Smoke tests + regression baseline ✓
- **tests/test_mesh_volume.py** — Gate G0.1 (volume conservation) ✓
- **tests/test_mesh_gradient.py** — Gate G0.2 (gradient recovery) ✓
- **tests/test_mesh_coloring.py** — Gate G0.3 (element coloring) ✓
- **tests/test_io_vtk.py** — Gate G0.4 (VTK round-trip) ✓
- **pyproject.toml** — Build metadata and dependencies
- **.copilot-instructions.md** — AI agent instructions for P0–P5

26/26 tests pass (`pytest tests/`). Three latent bugs found by manual code audit (not caught by
the existing suite, because nothing exercised these code paths) have been fixed:
- `mesh/metrics.py::build_face_adjacency` crashed under `@njit` (reflected-list dict values are
  not valid in numba nopython mode) — rewritten around a `numba.typed.Dict` keyed by sorted face,
  storing only the packed first-owner index instead of a growing list.
- `mesh/reader.py::write_mesh` silently dropped every named boundary group (`wall`, `farfield`,
  ...), writing only a legacy `"all_triangles"` block, and `.msh` was ambiguous between meshio's
  `ansys`/`gmsh` writers (defaulted to `ansys`, discarding all tag data). Now writes each boundary
  group as its own tagged `triangle` block plus `gmsh:physical`/`field_data`, explicitly via the
  `gmsh22` writer.
- `solve/picard.py::solve_laplace` reported `residual_norm` over *all* nodes, including Dirichlet
  (far-field) rows whose natural-BC flux imbalance is O(1) and never shrinks — swamping the actual
  free-dof residual (which was already converging to ~1e-10). Now restricted to free dofs and
  correctly nets out `body_source_rhs` when present (needed for the future MMS gate G1.1).

### 🔧 In progress, uncommitted to a gate yet (P1)
- **pyfp3d/kernels/residual.py** — Laplace residual (6.1) + SPD stiffness matrix (6.2) assembly
- **pyfp3d/solve/linear.py** — Dirichlet elimination (principal-submatrix reduction) + CG+PyAMG
- **pyfp3d/solve/picard.py** — `solve_laplace()` driver (P1's Picard loop degenerates to one solve)
- **pyfp3d/post/surface.py** — `nodal_gradient_recovery()`, volume-weighted average for surface Cp
- **tests/mesh_utils.py**, **cases/meshes/sphere_shell/** — structured-cube (MMS) and sphere-shell
  (G1.2) mesh generators; sphere-shell coarse/medium `.msh` + inspection PNGs are committed

**Known gaps before P1 can close:**
- No `tests/test_laplace_*.py` exists yet — G1.1 (MMS convergence), G1.2 (sphere Cp), and G1.3
  (mesh-independent CG iteration count) are all unimplemented as gates.
- Manually running the G1.2 setup against the committed sphere-shell meshes (uniform flow, exact
  Dirichlet BC at far field, natural BC at the wall) shows the φ field itself converges reasonably
  (~1% error), but surface Cp from `nodal_gradient_recovery` is far off the 2% target: ~26% max /
  9% mean error on the medium mesh (roughly halving coarse→medium, consistent with first-order
  boundary-recovery error). This is a numerical-accuracy gap in the boundary gradient recovery,
  not a logic bug — likely needs a better boundary evaluation (e.g. one-sided extrapolation or
  patch recovery) before G1.2 can be written and pass.

### ⏳ Next
- Write `tests/test_laplace_mms.py`, `test_laplace_sphere.py`, `test_laplace_cg_iterations.py`
  for G1.1–G1.3, after resolving the surface-Cp accuracy gap above.
- Generate NACA0012 mesh family (M0 gate) — not started.

## Quick Start

### 1. Install dependencies
```bash
cd /home/lrz/code/UP3D
pip install -e ".[dev]"
```

### 2. Run smoke tests
```bash
pytest tests/test_v0_freestream.py -xvs
```

Expected output:
```
tests/test_v0_freestream.py::test_import_pyfp3d PASSED
tests/test_v0_freestream.py::test_import_physics PASSED
tests/test_v0_freestream.py::test_isentropic_stagnation PASSED
tests/test_v0_freestream.py::test_isentropic_freestream PASSED
tests/test_v0_freestream.py::test_pressure_coefficient_bounds PASSED
```

### 3. Run physics module directly (self-test)
```bash
python -m pyfp3d.physics.isentropic
```

Expected output:
```
=== Isentropic Physics Self-Test ===

Freestream Mach M∞ = 0.5

At stagnation (q² = 0):
  ρ = 1.000000 (expected 1.0)
  Cp = 1.000000 (expected 1.0)

At freestream (q² = 1):
  ρ = 1.000000
  M = 0.500000 (expected 0.5)
  Cp = 0.000000 (expected 0.0)

Critical speed q*² = 0.923077 where M = 1.0:
  Computed M = 1.000000 (expected 1.0)

✓ All checks passed!
```

## Next Steps (P1 Completion)

P0 (mesh I/O, metrics, coloring, VTK writer, gates G0.1–G0.4) is done. To close P1:

1. **Resolve the surface-Cp accuracy gap** (see "Known gaps" above) in
   `post/surface.py` before G1.2 can be written meaningfully.

2. **Write `tests/test_laplace_mms.py`** (gate G1.1)
   - Manufactured φ_exact with a consistent FEM load vector for `body_source_rhs`
   - L² error vs h over 3 mesh levels (structured cube from `tests/mesh_utils.py`), slope ≥ 1.9

3. **Write `tests/test_laplace_sphere.py`** (gate G1.2)
   - Uniform flow past the sphere-shell mesh, Dirichlet at far field, natural at wall
   - Max |Cp − (1 − 9/4 sin²θ)| < 2% on the medium mesh

4. **Write `tests/test_laplace_cg_iterations.py`** (gate G1.3)
   - Run `solve_cg_amg` on coarse/medium/fine Laplace problems, iteration counts within 5%

5. **Create test meshes** (M0, parallel track)
   - Gmsh script: extruded NACA0012 + embedded wake surface
   - Generate coarse (30k), medium (150k), fine (700k) families
   - Validate topology with asserts

## References

- **Design & Theory:** [docs/design.md](../docs/design.md)
- **Roadmap & Gates:** [docs/roadmap.md](../docs/roadmap.md)
- **Agent Rules:** [docs/agent-rules.md](../docs/agent-rules.md)
- **These Instructions:** [.copilot-instructions.md](../.copilot-instructions.md)

---

**Last updated:** 2026-07-06  
**Status:** P0 closed (G0.1–G0.4 green, 26/26 tests). P1 in progress: Laplace residual/stiffness
assembly, CG+AMG linear solve, and the sphere-shell validation mesh are implemented, but the
G1.1–G1.3 gate tests are not yet written and surface-Cp accuracy needs work first (see
"Known gaps" above).
