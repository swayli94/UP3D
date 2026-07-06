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
    └── surface.py        # [P1] nodal_gradient_recovery() (volume-weighted, for interior fields)
                            #      and wall_tangential_gradient() (surface-only, for wall Cp --
                            #      the accurate one; see "Known gaps" for why it still isn't
                            #      accurate *enough* to close G1.2)

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
├── test_mesh_adjacency.py           # ✓ [P0] Regression test for build_face_adjacency fix
├── test_mesh_reader_roundtrip.py    # ✓ [P0] Regression test for write_mesh tag-loss fix
├── test_laplace_mms.py              # ✓ [P1] Gate G1.1 -- PASSES
├── test_laplace_cg_iterations.py    # ✓ [P1] Gate G1.3 -- PASSES
├── test_laplace_sphere.py           # ✓ [P1] Gate G1.2 -- strict xfail, see "Known gaps"
├── test_laplace_picard.py           # ✓ [P1] Regression test for solve_laplace residual_norm fix
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
- **tests/test_mesh_adjacency.py**, **tests/test_mesh_reader_roundtrip.py** — regression tests for
  two bugs found by manual audit (see below)
- **pyproject.toml** — Build metadata and dependencies
- **.copilot-instructions.md** — AI agent instructions for P0–P5

Three latent bugs found by manual code audit (not caught by the existing suite, because nothing
exercised these code paths) have been fixed, each now with a regression test:
- `mesh/metrics.py::build_face_adjacency` crashed under `@njit` (reflected-list dict values are
  not valid in numba nopython mode) — rewritten around a `numba.typed.Dict` keyed by sorted face,
  storing only the packed first-owner index instead of a growing list. (`test_mesh_adjacency.py`)
- `mesh/reader.py::write_mesh` silently dropped every named boundary group (`wall`, `farfield`,
  ...), writing only a legacy `"all_triangles"` block, and `.msh` was ambiguous between meshio's
  `ansys`/`gmsh` writers (defaulted to `ansys`, discarding all tag data). Now writes each boundary
  group as its own tagged `triangle` block plus `gmsh:physical`/`field_data`, explicitly via the
  `gmsh22` writer. (`test_mesh_reader_roundtrip.py`)
- `solve/picard.py::solve_laplace` reported `residual_norm` over *all* nodes, including Dirichlet
  (far-field) rows whose natural-BC flux imbalance is O(1) and never shrinks — swamping the actual
  free-dof residual (which was already converging to ~1e-10). Now restricted to free dofs and
  correctly nets out `body_source_rhs` when present. (`test_laplace_picard.py`)

### ✓ P1 gates G1.1 and G1.3 closed; G1.2 open (see below)
- **pyfp3d/kernels/residual.py** — Laplace residual (6.1) + SPD stiffness matrix (6.2) assembly
- **pyfp3d/solve/linear.py** — Dirichlet elimination (principal-submatrix reduction) + CG+PyAMG
- **pyfp3d/solve/picard.py** — `solve_laplace()` driver (P1's Picard loop degenerates to one solve)
- **pyfp3d/post/surface.py** — `nodal_gradient_recovery()` (volume-weighted, interior fields) and
  `wall_tangential_gradient()` (surface-only recovery, for wall Cp)
- **tests/mesh_utils.py**, **cases/meshes/sphere_shell/** — structured-cube (MMS) and sphere-shell
  (G1.2) mesh generators; sphere-shell coarse/medium `.msh` + inspection PNGs are committed
- **tests/test_laplace_mms.py** — Gate G1.1 (MMS convergence) ✓ — L2 slope ≈ 1.94–1.96 with a
  sin·cos manufactured solution and a proper 4-point quadrature-consistent load vector. (A
  harmonic-polynomial exact solution was tried first and rejected: this codebase's structured
  Kuhn-triangulated cube reproduces harmonic quadratics to machine precision at *every* h, giving
  zero convergence-order signal — the same reason central finite differences are exact for
  quadratics.)
- **tests/test_laplace_cg_iterations.py** — Gate G1.3 (CG+AMG mesh-independence) ✓ — iterations
  8→11→14 across an 8×/level node-count increase (n=8,16,32 cube), comfortably under a 2× cap.

**G1.2 (incompressible sphere Cp) is still open** — `tests/test_laplace_sphere.py::test_sphere_cp_medium_mesh`
is a `strict=True` xfail against the real <2% criterion, not a loosened threshold:
- The original `nodal_gradient_recovery` (volume-weighted average of the one-sided tets touching
  each wall node) gave ~26% max / 9% mean Cp error on the medium mesh — systematically low,
  because tangential velocity physically decays moving away from the wall and every incident tet
  sits on the inward side.
- A local least-squares ("SPR-style") extrapolation fix was tried and **rejected**: on a real
  graded mesh, some nodes' 1-ring element patches are nearly coplanar, and extrapolating a linear
  model through an ill-conditioned patch blew one node's Cp error up to 429 (found on the medium
  mesh). Don't resurrect this approach without a robust conditioning safeguard.
- The fix that shipped, `wall_tangential_gradient()`, computes the surface-tangential gradient
  directly from the wall's own P1 triangulation (physically exact for a natural-BC wall, where the
  normal component is zero by construction) instead of extrapolating the volume gradient. This
  improved medium-mesh max error to ~12%, and is a genuine, well-conditioned fix (no
  extrapolation, only interpolation/averaging) — but a mesh-refinement sweep (h_min = 0.08 → 0.04
  → 0.02 → 0.015, up to ~1.2M nodes / 7.4M tets) only reaches ~3.6% max error and visibly
  saturates rather than continuing to converge, so the remaining gap is *not* simply "refine the
  mesh more."
- Still open: whether the saturation is from near-wall graded-mesh element quality, the faceted
  (vs. true curved) wall geometry, or a genuine need for quadratic (not linear) surface patch
  recovery. Worth investigating before spending more mesh-refinement compute on it.

### ⏳ Next
- Investigate the G1.2 saturation (see above) — likely needs mesh-quality diagnostics on the
  sphere-shell graded region, or a quadratic surface recovery scheme, not more h-refinement.
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

P0 (mesh I/O, metrics, coloring, VTK writer, gates G0.1–G0.4) is done. G1.1 and G1.3 are done. To
close P1, only G1.2 remains:

1. **Investigate the G1.2 Cp-accuracy saturation** (see "Known gaps" above) — mesh refinement
   alone plateaus around ~3.6% max error, well above the 2% target. Candidate next steps: check
   near-wall element quality/aspect ratio in the graded region, or implement a quadratic (not
   linear) surface patch-recovery scheme in `post/surface.py`.

2. **Create test meshes** (M0, parallel track)
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
**Status:** P0 closed (G0.1–G0.4 green). P1: G1.1 (MMS) and G1.3 (CG+AMG mesh-independence) closed;
G1.2 (sphere Cp) open with a `strict=True` xfail tracking the real 2% criterion — see "Known gaps"
for what's been tried and ruled out. 41 tests total (40 passed, 1 xfailed), full suite ~9s.
