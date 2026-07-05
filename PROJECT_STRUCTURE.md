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
│   ├── residual.py       # [P1] Laplace → [P3] isentropic → [P4] Newton
│   └── upwind.py         # [P3] Artificial compressibility
├── solve/                # Linear and nonlinear solvers
│   ├── __init__.py
│   ├── linear.py         # [P1] CG + PyAMG preconditioner
│   └── newton.py         # [P4] Newton method
└── post/                 # Post-processing
    ├── __init__.py
    ├── vtk_out.py        # [P0] Write .vtu for ParaView
    └── artifacts.py      # [P0] Generate PNG + CSV for gate validation

cases/                     # Test cases and reference data
├── meshes/               # Mesh families (coarse/medium/fine)
│   ├── naca0012_2.5d/    # [M0] Extruded NACA0012 + wake surface
│   └── onera_m6/         # [M1] Swept wing
├── reference_data/       # Ground truth (DO NOT EDIT)
│   ├── sphere_incomp_cp.csv
│   ├── naca0012_cl.csv
│   └── ...
└── test_*.py             # [Deprecated] Integration tests (use tests/ now)

tests/                     # Unit and gate tests
├── conftest.py           # ✓ Pytest fixtures: artifacts_dir, mesh_dir, etc.
├── __init__.py
├── test_v0_freestream.py # ✓ [P0] Primary regression test
├── test_mesh_*.py        # [P0] Gates G0.1–G0.4
├── test_laplace_*.py     # [P1] Gates G1.1–G1.3
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
- **tests/conftest.py** — Pytest fixtures
- **tests/test_v0_freestream.py** — Smoke tests + regression baseline
- **pyproject.toml** — Build metadata and dependencies
- **.copilot-instructions.md** — AI agent instructions for P0–P5

### ⏳ Not Yet (P0)
- `pyfp3d/mesh/reader.py` — Mesh I/O (meshio)
- `pyfp3d/mesh/metrics.py` — Geometry: volumes, gradients, face adjacency
- `pyfp3d/mesh/coloring.py` — Element graph coloring for @prange
- `pyfp3d/post/vtk_out.py` — VTK writer
- `tests/test_mesh_*.py` — Gates G0.1–G0.4
- `cases/meshes/naca0012_2.5d/*` — Gmsh-generated mesh family

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

## Next Steps (P0 Completion)

1. **Implement `mesh/reader.py`** (gate G0.1)
   - Read .msh files via meshio
   - Convert to SoA (Structure-of-Arrays) format
   - Extract boundary tags

2. **Implement `mesh/metrics.py`** (gate G0.2)
   - Compute per-element volumes (tetrahedral geometry)
   - Compute element-wise gradients ∇φ
   - Build face adjacency lists

3. **Implement `mesh/coloring.py`** (gate G0.3)
   - Greedy graph coloring of element connectivity
   - Verify no same-color elements share a node

4. **Implement `post/vtk_out.py`** (gate G0.4)
   - Write nodal and cell fields to .vtu (ParaView format)
   - Round-trip test: read → write → read

5. **Create test meshes** (M0)
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
**Status:** P0 scaffolding complete; ready for mesh infrastructure implementation
