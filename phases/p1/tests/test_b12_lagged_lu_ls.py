"""
Track B / B12 gates G12.1 + G12.2: the lagged-LU direct-reuse path on the
level-set Newton solver (`solve_multivalued_newton`, precond=None).

Motivation (B11 measured, roadmap B12): on the fused 3D matrix the iterative
escapes fail beyond coarse -- ILU diverges at 2.5D medium lifting and
`factor_failed`s at M6 medium, AMG stalls throughout -- so at medium/M6 sizes
sparse-direct is the only converging tool and the cost driver is the NUMBER of
factorizations. `direct_refactor_every > 1` refactors the LU every k-th Newton
step and drives the intermediate steps with GMRES preconditioned by the stale
(exact) LU (the N6 lagged-LU mechanism ported from solve/newton.py, minus the
Woodbury -- the level-set system has no Gamma coupling).

  * G12.1 bit-identity: direct_refactor_every=1 (default) is the byte-identical
    per-step spsolve; the two new params default to 1 / 1e-8.
  * G12.2 equivalence: k=2 reaches the same converged gamma as the spsolve
    default to |dgamma| < 1e-8, 0 lim/flr, 0 GMRES stalls, and it actually
    reused the stale LU (n_refactor < n_newton and a reuse GMRES step ran).

All on the committed coarse 2.5D mesh (cheap). M6-scale timing is the demo's
job (cases/demo/b12_lagged_lu), not pytest.
"""

import inspect
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.solve.newton_ls import solve_multivalued_newton
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
ALPHA = 2.0

# seed enough that Newton converges but still takes several steps (so the
# lagged reuse path is exercised for k > 1)
COMMON = dict(m_inf=0.7, alpha_deg=ALPHA, farfield="neumann",
              n_seed=20, n_newton_max=25)


def _mesh(level="coarse"):
    path = M0_DIR / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    return read_mesh(path)


def _mvop(mesh, alpha=ALPHA):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


# ---------------------------------------------------------------------------
# G12.1 -- bit-identity convention.
# ---------------------------------------------------------------------------

def test_lagged_lu_param_defaults():
    p = inspect.signature(solve_multivalued_newton).parameters
    assert p["direct_refactor_every"].default == 1
    assert p["direct_reuse_rtol"].default == 1e-8


def test_default_bit_identical_to_spsolve():
    """direct_refactor_every=1 is exactly the spsolve default -- byte-identical
    solution, and no factorization-reuse bookkeeping engaged."""
    mesh = _mesh()
    a = solve_multivalued_newton(mesh=mesh, mvop=_mvop(mesh), **COMMON)
    b = solve_multivalued_newton(mesh=mesh, mvop=_mvop(mesh),
                                 direct_refactor_every=1, **COMMON)
    assert a["converged"] and b["converged"]
    assert np.array_equal(a["phi_ext"], b["phi_ext"])
    assert a["n_refactor"] == 0 and b["n_refactor"] == 0
    assert a["n_gmres_total"] == 0 and b["n_gmres_total"] == 0


# ---------------------------------------------------------------------------
# G12.2 -- lagged-LU == spsolve, and the LU is actually reused.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("k", [2, 1000])
def test_lagged_lu_matches_spsolve(k):
    mesh = _mesh()
    ref = solve_multivalued_newton(mesh=mesh, mvop=_mvop(mesh), **COMMON)
    got = solve_multivalued_newton(mesh=mesh, mvop=_mvop(mesh),
                                   direct_refactor_every=k, **COMMON)
    assert ref["converged"] and got["converged"]
    assert abs(got["gamma"] - ref["gamma"]) < 1e-8, (got["gamma"], ref["gamma"])
    assert got["n_limited"] == 0 and got["n_floored"] == 0
    assert got["n_gmres_stalled"] == 0
    # the stale LU was actually reused: fewer factorizations than Newton
    # iterations, and at least one reuse GMRES step ran
    assert got["n_refactor"] < got["n_newton"], (got["n_refactor"],
                                                 got["n_newton"])
    assert got["n_gmres_total"] > 0
    # k=1000 (once-per-solve) refactors only when a reuse fails; near a coarse
    # M0.70 solution the stale LU is near-exact, so it stays at a single factor
    if k == 1000:
        assert got["n_refactor"] <= 2, got["n_refactor"]
