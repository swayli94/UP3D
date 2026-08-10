"""
Track B / B13 gates GB13.1 + GB13.2: lagged-LU direct-reuse on the level-set
Picard OUTER loop (`solve_multivalued_lifting`, precond=None).

Motivation (roadmap B13): after B12, the M6-medium LS cost driver is the Picard
outer loop -- one 17.5 s spsolve per outer (seed 263 s, B11 lifting headline
455 s). B11 measured the iterative escapes to be coarse-only (ILU diverges at
2.5D medium, factor-fails at M6 medium; AMG stalls), so at medium/M6 sizes
sparse-direct is the only converging tool and the cost is the NUMBER of
factorizations. `direct_refactor_every > 1` applies the B12 lagged-LU mechanism
to the outer loop: refactor every k-th outer, stale-exact-LU-preconditioned
GMRES on the FRESH matrix in between (warm-started; the under-relaxed density
keeps matrix drift small).

  * GB13.1 bit-identity: direct_refactor_every=1 (default) is the byte-identical
    per-outer spsolve; the two new params default to 1 / 1e-8.
  * GB13.2 equivalence: k in {4, 1000} reaches the same converged gamma as the
    spsolve default to |dgamma| < 1e-8, converged, 0 GMRES stalls... note early
    outers with large density moves MAY legitimately stall-and-refactor, so the
    stall assertion is on the subsonic M0.5 case where drift is mild; the hard
    lock is n_refactor < n_outer (the LU is actually reused).

Gate IDs use the GB13.x prefix (Track V's GV precedent) -- NOT G13.x, which is
P13's heavily-cited namespace.

All on the committed coarse 2.5D mesh (cheap). M6-medium timing is the demo's
job (cases/demo/b13_lagged_picard), not pytest.
"""

import inspect
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.solve.picard_ls import solve_multivalued_lifting
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
ALPHA = 2.0
M_INF = 0.5


def _mesh(level="coarse"):
    path = M0_DIR / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    return read_mesh(path)


def _mvop(mesh):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


# ---------------------------------------------------------------------------
# GB13.1 -- bit-identity convention.
# ---------------------------------------------------------------------------

def test_lagged_lu_param_defaults():
    # direct_reuse_rtol is 1e-10 here, NOT B12's 1e-8: a Picard fixed point is
    # pinned only by its lag tolerances, so an inexact reuse step shifts the
    # stopping point (measured |dgamma| 8e-8 at rtol 1e-8).
    p = inspect.signature(solve_multivalued_lifting).parameters
    assert p["direct_refactor_every"].default == 1
    assert p["direct_reuse_rtol"].default == 1e-10


def test_default_bit_identical_to_spsolve():
    """direct_refactor_every=1 is exactly the per-outer spsolve default --
    byte-identical solution, no lagged-LU bookkeeping engaged."""
    mesh = _mesh()
    a = solve_multivalued_lifting(_mvop(mesh), mesh, M_INF, alpha_deg=ALPHA)
    b = solve_multivalued_lifting(_mvop(mesh), mesh, M_INF, alpha_deg=ALPHA,
                                  direct_refactor_every=1)
    assert a["converged"] and b["converged"]
    assert np.array_equal(a["phi_ext"], b["phi_ext"])
    assert a["n_refactor"] == 0 and b["n_refactor"] == 0
    assert a["n_gmres_total"] == 0 and b["n_gmres_total"] == 0


# ---------------------------------------------------------------------------
# GB13.2 -- lagged-LU == spsolve on the converged Picard fixed point, and the
# stale LU is actually reused.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("k", [4, 1000])
def test_lagged_lu_matches_spsolve(k):
    mesh = _mesh()
    ref = solve_multivalued_lifting(_mvop(mesh), mesh, M_INF, alpha_deg=ALPHA)
    got = solve_multivalued_lifting(_mvop(mesh), mesh, M_INF, alpha_deg=ALPHA,
                                    direct_refactor_every=k)
    assert ref["converged"] and got["converged"]
    assert abs(got["gamma"] - ref["gamma"]) < 1e-8, (got["gamma"], ref["gamma"])
    # the stale LU was actually reused: fewer factorizations than outers, and
    # at least one reuse GMRES step ran
    assert got["n_refactor"] < got["n_outer"], (got["n_refactor"],
                                                got["n_outer"])
    assert got["n_gmres_total"] > 0
    # subsonic M0.5: matrix drift per outer is mild (omega_rho=0.5), the stale
    # exact LU should carry reuse steps cleanly
    assert got["n_gmres_stalled"] == 0, got["n_gmres_stalled"]
    # k=1000 = refactor once, reuse for the whole solve
    if k == 1000:
        assert got["n_refactor"] == 1, got["n_refactor"]


def test_lagged_lu_farfield_neumann_also_works():
    """farfield="neumann" is the B6/B7 transonic + M6 3D recipe (and the M6
    headline demo condition) -- fixed Dirichlet set + fixed RHS, so the stale
    LU is a legal preconditioner there too. (The parametrized test above covers
    the default "vortex" far field, whose Dirichlet VALUES refresh each outer
    on the same node set.)"""
    mesh = _mesh()
    ref = solve_multivalued_lifting(_mvop(mesh), mesh, M_INF, alpha_deg=ALPHA,
                                    farfield="neumann")
    got = solve_multivalued_lifting(_mvop(mesh), mesh, M_INF, alpha_deg=ALPHA,
                                    farfield="neumann", direct_refactor_every=8)
    assert ref["converged"] and got["converged"]
    assert abs(got["gamma"] - ref["gamma"]) < 1e-8
    assert got["n_refactor"] < got["n_outer"]
