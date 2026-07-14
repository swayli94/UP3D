"""
Track B / B11 gate G11.3: the GMRES+AMG inner solver on the level-set path
reproduces the sparse-direct `spsolve` baseline, and `precond` defaults to None
(bit-identical).

The B11 escape from the splu wall wires `solve/linear.{build_amg_preconditioner,
build_ilu_preconditioner,solve_gmres}` into `solve_multivalued_laplace`,
`solve_multivalued_lifting`, and `solve_multivalued_newton`. These checks lock:
  * precond default None on every entry point (the bit-identity convention);
  * "ilu" reproduces spsolve on the quasi-2D coarse mesh to |Δ| < 1e-8 with 0
    GMRES stalls, on Laplace / lifting / Newton -- ILU is THE escape;
  * "amg" reproduces spsolve on the SPD Laplace (continuity-weld) system only;
  * the monitors (n_gmres_total, precond) are reported.

★ MEASURED BOUNDARY (design_track_b.md §5.3): AMG built on the SPD surrogate
converges as a GMRES preconditioner for the `continuity`-closure (Laplace)
system, but NOT for the `wake_ls`-closure LIFTING/transonic/Newton operator --
its g1+g2 wake-LS + nonlinear TE-Kutta rows are convection-like, not SPD, and
the surrogate cannot model them, so GMRES stalls at the restart cap
(measured: coarse M0.5 lifting, gamma 0.0033 vs 0.139, 80/80 outers stalled,
455 s). ILU factors the real fused matrix and converges (434 iters, exact).
So the tested/shipped iterative escape on the lifting path is ILU; AMG stays
wired for the SPD Laplace case and as the recorded knob.

All on the committed coarse 2.5D meshes (cheap). M6-scale timing is the demo's
job (cases/demo/b11_ls_infra), not pytest.
"""

import inspect
import os
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.constraints.dirichlet import freestream_phi
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.solve.newton_ls import solve_multivalued_newton
from pyfp3d.solve.picard_ls import (
    solve_multivalued_laplace,
    solve_multivalued_lifting,
    solve_multivalued_transonic,
)
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
ALPHA = 2.0
M_INF = 0.5

GATES = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"
run_gates = pytest.mark.skipif(
    not GATES, reason="transonic ramp is minutes; set PYFP3D_TRANSONIC_GATES=1")


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


def _all_boundary_nodes(mesh):
    b = [np.unique(f) for f in mesh.boundary_faces.values()]
    return np.unique(np.concatenate(b))


# ---------------------------------------------------------------------------
# The bit-identity convention pin.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("fn", [solve_multivalued_laplace,
                                solve_multivalued_lifting,
                                solve_multivalued_newton])
def test_precond_default_is_none(fn):
    assert inspect.signature(fn).parameters["precond"].default is None


# ---------------------------------------------------------------------------
# GMRES == spsolve.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("precond", ["ilu", "amg"])
def test_laplace_gmres_matches_spsolve(precond):
    mesh = _mesh()
    mvop = _mvop(mesh)
    b = _all_boundary_nodes(mesh)
    vals = freestream_phi(mesh.nodes[b], ALPHA)
    ref = solve_multivalued_laplace(mvop, b, vals)
    got = solve_multivalued_laplace(mvop, b, vals, precond=precond)
    assert np.max(np.abs(got["phi"] - ref["phi"])) < 1e-8


def test_lifting_gmres_matches_spsolve():
    """ILU is the shipped iterative escape on the lifting path (AMG stalls on
    the wake_ls operator -- see the module header)."""
    mesh = _mesh()
    mvop = _mvop(mesh)
    ref = solve_multivalued_lifting(mvop, mesh, M_INF, alpha_deg=ALPHA)
    got = solve_multivalued_lifting(mvop, mesh, M_INF, alpha_deg=ALPHA,
                                    precond="ilu")
    assert ref["converged"] and got["converged"]
    assert abs(got["gamma"] - ref["gamma"]) < 1e-8, (got["gamma"], ref["gamma"])
    assert got["n_gmres_stalled"] == 0
    assert got["precond"] == "ilu"
    assert got["n_gmres_total"] > 0
    # default path reports no GMRES activity
    assert ref["precond"] is None and ref["n_gmres_total"] == 0


def test_lifting_warm_start_iters_bounded():
    """Warm-started ILU-GMRES stays well under a full restart cycle per outer."""
    mesh = _mesh()
    mvop = _mvop(mesh)
    r = solve_multivalued_lifting(mvop, mesh, M_INF, alpha_deg=ALPHA,
                                  precond="ilu", gmres_restart=60)
    assert r["converged"]
    mean_iters = r["n_gmres_total"] / max(r["n_outer"], 1)
    assert mean_iters < 60, mean_iters


def test_newton_ls_gmres_matches_spsolve():
    """Coarse LS Newton at a mildly supercritical Mach (upwind active) -- the
    ILU-GMRES step reproduces the direct step's converged solution."""
    mesh = _mesh()
    mvop = _mvop(mesh)
    ref = solve_multivalued_newton(mesh=mesh, mvop=mvop, m_inf=0.7,
                                   alpha_deg=ALPHA, n_seed=25, n_newton_max=20)
    got = solve_multivalued_newton(mesh=mesh, mvop=mvop, m_inf=0.7,
                                   alpha_deg=ALPHA, n_seed=25, n_newton_max=20,
                                   precond="ilu", seed_precond="ilu")
    assert ref["converged"] and got["converged"]
    assert abs(got["gamma"] - ref["gamma"]) < 1e-8, (got["gamma"], ref["gamma"])
    assert got["precond"] == "ilu"


@run_gates
def test_transonic_forwarding_smoke():
    """precond= rides **kwargs through solve_multivalued_transonic; a short
    coarse ramp with precond='ilu' matches the default to |Δgamma| < 1e-6."""
    mesh = _mesh()
    mvop = _mvop(mesh)
    common = dict(alpha_deg=ALPHA, m_start=0.5, dm=0.05, upwind_c=1.5,
                  n_outer_seed=40, n_outer_level=120)
    ref = solve_multivalued_transonic(mvop, mesh, 0.7, **common)
    mvop2 = _mvop(mesh)
    got = solve_multivalued_transonic(mvop2, mesh, 0.7, precond="ilu", **common)
    assert abs(got["gamma"] - ref["gamma"]) < 1e-6, (got["gamma"], ref["gamma"])
