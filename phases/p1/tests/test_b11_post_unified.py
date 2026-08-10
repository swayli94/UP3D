"""
Track B / B11 gate G11.2: the unified post-processing dispatch
(`pyfp3d.post.unified`) is BIT-IDENTICAL to the legacy per-path functions on
both the conforming and the level-set paths.

B11 merges the shared cores of `post/surface.py` and `post/surface_ls.py`
(via `_cp_from_q2`, `_pressure_force`, `_wall_plane_crossings`,
`_section_curve_dict`, `_d11_wall_state`) under one keyword-dispatched upper
layer. These checks lock that the merge changed no number: every unified output
is `np.array_equal` (or exact-float-equal for scalars) to the function it
replaces. Solve-free -- synthetic states suffice for a numerical-identity check.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.post import unified
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.post.surface_ls import (
    cl_pressure_3d_levelset,
    section_cp_curve_levelset,
    wall_cp_levelset,
)
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"   # wake-embedded 2.5D
ALPHA = 2.0
M_INF = 0.5


def _mesh():
    path = M0_DIR / "coarse.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    return read_mesh(path)


def _conforming_state(mesh):
    """A non-trivial single-valued nodal field (freestream-like)."""
    return mesh.nodes[:, 0] + 0.3 * mesh.nodes[:, 1]


def _levelset_state(mesh):
    """mvop + a physically-scaled multivalued state: freestream on the main
    dofs plus a constant jump on the aux copies (the discrete shape of a
    circulation) -- same construction as the B7 D11 lock."""
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    phi = np.zeros(mvop.n_total)
    phi[: mvop.n_main] = mesh.nodes[:, 0]
    cut = np.flatnonzero(cm.ext_dof_of_node >= 0)
    phi[cm.ext_dof_of_node[cut]] = phi[cut] + 0.12
    return mvop, phi


def _zmid(mesh):
    z = mesh.nodes[:, 2]
    return 0.5 * (float(z.min()) + float(z.max()))


# ---------------------------------------------------------------------------
# Unified == legacy, both paths.
# ---------------------------------------------------------------------------

def test_unified_wall_cp_conforming_equals_legacy():
    mesh = _mesh()
    phi = _conforming_state(mesh)
    u = unified.wall_cp(mesh, phi=phi, m_inf=M_INF)
    leg = wall_force_coefficients(mesh.nodes, mesh.elements,
                                  mesh.boundary_faces["wall"], phi, m_inf=M_INF)
    assert np.array_equal(u["cp"], leg["cp_tri"])
    assert u["q2"].shape == u["cp"].shape
    assert u["upper"].dtype == bool


def test_unified_wall_cp_levelset_equals_legacy():
    mesh = _mesh()
    mvop, phi = _levelset_state(mesh)
    u = unified.wall_cp(mesh, mvop=mvop, phi_ext=phi, m_inf=M_INF)
    leg = wall_cp_levelset(mesh, mvop, phi, m_inf=M_INF)
    for k in ("x", "cp", "upper", "area", "n_out", "q2"):
        assert np.array_equal(u[k], leg[k]), k


def test_unified_forces_conforming_equals_legacy():
    mesh = _mesh()
    phi = _conforming_state(mesh)
    u = unified.wall_forces(mesh, phi=phi, alpha_deg=ALPHA, m_inf=M_INF,
                            s_ref=2.0)
    leg = wall_force_coefficients(mesh.nodes, mesh.elements,
                                  mesh.boundary_faces["wall"], phi,
                                  alpha_deg=ALPHA, m_inf=M_INF, s_ref=2.0)
    assert u["cl"] == leg["cl"]
    assert u["cd_pressure"] == leg["cd_pressure"]
    assert np.array_equal(u["cf"], leg["cf"])
    assert np.array_equal(u["cp_tri"], leg["cp_tri"])


def test_unified_forces_levelset_matches_cl_pressure_3d():
    mesh = _mesh()
    mvop, phi = _levelset_state(mesh)
    cp = wall_cp_levelset(mesh, mvop, phi, m_inf=M_INF)
    leg_cl = cl_pressure_3d_levelset(mesh, cp["cp"], cp["area"], cp["n_out"],
                                     ALPHA, s_ref=2.0)
    u = unified.wall_forces(mesh, mvop=mvop, phi_ext=phi, alpha_deg=ALPHA,
                            m_inf=M_INF, s_ref=2.0)
    assert u["cl"] == leg_cl


def test_unified_section_cp_conforming_equals_legacy():
    mesh = _mesh()
    phi = _conforming_state(mesh)
    z = _zmid(mesh)
    u = unified.section_cp(mesh, phi=phi, z=z, m_inf=M_INF)
    leg = section_cp_curve(mesh, phi, z=z, m_inf=M_INF)
    for k in ("x_upper", "cp_upper", "x_lower", "cp_lower", "chord", "x_le"):
        assert np.array_equal(np.asarray(u[k]), np.asarray(leg[k])), k


def test_unified_section_cp_levelset_equals_legacy():
    mesh = _mesh()
    mvop, phi = _levelset_state(mesh)
    z = _zmid(mesh)
    u = unified.section_cp(mesh, mvop=mvop, phi_ext=phi, z=z, m_inf=M_INF)
    leg = section_cp_curve_levelset(mesh, mvop, phi, z=z, m_inf=M_INF)
    for k in ("x_upper", "cp_upper", "x_lower", "cp_lower", "chord", "x_le"):
        assert np.array_equal(np.asarray(u[k]), np.asarray(leg[k])), k


def test_ls_section_smooth_passes_zero_is_default_bitwise():
    """The new opt-in smooth_passes kwarg (LS section) is bit-identical at 0."""
    mesh = _mesh()
    mvop, phi = _levelset_state(mesh)
    z = _zmid(mesh)
    a = section_cp_curve_levelset(mesh, mvop, phi, z=z, m_inf=M_INF)
    b = section_cp_curve_levelset(mesh, mvop, phi, z=z, m_inf=M_INF,
                                  smooth_passes=0)
    assert np.array_equal(a["cp_lower"], b["cp_lower"])
    assert np.array_equal(a["cp_upper"], b["cp_upper"])


# ---------------------------------------------------------------------------
# Dispatch guards.
# ---------------------------------------------------------------------------

def test_dispatch_requires_exactly_one_path():
    mesh = _mesh()
    phi = _conforming_state(mesh)
    with pytest.raises(ValueError, match="exactly one"):
        unified.wall_cp(mesh)                       # neither
    with pytest.raises(ValueError, match="exactly one"):
        unified.wall_cp(mesh, phi=phi, mvop=object())  # both


def test_dispatch_mvop_requires_phi_ext():
    mesh = _mesh()
    mvop, _ = _levelset_state(mesh)
    with pytest.raises(ValueError, match="phi_ext"):
        unified.wall_cp(mesh, mvop=mvop)            # mvop without phi_ext
