"""Track V V2 transpiration channel (binding: docs/roadmap/track_v.md
"V2 -- Transpiration channel through all three drivers", gates GV2.1(a)/(b);
design: docs/design_track_v.md; module under test:
pyfp3d/viscous/transpiration.py).

Covers:
  - assemble_transpiration_rhs exactness vs the consistent P1 mass-matrix
    load, sign convention b = -load(m_dot) included (module docstring;
    the end-to-end sign pin is the cylinder MMS below);
  - m_dot == 0 -> the EXACT zero vector (GV2.1(b) bit-identity at the
    assembly level);
  - surface_divergence_tri / surface_divergence_nodal exactness on linear
    fields + discrete conservation;
  - transpiration_from_delta_star on an analytic defect-flux field (the
    delta* -> m_dot operator's unit exercise; first LIVE use is V3);
  - edge_velocity_per_zone: default quadratic recovery vs the LE-band
    linear+smoothed branch (A4 per-zone discipline);
  - GV2.1(a) coarse-mesh regression lock: Fourier-mode blowing on the M0
    cylinder vs the analytic exterior Laplace solution (the full
    coarse/medium/fine order study lives in
    cases/analysis/v2_transpiration_channel/run.py);
  - GV2.1(b) Picard legs: solve_laplace / solve_subsonic /
    solve_subsonic_lifting with an explicit zero RHS are bit-identical to
    the channel-absent defaults, and a nonzero RHS is self-consistent
    (assemble_residual(phi) - b ~ 0 at convergence).

Runs in both lanes: default JIT and PYFP3D_NOJIT=1.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.picard import solve_laplace, solve_subsonic, solve_subsonic_lifting
from pyfp3d.viscous.surface_mesh import SurfaceMesh, structured_rectangle_surface
from pyfp3d.viscous.transpiration import (
    assemble_transpiration_rhs,
    edge_velocity_per_zone,
    surface_divergence_nodal,
    surface_divergence_tri,
    transpiration_from_delta_star,
)

from .mesh_utils import (
    cylinder_blowing_m_dot,
    cylinder_blowing_phi_exact,
    cylinder_phi_exact,
)

REPO_ROOT = Path(__file__).parent.parent
CYL_DIR = REPO_ROOT / "cases" / "meshes" / "cylinder_2.5d"
NACA_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"

V0, N_MODE = 0.1, 2


# ---------------------------------------------------------------------------
# assembly-level exactness (structured flat surface)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def flat_surface():
    xyz, tris = structured_rectangle_surface(0.0, 2.0, 0.0, 1.0, 8, 4)
    return xyz, tris, SurfaceMesh.from_wall_faces(xyz, tris)


def test_rhs_assembly_matches_consistent_load(flat_surface):
    """b = -<N_i, m_dot> with the EXACT P1 mass-matrix load (the 3-point
    edge-midpoint rule integrates N_i * m_dot exactly): b == -M @ m to
    roundoff."""
    xyz, tris, _ = flat_surface
    m = 1.0 + 2.0 * xyz[:, 0] + 3.0 * xyz[:, 2]
    rhs = assemble_transpiration_rhs(xyz, tris, m)
    load_exact = np.zeros(len(xyz))
    mass = (np.ones((3, 3)) + np.eye(3)) / 12.0
    for e in range(len(tris)):
        nn = tris[e]
        x0, x1, x2 = xyz[nn]
        area = 0.5 * np.linalg.norm(np.cross(x1 - x0, x2 - x0))
        np.add.at(load_exact, nn, area * (mass @ m[nn]))
    assert np.max(np.abs(rhs + load_exact)) < 1e-14


def test_rhs_zero_m_dot_is_exact_zero(flat_surface):
    """The GV2.1(b) assembly-level clause: m_dot == 0 assembles to the
    exact zero vector -- passing it through any driver is bit-identical
    to the channel being absent."""
    xyz, tris, _ = flat_surface
    rhs = assemble_transpiration_rhs(xyz, tris, np.zeros(len(xyz)))
    assert np.array_equal(rhs, np.zeros(len(xyz)))


def test_surface_divergence_linear_exact_and_conservative(flat_surface):
    """Strong-form surface divergence: exact per triangle for a P1 vector
    field; the lumped nodal projection is then exact too, and discretely
    conservative (sum node_area * div_nodal == sum area_tri * div_tri)."""
    _, _, sm = flat_surface
    vec = np.stack(
        [2.0 * sm.xyz[:, 0] + sm.xyz[:, 2],
         np.zeros(sm.n_node),
         sm.xyz[:, 0] - 3.0 * sm.xyz[:, 2]],
        axis=1,
    )  # div = 2 - 3 = -1 exactly
    div_tri = surface_divergence_tri(sm, vec)
    assert np.max(np.abs(div_tri + 1.0)) == 0.0
    div_nodal = surface_divergence_nodal(sm, vec)
    assert np.max(np.abs(div_nodal + 1.0)) == 0.0
    assert abs(
        np.sum(sm.node_area * div_nodal) - np.sum(sm.area_tri * div_tri)
    ) < 1e-15


def test_transpiration_from_delta_star_analytic(flat_surface):
    """delta* -> m_dot operator: rho = 1, u_e = x-hat, delta* = x gives
    Q = (x, 0, 0), div Q = 1 (the 2-D rule v_n = d(u_e delta*)/ds)."""
    _, _, sm = flat_surface
    ue = np.tile([1.0, 0.0, 0.0], (sm.n_node, 1))
    mdot = transpiration_from_delta_star(
        sm, np.ones(sm.n_node), ue, sm.xyz[:, 0]
    )
    assert np.max(np.abs(mdot - 1.0)) == 0.0


def test_edge_velocity_per_zone(flat_surface):
    """Default = quadratic recovery (exact for a quadratic field); the
    LE-band mask switches those nodes to the linear+smoothed branch
    (exact for a linear field, measurably NOT exact for a quadratic one
    -- the branch really switched)."""
    xyz, tris, _ = flat_surface
    # quadratic field: u_ex = (2x + 2z, 0, 2x - z)
    phi_q = xyz[:, 0] ** 2 + 2.0 * xyz[:, 0] * xyz[:, 2] - 0.5 * xyz[:, 2] ** 2
    ue_ex = np.stack(
        [2.0 * xyz[:, 0] + 2.0 * xyz[:, 2],
         np.zeros(len(xyz)),
         2.0 * xyz[:, 0] - xyz[:, 2]],
        axis=1,
    )
    ue = edge_velocity_per_zone(xyz, tris, phi_q)
    assert np.nanmax(np.abs(ue - ue_ex)) < 1e-12

    mask = xyz[:, 0] < 0.3
    ue_z = edge_velocity_per_zone(xyz, tris, phi_q, le_band_mask=mask)
    assert np.nanmax(np.abs(ue_z[~mask] - ue_ex[~mask])) < 1e-12
    assert np.nanmax(np.abs(ue_z[mask] - ue_ex[mask])) > 1e-6

    # linear field: exact through BOTH branches
    phi_l = 2.0 * xyz[:, 0] - xyz[:, 2]
    ue_l_ex = np.tile([2.0, 0.0, -1.0], (len(xyz), 1))
    ue_l = edge_velocity_per_zone(xyz, tris, phi_l, le_band_mask=mask)
    assert np.nanmax(np.abs(ue_l - ue_l_ex)) < 1e-12


# ---------------------------------------------------------------------------
# GV2.1(a) -- manufactured Fourier blowing on the M0 cylinder (coarse lock)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def cylinder_coarse():
    return read_mesh(CYL_DIR / "coarse.msh")


def _solve_cylinder_blowing(mesh):
    nodes, elements = mesh.nodes, mesh.elements
    wall_faces = mesh.boundary_faces["wall"]
    ff_nodes = np.unique(mesh.boundary_faces["farfield"])
    phi_ex = cylinder_blowing_phi_exact(nodes, v0=V0, n_mode=N_MODE)
    rhs = assemble_transpiration_rhs(
        nodes, wall_faces, cylinder_blowing_m_dot(nodes, v0=V0, n_mode=N_MODE)
    )
    result = solve_laplace(
        nodes, elements, ff_nodes, phi_ex[ff_nodes],
        body_source_rhs=rhs, rtol=1e-11, maxiter=3000,
    )
    return result, phi_ex


def test_cylinder_fourier_blowing_mms_coarse(cylinder_coarse):
    """GV2.1(a) coarse regression lock: rel max-norm phi error vs the
    analytic Fourier solution. Measured 2.1572e-02 at coarse (h_wall =
    0.10); the gate's coarse/medium/fine order study is the run.py
    evidence. The band doubles as the SIGN pin: a flipped transpiration
    sign lands at O(2), two orders outside the lock."""
    result, phi_ex = _solve_cylinder_blowing(cylinder_coarse)
    err = np.abs(result["phi"] - phi_ex)
    relmax = err.max() / np.abs(phi_ex).max()
    assert result["residual_norm"] < 1e-8
    assert relmax < 0.03, f"coarse MMS relmax {relmax:.4e} outside 0.03 lock"


# ---------------------------------------------------------------------------
# GV2.1(b) -- Picard legs: explicit zero RHS bit-identity + nonzero sanity
# ---------------------------------------------------------------------------


def test_solve_laplace_zero_rhs_bit_identical(cylinder_coarse):
    """solve_laplace: body_source_rhs=zeros is bit-identical to None."""
    mesh = cylinder_coarse
    nodes, elements = mesh.nodes, mesh.elements
    ff_nodes = np.unique(mesh.boundary_faces["farfield"])
    phi_ff = cylinder_phi_exact(nodes)[ff_nodes]
    a = solve_laplace(nodes, elements, ff_nodes, phi_ff,
                      rtol=1e-11, maxiter=3000)
    b = solve_laplace(nodes, elements, ff_nodes, phi_ff,
                      body_source_rhs=np.zeros(len(nodes)),
                      rtol=1e-11, maxiter=3000)
    assert np.array_equal(a["phi"], b["phi"])
    assert a["residual_norm"] == b["residual_norm"]


def test_solve_subsonic_channel(cylinder_coarse):
    """solve_subsonic (compressible Picard threading): zeros bit-identical
    to None; the nonzero transpiration load converges self-consistently
    (recorded residual is R - b by construction)."""
    mesh = cylinder_coarse
    nodes, elements = mesh.nodes, mesh.elements
    wall_faces = mesh.boundary_faces["wall"]
    ff_nodes = np.unique(mesh.boundary_faces["farfield"])
    phi_ff = nodes[ff_nodes, 0].copy()  # unit freestream far field

    a = solve_subsonic(nodes, elements, ff_nodes, phi_ff, m_inf=0.3)
    b = solve_subsonic(nodes, elements, ff_nodes, phi_ff, m_inf=0.3,
                       body_source_rhs=np.zeros(len(nodes)))
    assert np.array_equal(a["phi"], b["phi"])
    assert a["residual_history"] == b["residual_history"]

    rhs = assemble_transpiration_rhs(
        nodes, wall_faces, cylinder_blowing_m_dot(nodes, v0=0.02, n_mode=2)
    )
    c = solve_subsonic(nodes, elements, ff_nodes, phi_ff, m_inf=0.3,
                       body_source_rhs=rhs)
    assert c["converged"]
    assert c["residual_history"][-1] < 1e-8
    assert not np.array_equal(a["phi"], c["phi"])


@pytest.fixture(scope="module")
def naca_coarse_cut():
    mesh = read_mesh(NACA_DIR / "coarse.msh")
    return cut_wake(mesh)


def test_solve_subsonic_lifting_channel(naca_coarse_cut):
    """solve_subsonic_lifting (wake-cut Picard threading, the reduced_rhs
    T^T route): zeros bit-identical; the nonzero transpiration load
    converges with a shifted solution."""
    mc, wc = naca_coarse_cut
    n = len(mc.nodes)
    a = solve_subsonic_lifting(mc, wc, m_inf=0.3, alpha_deg=2.0)
    b = solve_subsonic_lifting(mc, wc, m_inf=0.3, alpha_deg=2.0,
                               body_source_rhs=np.zeros(n))
    assert np.array_equal(a["phi"], b["phi"])
    assert np.array_equal(a["gamma"], b["gamma"])
    assert a["residual_history"] == b["residual_history"]

    mdot = np.full(n, 0.0)
    mdot[np.unique(mc.boundary_faces["wall"])] = 0.01
    rhs = assemble_transpiration_rhs(mc.nodes, mc.boundary_faces["wall"], mdot)
    c = solve_subsonic_lifting(mc, wc, m_inf=0.3, alpha_deg=2.0,
                               body_source_rhs=rhs)
    assert c["converged"] and c["kutta_converged"]
    assert c["residual_history"][-1] < 1e-8
    assert not np.array_equal(a["phi"], c["phi"])
