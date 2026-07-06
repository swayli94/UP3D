"""
P4 upwinding unit tests + gate G4.2 (subcritical no-op).

G4.2 (roadmap): with M_inf = 0.5 the nu switch must be identically zero
and the results BIT-IDENTICAL to P3. The comparison runs the same driver
with upwind_c = 1.5 (P4 machinery active: upstream walk + rho_tilde
sweep in the loop) vs upwind_c = 0.0 (the literal P3 code path, sweep
bypassed): in the subcritical range nu == 0.0 exactly and
rho - 0.0*(rho - rho_u) == rho bitwise, so everything downstream of the
density sweep is unchanged bit-for-bit.
"""

import csv

import numpy as np
import pytest

from pyfp3d.kernels.jacobian import PicardOperator
from pyfp3d.kernels.upwind import UpwindOperator
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.physics.isentropic import density_field, q2_at_mach
from pyfp3d.solve.picard import solve_subsonic_lifting

from .mesh_utils import cube_boundary_mask, generate_structured_cube_mesh


@pytest.fixture(scope="module")
def cube():
    return generate_structured_cube_mesh(n=6, L=1.0)


class TestUpstreamWalk:
    def test_upstream_is_upstream(self, cube):
        """Uniform +x flow: every element's u(e) has a strictly negative
        streamwise centroid displacement (or is e itself at the inflow
        boundary), and the multi-hop reach covers ~the element extent."""
        nodes, elements = cube
        upw = UpwindOperator(nodes, elements)
        op = PicardOperator(nodes, elements)
        phi = nodes[:, 0].copy()
        grad, _ = op.velocities(phi)

        from pyfp3d.kernels.upwind import upstream_elements
        u = np.empty(op.n_tets, dtype=np.int64)
        upstream_elements(upw.face_neighbors, upw.centroids,
                          upw._nodes, upw._elements, grad, u)
        disp = (upw.centroids[u] - upw.centroids)[:, 0]  # flow is +x
        self_up = u == np.arange(op.n_tets)
        assert np.all(disp[~self_up] < 0.0)
        # self-upstream only at the upstream (x=0) boundary layer of cells
        assert np.all(upw.centroids[self_up, 0] < 0.25)

    def test_nu_exactly_zero_subcritical(self, cube):
        nodes, elements = cube
        upw = UpwindOperator(nodes, elements)
        op = PicardOperator(nodes, elements)
        phi = nodes[:, 0].copy()
        grad, q2 = op.velocities(phi)
        rho = density_field(q2, 0.5)
        rho_t = upw.rho_tilde(grad, q2, rho, 0.5, 1.5, 0.95)
        assert upw.nu_max == 0.0
        assert np.array_equal(rho_t, rho)

    def test_freestream_preservation_with_upwind(self, cube):
        """phi = x at M_inf = 0.8: q^2 = 1 is subcritical, nu == 0, and
        the upwinded residual stays machine-zero on interior nodes
        (design.md Sec 3 properties; hard rule 1 family)."""
        nodes, elements = cube
        upw = UpwindOperator(nodes, elements)
        op = PicardOperator(nodes, elements)
        phi = nodes[:, 0].copy()
        grad, q2 = op.velocities(phi)
        rho = density_field(q2, 0.8)
        rho_t = upw.rho_tilde(grad, q2, rho, 0.8, 1.5, 0.95)
        R = op.assemble_residual(phi, rho_t)
        interior = ~cube_boundary_mask(nodes, L=1.0)
        assert np.abs(R[interior]).max() < 1e-12

    def test_q2_limiter_inactive_below_cap(self):
        from pyfp3d.physics.isentropic import limit_q2_field

        q2 = np.linspace(0.0, 2.0, 11)
        assert np.array_equal(limit_q2_field(q2, 0.5, 3.0), q2)
        assert q2_at_mach(3.0, 0.0) == np.inf
        capped = limit_q2_field(np.array([50.0]), 0.8, 3.0)
        assert capped[0] == q2_at_mach(3.0, 0.8)


def test_g42_subcritical_noop_bitwise(artifacts_dir):
    """Gate G4.2: M_inf = 0.5 with the upwind machinery in the loop is
    bit-identical to the P3 path (upwind_c = 0), and nu == 0 everywhere."""
    from .conftest import REPO_ROOT

    mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh")
    mc, wc = cut_wake(mesh)

    r_p4 = solve_subsonic_lifting(mc, wc, m_inf=0.5, alpha_deg=2.0,
                                  upwind_c=1.5)
    r_p3 = solve_subsonic_lifting(mc, wc, m_inf=0.5, alpha_deg=2.0,
                                  upwind_c=0.0)

    bits_phi = np.array_equal(r_p4["phi"], r_p3["phi"])
    bits_gamma = np.array_equal(r_p4["gamma"], r_p3["gamma"])

    gate_dir = artifacts_dir / "G4.2"
    gate_dir.mkdir(parents=True, exist_ok=True)
    with open(gate_dir / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["check", "value"])
        w.writerow(["nu_max_at_M0.5", r_p4["nu_max"]])
        w.writerow(["n_nu_active", r_p4["n_nu_active"]])
        w.writerow(["phi_bitwise_identical", bits_phi])
        w.writerow(["gamma_bitwise_identical", bits_gamma])
        w.writerow(["max_abs_phi_diff",
                    float(np.max(np.abs(r_p4["phi"] - r_p3["phi"])))])
        w.writerow(["mach2_max", r_p4["mach2_max"]])

    assert r_p4["nu_max"] == 0.0
    assert r_p4["n_nu_active"] == 0
    assert bits_phi and bits_gamma
    assert r_p4["mach2_max"] < 0.95**2
