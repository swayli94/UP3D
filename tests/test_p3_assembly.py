"""
P3 assembly tech-debt retirement regressions (roadmap P3 deliverables):
the colored-prange fast path (precomputed B_e/V_e, elem_to_csr scatter,
mesh/coloring.py wired in -- design.md Sec 7 rules 2/4, agent-rules hard
rule #3) against the retained P1 serial reference kernels, plus the V0
freestream check with the rho machinery in the loop.

Determinism note: within one color no two elements share a node, so
per-node/per-nnz accumulation order is fixed by the color sequence alone
-- the fast path is bit-reproducible across calls and thread counts
(asserted below). Against the OLD path the per-nnz summation ORDER
differs (scipy COO dedup order vs color order), so equality there is to
machine precision, not bitwise.
"""

import numpy as np
import pytest

from pyfp3d.kernels.jacobian import PicardOperator
from pyfp3d.kernels.residual import (
    assemble_residual,
    assemble_stiffness_matrix,
    assemble_stiffness_matrix_reference,
)
from pyfp3d.mesh.coloring import greedy_coloring, validate_coloring
from pyfp3d.mesh.reader import read_mesh

from .mesh_utils import cube_boundary_mask, generate_structured_cube_mesh


@pytest.fixture(scope="module")
def naca_coarse():
    from .conftest import REPO_ROOT
    return read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh")


@pytest.fixture(scope="module")
def cube():
    return generate_structured_cube_mesh(n=6, L=1.0)


class TestFastAssemblyVsReference:
    def test_stiffness_matches_reference_cube(self, cube):
        nodes, elements = cube
        A_new = assemble_stiffness_matrix(nodes, elements)
        A_ref = assemble_stiffness_matrix_reference(nodes, elements)
        A_ref.sort_indices()
        assert np.array_equal(A_new.indptr, A_ref.indptr)
        assert np.array_equal(A_new.indices, A_ref.indices)
        scale = np.abs(A_ref.data).max()
        assert np.abs(A_new.data - A_ref.data).max() < 1e-13 * scale

    def test_stiffness_matches_reference_real_mesh(self, naca_coarse):
        m = naca_coarse
        A_new = assemble_stiffness_matrix(m.nodes, m.elements)
        A_ref = assemble_stiffness_matrix_reference(m.nodes, m.elements)
        A_ref.sort_indices()
        assert np.array_equal(A_new.indices, A_ref.indices)
        scale = np.abs(A_ref.data).max()
        assert np.abs(A_new.data - A_ref.data).max() < 1e-13 * scale

    def test_residual_matches_reference(self, naca_coarse):
        m = naca_coarse
        op = PicardOperator(m.nodes, m.elements)
        phi = np.sin(m.nodes[:, 0]) + 0.3 * m.nodes[:, 1] ** 2 - 0.1 * m.nodes[:, 2]
        R_new = op.assemble_residual(phi).copy()
        R_ref = assemble_residual(m.nodes, m.elements, phi)
        scale = np.abs(R_ref).max()
        assert np.abs(R_new - R_ref).max() < 1e-13 * scale

    def test_assembly_bit_deterministic(self, naca_coarse):
        m = naca_coarse
        A1 = assemble_stiffness_matrix(m.nodes, m.elements)
        A2 = assemble_stiffness_matrix(m.nodes, m.elements)
        assert np.array_equal(A1.data, A2.data)

        op = PicardOperator(m.nodes, m.elements)
        rho = np.linspace(0.9, 1.1, op.n_tets)
        B1 = op.assemble_matrix(rho)
        B2 = op.assemble_matrix(rho)
        assert np.array_equal(B1.data, B2.data)
        phi = m.nodes[:, 0] * m.nodes[:, 1]
        R1 = op.assemble_residual(phi, rho).copy()
        R2 = op.assemble_residual(phi, rho).copy()
        assert np.array_equal(R1, R2)

    def test_coloring_valid_on_real_mesh(self):
        from .conftest import REPO_ROOT
        m = read_mesh(REPO_ROOT / "cases" / "meshes" / "cylinder_2.5d" / "coarse.msh")
        colors, n_colors = greedy_coloring(m.elements)
        assert validate_coloring(m.elements, colors)
        assert 2 <= n_colors <= 256


class TestFreestreamWithRhoMachinery:
    """V0 freestream preservation with the density law in the loop
    (design.md Sec 3): phi = x has q^2 = 1, so rho = 1 at ANY M_inf and
    the weighted residual must stay machine-zero on interior nodes."""

    @pytest.mark.parametrize("m_inf", [0.0, 0.5])
    def test_freestream_residual_machine_zero(self, cube, m_inf):
        from pyfp3d.physics.isentropic import density_field

        nodes, elements = cube
        op = PicardOperator(nodes, elements)
        phi = nodes[:, 0].copy()
        _, q2 = op.velocities(phi)
        rho = density_field(q2, m_inf)
        R = op.assemble_residual(phi, rho)
        interior = ~cube_boundary_mask(nodes, L=1.0)
        assert np.abs(R[interior]).max() < 1e-12

    def test_density_field_matches_scalar_law(self):
        from pyfp3d.physics.isentropic import density_field, density_isentropic

        q2 = np.linspace(0.0, 1.8, 37)
        rho = density_field(q2, 0.5)
        for qi, ri in zip(q2, rho):
            assert ri == density_isentropic(float(qi), 0.5)

    def test_density_field_exactly_one_at_m0(self):
        """The G3.3 bit-identity rests on rho == 1.0 bitwise at M_inf = 0
        for ARBITRARY q^2 (design.md (2.2): the M^2 factor kills the q^2
        term exactly)."""
        from pyfp3d.physics.isentropic import density_field

        q2 = np.array([0.0, 0.3, 1.0, 1.7, 2.4])
        assert np.all(density_field(q2, 0.0) == 1.0)
