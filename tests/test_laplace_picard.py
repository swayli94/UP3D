"""
Regression tests for `solve.picard.solve_laplace`'s reported residual_norm.

It used to be computed over *all* nodes, including the Dirichlet far-field
rows. assemble_residual() returns the raw operator A@phi with no boundary
condition applied, so at a Dirichlet node the "residual" is really the
natural-BC flux imbalance -- an O(1) quantity that never shrinks, since that
row isn't part of the system actually being solved (it's overridden by the
prescribed value). That swamped the real free-dof residual (which was
already converging fine) and would make any future convergence diagnostic
built on it meaningless.
"""

import numpy as np
import pytest

from pyfp3d.kernels.residual import assemble_residual
from pyfp3d.solve.picard import solve_laplace

from .mesh_utils import cube_boundary_mask, generate_structured_cube_mesh


class TestResidualNormExcludesDirichletRows:
    def test_residual_norm_matches_free_dof_residual(self):
        nodes, elements = generate_structured_cube_mesh(n=6, L=1.0)
        boundary = cube_boundary_mask(nodes, L=1.0)
        dirichlet_nodes = np.where(boundary)[0]
        free = ~boundary

        # An exactly-representable linear field: converged phi should match
        # it almost exactly, and the *free*-dof residual should be tiny --
        # but the raw operator's boundary-row residual is O(1) regardless
        # (flux imbalance from the prescribed Dirichlet data), which is
        # exactly the value the old, buggy residual_norm reported.
        phi_exact = 2 * nodes[:, 0] + 3 * nodes[:, 1] - nodes[:, 2]
        result = solve_laplace(
            nodes, elements, dirichlet_nodes, phi_exact[dirichlet_nodes], rtol=1e-12, maxiter=2000,
        )

        boundary_row_residual = np.max(np.abs(assemble_residual(nodes, elements, result["phi"])[boundary]))
        assert boundary_row_residual > 1e-3, (
            "test assumption broken: boundary rows should have an O(1) flux "
            "imbalance for this setup, or this test isn't exercising the bug"
        )

        assert result["residual_norm"] < 1e-8, (
            f"residual_norm={result['residual_norm']:.3e} looks like it "
            f"includes the O(1) boundary-row imbalance ({boundary_row_residual:.3e}) "
            f"again, not just the free-dof residual"
        )

        free_dof_residual = np.max(np.abs(assemble_residual(nodes, elements, result["phi"])[free]))
        assert result["residual_norm"] == pytest.approx(free_dof_residual, rel=1e-9)

    def test_residual_norm_nets_out_body_source(self):
        """With a nonzero body_source_rhs (MMS use case), residual_norm must
        be b - A@phi restricted to free dofs, not bare A@phi -- otherwise a
        converged MMS solve would misreport an O(source) residual."""
        nodes, elements = generate_structured_cube_mesh(n=6, L=1.0)
        boundary = cube_boundary_mask(nodes, L=1.0)
        dirichlet_nodes = np.where(boundary)[0]

        phi_exact = nodes[:, 0] ** 2
        source = 2.0 * np.ones(len(nodes))  # crude but nonzero forcing

        result = solve_laplace(
            nodes, elements, dirichlet_nodes, phi_exact[dirichlet_nodes],
            body_source_rhs=source, rtol=1e-12, maxiter=3000,
        )

        assert result["residual_norm"] < 1e-8


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
