"""
Regression tests for the degenerate-element guard in `mesh.metrics.element_gradients`.

The old code silently returned *zero gradients* for any tet with
|det J| < 1e-20 -- an absolute threshold. Two failure modes, both locked in
here:

  1. A genuinely degenerate (coplanar) tet sailed straight into assembly and
     corrupted the stiffness matrix / residual with no error anywhere.
  2. The absolute threshold misfired the other way on perfectly well-shaped
     but *small* elements: a regular tet with edge ~1e-7 has det ~ 7e-22,
     below the old cutoff, so its (perfectly computable) gradients were
     zeroed too. The guard now scales with the element's own edge lengths
     (|det J| vs. 1e-12 x product of edge norms), so it is unit- and
     mesh-scale-invariant.
"""

import numpy as np
import pytest

from pyfp3d.kernels.residual import assemble_residual
from pyfp3d.mesh.metrics import element_gradients

from .mesh_utils import generate_structured_cube_mesh


class TestDegenerateElementGuard:
    def test_coplanar_tet_raises(self):
        nodes = np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], dtype=np.float64
        )  # all four nodes in the z = 0 plane
        elements = np.array([[0, 1, 2, 3]], dtype=np.int32)

        with pytest.raises(ValueError):
            element_gradients(nodes, elements, 0)

    def test_tiny_well_shaped_tet_keeps_exact_gradients(self):
        """A unit-shaped tet scaled to edge ~1e-7 (det ~ 1e-21, below the old
        absolute cutoff) must still recover a linear field's gradient exactly
        instead of being silently zeroed."""
        scale = 1e-7
        nodes = scale * np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64
        )
        elements = np.array([[0, 1, 2, 3]], dtype=np.int32)

        grads = element_gradients(nodes, elements, 0)
        a = np.array([2.0, -3.0, 0.5])
        phi = nodes @ a
        grad_phi = phi[elements[0]] @ grads

        assert np.allclose(grad_phi, a, rtol=1e-9), (
            f"tiny well-shaped tet: recovered gradient {grad_phi}, expected {a} "
            "(the old absolute threshold zeroed this element's gradients)"
        )

    def test_degenerate_tet_raises_through_assembly(self):
        """The guard must fire through the njit'd assembly call chain, not
        just on a direct call -- assembly is where silent corruption
        actually happened."""
        nodes, elements = generate_structured_cube_mesh(n=2, L=1.0)
        # Four corners of the z=0 face: a zero-volume "tet" appended to an
        # otherwise valid mesh.
        coplanar = np.array([[0, 18, 6, 24]], dtype=np.int32)
        elements_bad = np.vstack([elements, coplanar]).astype(np.int32)
        assert abs(np.linalg.det(
            (nodes[coplanar[0, 1:]] - nodes[coplanar[0, 0]]).T
        )) < 1e-14, "fixture assumption broken: appended tet should be coplanar"

        phi = nodes[:, 0].copy()
        with pytest.raises(ValueError):
            assemble_residual(nodes, elements_bad, phi)


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
