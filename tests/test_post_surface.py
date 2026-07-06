"""
Unit/regression tests for pyfp3d.post.surface wall gradient recovery.

`wall_tangential_gradient_quadratic` was added while investigating the G1.6
accuracy saturation (see PROJECT_STRUCTURE.md "Known gaps"): it is a
genuine, well-conditioned improvement over the linear `wall_tangential_gradient`
(exact for a locally quadratic field, vs. only exact for a locally linear one),
but investigation showed the recovery step is a minor contributor to the
gate's total error, not the dominant one -- these tests lock in the recovery
operator's own accuracy so that fact doesn't silently regress or get
re-litigated without evidence.
"""

import numpy as np

from pyfp3d.post.surface import wall_tangential_gradient, wall_tangential_gradient_quadratic


def _flat_quadratic_patch(n=9, spacing=0.1):
    """A flat (z=0) structured-triangle grid with a known quadratic phi(x, y)."""
    xs = np.arange(n) * spacing
    ys = np.arange(n) * spacing
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    nodes = np.column_stack([xx.ravel(), yy.ravel(), np.zeros(n * n)])

    def node_id(i, j):
        return i * n + j

    faces = []
    for i in range(n - 1):
        for j in range(n - 1):
            a, b, c, d = node_id(i, j), node_id(i + 1, j), node_id(i + 1, j + 1), node_id(i, j + 1)
            faces.append([a, b, c])
            faces.append([a, c, d])
    faces = np.array(faces, dtype=np.int64)

    a0, a1, a2, a3, a4, a5 = 1.0, 2.0, 3.0, 0.5, 0.7, -0.3
    x, y = nodes[:, 0], nodes[:, 1]
    phi = a0 + a1 * x + a2 * y + a3 * x**2 + a4 * x * y + a5 * y**2
    grad_exact = np.column_stack([a1 + 2 * a3 * x + a4 * y, a2 + a4 * x + 2 * a5 * y, np.zeros(n * n)])
    return nodes, faces, phi, grad_exact


def test_quadratic_recovery_exact_on_flat_quadratic_field():
    """A local quadratic fit must reproduce a globally quadratic field's
    gradient to near machine precision, at every node (interior or edge)."""
    nodes, faces, phi, grad_exact = _flat_quadratic_patch()
    grad = wall_tangential_gradient_quadratic(nodes, faces, phi)

    touched = np.unique(faces)
    err = np.linalg.norm(grad[touched] - grad_exact[touched], axis=1)
    assert np.max(err) < 1e-8, f"max gradient error {np.max(err):.3e}, expected near machine-zero"


def test_quadratic_recovery_beats_linear_on_quadratic_field():
    """The linear (per-triangle) recovery is only exact for locally linear
    fields, so it must show a real (non-machine-zero) error on the same
    quadratic field that the quadratic recovery reproduces exactly."""
    nodes, faces, phi, grad_exact = _flat_quadratic_patch()
    touched = np.unique(faces)

    grad_lin = wall_tangential_gradient(nodes, faces, phi)
    err_lin = np.linalg.norm(grad_lin[touched] - grad_exact[touched], axis=1)

    grad_quad = wall_tangential_gradient_quadratic(nodes, faces, phi)
    err_quad = np.linalg.norm(grad_quad[touched] - grad_exact[touched], axis=1)

    assert err_lin.max() > 1e-3, "expected the linear recovery to show a real bias on a curved field"
    assert err_quad.max() < 1e-8
    assert err_quad.max() < err_lin.max()


def test_quadratic_recovery_matches_linear_on_linear_field():
    """Both recovery schemes must be exact for a globally linear field
    (the common baseline case)."""
    nodes, faces, _, _ = _flat_quadratic_patch()
    x, y = nodes[:, 0], nodes[:, 1]
    phi = 1.0 + 2.0 * x - 3.0 * y
    grad_exact = np.tile(np.array([2.0, -3.0, 0.0]), (len(nodes), 1))
    touched = np.unique(faces)

    grad_lin = wall_tangential_gradient(nodes, faces, phi)
    grad_quad = wall_tangential_gradient_quadratic(nodes, faces, phi)

    assert np.max(np.abs(grad_lin[touched] - grad_exact[touched])) < 1e-10
    assert np.max(np.abs(grad_quad[touched] - grad_exact[touched])) < 1e-8


def test_inconsistent_wall_winding_raises():
    """Mixed triangle winding makes area-weighted vertex normals cancel; the
    tangent planes built from them are then garbage. This used to proceed
    silently -- now `_wall_vertex_normals` raises. The consistent-winding
    twin of the same patch must keep working (and still recover a linear
    field exactly through the 2-parameter fallback for tiny patches)."""
    import pytest

    nodes = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], dtype=np.float64
    )
    phi = 2.0 * nodes[:, 0] - 3.0 * nodes[:, 1]

    faces_consistent = np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int64)  # both +z
    grad = wall_tangential_gradient_quadratic(nodes, faces_consistent, phi)
    touched = np.unique(faces_consistent)
    assert np.allclose(grad[touched], [2.0, -3.0, 0.0], atol=1e-9)

    faces_mixed = np.array([[0, 1, 2], [1, 2, 3]], dtype=np.int64)  # second is -z
    with pytest.raises(ValueError, match="winding"):
        wall_tangential_gradient_quadratic(nodes, faces_mixed, phi)


def test_quadratic_recovery_sphere_medium_mesh_regression(mesh_dir):
    """Locks in the ~20x recovery-only improvement measured on the G1.6
    sphere case during the accuracy-saturation investigation (see
    PROJECT_STRUCTURE.md): feeding the *exact* analytic potential
    (no FEM solve at all) through both recovery schemes, quadratic recovery's
    own bias should be a small fraction of the linear scheme's."""
    from pyfp3d.mesh.reader import read_mesh

    mesh = read_mesh(mesh_dir / "sphere_shell" / "medium.msh")
    nodes = mesh.nodes
    wall_faces = mesh.boundary_faces["wall"]
    wall_nodes = np.unique(wall_faces)

    r = np.linalg.norm(nodes, axis=1)
    phi_exact = nodes[:, 0] * (1.0 + 0.5 / r**3)  # a=1 sphere, exact potential

    def cp_error(grad_wall):
        q2 = np.sum(grad_wall[wall_nodes] ** 2, axis=1)
        cos_t = nodes[wall_nodes, 0] / r[wall_nodes]
        cp_exact = 1.0 - 2.25 * (1.0 - cos_t**2)
        return np.abs((1.0 - q2) - cp_exact)

    err_lin = cp_error(wall_tangential_gradient(nodes, wall_faces, phi_exact))
    err_quad = cp_error(wall_tangential_gradient_quadratic(nodes, wall_faces, phi_exact))

    assert err_quad.max() < 0.1 * err_lin.max(), (
        f"quadratic recovery-only error ({err_quad.max():.4f}) should be well under "
        f"1/10 of the linear scheme's ({err_lin.max():.4f}) on an exact input"
    )
