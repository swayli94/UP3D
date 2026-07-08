"""P6 differentiable artificial density: correctness of the smooth
inflow-weighted upwind operator (design.md §3.2, Eq. 3.4).

Unit-level invariants on a structured cube -- no transonic solve. The
transonic shock-ladder / smoothness behaviour is validated separately by the
G4.1/G6.1 coarse checks.
"""

import numpy as np

from pyfp3d.kernels.jacobian import PicardOperator
from pyfp3d.kernels.upwind import UpwindOperator, compute_face_normals
from pyfp3d.physics.isentropic import density_field, mach_squared_field

from .mesh_utils import cube_boundary_mask, generate_structured_cube_mesh


def _cube():
    return generate_structured_cube_mesh(n=6, L=1.0)


def test_face_normals_unit_and_outward():
    nodes, elements = _cube()
    fn = np.empty((len(elements), 4, 3), dtype=np.float64)
    compute_face_normals(np.ascontiguousarray(nodes, dtype=np.float64),
                         np.ascontiguousarray(elements), fn)
    # unit length
    lens = np.linalg.norm(fn, axis=2)
    assert np.allclose(lens, 1.0, atol=1e-12)
    # outward: normal . (face_centroid - element_centroid) > 0
    fd = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])
    ec = nodes[elements].mean(axis=1)
    for f in range(4):
        face_c = nodes[elements[:, fd[f]]].mean(axis=1)
        out = np.einsum("ij,ij->i", fn[:, f, :], face_c - ec)
        assert np.all(out > 0.0)


def test_interior_faces_have_antiparallel_normals():
    """Two tets sharing a face see opposite outward normals across it."""
    nodes, elements = _cube()
    upw = UpwindOperator(nodes, elements)
    fn = np.empty((len(elements), 4, 3), dtype=np.float64)
    compute_face_normals(np.ascontiguousarray(nodes, dtype=np.float64),
                         np.ascontiguousarray(elements), fn)
    fnb = upw.face_neighbors
    n_checked = 0
    for e in range(len(elements)):
        for f in range(4):
            nb = fnb[e, f]
            if nb < 0:
                continue
            # find the shared face on nb (the local face whose neighbor is e)
            for g in range(4):
                if fnb[nb, g] == e:
                    dot = float(np.dot(fn[e, f], fn[nb, g]))
                    assert dot < -0.999   # anti-parallel
                    n_checked += 1
                    break
    assert n_checked > 0


def _kernel_op(nodes, elements):
    """The shipped P6 streamline-kernel operator (opt-in; the module default
    is the P4 walk until G6.1/G6.2 close)."""
    return UpwindOperator(nodes, elements, weighted=True, mode="kernel")


def test_weighted_noop_subcritical_bitwise():
    nodes, elements = _cube()
    upw = _kernel_op(nodes, elements)
    op = PicardOperator(nodes, elements)
    phi = nodes[:, 0].copy()                        # +x, q^2 = 1
    grad, q2 = op.velocities(phi)
    rho = density_field(q2, 0.5)
    rho_t = upw.rho_tilde(grad, q2, rho, 0.5, 1.5, 0.95)
    assert upw.nu_max == 0.0
    assert np.array_equal(rho_t, rho)               # bit-exact no-op


def test_uniform_density_is_invariant_even_supersonic():
    """rho_up is a convex blend of neighbour rho, so a spatially uniform rho
    gives rho_tilde == rho even where nu > 0 (rho_e - rho_up == 0)."""
    nodes, elements = _cube()
    upw = _kernel_op(nodes, elements)
    op = PicardOperator(nodes, elements)
    phi = 1.3 * nodes[:, 0].copy()                  # q^2 = 1.69 -> supersonic
    grad, q2 = op.velocities(phi)
    m2 = mach_squared_field(q2, 0.8)
    assert m2.max() > 0.95 ** 2                     # genuinely supersonic
    rho = np.full(op.n_tets, 0.7)                   # uniform density
    rho_t = upw.rho_tilde(grad, q2, rho, 0.8, 1.5, 0.95)
    assert upw.nu_max > 0.0                          # switch active
    assert np.allclose(rho_t, 0.7, atol=1e-12)      # yet no dissipation


def test_weighted_sweep_finite_and_deterministic():
    nodes, elements = _cube()
    upw = _kernel_op(nodes, elements)
    op = PicardOperator(nodes, elements)
    phi = 1.2 * nodes[:, 0] + 0.15 * nodes[:, 0] ** 2   # varying supersonic
    grad, q2 = op.velocities(phi)
    rho = density_field(np.minimum(q2, 3.0), 0.85)
    a = upw.rho_tilde(grad, q2, rho, 0.85, 1.5, 0.95).copy()
    b = upw.rho_tilde(grad, q2, rho, 0.85, 1.5, 0.95).copy()
    assert np.all(np.isfinite(a))
    assert np.array_equal(a, b)                     # prange-deterministic
    # dissipation only reduces density toward upstream, never below the floor
    assert a.min() >= 0.05 - 1e-15


def test_weighted_false_restores_walk_operator():
    """weighted=False must reproduce the P4 integer-walk rho_tilde exactly."""
    from pyfp3d.kernels.upwind import rho_tilde_sweep, upstream_elements

    nodes, elements = _cube()
    op = PicardOperator(nodes, elements)
    phi = 1.25 * nodes[:, 0].copy()
    grad, q2 = op.velocities(phi)
    rho = density_field(q2, 0.82)

    upw = UpwindOperator(nodes, elements, weighted=False)
    rho_t = upw.rho_tilde(grad, q2, rho, 0.82, 1.5, 0.95).copy()

    u = np.empty(op.n_tets, dtype=np.int64)
    upstream_elements(upw.face_neighbors, upw.centroids, upw._nodes,
                      upw._elements, grad, u)
    nu = np.empty(op.n_tets)
    ref = np.empty(op.n_tets)
    rho_tilde_sweep(q2, rho, u, 0.82, 0.95, 1.5, 1.4, 0.05, nu, ref)
    assert np.array_equal(rho_t, ref)
