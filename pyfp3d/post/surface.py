"""
Nodal post-processing derived from the element-constant P1 gradient.

Reference: design.md Sec 9 ("Nodal q^2, M, Cp via volume-weighted element
averages (or superconvergent patch recovery later)").
"""

import numpy as np

from pyfp3d.mesh.metrics import compute_tet_volumes, element_gradients


def _element_gradients_and_centroids(nodes: np.ndarray, elements: np.ndarray, phi: np.ndarray):
    n_tets = len(elements)
    grad_elem = np.empty((n_tets, 3), dtype=np.float64)
    for e in range(n_tets):
        grads = element_gradients(nodes, elements, e)
        grad_elem[e] = phi[elements[e]] @ grads
    centroids = nodes[elements].mean(axis=1)
    return grad_elem, centroids


def _node_to_element_patches(elements: np.ndarray, n_nodes: int):
    """CSR-style (offsets, element_ids) map from node -> incident elements."""
    n_tets = len(elements)
    flat_nodes = elements.reshape(-1)
    flat_elems = np.repeat(np.arange(n_tets), 4)

    order = np.argsort(flat_nodes, kind="stable")
    sorted_nodes = flat_nodes[order]
    sorted_elems = flat_elems[order]

    counts = np.bincount(sorted_nodes, minlength=n_nodes)
    offsets = np.zeros(n_nodes + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    return offsets, sorted_elems


def nodal_gradient_recovery(nodes: np.ndarray, elements: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """
    Recover a nodal gradient field from the element-constant P1 gradient by
    a volume-weighted average over each node's incident elements.

    This is exact for a globally linear field (all incident elements share
    the same constant gradient) but has an O(h) bias at domain boundaries:
    every incident tet extends *inward* from the node, so the average
    effectively samples the gradient somewhere in the tets' interior rather
    than at the node itself, and any nonuniformity across the patch (e.g.
    the tangential-velocity falloff away from a curved wall) is smoothed
    into an underestimate right at the boundary. `nodal_gradient_recovery_spr`
    corrects that bias for surface quantities like Cp; keep this function for
    interior-field use where the plain average is markedly cheaper.

    Args:
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) tetrahedral connectivity
        phi: (n_nodes,) nodal potential

    Returns:
        grad_nodal: (n_nodes, 3) recovered nodal gradient
    """
    n_nodes = len(nodes)
    volumes = compute_tet_volumes(nodes, elements)
    grad_elem, _ = _element_gradients_and_centroids(nodes, elements, phi)

    grad_sum = np.zeros((n_nodes, 3), dtype=np.float64)
    weight_sum = np.zeros(n_nodes, dtype=np.float64)
    for i in range(4):
        node_ids = elements[:, i]
        np.add.at(grad_sum, node_ids, volumes[:, None] * grad_elem)
        np.add.at(weight_sum, node_ids, volumes)

    return grad_sum / weight_sum[:, None]


def wall_tangential_gradient(nodes: np.ndarray, wall_faces: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """
    Surface velocity at a solid wall, computed directly from the wall's own
    triangulation instead of extrapolating the volume-element gradient.

    Rationale (see design.md Sec 5): the natural BC at a solid wall enforces
    rho * dphi/dn = 0, so the full 3D velocity there is exactly the *surface*
    (tangential) gradient of phi -- no normal component to reconstruct. This
    computes the standard constant-per-triangle P1 gradient in each wall
    triangle's own local tangent plane (project the triangle into 2D, solve
    the 2x2 linear system for [dphi/du, dphi/dv], map back to 3D), then does
    an area-weighted nodal average over each node's incident wall triangles.

    An earlier volume-based approach (averaging/extrapolating the tets'
    element-constant gradient at the wall) was tried and rejected: every
    incident tet sits to one side of the wall, and the tangential velocity
    physically decays moving away from the wall into the fluid, so any
    volume-based estimate is systematically biased low right at the surface.
    Naive local-least-squares extrapolation to correct that bias is also not
    a good fix -- on a real mesh, some nodes' 1-ring element patches are
    nearly coplanar, which makes the extrapolation's linear term ill-
    conditioned and can blow the recovered gradient up by orders of
    magnitude (observed: a single node's Cp error of 429 on an otherwise
    well-behaved mesh). Working entirely within the wall's own 2-manifold
    avoids that failure mode, since it only ever interpolates/averages
    (never extrapolates along a near-degenerate direction).

    Args:
        nodes: (n_nodes, 3) full mesh nodal coordinates
        wall_faces: (n_wall_tris, 3) triangle connectivity on the wall surface
        phi: (n_nodes,) nodal potential

    Returns:
        grad_wall: (n_nodes, 3) surface-tangential gradient, populated only
            at nodes touched by `wall_faces` (NaN elsewhere)
    """
    p0 = nodes[wall_faces[:, 0]]
    p1 = nodes[wall_faces[:, 1]]
    p2 = nodes[wall_faces[:, 2]]

    edge1 = p1 - p0
    area_vec = np.cross(edge1, p2 - p0)
    twice_area = np.linalg.norm(area_vec, axis=1)
    tri_normal = area_vec / twice_area[:, None]
    e1 = edge1 / np.linalg.norm(edge1, axis=1)[:, None]
    e2 = np.cross(tri_normal, e1)

    # Project the 3 vertices into each triangle's own local 2D (u, v) frame
    # (p0 is the local origin) and solve the linear-triangle gradient system
    # [[u1, v1], [u2, v2]] . [du, dv] = [f1 - f0, f2 - f0].
    u1, v1 = np.sum(edge1 * e1, axis=1), np.sum(edge1 * e2, axis=1)
    edge2 = p2 - p0
    u2, v2 = np.sum(edge2 * e1, axis=1), np.sum(edge2 * e2, axis=1)
    df1 = phi[wall_faces[:, 1]] - phi[wall_faces[:, 0]]
    df2 = phi[wall_faces[:, 2]] - phi[wall_faces[:, 0]]

    det = u1 * v2 - v1 * u2
    du = (v2 * df1 - v1 * df2) / det
    dv = (-u2 * df1 + u1 * df2) / det
    grad_tri = du[:, None] * e1 + dv[:, None] * e2

    area = 0.5 * twice_area
    n_nodes = len(nodes)
    grad_sum = np.zeros((n_nodes, 3), dtype=np.float64)
    weight_sum = np.zeros(n_nodes, dtype=np.float64)
    for i in range(3):
        np.add.at(grad_sum, wall_faces[:, i], area[:, None] * grad_tri)
        np.add.at(weight_sum, wall_faces[:, i], area)

    grad_wall = np.full((n_nodes, 3), np.nan, dtype=np.float64)
    touched = weight_sum > 0
    grad_wall[touched] = grad_sum[touched] / weight_sum[touched, None]
    return grad_wall
