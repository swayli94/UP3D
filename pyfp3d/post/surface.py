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
    into an underestimate right at the boundary. For wall quantities like Cp
    use `wall_tangential_gradient` / `wall_tangential_gradient_quadratic`,
    which avoid that bias by working in the wall's own 2-manifold; keep this
    function for interior-field use where the plain average is markedly
    cheaper.

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


def _wall_vertex_normals(nodes: np.ndarray, wall_faces: np.ndarray) -> np.ndarray:
    """Area-weighted vertex-normal field over the wall's own triangulation.

    Raises:
        ValueError: if any vertex's incident triangle normals nearly cancel
            (|sum of area-weighted normals| / sum of areas < 0.05). The
            averaged normal is then meaningless, and every tangent plane
            built from it downstream would be garbage. Two known causes:
            inconsistent triangle winding in the wall surface (a mesh/tagging
            defect -- fix the mesh), or a genuinely sharp crease such as a
            thin trailing edge (wedge angle below ~6 deg), where a vertex
            normal is ill-defined and the recovery scheme needs an explicit
            crease treatment instead of silent averaging.
    """
    p0, p1, p2 = nodes[wall_faces[:, 0]], nodes[wall_faces[:, 1]], nodes[wall_faces[:, 2]]
    area_vec = np.cross(p1 - p0, p2 - p0)
    twice_area = np.linalg.norm(area_vec, axis=1)
    tri_normal = area_vec / twice_area[:, None]
    area = 0.5 * twice_area

    n_nodes = len(nodes)
    normal_sum = np.zeros((n_nodes, 3), dtype=np.float64)
    weight_sum = np.zeros(n_nodes, dtype=np.float64)
    for i in range(3):
        np.add.at(normal_sum, wall_faces[:, i], area[:, None] * tri_normal)
        np.add.at(weight_sum, wall_faces[:, i], area)

    norms = np.linalg.norm(normal_sum, axis=1)
    touched = weight_sum > 0
    consistency = norms[touched] / weight_sum[touched]
    if np.any(consistency < 0.05):
        raise ValueError(
            "wall vertex normals nearly cancel at "
            f"{int(np.sum(consistency < 0.05))} node(s) "
            "(inconsistent wall-triangle winding, or a sharp crease like a "
            "thin trailing edge) -- see _wall_vertex_normals docstring"
        )

    vertex_normal = np.full((n_nodes, 3), np.nan, dtype=np.float64)
    touched_n = norms > 0
    vertex_normal[touched_n] = normal_sum[touched_n] / norms[touched_n, None]
    return vertex_normal


def _wall_node_rings(wall_faces: np.ndarray, n_nodes: int) -> list:
    """1-ring neighbor sets (via wall-triangle edges) per node, wall nodes only."""
    rings = [None] * n_nodes
    for tri in wall_faces:
        a, b, c = int(tri[0]), int(tri[1]), int(tri[2])
        for node, others in ((a, (b, c)), (b, (a, c)), (c, (a, b))):
            if rings[node] is None:
                rings[node] = set()
            rings[node].update(others)
    return rings


def wall_tangential_gradient_quadratic(
    nodes: np.ndarray, wall_faces: np.ndarray, phi: np.ndarray
) -> np.ndarray:
    """
    Quadratic (superconvergent-style) tangential-gradient recovery on a
    curved wall surface: at each wall node, fit a local quadratic model
    dphi(u, v) = c1*u + c2*v + c3*u^2 + c4*u*v + c5*v^2 (u, v in the node's
    own tangent plane, centered at the node so the constant term is zero by
    construction) through its 1-ring patch by least squares, then the
    recovered tangential gradient is just (c1, c2) mapped back to 3D -- no
    evaluation/extrapolation step needed since the fit is already centered
    at the node.

    This is exact (zero recovery bias) for any quadratic tangential field,
    vs. `wall_tangential_gradient`'s per-triangle *linear* fit + area-
    weighted average, which is only exact for a locally linear field. Empirically
    (see PROJECT_STRUCTURE.md "Known gaps"), this cuts the recovery
    operator's own bias by roughly 20x on an exact-input test (e.g. medium
    mesh: ~2.8% down to ~0.17% max Cp error with no FEM solve at all, phi
    sampled directly from the analytic solution) -- but only trims a few
    percent off the *total* gate-G1.6 pipeline error (e.g. ~12.0% down to
    ~11.6% on the medium mesh), because the recovery step was never the
    dominant error source: the dominant source is the volume PDE solve's own
    accuracy next to the curved-but-flat-faceted wall (a geometric/
    variational-crime consistency error from imposing the natural BC on the
    polyhedral wall approximation instead of the true curved surface), which
    this function cannot fix since it only ever sees the already-computed
    nodal phi values, not the operator that produced them.

    Conditioning: fit via `np.linalg.lstsq` (SVD-based minimum-norm solution,
    not a normal-equations solve or explicit matrix inverse), with an
    explicit rank check -- if the 1-ring patch is rank-deficient for the
    6-parameter quadratic (fewer than 5 independent tangential directions,
    e.g. a low-valence node), the patch is expanded to the 2-ring before
    refitting; if that is still degenerate, this falls back to a 2-parameter
    *linear* (planar) fit in the same tangent plane. Because every fit stays
    within the node's own well-conditioned 2D tangent-plane parametrization
    (never extrapolating along a near-null 3D direction) and uses a
    minimum-norm solver, this degrades gracefully instead of the catastrophic
    blowup (one node's Cp error of 429) hit by the earlier fully-3D
    volume-gradient least-squares extrapolation attempt -- see
    `wall_tangential_gradient`'s docstring for that history.

    Args:
        nodes: (n_nodes, 3) full mesh nodal coordinates
        wall_faces: (n_wall_tris, 3) triangle connectivity on the wall surface
        phi: (n_nodes,) nodal potential

    Returns:
        grad_wall: (n_nodes, 3) surface-tangential gradient, populated only
            at nodes touched by `wall_faces` (NaN elsewhere)
    """
    n_nodes = len(nodes)
    vertex_normal = _wall_vertex_normals(nodes, wall_faces)
    rings = _wall_node_rings(wall_faces, n_nodes)
    wall_nodes = np.unique(wall_faces)

    grad_wall = np.full((n_nodes, 3), np.nan, dtype=np.float64)
    for i in wall_nodes:
        n = vertex_normal[i]
        seed = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        e1 = seed - np.dot(seed, n) * n
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(n, e1)

        patch = np.fromiter(rings[i], dtype=np.int64)
        if len(patch) < 6:
            expanded = set(rings[i])
            for j in rings[i]:
                expanded.update(rings[j])
            expanded.discard(i)
            patch = np.fromiter(expanded, dtype=np.int64)

        d = nodes[patch] - nodes[i]
        u = d @ e1
        v = d @ e2
        dphi = phi[patch] - phi[i]

        M = np.column_stack([u, v, u * u, u * v, v * v])
        sol, _, rank, _ = np.linalg.lstsq(M, dphi, rcond=1e-6)
        if rank < 5:
            sol_lin, *_ = np.linalg.lstsq(M[:, :2], dphi, rcond=1e-6)
            c1, c2 = sol_lin
        else:
            c1, c2 = sol[0], sol[1]
        grad_wall[i] = c1 * e1 + c2 * e2

    return grad_wall
