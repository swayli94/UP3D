"""
Option A true-normal weak-flux correction for curved walls on flat-facet
meshes (design.md §5.1; roadmap gates G1.3/G1.4/G1.5).

After integration by parts the Galerkin wall term <v, grad(phi).n_facet> is
dropped as a natural BC, which enforces zero flux through the *flat facet*
(unit outward normal n_facet) instead of through the *true curved surface*
(unit normal n_true at the closest point). Decomposing

    n_facet = (n_facet . n_true) n_true + t,    t = n_facet - (n_facet . n_true) n_true

the physical condition grad(phi) . n_true = 0 kills the first term, so the
wall term is <v, grad(phi) . t> with |t| = O(h). Moving it to the RHS gives
the corrected weak form (stiffness matrix completely unchanged, SPD, AMG
reusable):

    integral( grad(phi) . grad(v) ) dV = <v, grad(phi_k) . t>_wall

where grad(phi_k) is either the exact analytic gradient (oracle mode, gates
G1.3/G1.4) or the lagged previous-iterate gradient of the adjacent element
(G1.5). The assembled vector plugs straight into
solve_laplace(body_source_rhs=...).

Note t is invariant under n_true -> -n_true, so only n_facet's orientation
matters: it must be DOMAIN-outward (divergence-theorem convention), which is
enforced here from the adjacent tet, not from the surface-mesh winding.

Geometry enters only through a replaceable closest_point_normal(points)
callback returning (projected_points, true_unit_normals); analytic
implementations for the cylinder and the sphere are provided below. All
callback results are precomputed into SoA arrays outside any hot loop
(agent-rules hard rule #3); the assembly itself is vectorized numpy for the
oracle pre-study and gets the njit treatment when G1.5 makes it hot.

The V0 freestream gate never touches this path: the correction is opt-in via
an explicitly assembled RHS, and plain solve_laplace() is unchanged.
"""

from typing import Callable, Dict, Tuple

import numpy as np

# 3-point edge-midpoint quadrature on a triangle: exact for quadratics.
# Point q is the midpoint of edge (i, j); the P1 shape functions there are
# N_i = N_j = 1/2, N_k = 0.
_EDGE_MIDPOINT_PAIRS = ((0, 1), (1, 2), (2, 0))


def cylinder_closest_point_normal(
    points: np.ndarray, radius: float = 1.0, axis_center: Tuple[float, float] = (0.0, 0.0)
) -> Tuple[np.ndarray, np.ndarray]:
    """Closest-point projection + true unit normal for a circular cylinder
    with axis along z (closed form: p(x) = a (x, y, 0)/r_xy, n = (x, y, 0)/r_xy)."""
    p = np.asarray(points, dtype=np.float64)
    xy = p[:, :2] - np.asarray(axis_center, dtype=np.float64)
    r = np.linalg.norm(xy, axis=1)
    if np.any(r < 1e-12 * radius):
        raise ValueError("query point on the cylinder axis: projection undefined")
    n = np.zeros_like(p)
    n[:, :2] = xy / r[:, None]
    proj = p.copy()
    proj[:, :2] = np.asarray(axis_center, dtype=np.float64) + radius * n[:, :2]
    return proj, n


def sphere_closest_point_normal(
    points: np.ndarray, radius: float = 1.0, center: Tuple[float, float, float] = (0.0, 0.0, 0.0)
) -> Tuple[np.ndarray, np.ndarray]:
    """Closest-point projection + true unit normal for a sphere
    (n = (x - c)/|x - c|); the G1.4 oracle experiment's callback."""
    p = np.asarray(points, dtype=np.float64)
    d = p - np.asarray(center, dtype=np.float64)
    r = np.linalg.norm(d, axis=1)
    if np.any(r < 1e-12 * radius):
        raise ValueError("query point at the sphere center: projection undefined")
    n = d / r[:, None]
    return np.asarray(center, dtype=np.float64) + radius * n, n


def wall_face_adjacent_tets(elements: np.ndarray, wall_faces: np.ndarray) -> np.ndarray:
    """Index of the unique tet owning each wall facet.

    Boundary facets have exactly one owner (asserted). Plain-python dict
    lookup keyed by the sorted node triple -- preprocessing, not a hot loop.
    """
    wanted = {tuple(sorted(f)): i for i, f in enumerate(np.asarray(wall_faces))}
    owners = np.full(len(wall_faces), -1, dtype=np.int64)
    face_defs = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
    for e, tet in enumerate(np.asarray(elements)):
        for fd in face_defs:
            key = tuple(sorted((tet[fd[0]], tet[fd[1]], tet[fd[2]])))
            i = wanted.get(key)
            if i is not None:
                if owners[i] != -1:
                    raise ValueError(
                        f"wall facet {i} owned by two tets ({owners[i]}, {e}): "
                        "not a boundary facet"
                    )
                owners[i] = e
    if np.any(owners < 0):
        raise ValueError(f"{np.sum(owners < 0)} wall facets have no owning tet")
    return owners


def wall_correction_geometry(
    nodes: np.ndarray,
    elements: np.ndarray,
    wall_faces: np.ndarray,
    closest_point_normal: Callable[[np.ndarray], Tuple[np.ndarray, np.ndarray]],
) -> Dict[str, np.ndarray]:
    """Precompute the SoA geometric data of the correction RHS.

    Returns dict with:
        owners   (F,)    adjacent (owning) tet of each wall facet
        area     (F,)    facet area
        n_facet  (F, 3)  domain-outward unit facet normal (oriented from the
                         owning tet, independent of surface-mesh winding)
        qp       (F, 3, 3)  the 3 edge-midpoint quadrature points
        n_true   (F, 3, 3)  true unit normal at each quadrature point's
                            closest-point projection
        t        (F, 3, 3)  tangential defect n_facet - (n_facet.n_true) n_true
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    wall_faces = np.asarray(wall_faces)
    tri = nodes[wall_faces]  # (F, 3, 3)

    area_vec = 0.5 * np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = np.linalg.norm(area_vec, axis=1)
    if np.any(area < 1e-30):
        raise ValueError("degenerate wall facet (zero area)")
    n_facet = area_vec / area[:, None]

    # Orient domain-outward: away from the owning tet's off-facet vertex.
    owners = wall_face_adjacent_tets(elements, wall_faces)
    tet_centroids = nodes[np.asarray(elements)[owners]].mean(axis=1)
    outward = np.einsum("fk,fk->f", n_facet, tri.mean(axis=1) - tet_centroids)
    n_facet[outward < 0.0] *= -1.0

    qp = 0.5 * (tri[:, [0, 1, 2]] + tri[:, [1, 2, 0]])  # (F, 3, 3)
    _, n_true = closest_point_normal(qp.reshape(-1, 3))
    n_true = n_true.reshape(qp.shape)

    ndotn = np.einsum("fk,fqk->fq", n_facet, n_true)
    t = n_facet[:, None, :] - ndotn[:, :, None] * n_true

    return {
        "owners": owners,
        "area": area,
        "n_facet": n_facet,
        "qp": qp,
        "n_true": n_true,
        "t": t,
    }


def assemble_wall_flux_correction_rhs(
    n_nodes: int,
    wall_faces: np.ndarray,
    geometry: Dict[str, np.ndarray],
    grad_at_qp: np.ndarray,
) -> np.ndarray:
    """Assemble b_i = <N_i, grad(phi) . t>_wall into a full (n_nodes,) vector
    (the body_source_rhs argument of solve_laplace).

    grad_at_qp: (F, 3, 3) gradient of phi at each quadrature point -- the
    exact analytic gradient in oracle mode (G1.3/G1.4), or the owning
    element's piecewise-constant gradient broadcast over the facet (G1.5).
    """
    wall_faces = np.asarray(wall_faces)
    flux_defect = np.einsum("fqk,fqk->fq", grad_at_qp, geometry["t"])  # (F, 3)
    # weight per quadrature point = area/3; shape-function value 1/2 at the
    # two edge vertices, 0 at the opposite one
    contrib = (geometry["area"][:, None] / 3.0) * 0.5 * flux_defect
    rhs = np.zeros(n_nodes, dtype=np.float64)
    for q, (i, j) in enumerate(_EDGE_MIDPOINT_PAIRS):
        np.add.at(rhs, wall_faces[:, i], contrib[:, q])
        np.add.at(rhs, wall_faces[:, j], contrib[:, q])
    return rhs


def facet_normal_deviation_angles(geometry: Dict[str, np.ndarray]) -> np.ndarray:
    """Deviation angle (radians) between each facet normal and the true
    surface normal at the facet's quadrature points, |angle(n_facet, n_true)|
    folded to [0, pi/2] (t is sign-invariant, so only misalignment matters).
    Returns (F, 3). The geometric-error self-check of the G1.3 artifacts."""
    ndotn = np.einsum("fk,fqk->fq", geometry["n_facet"], geometry["n_true"])
    return np.arccos(np.clip(np.abs(ndotn), 0.0, 1.0))
