"""
Nodal post-processing derived from the element-constant P1 gradient.

Reference: design.md Sec 9 ("Nodal q^2, M, Cp via volume-weighted element
averages (or superconvergent patch recovery later)").
"""

from typing import Tuple

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


def triangle_tangential_gradients(
    nodes: np.ndarray, wall_faces: np.ndarray, phi: np.ndarray
):
    """Constant P1 gradient of phi in each wall triangle's own tangent plane.

    At a solid wall the natural BC gives dphi/dn = 0, so this in-plane
    gradient IS the full wall velocity (see wall_tangential_gradient's
    docstring). Kept triangle-wise (no nodal averaging) so consumers like
    the force integral and the sectional-Cp curve stay well-defined across
    sharp creases (thin TE), where vertex averaging would mix the two
    sides' tangent planes.

    Returns:
        (grad_tri, area, tri_normal): (n_tris, 3) tangential gradient,
        (n_tris,) triangle areas, (n_tris, 3) unit normals of the STORED
        winding (no outward orientation implied).
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

    u1, v1 = np.sum(edge1 * e1, axis=1), np.sum(edge1 * e2, axis=1)
    edge2 = p2 - p0
    u2, v2 = np.sum(edge2 * e1, axis=1), np.sum(edge2 * e2, axis=1)
    df1 = phi[wall_faces[:, 1]] - phi[wall_faces[:, 0]]
    df2 = phi[wall_faces[:, 2]] - phi[wall_faces[:, 0]]

    det = u1 * v2 - v1 * u2
    du = (v2 * df1 - v1 * df2) / det
    dv = (-u2 * df1 + u1 * df2) / det
    grad_tri = du[:, None] * e1 + dv[:, None] * e2
    return grad_tri, 0.5 * twice_area, tri_normal


def wall_triangle_adjacency(wall_faces: np.ndarray) -> np.ndarray:
    """Edge-neighbour triangle of each wall triangle: (n_tris, 3) int, entry
    (t, k) is the triangle sharing edge k of t (edges (0,1),(1,2),(2,0)), or -1
    if that edge is a surface boundary. One-time geometric precompute for the
    P6 wall-gradient smoothing. A sharp trailing edge shares its TE edge between
    the upper and lower triangle -- still manifold (2 tris/edge), and the
    smoothing's normal gate (below) is what stops averaging across that fold."""
    from collections import defaultdict

    wf = np.asarray(wall_faces, dtype=np.int64)
    n = len(wf)
    adj = np.full((n, 3), -1, dtype=np.int64)
    emap = defaultdict(list)
    tri_edges = ((0, 1), (1, 2), (2, 0))
    for t in range(n):
        for k, (a, b) in enumerate(tri_edges):
            va, vb = wf[t, a], wf[t, b]
            emap[(min(va, vb), max(va, vb))].append((t, k))
    for occ in emap.values():
        if len(occ) == 2:
            (t0, k0), (t1, k1) = occ
            adj[t0, k0] = t1
            adj[t1, k1] = t0
    return adj


def smooth_wall_tangential_gradients(
    grad_tri: np.ndarray, out_normal: np.ndarray, area: np.ndarray,
    adj: np.ndarray, n_passes: int = 2, cos_thresh: float = 0.2,
) -> np.ndarray:
    """Area-weighted edge-neighbour smoothing of the per-triangle wall gradient,
    gated so it never averages across a sharp crease (roadmap P6 / gate G6.1).

    The per-triangle constant P1 gradient of a SMOOTH phi is O(1)-noisy on the
    sliver prism-split wall triangulation -- adjacent triangles oscillate,
    which the no-averaging wall-Cp extractor shows as a ~2-cell surface-Cp
    SAWTOOTH (root-caused 2026-07-08: nodal averaging on the same field drops
    the G6.1 metric ~330x, so the sawtooth is a RECOVERY artifact, not the
    artificial-density flux). A few Jacobi passes averaging each triangle with
    its edge neighbours removes it. The gate `out_normal[t] . out_normal[nb] >
    cos_thresh` (outward-oriented normals) only averages neighbours on the same
    smooth surface: across the sharp TE the two triangles' outward normals are
    nearly anti-parallel (dot ~ -1 < cos_thresh) so the crease is preserved --
    the exact reason `triangle_tangential_gradients` was kept per-triangle.

    Linear in `grad_tri` (a fixed smoothing operator at fixed geometry), so it
    is differentiable and reduces to identity where a triangle has no same-side
    neighbour. `n_passes=0` returns the input unchanged (bit-identical).
    """
    if n_passes <= 0:
        return grad_tri
    g = np.asarray(grad_tri, dtype=np.float64)
    n = len(g)
    for _ in range(n_passes):
        acc = g * area[:, None]
        wsum = area.copy()
        for k in range(3):
            nb = adj[:, k]
            valid = nb >= 0
            aligned = np.zeros(n, dtype=bool)
            aligned[valid] = (
                np.sum(out_normal[valid] * out_normal[nb[valid]], axis=1)
                > cos_thresh
            )
            m = valid & aligned
            acc[m] += g[nb[m]] * area[nb[m], None]
            wsum[m] += area[nb[m]]
        g = acc / wsum[:, None]
    return g


def wall_outward_normals(
    nodes: np.ndarray, elements: np.ndarray, wall_faces: np.ndarray
) -> np.ndarray:
    """Unit wall-triangle normals oriented OUT of the body (into the fluid).

    Orientation is taken from each face's owning tet (the fluid is where
    the tet is), so it needs no assumption about the stored winding or the
    body's shape.
    """
    el = np.asarray(elements, dtype=np.int64)
    tet_faces = np.concatenate(
        [el[:, [1, 2, 3]], el[:, [0, 2, 3]], el[:, [0, 1, 3]], el[:, [0, 1, 2]]]
    )
    tet_faces.sort(axis=1)
    tet_faces = np.ascontiguousarray(tet_faces)
    keys = tet_faces.view([("", tet_faces.dtype)] * 3).ravel()
    owners = np.tile(np.arange(len(el), dtype=np.int64), 4)
    order = np.argsort(keys)
    keys_sorted, owners_sorted = keys[order], owners[order]

    wf = np.ascontiguousarray(np.sort(np.asarray(wall_faces, dtype=np.int64), axis=1))
    wf_keys = wf.view([("", wf.dtype)] * 3).ravel()
    lo = np.searchsorted(keys_sorted, wf_keys, side="left")
    hi = np.searchsorted(keys_sorted, wf_keys, side="right")
    if np.any(hi - lo != 1):
        raise ValueError(
            f"{int(np.sum(hi - lo != 1))} wall face(s) not owned by exactly "
            "one tet"
        )
    owner = owners_sorted[lo]

    p0 = nodes[wall_faces[:, 0]]
    area_vec = np.cross(nodes[wall_faces[:, 1]] - p0, nodes[wall_faces[:, 2]] - p0)
    normal = area_vec / np.linalg.norm(area_vec, axis=1)[:, None]
    face_center = nodes[wall_faces].mean(axis=1)
    owner_centroid = nodes[np.asarray(elements)[owner]].mean(axis=1)
    # Fluid side = owner side; body-outward means pointing toward the fluid.
    flip = np.sum(normal * (owner_centroid - face_center), axis=1) < 0
    normal[flip] *= -1.0
    return normal


def wall_crease_angles(
    nodes: np.ndarray, elements: np.ndarray, wall_faces: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Turning angle of the wall surface across each interior edge of its
    triangulation: ``(angle_deg, edge_midpoint)``, one entry per edge shared by
    two wall triangles.

    The angle between the two OUTWARD normals, so it reads the geometry, not
    the winding. On a smooth surface it is O(h * curvature) and vanishes under
    refinement; on a sharp edge it converges to the edge's turning angle and
    stays there -- which is the discrete signature Track M's M5 gate uses to
    tell a rounded tip cap from a flat one, and which the sharp TE (~180 deg,
    by design -- it carries the Kutta condition) shows too.
    """
    normals = wall_outward_normals(nodes, elements, wall_faces)
    adj = wall_triangle_adjacency(wall_faces)
    wf = np.asarray(wall_faces, dtype=np.int64)

    tri = np.repeat(np.arange(len(wf), dtype=np.int64), 3)
    slot = np.tile(np.arange(3, dtype=np.int64), len(wf))
    nb = adj.ravel()
    # Interior edges only, and each edge once.
    keep = (nb >= 0) & (tri < nb)
    tri, slot, nb = tri[keep], slot[keep], nb[keep]

    cos = np.einsum("ij,ij->i", normals[tri], normals[nb])
    angle = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))
    edge = np.stack([wf[tri, slot], wf[tri, (slot + 1) % 3]], axis=1)
    return angle, nodes[edge].mean(axis=1)


def wall_force_coefficients(
    nodes: np.ndarray,
    elements: np.ndarray,
    wall_faces: np.ndarray,
    phi: np.ndarray,
    alpha_deg: float = 0.0,
    u_inf: float = 1.0,
    s_ref: float = 1.0,
    m_inf: float = 0.0,
    smooth_passes: int = 0,
) -> dict:
    """Pressure-integrated force coefficients on the wall (design.md Sec 9).

    Per wall triangle: velocity = in-plane tangential gradient (exactly the
    wall velocity under the natural BC), Cp from Bernoulli (m_inf = 0, P2)
    or the exact isentropic law (2.5) (m_inf > 0, P3), force
    dC_F = -Cp n_out dA / S_ref  with n_out the body-outward normal.

    `smooth_passes` > 0 applies the normal-gated edge-neighbour gradient
    smoothing (`smooth_wall_tangential_gradients`) before Cp; `0` is
    bit-identical to the original per-triangle integration. NOTE (G6.3, measured
    2026-07-08): smoothing is for the reported *Cp curve*, NOT the loads — on
    ONERA M6 coarse it moves CL_p ~1% further below the trustworthy CL_KJ
    (V6 2.40%→3.35%) because the ±sawtooth cancels in the integral and the
    averaging instead smears the LE suction peak. Keep `smooth_passes=0` for
    forces; the param stays opt-in for experiments only.

    Args:
        nodes, elements: CUT-mesh arrays (used to orient normals)
        wall_faces: (n, 3) wall triangles
        phi: (n_nodes,) potential on the cut mesh
        alpha_deg: incidence; lift is measured normal to the freestream
        u_inf: freestream speed
        s_ref: reference area (chord x span for the quasi-2D cases)
        m_inf: freestream Mach; 0.0 selects the incompressible Cp

    Returns:
        dict: cl, cd_pressure, cf (3-vector), cp_tri (per-triangle Cp)
    """
    from pyfp3d.physics.isentropic import (
        pressure_coefficient,
        pressure_coefficient_incompressible,
    )

    grad_tri, area, _ = triangle_tangential_gradients(nodes, wall_faces, phi)
    n_out = wall_outward_normals(nodes, elements, wall_faces)
    if smooth_passes > 0:
        adj = wall_triangle_adjacency(wall_faces)
        grad_tri = smooth_wall_tangential_gradients(
            grad_tri, n_out, area, adj, n_passes=smooth_passes)

    q2 = np.sum(grad_tri * grad_tri, axis=1) / u_inf**2
    if m_inf > 0.0:
        cp_tri = np.array([pressure_coefficient(q, m_inf) for q in q2])
    else:
        cp_tri = np.array([pressure_coefficient_incompressible(q) for q in q2])

    cf = -(cp_tri * area) @ n_out / s_ref
    a = np.deg2rad(alpha_deg)
    lift_dir = np.array([-np.sin(a), np.cos(a), 0.0])
    drag_dir = np.array([np.cos(a), np.sin(a), 0.0])
    return {
        "cl": float(cf @ lift_dir),
        "cd_pressure": float(cf @ drag_dir),
        "cf": cf,
        "cp_tri": cp_tri,
    }


def sectional_cl_from_gamma(
    gamma: np.ndarray, chord: float = 1.0, u_inf: float = 1.0
) -> np.ndarray:
    """Kutta-Joukowski sectional lift, cl_j = 2 Gamma_j / (U_inf c)
    (design.md Sec 9) -- the cross-check against pressure integration
    that gate G2.4 holds to < 1%. `chord` may be a scalar or a per-station
    array (broadcasts), so cl(eta) follows the taper on a 3D wing."""
    return 2.0 * np.atleast_1d(np.asarray(gamma, dtype=np.float64)) / (u_inf * chord)


def planform_area(
    nodes: np.ndarray,
    wall_faces: np.ndarray,
    thickness_axis: int = 1,
) -> float:
    """Half-wing reference area: the wall projected onto the chord-span plane.

    Both the upper and lower surfaces project onto the same planform, so the
    projected areas add; summing the absolute thickness-axis (lift-axis)
    component of every wall-face area vector and halving recovers the single
    planform once. The flat tip cap (normal ~ +span) has a near-zero
    thickness component and drops out automatically; the symmetry plane
    carries no wall faces, so this is exactly the half-model planform.

    Pairs with `wall_force_coefficients(s_ref=...)` and `cl_kj_3d(s_ref=...)`:
    the half-planform paired with the half-model pressure/circulation lift
    yields the true wing CL (no factor-of-2 correction). On the ONERA M6 half
    wing (chord x / lift y / span z) it reproduces the analytic
    0.5*(C_ROOT+C_TIP)*B_SEMI (design.md Sec 9).

    Args:
        nodes: (n_nodes, 3) coordinates.
        wall_faces: (n, 3) wall triangles.
        thickness_axis: coordinate the wing is thin in (lift axis; y=1 here).

    Returns:
        Half-wing planform area.
    """
    p = nodes[np.asarray(wall_faces, dtype=np.int64)]
    area_vec = 0.5 * np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    return 0.5 * float(np.sum(np.abs(area_vec[:, thickness_axis])))


def cl_kj_3d(
    gamma_stations: np.ndarray,
    station_z: np.ndarray,
    s_ref: float,
    b_semi: float,
    u_inf: float = 1.0,
) -> float:
    """3D Kutta-Joukowski lift coefficient from the spanwise circulation.

    Lift per unit span is rho_inf U_inf Gamma(z), so the half-wing lift is
    L = rho_inf U_inf integral_0^{b/2} Gamma(z) dz and, nondimensionally
    (rho_inf = 1, lengths in mesh units),

        CL_KJ = 2 * integral_0^{b/2} Gamma(z) dz / (U_inf * S_ref).

    This is the 3D analog of the 2.5D `cl = 2 Gamma` consistency check
    (design.md Sec 9): for an untapered wing (Gamma const, S = c*b/2) it
    collapses to 2 Gamma / c, the 2D section value.

    The integral is closed at both ends: the root (z=0, symmetry plane) has
    zero spanwise slope by even symmetry, and the tip carries Gamma = 0
    discretely (the tip TE corner is a single-valued free node excluded from
    the stations, see mesh/wake_cut.py). `station_z`/`gamma_stations` are the
    solver's per-station arrays (index-aligned); they are sorted here, so no
    pre-sorting is assumed.

    Args:
        gamma_stations: (n_st,) circulation per wake station.
        station_z: (n_st,) spanwise coordinate of each station.
        s_ref: half-wing reference area (use `planform_area`).
        b_semi: semi-span (tip z); the integral is closed to here with Gamma=0.
        u_inf: freestream speed.

    Returns:
        CL_KJ, comparable to `wall_force_coefficients(...)["cl"]`.
    """
    z = np.asarray(station_z, dtype=np.float64)
    g = np.asarray(gamma_stations, dtype=np.float64)
    order = np.argsort(z)
    z, g = z[order], g[order]
    tol = 1e-9 * b_semi
    if z[0] > tol:  # root not sampled: even symmetry -> flat extension to z=0
        z = np.concatenate(([0.0], z))
        g = np.concatenate(([g[0]], g))
    if z[-1] < b_semi - tol:  # close to the tip where Gamma is pinned to 0
        z = np.concatenate((z, [b_semi]))
        g = np.concatenate((g, [0.0]))
    return 2.0 * float(np.trapezoid(g, z)) / (u_inf * s_ref)


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
