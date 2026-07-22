"""Track V V2 -- transpiration channel (binding: docs/roadmap/track_v.md
"V2 -- Transpiration channel through all three drivers", 2026-07-22 re-spec).

The viscous-inviscid coupling BC: the boundary layer displaces the effective
body, modelled as a mass-flux blowing distribution on the wall (no mesh
motion, RHS-only). Two pieces live here:

1. ``assemble_transpiration_rhs`` -- the wall-RHS assembly: the Galerkin
   load of a nodal mass-flux field m_dot on the wall triangulation,
   scattered into a full (n_nodes,) volume RHS vector (the
   ``body_source_rhs`` argument of ``solve/picard.py::solve_laplace`` and
   the newly threaded compressible drivers). Structural template:
   ``solve/wall_correction.py::assemble_wall_flux_correction_rhs`` (the
   roadmap-pinned template; same 3-point edge-midpoint quadrature, exact
   for a P1 m_dot field, same np.add.at nodal scatter). Plain vectorized
   numpy like the template -- it gets the njit treatment when V3 makes it
   hot.

   SIGN CONVENTION (pinned by GV2.1(a) against the analytic Fourier-mode
   cylinder solution): m_dot > 0 is BLOWING, mass flux out of the body
   into the fluid, i.e. dphi/dn_body = m_dot / rho_e with n_body the
   body-outward normal. The weak FP form integrates by parts over the
   FLOW domain, whose outward normal at the wall points INTO the body
   (n_dom = -n_body), so the solver RHS is the NEGATIVE load:

       b_i = integral_Gamma N_i * rho_e * dphi/dn_dom dS
           = -integral_Gamma N_i * m_dot dS.

   With m_dot == 0 the returned vector is exactly zero (bit-identical to
   the channel being absent -- the GV2.1(b) discipline).

2. The delta* -> m_dot operator, m_dot = div_Gamma(rho_e * u_e * delta*):
   the surface divergence of the defect mass flux Q = rho_e * delta* *
   u_e (the 3-D generalization of the 2-D transpiration rule
   v_n = d(u_e delta*)/ds; positive divergence of the defect flux blows).
   Built on viscous/surface_mesh.py SurfaceMesh (gradN strong-form
   divergence + node_area lumped projection). GV2.1 exercises the channel
   with MANUFACTURED m_dot only; this operator is unit-tested here on
   analytic fields and gets its first live exercise in V3's coupled gates.

Also ``edge_velocity_per_zone``: the per-zone u_e extraction discipline of
A4 (LE/stagnation band linear+smoothed, elsewhere quadratic recovery;
cases/analysis/a4_ue_error_band/VERDICT.md). Zone masks are caller-supplied
(V3 wires the actual zones); the default is all-quadratic.
"""

from typing import Optional

import numpy as np

from pyfp3d.post.surface import (
    smooth_wall_tangential_gradients,
    triangle_tangential_gradients,
    wall_outward_normals,
    wall_tangential_gradient_quadratic,
    wall_triangle_adjacency,
)

# 3-point edge-midpoint quadrature on a triangle: exact for quadratics
# (solve/wall_correction.py precedent; P1 m_dot makes N_i * m_dot
# quadratic). Point q is the midpoint of edge (i, j); the P1 shape
# functions there are N_i = N_j = 1/2, N_k = 0.
_EDGE_MIDPOINT_PAIRS = ((0, 1), (1, 2), (2, 0))


def assemble_transpiration_rhs(
    nodes: np.ndarray, wall_faces: np.ndarray, m_dot: np.ndarray
) -> np.ndarray:
    """Assemble the transpiration wall RHS b = -<N_i, m_dot>_Gamma into a
    full (n_nodes,) vector (the body_source_rhs argument of the FP
    drivers; sign convention in the module docstring).

    Args:
        nodes: (n_nodes, 3) volume mesh nodal coordinates.
        wall_faces: (F, 3) wall triangle connectivity as VOLUME node ids
            (e.g. mesh.boundary_faces["wall"]; on a wake-cut mesh the
            cut-mesh node ids, TE-duplicated slaves included).
        m_dot: (n_nodes,) nodal mass flux [rho * velocity units]; only
            the wall-node entries are read. m_dot > 0 = blowing (out of
            the body into the fluid), m_dot < 0 = suction.

    Returns:
        rhs (n_nodes,) float64: the Galerkin load -integral_Gamma N_i
        m_dot dS, exact (mass-matrix consistent) for a P1 m_dot field.
        m_dot == 0 returns the exact zero vector.
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    wall_faces = np.asarray(wall_faces)
    m_dot = np.asarray(m_dot, dtype=np.float64)
    if m_dot.shape != (len(nodes),):
        raise ValueError(
            f"m_dot must be ({len(nodes)},), got {m_dot.shape}"
        )

    tri = nodes[wall_faces]  # (F, 3, 3)
    area_vec = 0.5 * np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = np.linalg.norm(area_vec, axis=1)
    if np.any(area < 1e-30):
        raise ValueError("degenerate wall facet (zero area)")

    mw = m_dot[wall_faces]  # (F, 3) nodal m_dot
    # P1 value at the midpoint of edge (i, j) is (m_i + m_j) / 2.
    m_qp = 0.5 * (mw[:, [0, 1, 2]] + mw[:, [1, 2, 0]])  # (F, 3)
    # weight per quadrature point = area/3; shape-function value 1/2 at
    # the two edge vertices, 0 at the opposite one
    contrib = (area[:, None] / 3.0) * 0.5 * m_qp
    load = np.zeros(len(nodes), dtype=np.float64)
    for q, (i, j) in enumerate(_EDGE_MIDPOINT_PAIRS):
        np.add.at(load, wall_faces[:, i], contrib[:, q])
        np.add.at(load, wall_faces[:, j], contrib[:, q])
    # Blowing is -dphi/dn_dom (module docstring): the solver RHS is the
    # negative load. The sign is pinned by GV2.1(a) against the analytic
    # exterior Laplace solution.
    return -load


def surface_divergence_tri(smesh, vec: np.ndarray) -> np.ndarray:
    """Per-triangle constant strong-form surface divergence of a nodal
    vector field: div_e = sum_b vec[nn[b]] . gradN[e, b] (the ibl3.py
    assembly idiom). Exact for a P1 (piecewise-linear) vec.

    Args:
        smesh: viscous/surface_mesh.py SurfaceMesh.
        vec: (smesh.n_node, 3) nodal vector field (global XYZ).

    Returns:
        div_tri (smesh.n_tri,) float64.
    """
    vec = np.asarray(vec, dtype=np.float64)
    if vec.shape != (smesh.n_node, 3):
        raise ValueError(
            f"vec must be ({smesh.n_node}, 3), got {vec.shape}"
        )
    return np.einsum("ebk,ebk->e", vec[smesh.triangles], smesh.gradN)


def surface_divergence_nodal(smesh, vec: np.ndarray) -> np.ndarray:
    """Nodal surface divergence via the lumped projection of
    surface_divergence_tri (integral div dS lumped area/3 per vertex,
    divided by smesh.node_area -- the mass-lumping idiom). Exact at a
    node whose incident triangles all carry the same constant divergence
    (e.g. a globally linear vec field)."""
    div_tri = surface_divergence_tri(smesh, vec)
    num = np.zeros(smesh.n_node, dtype=np.float64)
    np.add.at(
        num,
        smesh.triangles.reshape(-1),
        np.repeat(div_tri * smesh.area_tri / 3.0, 3),
    )
    return num / smesh.node_area


def defect_mass_flux(
    rho_e: np.ndarray, ue_tangent: np.ndarray, delta_star: np.ndarray
) -> np.ndarray:
    """Nodal defect mass flux Q = rho_e * delta* * u_e, (n_s, 3). All
    inputs are nodal surface-mesh arrays (u_e tangential by construction
    -- see edge_velocity_per_zone)."""
    rho_e = np.asarray(rho_e, dtype=np.float64)
    ue_tangent = np.asarray(ue_tangent, dtype=np.float64)
    delta_star = np.asarray(delta_star, dtype=np.float64)
    return (rho_e * delta_star)[:, None] * ue_tangent


def transpiration_from_delta_star(
    smesh, rho_e: np.ndarray, ue_tangent: np.ndarray, delta_star: np.ndarray
) -> np.ndarray:
    """Nodal transpiration mass flux m_dot = div_Gamma(rho_e u_e delta*)
    on the surface mesh, (n_s,) float64. Positive = blowing. First live
    exercise is V3's coupled gates (roadmap V2 deliverable note); here it
    is unit-tested on analytic fields only."""
    return surface_divergence_nodal(
        smesh, defect_mass_flux(rho_e, ue_tangent, delta_star)
    )


def edge_velocity_per_zone(
    nodes: np.ndarray,
    wall_faces: np.ndarray,
    phi: np.ndarray,
    elements: Optional[np.ndarray] = None,
    le_band_mask: Optional[np.ndarray] = None,
    n_smooth_passes: int = 2,
) -> np.ndarray:
    """Per-zone wall tangential velocity u_e = grad_Gamma(phi), (n_nodes, 3)
    with NaN off-wall (the wall_tangential_gradient_quadratic convention).

    A4 discipline (cases/analysis/a4_ue_error_band/VERDICT.md): quadratic
    recovery (post/surface.py::wall_tangential_gradient_quadratic)
    everywhere except the LE/stagnation band, which takes the
    linear+smoothed path (per-triangle P1 gradients + crease-gated
    area-weighted smoothing, P6 infra) followed by an area-weighted nodal
    average. ``le_band_mask`` is a boolean (n_nodes,) caller-supplied
    zone mask (V3 wires the actual zones); None = all-quadratic.

    ``elements`` orients the smoothing's crease gate with body-outward
    normals (wall_outward_normals); without it the stored winding is used
    (consistent-winding meshes only).
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    ue = wall_tangential_gradient_quadratic(nodes, wall_faces, phi)
    if le_band_mask is None or not np.any(le_band_mask):
        return ue

    grad_tri, area, tri_normal = triangle_tangential_gradients(
        nodes, wall_faces, phi
    )
    if elements is not None:
        out_normal = wall_outward_normals(nodes, elements, wall_faces)
    else:
        out_normal = tri_normal
    adj = wall_triangle_adjacency(wall_faces)
    grad_smooth = smooth_wall_tangential_gradients(
        grad_tri, out_normal, area, adj, n_passes=n_smooth_passes
    )
    # Area-weighted nodal average of the smoothed per-triangle gradients
    # (mass-lumping idiom; LE band is smooth, so no crease risk there).
    num = np.zeros((len(nodes), 3), dtype=np.float64)
    wsum = np.zeros(len(nodes), dtype=np.float64)
    contrib = grad_smooth * area[:, None]
    for i in range(3):
        np.add.at(num, wall_faces[:, i], contrib)
        np.add.at(wsum, wall_faces[:, i], area)
    ue_linear = np.full((len(nodes), 3), np.nan, dtype=np.float64)
    touched = wsum > 0
    ue_linear[touched] = num[touched] / wsum[touched, None]
    ue[le_band_mask] = ue_linear[le_band_mask]
    return ue
