"""Track V V5 Stages 1+2 -- fixed/linear sparse operators of the augmented
Newton system (binding: docs/roadmap/track_v.md GV5.1 + the 2026-07-22
pre-registered FD note; design: cases/analysis/v5_tight_coupling/
PRE_REGISTRATION.md).

The augmented residual couples the IBL state U into the inviscid reduced
residual through the transpiration chain of the loose loop
(viscous/coupling.py:710-835):

    R_red(phi, Gamma, U) = T^T (R_inv(phi, Gamma) - b(m(U))),
    m(U) = P . div_Gamma(rho_e u_e delta*(U))          (surface nodes),
    b    = assemble_transpiration_rhs(nodes, wall_faces, S m)  (cut mesh),

so the BL -> phi coupling block at a frozen outer state (phi, edge data)
is the assembled sparse product

    J_phi,BL = -T^T . W . S . P . L(rho_e, u_e) . D_ds*(U)

with every factor FIXED at that state (the PRE_REGISTRATION "all pieces
exist; the block is assembled, not re-derived" clause):

* W -- wall-load matrix of assemble_transpiration_rhs
  (viscous/transpiration.py:66-114), linear in m. Per wall triangle of
  area a the 3-point edge-midpoint quadrature (exact for a P1 m) gives
  rhs[vi] += -(a/6) m_i - (a/12) (m_j + m_k): each edge-midpoint point
  carries (a/3) * 0.5 * 0.5 * (m_i + m_j) onto its two edge vertices, so
  the diagonal collects twice (two incident edges).
* S -- 0/1 volume scatter m_vol = S m_surf via SurfaceMesh.volume_node_of
  (coupling.py:832-833).
* P -- diagonal 0/1 outflow/tail pin mask (coupling.py:822-831; identity
  on airfoil cases, where outflow_pin_surf is None).
* L -- the surface-divergence operator delta* -> m_surf at fixed
  (rho_e, u_e) (transpiration.py:117-174): per triangle e,
  div_e = sum_b rho_b delta*_b (u_e,b . gradN[e,b]), then mass-lumped
  m_i = (1/node_area_i) sum_{e in i} (a_e/3) div_e.
* D_ds* -- per-node 1x6 closure derivative row douts[:, OUT_DS1, :] of
  closures.closure_all (closures.py:870); the closure's internal DELTA_MIN
  floor masks d/d(delta) consistently (closures.py:448-454), which is what
  makes the FD gate well-defined at probe states.
* T -- the cut -> reduced wake-constraint transform inside
  NewtonWorkspace (constraints/wake.py:50-53; R_red = T^T R_vol at
  solve/newton.py:378), free-DOF rows selected by ws.free
  (solve/newton.py:379).

G (also built here): the phi -> u_e recovery map of
edge_velocity_per_zone (transpiration.py:177-229), LINEAR in phi at fixed
geometry/zones -- the GV5.1 FD-note decision: the zone partition is
state-independent (fixed geometry, x/c < le_band_x), so the recovery is a
fixed sparse linear map with no kink. Needed by BOTH phi-side blocks in
later stages. Quadratic zone: the per-node least-squares fit of
post/surface.py::wall_tangential_gradient_quadratic (590-678); LE band:
the P1 per-triangle gradient -> crease-gated smoothing -> area-weighted
nodal average chain (post/surface.py:80-118, 147-189;
transpiration.py:206-228).

Stage 2 (the BL <- phi block):

* D_ue -- the per-node (7, 3) Jacobian of the IBL3 edge-data packet
  (q, rho, mu, M, u_hat_1..3) w.r.t. u_e (edge_data_jacobian below): the
  coupling.py:723-726 chain q = |u_e|, rho = density_field(q^2),
  mu = 1/re_chord (constant), M = sqrt(mach_squared_field(q^2)),
  u_hat = u_e/q.
* J_BL,phi = J_e D_ue G (assemble_j_bl_phi): J_e is the IBL3 edge-data
  Jacobian of IBL3Solver.residual_edge_jacobian (ibl3.py:1398-1411; state
  fixed, veps/veps_s frozen by decision 5 -- the FD gate measures the
  omission).

Assembly happens once per case/outer state: plain numpy/scipy.sparse, no
numba (the viscous/coupling.py orchestration precedent).
"""

from typing import Optional

import numpy as np
import scipy.sparse as sp

from pyfp3d.post.surface import (
    _wall_node_rings,
    _wall_vertex_normals,
    wall_outward_normals,
    wall_triangle_adjacency,
)
from pyfp3d.physics.isentropic import (
    density_derivative_wrt_q_sq_field,
    mach_squared_derivative_wrt_q_sq_field,
    mach_squared_field,
)
from pyfp3d.viscous import closures as C

__all__ = [
    "wall_load_matrix",
    "surface_scatter_matrix",
    "pin_mask_matrix",
    "surface_divergence_delta_operator",
    "closure_ds1_jacobian",
    "edge_velocity_operator",
    "assemble_j_phi_bl",
    "edge_data_jacobian",
    "assemble_j_bl_phi",
]


def wall_load_matrix(nodes: np.ndarray, wall_faces: np.ndarray) -> sp.csr_matrix:
    """W: the matrix form of assemble_transpiration_rhs
    (viscous/transpiration.py:66-114), (n_cut, n_cut) CSR, nonzero columns
    (and rows) at wall nodes only.

    Per wall triangle (area a, vertices i, j, k) the edge-midpoint
    quadrature scatters (a/3) * 0.5 * 0.5 * (m_i + m_j) from each edge
    midpoint to its two vertices, so rhs_i = -(a/6) m_i - (a/12)(m_j+m_k)
    (the diagonal collects the two incident edges; the solver RHS is the
    NEGATIVE load -- the transpiration.py sign convention). Exact for a
    P1 m field: W @ m_vol == assemble_transpiration_rhs(nodes, wall_faces,
    m_vol) to machine precision (unit-tested in
    tests/test_v5_tight_jacobian.py).
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    wall_faces = np.asarray(wall_faces)
    tri = nodes[wall_faces]  # (F, 3, 3)
    area_vec = 0.5 * np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = np.linalg.norm(area_vec, axis=1)
    if np.any(area < 1e-30):
        raise ValueError("degenerate wall facet (zero area)")

    n_t = len(wall_faces)
    # per_tri[e, a, b] = coefficient of m at vertex b in rhs at vertex a.
    per_tri = np.empty((n_t, 3, 3), dtype=np.float64)
    per_tri[:] = (-area / 12.0)[:, None, None]
    idx = np.arange(3)
    per_tri[:, idx, idx] = (-area / 6.0)[:, None]
    rows = np.repeat(wall_faces, 3, axis=1)  # (F, 9): vertex a, thrice
    cols = np.tile(wall_faces, (1, 3))  # (F, 9): vertex b, cycling
    return sp.coo_matrix(
        (per_tri.reshape(-1), (rows.reshape(-1), cols.reshape(-1))),
        shape=(len(nodes), len(nodes)),
    ).tocsr()


def surface_scatter_matrix(smesh, n_volume: Optional[int] = None) -> sp.csr_matrix:
    """S: the 0/1 scatter m_vol = S @ m_surf with
    m_vol[sm.volume_node_of] = m_surf (coupling.py:832-833),
    (n_volume, n_s) CSR, exactly one 1 per column."""
    if n_volume is None:
        n_volume = len(smesh.node_map)
    return sp.coo_matrix(
        (
            np.ones(smesh.n_node, dtype=np.float64),
            (smesh.volume_node_of, np.arange(smesh.n_node, dtype=np.int64)),
        ),
        shape=(n_volume, smesh.n_node),
    ).tocsr()


def pin_mask_matrix(
    n_s: int, outflow_pin_surf: Optional[np.ndarray] = None
) -> sp.csr_matrix:
    """P: the diagonal 0/1 transpiration pin mask (coupling.py:822-831):
    boundary-data delta* at pinned (outflow/tail/tip-band) nodes generates
    no transpiration. None (airfoil cases) = identity."""
    diag = np.ones(n_s, dtype=np.float64)
    if outflow_pin_surf is not None:
        diag[np.asarray(outflow_pin_surf, dtype=bool)] = 0.0
    return sp.diags(diag, 0, shape=(n_s, n_s), format="csr")


def surface_divergence_delta_operator(
    smesh, rho_e: np.ndarray, ue_surf: np.ndarray
) -> sp.csr_matrix:
    """L: the matrix form of transpiration_from_delta_star at fixed
    (rho_e, u_e) (viscous/transpiration.py:117-174), (n_s, n_s) CSR:
    m_surf = L @ delta*.

    Per triangle e the strong-form P1 divergence is
    div_e = sum_b rho_b delta*_b (u_e,b . gradN[e,b]); the lumped nodal
    projection gives m_i = (1/node_area_i) sum_{e in i} (a_e/3) div_e, so
    per (e, local vertex a, local vertex b):

        L[tri[e,a], tri[e,b]] += (a_e/3) rho[n_b] (u_e[n_b] . gradN[e,b])
                                 / node_area[n_a].

    L @ ds == transpiration_from_delta_star(sm, rho_e, ue, ds) to machine
    precision (unit-tested).
    """
    rho_e = np.asarray(rho_e, dtype=np.float64)
    ue_surf = np.asarray(ue_surf, dtype=np.float64)
    tri = smesh.triangles  # (F, 3)
    # per (e, b): rho_b (u_e,b . gradN[e,b]) -- the div_e coefficient of
    # delta*_b before lumping
    dot = np.einsum("ebk,ebk->eb", ue_surf[tri], smesh.gradN)  # (F, 3)
    coeff = (smesh.area_tri[:, None] / 3.0) * rho_e[tri] * dot  # (F, 3)
    rows = np.repeat(tri, 3, axis=1)  # (F, 9): vertex a, thrice
    cols = np.tile(tri, (1, 3))  # (F, 9): vertex b, cycling
    vals = np.tile(coeff, (1, 3)) / smesh.node_area[rows]
    return sp.coo_matrix(
        (vals.reshape(-1), (rows.reshape(-1), cols.reshape(-1))),
        shape=(smesh.n_node, smesh.n_node),
    ).tocsr()


def closure_ds1_jacobian(douts: np.ndarray) -> sp.csr_matrix:
    """D_ds*: the per-node 1x6 delta* derivative row of the closure packet
    (douts[:, OUT_DS1, :], closures.py:870), as an (n_s, 6 n_s) CSR matrix
    against the IBL3 flat state layout (6 unknowns per node, ibl3.py
    layout)."""
    douts = np.asarray(douts, dtype=np.float64)
    n_s = douts.shape[0]
    rows = np.repeat(np.arange(n_s, dtype=np.int64), 6)
    cols = np.arange(6 * n_s, dtype=np.int64)
    return sp.coo_matrix(
        (douts[:, C.OUT_DS1, :].reshape(-1), (rows, cols)),
        shape=(n_s, 6 * n_s),
    ).tocsr()


def _quadratic_fit_weights(nodes, wall_faces, volume_node_of, skip_vol):
    """Per-wall-node quadratic-recovery weights, replicating
    post/surface.py::wall_tangential_gradient_quadratic (644-677) fit by
    fit: same vertex normals, same rings (the private helpers are imported
    so the patch sets are the reference function's own), same 2-ring
    expansion and rank fallback.

    The phi -> gradient map at node i is c = (c1, c2) with
    [c1; c2] = M^+ (phi[patch] - phi[i]) and grad = c1 e1 + c2 e2; the
    weights are the pseudoinverse rows, extracted by running the fit's own
    lstsq call on the identity RHS (same SVD path, same rcond=1e-6, so the
    rank/fallback decision is the reference's). Returns a list of
    (surf_row, col_nodes, weights) triples with weights (k, 3): the
    3-vector coefficient of phi[col] in u_e[node], plus the closing
    self-entry -sum(weights) at col = node itself.
    """
    vertex_normal = _wall_vertex_normals(nodes, wall_faces)
    rings = _wall_node_rings(wall_faces, len(nodes))
    out = []
    for si in range(len(volume_node_of)):
        i = int(volume_node_of[si])
        if skip_vol[i]:
            continue
        n = vertex_normal[i]
        seed = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array(
            [0.0, 1.0, 0.0]
        )
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
        M = np.column_stack([u, v, u * u, u * v, v * v])
        X, _, rank, _ = np.linalg.lstsq(M, np.eye(len(patch)), rcond=1e-6)
        if rank < 5:
            # the reference's degenerate-patch fallback: 2-parameter
            # linear fit in the same tangent plane
            X = np.linalg.lstsq(M[:, :2], np.eye(len(patch)), rcond=1e-6)[0]
        else:
            X = X[:2]
        # grad = c1 e1 + c2 e2 with [c1; c2] = X (phi[patch] - phi[i])
        a = X[0][:, None] * e1[None, :] + X[1][:, None] * e2[None, :]  # (k,3)
        out.append((si, patch, a))
    return out


def _le_band_chain(nodes, wall_faces, elements, n_smooth_passes):
    """The LE-band linear+smoothed recovery chain as one sparse scalar
    map: returns (SC, A1) with A1 (3F, n_cut) the per-triangle P1 gradient
    operator (replicating post/surface.py::triangle_tangential_gradients,
    lines 97-117) and SC (n_cut, F) the composed smoothing + nodal
    average (post/surface.py:147-189 + transpiration.py:206-228). The
    component map is kron(SC, I3) @ A1 in the interleaved (3e+c) layout;
    the smoothing weights are geometry-fixed, so the composition is exact.
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
    u1 = np.sum(edge1 * e1, axis=1)
    v1 = np.sum(edge1 * e2, axis=1)
    edge2 = p2 - p0
    u2 = np.sum(edge2 * e1, axis=1)
    v2 = np.sum(edge2 * e2, axis=1)
    det = u1 * v2 - v1 * u2
    # grad = du e1 + dv e2 with du = (v2 df1 - v1 df2)/det,
    # dv = (-u2 df1 + u1 df2)/det, df1 = phi[n1]-phi[n0], df2 likewise
    c1 = (v2 / det)[:, None] * e1 + (-u2 / det)[:, None] * e2  # df1 coeff
    c2 = (-v1 / det)[:, None] * e1 + (u1 / det)[:, None] * e2  # df2 coeff
    n_t = len(wall_faces)
    n_cut = len(nodes)
    slots = np.stack([-(c1 + c2), c1, c2], axis=1)  # (F, slot, xyz)
    rows_a1 = np.broadcast_to(
        3 * np.arange(n_t)[:, None, None] + np.arange(3)[None, None, :],
        (n_t, 3, 3),
    )
    cols_a1 = np.broadcast_to(wall_faces[:, :, None], (n_t, 3, 3))
    A1 = sp.coo_matrix(
        (slots.reshape(-1), (rows_a1.reshape(-1), cols_a1.reshape(-1))),
        shape=(3 * n_t, n_cut),
    ).tocsr()

    area = 0.5 * twice_area
    if elements is not None:
        out_normal = wall_outward_normals(nodes, elements, wall_faces)
    else:
        out_normal = tri_normal
    adj = wall_triangle_adjacency(wall_faces)
    # one-pass scalar smoothing map Sp (F, F), replicating
    # smooth_wall_tangential_gradients (post/surface.py:170-189) at fixed
    # geometry: g_t <- (a_t g_t + sum_nb a_nb g_nb) / (a_t + sum_nb a_nb)
    # over aligned (crease-gated) neighbours
    wsum = area.copy()
    sp_rows = [np.arange(n_t, dtype=np.int64)]
    sp_cols = [np.arange(n_t, dtype=np.int64)]
    sp_vals = [area.copy()]
    for k in range(3):
        nb = adj[:, k]
        valid = nb >= 0
        dots = np.zeros(n_t, dtype=np.float64)
        dots[valid] = np.sum(
            out_normal[valid] * out_normal[nb[valid]], axis=1
        )
        m = valid & (dots > 0.2)  # cos_thresh of the reference
        wsum[m] += area[nb[m]]
        sp_rows.append(np.where(m)[0])
        sp_cols.append(nb[m])
        sp_vals.append(area[nb[m]])
    Sp = sp.coo_matrix(
        (
            np.concatenate(sp_vals) / wsum[np.concatenate(sp_rows)],
            (np.concatenate(sp_rows), np.concatenate(sp_cols)),
        ),
        shape=(n_t, n_t),
    ).tocsr()
    # n_passes identical applications (weights fixed) -> Sp ** n_passes
    S_tot = sp.eye(n_t, format="csr")
    for _ in range(max(int(n_smooth_passes), 0)):
        S_tot = Sp @ S_tot

    # area-weighted nodal average (transpiration.py:218-227):
    # ue[node] = sum_{e in node} g_e a_e / sum_{e in node} a_e
    nw = np.zeros(n_cut, dtype=np.float64)
    np.add.at(nw, wall_faces.reshape(-1), np.repeat(area, 3))
    rows_a3 = wall_faces.reshape(-1)
    cols_a3 = np.repeat(np.arange(n_t, dtype=np.int64), 3)
    A3 = sp.coo_matrix(
        (np.repeat(area, 3) / nw[rows_a3], (rows_a3, cols_a3)),
        shape=(n_cut, n_t),
    ).tocsr()
    return A3 @ S_tot, A1


def edge_velocity_operator(
    nodes: np.ndarray,
    wall_faces: np.ndarray,
    volume_node_of: Optional[np.ndarray] = None,
    elements: Optional[np.ndarray] = None,
    le_band_mask: Optional[np.ndarray] = None,
    n_smooth_passes: int = 2,
) -> sp.csr_matrix:
    """G: the matrix form of edge_velocity_per_zone
    (viscous/transpiration.py:177-229) restricted to the wall rows actually
    consumed downstream, (3 n_s, n_cut) CSR with row 3*si+c = component c
    of u_e at surface node si (volume id volume_node_of[si]):
    (G @ phi).reshape(n_s, 3) == edge_velocity_per_zone(...)[volume_node_of]
    to machine precision (unit-tested on both zone paths).

    LINEAR in phi at fixed geometry/zones (the GV5.1 FD-note decision:
    zones are state-independent, so no kink rows exist). Quadratic-zone
    rows carry the node's lstsq fit weights (_quadratic_fit_weights); LE
    rows carry the composed P1-gradient -> crease-gated smoothing ->
    area-weighted nodal average chain (_le_band_chain).
    ``le_band_mask`` is the caller-supplied (n_cut,) volume bool of
    transpiration.py (None/empty = all-quadratic, the A4 default).
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    wall_faces = np.asarray(wall_faces)
    n_cut = len(nodes)
    if volume_node_of is None:
        volume_node_of = np.unique(wall_faces)
    volume_node_of = np.asarray(volume_node_of, dtype=np.int64)
    n_s = len(volume_node_of)
    le = np.zeros(n_cut, dtype=bool)
    if le_band_mask is not None:
        le = np.asarray(le_band_mask, dtype=bool)

    rows_l, cols_l, vals_l = [], [], []
    for si, patch, a in _quadratic_fit_weights(
        nodes, wall_faces, volume_node_of, le
    ):
        cols = np.concatenate([patch, [volume_node_of[si]]])
        vals = np.vstack([a, -a.sum(axis=0, keepdims=True)])  # (k+1, 3)
        for c in range(3):
            rows_l.append(np.full(len(cols), 3 * si + c, dtype=np.int64))
            cols_l.append(cols)
            vals_l.append(vals[:, c])

    if le.any():
        SC, A1 = _le_band_chain(nodes, wall_faces, elements, n_smooth_passes)
        G_le = sp.kron(SC, sp.eye(3, format="csr"), format="csr") @ A1
        le_surf = np.where(le[volume_node_of])[0]
        src = (3 * volume_node_of[le_surf, None] + np.arange(3)).reshape(-1)
        dst = (3 * le_surf[:, None] + np.arange(3)).reshape(-1)
        sub = G_le[src, :].tocoo()
        rows_l.append(dst[sub.row])
        cols_l.append(sub.col)
        vals_l.append(sub.data)

    return sp.coo_matrix(
        (
            np.concatenate(vals_l),
            (np.concatenate(rows_l), np.concatenate(cols_l)),
        ),
        shape=(3 * n_s, n_cut),
    ).tocsr()


def assemble_j_phi_bl(ws, W, S, P, L, D) -> sp.csr_matrix:
    """J_phi,BL = -(T^T W S P L D)[ws.free, :], the BL -> reduced-inviscid
    coupling block of the augmented Jacobian (module docstring;
    PRE_REGISTRATION "Jacobian blocks"). ``ws`` is a solve/newton.py
    NewtonWorkspace (con.T = the wake-constraint reduce transform,
    ws.free = the free reduced rows of R_red, solve/newton.py:378-379).

    Shape (ws.n_free, 6 n_s), CSR. All factors are fixed at the frozen
    outer state; the block is assembled, not re-derived.
    """
    J_cut = W @ S @ P @ L @ D  # (n_cut, 6 n_s): the wall RHS derivative
    J_red = ws.con.T.T @ J_cut  # (n_red, 6 n_s)
    return (-J_red[ws.free, :]).tocsr()


def edge_data_jacobian(
    ue_surf: np.ndarray, m_inf: float, gamma_air: float = 1.4
) -> sp.csr_matrix:
    """D_ue: the per-node (7, 3) Jacobian of the IBL3 edge-data packet
    (q, rho, mu, M, u_hat_1..3) w.r.t. the edge velocity u_e, as a
    block-diagonal (7 n_s, 3 n_s) CSR matrix (one dense 7x3 block per
    surface node).

    The packet map is the coupling.py:723-726 chain evaluated per node:
    q = |u_e|, rho = density_field(q^2), mu = 1/re_chord (constant -- zero
    row, coupling.py:681), M = sqrt(mach_squared_field(q^2)),
    u_hat = u_e/q. The block rows are

        q:     u_hat^T
        rho:   2 rho'(q^2) u_e^T         (density_derivative_wrt_q_sq)
        mu:    0
        M:     m2'(q^2) u_e^T / sqrt(m2) (mach_squared_derivative_wrt_q_sq;
               dM = dm2/(2 sqrt(m2)) dq^2, dq^2 = 2 u_e^T du_e)
        u_hat: (I - u_hat u_hat^T) / q

    Nodes with q <= 1e-12 (the _frames fallback, ibl3.py:105-112) get a
    zero block: the packet map is degenerate there and the FD gate masks
    them (tests/test_v5_tight_edge.py).
    """
    ue = np.asarray(ue_surf, dtype=np.float64)
    if ue.ndim != 2 or ue.shape[1] != 3:
        raise ValueError(f"ue_surf shape {ue.shape}, expected (n_s, 3)")
    n_s = ue.shape[0]
    q2 = np.sum(ue * ue, axis=1)
    q = np.sqrt(q2)
    live = q > 1.0e-12
    qs = np.where(live, q, 1.0)
    uhat = ue / qs[:, None]
    m2 = mach_squared_field(q2, m_inf, gamma_air)
    blk = np.zeros((n_s, 7, 3), dtype=np.float64)
    blk[:, 0, :] = uhat
    blk[:, 1, :] = (
        2.0
        * density_derivative_wrt_q_sq_field(q2, m_inf, gamma_air)[:, None]
        * ue
    )
    blk[:, 3, :] = (
        mach_squared_derivative_wrt_q_sq_field(q2, m_inf, gamma_air)[:, None]
        * ue
        / np.sqrt(np.maximum(m2, 1.0e-300))[:, None]
    )
    blk[:, 4:7, :] = (
        np.eye(3)[None, :, :]
        - uhat[:, :, None] * uhat[:, None, :]
    ) / qs[:, None, None]
    blk[~live] = 0.0
    ii = np.arange(n_s, dtype=np.int64)
    rows = np.broadcast_to(
        7 * ii[:, None, None] + np.arange(7, dtype=np.int64)[None, :, None],
        (n_s, 7, 3),
    ).reshape(-1)
    cols = np.broadcast_to(
        3 * ii[:, None, None] + np.arange(3, dtype=np.int64)[None, None, :],
        (n_s, 7, 3),
    ).reshape(-1)
    vals = blk.reshape(-1)
    nz = vals != 0.0
    return sp.coo_matrix(
        (vals[nz], (rows[nz], cols[nz])), shape=(7 * n_s, 3 * n_s)
    ).tocsr()


def assemble_j_bl_phi(J_e, D_ue, G) -> sp.csr_matrix:
    """J_BL,phi = J_e D_ue G: the phi -> IBL-residual block of the
    augmented Jacobian (GV5.1 Stage 2; the chain phi -> u_e -> edge packet
    -> R_BL with the IBL state fixed).

    J_e is the (6 n_s, 7 n_s) edge-data Jacobian of
    IBL3Solver.residual_edge_jacobian (ibl3.py:1398-1411; state fixed,
    veps/veps_s frozen by decision 5 -- the FD gate measures the
    omission), D_ue the (7 n_s, 3 n_s) packet Jacobian of
    edge_data_jacobian, G the (3 n_s, n_cut) phi -> u_e operator of
    edge_velocity_operator. Returns the (6 n_s, n_cut) CSR product."""
    return (J_e @ D_ue @ G).tocsr()
