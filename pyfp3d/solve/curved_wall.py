"""
Curved wall-adjacent elements (P11; design.md §5.1.3, the DP1 "> 5%" branch's
curved-element route for gate G1.6).

The G1.3/G1.4 oracles proved the sphere's wall-Cp error lives in the
*integration domain* (the flat-facet polyhedron vs the true curved surface),
not in the boundary-condition data — so the fix here changes the domain and
nothing else. Every tet containing at least one wall-surface edge gets a
quadratic (tet10-style) geometric map: the midpoints of wall edges are
projected onto the true surface through the replaceable
closest_point_normal(points) callback (solve/wall_correction.py), midpoints of
all other edges stay straight. Any tet sharing a projected edge contains a
wall edge by construction, so shared faces see identical quadratic geometry
from both sides — the curved layer is conforming.

The scalar field keeps its 4 vertex DOFs per tet (mapped P1): on a curved
element the basis is the reference-linear functions composed with the
quadratic map, and the element stiffness is integrated by quadrature,

    K_kl = sum_q w_q (J^{-T} g_k) . (J^{-T} g_l) |det J(xi_q)|

with g_k the constant reference gradients and J the (linear-in-xi) Jacobian
of the tet10 map. The assembled object is a sparse DELTA on the standard P1
matrix, dA = sum_curved (K(curved geometry) - K(straight geometry)), so
A_curved = A_P1 + dA keeps the CSR pattern, symmetry and (with Dirichlet
rows) SPD-ness, and the CG+AMG path is reused unchanged. Both K's go through
the SAME quadrature code, so when a projection moves nothing (planar
closest_point_normal on a planar wall) the delta is exactly zero, bitwise —
the G11.3 null test.

Pre-registered risk (design.md §5.1.3): the mapped-P1 space on a curved
element contains physical linears only up to the O(h^2) map perturbation
(superparametric consistency). Confined to this opt-in path; the V0
freestream gate's wall-free configs never build a delta.

OUTCOME (measured 2026-07-19, phase P11 close-out -- roadmap track_p.md §P11,
demo cases/demo/p11_curved_walls/): the risk FIRED and the route is a
recorded NEGATIVE for gate G1.6. The mapped-P1 gradient-consistency error is
O(h) (linear fields are not reproduced on curved elements: max deviation
0.138 on the coarse sphere), the same order as the O(h) facet-normal error
the curving removes, so the net medium-sphere Cp gain is 11.56% -> 11.33% --
identical to the G1.4 boundary-data oracle ceiling. The deeper measurement:
the wall geometric crime was never the dominant G1.6 error (structured-shell
control converges at ~2nd order with flat facets; the h_min-sweep order
collapse is the fixed-bulk-mesh pollution floor). The module stays as shipped
curved-geometry infrastructure + the evidence machinery for the locks in
tests/test_p11_curved_walls.py.

Vectorized numpy throughout (wall_correction.py precedent: the curved layer
is a thin O(surface) subset; njit only if it ever becomes hot).
"""

from typing import Callable, Dict, Tuple

import numpy as np
import scipy.sparse as sp

# Local edge numbering of the tet10 geometric map: geometry node 4 + e is the
# midpoint of vertex pair _TET_EDGES[e].
_TET_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

# Reference gradients of the barycentric coordinates lambda_a wrt xi
# (lambda_0 = 1 - xi1 - xi2 - xi3, lambda_a = xi_a).
_L = np.array(
    [[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    dtype=np.float64,
)

# Quadrature rules on the reference tet, given in barycentric coordinates;
# weights sum to the reference volume 1/6.
_A4 = 0.5854101966249685
_B4 = 0.1381966011250105
_RULES = {
    # 4-point, degree 2 (the default; design.md §5.1.3)
    "deg2": (
        np.array(
            [
                [_A4, _B4, _B4, _B4],
                [_B4, _A4, _B4, _B4],
                [_B4, _B4, _A4, _B4],
                [_B4, _B4, _B4, _A4],
            ],
            dtype=np.float64,
        ),
        np.full(4, 1.0 / 24.0, dtype=np.float64),
    ),
    # 5-point, degree 3 (the G11.4 quadrature A/B check; negative center
    # weight is fine for a diagnostic rule)
    "deg3": (
        np.array(
            [
                [0.25, 0.25, 0.25, 0.25],
                [0.5, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0],
                [1.0 / 6.0, 0.5, 1.0 / 6.0, 1.0 / 6.0],
                [1.0 / 6.0, 1.0 / 6.0, 0.5, 1.0 / 6.0],
                [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 0.5],
            ],
            dtype=np.float64,
        ),
        np.array([-2.0 / 15.0, 3.0 / 40.0, 3.0 / 40.0, 3.0 / 40.0, 3.0 / 40.0]),
    ),
}


def plane_closest_point_normal(
    origin: Tuple[float, float, float], normal: Tuple[float, float, float]
) -> Callable[[np.ndarray], Tuple[np.ndarray, np.ndarray]]:
    """closest_point_normal callback for a PLANE (the G11.3 null test: on a
    planar wall whose nodes lie in the plane, edge midpoints project onto
    themselves bitwise, so the stiffness delta must be exactly zero)."""
    o = np.asarray(origin, dtype=np.float64)
    n = np.asarray(normal, dtype=np.float64)
    n = n / np.linalg.norm(n)

    def _callback(points: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        p = np.asarray(points, dtype=np.float64)
        dist = (p - o) @ n
        proj = p - dist[:, None] * n
        normals = np.broadcast_to(n, p.shape).copy()
        return proj, normals

    return _callback


def _tet10_dshape(bary: np.ndarray) -> np.ndarray:
    """Reference derivatives dN_a/dxi_j of the 10 tet10 shape functions at
    the given barycentric points. Returns (10, Q, 3)."""
    q = np.asarray(bary, dtype=np.float64)  # (Q, 4)
    n_q = len(q)
    dN = np.zeros((10, n_q, 3), dtype=np.float64)
    for a in range(4):  # vertex functions N_a = lambda_a (2 lambda_a - 1)
        dN[a] = (4.0 * q[:, a] - 1.0)[:, None] * _L[a]
    for e, (a, b) in enumerate(_TET_EDGES):  # edge functions N = 4 la lb
        dN[4 + e] = 4.0 * (q[:, a][:, None] * _L[b] + q[:, b][:, None] * _L[a])
    return dN


def curved_wall_geometry(
    nodes: np.ndarray,
    elements: np.ndarray,
    wall_faces: np.ndarray,
    closest_point_normal: Callable[[np.ndarray], Tuple[np.ndarray, np.ndarray]],
) -> Dict[str, np.ndarray]:
    """Identify the wall-adjacent curved layer and build its tet10 geometry.

    Returns dict with:
        curved_tets   (C,)      indices of tets containing >= 1 wall edge
        geom10        (C,10,3)  curved geometry (wall-edge midpoints projected)
        geom10_flat   (C,10,3)  straight-midpoint geometry (same layout)
        edge_projected(C,6)     which local edges got projected
        max_offset, mean_offset   |projection - straight midpoint| stats over
                                  the projected midpoints (O(h^2) self-check)
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    elements = np.asarray(elements)
    wall_faces = np.asarray(wall_faces)
    n_nodes = len(nodes)

    # Unique wall-surface edges, encoded as min*n + max.
    wf = wall_faces.astype(np.int64)
    pairs = np.concatenate([wf[:, [0, 1]], wf[:, [1, 2]], wf[:, [2, 0]]], axis=0)
    codes = pairs.min(axis=1) * n_nodes + pairs.max(axis=1)
    wall_edge_codes = np.unique(codes)

    # Each tet's 6 edges, same encoding.
    el = elements.astype(np.int64)
    tet_edge_codes = np.empty((len(el), 6), dtype=np.int64)
    for e, (a, b) in enumerate(_TET_EDGES):
        lo = np.minimum(el[:, a], el[:, b])
        hi = np.maximum(el[:, a], el[:, b])
        tet_edge_codes[:, e] = lo * n_nodes + hi
    on_wall = np.isin(tet_edge_codes, wall_edge_codes)  # (T, 6)
    curved = np.flatnonzero(on_wall.any(axis=1))
    if len(curved) == 0:
        raise ValueError("no tet contains a wall edge: wall_faces empty or inconsistent")

    verts = nodes[el[curved]]  # (C, 4, 3)
    mids = np.empty((len(curved), 6, 3), dtype=np.float64)
    for e, (a, b) in enumerate(_TET_EDGES):
        mids[:, e] = 0.5 * (verts[:, a] + verts[:, b])

    geom10_flat = np.concatenate([verts, mids], axis=1)  # (C, 10, 3)
    geom10 = geom10_flat.copy()

    edge_projected = on_wall[curved]  # (C, 6)
    flat_mids = mids[edge_projected]  # (P, 3), duplicates across tets are
    # bitwise-identical (0.5*(a+b) is commutative), so per-tet projection is
    # conforming without a shared edge table.
    proj, _ = closest_point_normal(flat_mids)
    geom10[:, 4:][edge_projected] = proj

    offsets = np.linalg.norm(proj - flat_mids, axis=1)
    return {
        "curved_tets": curved,
        "geom10": geom10,
        "geom10_flat": geom10_flat,
        "edge_projected": edge_projected,
        "max_offset": float(offsets.max()),
        "mean_offset": float(offsets.mean()),
    }


def _jacobians(geom10: np.ndarray, rule: str) -> Tuple[np.ndarray, np.ndarray]:
    """J (C,Q,3,3) and |det J| (C,Q) of the tet10 map at the rule's points.

    Tet vertex ordering carries no orientation guarantee in this codebase
    (the P1 path uses V = |det|/6, metrics.py), so the element's own affine
    (vertex) map fixes the reference sign; the quadratic map must keep that
    sign at every quadrature point, else the projection folded the element.
    """
    bary, _ = _RULES[rule]
    dN = _tet10_dshape(bary)  # (10, Q, 3)
    J = np.einsum("cai,aqj->cqij", geom10, dN)
    detJ = np.linalg.det(J)
    det_affine = np.linalg.det(np.einsum("cai,aj->cij", geom10[:, :4], _L))
    signed = detJ * np.sign(det_affine)[:, None]
    if signed.min() <= 0.0:
        n_bad = int(np.sum(signed.min(axis=1) <= 0.0))
        raise ValueError(
            f"Jacobian sign flip in {n_bad} curved tet(s): the midpoint "
            "projection inverted the quadratic map (sliver next to a strongly "
            "curved wall?)"
        )
    return J, signed


def element_stiffness_tet10(geom10: np.ndarray, rule: str = "deg2") -> np.ndarray:
    """Mapped-P1 element stiffness (C,4,4) on tet10 geometry by quadrature.
    Exact (= the affine P1 stiffness) when the geometry is straight, since J
    is then constant and every rule integrates constants exactly."""
    _, w = _RULES[rule]
    J, detJ = _jacobians(geom10, rule)
    Jinv = np.linalg.inv(J)  # (C,Q,3,3)
    # Physical gradients of the 4 field basis functions: grad lambda_k = J^{-T} L_k
    G = np.einsum("cqji,kj->cqki", Jinv, _L)  # (C,Q,4,3)
    K = np.einsum("q,cq,cqki,cqli->ckl", w, detJ, G, G, optimize=True)
    return 0.5 * (K + K.transpose(0, 2, 1))


def curved_volumes(geom10: np.ndarray, rule: str = "deg2") -> np.ndarray:
    """Per-element volume sum_q w_q |det J| of the tet10 map (the G11.4
    volume-fidelity diagnostic)."""
    _, w = _RULES[rule]
    _, detJ = _jacobians(geom10, rule)
    return detJ @ w


def assemble_curved_stiffness_delta(
    n_nodes: int,
    elements: np.ndarray,
    geometry: Dict[str, np.ndarray],
    rule: str = "deg2",
) -> sp.csr_matrix:
    """Assemble dA = sum_curved (K_curved - K_flat) as a sparse (n,n) CSR.

    A_P1 + dA is the curved-wall stiffness matrix: symmetric, same sparsity
    pattern (the 4x4 blocks sit on existing node pairs), SPD after Dirichlet
    elimination. Plugs into solve_laplace(stiffness_delta=...).
    """
    K_curved = element_stiffness_tet10(geometry["geom10"], rule)
    K_flat = element_stiffness_tet10(geometry["geom10_flat"], rule)
    dK = K_curved - K_flat  # (C,4,4)

    tets = np.asarray(elements)[geometry["curved_tets"]].astype(np.int64)  # (C,4)
    rows = np.repeat(tets, 4, axis=1).ravel()  # i index varies slow
    cols = np.tile(tets, (1, 4)).ravel()  # j index varies fast
    dA = sp.coo_matrix(
        (dK.ravel(), (rows, cols)), shape=(n_nodes, n_nodes)
    ).tocsr()
    return dA
