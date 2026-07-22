"""
Tests for pyfp3d/viscous/surface_mesh.py (Track V V1: compact wall-surface
DOF numbering + P1 surface-FE geometry for the IBL3 integral boundary
layer, design_track_v.md §5/§6).

Covers: the standalone structured-rectangle generator, compact numbering
round-trip, P1 shape-gradient identities, adjacency/boundary extraction vs
the post pipeline, the per-node local basis (incl. crease fallback),
greedy coloring determinism/validity, the block-6 symbolic CSR Jacobian
pattern + elem_to_csr scatter map, and a cross-check against
tests/mesh_utils.py + pyfp3d/post/surface.py on a volume-mesh wall face.

Run with: pytest tests/test_v1_surface_mesh.py -xvs
"""

import numpy as np
import pytest
import scipy.sparse as sp

from pyfp3d.post.surface import (
    triangle_tangential_gradients,
    wall_triangle_adjacency,
)
from pyfp3d.viscous.surface_mesh import (
    QUAD_N,
    SurfaceMesh,
    structured_rectangle_surface,
)

from .mesh_utils import generate_structured_cube_mesh


def _plate(nx=6, nz=4):
    xyz, tris = structured_rectangle_surface(0.0, 1.2, -0.4, 0.4, nx, nz)
    return SurfaceMesh.from_wall_faces(xyz, tris)


# ---------------------------------------------------------------------------
# 1. structured_rectangle_surface
# ---------------------------------------------------------------------------


def test_rectangle_counts_normals_area():
    nx, nz = 5, 3
    x0, x1, z0, z1 = -0.2, 1.0, 0.3, 1.1
    xyz, tris = structured_rectangle_surface(x0, x1, z0, z1, nx, nz)
    assert xyz.shape == ((nx + 1) * (nz + 1), 3)
    assert tris.shape == (2 * nx * nz, 3)
    assert xyz.dtype == np.float64 and tris.dtype == np.int64
    # Row-major in x then z: node iz*(nx+1)+ix sits at (x[ix], 0, z[iz]).
    xs = np.linspace(x0, x1, nx + 1)
    zs = np.linspace(z0, z1, nz + 1)
    for iz in range(nz + 1):
        for ix in range(nx + 1):
            np.testing.assert_allclose(
                xyz[iz * (nx + 1) + ix], [xs[ix], 0.0, zs[iz]], atol=0.0
            )
    # Winding => all geometric unit normals exactly +y.
    x = xyz[tris]  # (F,3,3)
    n = np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0])
    n /= np.linalg.norm(n, axis=1)[:, None]
    np.testing.assert_allclose(n[:, 1], 1.0, atol=1e-15)
    np.testing.assert_allclose(n[:, [0, 2]], 0.0, atol=1e-15)
    # Total area exact.
    area = 0.5 * np.linalg.norm(
        np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]), axis=1
    )
    assert abs(area.sum() - (x1 - x0) * (z1 - z0)) < 1e-14


def test_rectangle_determinism():
    a = structured_rectangle_surface(0.0, 1.0, 0.0, 2.0, 4, 7)
    b = structured_rectangle_surface(0.0, 1.0, 0.0, 2.0, 4, 7)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])


# ---------------------------------------------------------------------------
# 2. Compact numbering round-trip
# ---------------------------------------------------------------------------


def test_compact_numbering_roundtrip():
    # Embed a plate patch into a larger, scrambled volume node set so the
    # compact map is exercised (volume ids are NOT 0..n_s-1).
    xyz, tris = structured_rectangle_surface(0.0, 1.0, 0.0, 1.0, 3, 2)
    rng = np.random.default_rng(42)
    n_vol = len(xyz) + 17  # extra off-surface volume nodes
    perm = rng.permutation(n_vol)[: len(xyz)]  # volume ids of surface nodes
    nodes = np.zeros((n_vol, 3))
    nodes[perm] = xyz
    wall_faces = perm[tris]

    sm = SurfaceMesh.from_wall_faces(nodes, wall_faces, name="wall")
    assert sm.n_node == len(xyz)
    assert sm.n_tri == len(tris)
    # node_map[volume_node_of] == arange(n_s)
    np.testing.assert_array_equal(
        sm.node_map[sm.volume_node_of], np.arange(sm.n_node)
    )
    # sorted-unique volume ids
    np.testing.assert_array_equal(sm.volume_node_of, np.sort(perm))
    # triangles index within range and invert the map
    assert sm.triangles.min() >= 0 and sm.triangles.max() < sm.n_node
    np.testing.assert_array_equal(sm.volume_node_of[sm.triangles], wall_faces)
    # xyz matches the source nodes
    np.testing.assert_array_equal(sm.xyz, nodes[sm.volume_node_of])
    np.testing.assert_array_equal(sm.xyz[sm.triangles], xyz[tris])
    # master-map hook present
    assert sm.volume_node_of.dtype == np.int64
    # repr smoke
    assert "n_node" in repr(sm) and "wall" in repr(sm)


# ---------------------------------------------------------------------------
# 3. gradN identities
# ---------------------------------------------------------------------------


def _check_gradN_identities(sm, tol_sum, tol_delta):
    x = sm.xyz[sm.triangles]  # (F,3,3)
    scale = max(1.0, float(np.abs(sm.gradN).max()))
    # sum_j gradN_j = 0 (partition of unity)
    assert np.abs(sm.gradN.sum(axis=1)).max() < tol_sum * scale
    # gradN_j . (x_k - x_l) = delta_jk - delta_jl
    for j in range(3):
        for k in range(3):
            for l in range(3):
                val = np.sum(sm.gradN[:, j, :] * (x[:, k, :] - x[:, l, :]), axis=1)
                expect = float(j == k) - float(j == l)
                assert np.abs(val - expect).max() < tol_delta


def test_gradN_on_structured_plate():
    _check_gradN_identities(_plate(8, 5), 1e-15, 1e-14)


def test_gradN_on_scrambled_triangles():
    # A few arbitrary non-uniform, non-coplanar triangles built by hand.
    nodes = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.3, 0.1, 0.2],
            [0.2, 0.9, 1.1],
            [1.1, 1.4, 0.7],
            [0.5, 0.4, 0.6],
            [-0.3, 0.8, 0.4],
        ]
    )
    faces = np.array(
        [[0, 1, 4], [1, 3, 4], [4, 3, 2], [0, 4, 5], [4, 2, 5]], dtype=np.int64
    )
    sm = SurfaceMesh.from_wall_faces(nodes, faces)
    _check_gradN_identities(sm, 1e-15, 1e-14)


# ---------------------------------------------------------------------------
# 4. Adjacency / boundary
# ---------------------------------------------------------------------------


def test_adjacency_matches_post_and_boundary_counts():
    nx, nz = 6, 4
    xyz, tris = structured_rectangle_surface(0.0, 1.0, 0.0, 1.0, nx, nz)
    sm = SurfaceMesh.from_wall_faces(xyz, tris)
    # Identical to the post pipeline on the same connectivity.
    np.testing.assert_array_equal(sm.adjacency, wall_triangle_adjacency(tris))
    # nx x nz plate boundary edge count = 2*nx + 2*nz.
    assert len(sm.boundary_edges) == 2 * nx + 2 * nz
    # Boundary edges are exactly the -1 adjacency entries.
    tri_edges = ((0, 1), (1, 2), (2, 0))
    expect = set()
    for t in range(sm.n_tri):
        for k, (a, b) in enumerate(tri_edges):
            if sm.adjacency[t, k] < 0:
                expect.add(tuple(sorted((sm.triangles[t, a], sm.triangles[t, b]))))
    got = {tuple(sorted(e)) for e in sm.boundary_edges}
    assert got == expect
    # boundary_node_mask consistent with boundary_edges.
    mask = np.zeros(sm.n_node, dtype=bool)
    mask[sm.boundary_edges.reshape(-1)] = True
    np.testing.assert_array_equal(sm.boundary_node_mask, mask)
    # Interior nodes of the plate: (nx-1)*(nz-1) of them, none boundary.
    assert sm.boundary_node_mask.sum() == sm.n_node - (nx - 1) * (nz - 1)


# ---------------------------------------------------------------------------
# 5. Local basis
# ---------------------------------------------------------------------------


def test_basis_on_flat_plate():
    sm = _plate(5, 5)
    # basis_y == +/- yhat (winding is +y on this generator, but allow sign).
    dots_y = np.abs(sm.basis_y[:, 1])
    assert np.all(dots_y > 1.0 - 1e-14)
    # basis_x == +/- Xhat (|dot| = 1).
    dots_x = np.abs(sm.basis_x[:, 0])
    np.testing.assert_allclose(dots_x, 1.0, atol=1e-14)
    # Orthonormal triple to 1e-14.
    assert np.abs(np.sum(sm.basis_x * sm.basis_y, axis=1)).max() < 1e-14
    assert np.abs(np.sum(sm.basis_x * sm.basis_z, axis=1)).max() < 1e-14
    assert np.abs(np.sum(sm.basis_y * sm.basis_z, axis=1)).max() < 1e-14
    np.testing.assert_allclose(
        np.linalg.norm(sm.basis_z, axis=1), 1.0, atol=1e-14
    )
    # Right-handed: z = y x x.
    np.testing.assert_allclose(
        sm.basis_z, np.cross(sm.basis_y, sm.basis_x), atol=1e-14
    )


def test_basis_crease_fallback_no_raise():
    # Shallow-V (tent) strip: ridge along x at z=0.5, y=0.05; ridge nodes
    # fan into triangles on both slopes. Must not raise, and basis_y must
    # stay unit length everywhere (fallback or average, both normalized).
    nodes = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 0.05, 0.5],
            [1.0, 0.05, 0.5],
        ]
    )
    faces = np.array(
        [[0, 1, 5], [0, 5, 4], [4, 5, 3], [4, 3, 2]], dtype=np.int64
    )
    sm = SurfaceMesh.from_wall_faces(nodes, faces)
    np.testing.assert_allclose(
        np.linalg.norm(sm.basis_y, axis=1), 1.0, atol=1e-14
    )
    np.testing.assert_allclose(
        np.linalg.norm(sm.basis_x, axis=1), 1.0, atol=1e-14
    )
    np.testing.assert_allclose(
        np.linalg.norm(sm.basis_z, axis=1), 1.0, atol=1e-14
    )
    # Ridge normals tilt away from +y but stay well-defined and orthogonal.
    assert np.abs(np.sum(sm.basis_x * sm.basis_y, axis=1)).max() < 1e-14


def test_basis_sharp_fold_fallback_path():
    # Extreme fold: two triangles back-to-back (near-antiparallel normals)
    # sharing both nodes of an edge => the area-weighted normal sum cancels
    # below the crease threshold at the shared nodes; fallback must kick in
    # (no raise, unit basis_y).
    nodes = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, -1.0, 1e-3],
        ]
    )
    faces = np.array([[0, 1, 2], [1, 0, 3]], dtype=np.int64)
    sm = SurfaceMesh.from_wall_faces(nodes, faces)
    np.testing.assert_allclose(
        np.linalg.norm(sm.basis_y, axis=1), 1.0, atol=1e-14
    )


# ---------------------------------------------------------------------------
# 6. Coloring
# ---------------------------------------------------------------------------


def _check_coloring_valid(sm):
    # No two same-color triangles share a node: per-node incident-triangle
    # color lists must have distinct entries.
    color_of = np.empty(sm.n_tri, dtype=np.int64)
    for c in range(sm.n_colors):
        seg = sm.color_elems[sm.color_offsets[c] : sm.color_offsets[c + 1]]
        color_of[seg] = c
    node_tris = [[] for _ in range(sm.n_node)]
    for t in range(sm.n_tri):
        for j in range(3):
            node_tris[sm.triangles[t, j]].append(t)
    for nd, tl in enumerate(node_tris):
        cs = [color_of[t] for t in tl]
        assert len(cs) == len(set(cs)), f"node {nd}: same-color neighbours"


def test_coloring_valid_and_deterministic():
    sm = _plate(8, 6)
    _check_coloring_valid(sm)
    # Partition covers all triangles exactly once, ascending within a color.
    np.testing.assert_array_equal(
        np.sort(sm.color_elems), np.arange(sm.n_tri)
    )
    for c in range(sm.n_colors):
        seg = sm.color_elems[sm.color_offsets[c] : sm.color_offsets[c + 1]]
        assert np.all(np.diff(seg) > 0)
    assert sm.color_offsets[0] == 0 and sm.color_offsets[-1] == sm.n_tri
    # Two builds bit-identical.
    sm2 = _plate(8, 6)
    np.testing.assert_array_equal(sm.color_offsets, sm2.color_offsets)
    np.testing.assert_array_equal(sm.color_elems, sm2.color_elems)
    # Structured plate: small color count (256 cap is the kernel guard).
    assert 2 <= sm.n_colors <= 8


# ---------------------------------------------------------------------------
# 7. Jacobian pattern (block 6 CSR + elem_to_csr)
# ---------------------------------------------------------------------------


def test_jacobian_pattern():
    sm = _plate(5, 4)
    pattern, elem_to_csr = sm.build_jacobian_pattern()
    n = sm.n_node
    assert pattern.shape == (6 * n, 6 * n)
    assert elem_to_csr.shape == (sm.n_tri, 3, 3, 6, 6)
    # Canonical CSR.
    assert pattern.has_sorted_indices and pattern.has_canonical_format
    assert np.all(pattern.data == 0.0)

    # Symmetric STRUCTURE: every (r, c) in the pattern has (c, r) present.
    coo = pattern.tocoo()
    pairs = set(zip(coo.row.tolist(), coo.col.tolist()))
    for r, c in pairs:
        assert (c, r) in pairs

    # Diagonal blocks present for every node (stored-entry set, not
    # nonzero(): the symbolic pattern's data is zeroed by construction).
    for A in range(n):
        for p in range(6):
            assert (6 * A + p, 6 * A + p) in pairs

    # nnz == 36 * (number of unique node pairs sharing a triangle).
    node_pairs = set()
    for t in range(sm.n_tri):
        for a in range(3):
            for b in range(3):
                node_pairs.add((sm.triangles[t, a], sm.triangles[t, b]))
    assert pattern.nnz == 36 * len(node_pairs)

    # elem_to_csr points at the correct (row, col) = (6A+p, 6B+q): explicit
    # index math via the CSR arrays, for every element/local pair (spot the
    # full map -- F is small here).
    indptr, indices = pattern.indptr, pattern.indices
    for e in range(sm.n_tri):
        for a in range(3):
            A = sm.triangles[e, a]
            for b in range(3):
                B = sm.triangles[e, b]
                for p in range(6):
                    row = 6 * A + p
                    seg = indices[indptr[row] : indptr[row + 1]]
                    for q in range(6):
                        col = 6 * B + q
                        expect = indptr[row] + int(np.searchsorted(seg, col))
                        idx = elem_to_csr[e, a, b, p, q]
                        assert idx == expect
                        assert seg[idx - indptr[row]] == col

    # Deterministic: second build identical.
    pattern2, elem2 = sm.build_jacobian_pattern()
    np.testing.assert_array_equal(pattern.indptr, pattern2.indptr)
    np.testing.assert_array_equal(pattern.indices, pattern2.indices)
    np.testing.assert_array_equal(elem_to_csr, elem2)

    # Colored scatter is conflict-free: within one color no two elements
    # touch the same CSR entry, so a colored assembly writing
    # data[elem_to_csr[e, a, b, p, q]] += ... is race-free.
    for c in range(sm.n_colors):
        seg = sm.color_elems[sm.color_offsets[c] : sm.color_offsets[c + 1]]
        used = elem_to_csr[seg].reshape(-1)
        assert len(used) == len(np.unique(used))


# ---------------------------------------------------------------------------
# 8. Cross-check against existing infra (volume cube + post pipeline)
# ---------------------------------------------------------------------------


def test_cross_check_with_post_pipeline():
    n = 4
    nodes, elements = generate_structured_cube_mesh(n, L=1.0)
    # Wall = the y=0 face of the cube: volume-tet faces whose 3 nodes all
    # have y == 0 (both triangles of each structured quad are present since
    # the Kuhn split is face-conforming).
    wall_faces = []
    tet_faces = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    seen = set()
    for e in range(len(elements)):
        tet = elements[e]
        for f in tet_faces:
            tri = tuple(int(tet[i]) for i in f)
            if all(abs(nodes[v][1]) < 1e-12 for v in tri):
                key = tuple(sorted(tri))
                if key not in seen:
                    seen.add(key)
                    wall_faces.append(tri)
    wall_faces = np.asarray(wall_faces, dtype=np.int64)
    assert len(wall_faces) == 2 * n * n

    sm = SurfaceMesh.from_wall_faces(nodes, wall_faces, elements=elements)

    # Areas agree with post/surface.py::triangle_tangential_gradients.
    _, area_post, _ = triangle_tangential_gradients(
        nodes, wall_faces, np.zeros(len(nodes))
    )
    np.testing.assert_allclose(sm.area_tri, area_post, atol=1e-15)

    # Adjacency agrees with wall_triangle_adjacency on the same input
    # (compact renumbering is a bijection, face order preserved).
    np.testing.assert_array_equal(
        sm.adjacency, wall_triangle_adjacency(wall_faces)
    )

    # No boundary on the closed-face patch interior: exactly the outer rim
    # of the y=0 face is boundary => 4*n boundary edges.
    assert len(sm.boundary_edges) == 4 * n


def test_quad_tables():
    sm = _plate(3, 2)
    # QUAD_N is the edge-midpoint shape table; weights are area/3 per point.
    np.testing.assert_array_equal(
        QUAD_N, np.array([[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]])
    )
    np.testing.assert_allclose(sm.quad_w, sm.area_tri / 3.0, atol=0.0)
    # The 3-point rule integrates the triangle area exactly:
    # sum_q w_q * sum_i N_i(q) = 3 * (area/3) = area  (sum_i N_i = 1 at
    # each Gauss point).
    np.testing.assert_allclose(
        (sm.quad_w[:, None, None] * QUAD_N[None, :, :]).sum(),
        sm.area_tri.sum(),
        atol=1e-15,
    )
    # node_area lumps to the total area.
    np.testing.assert_allclose(sm.node_area.sum(), sm.area_tri.sum(), atol=1e-15)


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
