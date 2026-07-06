"""
M0 single-layer quasi-2D extrusion (docs/roadmap.md Track M, gate link G2.5).

Takes a 2D triangulation and extrudes exactly ONE cell layer in z between two
parallel planes z0 and z0+dz. Every 2D triangle becomes one prism, subdivided
into 3 tets with a globally consistent rule: on every quad face, the diagonal
passes through the smallest global node index of that face (the "indirect"
prism subdivision of Dompierre, Labbe, Vallet & Camarero, "How to Subdivide
Pyramids, Prisms and Hexahedra into Tetrahedra", IMR 1999). Consistency
argument specialized to this extrusion: top-layer node indices are
bottom + n_2d, so the smallest global index on any lateral quad face is
always the smaller *bottom* node of its 2D edge -- both prisms sharing the
face (and any tagged boundary/interior quad on it) therefore pick the same
diagonal, with no communication needed.

Boundary tagging (M0 spec):
  - both z-planes           -> "symmetry"
  - lateral surfaces        -> inherit their 2D edge-group names
                               ("wall", "farfield")
  - tagged *interior* edges -> interior faces ("wake"), split by the same
                               diagonal rule so they coincide exactly with
                               tet faces (no hanging nodes / slivers)

Node duplication for the wake is NOT done here -- it stays in the solver
preprocessor (mesh/wake_cut.py, P2) per agent-rules hard rule 8.
"""

from typing import Dict, Iterable, Optional, Tuple

import numpy as np

from pyfp3d.mesh.reader import Mesh


def _enforce_ccw(points2d: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Return triangles reordered so every signed area is positive (CCW in xy).

    The prism split below assumes the bottom triangle is CCW seen from +z
    (so the top layer sits on the positive side); Gmsh usually delivers
    that, but this makes the extruder independent of the source's
    orientation convention.
    """
    tri = np.asarray(triangles, dtype=np.int64).copy()
    p0 = points2d[tri[:, 0]]
    p1 = points2d[tri[:, 1]]
    p2 = points2d[tri[:, 2]]
    signed2 = (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (
        p1[:, 1] - p0[:, 1]
    ) * (p2[:, 0] - p0[:, 0])
    if np.any(signed2 == 0.0):
        raise ValueError("degenerate (zero-area) triangle in 2D input mesh")
    flip = signed2 < 0.0
    tri[flip] = tri[flip][:, [0, 2, 1]]
    return tri


def _split_prisms(triangles_ccw: np.ndarray, n2: int) -> np.ndarray:
    """Subdivide one prism per CCW triangle into 3 tets, min-index diagonal rule.

    Vertices of the prism over triangle (a, b, c): bottom (a, b, c) at z0,
    top (a+n2, b+n2, c+n2) at z0+dz. The prism is first rotated cyclically so
    its smallest global index sits at local vertex 0 (the global minimum of
    the 6 prism vertices is always a bottom vertex, since top = bottom + n2);
    the two quad faces touching vertex 0 then automatically take their
    diagonal through it, and the third quad face (v1, v2, v2', v1') takes the
    diagonal through min(v1, v2).
    """
    a, b, c = triangles_ccw[:, 0], triangles_ccw[:, 1], triangles_ccw[:, 2]
    rot = np.argmin(np.stack([a, b, c]), axis=0)
    v0 = np.choose(rot, [a, b, c])
    v1 = np.choose(rot, [b, c, a])
    v2 = np.choose(rot, [c, a, b])
    v3, v4, v5 = v0 + n2, v1 + n2, v2 + n2

    n_tri = len(triangles_ccw)
    tets = np.empty((3 * n_tri, 4), dtype=np.int64)
    diag_a = v1 < v2  # diagonal v1 -- v5 on quad (v1, v2, v5, v4); else v2 -- v4

    tets[0::3] = np.where(
        diag_a[:, None],
        np.stack([v0, v1, v2, v5], axis=1),
        np.stack([v0, v1, v2, v4], axis=1),
    )
    tets[1::3] = np.where(
        diag_a[:, None],
        np.stack([v0, v1, v5, v4], axis=1),
        np.stack([v0, v4, v2, v5], axis=1),
    )
    tets[2::3] = np.stack([v0, v4, v5, v3], axis=1)
    return tets


def _split_lateral_quads(edges: np.ndarray, n2: int) -> np.ndarray:
    """Two triangles per extruded 2D edge, same min-index diagonal rule.

    Quad cycle over edge (e0, e1): (e0, e1, e1', e0') with X' = X + n2. The
    diagonal passes through min(e0, e1); both splits below preserve the
    cycle orientation, so the resulting triangle winding is consistent
    whenever the input edges are consistently oriented along their curve.
    """
    e0 = np.asarray(edges, dtype=np.int64)[:, 0]
    e1 = np.asarray(edges, dtype=np.int64)[:, 1]
    tris = np.empty((2 * len(e0), 3), dtype=np.int64)
    lo_first = e0 < e1  # diagonal e0 -- e1'; else e1 -- e0'
    tris[0::2] = np.where(
        lo_first[:, None],
        np.stack([e0, e1, e1 + n2], axis=1),
        np.stack([e0, e1, e0 + n2], axis=1),
    )
    tris[1::2] = np.where(
        lo_first[:, None],
        np.stack([e0, e1 + n2, e0 + n2], axis=1),
        np.stack([e1, e1 + n2, e0 + n2], axis=1),
    )
    return tris


def extrude_single_layer(
    points2d: np.ndarray,
    triangles: np.ndarray,
    edge_groups: Dict[str, np.ndarray],
    interior_edge_groups: Optional[Dict[str, np.ndarray]] = None,
    dz: float = 0.1,
    z0: float = 0.0,
    name: str = "quasi2d",
) -> Mesh:
    """
    Extrude a 2D triangulation into a single-layer quasi-2D tet mesh.

    Args:
        points2d: (n2, 2) planar node coordinates
        triangles: (n_tri, 3) 2D connectivity (any orientation; made CCW here)
        edge_groups: named *boundary* edge sets, e.g. {"wall": (m, 2), ...};
            edges index into points2d and must lie on the 2D boundary
        interior_edge_groups: named *interior* edge sets that become tagged
            interior face sheets, e.g. {"wake": (k, 2)}
        dz: layer thickness (> 0)
        z0: z of the bottom plane
        name: mesh name

    Returns:
        pyfp3d Mesh: nodes (2*n2, 3), tets (3*n_tri, 4), boundary_faces with
        the lateral groups, "symmetry" (both z-planes), and any interior
        groups (e.g. "wake" -- tagged faces, nodes NOT duplicated).
    """
    points2d = np.asarray(points2d, dtype=np.float64)
    if points2d.ndim != 2 or points2d.shape[1] != 2:
        raise ValueError(f"points2d must be (n, 2), got {points2d.shape}")
    if dz <= 0.0:
        raise ValueError(f"dz must be positive, got {dz}")
    n2 = len(points2d)

    tri = _enforce_ccw(points2d, triangles)
    used = np.zeros(n2, dtype=bool)
    used[tri.ravel()] = True
    if not used.all():
        raise ValueError(
            f"{int((~used).sum())} node(s) of points2d are referenced by no "
            "triangle -- compact the 2D mesh before extruding"
        )

    nodes = np.empty((2 * n2, 3), dtype=np.float64)
    nodes[:n2, :2] = points2d
    nodes[:n2, 2] = z0
    nodes[n2:, :2] = points2d
    nodes[n2:, 2] = z0 + dz

    tets = _split_prisms(tri, n2)

    boundary_faces: Dict[str, np.ndarray] = {}
    for tag, edges in edge_groups.items():
        boundary_faces[tag] = np.ascontiguousarray(
            _split_lateral_quads(edges, n2), dtype=np.int32
        )
    # Both z-planes in one "symmetry" group, wound outward (-z bottom, +z top).
    bottom = tri[:, [0, 2, 1]]
    top = tri + n2
    boundary_faces["symmetry"] = np.ascontiguousarray(
        np.vstack([bottom, top]), dtype=np.int32
    )
    for tag, edges in (interior_edge_groups or {}).items():
        boundary_faces[tag] = np.ascontiguousarray(
            _split_lateral_quads(edges, n2), dtype=np.int32
        )

    mesh = Mesh()
    mesh.nodes = nodes
    mesh.elements = np.ascontiguousarray(tets, dtype=np.int32)
    mesh.boundary_faces = boundary_faces
    mesh.element_tags = np.zeros(len(tets), dtype=np.int32)
    mesh.tag_names = ["fluid"]
    mesh.name = name
    mesh.validate()
    assert_quad_split_consistency(
        mesh, interior_groups=tuple(interior_edge_groups or ())
    )
    return mesh


def _face_counts(elements: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Sorted node triples of all tet faces and their occurrence counts."""
    el = np.asarray(elements, dtype=np.int64)
    faces = np.concatenate(
        [el[:, [1, 2, 3]], el[:, [0, 2, 3]], el[:, [0, 1, 3]], el[:, [0, 1, 2]]]
    )
    faces.sort(axis=1)
    return np.unique(faces, axis=0, return_counts=True)


def assert_quad_split_consistency(
    mesh: Mesh, interior_groups: Iterable[str] = ("wake",)
) -> None:
    """
    M0 preprocessor assert: every interior quad face was split identically
    from both sides (no hanging nodes / slivers), and every tagged face
    coincides with a tet face.

    Checks, on the tet mesh's face multiset:
      1. faces owned by exactly one tet  ==  the union of the *boundary*
         face groups (set equality, both directions). An inconsistently
         split interior quad would leave its 2x2 mismatched triangles
         single-owned in the interior, breaking this equality.
      2. every face in an *interior* group (e.g. "wake") is owned by exactly
         two tets -- i.e. the tagged sheet coincides with conforming
         interior tet faces.

    Raises AssertionError with a diagnostic on any violation.
    """
    interior_groups = set(interior_groups)
    unique_faces, counts = _face_counts(mesh.elements)

    def _key_view(arr: np.ndarray) -> np.ndarray:
        a = np.ascontiguousarray(np.sort(np.asarray(arr, dtype=np.int64), axis=1))
        return a.view([("", a.dtype)] * 3).ravel()

    face_keys = _key_view(unique_faces)
    order = np.argsort(face_keys)
    face_keys_sorted = face_keys[order]
    counts_sorted = counts[order]

    def _counts_for(tris: np.ndarray) -> np.ndarray:
        keys = _key_view(tris)
        pos = np.searchsorted(face_keys_sorted, keys)
        missing = (pos >= len(face_keys_sorted)) | (face_keys_sorted[np.minimum(pos, len(face_keys_sorted) - 1)] != keys)
        if np.any(missing):
            raise AssertionError(
                f"{int(missing.sum())} tagged face(s) are not faces of any tet "
                "(quad split mismatch between tagging and subdivision)"
            )
        return counts_sorted[pos]

    boundary_keys = []
    for tag, tris in mesh.boundary_faces.items():
        c = _counts_for(tris)
        if tag in interior_groups:
            if not np.all(c == 2):
                raise AssertionError(
                    f"interior group '{tag}': {int(np.sum(c != 2))} face(s) not "
                    "shared by exactly two tets (hanging node / sliver at the sheet)"
                )
        else:
            if not np.all(c == 1):
                raise AssertionError(
                    f"boundary group '{tag}': {int(np.sum(c != 1))} face(s) not "
                    "owned by exactly one tet"
                )
            boundary_keys.append(_key_view(tris))

    once_keys = face_keys_sorted[counts_sorted == 1]
    tagged_once = np.sort(np.concatenate(boundary_keys)) if boundary_keys else np.empty(0, once_keys.dtype)
    if len(tagged_once) != len(np.unique(tagged_once)):
        raise AssertionError("duplicate faces across boundary groups")
    if len(once_keys) != len(tagged_once) or not np.array_equal(once_keys, tagged_once):
        raise AssertionError(
            f"single-owner tet faces ({len(once_keys)}) != tagged boundary faces "
            f"({len(tagged_once)}): an interior quad was split inconsistently "
            "from its two sides, or the boundary tagging is incomplete"
        )
