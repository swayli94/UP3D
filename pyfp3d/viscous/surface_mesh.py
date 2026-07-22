"""
Compact wall-surface mesh for the Track V IBL3 integral boundary layer
(design_track_v.md §5.1/§5.2 discretization, §6 data-layout points).

``SurfaceMesh`` builds the compact surface DOF numbering that the rest of
``pyfp3d/viscous`` runs on: wall triangulations elsewhere in the repo
(``mesh.boundary_faces["wall"]``) index VOLUME nodes with no compact surface
numbering; this module creates it. On top of the numbering it precomputes the
P1 surface-FE geometry tables the IBL3 Galerkin kernels consume as SoA flat
arrays (design.md §7 rule 1):

- per-triangle area / unit normal / P1 shape gradients in GLOBAL XYZ;
- per-node lumped area and local Cartesian basis (design_track_v.md §2.1:
  ŷ = area-weighted vertex normal, x̂ = Gram-Schmidt against a global seed,
  ẑ = ŷ×x̂; residual equations are assembled in this local basis, which makes
  the scheme invariant under in-surface rotation so the TE kink needs no
  special equations);
- 3-point edge-midpoint quadrature tables (project precedent
  solve/wall_correction.py: exact for quadratics, N = 1/2 on the edge nodes);
- greedy 3-node element coloring (mesh/coloring.py algorithm generalized to
  ``elements.shape[1]`` nodes per element -- serial color loop + prange
  within a color gives bit-deterministic scatter order, design_track_v.md
  §5.2);
- symbolic CSR sparsity of the 6-dof/node Newton system plus an
  ``elem_to_csr`` scatter map (kernels/jacobian.py::build_csr_pattern /
  build_elem_to_csr idiom, block 6x6 per node pair).

Group awareness (design_track_v.md §6 point 1): the node table is
group-aware by construction -- one ``SurfaceMesh`` instance per boundary
group, carrying the group ``name``. V1 builds the single "wall" group only;
the V6 wake sheet reuses the same layout with wake closures. No multi-group
machinery is built here.

Master-map hook (design_track_v.md §6 point 2): ``volume_node_of`` keeps the
surface-id -> volume-node-id map so an IBL surface mesh built on the UNCUT
wall can later pull u_e from a cut-mesh (LS) solution. Keep this name.

Orientation: ``normal_tri`` follows the triangle winding AS STORED -- no
reorientation is attempted (the standalone V1 use has no volume mesh to
orient against). Outward orientation of the input faces is the caller's
responsibility.

Everything here is one-time preprocessing: plain NumPy everywhere except the
coloring and CSR scatter-map kernels, which are ``_njit`` (PYFP3D_NOJIT=1
falls back to pure Python, same numerics).
"""

import os

import numba
import numpy as np
import scipy.sparse as sp

from pyfp3d.post.surface import wall_triangle_adjacency

if os.environ.get("PYFP3D_NOJIT", "0") == "1":
    prange = range

    def _njit(*args, **kwargs):
        def _identity(func):
            return func

        return _identity
else:
    from numba import prange

    def _njit(*args, **kwargs):
        return numba.njit(*args, **kwargs)


# 3-point edge-midpoint quadrature on a triangle: exact for quadratics
# (solve/wall_correction.py precedent). Rows are the P1 shape-function
# values N_i at the midpoint of edge (0,1), (1,2), (2,0) respectively --
# the Gauss-point coordinates are uniform per triangle; only the weights
# (area/3 per point) vary, stored per-triangle as SurfaceMesh.quad_w.
QUAD_N = np.array(
    [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]], dtype=np.float64
)

# Crease fallback for the vertex normal (design_track_v.md §2.1 ŷ): if the
# area-weighted normal sum nearly cancels at a node (sharp fold), fall back
# to the unit normal of the largest-area incident triangle instead of
# raising -- the TE-kink case must survive.
_CREASE_RATIO_MIN = 0.05


def structured_rectangle_surface(x0, x1, z0, z1, nx, nz):
    """Structured rectangular surface in the y=0 plane (standalone V1
    verification mesh: no Gmsh, no volume mesh -- the tests/mesh_utils.py
    precedent of dependency-free synthetic meshes).

    Nodes are row-major in x then z: node id = iz*(nx+1) + ix at
    (x0 + ix*dx, 0, z0 + iz*dz). Each quad is split into 2 triangles along
    the SAME diagonal everywhere (SW corner to NE corner), wound CCW when
    viewed from +y so every triangle's geometric normal is +ŷ.

    Returns:
        (xyz, triangles): xyz ((nx+1)*(nz+1), 3) float64,
        triangles (2*nx*nz, 3) int64.
    """
    xs = np.linspace(x0, x1, nx + 1)
    zs = np.linspace(z0, z1, nz + 1)
    X, Z = np.meshgrid(xs, zs, indexing="xy")  # X[ix along cols], row = iz
    xyz = np.stack(
        [X.ravel(), np.zeros(X.size), Z.ravel()], axis=1
    ).astype(np.float64)

    def nid(ix, iz):
        return iz * (nx + 1) + ix

    triangles = np.empty((2 * nx * nz, 3), dtype=np.int64)
    t = 0
    for iz in range(nz):
        for ix in range(nx):
            a = nid(ix, iz)  # SW
            b = nid(ix + 1, iz)  # SE
            c = nid(ix + 1, iz + 1)  # NE
            d = nid(ix, iz + 1)  # NW
            # Diagonal a->c; (a, d, c) and (a, c, b) both wind CCW from +y.
            triangles[t] = (a, d, c)
            triangles[t + 1] = (a, c, b)
            t += 2
    return np.ascontiguousarray(xyz), triangles


def _node_to_element_csr(elements, n_nodes):
    """CSR node -> incident elements (offsets, element ids), generalized
    mesh/coloring.py::node_to_element_csr to elements.shape[1] nodes/elem."""
    n_per = elements.shape[1]
    flat_nodes = np.asarray(elements, dtype=np.int64).reshape(-1)
    flat_elems = np.repeat(np.arange(len(elements), dtype=np.int64), n_per)
    order = np.argsort(flat_nodes, kind="stable")
    counts = np.bincount(flat_nodes, minlength=n_nodes)
    offsets = np.zeros(n_nodes + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    return offsets, flat_elems[order]


@_njit(cache=True)
def _greedy_coloring_kernel(elements, node_offsets, node_elems):
    """mesh/coloring.py::_greedy_coloring_kernel generalized to
    elements.shape[1] nodes per element. Same visit order (elements in
    order, smallest available color) => identical deterministic assignment;
    serial loop, 256-color cap."""
    n_elems = len(elements)
    n_per = elements.shape[1]
    colors = np.full(n_elems, -1, dtype=np.int32)
    max_colors = 256
    mark = np.full(max_colors, -1, dtype=np.int64)

    for e in range(n_elems):
        for i in range(n_per):
            node = elements[e, i]
            for k in range(node_offsets[node], node_offsets[node + 1]):
                c = colors[node_elems[k]]
                if c >= 0:
                    mark[c] = e
        c = 0
        while c < max_colors and mark[c] == e:
            c += 1
        if c >= max_colors:
            raise ValueError(
                "greedy coloring exceeded 256 colors -- pathological mesh "
                "connectivity (a node shared by >255 elements)"
            )
        colors[e] = c
    return colors


def _color_partition_csr(elements, n_nodes):
    """(color_offsets, color_elems, n_colors) with the same output contract
    as mesh/coloring.py::color_partition_csr: elements of color c are
    color_elems[color_offsets[c]:color_offsets[c+1]] in ascending element
    order. A serial color loop + prange within a color is race-free (no two
    same-color elements share a node) and the scatter order is fixed by the
    color sequence alone => bit-deterministic assembly across runs and
    thread counts (design_track_v.md §5.2)."""
    node_offsets, node_elems = _node_to_element_csr(elements, n_nodes)
    colors = _greedy_coloring_kernel(
        np.ascontiguousarray(elements), node_offsets, node_elems
    )
    n_colors = int(np.max(colors) + 1)
    order = np.argsort(colors, kind="stable")
    counts = np.bincount(colors, minlength=n_colors)
    color_offsets = np.zeros(n_colors + 1, dtype=np.int64)
    np.cumsum(counts, out=color_offsets[1:])
    return color_offsets, order.astype(np.int64), n_colors


@_njit(cache=True)
def _build_elem_to_csr_block6(triangles, indptr, indices):
    """elem_to_csr[e, a, b, p, q] -> flat data index of
    (row = 6*triangles[e,a] + p, col = 6*triangles[e,b] + q) in the
    canonical CSR pattern (kernels/jacobian.py::build_elem_to_csr idiom,
    expanded to the 6x6 dof block per node pair). Binary search per entry;
    built once, one-time cost."""
    n_tris = len(triangles)
    emap = np.empty((n_tris, 3, 3, 6, 6), dtype=np.int64)
    for e in range(n_tris):
        for a in range(3):
            row_node = triangles[e, a]
            for p in range(6):
                row = 6 * row_node + p
                lo = indptr[row]
                hi = indptr[row + 1]
                for b in range(3):
                    col_node = triangles[e, b]
                    for q in range(6):
                        col = 6 * col_node + q
                        x = lo
                        y = hi
                        while x < y:
                            mid = (x + y) // 2
                            if indices[mid] < col:
                                x = mid + 1
                            else:
                                y = mid
                        emap[e, a, b, p, q] = x
    return emap


class SurfaceMesh:
    """Compact wall-surface DOF numbering + P1 surface-FE geometry for IBL3.

    All attributes are SoA, C-contiguous, int64/float64. Plain container:
    arrays are set in ``__init__`` by ``from_wall_faces`` (the only
    supported constructor path); no lazy properties.

    Attributes:
        name (str): boundary group name. V1 is single-group ("wall"); the
            V6 wake sheet reuses this layout as a second instance
            (design_track_v.md §6 point 1).
        n_node (int): number of compact surface nodes n_s.
        n_tri (int): number of surface triangles F.
        volume_node_of (n_s,) int64: surface id -> volume node id, sorted
            unique. MASTER-MAP HOOK (design_track_v.md §6 point 2): keeps
            the link to the volume mesh so a later coupling stage can pull
            u_e from a (cut-mesh) volume solution. Keep this name.
        node_map (N,) int64: volume node id -> surface id, or -1 for volume
            nodes not on this surface (N = len(nodes) passed to
            from_wall_faces).
        triangles (F,3) int64: compact connectivity (indices in [0, n_s)).
        xyz (n_s,3) f64: surface node coordinates.
        area_tri (F,) f64: triangle areas.
        normal_tri (F,3) f64: geometric unit normals from the triangle
            winding AS STORED (no reorientation; outward orientation is the
            caller's responsibility -- the standalone V1 use has no volume
            mesh to orient against).
        adjacency (F,3) int64: edge-neighbour triangle per edge
            (edges (0,1),(1,2),(2,0)), -1 on boundary edges; via
            post/surface.py::wall_triangle_adjacency on the compact
            connectivity.
        boundary_edges (n_be,2) int64: compact node pairs of boundary edges
            (from the -1 adjacency entries, in (triangle, edge) visit
            order).
        boundary_node_mask (n_s,) bool: nodes touched by a boundary edge.
        node_area (n_s,) f64: lumped nodal area (area/3 scattered per
            vertex; volume mass-lumping idiom np.add.at, solve/picard.py).
        basis_x, basis_y, basis_z (n_s,3) f64: per-node local Cartesian
            basis (design_track_v.md §2.1). basis_y = area-weighted vertex
            normal, normalized; crease fallback to the largest-area incident
            triangle's unit normal when |sum(a*n̂)|/sum(a) < 0.05 (no raise:
            the TE-kink case must survive). basis_x = normalize(seed -
            (seed.basis_y) basis_y), seed = global X̂ (Ŷ when
            |basis_y.X̂| > 0.9). basis_z = basis_y x basis_x. Orthonormal to
            round-off.
        gradN (F,3,3) f64: per-triangle P1 shape gradients in GLOBAL XYZ,
            gradN[e, i] = n̂ x (x_{i+2} - x_{i+1}) / (2A) (cyclic indices mod
            3; n̂ the triangle unit normal, A its area).
        quad_w (F,) f64: per-point weight of the 3-point edge-midpoint rule
            (= area_tri/3); point shape values are the module-level QUAD_N.
        color_offsets (n_colors+1,) int64, color_elems (F,) int64,
        n_colors (int): greedy colored partition for race-free,
            bit-deterministic prange assembly (see _color_partition_csr).
    """

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return (
            f"SurfaceMesh(name={self.name!r}, n_node={self.n_node}, "
            f"n_tri={self.n_tri}, n_colors={self.n_colors}, "
            f"n_boundary_edges={len(self.boundary_edges)})"
        )

    @classmethod
    def from_wall_faces(cls, nodes, wall_faces, elements=None, name="wall"):
        """Build from a wall triangulation referencing VOLUME node indices.

        Args:
            nodes: (N,3) float64 volume/surface node coordinates.
            wall_faces: (F,3) int array of triangle corners as indices into
                nodes (e.g. mesh.boundary_faces["wall"]). Triangle order is
                preserved; winding is used as stored.
            elements: optional volume connectivity, accepted for call-site
                symmetry with the volume pipeline and future cut-mesh
                master-map validation; V1 does not use it.
            name: boundary group name (V1: single group "wall").

        Returns:
            SurfaceMesh with all geometry/assembly tables precomputed.
        """
        nodes = np.ascontiguousarray(nodes, dtype=np.float64)
        wf = np.ascontiguousarray(wall_faces, dtype=np.int64)
        if wf.ndim != 2 or wf.shape[1] != 3:
            raise ValueError("wall_faces must be (F, 3)")

        self = cls(name)
        # Compact numbering: sorted-unique volume ids -> [0, n_s).
        volume_node_of = np.unique(wf)  # sorted unique
        n_surf = len(volume_node_of)
        node_map = np.full(len(nodes), -1, dtype=np.int64)
        node_map[volume_node_of] = np.arange(n_surf, dtype=np.int64)
        triangles = node_map[wf]

        self.n_node = n_surf
        self.n_tri = len(wf)
        self.volume_node_of = np.ascontiguousarray(volume_node_of)
        self.node_map = node_map
        self.triangles = np.ascontiguousarray(triangles)
        self.xyz = np.ascontiguousarray(nodes[volume_node_of])

        # --- per-triangle geometry ---------------------------------------
        x0 = self.xyz[triangles[:, 0]]
        x1 = self.xyz[triangles[:, 1]]
        x2 = self.xyz[triangles[:, 2]]
        area_vec = np.cross(x1 - x0, x2 - x0)
        twice_area = np.linalg.norm(area_vec, axis=1)
        if np.any(twice_area <= 0.0):
            raise ValueError("degenerate (zero-area) triangle in wall_faces")
        self.area_tri = np.ascontiguousarray(0.5 * twice_area)
        self.normal_tri = np.ascontiguousarray(area_vec / twice_area[:, None])

        # --- adjacency / boundary ----------------------------------------
        self.adjacency = wall_triangle_adjacency(self.triangles)
        tri_edges = ((0, 1), (1, 2), (2, 0))
        bedges = []
        for t in range(self.n_tri):
            for k, (ea, eb) in enumerate(tri_edges):
                if self.adjacency[t, k] < 0:
                    bedges.append((triangles[t, ea], triangles[t, eb]))
        self.boundary_edges = np.ascontiguousarray(
            np.asarray(bedges, dtype=np.int64).reshape(-1, 2)
        )
        boundary_node_mask = np.zeros(n_surf, dtype=bool)
        if len(bedges):
            boundary_node_mask[self.boundary_edges.reshape(-1)] = True
        self.boundary_node_mask = boundary_node_mask

        # --- lumped nodal area (mass-lumping idiom, solve/picard.py) -----
        node_area = np.zeros(n_surf, dtype=np.float64)
        np.add.at(
            node_area,
            triangles.reshape(-1),
            np.repeat(self.area_tri / 3.0, 3),
        )
        self.node_area = node_area

        # --- per-node local basis (design_track_v.md §2.1) ---------------
        contrib = self.area_tri[:, None] * self.normal_tri  # (F,3)
        normal_sum = np.zeros((n_surf, 3), dtype=np.float64)
        weight_sum = np.zeros(n_surf, dtype=np.float64)
        np.add.at(normal_sum, triangles.reshape(-1), np.repeat(contrib, 3, axis=0))
        np.add.at(weight_sum, triangles.reshape(-1), np.repeat(self.area_tri, 3))

        # Largest-area incident triangle per node (crease fallback source).
        best_area = np.full(n_surf, -1.0, dtype=np.float64)
        best_tri = np.zeros(n_surf, dtype=np.int64)
        for t in range(self.n_tri):
            a = self.area_tri[t]
            for j in range(3):
                nd = triangles[t, j]
                if a > best_area[nd]:
                    best_area[nd] = a
                    best_tri[nd] = t

        basis_y = np.empty((n_surf, 3), dtype=np.float64)
        sum_norm = np.linalg.norm(normal_sum, axis=1)
        ratio = sum_norm / weight_sum
        for nd in range(n_surf):
            if ratio[nd] < _CREASE_RATIO_MIN:
                # Sharp fold: the area-weighted sum nearly cancels -- take
                # the largest-area incident triangle's unit normal instead
                # of raising (TE-kink case must survive).
                basis_y[nd] = self.normal_tri[best_tri[nd]]
            else:
                basis_y[nd] = normal_sum[nd] / sum_norm[nd]

        seed = np.tile(np.array([1.0, 0.0, 0.0]), (n_surf, 1))
        flip = np.abs(basis_y[:, 0]) > 0.9
        seed[flip] = np.array([0.0, 1.0, 0.0])
        bx = seed - np.sum(seed * basis_y, axis=1)[:, None] * basis_y
        basis_x = bx / np.linalg.norm(bx, axis=1)[:, None]
        basis_z = np.cross(basis_y, basis_x)
        self.basis_x = np.ascontiguousarray(basis_x)
        self.basis_y = np.ascontiguousarray(basis_y)
        self.basis_z = np.ascontiguousarray(basis_z)

        # --- P1 shape gradients in GLOBAL XYZ ----------------------------
        # gradN[e, i] = n̂ x (x_{i+2} - x_{i+1}) / (2A), cyclic mod 3.
        xt = (x0, x1, x2)
        gradN = np.empty((self.n_tri, 3, 3), dtype=np.float64)
        for i in range(3):
            ip1 = (i + 1) % 3
            ip2 = (i + 2) % 3
            gradN[:, i, :] = (
                np.cross(self.normal_tri, xt[ip2] - xt[ip1])
                / twice_area[:, None]
            )
        self.gradN = np.ascontiguousarray(gradN)

        # --- quadrature weights (points are the module-level QUAD_N) -----
        self.quad_w = np.ascontiguousarray(self.area_tri / 3.0)

        # --- coloring for bit-deterministic prange assembly --------------
        (
            self.color_offsets,
            self.color_elems,
            self.n_colors,
        ) = _color_partition_csr(self.triangles, self.n_node)
        return self

    def build_jacobian_pattern(self):
        """Symbolic CSR sparsity of the 6-dof/node Newton system.

        Shape (6*n_s, 6*n_s): every node pair (A, B) sharing a triangle is
        expanded to its full 6x6 block (equations p, q = 0..5). Canonical
        CSR (sorted indices, duplicates summed, data zeroed), built once --
        kernels/jacobian.py::build_csr_pattern idiom, block 6.

        Returns:
            (pattern_csr, elem_to_csr): pattern_csr is the scipy csr_matrix;
            elem_to_csr is (F,3,3,6,6) int64 with the flat data index for
            each (element, local node a, local node b, eq p, eq q), so the
            colored assembly kernel writes
            ``data[elem_to_csr[e, a, b, p, q]] += ...`` with no search.
        """
        tri = self.triangles
        n_dof = 6 * self.n_node
        # COO of all (6A+p, 6B+q) pairs per element; scipy merges the
        # duplicates across (a,b) repetitions and shared elements.
        rows_node = np.repeat(tri, 3, axis=1).reshape(-1)  # (F*9,)
        cols_node = np.tile(tri, (1, 3)).reshape(-1)  # (F*9,)
        dof6 = np.arange(6)
        # Per node pair: rows cycle p (each repeated 6x), cols cycle q --
        # the zip covers the full 6x6 block (all 36 (p, q) combos).
        rows = np.repeat(6 * rows_node[:, None] + dof6[None, :], 6, axis=1).reshape(-1)
        cols = np.tile(6 * cols_node[:, None] + dof6[None, :], (1, 6)).reshape(-1)
        pattern = sp.coo_matrix(
            (np.ones(len(rows), dtype=np.float64), (rows, cols)),
            shape=(n_dof, n_dof),
        ).tocsr()
        pattern.sum_duplicates()
        pattern.sort_indices()
        pattern.data[:] = 0.0
        elem_to_csr = _build_elem_to_csr_block6(
            tri, pattern.indptr.astype(np.int64), pattern.indices.astype(np.int64)
        )
        return pattern, elem_to_csr
