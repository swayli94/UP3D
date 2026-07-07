"""
Wake-cut preprocessor: duplicate the nodes of the tagged interior "wake"
sheet so the potential can be multivalued around a lifting body
(design.md Sec 4; roadmap P2 deliverable).

Conventions (fixed here, relied on by constraints/wake.py and the Kutta
loop):

  - The sheet's two sides are labelled by `upper_hint` (default +y): the
    element star on the hint side of each wake face is the "+" (upper)
    side, the other the "-" (lower) side.
  - "+"-side copies are NEW nodes appended after the original ones
    (slaves); "-"-side elements keep the ORIGINAL node ids (masters).
    The jump is [phi] = phi(+) - phi(-) = Gamma  (design.md (4.1), (4.3)).
  - TE nodes (wall AND wake) ARE duplicated: in the continuum the jump at
    the trailing edge itself equals Gamma -- that is the Kutta condition
    -- so the sheet's jump must reach the wall. The originally specified
    single-valued TE (roadmap P2 topology asserts, pre-2026-07-06) was
    tried first and is QUANTITATIVELY WRONG: tapering the jump from Gamma
    to 0 over the first wake cell is equivalent to parking a point vortex
    of strength Gamma at the TE, whose wall suction spike integrates to a
    spurious force ~ Gamma^2/h that DIVERGES under refinement (measured on
    the coarse NACA0012 mesh: -0.27 of a 0.6 expected cl from the 6
    TE-adjacent triangles, peak |V| 4.6 U_inf). See roadmap P2 evidence
    note.
  - Wake nodes on any other boundary (symmetry planes, far field) ARE
    duplicated too; their boundary faces (including the wall faces at the
    TE) are re-pointed to the copy on the "+" side. A duplicated far-field
    node needs no special Dirichlet handling: the slave is eliminated to
    master + Gamma, which matches the vortex far-field correction's
    branch-cut jump across the wake (constraints/dirichlet.py places the
    cut on the wake sheet).
  - Sheet FREE edges (sheet boundary edges in the domain INTERIOR -- the
    M1 wing-tip edge, running from the tip TE corner downstream) are NOT
    duplicated: their element stars are not separated by the sheet, and
    physically the jump must vanish there (Gamma -> 0 at the tip; a
    trailing vortex line of the residual discrete Gamma_tip is the
    expected tip-vortex singularity, unlike the TE case where tapering a
    FINITE Gamma was measured to diverge). The tip TE corner node itself
    is a free-edge node: it stays single-valued and is excluded from the
    TE/Kutta stations, pinning Gamma(tip) = 0 discretely. Quasi-2D M0
    sheets have no free edges (both sheet sides lie on symmetry planes),
    so this path is exactly inert there.

Side classification is a per-node flood fill over the node's incident
elements (adjacency through non-wake faces only), seeded geometrically per
wake face -- no global planarity assumption, so a swept M1 wake works the
same way. Node duplication happens only here, never in Gmsh (agent-rules
hard rule 8).
"""

from typing import Dict, Tuple

import numpy as np

from pyfp3d.mesh.reader import Mesh


def _face_key_view(tris: np.ndarray) -> np.ndarray:
    """Sorted node triples as a structured view usable as dict/set keys."""
    a = np.ascontiguousarray(np.sort(np.asarray(tris, dtype=np.int64), axis=1))
    return a.view([("", a.dtype)] * 3).ravel()


_TET_FACES = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))


def _all_tet_faces(elements: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """(4*n_tets,) face keys and owning element ids, in tet-face order."""
    el = np.asarray(elements, dtype=np.int64)
    faces = np.concatenate([el[:, list(f)] for f in _TET_FACES])
    owners = np.tile(np.arange(len(el), dtype=np.int64), 4)
    return _face_key_view(faces), owners


def _node_to_elem_csr(elements: np.ndarray, n_nodes: int):
    flat_nodes = elements.reshape(-1).astype(np.int64)
    flat_elems = np.repeat(np.arange(len(elements), dtype=np.int64), 4)
    order = np.argsort(flat_nodes, kind="stable")
    counts = np.bincount(flat_nodes, minlength=n_nodes)
    offsets = np.zeros(n_nodes + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    return offsets, flat_elems[order]


class WakeCut:
    """Result of cut_wake(): duplication map, stations, Kutta probe nodes.

    Attributes:
        n_nodes_orig: node count before duplication
        master_nodes: (n_dup,) original ("-"-side) ids of duplicated nodes
        slave_nodes: (n_dup,) the new "+"-side ids (n_nodes_orig + k)
        free_nodes: (n_free,) sheet nodes on interior free edges (M1 tip
                  edge), kept single-valued (module docstring); empty on
                  quasi-2D sheets
        te_nodes: (n_te,) wall-and-wake nodes; duplicated like the rest of
                  the sheet (the jump reaches the wall -- module docstring)
        station_z: (n_st,) spanwise station coordinates (mean z of the
                      station's TE nodes). Stations group TE nodes by
                      their (x, y) position: on a quasi-2D extrusion all
                      TE nodes share one (x, y), so Gamma collapses to the
                      single scalar the M0 spec calls for; on a swept 3D
                      TE every node has its own (x, y) and stations are
                      per-node as design.md Sec 4 intends.
        node_station: (n_dup,) station index of each duplicated node
                      (wake lines inherit the Gamma of their TE station)
        te_station: (n_te,) station index of each TE node
        kutta_upper / kutta_lower: (n_te,) wall node ids one node off
                      EACH TE node on the +/- side; the Kutta target per
                      station averages phi[upper] - phi[lower] over the
                      station's TE nodes (constraints/wake.kutta_targets)
        wake_faces_minus / wake_faces_plus: (n_wake, 3) the sheet's two
                      copies in the cut mesh ("-": original ids, "+": with
                      duplicated ids); single-owner faces after the cut
    """

    def __init__(self):
        self.n_nodes_orig = 0
        self.master_nodes = np.empty(0, dtype=np.int64)
        self.slave_nodes = np.empty(0, dtype=np.int64)
        self.free_nodes = np.empty(0, dtype=np.int64)
        self.te_nodes = np.empty(0, dtype=np.int64)
        self.station_z = np.empty(0, dtype=np.float64)
        self.node_station = np.empty(0, dtype=np.int64)
        self.te_station = np.empty(0, dtype=np.int64)
        self.kutta_upper = np.empty(0, dtype=np.int64)
        self.kutta_lower = np.empty(0, dtype=np.int64)
        self.wake_faces_minus = np.empty((0, 3), dtype=np.int32)
        self.wake_faces_plus = np.empty((0, 3), dtype=np.int32)

    @property
    def n_stations(self) -> int:
        return len(self.station_z)

    def gamma_at_nodes(self, gamma_stations: np.ndarray) -> np.ndarray:
        """(n_dup,) jump value per duplicated node from per-station Gamma."""
        return np.asarray(gamma_stations, dtype=np.float64)[self.node_station]


def _edge_key_view(edges: np.ndarray) -> np.ndarray:
    """Sorted node pairs as a structured view usable for membership tests."""
    a = np.ascontiguousarray(np.sort(np.asarray(edges, dtype=np.int64), axis=1))
    return a.view([("", a.dtype)] * 2).ravel()


def _sheet_free_edge_nodes(mesh: Mesh, wake_faces: np.ndarray,
                           wake_tag: str) -> np.ndarray:
    """Nodes on the sheet's INTERIOR free edges (module docstring).

    A sheet boundary edge (used by exactly one wake face) whose edge is
    also an edge of some boundary triangle (wall / symmetry / far field)
    gets its node stars cut by that boundary; one that lies in the domain
    interior (M1 tip edge) does not, and its nodes must stay single-valued.
    """
    edges = np.concatenate([wake_faces[:, [0, 1]], wake_faces[:, [1, 2]],
                            wake_faces[:, [2, 0]]])
    keys = _edge_key_view(edges)
    uniq, counts = np.unique(keys, return_counts=True)
    sheet_boundary = uniq[counts == 1]
    if len(sheet_boundary) == 0:
        return np.empty(0, dtype=np.int64)

    bnd_edges = []
    for tag, tris in mesh.boundary_faces.items():
        if tag == wake_tag:
            continue
        t = np.asarray(tris, dtype=np.int64)
        bnd_edges.append(np.concatenate([t[:, [0, 1]], t[:, [1, 2]],
                                         t[:, [2, 0]]]))
    bnd_keys = (np.unique(_edge_key_view(np.concatenate(bnd_edges)))
                if bnd_edges else np.empty(0, dtype=keys.dtype))

    free = sheet_boundary[~np.isin(sheet_boundary, bnd_keys)]
    if len(free) == 0:
        return np.empty(0, dtype=np.int64)
    free_pairs = free.view(np.int64).reshape(-1, 2)
    return np.unique(free_pairs)


def _split_node_star(offsets, sorted_elems, elements64, node, wake_keys_set):
    """Partition the incident elements of `node` into the two sheet sides.

    Adjacency: two incident tets are connected iff they share a face that
    contains `node` and is not a wake face. Returns a list of components
    (each a list of element ids); a fully-cut star yields exactly 2.
    """
    elems = sorted_elems[offsets[node]:offsets[node + 1]]
    n_loc = len(elems)
    # The 3 faces of each incident tet that contain `node`.
    face_to_locals: Dict[bytes, list] = {}
    for li, e in enumerate(elems):
        tet = elements64[e]
        for f in _TET_FACES:
            tri = (tet[f[0]], tet[f[1]], tet[f[2]])
            if node not in tri:
                continue
            key = _face_key_view(np.array([tri]))[0].tobytes()
            if key in wake_keys_set:
                continue
            face_to_locals.setdefault(key, []).append(li)

    adj = [[] for _ in range(n_loc)]
    for locs in face_to_locals.values():
        if len(locs) == 2:
            adj[locs[0]].append(locs[1])
            adj[locs[1]].append(locs[0])

    comp = np.full(n_loc, -1, dtype=np.int64)
    n_comp = 0
    for start in range(n_loc):
        if comp[start] >= 0:
            continue
        stack = [start]
        comp[start] = n_comp
        while stack:
            cur = stack.pop()
            for nb in adj[cur]:
                if comp[nb] < 0:
                    comp[nb] = n_comp
                    stack.append(nb)
        n_comp += 1
    return elems, comp, n_comp


def cut_wake(
    mesh: Mesh,
    wake_tag: str = "wake",
    wall_tag: str = "wall",
    upper_hint=(0.0, 1.0, 0.0),
) -> Tuple[Mesh, WakeCut]:
    """Duplicate wake-sheet nodes and re-attach "+"-side elements/faces.

    Args:
        mesh: conforming mesh with an interior face group `wake_tag`
              (every wake face shared by exactly two tets) and a boundary
              group `wall_tag` (used to find TE nodes and Kutta probes)
        wake_tag, wall_tag: group names in mesh.boundary_faces
        upper_hint: direction whose side of the sheet is "+" (upper)

    Returns:
        (mesh_cut, wake_cut): a NEW Mesh (input untouched) whose
        boundary_faces has `wake_tag` replaced by "wake_minus"/"wake_plus",
        and the WakeCut map. Topology asserts (assert_wake_topology) are
        run before returning -- roadmap P2 requires them at preprocess
        time on every mesh.
    """
    if wake_tag not in mesh.boundary_faces:
        raise ValueError(f"mesh has no '{wake_tag}' face group")
    if wall_tag not in mesh.boundary_faces:
        raise ValueError(f"mesh has no '{wall_tag}' face group")

    nodes = mesh.nodes
    elements64 = np.asarray(mesh.elements, dtype=np.int64)
    n_nodes = len(nodes)
    hint = np.asarray(upper_hint, dtype=np.float64)
    hint /= np.linalg.norm(hint)

    wake_faces = np.asarray(mesh.boundary_faces[wake_tag], dtype=np.int64)
    wake_keys = _face_key_view(wake_faces)
    wake_keys_set = {k.tobytes() for k in wake_keys}

    # Owning tets of every face (wake faces must have exactly two).
    all_keys, all_owners = _all_tet_faces(elements64)
    order = np.argsort(all_keys)
    all_keys_sorted = all_keys[order]
    all_owners_sorted = all_owners[order]

    def _owners_of(tris: np.ndarray) -> list:
        keys = _face_key_view(tris)
        lo = np.searchsorted(all_keys_sorted, keys, side="left")
        hi = np.searchsorted(all_keys_sorted, keys, side="right")
        return [all_owners_sorted[a:b] for a, b in zip(lo, hi)]

    wake_owners = _owners_of(wake_faces)
    bad = [i for i, o in enumerate(wake_owners) if len(o) != 2]
    if bad:
        raise AssertionError(
            f"{len(bad)} wake face(s) not shared by exactly two tets "
            "(sheet does not conform to the volume mesh)"
        )

    # Geometric "+"-side owner per wake face (seed for the flood fill).
    centroids = nodes[elements64].mean(axis=1)
    face_centers = nodes[wake_faces].mean(axis=1)
    plus_owner = np.empty(len(wake_faces), dtype=np.int64)
    for i, (e0, e1) in enumerate(wake_owners):
        s0 = np.dot(centroids[e0] - face_centers[i], hint)
        s1 = np.dot(centroids[e1] - face_centers[i], hint)
        if s0 == s1:
            raise AssertionError(
                f"wake face {i}: cannot classify sides along upper_hint "
                f"{tuple(hint)} (both owners at offset {s0:g})"
            )
        plus_owner[i] = e0 if s0 > s1 else e1

    # Node sets. TE nodes (wall-and-wake) are duplicated ALONG WITH the
    # interior sheet nodes -- the jump reaches the wall (module docstring).
    # Free-edge nodes (M1 tip edge) are excluded from duplication AND from
    # the TE/Kutta stations: they stay single-valued, pinning the jump to
    # zero where the sheet ends inside the domain.
    wake_nodes = np.unique(wake_faces)
    wall_nodes = np.unique(np.asarray(mesh.boundary_faces[wall_tag], dtype=np.int64))
    is_wall = np.zeros(n_nodes, dtype=bool)
    is_wall[wall_nodes] = True
    free_nodes = _sheet_free_edge_nodes(mesh, wake_faces, wake_tag)
    is_free = np.zeros(n_nodes, dtype=bool)
    is_free[free_nodes] = True
    te_nodes = wake_nodes[is_wall[wake_nodes] & ~is_free[wake_nodes]]
    dup_nodes = wake_nodes[~is_free[wake_nodes]]
    if len(te_nodes) == 0:
        raise AssertionError("wake sheet does not touch the wall (no TE nodes)")

    # Spanwise stations: group TE nodes by (x, y) position (see the
    # WakeCut docstring -- a quasi-2D extrusion collapses to ONE station /
    # one scalar Gamma; a swept TE keeps one station per TE node). Wake
    # lines inherit the Gamma of the TE station they emanate from
    # (design.md Sec 4), assigned by nearest station z.
    z_extent = float(np.ptp(nodes[:, 2]))
    z_tol = 1e-9 * z_extent if z_extent > 0 else 1e-12
    bbox = float(np.max(np.ptp(nodes, axis=0)))
    xy_key = np.round(nodes[te_nodes, :2] / (1e-6 * bbox)).astype(np.int64)
    _, station_of_te = np.unique(xy_key, axis=0, return_inverse=True)
    n_st = int(station_of_te.max()) + 1
    te_station = station_of_te.astype(np.int64)
    station_z = np.array(
        [nodes[te_nodes[te_station == j], 2].mean() for j in range(n_st)]
    )

    if n_st == 1:
        node_station = np.zeros(len(dup_nodes), dtype=np.int64)
    else:
        z_vals = nodes[dup_nodes, 2]
        order_st = np.argsort(station_z)
        sz = station_z[order_st]
        idx = np.clip(np.searchsorted(sz, z_vals), 0, n_st - 1)
        left = np.clip(idx - 1, 0, n_st - 1)
        use_left = np.abs(z_vals - sz[left]) < np.abs(z_vals - sz[idx])
        node_station = order_st[np.where(use_left, left, idx)]

    # Per-node flood fill: split each duplicated node's element star.
    offsets, sorted_elems = _node_to_elem_csr(elements64, n_nodes)
    plus_of_face_by_node: Dict[int, Dict[int, int]] = {}
    for i, tri in enumerate(wake_faces):
        for nd in tri:
            plus_of_face_by_node.setdefault(int(nd), {})[i] = plus_owner[i]

    elements_cut = elements64.copy()
    slave_nodes = n_nodes + np.arange(len(dup_nodes), dtype=np.int64)
    slave_of = np.full(n_nodes, -1, dtype=np.int64)
    slave_of[dup_nodes] = slave_nodes
    plus_elems_of_node: Dict[int, np.ndarray] = {}

    for k, d in enumerate(dup_nodes):
        elems, comp, n_comp = _split_node_star(
            offsets, sorted_elems, elements64, int(d), wake_keys_set
        )
        if n_comp != 2:
            raise AssertionError(
                f"wake node {d}: element star splits into {n_comp} "
                "components (expected 2) -- sheet does not fully cut the "
                "star, or the mesh is non-manifold there"
            )
        # Which component is "+": every incident wake face names its
        # plus-side owner; all of them must agree on the same component.
        plus_comp = set()
        for fi, po in plus_of_face_by_node[int(d)].items():
            li = np.where(elems == po)[0]
            if len(li) != 1:
                raise AssertionError(
                    f"wake node {d}: plus-side owner of face {fi} is not "
                    "in the node's element star"
                )
            plus_comp.add(int(comp[li[0]]))
        if len(plus_comp) != 1:
            raise AssertionError(
                f"wake node {d}: inconsistent +/- side classification "
                "across its wake faces (check upper_hint vs sheet shape)"
            )
        pc = plus_comp.pop()
        plus_elems = elems[comp == pc]
        plus_elems_of_node[int(d)] = plus_elems
        for e in plus_elems:
            row = elements_cut[e]
            row[row == d] = slave_nodes[k]

    # Boundary faces: re-point faces owned by a "+"-side element.
    boundary_cut: Dict[str, np.ndarray] = {}
    for tag, tris in mesh.boundary_faces.items():
        if tag == wake_tag:
            continue
        tris64 = np.asarray(tris, dtype=np.int64).copy()
        touched = np.isin(tris64, dup_nodes).any(axis=1)
        if np.any(touched):
            owners = _owners_of(tris64[touched])
            for row_idx, own in zip(np.where(touched)[0], owners):
                if len(own) != 1:
                    raise AssertionError(
                        f"boundary face in '{tag}' owned by {len(own)} tets"
                    )
                e = own[0]
                tri = tris64[row_idx]
                for j in range(3):
                    d = tri[j]
                    if slave_of[d] >= 0 and e in plus_elems_of_node[int(d)]:
                        tri[j] = slave_of[d]
        boundary_cut[tag] = np.ascontiguousarray(tris64, dtype=np.int32)

    wake_minus = np.ascontiguousarray(wake_faces, dtype=np.int32)
    wake_plus64 = wake_faces.copy()
    mask = slave_of[wake_plus64] >= 0
    wake_plus64[mask] = slave_of[wake_plus64[mask]]
    wake_plus = np.ascontiguousarray(wake_plus64, dtype=np.int32)
    boundary_cut["wake_minus"] = wake_minus
    boundary_cut["wake_plus"] = wake_plus

    mesh_cut = Mesh()
    mesh_cut.nodes = np.vstack([nodes, nodes[dup_nodes]])
    mesh_cut.elements = np.ascontiguousarray(elements_cut, dtype=np.int32)
    mesh_cut.boundary_faces = boundary_cut
    mesh_cut.element_tags = mesh.element_tags.copy()
    mesh_cut.tag_names = list(mesh.tag_names)
    mesh_cut.name = f"{mesh.name}_cut"
    mesh_cut.validate()

    wc = WakeCut()
    wc.n_nodes_orig = n_nodes
    wc.master_nodes = dup_nodes
    wc.slave_nodes = slave_nodes
    wc.free_nodes = free_nodes
    wc.te_nodes = te_nodes
    wc.station_z = station_z
    wc.node_station = node_station
    wc.te_station = te_station
    wc.wake_faces_minus = wake_minus
    wc.wake_faces_plus = wake_plus
    wc.kutta_upper, wc.kutta_lower = _kutta_probe_nodes(
        mesh, wc, wall_nodes, wake_nodes, hint, z_tol
    )

    assert_wake_topology(mesh_cut, wc)
    return mesh_cut, wc


def _kutta_probe_nodes(mesh, wc, wall_nodes, wake_nodes, hint, z_tol):
    """Per TE node: the wall node one edge off the TE on each side.

    "One node off the TE" (design.md (4.4)): among the TE node's wall-edge
    neighbors, excluding wake/TE nodes, take the nearest node on the +hint
    side (upper) and on the -hint side (lower). On layered quasi-2D meshes
    candidates are first restricted to the TE node's own z-plane (the
    station's section plane); on an unstructured swept TE (M1) no wall
    neighbor shares the plane exactly, so a node whose strict pass comes
    up empty falls back to the unrestricted nearest +/- side neighbor --
    the probe is then off-plane by O(h), the same order as the "one node
    off the TE" approximation itself.
    """
    nodes = mesh.nodes
    wall_faces = np.asarray(mesh.boundary_faces["wall"], dtype=np.int64)
    is_wake = np.zeros(len(nodes), dtype=bool)
    is_wake[wake_nodes] = True

    neighbors: Dict[int, set] = {int(t): set() for t in wc.te_nodes}
    te_set = set(int(t) for t in wc.te_nodes)
    for tri in wall_faces:
        tri_i = [int(x) for x in tri]
        for a in tri_i:
            if a in te_set:
                for b in tri_i:
                    if b != a:
                        neighbors[a].add(b)

    n_te = len(wc.te_nodes)
    upper = np.full(n_te, -1, dtype=np.int64)
    lower = np.full(n_te, -1, dtype=np.int64)

    def _pick(t: int, same_plane: bool):
        up, lo = -1, -1
        d_up, d_lo = np.inf, np.inf
        for nb in neighbors[t]:
            if is_wake[nb] or nb in te_set:
                continue
            if same_plane and abs(nodes[nb, 2] - nodes[t, 2]) > 10 * z_tol + 1e-12:
                continue
            offset = nodes[nb] - nodes[t]
            side = float(np.dot(offset, hint))
            dist = float(np.linalg.norm(offset))
            if side > 0 and dist < d_up:
                up, d_up = nb, dist
            elif side < 0 and dist < d_lo:
                lo, d_lo = nb, dist
        return up, lo

    for k, t in enumerate(wc.te_nodes):
        up, lo = _pick(int(t), same_plane=True)
        if up < 0 or lo < 0:
            up, lo = _pick(int(t), same_plane=False)
        upper[k], lower[k] = up, lo

    if np.any(upper < 0) or np.any(lower < 0):
        missing = [int(wc.te_nodes[k]) for k in range(n_te)
                   if upper[k] < 0 or lower[k] < 0]
        raise AssertionError(
            f"TE nodes {missing}: no wall node found one edge off the TE "
            "on both sides (Kutta probes) -- check wall tagging / upper_hint"
        )
    return upper, lower


def assert_wake_topology(mesh_cut: Mesh, wc: WakeCut) -> None:
    """Roadmap P2 preprocess-time topology asserts, on the CUT mesh.

    1. each wake face has exactly one "+"-side and one "-"-side element
       (post-cut: every wake_plus/wake_minus face owned by exactly one tet,
       and the two owners were distinct pre-cut by construction);
    2. TE nodes ARE duplicated -- the sheet jump reaches the wall (this
       re-specs the original "TE nodes NOT duplicated" assert, which was
       shown to produce a divergent spurious TE force; module docstring
       and roadmap P2 evidence note);
    3. no orphan duplicated nodes (every slave referenced by some tet);
    4. duplicated-node count equals wake-sheet-node count minus the
       single-valued free-edge nodes (M1 tip edge; zero on quasi-2D
       sheets, so this is the original P2 assert there);
    5. wake-boundary edges handled: no element references BOTH a master and
       its own slave, and boundary faces reference the side-consistent copy
       (their owning tet references the same copy).
    """
    elements = np.asarray(mesh_cut.elements, dtype=np.int64)
    all_keys, all_owners = _all_tet_faces(elements)
    order = np.argsort(all_keys)
    all_keys_sorted = all_keys[order]
    all_owners_sorted = all_owners[order]

    def _owner_counts(tris):
        keys = _face_key_view(np.asarray(tris, dtype=np.int64))
        lo = np.searchsorted(all_keys_sorted, keys, side="left")
        hi = np.searchsorted(all_keys_sorted, keys, side="right")
        return hi - lo, lo

    # 1. one owner per side.
    for tag in ("wake_minus", "wake_plus"):
        counts, _ = _owner_counts(mesh_cut.boundary_faces[tag])
        assert np.all(counts == 1), (
            f"{tag}: {int(np.sum(counts != 1))} face(s) not owned by exactly "
            "one tet after the cut"
        )

    # 2. TE nodes ARE duplicated (jump reaches the wall).
    assert np.isin(wc.te_nodes, wc.master_nodes).all(), (
        "TE node(s) missing from the duplication map (the sheet jump must "
        "extend to the wall; see module docstring)"
    )

    # 3. no orphan slaves.
    referenced = np.zeros(len(mesh_cut.nodes), dtype=bool)
    referenced[elements.reshape(-1)] = True
    assert referenced[wc.slave_nodes].all(), (
        f"{int((~referenced[wc.slave_nodes]).sum())} duplicated node(s) "
        "referenced by no element (orphans)"
    )

    # 4. count equality: slaves == wake-sheet nodes (TE included) minus
    # the single-valued free-edge nodes.
    n_wake_nodes = len(np.unique(wc.wake_faces_minus))
    assert len(wc.slave_nodes) == n_wake_nodes - len(wc.free_nodes), (
        f"duplicated {len(wc.slave_nodes)} nodes but sheet has "
        f"{n_wake_nodes} nodes and {len(wc.free_nodes)} free-edge nodes"
    )

    # 5a. no element straddles the cut.
    slave_index = np.full(len(mesh_cut.nodes), -1, dtype=np.int64)
    slave_index[wc.slave_nodes] = np.arange(len(wc.slave_nodes))
    for e_nodes in (elements,):
        has_master = np.isin(e_nodes, wc.master_nodes)
        has_slave = slave_index[e_nodes] >= 0
        both = has_master.any(axis=1) & has_slave.any(axis=1)
        if np.any(both):
            for e in np.where(both)[0]:
                m_set = {int(x) for x in e_nodes[e][has_master[e]]}
                s_masters = {
                    int(wc.master_nodes[slave_index[x]])
                    for x in e_nodes[e][has_slave[e]]
                }
                assert not (m_set & s_masters), (
                    f"element {e} references a master and its own slave "
                    "(side classification leak across the sheet)"
                )

    # 5b. boundary faces agree with their owning tet's copy choice.
    for tag, tris in mesh_cut.boundary_faces.items():
        counts, lo = _owner_counts(tris)
        assert np.all(counts == 1), (
            f"boundary group '{tag}': {int(np.sum(counts != 1))} face(s) "
            "not owned by exactly one tet in the cut mesh"
        )
