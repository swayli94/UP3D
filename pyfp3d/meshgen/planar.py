"""
Vanilla-Gmsh 2D triangulations for the M0 quasi-2D cases (roadmap Track M).

Each builder meshes a planar domain with the Gmsh API (OCC/geo kernels, a
Distance+Threshold background size field graded from the wall -- same policy
as cases/meshes/sphere_shell/generate_sphere_shell.py) and returns plain
numpy arrays; the single-layer extrusion and prism->tet split happen in
extrude.py, never inside Gmsh. The wake is a Gmsh *embedded curve in
surface* (gmsh.model.mesh.embed) so triangle edges conform to it; node
duplication stays solver-side (agent-rules hard rule 8).

Gmsh is imported lazily so the solver test suite does not depend on it.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np


def _collect_2d(
    curve_groups: Dict[str, List[int]],
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]:
    """Extract (points2d, triangles, edge groups) from the current Gmsh model.

    Node tags are compacted to 0-based indices over the nodes actually
    referenced by triangles (Gmsh tags are 1-based and may have gaps).
    """
    import gmsh

    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
    index_of = np.full(int(node_tags.max()) + 1, -1, dtype=np.int64)

    tri_type = gmsh.model.mesh.getElementType("Triangle", 1)
    _, tri_nodes = gmsh.model.mesh.getElementsByType(tri_type)
    tri_tags_raw = np.asarray(tri_nodes, dtype=np.int64).reshape(-1, 3)

    used_tags = np.unique(tri_tags_raw)
    index_of[used_tags] = np.arange(len(used_tags))
    tag_to_row = {int(t): i for i, t in enumerate(node_tags)}
    points2d = coords[[tag_to_row[int(t)] for t in used_tags], :2]
    triangles = index_of[tri_tags_raw]

    line_type = gmsh.model.mesh.getElementType("Line", 1)
    edge_groups: Dict[str, np.ndarray] = {}
    for name, curves in curve_groups.items():
        edges = []
        for curve in curves:
            types, _, node_lists = gmsh.model.mesh.getElements(1, curve)
            for etype, nlist in zip(types, node_lists):
                if etype != line_type:
                    continue
                edges.append(np.asarray(nlist, dtype=np.int64).reshape(-1, 2))
        if not edges:
            raise RuntimeError(f"curve group '{name}' produced no line elements")
        stacked = np.vstack(edges)
        mapped = index_of[stacked]
        if np.any(mapped < 0):
            raise RuntimeError(
                f"curve group '{name}' references nodes not used by any triangle"
            )
        edge_groups[name] = mapped
    return points2d, triangles, edge_groups


def cylinder_annulus_2d(
    radius: float = 1.0,
    r_far: float = 20.0,
    h_wall: float = 0.1,
    h_far: float = 2.5,
    dist_min: float = 0.3,
    dist_max: float = 10.0,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]:
    """Circle-in-circle annular domain for the cylinder flow test case.

    Returns (points2d, triangles, edge_groups) with edge groups
    "wall" (inner circle) and "farfield" (outer circle).
    """
    import gmsh

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("cylinder_2d")

        outer = gmsh.model.occ.addDisk(0, 0, 0, r_far, r_far)
        inner = gmsh.model.occ.addDisk(0, 0, 0, radius, radius)
        gmsh.model.occ.cut([(2, outer)], [(2, inner)])
        gmsh.model.occ.synchronize()

        wall_curves, far_curves = [], []
        for dim, tag in gmsh.model.getEntities(1):
            bbox = gmsh.model.getBoundingBox(dim, tag)
            size = max(bbox[3] - bbox[0], bbox[4] - bbox[1])
            if abs(size / 2 - radius) < 0.25 * radius:
                wall_curves.append(tag)
            else:
                far_curves.append(tag)
        if not wall_curves or not far_curves:
            raise RuntimeError("could not identify wall/farfield curves")

        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", wall_curves)
        gmsh.model.mesh.field.setNumber(1, "Sampling", 200)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", h_wall)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", h_far)
        gmsh.model.mesh.field.setNumber(2, "DistMin", dist_min)
        gmsh.model.mesh.field.setNumber(2, "DistMax", dist_max)
        gmsh.model.mesh.field.setAsBackgroundMesh(2)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        gmsh.model.mesh.generate(2)
        return _collect_2d({"wall": wall_curves, "farfield": far_curves})
    finally:
        gmsh.finalize()


def naca0012_coordinates(n_half: int = 120) -> np.ndarray:
    """Closed-TE NACA0012 polyline, TE -> upper -> LE -> lower -> TE.

    Uses the closed-trailing-edge coefficient set (last coefficient -0.1036
    instead of -0.1015) so the TE is a single sharp point at (1, 0); the
    first and last rows are both exactly (1, 0) -- callers should drop one.
    Cosine clustering refines LE and TE.
    """
    beta = np.linspace(0.0, np.pi, n_half)
    x = 0.5 * (1.0 - np.cos(beta))  # 0 (LE) .. 1 (TE), cosine-clustered
    t = 0.12
    yt = 5.0 * t * (
        0.2969 * np.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x**2
        + 0.2843 * x**3
        - 0.1036 * x**4
    )
    upper = np.stack([x[::-1], yt[::-1]], axis=1)  # TE -> LE
    lower = np.stack([x[1:], -yt[1:]], axis=1)     # LE -> TE (skip repeated LE)
    coords = np.vstack([upper, lower])
    coords[0] = (1.0, 0.0)
    coords[-1] = (1.0, 0.0)
    return coords


def naca0012_wake_2d(
    r_far: float = 15.0,
    h_wall: float = 0.02,
    h_far: float = 3.0,
    h_wake: Optional[float] = None,
    dist_min: float = 0.1,
    dist_max: float = 8.0,
    wake_dist_max: float = 1.5,
    n_half: int = 120,
    embed_wake: bool = True,
    corridor_alpha_deg: Tuple[float, float] = (-6.0, 6.0),
    corridor_n_lines: int = 5,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """NACA0012 in a circular far field with an embedded wake line TE -> farfield.

    The far-field circle is centered at mid-chord (0.5, 0) and passes through
    the point (0.5 + r_far, 0) where the wake line ends, so the wake spans the
    whole distance from the sharp TE to the far-field boundary. The wake line
    is embedded in the surface (gmsh.model.mesh.embed): triangle edges conform
    to it, and its line elements are returned as the interior edge group
    "wake" for extrude.py to turn into the tagged interior face sheet.

    With ``embed_wake=False`` (the M3 wake-free family, roadmap Track M /
    design_track_b.md section 5.7) NO wake line is embedded -- the mesh
    topology knows nothing about the wake, which is Track B's level-set
    deliverable form. Because there is then no conforming sheet to attract
    refinement, a fan of ``corridor_n_lines`` size-field-only lines from the
    TE spanning ``corridor_alpha_deg`` (degrees, the intended alpha-sweep
    envelope; design_track_b.md D8) keeps element size ~h_wake in the wedge
    the level-set wake will sweep through. The fan lines are plain geometry
    used only in the Distance size field -- never embedded, so triangle
    edges do NOT conform to them. interior_edge_groups is then empty.

    Returns (points2d, triangles, edge_groups, interior_edge_groups):
    edge_groups has "wall" and "farfield"; interior_edge_groups has "wake"
    (embed_wake=True) or is {} (embed_wake=False).
    """
    import gmsh

    if h_wake is None:
        h_wake = 2.0 * h_wall

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("naca0012_2d")
        geo = gmsh.model.geo

        coords = naca0012_coordinates(n_half=n_half)[:-1]  # drop repeated TE
        pt = [geo.addPoint(x, y, 0.0) for x, y in coords]
        te = pt[0]
        i_le = int(np.argmin(coords[:, 0]))
        upper = geo.addSpline(pt[: i_le + 1])            # TE -> LE (upper)
        lower = geo.addSpline(pt[i_le:] + [te])          # LE -> TE (lower)
        foil_loop = geo.addCurveLoop([upper, lower])

        cx = 0.5
        center = geo.addPoint(cx, 0.0, 0.0)
        p_e = geo.addPoint(cx + r_far, 0.0, 0.0)
        p_n = geo.addPoint(cx, r_far, 0.0)
        p_w = geo.addPoint(cx - r_far, 0.0, 0.0)
        p_s = geo.addPoint(cx, -r_far, 0.0)
        arcs = [
            geo.addCircleArc(p_e, center, p_n),
            geo.addCircleArc(p_n, center, p_w),
            geo.addCircleArc(p_w, center, p_s),
            geo.addCircleArc(p_s, center, p_e),
        ]
        far_loop = geo.addCurveLoop(arcs)
        surface = geo.addPlaneSurface([far_loop, foil_loop])

        if embed_wake:
            wake_size_lines = [geo.addLine(te, p_e)]
            geo.synchronize()
            gmsh.model.mesh.embed(1, wake_size_lines, 2, surface)
        else:
            # M3 corridor fan: size-field-only lines from the TE covering the
            # alpha-sweep wedge. Kept strictly inside the far-field circle so
            # the standalone curves never touch the boundary; NOT embedded.
            a_lo, a_hi = np.radians(corridor_alpha_deg)
            fan_len = 0.97 * (r_far - 0.5)
            wake_size_lines = []
            for ang in np.linspace(a_lo, a_hi, corridor_n_lines):
                p_end = geo.addPoint(1.0 + fan_len * np.cos(ang),
                                     fan_len * np.sin(ang), 0.0)
                wake_size_lines.append(geo.addLine(te, p_end))
            geo.synchronize()

        # Graded size: h_wall near the airfoil, h_wake near the wake sheet,
        # h_far in the bulk; overall size = min of the two thresholds.
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", [upper, lower])
        gmsh.model.mesh.field.setNumber(1, "Sampling", 400)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", h_wall)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", h_far)
        gmsh.model.mesh.field.setNumber(2, "DistMin", dist_min)
        gmsh.model.mesh.field.setNumber(2, "DistMax", dist_max)

        gmsh.model.mesh.field.add("Distance", 3)
        gmsh.model.mesh.field.setNumbers(3, "CurvesList", wake_size_lines)
        gmsh.model.mesh.field.setNumber(3, "Sampling", 400)
        gmsh.model.mesh.field.add("Threshold", 4)
        gmsh.model.mesh.field.setNumber(4, "InField", 3)
        gmsh.model.mesh.field.setNumber(4, "SizeMin", h_wake)
        gmsh.model.mesh.field.setNumber(4, "SizeMax", h_far)
        gmsh.model.mesh.field.setNumber(4, "DistMin", dist_min)
        gmsh.model.mesh.field.setNumber(4, "DistMax", wake_dist_max)

        gmsh.model.mesh.field.add("Min", 5)
        gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2, 4])
        gmsh.model.mesh.field.setAsBackgroundMesh(5)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        gmsh.model.mesh.generate(2)
        curve_groups = {"wall": [upper, lower], "farfield": arcs}
        if embed_wake:
            curve_groups["wake"] = wake_size_lines
        points2d, triangles, groups = _collect_2d(curve_groups)
        interior = {"wake": groups.pop("wake")} if embed_wake else {}
        return points2d, triangles, groups, interior
    finally:
        gmsh.finalize()
