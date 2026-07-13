"""
M2 deliverable (part 2 of 2): ONERA M6 wing + simplified axisymmetric
fuselage, fused into a wing-body half model and meshed WAKE-FREE for the
level-set (Track B) solver path (roadmap Track M / M2; solver leg = B9).

Design choices (and why they differ from the plan sketch):

  * Wake-free only. M2's roadmap note says the wake-fuselage junction is the
    case a pre-embedded conforming sheet handles worst, and to schedule it
    with Track B's embedded (level-set) wake, which removes the need to embed
    a sheet at all. So this generator NEVER fragments/embeds a wake surface;
    the wake lives only as the source of a corridor size field, exactly like
    the M4 wake-free wing family (wing3d embed_wake=False). The solver builds
    its wake from the analytic TE polyline via WakeLevelSet (see te_polyline).

  * Walls split into "wall" (wing skin) and "fuselage" (fuselage skin), both
    natural no-penetration walls. The split is not needed to SOLVE (both are
    do-nothing boundaries) but it lets B9 post-process the fuselage lift
    separately (its gate: the axisymmetric body carries ~zero lift) and lets
    the caller keep fuselage nodes out of the wing TE control volumes at the
    junction. Surfaces are separated by testing whether a face lies on the
    fuselage's surface of revolution (y^2 + z^2 = R(x)^2).

  * Round tip cap by default (tip_cap="round"), so the only sharp wall
    features are the wing LE/TE that carry the physics -- the M5 geometry.

The wing solid, its round cap, the planform helpers and the mesh-extraction
and far-field machinery are reused verbatim from wing3d.py; this module only
adds the fuselage fuse and the wing/fuselage surface split.

Axis convention (fixed project-wide): chord +x, lift +y, span +z; symmetry
plane z = 0, wing tip z = B_SEMI. The wing root is buried in the fuselage; the
wing emerges at the junction near z = r_f, and the wing TE (y = 0) meets the
fuselage skin at z = r_f exactly -- the inboard end of the wake polyline.
"""

import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from pyfp3d.mesh.reader import Mesh
from pyfp3d.meshgen import wing3d
from pyfp3d.meshgen.wing3d import (
    B_SEMI,
    C_ROOT,
    MAC,
    TIP_CAP_RADIUS,
    TIP_CAPS,
    _add_round_tip_cap,
    _add_section_wire,
    _collect_3d,
    x_te,
)
from pyfp3d.meshgen.fuselage import FuselageParams, add_fuselage_solid, radius_at


def junction_z(p: FuselageParams) -> float:
    """Spanwise station where the wing TE (y = 0) meets the fuselage skin.

    On the chord plane the fuselage surface is at z = R(x) = r_f in the
    constant-radius section, and the M6 root chord (x in [0, 0.806]) lies
    wholly inside that section, so the TE-line junction is z = r_f exactly.
    """
    return p.r_f


def te_polyline(p: FuselageParams) -> np.ndarray:
    """Analytic TE polyline for WakeLevelSet: from the wing-fuselage junction
    (z = r_f) out to the tip TE corner (z = B_SEMI), in the chord plane y = 0.
    """
    z0 = junction_z(p)
    return np.array([
        [x_te(z0), 0.0, z0],
        [x_te(B_SEMI), 0.0, B_SEMI],
    ])


def _surface_on_fuselage(tag: int, p: FuselageParams, tol: float,
                         n: int = 6) -> bool:
    """True if 2D entity `tag` lies on the fuselage surface of revolution.

    Samples a parametric grid on the (synchronized) surface and tests the
    fraction of points satisfying sqrt(y^2 + z^2) == R(x). A wing face fails
    everywhere except right at the junction seam, so a 0.7 majority cleanly
    separates the two even though both are O(chord)-sized, non-symmetry,
    non-far surfaces (the bounding-box classifier cannot tell them apart).
    """
    import gmsh

    lo, hi = gmsh.model.getParametrizationBounds(2, tag)
    us = np.linspace(lo[0], hi[0], n)
    vs = np.linspace(lo[1], hi[1], n)
    uv: List[float] = []
    for u in us:
        for v in vs:
            uv.extend((float(u), float(v)))
    xyz = np.asarray(gmsh.model.getValue(2, tag, uv)).reshape(-1, 3)
    rr = np.hypot(xyz[:, 1], xyz[:, 2])
    R = np.array([radius_at(p, float(x)) for x in xyz[:, 0]])
    return float(np.mean(np.abs(rr - R) < tol)) > 0.7


def onera_m6_wingbody_mesh(
    h_wall: float,
    fuselage: Optional[FuselageParams] = None,
    h_far: Optional[float] = None,
    h_wake: Optional[float] = None,
    h_edge: Optional[float] = None,
    h_tip: Optional[float] = None,
    h_junction: Optional[float] = None,
    r_far: float = 15.0 * MAC,
    name: str = "onera_m6_wingbody",
    algorithm3d: int = 1,
    verbose: bool = False,
    tip_cap: str = "round",
    n_profile: int = 80,
) -> Mesh:
    """Generate the ONERA M6 wing-body half-model volume mesh (wake-free).

    One parameter (h_wall) sets the level; the rest default to scales of it,
    matching the wing3d / M5 roundtip policy (no h_far clamp, so a
    coarse/medium/fine ladder is self-similar).

    Args:
        h_wall: target element size on the wing + fuselage surface [m]
        fuselage: FuselageParams (default = the standard simplified body)
        h_far:  far-field size (default 120 h_wall -- NO clamp, self-similar)
        h_wake: wake-corridor size (default 3 h_wall)
        h_edge: wing LE/TE edge size (default 0.5 h_wall)
        h_tip:  round-cap size (default 0.25 h_wall; ignored if flat)
        h_junction: wing-fuselage junction-curve size (default 0.5 h_wall --
                the junction is a high-curvature intersection needing refinement)
        tip_cap: "round" (default, M5) or "flat"
        n_profile: fuselage meridian spline sample count

    Returns:
        Mesh with boundary groups "wall" (wing skin), "fuselage"
        (fuselage skin), "farfield", "symmetry". No "wake" group -- the wake
        is carried by the level set from te_polyline(fuselage).
    """
    import gmsh

    if tip_cap not in TIP_CAPS:
        raise ValueError(f"tip_cap must be one of {TIP_CAPS}, got {tip_cap!r}")
    p = fuselage if fuselage is not None else FuselageParams()
    if h_far is None:
        h_far = 120.0 * h_wall            # NO clamp (M1b defect); self-similar
    if h_wake is None:
        h_wake = 3.0 * h_wall
    if h_edge is None:
        h_edge = 0.5 * h_wall
    if h_tip is None:
        h_tip = 0.25 * h_wall
    if h_junction is None:
        h_junction = 0.5 * h_wall

    z_lo = -0.08 * C_ROOT
    z_junc = junction_z(p)
    # Far-field sphere centered on the whole body (nose..tail), not just the
    # wing -- the fuselage now sets the fore/aft extent.
    xc = 0.5 * (p.x_nose_tip + max(x_te(B_SEMI), p.x_tail_tip))
    x_down = xc + r_far + 0.5

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.model.add(name)
        occ = gmsh.model.occ

        # --- wing solid (loft + optional round cap), verbatim from wing3d ---
        root = _add_section_wire(occ, z_lo)
        tip = _add_section_wire(occ, B_SEMI)
        wing = occ.addThruSections(
            [root.loop, tip.loop], makeSolid=True, makeRuled=True
        )
        wing_vols = [dt for dt in wing if dt[0] == 3]
        if tip_cap == "round":
            wing_vols, _ = occ.fuse(wing_vols, _add_round_tip_cap(occ, tip))

        # --- fuse the fuselage on ------------------------------------------
        fus_vols = add_fuselage_solid(occ, p, n_profile=n_profile)
        body, _ = occ.fuse(wing_vols, fus_vols)

        # --- half-ball fluid domain, minus the wing-body -------------------
        ball = occ.addSphere(xc, 0.0, 0.0, r_far)
        below = occ.addBox(xc - 2 * r_far, -2 * r_far, -2 * r_far,
                           4 * r_far, 4 * r_far, 2 * r_far)
        half, _ = occ.cut([(3, ball)], [(3, below)])
        fluid, _ = occ.cut(half, body)
        occ.synchronize()

        vols = gmsh.model.getEntities(3)
        assert len(vols) == 1, f"expected one fluid volume, got {vols}"
        vol_tag = vols[0][1]

        # --- wake sheet (size-field source only, never embedded) -----------
        pa = occ.addPoint(x_te(z_junc), 0.0, z_junc)
        pb = occ.addPoint(x_te(B_SEMI), 0.0, B_SEMI)     # tip TE corner
        pc = occ.addPoint(x_down, 0.0, B_SEMI)
        pd = occ.addPoint(x_down, 0.0, z_junc)
        lines = [occ.addLine(pa, pb), occ.addLine(pb, pc),
                 occ.addLine(pc, pd), occ.addLine(pd, pa)]
        sheet = occ.addPlaneSurface([occ.addCurveLoop(lines)])
        occ.synchronize()

        # --- classify boundary surfaces ------------------------------------
        tol = 1e-3
        fus_tol = 5e-3          # spline vs piecewise-radius deviation margin
        groups: Dict[str, List[int]] = {"wall": [], "fuselage": [],
                                        "farfield": [], "symmetry": []}
        for dim, tag in gmsh.model.getEntities(2):
            if tag == sheet:
                continue                       # size-field-only, not a group
            bb = gmsh.model.getBoundingBox(dim, tag)
            extent = max(bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2])
            if (bb[5] - bb[2]) < tol and abs(bb[5]) < tol:
                groups["symmetry"].append(tag)
            elif extent > 1.2 * r_far:
                groups["farfield"].append(tag)
            elif _surface_on_fuselage(tag, p, fus_tol):
                groups["fuselage"].append(tag)
            else:
                groups["wall"].append(tag)
        for g, tags in groups.items():
            assert tags, f"surface classification found no '{g}' faces"

        # Round cap: the only wall geometry outboard of z = B_SEMI.
        tip_faces: List[int] = []
        if tip_cap == "round":
            z_cut = B_SEMI + 0.2 * TIP_CAP_RADIUS
            tip_faces = [t for t in groups["wall"]
                         if gmsh.model.getBoundingBox(2, t)[5] > z_cut]
            assert tip_faces, "rounded tip cap not found among the wall faces"

        # LE / TE edge curves for the edge-refinement field (straight,
        # thin-in-y, spanning z). After the fuse the exposed TE runs from the
        # junction to the tip; the buried inboard part is absorbed.
        edge_curves = _wing_edge_curves(gmsh, tol)

        # Wing-fuselage junction curves (intersection of wing & fuselage
        # skins) -- refined so the junction is well resolved.
        junction_curves = _junction_curves(gmsh, groups, p, fus_tol)

        # --- graded size fields (M0 policy: pure background field) ---------
        field = gmsh.model.mesh.field
        thresholds: List[int] = []

        def _dist_threshold(surfs=None, curves=None, size_min=h_wall,
                            dist_min=0.05, dist_max=0.55 * r_far, sampling=200):
            f = field.add("Distance")
            if surfs:
                field.setNumbers(f, "SurfacesList", surfs)
            if curves:
                field.setNumbers(f, "CurvesList", curves)
            field.setNumber(f, "Sampling", sampling)
            t = field.add("Threshold")
            field.setNumber(t, "InField", f)
            field.setNumber(t, "SizeMin", size_min)
            field.setNumber(t, "SizeMax", h_far)
            field.setNumber(t, "DistMin", dist_min)
            field.setNumber(t, "DistMax", dist_max)
            thresholds.append(t)
            return t

        _dist_threshold(surfs=groups["wall"] + groups["fuselage"],
                        size_min=h_wall)
        _dist_threshold(surfs=[sheet], size_min=h_wake,
                        dist_min=0.05, dist_max=1.2)
        if edge_curves:
            _dist_threshold(curves=edge_curves, size_min=h_edge,
                            dist_min=0.02, dist_max=0.3, sampling=400)
        if tip_faces:
            _dist_threshold(surfs=tip_faces, size_min=h_tip,
                            dist_min=0.02, dist_max=0.3)
        if junction_curves:
            _dist_threshold(curves=junction_curves, size_min=h_junction,
                            dist_min=0.02, dist_max=0.3, sampling=400)

        f_min = field.add("Min")
        field.setNumbers(f_min, "FieldsList", thresholds)
        field.setAsBackgroundMesh(f_min)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", algorithm3d)
        gmsh.option.setNumber("Mesh.Optimize", 1)

        gmsh.model.mesh.generate(3)

        mesh = _collect_3d(groups, name=name)
        _wingbody_asserts(mesh, p, r_far=r_far, xc=xc, tol=1e-7 * r_far,
                          tip_cap=tip_cap)
        return mesh
    finally:
        gmsh.finalize()


def _wing_edge_curves(gmsh, tol: float) -> List[int]:
    """Straight, thin-in-y curves spanning a good fraction of the span = the
    exposed wing LE and TE edges (the junction curve is thick in y)."""
    from pyfp3d.meshgen.wing3d import x_le
    found: List[int] = []
    for dim, tag in gmsh.model.getEntities(1):
        bb = gmsh.model.getBoundingBox(dim, tag)
        if (bb[4] - bb[1]) < tol and (bb[5] - bb[2]) > 0.05 * B_SEMI:
            # thin in y and spanning in z: an LE or TE edge (or the wake
            # sheet's own boundary lines -- but those were added after and
            # share the chord plane; screen by x being near x_le/x_te).
            zc = 0.5 * (bb[2] + bb[5])
            xmid = 0.5 * (bb[0] + bb[3])
            if abs(xmid - x_le(zc)) < 0.02 or abs(xmid - x_te(zc)) < 0.02:
                found.append(tag)
    return found


def _junction_curves(gmsh, groups: Dict[str, List[int]], p: FuselageParams,
                     tol: float) -> List[int]:
    """Curves shared between a wing ("wall") face and a fuselage face = the
    wing-fuselage junction. Identified as curves whose sampled points lie on
    the fuselage revolution surface AND that bound at least one wall face."""
    wall_set = set(groups["wall"])
    found: List[int] = []
    for dim, tag in gmsh.model.getEntities(1):
        # boundary faces of this curve
        up, _ = gmsh.model.getAdjacencies(1, tag)
        if not any(f in wall_set for f in up):
            continue
        if not any(f in set(groups["fuselage"]) for f in up):
            continue
        found.append(tag)
    return found


def _wingbody_asserts(mesh: Mesh, p: FuselageParams, r_far: float, xc: float,
                      tol: float, tip_cap: str) -> None:
    """Generation-time invariants for the wing-body mesh."""
    from pyfp3d.mesh.metrics import compute_tet_volumes

    vols = compute_tet_volumes(mesh.nodes, mesh.elements)
    assert vols.min() > 0.0, "non-positive tet volume"

    assert "wake" not in mesh.boundary_faces, \
        "wing-body mesh must be wake-free (level-set path)"

    sym_nodes = np.unique(mesh.boundary_faces["symmetry"])
    assert np.abs(mesh.nodes[sym_nodes, 2]).max() < tol, \
        "symmetry plane not at z = 0"

    far_nodes = np.unique(mesh.boundary_faces["farfield"])
    r = np.linalg.norm(mesh.nodes[far_nodes] - np.array([xc, 0.0, 0.0]), axis=1)
    assert np.all(np.abs(r - r_far) < 1e-3 * r_far), \
        "far-field nodes not on the sphere"

    # Fuselage skin genuinely lies on the revolution surface.
    fus_nodes = np.unique(mesh.boundary_faces["fuselage"])
    fx = mesh.nodes[fus_nodes, 0]
    frr = np.hypot(mesh.nodes[fus_nodes, 1], mesh.nodes[fus_nodes, 2])
    fR = np.array([radius_at(p, float(x)) for x in fx])
    # spline vs piecewise deviation ~ a few mm; allow 1% of r_f + 3 mm
    assert np.abs(frr - fR).max() < 0.01 * p.r_f + 3e-3, \
        "fuselage group contains faces off the revolution surface"

    # Tip TE corner present (wake attaches there), unmoved by the fuse.
    wall_nodes = np.unique(mesh.boundary_faces["wall"])
    corner = np.array([x_te(B_SEMI), 0.0, B_SEMI])
    d_corner = np.linalg.norm(mesh.nodes[wall_nodes] - corner, axis=1).min()
    assert d_corner < 1e-9, \
        f"no wall node at the tip TE corner {corner} (closest {d_corner:.2e})"

    # Wake polyline inboard end: the wing TE (y = 0) meets the fuselage skin
    # (z = r_f in the cylinder section) at an EXACT OCC vertex, so this is a
    # real wall node, not merely a nearby one (measured 1.5e-9 at coarse and
    # medium alike). The level-set wake's inboard end therefore lands on the
    # wall exactly, the way the tip TE corner does at the other end.
    z_junc = junction_z(p)
    te_junc = np.array([x_te(z_junc), 0.0, z_junc])
    d_junc = np.linalg.norm(mesh.nodes[wall_nodes] - te_junc, axis=1).min()
    assert d_junc < 1e-6, \
        (f"no wall node at the junction TE {te_junc} (closest {d_junc:.2e}); "
         "the wake polyline inboard end is unresolved")

    # Tip cap reach.
    z_wall_max = float(mesh.nodes[wall_nodes, 2].max())
    if tip_cap == "flat":
        assert z_wall_max <= B_SEMI + 1e-9, \
            f"flat cap reaches z = {z_wall_max}, past the tip {B_SEMI}"
    else:
        apex = B_SEMI + TIP_CAP_RADIUS
        assert B_SEMI < z_wall_max <= apex + 1e-9, \
            f"rounded cap reaches z = {z_wall_max}, expected ({B_SEMI}, {apex}]"
