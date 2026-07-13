"""
M1 deliverable: swept/tapered 3D half-wing tet mesh with an embedded wake
sheet (roadmap Track M, gate M1; feeds solver phase P5).

Case: ONERA M6 half wing. Geometry per the NASA GRC NPARC validation
archive (https://www.grc.nasa.gov/www/wind/valid/m6wing/m6wing.html):
ONERA D symmetric section, foilmod.txt coordinate set (the TE thickness is
"rounded off" to exactly zero, so the wake sheet attaches to a sharp TE
edge), straight-tapered planform.

Tip cap (``tip_cap``, Track M gates M1 and M5):

  "flat"  (default, M1)  the loft is closed by the planar tip section at
          z = B_SEMI. This is the standard FP-validation simplification --
          and it is WRONG in a way that matters: the flat cap meets the
          upper/lower surfaces at a SHARP CONVEX EDGE, which in potential
          flow is an edge singularity. P13/G13.3 measured it diverging
          under uniform refinement (tip-cap box peak Mach exponent
          p = +0.321) while the wake edge and the wing interior stayed
          bounded -- i.e. it, and not the wake, is what still blocks 3D
          grid convergence. Kept as the default because the P5 / P8-G8.2 /
          B7 / M1 regression locks are all anchored to it.

  "round" (M5)  the cap is the half body of revolution swept by rotating
          the tip section about its own chord line: the wing is closed by
          {sqrt(y^2 + (z - B_SEMI)^2) <= t(x)}, where t(x) is the tip
          section's local half thickness. There is no edge: the cap meets
          the wing surface tangentially (dy/dz = 0 on both sides at
          z = B_SEMI), and because t -> 0 at the LE and the TE the cap
          degenerates to a point exactly at each -- so the TE line, the
          wake sheet, the Kutta stations, B_SEMI and every solver-side
          semantic are UNCHANGED. Only the wall geometry near the tip
          moves, by at most TIP_CAP_RADIUS = 22 mm (1.9% of the semi-span)
          and only where the section has thickness.

Solver axis convention (fixed by mesh/wake_cut.py and the M0 cases):
chord along +x, lift along +y (wake_cut upper_hint), span along +z.
Root symmetry plane at z = 0, tip at z = B_SEMI. NASA planform units
(meters) are kept; the FP equations are scale-invariant.

Construction (all vanilla Gmsh/OCC -- node duplication stays solver-side,
agent-rules hard rule 8):
  1. Ruled loft between two section wires (root extended slightly below
     z = 0, tip). For a straight-tapered wing with a single section shape
     the two-section ruled loft is the EXACT planform surface, and its TE
     edge is exactly the straight segment between the section TE points.
  1b. tip_cap="round" only: the upper half of the tip section face is
     revolved a full turn about the tip chord line (an edge OF that face,
     so the revolution is the standard degenerate-at-the-axis kind, like a
     sphere from a half disc) and FUSED onto the loft. The revolved solid
     lies strictly inside the wing for z < B_SEMI -- verified analytically
     and asserted at generation: it protrudes by at most 4e-6 m, three
     orders below the finest edge size, and never reaches aft of the local
     TE, so it can neither cut the wake sheet nor spawn slivers. What it
     does is replace the flat cap's sharp edge with a tangent surface.
  2. Spherical far field of radius r_far (~15 reference chords) centered
     on the symmetry plane, cut to the z >= 0 half domain, minus the wing.
  3. Planar wake sheet in the chord plane y = 0, swept from the TE line
     (built from the same endpoints as the wing TE edge, so OCC fragment
     stitches them into one shared edge) downstream past the sphere and
     below z = 0; occ.fragment embeds it conformally in the fluid volume
     and trims it to: TE edge (wall), root edge (symmetry plane),
     downstream boundary (far-field sphere), tip edge (z = B_SEMI,
     interior FREE edge from the tip TE corner downstream).
  4. Boundary surfaces are classified geometrically and tagged
     wall / farfield / symmetry, the interior sheet is tagged wake.

Wake-tip closure (M1 gate): the sheet's tip edge starts exactly at the
wing tip TE corner (same OCC vertex after fragment) -- no crack, no
overlap. The tip edge itself lies in the domain interior; handling of its
non-splitting node stars is the solver preprocessor's job
(mesh/wake_cut.py free-edge support; design.md Sec 12 risk 1).

Gmsh is imported lazily so the solver test suite does not depend on it.
"""

import math
from typing import Dict, List, NamedTuple, Optional

import numpy as np

from pyfp3d.mesh.reader import Mesh

# ---------------------------------------------------------------------------
# ONERA D section, foilmod.txt (zero TE thickness). Columns: x/c, t/c.
# The section is symmetric: lower surface is t -> -t.
# ---------------------------------------------------------------------------
ONERA_D_UPPER = np.array([
    (0.0000000, 0.0000000), (0.0000165, 0.0006914), (0.0000696, 0.0014416),
    (0.0001675, 0.0022554), (0.0003232, 0.0031382), (0.0005508, 0.0040959),
    (0.0008657, 0.0051343), (0.0012868, 0.0062598), (0.0018364, 0.0074784),
    (0.0025441, 0.0087958), (0.0034428, 0.0102163), (0.0045704, 0.0117419),
    (0.0059751, 0.0133708), (0.0077112, 0.0150951), (0.0098413, 0.0168984),
    (0.0124479, 0.0187537), (0.0156171, 0.0206220), (0.0194609, 0.0224545),
    (0.0241067, 0.0242004), (0.0297008, 0.0258245), (0.0364261, 0.0273317),
    (0.0444852, 0.0287912), (0.0541248, 0.0303278), (0.0656303, 0.0320138),
    (0.0793366, 0.0338372), (0.0956354, 0.0357742), (0.1149796, 0.0377923),
    (0.1378963, 0.0398522), (0.1649976, 0.0419089), (0.1919327, 0.0436214),
    (0.2187096, 0.0450507), (0.2453310, 0.0462358), (0.2717978, 0.0471987),
    (0.2981113, 0.0479494), (0.3242726, 0.0484902), (0.3502830, 0.0488183),
    (0.3761446, 0.0489296), (0.4018567, 0.0488202), (0.4274223, 0.0484833),
    (0.4528441, 0.0479351), (0.4781197, 0.0471661), (0.5032514, 0.0461903),
    (0.5282426, 0.0450209), (0.5530937, 0.0436741), (0.5778043, 0.0421684),
    (0.6023757, 0.0405241), (0.6268104, 0.0387613), (0.6511093, 0.0368990),
    (0.6752726, 0.0349542), (0.6993027, 0.0329402), (0.7231995, 0.0308662),
    (0.7469658, 0.0287365), (0.7705998, 0.0265505), (0.7941055, 0.0243027),
    (0.8174828, 0.0219842), (0.8407324, 0.0195838), (0.8638564, 0.0170915),
    (0.8868235, 0.0145051), (0.9061905, 0.0121952), (0.9225336, 0.0101138),
    (0.9363346, 0.0083265), (0.9479946, 0.0068038), (0.9578511, 0.0055144),
    (0.9661860, 0.0044240), (0.9732361, 0.0035015), (0.9792020, 0.0027211),
    (0.9842508, 0.0020606), (0.9885252, 0.0015014), (0.9921438, 0.0010280),
    (0.9952080, 0.0006271), (0.9978030, 0.0002876), (1.0000000, 0.0000000),
])

# M6 planform (NASA GRC table 2), meters.
B_SEMI = 1.1963          # semi-span
TAPER = 0.562            # tip/root chord ratio
MAC = 0.64607            # mean aerodynamic chord (reference length)
LE_SWEEP_DEG = 30.0
C_ROOT = MAC / (2.0 / 3.0 * (1 + TAPER + TAPER**2) / (1 + TAPER))  # 0.8059
C_TIP = TAPER * C_ROOT

_TAN_SWEEP = math.tan(math.radians(LE_SWEEP_DEG))

#: Largest radius of the M5 rounded tip cap = the tip section's max half
#: thickness (the cap radius at chord station x is the section half thickness
#: there, so it tapers to zero at the LE and the TE). The rounded wing reaches
#: z = B_SEMI + TIP_CAP_RADIUS; the flat one stops at z = B_SEMI.
TIP_CAP_RADIUS = float(ONERA_D_UPPER[:, 1].max()) * C_TIP   # 0.02217 m

TIP_CAPS = ("flat", "round")


def chord_at(z: float) -> float:
    """Local chord of the straight-tapered planform (linear in z)."""
    return C_ROOT * (1.0 - (1.0 - TAPER) * (z / B_SEMI))


def x_le(z: float) -> float:
    """Leading-edge x of the swept planform."""
    return z * _TAN_SWEEP


def x_te(z: float) -> float:
    """Trailing-edge x (the wake sheet leading edge uses the SAME function
    as the wing sections so the OCC vertices coincide bitwise)."""
    return x_le(z) + chord_at(z)


class _Section(NamedTuple):
    """Tags of one section wire. `spl_up`, `p_le`, `p_te` are what the M5
    rounded cap revolves; the flat path uses only `loop`."""
    loop: int
    spl_up: int
    p_le: int
    p_te: int


def _add_section_wire(occ, z: float) -> _Section:
    """Closed two-spline section wire at span station z.

    Upper spline LE -> TE, lower spline TE -> LE; both share the exact LE
    and TE point tags so the loft produces distinct LE/TE BRep edges (the
    TE edge is what the wake sheet must stitch to).
    """
    c = chord_at(z)
    xle = x_le(z)
    upper_pts = [
        occ.addPoint(xle + xn * c, tn * c, z) for xn, tn in ONERA_D_UPPER
    ]
    # Lower surface: same x stations, negated thickness, TE -> LE,
    # reusing the shared TE and LE point tags at the ends.
    lower_inner = [
        occ.addPoint(xle + xn * c, -tn * c, z)
        for xn, tn in ONERA_D_UPPER[-2:0:-1]
    ]
    spl_up = occ.addSpline(upper_pts)
    spl_lo = occ.addSpline([upper_pts[-1]] + lower_inner + [upper_pts[0]])
    return _Section(occ.addCurveLoop([spl_up, spl_lo]), spl_up,
                    upper_pts[0], upper_pts[-1])


def _add_round_tip_cap(occ, tip: _Section) -> List[int]:
    """The M5 cap: revolve the upper half of the tip section face a full turn
    about the tip chord line, giving the solid of revolution
    {sqrt(y^2 + (z - B_SEMI)^2) <= t(x)}.

    The half face is bounded by the upper spline and the chord line, and the
    chord line IS the revolution axis -- so this is the degenerate-at-the-axis
    revolution OCC is built for (a sphere from a half disc), and it pinches to
    a point exactly at the LE and TE, where the section has zero thickness.
    """
    chord_line = occ.addLine(tip.p_te, tip.p_le)
    half_face = occ.addPlaneSurface([occ.addCurveLoop([tip.spl_up, chord_line])])
    cap = occ.revolve([(2, half_face)], 0.0, 0.0, B_SEMI, 1.0, 0.0, 0.0,
                      2.0 * math.pi)
    # The revolve leaves its base face behind as a free-standing surface. It
    # bounds no volume, but the geometric classification below would still see
    # it (chord-sized, thick in y, off the symmetry plane) and tag it "wall" --
    # whose triangles then reference nodes no tet uses. Drop it here.
    occ.remove([(2, half_face)], recursive=False)
    return [dt for dt in cap if dt[0] == 3]


def onera_m6_wing_mesh(
    h_wall: float,
    h_far: Optional[float] = None,
    h_wake: Optional[float] = None,
    h_edge: Optional[float] = None,
    r_far: float = 15.0 * MAC,
    name: str = "onera_m6",
    algorithm3d: int = 1,
    verbose: bool = False,
    embed_wake: bool = True,
    tip_cap: str = "flat",
    h_tip: Optional[float] = None,
) -> Mesh:
    """Generate the ONERA M6 half-wing volume mesh with embedded wake sheet.

    One parameter (h_wall) controls the level; everything else defaults to
    scales of it (same policy as the M0 family).

    With ``embed_wake=False`` (the M4 wake-free family, roadmap Track M /
    design_track_b.md section 5.7) the chord-plane sheet is still built,
    but ONLY as the source of the wake-corridor Distance size field: it is
    neither fragmented into the fluid BRep nor embedded, so the tet mesh
    does not conform to it and the mesh topology knows nothing about the
    wake -- Track B's deliverable form, now in 3D (swept TE, tip). The
    corridor sizing therefore matches the M1 family's (same h_wake, same
    Distance thresholds), which keeps the B7 A/B against the P5/P8
    baselines a controlled comparison. The returned Mesh then has NO
    "wake" group.

    Args:
        h_wall: target element size on the wing surface [m]
        h_far: size at the far field (default min(2.5, 120 h_wall))
        h_wake: size on the wake sheet (default 3 h_wall)
        h_edge: size on the LE/TE edges (default 0.5 h_wall; the TE value
                also controls the Kutta-probe neighborhood resolution)
        tip_cap: "flat" (M1, default -- planar cap, sharp convex edge) or
                "round" (M5 -- half body of revolution, no edge)
        h_tip: size on the rounded cap (default 0.25 h_wall, so it scales
                with the level and a refinement ladder stays self-similar).
                The cap radius is only 22 mm, so at h_wall the coarse cap
                would be about one element wide and would discretize back
                into a flat one -- the whole point of the geometry change
                would be lost. Ignored when tip_cap="flat".
        r_far: far-field sphere radius [m] (~15 MAC per roadmap M1)
        name: mesh name
        algorithm3d: gmsh Mesh.Algorithm3D (1 = Delaunay, robust with
                embedded surfaces)
        verbose: gmsh terminal output
        embed_wake: embed the wake sheet in the volume mesh (M1, default)
                or use it as a size field only (M4 wake-free)

    Returns:
        Mesh with boundary groups "wall", "farfield", "symmetry" and (when
        embed_wake) the interior sheet "wake"; tags follow the M0 naming so
        read_mesh / cut_wake ingest it unchanged.
    """
    import gmsh

    if tip_cap not in TIP_CAPS:
        raise ValueError(f"tip_cap must be one of {TIP_CAPS}, got {tip_cap!r}")
    if h_far is None:
        h_far = min(2.5, 120.0 * h_wall)
    if h_wake is None:
        h_wake = 3.0 * h_wall
    if h_edge is None:
        h_edge = 0.5 * h_wall
    if h_tip is None:
        h_tip = 0.25 * h_wall

    # Sphere center: mid planform on the symmetry plane.
    xc = 0.5 * (x_le(0.0) + x_te(B_SEMI))
    z_lo = -0.08 * C_ROOT       # root-side overhang, trimmed by the boolean
    x_down = xc + r_far + 0.5   # sheet overhang past the sphere

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.model.add(name)
        occ = gmsh.model.occ

        # --- wing solid: exact ruled loft root(extended) -> tip ----------
        root = _add_section_wire(occ, z_lo)
        tip = _add_section_wire(occ, B_SEMI)
        wing = occ.addThruSections(
            [root.loop, tip.loop], makeSolid=True, makeRuled=True
        )
        wing_vols = [dt for dt in wing if dt[0] == 3]

        if tip_cap == "round":
            # Fuse on the body of revolution; the flat cap face is strictly
            # interior to it and is absorbed, leaving a tangent surface.
            wing_vols, _ = occ.fuse(wing_vols, _add_round_tip_cap(occ, tip))

        # --- half-ball fluid domain --------------------------------------
        ball = occ.addSphere(xc, 0.0, 0.0, r_far)
        below = occ.addBox(xc - 2 * r_far, -2 * r_far, -2 * r_far,
                           4 * r_far, 4 * r_far, 2 * r_far)
        half, _ = occ.cut([(3, ball)], [(3, below)])
        fluid, _ = occ.cut(half, wing_vols)

        # --- wake sheet in the chord plane y = 0 -------------------------
        # Leading edge A->B reuses x_te() so it is the same line as the
        # wing TE edge (fragment stitches them); C, D overhang and are
        # trimmed against sphere / symmetry plane by the fragment.
        pa = occ.addPoint(x_te(z_lo), 0.0, z_lo)
        pb = occ.addPoint(x_te(B_SEMI), 0.0, B_SEMI)   # tip TE corner
        pc = occ.addPoint(x_down, 0.0, B_SEMI)
        pd = occ.addPoint(x_down, 0.0, z_lo)
        lines = [occ.addLine(pa, pb), occ.addLine(pb, pc),
                 occ.addLine(pc, pd), occ.addLine(pd, pa)]
        sheet = occ.addPlaneSurface([occ.addCurveLoop(lines)])

        tol = 1e-3
        if embed_wake:
            occ.fragment(fluid, [(2, sheet)])
            occ.synchronize()

            # --- drop sheet pieces trimmed OUTSIDE the fluid --------------
            # OCC bounding boxes are padded by ~1e-5, and the thinnest real
            # feature (a wing face's y-extent) is O(0.04 c), so 1e-3 m
            # cleanly separates "thin/at a plane" from "real extent".
            for dim, tag in gmsh.model.getEntities(2):
                bb = gmsh.model.getBoundingBox(dim, tag)
                thin_y = (bb[4] - bb[1]) < tol
                outside = bb[2] < -tol or bb[3] > xc + r_far + tol
                if thin_y and outside:
                    occ.remove([(dim, tag)], recursive=True)
            occ.synchronize()
            corridor_tags: List[int] = []
        else:
            # M4: the sheet stays a free-standing surface -- no fragment,
            # no embed. It is used ONLY by the wake Distance field below,
            # so the tet mesh never conforms to it (its own 2D mesh is
            # discarded by _collect_3d, which keeps only tet-referenced
            # nodes).
            occ.synchronize()
            corridor_tags = [sheet]

        vols = gmsh.model.getEntities(3)
        assert len(vols) == 1, f"expected one fluid volume, got {vols}"
        vol_tag = vols[0][1]

        # --- geometric classification of surfaces -------------------------
        groups: Dict[str, List[int]] = {"wall": [], "farfield": [],
                                        "symmetry": []}
        if embed_wake:
            groups["wake"] = []
        corridor_set = set(corridor_tags)
        for dim, tag in gmsh.model.getEntities(2):
            if tag in corridor_set:
                continue          # size-field-only sheet: not a mesh group
            bb = gmsh.model.getBoundingBox(dim, tag)
            extent = max(bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2])
            if (bb[5] - bb[2]) < tol and abs(bb[5]) < tol:
                groups["symmetry"].append(tag)
            elif embed_wake and (bb[4] - bb[1]) < tol:
                groups["wake"].append(tag)
            elif extent > 1.2 * r_far:
                # Only the far-field sphere patches span O(2 r_far); every
                # wing face is O(chord). (The symmetry disk also spans 2
                # r_far but is caught by the thin-z test above.)
                groups["farfield"].append(tag)
            else:
                groups["wall"].append(tag)
        for g, tags in groups.items():
            assert tags, f"surface classification found no '{g}' faces"

        # The rounded cap is the ONLY wall geometry outboard of z = B_SEMI --
        # the loft's own faces stop exactly there -- so a z_max test isolates
        # it. The threshold is a fraction of the cap radius: far above OCC's
        # ~1e-5 bounding-box padding, far below the cap's own 22 mm reach.
        tip_faces: List[int] = []
        if tip_cap == "round":
            z_cut = B_SEMI + 0.2 * TIP_CAP_RADIUS
            tip_faces = [t for t in groups["wall"]
                         if gmsh.model.getBoundingBox(2, t)[5] > z_cut]
            assert tip_faces, "rounded tip cap not found among the wall faces"

        if embed_wake:
            # The fragment stitched the sheet's boundary curves into the
            # fluid BRep (shared TE edge, symmetry/sphere trims) but the
            # sheet stays a free-standing surface; embedding makes the tet
            # mesh conform to it (wake faces shared by exactly two tets --
            # asserted by cut_wake at ingestion, and by _collect_3d node
            # conformity here).
            gmsh.model.mesh.embed(2, groups["wake"], 3, vol_tag)

        # LE / TE edges for the edge-refinement field, identified by their
        # exact segment endpoints (the planform edges are straight lines).
        def _curves_on_segment(p0, p1) -> List[int]:
            found = []
            for dim, tag in gmsh.model.getEntities(1):
                bb = gmsh.model.getBoundingBox(dim, tag)
                lo = np.minimum(p0, p1) - 1e-4
                hi = np.maximum(p0, p1) + 1e-4
                inside = (bb[0] >= lo[0] and bb[1] >= lo[1] and bb[2] >= lo[2]
                          and bb[3] <= hi[0] and bb[4] <= hi[1] and bb[5] <= hi[2])
                if not inside:
                    continue
                # thin in y and actually spanning in z
                if (bb[4] - bb[1]) < tol and (bb[5] - bb[2]) > 0.05 * B_SEMI:
                    found.append(tag)
            return found

        le_pts = (np.array([x_le(0.0), 0.0, 0.0]),
                  np.array([x_le(B_SEMI), 0.0, B_SEMI]))
        te_pts = (np.array([x_te(0.0), 0.0, 0.0]),
                  np.array([x_te(B_SEMI), 0.0, B_SEMI]))
        edge_curves = _curves_on_segment(*le_pts) + _curves_on_segment(*te_pts)

        # --- graded size fields (M0 policy: pure background field) --------
        field = gmsh.model.mesh.field
        f_wall = field.add("Distance")
        field.setNumbers(f_wall, "SurfacesList", groups["wall"])
        field.setNumber(f_wall, "Sampling", 200)
        t_wall = field.add("Threshold")
        field.setNumber(t_wall, "InField", f_wall)
        field.setNumber(t_wall, "SizeMin", h_wall)
        field.setNumber(t_wall, "SizeMax", h_far)
        field.setNumber(t_wall, "DistMin", 0.05)
        field.setNumber(t_wall, "DistMax", 0.55 * r_far)

        f_wake = field.add("Distance")
        field.setNumbers(
            f_wake, "SurfacesList",
            groups["wake"] if embed_wake else corridor_tags,
        )
        field.setNumber(f_wake, "Sampling", 200)
        t_wake = field.add("Threshold")
        field.setNumber(t_wake, "InField", f_wake)
        field.setNumber(t_wake, "SizeMin", h_wake)
        field.setNumber(t_wake, "SizeMax", h_far)
        field.setNumber(t_wake, "DistMin", 0.05)
        field.setNumber(t_wake, "DistMax", 1.2)

        thresholds = [t_wall, t_wake]
        if edge_curves:
            f_edge = field.add("Distance")
            field.setNumbers(f_edge, "CurvesList", edge_curves)
            field.setNumber(f_edge, "Sampling", 400)
            t_edge = field.add("Threshold")
            field.setNumber(t_edge, "InField", f_edge)
            field.setNumber(t_edge, "SizeMin", h_edge)
            field.setNumber(t_edge, "SizeMax", h_far)
            field.setNumber(t_edge, "DistMin", 0.02)
            field.setNumber(t_edge, "DistMax", 0.3)
            thresholds.append(t_edge)

        if tip_faces:
            f_tip = field.add("Distance")
            field.setNumbers(f_tip, "SurfacesList", tip_faces)
            field.setNumber(f_tip, "Sampling", 200)
            t_tip = field.add("Threshold")
            field.setNumber(t_tip, "InField", f_tip)
            field.setNumber(t_tip, "SizeMin", h_tip)
            field.setNumber(t_tip, "SizeMax", h_far)
            field.setNumber(t_tip, "DistMin", 0.02)
            field.setNumber(t_tip, "DistMax", 0.3)
            thresholds.append(t_tip)

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
        _builder_asserts(mesh, r_far=r_far, xc=xc, tol=1e-7 * r_far,
                         tip_cap=tip_cap)
        return mesh
    finally:
        gmsh.finalize()


def _collect_3d(surface_groups: Dict[str, List[int]], name: str) -> Mesh:
    """Extract (nodes, tets, named triangle groups) from the current Gmsh
    model into a solver Mesh. Node tags are compacted to 0-based indices
    over the nodes referenced by tets; every group triangle must reference
    only such nodes (i.e. the surface mesh conforms to the volume mesh)."""
    import gmsh

    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
    tag_to_row = np.full(int(node_tags.max()) + 1, -1, dtype=np.int64)
    tag_to_row[np.asarray(node_tags, dtype=np.int64)] = np.arange(len(node_tags))

    tet_type = gmsh.model.mesh.getElementType("Tetrahedron", 1)
    _, tet_nodes = gmsh.model.mesh.getElementsByType(tet_type)
    tets_raw = np.asarray(tet_nodes, dtype=np.int64).reshape(-1, 4)
    if len(tets_raw) == 0:
        raise RuntimeError("no tetrahedra generated")

    used_tags = np.unique(tets_raw)
    index_of = np.full(int(used_tags.max()) + 1, -1, dtype=np.int64)
    index_of[used_tags] = np.arange(len(used_tags))

    tri_type = gmsh.model.mesh.getElementType("Triangle", 1)
    boundary_faces: Dict[str, np.ndarray] = {}
    for gname, surf_tags in surface_groups.items():
        tris = []
        for s in surf_tags:
            types, _, node_lists = gmsh.model.mesh.getElements(2, s)
            for etype, nlist in zip(types, node_lists):
                if etype != tri_type:
                    continue
                tris.append(np.asarray(nlist, dtype=np.int64).reshape(-1, 3))
        if not tris:
            raise RuntimeError(f"surface group '{gname}' has no triangles")
        stacked = np.vstack(tris)
        if int(stacked.max()) >= len(index_of) or np.any(index_of[stacked] < 0):
            raise RuntimeError(
                f"group '{gname}' references nodes not used by any tet "
                "(surface mesh does not conform to the volume mesh)"
            )
        boundary_faces[gname] = np.ascontiguousarray(
            index_of[stacked], dtype=np.int32
        )

    mesh = Mesh()
    mesh.nodes = np.ascontiguousarray(coords[tag_to_row[used_tags]])
    mesh.elements = np.ascontiguousarray(index_of[tets_raw], dtype=np.int32)
    mesh.boundary_faces = boundary_faces
    mesh.element_tags = np.zeros(len(mesh.elements), dtype=np.int32)
    mesh.tag_names = ["fluid"]
    mesh.name = name
    mesh.validate()
    return mesh


def _builder_asserts(mesh: Mesh, r_far: float, xc: float, tol: float,
                     tip_cap: str = "flat") -> None:
    """Cheap generation-time invariants (full topology asserts run in the
    solver preprocessor, mesh/wake_cut.py)."""
    from pyfp3d.mesh.metrics import compute_tet_volumes

    vols = compute_tet_volumes(mesh.nodes, mesh.elements)
    assert vols.min() > 0.0, "non-positive tet volume"

    has_wake = "wake" in mesh.boundary_faces   # False for the M4 family
    if has_wake:
        wake_nodes = np.unique(mesh.boundary_faces["wake"])
        assert np.abs(mesh.nodes[wake_nodes, 1]).max() < tol, \
            "wake sheet not in the chord plane y = 0"

    sym_nodes = np.unique(mesh.boundary_faces["symmetry"])
    assert np.abs(mesh.nodes[sym_nodes, 2]).max() < tol, \
        "symmetry plane not at z = 0"

    far_nodes = np.unique(mesh.boundary_faces["farfield"])
    r = np.linalg.norm(mesh.nodes[far_nodes] - np.array([xc, 0.0, 0.0]),
                       axis=1)
    assert np.all(np.abs(r - r_far) < 1e-3 * r_far), \
        "far-field nodes not on the sphere"

    wall_nodes = np.unique(mesh.boundary_faces["wall"])
    if has_wake:
        te = np.intersect1d(wake_nodes, wall_nodes)
        assert len(te) >= 2, "wake sheet does not attach to the wing TE"

    # --- tip cap ----------------------------------------------------------
    # The wing must reach exactly z = B_SEMI when the cap is flat, and bulge
    # to B_SEMI + TIP_CAP_RADIUS when it is round -- with the tip TE corner
    # (where the cap radius vanishes and the wake attaches) unmoved either way.
    z_wall_max = float(mesh.nodes[wall_nodes, 2].max())
    if tip_cap == "flat":
        assert z_wall_max <= B_SEMI + 1e-9, (
            f"flat cap reaches z = {z_wall_max}, past the tip {B_SEMI}"
        )
    else:
        apex = B_SEMI + TIP_CAP_RADIUS
        assert B_SEMI < z_wall_max <= apex + 1e-9, (
            f"rounded cap reaches z = {z_wall_max}, expected in "
            f"({B_SEMI}, {apex}]"
        )
        # ... and it is actually resolved, not clipped to a couple of facets.
        assert z_wall_max > B_SEMI + 0.8 * TIP_CAP_RADIUS, (
            f"rounded cap only reaches z = {z_wall_max}: the mesh does not "
            f"resolve it (apex {apex}); lower h_tip"
        )
    corner = np.array([x_te(B_SEMI), 0.0, B_SEMI])
    d_corner = np.linalg.norm(mesh.nodes[wall_nodes] - corner, axis=1).min()
    assert d_corner < 1e-9, (
        f"no wall node at the tip TE corner {corner} (closest {d_corner:.2e}) "
        "-- the wake sheet attaches there"
    )
