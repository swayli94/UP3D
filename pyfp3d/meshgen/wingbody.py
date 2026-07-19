"""
M2 deliverable (part 2 of 2): ONERA M6 wing + simplified axisymmetric
fuselage, fused into a wing-body half model and meshed WAKE-FREE for the
level-set (Track B) solver path (roadmap Track M / M2; solver leg = B9).

Design choices (and why they differ from the plan sketch):

  * Wake-free by DEFAULT. M2's roadmap note says the wake-fuselage junction
    is the case a pre-embedded conforming sheet handles worst, and to schedule
    it with Track B's embedded (level-set) wake, which removes the need to
    embed a sheet at all. So by default this generator never fragments/embeds
    a wake surface; the wake lives only as the source of a corridor size
    field, exactly like the M4 wake-free wing family (wing3d
    embed_wake=False). The solver builds its wake from the analytic TE
    polyline via WakeLevelSet (see te_polyline).

  * embed_wake=True (B9 re-spec 2026-07-17, user-approved) additionally
    delivers the CONFORMING variant for the same geometry: the sheet is
    extended inboard to below the symmetry plane (z_lo, like wing3d) and
    downstream through the far-field sphere, then fragment+embed'ed. The
    fragment trims it automatically to: exposed wing TE (z in [r_f, B_SEMI])
    -> fuselage waterline (the y = 0 top meridian, z = R(x), junction TE ->
    tail tip) -> symmetry edge (z = 0 aft of the body) -> sphere arc ->
    free tip edge. Aft of the body the sheet MUST reach the symmetry plane
    (else the domain around body+wake is not simply connected and the branch
    cut fails), and it MUST reach the sphere (constraints/dirichlet.py
    prescribes branch-cut-multivalued far-field data; a sheet stopping short
    reproduces the P5 branch-ray artifact). cut_wake needs NO changes: the
    waterline edge lies in `fuselage` boundary triangles, so its nodes
    duplicate under the same boundary-edge rule as the wing-alone symmetry
    root edge, while TE/Kutta stations stay wing-only (wall∩wake).

  * Walls split into "wall" (wing skin) and "fuselage" (fuselage skin), both
    natural no-penetration walls. The split is not needed to SOLVE (both are
    do-nothing boundaries) but it lets B9 post-process the fuselage lift
    separately (its gate: the axisymmetric body carries ~zero lift) and lets
    the caller keep fuselage nodes out of the wing TE control volumes at the
    junction. Surfaces are separated by testing whether a face lies on the
    fuselage's surface of revolution (y^2 + z^2 = R(x)^2).

  * Round tip cap by default (tip_cap="round"), so the only sharp wall
    features are the wing LE/TE that carry the physics -- the M5 geometry.

  * The fuselage skin is sized SEPARATELY from the wing skin (h_body, default
    2 h_wall) -- user re-spec 2026-07-16. The body is now 5 root chords long
    with the wing in the middle (fuselage.py), so h_wall over the whole skin
    would spend most of the mesh resolving a body whose flow is a smooth
    displacement field. Near the wing the body still lands at h_wall, because
    the wing's own distance field wins the Min there; only the far body (nose,
    mid-body, afterbody) coarsens. The two tips are the exception and get
    their own refinement -- they are the only high-curvature parts of the body,
    and they are also the parts farthest from the wing.

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
from pyfp3d.meshgen.fuselage import (
    FuselageParams,
    add_fuselage_solid,
    add_fuselage_solid_split,
    radius_at,
)

#: Far-field sphere radius, 2026-07-16 (user-directed; was 15 MAC = 9.69, the
#: wing-ALONE convention inherited from wing3d). Once the body became 5 root
#: chords long the wing's MAC stopped being the scale that sets this: the old
#: sphere left the boundary only ~1.9 body lengths off the nose/tail tips, and
#: the body spanned ~42% of the far-field diameter. 25 MAC = 16.15 puts it
#: 3.5 body lengths clear of each tip and 25 MAC clear of the wing.
R_FAR = 25.0 * MAC

#: Far-field element size, in units of h_wall. Scaled WITH R_FAR (120 -> 200 =
#: 120 * 25/15) so that h_far / r_far is UNCHANGED. That ratio is what every
#: size field here is really written in: each Distance+Threshold grows from its
#: SizeMin to h_far across 0.55 * r_far, so holding h_far/r_far fixed holds
#: every near-body gradient — and therefore the whole mesh near the aircraft —
#: exactly where it was, while the domain grows. It is also why the bigger
#: domain is nearly free: the added shell's cells grow with it.
H_FAR_IN_H_WALL = 200.0

#: h_far of the design every FIXED refinement distance below was tuned against
#: (the 15-MAC / 120-h_wall family, through 2026-07-15). Those distances are
#: rescaled by h_far / (H_FAR_REF * h_wall) so that each field's GRADIENT --
#: SizeMin climbing to h_far across the distance -- is invariant when h_far
#: moves. Do not drop this: raising h_far alone (1.8 -> 3.0 at medium) steepens
#: every one of them by the same factor, and the wake corridor cannot take it.
#: Measured 2026-07-16 with the distances left at their old values: the corridor
#: under the sheet's inboard edge (z = 0.15, hanging over open fluid aft of the
#: body) demanded 0.3 m elements while the symmetry plane 1.2 m away sat at
#: 2.9 m, and the mesher answered with pancakes -- 72 tets under 5 deg, min
#: dihedral 0.17 deg vs the 2.0 deg bound. The 0.55 * r_far distances scale on
#: their own and were never part of the problem.
H_FAR_REF = 120.0


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
    h_body: Optional[float] = None,
    h_body_tip: Optional[float] = None,
    r_far: float = R_FAR,
    name: str = "onera_m6_wingbody",
    algorithm3d: int = 1,
    verbose: bool = False,
    tip_cap: str = "round",
    n_profile: int = 120,
    embed_wake: bool = False,
    junction_fillet: Optional[float] = None,
) -> Mesh:
    """Generate the ONERA M6 wing-body half-model volume mesh.

    Wake-free by default (the level-set / Track B form). embed_wake=True
    produces the CONFORMING variant (B9 re-spec 2026-07-17): the wake sheet
    is fragmented against the fluid and embedded, giving a "wake" boundary
    group whose inboard boundary is the fuselage waterline + the symmetry
    plane aft of the body, ready for mesh/wake_cut.py::cut_wake unchanged.
    The default path is bit-identical to the pre-embed_wake generator.

    One parameter (h_wall) sets the level; the rest default to scales of it,
    matching the wing3d / M5 roundtip policy (no h_far clamp, so a
    coarse/medium/fine ladder is self-similar).

    Args:
        h_wall: target element size on the WING surface [m]
        fuselage: FuselageParams (default = the standard simplified body)
        h_far:  far-field size (default 200 h_wall -- NO clamp, self-similar;
                scaled with r_far so h_far/r_far, and hence every gradient
                below, is fixed)
        h_wake: wake-corridor size (default 3 h_wall)
        h_edge: wing LE/TE edge size (default 0.5 h_wall)
        h_tip:  round-cap size (default 0.25 h_wall; ignored if flat)
        h_junction: wing-fuselage junction-curve size (default 0.5 h_wall --
                the junction is a high-curvature intersection needing refinement)
        h_body: fuselage skin size AWAY from the wing (default 2 h_wall). Near
                the wing the body still gets h_wall from the wing's own field.
        h_body_tip: fuselage nose/tail tip size (default 0.25 h_wall, the same
                policy as the M5 round tip cap) -- the tips are the body's only
                high-curvature regions, and the size ramps back up to h_body
                over the nose / afterbody length
        tip_cap: "round" (default, M5) or "flat"
        n_profile: fuselage meridian spline sample count
        junction_fillet: None (default, bit-identical) or a fairing RADIUS:
                the wing-fuselage reentrant corner is filled by a sphere-pipe
                fairing (ersatz tangent blend, _fillet_junction -- OCC's own
                fillet fails on the LE/TE cusps). A Ball-field strip at the
                sampled junction line keeps the blend resolved at h_junction.
                Wake-free path only (embed_wake is rejected). D3 leg of
                cases/analysis/b23_junction_discriminator.

    Returns:
        Mesh with boundary groups "wall" (wing skin), "fuselage"
        (fuselage skin), "farfield", "symmetry". Default: no "wake" group --
        the wake is carried by the level set from te_polyline(fuselage).
        embed_wake=True: plus the embedded interior "wake" group.
    """
    import gmsh

    if tip_cap not in TIP_CAPS:
        raise ValueError(f"tip_cap must be one of {TIP_CAPS}, got {tip_cap!r}")
    if junction_fillet is not None:
        if junction_fillet <= 0.0:
            raise ValueError("junction_fillet must be > 0 (a radius)")
        if embed_wake:
            raise ValueError("junction_fillet is supported on the wake-free "
                             "path only (D3 discriminator leg)")
    p = fuselage if fuselage is not None else FuselageParams()
    if h_far is None:
        # NO clamp (M1b defect); self-similar. 200, not 120: scaled with R_FAR.
        h_far = H_FAR_IN_H_WALL * h_wall
    if h_wake is None:
        h_wake = 3.0 * h_wall
    if h_edge is None:
        h_edge = 0.5 * h_wall
    if h_tip is None:
        h_tip = 0.25 * h_wall
    if h_junction is None:
        h_junction = 0.5 * h_wall
    if h_body is None:
        h_body = 2.0 * h_wall             # far body: coarser than the wing
    if h_body_tip is None:
        h_body_tip = 0.25 * h_wall        # as the M5 round tip cap (h_tip)

    z_lo = -0.08 * C_ROOT
    z_junc = junction_z(p)
    # Far-field sphere centered on the whole body (nose..tail), not just the
    # wing -- the fuselage now sets the fore/aft extent.
    xc = 0.5 * (p.x_nose_tip + max(x_te(B_SEMI), p.x_tail_tip))
    # Downstream end of the wake sheet. WAKE-FREE default: the sheet is a
    # SIZE-FIELD SOURCE only, so its extent is free -- and it must stop SHORT
    # of the far-field sphere, which `xc + r_far + 0.5` did not.
    # Measured 2026-07-16: a sheet that pierces the sphere drags its corridor
    # (h_wake = 3 h_wall within ~1.2 of the sheet) onto the OUTFLOW BOUNDARY,
    # which put 0.087 m triangles on a sphere sized h_far = 6.0 and produced
    # every bad tet in the mesh -- all 13 sat at the downstream pole, min
    # dihedral 3.03 deg vs 5.60 with the sheet clear. Nothing is served by
    # resolving a wake as it leaves a Dirichlet boundary. 0.8 r_far clears the
    # corridor by ~2.4 m and still leaves MORE downstream corridor (~16 MAC)
    # than the old 15-MAC sphere physically had room for.
    # EMBEDDED variant: the sheet MUST pierce the sphere (wing3d convention;
    # the fragment trims it to the sphere) -- the far-field Dirichlet data is
    # branch-cut-multivalued and needs duplicated wake nodes ON the sphere; a
    # sheet stopping short is the P5 branch-ray artifact. The downstream-pole
    # sliver cost recorded above is accepted and gated (Netgen is on).
    x_down = xc + (r_far + 0.5 if embed_wake else 0.8 * r_far)

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
        # embed_wake needs the skin split at the y = 0 meridians so the wake
        # sheet's waterline is a real seam edge (add_fuselage_solid_split);
        # the wake-free path keeps the un-split builder (bit-identical).
        fus_vols = (add_fuselage_solid_split(occ, p, n_profile=n_profile)
                    if embed_wake else
                    add_fuselage_solid(occ, p, n_profile=n_profile))
        body, _ = occ.fuse(wing_vols, fus_vols)

        # --- optional junction fillet (D3 discriminator leg) ----------------
        # Fillet the wing-fuselage intersection edges of the fused BODY solid,
        # before the fluid cut. junction_pts = the sampled pre-fillet junction
        # line, reused below as Ball-field centers so the blend stays resolved
        # at h_junction (the _junction_curves field finds no edge afterwards).
        junction_pts: Optional[np.ndarray] = None
        if junction_fillet is not None:
            body, junction_pts = _fillet_junction(gmsh, occ, body, p,
                                                  radius=junction_fillet)
            occ.synchronize()

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

        # --- wake sheet -----------------------------------------------------
        tol = 1e-3
        if embed_wake:
            # Conforming variant (B9). A rectangle in the chord plane y = 0,
            # spanning z = [z_lo, tip] (z_lo below the symmetry plane, like
            # wing3d) and x = [x_te(z), x_down] (LE collinear with the wing
            # TE line, which is straight, so it stitches the exposed TE),
            # PASSING THROUGH the body. occ.fragment trims it against the
            # fluid to the exact in-fluid cut: exposed wing TE (junction ->
            # tip) + the fuselage TOP WATERLINE (a real seam edge thanks to
            # add_fuselage_solid_split) + the z = 0 symmetry ray aft of the
            # body + the tip/downstream edges. The LE is pre-split at the
            # junction TE point so a sheet vertex lands exactly on the fused
            # body's junction-TE OCC vertex (d_junc < 1e-6).
            pa = occ.addPoint(x_te(z_lo), 0.0, z_lo)
            pj = occ.addPoint(x_te(z_junc), 0.0, z_junc)  # junction TE
            pb = occ.addPoint(x_te(B_SEMI), 0.0, B_SEMI)  # tip TE corner
            pc = occ.addPoint(x_down, 0.0, B_SEMI)
            pd = occ.addPoint(x_down, 0.0, z_lo)
            lines = [occ.addLine(pa, pj), occ.addLine(pj, pb),
                     occ.addLine(pb, pc), occ.addLine(pc, pd),
                     occ.addLine(pd, pa)]
            sheet = occ.addPlaneSurface([occ.addCurveLoop(lines)])
            occ.fragment(fluid, [(2, sheet)])
            occ.synchronize()

            # Drop sheet pieces trimmed OUTSIDE the fluid: below the symmetry
            # plane, beyond the sphere, or inside the body. OCC bboxes are
            # padded ~1e-5 and the thinnest real feature (a wing face's
            # y-extent) is O(0.04 c), so 1e-3 separates a sheet piece from
            # a real body/far surface.
            for dim, tag in gmsh.model.getEntities(2):
                bb = gmsh.model.getBoundingBox(dim, tag)
                if (bb[4] - bb[1]) >= tol:
                    continue                     # not a sheet piece
                outside = (bb[2] < -tol or bb[3] > xc + r_far + tol
                           or _sheet_piece_in_body(gmsh, tag, p))
                if outside:
                    occ.remove([(dim, tag)], recursive=True)
            occ.synchronize()

            # Re-query the volume: fragment/remove may renumber.
            vols = gmsh.model.getEntities(3)
            assert len(vols) == 1, f"expected one fluid volume, got {vols}"
            vol_tag = vols[0][1]
            corridor_tags: List[int] = []
        else:
            # Size-field source only, never embedded (the wake-free default).
            pa = occ.addPoint(x_te(z_junc), 0.0, z_junc)
            pb = occ.addPoint(x_te(B_SEMI), 0.0, B_SEMI)     # tip TE corner
            pc = occ.addPoint(x_down, 0.0, B_SEMI)
            pd = occ.addPoint(x_down, 0.0, z_junc)
            lines = [occ.addLine(pa, pb), occ.addLine(pb, pc),
                     occ.addLine(pc, pd), occ.addLine(pd, pa)]
            sheet = occ.addPlaneSurface([occ.addCurveLoop(lines)])
            occ.synchronize()
            corridor_tags = [sheet]

        # --- classify boundary surfaces ------------------------------------
        fus_tol = 5e-3          # spline vs piecewise-radius deviation margin
        groups: Dict[str, List[int]] = {"wall": [], "fuselage": [],
                                        "farfield": [], "symmetry": []}
        if embed_wake:
            groups["wake"] = []
        corridor_set = set(corridor_tags)
        for dim, tag in gmsh.model.getEntities(2):
            if tag in corridor_set:
                continue                       # size-field-only, not a group
            bb = gmsh.model.getBoundingBox(dim, tag)
            extent = max(bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2])
            if (bb[5] - bb[2]) < tol and abs(bb[5]) < tol:
                groups["symmetry"].append(tag)
            elif embed_wake and (bb[4] - bb[1]) < tol:
                # Thin-in-y = a surviving sheet piece (the wing skins are
                # O(0.04 c) thick in y; the split fuselage quarter-shells are
                # O(r_f); symmetry pieces were caught by thin-z above).
                groups["wake"].append(tag)
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

        if embed_wake:
            # The fragment stitched the sheet's boundary curves into the
            # fluid BRep (shared TE edge, waterline, symmetry/sphere trims)
            # but the sheet stays free-standing; embedding makes the tet mesh
            # conform to it (wake faces shared by exactly two tets --
            # asserted by cut_wake at ingestion). M1 trap: fragment alone
            # does NOT conform a sheet with an interior free edge (the tip).
            gmsh.model.mesh.embed(2, groups["wake"], 3, vol_tag)

        # --- graded size fields (M0 policy: pure background field) ---------
        field = gmsh.model.mesh.field
        thresholds: List[int] = []

        #: Scales every FIXED refinement distance below with h_far, holding each
        #: field's gradient at the H_FAR_REF design's. See H_FAR_REF.
        grad = h_far / (H_FAR_REF * h_wall)

        def _dist_threshold(surfs=None, curves=None, size_min=h_wall,
                            dist_min=0.05, dist_max=0.55 * r_far, sampling=200,
                            register=True):
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
            if register:
                thresholds.append(t)
            return t

        def _tip_ball(x_tip: float, ramp: float) -> int:
            """h_body_tip on the axis at a body tip, ramping to h_body over
            `ramp`. VOut is h_body, NOT h_far, because this field is read as a
            size FLOOR under Max, never under Min (see _fuselage_field).
            """
            f = field.add("Ball")
            field.setNumber(f, "XCenter", x_tip)
            field.setNumber(f, "YCenter", 0.0)
            field.setNumber(f, "ZCenter", 0.0)
            field.setNumber(f, "Radius", 0.2 * p.r_f)
            field.setNumber(f, "Thickness", ramp)
            field.setNumber(f, "VIn", h_body_tip)
            field.setNumber(f, "VOut", h_body)
            return f

        def _fuselage_field() -> int:
            """Body skin at h_body, EXCEPT near the two tips, where the size
            follows the body's local RADIUS (2026-07-16).

            Why the tips need their own law: the body is a surface of
            revolution, so its circumferential faceting angle is h / R(x), and
            R(x) collapses to zero at both ends -- to r_tail = 30 mm on the
            tail cap, and through the nose ellipsoid's 37 mm tip curvature
            radius. A single h_body over the whole skin therefore facets the
            afterbody into a polygonal needle however generous h_body looks
            against r_f (measured on the first cut of this re-spec: cone p99
            45.9 deg, 20 of the 26 worst edges on the aft cone, vs 20.7 deg on
            the cylinder). Ramping h_body_tip -> h_body over the afterbody
            length reproduces h ~ 0.4 R(x) along the cone almost exactly,
            because a cone's radius grows linearly with distance from its tip;
            over the nose it is conservative (R grows like sqrt there).

            Composed with Max, not Min: each of these two fields states a
            FLOOR on the size ("no finer than"), and the result is the finest
            size honouring both -- the tip balls floor the skin at h_body
            everywhere except the tips, and the distance ramp floors the size
            at h_far once away from the body. Min would instead read the tip
            balls' VOut = h_body as a CEILING and cap the whole domain, far
            field included, at 2 h_wall. The caller's Min over `thresholds`
            then applies the real ceilings (the wing, wake and edge fields).

            The witness that the far field escaped: n_tris_farfield in the
            stats CSV is 193 coarse / 732 medium, vs 193 / 730 on the pre-
            re-spec body -- untouched. Under Min it would have exploded.
            """
            ramp = _dist_threshold(surfs=groups["fuselage"],
                                   size_min=h_body_tip, register=False)
            tips = field.add("Min")
            field.setNumbers(tips, "FieldsList",
                             [_tip_ball(p.x_nose_tip, 0.5 * p.l_nose),
                              _tip_ball(p.x_tail_tip, p.l_tail + p.r_tail)])
            f = field.add("Max")
            field.setNumbers(f, "FieldsList", [ramp, tips])
            thresholds.append(f)
            return f

        # Wing skin and fuselage skin are sized separately (2026-07-16): the
        # body is 5 root chords long and only the part near the wing needs
        # h_wall. The wing's field wins the Min over the body near the
        # junction, so the transition needs no explicit blend region.
        _dist_threshold(surfs=groups["wall"], size_min=h_wall)
        _fuselage_field()
        _dist_threshold(surfs=(groups["wake"] if embed_wake else corridor_tags),
                        size_min=h_wake, dist_min=0.05, dist_max=1.2 * grad)
        if edge_curves:
            _dist_threshold(curves=edge_curves, size_min=h_edge,
                            dist_min=0.02, dist_max=0.3 * grad, sampling=400)
        if tip_faces:
            _dist_threshold(surfs=tip_faces, size_min=h_tip,
                            dist_min=0.02, dist_max=0.3 * grad)
        if junction_curves:
            _dist_threshold(curves=junction_curves, size_min=h_junction,
                            dist_min=0.02, dist_max=0.3 * grad, sampling=400)
        if junction_pts is not None and len(junction_pts):
            # D3 fillet leg: the crease is gone, so the junction-curve field
            # above found nothing. Keep the BLEND resolved with a chain of
            # Ball fields along the sampled pre-fillet junction line: VIn =
            # h_junction on the blend, ramping to h_wall over ~2 radii.
            balls: List[int] = []
            for pt in junction_pts:
                fb = field.add("Ball")
                field.setNumber(fb, "XCenter", float(pt[0]))
                field.setNumber(fb, "YCenter", float(pt[1]))
                field.setNumber(fb, "ZCenter", float(pt[2]))
                field.setNumber(fb, "Radius", 1.5 * junction_fillet)
                field.setNumber(fb, "Thickness", 2.0 * junction_fillet)
                field.setNumber(fb, "VIn", h_junction)
                field.setNumber(fb, "VOut", h_wall)
                balls.append(fb)
            f_balls = field.add("Min")
            field.setNumbers(f_balls, "FieldsList", balls)
            thresholds.append(f_balls)

        f_min = field.add("Min")
        field.setNumbers(f_min, "FieldsList", thresholds)
        field.setAsBackgroundMesh(f_min)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", algorithm3d)
        gmsh.option.setNumber("Mesh.Optimize", 1)
        # Netgen's tet optimizer on top of gmsh's own (2026-07-16). This family
        # has ONE fragile spot, and it is structural rather than a tuning
        # accident: the wake sheet's inboard edge sits at z = z_junc = 0.15,
        # hanging over open fluid everywhere aft of the body, so its corridor
        # prints a ~0.2 m-wide FINE RIBBON down the symmetry plane at y ~ 0
        # while the plane a metre away is at h_far. The mesher answers that
        # with pancakes. Measured on medium at r_far = 25 MAC, min dihedral by
        # h_far/h_wall: 120 -> 0.31 deg, 160 -> 4.80, 200 -> 2.63 -- NOT
        # monotone, i.e. a lottery, and the pre-2026-07-16 family was simply
        # winning it (5.88). Netgen removes the ribbon's slivers outright: the
        # same three cases give 0 tets under 5 deg and min dihedral 11-13.6,
        # BETTER than the 15-MAC family ever was and on ~9% FEWER tets.
        # It only touches the volume -- the fuselage crease p99 (10.07) and
        # every boundary invariant in _wingbody_asserts are unchanged, which is
        # what makes it a mesh-quality knob and not a change of geometry.
        # NOT set in wing3d.py: the M1/M4/M5 families stay bit-identical.
        #
        # embed_wake path: Netgen is OFF (the documented B9 fallback). Two
        # reasons. (1) The pathology Netgen fixed on the wake-FREE family --
        # the wake corridor's fine ribbon hanging over open fluid at
        # z = z_junc, slivering the symmetry plane -- is STRUCTURALLY REMOVED
        # by embedding: the sheet's inboard edge is now a real mesh boundary
        # riding the fuselage waterline and the aft symmetry plane, not a
        # size-field ghost over open fluid. (2) Netgen's re-tetting SEGFAULTS
        # on this embedded-sheet + split-body geometry (measured 2026-07-17,
        # SIGSEGV in generate(3) at the fine level); gmsh's own optimizer
        # (Mesh.Optimize=1) is kept. Quality is gated by the generator's
        # QUALITY_BOUNDS + the cut_wake ingest gate.
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 0 if embed_wake else 1)

        gmsh.model.mesh.generate(3)

        mesh = _collect_3d(groups, name=name)
        _wingbody_asserts(mesh, p, r_far=r_far, xc=xc, tol=1e-7 * r_far,
                          tip_cap=tip_cap, embed_wake=embed_wake,
                          filleted=junction_fillet is not None)
        return mesh
    finally:
        gmsh.finalize()


def _fillet_junction(gmsh, occ, body, p: FuselageParams,
                     radius: float, n_sample: int = 30,
                     fus_tol: float = 5e-3):
    """Fair the wing-fuselage reentrant corner with a straight CONE fairing
    pair (D3 discriminator leg).

    Why not occ.fillet / sphere-pipe / addPipe -- all measured failing or
    intractable 2026-07-19: OCC's blend dies on the LE/TE cusps ("Could not
    compute fillet", r = 0.02..0.05); a 60-sphere chain needs > 15 min of
    batched BOPs; addPipe loses the curve to OCC renumbering ("Unknown
    OpenCASCADE wire"). What works cheaply: TWO straight cone trains along
    the junction line. On the constant-radius section the junction is
    nearly straight, so a tube of radius `radius` centered at
    (y = +/-radius, z = r_f + radius) is tangent to BOTH the wing plane
    (y = 0) and the fuselage cylinder (sqrt(y^2+z^2) = r_f) -- an ersatz
    fairing, not a true fillet (the wing skin has thickness and the
    junction y wanders by O(0.04)), which is all the "kill the crease" A/B
    needs. Each side is taper-cone / cylinder / taper-cone so the ends die
    INSIDE the body. 6 tools, fused in batches of 3 (a ~60-tool single BOP
    silently merges nothing -- measured).

    STATUS 2026-07-19: geometry-side fuse succeeds, but the smoke mesh
    (canonical D4 case, h_wall = 0.04, rho = 0.045) put gmsh into a
    pathological state (> 15 min, RSS runaway past 4 GB where the unfilleted
    h04 mesh is 46k tets in seconds) -- the fused fairing leaves sliver
    surfaces that the surface mesher cannot close. Recorded LOW-COST
    INFEASIBLE per the D3 pre-registered fallback; the parameter stays as
    the documented dead end (default None keeps the default path
    bit-identical).

    Returns (new_body_dimtags, (n,3) junction-line points) -- the points
    feed the Ball-field strip that keeps the blend resolved (the
    crease-edge refinement field has no edge to hang on once the corner is
    filled).
    """
    occ.synchronize()
    vols3 = gmsh.model.getEntities(3)
    assert len(vols3) == 1, f"expected one body volume, got {vols3}"
    vol_tag = vols3[0][1]
    fus_faces = set()
    for dim, tag in gmsh.model.getEntities(2):
        if _surface_on_fuselage(tag, p, fus_tol):
            fus_faces.add(tag)
    jcurves: List[int] = []
    for dim, tag in gmsh.model.getEntities(1):
        up, _ = gmsh.model.getAdjacencies(1, tag)
        up = list(up)
        if not up:
            continue
        has_fus = any(f in fus_faces for f in up)
        has_other = any(f not in fus_faces for f in up)
        if has_fus and has_other:
            jcurves.append(tag)
    if not jcurves:
        raise RuntimeError("junction_fillet: no wing-fuselage junction edges "
                           "found on the fused body")
    pts: List[np.ndarray] = []
    for t in jcurves:
        lo, hi = gmsh.model.getParametrizationBounds(1, t)
        us = np.linspace(lo[0], hi[0], n_sample)
        uv: List[float] = []
        for u in us:
            uv.append(float(u))
        xyz = np.asarray(gmsh.model.getValue(1, t, uv)).reshape(-1, 3)
        pts.append(xyz)
    sampled = np.vstack(pts)

    # straight cone-train fairing along the junction chord extent
    x0 = float(sampled[:, 0].min()) - 0.02     # just fore of the LE cusp
    x1 = float(sampled[:, 0].max()) + 0.02     # just aft of the TE cusp
    taper = 0.12
    rho = radius
    z_c = p.r_f + rho
    tools: List[Tuple[int, int]] = []
    for sgn in (+1.0, -1.0):
        y_c = sgn * rho
        tools.append((3, occ.addCone(x0 - taper, y_c, z_c, taper, 0.0, 0.0,
                                     0.004, rho)))
        tools.append((3, occ.addCylinder(x0, y_c, z_c, x1 - x0, 0.0, 0.0,
                                         rho)))
        tools.append((3, occ.addCone(x1, y_c, z_c, taper, 0.0, 0.0,
                                     rho, 0.004)))
    vol = (3, vol_tag)
    for i in range(0, len(tools), 3):
        out, _ = occ.fuse([vol], tools[i:i + 3])
        vols_out = [dt for dt in out if dt[0] == 3]
        if len(vols_out) != 1:
            raise RuntimeError(
                f"junction_fillet: fairing batch {i // 3} produced "
                f"{len(vols_out)} volumes (expected 1)")
        vol = vols_out[0]
    return [vol], sampled


def _sheet_piece_in_body(gmsh, tag: int, p: FuselageParams,
                         n: int = 6) -> bool:
    """True if 2D entity `tag` (a thin-in-y sheet piece) sits INSIDE the
    fuselage: majority of a parametric sample grid has |z| < R(x). Kept fluid
    pieces lie at z > R(x) (or beyond the body's x-range, where R = 0)."""
    lo, hi = gmsh.model.getParametrizationBounds(2, tag)
    us = np.linspace(lo[0], hi[0], n)
    vs = np.linspace(lo[1], hi[1], n)
    uv: List[float] = []
    for u in us:
        for v in vs:
            uv.extend((float(u), float(v)))
    xyz = np.asarray(gmsh.model.getValue(2, tag, uv)).reshape(-1, 3)
    R = np.array([radius_at(p, float(x)) for x in xyz[:, 0]])
    inside = np.abs(xyz[:, 2]) < R - 1e-4
    return float(np.mean(inside)) > 0.5


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
                      tol: float, tip_cap: str,
                      embed_wake: bool = False,
                      filleted: bool = False) -> None:
    """Generation-time invariants for the wing-body mesh. `filleted` = the D3
    junction-fillet variant: the exact junction-TE vertex check is dropped
    (the blend legitimately moves the TE/skin meeting point off the analytic
    crease point by O(fillet radius)); every other invariant still applies."""
    from pyfp3d.mesh.metrics import compute_tet_volumes

    vols = compute_tet_volumes(mesh.nodes, mesh.elements)
    assert vols.min() > 0.0, "non-positive tet volume"

    if not embed_wake:
        assert "wake" not in mesh.boundary_faces, \
            "wing-body mesh must be wake-free (level-set path)"
    else:
        _wingbody_wake_asserts(mesh, p, r_far=r_far, xc=xc)

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
    # D3 fillet variant: the blend moves the TE/skin meeting point off the
    # analytic crease point by O(fillet radius) -- checked as a NEARBY node
    # instead of an exact vertex.
    z_junc = junction_z(p)
    te_junc = np.array([x_te(z_junc), 0.0, z_junc])
    d_junc = np.linalg.norm(mesh.nodes[wall_nodes] - te_junc, axis=1).min()
    assert d_junc < (0.05 if filleted else 1e-6), \
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


def _wingbody_wake_asserts(mesh: Mesh, p: FuselageParams, r_far: float,
                           xc: float) -> None:
    """Generation-time invariants specific to the EMBEDDED (conforming)
    variant. Fires before cut_wake ever sees the mesh, so a bad fragment is
    caught at the generator, with geometry still in scope."""
    assert "wake" in mesh.boundary_faces, "embed_wake=True but no wake group"
    wake_nodes = np.unique(mesh.boundary_faces["wake"])
    wn = mesh.nodes[wake_nodes]

    # Planar sheet in the chord plane.
    assert np.abs(wn[:, 1]).max() < 1e-7, "wake nodes off the y = 0 plane"

    # Both TE endpoints are shared wall∩wake nodes (exact vertices).
    wall_nodes = set(np.unique(mesh.boundary_faces["wall"]).tolist())
    shared = [n for n in wake_nodes.tolist() if n in wall_nodes]
    assert len(shared) >= 2, "wake sheet does not touch the wing wall"
    sh = mesh.nodes[shared]
    z_junc = junction_z(p)
    d_j = np.linalg.norm(sh - np.array([x_te(z_junc), 0.0, z_junc]),
                         axis=1).min()
    d_t = np.linalg.norm(sh - np.array([x_te(B_SEMI), 0.0, B_SEMI]),
                         axis=1).min()
    assert d_j < 1e-6, f"junction TE not a wall∩wake node (closest {d_j:.2e})"
    assert d_t < 1e-6, f"tip TE corner not a wall∩wake node (closest {d_t:.2e})"

    # The sheet's inboard boundary rides the fuselage waterline: wake nodes
    # shared with the fuselage group, on the revolution surface.
    fus_nodes = set(np.unique(mesh.boundary_faces["fuselage"]).tolist())
    waterline = [n for n in wake_nodes.tolist() if n in fus_nodes]
    assert waterline, "no wake∩fuselage (waterline) nodes -- sheet not stitched"
    wl = mesh.nodes[waterline]
    wlR = np.array([radius_at(p, float(x)) for x in wl[:, 0]])
    assert np.abs(np.abs(wl[:, 2]) - wlR).max() < 0.01 * p.r_f + 3e-3, \
        "waterline nodes off the revolution surface"

    # Aft of the body the sheet reaches the symmetry plane (simple
    # connectivity of the cut) ...
    aft_sym = (wn[:, 2] < 1e-6) & (wn[:, 0] > p.x_tail_tip - 1e-6)
    assert aft_sym.any(), \
        "no wake nodes on the symmetry plane aft of the body -- cut leaves " \
        "the domain multiply connected around body+wake"

    # ... and downstream it reaches the far-field sphere (branch-cut
    # Dirichlet data needs duplicated wake nodes ON the sphere).
    r = np.linalg.norm(wn - np.array([xc, 0.0, 0.0]), axis=1)
    assert r.max() > r_far - 1e-3 * r_far, \
        "wake sheet stops short of the far-field sphere"
