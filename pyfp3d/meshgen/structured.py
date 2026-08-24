"""Structured body-fitted O-grids for the 2.5-D cases (phase 3 task 2, route A).

WHY THIS EXISTS. Phase two measured that this project's unstructured generator has no
single-variable knob: `h_wall` moved the LE face count +41.7 %, `h_edge` sizes the LE *and*
the TE, `h_far` moved the near-body median cell size +24.04 %, and even the domain radius
`r_far` -- which is not a size field at all -- moved it +67.33 %. A conclusion of the form
"refining X controls Y" cannot be earned with knobs like that, and two conclusions had to be
rewritten because of it (docs/dev_phase_two/20260809-0600 §0b, 20260809-2100).

★★ The single-variable property here is a property of the CONSTRUCTION, not a lucky choice of
knob. Three decisions buy it:

  1. the radial grading is ANCHORED AT THE WALL with a FIXED growth rate, and the layer count
     is DERIVED to reach r_far -- so growing the domain APPENDS layers outward and leaves every
     near-body node bit-identical;
  2. the tangential distribution depends on `n_theta` alone, so wall-normal refinement cannot
     touch the wall triangle count or its tangential spacing;
  3. nothing is fed through a gmsh size field, which is what coupled everything before.

★ The tet split reuses `extrude_single_layer` (prism -> 3 tets, globally consistent
orientation, already locked by tests) rather than adding a new hex->tet splitter: each quad is
cut into two triangles by a fixed diagonal rule, so one hex becomes 6 tets. The STRUCTURE that
route (A) is being tested for -- an orthogonal near-wall layer and independently controllable
grading -- lives in the 2-D quad grid, not in the element type.

⚠ Registered risk (phases/p3/docs/dev_phase_three/20260811-0300-hex-mesh-prereg.md §2): splitting a hex
introduces DIAGONAL faces, so the face-neighbour graph is not grid-aligned and the upstream-donor
gain this route is meant to deliver may be partly lost in the split. That is measured, not assumed.
And an O-grid wall is still FLAT FACETS: the O(h) facet-normal error does not go away here.
"""

from typing import Dict, Optional, Tuple

import numpy as np


def wall_anchored_radii(
    r_wall: float,
    r_far: float,
    h_wall_normal: float,
    growth: float = 1.15,
) -> np.ndarray:
    """Radial node positions: first spacing `h_wall_normal` at the wall, each next spacing
    `growth` times the previous, layers appended until `r_far` is reached or passed.

    ★ This is the function that makes `r_far` a single-variable knob: the sequence is built
    OUTWARD FROM THE WALL, so a larger `r_far` only appends more entries and every earlier
    radius is bit-identical. The final spacing is not stretched to land exactly on `r_far` --
    doing that would make every interior radius depend on `r_far`, which is precisely the
    coupling this construction exists to avoid. The outer boundary therefore sits at the first
    node at or beyond `r_far`, and the achieved radius is reported by the caller.

    ⚠ NOT written as `r_wall + wall_anchored_distances(...)`, and that is deliberate:
    accumulating from `r_wall` and accumulating from 0 then shifting give results that differ in
    the last bits, because floating-point addition is not associative. Measured when the
    "cosmetic" refactor was attempted -- the cylinder's radial sequence stopped being
    bit-identical, which would have silently moved the COMMITTED cylinder meshes and the Q2
    numbers produced from them. Same family as the transport_sigma oracle lesson: the same
    factors in a different ORDER are not the same float. The two functions therefore share the
    construction as a DESCRIPTION, not as an implementation.
    """
    if not (r_far > r_wall > 0.0):
        raise ValueError(f"need r_far > r_wall > 0, got {r_far} and {r_wall}")
    if h_wall_normal <= 0.0 or growth < 1.0:
        raise ValueError(f"need h_wall_normal > 0 and growth >= 1, got "
                         f"{h_wall_normal} and {growth}")
    radii, dr = [r_wall], h_wall_normal
    while radii[-1] < r_far:
        radii.append(radii[-1] + dr)
        dr *= growth
        if len(radii) > 100_000:                     # runaway guard, not a tuning knob
            raise ValueError("radial layer count exceeded 100k -- check growth/h_wall_normal")
    return np.asarray(radii, dtype=np.float64)


def wall_anchored_distances(
    d_far: float,
    h_wall_normal: float,
    growth: float = 1.15,
) -> np.ndarray:
    """The same construction expressed as DISTANCES FROM THE WALL, starting at exactly 0.

    ★ Extracted so the cylinder (which wants radii) and the airfoil (which wants wall
    distances, and whose wall is not a circle so there is no radius to add) share ONE
    construction and therefore one prefix-stability guarantee. `wall_anchored_radii` is now a
    thin shift of this, and the refactor was accepted only after proving its output stays
    BIT-IDENTICAL for the cylinder -- the G0 evidence depends on that sequence.
    """
    if d_far <= 0.0:
        raise ValueError(f"need d_far > 0, got {d_far}")
    if h_wall_normal <= 0.0 or growth < 1.0:
        raise ValueError(f"need h_wall_normal > 0 and growth >= 1, got "
                         f"{h_wall_normal} and {growth}")
    dist, dr = [0.0], h_wall_normal
    while dist[-1] < d_far:
        dist.append(dist[-1] + dr)
        dr *= growth
        if len(dist) > 100_000:                      # runaway guard, not a tuning knob
            raise ValueError("radial layer count exceeded 100k -- check growth/h_wall_normal")
    return np.asarray(dist, dtype=np.float64)


def cylinder_o_grid_2d(
    radius: float = 1.0,
    r_far: float = 20.0,
    n_theta: int = 96,
    h_wall_normal: float = 0.02,
    growth: float = 1.15,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, float]]:
    """Structured O-grid annulus around a cylinder, returned as a TRIANGULATION.

    Knob scope, which is the whole point of this module:
        n_theta        -> wall tangential resolution ONLY
        h_wall_normal  -> wall-normal resolution ONLY (and the grading anchored on it)
        growth         -> grading rate ONLY
        r_far          -> domain size ONLY (appends outward layers)

    Returns (points2d, triangles, edge_groups, info) where `edge_groups` has "wall" (inner
    circle) and "farfield" (outer), matching `cylinder_annulus_2d`'s contract so the
    downstream pipeline is unchanged. `info` reports the achieved outer radius and the layer
    count, both DERIVED rather than requested.
    """
    radii = wall_anchored_radii(radius, r_far, h_wall_normal, growth)
    n_r = len(radii)
    theta = np.arange(n_theta, dtype=np.float64) * (2.0 * np.pi / n_theta)
    ct, st = np.cos(theta), np.sin(theta)

    #: node (i, j) = radii[j] * (cos t_i, sin t_i); index = j * n_theta + i, theta periodic
    pts = np.empty((n_r * n_theta, 2), dtype=np.float64)
    for j, r in enumerate(radii):
        pts[j * n_theta:(j + 1) * n_theta, 0] = r * ct
        pts[j * n_theta:(j + 1) * n_theta, 1] = r * st

    def idx(i, j):
        return j * n_theta + (i % n_theta)

    #: ★ each quad -> 2 triangles on a FIXED diagonal (i,j)-(i+1,j+1). A fixed rule (rather
    #: than a shorter-diagonal choice) keeps the split independent of the spacings, so the
    #: connectivity does not change when a resolution knob moves -- the G0 self-check reads
    #: bit-identity, and a spacing-dependent diagonal would break it for no benefit.
    tris = np.empty((2 * n_theta * (n_r - 1), 3), dtype=np.int64)
    k = 0
    for j in range(n_r - 1):
        for i in range(n_theta):
            a, b = idx(i, j), idx(i + 1, j)
            c, d = idx(i + 1, j + 1), idx(i, j + 1)
            tris[k] = (a, b, c); k += 1
            tris[k] = (a, c, d); k += 1

    wall = np.stack([np.arange(n_theta), (np.arange(n_theta) + 1) % n_theta], axis=1)
    top = (n_r - 1) * n_theta
    far = np.stack([top + np.arange(n_theta),
                    top + (np.arange(n_theta) + 1) % n_theta], axis=1)
    info = dict(n_theta=n_theta, n_radial=n_r, r_achieved=float(radii[-1]),
                h_wall_normal=float(h_wall_normal), growth=float(growth),
                h_wall_tangential=float(2.0 * np.pi * radius / n_theta))
    return pts, tris.astype(np.int64), {"wall": wall.astype(np.int64),
                                        "farfield": far.astype(np.int64)}, info


def _outward_normals_closed(pts: np.ndarray) -> np.ndarray:
    """Outward unit normals of a CLOSED polyline (n, 2), one per vertex.

    Averages the two adjacent edge normals, then fixes the global sign so the normals point
    AWAY from the centroid. Averaging (rather than taking one edge's normal) is what keeps the
    first offset layer smooth through the LE, where the edge directions turn fastest.
    """
    n = len(pts)
    e = np.roll(pts, -1, axis=0) - pts                       # edge i: pts[i] -> pts[i+1]
    en = np.stack([e[:, 1], -e[:, 0]], axis=1)               # right-hand normal per edge
    en /= np.linalg.norm(en, axis=1, keepdims=True)
    v = en + np.roll(en, 1, axis=0)                           # average with previous edge
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    c = pts.mean(axis=0)
    if np.mean(np.einsum("ij,ij->i", v, pts - c)) < 0.0:
        v = -v
    return v


def airfoil_o_grid_2d(
    surface: np.ndarray,
    r_far: float = 20.0,
    h_wall_normal: float = 0.004,
    growth: float = 1.15,
    circularize_from: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray],
             Dict[str, np.ndarray], Dict[str, float]]:
    """Structured O-grid around a CLOSED airfoil polyline, returned as a triangulation.

    Returns (points2d, triangles, edge_groups, interior_edge_groups, info) -- five values, one
    more than the cylinder, because a lifting airfoil needs the tagged wake sheet (see below).

    `surface` is (n_theta, 2), closed (first != last; the wrap is implicit) -- e.g.
    `planar.naca0012_coordinates(n_half)[:-1]`. Layers march outward along averaged vertex
    normals and are BLENDED toward a circle as they go, because a pure normal offset
    self-intersects in the concave region behind a sharp trailing edge. The blend weight is a
    function of the ACCUMULATED WALL DISTANCE only.

    ★ That last sentence is what preserves the single-variable property proved by G0: the radii
    come from `wall_anchored_radii`, which is prefix-stable, and the blend depends on the same
    prefix-stable accumulated distance -- so enlarging `r_far` appends layers and leaves every
    earlier node bit-identical, exactly as on the cylinder. Had the blend been written as
    "fraction of the way to r_far" it would have coupled every interior layer to the domain
    size, recreating the very defect this module exists to avoid.

    `circularize_from` is the wall distance (in chords) at which the blend reaches 1, i.e.
    where the layers have become circles.
    """
    surf = np.asarray(surface, dtype=np.float64)
    if surf.ndim != 2 or surf.shape[1] != 2:
        raise ValueError(f"surface must be (n, 2), got {surf.shape}")
    n_theta = len(surf)
    nrm = _outward_normals_closed(surf)
    centre = surf.mean(axis=0)
    #: wall-anchored radial spacings, reused verbatim from the cylinder construction
    r_seq = wall_anchored_distances(r_far, h_wall_normal, growth)
    n_r = len(r_seq)

    pts = np.empty((n_r * n_theta, 2), dtype=np.float64)
    ang = np.arctan2(surf[:, 1] - centre[1], surf[:, 0] - centre[0])
    for j, d in enumerate(r_seq):
        w = min(1.0, d / circularize_from) if circularize_from > 0 else 1.0
        offset = surf + nrm * d
        radius = d + 0.5                     # circle radius grows with the wall distance
        circle = centre + np.stack([radius * np.cos(ang), radius * np.sin(ang)], axis=1)
        pts[j * n_theta:(j + 1) * n_theta] = (1.0 - w) * offset + w * circle

    def idx(i, j):
        return j * n_theta + (i % n_theta)

    tris = np.empty((2 * n_theta * (n_r - 1), 3), dtype=np.int64)
    k = 0
    for j in range(n_r - 1):
        for i in range(n_theta):
            a, b = idx(i, j), idx(i + 1, j)
            c, d_ = idx(i + 1, j + 1), idx(i, j + 1)
            tris[k] = (a, b, c); k += 1
            tris[k] = (a, c, d_); k += 1

    wall = np.stack([np.arange(n_theta), (np.arange(n_theta) + 1) % n_theta], axis=1)
    top = (n_r - 1) * n_theta
    far = np.stack([top + np.arange(n_theta),
                    top + (np.arange(n_theta) + 1) % n_theta], axis=1)
    #: ★★ THE WAKE IS A GRID LINE. A lifting airfoil needs a tagged interior sheet for
    #: `cut_wake` to duplicate nodes across (the circulation branch cut). On an O-grid this has
    #: a natural home: `naca0012_coordinates` starts at the TE, so column i = 0 is the radial
    #: line running from the TE straight out to the far field -- the seam of the O-grid IS the
    #: wake. Nothing about the topology changes; this is TAGGING, the same thing
    #: `planar.naca0012_wake_2d` does for the unstructured family, except there the sheet has to
    #: be cut through arbitrary triangles while here it follows the grid.
    wake = np.stack([np.arange(n_r - 1) * n_theta,
                     (np.arange(n_r - 1) + 1) * n_theta], axis=1)
    seg = np.linalg.norm(np.roll(surf, -1, axis=0) - surf, axis=1)
    info = dict(n_theta=n_theta, n_radial=n_r, h_wall_normal=float(h_wall_normal),
                growth=float(growth), d_achieved=float(r_seq[-1]),
                h_wall_tangential_min=float(seg.min()),
                h_wall_tangential_max=float(seg.max()),
                circularize_from=float(circularize_from))
    return (pts, tris.astype(np.int64),
            {"wall": wall.astype(np.int64), "farfield": far.astype(np.int64)},
            {"wake": wake.astype(np.int64)}, info)

def airfoil_layer_block_2d(
    surface: np.ndarray,
    thickness: float = 0.08,
    h_wall_normal: float = 0.004,
    growth: float = 1.15,
    te_blend: float = 0.05,
):
    """A FINITE-HEIGHT structured layer block hugging a SHARP-TE airfoil (route A, v2).

    ★★ Replaces the v1 O-grid, which stretched ONE grid from the wall all the way to r_far and
    was measured bad: aspect ratio up to 1305, 574 cells above 100, min cell area 1.1e-07 behind
    the TE, and a polygonal far-field boundary (its outer ring inherited the airfoil's cosine
    clustering). The v1 orientation check said "consistent" and passed it -- orientation is not
    quality, and that guard was narrower than the conclusion it was used to support.

    Architecture (user-specified, and it is the production one): grow structured layers from the
    wall to a FINITE height, stop, and let an unstructured mesher fill the rest. A thin block
    keeps the aspect ratio controlled, needs no circularisation (so the LE is never distorted),
    and leaves the far field to a mesher that does it properly.

    ★★ SHARP-TE CONSTRAINT (user-specified). `naca0012_coordinates` uses the closed-TE
    coefficient set, so the TE is a single sharp point at (1, 0) with ZERO thickness. There the
    averaged vertex normal points along +x -- i.e. marching it would send the TE column straight
    DOWN THE WAKE LINE, which is degenerate, and the neighbouring upper/lower rays then cross
    just behind the TE. So the marching direction is blended toward +y on the upper side and -y
    on the lower side over the last `te_blend` of chord, reaching pure +/-y exactly at the TE,
    and the TE node is SPLIT into an upper and a lower copy. The block is therefore C-shaped and
    the wake line is its seam -- which is where the wake belongs anyway.

    Returns (points2d, triangles, outer_polyline_idx, wall_idx, info). The caller feeds
    `outer_polyline_idx` to the unstructured mesher as an inner boundary.
    """
    surf = np.asarray(surface, dtype=np.float64)
    if surf.ndim != 2 or surf.shape[1] != 2:
        raise ValueError(f"surface must be (n, 2), got {surf.shape}")
    #: open the closed polyline at the TE: index 0 is the TE on the UPPER side, and a duplicate
    #: TE node is appended for the LOWER side.
    open_surf = np.vstack([surf, surf[0]])
    n_s = len(open_surf)

    #: vertex normals on an OPEN polyline (no wrap): average the adjacent edge normals, and at
    #: the two ends use the single adjacent edge.
    e = np.diff(open_surf, axis=0)
    en = np.stack([e[:, 1], -e[:, 0]], axis=1)
    en /= np.linalg.norm(en, axis=1, keepdims=True)
    nrm = np.empty_like(open_surf)
    nrm[0] = en[0]
    nrm[-1] = en[-1]
    nrm[1:-1] = en[:-1] + en[1:]
    nrm[1:-1] /= np.linalg.norm(nrm[1:-1], axis=1, keepdims=True)
    c = surf.mean(axis=0)
    if np.mean(np.einsum("ij,ij->i", nrm, open_surf - c)) < 0.0:
        nrm = -nrm

    #: ★ TE blend by CHORDWISE distance, width `te_blend` (user ruling 2026-08-11). The side
    #: sign comes from the node ORDER, not from y: the TE nodes have y = 0 exactly, so a
    #: sign(y) rule would be undefined at precisely the station that matters most.
    x = open_surf[:, 0]
    x0 = x.max() - te_blend
    w = np.clip((x - x0) / te_blend, 0.0, 1.0) if te_blend > 0 else np.zeros(n_s)
    i_le = int(np.argmin(x))
    side = np.where(np.arange(n_s) <= i_le, 1.0, -1.0)     # first half = upper, second = lower
    target = np.stack([np.zeros(n_s), side], axis=1)
    d = (1.0 - w[:, None]) * nrm + w[:, None] * target
    d /= np.linalg.norm(d, axis=1, keepdims=True)

    dist = wall_anchored_distances(thickness, h_wall_normal, growth)
    n_r = len(dist)
    pts = (open_surf[None, :, :] + d[None, :, :] * dist[:, None, None]).reshape(-1, 2)

    def idx(i, j):
        return j * n_s + i

    tris = np.empty((2 * (n_s - 1) * (n_r - 1), 3), dtype=np.int64)
    k = 0
    for j in range(n_r - 1):
        for i in range(n_s - 1):
            a, b = idx(i, j), idx(i + 1, j)
            cc, dd = idx(i + 1, j + 1), idx(i, j + 1)
            tris[k] = (a, b, cc); k += 1
            tris[k] = (a, cc, dd); k += 1

    wall_idx = np.arange(n_s, dtype=np.int64)
    outer_idx = np.arange(n_s, dtype=np.int64) + (n_r - 1) * n_s
    seg = np.linalg.norm(np.diff(open_surf, axis=0), axis=1)
    info = dict(n_stations=n_s, n_layers=n_r, thickness_achieved=float(dist[-1]),
                h_wall_normal=float(h_wall_normal), growth=float(growth),
                te_blend=float(te_blend), te_split=True,
                h_wall_tangential_min=float(seg.min()),
                h_wall_tangential_max=float(seg.max()))
    return pts, tris, outer_idx, wall_idx, info


def layer_block_quality(pts: np.ndarray, tris: np.ndarray, wall_idx: np.ndarray,
                        n_stations: int, n_layers: int) -> Dict[str, float]:
    """The four checks v1 did NOT have, and whose absence let a bad mesh be measured on.

    ray_crossings   self-intersecting ("bow-tie") quads -- the failure the sharp-TE constraint
                    exists to prevent. Tested by same-sign quad triangles, which needs no
                    handedness convention; see the note at the test itself for the version that
                    got this wrong
    aspect_*        cell aspect ratios (v1 reached 1305)
    min_area        smallest cell (v1 hit 1.1e-07 behind the TE)
    orthogonality   deviation, in degrees, between the first marching step and the wall normal --
                    the price the TE blend charges, made visible instead of assumed small
    """
    P = pts.reshape(n_layers, n_stations, 2)
    step = P[1] - P[0]
    step_u = step / np.linalg.norm(step, axis=1, keepdims=True)
    tan = np.zeros_like(step_u)
    tan[:-1] = P[0, 1:] - P[0, :-1]
    tan[-1] = tan[-2]
    tan /= np.linalg.norm(tan, axis=1, keepdims=True)
    dev = np.degrees(np.abs(np.arcsin(np.clip(
        np.abs(np.einsum("ij,ij->i", step_u, tan)), -1.0, 1.0))))

    #: ★★ SELF-INTERSECTION, by an UNAMBIGUOUS test: for every quad (wall pair, outer pair) the
    #: two triangles of its fixed diagonal must have SAME-SIGN areas. A self-intersecting
    #: ("bow-tie") quad is exactly the case where they differ, and this needs no convention.
    #:
    #: The first version compared the tangential edge at the wall with the one at each layer and
    #: counted sign flips of their cross product. It reported 800 crossings on a block the
    #: unambiguous test finds CLEAN (0 of 1600): it was reading the LOWER surface's normal
    #: handedness as a crossing. A checker that cries wolf is worse than none -- it would have
    #: sent me to "fix" a mesh that was already fine -- so it was replaced and then validated on
    #: two cases whose answers were already known.
    n_cross = 0
    for j in range(n_layers - 1):
        a0 = P[j, :-1]; b0 = P[j, 1:]; c0 = P[j + 1, 1:]; d0 = P[j + 1, :-1]
        s1 = (b0[:, 0] - a0[:, 0]) * (c0[:, 1] - a0[:, 1]) - \
             (c0[:, 0] - a0[:, 0]) * (b0[:, 1] - a0[:, 1])
        s2 = (c0[:, 0] - a0[:, 0]) * (d0[:, 1] - a0[:, 1]) - \
             (d0[:, 0] - a0[:, 0]) * (c0[:, 1] - a0[:, 1])
        n_cross += int(np.count_nonzero(np.sign(s1) != np.sign(s2)))
    p = pts[tris]
    ee = np.stack([np.linalg.norm(p[:, 1] - p[:, 0], axis=1),
                   np.linalg.norm(p[:, 2] - p[:, 1], axis=1),
                   np.linalg.norm(p[:, 0] - p[:, 2], axis=1)], axis=1)
    ar = ee.max(axis=1) / np.maximum(ee.min(axis=1), 1e-300)
    a2 = np.abs((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
                - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1])) / 2.0
    return dict(ray_crossings=n_cross, aspect_median=float(np.median(ar)),
                aspect_p95=float(np.percentile(ar, 95)), aspect_max=float(ar.max()),
                min_area=float(a2.min()),
                orthogonality_dev_max_deg=float(dev.max()),
                orthogonality_dev_median_deg=float(np.median(dev)))

def n_layers_of(info):
    return info["n_layers"]


def airfoil_hybrid_2d(
    surface: np.ndarray,
    thickness: float = 0.08,
    h_wall_normal: float = 0.004,
    growth: float = 1.15,
    te_blend: float = 0.05,
    r_far: float = 20.0,
    h_far: float = 2.0,
    h_wake: Optional[float] = None,
    wake_dist_max: float = 1.5,
):
    """HYBRID 2-D airfoil mesh: structured wall layers + unstructured far field.

    The architecture the user specified, and the production one: grow structured layers from the
    wall to a FINITE height, stop, and let gmsh fill the rest. What it fixes relative to the v1
    O-grid, measured: aspect ratio max 1305 -> 47, min cell area 1.1e-07 -> 7.6e-07, and the
    polygonal far-field boundary is gone because the far field is no longer one stretched ring
    inheriting the airfoil's cosine clustering.

    ★★ NODE MATCHING is the part that has to be exact or the two regions do not join. Every
    segment of the block's outer boundary is added to gmsh as its own Line with
    `setTransfiniteCurve(line, 2)`, i.e. exactly two nodes -- its endpoints. gmsh therefore
    places boundary nodes at PRECISELY the block's points, and the merge is a coordinate lookup
    rather than an interpolation. Without the transfinite constraint gmsh would subdivide those
    lines wherever its size field asked for something finer, and the seam would be
    non-conforming with no error raised.

    ★ The wake line is embedded in the UNSTRUCTURED region (it runs downstream from the TE, which
    is where the block opens), so `cut_wake` gets its interior sheet exactly as on the existing
    unstructured family.

    Returns (points2d, triangles, edge_groups, interior_edge_groups, info).
    """
    import gmsh

    if h_wake is None:
        h_wake = 4.0 * h_wall_normal
    blk_pts, blk_tris, outer_idx, wall_idx, binfo = airfoil_layer_block_2d(
        surface, thickness=thickness, h_wall_normal=h_wall_normal, growth=growth,
        te_blend=te_blend)
    #: ★★ WELD THE TWO TE SURFACE NODES. The block SPLITS the TE so the upper and lower sides
    #: can march in +y and -y (the sharp-TE constraint), but after marching those two SURFACE
    #: nodes are at identical coordinates, and keeping them separate is a topology claim the
    #: geometry does not support. Left unwelded, gmsh has ONE TE point while the block has two,
    #: the merge can only attach the far field to one of them, and the seam is non-conforming --
    #: measured as 246 boundary edges against the 224 that wall+farfield account for. Splitting
    #: was a CONSTRUCTION need, not a topological one; `cut_wake` duplicates the TE again later
    #: for the wake itself.
    i_te_lower = binfo["n_stations"] - 1
    weld = np.arange(len(blk_pts), dtype=np.int64)
    weld[i_te_lower] = 0                                    # lower TE surface node -> upper
    keep = np.ones(len(blk_pts), dtype=bool); keep[i_te_lower] = False
    renum = np.cumsum(keep) - 1
    blk_pts = blk_pts[keep]
    blk_tris = renum[weld[blk_tris]]
    outer_idx = renum[weld[outer_idx]]
    wall_idx = renum[weld[wall_idx]]
    binfo = dict(binfo); binfo["te_welded"] = True

    outer = blk_pts[outer_idx]                     # open curve: TE-upper-outer -> nose -> TE-lower
    te = np.asarray(surface, dtype=np.float64)[0]  # the sharp TE point (1, 0)
    #: ★★ the two TE COLUMNS, node by node. The block subdivides each column into one edge per
    #: layer; the first version handed gmsh a single straight Line for each with two nodes, so
    #: the far field met the block along ONE edge where the block has ten. Measured as exactly 22
    #: extra boundary edges, all at x = 1.0 with |y| <= 0.0812 -- i.e. the two columns. The
    #: same-subdivision rule that was already applied to the outer curve has to apply here too,
    #: which is the whole content of "the seam must be conforming".
    n_s_w = binfo["n_stations"]
    col_up = renum[weld[np.arange(1, n_layers_of(binfo)) * n_s_w + 0]]
    col_lo = renum[weld[np.arange(1, n_layers_of(binfo)) * n_s_w + (n_s_w - 1)]]

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("airfoil_hybrid")
        geo = gmsh.model.geo

        #: the inner loop = upper TE column, the block's outer curve, lower TE column
        p_te = geo.addPoint(te[0], te[1], 0.0, h_wake)
        #: the upper column's interior nodes, then the outer curve, then the lower column's
        #: interior nodes reversed -- every segment its own Line so the subdivision MATCHES the
        #: block's (see the note where col_up/col_lo are built).
        p_col_up = [geo.addPoint(*blk_pts[i], 0.0) for i in col_up[:-1]]
        p_outer = [geo.addPoint(x, y, 0.0) for x, y in outer]
        p_col_lo = [geo.addPoint(*blk_pts[i], 0.0) for i in col_lo[:-1]]
        chain = [p_te] + p_col_up + p_outer + p_col_lo[::-1] + [p_te]
        lines = [geo.addLine(chain[i], chain[i + 1]) for i in range(len(chain) - 1)]
        inner_loop = geo.addCurveLoop(lines)

        cx = 0.5
        ctr = geo.addPoint(cx, 0.0, 0.0, h_far)
        pe = geo.addPoint(cx + r_far, 0.0, 0.0, h_far)
        pn = geo.addPoint(cx, r_far, 0.0, h_far)
        pw = geo.addPoint(cx - r_far, 0.0, 0.0, h_far)
        ps = geo.addPoint(cx, -r_far, 0.0, h_far)
        arcs = [geo.addCircleArc(pe, ctr, pn), geo.addCircleArc(pn, ctr, pw),
                geo.addCircleArc(pw, ctr, ps), geo.addCircleArc(ps, ctr, pe)]
        far_loop = geo.addCurveLoop(arcs)
        surf_tag = geo.addPlaneSurface([far_loop, inner_loop])

        #: wake line: TE -> the far-field boundary's east point, embedded so triangle edges
        #: conform to it. ★ It ENDS at `pe`, which already exists on the far-field loop; the
        #: first version also created a separate point at the same location and never used it,
        #: leaving a dangling gmsh point that showed up as 2 orphan nodes and 4 duplicate
        #: coordinates in the merged mesh -- found by the validity checks, not by inspection.
        wake_line = geo.addLine(p_te, pe)
        geo.synchronize()
        gmsh.model.mesh.embed(1, [wake_line], 2, surf_tag)

        #: ★ exactly two nodes per boundary segment -> gmsh nodes coincide with block nodes
        for ln in lines:
            gmsh.model.mesh.setTransfiniteCurve(ln, 2)

        fd = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(fd, "CurvesList", lines + [wake_line])
        ft = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(ft, "InField", fd)
        #: ★ the far-field size AT THE SEAM must continue the block's OUTER LAYER spacing, not a
        #: fraction of the block thickness. The first version used thickness*0.6 = 0.048 against
        #: an outer-layer spacing of 0.0141 -- a 3.4x jump visible as a size discontinuity right
        #: across the seam. The block knows its own last step; ask it.
        h_seam = float(np.diff(wall_anchored_distances(thickness, h_wall_normal, growth))[-1])
        gmsh.model.mesh.field.setNumber(ft, "SizeMin", max(h_wake, h_seam))
        gmsh.model.mesh.field.setNumber(ft, "SizeMax", h_far)
        gmsh.model.mesh.field.setNumber(ft, "DistMin", thickness)
        gmsh.model.mesh.field.setNumber(ft, "DistMax", wake_dist_max)
        gmsh.model.mesh.field.setAsBackgroundMesh(ft)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.model.mesh.generate(2)

        ntags, ncoord, _ = gmsh.model.mesh.getNodes()
        gpts = ncoord.reshape(-1, 3)[:, :2]
        gidx = {int(t): i for i, t in enumerate(ntags)}
        etypes, etags, enodes = gmsh.model.mesh.getElements(2, surf_tag)
        gtris = None
        for et, en in zip(etypes, enodes):
            if et == 2:
                gtris = np.array([gidx[int(t)] for t in en], dtype=np.int64).reshape(-1, 3)
        if gtris is None:
            raise RuntimeError("gmsh produced no triangles for the far field")
        far_edges = []
        for arc in arcs:
            _, _, en = gmsh.model.mesh.getElements(1, arc)
            if en:
                far_edges.append(np.array([gidx[int(t)] for t in en[0]],
                                          dtype=np.int64).reshape(-1, 2))
        _, _, wen = gmsh.model.mesh.getElements(1, wake_line)
        wake_edges = (np.array([gidx[int(t)] for t in wen[0]], dtype=np.int64).reshape(-1, 2)
                      if wen else np.zeros((0, 2), dtype=np.int64))
    finally:
        gmsh.finalize()

    #: ---- merge: block nodes first, then gmsh nodes minus the shared seam ------------------
    tol = 1e-9
    key = lambda q: (round(float(q[0]) / tol), round(float(q[1]) / tol))
    blk_key = {key(q): i for i, q in enumerate(blk_pts)}
    remap = np.full(len(gpts), -1, dtype=np.int64)
    extra = []
    for i, q in enumerate(gpts):
        j = blk_key.get(key(q))
        if j is not None:
            remap[i] = j
        else:
            remap[i] = len(blk_pts) + len(extra)
            extra.append(q)
    n_shared = int(np.count_nonzero(remap < len(blk_pts)))
    pts = np.vstack([blk_pts, np.asarray(extra).reshape(-1, 2)]) if extra else blk_pts
    tris = np.vstack([blk_tris, remap[gtris]])
    edge_groups = {"wall": np.stack([wall_idx[:-1], wall_idx[1:]], axis=1),
                   "farfield": remap[np.vstack(far_edges)]}
    interior = {"wake": remap[wake_edges]} if len(wake_edges) else {}

    #: ★ COMPACT: gmsh returns its geometry points too, and the far-field circle's CENTRE is one
    #: of them -- `addCircleArc` needs it, no element uses it. `extrude_single_layer` refuses a
    #: 2-D mesh with unreferenced nodes ("compact the 2D mesh before extruding"), which is the
    #: right guard in the right place; compacting belongs here, in the generator that produced
    #: the stray node.
    used = np.unique(tris)
    if len(used) != len(pts):
        keep2 = np.zeros(len(pts), dtype=bool); keep2[used] = True
        renum2 = np.cumsum(keep2) - 1
        pts = pts[keep2]
        tris = renum2[tris]
        edge_groups = {k: renum2[v] for k, v in edge_groups.items()}
        interior = {k: renum2[v] for k, v in interior.items()}
    info = dict(binfo)
    info.update(n_block_nodes=len(blk_pts), n_block_tris=len(blk_tris),
                n_far_nodes=len(extra), n_far_tris=len(gtris),
                n_seam_nodes_matched=n_shared, n_outer_boundary=len(outer_idx),
                n_nodes_after_compaction=len(pts), r_far=float(r_far), h_far=float(h_far))
    return pts, tris.astype(np.int64), edge_groups, interior, info

def airfoil_surface_distribution(
    coords: np.ndarray,
    n_stations: int,
    le_cluster: float = 0.8,
    local_window: Optional[Tuple[float, float]] = None,
    local_factor: float = 1.0,
) -> np.ndarray:
    """Resample a closed airfoil polyline: CLUSTER AT THE LE, COARSEN AT THE TE.

    ★★ Why this exists (user, 2026-08-11): `naca0012_coordinates` uses cosine clustering, which
    refines BOTH ends -- measured wall tangential spacing 0.000390 at the LE *and* at the TE.
    LE clustering is wanted (high curvature, suction peak). TE clustering is not: on a sharp-TE
    airfoil that region is geometrically almost straight, the Kutta condition is carried by the
    wake, and the clustering only buys 7.6e-07 cells, drives the aspect ratio, and -- worst --
    interacts with the TE +/-y lock, since every one of those densely packed last-5%-of-chord
    stations is rotating its marching direction at once.

    The mapping is `s(t) = t + a*sin(2*pi*t)/(2*pi)` over the closed arclength parameter, so
    `ds/dt = 1 + a*cos(2*pi*t)`: at the LE (t = 1/2) the spacing is (1 - a) x uniform, at the TE
    (t = 0, 1) it is (1 + a) x uniform. ONE knob, and it moves the two ends in OPPOSITE
    directions, which is exactly the asymmetry the geometry calls for.

    ⚠ This changes the wall point distribution relative to the unstructured family, which is
    generated straight from `naca0012_coordinates`. Any A/B between the families must say so:
    matching the wall DOF COUNT no longer implies matching their placement.

    ★★ `local_window = (x0, x1)` with `local_factor > 1` adds LOCAL chordwise refinement inside
    that chord band on BOTH surfaces, and is the one capability task 3 needs
    (phases/p3/docs/dev_phase_three/20260811-1500-task3-refinement-paradox-prereg.md): refine ONLY the
    shock band and leave everything else alone, so "the shock region's resolution" can be
    separated from "the global cell count". Phase two could never run that leg -- all four of its
    generator knobs were measured out of scope -- so it could only refine globally and could not
    attribute the result.

    The refinement is applied by warping the arclength parameter, so the station COUNT is
    unchanged: points are taken FROM outside the window and given TO it. That is deliberate --
    holding the count fixed means the leg cannot be confounded with a global DOF change, which is
    the whole point of the comparison. ⚠ The corollary, stated because it is a real limitation:
    outside the window the spacing gets COARSER, so this knob is not "window-only" in the
    bit-identical sense G0 demands of the others; what it guarantees is a fixed count and a
    monotone, measurable redistribution. The task-3 guard therefore checks the SIGN and the
    magnitude of the change inside and outside, not bit-identity, and says so.
    """
    c = np.asarray(coords, dtype=np.float64)
    if np.allclose(c[0], c[-1]):
        c = c[:-1]
    if not 0.0 <= le_cluster < 1.0:
        raise ValueError(f"le_cluster must be in [0, 1), got {le_cluster}")
    #: cumulative arclength of the closed polyline, normalised, starting at the TE
    seg = np.linalg.norm(np.diff(np.vstack([c, c[0]]), axis=0), axis=1)
    s_ref = np.concatenate([[0.0], np.cumsum(seg)])
    s_ref /= s_ref[-1]
    ref = np.vstack([c, c[0]])
    #: ★ anchor t = 1/2 on the ACTUAL LE (minimum x), not on half the arclength: the upper and
    #: lower runs are equal here but that is a property of a symmetric section, not a guarantee.
    i_le = int(np.argmin(c[:, 0]))
    s_le = s_ref[i_le]
    t = np.linspace(0.0, 1.0, n_stations + 1)[:-1]
    a = le_cluster
    u = t + a * np.sin(2.0 * np.pi * t) / (2.0 * np.pi)      # in [0, 1], clustered at u = 1/2
    #: stretch u so that u = 1/2 lands on the real LE arclength
    s = np.where(u <= 0.5, u * (s_le / 0.5), s_le + (u - 0.5) * ((1.0 - s_le) / 0.5))
    #: ★ local chordwise refinement: warp s toward the window's arclength interval. Implemented
    #: on the ARC parameter (not on x) so both surfaces are treated identically and the TE stays
    #: put; `local_factor` is the density multiplier inside the window.
    if local_window is not None and local_factor != 1.0:
        x0, x1 = float(local_window[0]), float(local_window[1])
        if not (x1 > x0):
            raise ValueError(f"local_window must be increasing, got {local_window}")
        if local_factor <= 0.0:
            raise ValueError(f"local_factor must be > 0, got {local_factor}")
        #: density rho(s) = 1 inside the window scaled by local_factor, 1 elsewhere; then s is
        #: replaced by the inverse of the normalised cumulative density, evaluated on a fine grid
        sg = np.linspace(0.0, 1.0, 20001)
        xg = np.interp(sg, s_ref, ref[:, 0])
        rho = np.where((xg >= x0) & (xg <= x1), local_factor, 1.0)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (rho[1:] + rho[:-1]) * np.diff(sg))])
        cum /= cum[-1]
        s = np.interp(s, cum, sg)                            # more points where rho is larger
    out = np.empty((n_stations, 2), dtype=np.float64)
    out[:, 0] = np.interp(s, s_ref, ref[:, 0])
    out[:, 1] = np.interp(s, s_ref, ref[:, 1])
    out[0] = c[0]                                            # keep the TE point exact
    return out

