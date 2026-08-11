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

⚠ Registered risk (docs/dev_phase_three/20260811-0300-hex-mesh-prereg.md §2): splitting a hex
introduces DIAGONAL faces, so the face-neighbour graph is not grid-aligned and the upstream-donor
gain this route is meant to deliver may be partly lost in the split. That is measured, not assumed.
And an O-grid wall is still FLAT FACETS: the O(h) facet-normal error does not go away here.
"""

from typing import Dict, Tuple

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
