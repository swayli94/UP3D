"""
M2 deliverable (part 1 of 2): a simplified axisymmetric fuselage, as an OCC
solid of revolution, for the ONERA M6 wing-body validation mesh
(roadmap Track M / M2; the solver leg is Track B / B9, level-set path).

Why a body of revolution built from ONE splined meridian rather than the
fused nose-ellipsoid + cylinder + cone + tail-sphere of the plan sketch:
fusing four primitives leaves C0 seams at every junction (the cylinder->cone
"shoulder" especially), and in potential flow a surface crease is a spurious
edge feature -- exactly the flat-tip-cap problem M5 removed. Splining a single
meridian profile through the same piecewise radius law and revolving it gives
a C2 skin with no seam creases, so the only sharp features on the whole
wing-body are the ones we *want* (the wing LE/TE, which carry the physics).

Axis convention (shared with wing3d.py / mesh.wake_cut, fixed project-wide):
chord along +x, lift along +y, span along +z. The fuselage is a body of
revolution about the +x axis, centered on the aircraft centerline
(y, z) = (0, 0). The aircraft symmetry plane is z = 0 (the wing root plane),
so the half model keeps z >= 0 and the fuselage contributes its z >= 0 half
skin to the fluid boundary, meeting the symmetry plane along its z = 0
meridian rims.

The default parameters place the wing (root chord x in [0, 0.806],
B_SEMI = 1.1963) so that the wing root is buried in the constant-radius
cylinder section and emerges at the wing-fuselage junction near z = r_f. The
wing TE is at y = 0 (the ONERA D section has zero TE thickness), so the TE
line meets the fuselage skin at z = r_f exactly -- this is where the level-set
wake polyline and the innermost Kutta station start (see wingbody.py).

Gmsh is imported lazily so the solver test suite does not depend on it; this
module's geometry helpers (FuselageParams, radius_at, profile_x) are pure
numpy and safe to import anywhere.
"""

import math
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np


@dataclass(frozen=True)
class FuselageParams:
    """Simplified fuselage: rounded nose (half ellipsoid), constant-radius
    body (cylinder), tapering afterbody (cone), rounded tail (sphere cap).

    All lengths in meters, on the wing3d planform scale (C_ROOT = 0.806,
    B_SEMI = 1.1963). Defaults are the plan-sketch values; r_f ~ 0.19 C_ROOT
    puts the wing-fuselage junction near z = r_f = 0.15.
    """
    r_f: float = 0.15          # body radius (constant cylinder section)
    l_nose: float = 0.30       # nose length (rounded fore-body)
    x_body_end: float = 1.20   # end of the constant-radius section
    l_tail: float = 0.60       # afterbody (cone) length
    r_tail: float = 0.03       # radius at the start of the tail sphere cap

    def __post_init__(self) -> None:
        for name in ("r_f", "l_nose", "x_body_end", "l_tail", "r_tail"):
            if getattr(self, name) <= 0.0:
                raise ValueError(f"FuselageParams.{name} must be > 0")
        if self.r_tail >= self.r_f:
            raise ValueError("r_tail must be smaller than r_f (afterbody tapers)")

    @property
    def x_nose_tip(self) -> float:
        return -self.l_nose

    @property
    def x_tail_start(self) -> float:
        """Where the cone ends and the tail sphere cap begins."""
        return self.x_body_end + self.l_tail

    @property
    def x_tail_tip(self) -> float:
        return self.x_tail_start + self.r_tail

    @property
    def length(self) -> float:
        return self.x_tail_tip - self.x_nose_tip


def radius_at(p: FuselageParams, x: float) -> float:
    """Meridian radius R(x) of the (pre-spline) piecewise profile.

    - nose  [-l_nose, 0]      : half ellipsoid, R = r_f * sqrt(1 - (x/l_nose)^2)
    - body  [0, x_body_end]   : cylinder, R = r_f
    - tail  [x_body_end, x_tail_start] : cone, r_f -> r_tail (linear)
    - cap   [x_tail_start, x_tail_tip] : sphere, R = sqrt(r_tail^2 - dx^2)
    R = 0 outside the body (on the axis at both tips).
    """
    if x <= p.x_nose_tip or x >= p.x_tail_tip:
        return 0.0
    if x < 0.0:
        return p.r_f * math.sqrt(max(0.0, 1.0 - (x / p.l_nose) ** 2))
    if x <= p.x_body_end:
        return p.r_f
    if x <= p.x_tail_start:
        t = (x - p.x_body_end) / p.l_tail
        return p.r_f * (1.0 - t) + p.r_tail * t
    dx = x - p.x_tail_start
    return math.sqrt(max(0.0, p.r_tail ** 2 - dx ** 2))


def profile_x(p: FuselageParams, n: int = 80) -> np.ndarray:
    """Meridian x sample stations for the profile spline.

    Cosine-clustered over the full body so points bunch at the nose and tail
    tips (where the curvature R->0 is highest), plus the three interior region
    boundaries as explicit knots so the spline honors the shoulders.
    """
    base = 0.5 * (1.0 - np.cos(np.linspace(0.0, math.pi, n)))
    xs = p.x_nose_tip + base * (p.x_tail_tip - p.x_nose_tip)
    knots = np.array([0.0, p.x_body_end, p.x_tail_start])
    return np.unique(np.concatenate([xs, knots]))


def profile_points(p: FuselageParams, n: int = 80) -> Tuple[np.ndarray, np.ndarray]:
    """(x, R) meridian samples; endpoints sit exactly on the axis (R = 0)."""
    xs = profile_x(p, n)
    rs = np.array([radius_at(p, float(x)) for x in xs])
    rs[0] = 0.0
    rs[-1] = 0.0
    return xs, rs


def add_fuselage_solid(occ, p: FuselageParams, n_profile: int = 80) -> List[Tuple[int, int]]:
    """Build the fuselage as a full body of revolution and return its
    volume (dim, tag) entities. `occ` is a live `gmsh.model.occ` handle.

    The meridian face lies in the z = 0 plane (points (x, R, 0), R >= 0),
    bounded by the profile spline and the x-axis segment closing it; a full
    2*pi revolve about the x-axis gives the solid. The revolve's leftover base
    face is dropped (it bounds no volume), mirroring wing3d._add_round_tip_cap.
    """
    xs, rs = profile_points(p, n_profile)
    pts = [occ.addPoint(float(x), float(r), 0.0) for x, r in zip(xs, rs)]
    spline = occ.addSpline(pts)                 # nose tip -> tail tip (R >= 0)
    axis = occ.addLine(pts[-1], pts[0])         # closes the loop along y = 0
    face = occ.addPlaneSurface([occ.addCurveLoop([spline, axis])])
    rev = occ.revolve([(2, face)], 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0 * math.pi)
    occ.remove([(2, face)], recursive=False)
    return [dt for dt in rev if dt[0] == 3]
