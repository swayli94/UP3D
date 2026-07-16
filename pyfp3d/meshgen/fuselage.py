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

PROPORTIONS -- user re-spec 2026-07-16. The body delivered 2026-07-13 was
2.13 long with its nose tip only 0.30 (= 0.37 C_ROOT) ahead of the wing root
LE, so the nose's own displacement flow sat on top of the wing -- not a
defensible subsonic wing-body for B9. The body is now sized and POSITIONED by
two rules instead of by absolute station numbers:

  * ``length = LENGTH_IN_ROOT_CHORDS * C_ROOT`` = 5 root chords = 4.030, and
    the body's MID-LENGTH station is ``x_center``, whose default is the wing
    root-chord midpoint -- i.e. the wing sits in the MIDDLE of the body, with
    2 root chords of body ahead of the root LE and 2 behind the root TE.
  * the nose half-ellipsoid's x semi-axis is ``l_nose = NOSE_DIAMETERS``
    body diameters = 0.60, i.e. fineness 4:1 against the r_f = 0.15 radius
    (was 0.30 = 1 diameter, a near-hemisphere).

Consequence worth knowing before reading a B9 result: the nose ellipsoid's
tip curvature radius is r_f^2 / l_nose = 37 mm, half of the old body's 75 mm,
so the nose tip needs its own mesh refinement -- see wingbody.py's tip balls.

Axis convention (shared with wing3d.py / mesh.wake_cut, fixed project-wide):
chord along +x, lift along +y, span along +z. The fuselage is a body of
revolution about the +x axis, centered on the aircraft centerline
(y, z) = (0, 0). The aircraft symmetry plane is z = 0 (the wing root plane),
so the half model keeps z >= 0 and the fuselage contributes its z >= 0 half
skin to the fluid boundary, meeting the symmetry plane along its z = 0
meridian rims.

The defaults place the wing (root chord x in [0, 0.806], B_SEMI = 1.1963) so
that the wing root is buried in the constant-radius cylinder section and
emerges at the wing-fuselage junction near z = r_f. The wing TE is at y = 0
(the ONERA D section has zero TE thickness), so the TE line meets the fuselage
skin at z = r_f exactly -- this is where the level-set wake polyline and the
innermost Kutta station start (see wingbody.py).

Gmsh is imported lazily so the solver test suite does not depend on it; this
module's geometry helpers (FuselageParams, radius_at, profile_x) are pure
numpy and safe to import anywhere.
"""

import math
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np

from pyfp3d.meshgen.wing3d import C_ROOT

#: Body radius. r_f ~ 0.19 C_ROOT puts the wing-fuselage junction at z = r_f.
R_F = 0.15

#: Nose length in body DIAMETERS: the nose half-ellipsoid's x semi-axis is
#: NOSE_DIAMETERS * (2 r_f). 2.0 => fineness 4:1 (user spec 2026-07-16).
NOSE_DIAMETERS = 2.0

#: Total body length in wing ROOT CHORDS (user spec 2026-07-16).
LENGTH_IN_ROOT_CHORDS = 5.0


@dataclass(frozen=True)
class FuselageParams:
    """Simplified fuselage: rounded nose (half ellipsoid), constant-radius
    body (cylinder), tapering afterbody (cone), rounded tail (sphere cap).

    All lengths in meters, on the wing3d planform scale (C_ROOT = 0.806,
    B_SEMI = 1.1963). The body is placed by its mid-length (`x_center`) and
    sized by its total `length`, NOT by absolute shoulder stations -- that is
    what makes "the wing sits in the middle of a 5-root-chord body" a rule the
    defaults enforce rather than an arithmetic coincidence.
    """
    r_f: float = R_F                                    # 0.15
    length: float = LENGTH_IN_ROOT_CHORDS * C_ROOT      # 4.0295 = 5 C_ROOT
    #: Body mid-length station. The default is the wing root-chord midpoint,
    #: so the wing is centered on the body (user spec 2026-07-16).
    x_center: float = 0.5 * C_ROOT                      # 0.4030
    l_nose: float = NOSE_DIAMETERS * 2.0 * R_F          # 0.60 = 2 diameters
    l_tail: float = 0.60       # afterbody (cone) length
    r_tail: float = 0.03       # radius at the start of the tail sphere cap

    def __post_init__(self) -> None:
        for name in ("r_f", "length", "l_nose", "l_tail", "r_tail"):
            if getattr(self, name) <= 0.0:
                raise ValueError(f"FuselageParams.{name} must be > 0")
        if self.r_tail >= self.r_f:
            raise ValueError("r_tail must be smaller than r_f (afterbody tapers)")
        if self.x_body_end <= self.x_nose_end:
            raise ValueError(
                "nose + afterbody leave no constant-radius section: "
                f"length={self.length} but l_nose + l_tail + r_tail = "
                f"{self.l_nose + self.l_tail + self.r_tail}"
            )

    @property
    def x_nose_tip(self) -> float:
        return self.x_center - 0.5 * self.length

    @property
    def x_nose_end(self) -> float:
        """Nose shoulder: where the ellipsoid reaches r_f and the cylinder
        begins. The wing root chord must lie aft of this (junction_z = r_f)."""
        return self.x_nose_tip + self.l_nose

    @property
    def x_tail_tip(self) -> float:
        return self.x_center + 0.5 * self.length

    @property
    def x_tail_start(self) -> float:
        """Where the cone ends and the tail sphere cap begins."""
        return self.x_tail_tip - self.r_tail

    @property
    def x_body_end(self) -> float:
        """End of the constant-radius section (cylinder -> cone)."""
        return self.x_tail_start - self.l_tail


def radius_at(p: FuselageParams, x: float) -> float:
    """Meridian radius R(x) of the (pre-spline) piecewise profile.

    - nose  [x_nose_tip, x_nose_end]   : half ellipsoid, x semi-axis l_nose,
                                         radial semi-axis r_f
    - body  [x_nose_end, x_body_end]   : cylinder, R = r_f
    - tail  [x_body_end, x_tail_start] : cone, r_f -> r_tail (linear)
    - cap   [x_tail_start, x_tail_tip] : sphere, R = sqrt(r_tail^2 - dx^2)
    R = 0 outside the body (on the axis at both tips).
    """
    if x <= p.x_nose_tip or x >= p.x_tail_tip:
        return 0.0
    if x < p.x_nose_end:
        # t = 1 at the nose tip, 0 at the shoulder
        t = (p.x_nose_end - x) / p.l_nose
        return p.r_f * math.sqrt(max(0.0, 1.0 - t * t))
    if x <= p.x_body_end:
        return p.r_f
    if x <= p.x_tail_start:
        t = (x - p.x_body_end) / p.l_tail
        return p.r_f * (1.0 - t) + p.r_tail * t
    dx = x - p.x_tail_start
    return math.sqrt(max(0.0, p.r_tail ** 2 - dx ** 2))


def profile_x(p: FuselageParams, n: int = 120) -> np.ndarray:
    """Meridian x sample stations for the profile spline.

    Cosine-clustered over the full body so points bunch at the nose and tail
    tips (where the curvature R->0 is highest), plus the three interior region
    boundaries as explicit knots so the spline honors the shoulders.
    """
    base = 0.5 * (1.0 - np.cos(np.linspace(0.0, math.pi, n)))
    xs = p.x_nose_tip + base * (p.x_tail_tip - p.x_nose_tip)
    knots = np.array([p.x_nose_end, p.x_body_end, p.x_tail_start])
    return np.unique(np.concatenate([xs, knots]))


def profile_points(p: FuselageParams, n: int = 120) -> Tuple[np.ndarray, np.ndarray]:
    """(x, R) meridian samples; endpoints sit exactly on the axis (R = 0)."""
    xs = profile_x(p, n)
    rs = np.array([radius_at(p, float(x)) for x in xs])
    rs[0] = 0.0
    rs[-1] = 0.0
    return xs, rs


def add_fuselage_solid(occ, p: FuselageParams, n_profile: int = 120) -> List[Tuple[int, int]]:
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
