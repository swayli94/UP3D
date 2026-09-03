'''
C-type structured grid generator for CFL3D, built on gmsh.

Generates a single-block C-grid around a sharp-trailing-edge airfoil and writes
the matching ``cfl3d.xyz`` (PLOT3D, stream binary, double precision) and
``cfl3d.inp`` files.

Topology (identical to the one produced by the reference ``cgrid`` tool, so the
resulting grids are drop-in replacements)::

    k = nk  +---------------------------------------------+   far field (1000)
            |                                             |
    k = 1   +------+=============+------+                     wall + wake cut
         j=1    j=jte1        j=jte2   j=nj

    j = 1      ... far field, downstream, LOWER side
    j = 1..jte1        lower wake cut          (1-1 blocking)
    j = jte1..jte2     airfoil surface, lower TE -> LE -> upper TE  (bc 2004)
    j = jte2..nj       upper wake cut           (1-1 blocking)
    k = 1              wall / wake cut
    k = nk             far field
    i = 1,2            the two symmetry planes of the 2-D case (z = 0 and 1)

How it is built
---------------

gmsh does the meshing.  The block is cut into a stack of zones in k; each zone
is a gmsh transfinite surface whose four sides are chains of straight two-node
segments, so the point distribution along every side is exactly the one
computed here rather than something gmsh invents.  gmsh fills each zone by
transfinite interpolation and the quad meshes are walked back into one (nj, nk)
lattice.

The zone boundaries matter.  A single transfinite zone spanning wall to far
field draws its k-lines as straight lines from the wall to the far field, which
leaves the wall at up to 45 degrees off the normal - useless for a RANS grid.
Instead each zone boundary is a *collar*: the previous one offset outwards
along its own normals.  A normal offset is orthogonal to the curve by
construction, so inside a zone the transfinite lines stay close to the wall
normal, and the wall itself comes out at 90 degrees to within a fraction of a
degree.

The collars are advanced in small sub-steps with the normals recomputed each
time, as in hyperbolic grid generation, and the step is halved whenever the
offset would fold.  The trailing edge needs extra care: the wake cut leaves the
surface at the trailing-edge angle, so the k = 1 line has a real slope
discontinuity there which a plain offset keeps re-sharpening until the front
tangles.  A small blend towards the neighbour average, applied only in a window
around the two trailing-edge points and ramped in with distance from the wall,
rounds it off.  It has to be local: the leading edge turns far more sharply
than the trailing edge, so any curvature-triggered smoothing hits the leading
edge first and destroys the orthogonality this construction exists to provide.

Beyond the collars a single transfinite zone reaches the far field, whose points
are spread uniformly - the wall clustering carried by the collars is meaningless
eighty chords out and handing it to the far field makes that zone's lines cross.

On the RAE2822 this produces a grid whose wall spacing is held to the requested
y+, whose k-lines leave the wall within 6 degrees of normal, and whose cell
volume ratio and interior angles match the reference elliptic generator; CFL3D
gives Cl and Cd within about 1% of the reference grid's answer.

A Winslow/Thomas-Middlecoff elliptic pass over the outer zones is available via
``n_smooth`` but is off by default - on these aspect ratios it degrades the
volume ratios more than it improves anything.

Usage::

    python cgrid_gmsh.py RAE2822.dat --out ./case --mach 0.725 --alpha 2.10

or from python::

    from cgrid_gmsh import CGridGmsh, read_airfoil

    x, y = read_airfoil('RAE2822.dat')
    cg   = CGridGmsh(n_foil=101, n_wake=41, n_grow=97, y_plus=0.9, r_far=80.0)
    cg.generate(x, y)
    cg.write('./case', mach=0.725, alpha=2.10, re_mil=10.5)
'''

from __future__ import annotations

import argparse
import os
import struct

import numpy as np
from scipy.optimize import brentq

try:
    import gmsh
except ImportError as exc:                                  # pragma: no cover
    raise ImportError('gmsh is required: pip install gmsh') from exc


# ---------------------------------------------------------------------------
#  airfoil input
# ---------------------------------------------------------------------------

def read_airfoil(fname):
    '''
    Read an airfoil coordinate file.

    Accepts the usual Selig layout (a title line followed by ``x y`` pairs that
    run from the trailing edge over the upper surface to the leading edge and
    back over the lower surface).  Any line that does not parse as two floats is
    skipped, so headers and blank lines are harmless.

    Returns the closed contour as two arrays ``x, y``.
    '''
    pts = []
    with open(fname, 'r', errors='replace') as f:
        for line in f:
            parts = line.replace(',', ' ').split()
            if len(parts) < 2:
                continue
            try:
                pts.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue

    if len(pts) < 20:
        raise ValueError(f'{fname}: only {len(pts)} coordinate pairs found')

    pts = np.asarray(pts)
    return pts[:, 0].copy(), pts[:, 1].copy()


def normalize_airfoil(x, y, tail_tol=1e-6):
    '''
    Put the airfoil into the canonical frame: leading edge at (0,0), trailing
    edge at (1,0), chord length 1.

    The trailing edge must be sharp - this generator builds a single-block
    C-grid, which has no room for a blunt base.  Returns ``xu, yu, xl, yl``,
    each running from the leading edge to the trailing edge.
    '''
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # the contour runs TE -> one surface -> LE -> other surface -> TE, so the
    # first and last points are the two trailing-edge points (identical when
    # the file closes the contour explicitly)
    te = np.array([0.5 * (x[0] + x[-1]), 0.5 * (y[0] + y[-1])])
    ile = int(np.argmax((x - te[0]) ** 2 + (y - te[1]) ** 2))
    le = np.array([x[ile], y[ile]])

    chord_vec = te - le
    chord = float(np.hypot(*chord_vec))
    cos_t, sin_t = chord_vec / chord

    xs = ((x - le[0]) * cos_t + (y - le[1]) * sin_t) / chord
    ys = (-(x - le[0]) * sin_t + (y - le[1]) * cos_t) / chord

    upper = slice(ile, None, -1)        # LE -> TE, following the original order
    lower = slice(ile, None, 1)

    xu, yu = xs[upper].copy(), ys[upper].copy()
    xl, yl = xs[lower].copy(), ys[lower].copy()

    # the file may run lower-surface-first; the upper side is the one with the
    # larger mean ordinate
    if np.mean(yu) < np.mean(yl):
        xu, yu, xl, yl = xl, yl, xu, yu

    # close the contour onto the exact trailing edge
    for arr in (xu, xl):
        arr[-1] = 1.0
    tail = abs(yu[-1] - yl[-1])
    if tail > tail_tol:
        raise ValueError(
            f'trailing edge thickness is {tail:.3e} chord - this generator '
            'only handles sharp trailing edges')
    yu[-1] = yl[-1] = 0.0

    return xu, yu, xl, yl


# ---------------------------------------------------------------------------
#  point distributions
# ---------------------------------------------------------------------------

def vinokur(n, s0, s1):
    '''
    Two-sided Vinokur stretching on [0, 1].

    ``n`` points, spacing ``s0`` at the first interval and ``s1`` at the last
    (both as a fraction of the total length).  Returns the monotone array of
    n normalized coordinates.
    '''
    if n < 3:
        return np.linspace(0.0, 1.0, n)

    a = np.sqrt(s1 / s0)
    b = 1.0 / ((n - 1) * np.sqrt(s0 * s1))
    xi = np.linspace(0.0, 1.0, n)

    if abs(b - 1.0) < 1e-8:
        u = xi
    elif b > 1.0:
        d = brentq(lambda t: np.sinh(t) / t - b, 1e-8, 50.0)
        u = 0.5 * (1.0 + np.tanh(d * (xi - 0.5)) / np.tanh(0.5 * d))
    else:
        d = brentq(lambda t: np.sin(t) / t - b, 1e-8, np.pi - 1e-8)
        u = 0.5 * (1.0 + np.tan(d * (xi - 0.5)) / np.tan(0.5 * d))

    t = u / (a + (1.0 - a) * u)
    t[0], t[-1] = 0.0, 1.0
    return t


def geometric(n, h1, length):
    '''
    One-sided geometric distribution on [0, length]: ``n`` points, first
    interval ``h1``.  Returns the coordinates and the growth ratio.
    '''
    ncell = n - 1
    if h1 * ncell >= length:
        raise ValueError('first spacing too large for the requested extent')

    ratio = brentq(lambda r: h1 * (r ** ncell - 1.0) / (r - 1.0) - length,
                   1.0 + 1e-12, 5.0)
    h = h1 * ratio ** np.arange(ncell)
    s = np.concatenate(([0.0], np.cumsum(h)))
    s *= length / s[-1]
    return s, ratio


def resample_curve(x, y, t_norm):
    '''
    Resample a polyline at the normalized arc-length positions ``t_norm``.
    Uses a cubic spline in arc length so that leading-edge curvature survives
    the redistribution.
    '''
    from scipy.interpolate import CubicSpline

    s = np.concatenate(([0.0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))))
    s /= s[-1]
    cs_x = CubicSpline(s, x)
    cs_y = CubicSpline(s, y)
    return cs_x(t_norm), cs_y(t_norm)


def y_plus_to_height(y_plus, re_chord):
    '''
    First cell height, in chords, for a target y+ at chord Reynolds number
    ``re_chord``.

    Flat-plate 1/7-power estimate, written in the only variable it actually
    depends on once everything is non-dimensionalised by the chord::

        cf = 0.0576 Re_c^-0.2
        y+ = (y/c) Re_c sqrt(cf/2)      ->      y/c = y+ / (Re_c sqrt(cf/2))

    Sizing this from a unit Reynolds number that has nothing to do with the
    case - which an earlier version of this file did - silently puts the wall
    spacing wherever it likes; here y+ comes out within about 10% of the
    target.
    '''
    cf = 0.0576 * re_chord ** (-0.2)
    return y_plus / (re_chord * np.sqrt(0.5 * cf))


# ---------------------------------------------------------------------------
#  the generator
# ---------------------------------------------------------------------------

class CGridGmsh:
    '''
    Single-block C-grid generator.

    Parameters
    ----------
    n_foil : points on the airfoil upper (or lower) surface, LE and TE included.
             Forced odd.  Total airfoil points ``nf = 2*n_foil - 1``.
    n_wake : points along one wake cut, TE and far field included.  Forced odd.
    n_grow : points from the wall to the far field.
    n_coarse : how many coarser levels CFL3D will build (mesh sequencing plus
             multigrid).  Every block edge is rounded down to
             ``2**n_coarse * m + 1`` so the block can actually be halved that
             many times; the default 2 matches ncg = 2 / mseq = 3 in the
             cfl3d.inp written here.
    y_plus : target y+ used to size the first cell height.
    h1     : first cell height in chords.  Overrides ``y_plus`` when given -
             use it for Euler grids, where wall clustering is pointless and a
             1e-3 chord first layer is plenty.
    r_far  : far-field radius, in chords.
    x_wake : downstream extent of the wake cut, in chords.
    re_chord : chord Reynolds number of the case, used to size the first
             cell from ``y_plus``.  Must match the run, or y+ will not.
    ds_le, ds_te : surface spacing at the leading and trailing edge, in chords.
             ``None`` picks ``le_factor`` and 1.05 times the mean surface
             spacing.  ds_te also sets the first wake interval, so the k = 1
             line has no spacing jump across the trailing edge.
    le_factor : leading-edge arc spacing as a fraction of the mean.  The nose
             is round, so a spacing that looks fine along the arc is coarse in
             x; this is the knob that resolves the suction peak.
    zone_layers : k layers per transfinite zone.  Fewer means more collars and
             lines that hug the wall normal further out, at the cost of more
             gmsh calls.
    collar_max : distance in chords out to which collars are marched.  Past
             roughly 0.1-0.4 chord (airfoil dependent) the offset front folds;
             the marcher reports that rather than producing a tangled grid.
    te_relax, te_start : strength of the trailing-edge rounding applied to the
             marching front, and the distance over which it is ramped in.
    relax_start : distance over which the collars' point spacing is relaxed
             towards uniform.  The default effectively disables it - relaxing
             inside the collar region slides wake points by whole chords and
             shears the cells.
    far_uniform : 1 spreads the far-field points uniformly (recommended),
             0 inherits the outermost collar's spacing.
    n_smooth : Winslow/Thomas-Middlecoff sweeps over the outer zones.  Off by
             default; see the module docstring.
    '''

    def __init__(self, n_foil=101, n_wake=41, n_grow=97, n_coarse=2,
                 y_plus=1.0, re_chord=6.5e6, h1=None,
                 r_far=80.0, x_wake=None,
                 ds_le=None, ds_te=None, le_factor=0.06,
                 zone_layers=4, collar_max=None,
                 relax_start=1e9, far_uniform=1.0, te_relax=0.3,
                 wake_ratio_max=1.5,
                 te_start=0.01, te_width=0.06, n_smooth=0):

        # CFL3D coarsens the block by a factor of two per multigrid / mesh
        # sequencing level, so every edge must have 2^n_coarse * m + 1 points.
        # Halving stops at a boundary-condition segment end as well, so the
        # trailing-edge indices have to land on the coarse grid too - otherwise
        # the coarse levels get a wall segment that does not line up.
        f = 2 ** int(n_coarse)
        self.mg_factor = f

        self.n_foil = int(np.floor((n_foil - 1) / (f // 2)) * (f // 2) + 1)
        self.n_wake = int(np.floor((n_wake - 1) / f) * f + 1)
        self.n_grow = int(np.floor((n_grow - 1) / f) * f + 1)

        self.nf = 2 * self.n_foil - 1
        self.nj = 2 * (self.n_foil + self.n_wake) - 3
        self.nk = self.n_grow
        self.jte1 = self.n_wake                     # 1-based j of the lower TE
        self.jte2 = self.nj - self.n_wake + 1       # 1-based j of the upper TE

        for name, n in (('jdim', self.nj), ('kdim', self.nk),
                        ('lower TE index', self.jte1),
                        ('upper TE index', self.jte2)):
            if (n - 1) % f:
                raise ValueError(
                    f'{name} = {n} is not {f}*m+1; with {n_coarse} coarser '
                    f'levels CFL3D cannot halve this block. Adjust n_foil / '
                    f'n_wake / n_grow, or lower n_coarse.')

        self.y_plus = y_plus
        self.h1 = h1
        self.r_far = r_far
        self.x_wake = r_far if x_wake is None else x_wake
        self.re_chord = re_chord
        self.ds_le, self.ds_te = ds_le, ds_te
        self.le_factor = le_factor
        self.zone_layers = int(zone_layers)
        self.collar_max = 0.03 if collar_max is None else collar_max
        self.relax_start = relax_start
        self.far_uniform = far_uniform
        self.wake_ratio_max = wake_ratio_max
        self.te_relax = te_relax
        self.te_start = te_start
        self.te_width = te_width
        self.n_smooth = n_smooth

        self.xy = None          # (nj, nk, 2) once generated
        self.info = {}

    # -- boundaries ---------------------------------------------------------

    def _inner_boundary(self, xu, yu, xl, yl):
        '''
        Build the k = 1 line: lower wake (far field -> TE), lower surface
        (TE -> LE), upper surface (LE -> TE), upper wake (TE -> far field).
        '''
        arc = np.sum(np.hypot(np.diff(xu), np.diff(yu)))
        ds_le = self.ds_le if self.ds_le else self.le_factor * arc / (self.n_foil - 1)
        ds_te = self.ds_te if self.ds_te else 1.05 * arc / (self.n_foil - 1)
        self.info['ds_le'], self.info['ds_te'] = ds_le, ds_te

        t = vinokur(self.n_foil, ds_le / arc, ds_te / arc)
        xs_u, ys_u = resample_curve(xu, yu, t)          # LE -> TE
        xs_l, ys_l = resample_curve(xl, yl, t)          # LE -> TE

        # wake cut: geometric growth away from the TE, starting at the surface
        # spacing there so the k = 1 line has no spacing jump
        sw, ratio_w = geometric(self.n_wake, ds_te, self.x_wake - 1.0)
        self.info['wake_ratio'] = ratio_w
        if ratio_w > self.wake_ratio_max:
            print(f'  WARNING: wake growth ratio {ratio_w:.3f} exceeds '
                  f'{self.wake_ratio_max:.2f}; raise n_wake or lower x_wake')
        xw = 1.0 + sw                                   # TE -> far field

        x = np.concatenate([xw[::-1][:-1], xs_l[::-1][:-1], xs_u, xw[1:]])
        y = np.concatenate([np.zeros(self.n_wake - 1), ys_l[::-1][:-1],
                            ys_u, np.zeros(self.n_wake - 1)])

        assert len(x) == self.nj, (len(x), self.nj)
        return np.stack([x, y], axis=1)

    def _outer_boundary(self, ref):
        '''
        Build the k = nk line: two horizontal lines at y = -/+ ``r_far`` joined
        upstream by a semicircle of the same radius.

        Points are placed by blending the normalized arc length of ``ref`` -
        the outermost collar curve - towards a uniform spread.  The collar
        still carries the strong leading-edge and wake clustering of the wall;
        handing that straight to a boundary eighty chords away makes the
        transfinite lines of the last zone cross.  ``far_uniform = 1`` (the
        default) drops the clustering entirely, which is what a far field
        wants anyway.
        '''
        r, xw = self.r_far, self.x_wake

        n_arc = 4000
        th = np.linspace(-0.5 * np.pi, -1.5 * np.pi, n_arc)
        path = np.vstack([
            np.stack([np.linspace(xw, 0.0, 400), np.full(400, -r)], axis=1),
            np.stack([r * np.cos(th), r * np.sin(th)], axis=1),
            np.stack([np.linspace(0.0, xw, 400), np.full(400, r)], axis=1)])

        s = np.concatenate(([0.0], np.cumsum(
            np.hypot(np.diff(path[:, 0]), np.diff(path[:, 1])))))
        s /= s[-1]

        t = np.concatenate(([0.0], np.cumsum(
            np.hypot(np.diff(ref[:, 0]), np.diff(ref[:, 1])))))
        t /= t[-1]
        w = self.far_uniform
        t = (1.0 - w) * t + w * np.linspace(0.0, 1.0, len(ref))

        return np.stack([np.interp(t, s, path[:, 0]),
                         np.interp(t, s, path[:, 1])], axis=1)

    # -- gmsh ---------------------------------------------------------------

    @staticmethod
    def _chain(geo, pts, first=None, last=None):
        '''
        Create gmsh points and two-node lines through ``pts``; reuse the given
        end point tags where the chain shares a corner with another chain.
        '''
        tags = []
        for i, (px, py) in enumerate(pts):
            if i == 0 and first is not None:
                tags.append(first)
            elif i == len(pts) - 1 and last is not None:
                tags.append(last)
            else:
                tags.append(geo.addPoint(float(px), float(py), 0.0))
        lines = [geo.addLine(tags[i], tags[i + 1]) for i in range(len(tags) - 1)]
        for t in lines:
            geo.mesh.setTransfiniteCurve(t, 2)
        return tags, lines

    @staticmethod
    def _normals(curve):
        '''Unit normals of a polyline, pointing to the left of the direction of
        travel.  With j running clockwise this is the outward direction.'''
        t = np.empty_like(curve)
        t[1:-1] = curve[2:] - curve[:-2]
        t[0] = curve[1] - curve[0]
        t[-1] = curve[-1] - curve[-2]
        t /= np.linalg.norm(t, axis=1)[:, None]
        return np.stack([-t[:, 1], t[:, 0]], axis=1)

    @staticmethod
    def _strip_area(a, b):
        '''Smallest signed cell area of the strip between two layer curves.'''
        u = b[1:] - a[:-1]
        v = b[:-1] - a[1:]
        return float(np.min(0.5 * (u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])))

    def _offset(self, curve, distance, max_smooth=60):
        '''
        Offset a layer curve outward by ``distance`` along its own normals.

        A plain offset is exactly orthogonal to the curve, which is what makes
        the transfinite lines inside the zone leave the wall along the normal.
        It folds if it is pushed further than the local radius of curvature, so
        the normal field is low-pass filtered just enough - and no more - to
        keep every cell of the strip positive.  Near the wall no filtering is
        needed at all, so the wall stays exactly orthogonal; the compromise is
        only ever made far away, where it costs nothing.
        '''
        n0 = self._normals(curve)
        n = n0.copy()
        for passes in range(max_smooth + 1):
            if passes:
                n[1:-1] = 0.25 * n[:-2] + 0.5 * n[1:-1] + 0.25 * n[2:]
                n /= np.linalg.norm(n, axis=1)[:, None]
            out = curve + distance * n
            out[0, 0] = out[-1, 0] = self.x_wake     # stay on the outflow line
            if self._strip_area(curve, out) > 0.0:
                return out, passes
        return None, max_smooth

    @staticmethod
    def _relax_spacing(curve, weight, passes=4):
        '''
        Even out the point spacing along a layer curve without moving the
        curve itself: the arc-length parameter of the points is smoothed and
        the curve is resampled there.

        A pure offset inherits the wall clustering, which piles points up in
        concave regions and eventually tangles the front.  ``weight`` ramps
        this correction in with distance, so the first layers - where a
        tangential shift would wreck cells that are 1e-6 chord tall - are left
        untouched.
        '''
        if weight <= 0.0:
            return curve

        u = np.concatenate(([0.0], np.cumsum(
            np.hypot(np.diff(curve[:, 0]), np.diff(curve[:, 1])))))
        u /= u[-1]

        v = u.copy()
        for _ in range(passes):
            v[1:-1] = 0.5 * v[1:-1] + 0.25 * (v[:-2] + v[2:])
        v = u + weight * (v - u)
        v = np.maximum.accumulate(v)

        return np.stack([np.interp(v, u, curve[:, 0]),
                         np.interp(v, u, curve[:, 1])], axis=1)

    def _te_taper(self, inner):
        '''
        Blend weight used to round the trailing-edge corner of a marching
        front.  Only the two trailing-edge kinks get smoothed: everywhere else
        the front is left exactly where the normal offset put it.

        A kink-detector keyed on the turning angle cannot be used instead - the
        leading edge turns far more sharply than the trailing edge does, and
        smoothing it destroys the very orthogonality this construction exists
        to provide.  The trailing edge, on the other hand, is a genuine slope
        discontinuity in the k = 1 line (the wake cut leaves the surface at the
        trailing-edge angle), and its position is known exactly.

        The window is measured in arc length, not in points: a fixed point
        count shrinks physically as the surface is refined, which sharpens the
        rounded corner again and shows up as a cell-volume jump at the collar
        on the finer grids.
        '''
        s = np.concatenate(([0.0], np.cumsum(
            np.hypot(np.diff(inner[:, 0]), np.diff(inner[:, 1])))))

        w = np.zeros(self.nj)
        for jte in (self.jte1 - 1, self.jte2 - 1):
            d = np.abs(s - s[jte]) / self.te_width
            w = np.maximum(w, np.where(d < 1.0,
                                       0.5 * (1.0 + np.cos(np.pi * d)), 0.0))
        w[0] = w[-1] = 0.0
        return w

    @staticmethod
    def _blend_towards_average(curve, w):
        '''Move each point a fraction ``w`` of the way to its neighbours.'''
        if not w.any():
            return curve
        out = curve.copy()
        avg = 0.5 * (curve[:-2] + curve[2:])
        out[1:-1] = ((1.0 - w[1:-1, None]) * curve[1:-1]
                     + w[1:-1, None] * avg)
        return out

    def _collar_ladder(self, inner, s_k, edges):
        '''
        Build the stack of layer curves that separate the meshing zones.  Each
        one is an offset of the previous, so the offset distance stays small
        compared with the local curvature and the normals are re-evaluated as
        the front moves out - the same idea as hyperbolic marching, but only
        one curve per zone is needed because gmsh fills the rest.
        '''
        curves = [inner]
        self._march_steps = 0
        self._te_w = self._te_taper(inner)
        for m in range(1, len(edges) - 1):
            curves.append(self._advance(curves[-1], s_k[edges[m - 1]],
                                        s_k[edges[m]]))
        self.info['march_steps'] = self._march_steps
        return curves

    def _advance(self, curve, s0, s1, step_frac=0.3, step_floor=0.02):
        '''
        Move the front from distance ``s0`` to ``s1``.

        The step is kept small compared with the distance already covered (the
        front's smallest features scale with it), and is halved whenever the
        offset would fold.  Recomputing the normals at every sub-step is what
        lets the front round off the trailing-edge kink instead of tangling on
        it, so the far-field zones stay valid.
        '''
        cur, s = curve, s0
        while s < s1 - 1e-14:
            step = min(s1 - s, max(step_frac * s, step_floor))
            while True:
                out, _ = self._offset(cur, step)
                if out is not None:
                    break
                step *= 0.5
                if step < (s1 - s0) * 1e-6:
                    raise RuntimeError(
                        f'marching front folded at {s:.4g} chords')
            s += step
            cur = self._relax_spacing(out, min(1.0, s / self.relax_start))
            ramp = min(1.0, s / self.te_start)
            cur = self._blend_towards_average(
                cur, self._te_w * self.te_relax * ramp)
            self._march_steps += 1
            if self._march_steps > 10000:
                raise RuntimeError('marching did not terminate')
        return cur

    def _mesh_zone(self, inner, outer, s_k):
        '''
        Mesh one zone: a transfinite surface between two layer curves.  All
        four sides are chains of two-node segments, so the point distributions
        are exactly the ones computed here rather than anything gmsh invents.
        '''
        nk = len(s_k)
        f = (s_k - s_k[0]) / (s_k[-1] - s_k[0])
        left = inner[0] + np.outer(f, outer[0] - inner[0])
        right = inner[-1] + np.outer(f, outer[-1] - inner[-1])

        gmsh.initialize()
        gmsh.option.setNumber('General.Terminal', 0)
        # the two sides of the wake cut are geometrically coincident; without
        # this gmsh would merge them and the C topology would collapse
        gmsh.option.setNumber('Geometry.AutoCoherence', 0)
        gmsh.model.add('cgrid_zone')
        geo = gmsh.model.geo

        p_in, l_in = self._chain(geo, inner)
        p_out, l_out = self._chain(geo, outer)
        _, l_lf = self._chain(geo, left, first=p_in[0], last=p_out[0])
        _, l_rt = self._chain(geo, right, first=p_in[-1], last=p_out[-1])

        loop = geo.addCurveLoop(l_lf + l_out + [-t for t in reversed(l_rt)]
                                + [-t for t in reversed(l_in)])
        surf = geo.addPlaneSurface([loop])
        geo.synchronize()

        gmsh.model.mesh.setTransfiniteSurface(
            surf, 'Left', [p_in[0], p_out[0], p_out[-1], p_in[-1]])
        gmsh.model.mesh.setRecombine(2, surf)
        gmsh.model.mesh.generate(2)

        xy = self._extract_lattice(p_in, self.nj, nk)
        gmsh.finalize()
        return xy

    def _extract_lattice(self, p_inner, nj, nk):
        '''
        Rebuild the (nj, nk, 2) lattice from the quad mesh.

        Every point of the k = 1 row was created explicitly, so that row is
        known exactly; the remaining rows follow by stepping through the quad
        that sits on top of each edge of the current row.
        '''
        tags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.asarray(coords).reshape(-1, 3)
        xyz = {int(t): coords[i, :2] for i, t in enumerate(tags)}

        etypes, _, enodes = gmsh.model.mesh.getElements(2)
        blocks = [np.asarray(en, dtype=np.int64).reshape(-1, 4)
                  for et, en in zip(etypes, enodes) if et == 3]
        quads = np.vstack(blocks) if blocks else None
        if quads is None or len(quads) != (nj - 1) * (nk - 1):
            raise RuntimeError(
                f'expected {(nj - 1) * (nk - 1)} quads, got '
                f'{0 if quads is None else len(quads)}')

        # edge -> [(quad index, position of the edge in the quad)]
        edge_map = {}
        for q, nodes in enumerate(quads):
            for e in range(4):
                a, b = nodes[e], nodes[(e + 1) % 4]
                edge_map.setdefault((min(a, b), max(a, b)), []).append((q, e))

        row = [int(gmsh.model.mesh.getNodes(0, p)[0][0]) for p in p_inner]
        rows = [row]
        used = set()

        for _ in range(nk - 1):
            nxt = [None] * nj
            for j in range(nj - 1):
                a, b = row[j], row[j + 1]
                cand = [qi for qi in edge_map[(min(a, b), max(a, b))]
                        if qi[0] not in used]
                if len(cand) != 1:
                    raise RuntimeError('lattice walk failed - ambiguous quad')
                q, e = cand[0]
                used.add(q)
                nodes = quads[q]
                # the quad is cyclic; the edge opposite (a,b) is (c,d)
                c, d = nodes[(e + 2) % 4], nodes[(e + 3) % 4]
                if nodes[e] == a:                   # quad runs a,b,c,d
                    up_a, up_b = d, c
                else:                               # quad runs b,a,d,c
                    up_a, up_b = c, d
                if nxt[j] is not None and nxt[j] != up_a:
                    raise RuntimeError('lattice walk failed - row mismatch')
                nxt[j], nxt[j + 1] = up_a, up_b
            row = [int(t) for t in nxt]
            rows.append(row)

        xy = np.zeros((nj, nk, 2))
        for k, r in enumerate(rows):
            for j, t in enumerate(r):
                xy[j, k] = xyz[t]
        return xy

    # -- elliptic smoothing -------------------------------------------------

    def _respace_outer(self, xy, k0, ratio):
        '''
        Re-space the k-lines above the collar so the wall-normal growth stays
        continuous.

        Inside the collar region every zone is a true normal offset, so the
        spacing is exactly the requested geometric distribution.  The last zone
        is not: transfinite interpolation lays its points out in proportion to
        the *length of that particular k-line*, and the far field is spread
        uniformly, so that length varies strongly with j.  The result is a
        spacing jump at the interface - measured at 1.47 on the RAE2822 grid,
        against 1.18 everywhere else.

        The fix is one-dimensional and leaves the k-lines exactly where they
        are: walk each line, start from the last collar interval times the
        global growth ratio, and re-solve the ratio so the geometric
        distribution lands on the far field.
        '''
        nj, nk = xy.shape[:2]
        out = xy.copy()
        worst = 1.0

        for j in range(nj):
            seg = xy[j, k0 - 1:]
            step = np.linalg.norm(np.diff(seg, axis=0), axis=1)
            s = np.concatenate(([0.0], np.cumsum(step)))
            length, ncell = s[-1], len(seg) - 1
            if ncell < 2 or length <= 0.0:
                continue

            h_prev = float(np.linalg.norm(xy[j, k0 - 1] - xy[j, k0 - 2]))
            h0 = min(h_prev * ratio, 0.9 * length / ncell)
            try:
                new, r = geometric(ncell + 1, h0, length)
            except ValueError:                  # cannot stretch: leave uniform
                continue
            worst = max(worst, r)

            out[j, k0 - 1:, 0] = np.interp(new, s, seg[:, 0])
            out[j, k0 - 1:, 1] = np.interp(new, s, seg[:, 1])

        self.info['k_ratio_outer'] = worst
        return out

    def _smooth(self, xy, n_iter, omega=1.5, k_freeze=1):
        '''
        Winslow smoothing with Thomas-Middlecoff control functions.

        Solves ``a (x_jj + phi x_j) - 2 b x_jk + c (x_kk + psi x_k) = 0`` by
        SOR, where ``phi`` comes from the point distribution on the k = const
        boundaries (blended in k) and ``psi`` from the j = const boundaries
        (blended in j).  Attaching the control functions to the metric
        coefficients rather than to J^2 keeps the system well scaled on grids
        with the 1e4 near-wall aspect ratios needed for y+ ~ 1.

        The effect is that the boundary point distributions - in particular the
        wall spacing that sets y+ - are held while the interior is pulled
        towards a smooth, near-orthogonal configuration.
        '''
        if n_iter <= 0:
            return xy

        x = xy[:, :, 0].copy()
        y = xy[:, :, 1].copy()
        nj, nk = x.shape

        def tm_source(a, b):
            '''-(a_t a_tt + b_t b_tt) / (a_t^2 + b_t^2), padded to full length'''
            at = 0.5 * (a[2:] - a[:-2])
            bt = 0.5 * (b[2:] - b[:-2])
            att = a[2:] - 2.0 * a[1:-1] + a[:-2]
            btt = b[2:] - 2.0 * b[1:-1] + b[:-2]
            den = at ** 2 + bt ** 2
            s = -(at * att + bt * btt) / np.where(den > 1e-30, den, 1e-30)
            return np.pad(s, 1, mode='edge')

        # phi from the k = const boundaries, linearly blended in k
        phi_wall = tm_source(x[:, 0], y[:, 0])
        phi_far = tm_source(x[:, -1], y[:, -1])
        eta = np.linspace(0.0, 1.0, nk)[None, :]
        phi = phi_wall[:, None] * (1.0 - eta) + phi_far[:, None] * eta

        # psi from the j = const boundaries, linearly blended in j
        psi_lf = tm_source(x[0, :], y[0, :])
        psi_rt = tm_source(x[-1, :], y[-1, :])
        xi = np.linspace(0.0, 1.0, nj)[:, None]
        psi = psi_lf[None, :] * (1.0 - xi) + psi_rt[None, :] * xi

        ph = phi[1:-1, 1:-1]
        ps = psi[1:-1, 1:-1]

        # red/black colouring so the sweeps are Gauss-Seidel (a plain
        # vectorised Jacobi update would diverge for omega > 1); the first
        # k_freeze rows come from the boundary-layer extrusion and are held
        jj, kk = np.meshgrid(np.arange(nj - 2), np.arange(nk - 2),
                             indexing='ij')
        free = kk >= (k_freeze - 1)
        colours = (((jj + kk) % 2 == 0) & free, ((jj + kk) % 2 == 1) & free)

        def residual_update(mask):
            xj = 0.5 * (x[2:, 1:-1] - x[:-2, 1:-1])
            yj = 0.5 * (y[2:, 1:-1] - y[:-2, 1:-1])
            xk = 0.5 * (x[1:-1, 2:] - x[1:-1, :-2])
            yk = 0.5 * (y[1:-1, 2:] - y[1:-1, :-2])

            a = xk ** 2 + yk ** 2                       # multiplies x_jj
            b = xj * xk + yj * yk
            c = xj ** 2 + yj ** 2                       # multiplies x_kk

            xjk = 0.25 * (x[2:, 2:] - x[2:, :-2] - x[:-2, 2:] + x[:-2, :-2])
            yjk = 0.25 * (y[2:, 2:] - y[2:, :-2] - y[:-2, 2:] + y[:-2, :-2])

            den = 2.0 * (a + c)
            den = np.where(den > 1e-300, den, 1e-300)

            rx = (a * (x[2:, 1:-1] + x[:-2, 1:-1]
                       + 0.5 * ph * (x[2:, 1:-1] - x[:-2, 1:-1]))
                  + c * (x[1:-1, 2:] + x[1:-1, :-2]
                         + 0.5 * ps * (x[1:-1, 2:] - x[1:-1, :-2]))
                  - 2.0 * b * xjk) / den
            ry = (a * (y[2:, 1:-1] + y[:-2, 1:-1]
                       + 0.5 * ph * (y[2:, 1:-1] - y[:-2, 1:-1]))
                  + c * (y[1:-1, 2:] + y[1:-1, :-2]
                         + 0.5 * ps * (y[1:-1, 2:] - y[1:-1, :-2]))
                  - 2.0 * b * yjk) / den

            dx = np.where(mask, omega * (rx - x[1:-1, 1:-1]), 0.0)
            dy = np.where(mask, omega * (ry - y[1:-1, 1:-1]), 0.0)
            x[1:-1, 1:-1] += dx
            y[1:-1, 1:-1] += dy
            return max(np.abs(dx).max(), np.abs(dy).max())

        for it in range(n_iter):
            shift = max(residual_update(c) for c in colours)
            if not np.isfinite(shift):
                raise RuntimeError('elliptic smoothing diverged - reduce omega')
            if shift < 1e-12 * self.r_far:
                self.info['smooth_iters'] = it + 1
                break
        else:
            self.info['smooth_iters'] = n_iter
        self.info['smooth_shift'] = float(shift)

        return np.stack([x, y], axis=2)

    # -- driver -------------------------------------------------------------

    def generate(self, x_foil, y_foil, verbose=True):
        '''
        Build the grid from a closed airfoil contour.  Stores the result in
        ``self.xy`` with shape (nj, nk, 2) and returns it.
        '''
        xu, yu, xl, yl = normalize_airfoil(x_foil, y_foil)

        if self.h1 is not None:
            h1 = self.h1                    # explicit spacing (Euler grids)
        else:
            h1 = y_plus_to_height(self.y_plus, self.re_chord)
        s_k, ratio_k = geometric(self.nk, h1, self.r_far)
        self.info['h1'] = h1
        self.info['k_ratio'] = ratio_k

        inner = self._inner_boundary(xu, yu, xl, yl)

        # k is split into zones; every zone is a transfinite surface between
        # two nearly parallel layer curves, so its interior lines stay close to
        # the wall normal instead of fanning straight out to the far field
        edges = [0]
        while True:
            nxt = edges[-1] + self.zone_layers
            if nxt >= self.nk - 1 or s_k[nxt] >= self.collar_max:
                break
            edges.append(nxt)
        edges.append(self.nk - 1)
        self.info['n_zones'] = len(edges) - 1
        self.info['zone_edges'] = edges
        self.info['collar'] = float(s_k[edges[-2]])

        curves = self._collar_ladder(inner, s_k, edges)
        curves.append(self._outer_boundary(curves[-1]))

        blocks = [self._mesh_zone(curves[m], curves[m + 1],
                                  s_k[edges[m]:edges[m + 1] + 1])
                  for m in range(len(edges) - 1)]
        xy = np.concatenate([blocks[0]] + [b[:, 1:] for b in blocks[1:]],
                            axis=1)
        xy = self._respace_outer(xy, edges[-2] + 1, ratio_k)
        xy = self._smooth(xy, self.n_smooth, k_freeze=edges[-2] + 1)

        self.xy = xy
        self.info.update(self.quality())

        if verbose:
            self.report()
        return xy

    # -- diagnostics --------------------------------------------------------

    def quality(self):
        '''Wall spacing, wall orthogonality and cell-area positivity.'''
        xy = self.xy
        j1, j2 = self.jte1 - 1, self.jte2 - 1

        h = np.hypot(xy[j1:j2 + 1, 1, 0] - xy[j1:j2 + 1, 0, 0],
                     xy[j1:j2 + 1, 1, 1] - xy[j1:j2 + 1, 0, 1])

        t = xy[j1 + 1:j2 + 1, 0] - xy[j1:j2, 0]
        t = 0.5 * (np.vstack([t[:1], t]) + np.vstack([t, t[-1:]]))
        n = xy[j1:j2 + 1, 1] - xy[j1:j2 + 1, 0]
        cos = np.sum(t * n, axis=1) / (np.linalg.norm(t, axis=1)
                                       * np.linalg.norm(n, axis=1))
        ang = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))

        a = xy[1:, 1:] - xy[:-1, :-1]
        b = xy[:-1, 1:] - xy[1:, :-1]
        area = 0.5 * (a[:, :, 0] * b[:, :, 1] - a[:, :, 1] * b[:, :, 0])

        rk = area[:, 1:] / area[:, :-1]
        rj = area[1:, :] / area[:-1, :]

        tj = xy[2:, 1:-1] - xy[:-2, 1:-1]
        tk = xy[1:-1, 2:] - xy[1:-1, :-2]
        cell = np.sum(tj * tk, axis=2) / (np.linalg.norm(tj, axis=2)
                                          * np.linalg.norm(tk, axis=2))
        cell = np.degrees(np.arccos(np.clip(cell, -1.0, 1.0)))

        dk = np.hypot(np.diff(xy[:, :, 0], axis=1), np.diff(xy[:, :, 1], axis=1))
        rk_line = dk[:, 1:] / dk[:, :-1]

        dj = np.hypot(np.diff(xy[:, 0, 0]), np.diff(xy[:, 0, 1]))
        i = self.jte1 - 1                       # 0-based lower TE
        dxw = abs(xy[i, 0, 0] - xy[i - 1, 0, 0])
        dxs = abs(xy[i + 1, 0, 0] - xy[i, 0, 0])
        jle = (self.jte1 + self.jte2) // 2 - 1

        return {'h1_min': float(h.min()), 'h1_max': float(h.max()),
                'k_ratio_max': float(np.maximum(rk_line, 1.0 / rk_line).max()),
                'te_match': float(dxw / dxs),
                'dx_le': float(abs(xy[jle + 1, 0, 0] - xy[jle, 0, 0])),
                'ortho_min': float(ang.min()), 'ortho_max': float(ang.max()),
                'vol_ratio': float(np.maximum(rk, 1.0 / rk).max()),
                'vol_ratio_j': float(np.maximum(rj, 1.0 / rj).max()),
                'angle_min': float(cell.min()), 'angle_max': float(cell.max()),
                'area_min': float(area.min()),
                'negative_cells': int(np.sum(area <= 0.0))}

    def report(self):
        i = self.info
        print(f'  grid            : {self.nj} x {self.nk}  '
              f'(i=2, airfoil j={self.jte1}..{self.jte2}), '
              f'all edges {self.mg_factor}m+1')
        if self.h1 is None:
            print(f'  target y+       : {self.y_plus:.2f}  ->  '
                  f'h1 = {i["h1"]:.3e}')
        else:
            print(f'  first layer     : h1 = {i["h1"]:.3e} (given)')
        print(f'  wall spacing    : {i["h1_min"]:.3e} .. {i["h1_max"]:.3e}')
        print(f'  k growth ratio  : {i["k_ratio"]:.4f} (collar), '
              f'{i.get("k_ratio_outer", float("nan")):.4f} (outer), '
              f'{i["k_ratio_max"]:.4f} (worst adjacent)')
        print(f'  wake growth     : {i["wake_ratio"]:.4f}   '
              f'TE spacing match: {i["te_match"]:.4f}')
        print(f'  LE spacing      : ds={i["ds_le"]:.3e}  '
              f'first dx={i["dx_le"]:.3e}')
        print(f'  zones / march   : {i["n_zones"]} transfinite zones, '
              f'collar to {i["collar"]:.3f} chords in {i["march_steps"]} steps')
        print(f'  wall orthogon.  : {i["ortho_min"]:.2f} .. '
              f'{i["ortho_max"]:.2f} deg (90 = orthogonal)')
        print(f'  cell vol ratio  : {i["vol_ratio"]:.3f} (k), '
              f'{i["vol_ratio_j"]:.3f} (j)')
        print(f'  interior angle  : {i["angle_min"]:.2f} .. '
              f'{i["angle_max"]:.2f} deg')
        print(f'  min cell area   : {i["area_min"]:.3e}, '
              f'negative cells: {i["negative_cells"]}')

    # -- output -------------------------------------------------------------

    def write_xyz(self, fname):
        '''
        PLOT3D grid, stream binary, double precision, single block, 2 planes in
        i.  Matches what this build of CFL3D expects (``form='unformatted',
        access='stream'``).
        '''
        nj, nk = self.nj, self.nk
        x = self.xy[:, :, 0]
        y = self.xy[:, :, 1]

        with open(fname, 'wb') as f:
            f.write(struct.pack('i', 1))
            f.write(struct.pack('iii', 2, nj, nk))
            for arr in (x, y):
                buf = np.empty((nk, nj, 2))
                buf[:, :, 0] = arr.T
                buf[:, :, 1] = arr.T
                f.write(buf.astype('<f8').tobytes())
            z = np.empty((nk, nj, 2))
            z[:, :, 0] = 0.0
            z[:, :, 1] = 1.0
            f.write(z.astype('<f8').tobytes())

    def write_inp(self, fname, mach=0.725, alpha=2.10, beta=0.0, re_mil=10.5,
                  tinf=460.0, cfl=1.0, cfl_max=5.0, iflagts=500, ncyc=1000,
                  mseq=3, ivisc=7, nitfo=0, nwrest=5000):
        '''
        Write the matching ``cfl3d.inp``.

        The block layout, the boundary-condition table and the 1-1 blocking of
        the wake cut are derived from the grid dimensions, so the file always
        agrees with the ``.xyz`` written alongside it.

        ``cfl``/``cfl_max``/``iflagts`` drive CFL3D's ramp: the CFL number
        starts at ``cfl`` and is raised to ``cfl*cfl_max`` over ``iflagts``
        cycles.  Starting cold at a high CFL blows this class of case up.

        ``ivisc = 0`` writes an Euler case: no turbulence model and the airfoil
        surface becomes an inviscid (flow-tangency) wall, bc 1005, instead of
        the adiabatic no-slip wall, bc 2004.

        ``nitfo`` runs that many cycles first order on each mesh-sequence
        level, which together with a small ``cfl`` is what gets a high-
        incidence case started.
        '''
        a1, a2 = self.nj, self.nk           # jdim, kdim
        a3 = self.n_wake                    # j index of the lower TE
        a5 = a1 - a3 + 1                    # j index of the upper TE

        L = []
        w = L.append
        w('FILES:            ')
        for name in ('cfl3d.xyz', 'plot3d_grid.xyz', 'plot3d_sol.bin',
                     'cfl3d.out', 'resid.out', 'cfl3d.turres', 'cfl3d.blomax',
                     'cfl3d.2out', 'cfl3d.prt', 'cfl3d.press', 'ovrlp.bin',
                     'patch.bin', 'cfl3d.restart'):
            w(f'{name:<18}')
        w('CFL3D V6 INPUT FILE GENERATED WITH cgrid_gmsh.py')
        w('     XMACH     ALPHA      BETA  REUE,MIL   TINF,DR     IALPH    IHSTRY')
        w(f'{mach:10.4f}{alpha:10.4f}{beta:10.4f}{re_mil:10.4f}{tinf:10.2f}'
          f'{1:10d}{0:10d}')
        w('      SREF      CREF      BREF       XMC       YMC       ZMC')
        w(f'{1.0:10.2f}{1.0:10.2f}{1.0:10.2f}{0.25:10.2f}{0.0:10.2f}{0.0:10.2f}')
        w('        DT     IREST   IFLAGTS      FMAX     IUNST    CFLTAU')
        w(f'{-abs(cfl):10.2f}{0:10d}{iflagts:10d}{cfl_max:10.2f}{0:10d}'
          f'{7.5:10.4f}')
        w('     NGRID   NPLOT3D    NPRINT    NWREST      ICHK       I2D    NTSTEP       ITA')
        w(f'{-1:10d}{1:10d}{-1:10d}{nwrest:10d}{0:10d}{1:10d}{1:10d}{1:10d}')
        w('       NCG       IEM  IADVANCE    IFORCE  IVISC(I)  IVISC(J)  IVISC(K)')
        w(f'{self.mg_factor.bit_length() - 1:10d}{0:10d}{0:10d}{333:10d}'
          f'{ivisc:10d}{ivisc:10d}{ivisc:10d}')
        w('      IDIM      JDIM      KDIM')
        w(f'{2:10d}{a1:10d}{a2:10d}')
        w('    ILAMLO    ILAMHI    JLAMLO    JLAMHI    KLAMLO    ILAMHI')
        w(f'{0:10d}' * 6)
        w('     INEWG    IGRIDC        IS        JS        KS        IE        JE        KE')
        w(f'{0:10d}' * 8)
        w('  IDIAG(I)  IDIAG(J)  IDIAG(K)  IFLIM(I)  IFLIM(J)  IFLIM(K)')
        w(f'{1:10d}{1:10d}{1:10d}{4:10d}{4:10d}{4:10d}')
        w('   IFDS(I)   IFDS(J)   IFDS(K)  RKAP0(I)  RKAP0(J)  RKAP0(K)')
        w(f'{1:10d}{1:10d}{1:10d}{0.3333:10.4f}{0.3333:10.4f}{0.3333:10.4f}')
        w('      GRID     NBCIO   NBCIDIM     NBCJO   NBCJDIM     NBCKO   NBCKDIM    IOVRLP')
        w(f'{1:10d}{1:10d}{1:10d}{1:10d}{1:10d}{3:10d}{1:10d}{0:10d}')
        w('I0:   GRID   SEGMENT    BCTYPE      JSTA      JEND      KSTA      KEND     NDATA')
        w(f'{1:10d}{1:10d}{1001:10d}{1:10d}{a1:10d}{1:10d}{a2:10d}{0:10d}')
        w('IDIM: GRID   SEGMENT    BCTYPE      JSTA      JEND      KSTA      KEND     NDATA')
        w(f'{1:10d}{1:10d}{1001:10d}{1:10d}{a1:10d}{1:10d}{a2:10d}{0:10d}')
        w('J0:   GRID   SEGMENT    BCTYPE      ISTA      IEND      KSTA      KEND     NDATA')
        w(f'{1:10d}{1:10d}{1000:10d}{1:10d}{2:10d}{1:10d}{a2:10d}{0:10d}')
        w('JDIM: GRID   SEGMENT    BCTYPE      ISTA      IEND      KSTA      KEND     NDATA')
        w(f'{1:10d}{1:10d}{1000:10d}{1:10d}{2:10d}{1:10d}{a2:10d}{0:10d}')
        w('K0:   GRID   SEGMENT    BCTYPE      ISTA      IEND      JSTA      JEND     NDATA')
        w(f'{1:10d}{1:10d}{0:10d}{1:10d}{2:10d}{1:10d}{a3:10d}{0:10d}')
        if ivisc == 0:
            # inviscid surface (flow tangency); no wall data to follow
            w(f'{1:10d}{2:10d}{1005:10d}{1:10d}{2:10d}{a3:10d}{a5:10d}{0:10d}')
        else:
            # viscous solid surface, adiabatic (Twtype = 0)
            w(f'{1:10d}{2:10d}{2004:10d}{1:10d}{2:10d}{a3:10d}{a5:10d}{2:10d}')
            w('    Twtype        Cq')
            w(f'{0.0:10.4f}{0.0:10.4f}')
        w(f'{1:10d}{3:10d}{0:10d}{1:10d}{2:10d}{a5:10d}{a1:10d}{0:10d}')
        w('KDIM: GRID   SEGMENT    BCTYPE      ISTA      IEND      JSTA      JEND     NDATA')
        w(f'{1:10d}{1:10d}{1000:10d}{1:10d}{2:10d}{1:10d}{a1:10d}{0:10d}')
        w('      MSEQ    MGFLAG    ICONSF       MTT      NGAM')
        w(f'{mseq:10d}{1:10d}{0:10d}{0:10d}{2:10d}')
        w('      ISSC EPSSSC(1) EPSSSC(2) EPSSSC(3)      ISSR EPSSSR(1) EPSSSR(2) EPSSSR(3)')
        w(f'{0:10d}{0.3:10.1f}{0.3:10.1f}{0.3:10.1f}{0:10d}'
          f'{0.3:10.1f}{0.3:10.1f}{0.3:10.1f}')
        w('      NCYC    MGLEVG     NEMGL     NITFO')
        for lev in range(1, mseq + 1):
            w(f'{ncyc:10d}{lev:10d}{0:10d}{nitfo:10d}')
        w('      MIT1      MIT2      MIT3      MIT4      MIT5      MIT6      MIT7      MIT8')
        for _ in range(mseq):
            w(f'{1:10d}' * 8)
        w('   1-1 BLOCKING DATA:')
        w('       NBLI')
        w(f'{1:10d}')
        w(' NUMBER   GRID     :    ISTA   JSTA   KSTA   IEND   JEND   KEND  ISVA1  ISVA2')
        w(f'{1:7d}{1:7d}      {1:8d}{1:7d}{1:7d}{2:7d}{a3:7d}{1:7d}{1:7d}{2:7d}')
        w(' NUMBER   GRID     :    ISTA   JSTA   KSTA   IEND   JEND   KEND  ISVA1  ISVA2')
        w(f'{1:7d}{1:7d}      {1:8d}{a1:7d}{1:7d}{2:7d}{a5:7d}{1:7d}{1:7d}{2:7d}')
        w('   PATCH SURFACE DATA:')
        w('    NINTER')
        w('     0')
        w('   PLOT3D OUTPUT:')
        w('    GRID IPTYPE ISTART   IEND   IINC JSTART   JEND   JINC KSTART   KEND   KINC')
        w(f'{1:6d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}')
        w(' IMOVIE')
        w('     0')
        w('   PRINT OUT:')
        w('    GRID IPTYPE ISTART   IEND   IINC JSTART   JEND   JINC KSTART   KEND   KINC')
        w(f'{1:6d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}{0:7d}')
        w('   CONTROL SURFACE:')
        w('   NCS')
        w('     0')
        w('    GRID ISTART   IEND   JSTART   JEND   KSTART   KEND  IWALL  INORM')

        with open(fname, 'w') as f:
            f.write('\n'.join(L) + '\n')

    def write(self, outdir, **inp_kwargs):
        '''Write ``cfl3d.xyz`` and ``cfl3d.inp`` into ``outdir``.'''
        os.makedirs(outdir, exist_ok=True)
        self.write_xyz(os.path.join(outdir, 'cfl3d.xyz'))
        self.write_inp(os.path.join(outdir, 'cfl3d.inp'), **inp_kwargs)
        return outdir

    def write_tecplot(self, fname):
        '''Dump the grid as an ASCII Tecplot/paraview-readable block.'''
        with open(fname, 'w') as f:
            f.write('VARIABLES = "X" "Y"\n')
            f.write(f'ZONE I={self.nj} J={self.nk} F=POINT\n')
            for k in range(self.nk):
                for j in range(self.nj):
                    f.write(f'{self.xy[j, k, 0]:20.10f}'
                            f'{self.xy[j, k, 1]:20.10f}\n')


# ---------------------------------------------------------------------------
#  command line
# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description='Generate a CFL3D C-grid for a sharp-TE airfoil with gmsh')
    p.add_argument('airfoil', help='airfoil coordinate file (Selig layout)')
    p.add_argument('--out', default='.', help='output directory')
    p.add_argument('--n-foil', type=int, default=101)
    p.add_argument('--n-wake', type=int, default=41)
    p.add_argument('--n-grow', type=int, default=97)
    p.add_argument('--r-far', type=float, default=80.0)
    p.add_argument('--y-plus', type=float, default=1.0)
    p.add_argument('--zone-layers', type=int, default=4,
                   help='k layers per transfinite zone (fewer = more collars, '
                        'closer to the wall normal, slower)')
    p.add_argument('--collar-max', type=float, default=0.03,
                   help='distance in chords out to which collars are marched; '
                        'beyond it one transfinite zone reaches the far field')
    p.add_argument('--te-relax', type=float, default=0.3,
                   help='trailing-edge rounding applied to the marching front')
    p.add_argument('--smooth', type=int, default=0,
                   help='Winslow/Thomas-Middlecoff sweeps over the outer zones '
                        '(off by default, see the module docstring)')
    p.add_argument('--mach', type=float, default=0.725)
    p.add_argument('--alpha', type=float, default=2.10)
    p.add_argument('--re-mil', type=float, default=6.5,
                   help='Reynolds number in millions, for cfl3d.inp')
    p.add_argument('--tinf', type=float, default=460.0)
    p.add_argument('--ncyc', type=int, default=1000)
    p.add_argument('--tecplot', action='store_true',
                   help='also write grid.dat for plotting')
    a = p.parse_args(argv)

    x, y = read_airfoil(a.airfoil)
    cg = CGridGmsh(n_foil=a.n_foil, n_wake=a.n_wake, n_grow=a.n_grow,
                   y_plus=a.y_plus, r_far=a.r_far,
                   re_chord=a.re_mil * 1e6,
                   zone_layers=a.zone_layers, collar_max=a.collar_max,
                   te_relax=a.te_relax, n_smooth=a.smooth)
    print(f'generating C-grid for {a.airfoil}')
    cg.generate(x, y)
    cg.write(a.out, mach=a.mach, alpha=a.alpha, re_mil=a.re_mil,
             tinf=a.tinf, ncyc=a.ncyc)
    if a.tecplot:
        cg.write_tecplot(os.path.join(a.out, 'grid.dat'))
    print(f'  written         : {os.path.join(a.out, "cfl3d.xyz")}, '
          f'{os.path.join(a.out, "cfl3d.inp")}')


if __name__ == '__main__':
    main()
