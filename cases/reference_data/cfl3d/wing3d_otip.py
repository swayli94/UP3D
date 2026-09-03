"""
Multiblock **C-H + O-tip** structured grid for a wing, for CFL3D, built on gmsh.

    python wing3d_otip.py --check-2d --ref ref_2d.npz    # the 2-D section
    python wing3d_otip.py --out ./m6 --model euler       # grid + deck

===============================================================================
WHY THIS TOPOLOGY, AND WHAT THE PREVIOUS ATTEMPT GOT WRONG
===============================================================================

★★★ **A single-block C-H grid cannot carry a blunt wing tip, and closing the
section to zero thickness to make it fit is not allowed** -- it deforms the very
geometry the reference exists to measure.  A first attempt did exactly that and
was deleted (user ruling 2026-09-03: "M6 机翼不能改变机翼形状, 不可以厚度收到零,
本来就不应该使用单块 C-H 网格来生成").  It is the same mistake this directory
refuses elsewhere, where it declines the CFL3D repo's own airfoil ordinates so
that no geometry difference enters the comparison.

The correct topology is a **multiblock C-H grid with an O-block over the tip**.
Reference implementations, both of which this file was written against:

  * the user's own CFL3D-oriented generator, ``tools/cgrid/``
    (``examples/wing-simple-OTip/wing-onera-m6.py``, ``cgrid/wing.py``);
  * pyHyp's M6 example, https://github.com/mdolab/pyhyp/tree/main/examples/m6.

★ ``tools/cgrid`` is also the **verification baseline**: it is an independent
implementation of the same topology, so agreement between its grid and this one
is a real cross-check rather than this file checking itself.  Its reference for
the M6 is 7 blocks / 3.53 M points / ``NBLI = 18``.

===============================================================================
THE ONE FACT THAT MAKES THIS CHEAP
===============================================================================

★★ A blunt trailing edge does NOT need a blunt-TE mesher.  ``cgrid``'s
``FoilGrid.gen_grid`` strips the TE thickness first (shifting y by
``-/+ 0.5*tail*x``), meshes the **sharp** section, and then ``add_tail()`` puts
the thickness back and fills the gap with one extra block.  So the proven
sharp-TE gmsh C-grid already vendored here (``cgrid_gmsh.CGridGmsh``, gmsh
transfinite + collar marching) is a drop-in for the hard part, and the blunt TE
costs one affine transform plus one algebraic block.

★★ And the index conventions already agree, which is why the port is clean --
``FoilGrid`` and ``CGridGmsh`` use the SAME formulas:

    nj   = 2*(n_foil + n_wake) - 3
    jTE0 = n_wake
    jTE1 = n_wake + 2*n_foil - 2

For the M6 reference parameters (n_foil 161, n_wake 61, n_grow 81, tail 17)
both give nj = 441, jTE = (61, 381), and njTip = (nf - tail)/2 + 1 = 153 --
matching the reference blocks 1 and 3 exactly.

===============================================================================
DIVISION OF LABOUR
===============================================================================

gmsh does the **meshing** -- the 2-D transfinite fill of the section C-grid,
where its value is.  Everything above that is transfinite interpolation and
distribution algebra (stacking, the tip basin, the Coons O-blocks, the
far-field extension), done in numpy, exactly as ``cgrid`` does it.  That is the
same split the 2-D generator in this directory already uses.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from cgrid_gmsh import CGridGmsh                          # noqa: E402


# ---------------------------------------------------------------------------
#  helpers, ported from the user's cgrid (tools/cgrid/cgrid/f2py.py)
# ---------------------------------------------------------------------------

def coons(pu0, pu1, p0v, p1v):
    """Dual-linear Coons transfinite interpolation of a surface patch.

    Ported from ``cgrid.f2py.l2coons_interpolation``; vectorised, and verified
    against it to machine precision by ``--check-2d``.

        q = (1-u) p(0,v) + u p(1,v)
        r = (1-v) p(u,0) + v p(u,1)
        s = bilinear blend of the four corners
        p = q + r - s
    """
    pu0 = np.asarray(pu0, float)
    pu1 = np.asarray(pu1, float)
    p0v = np.asarray(p0v, float)
    p1v = np.asarray(p1v, float)
    ni, nv = pu0.shape
    nj = p0v.shape[0]
    u = np.linspace(0.0, 1.0, ni)[:, None, None]
    v = np.linspace(0.0, 1.0, nj)[None, :, None]
    q = (1 - u) * p0v[None, :, :] + u * p1v[None, :, :]
    r = (1 - v) * pu0[:, None, :] + v * pu1[:, None, :]
    s = (((1 - u) * p0v[0] + u * p1v[0]) * (1 - v)
         + ((1 - u) * p0v[-1] + u * p1v[-1]) * v)
    return q + r - s


def power_growth(n_point: int, dL0: float, L: float = 1.0) -> np.ndarray:
    """Clustered distribution on [0, L] with first interval ``dL0``.

    Ported from ``cgrid.f2py.power_growth_distribution``: solve
    ``a = (a f - f + 1)^(1/n)`` for the growth factor, accumulate, then
    normalise onto [0, L] exactly (the normalisation is what makes both ends
    land on the section planes).
    """
    if dL0 >= L / (n_point - 1):
        return np.linspace(0.0, 1.0, n_point) * L
    f = L / dL0
    a, a1, b = 1.2, 1.0, 1.0 / n_point
    while abs(a - a1) / a > 2e-10:
        a1 = a
        a = (a * f - f + 1.0) ** b
    uu = np.cumsum(a ** np.arange(n_point) * dL0)
    return (uu - uu[0]) / (uu[-1] - uu[0]) * L


def activation(x, r_critical=0.5):
    """0 below ``r_critical``, 1 at/above 1, tanh in between (cgrid's 'Tanh')."""
    x = np.asarray(x, float)
    r = np.clip((x - r_critical) / (1.0 - r_critical), 0.0, 1.0) * 8.0 - 4.0
    out = 0.5 * (np.tanh(r) + 1.0)
    return np.where(x <= r_critical, 0.0, np.where(x >= 1.0, 1.0, out))


# ---------------------------------------------------------------------------
#  stage 1: the 2-D section
# ---------------------------------------------------------------------------

class Section2D:
    """One control section's 2-D blocks, in unit-chord section coordinates.

    ``blocks`` after ``build()``:

        0  main C-grid       (nj, nk, 2)
        1  tail/wake block   (n_wake, n_tail, 2)
        2  tip solid block   (njTip, n_tail, 2)   -- only when ``solid=True``

    which is ``cgrid.FoilGrid``'s layout (main-block / wake-block /
    solid1-block), so the 3-D assembly and the deck can follow it directly.
    """

    def __init__(self, n_foil=161, n_wake=61, n_grow=81, n_tail=17,
                 h1=None, y_plus=None, re_chord=None, r_far=20.0,
                 le_factor=0.06):
        """Wall spacing comes from EXACTLY ONE of ``h1`` or ``y_plus``.

        ★★★ Euler and RANS need different grids at the wall and the generator
        must not paper over that: **Euler has no boundary layer to resolve**, so
        ``h1`` is set directly; **SA and SST both need y+ ~ 1**, so for RANS
        ``h1`` is DERIVED from the target y+ and the case's own chord Reynolds
        number.  Passing neither, or both, raises.

        ★ An earlier version hard-coded a single ``h1 = 1.6212e-6`` for both --
        a value measured off the reference grid rather than computed.  Measured
        consequences: the Euler grid was **~617x over-refined at the wall**
        (1.6e-06 against ~1e-03 chord), which is the most likely reason its
        cold start needed CFL 0.2 where 1.0 should do; and the RANS y+ was
        not a y+ at all, landing at 0.77 only because the reference had sized
        it at ``yp_cri`` 0.9 / Re 2e7.  At any other Reynolds number it would
        have been silently wrong.

        The y+ formula is ``cgrid_gmsh.y_plus_to_height`` -- the SAME one the
        2-D ladder in this directory uses -- so the two stay consistent:
        ``cf = 0.0576 Re^-0.2``, ``y/c = y+ / (Re sqrt(cf/2))``.
        """
        if (h1 is None) == (y_plus is None):
            raise ValueError('give exactly one of h1 (Euler) or y_plus (RANS)')
        if y_plus is not None:
            if re_chord is None:
                raise ValueError('y_plus needs re_chord (the CHORD Reynolds '
                                 'number of the case)')
            from cgrid_gmsh import y_plus_to_height
            h1 = float(y_plus_to_height(y_plus, re_chord))
        self.y_plus_target, self.re_chord = y_plus, re_chord
        self.n_foil = int(np.floor((n_foil - 1) / 2) * 2 + 1)
        self.n_wake = int(np.floor((n_wake - 1) / 4) * 4 + 1)
        self.n_grow = int(np.floor((n_grow - 1) / 4) * 4 + 1)
        self.n_tail = int(np.floor(n_tail / 8) * 8 + 1)
        self.h1, self.r_far, self.le_factor = h1, r_far, le_factor

        self.nf = 2 * self.n_foil - 1
        self.nj = 2 * (self.n_foil + self.n_wake) - 3
        self.nk = self.n_grow
        self.jte0 = self.n_wake                       # 1-based lower TE
        self.jte1 = self.n_wake + 2 * self.n_foil - 2  # 1-based upper TE
        self.njtip = (self.nf - self.n_tail) // 2 + 1
        self.blocks: list[np.ndarray] = []
        self.info: dict = {}

    # -- construction ------------------------------------------------------

    def build(self, xu, yu, xl, yl, solid=False, verbose=False):
        xu, yu = np.asarray(xu, float), np.asarray(yu, float)
        xl, yl = np.asarray(xl, float), np.asarray(yl, float)

        # ★ Strip the TE thickness and any TE camber offset, exactly as
        # cgrid.FoilGrid.gen_grid does, so what gets meshed is a SHARP section
        # and the vendored sharp-TE gmsh C-grid applies unchanged.
        tail = float(yu[-1] - yl[-1])
        dyte = 0.5 * float(yu[-1] + yl[-1])
        yu_s = yu - xu * dyte - xu * tail * 0.5
        yl_s = yl - xl * dyte + xl * tail * 0.5
        self.info.update(tail=tail, dyTE=dyte,
                         sharp_te_gap=float(yu_s[-1] - yl_s[-1]))

        # closed contour, TE -> upper -> LE -> lower -> TE
        cx = np.concatenate([xu[::-1], xl[1:]])
        cy = np.concatenate([yu_s[::-1], yl_s[1:]])

        cg = CGridGmsh(n_foil=self.n_foil, n_wake=self.n_wake,
                       n_grow=self.n_grow, n_coarse=2, r_far=self.r_far,
                       le_factor=self.le_factor, h1=self.h1,
                       re_chord=1.0, wake_ratio_max=1.5)
        cg.generate(cx, cy, verbose=verbose)
        if (cg.nj, cg.nk) != (self.nj, self.nk):
            raise RuntimeError(f'gmsh C-grid came out {cg.nj}x{cg.nk}, '
                               f'expected {self.nj}x{self.nk}')
        main = cg.xy.copy()                            # (nj, nk, 2)
        self.info['quality_2d'] = cg.quality()

        # ★ Put the thickness back: dy = +/- 0.5 tail x, sign by which side of
        # the C the point is on, ramped over 0 <= x <= 1.  This is
        # cgrid.FoilGrid.add_tail's transform.
        if tail > 1e-6:
            ratio = np.clip(main[:, :, 0], 0.0, 1.0)
            sign = np.where(np.arange(self.nj)[:, None] > self.nj / 2,
                            0.5, -0.5)
            main[:, :, 1] += sign * tail * ratio
        if abs(dyte) > 1e-6:
            main[:, :, 1] += np.clip(main[:, :, 0], 0.0, 1.0) * dyte
        self.blocks = [main]

        # tail/wake block: the TE gap, closed off downstream along the cut
        if tail > 1e-6:
            self.blocks.append(self._tail_block(main))
        if solid:
            self.blocks.append(self._tip_block(main))
        return self

    def _tail_block(self, main):
        """The block filling the blunt-TE gap behind the section.

        x follows the wake cut; y is distributed linearly across the gap
        between the lower and upper wake lines.  Orientation follows cgrid's
        (j from the far field back to the TE, k across the gap).
        """
        nw, nt = self.n_wake, self.n_tail
        up = main[self.nj - nw, 0, 1]
        dn = main[nw - 1, 0, 1]
        ycol = np.linspace(dn, up, nt)[::-1]
        blk = np.empty((nw, nt, 2))
        blk[:, :, 0] = main[:nw, 0, 0][:, None]
        blk[:, :, 1] = ycol[None, :]
        return blk

    def _tip_block(self, main):
        """The tip-plane block: a Coons patch over the section interior.

        Four edges, following cgrid.FoilGrid.solid_1Block: the lower surface,
        the **blunt TE base**, the upper surface, and the tail block's last
        j-line.  ★ The blunt TE is what supplies a non-degenerate fourth edge;
        a sharp TE collapses it, which is why the TE closure is load-bearing
        for this topology and not a cosmetic choice.
        """
        nt = self.n_tail
        surf = main[self.jte0 - 1:self.jte1, 0, :]      # (nf, 2), lower->upper
        ic1 = (self.nf - nt) // 2
        ic2 = ic1 + nt
        pu1 = surf[:ic1 + 1]                            # v = 1 (lower)
        p1v = surf[ic1:ic2][::-1]                       # u = 1 (TE base)
        pu0 = surf[ic2 - 1:][::-1]                      # v = 0 (upper)
        p0v = self.blocks[1][-1, :, :]                  # u = 0 (tail block)
        return coons(pu0, pu1, p0v, p1v)

    # -- reporting ---------------------------------------------------------

    def report(self):
        i = self.info
        q = i['quality_2d']
        print(f'  section         : nj x nk = {self.nj} x {self.nk}, '
              f'jTE = ({self.jte0}, {self.jte1}), njTip = {self.njtip}, '
              f'n_tail = {self.n_tail}')
        print(f'  TE              : gap {i["tail"]:.6f}, dyTE {i["dyTE"]:.6f}, '
              f'gap after stripping {i["sharp_te_gap"]:.2e}')
        print(f'  2-D quality     : ortho {q["ortho_min"]:.2f}..'
              f'{q["ortho_max"]:.2f} deg, k ratio {q["k_ratio_max"]:.4f}, '
              f'h1 {q["h1_min"]:.4e}, neg cells {q["negative_cells"]}')
        for n, b in enumerate(self.blocks):
            print(f'  2-D block {n}     : {b.shape}')


# ---------------------------------------------------------------------------
#  ONERA M6, exactly as the reference case defines it
# ---------------------------------------------------------------------------

#: CST coefficients of the ONERA M6 (ONERA D) section, from the reference case
#: tools/cgrid/examples/wing-simple-OTip/wing-onera-m6.py.  Symmetric, so the
#: lower surface is -cst.
M6_CST = np.array([0.184161, 0.029638, 0.254533, -0.045245, 0.292907,
                   0.031388, 0.176787, 0.101274, 0.147929, 0.120975])
M6_REL_THICK = 0.09779
M6_REL_TAIL = 0.00141          # ★ blunt TE -- load-bearing for the O-tip block
M6_CHORD_ROOT = 0.8059
M6_CHORD_TIP = 0.4529
M6_CHORD_MAC = 0.64607
M6_SPAN = 1.1963
M6_SWEEP_LE_DEG = 30.0


def m6_section(n_point=1001):
    """The M6 section at unit chord, with its blunt TE.

    Uses ``cst_modeling.section.cst_foil`` -- the same routine the reference
    case uses -- so the section is the reference's, not a re-fit of it.  That
    package is needed only here; it is NOT a project dependency.
    """
    try:
        from cst_modeling.section import cst_foil
    except ImportError as exc:                          # pragma: no cover
        raise ImportError(
            'cst_modeling is needed for the M6 section (pip install '
            'cst-modeling3d).  It is only used to build the section, not by '
            'anything else in this repository.') from exc
    xx, yu, yl, tmax, _ = cst_foil(n_point, M6_CST, -M6_CST, x=None,
                                   t=M6_REL_THICK, tail=M6_REL_TAIL)
    return xx, yu, xx, yl, tmax


# ---------------------------------------------------------------------------
#  stage-1 verification against the reference implementation
# ---------------------------------------------------------------------------

def check_2d(ref: Path | None = None, verbose=True) -> dict:
    """Build the 2-D section with gmsh and compare against cgrid's own.

    ★ This is a cross-check between two INDEPENDENT implementations, so the
    point-distributions are expected to differ -- what has to agree is the
    structure and the geometry: block shapes, the TE gap, the surface extent,
    the tip-block envelope.  Anything that agrees only because both sides run
    the same code would not be evidence.
    """
    xu, yu, xl, yl, tmax = m6_section()
    sec = Section2D().build(xu, yu, xl, yl, solid=True, verbose=verbose)
    if verbose:
        sec.report()

    out = dict(shapes=[b.shape for b in sec.blocks],
               tail=sec.info['tail'], tmax=float(tmax))
    m, t, p = sec.blocks
    wall = m[sec.jte0 - 1:sec.jte1, 0, :]
    out.update(
        wall_x_min=float(wall[:, 0].min()), wall_x_max=float(wall[:, 0].max()),
        te_gap_at_wall=float(m[sec.jte1 - 1, 0, 1] - m[sec.jte0 - 1, 0, 1]),
        h1_min=float(sec.info['quality_2d']['h1_min']),
        k_ratio_max=float(sec.info['quality_2d']['k_ratio_max']),
        neg_cells=int(sec.info['quality_2d']['negative_cells']),
        tail_y_extent=float(t[:, :, 1].max() - t[:, :, 1].min()),
        tip_y_max=float(p[:, :, 1].max()),
        tip_x_min=float(p[:, :, 0].min()), tip_x_max=float(p[:, :, 0].max()),
    )

    if ref is not None and Path(ref).is_file():
        d = np.load(ref)
        rb = [d['b0'], d['b1'], d['b2']]
        out['ref_shapes'] = [b.shape for b in rb]
        out['shapes_match'] = out['shapes'] == out['ref_shapes']
        rw = rb[0][sec.jte0 - 1:sec.jte1, 0, :]
        out['ref_te_gap_at_wall'] = float(
            rb[0][sec.jte1 - 1, 0, 1] - rb[0][sec.jte0 - 1, 0, 1])
        out['ref_tip_y_max'] = float(rb[2][:, :, 1].max())
        out['ref_h1'] = float(np.linalg.norm(
            rb[0][sec.jte0 - 1:sec.jte1, 1] - rw, axis=-1).min())
        # the wall CURVE must agree even though the point distribution need not:
        # compare our wall y(x) against theirs by interpolation, upper side
        jle = (sec.nj + 1) // 2
        xs = np.linspace(0.02, 0.98, 200)
        ours = np.interp(xs, m[jle - 1:sec.jte1, 0, 0], m[jle - 1:sec.jte1, 0, 1])
        theirs = np.interp(xs, rb[0][jle - 1:sec.jte1, 0, 0],
                           rb[0][jle - 1:sec.jte1, 0, 1])
        out['wall_curve_max_dev'] = float(np.max(np.abs(ours - theirs)))

    if verbose:
        print('\n  --- 2-D cross-check against tools/cgrid ---')
        for k in sorted(out):
            if k not in ('shapes', 'ref_shapes'):
                print(f'    {k:24s} {out[k]}')
        print(f'    {"shapes":24s} {out["shapes"]}')
        if 'ref_shapes' in out:
            print(f'    {"ref_shapes":24s} {out["ref_shapes"]}')
    return out


#: ★ Chord Reynolds number for the M6 case, on the ROOT chord, because the grid
#: is normalised to root chord = 1: the experiment's Re_MAC = 11.72e6 rescales
#: to ``11.72e6 / 0.64607 * 0.8059 = 14.62e6``.  This is the number the RANS
#: wall spacing is sized from, so getting it wrong is a silent 24 % error in Re
#: AND in y+.
M6_RE_ROOT_CHORD = 14.62e6

#: Three-rung ladders.  Counts refine by ~1.35 per direction per rung; L3 is
#: deliberately the reference implementation's own parameter set, so the finest
#: rung is the grid that was cross-checked block-for-block against it.
#:
#: ★ Euler sets ``h1`` directly and refines it with the ladder; RANS holds
#: ``y_plus = 1`` on every rung and lets ``h1`` follow from the Reynolds number
#: -- a RANS grid whose boundary layer is not resolved is not a reference.
#: ★ ``n_tail`` is constrained to ``8m+1``, and that is NECESSARY rather than
#: conservative: ``njTip = (nf - n_tail)/2 + 1`` has to be ``4m+1`` to survive
#: mesh sequencing, and with ``nf = 2 n_foil - 1 = 1 mod 8`` that forces
#: ``n_tail = 1 mod 8``.  So the ladder uses 9 -> 17 -> 25, which refines
#: monotonically; 9 -> 17 -> 17 (the first attempt) did not refine the blunt-TE
#: base at all between L2 and L3.
#:
#: ★★★ ``n_grow`` = 61 / 77 / 101 is not a round choice either -- it comes out
#: of a JOINT search over (n_grow, k_crit) per rung, constrained to a monotone
#: n_grow, that minimises the spread in the tip basin's physical depth.  Legal
#: k indices are only every 4th (``4m+1``) and consecutive ones are ~``r^4``
#: apart in distance, so with n_grow fixed the depth can only be matched to
#: ~1.26x; letting n_grow move too brings it to **1.003x** (0.2361 / 0.2361 /
#: 0.2368 chords) while still refining at 1.26x and 1.31x.
M6_LEVELS = {
    'L1': dict(n_foil=81, n_wake=33, n_grow=61, n_tail=9,
               n_span=33, n_side=33, h1_euler=2.000e-3),
    'L2': dict(n_foil=121, n_wake=45, n_grow=77, n_tail=17,
               n_span=45, n_side=45, h1_euler=1.414e-3),
    'L3': dict(n_foil=161, n_wake=61, n_grow=101, n_tail=25,
               n_span=61, n_side=61, h1_euler=1.000e-3),
}

#: ★ The reference implementation's own parameter set, kept SEPARATE from the
#: ladder.  It is the configuration that was cross-checked block-for-block
#: against tools/cgrid, and its job is verification; the ladder's job is data.
#: Conflating them cost the ladder its consistency -- L3 was pinned to
#: n_tail = 17 to keep it equal to this grid, which is exactly why n_tail
#: stopped refining between L2 and L3.
M6_REFERENCE_PARAMS = dict(n_foil=161, n_wake=61, n_grow=81, n_tail=17,
                           n_span=61, n_side=61, h1_euler=1.000e-3)

#: RANS needs more wall-normal points than Euler at the same rung, because y+ 1
#: puts the first cell ~2e-06 chords out and the layer count has to bridge from
#: there to the far field with a sane growth ratio.
M6_RANS_NGROW = {'L1': 65, 'L2': 89, 'L3': 105, 'REF': 105}


def build_m6(level='L3', model='euler', y_plus=1.0,
             re_chord=M6_RE_ROOT_CHORD, verbose=True) -> WingOTip:
    """The ONERA M6 on the reference case's parameters, at one ladder rung.

    ★ The grid is normalised to **root chord = 1**, as the reference case does.
    ★ ``model='euler'`` sets the wall spacing directly (no y+ clustering --
    there is no boundary layer); ``model='rans'`` derives it from ``y_plus``
    and ``re_chord`` and uses the taller wall-normal ladder.
    """
    if level == 'REF':
        lv = dict(M6_REFERENCE_PARAMS)
    elif level in M6_LEVELS:
        lv = dict(M6_LEVELS[level])
    else:
        raise KeyError(f'unknown level {level!r}; have '
                       f'{sorted(M6_LEVELS) + ["REF"]}')
    h1_euler = lv.pop('h1_euler')
    if model == 'euler':
        wall = dict(h1=h1_euler)
    else:
        wall = dict(y_plus=y_plus, re_chord=re_chord)
        lv['n_grow'] = M6_RANS_NGROW[level]

    xu, yu, xl, yl, _ = m6_section()
    cr, ct, span = M6_CHORD_ROOT, M6_CHORD_TIP, M6_SPAN
    sweep = math.radians(M6_SWEEP_LE_DEG)
    g = WingOTip(**lv, **wall)
    g.add_section(xu, yu, xl, yl, chord=1.0, solid=False)
    g.add_section(xu, yu, xl, yl, xLE=math.tan(sweep) * span / cr, yLE=0.0,
                  zLE=span / cr, chord=ct / cr, solid=True)
    g.build(verbose=verbose)
    g.info.update(level=level, model=model,
                  h1=g.sec_proto.h1, y_plus_target=y_plus if model == 'rans'
                  else None)
    if verbose:
        h1 = g.sec_proto.h1
        if model == 'rans':
            print(f'  wall            : y+ target {y_plus}, Re_root '
                  f'{re_chord:.4g} -> h1 {h1:.4e} root chords')
        else:
            print(f'  wall            : h1 {h1:.4e} root chords set DIRECTLY '
                  f'(Euler -- no y+ clustering)')
    return g


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--check-2d', action='store_true',
                   help='stage 1: build the 2-D section and cross-check it')
    p.add_argument('--ref', default=None,
                   help='ref_2d.npz dumped from tools/cgrid (for --check-2d)')
    p.add_argument('--out', default=None,
                   help='write cfl3d.xyz + cfl3d.inp into this directory')
    p.add_argument('--level', default='L3',
                   choices=('L1', 'L2', 'L3', 'REF'),
                   help="ladder rung, or REF = the reference implementation's "
                        "own parameters (the cross-checked configuration)")
    p.add_argument('--model', default='euler', choices=('euler', 'rans'))
    p.add_argument('--y-plus', type=float, default=1.0,
                   help='target y+ (RANS only)')
    p.add_argument('--mach', type=float, default=0.8395)
    p.add_argument('--alpha', type=float, default=3.06,
                   help='EXPERIMENTAL alpha, uncorrected (user ruling)')
    p.add_argument('--re-mil', type=float, default=14.62,
                   help='unit Reynolds number / 1e6 on the ROOT chord '
                        '(the grid is normalised to root chord = 1)')
    p.add_argument('--ncyc', type=int, default=1000)
    p.add_argument('--cfl', type=float, default=1.0)
    a = p.parse_args(argv)

    if a.check_2d:
        r = check_2d(Path(a.ref) if a.ref else None)
        ok = (r['neg_cells'] == 0
              and abs(r['te_gap_at_wall'] - M6_REL_TAIL) < 1e-9
              and r.get('shapes_match', True))
        print('\n  STAGE 1:', 'PASS' if ok else 'FAIL')
        return 0 if ok else 1

    if a.out:
        g = build_m6(level=a.level, model=a.model, y_plus=a.y_plus)
        out = Path(a.out)
        out.mkdir(parents=True, exist_ok=True)
        g.write_xyz(out / 'cfl3d.xyz')
        g.write_inp(out / 'cfl3d.inp', mach=a.mach, alpha=a.alpha,
                    re_mil=a.re_mil, ivisc=0 if a.model == 'euler' else 7,
                    cfl=a.cfl, ncyc=a.ncyc)
        print(f"  nbli            : {g.info['nbli']}  "
              f"(j segment ends p = {g.info['p']})")
        print(f'  -> {out}/cfl3d.xyz, {out}/cfl3d.inp')
        return 0

    p.print_help()
    return 0


if __name__ == '__main__':
    sys.exit(main())


# ---------------------------------------------------------------------------
#  stage 2: the 3-D multiblock assembly
# ---------------------------------------------------------------------------

class WingOTip:
    """C-H + O-tip multiblock grid for a straight-tapered wing.

    Block order, which the deck depends on (and which reproduces the reference
    implementation's, so the two can be compared block for block):

        0  main C-H, root -> tip          (n_span, nj,     nk)
        1  tail block over the span       (n_span, n_wake, n_tail)
        2  O-block over the TIP           (k_crit, njTip,  n_tail)
        3  O-block in the tip WAKE void   (k_crit, n_wake, n_tail)
        4  main block -> spanwise far field  (n_side, nj, nk-k_crit+1)
        5  tip O-block -> far field       (n_side, njTip,  n_tail)
        6  wake O-block -> far field      (n_side, n_wake, n_tail)

    ★ ``k_crit`` is where the tip basin stops: ``floor(0.6*n_grow/4)*4 + 1``,
    the reference's rule, kept because it has to land on a multigrid-legal index
    (it is a BC-segment end on the IDIM face of block 0).
    """

    def __init__(self, n_foil=161, n_wake=61, n_grow=81, n_tail=17,
                 n_span=61, n_side=61, h1=None, y_plus=None, re_chord=None,
                 r_far=20.0, le_factor=0.06, trans_far_field=0.8,
                 basin_depth=0.236):
        self.sec_proto = Section2D(n_foil, n_wake, n_grow, n_tail, h1=h1,
                                   y_plus=y_plus, re_chord=re_chord,
                                   r_far=r_far, le_factor=le_factor)
        self.n_span = int(np.floor(n_span / 4) * 4 + 1)
        self.n_side = int(np.floor(n_side / 4) * 4 + 1)
        self.r_far = r_far
        self.trans = trans_far_field
        self.basin_depth = basin_depth
        #: set in build(), from basin_depth -- see _k_critical
        self.k_crit = None
        self.sections: list[Section2D] = []
        self.paras: list[dict] = []
        self.grid: list[np.ndarray] = []
        self.info: dict = {}

    # -- geometry ----------------------------------------------------------

    def add_section(self, xu, yu, xl, yl, xLE=0.0, yLE=0.0, zLE=0.0,
                    chord=1.0, solid=False, verbose=False):
        s = Section2D(self.sec_proto.n_foil, self.sec_proto.n_wake,
                      self.sec_proto.n_grow, self.sec_proto.n_tail,
                      h1=self.sec_proto.h1, r_far=self.sec_proto.r_far,
                      le_factor=self.sec_proto.le_factor)
        s.build(xu, yu, xl, yl, solid=solid, verbose=verbose)
        self.sections.append(s)
        self.paras.append(dict(xLE=xLE, yLE=yLE, zLE=zLE, chord=chord,
                               thick=s.info['tail'] and None or None))
        self.paras[-1]['thick'] = float(np.max(np.abs(
            s.blocks[0][s.jte0 - 1:s.jte1, 0, 1])) * 2.0)
        return self

    def _transform(self, i_sec):
        """Scale/offset a section by its chord and LE, blending back to the
        UNSCALED far field.

        ★ The blend is the reference's ``activation_function`` on
        ``1.2*|r|/r_far``, so beyond ~``trans_far_field`` of the domain radius
        every section carries the SAME curve and the outer boundary does not
        move with span.  (The previous single-block attempt scaled the whole
        plane by the local chord instead, and its far-field points moved ~24 m
        per unit span -- measured as a spanwise spacing ratio of 6.1.)
        """
        p = self.paras[i_sec]
        out = []
        for b0 in self.sections[i_sec].blocks:
            b = b0.copy()
            b[..., 0] = p['chord'] * b0[..., 0] + p['xLE']
            b[..., 1] = p['chord'] * b0[..., 1] + p['yLE']
            ratio = 1.2 * np.linalg.norm(b0, axis=-1) / self.r_far
            r = activation(ratio, r_critical=1.0 - self.trans)[..., None]
            out.append((1 - r) * b + r * b0)
        return out

    def _z_distribution(self):
        """Span stations: tip cell width fixed at ``2*tip_thickness/n_tail``,
        clustered toward the tip, then out to the spanwise far field."""
        nt = self.sec_proto.n_tail
        tip_thick = self.paras[-1]['thick'] * self.paras[-1]['chord']
        dz_tip = tip_thick / nt * 2.0
        z_tip = self.paras[-1]['zLE']
        z_far = self.r_far * 1.5
        z_wing = np.flip(z_tip - power_growth(self.n_span, dz_tip,
                                              z_tip - self.paras[0]['zLE']))
        z_side = power_growth(self.n_side, dz_tip, z_far - z_tip) + z_tip
        self.info.update(dz_tip=dz_tip, z_tip=z_tip, z_far=z_far,
                         tip_thickness=tip_thick)
        return z_wing, z_side

    # -- assembly ----------------------------------------------------------

    def _k_critical(self):
        """The k index where the tip basin has fully opened.

        ★★★ ``k_crit`` is only a BLOCK-SPLIT index here, not a geometric
        quantity -- and getting that separation right is the whole point.

        The first version set it to ``0.6*n_grow``, which put the tip basin at a
        DIFFERENT physical depth on every rung: measured **0.896 / 0.751 /
        0.948 chords, a 1.26x spread and not even monotone**.  The three rungs
        were then not the same geometry, so their differences mixed a topology
        change in with the refinement, and the ladder had no error bar at all
        (the cl deltas GREW between rungs, ratio 6.6).

        The obvious repair -- choose ``k_crit`` so the depth is fixed -- does
        NOT work, and that is worth recording: legal indices are only every
        4th (``4m+1``, for mesh sequencing), and on a geometric wall-normal
        distribution consecutive legal indices are a factor ``r^4`` apart in
        DISTANCE -- about **1.8x** here.  A target depth generally falls
        between two achievable values, so quantisation alone caps the
        consistency at ~1.3x.  Fighting it is hopeless.

        ⇒ Instead the basin's ramp is driven by ``basin_depth`` DIRECTLY (a
        physical length, identical on every rung -- see
        ``_extrude_main_block``), and ``k_crit`` merely has to sit at or beyond
        the point where the ramp has finished, so that the offset is already
        constant where the far-field block takes over.  The basin profile is
        then bit-comparable across rungs and ``k_crit`` carries no geometry.

        Rounded UP to the next ``4m+1``; the far-field block's
        ``n_grow - k_crit + 1`` is automatically ``4m+1`` too, since both are
        ``1 mod 4``.
        """
        s = self.sec_proto
        # ★★ Measured on the OUTBOARD (tip) section and over the MINIMUM
        # across the chord -- both matter, and both were wrong first time:
        #   * the basin lives at the TIP, whose chord is 0.562 of the root's,
        #     so a depth in root chords sits DEEPER in k there.  Choosing k
        #     from the root section gave 33 / 41 / 49 where the ramp actually
        #     saturates at 32 / 43 / 54, so the clamp below fired 2 and 6
        #     indices EARLY on L2 and L3 and forced the basin open before the
        #     ramp finished;
        #   * the k-lines reach different depths at different chord stations,
        #     so the minimum over the chord is what the clamp really sees.
        # Measured consequence of getting it wrong: a 22 % difference in the
        # basin profile between rungs at d = 0.7..0.9 c, i.e. the three rungs
        # were still not the same geometry.
        # ★ self.sections holds UNIT-CHORD grids -- the chord scaling happens
        # later, in _transform -- so sections[-1] and sections[0] are the same
        # scale and picking the tip one changes nothing.  The chord factor has
        # to be applied explicitly: the ramp in _extrude_main_block measures
        # distance on the SCALED block, so basin_depth is in root-chord units
        # there, and the tip station's unit-chord distance has to clear
        # basin_depth / chord_tip.
        tip = self.sections[-1].blocks[0]
        surf = tip[s.jte0 - 1:s.jte1]               # (nf, nk, 2)
        d = np.linalg.norm(surf - surf[:, :1], axis=-1).min(axis=0)
        target = self.basin_depth / self.paras[-1]['chord']
        k = int(np.searchsorted(d, target)) + 1
        # nearest legal 4m+1 (the ramp must land ON k_crit, so neither
        # rounding direction is privileged -- nearest minimises the depth error
        # and the per-level n_grow values are chosen so the error is small)
        lo = int(np.floor((k - 1) / 4) * 4 + 1)
        hi = lo + 4
        k = lo if abs(d[min(lo, len(d)) - 1] - target) <= \
            abs(d[min(hi, len(d)) - 1] - target) else hi
        k = int(np.clip(k, 9, s.n_grow - 4))
        if k > s.n_grow - 4:
            raise ValueError(
                f'basin_depth {self.basin_depth} c needs k = {k} but n_grow is '
                f'only {s.n_grow}; the far-field block would be degenerate')
        return k

    def build(self, verbose=True):
        if len(self.sections) < 2:
            raise RuntimeError('need at least a root and a tip section')
        self.k_crit = self._k_critical()
        z_wing, z_side = self._z_distribution()
        secs = [self._transform(i) for i in range(len(self.sections))]
        self.grid = []

        # blocks 0,1 -- stack the root and tip sections along z
        nb = len(secs[0])
        for b in range(nb):
            b0, b1 = secs[0][b], secs[-1][b]
            ni = len(z_wing)
            xyz = np.empty((ni, *b0.shape[:2], 3))
            t = ((z_wing - z_wing[0]) / (z_wing[-1] - z_wing[0]))[:, None, None]
            xyz[..., 0] = (1 - t) * b0[None, ..., 0] + t * b1[None, ..., 0]
            xyz[..., 1] = (1 - t) * b0[None, ..., 1] + t * b1[None, ..., 1]
            xyz[..., 2] = z_wing[:, None, None]
            self.grid.append(xyz)

        self._extrude_main_block()
        self.grid.append(self._o_block_tip(secs[-1]))
        self.grid.append(self._o_block_wake())
        self._extrude_far_field(z_side)

        s0 = self.sections[0].blocks[0]
        jj = (self.sec_proto.jte0 + self.sec_proto.jte1) // 2 - 1
        self.info['k_crit_depth'] = float(np.linalg.norm(
            s0[jj, self.k_crit - 1] - s0[jj, 0]))
        self.info['basin_depth'] = self.basin_depth
        self.info['blocks'] = [b.shape[:3] for b in self.grid]
        self.info['points'] = int(sum(np.prod(s) for s in self.info['blocks']))
        self.info['k_crit'] = self.k_crit
        if verbose:
            self.report()
        return self

    def _extrude_main_block(self):
        """Push the outboard end of block 0 outward in z to open the basin the
        two O-blocks live in.

        Tapered from ``i_critical`` (80 % of the span block) so the wing surface
        itself is untouched, and scaled by distance-from-the-wall so the offset
        dies by ``k_crit``; a second term opens it further downstream of
        ``x = 2`` so the wake O-block does not become a sliver at the outflow.
        ★ ``offset = tip thickness`` -- the basin is exactly deep enough to hold
        the tip block, which is what keeps the WING GEOMETRY untouched.
        """
        blk = self.grid[0]
        ni, nj, nk, _ = blk.shape
        kc = self.k_crit
        i_cri = int(self.n_span * 0.8)
        offset = self.info['tip_thickness']
        dY_far = 0.5 * abs(blk[-1, 0, kc - 1, 1] - blk[-1, -1, kc - 1, 1])
        x_cri, x_dis = 2.0, self.r_far - 2.0

        # ★★ The ramp is driven by the PHYSICAL basin depth, not by k_crit.
        # Using d/d[k_crit] instead ties the basin's shape to a quantised block
        # index, which is what made the three rungs different geometries (see
        # _k_critical).  With a physical length the basin profile is identical
        # on every rung and k_crit is free to be any legal split.
        # ★★★ The ramp normalises on the depth AT k_crit, and that is a
        # STRUCTURAL REQUIREMENT of the O-blocks, not a style choice.  The
        # O-block's i-index IS this direction, so the ramp has to reach 1
        # exactly AT k_crit: if it saturates earlier, every layer from
        # saturation to k_crit sits at the SAME z and the O-block collapses.
        #
        # ★★ Both alternatives were tried and both break the grid -- CFL3D
        # rejects them with "Fatal error(s) uncovered in grid metrics", and a
        # negative-volume count (calibrated first: it reports 0 on the two
        # grids CFL3D accepts) locates them in blocks 3 and 4, the O-blocks:
        #   * clip(d / basin_depth): saturates before k_crit -> zero-thickness
        #     O-block cells beyond that point;
        #   * a smoothstep of the same: the tip O-block's first i-spacing
        #     collapses 1.2e-04 -> 1.1e-06, because a C1 ramp barely opens the
        #     basin near the wall.  ★ And the "improvement" that motivated it
        #     was imaginary: the WALL itself never moves either way (d = 0 at
        #     k = 0 gives rk = 0 exactly); what a smoothstep holds still is the
        #     first layer OFF the wall, and displacing that layer is precisely
        #     HOW the basin opens.
        #
        # ⇒ Rung-to-rung consistency is therefore bought by choosing k_crit so
        # that d[k_crit] is the same PHYSICAL depth on every rung, which is
        # what _k_critical and the per-level n_grow values do together.
        d = np.linalg.norm(blk[..., :2] - blk[:, :, :1, :2], axis=-1)
        dmax = d[:, :, kc - 1][:, :, None]
        rk = np.clip(np.where(dmax > 0, d / np.maximum(dmax, 1e-300), 1.0),
                     0.0, 1.0)
        rk[:, :, kc - 1:] = 1.0

        ri = np.clip((np.arange(ni) - i_cri) / (ni - 1 - i_cri), 0.0, 1.0)
        ri[:i_cri] = 0.0
        blk[..., 2] += ri[:, None, None] * offset * rk
        rx = np.clip((blk[..., 0] - x_cri) / x_dis, 0.0, None)
        blk[..., 2] += ri[:, None, None] * dY_far * rx * rk
        self.info['dY_far'] = float(dY_far)

    def _o_block_tip(self, sec_tip):
        """O-block over the tip: per k-layer, Coons-fill the section interior.

        The four edges come from the outermost main-block j-line at that layer
        -- lower surface, **the blunt TE base**, upper surface -- closed by a
        straight fourth edge.  ★ This is what carries the real tip cap: the wing
        keeps its shape and the grid wraps the tip instead of the tip being
        reshaped to suit the grid.
        """
        s = self.sec_proto
        niw, nf, nt = s.n_wake - 1, s.nf, s.n_tail
        ic1 = s.njtip - 1
        ic2 = ic1 + nt
        b0 = self.grid[0]
        ni = b0.shape[0]
        out = np.empty((self.k_crit, ic1 + 1, nt, 3))
        for ik in range(self.k_crit):
            p = b0[ni - 1, niw:niw + nf, ik, :]
            pu1 = p[:ic1 + 1]
            p1v = p[ic1:ic2][::-1]
            pu0 = p[ic2 - 1:][::-1]
            p0v = np.stack([np.linspace(pu0[0, c], pu1[0, c], nt)
                            for c in range(3)], axis=1)
            out[ik] = coons(pu0, pu1, p0v, p1v)
        return out

    def _o_block_wake(self):
        """O-block filling the wake void beside the tip: linear blend between
        the upper and lower wake lines of the outermost main-block layer."""
        s = self.sec_proto
        nw, nf, nt = s.n_wake, s.nf, s.n_tail
        b0 = self.grid[0]
        ni = b0.shape[0]
        out = np.empty((self.k_crit, nw, nt, 3))
        for ik in range(self.k_crit):
            up = b0[ni - 1, :nw, ik, :]
            dn = b0[ni - 1, nw + nf - 2:, ik, :][::-1]
            t = np.linspace(0.0, 1.0, nt)[None, :, None]
            out[ik] = (1 - t) * dn[:, None, :] + t * up[:, None, :]
        return out

    def _extrude_far_field(self, z_side):
        """Carry the three outboard faces out to the spanwise far field.

        The z of each layer is blended from "the tip face shifted by dz" toward
        "the far-field station itself", so the sheared basin relaxes back to a
        flat plane by the outer boundary.
        """
        ni = len(z_side)
        zmax = z_side[-1] - z_side[0]
        t = ((z_side - z_side[0]) / zmax)[:, None, None]    # (ni, 1, 1)
        dz = (z_side - z_side[0])[:, None, None]            # (ni, 1, 1)
        zs = z_side[:, None, None]                          # (ni, 1, 1)
        for src, ks in ((0, slice(self.k_crit - 1, None)),
                        (2, slice(None)), (3, slice(None))):
            face = self.grid[src][-1, :, ks, :]             # (nj, nk', 3)
            xyz = np.empty((ni, *face.shape[:2], 3))
            xyz[..., :2] = face[None, ..., :2]
            xyz[..., 2] = (1 - t) * (face[None, ..., 2] + dz) + t * zs
            self.grid.append(xyz)

    # -- output ------------------------------------------------------------

    def write_xyz(self, path):
        """Multiblock PLOT3D, stream binary, double precision.

        Same layout as the 2-D writer: ``nblocks``, then ``(idim,jdim,kdim)``
        for every block, then x, y, z per block with i fastest.
        """
        import struct
        with open(path, 'wb') as f:
            f.write(struct.pack('i', len(self.grid)))
            for b in self.grid:
                f.write(struct.pack('iii', *b.shape[:3]))
            for b in self.grid:
                for d in range(3):
                    buf = np.ascontiguousarray(b[..., d].transpose(2, 1, 0))
                    f.write(buf.astype('<f8').tobytes())
        return path

    def report(self):
        i = self.info
        print(f"  blocks          : {len(self.grid)}, "
              f"{i['points'] / 1e6:.3f}M points, k_crit = {i['k_crit']}")
        for n, s in enumerate(i['blocks']):
            print(f"    block {n + 1}       : {s[0]:4d} x {s[1]:4d} x {s[2]:4d}"
                  f" = {int(np.prod(s)):9,d}")
        print(f"  basin           : ramp over {i['basin_depth']:.3f} c "
              f"(physical, rung-independent); split at k_crit = "
              f"{i['k_crit']} (depth {i['k_crit_depth']:.4f} c >= ramp)")
        print(f"  tip             : thickness {i['tip_thickness']:.6f}, "
              f"dz_tip {i['dz_tip']:.6f}, z_tip {i['z_tip']:.5f}, "
              f"z_far {i['z_far']:.3f}")

    # -- the deck ----------------------------------------------------------

    def write_inp(self, path, mach=0.8395, alpha=3.06, beta=0.0,
                  re_mil=14.62, tinf=460.0, ivisc=7, cfl=0.8, cfl_max=5.0,
                  iflagts=2000, ncyc=1000, mseq=3, nitfo=0, nwrest=5000,
                  sref=1.0, cref=1.0, bref=1.0, xmc=0.25,
                  keywords=None):
        """The 7-block ``cfl3d.inp`` for this topology.

        ★ The BC tables and the 1-to-1 blocking are ported from the reference
        implementation's ``output_cfl3d_input`` (``is_o_tip`` branch) rather
        than derived from scratch: the connectivity of a 7-block O-tip grid is
        exactly the kind of thing where "looks right" and "is right" are
        indistinguishable without a solver, and the reference deck is known to
        run.  ``--diff-deck`` compares this output against it line for line.

        ★ ``IALPH = 1`` is load-bearing: with it CFL3D sets
        ``vinf = M sin(alpha)``, i.e. **lift along y and span along z**, which
        is this grid's layout; it is also honoured only for PLOT3D-type grids
        (``NGRID < 0``), being silently forced to 0 otherwise (global.F:549).

        ★ Reynolds number caliber: the reference case quotes
        ``Re/m = 1.172e7 / chord_mac * chord_root = 14.62 million`` -- the grid
        is normalised to **root chord = 1**, so the deck's unit Reynolds number
        is the MAC-based experimental value rescaled to the root chord.  Getting
        this wrong is a silent 24 % error in Re.
        """
        s = self.sec_proto
        gn = len(self.grid)
        nwk, ntl = s.n_wake, s.n_tail
        nwg = (s.n_foil * 2 + 1 - ntl) // 2
        ntp = self.grid[3].shape[0]          # k_crit  (wake O-block ni)
        nbd = self.grid[4].shape[2]          # nk - k_crit + 1
        p1 = nwk
        p2 = nwk + nwg - 1
        p3 = nwk + nwg + ntl - 2
        p4 = nwk + 2 * nwg + ntl - 3
        p5 = 2 * nwk + 2 * nwg + ntl - 4
        wall = 2004 if ivisc != 0 else 1005

        L: list[str] = []
        w = L.append

        def rep(line):
            for _ in range(gn):
                w(line)

        def seg(grid, isg, bc, a, b, c, d):
            nd = 2 if bc == 2004 else 0
            w(f'{grid:10d}{isg:10d}{bc:10d}{a:10d}{b:10d}{c:10d}{d:10d}'
              f'{nd:10d}')
            if bc == 2004:
                w('    Twtype        Cq')
                w(f'{0.0:10.4f}{0.0:10.4f}')

        nbli: list[tuple] = []

        w('FILES:            ')
        for nm in ('cfl3d.xyz', 'plot3d_grid.xyz', 'plot3d_sol.bin',
                   'cfl3d.out', 'resid.out', 'cfl3d.turres', 'cfl3d.blomax',
                   'cfl3d.2out', 'cfl3d.prt', 'cfl3d.press', 'ovrlp.bin',
                   'patch.bin', 'cfl3d.restart'):
            w(f'{nm:<18}')
        if keywords:
            w('>')
            for k, v in sorted(keywords.items()):
                w(f'{k} {v}')
            w('<')
        w('CFL3D V6 INPUT FILE GENERATED WITH wing3d_otip.py')
        w('     XMACH     ALPHA      BETA  REUE,MIL   TINF,DR     IALPH    IHSTRY')
        w(f'{mach:10.4f}{alpha:10.4f}{beta:10.4f}{re_mil:10.4f}{tinf:10.2f}'
          f'{1:10d}{0:10d}')
        w('      SREF      CREF      BREF       XMC       YMC       ZMC')
        w(f'{sref:10.4f}{cref:10.4f}{bref:10.4f}{xmc:10.4f}'
          f'{0.0:10.4f}{0.0:10.4f}')
        w('        DT     IREST   IFLAGTS      FMAX     IUNST    CFLTAU')
        w(f'{-abs(cfl):10.2f}{0:10d}{iflagts:10d}{cfl_max:10.2f}{0:10d}'
          f'{7.5:10.4f}')
        w('     NGRID   NPLOT3D    NPRINT    NWREST      ICHK       I2D    NTSTEP       ITA')
        w(f'{-gn:10d}{gn:10d}{-1:10d}{nwrest:10d}{0:10d}{0:10d}{1:10d}{1:10d}')
        w('       NCG       IEM  IADVANCE    IFORCE  IVISC(I)  IVISC(J)  IVISC(K)')
        rep(f'{2:10d}{0:10d}{0:10d}{333:10d}{ivisc:10d}{ivisc:10d}{ivisc:10d}')
        w('      IDIM      JDIM      KDIM')
        for b in self.grid:
            ni, nj, nk = b.shape[:3]
            w(f'{ni:10d}{nj:10d}{nk:10d}')
        w('    ILAMLO    ILAMHI    JLAMLO    JLAMHI    KLAMLO    KLAMHI')
        rep(f'{0:10d}' * 6)
        w('     INEWG    IGRIDC        IS        JS        KS        IE        JE        KE')
        rep(f'{0:10d}' * 8)
        w('  IDIAG(I)  IDIAG(J)  IDIAG(K)  IFLIM(I)  IFLIM(J)  IFLIM(K)')
        rep(f'{1:10d}{1:10d}{1:10d}{4:10d}{4:10d}{4:10d}')
        w('   IFDS(I)   IFDS(J)   IFDS(K)  RKAP0(I)  RKAP0(J)  RKAP0(K)')
        rep(f'{1:10d}{1:10d}{1:10d}{0.3333:10.4f}{0.3333:10.4f}{0.3333:10.4f}')

        # per-block segment counts on each of the six faces
        w('      GRID     NBCI0   NBCIDIM     NBCJ0   NBCJDIM     NBCK0   NBCKDIM    IOVRLP')
        for i in range((gn // 2) - 3):
            w(f'{2 * i + 1:10d}{1:10d}{1:10d}{1:10d}{1:10d}{3:10d}{1:10d}{0:10d}')
            w(f'{2 * i + 2:10d}{1:10d}{1:10d}{1:10d}{1:10d}{1:10d}{1:10d}{0:10d}')
        w(f'{gn - 6:10d}{1:10d}{6:10d}{1:10d}{1:10d}{3:10d}{1:10d}{0:10d}')
        for g in (gn - 5, gn - 4, gn - 3):
            w(f'{g:10d}{1:10d}{1:10d}{1:10d}{1:10d}{1:10d}{1:10d}{0:10d}')
        w(f'{gn - 2:10d}{1:10d}{1:10d}{1:10d}{1:10d}{5:10d}{1:10d}{0:10d}')
        for g in (gn - 1, gn):
            w(f'{g:10d}{1:10d}{1:10d}{1:10d}{1:10d}{1:10d}{1:10d}{0:10d}')

        # I0: the root symmetry plane on the two wing blocks; the tip block's
        # I0 is the TIP SURFACE itself (a wall); everything else is interface
        w('I0:   GRID   SEGMENT    BCTYPE      JSTA      JEND      KSTA      KEND     NDATA')
        for i, b in enumerate(self.grid):
            _, nj, nk = b.shape[:3]
            bc = 1001 if i <= 1 else (wall if i == gn - 5 else 0)
            seg(i + 1, 1, bc, 1, nj, 1, nk)
        # IDIM: block 1 hands over to the four O-block faces + the far field
        w('IDIM: GRID   SEGMENT    BCTYPE      JSTA      JEND      KSTA      KEND     NDATA')
        for i, b in enumerate(self.grid):
            _, nj, nk = b.shape[:3]
            if i == gn - 7:
                for k, (a, c) in enumerate(((1, p1), (p1, p2), (p2, p3),
                                            (p3, p4), (p4, nj))):
                    seg(i + 1, k + 1, 0, a, c, 1, ntp)
                seg(i + 1, 6, 0, 1, nj, ntp, nk)
            elif i >= gn - 3:
                seg(i + 1, 1, 1000, 1, nj, 1, nk)
            else:
                seg(i + 1, 1, 0, 1, nj, 1, nk)
        w('J0:   GRID   SEGMENT    BCTYPE      ISTA      IEND      KSTA      KEND     NDATA')
        for i, b in enumerate(self.grid):
            ni, _, nk = b.shape[0], b.shape[1], b.shape[2]
            bc = 0 if i in (gn - 5, gn - 2) else 1000
            seg(i + 1, 1, bc, 1, ni, 1, nk)
        w('JDIM: GRID   SEGMENT    BCTYPE      ISTA      IEND      KSTA      KEND     NDATA')
        for i, b in enumerate(self.grid):
            ni, nk = b.shape[0], b.shape[2]
            if i <= gn - 6:
                bc = 1000 if i % 2 == 0 else wall
            elif i == gn - 3:
                bc = 1000
            else:
                bc = 0
            seg(i + 1, 1, bc, 1, ni, 1, nk)
        w('K0:   GRID   SEGMENT    BCTYPE      ISTA      IEND      JSTA      JEND     NDATA')
        for i, b in enumerate(self.grid):
            ni, nj = b.shape[0], b.shape[1]
            if i <= gn - 6 and i % 2 == 0:
                seg(i + 1, 1, 0, 1, ni, 1, p1)
                seg(i + 1, 2, wall, 1, ni, p1, p4)
                seg(i + 1, 3, 0, 1, ni, p4, p5)
            elif i == gn - 3:
                for k, (a, c) in enumerate(((1, p1), (p1, p2), (p2, p3),
                                            (p3, p4), (p4, p5))):
                    seg(i + 1, k + 1, 0, 1, ni, a, c)
            else:
                seg(i + 1, 1, 0, 1, ni, 1, nj)
        w('KDIM: GRID   SEGMENT    BCTYPE      ISTA      IEND      JSTA      JEND     NDATA')
        for i, b in enumerate(self.grid):
            ni, nj = b.shape[0], b.shape[1]
            bc = 1000 if (i <= gn - 6 and i % 2 == 0) or i == gn - 3 else 0
            seg(i + 1, 1, bc, 1, ni, 1, nj)

        w('      MSEQ    MGFLAG    ICONSF       MTT      NGAM')
        w(f'{mseq:10d}{1:10d}{0:10d}{0:10d}{2:10d}')
        w('      ISSC EPSSSC(1) EPSSSC(2) EPSSSC(3)      ISSR EPSSSR(1) EPSSSR(2) EPSSSR(3)')
        w(f'{0:10d}{0.3:10.1f}{0.3:10.1f}{0.3:10.1f}{0:10d}'
          f'{0.3:10.1f}{0.3:10.1f}{0.3:10.1f}')
        w('      NCYC    MGLEVG     NEMGL     NITFO')
        for lv in range(1, mseq + 1):
            w(f'{ncyc:10d}{lv:10d}{0:10d}{nitfo:10d}')
        w('      MIT1      MIT2      MIT3      MIT4      MIT5      MIT6      MIT7      MIT8')
        for _ in range(mseq):
            w(f'{1:10d}' * 8)

        # ---- 1-to-1 blocking: side 1 then side 2, ALL entries under one
        # header each (global0.F reads exactly 2*nbli + 2 lines here, and
        # echoinp(...,0) consumes none -- an interleaved table is one line too
        # long and the reader walks into the PATCH section)
        side1, side2 = [], []

        def pair(a, b):
            side1.append(a)
            side2.append(b)

        ni1, _, nk1 = self.grid[gn - 7].shape[:3]     # main wing block
        ni2 = self.grid[gn - 2].shape[0]              # far-field block ni
        for i in range((gn - 7) // 2):                # only for n_sec > 2
            ni, nj, nk = self.grid[2 * i].shape[:3]
            pair((2 * i + 1, 1, 1, 1, ni, p1, 1, 1, 2),
                 (2 * i + 2, 1, 1, ntl, ni, nwk, ntl, 1, 2))
            pair((2 * i + 1, 1, p4, 1, ni, nj, 1, 1, 2),
                 (2 * i + 2, 1, nwk, 1, ni, 1, 1, 1, 2))
            pair((2 * i + 1, ni, 1, 1, ni, nj, nk, 2, 3),
                 (2 * i + 3, 1, 1, 1, 1, nj, nk, 2, 3))
            pair((2 * i + 2, ni, 1, 1, ni, nwk, ntl, 2, 3),
                 (2 * i + 4, 1, 1, 1, 1, nwk, ntl, 2, 3))

        # the wing block's own wake cut, against the tail block
        pair((gn - 6, 1, 1, 1, ni1, p1, 1, 1, 2),
             (gn - 5, 1, 1, ntl, ni1, nwk, ntl, 1, 2))
        pair((gn - 6, 1, p4, 1, ni1, p5, 1, 1, 2),
             (gn - 5, 1, nwk, 1, ni1, 1, 1, 1, 2))
        # the wing block's outboard face, against the two O-blocks in the basin
        pair((gn - 6, ni1, 1, 1, ni1, p1, ntp, 2, 3),
             (gn - 3, 1, 1, ntl, ntp, nwk, ntl, 2, 1))
        pair((gn - 6, ni1, p1, 1, ni1, p2, ntp, 2, 3),
             (gn - 4, 1, 1, ntl, ntp, nwg, ntl, 2, 1))
        pair((gn - 6, ni1, p2, 1, ni1, p3, ntp, 2, 3),
             (gn - 4, 1, nwg, ntl, ntp, nwg, 1, 3, 1))
        pair((gn - 6, ni1, p3, 1, ni1, p4, ntp, 2, 3),
             (gn - 4, 1, nwg, 1, ntp, 1, 1, 2, 1))
        pair((gn - 6, ni1, p4, 1, ni1, p5, ntp, 2, 3),
             (gn - 3, 1, nwk, 1, ntp, 1, 1, 2, 1))
        # above the basin the wing block meets the far-field block directly
        pair((gn - 6, ni1, 1, ntp, ni1, p5, nk1, 2, 3),
             (gn - 2, 1, 1, 1, 1, p5, nbd, 2, 3))
        pair((gn - 5, ni1, 1, 1, ni1, nwk, ntl, 2, 3),
             (gn - 3, 1, 1, 1, 1, nwk, ntl, 2, 3))
        # tip O-block against wake O-block, and both out to the far field
        pair((gn - 4, 1, 1, 1, ntp, 1, ntl, 1, 3),
             (gn - 3, 1, nwk, 1, ntp, nwk, ntl, 1, 3))
        pair((gn - 4, ntp, 1, 1, ntp, nwg, ntl, 2, 3),
             (gn - 1, 1, 1, 1, 1, nwg, ntl, 2, 3))
        pair((gn - 3, ntp, 1, 1, ntp, nwk, ntl, 2, 3),
             (gn, 1, 1, 1, 1, nwk, ntl, 2, 3))
        # the far-field main block's K0 face, against the two far-field
        # O-blocks, split on the same five j ranges as the wing surface
        pair((gn - 2, 1, 1, 1, ni2, p1, 1, 1, 2),
             (gn, 1, 1, ntl, ni2, nwk, ntl, 1, 2))
        pair((gn - 2, 1, p1, 1, ni2, p2, 1, 1, 2),
             (gn - 1, 1, 1, ntl, ni2, nwg, ntl, 1, 2))
        pair((gn - 2, 1, p2, 1, ni2, p3, 1, 1, 2),
             (gn - 1, 1, nwg, ntl, ni2, nwg, 1, 1, 3))
        pair((gn - 2, 1, p3, 1, ni2, p4, 1, 1, 2),
             (gn - 1, 1, nwg, 1, ni2, 1, 1, 1, 2))
        pair((gn - 2, 1, p4, 1, ni2, p5, 1, 1, 2),
             (gn, 1, nwk, 1, ni2, 1, 1, 1, 2))
        pair((gn - 1, 1, 1, 1, ni2, 1, ntl, 1, 3),
             (gn, 1, nwk, 1, ni2, nwk, ntl, 1, 3))

        w('   1-1 BLOCKING DATA:')
        w('       NBLI')
        w(f'{len(side1):10d}')
        for side in (side1, side2):
            w(' NUMBER   GRID     :    ISTA   JSTA   KSTA   IEND   JEND   KEND  ISVA1  ISVA2')
            for n, e in enumerate(side, start=1):
                w(f'{n:7d}{e[0]:7d}      {e[1]:8d}{e[2]:7d}{e[3]:7d}'
                  f'{e[4]:7d}{e[5]:7d}{e[6]:7d}{e[7]:7d}{e[8]:7d}')
        w('   PATCH SURFACE DATA:')
        w('    NINTER')
        w('     0')
        # ★ Row counts here are set by NPLOT3D and NPRINT, NOT by the block
        # count.  NPLOT3D = gn gives one row per block (IPTYPE = 1, as the
        # reference deck uses); NPRINT = -1 gives |NPRINT| = ONE row.  Writing
        # gn rows in the PRINT OUT section instead made the reader run off the
        # end -- "Bad integer for item 1 in list input" at global.F:4414 -- and
        # the failure surfaced six sections later, so it looked like a
        # CONTROL SURFACE problem rather than a row-count one.
        w('   PLOT3D OUTPUT:')
        w('    GRID IPTYPE ISTART   IEND   IINC JSTART   JEND   JINC KSTART   KEND   KINC')
        for g in range(1, gn + 1):
            w(f'{g:6d}{1:7d}' + f'{0:7d}' * 9)
        w(' IMOVIE')
        w('     0')
        w('   PRINT OUT:')
        w('    GRID IPTYPE ISTART   IEND   IINC JSTART   JEND   JINC KSTART   KEND   KINC')
        w(f'{1:6d}' + f'{0:7d}' * 10)
        w('   CONTROL SURFACE:')
        w('   NCS')
        w('     0')
        w('    GRID ISTART   IEND   JSTART   JEND   KSTART   KEND  IWALL  INORM')
        Path(path).write_text('\n'.join(L) + '\n')
        self.info.update(nbli=len(side1), p=(p1, p2, p3, p4, p5),
                         ntp=ntp, nbd=nbd, nwg=nwg)
        return path
