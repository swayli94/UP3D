"""
Shared machinery for the CFL3D reference datasets (D05/D06/D08/D09 …).

This module builds the grid, writes the deck, runs the solver and parses the
output.  The per-dataset case lists and the CSV writing live in
``generate_cfl3d_reference.py``; everything here is case-agnostic.

External tool (NOT a repo dependency, both gitignored under ``tools/``):
  - the solver binary ``tools/cfl3d_seq`` (Aerolab CFL3D 6.7);
  - nothing else -- the grid generator ``cgrid_gmsh.py`` is VENDORED next to
    this file (see README.md for its source repo, commit and md5).

Geometry provenance -- the one thing that must NOT come from the CFL3D side:
the airfoil coordinates are written from **our own** generators
(``pyfp3d.meshgen.planar``), so the reference and the pyFP3D solution stand on
a bit-identical section.  Taking the CFL3D repo's own ``NACA0012.dat`` /
``RAE2822.dat`` would put a geometry difference inside every comparison --
the "are the two numbers the same thing?" family of error this project has
logged six times.
"""

from __future__ import annotations

import os
import re
import struct
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))          # the vendored cgrid_gmsh
sys.path.insert(0, str(REPO_ROOT))     # pyfp3d

from cgrid_gmsh import CGridGmsh, read_airfoil           # noqa: E402

TITLE_LINE = 'CFL3D V6 INPUT FILE GENERATED WITH cgrid_gmsh.py'
ILAM_HEADER = '    ILAMLO    ILAMHI    JLAMLO    JLAMHI    KLAMLO    ILAMHI'


# ---------------------------------------------------------------------------
#  the external binary
# ---------------------------------------------------------------------------

def find_solver(name: str = 'cfl3d_seq') -> Path:
    """Locate the gitignored CFL3D binary under tools/."""
    cand = REPO_ROOT / 'tools' / name
    if cand.is_file() and os.access(cand, os.X_OK):
        return cand
    raise FileNotFoundError(
        f'{cand} not found or not executable. The CFL3D binaries are an '
        'external tool, not a repo dependency: build them from the Aerolab '
        'CFL3D 6.7 tree (build/makefile_linux_gfortran, target cfl3d_seq) and '
        'drop the executable at tools/cfl3d_seq. See README.md here.')


# ---------------------------------------------------------------------------
#  geometry -- written from OUR generators, never from the CFL3D repo
# ---------------------------------------------------------------------------

#: Source contour density.  cgrid_gmsh re-splines the contour in arc length, so
#: the source must be finer than the grid spacing it will ask for at the nose
#: (~9e-6 chord).  n_half = 1001 gives a first LE interval of 2.5e-6.
N_HALF_SOURCE = 1001

RAE2822_ORDINATES = REPO_ROOT / 'cases' / 'meshes' / 'rae2822_2.5d' / 'rae2822.dat'


def write_geometry(kind: str, path: Path) -> Path:
    """Write a Selig-layout coordinate file for ``kind`` from our own library."""
    from pyfp3d.meshgen.planar import (
        load_airfoil_ordinates,
        naca0012_coordinates,
        pointset_airfoil_coordinates,
    )

    if kind == 'naca0012':
        coords = naca0012_coordinates(n_half=N_HALF_SOURCE)
        title = ('NACA0012 (closed-TE coefficient set, -0.1036) from '
                 'pyfp3d.meshgen.planar.naca0012_coordinates')
    elif kind == 'rae2822':
        x, z_lo, z_up = load_airfoil_ordinates(RAE2822_ORDINATES)
        # Cook Table 6.1 tabulates the lower ordinate positive-DOWN; the
        # physical lower z is its negative (see load_airfoil_ordinates).
        coords = pointset_airfoil_coordinates(x, -z_lo, z_up,
                                              n_half=N_HALF_SOURCE)
        title = ('RAE2822 (Cook Table 6.1 ordinates, PCHIP) from '
                 'pyfp3d.meshgen.planar.pointset_airfoil_coordinates')
    else:
        raise ValueError(f'unknown geometry {kind!r}')

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as f:
        f.write(title + '\n')
        for xx, yy in coords:
            f.write(f'{xx:20.12f}{yy:20.12f}\n')
    return path


# ---------------------------------------------------------------------------
#  mesh ladder
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Level:
    """One rung of the grid-convergence ladder."""
    name: str
    n_foil: int
    n_wake: int
    n_grow: int
    h1: float | None = None        # Euler: first cell height, chords
    y_plus: float | None = None    # RANS: target y+ (sized per case Re)


#: ONE tangential/wake ladder for both equation sets, refining by ~sqrt(2) per
#: rung in every direction.  The Euler and RANS ladders differ ONLY in the
#: wall-normal direction, so an Euler-vs-RANS read at the same condition
#: differs in the wall treatment and not in the surface discretisation -- which
#: is what makes the FP-Euler / (FP+IBL)-RANS error decomposition clean.
#:
#: ★ n_wake = 41/57/81 and not 29/41/57, and that is a MEASUREMENT: the wake
#: cut reaches 80 chords geometrically from the trailing-edge spacing, so 29
#: points make its growth ratio large enough to destabilise the transonic RANS
#: cold start.  Isolated at L1 on M 0.778 / alpha 2.03 / Re 6e6 / SST, one
#: variable at a time:
#:
#:     n_wake 29, CFL 1                       diverged  NaN block 3 cycle 8
#:     n_wake 29, CFL 1, NITFO 500            diverged  NaN block 3 cycle 10
#:     n_wake 29, CFL 0.5 ramp 10 NITFO 500   ok        cl 0.333826
#:     n_wake 41, CFL 1                       ok        cl 0.332273
#:     n_wake 29, CFL 1, UPSTREAM geometry    diverged  NaN block 3 cycle 17
#:
#: The last leg is the one that matters for provenance: swapping in the CFL3D
#: repo's own NACA0012 ordinates changes nothing, so the failure is the wake
#: grading and NOT our section.
TANGENTIAL = ((101, 41), (141, 57), (201, 81))       # (n_foil, n_wake)

#: Euler ladder: normal count and first cell height both refine by ~sqrt(2).
EULER_LEVELS = (
    Level('L1', n_foil=101, n_wake=41, n_grow=49, h1=2.000e-3),
    Level('L2', n_foil=141, n_wake=57, n_grow=69, h1=1.414e-3),
    Level('L3', n_foil=201, n_wake=81, n_grow=97, h1=1.000e-3),
)

#: RANS ladder.  y+ is held at 1 on every rung -- a RANS result whose boundary
#: layer is not resolved is not a reference (data request 4.7 item 3) -- so the
#: normal ladder refines the OUTER boundary layer and wake, and the wall-normal
#: growth ratio comes down ~1.19 -> 1.14 -> 1.11.  Starting the normal count as
#: low as the Euler ladder's 49 would put that ratio at ~1.42, which is not a
#: RANS grid.
RANS_LEVELS = (
    Level('L1', n_foil=101, n_wake=41, n_grow=97, y_plus=1.0),
    Level('L2', n_foil=141, n_wake=57, n_grow=129, y_plus=1.0),
    Level('L3', n_foil=201, n_wake=81, n_grow=161, y_plus=1.0),
)


# ---------------------------------------------------------------------------
#  case description
# ---------------------------------------------------------------------------

@dataclass
class Case:
    tag: str
    geometry: str                  # 'naca0012' | 'rae2822'
    mach: float
    alpha: float                   # degrees, EXPERIMENTAL alpha, uncorrected
    re_mil: float                  # chord Reynolds number / 1e6
    ivisc: int                     # 0 = Euler, 6 = SA, 7 = k-omega SST
    x_tr: float | None = None      # forced transition x/c, None = fully turbulent
    tinf: float = 460.0
    request: str = ''             # data-request id (2D-1 ... R-7)
    note: str = ''
    keywords: dict = field(default_factory=dict)   # CFL3D keyword block
    solver: dict = field(default_factory=dict)   # NCYC/CFL_START/... overrides

    @property
    def model(self) -> str:
        return 'euler' if self.ivisc == 0 else 'rans'

    @property
    def turb_model(self) -> str:
        return {0: 'none', 4: 'baldwin-lomax', 5: 'k-omega',
                6: 'sa', 7: 'sst'}.get(self.ivisc, f'ivisc{self.ivisc}')


#: Startup recipes.  ★ A startup recipe is a CALIBRATION of the solver on
#: THESE grids, not a statement about the flow -- same status as this project's
#: EW forcing, taper r_c and descent10 threshold -- so it is chosen by
#: measurement and quoted with what it was measured on.
#:
#: Euler: cold start at CFL 1 with no first-order phase reaches 7.7-9.9 decades
#: of residual drop on all 24 runs, so it is left alone.
EULER_SOLVER = dict(NCYC=1000, CFL_START=1.0, CFL_RAMP=5.0,
                    IFLAGTS=500, NITFO=0, MSEQ=3)

#: RANS: measured at L1 on M 0.778 / alpha 2.03 / Re 6e6 / SST (the case that
#: failed first), one variable at a time, on the 41-point wake:
#:
#:     CFL 1, NITFO 0,   NCYC 2000  ok        cl 0.332329  drop 2.99  121.7 s
#:     CFL 1, NITFO 500, NCYC 2000  ok        cl 0.332325  drop 3.42   68.5 s
#:     CFL 0.5 ramp 10, NITFO 500   diverged  NaN block 3 cycle 24
#:
#: ★ The middle row is adopted: it is FASTER and better converged than the
#: plain cold start, and it agrees with it to FIVE digits in cl -- which is the
#: property that matters for a reference, because it says the answer does not
#: depend on the startup.  ★ And the third row is why the recipe was measured
#: rather than reasoned: ramping to CFL 5 diverges, so "gentler start" is not
#: monotone in the knobs and a plausible-looking choice would have been wrong.
RANS_SOLVER = dict(NCYC=2000, CFL_START=1.0, CFL_RAMP=5.0,
                   IFLAGTS=500, NITFO=500, MSEQ=3)

#: ★★ Fallback for a diverged RANS leg, retried ONCE -- the same shape as this
#: project's own ``_SEED_FALLBACK``: it fires only where the first attempt
#: already failed, so it cannot move a result a caller would have used, and the
#: recipe actually used is recorded per row in the ``recipe`` column.
#:
#: Terminal CFL is 1.0 here (0.05 x 20), i.e. it ramps to a LOWER terminal CFL
#: than the default's 5.0 and starts 20x below it, with a long first-order
#: phase.  Measured on the two failing legs:
#:
#:     CFL_START 0.1  x10 (term 1.0)  diverged  cycle 19  (both legs)
#:     CFL_START 0.05 x20 (term 1.0)  OK        both legs
#:     CFL_START 0.02 x50 (term 1.0)  OK        both legs
#:
#: ★★★ But read what it buys with the convergence caliber attached, because the
#: two recipes that work do NOT always agree:
#:
#:     M0.5 xtr005 L2:  0.05 -> cl 0.248883 | 0.02 -> cl 0.249002   (0.05 %)
#:     M0.778 L3:       0.05 -> cl 0.373755 | 0.02 -> cl 0.356217   (4.9 %)
#:
#: and the residual drop is ~1.1-1.6 decades in all four.  So the L2 row is
#: startup-INDEPENDENT and usable, while the M0.778 L3 row is startup-DEPENDENT
#: at 4.9 % and is NOT a converged solution.  ⇒ a row rescued by this fallback
#: must be read with its ``resid_decades``; the README states the floor.
RANS_FALLBACK = dict(NCYC=2000, CFL_START=0.05, CFL_RAMP=20.0,
                     IFLAGTS=500, NITFO=1500, MSEQ=3)


def default_solver(case: 'Case') -> dict:
    return dict(EULER_SOLVER if case.ivisc == 0 else RANS_SOLVER)


#: ★★★ Eddy-viscosity cap, on EVERY RANS row.  CFL3D's default is
#: ``edvislim = 1e10``, i.e. effectively unbounded.
#:
#: Why it is here.  Measured by factorial on M0.778 / alpha 2.03 / Re6e6, one
#: count at a time, with everything else fixed:
#:
#:     (n_foil 101, n_wake 41, n_grow  97) SST  281x 97  ok, 3.42 dec
#:     (n_foil 201, n_wake 41, n_grow  97) SST  481x 97  ok, 0.79 dec
#:     (n_foil 101, n_wake 81, n_grow  97) SST  361x 97  ok, 3.39 dec
#:     (n_foil 101, n_wake 41, n_grow 161) SST  281x161  DIVERGED cycle 9
#:     (n_foil 201, n_wake 81, n_grow 161) SA   561x161  ok, 3.66 dec
#:     (n_foil 201, n_wake 81, n_grow 161) SST  561x161  DIVERGED cycle 11
#:
#: ⇒ the trigger is the WALL-NORMAL COUNT and it is SST-specific: SA converges
#: on the identical 561x161 grid.  With ``edvislim 1.e05`` that same SST/L3 case
#: converges to |R| 2.54e-09 (3.44 decades) and continues its own ladder:
#: cl 0.332325 (L1) -> 0.339302 (L2) -> 0.342969 (L3), monotone, ~3.4 decades
#: at every rung.
#:
#: ★★ It is applied at EVERY rung, not just the one that needed it: a limiter
#: switched on for one level only would make the three rungs THREE DIFFERENT
#: MODELS, which is the cross-provenance error this dataset exists to avoid.
#: That obliges measuring what it does where the un-capped run already
#: converged, so it was A/B'd (M0.778, L1, one variable):
#:
#:     SST  default(1e10)  cl 0.332325  cd 0.025884  cdv 0.005371  |R| 6.714991e-09
#:     SST  edvislim 1.e05 cl 0.332325  cd 0.025884  cdv 0.005371  |R| 6.649566e-09
#:     SA   default(1e10)  cl 0.413929  cd 0.030476  cdv 0.006301  |R| 1.011592e-08
#:     SA   edvislim 1.e05 cl 0.413929  cd 0.030476  cdv 0.006301  |R| 1.011592e-08
#:
#: ★ SST agrees to six digits in cl, cd AND cdv while |R| moves in the third
#: digit -- so the cap DID act during the transient and left the converged
#: answer alone, which is what makes it admissible.  ★ SA is BIT-identical
#: including |R|, and that is structural rather than lucky: ``edvislim`` appears
#: in twoeqn.F / threeeqn.F / foureqn.F and NOT in spalart.F, so it cannot
#: touch the SA arm at all.
#:
#: ★ Compare the alternative that was rejected: the gentler-startup rescue also
#: gets M0.778/L3 to run, but only to ~1.1 decades and with a 4.9 % cl
#: difference between two startup recipes -- a snapshot, not a solution.  See
#: RANS_FALLBACK.
RANS_KEYWORDS = {'edvislim': '1.e05'}


def default_keywords(case: 'Case') -> dict:
    return {} if case.ivisc == 0 else dict(RANS_KEYWORDS)


# ---------------------------------------------------------------------------
#  build
# ---------------------------------------------------------------------------

def _laminar_j_range(xy, jte1, jte2, x_tr):
    """1-based j indices bracketing the wall stations with x/c <= x_tr.

    The C-grid's j runs lower far field -> lower TE (jte1) -> LE -> upper TE
    (jte2) -> upper far field, so "x/c <= x_tr on BOTH surfaces" is ONE
    contiguous j interval straddling the leading edge.  CFL3D's laminar test is
    ``j >= jlamlo .and. j < jlamhi`` on cell indices and it halves the pair for
    every coarser mesh-sequence level itself (global.F:1361), so the card is
    written for the finest level.
    """
    j = np.arange(jte1, jte2 + 1)                       # 1-based, wall segment
    x = xy[jte1 - 1:jte2, 0, 0]                         # k = 1 wall line
    inside = x <= x_tr
    if not inside.any():
        raise ValueError(f'x_tr = {x_tr} lands ahead of the first wall station')
    return int(j[inside].min()), int(j[inside].max()) + 1


def build_case(case: Case, level: Level, workdir: Path, geom_file: Path,
               verbose: bool = False, recipe: dict | None = None) -> dict:
    """Write cfl3d.xyz + cfl3d.inp into ``workdir``; return the grid record.

    ``recipe`` overrides the per-equation-set default (used by the diverged-leg
    fallback); ``case.solver`` still wins over it, so a per-case override such
    as the near-stall startup is never silently discarded.
    """
    s = dict(recipe if recipe is not None else default_solver(case),
             **case.solver)

    if case.ivisc == 0:
        wall = dict(h1=level.h1, re_chord=case.re_mil * 1e6)
    else:
        wall = dict(y_plus=level.y_plus, re_chord=case.re_mil * 1e6)

    cg = CGridGmsh(n_foil=level.n_foil, n_wake=level.n_wake,
                   n_grow=level.n_grow, n_coarse=2, r_far=80.0,
                   le_factor=0.06, wake_ratio_max=1.5, **wall)
    x, y = read_airfoil(geom_file)
    cg.generate(x, y, verbose=verbose)

    workdir.mkdir(parents=True, exist_ok=True)
    cg.write(str(workdir),
             mach=case.mach, alpha=case.alpha, re_mil=case.re_mil,
             tinf=case.tinf, cfl=s['CFL_START'], cfl_max=s['CFL_RAMP'],
             iflagts=s['IFLAGTS'], ncyc=s['NCYC'], nitfo=s['NITFO'],
             mseq=s['MSEQ'], ivisc=case.ivisc)

    rec = dict(nj=cg.nj, nk=cg.nk, jte1=cg.jte1, jte2=cg.jte2,
               n_foil=cg.n_foil, n_wake=cg.n_wake, n_grow=cg.n_grow,
               jlamlo=0, jlamhi=0, x_tr_actual='')
    rec.update({k: v for k, v in cg.info.items()
                if k in ('h1', 'wall_growth_max', 'wake_ratio',
                         'orthogonality_min_deg', 'orthogonality_max_deg',
                         'ds_le', 'negative_volumes')})
    rec['h1'] = float(cg.info.get('h1', wall.get('h1') or 0.0))

    if case.x_tr is not None:
        jlo, jhi = _laminar_j_range(cg.xy, cg.jte1, cg.jte2, case.x_tr)
        _patch_laminar_card(workdir / 'cfl3d.inp', jlamlo=jlo, jlamhi=jhi,
                            kdim=cg.nk, idim=2)
        rec['jlamlo'], rec['jlamhi'] = jlo, jhi
        xw = cg.xy[jlo - 1:jhi, 0, 0]
        rec['x_tr_actual'] = f'{float(xw.max()):.6f}'
    if case.keywords:
        _insert_keywords(workdir / 'cfl3d.inp', case.keywords)
    rec['keywords'] = ' '.join(f'{k}={v}' for k, v in
                               sorted(case.keywords.items()))
    return rec


def _insert_keywords(path: Path, keywords: dict):
    """Prepend a CFL3D keyword block.

    ``readkey`` is called before the title line is read: it reads one line, and
    if it starts with ``>`` it consumes keyword lines until one starts with
    ``<``, otherwise it treats that line AS the title (readkey.F:1071-1093).
    So the block goes immediately before the title, i.e. after the FILES card.
    """
    lines = path.read_text().split('\n')
    try:
        t = lines.index(TITLE_LINE)
    except ValueError:
        raise RuntimeError(f'{path}: title line not found; cannot place the '
                           'keyword block')
    if lines[t - 1].startswith('<'):
        raise RuntimeError(f'{path}: a keyword block is already present')
    block = ['>'] + [f'{k} {v}' for k, v in sorted(keywords.items())] + ['<']
    lines[t:t] = block
    path.write_text('\n'.join(lines))


def _patch_laminar_card(path: Path, jlamlo: int, jlamhi: int, kdim: int, idim: int):
    """Fill in the laminar-region card: forced transition at a fixed x/c.

    The vendored generator always writes the card as six zeros (= no laminar
    region, fully turbulent).  The substitution is ASSERTED rather than
    attempted: a silently unpatched deck would run FULLY TURBULENT while the
    CSV claimed a transition location the solver never used.

    Inside the index box CFL3D zeroes the turbulence PRODUCTION term
    (``cutoff = 0``, twoeqn.F:2411 for k-omega/SST, spalart.F for SA), so no
    turbulence is generated ahead of the trip and mu_t stays at freestream
    level -- the standard CFL3D fixed-transition mechanism.

    ★ ``i_lam_forcezero 1`` -- which additionally forces ``vist3d = 0`` inside
    the box -- is NOT used, and that is a MEASUREMENT, not a preference.  With
    it on, every tripped case NaNs at cycle 4 on the coarsest mesh-sequence
    level; with it off and everything else identical the same case converges to
    |R| 4.6e-09.  Isolated one variable at a time at L1 (M 0.5, alpha 2,
    Re 3e6, SST):

        no trip                        ok        cl 0.248504  |R| 4.65e-09
        trip 0.05, no forcezero        ok        cl 0.250708  |R| 4.65e-09
        trip 0.05, forcezero, k box 40 diverged  NaN, block 3 cycle 4
        trip 0.05, forcezero, full box diverged  NaN, block 3 cycle 4

    The k extent of the box makes no difference, so it is the hard zeroing
    itself and not the region size.  The ``forcezero_echo`` column therefore
    reads 0 on every row: it records that the flag was not used.
    """
    lines = path.read_text().split('\n')

    try:
        i = lines.index(ILAM_HEADER)
    except ValueError:
        raise RuntimeError(f'{path}: laminar-region card header not found; '
                           'the vendored cgrid_gmsh deck layout has changed')
    if lines[i + 1].split() != ['0'] * 6:
        raise RuntimeError(f'{path}: laminar-region card is not six zeros '
                           f'({lines[i + 1]!r})')
    lines[i + 1] = (f'{1:10d}{idim:10d}{jlamlo:10d}{jlamhi:10d}'
                    f'{1:10d}{kdim:10d}')
    path.write_text('\n'.join(lines))


# ---------------------------------------------------------------------------
#  run
# ---------------------------------------------------------------------------

def run_case(workdir: Path, solver: Path, timeout: float = 7200.0) -> dict:
    """Run the solver in ``workdir``.  Returns status + wall time."""
    for junk in ('cfl3d.restart', 'cfl3d.error', 'cfl3d.out',
                 'clcd_total.dat', 'cfl3d.prt'):
        p = workdir / junk
        if p.exists():
            p.unlink()

    t0 = time.perf_counter()
    with open(workdir / 'run.log', 'wb') as log:
        try:
            rc = subprocess.call([str(solver)], cwd=str(workdir),
                                 stdout=log, stderr=subprocess.STDOUT,
                                 timeout=timeout)
        except subprocess.TimeoutExpired:
            return dict(status='timeout', wall_s=time.perf_counter() - t0)
    wall = time.perf_counter() - t0

    status = 'ok'
    err = workdir / 'cfl3d.error'
    if rc != 0:
        status = f'exit{rc}'
    elif not err.is_file():
        status = 'no-error-file'
    elif 'terminated normally' not in err.read_text(errors='replace'):
        status = 'diverged'
    return dict(status=status, wall_s=wall)


# ---------------------------------------------------------------------------
#  parse
# ---------------------------------------------------------------------------

FORCE_KEYS = ('CL', 'CD', 'CDp', 'CDv', 'CZ', 'CY', 'CX', 'wetted_area',
              'CMX', 'CMY', 'CMZ')


def read_forces(out: Path) -> dict:
    """Final all-blocks force/moment summary and the y+ statistics."""
    txt = out.read_text(errors='replace')
    blocks = txt.split('SUMMARY OF FORCES AND MOMENTS - ALL GLOBAL BLOCKS')
    res = {}
    if len(blocks) > 1:
        nums = re.findall(r'[-+]?\d\.\d+E[-+]\d+', blocks[-1][:600])
        for key, val in zip(FORCE_KEYS, nums):
            res[key] = float(val)

    m = re.search(r'Y\+ AVG\s+Y\+ STD DEV\s+NY\+ > 5\s+NPTS\s*\n\s*'
                  r'([-+.\dEe]+)\s+([-+.\dEe]+)', txt)
    if m:
        res['yplus_avg'] = float(m.group(1))
    m = re.search(r'Y\+ MAX\s+JLOC\s+ILOC\s+Y\+ MIN\s+JLOC\s+ILOC\s*\n\s*'
                  r'([-+.\dEe]+)\s+\d+\s+\d+\s+([-+.\dEe]+)', txt)
    if m:
        res['yplus_max'] = float(m.group(1))
        res['yplus_min'] = float(m.group(2))

    # The solver echoes the laminar region it actually used.  This is the only
    # independent confirmation that the patched card took effect, so it is read
    # back rather than trusted.
    m = re.search(r'laminar region is:\s*\n\s*i=\s*(\d+)\s+to\s*(\d+),\s*'
                  r'j=\s*(\d+)\s+to\s*(\d+),\s*k=\s*(\d+)\s+to\s*(\d+)', txt)
    res['laminar_echo'] = ('' if m else 'none')
    if m:
        res['laminar_echo'] = (f'i{m.group(1)}-{m.group(2)}_'
                               f'j{m.group(3)}-{m.group(4)}_'
                               f'k{m.group(5)}-{m.group(6)}')
    res['forcezero_echo'] = int('forcing vist3d=0' in txt)

    # The keyword block is echoed verbatim between the two banner lines, so
    # what the solver actually READ can be recovered rather than trusted.
    m = re.search(r'>-+ begin keyword-driven input section -+>\n(.*?)'
                  r'<-+ end keyword-driven input section', txt, re.S)
    res['keyword_echo'] = ('' if not m else
                           ' '.join(m.group(1).split()))
    return res


def read_history(path: Path):
    """clcd_total.dat: level, block, iter, res, tres, cl, cd, cdp, cdf.

    ★ The file's own header names NINE columns but every row carries TEN
    values -- there is one unnamed trailing column (side force).  Verified
    against the cfl3d.out summary on a converged Euler run: columns 5..8 are
    cl / cd / cdp / cdv and reproduce it exactly (cl 0.34823944,
    cd = cdp 0.022600953, cdv 0), so the named mapping holds and only the tail
    is extra.  Read as p[3:9] and do not widen it without re-checking.
    """
    rows = []
    for line in open(path, errors='replace'):
        p = line.split()
        if len(p) < 9 or not p[0].isdigit():
            continue
        try:
            rows.append([int(p[0]), int(p[1]), int(p[2])]
                        + [float(t) for t in p[3:9]])
        except ValueError:
            continue
    return np.array(rows) if rows else None


def convergence(hist) -> dict:
    """Residual caliber: final value and decades dropped on the finest level."""
    if hist is None:
        return dict(resid_final='', resid_decades='', ncyc_total='')
    lev = hist[:, 0]
    fine = hist[lev == lev.max()]
    r = fine[:, 3]
    r = r[r > 0]
    if r.size == 0:
        return dict(resid_final='', resid_decades='', ncyc_total=len(hist))
    return dict(resid_final=f'{r[-1]:.6e}',
                resid_decades=f'{np.log10(r[0] / r[-1]):.2f}',
                ncyc_total=int(len(hist)))


def read_wall_cp(prt: Path):
    """k = 1 wall line from cfl3d.prt, ordered lower TE -> LE -> upper TE.

    One row per cell centre: I J K X Y Z U V W P T MACH cp mut.  Rows outside
    0 <= x <= 1 are the wake cut and are dropped.
    """
    x, y, cp, mach = [], [], [], []
    with open(prt, 'r', errors='replace') as f:
        for line in f:
            p = line.split()
            if len(p) != 14:
                continue
            try:
                v = [float(t) for t in p]
            except ValueError:
                continue
            if int(v[2]) != 1:
                continue
            if -1e-9 <= v[3] <= 1.0 + 1e-9:
                x.append(v[3])
                y.append(v[4])
                cp.append(v[12])
                mach.append(v[11])
    if not x:
        raise RuntimeError(f'no wall data in {prt}')
    return (np.array(x), np.array(y), np.array(cp), np.array(mach))


def split_surfaces(x, y, cp, mach):
    """Split a wall line running lower TE -> LE -> upper TE."""
    ile = int(np.argmin(x))
    up = slice(ile, None)
    lo = slice(0, ile + 1)
    return (dict(x=x[up], y=y[up], cp=cp[up], mach=mach[up]),
            dict(x=x[lo], y=y[lo], cp=cp[lo], mach=mach[lo]))
