"""
Reference-data generator: CFL3D Euler on the ONERA M6 wing, for gate D07
(request items 3D-1 / 3D-2 of the phase-six data request).

    python generate_m6_reference.py --from-runs <ladder dir>

Reads a finished ladder produced by ``wing3d_otip.py`` + ``tools/cfl3d_seq``
and writes ``euler_onera_m6/``.  It performs no solve; the runs live in the
gitignored work directory and the CSVs here are the evidence.

★★★ Two calibration facts that any reader of these files needs, both of them
measured rather than assumed:

1. **The deck carries ``SREF = 1.00``** (as the reference implementation's deck
   does), so the CL/CD that CFL3D reports are NOT normalised by the wing area.
   This file publishes the wing coefficients, obtained by dividing by the
   exposed half-planform area in root-chord units,
   ``0.5 (1 + ct/cr) (b/cr) = 1.15932`` -- and keeps the raw values in their
   own columns so the conversion is auditable.  Unconverted, L3's CL reads
   0.3321 where the M6 inviscid value at this condition is ~0.28-0.30; the
   converted 0.2864 lands in that band.
2. **Force coefficients come from the ``cfl3d.out`` SUMMARY, never from
   ``clcd_total.dat``** -- on a multiblock grid the final iteration writes one
   row per block there, and the ``blk = 1`` row was 1.6 % off in cd.

★★ And the honest limitation, stated here because it decides what may be
gated: **cl is NOT in the asymptotic range on this ladder** (its rung-to-rung
differences still grow, ratio 1.74), while **cd is** (ratio 0.74).  So cd and
the shock positions carry an error bar and cl does not; cl is RECORDED.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import re

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))

from wing3d_otip import (                                   # noqa: E402
    M6_CHORD_MAC, M6_CHORD_ROOT, M6_CHORD_TIP, M6_RE_ROOT_CHORD, M6_SPAN,
    M6_SWEEP_LE_DEG, build_m6,
)
from pyfp3d.post.shock import cp_critical, shock_metrics    # noqa: E402

#: The seven measured stations of AGARD AR-138 TEST 2308, read from
#: cases/reference_data/onera_m6_experiment/experiment-Cp.dat.
EXP_STATIONS = (0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99)

#: Exposed half-planform area, in ROOT-CHORD units (the grid is normalised to
#: root chord = 1).  This is the divisor that turns CFL3D's SREF = 1 output
#: into wing coefficients.
S_HALF = 0.5 * (1.0 + M6_CHORD_TIP / M6_CHORD_ROOT) * (M6_SPAN / M6_CHORD_ROOT)

MACH, ALPHA = 0.8395, 3.06


# ---------------------------------------------------------------------------
#  readers
# ---------------------------------------------------------------------------

def deck_ncyc(d: Path):
    """-> (ncyc, mseq) from the deck.

    ★ mseq must come from the DECK, not from the highest level seen in the
    history file: a run that died before reaching the finest sequence level
    writes no rows for it, so a file-derived mseq makes the run look complete
    one level early.
    """
    t = (d / 'cfl3d.inp').read_text().splitlines()
    for i, l in enumerate(t):
        if 'NCYC' in l and 'NITFO' in l:
            rows = []
            for k in range(i + 1, len(t)):
                p = t[k].split()
                if len(p) != 4 or not p[0].lstrip('-').isdigit():
                    break
                rows.append(int(p[0]))
            return rows[0], len(rows)
    raise RuntimeError('no NCYC block')


def run_status(d: Path):
    """-> (ok: bool, detail: str, last_it: int, want_it: int)

    ★★ ``cfl3d.error`` is written on BOTH outcomes and is the authoritative
    field; the normal-termination banner appears nowhere else (not in
    cfl3d.out, not on stdout).  Success is `error code: 0`, failure is
    `error code: -1` followed by the message.  Three wrong versions of this
    function preceded the right one -- "NaN in cfl3d.out" (the message is not
    there), "a SUMMARY exists" (CFL3D writes one while dying), and "the error
    file is non-empty" (so is a successful run's).  Each was caught only by
    running it against cases whose answer was already known.
    """
    d = Path(d)
    rows = ([l.split() for l in (d / 'clcd_total.dat').read_text().splitlines()
             if l.split() and l.split()[0].isdigit()]
            if (d / 'clcd_total.dat').is_file() else [])
    last_it = max((int(r[2]) for r in rows), default=0)
    nlv = max((int(r[0]) for r in rows), default=0)
    ncyc, mseq = deck_ncyc(d)
    want = mseq * ncyc

    ef = d / 'cfl3d.error'
    if not (ef.is_file() and ef.stat().st_size):
        return False, f'STILL RUNNING at IT {last_it}/{want}', last_it, want
    txt = ef.read_text(errors='replace')
    m = re.search(r'error code:\s*(-?\d+)', txt)
    code = int(m.group(1)) if m else None
    if code is None:
        return False, 'cfl3d.error has no error code', last_it, want
    if code != 0:
        lines = [l.strip() for l in txt.splitlines() if l.strip()]
        return False, f'DIVERGED (code {code}): {lines[-1]}', last_it, want
    if want and last_it < want:
        return False, f'code 0 but stopped at IT {last_it}/{want}', last_it, want
    return True, f'ok, IT {last_it}/{want}', last_it, want



def read_summary(out: Path) -> dict:
    """cl / cd / cdp / cdv / wetted area from the cfl3d.out SUMMARY.

    ★ The authoritative source -- see the module docstring.
    """
    txt = out.read_text(errors='replace')
    blk = txt.split('SUMMARY OF FORCES AND MOMENTS - ALL GLOBAL BLOCKS')
    if len(blk) < 2:
        raise RuntimeError(f'{out}: no force summary (did the run finish?)')
    n = re.findall(r'[-+]?\d\.\d+E[-+]\d+', blk[-1][:600])
    keys = ('CL', 'CD', 'CDp', 'CDv', 'CZ', 'CY', 'CX', 'wetted_area',
            'CMX', 'CMY', 'CMZ')
    return {k: float(v) for k, v in zip(keys, n)}


def read_residual(hist: Path) -> dict:
    """Final residual and decades dropped on the finest mesh-sequence level.

    ★ Only column 3 is read -- the force columns of this file are not safe,
    see the module docstring.
    """
    rows = [l.split() for l in open(hist, errors='replace')
            if l.split() and l.split()[0].isdigit()]
    if not rows:
        return dict(resid_final='', resid_decades='', ncyc='')
    lv = max(int(r[0]) for r in rows)
    fine = [r for r in rows if int(r[0]) == lv]
    r = np.array([float(x[3]) for x in fine])
    r = r[r > 0]
    return dict(resid_final=f'{r[-1]:.6e}',
                resid_decades=f'{math.log10(r[0] / r[-1]):.2f}',
                ncyc=len(fine))


def yplus_patches(out: Path):
    """Every Y+ block in cfl3d.out, keyed by the wall patch it belongs to.

    ★★★ THERE ARE FOUR BLOCKS, NOT ONE, and reading only the first is a
    reporting defect that hides an unresolved patch.  The M6 deck marks three
    viscous walls, and CFL3D reports each separately plus an aggregate:

        JLOC/ILOC  the WING SURFACE      y+ 0.65-0.69, max 1.26, 0 above 5
        KLOC/ILOC  the BLUNT TE BASE     y+ 231-321,  ALL points above 5
        JLOC/KLOC  the TIP O-BLOCK       y+ 0.032-0.036, 0 above 5
        BLOCK/GRID the AGGREGATE         y+ 12.2-13.2

    ★★ Neither the first block nor the aggregate may be quoted as "the y+".
    The first hides the TE base; the aggregate is a POINT average over patches
    whose spacing differs by four orders of magnitude, and it is dominated by a
    patch carrying 0.071 % of the wetted area (3.8-5.1 % of the wall points).

    ★ Why the TE base is coarse and why it does not contaminate cd_friction:
    its wall-normal direction is the TE-THICKNESS direction, whose spacing is
    set by ``n_tail`` and not by the y+ clustering.  Measured on L1: it carries
    Cdv = 0.0000000 EXACTLY -- its face normal is streamwise, so wall shear
    lies in the face and contributes nothing to Cd -- but Cdp = -0.0004556,
    which is **-3.0 % of the total pressure drag**.  That share is the exposed
    quantity.

    ★★ And it does NOT shrink with refinement, which is the more useful fact.
    Measured across the ladder as n_tail goes 9 -> 17, the base's y+ drops
    303 -> 231 (a factor 1.3) while its cdp SHARE stays flat: -3.027 % ->
    -3.057 % (SST), -3.229 % -> -3.274 % (SA).  ⇒ refining n_tail improves the
    base's wall resolution without moving its drag share, because the base
    pressure is set by the separated-flow topology rather than by the wall
    layer.  A first draft of this docstring said the share shrinks 1.3x per
    step -- that conflated the y+ with the share.
    """
    txt = out.read_text(errors='replace')
    pats = {}
    # ★ the header is `Y+ MAX <i> <j>   Y+ MIN <i> <j>` -- six tokens, not two.
    #   An earlier pattern required the line to END after the two index labels
    #   and so matched nothing, returning None for every patch SILENTLY.
    for m in re.finditer(r'Y\+ MAX\s+(\w+)\s+(\w+)\s+Y\+ MIN[^\n]*\n'
                         r'\s*([-\d.E+]+)'
                         r'(.*?)Y\+ AVG\s+Y\+ STD DEV\s+NY\+ > 5\s+NPTS'
                         r'\s*\n\s*([-\d.E+]+)\s+([-\d.E+]+)\s+(\d+)'
                         r'\s+(\d+)', txt, re.S):
        a, b, mx, _mid, avg, _sd, ngt5, npts = m.groups()
        pats[f'{a}/{b}'] = dict(avg=float(avg), max=float(mx),
                                ngt5=int(ngt5), npts=int(npts))
    agg = re.search(r'DN MAX\s+ILOC\s+JLOC\s+KLOC\s+BLOCK\s+GRID.*?'
                    r'Y\+ AVG\s+Y\+ STD DEV\s+NY\+ > 5\s+NPTS\s*\n\s*'
                    r'([-\d.E+]+)\s+([-\d.E+]+)\s+(\d+)\s+(\d+)', txt, re.S)
    out_d = {
        'wing': pats.get('JLOC/ILOC'),
        'te_base': pats.get('KLOC/ILOC'),
        'tip': pats.get('JLOC/KLOC'),
        'aggregate': (dict(avg=float(agg.group(1)), ngt5=int(agg.group(3)),
                           npts=int(agg.group(4))) if agg else None)}
    # ★ A MISSING PATCH MUST BE LOUD.  Returning None silently is how an empty
    #   y+ column reaches a published dataset -- the same class as `.get`-ing a
    #   ramp-honesty field with a default.  A viscous run has all four.
    missing = [k for k, v in out_d.items() if v is None]
    if missing and out_d['aggregate'] is not None:
        raise RuntimeError(
            f'{out}: y+ patches {missing} not found, but an aggregate block '
            f'is present -- the parser is out of step with the output format, '
            f'not the run.  Found keys: {sorted(pats)}')
    return out_d


def blockforce_shares(d: Path):
    """Per-block Cdp / Cdv / wetted area from blockforce.dat, finest level.

    ★ This is how the TE base's exposure was BOUNDED rather than argued from
    its area: the file gives Cdv exactly 0.0000000 on that block (geometry --
    a streamwise face normal) against Cdp = -3.0 % of the total.
    """
    f = d / 'blockforce.dat'
    if not f.is_file():
        return {}
    rows = [l.split() for l in f.read_text().splitlines()
            if l.split() and l.split()[0].isdigit()]
    if not rows:
        return {}
    it = max(int(r[1]) for r in rows)
    last = {int(r[0]): r for r in rows if int(r[1]) == it}
    tot = max(last)                      # the aggregate row has the top index
    def g(b, k):
        return float(last[b][k]) if b in last else float('nan')
    return dict(total_cdp=g(tot, 4), total_cdv=g(tot, 5),
                total_area=g(tot, 9),
                te_cdp=g(4, 4), te_cdv=g(4, 5), te_area=g(4, 9),
                tip_cdp=g(7, 4), tip_cdv=g(7, 5), tip_area=g(7, 9),
                wing_cdp=g(1, 4), wing_cdv=g(1, 5), wing_area=g(1, 9))


def read_wall_prt(prt: Path, grid: int = 1):
    """Wall points of one GRID from cfl3d.prt.

    The printout is blocked by ``BLOCK n (GRID g) IDIM,JDIM,KDIM=...`` headers
    and only the grids that actually carry a wall appear.  Grid 1 is the main
    C-H block, whose K0 wall segment is the wing surface; grid 2 is the blunt-TE
    tail block and grid 3 the tip O-block, and neither is needed for sectional
    Cp.  Each data row is ``I J K X Y Z U V W P T MACH cp mut``.
    """
    want, rows = False, []
    for line in open(prt, errors='replace'):
        m = re.match(r'\s*BLOCK\s+\d+\s+\(GRID\s+(\d+)\)', line)
        if m:
            want = int(m.group(1)) == grid
            continue
        p = line.split()
        if not want or len(p) != 14:
            continue
        try:
            v = [float(t) for t in p]
        except ValueError:
            continue
        rows.append((int(v[0]), int(v[1]), int(v[2]), v[3], v[4], v[5], v[12]))
    if not rows:
        raise RuntimeError(f'{prt}: no wall rows for grid {grid}')
    return rows


# ---------------------------------------------------------------------------
#  sectional Cp at the experimental stations
# ---------------------------------------------------------------------------

def station_cp(rows, jte0, jte1, z_tip):
    """Cp at each experimental station, interpolated spanwise.

    The section grid is one template mapped to every span station, so a
    chordwise index ``j`` is the same relative position at every station and
    the interpolation can be done index-wise in the spanwise direction -- no
    resampling in x, which would smear the leading-edge peak.
    """
    surf = [r for r in rows if r[2] == 1 and jte0 <= r[1] <= jte1]
    ii = sorted({r[0] for r in surf})
    jj = sorted({r[1] for r in surf})
    by = {(r[0], r[1]): r for r in surf}
    z_of = np.array([by[(i, jj[0])][5] for i in ii])
    eta_of = z_of / z_tip

    out = {}
    for eta in EXP_STATIONS:
        k = int(np.clip(np.searchsorted(eta_of, eta), 1, len(ii) - 1))
        i0, i1 = ii[k - 1], ii[k]
        e0, e1 = eta_of[k - 1], eta_of[k]
        w = 0.0 if e1 == e0 else (eta - e0) / (e1 - e0)
        xs, ys, cps = [], [], []
        for j in jj:
            a, b = by[(i0, j)], by[(i1, j)]
            xs.append((1 - w) * a[3] + w * b[3])
            ys.append((1 - w) * a[4] + w * b[4])
            cps.append((1 - w) * a[6] + w * b[6])
        x = np.array(xs)
        # local chord fraction from the planform, in root-chord units
        z = eta * z_tip
        xle = math.tan(math.radians(M6_SWEEP_LE_DEG)) * z
        c = 1.0 - (1.0 - M6_CHORD_TIP / M6_CHORD_ROOT) * (z / z_tip)
        out[eta] = dict(x_c=(x - xle) / c, y_c=np.array(ys) / c,
                        cp=np.array(cps), eta_grid=(1 - w) * e0 + w * e1)
    return out


def upstream_supersonic_depth(x, cp, x_shock, window=0.08):
    """How far BELOW Cp* the flow gets in the ``window`` just upstream of the
    detected crossing.

    ★★★ THIS IS A PREMISE CHECK ON THE DETECTOR, not a physical quantity.
    ``shock_metrics`` returns, by its own documented contract, the **LAST**
    supersonic->subsonic crossing of Cp*.  That is the terminating shock only
    when the flow upstream of it is DECISIVELY supersonic.  Where a section
    instead carries a long marginally-sonic plateau -- Cp hovering within a
    few hundredths of Cp* -- the "last crossing" is wherever the plateau
    finally drifts above Cp*, which can sit a quarter of a chord downstream of
    the actual compression.  Measured on the M6 tip station: the compression
    is at x/c 0.236 and the last crossing is at 0.469.

    Returning a small number therefore means **the detector's premise fails
    here and its output is not a shock position** -- not that the shock is
    weak.  The reference dataset publishes this column so the artifact
    declares itself, and ``grid_convergence`` refuses an error bar for any
    station where it fails on any rung.

    ★ The fix is deliberately NOT "use a different rule at that station".
    The pyFP3D side of every D07 comparison is read with this SAME
    ``shock_metrics``; a reference side read with a different rule would make
    the two numbers different things, which is the criterion defect this
    project has now hit six times.
    """
    x = np.asarray(x, float)
    cp = np.asarray(cp, float)
    m = (x > x_shock - window) & (x < x_shock)
    if not m.any():
        return float('nan')
    return float(cp_critical(MACH) - cp[m].min())


#: below this the upstream flow is not decisively supersonic and the detected
#: crossing is a Cp*-grazing artifact.  Calibrated on the measured spread, not
#: chosen by hand: the five stations with a genuine terminating shock read
#: 0.175-0.438 on every rung, the tip station reads 0.045/0.017/0.005, and the
#: experiment's own curves read 0.349-0.808 at the six outboard stations.  Any
#: cut inside 0.05-0.17 separates the same two groups, so the number is a
#: CALIBRATION of a 4x gap, with the same status as the EW forcing and the
#: taper r_c -- it is not a physical threshold.
DEPTH_MIN = 0.05


def split_surfaces(st):
    """Split one station's contour at the leading edge (minimum x/c)."""
    i = int(np.argmin(st['x_c']))
    lo = slice(0, i + 1)
    up = slice(i, None)
    return ({k: st[k][up] for k in ('x_c', 'y_c', 'cp')},
            {k: st[k][lo] for k in ('x_c', 'y_c', 'cp')})


# ---------------------------------------------------------------------------
#  output
# ---------------------------------------------------------------------------

#: ★★ ONE SCHEMA for both 3-D datasets.  The 2-D round shipped Euler with 30
#: columns and RANS with 32, which `csv.DictReader` reads without complaint --
#: a reader written against one file silently gets empty strings from the
#: other.  Euler rows here carry ``turb = none`` and empty y+ / TE-share
#: columns; every column exists in both files.
FORCE_COLUMNS = [
    'level', 'turb', 'points', 'nj', 'nk', 'ni', 'n_tail', 'k_crit',
    'basin_depth',
    'mach', 'alpha_deg', 're_root_chord', 'h1_wall',
    'cl', 'cd', 'cd_pressure', 'cd_friction', 'cm_quarter_chord',
    'cl_raw_sref1', 'cd_raw_sref1', 's_ref_half_planform', 'wetted_area',
    'resid_final', 'resid_decades', 'ncyc_fine', 'wall_s',
    # ★ y+ PER WALL PATCH -- see yplus_patches().  Neither the wing figure nor
    #   the aggregate may be quoted alone as "the y+".
    'yplus_wing_avg', 'yplus_wing_max', 'n_yplus_gt5_wing', 'n_pts_wing',
    'yplus_te_base_avg', 'yplus_te_base_max', 'n_pts_te_base',
    'yplus_tip_avg', 'yplus_aggregate_avg', 'n_yplus_gt5_all', 'n_pts_all',
    # ★ the TE base's measured exposure, which BOUNDS the coarse patch
    'te_base_area_frac', 'te_base_cdv_share', 'te_base_cdp_share',
    'note',
]
CP_COLUMNS = ['level', 'turb', 'eta_requested', 'eta_grid', 'x_c', 'y_c',
              'cp', 'surface']
SHOCK_COLUMNS = ['level', 'turb', 'eta_requested', 'eta_grid', 'surface',
                 'has_shock',
                 'x_shock', 'n_cells', 'monotone', 'cp_min', 'cp_pre_shock',
                 'cp_post_shock', 'cp_critical', 'upstream_depth',
                 'detector_premise']
#: built dynamically in main() -- the rung list is not fixed at three
GC_BASE = ['quantity', 'turb']


def implied_order(h, vals, i):
    """The order p that reproduces the observed delta ratio at triple i.

    ★★★ THE RATIO NEEDS A CALIBRATION AND FOR A LONG TIME HAD NONE.  Comparing
    it against 1.0 asks only "did the deltas stop shrinking"; it does not ask
    whether they shrank at the rate a converging scheme produces.  For a
    p-order quantity the ratio is ``|h3^p - h2^p| / |h2^p - h1^p|`` -- a
    property of the LADDER as well as of p.  On this ladder (h ~ N^(-1/3) on
    the total point count, near-uniform at r = 1.376 / 1.340) that is 0.496 at
    p = 2 and 0.675 at p = 1, so "ratio < 1" admits order-0.2 convergence.

    ★ Reporting the implied order is also what lets a reader price an error
    bar: at p = 0.68 halving it needs a 2.8x point-count increase, not 1.4x.

    ★★ This is how the first-order defect was found -- cd's ratio of 0.744
    implies p = 0.68, which a kappa = 1/3 scheme should not give, and that
    inconsistency pointed at the scheme rather than at any single number.
    """
    h1, h2, h3 = h[i], h[i + 1], h[i + 2]
    d12, d23 = vals[i + 1] - vals[i], vals[i + 2] - vals[i + 1]
    if d12 == 0:
        return float('inf'), float('nan')
    ratio = abs(d23) / abs(d12)
    ps = np.linspace(0.01, 6.0, 6000)
    rr = np.abs(h3 ** ps - h2 ** ps) / np.abs(h2 ** ps - h1 ** ps)
    j = int(np.argmin(np.abs(rr - ratio)))
    p = float(ps[j]) if abs(rr[j] - ratio) < 5e-3 else float('nan')
    return ratio, p


def grid_convergence(rows_by_level, shock_by_level, order):
    """Per quantity: every rung, every delta, and the ratio + implied order on
    EVERY consecutive triple.

    ★★ With four rungs the two triples answer different questions: the finest
    one gives the error bar, and the pair together says whether the quantity is
    ENTERING the asymptotic range (ratio falling) or leaving it.  Three rungs
    could only ever give a single number with nothing to compare it to.

    The error bar comes from the FINEST triple, and only when that triple is
    asymptotic AND the quantity's detector premise held on every rung.
    """
    lv = [l for l in order if l in rows_by_level]
    N = np.array([float(rows_by_level[l]['points']) for l in lv])
    h = N ** (-1.0 / 3.0)
    h = h / h[0]

    q, bad = {}, set()
    for l in lv:
        for k in ('cl', 'cd', 'cd_pressure', 'cm_quarter_chord'):
            q.setdefault(k, {})[l] = float(rows_by_level[l][k])
    for l, sh in shock_by_level.items():
        if l not in lv:
            continue
        for sr in sh:
            if sr['surface'] != 'upper' or sr['x_shock'] == '':
                continue
            name = f"x_shock_upper_eta{sr['eta_requested']}"
            q.setdefault(name, {})[l] = float(sr['x_shock'])
            q.setdefault(f"cp_min_upper_eta{sr['eta_requested']}",
                         {})[l] = float(sr['cp_min'])
            # ★ one bad rung disqualifies the station: the deltas would be
            #   differences between DIFFERENT FEATURES.
            if sr['detector_premise'].startswith('FAILS'):
                bad.add(name)

    out = []
    for name, v in q.items():
        if not all(l in v for l in lv):
            continue
        vals = np.array([v[l] for l in lv])
        rec = {'quantity': name}
        for l in lv:
            rec[l] = f'{v[l]:.6f}'
        for i in range(len(lv) - 1):
            rec[f'delta_{lv[i]}_{lv[i+1]}'] = f'{vals[i+1] - vals[i]:+.6f}'
        last = None
        for i in range(len(lv) - 2):
            r, p = implied_order(h, vals, i)
            tag = f'{lv[i]}{lv[i+1]}{lv[i+2]}'
            rec[f'ratio_{tag}'] = f'{r:.3f}'
            rec[f'order_{tag}'] = ('' if p != p else f'{p:.2f}')
            last = (r, tag)
        asym = last is not None and last[0] < 1.0 and name not in bad
        rec['asymptotic'] = 'yes' if asym else 'no'
        rec['basis'] = last[1] if last else ''
        rec['error_bar'] = (f'{abs(vals[-1] - vals[-2]):.6f}' if asym else
                            ('NONE (detector premise fails -- Cp*-grazing)'
                             if name in bad else 'NONE (not asymptotic)'))
        out.append(rec)
    cols = (GC_BASE + list(lv)
            + [f'delta_{lv[i]}_{lv[i+1]}' for i in range(len(lv) - 1)]
            + sum([[f'ratio_{lv[i]}{lv[i+1]}{lv[i+2]}',
                    f'order_{lv[i]}{lv[i+1]}{lv[i+2]}']
                   for i in range(len(lv) - 2)], [])
            + ['asymptotic', 'basis', 'error_bar'])
    return out, cols


def write_csv(path: Path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(columns), extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    # ★ relative_to raises for a path outside the repo (e.g. --out /tmp/...).
    #   This exact defect was fixed once already in cfl3d_runner.py's _write_csv
    #   and was not backported here -- discipline #9, two modules one fix.
    try:
        shown = path.relative_to(HERE.parents[2])
    except ValueError:
        shown = path
    print(f'  -> {shown}  ({len(rows)} rows)')


def turbulence_spread(forces, levels, turbs=('sst', 'sa')):
    """|SST - SA| per level and quantity -- the D10 resolution floor.

    ★★ Measured L1/L2: the model spread in cl is 3.3-3.7 %, against a
    rung-to-rung grid delta of ~1.5 %.  ⇒ **the turbulence model, not the
    grid, is the dominant uncertainty in this dataset**, so a gate built on
    the grid error bar alone would claim a precision the data does not have.
    """
    by = {(r['level'], r['turb']): r for r in forces}
    rows = []
    for lv in levels:
        a, b = by.get((lv, turbs[0])), by.get((lv, turbs[1]))
        if not (a and b):
            continue
        for k in ('cl', 'cd', 'cd_pressure', 'cd_friction',
                  'cm_quarter_chord'):
            va, vb = float(a[k]), float(b[k])
            rows.append(dict(
                level=lv, quantity=k, **{turbs[0]: f'{va:.6f}',
                                         turbs[1]: f'{vb:.6f}'},
                abs_spread=f'{abs(vb - va):.6f}',
                rel_spread_pct=('' if va == 0 else
                                f'{abs(vb - va) / abs(va) * 100:.3f}')))
    cols = ['level', 'quantity', turbs[0], turbs[1], 'abs_spread',
            'rel_spread_pct']
    return rows, cols


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--from-runs', required=True,
                   help='directory holding the run directories')
    p.add_argument('--model', default='euler', choices=('euler', 'rans'))
    p.add_argument('--turb', default='sst,sa',
                   help='rans only: comma-separated turbulence models')
    p.add_argument('--levels', nargs='*',
                   default=['L1', 'L2', 'L3', 'L4', 'REF'])
    p.add_argument('--suffix', default='',
                   help="run-dir suffix, e.g. '_safe' (euler only)")
    p.add_argument('--out', default=None)
    a = p.parse_args(argv)

    root = Path(a.from_runs)
    default_out = f'{a.model}_onera_m6'
    out = Path(a.out) if a.out else HERE / default_out
    turbs = tuple(a.turb.split(',')) if a.model == 'rans' else ('none',)

    forces, cps, shocks = [], [], []
    by_lt, sh_by_lt = {}, {}          # keyed (level, turb)

    for lv in a.levels:
        for tb in turbs:
            d = (root / f'{lv}{a.suffix}' if a.model == 'euler'
                 else root / f'{lv}_{tb}')
            if not (d / 'cfl3d.out').is_file():
                print(f'  ! {lv}/{tb}: no run in {d}')
                continue
            # ★★★ A DIVERGED RUN STILL WRITES A SUMMARY -- see run_status.
            ok, why, last_it, want_it = run_status(d)
            if not ok:
                raise RuntimeError(
                    f'{lv}/{tb}: refusing to publish -- {why}.  A CFL3D '
                    f'summary exists for failed runs too, so it is not '
                    f'evidence of convergence.')
            g = build_m6(level=lv, model=a.model, verbose=False)
            sec = g.sec_proto
            smry = read_summary(d / 'cfl3d.out')
            res = read_residual(d / 'clcd_total.dat')
            yp = yplus_patches(d / 'cfl3d.out')
            bf = blockforce_shares(d)

            def yv(patch, key, fmt='{:.4f}'):
                v = yp.get(patch)
                return '' if not v or key not in v else fmt.format(v[key])

            def share(num, den):
                if not bf or den not in bf or bf[den] == 0:
                    return ''
                return f'{bf[num] / bf[den] * 100:+.3f}'

            row = dict(
                level=lv, turb=tb, points=g.info['points'],
                nj=sec.nj, nk=sec.nk, ni=g.grid[0].shape[0],
                n_tail=sec.n_tail, k_crit=g.k_crit,
                basin_depth=f'{g.basin_depth:.4f}',
                mach=f'{MACH:.4f}', alpha_deg=f'{ALPHA:.2f}',
                re_root_chord=f'{M6_RE_ROOT_CHORD:.4e}',
                h1_wall=f'{sec.h1:.4e}',
                cl=f"{smry['CL'] / S_HALF:.6f}",
                cd=f"{smry['CD'] / S_HALF:.6f}",
                cd_pressure=f"{smry['CDp'] / S_HALF:.6f}",
                cd_friction=f"{smry['CDv'] / S_HALF:.6f}",
                cm_quarter_chord=f"{smry['CMZ'] / S_HALF:.6f}",
                cl_raw_sref1=f"{smry['CL']:.6f}",
                cd_raw_sref1=f"{smry['CD']:.6f}",
                s_ref_half_planform=f'{S_HALF:.5f}',
                wetted_area=f"{smry['wetted_area']:.6f}",
                yplus_wing_avg=yv('wing', 'avg'),
                yplus_wing_max=yv('wing', 'max'),
                n_yplus_gt5_wing=yv('wing', 'ngt5', '{:d}'),
                n_pts_wing=yv('wing', 'npts', '{:d}'),
                yplus_te_base_avg=yv('te_base', 'avg', '{:.2f}'),
                yplus_te_base_max=yv('te_base', 'max', '{:.2f}'),
                n_pts_te_base=yv('te_base', 'npts', '{:d}'),
                yplus_tip_avg=yv('tip', 'avg', '{:.5f}'),
                yplus_aggregate_avg=yv('aggregate', 'avg', '{:.3f}'),
                n_yplus_gt5_all=yv('aggregate', 'ngt5', '{:d}'),
                n_pts_all=yv('aggregate', 'npts', '{:d}'),
                te_base_area_frac=share('te_area', 'total_area'),
                te_base_cdv_share=share('te_cdv', 'total_cdv'),
                te_base_cdp_share=share('te_cdp', 'total_cdp'),
                note=(f'{a.model}, alpha = experimental UNCORRECTED; grid '
                      f'normalised to root chord = 1; coefficients divided by '
                      f'the half-planform area (deck SREF = 1)'),
                **res)
            row['ncyc_fine'] = res['ncyc']
            forces.append(row)
            by_lt[(lv, tb)] = row

            st = station_cp(read_wall_prt(d / 'cfl3d.prt', grid=1),
                            sec.jte0, sec.jte1, g.info['z_tip'])
            lv_sh = []
            for eta in EXP_STATIONS:
                up, lo = split_surfaces(st[eta])
                for name, side in (('upper', up), ('lower', lo)):
                    order = np.argsort(side['x_c'])
                    for i2 in order:
                        cps.append(dict(
                            level=lv, turb=tb, eta_requested=f'{eta:.2f}',
                            eta_grid=f"{st[eta]['eta_grid']:.4f}",
                            x_c=f"{side['x_c'][i2]:.6f}",
                            y_c=f"{side['y_c'][i2]:.6f}",
                            cp=f"{side['cp'][i2]:.6f}", surface=name))
                    m = shock_metrics(side['x_c'], side['cp'], MACH)
                    pre = post = depth = premise = ''
                    if m['has_shock']:
                        xs = m['x_shock']
                        xa = side['x_c'][order]
                        ca = side['cp'][order]
                        pre = f'{float(np.interp(xs - 0.05, xa, ca)):.6f}'
                        if xs + 0.05 <= xa.max():
                            post = f'{float(np.interp(xs + 0.05, xa, ca)):.6f}'
                        dep = upstream_supersonic_depth(xa, ca, xs)
                        depth = f'{dep:.6f}'
                        premise = ('ok' if dep >= DEPTH_MIN else
                                   'FAILS -- Cp*-grazing, not a shock '
                                   'position')
                    rec = dict(
                        level=lv, turb=tb, eta_requested=f'{eta:.2f}',
                        eta_grid=f"{st[eta]['eta_grid']:.4f}", surface=name,
                        has_shock=int(m['has_shock']),
                        x_shock='' if not m['has_shock']
                        else f"{m['x_shock']:.6f}",
                        n_cells=m['n_cells'], monotone=int(m['monotone']),
                        cp_min=f"{m['cp_min']:.6f}", cp_pre_shock=pre,
                        cp_post_shock=post,
                        cp_critical=f'{cp_critical(MACH):.6f}',
                        upstream_depth=depth, detector_premise=premise)
                    shocks.append(rec)
                    lv_sh.append(rec)
            sh_by_lt[(lv, tb)] = lv_sh
            print(f'  {lv}/{tb}: cl {row["cl"]} cd {row["cd"]} '
                  f'cdv {row["cd_friction"]} '
                  f'y+wing {row["yplus_wing_avg"] or "n/a"} '
                  f'|R| {row["resid_final"]} ({row["resid_decades"]} dec)')

    write_csv(out / 'forces.csv', FORCE_COLUMNS, forces)
    write_csv(out / 'cp_stations.csv', CP_COLUMNS, cps)
    write_csv(out / 'shock.csv', SHOCK_COLUMNS, shocks)

    ladder = [l for l in a.levels if l != 'REF']
    gc_all, gc_cols = [], None
    for tb in turbs:
        rb = {l: by_lt[(l, tb)] for l in ladder if (l, tb) in by_lt}
        sb = {l: sh_by_lt[(l, tb)] for l in ladder if (l, tb) in sh_by_lt}
        if len(rb) < 3:
            print(f'  ! {tb}: only {len(rb)} rungs -- no convergence table')
            continue
        gc, gc_cols = grid_convergence(rb, sb, [l for l in ladder if l in rb])
        for r in gc:
            r['turb'] = tb
        gc_all += gc
    if gc_all:
        write_csv(out / 'grid_convergence.csv', gc_cols, gc_all)
        n_as = sum(1 for r in gc_all if r['asymptotic'] == 'yes')
        print(f'\n  asymptotic: {n_as} of {len(gc_all)} quantities have an '
              f'error bar')
        for r in gc_all:
            if r['asymptotic'] != 'yes':
                key = f"ratio_{r['basis']}"
                print(f"    RECORDED only: {r['quantity']} [{r['turb']}] "
                      f"({key} {r.get(key, '?')})")

    if a.model == 'rans' and len(turbs) == 2:
        ts, ts_cols = turbulence_spread(forces, ladder, turbs)
        if ts:
            write_csv(out / 'turbulence_spread.csv', ts_cols, ts)
            cl = [float(r['rel_spread_pct']) for r in ts
                  if r['quantity'] == 'cl' and r['rel_spread_pct']]
            if cl:
                print(f'\n  |{turbs[0].upper()}-{turbs[1].upper()}| in cl: '
                      f'{min(cl):.2f}-{max(cl):.2f} % '
                      f'-- the model spread, and the D10 resolution floor')
    return 0


if __name__ == '__main__':
    sys.exit(main())
