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

FORCE_COLUMNS = [
    'level', 'points', 'nj', 'nk', 'ni', 'n_tail', 'k_crit', 'basin_depth',
    'mach', 'alpha_deg', 're_root_chord', 'h1_wall',
    'cl', 'cd', 'cd_pressure', 'cd_friction', 'cm_quarter_chord',
    'cl_raw_sref1', 'cd_raw_sref1', 's_ref_half_planform', 'wetted_area',
    'resid_final', 'resid_decades', 'ncyc_fine', 'wall_s', 'note',
]
CP_COLUMNS = ['level', 'eta_requested', 'eta_grid', 'x_c', 'y_c', 'cp',
              'surface']
SHOCK_COLUMNS = ['level', 'eta_requested', 'eta_grid', 'surface', 'has_shock',
                 'x_shock', 'n_cells', 'monotone', 'cp_min', 'cp_pre_shock',
                 'cp_post_shock', 'cp_critical', 'upstream_depth',
                 'detector_premise']
#: built dynamically in main() -- the rung list is not fixed at three
GC_BASE = ['quantity']


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
    print(f'  -> {path.relative_to(HERE.parents[2])}  ({len(rows)} rows)')


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--from-runs', required=True,
                   help='directory holding <LEVEL>_safe/ run directories')
    p.add_argument('--levels', nargs='*',
                   default=['L1', 'L2', 'L3', 'L4', 'REF'])
    p.add_argument('--suffix', default='',
                   help="run-dir suffix, e.g. '_safe'")
    p.add_argument('--out', default=None)
    a = p.parse_args(argv)

    root = Path(a.from_runs)
    out = Path(a.out) if a.out else HERE / 'euler_onera_m6'
    forces, cps, shocks = [], [], []
    by_level, sh_by_level = {}, {}

    for lv in a.levels:
        d = root / f'{lv}{a.suffix}'
        if not (d / 'cfl3d.out').is_file():
            print(f'  ! {lv}: no run in {d}')
            continue
        # ★★★ A DIVERGED RUN STILL WRITES A SUMMARY.  CFL3D emits
        # "SUMMARY OF FORCES" on its way out of a failed solve, so reading the
        # summary is NOT evidence the solve finished.  Measured: the M6 L2 rung
        # died with `NaN ... block 1 cycle 540` and its summary said
        # cl = +0.339608, which a runner checking only for the summary reported
        # as "ok" and would have published.  Success is established from
        # cfl3d.error's error code -- the only field that distinguishes the two,
        # and the only place the normal-termination banner is written.
        ok, why, last_it, want_it = run_status(d)
        if not ok:
            raise RuntimeError(
                f'{lv}: refusing to publish -- {why}.  A CFL3D summary exists '
                f'for failed runs too, so it is not evidence of convergence.')
        g = build_m6(level=lv, model='euler', verbose=False)
        s = g.sec_proto
        smry = read_summary(d / 'cfl3d.out')
        res = read_residual(d / 'clcd_total.dat')
        row = dict(
            level=lv, points=g.info['points'], nj=s.nj, nk=s.nk,
            ni=g.grid[0].shape[0], n_tail=s.n_tail, k_crit=g.k_crit,
            basin_depth=f"{g.basin_depth:.4f}",
            mach=f'{MACH:.4f}', alpha_deg=f'{ALPHA:.2f}',
            re_root_chord=f'{M6_RE_ROOT_CHORD:.4e}',
            h1_wall=f'{s.h1:.4e}',
            cl=f"{smry['CL'] / S_HALF:.6f}",
            cd=f"{smry['CD'] / S_HALF:.6f}",
            cd_pressure=f"{smry['CDp'] / S_HALF:.6f}",
            cd_friction=f"{smry['CDv'] / S_HALF:.6f}",
            cm_quarter_chord=f"{smry['CMZ'] / S_HALF:.6f}",
            cl_raw_sref1=f"{smry['CL']:.6f}",
            cd_raw_sref1=f"{smry['CD']:.6f}",
            s_ref_half_planform=f'{S_HALF:.5f}',
            wetted_area=f"{smry['wetted_area']:.6f}",
            note=('Euler, alpha = experimental UNCORRECTED; grid normalised to '
                  'root chord = 1; coefficients divided by the half-planform '
                  'area (deck SREF = 1)'),
            **res)
        row['ncyc_fine'] = res['ncyc']
        forces.append(row)
        by_level[lv] = row

        st = station_cp(read_wall_prt(d / 'cfl3d.prt', grid=1),
                        s.jte0, s.jte1, g.info['z_tip'])
        lv_sh = []
        for eta in EXP_STATIONS:
            up, lo = split_surfaces(st[eta])
            for name, side in (('upper', up), ('lower', lo)):
                order = np.argsort(side['x_c'])
                for i in order:
                    cps.append(dict(level=lv, eta_requested=f'{eta:.2f}',
                                    eta_grid=f"{st[eta]['eta_grid']:.4f}",
                                    x_c=f"{side['x_c'][i]:.6f}",
                                    y_c=f"{side['y_c'][i]:.6f}",
                                    cp=f"{side['cp'][i]:.6f}",
                                    surface=name))
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
                               'FAILS -- Cp*-grazing, not a shock position')
                rec = dict(level=lv, eta_requested=f'{eta:.2f}',
                           eta_grid=f"{st[eta]['eta_grid']:.4f}",
                           surface=name, has_shock=int(m['has_shock']),
                           x_shock='' if not m['has_shock']
                           else f"{m['x_shock']:.6f}",
                           n_cells=m['n_cells'], monotone=int(m['monotone']),
                           cp_min=f"{m['cp_min']:.6f}", cp_pre_shock=pre,
                           cp_post_shock=post,
                           cp_critical=f'{cp_critical(MACH):.6f}',
                           upstream_depth=depth, detector_premise=premise)
                shocks.append(rec)
                lv_sh.append(rec)
        sh_by_level[lv] = lv_sh
        print(f'  {lv}: cl {row["cl"]} cd {row["cd"]} '
              f'|R| {row["resid_final"]} ({row["resid_decades"]} dec)')

    write_csv(out / 'forces.csv', FORCE_COLUMNS, forces)
    write_csv(out / 'cp_stations.csv', CP_COLUMNS, cps)
    write_csv(out / 'shock.csv', SHOCK_COLUMNS, shocks)
    ladder = [l for l in a.levels if l != 'REF']
    gc, gc_cols = grid_convergence(
        {k: v for k, v in by_level.items() if k != 'REF'},
        {k: v for k, v in sh_by_level.items() if k != 'REF'}, ladder)
    write_csv(out / 'grid_convergence.csv', gc_cols, gc)
    n_as = sum(1 for r in gc if r['asymptotic'] == 'yes')
    print(f'\n  asymptotic: {n_as} of {len(gc)} quantities have an error bar')
    for r in gc:
        if r['asymptotic'] != 'yes':
            print(f"    RECORDED only: {r['quantity']} (ratio {r['ratio']})")
    return 0


if __name__ == '__main__':
    sys.exit(main())
