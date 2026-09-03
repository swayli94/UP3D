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
                 'cp_post_shock', 'cp_critical']
GC_COLUMNS = ['quantity', 'L1', 'L2', 'L3', 'delta_L1_L2', 'delta_L2_L3',
              'ratio', 'asymptotic', 'error_bar']


def grid_convergence(rows_by_level, shock_by_level):
    """Per quantity: the three rungs, both deltas, and whether the deltas are
    SHRINKING.

    ★★ ``ratio = |L2->L3| / |L1->L2|``.  Below 1 the rung-to-rung differences
    are shrinking and ``delta_L2_L3`` is usable as an error bar; at or above 1
    the ladder has not reached the asymptotic range and the quantity is
    RECORDED with no error bar, however small the delta happens to look.  This
    is the distinction the previous, inconsistent ladder could not make.
    """
    out = []
    q = {}
    for lv, r in rows_by_level.items():
        for k in ('cl', 'cd', 'cd_pressure', 'cm_quarter_chord'):
            q.setdefault(k, {})[lv] = float(r[k])
    for lv, sh in shock_by_level.items():
        for s in sh:
            if s['x_shock'] == '' or s['surface'] != 'upper':
                continue
            q.setdefault(f"x_shock_upper_eta{s['eta_requested']}", {})[lv] = \
                float(s['x_shock'])
    for name, v in q.items():
        if not all(l in v for l in ('L1', 'L2', 'L3')):
            continue
        d12, d23 = v['L2'] - v['L1'], v['L3'] - v['L2']
        ratio = abs(d23) / abs(d12) if d12 else float('inf')
        asym = ratio < 1.0
        out.append(dict(
            quantity=name,
            L1=f"{v['L1']:.6f}", L2=f"{v['L2']:.6f}", L3=f"{v['L3']:.6f}",
            delta_L1_L2=f'{d12:+.6f}', delta_L2_L3=f'{d23:+.6f}',
            ratio=f'{ratio:.3f}', asymptotic='yes' if asym else 'no',
            error_bar=f'{abs(d23):.6f}' if asym else 'NONE (not asymptotic)'))
    return out


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
    p.add_argument('--levels', nargs='*', default=['L1', 'L2', 'L3', 'REF'])
    p.add_argument('--out', default=None)
    a = p.parse_args(argv)

    root = Path(a.from_runs)
    out = Path(a.out) if a.out else HERE / 'euler_onera_m6'
    forces, cps, shocks = [], [], []
    by_level, sh_by_level = {}, {}

    for lv in a.levels:
        d = root / f'{lv}_safe'
        if not (d / 'cfl3d.out').is_file():
            print(f'  ! {lv}: no run in {d}')
            continue
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
                pre = post = ''
                if m['has_shock']:
                    xs = m['x_shock']
                    xa = side['x_c'][order]
                    ca = side['cp'][order]
                    pre = f'{float(np.interp(xs - 0.05, xa, ca)):.6f}'
                    if xs + 0.05 <= xa.max():
                        post = f'{float(np.interp(xs + 0.05, xa, ca)):.6f}'
                rec = dict(level=lv, eta_requested=f'{eta:.2f}',
                           eta_grid=f"{st[eta]['eta_grid']:.4f}",
                           surface=name, has_shock=int(m['has_shock']),
                           x_shock='' if not m['has_shock']
                           else f"{m['x_shock']:.6f}",
                           n_cells=m['n_cells'], monotone=int(m['monotone']),
                           cp_min=f"{m['cp_min']:.6f}", cp_pre_shock=pre,
                           cp_post_shock=post,
                           cp_critical=f'{cp_critical(MACH):.6f}')
                shocks.append(rec)
                lv_sh.append(rec)
        sh_by_level[lv] = lv_sh
        print(f'  {lv}: cl {row["cl"]} cd {row["cd"]} '
              f'|R| {row["resid_final"]} ({row["resid_decades"]} dec)')

    write_csv(out / 'forces.csv', FORCE_COLUMNS, forces)
    write_csv(out / 'cp_stations.csv', CP_COLUMNS, cps)
    write_csv(out / 'shock.csv', SHOCK_COLUMNS, shocks)
    gc = grid_convergence({k: v for k, v in by_level.items() if k != 'REF'},
                          {k: v for k, v in sh_by_level.items() if k != 'REF'})
    write_csv(out / 'grid_convergence.csv', GC_COLUMNS, gc)
    n_as = sum(1 for r in gc if r['asymptotic'] == 'yes')
    print(f'\n  asymptotic: {n_as} of {len(gc)} quantities have an error bar')
    for r in gc:
        if r['asymptotic'] != 'yes':
            print(f"    RECORDED only: {r['quantity']} (ratio {r['ratio']})")
    return 0


if __name__ == '__main__':
    sys.exit(main())
