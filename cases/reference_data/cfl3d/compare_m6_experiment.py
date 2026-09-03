"""
Bias of the D07 Euler reference against the committed seven-station M6
experiment, plus the Cp comparison figure.

    python compare_m6_experiment.py

★★★ This is a RECORDED BIAS, not a gate, and the reason is ruling 3 of the
phase-six data request: the experiment sits at a different MODEL LEVEL from an
Euler solution.  Its bias has a KNOWN DIRECTION, written down in advance so the
first reader does not mistake it for a defect:

    the experiment is VISCOUS -- boundary-layer displacement reduces the
    effective camber and pushes the shock UPSTREAM -- so an inviscid solution
    at the SAME angle of attack must carry MORE lift and a shock FURTHER AFT.

★★ Per the criterion-defect list, a single-sided expectation has to say where
the OPPOSITE outcome lands: a Euler shock UPSTREAM of the measured one is
**not "better agreement"**, it contradicts the mechanism and means stopping to
look for a setup error.  ``direction`` in the CSV is that test, and it is
premise-gated -- see below.

★★★ CORRECTION, measured 2026-09-03.  An earlier version of this file put the
SUCTION PEAK in the same prediction ("inviscid must also peak DEEPER") and
duly reported all six valid stations as "half the mechanism missing".  That
was a defect in the prediction, not a finding about the solution: the peak has
a second mechanism of the opposite sign -- a finite grid always under-resolves
a leading-edge suction peak, making it shallower -- and on this ladder that
one DOMINATES.  Measured: ``cp_min`` deepens MONOTONICALLY across all three
rungs at every station, and the L2->L3 step alone covers 37-66 % of the whole
remaining gap to the experiment (grid_convergence.csv gives every
``cp_min_upper_eta*`` as non-asymptotic).  A quantity still moving by half its
own residual gap cannot test a model difference.  ⇒ the peak comparison is
RECORDED with no predicted sign; only the shock position carries one.

Writes euler_onera_m6/experiment_bias.csv and euler_onera_m6/cp_vs_experiment.png.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO))

from pyfp3d.post.shock import shock_metrics                # noqa: E402
from generate_m6_reference import (                        # noqa: E402
    DEPTH_MIN, upstream_supersonic_depth)

EXP = REPO / 'cases' / 'reference_data' / 'onera_m6_experiment' / 'experiment-Cp.dat'
DATA = HERE / 'euler_onera_m6'
MACH = 0.8395
STATIONS = (0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99)


def read_experiment():
    """Seven zones of NP, X/L, Y/b, Z/L, Cp.

    ★ Upper/lower is taken from the sign of Z/L, the surface ordinate.  That is
    valid HERE because the ONERA D section is symmetric; on a cambered section
    it is exactly the rule that mis-sorted 37 RAE2822 points (the D12 erratum),
    so it must not be copied to a cambered case.
    """
    zones, cur = [], None
    for ln in open(EXP, errors='replace'):
        if 'ZONE' in ln.upper():
            cur = []
            zones.append(cur)
            continue
        p = ln.split()
        if cur is None or len(p) != 5:
            continue
        try:
            cur.append([float(x) for x in p])
        except ValueError:
            pass
    out = {}
    for rows in zones:
        a = np.array(rows)
        eta = round(float(a[0, 2]), 3)
        up = a[a[:, 3] > 0]
        lo = a[a[:, 3] < 0]
        out[eta] = dict(
            upper=(up[np.argsort(up[:, 1]), 1], up[np.argsort(up[:, 1]), 4]),
            lower=(lo[np.argsort(lo[:, 1]), 1], lo[np.argsort(lo[:, 1]), 4]))
    return out


def read_cfl3d(level='L3'):
    rows = [r for r in csv.DictReader(open(DATA / 'cp_stations.csv'))
            if r['level'] == level]
    out = {}
    for r in rows:
        k = (float(r['eta_requested']), r['surface'])
        out.setdefault(k, []).append((float(r['x_c']), float(r['cp'])))
    res = {}
    for (eta, side), v in out.items():
        v.sort()
        res.setdefault(eta, {})[side] = (np.array([a for a, _ in v]),
                                         np.array([b for _, b in v]))
    return res


def main():
    exp = read_experiment()
    lev = {l: read_cfl3d(l) for l in ('L1', 'L2', 'L3')}
    cfd = lev['L3']

    rows = []
    for eta in STATIONS:
        ekey = min(exp, key=lambda z: abs(z - eta))
        rec = dict(eta=f'{eta:.2f}', eta_experiment=f'{ekey:.3f}')
        pooled = []
        for side in ('upper', 'lower'):
            xe, cpe = exp[ekey][side]
            xc, cpc = cfd[eta][side]
            m = (xe >= xc.min()) & (xe <= xc.max())
            interp = np.interp(xe[m], xc, cpc)
            d = interp - cpe[m]
            pooled.append(d)
            rec[f'rms_{side}'] = f'{np.sqrt(np.mean(d ** 2)):.6f}'
            rec[f'n_{side}'] = int(m.sum())
        allд = np.concatenate(pooled)
        rec['rms_pooled'] = f'{np.sqrt(np.mean(allд ** 2)):.6f}'

        # the two quantities whose direction is predicted
        se = shock_metrics(*exp[ekey]['upper'], MACH)
        sc = shock_metrics(*cfd[eta]['upper'], MACH)
        rec['x_shock_experiment'] = ('' if not se['has_shock']
                                     else f"{se['x_shock']:.6f}")
        rec['x_shock_euler'] = ('' if not sc['has_shock']
                                else f"{sc['x_shock']:.6f}")
        rec['cp_min_experiment'] = f"{se['cp_min']:.6f}"
        rec['cp_min_euler'] = f"{sc['cp_min']:.6f}"
        de = (upstream_supersonic_depth(*exp[ekey]['upper'], se['x_shock'])
              if se['has_shock'] else float('nan'))
        dc = (upstream_supersonic_depth(*cfd[eta]['upper'], sc['x_shock'])
              if sc['has_shock'] else float('nan'))
        rec['depth_experiment'] = f'{de:.6f}'
        rec['depth_euler'] = f'{dc:.6f}'

        if not (se['has_shock'] and sc['has_shock']):
            rec['d_x_shock'] = rec['d_cp_min'] = ''
            rec['direction'] = 'n/a (no shock on one side)'
        elif dc < DEPTH_MIN or de < DEPTH_MIN:
            # ★★ THE PREMISE GATE COMES FIRST.  Without it the tip station
            # passes a one-sided "is it downstream?" test for entirely the
            # wrong reason: its dx = +0.255 is FOUR TIMES the largest genuine
            # bias, and it is downstream because the detector returned a
            # Cp*-grazing point, not because the flow is inviscid.  A
            # one-sided band cannot tell a big right answer from a wrong one.
            rec['d_x_shock'] = f"{sc['x_shock'] - se['x_shock']:+.6f}"
            rec['d_cp_min'] = f"{sc['cp_min'] - se['cp_min']:+.6f}"
            rec['direction'] = ('WITHDRAWN -- detector premise fails '
                                '(Cp*-grazing); dx is not a shock bias')
        else:
            dx = sc['x_shock'] - se['x_shock']
            dpk = sc['cp_min'] - se['cp_min']
            rec['d_x_shock'] = f'{dx:+.6f}'
            rec['d_cp_min'] = f'{dpk:+.6f}'
            # inviscid at the same alpha: shock DOWNSTREAM.  Both signs are
            # stated so the opposite outcome has a landing spot.
            if dx <= 0:
                rec['direction'] = ('ANOMALY -- shock UPSTREAM of experiment, '
                                    'contradicts the viscous mechanism')
            else:
                # peak sign deliberately NOT tested -- see the correction in
                # the module docstring; d_cp_min is RECORDED beside it.
                rec['direction'] = 'as predicted (shock aft of experiment)'
        rows.append(rec)

    cols = ['eta', 'eta_experiment', 'rms_upper', 'n_upper', 'rms_lower',
            'n_lower', 'rms_pooled', 'x_shock_experiment', 'x_shock_euler',
            'd_x_shock', 'cp_min_experiment', 'cp_min_euler', 'd_cp_min',
            'depth_experiment', 'depth_euler', 'direction']
    with open(DATA / 'experiment_bias.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    print(f'  -> {(DATA / "experiment_bias.csv").relative_to(REPO)}')

    print(f'\n  {"y/b":>5s} {"RMS_up":>8s} {"RMS_lo":>8s} {"RMS":>8s} '
          f'{"x_s exp":>8s} {"x_s Eul":>8s} {"dx":>9s}  direction')
    for r in rows:
        print(f'  {r["eta"]:>5s} {r["rms_upper"]:>8s} {r["rms_lower"]:>8s} '
              f'{r["rms_pooled"]:>8s} {r["x_shock_experiment"] or "--":>8s} '
              f'{r["x_shock_euler"] or "--":>8s} {r["d_x_shock"] or "--":>9s}'
              f'  {r["direction"]}')
    p = [float(r['rms_pooled']) for r in rows]
    print(f'\n  pooled RMS {min(p):.4f} .. {max(p):.4f}  mean {np.mean(p):.4f}')
    ok = [r for r in rows if r['direction'].startswith('as predicted')]
    print(f'  direction as predicted at {len(ok)} of {len(rows)} stations; '
          f'{len(rows) - len(ok)} withdrawn on the detector premise')

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use('Agg')                    # headless, per the project rule
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(17.5, 8.2), sharex=True)
    for ax, eta in zip(axes.ravel(), STATIONS):
        ekey = min(exp, key=lambda z: abs(z - eta))
        for lv, sty in (('L1', dict(lw=0.8, ls=':', color='#9aa0a6')),
                        ('L2', dict(lw=0.9, ls='--', color='#5b8def')),
                        ('L3', dict(lw=1.7, color='#1a56db'))):
            for side in ('upper', 'lower'):
                x, cp = lev[lv][eta][side]
                ax.plot(x, cp, label=f'CFL3D Euler {lv}'
                        if side == 'upper' else None, **sty)
        for side, mk in (('upper', 'o'), ('lower', 's')):
            xe, cpe = exp[ekey][side]
            ax.plot(xe, cpe, mk, ms=4.2, mfc='none', mew=1.1, color='#c0392b',
                    label='experiment (TEST 2308)' if side == 'upper' else None)
        r = next(q for q in rows if q['eta'] == f'{eta:.2f}')
        bad = r['direction'].startswith('WITHDRAWN')
        ax.set_title(f'y/b = {eta:.2f}   pooled RMS {r["rms_pooled"]}'
                     + ('   [shock WITHDRAWN]' if bad else ''),
                     fontsize=10, color='#c0392b' if bad else 'black')
        if r['x_shock_euler']:
            ax.axvline(float(r['x_shock_euler']), color='#1a56db', lw=0.8,
                       ls='-.', alpha=0.55 if not bad else 0.9)
        if r['x_shock_experiment']:
            ax.axvline(float(r['x_shock_experiment']), color='#c0392b',
                       lw=0.8, ls='-.', alpha=0.55)
        ax.axhline(-0.3282, color='#2e7d32', lw=0.7, ls=(0, (4, 3)),
                   alpha=0.7)
        ax.invert_yaxis()
        ax.grid(alpha=.25, lw=.5)
        ax.set_xlim(-0.02, 1.02)
        if ax is axes.ravel()[0]:
            ax.legend(fontsize=8, loc='lower right')
    for ax in axes[1]:
        ax.set_xlabel('x/c')
    for ax in axes[:, 0]:
        ax.set_ylabel('$C_p$')
    axes.ravel()[-1].axis('off')
    axes.ravel()[-1].text(
        0.02, 0.95,
        'ONERA M6, AGARD AR-138 TEST 2308\n'
        'M 0.8395, alpha 3.06 (experimental, UNCORRECTED)\n\n'
        'CFL3D EULER vs VISCOUS experiment: the bias\n'
        'direction is KNOWN IN ADVANCE -- inviscid at the\n'
        'same alpha carries more lift and its shock sits\n'
        'FURTHER AFT.  A shock further FORWARD than the\n'
        'measurement would contradict that and would be\n'
        'an anomaly, not better agreement.\n\n'
        'RECORDED bias, not a gate (ruling 3: the gate\n'
        'against experiment belongs to FP+IBL, not to a\n'
        'model one level away).\n\n'
        'Only the SHOCK POSITION has a predicted sign.\n'
        'The suction peak does not: it is still deepening\n'
        'monotonically with refinement, the last rung\n'
        'covering 37-66% of the gap to experiment, so\n'
        'grid truncation dominates it.\n\n'
        'Dash-dot = detected shock (blue Euler, red exp).\n'
        'Green dashes = Cp* = -0.328.\n\n'
        'y/b = 0.99 is WITHDRAWN: Cp grazes Cp* over a\n'
        'long plateau there, so the detector\'s "last\n'
        'sonic crossing" lands 0.23c aft of the actual\n'
        'compression.  Its delta ratio was 0.337 -- it\n'
        'would have PASSED the convergence test, on a\n'
        'non-feature.',
        va='top', ha='left', fontsize=8.5, family='monospace')
    fig.suptitle('D07 Euler reference vs the committed seven-station '
                 'experiment', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = DATA / 'cp_vs_experiment.png'
    fig.savefig(out, dpi=125)
    print(f'  -> {out.relative_to(REPO)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
