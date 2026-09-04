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

import argparse
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
#: set by main() from --dataset; both 3-D datasets share the cp_stations schema
DATA = HERE / 'euler_onera_m6'
#: ★ True only for the Euler dataset.  Gates the inviscid-vs-viscous direction
#: prediction, which has no basis when both sides are viscous.
INVISCID = True
MACH = 0.8395
STATIONS = (0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99)


def experiment_bracket(x, x_shock):
    """Width of the experimental interval containing the Cp* crossing.

    ★★★ A shock position INFERRED from discretely sampled Cp cannot be located
    more precisely than the sampling interval that brackets the crossing --
    roughly half of it, geometrically.  The M6 experiment samples every
    0.040-0.050 c, so a computed-vs-measured displacement smaller than
    ~0.020 c carries NO information about the shock's true position.

    ★★ Measured, and it changes the verdict: at SECOND order every station's
    displacement is 0.11-0.62x the local bracket, i.e. the direction test has
    NO POWER against this experiment at all.  The FIRST-order dataset's
    displacements (+0.025..+0.072 c, up to 1.4x the bracket) only looked
    convincing because the solution was wrong by more than the reference could
    resolve -- fixing the scheme removed the "signal".

    ⇒ This is the criterion-defect family again, in its own distinctive form:
    the criterion compared a quantity against a reference whose RESOLUTION is
    coarser than every difference being reported.
    """
    x = np.asarray(x, float)
    i = int(np.searchsorted(x, x_shock))
    return float(x[min(i, len(x) - 1)] - x[max(i - 1, 0)])


#: a displacement below this fraction of the experimental bracket is
#: UNRESOLVABLE.  0.5 is the geometric half-interval, i.e. a CALIBRATION of the
#: reference's sampling, with the same status as the EW forcing and DEPTH_MIN.
RESOLVE_FRAC = 0.5


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


def read_cfl3d(level='L4', turb=None):
    """Cp at the stations, for one rung and (for RANS) one turbulence model."""
    rows = [r for r in csv.DictReader(open(DATA / 'cp_stations.csv'))
            if r['level'] == level
            and (turb is None or r.get('turb', 'none') == turb)]
    if not rows:
        raise RuntimeError(
            f'{DATA}/cp_stations.csv has no rows for level={level!r} '
            f'turb={turb!r} -- check the ladder actually reached that rung '
            f'rather than publishing an empty comparison')
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


def main(argv=None):
    global DATA
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dataset', default='euler_onera_m6')
    ap.add_argument('--level', default=None,
                    help='rung to read the bias on; default = the finest '
                         'present in the dataset')
    ap.add_argument('--turb', default=None, help='rans only, e.g. sst')
    ap.add_argument('--levels-plot', default=None,
                    help='comma-separated rungs to draw; default = all present')
    a = ap.parse_args(argv)
    global INVISCID
    DATA = HERE / a.dataset
    INVISCID = 'euler' in a.dataset

    present = sorted({r['level'] for r in
                      csv.DictReader(open(DATA / 'cp_stations.csv'))
                      if r['level'] != 'REF'})
    # ★ the FINEST rung, taken from the file rather than hard-coded: reading a
    #   bias on L3 while an L4 exists would compare the experiment against a
    #   rung the dataset itself supersedes.
    level = a.level or sorted(present)[-1]
    plot_levels = (a.levels_plot.split(',') if a.levels_plot else present)
    turb = a.turb
    exp = read_experiment()
    lev = {l: read_cfl3d(l, turb) for l in plot_levels}
    cfd = lev[level]         # the bias is read on the FINEST rung present
    tag = f' [{turb}]' if turb else ''
    print(f'  dataset {a.dataset}{tag}, bias read on {level} '
          f'(rungs present: {", ".join(present)})')

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
            brk = experiment_bracket(exp[ekey]['upper'][0], se['x_shock'])
            rec['exp_bracket'] = f'{brk:.4f}'
            rec['dx_over_bracket'] = f'{abs(dx) / brk:.2f}'
            if abs(dx) < RESOLVE_FRAC * brk:
                # ★ THE RESOLUTION GATE COMES BEFORE THE SIGN.  Reading a sign
                #   off a difference finer than the reference's own sampling is
                #   reading noise -- in either direction.
                rec['direction'] = (
                    f'UNRESOLVABLE -- |dx| is {abs(dx)/brk:.2f}x the '
                    f'experiment\'s own {brk:.4f} c sampling interval')
            elif not INVISCID:
                # ★★★ THE PREDICTION IS EULER-ONLY AND MUST NOT BE
                #   TRANSPLANTED.  Its whole basis is that an INVISCID
                #   solution at the same alpha lacks the boundary-layer
                #   displacement that moves the measured shock upstream.  A
                #   RANS solution HAS that displacement, so there is no
                #   inviscid/viscous gap to predict a sign from: what remains
                #   between RANS and the experiment is turbulence modelling,
                #   transition placement and wind-tunnel corrections, none of
                #   which has a pre-registered direction here.
                #   Applying the Euler clause to RANS would have reported
                #   "as predicted" and "ANOMALY" verdicts with no premise
                #   behind either.
                rec['direction'] = (f'RECORDED, no predicted sign (viscous vs '
                                    f'viscous); dx {dx:+.4f} = '
                                    f'{abs(dx) / brk:.2f}x the bracket')
            elif dx <= 0:
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
            'depth_experiment', 'depth_euler', 'exp_bracket',
            'dx_over_bracket', 'direction']
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
    # ★ the previous version counted every non-"as predicted" row as
    #   "withdrawn on the detector premise", which mislabelled an ANOMALY and
    #   an UNRESOLVABLE as a detector problem.  Count each class.
    from collections import Counter
    cls = Counter(r['direction'].split(' --')[0].split(' (')[0]
                  for r in rows)
    print('  direction verdicts: ' + ', '.join(
        f'{k} {v}' for k, v in sorted(cls.items(), key=lambda kv: -kv[1])))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use('Agg')                    # headless, per the project rule
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(17.5, 8.2), sharex=True)
    for ax, eta in zip(axes.ravel(), STATIONS):
        ekey = min(exp, key=lambda z: abs(z - eta))
        styles = {'L1': dict(lw=0.7, ls=':', color='#b0b6bd'),
                  'L2': dict(lw=0.8, ls=':', color='#9aa0a6'),
                  'L3': dict(lw=0.9, ls='--', color='#5b8def'),
                  'L4': dict(lw=1.7, color='#1a56db')}
        for lv in plot_levels:
            sty = dict(styles.get(lv, dict(lw=1.2, color='#1a56db')))
            if lv == level:
                sty.update(lw=1.7, ls='-', color='#1a56db')
            for side in ('upper', 'lower'):
                x, cp = lev[lv][eta][side]
                ax.plot(x, cp, label=f'CFL3D {lv}{tag}'
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
    # ★★★ THE CAPTION IS COMPUTED FROM `rows`, NOT WRITTEN AS PROSE.
    #   The first RANS figure inherited the Euler caption verbatim and so
    #   asserted "only y/b = 0.20 resolves" (it is y/b = 0.99 for RANS), the
    #   Euler first-order displacement range, and a premise failing "on L1 AND
    #   L4" when this dataset has no L4.  Stale text on a regenerated artifact
    #   is the erratum-checklist class; deriving it removes the failure mode.
    res = [r for r in rows if r['direction'].startswith(('as predicted',
                                                         'RECORDED, no'))]
    unres = [r for r in rows if r['direction'].startswith('UNRESOLVABLE')]
    wd = [r for r in rows if r['direction'].startswith('WITHDRAWN')]
    brks = [float(r['exp_bracket']) for r in rows if r.get('exp_bracket')]
    fr = [float(r['dx_over_bracket']) for r in unres
          if r.get('dx_over_bracket')]
    pooled = [float(r['rms_pooled']) for r in rows]
    lines = [
        'ONERA M6, AGARD AR-138 TEST 2308',
        'M 0.8395, alpha 3.06 (experimental, UNCORRECTED)',
        f'CFL3D {"Euler" if INVISCID else "RANS"}, SECOND ORDER '
        f'(RKAP0 1/3, ICHK 2){tag}',
        f'bias read on {level}; rungs drawn: {", ".join(plot_levels)}',
        '',
        'RECORDED bias, not a gate.',
    ]
    lines += ([
        '(ruling 3: the gate against experiment belongs',
        'to FP+IBL, not to a model one level away.)',
        '',
        'PREDICTED SIGN: inviscid at the same alpha lacks',
        'the boundary-layer displacement that moves the',
        'measured shock upstream, so its shock must sit',
        'AFT.  A shock forward of the measurement would',
        'be an ANOMALY, not better agreement.',
    ] if INVISCID else [
        '',
        '*** NO PREDICTED SIGN.  The Euler prediction',
        '(inviscid lacks the displacement that moves the',
        'measured shock upstream) has NO basis here --',
        'RANS HAS that displacement.  What is left is',
        'turbulence modelling, transition and tunnel',
        'corrections, none with a pre-registered',
        'direction.  So the shock rows are RECORDED.',
    ])
    lines += [
        '',
        '*** RESOLUTION GATE.  The experiment samples Cp',
        f'every {min(brks):.3f}-{max(brks):.3f} c, so a shock position',
        'inferred from it is good to about half that.',
        f'{len(unres)} of {len(rows)} displacements are '
        + (f'{min(fr):.2f}-{max(fr):.2f}x' if fr else 'n/a'),
        'that interval => UNRESOLVABLE, in either sign.',
    ]
    if res:
        lines.append('Resolvable: ' + ', '.join(
            f"y/b {r['eta']} ({r['dx_over_bracket']}x)" for r in res))
    lines += [
        '',
        f'What survives: pooled Cp RMS {min(pooled):.4f}-{max(pooled):.4f},',
        'comparing Cp AT the measured points and',
        'inferring no position.',
        '',
        'Dash-dot = detected shock (blue CFD, red exp).',
        'Green dashes = Cp* = -0.328.',
    ]
    if wd:
        lines += ['', 'WITHDRAWN (Cp*-grazing, the detector\'s "last sonic',
                  'crossing" is not the compression): '
                  + ', '.join(f"y/b {r['eta']}" for r in wd)]
    axes.ravel()[-1].text(0.02, 0.98, '\n'.join(lines), va='top', ha='left',
                          fontsize=8.2, family='monospace')
    fig.suptitle(f'CFL3D {"Euler" if INVISCID else "RANS"}{tag} reference vs '
                 f'the committed seven-station experiment', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = DATA / 'cp_vs_experiment.png'
    fig.savefig(out, dpi=125)
    print(f'  -> {out.relative_to(REPO)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
