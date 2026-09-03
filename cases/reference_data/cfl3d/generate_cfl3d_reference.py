"""
Reference-data generator: CFL3D 6.7 Euler and RANS solutions for the
2-D gates D05 / D06 (Euler) and D08 / D09 (RANS).

    python generate_cfl3d_reference.py --set euler          # D05 + D06
    python generate_cfl3d_reference.py --set rans           # D08 + D09
    python generate_cfl3d_reference.py --set all --jobs 8
    python generate_cfl3d_reference.py --set euler --levels L3   # one rung
    python generate_cfl3d_reference.py --derive-only        # rebuild the
        # derived CSVs (shock.csv, grid_convergence.csv, turbulence_spread.csv)
        # from the committed cp_*.csv + forces.csv, with NO solver run

★ What this is and is NOT (data request note, docs/dev_phase_six/
  20260824-0400-gate-taxonomy-analysis.md section 4):
  **a CFL3D solution is another NUMERICAL solution, not truth.** Every gate
  built on it must read "difference from a recognised solver's Euler/RANS
  solution at a stated grid and convergence caliber", never "error".
  That is why every row here carries its grid level, its residual drop and
  its y+, and why every case is run on a three-rung ladder: a single grid is
  one number with no error bar.

★ alpha caliber (user ruling 2026-08-24): every case uses the **experimental
  angle of attack, uncorrected** -- the same caliber the pyFP3D side uses --
  so the reference carries no correction constant of its own.

Outputs (four dataset directories beside this file):
    euler_naca0012/   D05      euler_rae2822/   D06
    rans_naca0012/    D08      rans_rae2822/    D09
each with
    cp_<tag>_<model>_<level>.csv   x_c, y_c, cp, mach_local, surface
    forces.csv                     integrals + convergence + grid caliber
    shock.csv                      shock position / pre- and post-shock Cp
    grid_convergence.csv           value per rung + the finest-pair delta
    turbulence_spread.csv          |SST - SA| per quantity  (RANS only)
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))

import cfl3d_runner as R                                        # noqa: E402
from pyfp3d.post.shock import cp_critical, shock_metrics        # noqa: E402


# ---------------------------------------------------------------------------
#  case lists
# ---------------------------------------------------------------------------
#
# Tags encode the condition, not the request number, so a row stays readable
# without the plan document open.  The request id (2D-1 ... R-7) is carried in
# the `request` column of forces.csv.

def _euler(tag, geom, mach, alpha, re_mil, request, note=''):
    return R.Case(tag=tag, geometry=geom, mach=mach, alpha=alpha,
                  re_mil=re_mil, ivisc=0, request=request, note=note)


def _rans(tag, geom, mach, alpha, re_mil, ivisc, request, x_tr=None,
          note='', solver=None):
    c = R.Case(tag=tag, geometry=geom, mach=mach, alpha=alpha,
               re_mil=re_mil, ivisc=ivisc, x_tr=x_tr, request=request,
               note=note, solver=solver or {})
    c.keywords = R.default_keywords(c)      # the eddy-viscosity cap
    return c


#: Reynolds number on an Euler run.  CFL3D reads REUE,MIL from the deck but
#: ivisc = 0 evaluates no viscous term, so the value is inert; it is set to the
#: matching experiment's Re where one exists purely so the row is readable.
#: `--verify-euler-re` measures the inertness rather than asserting it.
EULER_CASES = [
    _euler('n0012_m0500_a2.00',  'naca0012', 0.500,  2.00, 6.0, '2D-1',
           'subsonic core + Kutta; second opinion on our own naca0012_m05'),
    _euler('n0012_m0800_a1.25',  'naca0012', 0.800,  1.25, 6.0, '2D-2',
           'THE M1 target condition (naca0012_m080 shock reference)'),
    _euler('n0012_m0720_a0.00',  'naca0012', 0.720,  0.00, 6.0, '2D-3',
           'subcritical->supercritical, alpha=0 removes lift/wake coupling'),
    _euler('n0012_m0750_a0.00',  'naca0012', 0.750,  0.00, 6.0, '2D-3',
           'as 2D-3 lower Mach; single-variable read on artificial dissipation'),
    _euler('n0012_m0778_a2.03',  'naca0012', 0.778,  2.03, 6.0, '2D-4',
           'same condition as naca0012_experiment M0.778'),
    _euler('n0012_m0803_am0.10', 'naca0012', 0.803, -0.10, 6.5, '2D-4',
           'same condition as naca0012_experiment M0.803'),
    _euler('rae2822_m0725_a2.55', 'rae2822', 0.725,  2.55, 6.5, '2D-5',
           'RAE2822 case 7; same condition as rae2822_experiment ExpCase7'),
    _euler('rae2822_m0730_a3.19', 'rae2822', 0.730,  3.19, 6.5, '2D-6',
           'RAE2822 case 9/10; same condition as rae2822_experiment M0.73'),
]

#: Near-stall startup: measured by the upstream verification harness, a cold
#: start at CFL 1 blows up in ~15 cycles at 12.86 deg.  CFL_START x CFL_RAMP
#: = 1.0 here, i.e. it ramps to the same terminal CFL as every other case but
#: starts 20x below it.
_STALL_START = dict(CFL_START=0.05, CFL_RAMP=20.0, NITFO=1500, NCYC=2000)

#: RANS conditions.  Two turbulence models per condition is a REQUIREMENT, not
#: a bonus: their disagreement is the resolution noise floor of any gate built
#: on this dataset (data request section 4.7), and a tolerance tighter than
#: that floor cannot discriminate in principle.
_RANS_CONDITIONS = [
    # tag                          geom        M      alpha   Re/1e6  x_tr  request
    ('n0012_m0500_a2.00_xtr005', 'naca0012', 0.500,  2.00,  3.0, 0.05, 'R-1',
     'trip 0.05 -- pairs with naca0012_viscous_xfoil xtr005 (XFOIL/RANS/us)',
     None),
    ('n0012_m0500_a2.00_xtr030', 'naca0012', 0.500,  2.00,  3.0, 0.30, 'R-1',
     'trip 0.30 -- pairs with naca0012_viscous_xfoil xtr030', None),
    ('n0012_m0778_a2.03',        'naca0012', 0.778,  2.03,  6.0, None, 'R-2',
     'fully turbulent; same condition as naca0012_experiment M0.778', None),
    ('n0012_m0803_am0.10',       'naca0012', 0.803, -0.10,  6.5, None, 'R-3',
     'fully turbulent; same condition as naca0012_experiment M0.803', None),
    ('n0012_m0352_a12.86',       'naca0012', 0.352, 12.86,  3.0, None, 'R-4',
     'NEAR STALL -- widest model spread, deliberately NOT gate material',
     _STALL_START),
    ('rae2822_m0725_a2.55',      'rae2822',  0.725,  2.55,  6.5, 0.03, 'R-5',
     'trip 0.03 = the Track V x_tr for RAE2822', None),
    ('rae2822_m0730_a3.19',      'rae2822',  0.730,  3.19,  6.5, 0.03, 'R-6',
     'trip 0.03 = the Track V x_tr for RAE2822', None),
]

TURB_MODELS = ((7, 'sst'), (6, 'sa'))

#: one Case per (condition, turbulence model)
RANS_CASES = []
for (tag, geom, mach, alpha, re_mil, x_tr, request, note,
     solver) in _RANS_CONDITIONS:
    for ivisc, _name in TURB_MODELS:
        RANS_CASES.append(_rans(tag, geom, mach, alpha, re_mil, ivisc,
                                request, x_tr=x_tr, note=note,
                                solver=solver))


def dataset_dir(case: R.Case) -> Path:
    geom = 'naca0012' if case.geometry == 'naca0012' else 'rae2822'
    return HERE / f'{case.model}_{geom}'


# ---------------------------------------------------------------------------
#  one run
# ---------------------------------------------------------------------------

FORCE_COLUMNS = [
    'request', 'case', 'geometry', 'model', 'turb_model', 'level',
    'mach', 'alpha_deg', 're_chord', 'x_tr', 'x_tr_actual',
    'nj', 'nk', 'n_surface_cells', 'h1_wall', 'yplus_avg', 'yplus_max',
    'cl', 'cd', 'cd_pressure', 'cd_friction', 'cm_quarter_chord',
    'resid_final', 'resid_decades', 'ncyc_total', 'recipe',
    'status', 'wall_s',
    'keywords', 'laminar_echo', 'forcezero_echo', 'note',
]


#: set by main() before the solver phase; the retry path needs it and the
#: thread pool's payload does not carry it.
geom_files_global: dict = {}


def case_dir(case: R.Case, level: R.Level, workroot: Path) -> Path:
    return workroot / f'{case.model}_{case.turb_model}' / case.tag / level.name


def _recipe_tag(case: R.Case, recipe: dict, rescued: bool = False) -> str:
    """A tag that identifies the startup actually used.

    ★ The ``_fallback`` suffix is load-bearing, not decoration.  ``R-4``'s
    designed near-stall override and ``RANS_FALLBACK`` carry the SAME knob
    values (CFL 0.05 x 20, NITFO 1500, NCYC 2000), so a tag built from the
    knobs alone cannot tell "this startup was chosen for this case" from "this
    leg diverged and was retried" -- and the README promises the recipe column
    distinguishes them.  Measured on this dataset the fallback fired ZERO
    times, so nothing is mislabelled today; the suffix is what keeps that
    claim true when it does fire.
    """
    base = dict(recipe, **case.solver)
    tag = (f"cfl{base['CFL_START']:g}x{base['CFL_RAMP']:g}"
           f"_nitfo{base['NITFO']}_ncyc{base['NCYC']}_mseq{base['MSEQ']}")
    return tag + '_fallback' if rescued else tag


def build_one(case: R.Case, level: R.Level, workroot: Path, geom_files: dict,
              rebuild: bool):
    """Grid + deck for one run.  Returns (grid record, needs_solve).

    ★ MUST run on the main thread: gmsh installs a SIGINT handler in
    ``initialize()``, and ``signal.signal`` raises outside the main thread.
    That is why building and solving are two phases here instead of one
    parallel map -- only the solver runs concurrently.
    """
    d = case_dir(case, level, workroot)
    done = ((d / 'cfl3d.out').is_file() and (d / 'cfl3d.prt').is_file()
            and (d / 'grid_record.csv').is_file())
    if done and not rebuild:
        return _reload_grid_record(d), False
    rec = R.build_case(case, level, d, geom_files[case.geometry])
    rec['recipe'] = _recipe_tag(case, R.default_solver(case))
    _store_grid_record(d, rec)
    return rec, True


def run_one(case: R.Case, level: R.Level, workroot: Path, rec: dict,
            needs_solve: bool, solver: Path) -> dict:
    d = case_dir(case, level, workroot)
    tag = f'{case.tag}/{case.turb_model}/{level.name}'
    run = (R.run_case(d, solver) if needs_solve
           else dict(status='cached', wall_s=0.0))

    # ★ Retry ONCE with the gentler startup.  It fires only where the first
    # attempt already diverged, so it cannot move a result that would have been
    # used; the recipe actually used lands in the `recipe` column, and the row's
    # `resid_decades` is what says whether the rescue produced a converged
    # solution or only a snapshot (see RANS_FALLBACK's note).
    if run['status'] == 'diverged' and case.ivisc != 0:
        print(f'    {tag:52s} diverged -> retrying with the gentle startup')
        rec = R.build_case(case, level, d, geom_files_global[case.geometry],
                           recipe=R.RANS_FALLBACK)
        rec['recipe'] = _recipe_tag(case, R.RANS_FALLBACK,
                                    rescued=True)
        _store_grid_record(d, rec)
        retry = R.run_case(d, solver)
        retry['wall_s'] += run['wall_s']
        run = retry

    row = dict.fromkeys(FORCE_COLUMNS, '')
    row.update(
        request=case.request, case=case.tag,
        geometry=case.geometry, model=case.model,
        turb_model=case.turb_model, level=level.name,
        mach=f'{case.mach:.4f}', alpha_deg=f'{case.alpha:.4f}',
        re_chord=f'{case.re_mil * 1e6:.4e}',
        x_tr='' if case.x_tr is None else f'{case.x_tr:.4f}',
        x_tr_actual=rec.get('x_tr_actual', ''),
        nj=rec['nj'], nk=rec['nk'],
        n_surface_cells=rec['jte2'] - rec['jte1'],
        h1_wall=f'{rec["h1"]:.4e}',
        recipe=rec.get('recipe', ''),
        status=run['status'], wall_s=f'{run["wall_s"]:.1f}',
        note=case.note,
    )
    if run['status'] not in ('ok', 'cached'):
        print(f'  ! {tag:52s} {run["status"]}')
        return row

    f = R.read_forces(d / 'cfl3d.out')
    h = R.read_history(d / 'clcd_total.dat')
    row.update(
        cl=f'{f.get("CL", float("nan")):.6f}',
        cd=f'{f.get("CD", float("nan")):.6f}',
        cd_pressure=f'{f.get("CDp", float("nan")):.6f}',
        cd_friction=f'{f.get("CDv", float("nan")):.6f}',
        cm_quarter_chord=f'{f.get("CMZ", float("nan")):.6f}',
        yplus_avg='' if 'yplus_avg' not in f else f'{f["yplus_avg"]:.3f}',
        yplus_max='' if 'yplus_max' not in f else f'{f["yplus_max"]:.3f}',
        laminar_echo=f.get('laminar_echo', ''),
        forcezero_echo=f.get('forcezero_echo', ''),
        keywords=f.get('keyword_echo', ''),
    )
    row.update(R.convergence(h))

    # ★ Every requested keyword must appear in the solver's own echo of the
    # keyword block.  A keyword the build silently dropped -- or one this build
    # does not recognise -- would otherwise be recorded in the CSV as if it had
    # been in force.  Same "mentioned != used" family as the transition check
    # below and as F06's shock reference.
    for k in case.keywords:
        if k not in row['keywords']:
            raise RuntimeError(
                f'{tag}: keyword {k!r} was requested but does not appear in '
                f'the solver echo ({row["keywords"]!r})')

    # transition caliber must be verified, not trusted: a silently unpatched
    # deck runs fully turbulent while the CSV claims a trip location.
    if case.x_tr is not None and row['laminar_echo'] in ('', 'none'):
        raise RuntimeError(
            f'{tag}: x_tr = {case.x_tr} was requested but cfl3d.out echoes no '
            'laminar region -- the deck patch did not take effect')
    if case.x_tr is None and row['laminar_echo'] not in ('', 'none'):
        raise RuntimeError(f'{tag}: fully turbulent requested but the solver '
                           f'echoes a laminar region {row["laminar_echo"]}')

    up, lo = R.split_surfaces(*R.read_wall_cp(d / 'cfl3d.prt'))
    out = dataset_dir(case)
    out.mkdir(parents=True, exist_ok=True)
    _write_cp(out / f'cp_{case.tag}_{case.turb_model}_{level.name}.csv', up, lo)

    print(f'    {tag:52s} cl={row["cl"]:>9s} cd={row["cd"]:>9s} '
          f'|R|={row["resid_final"] or "n/a":>12s} {run["wall_s"]:6.1f} s')
    return row


def _store_grid_record(d: Path, rec: dict):
    with open(d / 'grid_record.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(rec.keys())
        w.writerow(rec.values())


def _reload_grid_record(d: Path) -> dict:
    with open(d / 'grid_record.csv', newline='') as f:
        rows = list(csv.DictReader(f))
    r = rows[0]
    for k in ('nj', 'nk', 'jte1', 'jte2', 'jlamlo', 'jlamhi'):
        r[k] = int(r[k])
    r['h1'] = float(r['h1'])
    return r


def _write_cp(path: Path, up: dict, lo: dict):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['x_c', 'y_c', 'cp', 'mach_local', 'surface'])
        for side, s in (('upper', up), ('lower', lo)):
            order = np.argsort(s['x'])
            for i in order:
                w.writerow([f'{s["x"][i]:.6f}', f'{s["y"][i]:.6f}',
                            f'{s["cp"][i]:.6f}', f'{s["mach"][i]:.6f}',
                            side])


# ---------------------------------------------------------------------------
#  derived tables
# ---------------------------------------------------------------------------

def _read_cp(path: Path):
    x = {'upper': [], 'lower': []}
    cp = {'upper': [], 'lower': []}
    with open(path, newline='') as f:
        for r in csv.DictReader(f):
            x[r['surface']].append(float(r['x_c']))
            cp[r['surface']].append(float(r['cp']))
    return ({k: np.asarray(v) for k, v in x.items()},
            {k: np.asarray(v) for k, v in cp.items()})


SHOCK_COLUMNS = ['case', 'model', 'turb_model', 'level', 'mach', 'surface',
                 'has_shock', 'x_shock', 'n_cells', 'monotone',
                 'cp_min', 'cp_pre_shock', 'cp_post_shock', 'cp_critical']


def derive_shock(dsdir: Path, rows: list[dict]) -> list[dict]:
    """Shock table, computed from the committed cp CSVs with OUR OWN operator.

    ★ The extraction operator is `pyfp3d.post.shock.shock_metrics`, the very
    function the pyFP3D side is read with.  Using CFL3D's own shock definition
    on one side and ours on the other would compare two different quantities --
    the error family this project has logged six times.  The consequence to
    keep in mind: these columns are DERIVED, so they move if that operator
    changes; the primary data are the cp_*.csv curves, and `--derive-only`
    rebuilds this table from them without touching the solver.

    cp_pre_shock / cp_post_shock are sampled 0.05 c either side of the
    crossing.  The post-shock column exists to be RECORDED, not gated: full
    potential is isentropic and irrotational, so a real model difference is
    EXPECTED behind the shock (data request ruling 1).
    """
    out = []
    for row in rows:
        if row['status'] not in ('ok', 'cached') or not row['cl']:
            continue
        p = dsdir / f'cp_{row["case"]}_{row["turb_model"]}_{row["level"]}.csv'
        if not p.is_file():
            continue
        xs, cps = _read_cp(p)
        m_inf = float(row['mach'])
        for surface in ('upper', 'lower'):
            x, cp = xs[surface], cps[surface]
            sm = shock_metrics(x, cp, m_inf)
            pre = post = ''
            if sm['has_shock']:
                xsh = sm['x_shock']
                pre = f'{float(np.interp(xsh - 0.05, x, cp)):.6f}'
                if xsh + 0.05 <= x.max():
                    post = f'{float(np.interp(xsh + 0.05, x, cp)):.6f}'
            out.append(dict(
                case=row['case'], model=row['model'],
                turb_model=row['turb_model'], level=row['level'],
                mach=row['mach'], surface=surface,
                has_shock=int(sm['has_shock']),
                x_shock='' if not sm['has_shock'] else f'{sm["x_shock"]:.6f}',
                n_cells=sm['n_cells'], monotone=int(sm['monotone']),
                cp_min=f'{sm["cp_min"]:.6f}',
                cp_pre_shock=pre, cp_post_shock=post,
                cp_critical=f'{cp_critical(m_inf):.6f}'))
    return out


GC_QUANTITIES = ('cl', 'cd', 'cd_pressure', 'cd_friction',
                 'cm_quarter_chord')
GC_COLUMNS = (['case', 'turb_model', 'quantity']
              + ['L1', 'L2', 'L3']
              + ['delta_L2_L3', 'rel_delta_L2_L3'])


def derive_grid_convergence(rows: list[dict], shock: list[dict]) -> list[dict]:
    """Per case: the value on each rung and the delta between the two finest.

    ★ That delta IS the error bar the data request asks for (item 3): a single
    grid is one number with no uncertainty, and a reference without an
    uncertainty cannot carry a gate.  It is reported as an interval, not
    extrapolated -- three rungs at ratio ~sqrt(2) in every direction do not by
    themselves establish an asymptotic order on a shocked solution.
    """
    by_case: dict[tuple, dict] = {}
    for row in rows:
        if row['status'] not in ('ok', 'cached'):
            continue
        for q in GC_QUANTITIES:
            if row[q] == '':
                continue
            by_case.setdefault((row['case'], row['turb_model'], q), {})[
                row['level']] = float(row[q])
    for s in shock:
        if s['x_shock'] == '':
            continue
        key = (s['case'], s['turb_model'], f'x_shock_{s["surface"]}')
        by_case.setdefault(key, {})[s['level']] = float(s['x_shock'])

    out = []
    for (case, turb, q), vals in sorted(by_case.items()):
        r = dict(case=case, turb_model=turb, quantity=q,
                 L1='', L2='', L3='', delta_L2_L3='', rel_delta_L2_L3='')
        for lev, v in vals.items():
            r[lev] = f'{v:.6f}'
        if 'L2' in vals and 'L3' in vals:
            d = vals['L3'] - vals['L2']
            r['delta_L2_L3'] = f'{d:.6f}'
            scale = max(abs(vals['L3']), 1e-12)
            r['rel_delta_L2_L3'] = f'{d / scale:.4%}'
        out.append(r)
    return out


TS_COLUMNS = ['case', 'level', 'quantity', 'sst', 'sa', 'spread',
              'rel_spread', 'coverage']


def derive_turbulence_spread(rows: list[dict],
                             shock: list[dict]) -> list[dict]:
    """|SST - SA| per quantity = the NOISE FLOOR of any gate on this dataset.

    ★ Data request section 4.7: the turbulence model is a free choice, and the
    disagreement between two accepted models can exceed the quantity a gate
    wants to bound.  A criterion tighter than this spread is UNDEFINED, not
    failed -- the same logic as the A4 2.5 % input band, and the reason GS4.1
    round 15 demoted E-CF from a gate to RECORDED.
    """
    vals: dict[tuple, dict] = {}
    for row in rows:
        if row['status'] not in ('ok', 'cached'):
            continue
        for q in GC_QUANTITIES:
            if row[q] == '':
                continue
            vals.setdefault((row['case'], row['level'], q), {})[
                row['turb_model']] = float(row[q])
    for s in shock:
        if s['x_shock'] == '':
            continue
        vals.setdefault((s['case'], s['level'],
                         f'x_shock_{s["surface"]}'), {})[
            s['turb_model']] = float(s['x_shock'])

    out = []
    for (case, level, q), v in sorted(vals.items()):
        row = dict(case=case, level=level, quantity=q,
                   sst='' if 'sst' not in v else f'{v["sst"]:.6f}',
                   sa='' if 'sa' not in v else f'{v["sa"]:.6f}',
                   spread='', rel_spread='', coverage='both')
        if 'sst' not in v or 'sa' not in v:
            # ★ Emit the row anyway with coverage = which model is MISSING.
            # A spread that is simply absent from the file reads as "no
            # disagreement here"; a spread that says UNDEFINED reads as "this
            # gate has no measured noise floor at this rung".  The project has
            # paid for that distinction once already -- M1c's seed-0 leg failed
            # on COVERAGE, not on tolerance, and no tolerance could cure it.
            row['coverage'] = ('sst_only' if 'sa' not in v else 'sa_only')
            row['spread'] = 'UNDEFINED'
            out.append(row)
            continue
        spread = abs(v['sst'] - v['sa'])
        scale = max(abs(v['sst']), abs(v['sa']), 1e-12)
        row['spread'] = f'{spread:.6f}'
        row['rel_spread'] = f'{spread / scale:.4%}'
        out.append(row)
    return out


def _write_csv(path: Path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(columns), extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    try:
        shown = path.relative_to(HERE.parents[2])
    except ValueError:          # a work dir outside the repo
        shown = path
    print(f'  -> {shown}  ({len(rows)} rows)')


def _read_forces(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------------------
#  driver
# ---------------------------------------------------------------------------

def emit(cases: list[R.Case], rows: list[dict]):
    """Write forces.csv + the derived tables into each dataset directory."""
    for dsdir in sorted({dataset_dir(c) for c in cases}):
        mine = [r for r in rows
                if dataset_dir_of_row(r) == dsdir]
        if not mine:
            continue
        mine.sort(key=lambda r: (r['case'], r['turb_model'], r['level']))
        _write_csv(dsdir / 'forces.csv', FORCE_COLUMNS, mine)
        shock = derive_shock(dsdir, mine)
        _write_csv(dsdir / 'shock.csv', SHOCK_COLUMNS, shock)
        _write_csv(dsdir / 'grid_convergence.csv', GC_COLUMNS,
                   derive_grid_convergence(mine, shock))
        if mine[0]['model'] == 'rans':
            _write_csv(dsdir / 'turbulence_spread.csv', TS_COLUMNS,
                       derive_turbulence_spread(mine, shock))


def dataset_dir_of_row(row: dict) -> Path:
    geom = 'naca0012' if row['geometry'] == 'naca0012' else 'rae2822'
    return HERE / f'{row["model"]}_{geom}'


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--set', choices=('euler', 'rans', 'all'), default='all')
    p.add_argument('--levels', nargs='*', default=None,
                   help='subset of L1 L2 L3 (default: all three)')
    p.add_argument('--cases', nargs='*', default=None,
                   help='subset of case tags')
    p.add_argument('--jobs', type=int, default=4)
    p.add_argument('--workdir', default=None,
                   help='run directory (default tools/cfl3d_work, gitignored)')
    p.add_argument('--rebuild', action='store_true',
                   help='re-run cases that already have output')
    p.add_argument('--derive-only', action='store_true',
                   help='rebuild the derived CSVs from the committed data')
    a = p.parse_args(argv)

    cases = ((EULER_CASES if a.set in ('euler', 'all') else [])
             + (RANS_CASES if a.set in ('rans', 'all') else []))
    if a.cases:
        cases = [c for c in cases if c.tag in a.cases]
    if not cases:
        sys.exit('no cases selected')

    if a.derive_only:
        rows = []
        for dsdir in sorted({dataset_dir(c) for c in cases}):
            rows += _read_forces(dsdir / 'forces.csv')
        if not rows:
            sys.exit('--derive-only needs a committed forces.csv')
        emit(cases, rows)
        return 0

    solver = R.find_solver()
    workroot = Path(a.workdir) if a.workdir else (
        HERE.parents[2] / 'tools' / 'cfl3d_work')
    print(f'solver   : {solver}')
    print(f'work dir : {workroot}')
    print(f'load avg : {os.getloadavg()}')

    global geom_files_global
    geom_files = {}
    for kind in sorted({c.geometry for c in cases}):
        geom_files[kind] = R.write_geometry(
            kind, workroot / 'geometry' / f'{kind}.dat')
        print(f'geometry : {kind} <- pyfp3d.meshgen.planar '
              f'({geom_files[kind]})')
    geom_files_global = geom_files

    work = []
    for c in cases:
        ladder = R.EULER_LEVELS if c.ivisc == 0 else R.RANS_LEVELS
        for lev in ladder:
            if a.levels and lev.name not in a.levels:
                continue
            work.append((c, lev))

    # phase 1 -- grids, SERIAL on the main thread (gmsh + signal handlers)
    t0 = time.perf_counter()
    built = []
    for c, lev in work:
        rec, needs = build_one(c, lev, workroot, geom_files, a.rebuild)
        built.append((c, lev, rec, needs))
    n_solve = sum(1 for *_, needs in built if needs)
    print(f'\n{len(work)} run(s): {n_solve} to solve, '
          f'{len(work) - n_solve} cached; grids built in '
          f'{time.perf_counter() - t0:.0f} s\n')

    # phase 2 -- solver runs, parallel (cfl3d_seq is single-threaded)
    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        rows = list(pool.map(
            lambda t: run_one(t[0], t[1], workroot, t[2], t[3], solver),
            built))
    print(f'\nsolver wall {time.perf_counter() - t0:.0f} s\n')

    # append/merge with whatever is already committed so a partial run does
    # not silently truncate the dataset
    merged: dict[tuple, dict] = {}
    for dsdir in sorted({dataset_dir(c) for c in cases}):
        for r in _read_forces(dsdir / 'forces.csv'):
            merged[(r['case'], r['turb_model'], r['level'])] = r
    for r in rows:
        merged[(r['case'], r['turb_model'], r['level'])] = r
    emit(cases, list(merged.values()))

    bad = [r for r in rows if r['status'] not in ('ok', 'cached')]
    if bad:
        print(f'{len(bad)} run(s) did not finish cleanly:')
        for r in bad:
            print(f"  {r['case']}/{r['turb_model']}/{r['level']}: "
                  f"{r['status']}")
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
