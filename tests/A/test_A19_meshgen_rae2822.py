"""GV5.2 meshgen point-set path tests (pre-registration Sec.7).

Pure-python, no gmsh: the committed RAE2822 ordinate load, the PCHIP
resample convention, the TE wedge measure, and the Cp-compare helpers
on BOTH committed experiment layouts (two-zone P1 / wrapped P2). Mesh
generation itself is validated by the gate's band (a) + committed stats.
"""
import importlib.util
import os

import numpy as np
import pytest

from pyfp3d.meshgen.planar import (
    load_airfoil_ordinates,
    pointset_airfoil_coordinates,
    te_wedge_angle_deg,
)

from tests.conftest import REPO_ROOT
ROOT = str(REPO_ROOT)   #: ★ 与目录深度无关（2026-08-24 重编号）
DAT = os.path.join(ROOT, "cases", "meshes", "rae2822_2.5d", "rae2822.dat")
EXP_DIR = os.path.join(ROOT, "cases", "reference_data", "rae2822_experiment")
EXP_P1 = os.path.join(
    EXP_DIR, "ExpCase7_RAE2822_M0.725_AoA2.55_Rec6.5e6.dat")
EXP_P2 = os.path.join(
    EXP_DIR, "Expe_RAE2822_M0.73_AoA3.19_Rec6.5e6.dat")
RUN = os.path.join(ROOT, "bench", "studies", "v5_2_rae2822", "run.py")


def _load_runner():
    """The GV5.2 runner as a module (the v5-gate runner-import idiom)."""
    if not os.path.exists(RUN):
        pytest.skip("GV5.2 runner not yet present")
    spec = importlib.util.spec_from_file_location("gv5_2_run", RUN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_ordinates_load():
    x, z_lo, z_up = load_airfoil_ordinates(DAT)
    assert x.size == 65
    assert x[0] == 0.0 and x[-1] == 1.0
    assert np.all(np.diff(x) > 0.0)
    # Cook Table 6.1 spot values: max upper ordinate at x = 0.45099
    i = int(np.argmax(z_up))
    assert abs(x[i] - 0.45099) < 1e-12 and abs(z_up[i] - 0.06286) < 1e-12
    assert z_up[0] == 0.0 and z_lo[0] == 0.0
    assert z_up[-1] == 0.0 and z_lo[-1] == 0.0
    # the positive-down lower convention (physical z_lower = -tabulated):
    # thickness 12.1% @ x/c 0.379, camber 1.3% @ 0.757 (airfoiltools sig.)
    thick = z_up + z_lo          # = z_up - (-z_lo)
    camber = 0.5 * (z_up - z_lo)  # physical: (z_up + z_lo_phys)/2
    assert abs(thick[int(np.argmax(thick))] - 0.121) < 0.002
    assert abs(camber[int(np.argmax(camber))] - 0.013) < 0.001
    # physical section: strictly positive thickness inside (LE/TE sharp)
    assert np.all(thick[1:-1] > 0.0)


def test_resample_convention_and_bounds():
    x, z_lo, z_up = load_airfoil_ordinates(DAT)
    z_lo = -z_lo  # physical z (the generator's convention fix)
    coords = pointset_airfoil_coordinates(x, z_lo, z_up, n_half=120)
    # the closed-polyline convention: first/last exactly (1, 0)
    assert tuple(coords[0]) == (1.0, 0.0)
    assert tuple(coords[-1]) == (1.0, 0.0)
    # unique LE point
    assert int(np.argmin(coords[:, 0])) > 0
    assert coords[:, 0].min() == 0.0 and coords[:, 0].max() == 1.0
    # PCHIP shape preservation: no overshoot beyond the data extrema
    z = coords[1:-1, 1]
    assert z.max() <= max(z_up.max(), abs(z_lo.min())) + 1e-12
    assert z.min() >= z_lo.min() - 1e-12
    # the physical contour does NOT self-intersect (the 2026-07-24
    # convention erratum: the un-negated lower folds over the upper and
    # the 2-D gmsh build never terminates)
    upper = coords[: int(np.argmin(coords[:, 0])) + 1]
    lower = coords[int(np.argmin(coords[:, 0])):]
    assert upper[:, 1].min() >= -1e-12  # upper stays at/above datum-ish
    assert lower[:, 1].max() <= 2.0e-3  # lower stays below save the TE poke
    # cosine clustering: the end spacing is much finer than mid-chord
    xs = np.unique(coords[:, 0])
    d = np.diff(xs)
    assert d.min() < 0.05 * d[len(d) // 2]


def test_te_wedge_angle():
    x, z_lo, z_up = load_airfoil_ordinates(DAT)
    wedge = te_wedge_angle_deg(x, z_lo, z_up)
    # the GV5.2 prior (~13 deg); the A4 guard (~6 deg) must stay quiet
    assert 10.0 <= wedge <= 16.0
    assert wedge > 6.0


def test_cp_compare_helpers_both_layouts():
    mod = _load_runner()
    up1, lo1 = mod.load_experiment_cp(EXP_P1)
    assert up1[0].size == 53 and lo1[0].size == 46
    # the P1 upper shock bracket (committed): the windowed compression
    # max lands at the bracket edge station
    xs, cp = up1
    assert 0.5 <= mod.shock_x(xs, cp) <= 0.575
    up2, lo2 = mod.load_experiment_cp(EXP_P2)
    assert up2[0].size + lo2[0].size == 67
    assert 0.5 <= mod.shock_x(up2[0], up2[1]) <= 0.6


def test_cp_rms_and_shock_helpers():
    mod = _load_runner()
    # identical curves -> zero RMS; a unit offset -> RMS 1
    xe = np.linspace(0.0, 1.0, 11)
    cf = np.linspace(0.0, 1.0, 101)
    assert mod.cp_rms(cf, np.zeros_like(cf), xe, np.zeros_like(xe)) == 0.0
    assert mod.cp_rms(cf, np.ones_like(cf), xe,
                      np.zeros_like(xe)) == pytest.approx(1.0)
    # shock_x picks the windowed compression-jump station
    xs = np.linspace(0.0, 1.0, 201)
    cp = -1.0 + 0.8 / (1.0 + np.exp(-(xs - 0.53) / 0.005))
    assert abs(mod.shock_x(xs, cp) - 0.53) < 0.01
    # an expansion jump at the same place is NOT the shock branch
    cp_exp = -0.2 - 0.8 / (1.0 + np.exp(-(xs - 0.53) / 0.005))
    assert abs(mod.shock_x(xs, cp_exp) - 0.53) > 0.01
