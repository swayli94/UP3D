"""
Track B / B18: wing-body transonic (M0.84) -- the conforming path reaches it,
the level-set path is junction-limited.

The heavy wing-body M0.5->0.84 ramps live in the B18 demo
(cases/demo/b18_wingbody_transonic) under PYFP3D_TRANSONIC_GATES. These ungated
tests lock the WIRING that the demo relies on, cheaply on the committed 2.5D NACA
mesh: both Mach-continuation drivers exist and pass the transonic knobs through,
and the level-set ramp returns the HONEST ceiling fields (target_reached /
m_last_converged / m_final) that the demo uses to record a non-converged ramp
rather than censusing a state that never reached the target (the P13/G13.3
erratum -- never census a state whose ramp did not reach the target).
"""

import inspect
import os
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.solve.newton import solve_newton_transonic
from pyfp3d.solve.newton_ls import solve_multivalued_newton_transonic
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
ALPHA = 2.0


def _naca(level="coarse"):
    path = M0_DIR / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    return read_mesh(path)


def _naca_mvop(mesh):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
                       direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


# ---------------------------------------------------------------------------
# 1. both transonic drivers exist with the target knob (conforming target is
#    `m_inf`, the level-set target is `m_target`)
# ---------------------------------------------------------------------------
def test_both_transonic_drivers_present():
    pc = inspect.signature(solve_newton_transonic).parameters
    pl = inspect.signature(solve_multivalued_newton_transonic).parameters
    assert "m_inf" in pc and "m_start" in pc
    assert "m_target" in pl and "m_start" in pl
    # the LS ramp forwards per-level Newton kwargs (farfield/farfield_aux) via
    # **newton_kw -- the demo relies on farfield="freestream", farfield_aux=
    # "pin_gamma" reaching each level solve
    assert any(p.kind == p.VAR_KEYWORD for p in pl.values())


# ---------------------------------------------------------------------------
# 2. the LS ramp returns the honest ceiling fields on a REACHED target
# ---------------------------------------------------------------------------
def test_ls_transonic_honest_fields_reached():
    mesh = _naca()
    r = solve_multivalued_newton_transonic(
        _naca_mvop(mesh), mesh, m_target=0.72, alpha_deg=ALPHA, m_start=0.70,
        dm=0.05, farfield="vortex", freeze_tol=1e-6, n_seed=20, n_newton_max=40)
    assert set(("target_reached", "m_last_converged", "m_final", "levels")) <= set(r)
    if r["target_reached"]:
        assert abs(r["m_final"] - 0.72) < 1e-9
        assert abs(r["m_last_converged"] - 0.72) < 1e-9
    # levels is a per-level record list; the schedule ends at the target
    assert r["mach_schedule"][-1] == pytest.approx(0.72)


# ---------------------------------------------------------------------------
# 3. honest fields on a NON-converged ramp: a 1-step budget cannot reach a
#    transonic target => target_reached is False and we must NOT claim the target
# ---------------------------------------------------------------------------
def test_ls_transonic_honest_fields_not_reached():
    mesh = _naca()
    r = solve_multivalued_newton_transonic(
        _naca_mvop(mesh), mesh, m_target=0.80, alpha_deg=ALPHA, m_start=0.75,
        dm=0.05, dm_min=0.02, farfield="vortex", freeze_tol=1e-6,
        n_seed=1, n_newton_max=1)
    assert r["target_reached"] is False
    # the returned state is at m_final (the level it died on), NOT the target;
    # m_last_converged is None (nothing converged) or < target
    assert r["m_final"] <= 0.80
    assert r["m_last_converged"] is None or r["m_last_converged"] < 0.80


# ---------------------------------------------------------------------------
# 4. farfield_aux="pin_gamma" is accepted by the LS ramp (inert on vortex; the
#    demo uses it on freestream) -- the knob reaches the per-level solve
# ---------------------------------------------------------------------------
def test_ls_transonic_accepts_pin_gamma():
    mesh = _naca()
    r = solve_multivalued_newton_transonic(
        _naca_mvop(mesh), mesh, m_target=0.72, alpha_deg=ALPHA, m_start=0.70,
        dm=0.05, farfield="vortex", farfield_aux="pin_gamma", freeze_tol=1e-6,
        n_seed=20, n_newton_max=40)
    assert np.all(np.isfinite(r["phi_ext"]))
