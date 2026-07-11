"""P10/G10.2: level-adaptive intermediate continuation tolerance
(solve/newton.py::solve_newton_transonic `intermediate_tol`).

Suite-cheap coverage only (the A/B acceptance runs on the G8.1/G8.2
recipes are one-shot gated evidence, cases/demo/p10_newton_usability/):

  - default path instrumentation: `accept_reason` == "tol" on a normally
    converged solve (additive key, default behaviour bit-identical -- the
    whole existing newton test set is the numeric lock).
  - the opt-in knob on a subsonic Mach ramp: intermediate levels accept
    loosely (accept_reason recorded, no freeze machinery), the FINAL
    level keeps tol 1e-10 and converges to the same discrete solution as
    the default run (unique subsonic solution => Gamma matches to
    ~roundoff-amplified level), and the loose run spends no more Newton
    steps than the default.
"""

import numpy as np

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.newton import solve_newton_lifting, solve_newton_transonic

from .conftest import REPO_ROOT

UPWIND_C = 1.5
M_CRIT = 0.95


def _coarse():
    mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
                     / "coarse.msh")
    return cut_wake(mesh)


def test_accept_reason_default_tol():
    mc, wc = _coarse()
    r = solve_newton_lifting(mc, wc, m_inf=0.5, alpha_deg=2.0,
                             upwind_c=UPWIND_C, m_crit=M_CRIT)
    assert r["converged"]
    assert r["accept_reason"] == "tol"


def test_g102_intermediate_tol_subsonic_ramp():
    mc, wc = _coarse()
    kw = dict(m_inf=0.5, alpha_deg=2.0, m_start=0.30, dm=0.10,
              upwind_c=UPWIND_C, m_crit=M_CRIT,
              newton_kw=dict(precond="direct"))
    r_ref = solve_newton_transonic(mc, wc, **kw)
    r_ad = solve_newton_transonic(mc, wc, intermediate_tol=1e-5, **kw)

    for r in (r_ref, r_ad):
        assert r["converged"]
        assert r["residual_history"][-1] < 1e-10
        assert r["n_limited"] == 0 and r["n_floored"] == 0
    # final level always runs the strict criterion
    assert r_ad["level_results"][-1]["accept_reason"] == "tol"
    # every accepted level records its criterion
    assert all(lr["accept_reason"] is not None
               for lr in r_ad["level_results"])
    # same discrete solution (the final level re-tightens to 1e-10)
    assert np.max(np.abs(r_ad["gamma"] - r_ref["gamma"])) < 1e-6
    # the loose intermediate levels cannot cost MORE steps in total
    steps = lambda r: sum(lr["n_newton"] for lr in r["level_results"])
    assert steps(r_ad) <= steps(r_ref)
