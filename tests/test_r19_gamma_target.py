"""R19 -- prescribed circulation on the Newton path (`gamma_target`).

Binding text: docs/dev_phase_five/20260823-1700-r19-prereg.md.

The Picard path has had `gamma_fixed` since phase 1 and `newton.py` had no equivalent
(discipline #9). The change generalises the B31 pin's target from zero to a prescribed
value on BOTH Kutta rows -- the probe row (the default estimator) and the pressure-path
blend. The Jacobian is unchanged in both, because the target is a constant, so the
existing FD locks still cover it; these tests lock the parts FD cannot see.

★ Semantics worth stating: with tip_taper = 0 the Kutta condition is NOT enforced, so a
pinned state is a DIAGNOSTIC probe (A2's fixed-Gamma discriminator), not a physical
solution. That is the intended use.
"""
import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting

MESH = "cases/meshes/naca0012_2.5d/xcoarse.msh"
#: the M1 gate recipe, imported in spirit but written here so the test does not
#: depend on a bench module
KW = dict(upwind_c=1.0, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
          precond="direct", direct_refactor_every=4, n_newton_max=80,
          n_picard_seed=0)
M_INF, ALPHA = 0.80, 1.25


@pytest.fixture(scope="module")
def case():
    mc, wc = cut_wake(read_mesh(MESH))
    return mc, wc, np.zeros(wc.n_stations)


def _cl(mc, phi):
    dz = float(np.ptp(mc.nodes[:, 2]))
    return float(wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], np.asarray(phi, float),
        alpha_deg=ALPHA, m_inf=M_INF, s_ref=dz)["cl"])


def test_default_is_the_same_code_path(case):
    """gamma_target=None must be BIT-IDENTICAL to omitting it -- the guard is written
    as `if gamma_target is None` precisely so this needs no floating-point argument."""
    mc, wc, _ = case
    a = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **KW)
    b = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             gamma_target=None, **KW)
    assert np.array_equal(np.asarray(a["phi"], float), np.asarray(b["phi"], float))


def test_zero_target_is_bit_identical_to_no_target(case):
    """★ This one exercises the CHANGED expression: with tip_taper = 0 the pin is live,
    and a target of exactly zero must reproduce the pre-R19 pin bit-for-bit."""
    mc, wc, z = case
    a = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, tip_taper=z, **KW)
    b = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, tip_taper=z,
                             gamma_target=z, **KW)
    assert np.array_equal(np.asarray(a["phi"], float), np.asarray(b["phi"], float))
    assert float(np.asarray(a["gamma"])[0]) == float(np.asarray(b["gamma"])[0])


@pytest.mark.parametrize("target", [0.10, 0.25])
def test_gamma_hits_the_prescribed_value(case, target):
    """The pin is a unit-slope row, so the converged Gamma must BE the target. Measured
    bit-exact on this case; the band is 1e-10 so a solver change cannot silently
    loosen it into 'approximately pinned'."""
    mc, wc, z = case
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, tip_taper=z,
                             gamma_target=np.full(wc.n_stations, target), **KW)
    assert bool(r["converged"])
    assert abs(float(np.asarray(r["gamma"])[0]) - target) <= 1e-10


def test_gamma_target_changes_the_lift_monotonically(case):
    """A prescribed circulation must move cl_p, and in the right direction -- otherwise
    the pin is wired to something that does not reach the flow."""
    mc, wc, z = case
    cls = []
    for g in (0.05, 0.15, 0.25):
        r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, tip_taper=z,
                                 gamma_target=np.full(wc.n_stations, g), **KW)
        cls.append(_cl(mc, r["phi"]))
    assert cls[0] < cls[1] < cls[2], cls


def test_shape_is_validated(case):
    mc, wc, _ = case
    with pytest.raises(ValueError, match="gamma_target must be"):
        solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             tip_taper=np.zeros(wc.n_stations),
                             gamma_target=np.zeros(wc.n_stations + 1), **KW)
