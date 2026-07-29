"""GS1b.3 always-on locks for the entropy-corrected density.

Pre-registered criteria (docs/dev_phase_two/20260729-0700-s1b-entropy-
implementation.md sec 2.1):

  A  the 0-D identity: with rho_s = sigma*rho_isen and sigma = p02/p01 the
     full-potential mass-conservation jump reproduces Rankine-Hugoniot
  C  the Jacobian is FD-exact for the FROZEN-sigma system, WITH the epsilon
     discriminator (an epsilon-independent error is a missing term, not FD
     noise -- the phase-one B19 lesson)
  D  default off is bit-identical

plus the kernel's own contracts (chain product, smeared-shock peak detection,
donor-cycle detection) -- each of the last two caught a real bug during
implementation, which is why they are locked rather than trusted.
"""

import numpy as np
import pytest

from pyfp3d.kernels.entropy import (EntropyOperator, total_pressure_ratio,
                                    transport_sigma)
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.physics.isentropic import (GAMMA, critical_speed_squared,
                                       density_isentropic, q2_at_mach)
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting

from ._tol import assert_rel_close

MESH = "cases/meshes/naca0012_2.5d/coarse.msh"
M_INF, ALPHA = 0.7875, 1.25


def _rh(m1, gamma=GAMMA):
    m2 = m1 * m1
    return ((gamma + 1.0) * m2) / ((gamma - 1.0) * m2 + 2.0)


def _fp_jump(m1, sigma, m_inf=0.80, gamma=GAMMA):
    """Downstream density ratio from FP mass conservation with the corrected
    density rho2 = sigma*rho_isen(q2)."""
    u1 = float(np.sqrt(q2_at_mach(m1, m_inf, gamma)))
    rho1 = float(density_isentropic(u1 * u1, m_inf, gamma))
    mdot = rho1 * u1
    lo, hi = 1e-9, u1
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if sigma * float(density_isentropic(mid * mid, m_inf, gamma)) * mid \
                < mdot:
            lo = mid
        else:
            hi = mid
    u2 = 0.5 * (lo + hi)
    return sigma * float(density_isentropic(u2 * u2, m_inf, gamma)) / rho1


# --------------------------------------------------------------- criterion A
@pytest.mark.parametrize("m1", [1.15, 1.20, 1.30, 1.35, 1.40, 1.50, 1.60])
def test_a_corrected_jump_is_rankine_hugoniot(m1):
    """A: the entropy-corrected FP jump IS R-H, to better than 1e-6 relative.

    This is the whole justification for sigma = p02/p01 with NO further
    exponent: rho0 = p0/(R*T0) and T0 is preserved across the shock, so
    rho02/rho01 = p02/p01 exactly. GS1b.2 first wrote sigma^(1/(gamma-1))
    (the exponent applied twice) -- this test is what makes that
    unrepeatable.
    """
    sigma = total_pressure_ratio(m1)
    assert_rel_close(_fp_jump(m1, sigma), _rh(m1), rtol=1e-6)


@pytest.mark.parametrize("m1", [1.30, 1.40])
def test_a_wrong_exponents_do_not_reproduce_rh(m1):
    """A (negative control): the exponents that LOOK plausible are wrong by
    percents, so criterion A actually discriminates."""
    sigma = total_pressure_ratio(m1)
    target = _rh(m1)
    for bad in (sigma ** (1.0 / (GAMMA - 1.0)), sigma ** (1.0 / GAMMA)):
        rel = abs(_fp_jump(m1, bad) / target - 1.0)
        assert rel > 1e-2, f"exponent {bad} is not discriminated (rel {rel})"


def test_subsonic_sigma_is_exactly_one():
    """No shock => no entropy => sigma == 1.0 bitwise (so a subcritical solve
    is bit-identical whether or not the correction is enabled)."""
    assert total_pressure_ratio(0.99) == 1.0
    assert total_pressure_ratio(1.0) == 1.0


# ------------------------------------------------------- kernel contracts
def test_sigma_is_the_product_along_the_donor_chain():
    mach = np.array([1.30, 1.30, 0.80, 1.25, 0.70])       # two shocks in series
    q2 = np.array([float(q2_at_mach(m, 0.8)) for m in mach])
    ent = EntropyOperator(5)
    sig = ent.sigma(q2, np.array([0, 0, 1, 2, 3]), 0.8).copy()
    want = total_pressure_ratio(1.30) * total_pressure_ratio(1.25)
    assert_rel_close(sig[4], want)
    assert sig[0] == 1.0 and sig[1] == 1.0


def test_peak_mach_is_read_through_a_smeared_shock():
    """The detector must report the PEAK pre-shock Mach, not the last
    supersonic cell's.

    Regression for a measured bug: the first implementation used the donor's own
    Mach and reported sigma_min 0.9968 (M1 ~ 1.14) on a state whose M_max was
    1.366 -- it was reading the artificial-density smearing, not the shock.
    """
    mach = np.array([1.05, 1.37, 1.25, 1.10, 0.95, 0.85])
    q2 = np.array([float(q2_at_mach(m, 0.8)) for m in mach])
    ent = EntropyOperator(6)
    sig = ent.sigma(q2, np.array([0, 0, 1, 2, 3, 4]), 0.8).copy()
    assert_rel_close(ent.m1_max, 1.37)
    assert_rel_close(sig[-1], total_pressure_ratio(1.37))


def test_donor_cycle_is_detected_not_silently_wrong():
    """A donor cycle must report converged=False.

    Regression for a measured bug: with the weaker test "the ancestor pointers
    stopped moving", a 2-cycle settles at A[0]=0, A[1]=1 -- each element becomes
    its own ancestor and looks exactly like a chain root -- while the product
    keeps squaring every round. It reported converged with a corrupted sigma.
    """
    mach = np.array([1.30, 1.30, 0.80, 0.75, 0.70])
    q2 = np.array([float(q2_at_mach(m, 0.8)) for m in mach])
    ent = EntropyOperator(5)
    ent.sigma(q2, np.array([1, 0, 1, 2, 3]), 0.8)        # u(0)=1, u(1)=0
    assert ent.converged is False
    assert ent.n_rounds == ent.n_round


def test_transport_is_thread_order_independent():
    """Same input twice must give bit-identical sigma (the transport writes into
    separate per-round buffers precisely so threading cannot reorder it)."""
    rng = np.random.default_rng(0)
    n = 4000
    s = 1.0 - 0.05 * rng.random(n)
    up = np.maximum(np.arange(n) - 1, 0)
    out = [np.empty(n) for _ in range(2)]
    for o in out:
        transport_sigma(s, up, 24, o, np.empty(n, dtype=np.int64),
                        np.empty(n, dtype=np.int64), np.empty(n))
    assert np.array_equal(out[0], out[1])


# ------------------------------------------ criterion D / the scope boundary
def test_scope_boundary_no_shock_means_bit_identical():
    """★ The SCOPE of the correction, proved rather than surveyed: with no
    supersonic cell anywhere, ON and OFF are bit-identical.

    sigma = p02/p01 is exactly 1.0 at and below M = 1 and the detector only fires
    on a supersonic->subsonic donor transition, so a subcritical solve cannot see
    the correction at all. That is what bounds the re-baselining: only committed
    evidence with a supersonic zone can move.
    """
    mc, wc = cut_wake(read_mesh(MESH))
    kw = dict(m_inf=0.50, alpha_deg=ALPHA, upwind_c=1.5, m_crit=0.95,
              precond="direct", n_newton_max=40)
    off = solve_newton_lifting(mc, wc, entropy_correction=False, **kw)
    on = solve_newton_lifting(mc, wc, entropy_correction=True, **kw)
    assert float(np.sqrt(off["mach2_max"])) < 1.0, \
        "the probe is supercritical -- it cannot test the no-shock null"
    assert np.array_equal(off["phi"], on["phi"])
    assert np.array_equal(off["gamma"], on["gamma"])
    assert on["sigma_min"] == 1.0 and on["n_shock_cells"] == 0


def test_the_switch_actually_switches_at_transonic():
    """The companion to the null: where there IS a shock, ON and OFF must give
    measurably different states (otherwise the flag is decorative)."""
    mc, wc = cut_wake(read_mesh(MESH))
    kw = dict(m_inf=M_INF, alpha_deg=ALPHA, upwind_c=1.5, m_crit=0.95,
              freeze_tol=1e-6, freeze_refresh_max=8, precond="direct",
              direct_refactor_every=4, n_newton_max=80)
    off = solve_newton_lifting(mc, wc, entropy_correction=False, **kw)
    on = solve_newton_lifting(mc, wc, entropy_correction=True, **kw)
    assert off["sigma_min"] is None and on["sigma_min"] < 1.0
    assert abs(float(on["gamma"][0]) / float(off["gamma"][0]) - 1.0) > 1e-3


def test_off_path_is_deterministic():
    """The OFF path must be bit-reproducible run to run, so it stays usable as
    the isentropic reference leg of every ON/OFF comparison."""
    mc, wc = cut_wake(read_mesh(MESH))
    kw = dict(m_inf=M_INF, alpha_deg=ALPHA, upwind_c=1.5, m_crit=0.95,
              freeze_tol=1e-6, freeze_refresh_max=8, precond="direct",
              direct_refactor_every=4, n_newton_max=80,
              entropy_correction=False)
    a = solve_newton_lifting(mc, wc, **kw)
    b = solve_newton_lifting(mc, wc, **kw)
    assert np.array_equal(a["phi"], b["phi"])
    assert a["n_newton"] == b["n_newton"]


def test_upstream_map_must_be_built_before_use():
    """Regression for a SEGMENTATION FAULT (2026-07-29, found by flipping the
    default ON): the Picard driver read UpwindOperator._upstream before any walk
    had filled it -- uninitialised int64 garbage, which numba then used as array
    indices. The operator now validates the map instead of crashing, and
    `upstream_map(grad)` is the way to obtain one.
    """
    ent = EntropyOperator(4)
    q2 = np.full(4, float(q2_at_mach(1.2, 0.8)))
    with pytest.raises(ValueError, match="out of range"):
        ent.sigma(q2, np.array([0, 1, 999999, 3]), 0.8)
    with pytest.raises(ValueError, match="expected"):
        ent.sigma(q2, np.array([0, 1]), 0.8)


# --------------------------------------------------------------- criterion C
def test_c_jacobian_is_fd_exact_for_the_frozen_sigma_system():
    """C: the Jacobian of the FROZEN-sigma system is the derivative of its
    residual -- asserted the way the pre-registration wrote it, in TWO parts:

      (i) at the sweet-spot epsilon the relative error is < 1e-6;
      (ii) the error SCALES like 1/epsilon across three decades.

    (ii) is the discriminator, not decoration (phase-one B19): a MISSING TERM
    gives an epsilon-INDEPENDENT relative error (B19 measured 1.532e-01 at every
    epsilon, spread 1.00), while FD roundoff grows like 1/epsilon. Asserting only
    (i) at a single epsilon cannot tell those apart, and asserting (i) at every
    epsilon is simply wrong -- at eps = 1e-8 FD roundoff alone is ~3e-6.
    """
    mc, wc = cut_wake(read_mesh(MESH))
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA)
    ws.set_mach(M_INF)
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             upwind_c=1.5, m_crit=0.95, precond="direct",
                             n_newton_max=6, entropy_correction=True)
    phi_free = np.asarray(r["phi"], dtype=np.float64)[:ws.n_red][ws.free].copy()
    gamma = np.asarray(r["gamma"], dtype=np.float64).copy()

    # freeze sigma exactly as the driver does, then never touch it again
    _, _, st0 = ws.eval_residual(phi_free, gamma, 1.5, 0.95, 3.0, 0.05)
    ws.refresh_sigma(st0, frozen=None)
    _, _, st0 = ws.eval_residual(phi_free, gamma, 1.5, 0.95, 3.0, 0.05)
    assert ws.sigma_frozen is not None and ws.sigma_frozen.min() < 1.0, \
        "the probe state carries no shock -- the test would be vacuous"

    J_ff, _ = ws.assemble_coupled(st0, 1.5, 0.95, 0.05)
    rng = np.random.default_rng(7)
    v = rng.standard_normal(len(phi_free))
    v /= np.linalg.norm(v)
    jv = J_ff @ v
    rels = {}
    for eps in (1e-6, 1e-7, 1e-8):
        Rp, _, _ = ws.eval_residual(phi_free + eps * v, gamma, 1.5, 0.95,
                                    3.0, 0.05)
        Rm, _, _ = ws.eval_residual(phi_free - eps * v, gamma, 1.5, 0.95,
                                    3.0, 0.05)
        fd = (Rp - Rm) / (2.0 * eps)
        rels[eps] = float(np.linalg.norm(jv - fd) / np.linalg.norm(fd))
    spread = rels[1e-8] / rels[1e-6]
    print("\n  ||Jv-FD||/||FD||: "
          + "  ".join(f"eps={e:g}: {v:.3e}" for e, v in rels.items())
          + f"   spread(1e-8/1e-6) = {spread:.1f}")
    assert min(rels.values()) < 1e-6, f"not FD-exact at any eps: {rels}"
    assert spread > 10.0, (
        f"the error does not scale with eps (spread {spread:.2f}) -- that is a "
        f"MISSING TERM, not FD roundoff: {rels}")
