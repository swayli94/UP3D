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


# ---------------------------------------------------------------------------
# the m_cap guard (pre-registered
# docs/dev_phase_two/20260731-2000-entropy-mcap-prereg.md; the defect it fixes
# was measured on M6 medium: m1_max = m_cap exactly, sigma_min = 0.0, 57 floored
# cells, |R| 2.49e-06 un-converged -- the G8.2 signature). All four were verified
# FAILING against the pre-fix kernel before being committed.
# ---------------------------------------------------------------------------

M_CAP = 3.0                    # NewtonWorkspace.eval_residual's default


def _q2(mach, m_inf=0.8):
    return np.array([float(q2_at_mach(m, m_inf)) for m in mach])


def test_mcap_limited_donor_gets_no_correction():
    """A post-shock cell whose donor's speed was LIMITED must not be corrected.

    Pre-fix this read the cap as a physical pre-shock Mach and returned
    sigma_RH(3.0) ~ 0.328 -- a 67 % density cut invented out of a limiter.
    """
    q2 = _q2([M_CAP, 0.80])
    lim = np.array([False, True])          # the donor was limited
    ent = EntropyOperator(2)
    sig = ent.sigma(q2, np.array([0, 0]), 0.8, lim=lim).copy()
    assert sig[1] == 1.0, "a limited donor must produce no entropy correction"
    assert ent.m1_max == 0.0
    assert ent.sigma_min == 1.0


def test_mcap_walk_stops_at_a_limited_cell():
    """The knee walk must stop BEFORE reading a limited cell's Mach.

    Chain: 3.00 (limited) -> 1.30 -> 1.20 -> 0.85. M1 must be 1.30, the last
    physical value, not the cap.
    """
    q2 = _q2([M_CAP, 1.30, 1.20, 0.85])
    lim = np.array([False, True, True, True])
    ent = EntropyOperator(4)
    sig = ent.sigma(q2, np.array([0, 0, 1, 2]), 0.8, lim=lim).copy()
    assert_rel_close(ent.m1_max, 1.30)
    assert_rel_close(sig[3], total_pressure_ratio(1.30))


def test_mcap_chain_product_does_not_collapse():
    """A run of capped cells feeding a subsonic cell must produce no correction.

    Note what this construction does and does not contain: the detector fires only
    where a SUBSONIC cell has a SUPERSONIC donor, so eleven capped cells in a row
    are one shock, worth a single sigma_RH(3.0) = 0.32834 factor pre-fix -- not
    0.32834**11. My first version of this test asserted the latter and failed,
    which is how the mechanism got measured properly; the two ways sigma actually
    reaches 0 are locked in test_mcap_the_mechanism_itself_is_locked.
    """
    n = 12
    mach = np.array([M_CAP] * (n - 1) + [0.85])
    lim = np.array([False] * (n - 1) + [True])
    ent = EntropyOperator(n)
    sig = ent.sigma(_q2(mach), np.arange(-1, n - 1).clip(0), 0.8,
                    lim=lim).copy()
    assert sig.min() == 1.0, f"sigma collapsed to {sig.min()!r}"


def test_mcap_lim_none_is_the_all_true_control():
    """lim=None must reproduce an explicit all-true mask BITWISE (it is the
    control leg of the round's A/B, so it may not drift from it)."""
    q2 = _q2([1.05, 1.37, 1.25, 1.10, 0.95, 0.85])
    up = np.array([0, 0, 1, 2, 3, 4])
    a = EntropyOperator(6).sigma(q2, up, 0.8).copy()
    b = EntropyOperator(6).sigma(q2, up, 0.8,
                                 lim=np.ones(6, dtype=bool)).copy()
    assert np.array_equal(a, b)


def test_mcap_lim_length_is_validated():
    ent = EntropyOperator(4)
    with pytest.raises(ValueError, match="lim mask has 3 entries"):
        ent.sigma(_q2([1.3, 1.2, 0.9, 0.8]), np.array([0, 0, 1, 2]), 0.8,
                  lim=np.ones(3, dtype=bool))


def test_mcap_the_mechanism_itself_is_locked():
    """Lock the two MEASURED routes from a limited cell to a collapsed sigma.

    The other mcap tests fail against the pre-fix library with a TypeError (the
    argument did not exist), which demonstrates less than it looks like. This one
    drives both routes through the SHIPPED API, since lim=None reproduces the
    pre-fix path exactly:

      (A) separate capped shocks along one ACYCLIC chain multiply, one
          sigma_RH(3.0) = 0.32834 factor each: twelve of them give 1.570e-06 and
          the transport still reports converged.
      (B) a capped shock inside a DONOR CYCLE gives exactly 0.0 with
          converged=False -- pointer doubling squares the accumulated product
          every round, so 0.32834 underflows.

    M6 medium pre-fix read sigma_min exactly 0.0, so it was route (B); m1_max
    there was the cap itself, which is what the guard removes at the root. That
    attribution is why "the G8.2 donor-cycle signature" is the right name for it.
    """
    # (A) acyclic: 12 separate capped shocks
    k = 12
    mach, up = [], []
    for i in range(k):
        mach += [M_CAP, 0.85]
        up += [max(2 * i - 1, 0), 2 * i]
    ent_a = EntropyOperator(len(mach))
    a = ent_a.sigma(_q2(np.array(mach)), np.array(up), 0.8).copy()
    assert_rel_close(a.min(), total_pressure_ratio(M_CAP) ** k, rtol=1e-9)
    assert ent_a.converged, "an acyclic chain must not exhaust the transport"

    # (B) the same shock inside a 3-cycle: exact zero, and NOT reported converged
    q2_b, up_b = _q2(np.array([M_CAP, 0.85, 0.90])), np.array([2, 0, 1])
    ent_b = EntropyOperator(3)
    b = ent_b.sigma(q2_b, up_b, 0.8).copy()
    assert b.min() == 0.0
    assert not ent_b.converged, (
        "a donor cycle must be reported, not silently returned as a solution")

    # and with the guard, neither route fires at all
    guarded = EntropyOperator(3).sigma(
        q2_b, up_b, 0.8, lim=np.array([False, True, True])).copy()
    assert guarded.min() == 1.0
