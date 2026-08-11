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

from pyfp3d.kernels.entropy import (EntropyOperator, shock_factor_sweep,
                                    total_pressure_ratio, transport_sigma)
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.physics.isentropic import (GAMMA, critical_speed_squared,
                                       density_isentropic, q2_at_mach)
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting

from ._tol import assert_rel_close

MESH = "cases/meshes/naca0012_2.5d/coarse.msh"
M_INF, ALPHA = 0.7875, 1.25
#: NewtonWorkspace.eval_residual's default. Module-level rather than beside the mcap
#: tests because a @parametrize case needs it at import time.
M_CAP = 3.0


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


def _sigma_by_walking(s, up):
    """The answer this kernel exists to compute, by definition: walk each chain
    hop by hop and multiply, stopping when a node repeats. Independent of pointer
    doubling, so it is an oracle rather than a second opinion.

    Returns (sigma, feeds_shocked_cycle). The flag matters because sigma is only
    DEFINED where the walk terminates: a chain reaching a cycle that carries
    s < 1 has no finite product (going round multiplies again every lap), which is
    exactly the case the kernel must refuse rather than answer.
    """
    n = len(up)
    out = np.ones(n)
    bad = np.zeros(n, dtype=bool)
    for e in range(n):
        c, seen, order = e, set(), []
        while c not in seen:
            seen.add(c)
            order.append(c)
            out[e] *= s[c]
            if up[c] == c:
                break
            c = up[c]
        else:                                   # loop exited via the condition
            cyc = order[order.index(c):]
            bad[e] = any(s[k] != 1.0 for k in cyc)
    return out, bad


def test_harmless_donor_cycle_converges_with_the_correct_sigma():
    """A donor cycle carrying NO shock must converge, with sigma right.

    ★ RE-SPECIFIED 2026-08-05, and the old assertion was backwards on this data.
    This case used to assert converged is False, under the name
    "test_donor_cycle_is_detected_not_silently_wrong" and the reasoning "a 2-cycle
    settles while the product keeps squaring every round, reporting converged with
    a corrupted sigma". Measured, on this test's own numbers: the two cycle
    elements are the SUPERSONIC pair, the detector only fires on a
    supersonic->subsonic donor transition, so s == 1 on both and the product
    squares 1 -> 1 forever. Nothing is corrupted; sigma here is exactly right,
    and it is asserted below against the walking oracle.

    So the old lock did not test what its name said. It tested the old
    termination condition ("every ancestor is a genuine root"), which a cycle can
    never satisfy -- and that is a real defect, not a safeguard: measured on two
    solves, exactly 2 of 68624 and 2 of 90099 elements sat on such a cycle with
    s == 1 on both, and a solve at |R| = 8.85e-15 with zero clamped and zero shock
    cells was reported FAILED because of them
    (docs/dev_phase_two/20260805-0200-sigma-transport-root-cause.md).

    The mechanism the old docstring described is real, but it lives in the
    SHOCKED-cycle case, which that data did not contain and which the next test
    and test_mcap_the_mechanism_itself_is_locked route (B) now lock explicitly.
    """
    mach = np.array([1.30, 1.30, 0.80, 0.75, 0.70])
    q2 = np.array([float(q2_at_mach(m, 0.8)) for m in mach])
    up = np.array([1, 0, 1, 2, 3])                       # u(0)=1, u(1)=0
    s = np.empty(5)
    shock_factor_sweep(q2, up, 0.8, GAMMA, EntropyOperator.MAX_WALK_DEFAULT,
                       EntropyOperator.KNEE_FRAC_DEFAULT, np.ones(5, bool),
                       s, np.empty(5))
    assert s[0] == 1.0 and s[1] == 1.0, \
        "premise: the cycle must carry no shock, else this is the other test"
    ent = EntropyOperator(5)
    out = ent.sigma(q2, up, 0.8).copy()
    assert ent.converged is True
    assert ent.n_rounds < ent.n_round
    ref, bad = _sigma_by_walking(s, up)
    assert not bad.any(), "premise: nothing may feed a shocked cycle here"
    assert np.array_equal(out, ref)


@pytest.mark.parametrize("mach,up,cyc", [
    ([1.30, 0.80], [1, 0], "2-cycle"),
    ([M_CAP, 0.85, 0.90], [2, 0, 1], "3-cycle, capped shock"),
    ([1.45, 0.82, 0.88, 0.80, 0.75, 0.70], [2, 0, 1, 2, 3, 4],
     "3-cycle + downstream chain"),
])
def test_shocked_donor_cycle_is_still_refused(mach, up, cyc):
    """A donor cycle that CARRIES a shock must still report converged=False.

    This is the case the previous test's name claimed and its data did not have:
    with s < 1 somewhere on the cycle, pointer doubling squares the accumulated
    product every round and it runs to exactly 0 -- a zero density, i.e. garbage
    -- so the transport must exhaust its rounds and the caller must refuse.

    ★ It is the criterion that decides this, and one candidate criterion FAILED
    here and was reverted the same round: "terminate when the product stops
    moving" accepts these, because a product squaring toward zero REACHES zero
    and then stops moving. Zero is the most stable fixed point there is. The
    shipped criterion asks whether the ancestor's product is exactly 1 --
    "contributes nothing" -- which 0 is not, so the refusal survives.
    """
    q2 = np.array([float(q2_at_mach(m, 0.8)) for m in mach])
    ent = EntropyOperator(len(mach))
    out = ent.sigma(q2, np.asarray(up), 0.8).copy()
    assert ent.converged is False, f"{cyc}: a shocked cycle must be refused"
    assert ent.n_rounds == ent.n_round
    assert out.min() == 0.0, "the collapse this refusal exists to catch"


@pytest.mark.parametrize("seed,shocked_cycle", [(0, False), (0, True),
                                               (1, False), (1, True),
                                               (2, True), (3, True), (4, True)])
def test_transport_equals_the_walking_oracle_on_random_graphs(seed,
                                                              shocked_cycle):
    """★ COMPLETENESS, absolutely: on random donor graphs the transport's sigma
    must equal the hop-by-hop walk.

    This is the property the termination criterion has to preserve, and the reason
    it is tested here rather than on a live solve: completeness is a KERNEL
    property, and a live solve can only show that some criterion fired, not that
    the accumulated product was finished. Pointer doubling and the walking oracle
    share no code, so agreement is evidence rather than a second opinion.

    ★ It earns its place: it FAILED, on 4 of 5 graphs, against a candidate
    criterion ("stop when the ancestor's own product is 1") that passed every
    hand-built case including the collapse locks. That criterion drops a shock
    sitting further upstream than the segment it inspects. Neither the hand-built
    tests nor a live solve caught it.

    Each graph carries all four cases the criterion must separate: genuine roots,
    long chains, a HARMLESS cycle, and -- when `shocked_cycle` -- a cycle carrying
    a shock, which must make the transport exhaust its rounds so the caller
    refuses. Both are asserted as premises rather than hoped for.

    Two properties of the comparison are deliberate:

      the tolerance is rtol 1e-12, not bit-equality. Pointer doubling multiplies
      the same factors in a different ORDER than a sequential walk and floating
      point multiplication is not associative, so the last bits must differ. The
      first draft asserted array_equal and failed on graphs where the two answers
      agreed to every printed digit -- a test defect, not a kernel one.

      s == 1 is imposed at every ROOT. That is not a convenience: the kernel
      documents it as a precondition and shock_factor_sweep guarantees it by
      skipping u == e, so a root carrying s < 1 has no defined answer (doubling
      squares it toward zero, which is the documented behaviour). The first draft
      violated it and read the consequence as a kernel bug.
    """
    rng = np.random.default_rng(seed)
    n = 3000
    # a forest of chains: each element points to a random LOWER index (acyclic),
    # plus a few long-range edges, then roots and cycles are stitched in
    up = np.maximum(np.arange(n) - rng.integers(1, 40, n), 0).astype(np.int64)
    up[rng.choice(n, 12, replace=False)] = rng.choice(n, 12, replace=False)
    for r in rng.choice(n, 20, replace=False):
        up[r] = r
    h = rng.choice(n, 2, replace=False)                  # the harmless 2-cycle
    up[h[0]], up[h[1]] = h[1], h[0]
    c = rng.choice(n, 3, replace=False)                  # and a 3-cycle
    up[c[0]], up[c[1]], up[c[2]] = c[1], c[2], c[0]

    s = np.ones(n)
    s[rng.choice(n, 200, replace=False)] = 1.0 - 0.3 * rng.random(200)
    s[up == np.arange(n)] = 1.0                          # the ROOT precondition
    s[h] = 1.0                                           # keep this cycle harmless
    if shocked_cycle:
        s[c[0]] = 0.83
    else:
        s[c] = 1.0
    assert np.all(s[up == np.arange(n)] == 1.0), "root precondition violated"
    assert s[h[0]] == 1.0 and s[h[1]] == 1.0, "the harmless cycle must be harmless"

    out = np.empty(n)
    rounds = transport_sigma(s, up, 24, out, np.empty(n, dtype=np.int64),
                            np.empty(n, dtype=np.int64), np.empty(n))
    ref, bad = _sigma_by_walking(s, up)
    #: `bad` comes from the ORACLE, not from the kernel's output. Excluding
    #: "whatever came out as 0" would let a wrong answer hide behind the
    #: exclusion; excluding "what has no defined answer" cannot.
    assert np.allclose(out[~bad], ref[~bad], rtol=1e-12, atol=0.0), (
        f"seed {seed}: transport disagrees with the walk on "
        f"{np.count_nonzero(~np.isclose(out[~bad], ref[~bad], rtol=1e-12))} of "
        f"{np.count_nonzero(~bad)} elements with a defined sigma")
    if shocked_cycle:
        assert bad.any(), "premise: a shocked cycle must be present"
        assert rounds == 24, ("a shocked cycle is present, so the transport must "
                              "exhaust its rounds and let the caller refuse")
    else:
        assert not bad.any(), "premise: no element may lack a defined sigma"
        assert rounds < 24, ("every cycle here is harmless, so the transport must "
                             "settle -- this is the false failure being fixed")


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


# --- the TEMPORARY sigma_scale instrument (phase 3, 2026-08-12) --------------------------------
# ★ These two tests exist only while the knob does. It was pre-registered as an INSTRUMENT with
# its removal criterion fixed in advance (docs/dev_phase_three/20260812-0300-sigma-strength-prereg
# .md), per the B20 precedent where a knob built solely to make an A/B measurable was deleted on
# adoption. When sigma_scale goes, these go with it.
#
# What they lock is the property the whole sweep rests on: at the default the instrument must not
# touch anything. That was verified once by a 290 s bench run against the committed M1 anchors
# (12 rows x 14 fields, 0 differences); a unit test makes it reproducible in milliseconds.

def _two_element_shock():
    """One supersonic donor feeding one subsonic element -- the minimal post-shock pair."""
    q2 = np.array([1.5, 0.6])
    up = np.array([0, 0])
    return q2, up


def test_sigma_scale_one_is_bit_identical_to_the_unscaled_path():
    """★ BY CONSTRUCTION, not by the algebra happening to be exact: `1 - 1.0*(1 - s)` is NOT
    bit-identical to `s` in floating point (0.9812 -> 0.9811999999999999), so a branchless blend
    would have silently moved every committed number at the DEFAULT setting. The kernel therefore
    short-circuits on theta == 1.0, and this test is what keeps that short-circuit in place."""
    q2, up = _two_element_shock()
    base = EntropyOperator(2).sigma(q2, up, 0.8).copy()
    scaled = EntropyOperator(2, sigma_scale=1.0).sigma(q2, up, 0.8).copy()
    assert np.array_equal(base, scaled), "theta = 1 must not perturb the array at all"
    assert base.min() < 1.0, "premise: this pair must actually carry a correction"


def test_sigma_scale_zero_gives_exactly_one_and_the_dial_is_monotone():
    """theta = 0 must give sigma identically 1.0 -- that is what makes agreement with
    entropy_correction=False checkable -- and intermediate values must interpolate monotonically,
    since the dose-response reading is meaningless if the dial itself is not ordered."""
    q2, up = _two_element_shock()
    off = EntropyOperator(2, sigma_scale=0.0).sigma(q2, up, 0.8).copy()
    assert np.all(off == 1.0)

    mins = [EntropyOperator(2, sigma_scale=t).sigma(q2, up, 0.8).min()
            for t in (0.0, 0.25, 0.5, 0.75, 1.0)]
    assert all(b <= a for a, b in zip(mins, mins[1:])), f"dial not monotone: {mins}"
    #: the blend is exact at the endpoints and linear between them
    full = mins[-1]
    assert mins[2] == pytest.approx(1.0 - 0.5 * (1.0 - full), rel=1e-15)


@pytest.mark.parametrize("bad", [-0.1, 1.1])
def test_sigma_scale_out_of_range_raises(bad):
    with pytest.raises(ValueError, match="sigma_scale"):
        EntropyOperator(2, sigma_scale=bad)
