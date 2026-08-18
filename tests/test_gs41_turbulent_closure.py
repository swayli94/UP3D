"""GS4.1 round 5: locks for the turbulent branch of `closures_2d.py`.

Ungated, no solve. The anchors here are the ones that would be silently lost if
the closure were edited: the sourced constants, the dissipation identity that a
memory-written version got wrong, and the fact that a flat plate settles onto a
physical turbulent shape factor instead of collapsing.

Verdict: docs/dev_phase_four/20260819-1500-turbulent-closure-verdict.md
"""

import numpy as np
import pytest

from pyfp3d.viscous import closures_2d as C2
from pyfp3d.viscous import strip2d as S

RHO, MU, U = 1.0, 1.0e-5, 1.0


class TestSourcedConstants:
    """Each of these is a transcription from a cited line; a silent edit to any
    of them changes the model."""

    def test_constants_match_their_sources(self):
        assert C2.GACON == 6.70          # xbl.f:1559
        assert C2.GBCON == 0.75          # xbl.f:1560
        assert C2.GCCON == 18.0          # xbl.f:1561
        assert C2.HSMIN == 1.500         # xblsys.f:2394
        assert C2.DHSINF == 0.015        # xblsys.f:2394

    def test_ctcon_is_derived_not_typed(self):
        """xbl.f:1569 computes CTCON from GACON and GBCON. Typing the number
        instead would let the three drift apart."""
        assert C2.CTCON == pytest.approx(0.5 / (C2.GACON ** 2 * C2.GBCON),
                                         rel=1e-15)

    def test_kinematic_shape_is_the_identity_incompressible(self):
        for H in (1.3, 1.4, 2.591):
            assert C2.h_kinematic(H, 0.0) == pytest.approx(H, rel=1e-15)


class TestDissipationIdentity:
    """The single error that drove the memory-written version to H = 0.60."""

    def test_di_is_two_cd_over_hstar(self):
        for H, ret in ((1.3, 800.0), (1.4, 2000.0), (2.5, 1.0e4)):
            p = C2.packet_turb(ret * MU / (RHO * U), H, U, rho=RHO, mu=MU)
            assert C2.dissipation_identity(p["cD"], p["H_star"]) == \
                pytest.approx(2.0 * p["cD"] / p["H_star"], rel=1e-14)

    def test_cd_is_not_di(self):
        """★ Guards the actual mistake: c_D and DI differ by H*/2, and H* is
        near 1.75 here, so confusing them is a ~13 % error in dissipation."""
        p = C2.packet_turb(2.0e-2, 1.4, U, rho=RHO, mu=MU)
        di = C2.dissipation_identity(p["cD"], p["H_star"])
        assert abs(di / p["cD"] - 1.0) > 0.1


class TestRange:
    def test_unphysical_H_is_reported_not_clamped(self):
        with pytest.raises(C2.ClosureRangeError):
            C2.packet_turb(1.0e-3, 0.60, U)      # the memory attempt's value
        with pytest.raises(C2.ClosureRangeError):
            C2.packet_turb(1.0e-3, 5.0, U)

    def test_both_h_star_arms_are_close(self):
        """The XFOIL-current and the commented-out older correlation differ by
        under a percent over the measured band -- recorded, not adjudicated."""
        for ret in (600.0, 2000.0, 1.0e4):
            a = C2.h_star_turb(1.4, ret, arm="new")
            b = C2.h_star_turb(1.4, ret, arm="old")
            assert abs(a / b - 1.0) < 0.02


class TestFlatPlate:
    def test_turbulent_plate_settles_on_a_physical_shape_factor(self):
        """The round's physical result: no collapse, and H lands where a ZPG
        turbulent boundary layer belongs."""
        y0 = C2.blasius_state(0.05, ue=U, rho=RHO, mu=MU, H=2.591100)
        st = S.march_correlation(np.geomspace(5.1, 400.0, 30), y0, 0.05,
                                 S.flat_plate_ue(U), rho=RHO, mu=MU,
                                 n_substep=4000, x_tr=5.0)
        m = st.re_theta >= 800.0
        assert np.all(st.H > 1.05) and np.all(st.H < 4.0)
        assert 1.25 < st.H[m].min() and st.H[m].max() < 1.50

    def test_H_decreases_with_reynolds_number(self):
        """★ Recorded because it is what T-EQUIL got wrong: a ZPG turbulent H
        genuinely DRIFTS with Re_theta (about 5 %/decade), so a criterion
        demanding a constant H is asking for something physically false. See
        the verdict, section 3."""
        y0 = C2.blasius_state(0.05, ue=U, rho=RHO, mu=MU, H=2.591100)
        st = S.march_correlation(np.geomspace(5.1, 400.0, 30), y0, 0.05,
                                 S.flat_plate_ue(U), rho=RHO, mu=MU,
                                 n_substep=4000, x_tr=5.0)
        m = st.re_theta >= 800.0
        H, ret = st.H[m], st.re_theta[m]
        assert H[-1] < H[0]
        slope = abs(np.log(H[-1]/H[0]) / np.log10(ret[-1]/ret[0]))
        assert 0.02 < slope < 0.12, f"drift {slope:.3f}/decade"

    def test_laminar_path_is_untouched(self):
        """G-LEGACY as a suite lock: x_tr=None still reproduces round 3."""
        assert C2.zpg_fixed_point() == pytest.approx(2.590433, rel=1e-6)
        y0 = C2.blasius_state(0.01, H=2.591100)
        st = S.march_correlation(np.array([1.0, 100.0]), y0, 0.01,
                                 S.flat_plate_ue(U), n_substep=2000)
        assert st.H[-1] == pytest.approx(2.590433, abs=5e-4)


class TestAttraction:
    """GS4.1 round 6 (E-ATTRACT): the turbulent branch forgets its initial
    condition. This is what "equilibrium" means -- round 5's T-EQUIL wrongly
    tested it as a constant H, which a ZPG turbulent layer does not have."""

    def test_two_different_seeds_collapse_onto_one_curve(self):
        """★ The seed is DERIVED here, not hardcoded.

        An earlier version pinned (4.793e-3, 1.8080), taken from a march at the
        time it was written. Round 8's CtauEQ correction moved that state to
        (4.7898e-3, 1.8335) -- H by 1.41 % -- leaving a stale pre-fix number
        sitting in the suite with no provenance attached. It never failed,
        because the seed is an INPUT and the assertion is about convergence, so
        nothing would ever have flagged it. Computing it removes the whole class
        of problem instead of re-pinning a value that will drift again.
        """
        lead = S.march_correlation(np.array([5.1]),
                                   C2.blasius_state(0.05, ue=U, rho=RHO, mu=MU,
                                                    H=2.591100),
                                   0.05, S.flat_plate_ue(U), rho=RHO, mu=MU,
                                   n_substep=8000, x_tr=5.0)
        th, H0 = float(lead.theta[0]), float(lead.H[0])
        stations = np.geomspace(5.3, 200.0, 40)
        kw = dict(rho=RHO, mu=MU, n_substep=4000, x_tr=5.0)
        a = S.march_correlation(stations, (th, H0), 5.1,
                                S.flat_plate_ue(U), **kw)
        b = S.march_correlation(stations, (th, H0 * 1.15), 5.1,
                                S.flat_plate_ue(U), **kw)
        assert 0.15 >= 0.10                             # seeds really differ
        sep = np.abs(b.H - a.H) / a.H
        assert sep[-1] < 1e-4, f"seeds did not converge: {sep[-1]:.2e}"
        assert sep.argmax() == 0                        # and it only shrinks


class TestPostTransitionRelaxation:
    """GS4.1 round 8. ★ These exist because the round-8 GCC fix moved the
    turbulent readings most strongly at Re_theta < 800 -- and every lock above
    windows at Re_theta >= 800, so the suite did NOT notice. A lock whose window
    excludes the region a change acts on is not covering that change.
    """

    def _plate(self):
        y0 = C2.blasius_state(0.05, ue=U, rho=RHO, mu=MU, H=2.591100)
        return S.march_correlation(np.geomspace(5.1, 400.0, 60), y0, 0.05,
                                   S.flat_plate_ue(U), rho=RHO, mu=MU,
                                   n_substep=8000, x_tr=5.0)

    def test_H_in_the_relaxation_region_is_anchored(self):
        """H at FIXED Re_theta, interpolated -- station-layout independent.

        ★ An earlier version anchored the maximum H below Re_theta 800, which is
        not well defined: it depends on how close the first station lands to
        x_tr, and it read 1.83 here against 1.5213 in the round's own script.
        Two different station sets, which is question 5 again. Interpolating to
        a fixed Re_theta removes the dependence -- verified stable to 0.005
        across 40, 60 and 90 stations.

        ★ Re-pinned in round 9 leg A, once the five missing turbulent terms were
        transcribed: the anchors moved -0.21 %, +0.29 % and +0.85 %, and only the
        Re_theta 3000 one went red. The other two tolerances are set by the
        station-layout spread (0.0055 at Re_theta 600), which is wider than the
        change -- so a lock whose tolerance is dominated by a nuisance
        sensitivity cannot see a real move of comparable size. Recorded rather
        than tightened, because that spread is real.
        """
        st = self._plate()
        o = np.argsort(st.re_theta)
        ret, H = st.re_theta[o], st.H[o]
        assert np.interp(600.0, ret, H) == pytest.approx(1.5088, abs=0.01)
        assert np.interp(1000.0, ret, H) == pytest.approx(1.4558, abs=0.005)
        assert np.interp(3000.0, ret, H) == pytest.approx(1.3926, abs=0.003)

    def test_the_relaxation_decays(self):
        """H must fall monotonically away from that peak -- the transition
        transient relaxes rather than persisting."""
        st = self._plate()
        i = int(np.argmax(st.H))
        tail = st.H[i:]
        assert np.all(np.diff(tail) <= 1e-9), "H does not decay after its peak"


class TestFiveMissingTerms:
    """GS4.1 round 9 leg A. Five XFOIL terms that were absent from this module
    until round 9's dry run read the turbulent block of `xblsys.f` whole instead
    of reading the lines its constants appear on.

    ★★ None of them could have been caught by a guard over the constants already
    written down: 0.995 was a number never typed, DFAC a function that did not
    exist, and max(CFT, CFL) a branch. That is why these locks assert the TERMS,
    not just the constants.
    """

    def test_the_new_constants_match_their_sources(self):
        assert C2.CD_OUT_US == 0.995          # xblsys.f:1014, 1029 -- not 1.0
        assert C2.CD_LAMSTRESS == 0.15        # xblsys.f:1029
        assert C2.US_CLAMP_TRIG == 0.95       # xblsys.f:836
        assert C2.US_CLAMP_VAL == 0.98        # xblsys.f:838
        assert C2.DFAC_C == 2.1               # xblsys.f:970

    def test_cf_at_a_turbulent_station_is_the_max_of_two(self):
        """`xblsys.f:913-921`. Fires near separation, where the laminar
        correlation exceeds the turbulent one."""
        assert C2.cf_turb(1.45, 2000.0) == C2.cf_turb_wall(1.45, 2000.0)
        assert C2.cf_lam_xfoil(2.90, 600.0) > C2.cf_turb_wall(2.90, 600.0)
        assert C2.cf_turb(2.90, 600.0) == C2.cf_lam_xfoil(2.90, 600.0)

    def test_the_wall_dissipation_uses_the_raw_cft(self):
        """★ The asymmetry is XFOIL's: the momentum equation gets
        `max(CFT, CFL)` while the wall dissipation keeps the raw CFT
        (`xblsys.f:959` uses CF2T). Round 9's A-USED probe first read c_D and so
        reported the max as unreachable -- it does not enter c_D at all."""
        p = C2.packet_turb(6.0e-3, 2.90, 1.0, rho=RHO, mu=MU)
        assert p["cf"] > p["cf_wall"]                 # the max fired
        assert p["cf_wall"] == C2.cf_turb_wall(2.90, p["re_theta"])
        cd = C2.cd_turb(p["cf_wall"], p["Us"], p["Ctau_eq"], p["H_star"],
                        p["DFAC"], p["re_theta"])
        assert p["cD"] == pytest.approx(cd, rel=1e-15)

    def test_dfac_fades_the_wall_term_towards_hk_one(self):
        """`xblsys.f:965-985`: `DFAC = (1 + tanh((Hk-1)/(Hmin-1)))/2` with
        `Hmin = 1 + 2.1/ln(Re_theta)`.

        ★ My first version of this lock asserted DFAC = 1/2 at Hk = Hmin, which
        is wrong: that argument is 1, not 0, so DFAC there is tanh(1) shifted =
        0.8808. The half point sits at Hk = 1, where no wake layer exists at all.
        The test was wrong and the transcription was right -- which is the way
        round it should be caught.
        """
        for ret in (600.0, 3000.0, 1.0e5):
            hmin = 1.0 + C2.DFAC_C / np.log(ret)
            assert C2.dfac_low_hk(1.0, ret) == pytest.approx(0.5, rel=1e-12)
            assert C2.dfac_low_hk(hmin, ret) == pytest.approx(
                0.5 + 0.5 * np.tanh(1.0), rel=1e-12)
            assert 0.99 < C2.dfac_low_hk(3.0, ret) <= 1.0
            assert 0.85 < C2.dfac_low_hk(1.45, ret) < 1.0
        # strongest at low Re_theta, where Hmin is largest
        assert C2.dfac_low_hk(1.45, 600.0) < C2.dfac_low_hk(1.45, 30000.0)

    def test_cd_carries_all_three_contributions(self):
        """Wall + outer turbulent + outer laminar stress. Removing any one
        changes c_D, which is what rounds 5-8 were silently doing to two of
        them."""
        p = C2.packet_turb(2.0e-2, 1.40, 1.0, rho=RHO, mu=MU)
        us, ret = p["Us"], p["re_theta"]
        wall = 0.5 * p["cf_wall"] * us * p["DFAC"]
        outer = p["Ctau_eq"] * (C2.CD_OUT_US - us)
        lam = C2.CD_LAMSTRESS * (C2.CD_OUT_US - us) ** 2 / ret
        assert p["cD"] == pytest.approx(wall + outer + lam, rel=1e-14)
        assert lam / p["cD"] > 0.005          # not negligible at Re_theta 2000
        pre_round9 = 0.5 * p["cf"] * us + p["Ctau_eq"] * (1.0 - us)
        assert abs(p["cD"] / pre_round9 - 1.0) > 0.005

    def test_the_us_clamp_is_wired_but_unreachable(self):
        """★★ Round 9 addendum #2. The clamp is transcribed faithfully and IS
        wired into `slip_velocity`, but `H*` is bounded by 2 as Hk -> 1 while
        `H_TURB_LO` holds H at 1.05, so raw Us never exceeds 0.9221 and the 0.95
        trigger cannot fire through `packet_turb`.

        A-USED as first registered would have called that a FAIL, which is the
        one-sided-criterion mistake: "faithful but unreachable" and "never wired
        in" need opposite responses. This lock records BOTH halves, so lowering
        H_TURB_LO -- the change that would make it reachable -- lands here.
        """
        assert C2.slip_velocity(2.4, 1.05) == C2.US_CLAMP_VAL      # wired
        hi = max(0.5 * C2.h_star_turb(H, ret)
                 * (1.0 - (H - 1.0) / (C2.GBCON * H))
                 for H in np.linspace(C2.H_TURB_LO, C2.H_TURB_HI, 60)
                 for ret in np.geomspace(200.0, 1.0e8, 40))
        assert hi == pytest.approx(0.922068, abs=1e-5)             # unreachable
        assert hi < C2.US_CLAMP_TRIG


class TestLagEquation:
    """GS4.1 round 9 leg B. The shear-stress lag, transcribed from
    `xblsys.f:1769-1771`'s two-point residual in its continuous limit.

    ★ Rounds 5-8 deliberately did NOT carry it, because a zero-pressure-gradient
    plate cannot test a lag: equilibrium IS Ctau = CtauEQ there, so the source
    term vanishes identically. These locks assert the machinery and the ONE
    thing a plate CAN show -- that the two arms converge downstream.
    """

    def test_the_lag_constants_match_their_sources(self):
        assert C2.SCCON == 5.6                # xbl.f:1558
        assert C2.DUXCON == 1.0               # xbl.f:1567
        assert C2.DLCON == 0.9                # xbl.f:1562 -- wake only
        assert C2.HDMAX == 12.0               # xblsys.f:1112
        assert (C2.DE_A, C2.DE_B) == (3.15, 1.72)     # xblsys.f:1103
        assert (C2.CTRCON, C2.CTRCEX) == (1.8, 3.3)   # xbl.f:1564-1565

    def test_delta_is_the_thickness_not_delta_star(self):
        """★ The lag's relaxation LENGTH is Delta, about eight times delta*.
        Confusing them puts the relaxation rate out by that factor, which is the
        registered way leg B's first prediction can be wrong."""
        theta, H = 2.0e-2, 1.45
        de = C2.bl_thickness(theta, H)
        assert de == pytest.approx((C2.DE_A + C2.DE_B / (H - 1.0)) * theta
                                   + theta * H, rel=1e-14)
        assert de / (theta * H) > 5.0                       # >> delta*
        # the cap is the source's own guard as Hk -> 1
        assert C2.bl_thickness(theta, 1.0001) == C2.HDMAX * theta
        assert C2.bl_thickness(theta, 1.0) == C2.HDMAX * theta

    def test_the_transition_seed_is_three_decades_below_equilibrium(self):
        """`xblsys.f:1393, 1403`: S_tr = CTRCON exp(-CTRCEX/(Hk-1)) sqrt(CtauEQ).
        The prefactor is what makes the lag arm distinguishable at all."""
        ct_eq = 1.5e-3
        for H, lo, hi in ((1.5, 1e-4, 1e-2), (2.0, 1e-2, 0.2)):
            r = C2.s_tau_at_transition(H, ct_eq) / np.sqrt(ct_eq)
            assert lo < r < hi, f"H={H} prefactor {r:.3e}"
        assert C2.s_tau_at_transition(1.5, ct_eq) < 0.01 * np.sqrt(ct_eq)

    def test_the_lag_source_vanishes_at_equilibrium(self):
        """The defining property: with Ctau = CtauEQ and no pressure gradient
        beyond UQ, the relaxation term is exactly zero."""
        assert C2.lag_rate(0.04, 0.04, 0.55, 0.1, 0.0, 0.0) == 0.0
        assert C2.lag_rate(0.02, 0.04, 0.55, 0.1, 0.0, 0.0) > 0.0   # climbs
        assert C2.lag_rate(0.06, 0.04, 0.55, 0.1, 0.0, 0.0) < 0.0   # decays

    def test_rhs_turb_reduces_to_the_equilibrium_arm(self):
        """★ G-LEGACY at the equation level: handing the lag arm exactly
        sqrt(CtauEQ) must reproduce the equilibrium arm's first two components,
        so the two arms differ only through the transported Ctau."""
        theta, H, ue = 2.0e-2, 1.45, 1.0
        p = C2.packet_turb(theta, H, ue, rho=RHO, mu=MU)
        eq = C2.rhs_turb(theta, H, ue, 0.0, rho=RHO, mu=MU)
        lg = C2.rhs_turb(theta, H, ue, 0.0, rho=RHO, mu=MU,
                         s_tau=np.sqrt(p["Ctau_eq"]))
        assert len(eq) == 2 and len(lg) == 3
        assert lg[0] == pytest.approx(eq[0], rel=1e-15)
        assert lg[1] == pytest.approx(eq[1], rel=1e-15)
        # ★ and the third component is NOT zero here. My first version asserted
        # it was, on the intuition that "Ctau = CtauEQ means equilibrium" -- but
        # the source says otherwise: with the relaxation term vanishing, the
        # remaining term is DUXCON (UQ - u_e'/u_e), and UQ is the pressure
        # gradient an equilibrium layer WANTS, which a flat plate does not
        # supply. So a ZPG plate is not an equilibrium state for the lag, and
        # dS/dxi = S UQ exactly. Asserting the transcription, not the intuition.
        uq = C2.uq_equilibrium(p["cf_wall"], H, p["re_theta"], theta)
        assert lg[2] == pytest.approx(np.sqrt(p["Ctau_eq"]) * uq, rel=1e-12)
        assert uq < 0.0                      # measured; pulls Ctau below CtauEQ

    def test_the_two_arms_converge_on_a_flat_plate(self):
        """★★ The plate's one honest lag reading, and it is a NEGATIVE result:
        the arms agree to ~1 % far downstream, which is exactly why rounds 5-8
        could not have tested the lag. Near transition they differ by 4x."""
        y0 = C2.blasius_state(0.05, ue=U, rho=RHO, mu=MU, H=2.591100)
        st = np.geomspace(5.1, 400.0, 30)
        kw = dict(rho=RHO, mu=MU, n_substep=4000, x_tr=5.0)
        a = S.march_correlation(st, y0, 0.05, S.flat_plate_ue(U), **kw)
        b = S.march_correlation(st, y0, 0.05, S.flat_plate_ue(U), lag=True, **kw)
        assert b.ctau[0] / a.ctau[0] < 0.5            # still climbing
        assert abs(b.ctau[-1] / a.ctau[-1] - 1.0) < 0.05   # converged
        assert b.H[0] > a.H[0]                        # less stress, fuller H

    def test_lag_off_is_bit_identical(self):
        """G-LEGACY: the default path did not move when the lag arm landed."""
        y0 = C2.blasius_state(0.05, ue=U, rho=RHO, mu=MU, H=2.591100)
        st = S.march_correlation(np.geomspace(5.1, 400.0, 30), y0, 0.05,
                                 S.flat_plate_ue(U), rho=RHO, mu=MU,
                                 n_substep=4000, x_tr=5.0)
        assert st.H[-1] == pytest.approx(1.2932778384340817, rel=1e-14)

    def test_ald_is_explicit_so_the_wall_assumption_is_visible(self):
        """★ `ALD` is 1 on a wall layer and DLCON in the wake
        (`xblsys.f:1701-1705`, verified from the source). It is an argument, not
        a hardcoded 1, because round 5 annotated GCCON "wake only, unused here"
        on exactly this pattern and was wrong -- so an out-of-scope constant gets
        a visible use site rather than a comment."""
        base = dict(s_tau=0.04, s_tau_eq=0.04, us=0.55, delta=0.1, uq=0.0,
                    due_over_ue=0.0)
        assert C2.lag_rate(**base) == 0.0                     # wall, ALD = 1
        assert C2.lag_rate(ald=C2.DLCON, **base) != 0.0       # wake, ALD = 0.9


class TestSeedAgainstXfoil:
    """GS4.1 round 17 (G20). ★★★ The first EXTERNALLY validated anchor for the
    seed chain, and the only one in this file whose expected values come from
    outside the project.

    The design is what makes it clean. XFOIL's forced transition point is placed
    exactly ON a station (`XTR` = that station's x/c), and `xblsys.f:435` tests
    `XIFORC .LE. X2` while `:451` assigns `XT = XIFORC` -- so the transition
    interval's TURBULENT part has zero length. That station therefore keeps its
    LAMINAR theta and H, and its stored CTAU is exactly XFOIL's own
    `ST = CTR * CQ` at that state. No march, no discretisation, no interpolation
    stands between the two numbers.

    Measured 2026-08-21 with the locally rebuilt XFOIL 6.99 on NACA 0012,
    Re 3e6, M 0, alpha 2 deg, 280 panels, XTR 0.05370 both surfaces; verified by
    the polar's Top_Xtr = Bot_Xtr = 0.05370 and by H staying on the laminar
    plateau at that station and collapsing at the next.

    ★ What this does NOT establish: in general use XT falls INSIDE an interval and
    the state fed to the seed is our own laminar march's, not XFOIL's. This
    verifies the seed FORMULA, not the seed's input in general use.
    """

    #: (theta, H, ue) at the on-station transition point, and XFOIL's own CT there
    XFOIL_ST = (
        (7.900000e-05, 2.5660, 1.32970, 1.661500e-02),   # upper
        (6.800000e-05, 2.3980, 1.01430, 1.203100e-02),   # lower
    )

    def test_seed_chain_reproduces_xfoils_own_ST(self):
        for theta, H, ue, ct_xfoil in self.XFOIL_ST:
            p = C2.packet_turb(theta, H, ue, rho=1.0, mu=1.0 / 3.0e6)
            ours = C2.s_tau_at_transition(H, p["Ctau_eq"])
            assert ours == pytest.approx(ct_xfoil, rel=0.02), \
                f"H={H}: ours {ours:.6e} vs XFOIL {ct_xfoil:.6e}"
            # measured 0.12 % and 0.15 % -- pinned an order tighter than the band
            assert abs(ours / ct_xfoil - 1.0) < 0.005

    def test_the_implied_CTR_matches_the_formula(self):
        """★ Splits the chain: CTR alone, against what XFOIL's ST implies."""
        for theta, H, ue, ct_xfoil in self.XFOIL_ST:
            p = C2.packet_turb(theta, H, ue, rho=1.0, mu=1.0 / 3.0e6)
            ctr_formula = C2.CTRCON * np.exp(-C2.CTRCEX / (H - 1.0))
            ctr_implied = ct_xfoil / np.sqrt(p["Ctau_eq"])
            assert ctr_formula == pytest.approx(ctr_implied, rel=0.005)
