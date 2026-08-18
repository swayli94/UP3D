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
        stations = np.geomspace(5.3, 200.0, 40)
        kw = dict(rho=RHO, mu=MU, n_substep=4000, x_tr=5.0)
        a = S.march_correlation(stations, (4.793e-3, 1.8080), 5.1,
                                S.flat_plate_ue(U), **kw)
        b = S.march_correlation(stations, (4.793e-3, 1.8080 * 1.15), 5.1,
                                S.flat_plate_ue(U), **kw)
        assert abs(1.8080 * 0.15) / 1.8080 >= 0.10      # seeds really differ
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
        """
        st = self._plate()
        o = np.argsort(st.re_theta)
        ret, H = st.re_theta[o], st.H[o]
        assert np.interp(600.0, ret, H) == pytest.approx(1.512, abs=0.01)
        assert np.interp(1000.0, ret, H) == pytest.approx(1.4516, abs=0.005)
        assert np.interp(3000.0, ret, H) == pytest.approx(1.3809, abs=0.003)

    def test_the_relaxation_decays(self):
        """H must fall monotonically away from that peak -- the transition
        transient relaxes rather than persisting."""
        st = self._plate()
        i = int(np.argmax(st.H))
        tail = st.H[i:]
        assert np.all(np.diff(tail) <= 1e-9), "H does not decay after its peak"
