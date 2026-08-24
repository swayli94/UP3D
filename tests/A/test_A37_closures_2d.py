"""GS4.1 round 3: locks for the correlation closure (`closures_2d.py`).

Ungated, no solve. Same reasoning as round 1's locks: a library module whose
only assertions live in a bench script runs on no cadence at all.

★ The anchors below come from an INDEPENDENT Falkner-Skan ODE integration
(`bench/studies/gs41_a2_correlation/results/g_oracle.csv`), not from the
correlations themselves. That is what makes them a transcription lock rather
than a tautology.
"""

import numpy as np
import pytest

from pyfp3d.viscous import closures_2d as C2

#: Blasius, from the ODE oracle (round 1's falkner_skan, committed in g_oracle.csv)
H_BLASIUS = 2.591100
HS_ODE_BLASIUS = 1.572580        # theta*/theta
CF_ODE_BLASIUS = 0.220520        # Re_theta * c_f/2
CD_ODE_BLASIUS = 0.220520        # Re_theta * 2 c_D / H*


class TestTranscription:
    """The correlations against the ODE oracle -- these catch a typo."""

    def test_h_star_at_blasius(self):
        assert C2.h_star(H_BLASIUS) == pytest.approx(HS_ODE_BLASIUS, rel=1e-3)

    def test_cf_at_blasius(self):
        assert C2.re_theta_cf_half(H_BLASIUS) == pytest.approx(
            CF_ODE_BLASIUS, rel=1e-3)

    def test_cd_at_blasius(self):
        assert C2.re_theta_2cd_over_hstar(H_BLASIUS) == pytest.approx(
            CD_ODE_BLASIUS, rel=1e-3)

    @pytest.mark.parametrize("H,hs_ode", [(3.22001, 1.529160),
                                          (2.81817, 1.551960),
                                          (2.49566, 1.583500),
                                          (2.21623, 1.625750)])
    def test_h_star_across_the_wedge_family(self, H, hs_ode):
        assert C2.h_star(H) == pytest.approx(hs_ode, rel=2e-3)

    def test_dh_star_dh_matches_finite_difference(self):
        for H in (2.2, 2.591100, 3.0, 3.5):
            fd = (C2.h_star(H + 1e-7) - C2.h_star(H - 1e-7)) / 2e-7
            assert C2.dh_star_dh(H) == pytest.approx(fd, rel=1e-5)


class TestStructure:
    """The two facts the module is built around."""

    def test_blasius_is_a_fixed_point_by_construction(self):
        """The whole content of route (a2), stated as the quantity that is
        actually well posed.

        An earlier version of this test bounded `dH/dxi` directly, which is a
        RATE (units 1/length) compared against a dimensionless number -- it
        could be made to pass or fail by choosing units. The fixed point itself
        is dimensionless and Reynolds-independent, so that is what gets locked.

        Profile family (round 1): H = 2.708292, +4.52 % off Blasius.
        This family: within 0.03 %.
        """
        H_fix = C2.zpg_fixed_point()
        assert abs(H_fix / H_BLASIUS - 1.0) < 3e-4, H_fix

    def test_the_fixed_point_is_reynolds_independent(self):
        """Both source terms carry one factor 1/Re_theta, so the similar state
        cannot depend on Reynolds number. If it did, the collapse would be
        wrong."""
        H_fix = C2.zpg_fixed_point()
        for theta, ue in ((1e-4, 1.0), (1e-2, 1.0), (1e-3, 40.0)):
            _, dH = C2.rhs(theta, H_fix, ue, 0.0)
            assert abs(dH) < 1e-6 * abs(
                C2.rhs(theta, H_fix * 1.05, ue, 0.0)[1]) + 1e-9

    def test_direct_form_is_singular_at_separation(self):
        """dH*/dH vanishes at H = 4. Known property, GS4.2's motivation."""
        assert C2.dh_star_dh(4.0) == pytest.approx(0.0, abs=1e-15)
        assert abs(C2.dh_star_dh(3.99)) < 1e-3
        assert abs(C2.dh_star_dh(H_BLASIUS)) > 0.1

    def test_range_errors_report_rather_than_clamp(self):
        with pytest.raises(C2.ClosureRangeError):
            C2.packet(1e-3, 1.0, 1.0)
        with pytest.raises(C2.ClosureRangeError):
            C2.packet(1e-3, float("nan"), 1.0)


class TestAuthority:
    """The split fixed before the module was written."""

    def test_the_two_closure_families_do_not_import_each_other(self):
        import ast

        def imports(path):
            out = set()
            for node in ast.walk(ast.parse(open(path).read())):
                if isinstance(node, ast.Import):
                    out.update(a.name for a in node.names)
                elif isinstance(node, ast.ImportFrom):
                    out.add(node.module or "")
            return out

        from pyfp3d.viscous import closures as C1
        c2 = imports(C2.__file__)
        c1 = imports(C1.__file__)
        assert not [m for m in c2 if "closures" in m and "closures_2d" not in m]
        assert not [m for m in c1 if "closures_2d" in m]


class TestMarch:
    def test_correlation_march_holds_blasius_over_three_decades(self):
        from pyfp3d.viscous import strip2d as S
        y0 = C2.blasius_state(0.01, H=H_BLASIUS)
        st = S.march_correlation(np.geomspace(0.1, 100.0, 13), y0, 0.01,
                                 S.flat_plate_ue(1.0), n_substep=2000)
        assert np.max(np.abs(st.H - H_BLASIUS)) / H_BLASIUS < 5e-4
        cf_r = st.cf * np.sqrt(st.re_x)
        assert np.max(np.abs(cf_r - 0.664115)) / 0.664115 < 5e-4

    def test_profile_march_is_untouched_by_the_new_route(self):
        """G-LEGACY as a suite lock: adding (a2) moved no round-1 reading."""
        from pyfp3d.viscous import strip2d as S
        A, H, cf = S.similarity_fixed_point(m=0.0, rho=1.0, mu=1.0e-5)
        # ★ Re-pinned in phase-5 round 2 when closures.py's laminar quadrature went
        # 8 -> 24 points. The bands are UNCHANGED (rel=1e-6); only the values move:
        # A by -5.0e-05, H by -1.3e-04, cf by +9.3e-05.
        # ★★ An earlier version of this comment said "cf did not move at this
        # precision" -- INFERRED from cf not appearing in the failure list rather
        # than measured, and wrong. It moved; its own assertion simply had not been
        # reached, because the A assertion above it fails first. Same trap as b7/b9:
        # the first failing assert hides the rest.
        assert A == pytest.approx(8.028406562, rel=1e-6)
        assert H == pytest.approx(2.707931, rel=1e-6)
        assert cf == pytest.approx(0.710301160, rel=1e-6)

    def test_separation_stops_the_leg_rather_than_limping_on(self):
        from pyfp3d.viscous import strip2d as S
        # a strong adverse gradient drives H up into the guard
        with pytest.raises(C2.ClosureRangeError):
            S.march_correlation(np.geomspace(0.2, 50.0, 10),
                                C2.blasius_state(0.1, H=H_BLASIUS), 0.1,
                                S.falkner_skan_ue(-0.09), n_substep=4000)
