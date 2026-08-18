"""GS4.1 round 4: lock the ROOT CAUSE of the audit's one substantive finding.

The audit (docs/dev_phase_four/20260819-1100-closures-source-audit-verdict.md)
found that `closures.py`'s laminar Gauss rule is too coarse for the
kinetic-energy thicknesses, and that the comment justifying the rule is wrong
about why 8 points suffice.

★ This file locks the ROOT CAUSE -- the integrand's polynomial degree, hence the
number of points needed -- and NOT the error magnitude. That choice matters: the
degree is a property of D13's profile family and stays true whether or not the
quadrature is ever fixed, so this test cannot rot into a lock on a defect. If
someone changes ETA_LAM, the last assertion points them at the verdict.

No solve, no library change. The finding is reported, not fixed: raising the
rule would move every committed Track V number and needs its own round with a
re-baseline errata list.
"""

import numpy as np
import pytest

from pyfp3d.viscous import closures as C


def _profile(e, A=7.9, B=0.6, Psi=-0.25):
    """D13 (42)(43), independently of closures.py -- degree bookkeeping only."""
    f0 = 6*e**2 - 8*e**3 + 3*e**4                       # degree 4
    f1 = e - 3*e**2 + 3*e**3 - e**4                     # degree 4
    f2 = (e - 4*e**2 + 6*e**3 - 4*e**4 + e**5)*(1 - e)**2   # degree 7
    f3 = (e**2 - 3*e**3 + 3*e**4 - e**5)*(1 - e)**2         # degree 7
    U = A*(1 - 0.6*(A - 3)*e**3)*f1 + f0                # degree 3+4 = 7
    W = B*f2 + Psi*f3                                   # degree 7
    return U, W


def _gauss(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5*(x + 1.0), 0.5*w


def _phi_star_1(n):
    """D13 (60) phi*_1 integrand 1 - U(U^2+W^2): CUBIC in a degree-7 profile."""
    e, w = _gauss(n)
    U, W = _profile(e)
    return float(np.sum(w*(1 - U*(U*U + W*W))))


def _phi_11(n):
    """D13 (60) phi_11 integrand 1 - U^2: QUADRATIC, degree 14."""
    e, w = _gauss(n)
    U, _ = _profile(e)
    return float(np.sum(w*(1 - U*U)))


class TestQuadratureDegree:
    """Gauss with n points is exact to degree 2n-1."""

    def test_kinetic_energy_integrand_is_degree_21(self):
        """Exact at n = 11 (2n-1 = 21) and not before -- so it IS degree 21."""
        ref = _phi_star_1(200)
        assert abs(_phi_star_1(10)/ref - 1) > 1e-9, "would not be degree 21"
        assert abs(_phi_star_1(11)/ref - 1) < 1e-13
        assert abs(_phi_star_1(12)/ref - 1) < 1e-13

    def test_quadratic_thickness_is_degree_14(self):
        """Control: this one IS covered by 8 points, which is why the library's
        rule is right for the quadratic thicknesses and wrong only for the
        cubic ones."""
        ref = _phi_11(200)
        assert abs(_phi_11(7)/ref - 1) > 1e-9
        assert abs(_phi_11(8)/ref - 1) < 1e-13

    def test_eight_points_are_insufficient_for_the_cubic_thicknesses(self):
        """The finding, as a number, from the degree fact alone."""
        err = abs(_phi_star_1(8)/_phi_star_1(200) - 1)
        assert 1e-5 < err < 1e-4, f"quadrature error moved: {err:.3e}"


class TestLibraryRuleUnchanged:
    def test_laminar_rule_is_still_eight_points(self):
        """★ If this fails, someone changed the laminar quadrature. That is the
        registered follow-up F1 and it moves every committed Track V number --
        see docs/dev_phase_four/20260819-1100-closures-source-audit-verdict.md
        section 2, and bring a re-baseline errata list.
        """
        assert C.ETA_LAM.size == 8, (
            f"ETA_LAM is now {C.ETA_LAM.size} points, was 8. If this is the "
            "F1 fix, the verdict's section 2.3 requires a re-baseline errata "
            "list in the same commit.")

    def test_the_rule_is_a_gauss_legendre_map_onto_unit_interval(self):
        x, w = np.polynomial.legendre.leggauss(8)
        assert C.ETA_LAM == pytest.approx(0.5*(x + 1.0), rel=1e-15)
        assert C.W_LAM == pytest.approx(0.5*w, rel=1e-15)
        assert float(np.sum(C.W_LAM)) == pytest.approx(1.0, rel=1e-15)
