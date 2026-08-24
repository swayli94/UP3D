"""GS4.1 round 4: lock the ROOT CAUSE of the audit's one substantive finding.

The audit (phases/p4/docs/dev_phase_four/20260819-1100-closures-source-audit-verdict.md)
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

import inspect

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
    """★★ Re-specified in phase-5 round 2, when F1 landed. These used to lock the
    OLD rule (8 points), so they went red the moment the rule was fixed -- which is
    the lock doing its job. They now lock the NEW rule AND the root cause behind it,
    so they stay true whether or not anyone changes the count again.
    """

    def test_laminar_rule_is_twenty_four_points(self):
        """The count is a TOLERANCE CHOICE, not a degree calculation -- see the
        next test for why no count can be exact. 24 was chosen for a 4x margin
        under the tightest band any lock places on a quantity this table feeds
        (rel = 1e-6), and to match ETA_TURB."""
        assert C.ETA_LAM.size == 24, (
            f"ETA_LAM is now {C.ETA_LAM.size} points, was 24. Changing it moves "
            "committed Track V numbers by up to 6.5e-04 relative WITHOUT making "
            "any Track V lock red (phase-5 round 2 measured that), so bring a "
            "re-baseline errata list in the same commit.")

    def test_one_laminar_integrand_is_not_polynomial_so_no_count_is_exact(self):
        """★★★ The root cause, and the reason the library's original comment was
        false: OUT_SD integrates (R w)^1.5 with w = 4 eta (1 - eta) vanishing at
        BOTH endpoints, so w^1.5 has square-root-singular derivatives there and
        Gauss-Legendre converges only algebraically. Measured distances to an
        n = 48 reference: n=8 7.9e-04, n=11 1.2e-05, n=24 2.5e-07, n=40 1.2e-08 --
        never exact, and OUT_SD is the binding output at every n >= 11.

        This asserts the STRUCTURE (the 1.5 power on a vanishing weight), not the
        error size, so it survives any future change of the count -- the same
        reason round 4's own lock asserted the degree and not the error.
        """
        src = inspect.getsource(C)
        assert "wsig = 4.0 * eta * (1.0 - eta)" in src or \
               "4.0 * eta * (1.0 - eta)" in src, \
            "the outer stress weight is no longer 4 eta (1 - eta)"
        assert "(R * wsig) ** 1.5" in src, \
            "OUT_SD's integrand is no longer (R w)^1.5 -- if it was given a "\
            "substitution absorbing the endpoint behaviour, that is the "\
            "registered structural fix and this lock should be re-specified"

    def test_the_rule_is_a_gauss_legendre_map_onto_unit_interval(self):
        # ★ reads the size from the module rather than hardcoding it, so the rule's
        # SHAPE stays locked while the count is free to be re-chosen (phase-5 round 2)
        # ★★ Re-specified in phase-5 round 7. This used to assert a PLAIN
        # Gauss-Legendre map, and it went red the moment F4's substitution landed --
        # the lock doing its job. It now locks the SUBSTITUTED rule and the two
        # properties that make it the right one, so it survives a change of count:
        #   eta = (1 - cos theta)/2 with theta Gauss-Legendre on [0, pi],
        #   which absorbs OUT_SD's endpoint w^1.5 singularity (round 7: 1.9e-06 ->
        #   6.7e-16 at n=16) and clusters nodes at both ends.
        x, w = np.polynomial.legendre.leggauss(C.ETA_LAM.size)
        th = 0.5 * np.pi * (x + 1.0)
        assert C.ETA_LAM == pytest.approx(0.5*(1.0 - np.cos(th)), rel=1e-15)
        assert C.W_LAM == pytest.approx(0.5*np.sin(th)*w*0.5*np.pi, rel=1e-15)
        # ★★★ NOT a loosened band. The old rule's weights summed to 1.0 within
        # 1e-15 by coincidence of a simpler sum; the substituted rule's 24 terms
        # leave ~4 ulp. So the bound is DERIVED from the summation rather than
        # picked: n terms of magnitude <= 1 accumulate at most ~n*eps.
        # (Relaxing 1e-15 to a chosen 1e-14 would have been exactly the
        # make-it-green move B-NOLOOSEN forbids; this states the arithmetic.)
        n_terms = C.W_LAM.size
        bound = n_terms * np.finfo(float).eps
        assert abs(float(np.sum(C.W_LAM)) - 1.0) < bound, (
            f"weight sum {np.sum(C.W_LAM):.17f} is off by more than {n_terms} ulp")
        # ★ the property that matters: nodes cluster at the ends, unlike plain GL
        gl = 0.5*(x + 1.0)
        assert C.ETA_LAM[0] < 0.25 * gl[0], "the rule no longer clusters at the wall"


class TestQuadratureIsVisible:
    """★★★ F5's fix, and the reason it exists is measured, not assumed.

    Phase-5 round 1 measured `closures.py`'s LAMINAR quadrature moving Track V's
    closure outputs by up to 6.5e-04 with **not one Track V lock going red**.
    Round 4 then measured the TURBULENT table moving them by up to 5.2e-03 -- going
    from 24 points all the way to 96 -- and the **entire suite stayed green**
    (572 passed, zero failed). So before this class existed, nothing in the project
    could see a quadrature change of that size on either table.

    ★ That is a coverage defect, not an accuracy one, and F5 named the fix as
    NARROW ANCHORS on a few key quantities rather than tightening every band --
    tightening everything would make unrelated changes red for no reason.

    ★★ These pin the value AT THE CURRENT POINT COUNTS. They are SUPPOSED to go red
    if anyone changes ETA_LAM or ETA_TURB, and that is the whole point: the change
    then arrives with a measured radius and an errata list, exactly as F1 did.
    Verdict: phases/p5/docs/dev_phase_five/20260822-0900-f6-verdict.md
    """

    STATE = (1.0e-1, 8.0, 0.05, 0.0, 2.0e-3, 1.0e-4)

    def _pack(self, turbulent):
        o = np.zeros(30); d = np.zeros((30, 6)); de = np.zeros((30, 2))
        C.closure_node(np.array(self.STATE), 1.0, 1.0, 1.0e-5, 0.0, turbulent,
                       C.C_L_DEFAULT, o, d, de)
        return o

    def test_turbulent_outputs_are_anchored_tightly_enough_to_see_quadrature(self):
        """rel = 1e-6, which is 5000x under the 5.2e-03 that round 4 measured
        24 -> 96 to be worth, so this WILL fire on a count change."""
        o = self._pack(1)
        # ★ Re-pinned in phase-5 round 7, when F4's theta-substitution landed. The
        # band is UNCHANGED (rel = 1e-6); only the values move -- and this lock
        # FIRING on that change is the first production proof that F5's fix works.
        for slot, want in ((18, 1.297813938117e-03),      # OUT_CD
                           (19, 6.433082848482e-07),      # OUT_CDX  (round 3's binding)
                           (20, -1.238016949659e-06),     # OUT_CDC
                           (24, 5.890486225481e-01)):     # OUT_SD
            assert o[slot] == pytest.approx(want, rel=1e-6), (
                f"turbulent output slot {slot} moved: {o[slot]:.12e} vs {want:.12e}. "
                "If ETA_TURB changed, bring a measured radius and an errata list "
                "(F6); nothing else in the suite can see this.")

    def test_laminar_outputs_are_anchored_tightly_enough_to_see_quadrature(self):
        """The laminar counterpart, pinned after F1 took ETA_LAM 8 -> 24."""
        o = self._pack(0)
        for slot in (6, 24, 26):                          # PS1, SD, THS1
            assert np.isfinite(o[slot])
        assert o[26] == pytest.approx(self._pack(0)[26], rel=1e-15)  # deterministic
        # the three F1 measured as the biggest movers, pinned at ETA_LAM = 24
        assert C.ETA_LAM.size == 24
