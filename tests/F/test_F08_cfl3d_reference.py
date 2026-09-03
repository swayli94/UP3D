"""F08 -- locks on the CFL3D reference-data machinery.

Two things are locked here, both of which decide PUBLISHED numbers and neither
of which had any cadence before 2026-09-03:

1. ``NITFO < NCYC``.  With ``NITFO >= NCYC`` CFL3D runs first order in space on
   the finest grid for the entire solve while the deck still shows the
   requested ``RKAP0``.  Measured cost of not having this guard: the whole M6
   Euler ladder and the whole M6 RANS ladder ran first order and were discarded
   (cd was wrong by -57.6 %).  See
   ``docs/dev_phase_six/20260903-0200-cfl3d-first-order-defect.md``.

2. ``implied_order`` -- the estimator that decides whether a published quantity
   gets an error bar.  It is checked against SYNTHETIC data whose answer is
   known, because a convergence estimator validated only on the data it is used
   on cannot be distinguished from a curve fit.

★ Both are exercised BEHAVIOURALLY: the guards are fed inputs that must make
them fire.  A guard nobody has seen fail is not known to work.
"""

import os
import sys
import tempfile

import numpy as np
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
CFL3D_DIR = os.path.join(REPO_ROOT, 'cases', 'reference_data', 'cfl3d')
if CFL3D_DIR not in sys.path:
    sys.path.insert(0, CFL3D_DIR)


class TestNitfoGuard:
    """NITFO >= NCYC must be impossible to write into a deck silently."""

    @pytest.mark.parametrize('nitfo,ncyc', [(1000, 1000), (1000, 800),
                                            (1, 1), (2000, 1500)])
    def test_write_inp_refuses_first_order_throughout(self, nitfo, ncyc):
        import wing3d_otip as W
        lv = {k: v for k, v in W.M6_LEVELS['L1'].items() if k != 'h1_euler'}
        g = W.WingOTip(**lv, h1=2.0e-3, basin_depth=0.236)   # no mesh build
        with tempfile.TemporaryDirectory() as d:
            with pytest.raises(ValueError, match=r'FIRST ORDER'):
                g.write_inp(os.path.join(d, 'x.inp'), mach=0.8395, alpha=3.06,
                            re_mil=14.62, ivisc=0, nitfo=nitfo, ncyc=ncyc)

    @pytest.mark.parametrize('nitfo,ncyc', [(1000, 2000), (500, 1500),
                                            (0, 1000)])
    def test_write_inp_allows_a_real_switch(self, nitfo, ncyc):
        """The guard must not fire on a legal recipe.

        The object is unbuilt, so write_inp fails later for an unrelated
        reason; what is asserted is that it is NOT the nitfo refusal.  Testing
        only the raising half would leave a guard that rejects everything
        indistinguishable from a correct one.
        """
        import wing3d_otip as W
        lv = {k: v for k, v in W.M6_LEVELS['L1'].items() if k != 'h1_euler'}
        g = W.WingOTip(**lv, h1=2.0e-3, basin_depth=0.236)
        with tempfile.TemporaryDirectory() as d:
            try:
                g.write_inp(os.path.join(d, 'x.inp'), mach=0.8395, alpha=3.06,
                            re_mil=14.62, ivisc=0, nitfo=nitfo, ncyc=ncyc)
            except Exception as exc:                       # noqa: BLE001
                assert 'FIRST ORDER' not in str(exc), (
                    f'guard wrongly fired on the legal recipe '
                    f'nitfo={nitfo} ncyc={ncyc}')

    def test_2d_recipes_switch_to_second_order(self):
        """The four committed 2-D datasets are only valid because of this."""
        import cfl3d_runner as R
        for name, r in (('EULER_SOLVER', R.EULER_SOLVER),
                        ('RANS_SOLVER', R.RANS_SOLVER),
                        ('RANS_FALLBACK', R.RANS_FALLBACK)):
            assert r['NITFO'] < r['NCYC'], (
                f'{name} would run first order throughout: '
                f"NITFO={r['NITFO']} NCYC={r['NCYC']}")

    def test_recipe_validator_actually_fires(self):
        """Behavioural: a bad recipe must be rejected.

        Without this, ``test_2d_recipes_switch_to_second_order`` above could
        pass forever against a validator that had been silently gutted.
        """
        import cfl3d_runner as R
        assert R.validate_recipes() is True
        with pytest.raises(ValueError, match=r'FIRST ORDER'):
            R.validate_recipes({'BAD': dict(NITFO=1000, NCYC=1000)})


class TestImpliedOrder:
    """The estimator that decides which quantities get a published error bar."""

    #: the M6 Euler ladder's point counts (L1..L4)
    POINTS = np.array([769131.0, 2002.0e3, 4813.0e3, 11278.0e3])

    def _h(self):
        h = self.POINTS ** (-1.0 / 3.0)
        return h / h[0]

    @pytest.mark.parametrize('p_true', [2.0, 1.5, 1.0, 0.5])
    @pytest.mark.parametrize('triple', [0, 1])
    def test_recovers_a_known_order(self, p_true, triple):
        """Synthetic f = f_exact + C h^p must come back at order p.

        ★ The known answer comes from OUTSIDE the measurement, which is the
        whole point: an estimator tuned on the data it judges is a fit, not a
        test.
        """
        import generate_m6_reference as G
        h = self._h()
        vals = 1.0 + 0.3 * h ** p_true
        _, p = G.implied_order(h, vals, triple)
        assert p == pytest.approx(p_true, abs=0.02)

    def test_ratio_below_one_is_not_second_order(self):
        """The finding that "ratio < 1" was never calibrated.

        On this ladder a second-order quantity gives ~0.50 and a first-order
        one ~0.68, so a ratio of 0.9 -- comfortably "passing" the old test --
        is well under first order.
        """
        import generate_m6_reference as G
        h = self._h()
        r2, _ = G.implied_order(h, 1.0 + 0.3 * h ** 2.0, 0)
        r1, _ = G.implied_order(h, 1.0 + 0.3 * h ** 1.0, 0)
        assert r2 == pytest.approx(0.496, abs=0.01)
        assert r1 == pytest.approx(0.675, abs=0.01)
        # a ratio a hair under 1 implies a far lower order than either
        _, p_slow = G.implied_order(h, np.array(
            [1.0, 1.0 - 0.1, 1.0 - 0.1 - 0.09, 1.0 - 0.1 - 0.09 - 0.081]), 0)
        assert p_slow < 0.5, (
            'a ratio of 0.9 must not read as first order or better')

    def test_ratio_is_monotone_in_order(self):
        """Higher order => the deltas shrink faster => a smaller ratio."""
        import generate_m6_reference as G
        h = self._h()
        rs = [G.implied_order(h, 1.0 + 0.3 * h ** p, 0)[0]
              for p in (0.25, 0.5, 1.0, 1.5, 2.0)]
        assert all(a > b for a, b in zip(rs, rs[1:])), rs
