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


class TestRunStatus:
    """A CFL3D summary is NOT evidence that the solve finished.

    ★★★ Measured: the M6 L2 rung died with `NaN detected after residual
    evaluation, block 1 cycle 540` and still wrote a SUMMARY OF FORCES giving
    cl = +0.339608.  A runner that checked for the summary reported it "ok" and
    would have published a diverging state as reference data.

    ★★ Three wrong detectors preceded the right one, and each was caught only
    by running it against cases whose outcome was already known:
      1. "NaN appears in cfl3d.out"  -- the message is written to cfl3d.error;
      2. "a SUMMARY block exists"    -- failed runs write one too;
      3. "cfl3d.error is non-empty"  -- SUCCESSFUL runs write to it as well,
         and in fact the normal-termination banner appears nowhere else.
    The authoritative field is the error CODE.  These synthetic cases lock all
    four behaviours so no future simplification can reintroduce any of them.
    """

    OK_ERR = ' error code:\n   0\n\n execution terminated normally\n'
    BAD_ERR = (' error code:\n  -1\n\n abnormal termination due to cfl3d '
               'error check\n (error message follows)\n\n dump of unit 11 '
               '(main output) buffer:\n\n NaN detected after residual '
               'evaluation, block   1 cycle  540\n')
    INP = ('      NCYC    MGLEVG     NEMGL     NITFO\n'
           '      1500         1         0       500\n'
           '      1500         2         0       500\n'
           '      1500         3         0       500\n')

    def _case(self, tmp_path, err, last_it, with_summary=True):
        d = tmp_path
        (d / 'cfl3d.inp').write_text(self.INP)
        rows = ['Variables = LV BLK IT res']
        for lv in (1, 2, 3):
            for it in range(1, 1501):
                g = (lv - 1) * 1500 + it
                if g > last_it:
                    break
                rows.append(f'{lv:6d}{1:4d}{g:7d}  0.1E-06  0.1E-06  '
                            f'0.33E+00  0.15E-01  0.15E-01  0.0E+00  0.0E+00')
        (d / 'clcd_total.dat').write_text('\n'.join(rows) + '\n')
        if err is not None:
            (d / 'cfl3d.error').write_text(err)
        out = 'blah\n'
        if with_summary:
            # ★ verbatim the shape CFL3D writes, so read_summary really parses
            #   it -- the trap is only live if the reader would have returned a
            #   number for this dead run.
            out += ('SUMMARY OF FORCES AND MOMENTS - ALL GLOBAL BLOCKS\n\n'
                    '          CL                CD               CDp'
                    '               CDv\n'
                    '  0.33960762016E+00  0.14795440792E-01'
                    '  0.14795440792E-01  0.00000000000E+00\n'
                    '          CZ                CY               CX'
                    '            wetted area\n'
                    '  0.1E-01  0.3E+00  -0.2E-02  0.15E+01\n')
        (d / 'cfl3d.out').write_text(out)
        return d

    def test_completed_run_is_ok(self, tmp_path):
        import generate_m6_reference as G
        self._case(tmp_path, self.OK_ERR, 4500)
        ok, why, last, want = G.run_status(tmp_path)
        assert ok, why
        assert (last, want) == (4500, 4500)

    def test_diverged_run_is_rejected_despite_its_summary(self, tmp_path):
        """The exact L2 failure: died at 3539 of 4500, summary present."""
        import generate_m6_reference as G
        self._case(tmp_path, self.BAD_ERR, 3539, with_summary=True)
        ok, why, last, want = G.run_status(tmp_path)
        assert not ok
        assert 'DIVERGED' in why and 'cycle  540' in why
        # and the summary really is there, i.e. the trap is live
        got = G.read_summary(tmp_path / 'cfl3d.out')
        assert got['CL'] == pytest.approx(0.33960762016), (
            'the trap must be live: the reader returns a number for this '
            'DEAD run, which is exactly why run_status is required')

    def test_nonempty_error_file_is_not_failure(self, tmp_path):
        """Guards against the third wrong detector."""
        import generate_m6_reference as G
        self._case(tmp_path, self.OK_ERR, 4500)
        assert (tmp_path / 'cfl3d.error').stat().st_size > 0
        assert G.run_status(tmp_path)[0] is True

    def test_short_run_with_code_zero_is_rejected(self, tmp_path):
        """Code 0 but the counter never reached mseq*ncyc."""
        import generate_m6_reference as G
        self._case(tmp_path, self.OK_ERR, 3000)
        ok, why, last, want = G.run_status(tmp_path)
        assert not ok, why
        assert (last, want) == (3000, 4500), (last, want)

    def test_still_running_is_not_reported_as_success(self, tmp_path):
        import generate_m6_reference as G
        self._case(tmp_path, None, 3490)
        ok, why, _, _ = G.run_status(tmp_path)
        assert not ok and 'STILL RUNNING' in why
