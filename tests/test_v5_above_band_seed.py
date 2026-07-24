"""GV5.1c synthetic seed-helper tests (a band-(a) suite leg).

Pure-map tests of the GV5.1c case runner's helper functions -- no mesh,
no heavy compute: the delta-perturbation mask, the calibration
bisection, the above-band triple filter (the GV5.1b _p_series,
imported), the regression slope, the contraction factors, and the
pooled window-verdict logic. The runner module is loaded by path (the
case-runner loader precedent). Binding text:
cases/analysis/v5_1c_above_band_window/PRE_REGISTRATION.md.
"""

import importlib.util
import os

import numpy as np

RUN = os.path.join(
    os.path.dirname(__file__), "..", "cases", "analysis",
    "v5_1c_above_band_window", "run.py")
spec = importlib.util.spec_from_file_location("gv51c_run", RUN)
gv51c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gv51c)


def _state_x(n_free=3, n_st=2, n_s=4):
    """A synthetic x = (phi_free, gamma, U.ravel())."""
    rng = np.arange(n_free + n_st + 6 * n_s, dtype=np.float64) + 1.0
    return rng, n_free, n_st, n_s


# ---------------------------------------------------------------------------
# perturb_delta
# ---------------------------------------------------------------------------


def test_perturb_scales_free_delta_only():
    x, n_free, n_st, n_s = _state_x()
    inflow = np.array([True, False, False, True])
    out = gv51c.perturb_delta(x, n_free, n_st, inflow, 0.5)
    U_in = x[n_free + n_st:].reshape(n_s, 6)
    U_out = out[n_free + n_st:].reshape(n_s, 6)
    # free rows: delta x1.5; pinned rows: bit-identical
    assert np.allclose(U_out[1, 0], 1.5 * U_in[1, 0])
    assert np.allclose(U_out[2, 0], 1.5 * U_in[2, 0])
    assert np.array_equal(U_out[0], U_in[0])
    assert np.array_equal(U_out[3], U_in[3])
    # every other column and the phi/gamma blocks untouched
    assert np.array_equal(U_out[:, 1:], U_in[:, 1:])
    assert np.array_equal(out[:n_free + n_st], x[:n_free + n_st])
    # the input is not mutated
    assert np.array_equal(
        x, np.arange(n_free + n_st + 6 * n_s, dtype=np.float64) + 1.0)


def test_perturb_zero_eps_is_identity():
    x, n_free, n_st, _ = _state_x()
    out = gv51c.perturb_delta(x, n_free, n_st,
                              np.array([True, False, True, False]), 0.0)
    assert np.array_equal(out, x)


# ---------------------------------------------------------------------------
# calibrate_eps
# ---------------------------------------------------------------------------


def test_calibrate_hits_window_on_monotone_map():
    # fbl(eps) = 1e-6 + 1e-3 eps^2: window [5e-2, 5e-1] reached at
    # eps in [~7.07, ~22.36]
    fn = lambda eps: 1.0e-6 + 1.0e-3 * eps ** 2
    eps, f, trace, ok, why = gv51c.calibrate_eps(fn, (5.0e-2, 5.0e-1))
    assert ok
    assert 5.0e-2 <= f <= 5.0e-1
    assert abs(f - (1.0e-6 + 1.0e-3 * eps ** 2)) < 1e-15
    assert len(trace) <= gv51c.MAX_CALIB_EVALS
    # determinism: same call, same eps
    eps2 = gv51c.calibrate_eps(fn, (5.0e-2, 5.0e-1))[0]
    assert eps == eps2


def test_calibrate_unreachable_window_aborts():
    eps, f, trace, ok, why = gv51c.calibrate_eps(lambda eps: 1.0e-6,
                                                 (5.0e-2, 5.0e-1))
    assert not ok
    assert eps is None
    assert "unreachable" in why


def test_calibrate_inf_steers_down():
    # inf above eps=1 (a throwing state): the bracket must steer below it
    def fn(eps):
        return float("inf") if eps > 1.0 else 1.0e-6 + 0.2 * eps
    eps, f, trace, ok, why = gv51c.calibrate_eps(fn, (5.0e-2, 5.0e-1))
    assert ok
    assert eps <= 1.0
    assert 5.0e-2 <= f <= 5.0e-1


# ---------------------------------------------------------------------------
# the above-band triple filter (the committed GV5.1b _p_series, imported)
# ---------------------------------------------------------------------------


def test_p_series_band_filter():
    e = [1.0e-1, 1.0e-2, 1.0e-4, 1.0e-8]
    band = 3.16e-5
    ps = gv51c.gb._p_series(e, band=band)
    assert [i for i, _ in ps] == [2, 3]
    assert all(abs(p - 2.0) < 1e-9 for _, p in ps)
    # a non-contraction inside a triple excludes that triple only:
    # (i=2: 2e-1 > 1e-1 breaks the triple; i=3 is a clean contraction)
    e2 = [1.0e-1, 2.0e-1, 1.0e-2, 1.0e-4]
    ps2 = gv51c.gb._p_series(e2, band=band)
    assert [i for i, _ in ps2] == [3]
    # triples whose predecessors sit below the band are excluded
    e3 = [1.0e-1, 1.0e-4, 1.0e-9, 1.0e-12]
    ps3 = gv51c.gb._p_series(e3, band=band)
    assert [i for i, _ in ps3] == [2]


# ---------------------------------------------------------------------------
# regression_slope / contraction_factors
# ---------------------------------------------------------------------------


def test_regression_slope_exact_on_quadratic():
    band = 1.0e-5
    e = [0.5]
    for _ in range(5):
        e.append(1.5 * e[-1] ** 2)
    slope, n = gv51c.regression_slope(e, band)
    assert n >= 3
    assert abs(slope - 2.0) < 1e-9
    # fewer than 3 above-band pairs -> None
    slope2, n2 = gv51c.regression_slope([1.0e-1, 1.0e-2, 1.0e-9], band)
    assert slope2 is None and n2 == 2


def test_contraction_factors():
    band = 1.0e-5
    e = [1.0e-1, 1.0e-2, 1.0e-3, 2.0e-3, 1.0e-9]
    cf = gv51c.contraction_factors(e, band)
    # i=1,2 contract above the band (1 decade each); i=3 expands (skip);
    # i=4 contracts and its predecessor 2e-3 > band (counts, ~6.3 decades)
    assert [i for i, _ in cf] == [1, 2, 4]
    assert all(abs(d - 1.0) < 1e-12 for _, d in cf[:2])


# ---------------------------------------------------------------------------
# window_verdict (the pooled band-(b) read)
# ---------------------------------------------------------------------------


def test_window_verdict_branches():
    # converged outright -> PASS on both levels
    v, _ = gv51c.window_verdict([], True, "medium")
    assert v == "PASS"
    # >= 3 triples, median in [1.5, 2.5]: medium PASS / FAIL binding,
    # coarse always RECORDED
    good = [(2, 2.0), (3, 1.9), (4, 2.1)]
    assert gv51c.window_verdict(good, False, "medium")[0] == "PASS"
    bad = [(2, 1.0), (3, 1.1), (4, 0.9)]
    assert gv51c.window_verdict(bad, False, "medium")[0] == "FAIL"
    assert gv51c.window_verdict(good, False, "coarse")[0] == "RECORDED"
    assert gv51c.window_verdict(bad, False, "coarse")[0] == "RECORDED"
    # < 3 triples -> the RECORDED fallback
    v, key = gv51c.window_verdict([(2, 2.0)], False, "medium")
    assert v == "RECORDED" and "fallback" in key
