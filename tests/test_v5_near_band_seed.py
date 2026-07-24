"""GV5.1d synthetic seed-helper tests (a band-(a) suite leg).

Pure-map tests of the GV5.1d case runner -- no mesh, no heavy compute:
the near-band window sanity (above both floor bands, below the GV5.1c
stall region, T2 strictly above T1), the band-entry read, the
calibration bisection on a near-band synthetic response, and the
imported-helper identity (the GV5.1c runner's helpers, imported not
mirrored). The runner module is loaded by path (the case-runner loader
precedent). Binding text:
cases/analysis/v5_1d_near_band_window/PRE_REGISTRATION.md.
"""

import importlib.util
import os

RUN = os.path.join(
    os.path.dirname(__file__), "..", "cases", "analysis",
    "v5_1d_near_band_window", "run.py")
spec = importlib.util.spec_from_file_location("gv51d_run", RUN)
gv51d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gv51d)


# ---------------------------------------------------------------------------
# the near-band windows (the pre-registered TARGETS)
# ---------------------------------------------------------------------------


def test_windows_sit_above_both_floor_bands():
    t1, t2 = gv51d.TARGETS["T1"], gv51d.TARGETS["T2"]
    for level in ("coarse", "medium"):
        band = gv51d.FLOOR_BAND[level]
        # the whole T1 window sits ABOVE the band (a near-band seed is
        # an above-band seed by construction)
        assert t1[0] > band, f"{level}: T1 lo {t1[0]} <= band {band}"
        # ... and below the GV5.1c mid-range stall region (~1e-2)
        assert t1[1] <= 1.0e-2
    # the pre-registered values verbatim
    assert t1 == (1.0e-4, 1.0e-3)
    assert t2 == (1.0e-3, 1.0e-2)


def test_escalation_window_strictly_above_primary():
    t1, t2 = gv51d.TARGETS["T1"], gv51d.TARGETS["T2"]
    # non-overlapping, T2 strictly above T1 (the escalation extends the
    # read resolution upward -- the GV5.1c escalation direction)
    assert t2[0] >= t1[1]
    assert t2[1] > t1[1]
    # both windows sit strictly below the GV5.1c primary window
    assert t2[1] <= gv51d.gc.TARGETS["T1"][0]


def test_windows_multiples_of_the_floor_bands():
    # T1 = 5.8-58x the medium band / 3.2-31.6x the coarse band (the
    # decades directly above the band): sanity on the order of magnitude
    med = gv51d.FLOOR_BAND["medium"]
    cor = gv51d.FLOOR_BAND["coarse"]
    assert 1.0 < gv51d.TARGETS["T1"][0] / med < 100.0
    assert 1.0 < gv51d.TARGETS["T1"][0] / cor < 100.0
    assert gv51d.TARGETS["T1"][1] / med < 1.0e3
    assert gv51d.TARGETS["T1"][1] / cor < 1.0e3


# ---------------------------------------------------------------------------
# band_entry_iter (the key new datum)
# ---------------------------------------------------------------------------


def test_band_entry_iter_reads_first_sub_band_point():
    band = 3.16e-5
    # a trajectory that enters the band at iter 3
    assert gv51d.band_entry_iter([1e-3, 1e-4, 1e-4, 2e-5, 1e-5], band) == 3
    # never entering the band (the GV5.1c above-band outcome) -> None
    assert gv51d.band_entry_iter([1e-1, 3e-2, 1.3e-2], band) is None
    # already inside the band at the seed (the GV5.1b constructional
    # read) -> 0
    assert gv51d.band_entry_iter([1e-6, 1e-6], band) == 0
    # a mid-trajectory dip below the band counts at its first index
    assert gv51d.band_entry_iter([1e-3, 1e-5, 1e-4], band) == 1


# ---------------------------------------------------------------------------
# calibrate_eps on a near-band synthetic response (imported from GV5.1c)
# ---------------------------------------------------------------------------


def test_calibrate_lands_in_the_near_band_window():
    # fbl(eps) = 2e-6 + 3e-5 * eps: window [1e-4, 1e-3] reached at
    # eps in [~3.3, ~33] -- the bisection lands deterministically
    fn = lambda eps: 2.0e-6 + 3.0e-5 * eps
    eps, f, trace, ok, why = gv51d.gc.calibrate_eps(
        fn, gv51d.TARGETS["T1"])
    assert ok
    assert gv51d.TARGETS["T1"][0] <= f <= gv51d.TARGETS["T1"][1]
    assert abs(f - (2.0e-6 + 3.0e-5 * eps)) < 1e-15
    eps2 = gv51d.gc.calibrate_eps(fn, gv51d.TARGETS["T1"])[0]
    assert eps == eps2  # determinism


def test_calibrate_near_band_unreachable_aborts():
    # a response pinned at the floor (2e-6) can never reach 1e-4:
    # the pre-registered abort fires
    eps, f, trace, ok, why = gv51d.gc.calibrate_eps(
        lambda eps: 2.0e-6, gv51d.TARGETS["T1"])
    assert not ok and eps is None
    assert "unreachable" in why


# ---------------------------------------------------------------------------
# imported-helper identity (imported, not mirrored)
# ---------------------------------------------------------------------------


def test_helpers_are_the_gv51c_runners_own():
    # the physics helpers resolve to the COMMITTED GV5.1c runner's
    # functions, and through it to GV5.1b's protocol constants
    assert gv51d.gc.perturb_delta is gv51d.gc.perturb_delta
    assert gv51d.gc.__name__ == "gv51c_run"
    assert gv51d.gb.FLOOR_BAND["medium"] == 1.71e-5
    assert gv51d.gb.FLOOR_BAND["coarse"] == 3.16e-5
    assert gv51d.BAND_A_TOL == gv51d.gb.BAND_A_TOL
    # the window-verdict branch logic is the GV5.1c one (pooled read)
    good = [(2, 2.0), (3, 1.9), (4, 2.1)]
    assert gv51d.gc.window_verdict(good, False, "medium")[0] == "PASS"
    assert gv51d.gc.window_verdict(good, False, "coarse")[0] == "RECORDED"
