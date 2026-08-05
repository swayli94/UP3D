"""The cold-start Picard-seed fallback (user ruling 2026-08-05).

Mechanism it exists for, MEASURED in docs/dev_phase_two/20260805-2200-seed-exposure.md:
the 2026-08-02 flip of `n_picard_seed`'s default from 5 to 0 does not fail on its own.
It fails on a conjunction -- no seed, AND a cold start directly at a supercritical
M_inf, AND a mesh fine enough to resolve the supersonic pocket. NACA0012 M0.80 medium
died with M_max exactly at m_cap from Newton step 3; the same case through a Mach ramp
at the same seed converges, because the previous level's converged solution does the
seed's job.

So the fallback retries ONCE with a seed when, and only when, neither a seed nor a
warm start is present and the attempt ended clamped.

The gating is what these tests lock, and they lock it CHEAPLY by lowering `m_cap` on
the coarse mesh instead of paying for the medium transonic case: a low cap makes a
1-second solve clamp on demand, which is what the three gate conditions need in order
to be exercised at all. The expensive real recovery is gated separately at the bottom.
"""

import os

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.newton import _SEED_FALLBACK, solve_newton_lifting

from .conftest import REPO_ROOT

MESH = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh"
#: enough dissipation and Newton budget to be a normal solve, small enough to be cheap
KW = dict(alpha_deg=1.25, upwind_c=1.5, m_crit=0.95, precond="direct",
          n_newton_max=30)
#: a cap this low clamps on the coarse mesh at M0.80 within a couple of steps. It is a
#: TEST INSTRUMENT for reaching the clamped branch, not a physical setting -- the
#: production cap is 3.0 and the failure this fallback addresses happens there.
M_CAP_LOW = 1.02


@pytest.fixture(scope="module")
def cut():
    return cut_wake(read_mesh(MESH))


def test_no_fallback_on_a_clean_solve(cut):
    """The success path is untouched: a solve that converges unclamped must not fire
    the fallback. This is guaranteed by construction -- the trigger requires clamped
    cells and GS1.4 already refuses to report a clamped state as converged -- and
    asserted here so a future change to either half is caught."""
    mc, wc = cut
    r = solve_newton_lifting(mc, wc, m_inf=0.50, **KW)
    assert r["converged"] and r["n_limited"] == 0 and r["n_floored"] == 0
    assert r["seed_fallback"]["fired"] is False
    assert r["seed_fallback"]["seed"] is None


def test_fallback_fires_on_a_clamped_cold_start(cut):
    """The trigger itself. With the cap lowered the cold seed-0 attempt clamps, so the
    retry must happen -- and the record must carry the FIRST attempt's numbers, not
    only the retry's, because a fallback that hides what it replaced cannot be
    audited."""
    mc, wc = cut
    r = solve_newton_lifting(mc, wc, m_inf=0.80, m_cap=M_CAP_LOW, **KW)
    fb = r["seed_fallback"]
    assert fb["fired"] is True, (
        "premise: this configuration must clamp on the first attempt, else the test "
        "is not exercising the trigger")
    assert fb["seed"] == _SEED_FALLBACK
    if fb["accepted"]:
        assert r["converged"]
        assert fb["first_n_limited"] > 0 or fb["first_n_floored"] > 0
    else:
        # the retry failed too: the ORIGINAL must be returned, so the result's own
        # diagnostics describe the default path rather than the fallback's
        assert not r["converged"]
        assert "retry_accept_reason" in fb


def test_a_warm_start_excludes_the_fallback(cut):
    """A ramp's intermediate levels warm-start from the last converged level, and
    replacing that with a Picard seed would make them WORSE. So `phi_init` must
    suppress the fallback even when the attempt clamps."""
    mc, wc = cut
    seed = solve_newton_lifting(mc, wc, m_inf=0.50, **KW)
    assert seed["converged"], "the warm start itself must be a solution"
    r = solve_newton_lifting(mc, wc, m_inf=0.80, m_cap=M_CAP_LOW,
                             phi_init=seed["phi"], gamma_init=seed["gamma"], **KW)
    assert r["n_limited"] > 0 or r["n_floored"] > 0, (
        "premise: this leg must clamp, else it cannot show the suppression")
    assert r["seed_fallback"]["fired"] is False


def test_an_explicit_seed_excludes_the_fallback(cut):
    """A caller who asked for a seed already has one; retrying would be pointless work
    and would hide that their seed was insufficient."""
    mc, wc = cut
    r = solve_newton_lifting(mc, wc, m_inf=0.80, m_cap=M_CAP_LOW,
                             n_picard_seed=3, **KW)
    assert r["n_limited"] > 0 or r["n_floored"] > 0, "premise: must clamp"
    assert r["seed_fallback"]["fired"] is False


def test_the_retry_cannot_recurse(cut):
    """The retry runs with `_seed_retry=True`, so if it also clamps it must not spawn
    another. Asserted through the public result: a fired-and-rejected record proves
    the retry returned rather than recursing (an infinite recursion would raise)."""
    mc, wc = cut
    r = solve_newton_lifting(mc, wc, m_inf=0.80, m_cap=M_CAP_LOW, **KW)
    assert isinstance(r["seed_fallback"]["fired"], bool)


@pytest.mark.skipif(not os.environ.get("PYFP3D_TRANSONIC_GATES"),
                    reason="heavy: NACA medium transonic (~50 s)")
def test_the_real_case_the_fallback_exists_for():
    """★ The measured failure itself: NACA0012 M0.80 alpha 1.25 on the MEDIUM mesh,
    cold, seed 0 -- the M1 gate's own condition.

    Before the fallback this returned |R| = 3.29e-02 with M_max exactly at m_cap
    (7265 limited / 758 floored) and cl going negative. With it, the recovered state
    must match what an explicit seed 5 produces, because that is what the fallback
    does -- the recovery is not a new answer, it is the pre-2026-08-02 answer.
    """
    p = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "medium.msh"
    if not p.exists():
        pytest.skip("naca0012_2.5d/medium.msh not generated")
    mc, wc = cut_wake(read_mesh(p))
    kw = dict(alpha_deg=1.25, upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6,
              freeze_refresh_max=8, precond="direct", direct_refactor_every=4,
              n_newton_max=80)
    r = solve_newton_lifting(mc, wc, m_inf=0.80, **kw)
    fb = r["seed_fallback"]
    assert fb["fired"] and fb["accepted"], f"fallback did not recover: {fb}"
    assert r["converged"] and r["n_limited"] == 0 and r["n_floored"] == 0
    assert float(r["residual_history"][-1]) < 1e-9
    ref = solve_newton_lifting(mc, wc, m_inf=0.80, n_picard_seed=_SEED_FALLBACK,
                               **kw)
    assert ref["converged"]
    assert np.allclose(np.asarray(r["phi"]), np.asarray(ref["phi"]),
                       rtol=1e-12, atol=0.0), (
        "the fallback must reproduce the explicit-seed answer -- it is the same "
        "solve, so anything else means the retry is not passing the settings through")
