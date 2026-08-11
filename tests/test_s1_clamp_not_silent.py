"""GS1.4: a clamped state is never reported as converged.

Why this exists
---------------
GS1.1 measured that the artificial-density floor can host a SPURIOUS machine-zero
solution: on one nozzle boundary-value problem with a provably unique answer, the
Newton converged to |R| = 9.2e-13 with the shock 40 cells away from the exact
position, and that state had 72 elements pinned at exactly `rho_floor`. A
four-decade A/B showed the physical branch is completely insensitive to the floor
(n_floored = 0 throughout) while the spurious branch only EXISTS for
floor >= 0.02 (docs/dev_phase_two/20260728-1640-s1-shock-bench.md 4.4).

A clamp is an algebraic override: where it binds, the element no longer
discretises the full-potential equation, so "the residual is zero" says nothing
about the flow. Every driver must therefore refuse to call such a state
converged. `solve_newton_lifting` and `solve_transonic_lifting` already did;
`solve_subsonic_lifting` -- the Picard SEED the other paths start from -- did
not, and GS1.4 fixed that. (The list also named the level-set
`solve_multivalued_newton`; that route was abandoned by ruling D5 and deleted in
phase 3, so it is dropped from the sentence -- the CONTRACT is unchanged.)

These tests lock the contract, not a number.
"""

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.newton import solve_newton_lifting
from pyfp3d.solve.picard import solve_subsonic_lifting


@pytest.fixture(scope="module")
def coarse(mesh_dir_2p5d):
    return cut_wake(read_mesh(mesh_dir_2p5d / "coarse.msh"))


@pytest.fixture(scope="module")
def mesh_dir_2p5d():
    from .conftest import REPO_ROOT
    d = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
    if not (d / "coarse.msh").exists():
        pytest.skip("naca0012_2.5d/coarse.msh not generated")
    return d


def test_picard_lifting_reports_clamped_flag(coarse):
    """The one-flag contract: `clamped` is present and agrees with the two
    counters it summarises."""
    mc, wc = coarse
    r = solve_subsonic_lifting(mc, wc, m_inf=0.5, alpha_deg=1.25,
                               n_picard_max=20)
    assert "clamped" in r, "solve_subsonic_lifting must report `clamped`"
    assert r["clamped"] == (r["n_limited"] > 0 or r["n_floored"] > 0)
    # a benign subsonic case must be clamp-free and converged
    assert not r["clamped"]
    assert r["converged"]


def test_newton_lifting_reports_clamped_flag(coarse):
    mc, wc = coarse
    r = solve_newton_lifting(mc, wc, m_inf=0.5, alpha_deg=1.25)
    assert "clamped" in r
    assert r["clamped"] == (r["n_limited"] > 0 or r["n_floored"] > 0)
    assert not r["clamped"] and r["converged"]


def test_picard_lifting_refuses_convergence_while_clamped(coarse):
    """Force the clamp to bind by lowering the speed limiter far below the
    flow's own peak: the driver must then NEVER report converged, however
    quiet the density lag gets.

    m_cap = 0.55 caps q^2 below the M 0.5 suction peak (local
    M ~ 0.68), so a wide patch is limited -- an unmistakably clamped state.
    """
    mc, wc = coarse
    r = solve_subsonic_lifting(mc, wc, m_inf=0.5, alpha_deg=1.25,
                               m_cap=0.55, n_picard_max=30)
    assert r["n_limited"] > 0, (
        "the m_cap probe did not actually clamp anything -- the test would "
        "be vacuous")
    assert r["clamped"]
    assert not r["converged"], (
        f"a clamped state was reported as converged "
        f"({r['n_limited']} limited / {r['n_floored']} floored)")


def test_clamped_seed_is_visible_to_callers(coarse):
    """A caller that only looks at `converged` must be safe: the clamped run
    above must not be mistakable for the clean one."""
    mc, wc = coarse
    clean = solve_subsonic_lifting(mc, wc, m_inf=0.5, alpha_deg=1.25,
                                   n_picard_max=20)
    clamped = solve_subsonic_lifting(mc, wc, m_inf=0.5, alpha_deg=1.25,
                                     m_cap=0.55, n_picard_max=30)
    assert clean["converged"] and not clean["clamped"]
    assert clamped["clamped"] and not clamped["converged"]
    # and the two states really are different flows, so the flag is not
    # decorative
    assert abs(float(clean["gamma"][0]) - float(clamped["gamma"][0])) > 1e-4
    assert not np.allclose(clean["phi"], clamped["phi"], rtol=1e-6)
