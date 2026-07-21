"""
B31 / GB31.2b tests: the tip_taper row blend for the PRESSURE Kutta
estimator (solve/newton.py, NewtonWorkspace class docstring).

The pressure row is homogeneous (F_j = P_j = station-mean(|q_u|^2 -
|q_l|^2) = 0), so the probe path's taper scaling (F_j = taper_j *
target_j - Gamma_j) does not port: scaling a homogeneous row is a no-op.
B31 blends the row with a Gamma pin instead (the B8 row-blend form, but
welding an EXPLICIT Newton unknown -- no LS-style re-leveling channel):

    F_j = taper_j * sigma_j * F_raw_j + (1 - taper_j) * s_j * Gamma_j

with sigma_j the SAME frozen P14 scaling and s_j = sign(diag D0)_j the
SIGN of the row's own Gamma-sensitivity, both frozen at the first
residual evaluation of the workspace (recorded as kutta_sigma /
kutta_weld_sign). The pin carries the row's own orientation because the
naive unsigned weld -(1-t)*Gamma is only correct under a dF/dGamma ~ -I
orientation, while the conforming meshes carry diag D > 0 (measured +80
on the NACA 2.5D coarse): there it amplifies mid-taper loading and the
blended Gamma block crosses zero at t = 1/(1+D). With the row's own
sign the blend unloads monotonically toward the pinned endpoint (at the
fixed-phi linearization Gamma_b = t * Gamma*, the probe taper's
Gamma_eff = taper * Gamma_Kutta semantics; the re-equilibrated R = 0
manifold unloads HARDER -- the pressure row's Gamma-slope along the
manifold is the Schur-complement slope, much smaller than the frozen
direct slope the weld carries -- measured kappa(0.7) ~ 0.14 on the NACA
2.5D coarse) and the raw Gamma block stays diag(D0)-conditioned at
every t. In raw units (F/sigma) the weld slope is s_j/sigma_j =
sign*max(|dj|, 0.1 median), the floored frozen dF_raw/dGamma diagonal
-- the B8 weld analog, dimensionally consistent.

Fast (NACA0012 coarse, one station -- the blend is per-station, so the
endpoints are separate one-station workspaces). Covers:
  (a) tip_taper=None / all-ones bit-identity for BOTH estimators;
  (b) blend algebra: value identity, the taper=1 / taper=0 endpoints,
      the exact elimination identity, and FD checks of dF/dphi_free and
      dF/dGamma through eval_residual against the analytic blocks
      (mirroring test_p14_te_pressure's FD protocol);
  (c) sigma freeze/record semantics under the blend.

Plus, on the gitignored M6 wing-body conforming coarse mesh (skip-
guarded), the REAL vanish_smooth taper profile: inboard (taper == 1)
stations bitwise equal the production row while the tip stations unload.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.constraints.wake import tip_taper_factors
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting

ALPHA = 2.0
M_INF = 0.5
UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR = 1.5, 0.95, 3.0, 0.05
T_BLEND = 0.7

REPO_ROOT = Path(__file__).parent.parent
M6_CONF_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6_wingbody_conforming"


@pytest.fixture(scope="module")
def naca_case():
    mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
                     / "coarse.msh")
    mc, wc = cut_wake(mesh)
    return mc, wc


def _solve(mc, wc, estimator, taper="default"):
    kw = {}
    if taper != "default":
        kw["tip_taper"] = taper
    return solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                kutta_estimator=estimator, precond="direct",
                                **kw)


@pytest.fixture(scope="module")
def blend_case(naca_case):
    """Converged M0.5 pressure+taper solve (taper = 0.7 on the single
    station) plus its workspace -- the state the algebra/FD checks run at."""
    mc, wc = naca_case
    taper = np.full(wc.n_stations, T_BLEND)
    r = _solve(mc, wc, "pressure", taper)
    return mc, wc, r


def _state(ws, r):
    phi_free = np.asarray(r["phi"])[:ws.n_red][ws.free].copy()
    gamma = np.asarray(r["gamma"], dtype=np.float64).copy()
    _, F, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                   M_CAP, RHO_FLOOR)
    return phi_free, gamma, F, state


# ------------------------------------------------- (a) default bit-identity


@pytest.mark.parametrize("estimator", ["probe", "pressure"])
def test_default_bit_identical(naca_case, estimator):
    """tip_taper omitted vs explicit all-ones: byte-identical phi, gamma
    (and the frozen sigma on the pressure path) -- the default-off gate
    for BOTH estimators."""
    mc, wc = naca_case
    r_def = _solve(mc, wc, estimator)
    r_ones = _solve(mc, wc, estimator, np.ones(wc.n_stations))
    assert r_def["converged"] and r_ones["converged"]
    assert np.array_equal(r_def["phi"], r_ones["phi"])
    assert np.array_equal(r_def["gamma"], r_ones["gamma"])
    if estimator == "pressure":
        assert np.array_equal(r_def["kutta_sigma"], r_ones["kutta_sigma"])
        assert not r_ones["workspace"]._taper_active


# ------------------------------------------------------- (b) blend algebra


def test_blend_value_algebra(blend_case):
    """F == sigma * taper * F_raw + (1 - taper) * s * Gamma with F_raw the
    unscaled pressure residual and s the frozen weld sign -- the
    eval_residual expression verified against the TEControlVolumes oracle
    at the converged state."""
    _, _, r = blend_case
    ws = r["workspace"]
    assert ws._taper_active
    phi_free, gamma, F, state = _state(ws, r)
    F_raw = ws.cvs.residual_stations(state["phi_cut"])
    expect = (ws.kutta_sigma * (ws.tip_taper * F_raw)
              + (1.0 - ws.tip_taper) * (ws.kutta_weld_sign * gamma))
    np.testing.assert_allclose(F, expect, rtol=1e-14, atol=1e-16)


def test_blend_unloads_gamma(blend_case, naca_case):
    """The blend unloads Gamma through the pressure row, monotonically as
    taper -> 0, preserving the sign: 0 < Gamma(t=0.3) < Gamma(t=0.7) <
    Gamma*. NOTE the quasi-2D single-station limit unloads HARDER than
    the probe taper's exact Gamma_b = t * Gamma* (measured kappa(0.7) ~
    0.14): along the R = 0 manifold the pressure row's Gamma-slope is
    the Schur-complement slope, much smaller than the frozen direct
    slope the weld carries, so mid-blend the pin dominates. The endpoint
    semantics (t = 1 production, t -> 0 pinned) are exact either way;
    the realistic multi-station comparison is the M6 coarse gate."""
    mc, wc, r = blend_case
    assert r["converged"]
    r_un = _solve(mc, wc, "pressure")
    r_hard = _solve(mc, wc, "pressure", np.full(wc.n_stations, 0.3))
    assert r_hard["converged"]
    g_un, g_bl, g_hard = (r_un["gamma"][0], r["gamma"][0],
                          r_hard["gamma"][0])
    assert g_un > 0.0
    assert 0.0 < g_hard < g_bl < g_un, (
        f"not a monotone signed unload: {g_hard:.4g} < {g_bl:.4g} < "
        f"{g_un:.4g} failed")
    assert g_bl < 0.9 * g_un, "the t = 0.7 blend must visibly unload"


def test_blend_endpoint_taper_zero(naca_case):
    """taper = 0: F = s * Gamma (pure pin at the row's own orientation)
    and the eliminated step is dGamma = -Gamma exactly (K_tilde == 0,
    F_tilde == -Gamma)."""
    mc, wc = naca_case
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA,
                         tip_taper=np.zeros(wc.n_stations),
                         kutta_estimator="pressure")
    ws.set_mach(M_INF)
    rng = np.random.default_rng(3)
    phi_free = rng.standard_normal(ws.n_free) * 0.01
    gamma = np.full(ws.n_st, 0.123)         # generic nonzero circulation
    _, F, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                   M_CAP, RHO_FLOOR)
    np.testing.assert_allclose(F, ws.kutta_weld_sign * gamma,
                               rtol=1e-14, atol=1e-16)
    K_tilde, F_tilde = ws.kutta_blocks(state, F)
    x = rng.standard_normal(ws.n_free)
    np.testing.assert_allclose(K_tilde @ x, np.zeros(ws.n_st),
                               rtol=0, atol=1e-15)
    np.testing.assert_allclose(F_tilde, -gamma, rtol=1e-12, atol=1e-14)


def test_blend_endpoint_taper_one_is_production(naca_case):
    """taper = 1 (all-ones): the blend never activates -- F is the
    production sigma * F_raw row bitwise (the inactive path; the
    mixed-profile inboard clause is the M6 test below)."""
    mc, wc = naca_case
    taper = np.ones(wc.n_stations)
    r = _solve(mc, wc, "pressure", taper)
    ws = r["workspace"]
    phi_free, gamma, F, state = _state(ws, r)
    F_raw = ws.cvs.residual_stations(state["phi_cut"])
    assert np.array_equal(F, ws.kutta_sigma * F_raw)


def test_blend_elimination_identity(blend_case):
    """The eliminated blocks satisfy the exact block-row equation
    K dphi + D dgamma = -F with the ANALYTIC blended blocks
    K = sigma * taper * Kp_free,
    D = sigma * taper * D_raw + diag((1 - taper) * s)
    and dgamma = K_tilde @ dphi + F_tilde -- LU-exact, no FD noise."""
    _, _, r = blend_case
    ws = r["workspace"]
    phi_free, gamma, F, state = _state(ws, r)
    t = ws.tip_taper
    Kp_cut, D_raw = ws.cvs.newton_rows(state["phi_cut"])
    Kp_free = (Kp_cut @ ws.con.T).tocsr()[:, ws.free]
    K_code = ws.kutta_sigma[:, None] * (t[:, None] * Kp_free.toarray())
    D_code = (ws.kutta_sigma[:, None] * (t[:, None] * D_raw)
              + np.diag((1.0 - t) * ws.kutta_weld_sign))
    K_tilde, F_tilde = ws.kutta_blocks(state, F)
    rng = np.random.default_rng(5)
    for _ in range(3):
        dphi = rng.standard_normal(ws.n_free)
        dgamma = K_tilde @ dphi + F_tilde
        row = K_code @ dphi + D_code @ dgamma + F
        assert np.abs(row).max() < 1e-10


def test_blend_jacobian_fd_phi(blend_case):
    """dF/dphi_free = sigma * taper * Kp_free vs central FD of F through
    eval_residual (F exactly quadratic -- the tapered analog of p14's
    test_newton_pressure_kutta_row_fd)."""
    _, _, r = blend_case
    ws = r["workspace"]
    phi_free, gamma, _, state = _state(ws, r)
    Kp_cut, _ = ws.cvs.newton_rows(state["phi_cut"])
    Kp_free = (Kp_cut @ ws.con.T).tocsr()[:, ws.free]
    coef = ws.kutta_sigma * ws.tip_taper

    rng = np.random.default_rng(7)
    eps = 1e-4
    for _ in range(3):
        delta = rng.standard_normal(ws.n_free)
        delta /= np.abs(delta).max()
        _, F_p, _ = ws.eval_residual(phi_free + eps * delta, gamma,
                                     UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
        _, F_m, _ = ws.eval_residual(phi_free - eps * delta, gamma,
                                     UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
        fd = (F_p - F_m) / (2.0 * eps)
        exact = coef * (Kp_free @ delta)
        scale = max(np.abs(fd).max(), 1e-30)
        assert np.abs(exact - fd).max() / scale < 1e-8


def test_blend_jacobian_fd_gamma(blend_case):
    """dF/dGamma = sigma * taper * D_raw + diag((1 - taper) * s) vs central
    FD through eval_residual (the tapered analog of p14's gamma-column
    FD: the explicit frozen-slope pin derivative plus the slave-jump
    chain, no far-field term)."""
    _, _, r = blend_case
    ws = r["workspace"]
    phi_free, gamma, _, state = _state(ws, r)
    _, D_raw = ws.cvs.newton_rows(state["phi_cut"])
    t = ws.tip_taper
    D_code = (ws.kutta_sigma[:, None] * (t[:, None] * D_raw)
              + np.diag((1.0 - t) * ws.kutta_weld_sign))

    eps = 1e-6
    for j in range(ws.n_st):
        dg = np.zeros(ws.n_st)
        dg[j] = eps
        _, F_p, _ = ws.eval_residual(phi_free, gamma + dg, UPWIND_C,
                                     M_CRIT, M_CAP, RHO_FLOOR)
        _, F_m, _ = ws.eval_residual(phi_free, gamma - dg, UPWIND_C,
                                     M_CRIT, M_CAP, RHO_FLOOR)
        fd = (F_p - F_m) / (2.0 * eps)
        scale = max(np.abs(fd).max(), 1e-30)
        assert np.abs(D_code[:, j] - fd).max() / scale < 1e-7


def test_blend_weld_slope_is_frozen_row_slope(naca_case):
    """The raw-units weld slope s_j/sigma_j equals the row's own floored
    frozen Gamma-slope sign(dj)*max(|dj|, 0.1 median) at the seed state
    (the B8-analog unit bridge). NOTE: unlike the untapered pressure
    path, post-hoc sigma rescaling is NOT a symmetry of the blend -- the
    pressure:weld ratio is physical, so sigma must be (and is) frozen."""
    mc, wc = naca_case
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA,
                         tip_taper=np.full(wc.n_stations, T_BLEND),
                         kutta_estimator="pressure")
    ws.set_mach(M_INF)
    rng = np.random.default_rng(23)
    phi_free = rng.standard_normal(ws.n_free) * 0.01
    gamma = np.full(ws.n_st, 0.05)
    _, _, state0 = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                    M_CAP, RHO_FLOOR)
    D0 = ws.cvs.gamma_jacobian(state0["phi_cut"], mode="exact")
    dj = np.diag(D0)
    adj = np.abs(dj)
    floored = np.sign(dj) * np.maximum(adj, 0.1 * np.median(adj))
    assert np.all(ws.kutta_weld_sign == np.sign(dj))
    np.testing.assert_allclose(ws.kutta_weld_sign / ws.kutta_sigma,
                               floored, rtol=1e-14, atol=0.0)


def test_blend_taper_zero_solve_pins_gamma(naca_case):
    """End-to-end: a taper = 0 pressure solve converges and pins the
    circulation to ~0 (the homogeneous-row analog of the probe path's
    Gamma = 0 station)."""
    mc, wc = naca_case
    r = _solve(mc, wc, "pressure", np.zeros(wc.n_stations))
    assert r["converged"]
    assert np.abs(r["gamma"]).max() < 1e-6


# -------------------------------------------------- (c) sigma freeze record


def test_sigma_frozen_once_and_recorded(naca_case):
    """Under the blend sigma is STILL frozen at the first residual
    evaluation of the workspace (the seed state), follows the P14 freeze
    formula 1/max(|D_jj|, 0.1 median), is not refit at later states, and
    is recorded on the result dict with the flip count."""
    mc, wc = naca_case
    taper = np.full(wc.n_stations, T_BLEND)
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA, tip_taper=taper,
                         kutta_estimator="pressure")
    ws.set_mach(M_INF)
    assert ws.kutta_sigma is None
    rng = np.random.default_rng(17)
    phi_free = rng.standard_normal(ws.n_free) * 0.01
    gamma = np.full(ws.n_st, 0.05)
    _, _, state0 = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                    M_CAP, RHO_FLOOR)
    sigma0 = ws.kutta_sigma.copy()
    sign0 = ws.kutta_weld_sign.copy()
    assert np.all(np.abs(sign0) == 1.0)
    # freeze formula against the FIRST state's exact D (newton.py)
    D0 = ws.cvs.gamma_jacobian(state0["phi_cut"], mode="exact")
    adj = np.abs(np.diag(D0))
    np.testing.assert_allclose(
        sigma0, 1.0 / np.maximum(adj, 0.1 * np.median(adj)),
        rtol=1e-14, atol=0.0)
    # later states refit NEITHER sigma nor the weld sign
    for _ in range(2):
        ws.eval_residual(phi_free + 0.01 * rng.standard_normal(ws.n_free),
                         gamma + 0.01, UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
        assert np.array_equal(ws.kutta_sigma, sigma0)
        assert np.array_equal(ws.kutta_weld_sign, sign0)
    # recorded on the solve result
    r = _solve(mc, wc, "pressure", taper)
    assert np.array_equal(r["kutta_sigma"], r["workspace"].kutta_sigma)
    assert r["kutta_sigma"] is not None
    assert isinstance(r["kutta_sigma_sign_flips"], int)
    assert np.array_equal(r["kutta_weld_sign"],
                          r["workspace"].kutta_weld_sign)


def test_unknown_estimator_still_raises(naca_case):
    """The NotImplementedError is gone but the estimator ValueError stays."""
    mc, wc = naca_case
    with pytest.raises(ValueError, match="kutta_estimator"):
        NewtonWorkspace(mc, wc, alpha_deg=ALPHA,
                        tip_taper=np.full(wc.n_stations, 0.5),
                        kutta_estimator="typo")


# ----------------------- M6 wing-body conforming coarse (gitignored mesh):
# the real multi-station vanish_smooth profile -- inboard stations bitwise
# production, tip stations blended, tip-most taper ~ 0.


@pytest.fixture(scope="module")
def m6_conf():
    path = M6_CONF_DIR / "coarse.msh"
    if not path.exists():
        pytest.skip(f"{path} missing (run cases/meshes/"
                    "onera_m6_wingbody_conforming/"
                    "generate_onera_m6_wingbody_conforming.py)")
    return cut_wake(read_mesh(str(path)))


def test_m6_real_taper_profile_algebra(m6_conf):
    """GB31.2b production taper (vanish_smooth, r_c = 0.05 * B_SEMI) on the
    real geometry, evaluated at a generic state (no solve): stations with
    taper == 1 carry the production row BITWISE, tapered stations match
    the blend formula, and the tip-most station sits at taper ~ 0 (the
    compact ~5%-span unload)."""
    from pyfp3d.constraints.dirichlet import freestream_phi
    from pyfp3d.meshgen.wing3d import B_SEMI

    mc, wc = m6_conf
    taper = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth",
                              0.05 * B_SEMI)
    assert np.any(taper == 1.0) and taper[-1] < 0.05
    ws = NewtonWorkspace(mc, wc, alpha_deg=3.06,
                         farfield_spanwise_gamma=True, tip_taper=taper,
                         kutta_estimator="pressure")
    ws.set_mach(M_INF)
    phi_free = freestream_phi(mc.nodes[:ws.n_red], 3.06, 1.0)[ws.free]
    gamma = np.full(ws.n_st, 0.05)
    _, F, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                   M_CAP, RHO_FLOOR)
    F_raw = ws.cvs.residual_stations(state["phi_cut"])
    expect = (ws.kutta_sigma * (taper * F_raw)
              + (1.0 - taper) * (ws.kutta_weld_sign * gamma))
    np.testing.assert_allclose(F, expect, rtol=1e-14, atol=1e-16)
    inb = taper == 1.0
    assert np.array_equal(F[inb], (ws.kutta_sigma * F_raw)[inb])
