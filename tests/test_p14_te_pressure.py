"""
P14 stage-1 tests: the wall-adjacent-CV pressure-equality Kutta estimator
(constraints/te_pressure.py) and its solve_laplace_lifting wiring.

Fast (NACA0012 coarse only, two Laplace solves in the module fixture).
Covers the G14.1 construction/non-degeneracy clauses, the exact-factorization
identity, an FD guard of dF/dGamma (exact to roundoff -- F is exactly
quadratic in (phi, Gamma)), the G14.3 NACA lift-preservation band, and the
G14.4 probe-default bit-identity at the Laplace driver.
"""

import numpy as np
import pytest

from pyfp3d.constraints.dirichlet import farfield_dirichlet
from pyfp3d.constraints.te_pressure import TEControlVolumes
from pyfp3d.constraints.wake import WakeConstraint, kutta_targets
from pyfp3d.kernels.residual import assemble_stiffness_matrix
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting
from pyfp3d.solve.picard import solve_laplace_lifting

ALPHA = 2.0
M_INF = 0.5
UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR = 1.5, 0.95, 3.0, 0.05


@pytest.fixture(scope="module")
def naca_case():
    from .conftest import REPO_ROOT

    mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
                     / "coarse.msh")
    mc, wc = cut_wake(mesh)
    dir_nodes, _ = farfield_dirichlet(mc, wc, ALPHA, np.zeros(wc.n_stations))
    cvs = TEControlVolumes(mc, wc, dirichlet_nodes=dir_nodes)
    con = WakeConstraint(assemble_stiffness_matrix(mc.nodes, mc.elements), wc)
    r_probe = solve_laplace_lifting(mc, wc, alpha_deg=ALPHA)
    return mc, wc, cvs, con, r_probe


def test_cv_construction_invariants(naca_case):
    """G14.1 construction clause: the constructor's asserts all passed
    (exact wall-face ownership, probe-membership side identity, zero
    Dirichlet contact -- a raise would have failed the fixture) and the
    fans are two-sided non-empty."""
    _, wc, cvs, _, _ = naca_case
    fs = cvs.fan_stats()
    assert fs["n_te"] == len(wc.te_nodes) == 2       # quasi-2D: 2 TE nodes
    assert fs["n_st"] == wc.n_stations == 1          # ... in 1 station
    assert fs["fan_u_min"] >= 1 and fs["fan_l_min"] >= 1


def test_factorization_identity(naca_case):
    """s_bar.(q_u - q_l) == |q_u|^2 - |q_l|^2 exactly (the te_kutta_coo
    factorization) -- so the frozen-mean implied target has the nonlinear
    pressure equality as its exact fixed point."""
    _, wc, cvs, _, r = naca_case
    phi = r["phi"]
    qu, ql = cvs.te_velocities(phi)
    lhs_node = np.einsum("kd,kd->k", qu + ql, qu - ql)
    rhs_node = np.einsum("kd,kd->k", qu, qu) - np.einsum("kd,kd->k", ql, ql)
    np.testing.assert_allclose(lhs_node, rhs_node, rtol=0, atol=1e-14)
    counts = np.bincount(wc.te_station, minlength=wc.n_stations)
    F_direct = np.bincount(wc.te_station, weights=rhs_node,
                           minlength=wc.n_stations) / counts
    np.testing.assert_allclose(cvs.residual_stations(phi), F_direct,
                               rtol=1e-14, atol=1e-16)


def test_gamma_jacobian_fd_exact(naca_case):
    """Central FD of the station residual w.r.t. Gamma at FIXED phi_red vs
    the dense exact D: F is exactly quadratic in Gamma, so this is exact to
    roundoff and validates the whole _rows_cut @ G chain."""
    _, wc, cvs, con, r = naca_case
    gamma = r["gamma"]
    phi_red = r["phi"][: wc.n_nodes_orig]
    D = cvs.gamma_jacobian(con.expand(phi_red, gamma), mode="exact")
    eps = 1e-6
    n_st = wc.n_stations
    for j in range(n_st):
        e = np.zeros(n_st)
        e[j] = eps
        Fp = cvs.residual_stations(con.expand(phi_red, gamma + e))
        Fm = cvs.residual_stations(con.expand(phi_red, gamma - e))
        fd = (Fp - Fm) / (2 * eps)
        np.testing.assert_allclose(D[:, j], fd, rtol=1e-8, atol=1e-9)


def test_gamma_nondegenerate(naca_case):
    """G14.1 non-degeneracy clause: dF/dGamma bounded away from zero at a
    converged state (the recovered two-sided velocity responds to Gamma)."""
    _, wc, cvs, con, r = naca_case
    phi_cut = con.expand(r["phi"][: wc.n_nodes_orig], r["gamma"])
    D = cvs.gamma_jacobian(phi_cut, mode="exact")
    assert np.abs(np.diag(D)).min() > 1.0
    M = cvs.gamma_jacobian(phi_cut, mode="frozen")
    assert np.linalg.cond(M) < 100.0


def test_pressure_laplace_converges_and_preserves_lift(naca_case):
    """G14.3 NACA leg: the pressure-estimator Laplace solve converges in
    the existing update budget and lands within 2% of the probe-path Gamma
    (B4 measured the two closures 0.2%/0.8% apart on this family)."""
    mc, wc, _, _, r_probe = naca_case
    r = solve_laplace_lifting(mc, wc, alpha_deg=ALPHA,
                              kutta_estimator="pressure")
    assert r["kutta_converged"]
    assert r["n_kutta_updates"] < 20                 # the G2.3 budget
    rel = abs(r["gamma"][0] - r_probe["gamma"][0]) / abs(r_probe["gamma"][0])
    assert rel < 0.02, f"pressure vs probe Gamma differs {100 * rel:.2f}%"


def test_implied_target_fixed_point(naca_case):
    """At the converged pressure solve, the implied target's fixed point IS
    the nonlinear pressure equality: Gamma* == Gamma and F ~ 0."""
    mc, wc, cvs, con, _ = naca_case
    r = solve_laplace_lifting(mc, wc, alpha_deg=ALPHA,
                              kutta_estimator="pressure")
    phi_cut = con.expand(r["phi"][: wc.n_nodes_orig], r["gamma"])
    g_star = cvs.implied_targets(phi_cut, r["gamma"])
    assert np.abs(g_star - r["gamma"]).max() < 1e-7
    assert np.abs(cvs.residual_stations(phi_cut)).max() < 1e-5


def test_probe_default_bit_identical(naca_case):
    """G14.4 at the Laplace driver: the default and the explicit "probe"
    request produce byte-identical Gamma and phi."""
    mc, wc, _, _, r_default = naca_case
    r_probe = solve_laplace_lifting(mc, wc, alpha_deg=ALPHA,
                                    kutta_estimator="probe")
    assert np.array_equal(r_default["gamma"], r_probe["gamma"])
    assert np.array_equal(r_default["phi"], r_probe["phi"])


def test_unknown_estimator_raises(naca_case):
    mc, wc, _, _, _ = naca_case
    with pytest.raises(ValueError, match="kutta_estimator"):
        solve_laplace_lifting(mc, wc, alpha_deg=ALPHA,
                              kutta_estimator="typo")


# ---------------------------------------------------- stage 2: Newton path


@pytest.fixture(scope="module")
def newton_pressure_case(naca_case):
    mc, wc, _, _, _ = naca_case
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             kutta_estimator="pressure", precond="direct")
    return mc, wc, r


def test_newton_pressure_converges_matches_probe(newton_pressure_case):
    """Coupled Newton with the pressure estimator converges at M0.5 and
    lands within 1% of the probe path's Gamma (G14.3 NACA leg)."""
    mc, wc, r = newton_pressure_case
    assert r["converged"]
    assert r["kutta_estimator"] == "pressure"
    assert r["kutta_sigma"] is not None
    r_probe = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                   precond="direct")
    assert r_probe["kutta_sigma"] is None
    rel = abs(r["gamma"][0] - r_probe["gamma"][0]) / abs(r_probe["gamma"][0])
    assert rel < 0.01, f"pressure vs probe Gamma differs {100 * rel:.2f}%"


def test_newton_pressure_kutta_row_fd(newton_pressure_case):
    """sigma-scaled K_p rows vs central FD of F w.r.t. phi_free through
    eval_residual: F is exactly quadratic, so this is roundoff-exact (the
    pressure analog of test_kutta_row_exact)."""
    _, _, r = newton_pressure_case
    ws = r["workspace"]
    phi_free = np.asarray(r["phi"])[:ws.n_red][ws.free].copy()
    gamma = np.asarray(r["gamma"], dtype=np.float64).copy()
    _, _, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                   M_CAP, RHO_FLOOR)
    Kp_cut, _ = ws.cvs.newton_rows(state["phi_cut"])
    Kp_free = (Kp_cut @ ws.con.T).tocsr()[:, ws.free]

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
        exact = ws.kutta_sigma * (Kp_free @ delta)
        scale = max(np.abs(fd).max(), 1e-30)
        assert np.abs(exact - fd).max() / scale < 1e-8


def test_newton_pressure_gamma_jacobian_fd(newton_pressure_case):
    """Central FD of F w.r.t. Gamma THROUGH eval_residual (which also moves
    the far-field Dirichlet values by V_red) vs sigma * D from newton_rows:
    equality end-to-end proves the CV Dirichlet-contact assert's promise
    that dF/dGamma carries only the slave-jump chain (the pressure analog
    of test_gamma_column_fd's dF/dGamma = -I clause)."""
    _, _, r = newton_pressure_case
    ws = r["workspace"]
    phi_free = np.asarray(r["phi"])[:ws.n_red][ws.free].copy()
    gamma = np.asarray(r["gamma"], dtype=np.float64).copy()
    _, _, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                   M_CAP, RHO_FLOOR)
    _, D = ws.cvs.newton_rows(state["phi_cut"])
    D_s = ws.kutta_sigma[:, None] * D

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
        assert np.abs(D_s[:, j] - fd).max() / scale < 1e-7


def test_newton_pressure_sigma_independence(newton_pressure_case):
    """sigma cancels in the eliminated blocks (-(sigma D)^{-1}(sigma F)):
    scaling the frozen sigma leaves K_tilde @ x and F_tilde unchanged, so
    the Newton iterates are sigma-independent."""
    _, _, r = newton_pressure_case
    ws = r["workspace"]
    phi_free = np.asarray(r["phi"])[:ws.n_red][ws.free].copy()
    gamma = np.asarray(r["gamma"], dtype=np.float64).copy()
    _, F, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                   M_CAP, RHO_FLOOR)
    K1, F1 = ws.kutta_blocks(state, F)
    sigma_old = ws.kutta_sigma
    try:
        ws.kutta_sigma = 7.0 * sigma_old
        K2, F2 = ws.kutta_blocks(state, 7.0 * F)
    finally:
        ws.kutta_sigma = sigma_old
    rng = np.random.default_rng(11)
    x = rng.standard_normal(ws.n_free)
    np.testing.assert_allclose(K1 @ x, K2 @ x, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(F1, F2, rtol=1e-12, atol=1e-16)


def test_newton_m0_matches_laplace_pressure(naca_case):
    """m_inf = 0 pressure Newton (cold start) and the pressure Laplace
    driver converge to the same circulation -- two independent closures of
    the same nonlinear pressure equality on the linear PDE."""
    mc, wc, _, _, _ = naca_case
    r_n = solve_newton_lifting(mc, wc, m_inf=0.0, alpha_deg=ALPHA,
                               kutta_estimator="pressure", n_picard_seed=0,
                               precond="direct")
    r_l = solve_laplace_lifting(mc, wc, alpha_deg=ALPHA,
                                kutta_estimator="pressure")
    assert r_n["converged"]
    assert abs(r_n["gamma"][0] - r_l["gamma"][0]) < 1e-6


def test_newton_workspace_estimator_mismatch_raises(newton_pressure_case):
    mc, wc, r = newton_pressure_case
    with pytest.raises(ValueError, match="kutta_estimator"):
        solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             workspace=r["workspace"])   # default = probe


def test_pressure_tip_taper_not_implemented(naca_case):
    mc, wc, _, _, _ = naca_case
    with pytest.raises(NotImplementedError, match="tip_taper"):
        NewtonWorkspace(mc, wc, alpha_deg=ALPHA,
                        tip_taper=np.full(wc.n_stations, 0.5),
                        kutta_estimator="pressure")
