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
from pyfp3d.solve.picard import solve_laplace_lifting

ALPHA = 2.0


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
