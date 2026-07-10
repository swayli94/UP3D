"""P8/N3+N4: the fully-coupled (phi_red, Gamma) Newton driver
(solve/newton.py, design.md Sec 8.1) -- subsonic milestone tests.

Covers, on the NACA0012 coarse 2.5D case at M0.5 / alpha 2:

  - the Gamma-Jacobian column B = J_red[free, dir] @ V_red + H_J[free, :]
    against a central difference of the coupled residual (subsonic, so
    the walk selection is inert and the raw residual is FD-safe). This is
    the far-field-vortex-column trap detector (roadmap P8: the column is
    folded silently into the Picard RHS and is easy to omit from the
    Newton Gamma block).
  - the Kutta row K = dF/dphi_free (affine -> FD exact to roundoff).
  - the far-field linearity assumption vals(Gamma) = vals0 + V_red Gamma
    that set_mach's unit-Gamma probing relies on.
  - GMRES+AMG on the assembled Newton J_ff vs a direct sparse solve, and
    the supersonic nonsymmetry of the exact Jacobian (why GMRES, not CG).
  - N4 acceptance: coupled Newton matches the P3 Picard solution
    (|dcl/cl| < 5e-3, ||dGamma||_inf < 1e-6), converges in a handful of
    steps to ||R||_inf < 1e-10, with terminal-quadratic order.

Gated (PYFP3D_TRANSONIC_GATES=1): the G8.1 terminal-quadratic runs on the
RE-SPECCED case set (user-approved 2026-07-11): coarse M0.80/alpha 1.25
and medium M0.7875/alpha 1.25. The original "G4.1 case" (M0.80 medium)
sits at the edge of the FP non-uniqueness fold on the medium mesh --
measured dcl/dM ~ 10 between M0.775 (cl 0.396) and M0.7875 (cl 0.523),
no reachable isolated solution at M0.80 -- so it is recorded as a
model-validity finding (design.md Sec 12 risks 2/3), not a solver target.
The physics bands asserted here are REGRESSION LOCKS around the measured
Newton solutions (the Euler-anchored G4.1 band does not bind conservative
FP at this shock strength; the P4 Picard states are stall artifacts with
Newton residual 2.2e-4 -- see the roadmap P4 erratum).
"""

import os

import numpy as np
import pytest
import scipy.sparse.linalg as spla

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting
from pyfp3d.solve.picard import solve_subsonic_lifting

M_INF = 0.5
ALPHA = 2.0
UPWIND_C = 1.5
M_CRIT = 0.95
M_CAP = 3.0
RHO_FLOOR = 0.05

run_gates = pytest.mark.skipif(
    os.environ.get("PYFP3D_TRANSONIC_GATES", "0") != "1",
    reason="transonic Newton gate runs take minutes; "
           "set PYFP3D_TRANSONIC_GATES=1 for the gate-closure run",
)

#: the N5 transonic Newton recipe (measured on the 2.5D family):
#: direct exact steps (the shock-position soft mode leaves Krylov-eta
#: steps stalled), stall-adaptive freeze with active-set refresh, fine
#: Mach steps near the top of the ramp.
NEWTON_TRANSONIC_RECIPE = dict(
    dm=0.025, dm_min=0.003, freeze_tol=1e-6,
    newton_kw=dict(freeze_refresh_max=8, precond="direct",
                   n_newton_max=60),
)


def _case_args():
    return dict(upwind_c=UPWIND_C, m_crit=M_CRIT, m_cap=M_CAP,
                rho_floor=RHO_FLOOR)


@pytest.fixture(scope="module")
def coarse_mesh():
    from .conftest import REPO_ROOT

    mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
                     / "coarse.msh")
    return cut_wake(mesh)


@pytest.fixture(scope="module")
def newton_case(coarse_mesh):
    mc, wc = coarse_mesh
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             **_case_args())
    assert r["converged"]
    return mc, wc, r


@pytest.fixture(scope="module")
def picard_case(coarse_mesh):
    mc, wc = coarse_mesh
    r = solve_subsonic_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA)
    assert r["converged"] and r["kutta_converged"]
    return r


def _cl(mc, phi):
    dz = float(np.ptp(mc.nodes[:, 2]))
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
        alpha_deg=ALPHA, s_ref=1.0 * dz, m_inf=M_INF,
    )
    return float(forces["cl"])


# ------------------------------------------------- Gamma-Jacobian blocks


def test_farfield_values_linear_in_gamma(coarse_mesh):
    """vals_red(Gamma) = vals0_red + V_red @ Gamma to machine precision --
    the linearity set_mach's unit-Gamma probing relies on (guards against
    a later nonlinear far-field change silently breaking the Newton
    Gamma column)."""
    mc, wc = coarse_mesh
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA)
    ws.set_mach(M_INF)
    rng = np.random.default_rng(3)
    from pyfp3d.constraints.dirichlet import farfield_dirichlet

    for _ in range(3):
        g = rng.standard_normal(ws.n_st)
        _, vals = farfield_dirichlet(
            mc, wc, ALPHA, g, 1.0, ws.vortex_center, beta=ws.beta,
            spanwise_gamma=ws.spanwise,
        )
        vals_red = ws._reduce_ff_values(vals)
        affine = ws.vals0_red + ws.V_red @ g
        assert np.max(np.abs(vals_red - affine)) < 1e-13


def test_gamma_column_fd(newton_case):
    """B[:, j] vs central FD of the coupled residual w.r.t. Gamma_j at the
    converged state. An omitted far-field vortex column (or a stale
    Picard-level h_j in place of the exact T^T J g_j) fails this at O(1)."""
    mc, wc, r = newton_case
    ws = r["workspace"]
    phi_free = np.asarray(r["phi"])[:ws.n_red][ws.free].copy()
    gamma = np.asarray(r["gamma"], dtype=np.float64).copy()

    _, _, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                   M_CAP, RHO_FLOOR)
    _, B = ws.assemble_coupled(state, UPWIND_C, M_CRIT, RHO_FLOOR)

    eps = 1e-5
    for j in range(ws.n_st):
        dg = np.zeros(ws.n_st)
        dg[j] = eps
        R_p, F_p, _ = ws.eval_residual(phi_free, gamma + dg, UPWIND_C,
                                       M_CRIT, M_CAP, RHO_FLOOR)
        R_m, F_m, _ = ws.eval_residual(phi_free, gamma - dg, UPWIND_C,
                                       M_CRIT, M_CAP, RHO_FLOOR)
        fd = (R_p - R_m) / (2.0 * eps)
        col = np.asarray(B[:, j].todense()).ravel()
        scale = np.abs(fd).max()
        assert scale > 0.0
        rel = np.abs(col - fd).max() / scale
        assert rel < 1e-6, f"Gamma column {j}: rel err {rel:.3e}"
        # dF/dGamma = -I exactly (kutta_targets does not read Gamma)
        fd_F = (F_p - F_m) / (2.0 * eps)
        expected = np.zeros(ws.n_st)
        expected[j] = -1.0
        assert np.max(np.abs(fd_F - expected)) < 1e-9


def test_kutta_row_exact(newton_case):
    """K @ delta vs FD of F w.r.t. phi_free: kutta_targets is affine, so
    the central difference is exact to roundoff -- including the
    shared-probe rows (adjacent stations reusing a TE probe node)."""
    mc, wc, r = newton_case
    ws = r["workspace"]
    phi_free = np.asarray(r["phi"])[:ws.n_red][ws.free].copy()
    gamma = np.asarray(r["gamma"], dtype=np.float64).copy()

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
        assert np.max(np.abs(ws.K @ delta - fd)) < 1e-10


# ------------------------------------------------------- linear solve (N3)


def test_gmres_amg_solves_newton_jacobian(newton_case):
    """GMRES preconditioned by AMG on the SPD Picard block solves the
    assembled Newton J_ff to direct-solve accuracy."""
    from pyfp3d.solve.linear import build_amg_preconditioner, solve_gmres

    mc, wc, r = newton_case
    ws = r["workspace"]
    phi_free = np.asarray(r["phi"])[:ws.n_red][ws.free].copy()
    gamma = np.asarray(r["gamma"], dtype=np.float64).copy()
    _, _, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                   M_CAP, RHO_FLOOR)
    J_ff, _ = ws.assemble_coupled(state, UPWIND_C, M_CRIT, RHO_FLOOR)

    A_pic = ws.op.assemble_matrix(state["rho_t"])
    A_ff = (ws.con.T.T @ (A_pic @ ws.con.T)).tocsr()[ws.free][:, ws.free]
    _, M_pre = build_amg_preconditioner(A_ff.tocsr())

    rng = np.random.default_rng(11)
    b = rng.standard_normal(ws.n_free)
    x_ref = spla.spsolve(J_ff.tocsc(), b)
    x, n_it, _ = solve_gmres(J_ff, b, M=M_pre, rtol=1e-12)
    rel = np.max(np.abs(x - x_ref)) / np.max(np.abs(x_ref))
    assert rel < 1e-8, f"GMRES vs direct rel err {rel:.3e} ({n_it} iters)"


def test_supersonic_jacobian_is_nonsymmetric():
    """The Term-3 upstream coupling has no transpose partner: the exact
    Jacobian is nonsymmetric wherever the pocket is active (why the Newton
    linear solve is GMRES, not CG -- design.md Sec 6.3)."""
    from pyfp3d.kernels.jacobian import PicardOperator
    from pyfp3d.kernels.upwind import UpwindOperator
    from pyfp3d.physics.isentropic import density_field

    from .mesh_utils import generate_structured_cube_mesh
    from .test_p8_jacobian import _degeneracy_breaker

    nodes, elements = generate_structured_cube_mesh(n=6, L=1.0)
    op = PicardOperator(nodes, elements)
    upw = UpwindOperator(nodes, elements, weighted=False)
    x = nodes[:, 0]
    phi = x + 0.3 * x ** 2 + _degeneracy_breaker(len(nodes))
    grad, q2 = op.velocities(phi)
    grad, q2 = grad.copy(), q2.copy()
    rho = density_field(q2, 0.8)
    rho_t = upw.rho_tilde(grad, q2, rho, 0.8, UPWIND_C, M_CRIT).copy()
    s_e, s_u, upstream = upw.rho_tilde_sensitivities(
        grad, q2, rho, 0.8, UPWIND_C, M_CRIT)
    J = op.assemble_newton_jacobian(phi, rho_t, s_e.copy(), s_u.copy(),
                                    upstream.copy())
    asym = J - J.T
    asym.eliminate_zeros()
    assert op.n_term3_active > 0
    assert asym.nnz > 0
    assert np.max(np.abs(asym.data)) > 1e-8


# --------------------------------------------------- N4 acceptance ladder


def test_newton_subsonic_matches_p3(newton_case, picard_case):
    """The coupled Newton lands on the SAME discrete solution as the P3
    nested Picard/secant: identical discretization, so the agreement is
    tight (gate wording < 0.5% on cl; measured ~1e-7)."""
    mc, wc, r_newton = newton_case
    r_picard = picard_case

    assert r_newton["n_newton"] <= 10
    assert r_newton["residual_history"][-1] < 1e-10
    assert r_newton["n_limited"] == 0 and r_newton["n_floored"] == 0
    assert np.max(np.abs(np.asarray(r_newton["gamma"])
                         - np.asarray(r_picard["gamma"]))) < 1e-6

    cl_n = _cl(mc, r_newton["phi"])
    cl_p = _cl(mc, r_picard["phi"])
    assert abs(cl_n / cl_p - 1.0) < 5e-3, (
        f"cl mismatch: newton {cl_n:.6f} vs picard {cl_p:.6f}")


def test_newton_subsonic_terminal_order(coarse_mesh):
    """Terminal quadratic convergence from a freestream cold start (no
    Picard seed): observed order p_k reaches ~2 and the final step is a
    super-linear residual collapse. The Eisenstat-Walker forcing makes
    the FIRST steps inexact by design, so the assertion is on the
    terminal behaviour (G8.1's protocol, measured here subsonic)."""
    mc, wc = coarse_mesh
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             n_picard_seed=0, **_case_args())
    assert r["converged"]
    assert r["n_newton"] <= 8
    h = r["residual_history"]
    assert h[-1] < 1e-10
    # terminal super-linear collapse: last step gains >= 3 digits
    assert h[-1] / h[-2] < 1e-3
    orders = r["newton_orders"]
    assert len(orders) >= 1
    assert max(orders) > 1.8, f"observed orders {orders}"
    assert orders[-1] > 1.5, f"terminal order {orders[-1]:.2f}"


# ------------------------------------------------- gated: G8.1 transonic


def _assert_terminal_quadratic(r):
    """G8.1 protocol: the final level must end below 1e-9 AND contain a
    consecutive PAIR of >= 1.5-digit residual collapses (only a quadratic
    tail does that at these levels; cf. Lopez Table 4.9). The pair is
    searched over the level history rather than read off the last two
    entries because the freeze-refresh honesty re-evaluations interleave
    live-residual jumps with the (quadratic) frozen phases -- a final
    frozen phase can legitimately converge in ONE step after a refresh."""
    h = r["residual_history"]
    assert h[-1] < 1e-9, f"final residual {h[-1]:.3e}"
    drops = [h[i + 1] / h[i] for i in range(len(h) - 1)]
    pairs = [(drops[i], drops[i + 1]) for i in range(len(drops) - 1)]
    assert any(a < 3e-2 and b < 3e-2 for a, b in pairs), (
        f"no consecutive quadratic-like drop pair in {drops}")


def _transonic_case(mesh_file, m_inf):
    from .conftest import REPO_ROOT
    from pyfp3d.post.section_cut import wall_cp_curve
    from pyfp3d.post.shock import shock_report
    from pyfp3d.solve.newton import solve_newton_transonic

    mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
                     / mesh_file)
    mc, wc = cut_wake(mesh)
    r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=1.25,
                               **NEWTON_TRANSONIC_RECIPE)
    dz = float(np.ptp(mc.nodes[:, 2]))
    rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=m_inf),
                       m_inf)
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
        alpha_deg=1.25, s_ref=dz, m_inf=m_inf)
    return r, rep, forces


@run_gates
def test_g81_terminal_quadratic_coarse_m080():
    """G8.1 (re-specced): coarse M0.80/alpha1.25 -- terminal quadratic
    convergence to the TRUE discrete solution (regression lock: shock
    0.658, cl 0.459, M_max 1.408 -- the Newton answer, dissipation-scan
    robust and continuation-path independent; NOT the P4 Picard stall
    state 0.604/0.334, whose Newton residual is 2.2e-4)."""
    r, rep, forces = _transonic_case("coarse.msh", 0.80)
    assert r["converged"]
    _assert_terminal_quadratic(r)
    assert r["n_limited"] == 0 and r["n_floored"] == 0
    assert r["F_history"][-1] < 1e-12                   # coupled Kutta
    assert rep["upper"]["has_shock"] and rep["upper"]["monotone"]
    assert abs(rep["upper"]["x_shock"] - 0.658) < 0.012
    assert abs(forces["cl"] - 0.459) < 0.01
    assert abs(float(np.sqrt(r["mach2_max"])) - 1.408) < 0.02


@run_gates
def test_g81_terminal_quadratic_medium_m07875():
    """G8.1 (re-specced): medium M0.7875/alpha1.25 -- the strongest
    medium-mesh condition with a reachable isolated solution (M0.80 sits
    at the fold, recorded finding). Regression lock around the measured
    Newton solution: shock 0.674, cl 0.523, M_max 1.404; the frozen
    finish reaches ~8e-11 with the live assignment-discontinuity floor
    reported (~1.3e-7 -- the C0 walk flux's intrinsic floor on this
    tie-degenerate prism mesh)."""
    r, rep, forces = _transonic_case("medium.msh", 0.7875)
    assert r["converged"]
    _assert_terminal_quadratic(r)
    assert r["n_limited"] == 0 and r["n_floored"] == 0
    assert rep["upper"]["has_shock"] and rep["upper"]["monotone"]
    assert abs(rep["upper"]["x_shock"] - 0.674) < 0.012
    assert abs(forces["cl"] - 0.523) < 0.01
    assert abs(float(np.sqrt(r["mach2_max"])) - 1.404) < 0.02
    if r["residual_unfrozen"] is not None:
        assert r["residual_unfrozen"] < 1e-5


def test_newton_incompressible_single_step(coarse_mesh):
    """m_inf = 0: the problem is linear (rho == 1), so Newton from
    freestream converges in ONE step to the P2 lifting Laplace solution --
    provided the inner solve is exact (the Eisenstat-Walker default
    eta_0 = 1e-2 would trade this for a few cheap inexact steps)."""
    from pyfp3d.solve.picard import solve_laplace_lifting

    mc, wc = coarse_mesh
    r = solve_newton_lifting(mc, wc, m_inf=0.0, alpha_deg=ALPHA,
                             n_picard_seed=0, ew_eta0=1e-10,
                             ew_eta_max=1e-10, **_case_args())
    assert r["converged"]
    assert r["n_newton"] <= 2                 # 1 linear step (+ roundoff)
    r_p2 = solve_laplace_lifting(mc, wc, alpha_deg=ALPHA)
    assert np.max(np.abs(np.asarray(r["gamma"])
                         - np.asarray(r_p2["gamma"]))) < 1e-7
