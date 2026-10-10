"""A52 —— 共形 Newton 驱动的模块性质（A 类）。

★ 2026-08-24 从 `test_p8_newton.py` 拆出。断言的是 **FD 精确性、Kutta 行、GMRES+AMG、
Term 3 的非对称结构、Newton 二次收敛阶、与 Picard 路径的一致性** —— **不需要任何真值**。
★ 已提交物理值的那三条（G8.1 coarse/medium、G8.2 M6）在 `tests/E/test_E01_p8_newton_anchors.py`。
★ 实测：两半**共享模块级名字 0 个** ⇒ 干净分割，不需要共享模块。
"""
import os
import time
from pathlib import Path
import numpy as np
import pytest
import scipy.sparse.linalg as spla
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting
from pyfp3d.solve.picard import solve_subsonic_lifting
from bench.recipes import NEWTON_M6_RECIPE  # noqa: E402  (裁决 a, 2026-08-24)


M_INF = 0.5


ALPHA = 2.0


UPWIND_C = 1.5


M_CRIT = 0.95


M_CAP = 3.0


RHO_FLOOR = 0.05


def _case_args():
    return dict(upwind_c=UPWIND_C, m_crit=M_CRIT, m_cap=M_CAP,
                rho_floor=RHO_FLOOR)


@pytest.fixture(scope="module")
def coarse_mesh():
    from tests.conftest import REPO_ROOT

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

    from tests.mesh_utils import generate_structured_cube_mesh
    from tests.A.test_A21_newton_jacobian_fd import _degeneracy_breaker

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
    # per-step Gamma instrumentation (capability demo): 1:1 with the
    # residual history, ends at the returned circulation
    assert len(r["gamma_history"]) == len(h)
    assert np.allclose(r["gamma_history"][-1], r["gamma"], atol=1e-12)


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


class TestPreconditionersAreAccelerationNotPhysics:
    r"""★★★ W3.7 / H20（2026-09-06）：**三个预条件子必须给同一个答案。**

    审计实测:`precond="ilu"` 在 `newton.py:1346` 实装,而 **`tests/` 里零引用**
    —— 一条上线的代码路径,没有任何东西跑过它。
    ⇒ 处置本来是二选一(补冒烟 / 声明弃用)。**先量再定**:实测它跑得通,
    且与 `direct` 一致到 **1.63e-13** ⇒ 不该弃用,该补锁。

    判据形状不是「ilu 能跑」,而是 **预条件子是加速手段,不该改变答案**:
    三条腿的 cl 必须一致到远严于任何物理判据的带。
    ★ 带取 `1e-9`(实测最坏 1.6e-13,余量 ~6000x)—— 它**不是**机器零:
    三者走不同的 Krylov 迭代序列,终止在 `tol_residual` 附近的不同点上,
    所以位相同在原理上就不该期待。
    ★ G-DOMAIN:相反结果 = 某个预条件子给出不同的 cl ⇒ 那是**求解器缺陷**
    (预条件子改变了解),不是容差问题,**不要放宽这条带**。
    """

    #: 实测 2026-09-06(coarse, M0.50/α2.0):
    #:   direct cl=0.277589780302  amg cl=0.277589780299  ilu cl=0.277589780302
    PRECOND_CL_TOL = 1e-9

    def test_direct_amg_ilu_agree(self, coarse_mesh):
        mc, wc = coarse_mesh
        dz = float(np.ptp(mc.nodes[:, 2]))
        cl = {}
        for pc in ("direct", "amg", "ilu"):
            r = solve_newton_lifting(mc, wc, m_inf=0.50, alpha_deg=2.0,
                                     upwind_c=1.5, m_crit=0.95, precond=pc,
                                     n_picard_seed=5, n_newton_max=400)
            assert r["converged"], f"precond={pc} 不收敛 —— 前提变了"
            f = wall_force_coefficients(mc.nodes, mc.elements,
                                        mc.boundary_faces["wall"],
                                        np.asarray(r["phi"]), alpha_deg=2.0,
                                        u_inf=1.0, s_ref=dz, m_inf=0.50)
            cl[pc] = float(f["cl"])
        base = cl["direct"]
        for pc in ("amg", "ilu"):
            rel = abs(cl[pc] - base) / abs(base)
            assert rel <= self.PRECOND_CL_TOL, (
                f"precond={pc} 的 cl {cl[pc]:.12f} 与 direct {base:.12f} "
                f"差 {rel:.2e} > {self.PRECOND_CL_TOL:.0e} —— **预条件子改变了解**,"
                "这是求解器缺陷,不是容差问题")

    def test_the_ilu_path_is_actually_reached(self, coarse_mesh):
        """前提断言:`precond="ilu"` 真的走到了 ILU 分支,不是被静默回退。

        ★ 少了它,上一条可能只是「三次都跑了 direct」⇒ 按构造恒真。

        ★★ **补丁必须打在消费者模块上,不是定义模块上**(2026-09-06 实测踩到)。
        `newton.py:71` 写的是 `from pyfp3d.solve.linear import
        build_ilu_preconditioner`,名字在 **import 时**就绑进了 `newton` 的命名空间;
        去 patch `pyfp3d.solve.linear` 的属性**打不到它**,于是这条前提断言会
        报一个**假红**(我第一版就是这么写的,被它自己红出来的)。
        """
        import pyfp3d.solve.newton as N
        seen = []
        orig = N.build_ilu_preconditioner
        N.build_ilu_preconditioner = lambda *a, **k: (seen.append(1), orig(*a, **k))[1]
        try:
            mc, wc = coarse_mesh
            solve_newton_lifting(mc, wc, m_inf=0.50, alpha_deg=2.0, upwind_c=1.5,
                                 m_crit=0.95, precond="ilu", n_picard_seed=5,
                                 n_newton_max=400)
        finally:
            N.build_ilu_preconditioner = orig
        assert seen, "precond='ilu' 没有调用 build_ilu_preconditioner —— 被静默回退了"

