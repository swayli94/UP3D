"""A55 -- G10.1:生产 Newton 驱动器的**非升力入口**(A 类:机制/一致性,无外部真值)。

约束文本(判据事前写死):docs/dev_phase_six/20261011-1100-w22a-prereg.md (W2.2a),
判据原文 = phases/p1/docs/roadmap/track_p.md:1066 (G10.1)。

`solve_newton_nonlifting` 不是第二个驱动器:它以 `wc=None` 调 `solve_newton_lifting`
—— 同一个 NewtonWorkspace(泛化:恒等约化、n_st=0、无 Kutta 行、无 Γ 未知数、
远场无涡基)、同一个 Newton 循环。因此 m_cap/rho_floor 与钳制记账
(n_limited/n_floored/clamped + 钳死不许报收敛)、熵修正(默认 ON)、accept_reason、
residual_history、newton_orders、timings 全部在位 —— 这正是它相对重写驱动器
(`tests/_nozzle_case.py::DuctSystem`)的意义。

五条判据(预注册 §3,断言消息里写清「红了意味着什么」):
  C1  sphere_shell medium M0.3:新入口 vs solve_subsonic(同一离散),
      |Δφ|∞ < 1e-8 且 peak-Cp 差 < 1e-6。
  C2  m_inf=0:一步收敛,且与 solve_laplace 机器零一致。
  C3  终端二次收敛阶(newton_orders 尾部 ≈ 2)。
  C4  位相同守卫(A27 模板):lifting 算例重复求解 array_equal。
  C5  熵修正在新入口可用且默认 ON:熵 ON 的 sphere 腿收敛,且亚声速下
      sigma ≡ 1(与 OFF 位相同 —— 惰性是逐位的,不是近似的)。

W2.2b(20261011-1300 预注册)扩两条机制腿:处方 Dirichlet 模式
(dirichlet_nodes/dirichlet_values 取代 farfield 默认)的销钉读回 +
「自由流等价值 ⇔ farfield 默认」逐位一致(§4 第二条 kill 的适定形式),
以及处方模式与 lifting 形状/复用 workspace 的互斥报错。

实测(2026-10-11,16 线程,sphere medium 16793 节点):
  C1 |Δφ|∞ = 6.0e-12(带 1e-8),peak-Cp 差 7.3e-11(带 1e-6);
  C2 n_newton = 1,|Δφ| vs solve_laplace(rtol=1e-13) = 5.8e-12;
  C3 newton_orders = [2.86, 1.87](末值被装配舍入地板压低,典型形态);
  C5 n_shock = 0、sigma_min = 1.0、ON/OFF array_equal(phi) = True。
"""
import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.physics.isentropic import pressure_coefficient
from pyfp3d.post.surface import wall_tangential_gradient_quadratic
from pyfp3d.solve.newton import (
    NewtonWorkspace,
    solve_newton_lifting,
    solve_newton_nonlifting,
)
from pyfp3d.solve.picard import solve_laplace, solve_subsonic

M_INF = 0.3
#: C1/C3/C5 共享腿用默认 precond("amg" —— 生产默认,W3.7 教训:默认路径必须有锁);
#: 只有 C2 的「一步 ≡ Laplace」需要精确内解(precond="direct",cf. A52 的同注:
#: EW 默认 forcing 会用几个便宜的非精确步换掉一步到位)。
TOL_SHARP = 1e-12   # 双方都比生产默认收得更紧,让 C1 的比值逼近机器零


@pytest.fixture(scope="module")
def sphere():
    from tests.conftest import REPO_ROOT

    mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "sphere_shell"
                     / "medium.msh")
    ff = np.unique(mesh.boundary_faces["farfield"])
    # alpha = 0, u_inf = 1 的自由流 Dirichlet —— 与新入口内部
    # _farfield_dirichlet_nonlifting 的 freestream_phi 逐位一致。
    phi_ff = mesh.nodes[ff, 0]
    return mesh, ff, phi_ff


@pytest.fixture(scope="module")
def newton_m03(sphere):
    """C1/C3/C5 共享的新入口腿:默认(amg)预条件子,熵修正默认 ON。"""
    mesh, _, _ = sphere
    r = solve_newton_nonlifting(mesh, m_inf=M_INF, tol_residual=TOL_SHARP)
    assert r["converged"], "前提失败:新入口在 sphere M0.3 不收敛"
    return r


@pytest.fixture(scope="module")
def picard_m03(sphere):
    mesh, ff, phi_ff = sphere
    r = solve_subsonic(mesh.nodes, mesh.elements, ff, phi_ff, m_inf=M_INF,
                       phi_init=mesh.nodes[:, 0].copy(), rtol=1e-13,
                       tol_rho=1e-13, n_picard_max=200)
    assert r["converged"], "前提失败:solve_subsonic 在 sphere M0.3 不收敛"
    return r


def _wall_cp(mesh, phi):
    wall = mesh.boundary_faces["wall"]
    wn = np.unique(wall)
    g = wall_tangential_gradient_quadratic(mesh.nodes, wall,
                                           np.asarray(phi, float))
    q2 = np.sum(g[wn] ** 2, axis=1)
    return np.array([pressure_coefficient(q, M_INF) for q in q2])


def test_c1_sphere_m03_matches_solve_subsonic(sphere, newton_m03, picard_m03):
    mesh, _, _ = sphere
    r_new, r_pic = newton_m03, picard_m03
    # 钳制记账在位且必须干净(不干净时驱动器本身拒绝报收敛 —— GS1.4 语义原样继承)
    assert r_new["n_limited"] == 0 and r_new["n_floored"] == 0
    assert not r_new["clamped"]
    # 无意义键的空/惰形态(docstring 写明的那部分契约)
    assert np.asarray(r_new["gamma"]).shape == (0,)
    assert all(f == 0.0 for f in r_new["F_history"])
    assert r_new["kutta_sigma"] is None and r_new["kutta_weld_sign"] is None

    dphi = float(np.max(np.abs(np.asarray(r_new["phi"])
                               - np.asarray(r_pic["phi"]))))
    assert dphi < 1e-8, (
        f"|Δφ|∞ = {dphi:.3e} > 1e-8 —— 同一离散的两条求解路径(生产 Newton vs "
        "Picard solve_subsonic)真分歧。按预注册 §4 这是 kill criterion:先归因"
        "(装配/远场/边界),不许放宽带。实测基线 6.0e-12。")

    cp_new, cp_pic = _wall_cp(mesh, r_new["phi"]), _wall_cp(mesh, r_pic["phi"])
    i_peak = int(np.argmin(cp_pic))          # 吸力峰 = 球的赤道点
    d_peak = float(abs(cp_new[i_peak] - cp_pic[i_peak]))
    d_all = float(np.max(np.abs(cp_new - cp_pic)))
    assert d_peak < 1e-6, (
        f"peak-Cp 差 = {d_peak:.3e} > 1e-6(全场 max|ΔCp| = {d_all:.3e})—— "
        "壁面 Cp 恢复对新旧路径分歧敏感,红了先归因同上。实测基线 7.3e-11。")


def test_c1_precond_is_acceleration_not_physics(sphere, newton_m03):
    """附属锁(W3.7 形状):direct 腿必须与默认 amg 腿落在同一离散解上。"""
    mesh, _, _ = sphere
    r_dir = solve_newton_nonlifting(mesh, m_inf=M_INF, precond="direct",
                                    tol_residual=TOL_SHARP)
    assert r_dir["converged"]
    dphi = float(np.max(np.abs(np.asarray(r_dir["phi"])
                               - np.asarray(newton_m03["phi"]))))
    assert dphi < 1e-8, (
        f"direct vs amg |Δφ|∞ = {dphi:.3e} —— 预条件子改变了解,这是求解器缺陷,"
        "不是容差问题。实测基线 2.6e-13。")


def test_c2_m0_one_step_equals_solve_laplace(sphere):
    """m_inf=0:rho ≡ 1(逐位,isentropic.py 的 G3.3 保证),问题线性 ⇒
    精确内解(precond='direct')下 Newton 必须**一步**收敛到 P1 Laplace 离散解。"""
    mesh, ff, phi_ff = sphere
    r = solve_newton_nonlifting(mesh, m_inf=0.0, precond="direct")
    assert r["converged"], (
        "m_inf=0 的线性问题没收敛 —— 初值/线性化链路断了(预注册 §4 kill)")
    assert r["n_newton"] == 1, (
        f"n_newton = {r['n_newton']} != 1 —— 线性问题多步才收敛 = 精确内解或"
        "残差/Jacobian 一致性被破坏(预注册 §4 kill)。residual_history = "
        f"{r['residual_history']}")
    r_lap = solve_laplace(mesh.nodes, mesh.elements, ff, phi_ff,
                          rtol=1e-13, maxiter=10000)
    dphi = float(np.max(np.abs(np.asarray(r["phi"])
                               - np.asarray(r_lap["phi"]))))
    assert dphi < 1e-9, (
        f"|Δφ| vs solve_laplace = {dphi:.3e} —— 超出机器零:新入口的 m=0 系统"
        "不再是 P1 Laplace 离散(LU 精确步 vs CG@1e-13 的差只能在这个量级)。"
        "实测基线 5.8e-12。")


def test_c3_terminal_quadratic_order(newton_m03):
    """终端二次阶(G10.1 原文;阈值取 A52 同型:max > 1.8,末值 > 1.5 —
    末值被装配舍入地板压低是典型形态,不是降阶)。"""
    orders = newton_m03["newton_orders"]
    assert len(orders) >= 1, "newton_orders 为空 —— 收敛太快或历史没记上"
    assert max(orders) > 1.8 and orders[-1] > 1.5, (
        f"newton_orders = {orders} —— 终端二次阶消失 = Jacobian 与残差不再一致"
        "(精确性坏了)或线性求解太松。实测基线 [2.86, 1.87]。")


def test_c4_lifting_path_bit_identical():
    """位相同守卫(A27 模板:同会话两次求解 array_equal)。

    ★ 本提交对 lifting 路径是**纯新增分支**(wc=None),位相同靠构造保持;
    这条锁的是「新分支没有经共享状态泄漏进 lifting 路径」这一可测形态。
    提交前/后的 phi 逐位证据在 W2.2a 报告里(改动前缓存 vs 改动后重算)。
    红了 = lifting 路径的数值被这次改动动了 —— 预注册 §1 的 FAIL 形态。
    """
    from tests.conftest import REPO_ROOT

    # A27 的算例与 recipe(M1 门配方),两条腿:默认 + B31/R19 钉扎
    KW = dict(upwind_c=1.0, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80,
              n_picard_seed=0)
    mc, wc = cut_wake(read_mesh(REPO_ROOT / "cases" / "meshes"
                                / "naca0012_2.5d" / "xcoarse.msh"))
    z = np.zeros(wc.n_stations)
    for kw_extra in (dict(), dict(tip_taper=z, gamma_target=z)):
        a = solve_newton_lifting(mc, wc, m_inf=0.80, alpha_deg=1.25,
                                 **KW, **kw_extra)
        b = solve_newton_lifting(mc, wc, m_inf=0.80, alpha_deg=1.25,
                                 **KW, **kw_extra)
        assert np.array_equal(np.asarray(a["phi"]), np.asarray(b["phi"])), (
            f" lifting 路径(phi)不再逐位确定,kw_extra={kw_extra} —— "
            "wc=None 分支泄漏进了 lifting 路径或共享状态被破坏")
        assert np.array_equal(np.asarray(a["gamma"]), np.asarray(b["gamma"]))


def test_c5_entropy_correction_on_by_default_and_converges(sphere, newton_m03):
    """熵修正在新入口**可用且默认 ON**(生产配置):共享腿本身就是熵 ON 的
    收敛腿;亚声速状态下 sigma ≡ 1(无激波 ⇒ 修正必须逐位惰性)。"""
    mesh, _, _ = sphere
    r_on = newton_m03
    assert r_on["entropy_correction"] is True, (
        "熵修正不再是默认 ON —— 生产配置被改了(GS1b.11 裁决)")
    assert r_on["sigma_transport_converged"], (
        "熵输运不收敛 —— 按 GS1.4 合同驱动器不许报收敛,这里却收敛了/或真失败")
    assert r_on["n_shock_cells"] == 0 and r_on["sigma_min"] == 1.0, (
        f"亚声速 sphere 出现激波单元/非 1 sigma(n_shock="
        f"{r_on['n_shock_cells']}, sigma_min={r_on['sigma_min']})—— "
        "熵修正在无激波流场不再惰性")
    r_off = solve_newton_nonlifting(mesh, m_inf=M_INF, tol_residual=TOL_SHARP,
                                    entropy_correction=False)
    assert r_off["converged"]
    assert np.array_equal(np.asarray(r_on["phi"]), np.asarray(r_off["phi"])), (
        "熵 ON/OFF 在亚声速 sphere 上不逐位相同 —— sigma ≡ 1 的惰性是逐位的,"
        "红了意味着修正泄漏进了无激波流场(或 ON 腿的额外残差求值污染了迭代)")


def test_kutta_options_rejected_on_raw_mesh(sphere):
    """wc=None 时 Kutta 概念(tip_taper / gamma_target / 压力估计子)必须
    显式报错,而不是被静默忽略(静默接受无意义参数是上次审计的同型缺陷)。"""
    mesh, _, _ = sphere
    with pytest.raises(ValueError, match="wake cut"):
        NewtonWorkspace(mesh, None, tip_taper=np.ones(1))
    with pytest.raises(ValueError, match="wake cut"):
        NewtonWorkspace(mesh, None, gamma_target=np.zeros(1))
    with pytest.raises(ValueError, match="wake cut"):
        NewtonWorkspace(mesh, None, kutta_estimator="pressure")


class TestPrescribedDirichlet:
    """W2.2b:处方 Dirichlet 模式(C05/C08 的入口形状)的机制锁 —— sphere
    coarse 上跑,与上面的 farfield 默认腿同一个网格族。

    ★ 等价位同腿就是预注册 §4 第二条 kill criterion 的**适定形式**:把远场
    节点上的自由流值原样作为处方值给出,与默认 farfield 路径**本该等价**
    ⇒ 必须逐位相同(同一批 dir_red/vals0_red ⇒ 同一条求解链路)。红了 =
    处方链路的实现错了,不是容差问题。
    """

    def test_pinned_bitexact_and_farfield_equivalent(self):
        from tests.conftest import REPO_ROOT

        mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "sphere_shell"
                         / "coarse.msh")
        ff = np.unique(mesh.boundary_faces["farfield"])
        vals = mesh.nodes[ff, 0].copy()     # alpha=0, u_inf=1 自由流,逐位等于默认
        a = solve_newton_nonlifting(mesh, m_inf=M_INF)
        b = solve_newton_nonlifting(mesh, m_inf=M_INF, dirichlet_nodes=ff,
                                    dirichlet_values=vals)
        assert a["converged"] and b["converged"], (
            "处方腿/默认腿不收敛 —— 前提变了")
        assert np.array_equal(np.asarray(b["phi"])[ff], vals), (
            "收敛态的处方节点值不等于处方值 —— Dirichlet 没压住(读回断言)")
        assert np.array_equal(np.asarray(a["phi"]), np.asarray(b["phi"])), (
            "处方(自由流等价值)与 farfield 默认不逐位相同 —— 预注册 §4 kill:"
            "本该等价的两条链路给出了不同答案,实现错了")

    def test_validation_and_mutual_exclusion(self):
        from tests.conftest import REPO_ROOT

        mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "sphere_shell"
                         / "coarse.msh")
        ff = np.unique(mesh.boundary_faces["farfield"])
        vals = mesh.nodes[ff, 0].copy()
        # 只给一个、长度不齐、节点越界 —— 全部必须显式报错
        with pytest.raises(ValueError, match="together"):
            solve_newton_nonlifting(mesh, m_inf=M_INF, dirichlet_nodes=ff)
        with pytest.raises(ValueError, match="together"):
            solve_newton_nonlifting(mesh, m_inf=M_INF, dirichlet_values=vals)
        with pytest.raises(ValueError, match="equal length"):
            solve_newton_nonlifting(mesh, m_inf=M_INF,
                                    dirichlet_nodes=ff[:3],
                                    dirichlet_values=vals)
        with pytest.raises(ValueError, match="out of range"):
            solve_newton_nonlifting(
                mesh, m_inf=M_INF,
                dirichlet_nodes=np.append(ff, len(mesh.nodes)),
                dirichlet_values=np.append(vals, 0.0))
        # 与 lifting 形状互斥(wc 给定时不许带处方)
        from pyfp3d.mesh.wake_cut import cut_wake
        mc, wc = cut_wake(read_mesh(REPO_ROOT / "cases" / "meshes"
                                    / "naca0012_2.5d" / "xcoarse.msh"))
        with pytest.raises(ValueError, match="wc=None"):
            NewtonWorkspace(mc, wc, dirichlet_nodes=np.array([0]),
                            dirichlet_values=np.array([0.0]))
        # 与复用 workspace 互斥(处方只能进 NewtonWorkspace 构造器)
        ws = NewtonWorkspace(mesh, None)
        with pytest.raises(ValueError, match="workspace"):
            solve_newton_nonlifting(mesh, m_inf=M_INF, workspace=ws,
                                    dirichlet_nodes=ff, dirichlet_values=vals)
