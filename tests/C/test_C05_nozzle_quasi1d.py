r"""C05 — 准一维 Laval 喷管：**项目里唯一带激波的精确解**（C 类：解析解 + 图）。

载体在 `tests/_nozzle_case.py`（phase 2 的 GS1.1 试验台按 AST 抽回活树；那份文件的
docstring 说明了它护住 `pyfp3d/kernels/` 的人工密度/迎风/熵修正，而**护不到**生产驱动
`solve/newton.py` 的全局化那一层）。

★★★ **判据立在「位置」上，不是残差上 —— 而这个许可是本门自己证的，不是引来的。**
phase 5 的侦察记录留下一条禁令：「**不得**引用 h=0.05 的喷管腿 —— 全部 `reason=cap`
非收敛」。那条禁令针对的是**收敛性**主张。本门实测（`TestCriterionIsOnPosition`）：

| `n_max` | `\|R\|` | `x_shock` |
|---|---|---|
| 20 | 4.45e-08 | 11.96160764 |
| 40 | 4.04e-08 | 11.96160687 |
| 80 | 1.89e-08 | 11.96160818 |

迭代数翻两番，**残差只在 ~2e-8 的地板上晃**，而**位置只动 1.3e-06 = 2.6e-05 个单元**。
⇒ 那条腿对**位置**可用、对**残差**仍不可用。**把禁令变成一个带自检的条件许可**，
而不是照抄或照违。

★★★ **本门第二半更值钱：一个残差到机器零的解可以错 40 个单元。**
从 `x_s = 8` 的精确解出发（进出口 Dirichlet 仍是 `x_s = 12` 的真值），nx = 200 的
Newton 收敛到 **`converged=True`, `reason=tol`, |R| = 1.30e-15**，而激波停在
**x = 7.98 而不是 12.0**，即 **−40.17 个单元 / 喷管全长的 20 %**。
★★ 而 `verify_uniqueness` 证明**连续问题的解是唯一的**（Δφ_shocked 19.3186 >
Δφ_subsonic_max 17.6307）⇒ 那个伪根是**离散算子**的，连续唯一性挡不住它。
⇒ **在一个有真值的案例上证明了「`converged` ∧ 残差机器零」不是正确性证书** ——
`bench/usability.py` 在翼型上只能拿「共识」当锚点，这里锚点是**解析解**。

★ nx = 100 的同一条腿更极端：收敛到 |R| = 1.32e-12，而 `shock_from_profile` 返回
**nan** —— 整场亚声速，激波**根本不存在**了。
"""
import os

import numpy as np
import pytest

from tests import _nozzle_case as NZ
from tests._gate_evidence import assert_matches_committed, fmt
from tests.conftest import REPO_ROOT, gate_figures_enabled

M_INF = 0.80
UPWIND_C = 1.5
#: 阶梯：h = 20/nx。nx = 400 只用 n_max = 20（位置早已冻住，见 docstring 的表），
#: 于是本门总计 ~11 s 而不是 ~30 s。
LEVELS = ((100, 80), (200, 80), (400, 20))

#: —— 实测（2026-08-26，本文件）——
#: err_x     -0.13015 / -0.06845 / -0.03839
#: err_cells -0.6507  / -0.6845  / -0.7678
#: 阶        0.927 (h 0.2->0.1) / 0.834 (0.1->0.05)
ORDER_MIN = 0.70          # 判据，非实测值：实测最小 0.834，留约 20 % 余量
CELLS_MAX = 1.00          # 亚单元精度：实测最大 |err_cells| = 0.768
POS_INVARIANCE_CELLS = 1e-3   # 实测 2.6e-05 个单元，判据留 38x 余量
WRONG_CELLS_MIN = 10.0
#: —— W1.2 / H7 的三条（标定，不是实测值）——
MDOT_IMBALANCE_MAX = 1e-10    # 结构守卫；实测 1.3e-16 / 3.5e-15 / 4.7e-15
MDOT_REL_ERR_COARSEST = 0.005  # nx=100 对精确 mdot：实测 0.001986（2.5x 余量）
MDOT_CONTRACTION = 1.5         # 逐级收缩；实测 2.15x / 2.07x（阶 ~1.05）    # 反例腿：实测 -40.17


#: ★★ 记忆化解析解。**实测它才是本门的成本**，不是求解：`exact_solution` 在
#: `n_quad = 40001` 个点上做 Python 循环求根（每点两个根），一次 ~8 s，而每条腿调它
#: 两次。加缓存后本门 75 s -> ~40 s，**数字一位不变**（纯函数，只有 m_inf 与 x_s 两个键）。
#: ★ 缓存放在这里而不是 `tests/_nozzle_case.py`，是为了让那份文件保持与归档
#: **逐字一致的 AST 抽取结果** —— 它的 docstring 是这么声明的。
_EXACT = {}


def _exact(x_s=None):
    key = (M_INF, x_s)
    if key not in _EXACT:
        _EXACT[key] = (NZ.exact_solution(M_INF) if x_s is None
                       else NZ.exact_solution(M_INF, x_s=x_s))
    return _EXACT[key]


def _solve(nx, n_max, x_s_init=None, upwind_c=UPWIND_C):
    """一条腿。Dirichlet 数据**永远**取自 x_s = 12 的精确解；只有初值变。"""
    ex = _exact()
    ny = max(6, nx // 16)
    mesh = NZ.nozzle_mesh(nx, ny)
    sysd = NZ.DuctSystem(mesh, m_inf=M_INF, upwind_c=upwind_c)
    phi_bc = ex["phi_of_x"](mesh.nodes[:, 0])
    if x_s_init is None or x_s_init == ex["x_s"]:
        phi0 = phi_bc.copy()
    else:
        phi0 = _exact(x_s_init)["phi_of_x"](mesh.nodes[:, 0])
        phi0[sysd.dir_nodes] = phi_bc[sysd.dir_nodes]      # 真值边界不动
    phi, info = sysd.newton(phi0, n_max=n_max, tol=1e-11)
    #: ★★ W1.2 / H7: 残差导出的边界质量流。R_i = ∫ rho~ grad(phi).grad(N_i) dV，
    #: 分部积分给 R_i = ∮ N_i rho dphi/dn dS（内点的体积项在收敛解上为零），
    #: 而一个边界组上 sum_i N_i = 1 ⇒ **该组节点的 R 之和就是穿过它的质量流**。
    #: 零额外求解，零新机器。
    R, _ = sysd.residual(phi)
    inlet = np.unique(mesh.boundary_faces["inlet"].ravel())
    outlet = np.unique(mesh.boundary_faces["outlet"].ravel())
    dz = float(np.ptp(mesh.nodes[:, 2]))
    mdot_in, mdot_out = float(R[inlet].sum()), float(R[outlet].sum())
    mdot_exact = float(ex["mdot"]) * dz          # 精确解是二维的，网格是一层拉伸
    xc, ux = NZ.element_u(sysd, phi)
    x_shock, n_sup, xb, ub = NZ.shock_from_profile(xc, ux, ex["u_star"], nx)
    h = NZ.LENGTH / nx
    err = (x_shock - ex["x_s"]) if np.isfinite(x_shock) else float("nan")
    return dict(nx=nx, h=h, n_max=n_max, x_shock=float(x_shock),
                err_x=float(err), err_cells=float(err / h),
                converged=bool(info["converged"]), reason=info["reason"],
                n_newton=int(info["n_newton"]),
                residual=float(info["residual_history"][-1]),
                x_s_exact=float(ex["x_s"]), u_star=float(ex["u_star"]),
                mdot_in=mdot_in, mdot_out=mdot_out, mdot_exact=mdot_exact,
                #: 守恒：进出之差（相对精确 mdot）；精度：量值对精确解的相对偏差
                mdot_imbalance=abs(mdot_in + mdot_out) / mdot_exact,
                mdot_rel_err=abs(abs(mdot_out) - mdot_exact) / mdot_exact,
                xb=xb, ub=ub)


@pytest.fixture(scope="module")
def sweep():
    return {nx: _solve(nx, n_max) for nx, n_max in LEVELS}


class TestExactSolutionIsWellPosed:
    """先证明参照配得上「精确解」这个名字，再拿它当尺子。"""

    def test_the_shocked_solution_is_the_only_one_meeting_the_imposed_phi(self):
        """★★ 喉部固定质量流 ⇒ 解族只由 x_s 参数化，Δφ 对 x_s **严格单调**；
        全亚声速（未壅塞）族的 Δφ 落在**不相交的、更低的**区间。

        ★ 这不是引文而是可算的：`verify_uniqueness` 实测
        Δφ_shocked **19.3186** > Δφ_subsonic_max **17.6307**，余量 1.688。
        ★★★ 而**唯一性只在连续层面成立** —— 离散算子照样有伪根，见
        `TestConvergedIsNotCorrect`。这两条必须一起读，否则会把连续唯一性
        误当成「求解器不可能算错」。
        """
        dp, dsub, unique = NZ.verify_uniqueness(M_INF)
        assert unique, (
            f"Δφ_shocked {dp:.6f} 未高于全亚声速族的上确界 {dsub:.6f} —— "
            "边界数据不再唯一确定激波位置，本门的整把尺子失效")
        assert dp - dsub > 1.0, f"唯一性余量只有 {dp - dsub:.4f}（实测 1.688）"


class TestShockPositionConvergence:
    """裁决③：C 类判**收敛阶** ∧ **绝对误差不过大**。"""

    def test_position_order(self, sweep):
        """★ 相邻两级的**物理位置**误差阶（实测 0.927 / 0.834）。

        ★★ 亚一阶是**预期而非缺陷**：激波被捕捉在大致固定的**单元数**内
        （`err_cells` 0.651 → 0.685 → 0.768，缓慢增长），于是物理误差 ~ h 但带一个
        略慢于 1 的修正。**反向落点**：若人工密度或迎风被改坏，阶会掉向 0 且
        `err_cells` 会开始随 h 增长 —— 下一条正是判那个。
        """
        ns = [nx for nx, _ in LEVELS]
        for i in range(1, len(ns)):
            a, b = sweep[ns[i - 1]], sweep[ns[i]]
            p = float(np.log(abs(a["err_x"]) / abs(b["err_x"]))
                      / np.log(a["h"] / b["h"]))
            assert p >= ORDER_MIN, (
                f"位置误差阶 h={a['h']:.3f} -> {b['h']:.3f} = {p:.3f} < {ORDER_MIN}"
                f"（实测基线 0.927 / 0.834）")

    def test_shock_lands_within_one_cell_at_every_level(self, sweep):
        """★★ 绝对判据用**单元**计，不用长度计 —— 长度判据会随 h 自动变松，
        而「激波落在几个单元内」才是这个数值方法真正的能力陈述。
        实测 0.651 / 0.685 / 0.768，全部亚单元。
        """
        for nx, _ in LEVELS:
            d = sweep[nx]
            assert abs(d["err_cells"]) <= CELLS_MAX, (
                f"nx={nx}: 激波偏离 {d['err_cells']:.3f} 个单元 > {CELLS_MAX}"
                f"（实测 0.651/0.685/0.768）")

    def test_error_sign_and_monotonicity(self, sweep):
        """★ 三级的物理误差必须逐级下降，且**始终在精确位置上游**（负号）。
        符号是内容：捕捉格式把激波抹在它上游的若干单元里，符号翻转说明
        抹的方向变了，那是迎风方向出了问题。
        """
        e = [sweep[nx]["err_x"] for nx, _ in LEVELS]
        assert all(x < 0 for x in e), f"误差符号不全为负（激波应在精确位置上游）：{e}"
        assert abs(e[0]) > abs(e[1]) > abs(e[2]), f"|err_x| 非单调：{e}"


class TestCriterionIsOnPosition:
    """★★★ 本门为什么可以用一条 `reason=cap` 的腿 —— 许可在这里自证。"""

    def test_position_is_invariant_to_the_iteration_cap(self):
        """残差有地板，位置没有。实测 n_max 20 -> 40：|R| 4.45e-08 -> 4.04e-08
        （只降 1.1×，是地板不是收敛），而 x_shock 只动 **7.7e-07 = 1.5e-05 个单元**。

        ★★ 这条一旦红，意味着 h=0.05 那条腿的**位置读数也不再可信**，
        于是本门余下的判据必须退回两级阶梯 —— 所以它排在收敛判据之外单列。
        """
        a = _solve(400, n_max=20)
        b = _solve(400, n_max=40)
        moved = abs(b["x_shock"] - a["x_shock"]) / a["h"]
        assert moved < POS_INVARIANCE_CELLS, (
            f"n_max 20 -> 40 让激波位置移动 {moved:.3e} 个单元 "
            f"> {POS_INVARIANCE_CELLS}（实测 1.5e-05）—— "
            "位置不再对迭代数不变，这条腿对位置判据也失去许可")
        #: ★ 同时记下残差**没有**收敛，免得有人把本条读成「其实收敛了」
        assert not a["converged"] and a["reason"] == "cap", (
            f"nx=400/n_max=20 现在报 converged={a['converged']} "
            f"reason={a['reason']} —— 若它真的收敛了，本条的前提变了，"
            "请重新测量并更新 docstring 的表")


class TestConvergedIsNotCorrect:
    """★★★ 本门最值钱的一半：在**有真值**的案例上，`converged` 不是正确性证书。"""

    def test_a_machine_zero_residual_can_be_wrong_by_forty_cells(self):
        """从 x_s = 8 的精确解出发（**边界数据仍是 x_s = 12 的真值**），
        nx = 200 收敛到 |R| = 1.30e-15，而激波停在 x ≈ 7.98。

        ★★ 本条断言的是**缺陷仍然存在**。若它变红，那是好消息，但**必须走勘误**：
        `bench/usability.py` 的整套论证、以及 R23 的判据设计，都以「converged ∧ 0 钳制
        不是充分可用性判据」为前提。前提若变了，那批结论要重读。
        """
        d = _solve(200, n_max=80, x_s_init=8.0)
        assert d["converged"] and d["residual"] < 1e-11, (
            f"反例腿不再收敛（converged={d['converged']}, |R|={d['residual']:.2e}）—— "
            "本条要证的是「收敛且错」，不收敛的腿证不了")
        assert abs(d["err_cells"]) > WRONG_CELLS_MIN, (
            f"反例腿只错 {d['err_cells']:.2f} 个单元 < {WRONG_CELLS_MIN}"
            f"（实测 -40.17）。★ 若求解器真的被修好了这是好消息，但"
            "请按纪律 11 走勘误：usability/R23 的论证前提变了")

    def test_the_coarse_leg_loses_the_shock_entirely(self):
        """★ 同一条腿在 nx = 100 上更极端：收敛到 |R| ~ 1e-12，而声速穿越
        **根本不存在**（`shock_from_profile` 返回 nan）—— 整场亚声速。

        ★★ 记它是因为两种失败**形状不同**：nx=200 是「激波在错的位置」，
        nx=100 是「没有激波」。一个只看位置误差的判据**看不见后者**
        （nan 不参与比较），所以它需要自己一条。
        """
        d = _solve(100, n_max=80, x_s_init=8.0)
        assert d["converged"], f"前提变了：该腿不再收敛（{d['reason']}）"
        assert not np.isfinite(d["x_shock"]), (
            f"nx=100/perturbed 现在有激波了（x={d['x_shock']:.4f}）—— "
            "本条记录的是「收敛到一个无激波解」这种失败形状，前提变了请重新测量")


class TestGlobalMassConservation:
    r"""★★★ W1.2 / H7（2026-09-06）：**全局质量守恒**，此前全仓一条都没有。

    求解的是守恒形式 ∇·(ρ∇φ) = 0，而在这一轮之前**没有任何一道门断言过
    「进出边界的质量流」**（A11 守的是单元**几何体积**，是另一回事）。
    喷管是立这条门最干净的地方：进出口 Dirichlet、壁面自然 BC、
    **而且 `mdot` 有精确值**。

    ★★ 两条判据的强度**不一样，必须分开说**：

    | 判据 | 实测 | 它有多强 |
    |---|---|---|
    | (a) 进 + 出 = 0 | 1.3e-16 … 4.7e-15 | ★ **接近恒等式**：Σ_i R_i ≡ 0 是装配自带的（Σ_i N_i ≡ 1），收敛解的自由节点残差又 ~0。**它护不住离散精度**，护的是**结构**：边界组取错、Dirichlet 行没被覆盖、节点集重叠——都会让它跳好几个数量级 |
    | (b) 量值对精确 mdot | **+0.1986 % → +0.0922 % → +0.0445 %** | ★★★ **真判据**：对**精确解**的误差，逐级收缩 2.15x / 2.07x，h 每次减半 ⇒ 观察阶 ≈ **1.05**（一阶，与捕获激波一致） |

    ⇒ **(b) 才是这道门的价值**，(a) 是廉价的结构守卫，两者都记录。

    ★★★ **G-TEETH 是量出来的，不是"验证过会红"**（2026-09-06，W3.1 的做法提前用）：

    | 扰动 | (b) mdot_rel_err | (a) imbalance |
    |---|---|---|
    | 无 | 0.001986 green | 1.3e-16 green |
    | ρ̃ ×(1+1e-6) | 0.001987 green | 1.2e-15 green |
    | ρ̃ ×(1+1e-3) | 0.002988 green | 7.2e-15 green |
    | ρ̃ ×(1+3e-3) | 0.004992 green（**擦线**） | 4.3e-13 green |
    | ρ̃ ×(1+5e-3) | 0.006996 **RED** | 3.6e-12 green |
    | 进口节点集少一个 | — | **5.8e-02 RED** |
    | 进出口取成同一组 | — | **2.0e+00 RED** |

    ⇒ **(b) 的最小可见退化 ≈ ρ̃ 上 0.3 % 的偏置**。计划书里建议的"1e-6 常数偏置
    必红"**实测差了四个数量级** —— 记在这里，免得下次照抄那句话。
    ⇒ **(a) 对 ρ̃ 的均匀偏置完全不动**（进出两侧同比例缩放，差仍为零），
    但对**结构性**错误跳 14 个数量级。**两条判据量的不是同一件事**，
    这张表就是证据。
    ★ G-DOMAIN：相反结果 = mdot 偏低（`mdot_out < mdot_exact`）也照样判 ——
    判据取的是**绝对值**，方向另行 RECORDED（实测三级都偏高，见下）。
    """

    def test_flux_sign_convention_is_the_measured_one(self, sweep):
        """进口为负、出口为正 —— **量出来的**，不是从公式抄的。

        ★ 它是接线检查：把两个边界组弄反、或把某一组的节点集算错，这条先红。
        """
        for nx, _ in LEVELS:
            d = sweep[nx]
            assert d["mdot_in"] < 0.0 < d["mdot_out"], (
                f"nx={nx}: sign convention moved "
                f"(in {d['mdot_in']:+.4e}, out {d['mdot_out']:+.4e})")

    def test_mass_is_conserved_between_inlet_and_outlet(self, sweep):
        """(a) 结构守卫：|进 + 出| / mdot_exact ≤ 1e-10（实测 ≤ 4.7e-15）。

        ★ 阈值离实测有 ~2e4 倍余量，**这是有意的**：这个量是机器精度量级的
        近恒等式，能让它动的是**结构性**错误（边界组、Dirichlet 覆盖、节点集），
        而那种错误会把它抬高好几个数量级，不是抬高 2 倍。
        """
        for nx, _ in LEVELS:
            d = sweep[nx]
            assert d["mdot_imbalance"] <= MDOT_IMBALANCE_MAX, (
                f"nx={nx}: inlet+outlet imbalance "
                f"{d['mdot_imbalance']:.3e} > {MDOT_IMBALANCE_MAX:.0e}")

    def test_mass_flux_matches_the_exact_solution(self, sweep):
        """(b) 真判据：最粗一级对**精确 mdot** 的相对误差 ≤ 带。"""
        d = sweep[LEVELS[0][0]]
        assert d["mdot_rel_err"] <= MDOT_REL_ERR_COARSEST, (
            f"nx={LEVELS[0][0]}: |mdot| off exact by "
            f"{100 * d['mdot_rel_err']:.4f} % > "
            f"{100 * MDOT_REL_ERR_COARSEST:.2f} %")

    def test_mass_flux_error_contracts_with_refinement(self, sweep):
        """(b) 续：逐级收缩 —— 一个不收缩的守恒误差不是离散误差。"""
        e = [sweep[nx]["mdot_rel_err"] for nx, _ in LEVELS]
        for a, b, (na, _), (nb, _) in zip(e[:-1], e[1:], LEVELS[:-1], LEVELS[1:]):
            assert a / b >= MDOT_CONTRACTION, (
                f"nx {na}->{nb}: mdot error contracts only {a / b:.2f}x "
                f"(< {MDOT_CONTRACTION}) -- {a:.3e} -> {b:.3e}")

    def test_the_sign_of_the_discretisation_bias_is_recorded(self, sweep):
        """RECORDED，不设方向门：实测三级 `|mdot|` 都**偏高**于精确值。

        ★ 只断言「三级同号」——**一致性**是可判的，而"应该偏高还是偏低"
        本项目没有可引用的出处（同 D05 拒绝给激波偏移登记符号）。
        """
        signs = [np.sign(abs(sweep[nx]["mdot_out"]) - sweep[nx]["mdot_exact"])
                 for nx, _ in LEVELS]
        assert len(set(signs)) == 1, (
            f"mdot bias changes sign across the ladder: {signs} -- "
            "that is a finding, not a tolerance question")


class TestCommittedEvidenceIsLoadBearing:
    r"""★★★ 新鲜计算 vs 已提交 `summary.csv`。设计与三个坑见 `tests/_gate_evidence.py`。
    ★ 零额外计算：用 `sweep` 已算好的结果（**反例腿不在 sweep 里**，它由 §
    `TestConvergedIsNotCorrect` 单独求解，因此不进本锁）。"""

    def test_fresh_run_reproduces_the_committed_summary(self, sweep):
        fresh = {("ladder", str(nx)): dict(x_shock=sweep[nx]["x_shock"],
                                  err_x=sweep[nx]["err_x"],
                                  err_cells=sweep[nx]["err_cells"],
                                  #: ★ W1.2 加的三列 —— 一道门的数字必须被锁，
                                  #: 否则它就是"只写不读"那一族（2026-08-28 实测
                                  #: 19 个证据 CSV 里 16 个无人读）。
                                  mdot_out=sweep[nx]["mdot_out"],
                                  mdot_exact=sweep[nx]["mdot_exact"],
                                  mdot_rel_err=sweep[nx]["mdot_rel_err"])
                 for nx, _ in LEVELS}
        n = assert_matches_committed(
            os.path.join(str(REPO_ROOT), "cases", "gates", "C05_nozzle_quasi1d"),
            fresh, ("x_shock", "err_x", "err_cells",
                    "mdot_out", "mdot_exact", "mdot_rel_err"),
            key_of=lambda r: (r["leg"], r["nx"]),
            refresh_hint="PYFP3D_GATE_FIGURES=1 pytest tests/C/test_C05_nozzle_quasi1d.py")
        assert n >= 18, f"只比了 {n} 个数（3 级 x 6 列 = 18）"

@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_nozzle_figure(sweep, gate_evidence_dir):
    import csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ex = _exact()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(ex["x"], ex["u"], "k-", lw=1.5, label="exact (quasi-1D)")
    for nx, _ in LEVELS:
        d = sweep[nx]
        ax[0].plot(d["xb"], d["ub"], ".", ms=3, label=f"nx={nx} (h={d['h']:.3f})")
    ax[0].axhline(ex["u_star"], color="0.6", ls=":", lw=1, label="u*")
    ax[0].axvline(ex["x_s"], color="r", ls="--", lw=1, label="exact x_s")
    ax[0].set_xlabel("x"), ax[0].set_ylabel("u")
    ax[0].set_title(f"C05 Laval nozzle, M_inf = {M_INF}, C = {UPWIND_C}")
    ax[0].legend(fontsize=8), ax[0].grid(alpha=.3)

    hs = np.array([sweep[nx]["h"] for nx, _ in LEVELS])
    ex_err = np.array([abs(sweep[nx]["err_x"]) for nx, _ in LEVELS])
    ax[1].loglog(hs, ex_err, "o-", label="|err_x| (length)")
    ax[1].loglog(hs, [abs(sweep[nx]["err_cells"]) for nx, _ in LEVELS], "s-",
                 label="|err_cells|")
    ax[1].loglog(hs, ex_err[0] * hs / hs[0], "k--", lw=.8, label="order 1")
    ax[1].set_xlabel("h"), ax[1].set_ylabel("shock-position error")
    ax[1].set_title("position converges; cells stay sub-unit")
    ax[1].legend(fontsize=8), ax[1].grid(alpha=.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "c05_nozzle_shock.png"), dpi=130)
    plt.close(fig)

    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        #: ★ `leg` 列是 2026-08-28 加的：阶梯腿与反例腿都是 nx = 200，
        #: 没有它键就不唯一，证据锁会比错行（`tests/_gate_evidence.py` 的重复键守卫）。
        w.writerow(["leg", "nx", "h", "n_max", "converged", "reason", "n_newton",
                    "residual", "x_shock", "x_s_exact", "err_x", "err_cells",
                    "mdot_in", "mdot_out", "mdot_exact", "mdot_imbalance",
                    "mdot_rel_err"])

        def _row(tag, nx, n_max, d):
            return [tag, nx, f"{d['h']:.5f}", n_max, d["converged"], d["reason"],
                    d["n_newton"], f"{d['residual']:.3e}", fmt(d["x_shock"]),
                    d["x_s_exact"], fmt(d["err_x"]), fmt(d["err_cells"]),
                    fmt(d["mdot_in"]), fmt(d["mdot_out"]), fmt(d["mdot_exact"]),
                    fmt(d["mdot_imbalance"]), fmt(d["mdot_rel_err"])]

        for nx, n_max in LEVELS:
            w.writerow(_row("ladder", nx, n_max, sweep[nx]))
        w.writerow(_row("counterexample", 200, 80,
                        _solve(200, n_max=80, x_s_init=8.0)))
