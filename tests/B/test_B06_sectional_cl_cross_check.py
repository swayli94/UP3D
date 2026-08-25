"""B06 — 截面升力系数：**Cp 积分 vs Kutta-Joukowski** 的交叉校核（B 类）。

★★ 2026-08-25 使用者要求：加入「由截面压力分布积分出的截面 cl」，并与
`sectional_cl_from_gamma` 对比。本文件是那条对比。

## 这条门在检验什么

**两个独立的算法必须给出同一个物理量。**
- `sectional_cl_from_gamma`：`2Γ/(U c)` —— 走**环量**，来自 Kutta 条件与尾迹跳跃；
- `sectional_cl_from_cp`：`∮ Cp` —— 走**壁面压力**，来自 φ 的壁面切向梯度。

⇒ 它们共享的只有 φ。**一致 ⇒ Kutta 行、尾迹跳跃、壁面恢复、Cp 定律四者自洽；
不一致 ⇒ 其中至少一处错了，而单看任一个都看不出来。**
★ 全翼层面已有先例（G2.4 把 `cl_p` 与 `cl_KJ` 卡在 < 1 %）；**截面层面此前没有**。

## ★★★ 为什么判据只立在 2.5-D —— 使用者 2026-08-25 的方法论裁决

我最初把这条门建在 M6 上并想断言「平滑让两法更一致」。**两处都错**：

1. **M6 不能当判据的载体。** 非结构三角形网格**展向分布不均匀**，截面切割时
   **上下表面并不对称** ⇒ 那里的差**混着网格效应**，无法归给求解器。
   实测 M6 coarse 五站 RMS 差 **0.0116（约 4 %）**，而它**不可归因**。
2. **我给的机理是编的。** 我写过「积分本身把 ±配对恢复了」—— 未验证；
   使用者指出上下表面不对称时抵消无从保证。★ **我把解释当成了定位**
   （这句本项目早就写过，我又犯了一次）。

⇒ **唯一合理的判据要建立在排除展向不均匀的网格上**，即 **2.5-D 拉伸网格**。
那里的差**只能**归给求解器/后处理本身。

## 实测（2.5-D，α 2°，2026-08-25）

| 网格 | M | `cl_γ` | Cp 积分 | 相对差 | 全翼 `cl_p` |
|---|---|---|---|---|---|
| coarse | 0.00 | 0.23683 | 0.23450 | **0.987 %** | 0.23443 |
| coarse | 0.50 | 0.28019 | 0.27788 | **0.826 %** | 0.27759 |
| medium | 0.00 | 0.24073 | 0.24030 | ★ **0.178 %** | 0.23994 |
| medium | 0.50 | 0.28528 | 0.28485 | ★ **0.150 %** | 0.28437 |

★★ **coarse ~0.9 % → medium ~0.16 %，约 5.5x 收缩** ⇒ 这是**离散误差在收敛**，
不是内部错误。**收缩本身才是证据**，单看一个绝对值说明不了。
★ 而全翼 `cl_p` 与截面 Cp 积分**几乎重合**（0.23450 / 0.23443）——
两条独立的压力积分路径互相印证。

## RECORDED（不设门）

- **平滑在 cl 这一项上让两法更不一致**（medium M0.5 0.150 % → 0.458 %；
  coarse 0.826 % → 1.421 %，四组同向），方向与 G6.3 在全翼上的一致。
  ★ **只记方向，不写机理** —— 机理未定。
- **M6 上 4 % 的差里，展向不均匀占多少、求解器占多少：未知。**
  2.5-D 分不开（没有三维效应），M6 分不开（有混杂）⇒ 需要一个
  **排除展向混杂的三维载体**。登记为开放问题。
"""
import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.surface import (sectional_cl_from_cp, sectional_cl_from_gamma,
                                 wall_force_coefficients)
from pyfp3d.solve.newton import solve_newton_lifting
from tests.conftest import REPO_ROOT

M_INF, ALPHA = 0.50, 2.0
#: medium 实测 0.150 %；带取 0.5 %，余量 3.3x。
AGREE_MAX_MEDIUM = 0.005
#: coarse 0.826 % -> medium 0.150 % = 5.5x；带取 2.0x，余量 2.8x。
CONTRACTION_MIN = 2.0


def test_flat_plate_constant_dcp_is_exact():
    """解析自检：常压差 `dCp = 1` 的平板给 `cn = 1`。零求解。

    ★ 先钉住积分器本身，否则下面用真实解做的对比**分不清是求解器错还是积分器错**。
    """
    c = dict(x_upper=np.array([0.0, 1.0]), cp_upper=np.array([-0.5, -0.5]),
             x_lower=np.array([0.0, 1.0]), cp_lower=np.array([0.5, 0.5]))
    assert abs(sectional_cl_from_cp(c, 0.0) - 1.0) < 1e-12


def test_axial_term_is_not_silently_dropped():
    """★★ 没有几何时 `alpha != 0` **必须报错**，不能悄悄只返回 cn。

    静默丢掉轴向项会让大迎角下的 cl **无声偏小** —— 安静的错误比响的错误糟。
    """
    c = dict(x_upper=np.array([0.0, 1.0]), cp_upper=np.array([-0.5, -0.5]),
             x_lower=np.array([0.0, 1.0]), cp_lower=np.array([0.5, 0.5]))
    assert abs(sectional_cl_from_cp(c, 0.0) - 1.0) < 1e-12
    with pytest.raises(ValueError, match="y_upper"):
        sectional_cl_from_cp(c, 5.0)


def _case(level, alpha=ALPHA):
    p = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / ("%s.msh" % level)
    if not p.exists():
        pytest.skip("naca0012_2.5d/%s.msh 未生成" % level)
    mc, wc = cut_wake(read_mesh(str(p)))
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=alpha,
                             precond="direct", n_newton_max=60)
    assert r["converged"] and r["n_limited"] == 0 and r["n_floored"] == 0, (
        "%s 未给出干净的收敛态（conv=%s lim=%s flr=%s）—— 交叉校核建立在它之上，"
        "所以这是前置条件" % (level, r["converged"], r["n_limited"], r["n_floored"]))
    dz = float(np.ptp(mc.nodes[:, 2]))
    c = section_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=M_INF, u_inf=1.0)
    cl_g = float(np.median(sectional_cl_from_gamma(
        np.asarray(r["gamma"], dtype=np.float64), chord=c["chord"], u_inf=1.0)))
    cl_cp = sectional_cl_from_cp(c, alpha)
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"],
                                r["phi"], alpha_deg=alpha, s_ref=dz, m_inf=M_INF)
    #: 只算 cn 的版本 —— 用来把**轴向分量**单独拎出来（见 test_the_axial_term_is_covered）
    c_no_geom = dict(c); c_no_geom["y_upper"] = None; c_no_geom["y_lower"] = None
    cn_only = sectional_cl_from_cp(c_no_geom, 0.0)
    return (cl_g, cl_cp, float(f["cl"]), abs(cl_cp - cl_g) / abs(cl_g), cn_only)


@pytest.fixture(scope="module")
def ladder():
    return {lv: _case(lv) for lv in ("coarse", "medium")}


#: ★★ 轴向项在生产迎角 alpha 2 上只占 cl 的 **0.16 %** —— G-TEETH 实测：把它的符号
#: 翻转，本文件**全绿**，即 `sectional_cl_from_cp` 的**轴向那一半没有被检验**。
#: 这与本季记的「窗口覆盖不到被改动的区域」同族。
#: 实测 alpha 6：轴向项 **1.01 %** of cl；alpha 12.86 在 M0.5 medium 上**不收敛**
#: （2925 limited / 727 floored）⇒ 6 度是可用的最大值。
ALPHA_AXIAL = 6.0


@pytest.fixture(scope="module")
def high_alpha():
    return _case("medium", alpha=ALPHA_AXIAL)


def test_the_section_curve_carries_its_geometry():
    """★ `y_upper/y_lower` 出来了，且与 `x_*` 等长。

    这是 2026-08-25 那次贯通（`_wall_section_points` 返回元数 3 → 4）的直接锁：
    `ys` 一旦又被丢掉，`sectional_cl_from_cp` 会退回「只有 cn」并在 alpha != 0 时报错。
    """
    _ = _case  # noqa: F841  (保持与上面同一个装配路径)
    p = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh"
    if not p.exists():
        pytest.skip("naca0012_2.5d/coarse.msh 未生成")
    mc, wc = cut_wake(read_mesh(str(p)))
    dz = float(np.ptp(mc.nodes[:, 2]))
    c = section_cp_curve(mc, np.zeros(len(mc.nodes)), z=0.5 * dz, m_inf=0.0)
    for side in ("upper", "lower"):
        y = c["y_%s" % side]
        assert y is not None, "y_%s 是 None —— 截面几何又被丢掉了" % side
        assert len(y) == len(c["x_%s" % side]), "y_%s 与 x_%s 不等长" % (side, side)


def test_cp_integration_agrees_with_kutta_joukowski_on_a_clean_mesh(ladder):
    """★★★ 干净网格（2.5-D，无展向混杂）上两法一致到 < 0.5 %。

    ★ 只立在 2.5-D：M6 上的差**混着展向不均匀**，不可归因（见模块 docstring）。
    """
    cl_g, cl_cp, cl_p, rel, _ = ladder["medium"]
    assert rel < AGREE_MAX_MEDIUM, (
        "2.5-D medium 上截面 cl 两法不一致：Cp 积分 %.6f vs KJ %.6f（%.3f %%，"
        "要求 < %.1f %%）。\n"
        "  ★ 它们只共享 φ ⇒ 不一致意味着 Kutta 行 / 尾迹跳跃 / 壁面恢复 / Cp 定律\n"
        "  至少有一处错了。★ 这张网格**排除了展向不均匀**，所以这个差不能推给网格。"
        % (cl_cp, cl_g, 100 * rel, 100 * AGREE_MAX_MEDIUM))


def test_the_disagreement_contracts_under_refinement(ladder):
    """★★ **收缩才是证据** —— 单看一个绝对值分不清「离散误差」与「内部错误」。

    实测 coarse 0.826 % -> medium 0.150 % = 5.5x。一个**不收缩**的差意味着
    某处有真错误，而不是分辨率不够。
    """
    rel_c, rel_m = ladder["coarse"][3], ladder["medium"][3]
    assert rel_m > 0.0, "medium 上两法完全相同 —— 检查是不是两条路退化成了同一条"
    ratio = rel_c / rel_m
    assert ratio >= CONTRACTION_MIN, (
        "两法之差没有随加密收缩（coarse %.3f %% -> medium %.3f %%，仅 %.2fx，"
        "要求 >= %.1fx）—— 不收缩的差是**错误**不是分辨率"
        % (100 * rel_c, 100 * rel_m, ratio, CONTRACTION_MIN))


def test_the_two_pressure_integrations_agree(ladder):
    """★ 全翼 `cl_p` 与截面 Cp 积分几乎重合（实测 0.28485 vs 0.28437）。

    两条**独立的压力积分路径**（全表面面积分 vs 截面围线积分）互相印证 ——
    它们与上面那条 KJ 对比是**不同的**校核：这条只查压力侧，不涉及环量。
    """
    for level in ("coarse", "medium"):
        _, cl_cp, cl_p, _, _ = ladder[level]
        rel = abs(cl_cp - cl_p) / abs(cl_p)
        assert rel < 0.01, (
            "%s：截面 Cp 积分 %.6f 与全翼 cl_p %.6f 差 %.3f %% —— "
            "两条压力积分路径应当几乎重合" % (level, cl_cp, cl_p, 100 * rel))


def test_the_axial_term_is_covered(high_alpha):
    """★★★ 覆盖 `sectional_cl_from_cp` 的**轴向那一半**。

    ★ 这条是 **G-TEETH 抓出来的空缺**，不是设计时想到的：在生产迎角 alpha 2 上把
    轴向项的符号翻转，本文件**全绿** —— 因为 NACA0012 对称且 `sin 2 deg = 0.035`，
    轴向项只占 cl 的 **0.16 %**。⇒ 那半**没有被检验**。

    实测 alpha 6：轴向项 **1.01 %** of cl，符号翻转会让 cl 动 ~2 %，越过下面的带。
    ★ alpha 12.86（reference_data 里有对应实验工况）**不可用** —— M0.5 medium 上
    不收敛（2925 limited / 727 floored），所以 6 度是可用的最大值。**这是覆盖的边界，
    不是选择**：轴向项在 alpha 2 的生产工况上仍然测不到。
    """
    cl_g, cl_cp, cl_p, rel, cn_only = high_alpha
    axial = cl_cp - cn_only * np.cos(np.deg2rad(ALPHA_AXIAL))
    frac = abs(axial) / abs(cl_cp)
    assert frac > 0.005, (
        "alpha %.1f 上轴向项只占 cl 的 %.3f %% —— 低于 0.5 %% 时符号错误检测不出来，"
        "这条门就覆盖不到它声称覆盖的那一半" % (ALPHA_AXIAL, 100 * frac))
    assert rel < AGREE_MAX_MEDIUM, (
        "alpha %.1f 上两法差 %.3f %%（要求 < %.1f %%）：Cp 积分 %.6f vs KJ %.6f。"
        "★ 这个迎角下轴向项占 %.2f %%，所以本条同时覆盖它"
        % (ALPHA_AXIAL, 100 * rel, 100 * AGREE_MAX_MEDIUM, cl_cp, cl_g, 100 * frac))
