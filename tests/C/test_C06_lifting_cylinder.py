"""C06 — 带环量圆柱：**最便宜的精确 Kutta 检验**（C 类：解析解 + 图）。

★★★ **这道门为什么值得存在。** 圆柱**没有尖尾缘**，所以 Kutta 条件选不出环量 ——
Γ 必须**被规定**。于是这道门能做一件别处做不到的事：**给进去的是 Γ，量出来的是
表面压力积分的升力**，两条完全不同的路，中间**没有尾缘奇点**污染。
Kármán–Trefftz（C07）能检验的是同一件事**加上**有限尾缘角；本门是它之前的台阶。

精确解（不可压位势，半径 a、来流 U、环量 Γ）::

    u_theta(theta) = -2 U sin(theta) - Gamma / (2 pi a)
    Cp(theta)      = 1 - (2 sin(theta) + Gamma / (2 pi a U))**2
    L = rho U Gamma   =>   cl = 2 Gamma / (U c),  c = 2a

★★ **摆法与翼型约定同构**：半径 0.5、圆心 (0.5, 0) ⇒ 弦长 1、尾缘 (1, 0)、前缘 (0, 0)。
这不是美学 —— `airfoil_wake_2d` 取 `pt[0]` 为尾缘、`argmin(x)` 为前缘、远场圆心
**硬编码在半弦 (0.5, 0)**，对齐之后**一行新的网格生成代码都不需要**。

★★ **几何误差已实测，不是假设**：折线经 `addSpline` 后壁面节点到精确圆的偏差
**max |r - a| = 1.74e-09（相对 3.5e-09）**，比本门要测的网格误差低**九个量级** ⇒
不是混淆项。

★★★ **判据必须全尺度加密，而这是实测逼出来的。** 第一版阶梯沿用了 NACA 生成器的
`h_far = min(3.0, 150*h)`，三级全被钳在 3.0 ⇒ **只加密壁面、体网格没细**，实测：

| 阶梯 | Cp 阶 (c->m) | cl 误差 |
|---|---|---|
| `h_far` 钳在 3.0（只细壁面） | **0.325** | 1.124 % -> 1.136 %，**不动** |
| `h_far = 150*h`（全尺度） | **1.686** | 3.477 % -> 1.136 % |

这正是 P11 在球上测到的**固定体网格污染地板**，同一个形状。⇒ 本门的阶梯**必须**让
`h_far` 跟着 `h_wall` 走；若有人把它钳回去，收敛阶会塌到 0.3 而绝对误差看起来还行。

★ 第四级 h = 0.01（25652 节点、30 s）**只用来定渐近率、不进门** —— 项目约定
`fine` 不在阶梯里。它测得 Cp 阶 1.329 / cl 阶 1.600，确认 c->m 的 1.686/1.614
不是偶然。
"""
import os

import numpy as np
import pytest

from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.extrude import extrude_single_layer
from pyfp3d.meshgen.planar import airfoil_wake_2d
from pyfp3d.post.surface import (wall_force_coefficients,
                                 wall_tangential_gradient_quadratic)
from pyfp3d.solve.picard import solve_laplace_lifting
from tests._gate_evidence import assert_matches_committed, fmt
from tests.conftest import REPO_ROOT, gate_figures_enabled

#: 半径与圆心 —— 对齐 `airfoil_wake_2d` 的翼型约定（弦长 1，远场圆心硬编码在 (0.5, 0)）
A, XC = 0.5, 0.5
U_INF = 1.0
#: 规定环量。0.6 给出 cl = 1.2，Cp 的驻点结构清晰而不至于把上表面推到强负压。
GAMMA = 0.6
#: 正则阶梯（`fine` 按项目约定不在门内）
LEVELS = (("xcoarse", 0.08), ("coarse", 0.04), ("medium", 0.02))
R_FAR = 15.0

#: —— 实测值（2026-08-26，本文件的 sweep）——
#: max|Cp|  8.6698e-01 / 3.2466e-01 / 1.0093e-01
#: rms Cp   4.0400e-01 / 1.5480e-01 / 5.1721e-02
#: cl_p     1.159284   / 1.158277   / 1.186373   （精确 1.2）
#: 阶 Cp    1.417 (xc->c) / 1.686 (c->m)；阶 cl 误差 -0.035 / 1.614
ORDER_MIN = 1.0          # ★ 判据，不是实测值：实测最小 1.417，留 40 % 余量
CP_MAX_MEDIUM = 0.15     # 实测 0.1009
CP_RMS_MEDIUM = 0.08     # 实测 0.0517
CL_REL_MEDIUM = 0.02     # 实测 0.01136
CL_CONTRACTION = 2.0     # coarse -> medium 至少收缩 2x；实测 3.477/1.136 = 3.06


def _exact_cp(theta):
    """Cp(theta) = 1 - (2 sin θ + Γ/(2πaU))**2 —— 符号约定**实测**确定（见下）。

    ★ 反号的话 max|ΔCp| 从 0.10 跳到 **1.60**，相差 16 倍，毫无歧义。
    这个约定是量出来的，不是从公式书上抄的。
    """
    return 1.0 - (2.0 * np.sin(theta) + GAMMA / (2.0 * np.pi * A * U_INF)) ** 2


def _circle_polyline(n_half):
    """闭合折线，与 `naca0012_coordinates` 同约定：尾缘 -> 上 -> 前缘 -> 下 -> 尾缘。"""
    th = np.linspace(0.0, np.pi, n_half + 1)
    upper = np.column_stack([XC + A * np.cos(th), A * np.sin(th)])
    th2 = np.linspace(np.pi, 2.0 * np.pi, n_half + 1)[1:]
    lower = np.column_stack([XC + A * np.cos(th2), A * np.sin(th2)])
    return np.vstack([upper, lower])


def _one_level(h, gamma=GAMMA):
    p2, tri, eg, ig = airfoil_wake_2d(
        _circle_polyline(max(80, int(round(2.0 / h)))), model_name="c06_cylinder",
        r_far=R_FAR, h_wall=h,
        #: ★★ 全尺度加密 —— 见 docstring 的表：钳住它会把阶塌到 0.325
        h_far=150.0 * h, h_wake=3.0 * h,
        dist_min=0.1, dist_max=6.0, wake_dist_max=1.5)
    dz = 2.0 * h
    mesh = extrude_single_layer(p2, tri, eg, ig, dz=dz, name="c06")
    mc, wc = cut_wake(mesh)
    n_station = int(np.unique(wc.te_station).size)
    r = solve_laplace_lifting(mc, wc, alpha_deg=0.0, u_inf=U_INF,
                              gamma_fixed=np.full(n_station, gamma))
    phi = np.asarray(r["phi"])
    wf = mc.boundary_faces["wall"]
    wall = np.unique(wf)
    grad = wall_tangential_gradient_quadratic(mc.nodes, wf, phi)
    cp = 1.0 - np.sum(grad[wall] ** 2, axis=1) / U_INF ** 2
    theta = np.arctan2(mc.nodes[wall, 1], mc.nodes[wall, 0] - XC)
    err = cp - (1.0 - (2.0 * np.sin(theta) + gamma / (2.0 * np.pi * A * U_INF)) ** 2)
    f = wall_force_coefficients(mc.nodes, mc.elements, wf, phi,
                                alpha_deg=0.0, u_inf=U_INF, s_ref=dz, m_inf=0.0)
    cl_exact = 2.0 * gamma / (U_INF * 2.0 * A)
    return dict(
        h=h, n_nodes=int(mc.nodes.shape[0]), n_wall=int(wall.size),
        theta=theta, cp=cp, err=err,
        cp_max=float(np.max(np.abs(err))), cp_rms=float(np.sqrt(np.mean(err ** 2))),
        cl=float(f["cl"]), cl_exact=cl_exact,
        cl_rel=float(abs(f["cl"] - cl_exact) / cl_exact),
        residual=float(r["residual_norm"]), n_kutta=int(r["n_kutta_updates"]),
        r_wall_dev=float(np.max(np.abs(
            np.hypot(mc.nodes[wall, 0] - XC, mc.nodes[wall, 1]) - A))),
    )


@pytest.fixture(scope="module")
def sweep():
    """三级阶梯，一次求解一级 —— 实测 ~3 s 总计，属快层。"""
    return {name: _one_level(h) for name, h in LEVELS}


def _order(a, b, ha, hb):
    return float(np.log(a / b) / np.log(ha / hb))


class TestGeometryAndWiring:
    """先证明载体是它自称的东西，再谈误差 —— 否则误差在量什么都不知道。"""

    def test_wall_nodes_lie_on_the_exact_circle(self, sweep):
        """★ 几何误差必须远低于要测的离散误差，否则收敛阶量的是样条不是求解器。"""
        for name, d in sweep.items():
            assert d["r_wall_dev"] < 1e-6, (
                f"{name}: 壁面节点到精确圆偏差 {d['r_wall_dev']:.3e} —— "
                "样条几何误差已进入可测范围，本门的收敛阶不再只反映离散误差")

    def test_prescribed_gamma_skips_the_kutta_loop(self, sweep):
        """★★ 圆柱没有尖尾缘 ⇒ Kutta 选不出 Γ。本门**规定** Γ，
        `solve_laplace_lifting` 必须真的跳过 Kutta 外循环（`n_kutta_updates == 0`）。
        若它跑了循环，那 Γ 就不是我们给的那个，整道门在量别的东西。"""
        for name, d in sweep.items():
            assert d["n_kutta"] == 0, (
                f"{name}: n_kutta_updates = {d['n_kutta']} != 0 —— "
                "gamma_fixed 没有短路 Kutta 环，Γ 不再是规定值")
            assert d["residual"] < 1e-8, f"{name}: |R| = {d['residual']:.2e} 未收敛"

    def test_sign_convention_is_the_measured_one(self, sweep):
        """★★★ Γ 的符号约定是**量出来**的，不是抄来的：反号时 max|ΔCp| 从 0.10
        跳到 1.60（16 倍）。本条把那个差距钉住，于是将来若约定被改，门会红而不是
        悄悄给出一个「看起来也还行」的误差。"""
        d = sweep["medium"]
        wrong = d["cp"] - (1.0 - (2.0 * np.sin(d["theta"])
                                  - GAMMA / (2.0 * np.pi * A * U_INF)) ** 2)
        ratio = float(np.max(np.abs(wrong))) / d["cp_max"]
        assert ratio > 5.0, (
            f"反号误差只有正号的 {ratio:.2f} 倍（实测 15.9）—— 两个符号已难以区分，"
            "说明 Γ 没有真正进入解，或者环量太小以致这道门测不到它")


class TestConvergence:
    """裁决③：C 类判**收敛阶** ∧ **绝对误差不过大**。两条都要。"""

    def test_cp_convergence_order(self, sweep):
        """★ 相邻两级都要达到阶 >= 1.0（实测 1.417 / 1.686）。

        ★★ **反向落点**（判据不可单侧）：若求解器退化，阶会掉向 0；
        若有人把 `h_far` 钳回固定值，实测阶塌到 **0.325** —— 本条正是为那种情形红的。
        """
        names = [n for n, _ in LEVELS]
        for i in range(1, len(names)):
            a, b = sweep[names[i - 1]], sweep[names[i]]
            for key, label in (("cp_max", "max|Cp|"), ("cp_rms", "rms Cp")):
                p = _order(a[key], b[key], a["h"], b["h"])
                assert p >= ORDER_MIN, (
                    f"{label} 阶 {names[i-1]} -> {names[i]} = {p:.3f} < {ORDER_MIN}"
                    f"（实测基线 1.417 / 1.686）。★ 最常见的原因不是求解器，是阶梯 ——"
                    " 检查 h_far 是否跟着 h_wall 一起细（钳住它实测塌到 0.325）")

    def test_absolute_error_at_medium(self, sweep):
        """★ 阶对了但绝对误差可以仍然很大 —— 裁决③要求两条都判。"""
        d = sweep["medium"]
        assert d["cp_max"] <= CP_MAX_MEDIUM, (
            f"medium max|ΔCp| = {d['cp_max']:.4f} > {CP_MAX_MEDIUM}（实测 0.1009）")
        assert d["cp_rms"] <= CP_RMS_MEDIUM, (
            f"medium rms ΔCp = {d['cp_rms']:.4f} > {CP_RMS_MEDIUM}（实测 0.0517）")

    def test_error_is_monotone_in_h(self, sweep):
        """★ 单调性是比阶更弱、但比阶更难作弊的条件：三级误差必须逐级下降。"""
        vals = [sweep[n]["cp_rms"] for n, _ in LEVELS]
        assert vals[0] > vals[1] > vals[2], f"rms ΔCp 非单调: {vals}"


class TestKuttaJoukowski:
    """★★★ 本门的核心：**给进 Γ，量出 cl**，两条完全不同的路。"""

    def test_lift_matches_the_prescribed_circulation(self, sweep):
        """cl_p（表面压力积分）对 2Γ/(U c)（Kutta–Joukowski）。

        ★ 这不是一个恒等式的重述：左边走的是壁面 Cp 恢复 + 面积分，右边只是我们
        输入的那个数。两者一致，说明**环量的输入、势场的求解、压力的恢复、力的积分**
        这四步整条链是自洽的。
        """
        d = sweep["medium"]
        assert d["cl_rel"] <= CL_REL_MEDIUM, (
            f"medium cl_p = {d['cl']:.6f} vs 精确 {d['cl_exact']:.6f}"
            f"（相对 {100*d['cl_rel']:.3f} % > {100*CL_REL_MEDIUM:.1f} %，实测 1.136 %）")

    def test_lift_error_contracts_with_refinement(self, sweep):
        """★★ 只判 medium 的绝对值会被一个**不收敛**的实现蒙混过去，所以还要判收缩。

        ★ 但收缩只从 **coarse -> medium** 判，因为 **xcoarse -> coarse 实测是平的**
        （3.393 % -> 3.477 %，阶 -0.035）—— 那一级只有 598 节点、h_far = 12 而域半径
        才 15，仍在前渐近区。**如实 RECORDED，不塞进判据**：把一个已知处于前渐近区的
        点纳入收敛判据，是本项目记过的「判据覆盖不了被比较量的定义域」。
        """
        c, m = sweep["coarse"], sweep["medium"]
        ratio = c["cl_rel"] / m["cl_rel"]
        assert ratio >= CL_CONTRACTION, (
            f"cl 误差 coarse -> medium 只收缩 {ratio:.2f}x < {CL_CONTRACTION}x"
            f"（实测 3.06x：{100*c['cl_rel']:.3f} % -> {100*m['cl_rel']:.3f} %）")

    def test_a_second_circulation_behaves_the_same(self):
        """★ 一个 Γ 值可能碰巧对上。换一个（0.6 -> 1.0）再看相对误差是否同量级。

        ★★ 判据是**相对**的（对各自的精确值），不是绝对的 —— 两个 Γ 的精确 cl 差 1.67 倍，
        用绝对阈值会让大 Γ 那腿自动更容易通过。这是四次判据缺陷里的第 3 条。
        """
        d = _one_level(0.02, gamma=1.0)
        assert d["cl_rel"] <= CL_REL_MEDIUM, (
            f"Γ = 1.0 时 cl_p = {d['cl']:.6f} vs 精确 {d['cl_exact']:.6f} "
            f"（相对 {100*d['cl_rel']:.3f} %）")


class TestCommittedEvidenceIsLoadBearing:
    r"""★★★ 新鲜计算 vs 已提交 `summary.csv` —— 把 `conftest` 的承诺变成事实。
    设计与三个坑见 `tests/_gate_evidence.py`（不能「值 == 它自己的来源」·
    判据阈值不进 CSV · 证据按真值精度存）。★ 零额外计算：用的是 `sweep` 已算好的结果。
    """

    def test_fresh_run_reproduces_the_committed_summary(self, sweep):
        fresh = {(name,): dict(cp_max=sweep[name]["cp_max"],
                               cp_rms=sweep[name]["cp_rms"],
                               cl_p=sweep[name]["cl"],
                               cl_rel=sweep[name]["cl_rel"])
                 for name, _ in LEVELS}
        n = assert_matches_committed(
            os.path.join(str(REPO_ROOT), "cases", "gates", "C06_lifting_cylinder"),
            fresh, ("cp_max", "cp_rms", "cl_p", "cl_rel"),
            key_of=lambda r: (r["level"],),
            refresh_hint="PYFP3D_GATE_FIGURES=1 pytest tests/C/test_C06_lifting_cylinder.py")
        assert n >= 12, f"只比了 {n} 个数（3 级 x 4 列 = 12）"


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_lifting_cylinder_figure(sweep, gate_evidence_dir):
    """把 Cp(θ) 对精确解、以及收敛曲线，写进被跟踪的证据目录。"""
    import matplotlib
    matplotlib.use("Agg")
    import csv

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    th = np.linspace(-np.pi, np.pi, 721)
    ax[0].plot(th, _exact_cp(th), "k-", lw=1.6, label="exact")
    for name, _ in LEVELS:
        d = sweep[name]
        o = np.argsort(d["theta"])
        ax[0].plot(d["theta"][o], d["cp"][o], ".", ms=3,
                   label=f"{name} (n={d['n_nodes']})")
    ax[0].set_xlabel("theta [rad]"), ax[0].set_ylabel("Cp")
    ax[0].set_title(f"C06 lifting cylinder, Gamma = {GAMMA}")
    ax[0].invert_yaxis(), ax[0].legend(fontsize=8), ax[0].grid(alpha=.3)

    hs = np.array([sweep[n]["h"] for n, _ in LEVELS])
    for key, mk, lb in (("cp_max", "o-", "max|dCp|"), ("cp_rms", "s-", "rms dCp"),
                        ("cl_rel", "^-", "cl rel err")):
        ax[1].loglog(hs, [sweep[n][key] for n, _ in LEVELS], mk, label=lb)
    ax[1].loglog(hs, 0.5 * hs ** 1.0 / hs[0] ** 1.0 * sweep["xcoarse"]["cp_max"],
                 "k--", lw=.8, label="order 1")
    ax[1].set_xlabel("h_wall"), ax[1].set_ylabel("error")
    ax[1].set_title("all scales refined together"), ax[1].legend(fontsize=8)
    ax[1].grid(alpha=.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "c06_lifting_cylinder.png"), dpi=130)
    plt.close(fig)

    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["level", "h_wall", "n_nodes", "n_wall", "cp_max", "cp_rms",
                    "cl_p", "cl_exact", "cl_rel", "residual", "n_kutta"])
        for name, _ in LEVELS:
            d = sweep[name]
            w.writerow([name, d["h"], d["n_nodes"], d["n_wall"],
                        #: ★ 实测列一律 `.9e`（`_gate_evidence.fmt`）—— 证据要按真值精度存，
                        #: 不是显示精度；理由见 `tests/_gate_evidence.py` 的第 3 条。
                        fmt(d["cp_max"]), fmt(d["cp_rms"]),
                        fmt(d["cl"]), fmt(d["cl_exact"]),
                        fmt(d["cl_rel"]), f"{d['residual']:.3e}", d["n_kutta"]])
