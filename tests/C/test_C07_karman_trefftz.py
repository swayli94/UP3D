r"""C07 — Kármán–Trefftz：**唯一能把 Kutta 条件对精确 Γ 检验**的门（C 类）。

精确解在 `tests/_kt_case.py`（那份文件的 docstring 记着验它用的**四个独立 oracle**）。

★★★ **C06 与 C07 的分工**：圆柱**没有尖尾缘** ⇒ Kutta 选不出环量，C06 里 Γ 是
**规定**的（`gamma_fixed`），检验的是「给进 Γ，量出 cl」。KT 有**有限尾缘角** ⇒
Kutta **必须自己选出** Γ，而精确值有解析式 ⇒ 本门检验的是**选择本身**。

★★★ **本轮最值钱的读数：积分量收敛，前缘点值不收敛 —— 而且是在同一个解上。**
（τ=10°、α=5°、全尺度加密，实测 2026-08-26）

| 量 | h 0.08→0.04 | h 0.04→0.02 | medium 绝对值 |
|---|---|---|---|
| **Γ**（Kutta 选出的） | **2.109** | **1.493** | 2.803 % |
| **cl**（表面压力积分） | **1.957** | **1.468** | 3.045 % |
| Cp **MID** 带 (0.05–0.90) | **2.893** | **1.761** | 0.0108 |
| Cp **TE** 带 (0.90–1.01) | 0.532 | 0.645 | 0.1778 |
| Cp **LE** 带 (x<0.05) | ★ **无样本** | **0.085** | 0.1610 |

★★ **两个慢带都有已定位的成因，都是能力陈述而不是缺陷**：
· **LE** —— 这是本项目 G1.6 那一族（球上实测 LE 带阶 0.37/0.87、「medium 的 11.6 %
≈ h=0.08 处 P1 场的固有 max-norm 能力」），**现在第一次在一个有精确解的二维翼型上复现**。
排除过两个混淆项：把折线几何点数 ×20（100→2000）只让 LE 阶从 0.241 动到 0.287；
把 oracle 采样加密 5×（匹配距 3.5e-05→6e-06）只让 LE RMS 从 0.1594 动到 0.1574。
· **TE** —— 有限楔角本身：精确解的 \|q\| 以 (ζ-b)^(τ/π) 趋零，τ=10° 时指数只有 **0.0556**，
θ=1e-5 处精确 \|q\| 还有 0.49 ⇒ **精确解自己就没有干净的驻点**，数值解自然更没有。

★★★ **所以 LE / TE 两带 RECORDED、不设判据，而且 LE 有一个更硬的理由**：
h=0.08 时 x<0.05 内**一个壁面节点都没有**（RMS 为 nan）—— **自变量改变了样本是否存在**，
正是四条判据缺陷里的第 2 条。在一个会随 h 出现/消失的样本集上比较「谱宽」是没有定义的。
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
from tests._kt_case import KarmanTrefftz
from tests.conftest import gate_figures_enabled

LEVELS = (("xcoarse", 0.08), ("coarse", 0.04), ("medium", 0.02))
BANDS = (("LE", 0.0, 0.05), ("MID", 0.05, 0.90), ("TE", 0.90, 1.01))

#: 判据（非实测值；实测基线写在 docstring 的表里）
ORDER_MIN = 1.20        # Γ / cl / MID：实测最小 1.468（cl，0.04->0.02），余量 22 %
REL_MAX_MEDIUM = 0.035  # Γ 2.803 % / cl 3.045 %
#: ★ cl 对这条判据的余量只有 **15 %** —— 如实记下来。它是本门最紧的一条，
#: 而解是确定性的 Laplace 解（无线程/随机来源），所以余量小不等于会抖。
MID_ABS_MEDIUM = 0.03   # 实测 0.0108
#: Γ 与 cl 两条独立路线的误差必须彼此吻合：实测 2.72 % vs 2.86 %，相差 5 %
ROUTES_AGREE = 0.25
#: MID 阶必须显著高于 LE 阶 —— 这条**把本轮的发现本身变成门**：实测 1.704 / 0.118 = 14.4x
MID_OVER_LE = 5.0


@pytest.fixture(scope="module")
def kt():
    return KarmanTrefftz()


def _one_level(kt, h):
    p2, tri, eg, ig = airfoil_wake_2d(
        kt.polyline(400), model_name="c07_kt", r_far=15.0, h_wall=h,
        #: ★★ 全尺度加密 —— C06 实测把 h_far 钳成固定值会让阶塌到 0.325（P11 的体网格污染地板）
        h_far=150.0 * h, h_wake=3.0 * h,
        dist_min=0.1, dist_max=6.0, wake_dist_max=1.5)
    dz = 2.0 * h
    mesh = extrude_single_layer(p2, tri, eg, ig, dz=dz, name="c07")
    mc, wc = cut_wake(mesh)
    #: ★ 不传 gamma_fixed ⇒ Kutta 自己选 Γ，这是本门的全部意义
    r = solve_laplace_lifting(mc, wc, alpha_deg=np.rad2deg(kt.alpha), u_inf=kt.u)
    phi = np.asarray(r["phi"])
    wf = mc.boundary_faces["wall"]
    wall = np.unique(wf)
    gamma = float(np.atleast_1d(r["gamma"])[0])
    f = wall_force_coefficients(mc.nodes, mc.elements, wf, phi,
                                alpha_deg=np.rad2deg(kt.alpha), u_inf=kt.u,
                                s_ref=dz, m_inf=0.0)
    grad = wall_tangential_gradient_quadratic(mc.nodes, wf, phi)
    cp = 1.0 - np.sum(grad[wall] ** 2, axis=1) / kt.u ** 2
    cp_ex, dist = kt.cp_at(mc.nodes[wall, 0] + 1j * mc.nodes[wall, 1])
    err = cp - cp_ex
    x = mc.nodes[wall, 0]
    bands, counts = {}, {}
    for name, a, b in BANDS:
        m = (x >= a) & (x < b)
        counts[name] = int(m.sum())
        bands[name] = float(np.sqrt(np.mean(err[m] ** 2))) if m.any() else float("nan")
    return dict(
        h=h, n_nodes=int(mc.nodes.shape[0]), n_wall=int(wall.size),
        gamma=gamma, gamma_rel=abs(gamma - kt.gamma) / abs(kt.gamma),
        cl=float(f["cl"]), cl_rel=abs(f["cl"] - kt.cl) / abs(kt.cl),
        cp_rms=float(np.sqrt(np.mean(err ** 2))), bands=bands, counts=counts,
        match_max=float(dist.max()), n_kutta=int(r["n_kutta_updates"]),
        kutta_converged=bool(r["kutta_converged"]),
        x=x, cp=cp, cp_exact=cp_ex)


@pytest.fixture(scope="module")
def sweep(kt):
    return {name: _one_level(kt, h) for name, h in LEVELS}


def _order(a, b, ha, hb):
    return float(np.log(a / b) / np.log(ha / hb))


class TestExactSolutionOracles:
    """★★★ 先证明精确解配得上「精确」，而且**不用它自己的公式来证**。"""

    def test_circulation_from_a_contour_integral(self, kt):
        r"""oracle ①：沿**物理平面**翼型围线数值积分速度，对解析 Γ。

        ★ 这条与 `gamma_raw = 4πUR sin(α+β)` 那行**共享零行代码** —— 它走的是
        映射后的速度场。实测 \|Γ\| 相对差 4.1e-07（4e4 点围线求和的离散误差）。
        """
        from tests._kt_case import _dzdzeta
        th = (np.arange(kt.z.size) + 0.5) * (2.0 * np.pi / kt.z.size)
        #: 复速度 = u - i v
        wv = kt._w_prime(kt._circle(th)) / _dzdzeta(kt._circle(th), kt.b, kt.n)
        dz = 0.5 * (np.roll(kt.z, -1) - np.roll(kt.z, 1)) * kt.chord
        circ = float(np.sum(wv.real * dz.real + (-wv.imag) * dz.imag))
        rel = abs(abs(circ) - abs(kt.gamma_raw)) / abs(kt.gamma_raw)
        assert rel < 1e-5, (
            f"围线环流 {circ:.8f} 与解析 Γ {kt.gamma_raw:.8f} 不符（相对 {rel:.2e}，"
            "实测 4.1e-07）—— 映射、复势或 Kutta 公式三者之一错了")

    def test_trailing_edge_angle_is_the_tau_parameter(self, kt):
        """oracle ②：从 z(θ) 的**几何**切线量尾缘夹角，对参数 τ。实测 10.207° vs 10.000°。"""
        i = int(np.argmax(kt.z.real))
        k = max(3, kt.z.size // 400)
        m = kt.z.size
        a_up = np.angle(kt.z[(i + k) % m] - kt.z[i])
        a_lo = np.angle(kt.z[(i - k) % m] - kt.z[i])
        inc = np.rad2deg(abs(np.angle(np.exp(1j * (a_up - a_lo)))))
        assert abs(inc - kt.tau_deg) < 1.0, (
            f"几何量得的尾缘夹角 {inc:.4f}° 与参数 τ = {kt.tau_deg}° 不符（实测 10.207°）")

    def test_kutta_puts_a_stagnation_point_exactly_at_the_te_preimage(self, kt):
        """oracle ③：Kutta 条件 = 圆平面复速度在 **ζ = b** 处为零。

        ★★ **判据在 ζ = b 上解析求值，不去找采样点的最小值** —— 第一版正是那么写的，
        于是它红了：翼型有**两个**驻点（实测前驻点 θ ≈ 190.00°、后驻点 θ ≈ 359.955°），
        而离散 `argmin` 取哪一个**纯看半格偏移离谁更近**。断言「全局最小落在尾缘」
        是**判据没有覆盖被比较量的定义域** —— 本项目记过的那一族，这次落在
        「一个集合里的某个元素」被写成了「那个元素」。
        """
        w_te = kt._w_prime(complex(kt.b, 0.0))
        scale = kt.u
        assert abs(w_te) / scale < 1e-12, (
            f"圆平面复速度在 ζ = b 处不为零（|dw/dζ| = {abs(w_te):.3e}）—— "
            "Kutta 条件没有被精确解满足，本门的 Γ 参照失效")
        #: ★ 并且它必须是**后**驻点：前驻点在别处，两者不可混
        th = (np.arange(kt.z.size) + 0.5) * (2.0 * np.pi / kt.z.size)
        q = np.abs(kt._w_prime(kt._circle(th)))
        two = np.rad2deg(np.sort(th[np.argsort(q)[:2]]))
        near_te = min(min(t, 360.0 - t) for t in two)
        assert near_te < 0.5, (
            f"两个驻点在 θ = {two} °，没有一个落在尾缘原像（θ = 0）附近"
            f"（最近 {near_te:.3f}°，实测 0.045°）")

    def test_te_speed_decays_at_the_map_exponent(self, kt):
        r"""★★ oracle ④，四个里最强的一个：θ→0 时 \|q\| 每十倍的比值必须等于
        **0.1^(τ/π)** —— 它把观察到的尾缘行为**通过映射的指数**接回 τ 参数。

        实测 0.879，理论 0.1^(10/180) = 0.8797。
        ★ 配套对照写在载体的 docstring 里：**τ = 0（Joukowski 尖点）时 \|q\| → 0.9056
        常数非零**，而 τ > 0 时 → 0。两个情形干净分开是实现正确最有力的证据。
        """
        q = kt.surface_speed_at(np.array([1e-3, 1e-4, 1e-5]))
        ratios = q[1:] / q[:-1]
        want = 0.1 ** (kt.tau_deg / 180.0)
        assert np.allclose(ratios, want, rtol=5e-3), (
            f"尾缘速度每十倍衰减比 {ratios} 与映射指数给的 {want:.6f} 不符 "
            "—— 尾缘的解析行为与 τ 参数脱钩了")


class TestKuttaSelectsTheExactCirculation:
    """★★★ 本门的核心：Kutta **自己选出**的 Γ，对解析 Γ。"""

    def test_the_kutta_loop_actually_ran(self, sweep):
        """★ 前提：本门**不传** `gamma_fixed`，所以 Kutta 外循环必须真的跑过。
        若它没跑（`n_kutta_updates == 0`），那 Γ 是初值而不是被选出来的，
        整道门在量别的东西 —— 这与 C06 的前提**恰好相反**，两处对照着读。
        """
        for name, _ in LEVELS:
            d = sweep[name]
            assert d["n_kutta"] > 0 and d["kutta_converged"], (
                f"{name}: n_kutta_updates = {d['n_kutta']}, "
                f"converged = {d['kutta_converged']} —— Γ 不是 Kutta 选出来的")

    def test_gamma_converges_to_the_exact_value(self, sweep, kt):
        """阶 + 绝对误差（裁决③两条都判）。实测阶 2.109 / 1.536，medium 2.72 %。"""
        names = [n for n, _ in LEVELS]
        for i in range(1, len(names)):
            a, b = sweep[names[i - 1]], sweep[names[i]]
            p = _order(a["gamma_rel"], b["gamma_rel"], a["h"], b["h"])
            assert p >= ORDER_MIN, (
                f"Γ 误差阶 {names[i-1]} -> {names[i]} = {p:.3f} < {ORDER_MIN}"
                f"（实测 2.109 / 1.493）")
        d = sweep["medium"]
        assert d["gamma_rel"] <= REL_MAX_MEDIUM, (
            f"medium Γ = {d['gamma']:.6f} vs 精确 {kt.gamma:.6f}"
            f"（{100*d['gamma_rel']:.3f} % > {100*REL_MAX_MEDIUM:.1f} %）")

    def test_lift_converges_and_agrees_with_the_circulation_route(self, sweep, kt):
        """★★ cl 走的是**表面 Cp 恢复 + 面积分**，Γ 走的是 Kutta 行 —— 两条独立路线。
        它们的误差必须彼此吻合，否则「两个数都对」可能是两个不同的巧合。
        实测 2.803 % vs 3.045 %，相差 8 %。
        """
        names = [n for n, _ in LEVELS]
        for i in range(1, len(names)):
            a, b = sweep[names[i - 1]], sweep[names[i]]
            p = _order(a["cl_rel"], b["cl_rel"], a["h"], b["h"])
            assert p >= ORDER_MIN, (
                f"cl 误差阶 {names[i-1]} -> {names[i]} = {p:.3f} < {ORDER_MIN}"
                f"（实测 1.957 / 1.468）")
        d = sweep["medium"]
        assert d["cl_rel"] <= REL_MAX_MEDIUM, (
            f"medium cl = {d['cl']:.6f} vs 精确 {kt.cl:.6f}"
            f"（{100*d['cl_rel']:.3f} %）")
        spread = abs(d["cl_rel"] - d["gamma_rel"]) / max(d["cl_rel"], d["gamma_rel"])
        assert spread <= ROUTES_AGREE, (
            f"两条路线的误差不吻合：Γ {100*d['gamma_rel']:.3f} % vs "
            f"cl {100*d['cl_rel']:.3f} %（相差 {100*spread:.1f} % > "
            f"{100*ROUTES_AGREE:.0f} %，实测 8 %）")


class TestPointwiseVersusIntegrated:
    """★★★ 把本轮的发现本身做成门：**积分量收敛，前缘点值不收敛**。"""

    def test_mid_band_converges(self, sweep):
        """中段 Cp（0.05–0.90 c）：实测阶 2.893 / 1.704，medium RMS 0.0112。"""
        names = [n for n, _ in LEVELS]
        for i in range(1, len(names)):
            a, b = sweep[names[i - 1]], sweep[names[i]]
            p = _order(a["bands"]["MID"], b["bands"]["MID"], a["h"], b["h"])
            assert p >= ORDER_MIN, (
                f"MID 带阶 {names[i-1]} -> {names[i]} = {p:.3f} < {ORDER_MIN}"
                f"（实测 2.893 / 1.761）")
        assert sweep["medium"]["bands"]["MID"] <= MID_ABS_MEDIUM, (
            f"medium MID 带 RMS {sweep['medium']['bands']['MID']:.4f} > {MID_ABS_MEDIUM}")

    def test_le_band_is_the_slow_one_and_that_is_the_finding(self, sweep):
        """★★ LE 带的阶必须**显著低于** MID 带 —— 这不是在容忍缺陷，是在**锁住一个结论**：
        G1.6 那一族（P1 场在前缘吸力峰上的固有能力）在一个**有精确解的二维翼型**上复现。

        实测 MID 1.761 / LE 0.085 = **20.6×**。
        ★★ **反向落点**：若这条红了（LE 追上 MID），那是**好消息**，但意味着 G1.6 的
        能力陈述变了 —— 请按纪律 11 走勘误，`PROJECT_STRUCTURE.md` 的「Known gaps」
        与 P11 的再归因都以它为前提。
        ★ 排除过两个混淆项：折线几何点数 ×20 只让 LE 阶 0.241→0.287；oracle 采样
        加密 5× 只让 LE RMS 0.1594→0.1574。
        """
        c, m = sweep["coarse"], sweep["medium"]
        p_mid = _order(c["bands"]["MID"], m["bands"]["MID"], c["h"], m["h"])
        p_le = _order(c["bands"]["LE"], m["bands"]["LE"], c["h"], m["h"])
        assert p_mid / max(p_le, 1e-3) >= MID_OVER_LE, (
            f"MID 阶 {p_mid:.3f} 不再显著高于 LE 阶 {p_le:.3f}"
            f"（比值 {p_mid/max(p_le,1e-3):.1f} < {MID_OVER_LE}，实测 20.6）。"
            "★ 若 LE 真的追上来了这是好消息，但 G1.6 的能力陈述变了，请走纪律 11 勘误")

    def test_the_le_band_sample_set_is_not_comparable_at_the_coarsest_level(self, sweep):
        """★★★ 为什么 LE 带只判「比 MID 慢」而不判它自己的阶：**xcoarse 上
        x < 0.05 内一个壁面节点都没有**（RMS 为 nan）。

        **自变量改变了样本是否存在** —— 四条判据缺陷里的第 2 条。在一个会随 h
        出现/消失的样本集上比较谱宽是**没有定义的**，所以本条把这个事实**断言下来**，
        免得将来有人看到 LE 那一列就顺手补一个三级的阶。
        """
        assert sweep["xcoarse"]["counts"]["LE"] == 0, (
            f"xcoarse 的 LE 带现在有 {sweep['xcoarse']['counts']['LE']} 个节点 —— "
            "样本集变了，可以重新考虑给 LE 带一个三级判据（并更新 docstring 的表）")
        for name in ("coarse", "medium"):
            assert sweep[name]["counts"]["LE"] >= 3, (
                f"{name} 的 LE 带只有 {sweep[name]['counts']['LE']} 个节点，"
                "RMS 不具代表性")


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_kt_figure(sweep, kt, gate_evidence_dir):
    import csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    o = np.argsort(kt.z.real)
    ax[0].plot(kt.z.real, kt.cp, "k.", ms=0.6, label="exact (conformal)")
    for name, _ in LEVELS:
        d = sweep[name]
        ax[0].plot(d["x"], d["cp"], ".", ms=3, label=f"{name} (n={d['n_nodes']})")
    ax[0].set_xlabel("x/c"), ax[0].set_ylabel("Cp"), ax[0].invert_yaxis()
    ax[0].set_title(f"C07 Karman-Trefftz  tau={kt.tau_deg}deg  alpha=5deg")
    ax[0].legend(fontsize=8), ax[0].grid(alpha=.3)

    hs = np.array([sweep[n]["h"] for n, _ in LEVELS])
    for key, mk, lb in (("gamma_rel", "o-", "Gamma rel err"),
                        ("cl_rel", "^-", "cl rel err")):
        ax[1].loglog(hs, [sweep[n][key] for n, _ in LEVELS], mk, label=lb)
    for bn, mk in (("MID", "s-"), ("TE", "d-"), ("LE", "v-")):
        v = [sweep[n]["bands"][bn] for n, _ in LEVELS]
        ax[1].loglog(hs, v, mk, ms=4, label=f"Cp {bn} band")
    ax[1].loglog(hs, sweep["xcoarse"]["gamma_rel"] * hs / hs[0], "k--", lw=.8,
                 label="order 1")
    ax[1].set_xlabel("h_wall"), ax[1].set_ylabel("error")
    ax[1].set_title("integrated converges; LE band does not")
    ax[1].legend(fontsize=7), ax[1].grid(alpha=.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "c07_karman_trefftz.png"), dpi=130)
    plt.close(fig)

    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["level", "h_wall", "n_nodes", "n_wall", "n_kutta",
                    "gamma", "gamma_exact", "gamma_rel", "cl", "cl_exact", "cl_rel",
                    "cp_rms", "band_LE", "band_MID", "band_TE",
                    "n_LE", "n_MID", "n_TE", "match_max"])
        for name, _ in LEVELS:
            d = sweep[name]
            w.writerow([name, d["h"], d["n_nodes"], d["n_wall"], d["n_kutta"],
                        f"{d['gamma']:.6f}", f"{kt.gamma:.6f}", f"{d['gamma_rel']:.6e}",
                        f"{d['cl']:.6f}", f"{kt.cl:.6f}", f"{d['cl_rel']:.6e}",
                        f"{d['cp_rms']:.6f}",
                        *[f"{d['bands'][b]:.6f}" for b, _, _ in BANDS],
                        *[d["counts"][b] for b, _, _ in BANDS],
                        f"{d['match_max']:.3e}"])
