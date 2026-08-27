r"""D12 — RAE2822 Case 7 / Case 9 实验 Cp：**无粘模型的过升偏置**（D 类）。

★★★ 与 D11 同一形状、同一裁决：**对无粘用 Euler 设门、对实验记录偏置**（裁决三）。
本门的**主产物是偏置**；能设门的只有偏置的**方向与结构**。

★★ **数据在 2026-08-27 之前只被 `tests/A/test_A19_meshgen_rae2822.py` 读来建网格**，
没有任何比对 —— 与 D11 那一批同批发现，同一个包住 `open` 的运行时探针。

★★★ **实测（medium，`smooth_passes=1`，生产配方，实验攻角不修正）**：

| 工况 | 上表面 Cp RMS | 激波位置 计算 − 实验 | **cn**（法向力） | 对实验偏置 |
|---|---|---|---|---|
| **Case 7** M0.725 α2.55° | 0.3257 | **+0.1243 c** | 1.0174 vs 0.6576 | **+54.7 %** |
| **Case 9** M0.73 α3.19° | 0.2454 | **+0.1194 c** | 1.1150 vs 0.8077 | **+38.0 %** |

⇒ **无粘全速势把法向力高估 38–55 %，激波压后 ~0.12 c** —— 两者是同一件事的两面：
边界层位移把真实激波前移、并使翼型有效弯度下降，无粘模型两样都没有。
★ 与 Track V 的粘性耦合路径对照：GV5.2 在**有粘**路径上测到激波仍偏后 **0.06–0.10 c**
（`phases/p2/bench/studies/v5_2_rae2822/`）—— 本门的 0.12 c 是**无粘**的量级，两者不可混。

★★ **cn 而不是 cl，是为了同尺**：实验数据只有 (x, Cp)，**没有 y** ⇒ 轴向项不可得。
两边都只取 `cn = ∫(cp_lower − cp_upper) d(x/c)` 才是苹果对苹果。
（α 这个量级下轴向项约占 1 %，B06 在 α=6° 实测 1.01 %。）

★★★ **Case 9 的 coarse 腿是 `clamping`，而它与 D11 的失败方向恰好相反**：
13 limited / **2857 floored**、**M_max 撞上 `m_cap = 3.0`**（分类器：
`clamping, clamps 13/2857`）。按模式表，`clamping` 指向 **physics/geometry** ——
coarse 网格分辨不了 α=3.19° 的流动。
★ 而 D11 里是 **coarse 成功、medium 极限环**；这里是 **coarse 钳制、medium 干净**。
**两个方向都出现过 ⇒ "细网格更可靠"不是一条可以默认的假设**，每个工况都要自己量。

★ 六条腿全部带 `sigma-freeze: frozen_in_transient`，RECORDED（与 D11 同）。
★ 攻角说明：文件名给的 2.55° / 3.19° 与本项目 GV5.2 用的一致。RAE2822 的
AGARD 表里 Case 7 的名义攻角是 2.92°、2.55° 是风洞修正值 —— 也就是说这里的
"不修正"指的是**照抄数据文件给的角度、不再做二次修正**，而不是"名义角"。
"""
import os
import re

import numpy as np
import pytest

from bench.failure_modes import classify_failure
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.shock import shock_report
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting
from tests.conftest import REPO_ROOT, gate_figures_enabled

CASES = (
    ("Case7_M0.725_a2.55", 0.725, 2.55, "ExpCase7_RAE2822_M0.725_AoA2.55_Rec6.5e6.dat"),
    ("Case9_M0.73_a3.19", 0.73, 3.19, "Expe_RAE2822_M0.73_AoA3.19_Rec6.5e6.dat"),
)
LEVELS = ("coarse", "medium")
RECIPE = dict(upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80,
              n_picard_seed=5)

#: —— 判据（非实测值；实测基线见 docstring 的表）——
OVERLIFT_MAX = 0.80      # cn 相对偏置上界：实测 +0.526 / +0.363
SHOCK_AFT_MIN = 0.05     # 激波后移下界：实测 +0.1243 / +0.1194，远超分辨散布
CN_EXP = {"Case7_M0.725_a2.55": 0.6576, "Case9_M0.73_a3.19": 0.8077}


def _load_experiment(fn):
    """★ Case 7 带显式 `ZONE T = "Exp Upper"/"Exp Lower"`；Case 9 是单分区，按 x 转折切，
    并以「Cp 更负的那一段是上表面」定侧 —— 不靠出现顺序。"""
    z, cur = {}, None
    with open(os.path.join(str(REPO_ROOT), "cases", "reference_data",
                           "rae2822_experiment", fn), encoding="utf-8") as fh:
        for ln in fh:
            if ln.startswith(("#", "VARIABLES")):
                continue
            m = re.match(r'ZONE T = "(.*)"', ln.strip())
            if m:
                cur = m.group(1); z[cur] = []; continue
            p = ln.split()
            if len(p) == 2 and cur is not None:
                try:
                    z[cur].append((float(p[0]), float(p[1])))
                except ValueError:
                    pass
    z = {k: np.array(v) for k, v in z.items() if v}
    if len(z) == 2:
        u = [v for k, v in z.items() if "Upper" in k][0]
        l = [v for k, v in z.items() if "Lower" in k][0]
        return {"upper": (u[:, 0], u[:, 1]), "lower": (l[:, 0], l[:, 1])}
    v = list(z.values())[0]
    x, cp = v[:, 0], v[:, 1]
    dx = np.diff(x)
    k = int(np.nonzero(np.sign(dx[:-1]) != np.sign(dx[1:]))[0][0]) + 1
    a, b = (x[:k + 1], cp[:k + 1]), (x[k + 1:], cp[k + 1:])
    return {"upper": a, "lower": b} if a[1].min() < b[1].min() else {"upper": b, "lower": a}


def _cn(xu, cu, xl, cl):
    """`cn = ∫(cp_l − cp_u) d(x/c)`，两边同一把尺（实验无 y ⇒ 不含轴向项）。"""
    g = np.linspace(0.0, 1.0, 801)
    su, sl = np.argsort(xu), np.argsort(xl)
    return float(np.trapezoid(np.interp(g, np.asarray(xl)[sl], np.asarray(cl)[sl])
                              - np.interp(g, np.asarray(xu)[su], np.asarray(cu)[su]), g))


def _experimental_shock_x(exp):
    x, cp = exp["upper"]
    o = np.argsort(x)
    x, cp = x[o], cp[o]
    d = np.diff(cp)
    j = int(np.argmax(d))
    return 0.5 * (x[j] + x[j + 1]) if d[j] > 0.15 else None


def _one(level, m_inf, alpha, fn):
    mc, wc = cut_wake(read_mesh(os.path.join(
        str(REPO_ROOT), "cases", "meshes", "rae2822_2.5d", f"{level}.msh")))
    dz = float(np.ptp(mc.nodes[:, 2]))
    r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=alpha, **RECIPE)
    phi = np.asarray(r["phi"])
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
                                alpha_deg=alpha, u_inf=1.0, s_ref=dz, m_inf=m_inf)
    cur = section_cp_curve(mc, phi, z=float(np.mean(mc.nodes[:, 2])),
                           smooth_passes=1, m_inf=m_inf)
    exp = _load_experiment(fn)
    rms = {s: float(np.sqrt(np.mean(
        (np.interp(exp[s][0], cur[f"x_{s}"], cur[f"cp_{s}"]) - exp[s][1]) ** 2)))
        for s in ("upper", "lower")}
    h = np.asarray(r["residual_history"], float)
    #: ★ 分类器只用于**未收敛**腿（D11 的教训：收敛腿残差尾巴天然平，会被读成 limit_cycle）
    if r.get("converged"):
        mode = ("converged", f"|R|={h[-1]:.2e} in {len(h)} steps")
    else:
        mode = classify_failure(
            h, np.asarray(r.get("clamp_history", [0] * len(h)), float),
            np.asarray(r.get("F_history", h), float), r.get("n_gmres_stalled", 0),
            r.get("accept_reason"), r.get("n_limited", 0), r.get("n_floored", 0))
    return dict(
        converged=bool(r.get("converged")), mode=mode[0], mode_detail=mode[1],
        residual=float(h[-1]), n_newton=int(len(h)),
        n_limited=int(r.get("n_limited", 0)), n_floored=int(r.get("n_floored", 0)),
        mach_max=float(np.sqrt(r["mach2_max"])), cl=float(f["cl"]),
        cn=_cn(cur["x_upper"], cur["cp_upper"], cur["x_lower"], cur["cp_lower"]),
        cn_exp=_cn(*exp["upper"], *exp["lower"]),
        rms_upper=rms["upper"], rms_lower=rms["lower"],
        x_shock=shock_report(cur, m_inf)["upper"].get("x_shock"),
        x_shock_exp=_experimental_shock_x(exp), curve=cur, exp=exp)


@pytest.fixture(scope="module")
def runs():
    return {(n, lv): _one(lv, m, a, fn) for n, m, a, fn in CASES for lv in LEVELS}


class TestWhatIsGateable:
    r"""★ 裁决三：对实验只记录偏置 ⇒ 能设门的是**方向、结构与量级上界**。"""

    def test_the_experimental_normal_force_is_the_committed_one(self, runs):
        r"""★ 先锁尺子：从已提交的实验 Cp 积出的 `cn` 不许静默漂移。
        实测 Case 7 **0.6576**、Case 9 **0.8077**。"""
        for name, _, _, _ in CASES:
            got = runs[(name, "medium")]["cn_exp"]
            assert abs(got - CN_EXP[name]) < 5e-3, (
                f"{name} 的实验 cn 从 {CN_EXP[name]} 变成 {got:.4f} —— "
                "实验文件或积分口径变了，请按纪律 11 走勘误")

    def test_the_inviscid_solution_overpredicts_the_normal_force(self, runs):
        r"""★★★ **方向是可判的，而且有物理内容**：边界层位移使翼型**有效弯度下降**，
        所以无粘解必然**高估**法向力。实测 medium：Case 7 **+54.7 %**、Case 9 **+38.0 %**。

        ★ 只判**方向**与一个宽松上界（≤ 80 %）—— 量级本身是 RECORDED，因为它度量的是
        我们刻意不建模的物理。**反向落点**：若无粘解开始**低估**，那不是精度问题，是机理反了。
        """
        for name, _, _, _ in CASES:
            d = runs[(name, "medium")]
            assert d["converged"] and d["n_limited"] == 0 and d["n_floored"] == 0, (
                f"{name} medium 不干净：{d['mode']} ({d['mode_detail']})")
            rel = d["cn"] / d["cn_exp"] - 1.0
            assert rel > 0, (
                f"{name}: cn {d['cn']:.4f} 未高于实验 {d['cn_exp']:.4f}"
                f"（{100*rel:+.1f} %）—— 边界层减弯的机理反了")
            assert rel <= OVERLIFT_MAX, (
                f"{name}: 过升 {100*rel:+.1f} % 超过 {100*OVERLIFT_MAX:.0f} %"
                f"（实测 +54.7 % / +38.0 %）")

    def test_the_inviscid_shock_sits_well_downstream(self, runs):
        r"""★★ 激波后移，**两个工况都远超分辨散布**（实测 +0.1243 / +0.1194 c）——
        与 D11 里 M0.803 那种 ±0.015 跨零的情形不同，这里可以设门。

        ★ 对照：GV5.2 在**有粘耦合**路径上测到仍偏后 **0.06–0.10 c**。本门的 0.12 c
        是**无粘**量级，两者不可混为一谈。
        """
        for name, _, _, _ in CASES:
            d = runs[(name, "medium")]
            assert d["x_shock"] is not None and d["x_shock_exp"] is not None, name
            off = d["x_shock"] - d["x_shock_exp"]
            assert off >= SHOCK_AFT_MIN, (
                f"{name}: 激波偏置 {off:+.4f} c < {SHOCK_AFT_MIN}"
                f"（实测 +0.1243 / +0.1194）")


class TestWhatIsOnlyRecorded:
    r"""★ 把"只记录"的事实**锁成断言**，这样它们一旦改变会响而不是被悄悄重读。"""

    def test_the_bias_grows_with_refinement(self, runs):
        r"""★★ Case 7 的过升**随加密变大**：+39.5 % → **+54.7 %**（cn 0.9170 → 1.0174，
        而实验 0.6576）。★ 这不是"加密使解变差"，是**加密让无粘解更接近它自己的极限**，
        而那个极限离实验更远 —— 缺的物理不会因为网格变细而出现。
        ★★ 本条断言这个事实仍然成立。若它反过来（加密使偏置变小），那说明粗网格的
        误差此前在**抵消**模型误差，那是一个必须记录的巧合。
        """
        c, m = runs[("Case7_M0.725_a2.55", "coarse")], runs[("Case7_M0.725_a2.55", "medium")]
        assert c["converged"] and m["converged"]
        rc, rm = c["cn"] / c["cn_exp"] - 1.0, m["cn"] / m["cn_exp"] - 1.0
        assert rm > rc, (
            f"Case 7 的过升不再随加密变大（{100*rc:+.1f} % → {100*rm:+.1f} %）—— "
            "粗网格误差可能在抵消模型误差，请记录")

    def test_the_coarse_case9_leg_is_clamping_and_is_not_quoted(self, runs):
        r"""★★★ Case 9 的 **coarse** 腿是 `clamping`：**13 limited / 2857 floored**，
        **M_max 撞上 `m_cap = 3.0`**。按模式表，`clamping` 指向 **physics/geometry** ——
        coarse 网格分辨不了 α=3.19° 的流动，它的每一个数都不作数。

        ★★★ **而它与 D11 的失败方向恰好相反**：D11 里是 **coarse 成功、medium 极限环**，
        这里是 **coarse 钳制、medium 干净**。**两个方向都出现过 ⇒「细网格更可靠」
        不是一条可以默认的假设**，每个工况都要自己量。本条把这个反转锁住。
        """
        c = runs[("Case9_M0.73_a3.19", "coarse")]
        m = runs[("Case9_M0.73_a3.19", "medium")]
        assert not c["converged"] and c["mode"] == "clamping", (
            f"Case 9 coarse 现在是 {c['mode']}（{c['mode_detail']}）而不是 clamping —— "
            "失败模式变了，修法也就不同（见 CLAUDE.md 的模式表）")
        assert c["n_floored"] > 500 and c["mach_max"] >= 2.99, (
            f"钳制规模变了：{c['n_limited']}/{c['n_floored']}，M_max {c['mach_max']:.4f}"
            "（实测 13/2857，M_max 撞 m_cap=3.0）")
        assert m["converged"] and m["n_limited"] == 0 and m["n_floored"] == 0, (
            "Case 9 medium 不再干净 —— 那个反转（coarse 坏 / medium 好）没了")

    def test_the_two_cases_have_comparable_bias(self, runs):
        r"""★ 两个工况的过升在同一量级（+54.7 % / +38.0 %，比值 1.44）——
        RECORDED：它说明这个偏置是**模型的性质**，不是某一个工况的巧合。"""
        r7 = runs[("Case7_M0.725_a2.55", "medium")]
        r9 = runs[("Case9_M0.73_a3.19", "medium")]
        a = r7["cn"] / r7["cn_exp"] - 1.0
        b = r9["cn"] / r9["cn_exp"] - 1.0
        assert 1.0 / 3.0 < a / b < 3.0, (
            f"两工况的过升不再同量级（{100*a:+.1f} % vs {100*b:+.1f} %，比 {a/b:.2f}）")


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_rae2822_bias_figure(runs, gate_evidence_dir):
    import csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    for k, (name, m_inf, alpha, _) in enumerate(CASES):
        d = runs[(name, "medium")]
        for side, mk in (("upper", "-"), ("lower", "--")):
            xe, ce = d["exp"][side]
            ax[k].plot(xe, ce, "ko", ms=3, label="experiment" if side == "upper" else None)
            ax[k].plot(d["curve"][f"x_{side}"], d["curve"][f"cp_{side}"], mk, lw=1.2,
                       color="C0", label="inviscid FP (medium)" if side == "upper" else None)
        if d["x_shock_exp"] is not None:
            ax[k].axvline(d["x_shock_exp"], color="k", ls=":", lw=1)
        if d["x_shock"] is not None:
            ax[k].axvline(d["x_shock"], color="C0", ls=":", lw=1)
        rel = 100 * (d["cn"] / d["cn_exp"] - 1.0)
        ax[k].set_title(f"RAE2822  M {m_inf}  alpha {alpha}\n"
                        f"cn {d['cn']:.4f} vs exp {d['cn_exp']:.4f}  ({rel:+.1f} %)"
                        f"   shock {d['x_shock'] - d['x_shock_exp']:+.4f} c", fontsize=9)
        ax[k].set_xlabel("x/c"), ax[k].invert_yaxis(), ax[k].grid(alpha=.3)
        ax[k].legend(fontsize=7)
    ax[0].set_ylabel("Cp")
    fig.suptitle("D12  inviscid over-lift against experiment -- the BIAS is the product",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "d12_rae2822_bias.png"), dpi=130)
    plt.close(fig)

    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "level", "converged", "mode", "residual", "n_newton",
                    "n_limited", "n_floored", "mach_max", "cl_p", "cn", "cn_exp",
                    "cn_rel", "rms_upper", "rms_lower", "x_shock", "x_shock_exp",
                    "shock_offset"])
        for name, _, _, _ in CASES:
            for lv in LEVELS:
                d = runs[(name, lv)]
                off = ("" if d["x_shock"] is None or d["x_shock_exp"] is None
                       else f"{d['x_shock'] - d['x_shock_exp']:+.4f}")
                w.writerow([name, lv, d["converged"], d["mode"], f"{d['residual']:.3e}",
                            d["n_newton"], d["n_limited"], d["n_floored"],
                            f"{d['mach_max']:.4f}", f"{d['cl']:.6f}", f"{d['cn']:.4f}",
                            f"{d['cn_exp']:.4f}", f"{d['cn']/d['cn_exp']-1.0:+.4f}",
                            f"{d['rms_upper']:.4f}", f"{d['rms_lower']:.4f}",
                            "" if d["x_shock"] is None else f"{d['x_shock']:.4f}",
                            "" if d["x_shock_exp"] is None else f"{d['x_shock_exp']:.4f}",
                            off])
