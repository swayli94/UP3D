r"""D11 — NACA0012 三工况实验 Cp：**记录无粘模型对实验的偏置**（D 类）。

★★★ **为什么这道门"主要不设门"，而这是裁决三的直接后果。**
裁决三（2026-08-24）：**对无粘用 Euler 设门 · 对有粘耦合用 RANS 设门 · 对实验记录偏置**。
本求解器是**无粘全速势**，所以拿它对实验做 pass/fail 会把**我们刻意不建模的物理**
（边界层位移、激波-边界层干涉、分离）算成缺陷。⇒ 本门的**主产物是偏置本身**。

★★ **数据在 2026-08-27 之前从未被任何代码读过** —— 提交了却没有读者。发现它的不是 grep
（grep 只看得见"提到"），是一个包住 `open` 跑全套的**运行时探针**。

★★★ **三个工况按"粘性物理有多重要"干净分层，这就是本门的读数**（medium，`smooth_passes=1`）：

| 工况 | 上表面 Cp RMS | 激波位置 计算 − 实验 | cl_p | 备注 |
|---|---|---|---|---|
| **M0.803 α = −0.10°**（近零升） | **0.0779** | **+0.0129 c** | −0.0423 | 无粘几乎够用 |
| **M0.778 α = 2.03°** | 0.3944 | **+0.2281 c** | 0.7062 | ★ 该腿是 **limit_cycle**，不作数 |
| 同上，coarse（收敛） | 0.2306 | **+0.1219 c** | 0.5751 | 经典激波后移偏置 |
| **M0.352 α = 12.86°**（近失速） | **1.5235** | 无真实激波 | 1.6764 | 无粘不可用 |

★★★ **激波前 / 激波后的分离，正好印证裁决一把 (b) 设门、(d) 只记录的形状**
（用**实验**激波位置切，两边同一把尺）：

| 工况 | 激波**前** RMS | 激波**后** RMS |
|---|---|---|
| M0.803 coarse → medium | **0.0926 → 0.0369**（**改善**） | 0.0609 → **0.0807**（**恶化**） |
| M0.778 coarse → medium | 0.2081 → 0.1994 | 0.1943 → **0.4777** |

⇒ **激波前是无粘模型能负责的地方，激波后不是。** 本门因此只在**激波前**与**近零升**
那一格上设门，其余 RECORDED。

★★ **M0.778 medium 不是"conv=False"，是 `limit_cycle`**（`bench/failure_modes.classify_failure`
判定，`descent10 = 0.599`，残差尾巴 **2.29e-07 → 3.82e-07 → 2.61e-07 周期 3**，三位有效数字
重复）。纪律要求**分类而不是报一个 conv=False**，而分类的后果是：**那一腿的数字不被任何
判据引用**，只作 RECORDED。

★ 全部六腿都带 `sigma-freeze WARNING: frozen_in_transient` —— 熵修正冻结在一个场已经离开的
状态上。RECORDED：本门不因此拒绝读数（GS1.4 的钳制契约管的是钳制），但它是解释
"激波后为何恶化"时必须一并摆上的一项。
"""
import os

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

#: 实验工况**逐字照抄**文件名，且按裁决二使用**实验攻角、不修正**
CASES = (
    ("M0.803_a-0.10", 0.803, -0.10, "Expe_NACA0012_M0.803_AoA-0.1_Rec6.5e6.dat"),
    ("M0.778_a2.03", 0.778, 2.03, "Expe_NACA0012_M0.778_AoA2.03_Rec6e6.dat"),
    ("M0.352_a12.86", 0.352, 12.86, "Expe_NACA0012_M0.352_Rec3e6_AoA12.86.dat"),
)
LEVELS = ("coarse", "medium")
#: 生产配方，与 `bench/run_m1_gate.py` 的调用逐字一致
RECIPE = dict(upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80,
              n_picard_seed=5)

#: —— 判据（非实测值；实测基线见 docstring 的两张表）——
PRE_SHOCK_MAX = 0.06        # M0.803 medium 激波前 RMS：实测 0.0369
PRE_SHOCK_CONTRACT = 1.5    # coarse -> medium 收缩：实测 0.0926/0.0369 = 2.51x
CL_NEAR_ZERO = 0.10         # 对称翼型 α = −0.10°：实测 |cl_p| = 0.0423
STRATIFY = 1.5              # 相邻层级的 RMS 比：实测 0.0633 -> 0.3986 -> 1.5477
SHOCK_AFT_MIN = 0.0         # 无粘激波必须在实验激波**下游**：实测 +0.013 / +0.122


def _load_experiment(fn):
    """★ 数据只有 (X, Cp)，没有上下表面标志 —— 按 x 的单调转折切。
    首段 Cp 强负 ⇒ 上表面（α > 0 时上表面是吸力面；实测三例皆然）。
    ★ M0.803 那一份**只有单面 24 点** ⇒ 它**算不出 cl**，这是数据本身的限制。
    """
    d = np.loadtxt(os.path.join(str(REPO_ROOT), "cases", "reference_data",
                                "naca0012_experiment", fn), skiprows=2)
    x, cp = d[:, 0], d[:, 1]
    dx = np.diff(x)
    turn = np.nonzero(np.sign(dx[:-1]) != np.sign(dx[1:]))[0]
    if turn.size == 0:
        return {"upper": (x, cp), "lower": None}
    k = int(turn[0]) + 1
    return {"upper": (x[:k + 1], cp[:k + 1]), "lower": (x[k + 1:], cp[k + 1:])}


def _experimental_shock_x(exp):
    """实验激波位置：上表面 Cp 沿 x 的最大**正向**跳。跳幅 < 0.15 视为无激波。"""
    x, cp = exp["upper"]
    o = np.argsort(x)
    x, cp = x[o], cp[o]
    d = np.diff(cp)
    j = int(np.argmax(d))
    return 0.5 * (x[j] + x[j + 1]) if d[j] > 0.15 else None


def _one(level, m_inf, alpha, fn):
    mc, wc = cut_wake(read_mesh(os.path.join(
        str(REPO_ROOT), "cases", "meshes", "naca0012_2.5d", f"{level}.msh")))
    dz = float(np.ptp(mc.nodes[:, 2]))
    r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=alpha, **RECIPE)
    phi = np.asarray(r["phi"])
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
                                alpha_deg=alpha, u_inf=1.0, s_ref=dz, m_inf=m_inf)
    #: ★ 截面 Cp 的**提取与对比**要平滑（使用者裁决 2026-08-25）；载荷不平滑
    cur = section_cp_curve(mc, phi, z=float(np.mean(mc.nodes[:, 2])),
                           smooth_passes=1, m_inf=m_inf)
    exp = _load_experiment(fn)
    rms = {}
    for side in ("upper", "lower"):
        if exp[side] is None:
            rms[side] = float("nan"); continue
        xe, ce = exp[side]
        ci = np.interp(xe, cur[f"x_{side}"], cur[f"cp_{side}"])
        rms[side] = float(np.sqrt(np.mean((ci - ce) ** 2)))
    x_exp = _experimental_shock_x(exp)
    pre = post = float("nan")
    if x_exp is not None:
        xe, ce = exp["upper"]
        ci = np.interp(xe, cur["x_upper"], cur["cp_upper"])
        a, b = xe < x_exp - 0.03, xe > x_exp + 0.03
        if a.sum() > 2: pre = float(np.sqrt(np.mean((ci[a] - ce[a]) ** 2)))
        if b.sum() > 2: post = float(np.sqrt(np.mean((ci[b] - ce[b]) ** 2)))
    h = np.asarray(r["residual_history"], float)
    #: ★★★ `classify_failure` 是给**失败**腿用的：一条收敛腿的残差尾巴天然是平的，
    #: 分类器会把它读成 `limit_cycle`。第一版对每条腿都调它，于是
    #: "mode == limit_cycle" 这个断言**几乎对所有腿都成立、根本不判别** ——
    #: 同一族判据缺陷。⇒ 只对**未收敛**的腿分类。
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
        rms_upper=rms["upper"], rms_lower=rms["lower"],
        pre=pre, post=post, x_shock=shock_report(cur, m_inf)["upper"].get("x_shock"),
        x_shock_exp=x_exp, curve=cur, exp=exp)


@pytest.fixture(scope="module")
def runs():
    return {(name, lv): _one(lv, m, a, fn)
            for name, m, a, fn in CASES for lv in LEVELS}


class TestWhatIsGateableAgainstExperiment:
    r"""★ 裁决三：对实验只记录偏置。**能设门的，只有偏置的方向与结构**，加上
    无粘模型确实负责的那一格（激波前、近零升）。"""

    def test_pre_shock_cp_at_near_zero_lift(self, runs):
        r"""★★ **激波前**、**近零升**：无粘全速势该负责的那一格。
        实测 M0.803 α=−0.10°：coarse 0.0926 → medium **0.0369**（收缩 2.51×）。

        ★ **反向落点**：若这一格也开始随加密恶化，那就不再是"实验含粘性"能解释的了 ——
        问题回到求解器本身。
        """
        c, m = runs[("M0.803_a-0.10", "coarse")], runs[("M0.803_a-0.10", "medium")]
        assert m["converged"] and m["n_limited"] == 0 and m["n_floored"] == 0, (
            f"medium 腿不干净：conv={m['converged']} lim={m['n_limited']} flr={m['n_floored']}")
        assert m["pre"] <= PRE_SHOCK_MAX, (
            f"激波前 RMS = {m['pre']:.4f} > {PRE_SHOCK_MAX}（实测 0.0369）")
        assert c["pre"] / m["pre"] >= PRE_SHOCK_CONTRACT, (
            f"激波前 RMS 只收缩 {c['pre']/m['pre']:.2f}× < {PRE_SHOCK_CONTRACT}"
            f"（实测 2.51×：{c['pre']:.4f} → {m['pre']:.4f}）")

    def test_symmetric_airfoil_at_zero_alpha_has_almost_no_lift(self, runs):
        r"""★ 对称翼型在 α = −0.10°：|cl_p| 必须很小。实测 **0.0423**（medium）。
        这一条几乎不含粘性成分，所以对实验也是可判的。"""
        m = runs[("M0.803_a-0.10", "medium")]
        assert abs(m["cl"]) <= CL_NEAR_ZERO, (
            f"|cl_p| = {abs(m['cl']):.4f} > {CL_NEAR_ZERO}（实测 0.0423）")
        assert m["cl"] < 0, f"α<0 而 cl_p = {m['cl']:+.4f} 不为负 —— 符号错了"

    def test_the_inviscid_shock_sits_downstream_of_the_experimental_one(self, runs):
        r"""★★★ **偏置的方向是可判的，而且它有物理内容**：真实边界层的位移厚度把激波
        **往前推**，所以无粘解的激波必然在实验激波**下游**。

        ★★★ **但这条只挂在偏置远超分辨散布的那一格上，而这是实测逼出来的**：
        M0.778 coarse **+0.1219 c**（设门）；而 **M0.803 的偏置是 coarse −0.0151 c、
        medium +0.0129 c —— 跨越了零**。那一格的激波弱、位置分辨不到 ±0.015 以内，
        **符号在那里没有意义** ⇒ **RECORDED，不设门**。
        ★ 我第一版把这条也挂在 M0.803 medium 上（+0.0129），那是**拿一个小于自身散布的
        数当判据** —— 本项目记过的判据缺陷同族。
        ★ M0.778 的 **medium 腿不参与判据** —— 它是 `limit_cycle`（见下一类）。
        ★★ **反向落点**：若 M0.778 那个 +0.12 c 变成负的，那不是精度问题，是机理反了。
        """
        rec = runs[("M0.803_a-0.10", "coarse")], runs[("M0.803_a-0.10", "medium")]
        offs = [d["x_shock"] - d["x_shock_exp"] for d in rec]
        assert min(offs) < 0 < max(offs), (
            f"M0.803 的激波偏置不再跨越零（{offs}）—— 若它现在稳定为一个符号，"
            "那一格就可以升格为判据了，请重新测量并更新本条")
        for key in (("M0.778_a2.03", "coarse"),):
            d = runs[key]
            assert d["converged"], f"{key} 未收敛，不能引用"
            assert d["x_shock"] is not None and d["x_shock_exp"] is not None, key
            off = d["x_shock"] - d["x_shock_exp"]
            assert off > SHOCK_AFT_MIN, (
                f"{key[0]} {key[1]}: 计算激波 {d['x_shock']:.4f} 不在实验激波 "
                f"{d['x_shock_exp']:.4f} 下游（偏置 {off:+.4f} c）—— 边界层位移的机理反了")

    def test_the_three_cases_stratify_by_how_viscous_they_are(self, runs):
        r"""★★★ **本门最结构性的一条**：无粘模型的误差必须随"该工况多依赖粘性物理"单调增长。
        实测 medium 上表面 RMS：**0.0633（近零升） < 0.3986（α=2°） < 1.5477（近失速）**。

        ★ 这不是精度判据，是**适用范围**判据：它一旦破了，说明三个工况的相对难度变了，
        而本门的整套"记录偏置"叙述就要重读。
        """
        v = [runs[(n, "medium")]["rms_upper"] for n, _, _, _ in CASES]
        for a, b in zip(v, v[1:]):
            assert b / a >= STRATIFY, (
                f"分层破了：{['%.4f' % q for q in v]} —— 相邻比 {b/a:.2f} < {STRATIFY}")


class TestWhatIsOnlyRecorded:
    r"""★ 裁决一：(d) 激波后 Cp 只记录。本类把"只记录"的那些**锁成事实**，
    这样它们一旦改变会响，而不是被悄悄当成改进或退化。"""

    def test_post_shock_does_not_improve_with_refinement(self, runs):
        r"""★★ 激波后 RMS **不随加密改善**（实测 0.0609 → **0.0807**，变差）——
        那正是"实验含激波-边界层干涉、无粘模型不含"的直接表现，也正是裁决一
        把它排除在判据外的理由。

        ★★ 本条断言的是**这个事实仍然成立**。若它变成"改善"，那是好消息，
        但**必须走纪律 11 勘误**：本门"只记录激波后"的整个理由建立在它之上。
        """
        c, m = runs[("M0.803_a-0.10", "coarse")], runs[("M0.803_a-0.10", "medium")]
        assert m["post"] > c["post"], (
            f"激波后 RMS 现在随加密改善了（{c['post']:.4f} → {m['post']:.4f}）—— "
            "好消息，但本门把 (d) 排除在判据外的理由变了，请走纪律 11 勘误")

    def test_near_stall_is_out_of_the_inviscid_envelope(self, runs):
        r"""★ M0.352 α = 12.86°：上表面 RMS **1.0451 → 1.5477**，且**加密使其恶化**。
        近失速的分离流不在无粘全速势的包络内 —— **记录这个事实，而不是给它一个容差**。
        ★★ 该工况**没有可比的激波结构**，而这一点要按实测说，不能靠推断：检波器确实
        在实验数据上报出一个"最大正跳"，位置 **x = 0.0237**（本条断言它 < 0.05）——
        那是 α = 12.86° 下从 Cp = −6.81 起的**前缘吸力恢复**，不是激波。
        ★ 我第一版写的是 `x_shock_exp is None`，那是从 `pre = nan` **推断**来的；
        而 nan 的真实原因是"x = 0.024 之前不足 2 个点"。**又一次没有引用管线实际跑出的数。**
        ⇒ 后果也一并锁住：**激波前 RMS 必须是 nan**，即这一格没有"激波前"可判。
        """
        c, m = runs[("M0.352_a12.86", "coarse")], runs[("M0.352_a12.86", "medium")]
        assert m["rms_upper"] > c["rms_upper"] > 0.8, (
            f"近失速工况的 RMS 不再又大又随加密恶化（{c['rms_upper']:.4f} → "
            f"{m['rms_upper']:.4f}）—— 无粘包络的边界变了，请重新测量并更新 docstring")
        assert m["x_shock_exp"] < 0.05, (
            f"检波器在该工况报出的最大正跳位置 {m['x_shock_exp']:.4f} 不再落在前缘 "
            "(<0.05) —— 它可能真的检出了一个激波，本条的前提变了")
        assert np.isnan(m["pre"]), (
            f"该工况现在有可判的激波前区段（pre = {m['pre']:.4f}）—— 前提变了")

    def test_the_failing_leg_is_classified_not_just_reported(self, runs):
        r"""★★★ 纪律：**分类它，永远不要只报 conv=False**。
        M0.778 medium 实测 `limit_cycle`（`descent10 = 0.599`），残差尾巴是**周期 3**：
        2.29e-07 → 3.82e-07 → 2.61e-07，三位有效数字重复。

        ★★ 而分类的**后果**是本条真正要锁的：**那一腿的数字不被任何判据引用**。
        上面每一条用到 M0.778 的地方都显式取 coarse。
        """
        d = runs[("M0.778_a2.03", "medium")]
        assert not d["converged"], "M0.778 medium 现在收敛了 —— 前提变了，可以重新考虑用它"
        assert d["mode"] == "limit_cycle", (
            f"分类从 limit_cycle 变成 {d['mode']}（{d['mode_detail']}）—— "
            "失败模式变了，修法也就不同（见 CLAUDE.md 的模式表）")
        assert d["n_limited"] == 0 and d["n_floored"] == 0, (
            "出现钳制了 —— 那是另一种模式（physics/geometry），不是 limit_cycle")
        #: ★★ 判别性：收敛的腿**不得**带失败模式标签。第一版没有这一条，
        #: 于是分类器把每条收敛腿也标成 limit_cycle 而断言照样绿。
        for key in (("M0.803_a-0.10", "medium"), ("M0.778_a2.03", "coarse")):
            assert runs[key]["mode"] == "converged", (
                f"{key} 收敛却被标成 {runs[key]['mode']} —— 分类器被用在了它不适用的腿上")


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_experiment_bias_figure(runs, gate_evidence_dir):
    import csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for k, (name, m_inf, alpha, _) in enumerate(CASES):
        d = runs[(name, "medium")]
        for side, mk in (("upper", "-"), ("lower", "--")):
            if d["exp"][side] is None:
                continue
            xe, ce = d["exp"][side]
            ax[k].plot(xe, ce, "ko", ms=3,
                       label="experiment" if side == "upper" else None)
            ax[k].plot(d["curve"][f"x_{side}"], d["curve"][f"cp_{side}"], mk,
                       lw=1.2, color="C0",
                       label="inviscid FP (medium)" if side == "upper" else None)
        if d["x_shock_exp"] is not None:
            ax[k].axvline(d["x_shock_exp"], color="k", ls=":", lw=1)
        if d["x_shock"] is not None:
            ax[k].axvline(d["x_shock"], color="C0", ls=":", lw=1)
        ax[k].set_title(f"M {m_inf}  alpha {alpha}\nupper RMS {d['rms_upper']:.4f}"
                        + ("  (limit_cycle)" if not d["converged"] else ""), fontsize=9)
        ax[k].set_xlabel("x/c"), ax[k].invert_yaxis(), ax[k].grid(alpha=.3)
        ax[k].legend(fontsize=7)
    ax[0].set_ylabel("Cp")
    fig.suptitle("D11  inviscid full potential vs experiment -- the BIAS is the product",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "d11_experiment_bias.png"), dpi=130)
    plt.close(fig)

    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "level", "converged", "mode", "residual", "n_newton",
                    "n_limited", "n_floored", "mach_max", "cl_p", "rms_upper",
                    "rms_lower", "pre_shock_rms", "post_shock_rms",
                    "x_shock", "x_shock_exp", "shock_offset"])
        for name, _, _, _ in CASES:
            for lv in LEVELS:
                d = runs[(name, lv)]
                off = ("" if d["x_shock"] is None or d["x_shock_exp"] is None
                       else f"{d['x_shock'] - d['x_shock_exp']:+.4f}")
                w.writerow([name, lv, d["converged"], d["mode"], f"{d['residual']:.3e}",
                            d["n_newton"], d["n_limited"], d["n_floored"],
                            f"{d['mach_max']:.4f}", f"{d['cl']:.6f}",
                            f"{d['rms_upper']:.4f}", f"{d['rms_lower']:.4f}",
                            f"{d['pre']:.4f}", f"{d['post']:.4f}",
                            "" if d["x_shock"] is None else f"{d['x_shock']:.4f}",
                            "" if d["x_shock_exp"] is None else f"{d['x_shock_exp']:.4f}",
                            off])
