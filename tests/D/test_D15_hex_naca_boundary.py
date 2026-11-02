r"""D15 — **hex 族 NACA 能力边界门：「缺陷仍在」式**（phase 7 第 3 轮，S2）。

★★★ **它断言的是「缺陷还在」，所以红 = 好消息** —— 形状照抄 D09/D10/D14。

## 缺陷内容（p6 第 47 轮归因闭合；证据 `phases/p6/bench/studies/hex_cl_attribution/`）

hex 族（`cases/meshes/naca0012_hex_2.5d/`）在 M1 靶子工况 M0.80/α1.25 上收敛到
tet 族**从未出现**的上翼面中部气泡/双激波结构：

| 量（medium） | hex | tet | 出处 |
|---|---|---|---|
| cl | **0.258223** | **0.341340** | 加细不收缩：fine 0.256059（等 h 仍差 −25 %） |
| 气泡计（上翼面 [0.2,0.55] max Cp） | **+0.199 @x=0.457** | −0.891 | coarse −0.297 → medium +0.199 → fine +0.383：**涌现、增强、前移** |
| 逐段 Δcl（hex−tet，总 −0.0827） | front −0.0341 · shock −0.0509 | — | LE −0.0003 · TE +0.0026：**亏损 100 % 在前中段+激波区** |
| x_shock | 0.6094 | 0.5991 | 差 0.010 c：**缺陷在升力/结构，不在激波定位** |

**残余未分开**（生成器缺陷 vs 熵冻结交互：hex 各腿 sigma 冻结搅动 |dσ| 0.19–0.21
vs tet 0.067）—— 本轮的拆分探针在 `bench/studies/hex_residual_split/`，
分开与否如实登记，**门不依赖探针的结论**。

## ★ 红了怎么办（写在这里，免得有人改回去换绿）

气泡消失 / 亏损闭合 = **缺陷被修好，是好消息**：走纪律 10/11 勘误 —— 重读本文件
全部判据、`cases/gates/INDEX.md` 条目、能力边界 §9 的移交句，并把本门**重新定型**
为一条正常回归门。**不许把判据改宽去换绿。**

## G-CADENCE（周期登记）

- **默认套件**：A1–A3、A5、A6、G-TEETH(T2)、证据锁（hex coarse/medium + tet medium 锚）。
- **门控全集**（`PYFP3D_TRANSONIC_GATES=1`）：加 A4 的 fine 腿（它撞 400 步帽，
  是最贵的一条）与导图腿。
- fine 腿**进 CSV 但不进值锁**：它是 limit-cycle 停滞腿，cl 对停滞水平的散布
  有 D05 线程实验标定（0.15–0.36 %）—— 对它只断言不等式（指纹增强），
  这正是 D14「只锁跨线程稳定的腿」的做法。

## 线程稳定性（A/B 实测，决定锁名单）

2026-10-12 第 3 轮实测（同一代码、同一网格，@8t vs @16t committed 证据）：
三条默认腿的四列被锁值相对差全部 **1e-11 … 2e-10**（最大 = hex coarse 气泡计
1.5e-10）⇒ **cl / x_shock / bubble_max_cp / bubble_x_at_max 全部可锁**，
锁名单 = 三条默认腿全列。hex fine 不进值锁（停滞腿，散布有 D05 标定 0.15–0.36 %），
对它只断言 A4b 的不等式。
"""
import os

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.shock import shock_report
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting
from tests._gate_evidence import assert_matches_committed, blas_thread_count, fmt
from tests.conftest import REPO_ROOT, gate_figures_enabled

#: ★ 判据配方**从 D05 导入**，不重打一遍（D14 同款）：抄一份常数就会让两边悄悄分叉。
from tests.D.test_D05_euler_naca0012 import RECIPE

M_INF, ALPHA = 0.80, 1.25
MESH_DIRS = {"tet": "naca0012_2.5d", "hex": "naca0012_hex_2.5d"}
#: 默认套件的腿；hex fine 撞 400 步帽、最贵 ⇒ 门控全集才跑（G-CADENCE 见 docstring）
LEGS = (("hex", "coarse"), ("hex", "medium"), ("tet", "medium"))
GATED_LEGS = (("hex", "fine"),)

#: ★ 气泡窗口与分段边界**事前写死**（与研究逐字一致，看数之前的约定不许事后移动）
BUBBLE_X0, BUBBLE_X1 = 0.20, 0.55
SEGMENTS = (("LE_peak", 0.00, 0.05), ("front", 0.05, 0.45),
            ("shock", 0.45, 0.75), ("TE_unload", 0.75, 1.00))

#: —— 判据（标定；实测基线见模块 docstring 的表）——
BUBBLE_PRESENT_MIN = 0.0    # ★ 现象定义本身（符号翻转 = 气泡存在），不是容差
DEFICIT_MIN = 0.05          # cl 亏损下界：实测 0.083，余量 ~40 %
XSHOCK_AGREE_MAX = 0.02     # 两族激波位置一致：实测 0.010 c，余量 2×
FRONT_SHOCK_SHARE_MIN = 0.90    # front+shock 占 |总亏损| 的份额：实测 1.03
EDGE_SHARE_MAX = 0.05       # LE / TE 各占份额上界：实测 0.004 / 0.031

TRANSONIC = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"


def _one(family, level):
    """一腿：求解 + 力系数 + 截面 Cp + 激波 + 气泡计 + 分段 cl。"""
    mc, wc = cut_wake(read_mesh(os.path.join(
        str(REPO_ROOT), "cases", "meshes", MESH_DIRS[family], f"{level}.msh")))
    dz = float(np.ptp(mc.nodes[:, 2]))
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **RECIPE)
    phi = np.asarray(r["phi"])
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], phi,
                                alpha_deg=ALPHA, u_inf=1.0, s_ref=dz,
                                m_inf=M_INF)
    cur = section_cp_curve(mc, phi, z=float(np.mean(mc.nodes[:, 2])),
                           smooth_passes=1, m_inf=M_INF)
    up = shock_report(cur, M_INF)["upper"]
    frz = r.get("sigma_freeze_report") or {}
    hist = np.asarray(r["residual_history"], float)
    bubble_cp, bubble_x = _bubble(cur)
    return dict(family=family, level=level,
                n_nodes=int(mc.nodes.shape[0]), n_tets=int(mc.elements.shape[0]),
                cl=float(f["cl"]), cd=float(f["cd_pressure"]),
                x_shock=(float(up["x_shock"]) if up.get("has_shock")
                         else float("nan")),
                bubble_max_cp=bubble_cp, bubble_x_at_max=bubble_x,
                converged=bool(r.get("converged")),
                residual=float(hist[-1]), n_newton=int(len(hist)),
                n_limited=int(r.get("n_limited", 0)),
                n_floored=int(r.get("n_floored", 0)),
                sigma_frozen=bool(frz.get("frozen_in_transient", False)),
                sigma_churn=bool(frz.get("selection_churn", False)),
                curve=cur)


def _sorted_unique(x, y):
    o = np.argsort(x)
    xs, ys = np.asarray(x)[o], np.asarray(y)[o]
    keep = np.concatenate(([True], np.diff(xs) > 0))
    return xs[keep], ys[keep]


def _bubble(cur):
    """气泡计：上翼面 x∈[0.2,0.55] 的 max Cp 及其位置（正 = 中部减速气泡存在）。"""
    xs, ys = _sorted_unique(cur["x_upper"], cur["cp_upper"])
    m = (xs >= BUBBLE_X0) & (xs <= BUBBLE_X1)
    i = int(np.argmax(ys[m]))
    return float(ys[m][i]), float(xs[m][i])


def _segment_cl(cur, x0, x1, n=4001):
    """一段对 cl 的贡献：∫(cp_lower − cp_upper) dx（2-D 截面，弦长 1）。"""
    xg = np.linspace(x0, x1, n)
    tot = 0.0
    for side, sgn in (("lower", +1.0), ("upper", -1.0)):
        xs, ys = _sorted_unique(cur[f"x_{side}"], cur[f"cp_{side}"])
        tot += sgn * float(np.trapezoid(np.interp(xg, xs, ys), xg))
    return tot


def _defect_present(bubble_max_cp, deficit):
    """A1 ∧ A3 的判别函数。★ G-TEETH T2 拿 tet 对照喂它，**必须返回 False**
    —— 否则这两条判据是恒真的装饰（`TestTheCriterionIsNotVacuous`）。"""
    return bubble_max_cp > BUBBLE_PRESENT_MIN and deficit > DEFICIT_MIN


@pytest.fixture(scope="module")
def runs():
    return {(fam, lv): _one(fam, lv) for fam, lv in LEGS}


@pytest.fixture(scope="module")
def gated_runs():
    """门控腿（hex fine）。只在 PYFP3D_TRANSONIC_GATES=1 下被消费 ⇒
    默认套件里这个 fixture 不会被建立。"""
    return {(fam, lv): _one(fam, lv) for fam, lv in GATED_LEGS}


class TestTheDefectIsStillThere:
    r"""★★★ 断言「缺陷仍在」。**红 = 修好了 = 好消息**，处置见模块 docstring ——
    红了走纪律 10/11 勘误，**不许把判据改宽去换绿**。"""

    def test_the_bubble_is_still_on_hex_medium(self, runs):
        """A1：hex medium 上翼面中部的减速气泡仍在（max Cp > 0）。
        研究实测 +0.199 @x=0.457；tet 对照 −0.891。"""
        d = runs[("hex", "medium")]
        assert d["bubble_max_cp"] > BUBBLE_PRESENT_MIN, (
            f"hex medium 的气泡消失了（max Cp {d['bubble_max_cp']:+.4f} "
            f"@x={d['bubble_x_at_max']:.3f}）—— **this is an improvement, "
            f"not a regression**: 缺陷被修好或绕开了。走纪律 10/11 勘误并重定型本门，"
            "不许调 BUBBLE_PRESENT_MIN")

    def test_the_tet_control_has_no_bubble(self, runs):
        """A2：tet medium 同窗**没有**气泡（研究实测 −0.891）。
        ★ 若 tet 也长出气泡 ⇒ 现象不再是 hex 专属，归因形状变了 ⇒ 重读。"""
        d = runs[("tet", "medium")]
        assert d["bubble_max_cp"] < -BUBBLE_PRESENT_MIN, (
            f"tet medium 也长出了气泡（max Cp {d['bubble_max_cp']:+.4f} "
            f"@x={d['bubble_x_at_max']:.3f}）—— 现象不再是 hex 族专属，"
            "「hex 收敛到 tet 没有的结构」这个前提变了，重新归因")

    def test_the_lift_deficit_is_still_there(self, runs):
        """A3：cl 亏损仍在（cl_tet − cl_hex > 0.05；研究实测 0.083）。"""
        deficit = runs[("tet", "medium")]["cl"] - runs[("hex", "medium")]["cl"]
        assert deficit > DEFICIT_MIN, (
            f"hex medium 的 cl 亏损闭合到 {deficit:+.4f}（带 {DEFICIT_MIN}）—— "
            f"hex {runs[('hex','medium')]['cl']:.6f} vs tet "
            f"{runs[('tet','medium')]['cl']:.6f}。**this is an improvement**: "
            "走勘误并把本门重定型为正常回归门；不许调 DEFICIT_MIN")

    def test_the_fingerprint_strengthens_coarse_to_medium(self, runs):
        """A4a：指纹 = 气泡随加细**涌现**。hex coarse → medium 气泡计上升
        （研究实测 −0.297 → +0.199）。变平/反转 ⇒ 结构性质变了 ⇒ 红，重读。"""
        c, m = runs[("hex", "coarse")], runs[("hex", "medium")]
        assert m["bubble_max_cp"] > c["bubble_max_cp"], (
            f"气泡计不再随加细增强：coarse {c['bubble_max_cp']:+.4f} → "
            f"medium {m['bubble_max_cp']:+.4f} —— 指纹变了，重读归因")

    def test_shock_position_agrees_across_families(self, runs):
        """A5：两族激波位置一致（≤ 0.02 c；研究实测 0.010 c）——
        把缺陷约束在 cl/结构上；超带 ⇒ 缺陷形状变了 ⇒ 红，重读。"""
        d = abs(runs[("hex", "medium")]["x_shock"]
                - runs[("tet", "medium")]["x_shock"])
        assert d <= XSHOCK_AGREE_MAX, (
            f"两族激波位置分家到 {d:.4f} c（带 {XSHOCK_AGREE_MAX}）—— "
            "缺陷不再只是升力/结构，激波定位也卷进来了，重读归因")

    def test_the_deficit_lives_in_front_and_shock_segments(self, runs):
        """A6：亏损定位 —— front+shock 两段占 |总亏损| ≥ 90 %，LE/TE 各占 ≤ 5 %
        （研究实测 1.03 / 0.004 / 0.031）。亏损搬到 LE/TE ⇒ 归因形状变了 ⇒ 红。"""
        hm, tm = runs[("hex", "medium")]["curve"], runs[("tet", "medium")]["curve"]
        d = {name: _segment_cl(hm, x0, x1) - _segment_cl(tm, x0, x1)
             for name, x0, x1 in SEGMENTS}
        total = abs(sum(d.values()))
        share_fs = (abs(d["front"]) + abs(d["shock"])) / total
        assert share_fs >= FRONT_SHOCK_SHARE_MIN, (
            f"front+shock 份额降到 {share_fs:.3f}（带 {FRONT_SHOCK_SHARE_MIN}）—— "
            f"逐段差 {d}；亏损的定位变了，重读归因")
        for edge in ("LE_peak", "TE_unload"):
            share = abs(d[edge]) / total
            assert share <= EDGE_SHARE_MAX, (
                f"{edge} 份额涨到 {share:.3f}（带 {EDGE_SHARE_MAX}）—— "
                f"逐段差 {d}；亏损不再局限在前中段+激波区，重读归因")


@pytest.mark.skipif(not TRANSONIC,
                    reason="hex fine 撞 400 步帽、最贵 —— set PYFP3D_TRANSONIC_GATES=1")
class TestTheFingerprintAtFine:
    r"""A4b：fine 腿按 **physical 口径**读（\|R\| ≤ 1e-6；converged 旗标是
    线程依赖噪声，不对它设断言 —— D05 的教训）。fine 进 CSV 但**不进值锁**
    （停滞腿的 cl 散布有 D05 标定 0.15–0.36 %）。"""

    def test_the_fingerprint_strengthens_medium_to_fine(self, runs, gated_runs):
        """研究实测：medium +0.199 → fine +0.383（增强并前移，x 0.457 → 0.376）。"""
        f = gated_runs[("hex", "fine")]
        assert f["residual"] <= 1e-6, (
            f"hex fine 连 physical 口径都不到了（|R| {f['residual']:.2e} > 1e-6）—— "
            "它的读数不再可用，重测")
        m = runs[("hex", "medium")]
        assert f["bubble_max_cp"] > m["bubble_max_cp"], (
            f"气泡计 medium→fine 不再增强：{m['bubble_max_cp']:+.4f} → "
            f"{f['bubble_max_cp']:+.4f} —— 指纹变了，重读归因")


class TestTheCriterionIsNotVacuous:
    r"""★★ G-TEETH T2（常驻）：把 **tet 对照**的新鲜数据喂给 A1∧A3 的判别函数，
    它必须返回「缺陷不存在」—— 否则那两条判据是恒真的装饰，测不出修好。"""

    def test_the_predicate_absolves_the_tet_control(self, runs):
        tm = runs[("tet", "medium")]
        deficit_against_itself = tm["cl"] - tm["cl"]
        assert not _defect_present(tm["bubble_max_cp"], deficit_against_itself), (
            f"判别函数对 tet 对照也报「缺陷在」（bubble {tm['bubble_max_cp']:+.4f}）"
            "—— 判据恒真，G-TEETH 失败")
        #: 而它对自己该报「在」的数据（hex medium）必须报「在」—— 双向
        hm = runs[("hex", "medium")]
        deficit = tm["cl"] - hm["cl"]
        assert _defect_present(hm["bubble_max_cp"], deficit), (
            "判别函数对 hex medium 报「缺陷不在」—— 与本文件头部的实测表矛盾，"
            "要么缺陷被修好（好消息，走勘误），要么判别函数写错了")


class TestCommittedEvidenceIsLoadBearing:
    r"""新鲜计算 vs 已提交 `summary.csv`。设计见 `tests/_gate_evidence.py`。

    ★ **只锁三条默认腿**；hex fine 进 CSV 但不进值锁（停滞腿，散布有 D05 标定）。
    ★ `residual` 不进值锁（D14 已登记：机器零附近比的是舍入对舍入）。
    ★ 锁名单由本轮 8t vs 16t 的 A/B 实测决定 —— 结果写在模块 docstring 的表里。
    """

    MEASURED = ("cl", "x_shock", "bubble_max_cp", "bubble_x_at_max", "blas_threads")

    def test_matches_committed_summary(self, runs, gate_evidence_dir):
        fresh = {f"{fam}|{lv}": {**{k: fmt(runs[(fam, lv)][k])
                                    for k in self.MEASURED if k != "blas_threads"},
                                 "blas_threads": fmt(blas_thread_count())}
                 for fam, lv in LEGS}
        n = assert_matches_committed(
            gate_evidence_dir, fresh, self.MEASURED,
            key_of=lambda r: r["leg"],
            refresh_hint="PYFP3D_TRANSONIC_GATES=1 PYFP3D_GATE_FIGURES=1 pytest "
                         "tests/D/test_D15_hex_naca_boundary.py")
        assert n >= 15, f"只比了 {n} 个数（3 条默认腿 x 5 列 = 15）"


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图/CSV 证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_evidence(runs, gated_runs, gate_evidence_dir):
    """写 `summary.csv` + Cp 对照图。

    ★ 刷新是**两遍**（`tests/_gate_evidence.py` 第 3 条）：先带
    `PYFP3D_TRANSONIC_GATES=1 PYFP3D_GATE_FIGURES=1` 跑一遍刷新（fine 腿在这一遍
    进 CSV），再不带标志跑一遍验证。
    """
    import csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    all_legs = dict(runs)
    if TRANSONIC:
        all_legs.update(gated_runs)

    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        #: ★ 出处列：n_threads 是 numba 旋钮，blas_threads 读库（H30 的教训）
        w.writerow(["leg", "n_nodes", "blas_threads", "n_threads",
                    "cl", "cd", "x_shock", "bubble_max_cp", "bubble_x_at_max",
                    "converged", "residual", "n_newton", "n_limited",
                    "n_floored", "sigma_frozen", "sigma_churn"])
        for (fam, lv), d in all_legs.items():
            w.writerow([f"{fam}|{lv}", d["n_nodes"], blas_thread_count(),
                        os.environ.get("NUMBA_NUM_THREADS", "unset"),
                        fmt(d["cl"]), fmt(d["cd"]), fmt(d["x_shock"]),
                        fmt(d["bubble_max_cp"]), fmt(d["bubble_x_at_max"]),
                        int(d["converged"]), fmt(d["residual"]), d["n_newton"],
                        d["n_limited"], d["n_floored"],
                        int(d["sigma_frozen"]), int(d["sigma_churn"])])

    #: ★ 主图是 Cp 分布（使用者裁决 2026-09-05）：气泡在 Cp 上肉眼可见，
    #:   cl 柱状图藏住它。
    fig, ax = plt.subplots(figsize=(10, 6.0))
    styles = {("hex", "coarse"): dict(lw=0.9, ls=":", color="#9aa0a6"),
              ("hex", "medium"): dict(lw=1.3, ls="--", color="#e08a00"),
              ("hex", "fine"): dict(lw=1.5, color="#2a9d3f"),
              ("tet", "medium"): dict(lw=1.3, ls="--", color="#1a56db")}
    for key, sty in styles.items():
        d = all_legs.get(key)
        if d is None:
            continue
        for side in ("upper", "lower"):
            xs, ys = _sorted_unique(d["curve"][f"x_{side}"], d["curve"][f"cp_{side}"])
            ax.plot(xs, ys, label=" ".join(key) if side == "upper" else None, **sty)
    ax.axvspan(BUBBLE_X0, BUBBLE_X1, alpha=0.08, color="#c0392b",
               label="bubble window")
    ax.invert_yaxis(), ax.grid(alpha=.3)
    ax.set_xlabel("x/c"), ax.set_ylabel("$C_p$")
    hm, tm = runs[("hex", "medium")], runs[("tet", "medium")]
    ax.set_title(f"D15 hex NACA bubble, M{M_INF}/α{ALPHA} — "
                 f"cl: hex med {hm['cl']:.4f} / tet med {tm['cl']:.4f}; "
                 f"bubble hex med {hm['bubble_max_cp']:+.3f} "
                 f"(tet {tm['bubble_max_cp']:+.3f})", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "d15_hex_naca_boundary.png"),
                dpi=130)
    plt.close(fig)
