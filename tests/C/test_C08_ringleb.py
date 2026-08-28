r"""C08 — Ringleb：**唯一 2-D 跨声速精确解**，且人工耗散第一次能对真值定价（C 类）。

精确解在 `tests/_ringleb_case.py`（那份文件的 docstring 记着验它的**七个 oracle**，
以及被 oracle 抓出来的**两个真实符号错误**）。求解走 `tests/_nozzle_case.DuctSystem`
（进出口 Dirichlet、壁面自然）—— 与 C05 同一个驱动，因此**同样护住 `pyfp3d/kernels/`
的人工密度/迎风/熵修正，同样护不到 `solve/newton.py` 的全局化那一层**。

★★★ **C08 与 C05 的分工，是本门存在的全部理由。**
喷管里**有激波** ⇒ 数值解**本该**与等熵精确解在激波处不同，所以那里只能判**激波位置**。
Ringleb 的超声速口袋是**无激波、光滑**的 ⇒ **人工密度在那里加的每一分耗散都是纯误差**。

★★★ **实测（h = 0.08，1080 节点，469/2730 个超声速单元，精确 M_max = 1.2464）：**

| `upwind_c` | `m_crit` | 收敛 | 超声速区 RMS | 亚声速区 RMS | 数值 M_max |
|---|---|---|---|---|---|
| **1.50（生产）** | 0.95 | ✓ 11 步 | **1.00e-01** | 4.68e-02 | 1.2182 |
| 0.50 | 0.95 | ✓ 6 步 | 4.77e-02 | 1.99e-02 | 1.2676 |
| 0.25 | 0.95 | ✓ 8 步 | 3.76e-02 | 1.47e-02 | 1.2974 |
| **0.00** | 0.95 | ★ **不收敛** | 4.05e-02 | 1.02e-02 | 1.3102 |
| 1.50 | **1.20** | ✓ 37 步 | **2.93e-02** | 9.86e-03 | 1.2747 |

⇒ 生产配置在无激波的超声速区造成 **10.0 % 的速度误差**；把 C 降到 0.25 减到 3.8 %
（**2.7 倍**）；而 **C = 0 直接不收敛** ⇒ **那份耗散买的是稳健性，代价现在有价格。**

★★★ **判据的形状必须说清楚，否则会被误读成缺陷。** 裁决③要「收敛阶 ∧ 绝对误差」：
· **阶**判在 **`R(exact)`** —— 把精确解代进离散残差，实测 **2.32 / 2.28 阶**收敛
  ⇒ **离散化本身是相容的，而且是二阶的**。这一半**不需要求解**，因此也便宜。
· **绝对值**判在**解**上。而**解的误差不随 h 收敛**（实测 RMS 6.4e-2 → 5.9e-2 → 7.3e-2）
  —— **那不是缺陷也不是 bug，是模型误差**：人工密度不随 h 消失。
  **证明它是模型误差的不是论证，是上表的剂量-响应**：误差随 `upwind_c` 单调下降。
"""
import os

import numpy as np
import pytest

from pyfp3d.physics.isentropic import mach_number_squared
from tests import _ringleb_case as RB
from tests._nozzle_case import DuctSystem
from tests._gate_evidence import assert_matches_committed, fmt
from tests.conftest import REPO_ROOT, gate_figures_enabled

LEVELS = (0.12, 0.08, 0.055)
H_STUDY = 0.08                      # 耗散扫描所在的级
C_SWEEP = (1.5, 0.5, 0.25, 0.0)     # 生产 -> 0
UPWIND_PROD, MCRIT_PROD = 1.5, 0.95

#: 判据（非实测值；实测基线见 docstring 的表）
ORDER_MIN = 1.50        # R(exact) 的阶：实测 2.32 / 2.28，余量 52 %
SUP_RMS_MAX = 0.13      # 生产配置在超声速区的速度 RMS：实测 0.1001
SUP_OVER_SUB = 1.5      # 超声速区误差必须显著高于亚声速区：实测 2.14x
DOSE_MIN = 2.0          # C 1.5 -> 0.25 的误差收缩：实测 2.66x
MACH_REL_MAX = 0.05     # 生产 M_max 对精确的相对偏差：实测 -2.26 %


def _solve(case, C=UPWIND_PROD, m_crit=MCRIT_PROD, n_max=60):
    s = DuctSystem(case["mesh"], m_inf=RB.M_INF, upwind_c=C, m_crit=m_crit)
    #: ★ `assemble_residual` 返回**共享缓冲区的视图** —— 不 copy 会被 Newton 覆写。
    #: 本轮实测踩过：R(exact) 读到 1.2e-14（其实是最终残差），据此算出的"阶 2.1/2.8"是假数。
    r0, _ = s.residual(case["phi_exact"])
    r_exact = float(np.max(np.abs(r0.copy()[s.free])))
    phi, info = s.newton(case["phi_exact"].copy(), n_max=n_max, tol=1e-11)
    _, q2 = s.op.velocities(phi)
    q2 = q2.copy()                                   # ★ 同上，velocities 也是视图
    e = np.sqrt(q2) - case["qhat_cell"]
    sup = case["mach_cell"] > 1.0
    m_num = np.sqrt(np.maximum(mach_number_squared(q2, RB.M_INF, RB.G), 0.0))
    return dict(
        r_exact=r_exact, converged=bool(info["converged"]), reason=info["reason"],
        n_newton=int(info["n_newton"]), res=float(info["residual_history"][-1]),
        rms=float(np.sqrt(np.mean(e**2))), max=float(np.abs(e).max()),
        rms_sup=float(np.sqrt(np.mean(e[sup]**2))),
        rms_sub=float(np.sqrt(np.mean(e[~sup]**2))),
        n_sup=int(sup.sum()), n_cell=int(sup.size),
        m_num_max=float(m_num.max()), m_ex_max=float(case["mach_cell"].max()),
        q2=q2, err=e)


@pytest.fixture(scope="module")
def cases():
    return {h: RB.build_case(h) for h in LEVELS}


@pytest.fixture(scope="module")
def prod(cases):
    return _solve(cases[H_STUDY])


@pytest.fixture(scope="module")
def sweep(cases):
    c = cases[H_STUDY]
    out = {C: _solve(c, C=C) for C in C_SWEEP}
    out["mcrit120"] = _solve(c, C=UPWIND_PROD, m_crit=1.20)
    return out


class TestTheExactSolutionIsAConsistentDiscreteSolution:
    r"""★ 先证明尺子配得上 —— 阶判在这里，而且**一次求解都不用**。"""

    def test_residual_of_the_exact_solution_converges(self, cases):
        r"""把精确解代进离散残差，它必须随 h 以二阶量级趋零。实测 2.32 / 2.28。

        ★★ 这一条是本门**唯一**的收敛阶判据，理由写在 docstring：**解**的误差
        由人工密度主导、不随 h 收敛，所以拿它判阶会把一个模型误差误报成离散化缺陷。
        ★ 反向落点：若离散化本身坏了（装配/度量），这一条会掉向 0 而下面的耗散
        剂量-响应仍可能正常 —— 两者分工不同，都要有。
        """
        rs = []
        for h in LEVELS:
            s = DuctSystem(cases[h]["mesh"], m_inf=RB.M_INF,
                           upwind_c=UPWIND_PROD, m_crit=MCRIT_PROD)
            r0, _ = s.residual(cases[h]["phi_exact"])
            rs.append(float(np.max(np.abs(r0.copy()[s.free]))))
        for i in range(1, len(LEVELS)):
            p = float(np.log(rs[i-1]/rs[i]) / np.log(LEVELS[i-1]/LEVELS[i]))
            assert p >= ORDER_MIN, (
                f"R(exact) 阶 h {LEVELS[i-1]} -> {LEVELS[i]} = {p:.3f} < {ORDER_MIN}"
                f"（实测 2.32 / 2.28；三级读数 {['%.2e' % r for r in rs]}）")

    def test_the_domain_really_has_a_supersonic_pocket(self, cases):
        r"""★ 前提：没有超声速区，本门的全部内容都不存在。
        实测 h=0.08 上 469/2730 个单元超声速，精确 M_max = 1.2464。"""
        c = cases[H_STUDY]
        n_sup = int((c["mach_cell"] > 1.0).sum())
        assert n_sup > 100, f"超声速单元只有 {n_sup} 个（实测 469）"
        assert 1.15 < c["mach_cell"].max() < 1.35, (
            f"精确 M_max = {c['mach_cell'].max():.4f} 离开了 1.15~1.35 的设计带 —— "
            "几何参数被改了？峰值 Mach 是刻意压在真实翼型激波量级的")


class TestArtificialDissipationPricedAgainstTruth:
    r"""★★★ 本门的核心：在**无激波**的超声速区，人工密度加的每一分都是纯误差。"""

    def test_production_absolute_error(self, prod):
        r"""生产配置（C = 1.5, m_crit = 0.95）在超声速区的速度 RMS。实测 0.1001。"""
        assert prod["converged"], f"生产腿不收敛（{prod['reason']}）"
        assert prod["rms_sup"] <= SUP_RMS_MAX, (
            f"超声速区速度 RMS = {prod['rms_sup']:.4e} > {SUP_RMS_MAX}（实测 0.1001）")

    def test_the_error_lives_in_the_supersonic_pocket(self, prod):
        r"""★★ 超声速区误差必须显著高于亚声速区 —— 这是**归因**，不是精度。
        实测 1.00e-01 vs 4.68e-02 = **2.14x**。若两者拉平，误差就不是人工密度了，
        本门余下的剂量-响应也就失去前提。"""
        ratio = prod["rms_sup"] / prod["rms_sub"]
        assert ratio >= SUP_OVER_SUB, (
            f"超声速/亚声速 误差比 = {ratio:.2f} < {SUP_OVER_SUB}（实测 2.14）")

    def test_dose_response_in_upwind_c(self, sweep):
        r"""★★★ **剂量-响应**：误差随 `upwind_c` 单调下降。

        实测超声速区 RMS：C=1.5 **1.00e-01** → 0.5 **4.77e-02** → 0.25 **3.76e-02**
        （1.5 → 0.25 收缩 **2.66x**）。
        ★ 这一条**证明**那 10 % 是人工密度而不是离散化 —— 论证做不到，剂量-响应做得到。
        ★ C = 0 不进单调判据：它**不收敛**（实测 60 步、\|R\| 2.5e-07），
        一个非收敛腿的 RMS 不与收敛腿可比（本项目记过的判据缺陷第 2 条的同族）。
        """
        vals = [sweep[C]["rms_sup"] for C in C_SWEEP if C > 0]
        assert all(a > b for a, b in zip(vals, vals[1:])), (
            f"超声速区误差对 upwind_c 非单调：{['%.3e' % v for v in vals]}")
        shrink = vals[0] / vals[-1]
        assert shrink >= DOSE_MIN, (
            f"C {C_SWEEP[0]} -> {[c for c in C_SWEEP if c>0][-1]} 只收缩 {shrink:.2f}x "
            f"< {DOSE_MIN}（实测 2.66）")

    def test_zero_dissipation_does_not_converge(self, sweep):
        r"""★★ 那份耗散买的是**稳健性** —— 现在它有价格。

        实测：C = 1.5 收敛 11 步；**C = 0 撞 60 步上限、\|R\| = 2.5e-07**。
        ★★ 本条断言的是**权衡仍然存在**。若 C = 0 将来收敛了，那是好消息，
        但**必须走纪律 11 勘误**：`upwind_c = 1.5` 的生产默认（A53 锁着）以这个权衡为理由。
        """
        assert not sweep[0.0]["converged"], (
            f"C = 0 现在收敛了（{sweep[0.0]['reason']}, {sweep[0.0]['n_newton']} 步）—— "
            "权衡变了，upwind_c 的生产默认需要重新论证，请走纪律 11")
        assert sweep[UPWIND_PROD]["converged"], "生产 C 反而不收敛"

    def test_m_crit_is_the_cheaper_lever(self, sweep, prod):
        r"""★ 把 `m_crit` 从 0.95 抬到 1.20（只在 M > 1.2 才上人工密度）：
        超声速区 RMS **1.00e-01 → 2.93e-02**（3.4x），且**仍然收敛**（37 步 vs 11 步）。

        ★★ RECORDED 而非建议改默认：代价是 Newton 步数 3.4 倍，且本门只测了**一个**
        无激波算例 —— 有激波时 `m_crit` 抬高会把激波前的捕捉削弱，那不在本门的域内。
        """
        a = sweep["mcrit120"]
        assert a["converged"], f"m_crit=1.20 腿不收敛（{a['reason']}）"
        assert a["rms_sup"] < prod["rms_sup"], (
            f"m_crit=1.20 的超声速误差 {a['rms_sup']:.3e} 未低于生产 {prod['rms_sup']:.3e}")

    def test_peak_mach_is_undershot_by_the_dissipation(self, prod, sweep):
        r"""★ 生产配置**低估**峰值 Mach（1.2182 vs 精确 1.2464 = −2.26 %），
        而 C = 0 **高估**（1.3102 = +5.1 %）—— 方向相反，这是耗散抹峰的直接读数。"""
        rel = (prod["m_num_max"] - prod["m_ex_max"]) / prod["m_ex_max"]
        assert abs(rel) <= MACH_REL_MAX, (
            f"生产 M_max = {prod['m_num_max']:.4f} vs 精确 {prod['m_ex_max']:.4f}"
            f"（{100*rel:+.2f} % > ±{100*MACH_REL_MAX:.0f} %）")
        assert rel < 0, f"生产配置不再低估峰值 Mach（{100*rel:+.2f} %）—— 耗散的方向变了"


class TestCommittedEvidenceIsLoadBearing:
    r"""★★★ 新鲜计算 vs 已提交 `summary.csv`。设计与三个坑见 `tests/_gate_evidence.py`。
    ★ 零额外计算：用 `sweep` 已算好的五条腿（`upwind_c` 扫描 + `m_crit` 腿）。"""

    def test_fresh_run_reproduces_the_committed_summary(self, sweep):
        fresh = {}
        for C in C_SWEEP:
            d = sweep[C]
            fresh[(f"C={C}",)] = dict(r_exact=d["r_exact"], rms=d["rms"],
                                      rms_sup=d["rms_sup"], rms_sub=d["rms_sub"],
                                      m_num_max=d["m_num_max"])
        d = sweep["mcrit120"]
        fresh[("m_crit=1.20",)] = dict(r_exact=d["r_exact"], rms=d["rms"],
                                       rms_sup=d["rms_sup"], rms_sub=d["rms_sub"],
                                       m_num_max=d["m_num_max"])
        n = assert_matches_committed(
            os.path.join(str(REPO_ROOT), "cases", "gates", "C08_ringleb"), fresh,
            ("r_exact", "rms", "rms_sup", "rms_sub", "m_num_max"),
            key_of=lambda r: (r["leg"],),
            refresh_hint="PYFP3D_GATE_FIGURES=1 pytest tests/C/test_C08_ringleb.py")
        assert n >= 25, f"只比了 {n} 个数（5 腿 x 5 列 = 25）"

@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_ringleb_figure(cases, prod, sweep, gate_evidence_dir):
    import csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c = cases[H_STUDY]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    sc = ax[0].scatter(c["ctr"][:, 0], c["ctr"][:, 1], c=prod["err"], s=4,
                       cmap="RdBu_r", vmin=-0.3, vmax=0.3)
    ax[0].tricontour(c["ctr"][:, 0], c["ctr"][:, 1], c["mach_cell"], levels=[1.0],
                     colors="k", linewidths=1.2)
    fig.colorbar(sc, ax=ax[0], label="speed error (q_num - q_exact)")
    ax[0].set_xlabel("x"), ax[0].set_ylabel("y"), ax[0].set_aspect("equal")
    ax[0].set_title(f"C08 Ringleb, h={H_STUDY}  (black = exact sonic line)")

    cs = [C for C in C_SWEEP if C > 0]
    ax[1].plot(cs, [sweep[C]["rms_sup"] for C in cs], "o-", label="supersonic RMS")
    ax[1].plot(cs, [sweep[C]["rms_sub"] for C in cs], "s-", label="subsonic RMS")
    ax[1].plot([0.0], [sweep[0.0]["rms_sup"]], "rx", ms=9,
               label="C=0 (NOT converged)")
    ax[1].plot([UPWIND_PROD], [sweep["mcrit120"]["rms_sup"]], "g^", ms=8,
               label="m_crit=1.20")
    ax[1].set_xlabel("upwind_c"), ax[1].set_ylabel("speed RMS")
    ax[1].set_title("artificial dissipation, priced against truth")
    ax[1].legend(fontsize=8), ax[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "c08_ringleb.png"), dpi=130)
    plt.close(fig)

    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["leg", "h", "upwind_c", "m_crit", "converged", "reason", "n_newton",
                    "residual", "r_exact", "rms", "rms_sup", "rms_sub",
                    "n_sup", "n_cell", "m_num_max", "m_ex_max"])
        for C in C_SWEEP:
            d = sweep[C]
            w.writerow([f"C={C}", H_STUDY, C, MCRIT_PROD, d["converged"], d["reason"],
                        d["n_newton"], f"{d['res']:.3e}", fmt(d["r_exact"]),
                        fmt(d["rms"]), fmt(d["rms_sup"]), fmt(d["rms_sub"]),
                        d["n_sup"], d["n_cell"], fmt(d["m_num_max"]),
                        fmt(d["m_ex_max"])])
        d = sweep["mcrit120"]
        w.writerow(["m_crit=1.20", H_STUDY, UPWIND_PROD, 1.20, d["converged"], d["reason"],
                    d["n_newton"], f"{d['res']:.3e}", fmt(d["r_exact"]),
                    fmt(d["rms"]), fmt(d["rms_sup"]), fmt(d["rms_sub"]),
                    d["n_sup"], d["n_cell"], fmt(d["m_num_max"]), fmt(d["m_ex_max"])])
