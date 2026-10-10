r"""D14 — **α = 0 跨声速非唯一性：一条能力边界门**（W1.5 / H11，2026-09-06）。

★★★ **它断言的是「仍然坏着」，所以红 = 好消息** —— 形状照抄 D09 / D10。

## 为什么是这道门

`docs/design.md` **§2「What FP cannot do」与 §12 risk 2 都写着**：保守型全速势
在 **M∞ ≈ 0.82–0.85、低升力**上有 Steinhoff–Jameson 非唯一性，登记的缓解措施
逐字是 *"always ramp Mach upward from a subcritical converged state"*。

而 2026-09-05 的门审计量出三件事：

1. **声明包线（M∞ 0.3–0.87）包含那条已登记的带**，而 **D05 的 α = 0 腿只到
   M 0.75** —— **门恰好停在风险开始的地方**；
2. 那条缓解措施**没有被任何代码或门强制**（D05/D11/D12 正是冷启动跑
   M 0.778/0.80/0.803 的）；
3. **实测**：固定网格、固定配方、固定 α = 0、固定 M 0.86，**只改初值**
   （`n_picard_seed`），coarse 上得到**四个各自收敛到 |R| ≈ 3e-13、零钳制的
   不同解**，cl 跨越 **−0.005919 … +0.435982**。

⇒ **判据的形状早就写好了** —— D05 的 `test_zero_incidence_carries_no_lift`
（`|cl| ≤ ZERO_LIFT_CL_MAX = 0.005`）。本门**不改那条判据**，只把它的工况点
往上延伸，并按能力边界的方式记录结果。

## 实测（`bench/studies/gate_audit_20260905/`，α = 0，冷启动，生产配方）

| M∞ | medium cl | \|cl\| / 0.005 | 收敛 |
|---|---|---|---|
| 0.78 | +0.000061 | **0.012** | ✓ |
| 0.80 | −0.006669 | **1.33** | ✓ |
| 0.82 | −0.023454 | **4.69** | ✓ |
| 0.84 | +0.226627 | 45 | ✗ |
| 0.86 | −0.520903 | **104** | ✓（\|R\| 2.6e-13，零钳制） |

⇒ **可信上界落在 M 0.78 与 0.80 之间**，而不是声明的 0.87。
★★★ 最硬的一句不需要「多解」就成立：**medium / M 0.86 / α = 0 上有一个收敛到
|R| = 2.6e-13、零钳制的解，它在一个对称翼型的零迎角上带 cl = −0.52。**
物理上不可能，**而每一个可用性指标都说它是好的** —— GS1.4 的「钳制不静默」
合约在原理上看不见这一类。

## 口径：**哪些数字稳到可以锁，是量出来的**

跨线程 A/B（@8t vs @16t，同一份代码、同一张网格）：

| 腿 | cl 相对差 | 收敛旗标 | 能锁数值? |
|---|---|---|---|
| coarse 全部 5 档 | 7e-12 … 2e-10 | 一致 | ✓ |
| seed 全部 5 条 | 1e-11 … 2e-10 | 一致 | ✓ |
| medium M 0.78 | 5.9e-10 | 一致 | ✓ |
| **medium M 0.80** | **1.1e-01** | **1 → 0（翻转）** | ✗ |
| **medium M 0.82** | 1.8e-05 | 一致 | ✗ |
| **medium M 0.84** | 1.3e-06 | 一致 | ✗ |
| **medium M 0.86** | 1.6e-07 | 一致 | ✗ |

★★★ **不稳的恰好就是落在非唯一带里的那四条 medium 腿** —— 这与分岔的读法
自洽：转折点邻域上 Jacobian 趋于奇异，**线程归约顺序的舍入就足以选中不同的分支
或不同的停滞点**。⇒ 这不是噪声，是**同一个现象的又一个表现**。

⇒ **证据锁只锁稳的那 11 条腿**；不稳的 4 条**写进 CSV 但不进数值锁**，
对它们只断言**不等式**（`|cl|/band > 1`），而那个判定跨线程是稳的
（1.33/1.50 · 4.69/4.69 · 45/45 · 104/104）。
★ 这正是 D05 自己定下的原则 —— **「设门/记录按分歧的大小划分（稳定），
不按旗标（不稳定）」** —— 在本门上的应用；也是 H29（D05 的证据锁继承了线程不稳）
不再重演的做法。
★ **红了怎么办**：不是调阈值。要么**缩包线**（把 M∞ 上界写成实测值），
要么**改激波处理**（Salas–Jameson–Melnik 把根因指向保守位势模型的激波近似）。
★ 熵修正 ON/OFF 在 coarse 与 medium 上**方向相反** ⇒ **不登记方向**。
"""
import os

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting
from tests._gate_evidence import assert_matches_committed, fmt
from tests.conftest import REPO_ROOT, gate_figures_enabled
#: ★ 判据与配方**都从 D05 导入**，不重打一遍 —— 本门的全部意义是「同一条判据，
#: 更高的 M」，抄一份常数就会让两边悄悄分叉。
from tests.D.test_D05_euler_naca0012 import RECIPE, ZERO_LIFT_CL_MAX

MACHS = (0.78, 0.80, 0.82, 0.84, 0.86)
LEVELS = ("coarse", "medium")
SEEDS = (0, 2, 5, 8, 12)
SEED_MACH = 0.86

#: —— 判据（标定，不是实测值）——
LAST_GOOD_MACH = 0.78          # 这一档必须仍然守住 D05 的带；实测 0.012x
FIRST_BAD_MACH = 0.82          # 这一档必须仍然破；实测 4.69x（medium）
BRANCH_SPREAD_MIN = 0.05       # M0.86 coarse 多解的 cl 极差；实测 0.4419
BRANCH_LEGS_MIN = 3            # 收敛且零钳制的腿数；实测 4


def _mesh(level):
    p = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / f"{level}.msh"
    if not p.exists():                       # W0.1 的约定
        pytest.skip(f"naca0012_2.5d/{level}.msh not generated")
    return cut_wake(read_mesh(p))


def _solve(mc, wc, m_inf, seed=None):
    kw = dict(RECIPE)
    if seed is not None:
        kw["n_picard_seed"] = seed
    r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=0.0, **kw)
    dz = float(np.ptp(mc.nodes[:, 2]))
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"],
                                np.asarray(r["phi"]), alpha_deg=0.0, u_inf=1.0,
                                s_ref=dz, m_inf=m_inf)
    return dict(mach=m_inf, cl=float(f["cl"]),
                cl_over_band=abs(float(f["cl"])) / ZERO_LIFT_CL_MAX,
                converged=int(bool(r.get("converged"))),
                residual=float(np.asarray(r["residual_history"], float)[-1]),
                n_limited=int(r.get("n_limited", 0)),
                n_floored=int(r.get("n_floored", 0)))


@pytest.fixture(scope="module")
def sweep():
    out = {}
    for lv in LEVELS:
        mc, wc = _mesh(lv)
        for m in MACHS:
            out[(lv, m)] = _solve(mc, wc, m)
    return out


@pytest.fixture(scope="module")
def seeds():
    mc, wc = _mesh("coarse")
    return {s: _solve(mc, wc, SEED_MACH, seed=s) for s in SEEDS}


class TestTheCapabilityLimitIsRecorded:
    r"""断言「仍然坏着」。**红 = 跑通了 = 好消息** —— 见模块 docstring 的处置。"""

    def test_the_alpha0_band_still_holds_at_the_last_good_mach(self, sweep):
        """M 0.78 上 D05 的带**仍然守得住** —— 这是包线的实测下沿。

        ★ 它红了才是真的坏消息：说明连 0.78 都保不住。
        """
        for lv in LEVELS:
            d = sweep[(lv, LAST_GOOD_MACH)]
            assert d["converged"], f"{lv} M{LAST_GOOD_MACH} 不收敛，前提变了"
            assert d["cl_over_band"] < 1.0, (
                f"{lv} M{LAST_GOOD_MACH}: |cl| = {abs(d['cl']):.6f} 已越过 D05 的带 "
                f"{ZERO_LIFT_CL_MAX} —— 可信上界比 0.78 还低了")

    def test_the_alpha0_band_is_still_violated_above(self, sweep):
        """M ≥ 0.82 上带**仍然破着**（medium）。

        ★★ **红了意味着改进**：断言消息写明下一步 —— 这不是回归，
        是「非唯一性被治好了或被绕开了」，那时本门要**重新定型为一条正常的
        零升力门**（判据现成，就是 D05 的那条）。
        """
        for m in [x for x in MACHS if x >= FIRST_BAD_MACH]:
            d = sweep[("medium", m)]
            assert d["cl_over_band"] > 1.0, (
                f"medium M{m}: |cl| = {abs(d['cl']):.6f} 现在落回 D05 的带 "
                f"{ZERO_LIFT_CL_MAX} 之内 —— **this is an improvement, not a "
                f"regression**: re-specify D14 as a plain zero-lift gate "
                f"(D05's criterion, extended in Mach) and shrink this file to "
                f"the legs that still fail")

    def test_a_converged_zero_clamp_solution_carries_impossible_lift(self, sweep):
        """最硬的一条：M 0.86 medium 上有一个**干净收敛**的解带着 |cl| ~ 0.5。

        ★ 「干净」= converged 且 `n_limited = n_floored = 0` ⇒ GS1.4 的
        clamp-not-silent 合约**看不见它**，这正是本门存在的理由。
        """
        d = sweep[("medium", 0.86)]
        assert d["converged"] and d["n_limited"] == 0 and d["n_floored"] == 0, (
            f"前提变了：M0.86 medium 不再是干净收敛 "
            f"(conv={d['converged']}, lim/flr={d['n_limited']}/{d['n_floored']})")
        assert d["residual"] < 1e-10, f"|R| = {d['residual']:.2e} 不再是机器零"
        assert abs(d["cl"]) > 0.1, (
            f"零迎角对称翼型的 |cl| 降到 {abs(d['cl']):.4f} —— 好消息，重新定型本门")

    def test_the_solution_is_not_unique(self, seeds):
        """**只改初值**就得到不同的解 —— 非唯一性本身。

        ★ 只数**收敛且零钳制**的腿；未收敛的腿不参与（本门不读 `accept_reason`，
        所以不给它们命名模式）。
        ★ G-DOMAIN：相反结果 = 各腿一致 ⇒ 该 M 上解唯一 ⇒ 红，且**是好消息**。
        """
        clean = [d for d in seeds.values()
                 if d["converged"] and d["n_limited"] == 0 == d["n_floored"]]
        assert len(clean) >= BRANCH_LEGS_MIN, (
            f"只有 {len(clean)} 条干净收敛的腿（需要 {BRANCH_LEGS_MIN}）—— "
            "前提变了，不能据此判唯一性")
        cl = np.array([d["cl"] for d in clean])
        spread = float(cl.max() - cl.min())
        assert spread >= BRANCH_SPREAD_MIN, (
            f"{len(clean)} 条干净收敛的腿 cl 极差只有 {spread:.6f} "
            f"(< {BRANCH_SPREAD_MIN}) —— 解看起来唯一了。**this is an "
            f"improvement**: 处置是缩包线或改激波处理带来的效果，"
            f"重新定型本门；**不要调这个阈值**")


#: ★★★ 只锁**跨线程稳定**的腿（模块 docstring 的表）。medium M >= 0.80 的四条
#: 落在非唯一带里，线程归约顺序就能移动它们 —— 锁它们等于把环境焊进门里。
_LOCKABLE_SWEEP = [("coarse", m) for m in MACHS] + [("medium", 0.78)]
_UNSTABLE_SWEEP = [("medium", m) for m in MACHS if m >= 0.80]


class TestCommittedEvidenceIsLoadBearing:
    r"""新鲜计算 vs 已提交 `summary.csv`，**只对稳定腿**（口径见模块 docstring）。"""

    COLS = ("cl", "cl_over_band", "residual")

    def test_the_unstable_legs_are_excluded_for_a_measured_reason(self):
        """前提断言：被排除的正好是那四条，且它们**确实**在非唯一带里。

        ★ 它红了说明有人改了排除名单而没重做跨线程 A/B。
        """
        assert _UNSTABLE_SWEEP == [("medium", m) for m in (0.80, 0.82, 0.84, 0.86)]
        assert all(m >= FIRST_BAD_MACH - 0.02 for _lv, m in _UNSTABLE_SWEEP)
        assert set(_LOCKABLE_SWEEP) | set(_UNSTABLE_SWEEP) == {
            (lv, m) for lv in LEVELS for m in MACHS}, "腿的划分不完整"

    def test_matches_committed_summary(self, sweep, seeds, gate_evidence_dir):
        fresh = {f"sweep|{lv}|{m}": {k: fmt(sweep[(lv, m)][k]) for k in self.COLS}
                 for lv, m in _LOCKABLE_SWEEP}
        fresh.update({f"seed|{s}|{SEED_MACH}":
                      {k: fmt(seeds[s][k]) for k in self.COLS} for s in SEEDS})
        n = assert_matches_committed(
            gate_evidence_dir, fresh, self.COLS,
            key_of=lambda r: f"{r['arm']}|{r['key']}|{r['mach']}",
            refresh_hint="PYFP3D_GATE_FIGURES=1 pytest "
                         "tests/D/test_D14_alpha0_transonic_nonuniqueness.py")
        assert n >= 33, f"只比了 {n} 个数（11 条稳定腿 x 3 列 = 33）"

    def test_the_unstable_legs_still_give_a_stable_VERDICT(self, sweep):
        """不锁数值，但**判定**必须稳：那四条跨线程都远在带外（1.3x … 104x）。"""
        for lv, m in _UNSTABLE_SWEEP:
            assert sweep[(lv, m)]["cl_over_band"] > 1.0, (
                f"{lv} M{m}: |cl|/band = {sweep[(lv, m)]['cl_over_band']:.2f} "
                "—— 判定本身翻了，见 test_the_alpha0_band_is_still_violated_above")


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图/CSV 证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_evidence(sweep, seeds, gate_evidence_dir):
    """★ D09/D10 **不写证据 CSV**（审计的 H22）。本门是新建的，所以不重复那个形状。"""
    import csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(str(gate_evidence_dir), exist_ok=True)
    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["arm", "key", "mach", "cl", "cl_over_band", "converged",
                    "residual", "n_limited", "n_floored"])
        for lv in LEVELS:
            for m in MACHS:
                d = sweep[(lv, m)]
                w.writerow(["sweep", lv, m, fmt(d["cl"]), fmt(d["cl_over_band"]),
                            d["converged"], fmt(d["residual"]),
                            d["n_limited"], d["n_floored"]])
        for s in SEEDS:
            d = seeds[s]
            w.writerow(["seed", s, SEED_MACH, fmt(d["cl"]),
                        fmt(d["cl_over_band"]), d["converged"],
                        fmt(d["residual"]), d["n_limited"], d["n_floored"]])

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for lv, mk in zip(LEVELS, ("o-", "s-")):
        ax[0].semilogy([m for m in MACHS],
                       [sweep[(lv, m)]["cl_over_band"] for m in MACHS], mk,
                       label=lv)
    ax[0].axhline(1.0, color="r", ls="--", lw=1, label="D05 band |cl| = 0.005")
    ax[0].set_xlabel("M_inf"), ax[0].set_ylabel("|cl| / band")
    ax[0].set_title("D14  alpha = 0, symmetric section: spurious lift")
    ax[0].legend(fontsize=8), ax[0].grid(alpha=.3, which="both")

    clean = [(s, seeds[s]) for s in SEEDS
             if seeds[s]["converged"] and seeds[s]["n_limited"] == 0
             == seeds[s]["n_floored"]]
    ax[1].bar([str(s) for s, _ in clean], [d["cl"] for _, d in clean],
              color="tab:red")
    ax[1].axhline(0.0, color="k", lw=.8)
    ax[1].set_xlabel("n_picard_seed  (ONLY the initial guess changes)")
    ax[1].set_ylabel("cl")
    ax[1].set_title(f"coarse, M {SEED_MACH}, alpha = 0: converged, zero-clamp\n"
                    f"solutions of the SAME discrete problem")
    ax[1].grid(alpha=.3, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir),
                             "d14_alpha0_nonuniqueness.png"), dpi=130)
    plt.close(fig)
