r"""D05 — pyFP3D 无粘 vs CFL3D Euler — NACA0012（2D-1..2D-4）。

工况（`cases/reference_data/cfl3d/euler_naca0012/`，6 工况 × 3 档）：
M0.50/α2.0 · **M0.80/α1.25 = M1 靶子本身** · M0.72 与 M0.75/α0 ·
M0.778/α2.03 · M0.803/α−0.1。参考侧的档 **L1–L4**，本门读**最细的 L4**
（2026-09-05 参考补第四档；★ 力系数与 Cp 必须读**同一档**，否则是跨档比较）。

---

## ★★★ 判据是在**测量之后**登记的 —— 而第一次测量是错的

先测再登记。但第一次测量在 **8 线程**上做，而项目纪律 #1 是**上限 16，含
BLAS/OMP**。由此得出的"跨声速举力腿不能设门"是**测量条件的产物**：

| 算例（medium） | pyFP3D cl | CFL3D cl | Δcl | @8t 收敛 | **@16t 收敛** |
|---|---|---|---|---|---|
| M0.50/α2.0 | 0.284372 | 0.284256 | **+0.04 %** | ✓ | ✓ |
| M0.72/α0 | 0.001164 | 0.000000 | — | ✓ | ✓ |
| M0.75/α0 | 0.001273 | 0.000000 | — | ✓ | ✓ |
| **M0.80/α1.25** | 0.341340 | 0.348239 | **−1.98 %** | ✗ | **✓** |
| M0.778/α2.03 | 0.706197 | 0.490875 | +43.9 % | ✗ | ✗ |
| **M0.803/α−0.10** | −0.042311 | −0.028789 | −47 % | **✓** | **✗** |

★★★ **收敛旗标在两个方向上都随线程数翻转，而 cl 五位不变。**
M0.80 在 8/12 线程 80 步封顶（|R| ~1e-06）、16 线程收敛到 |R| 2.9e-13，cl 差
0.36 %；M0.803 反过来，8 线程收敛而 16 线程不收敛，cl 五位相同
（−0.042314 / −0.042311）。⇒ **旗标是稳定答案周围的线程依赖噪声。**

⇒ **设门/记录的划分依据必须是分歧的大小（稳定），不是收敛旗标（不稳定）。**
本门**不对 `converged` 设任何断言**，只把它连同线程数写进证据。
（项目里 `test_p4_transonic` 正因线程依赖才设成非严格 xfail —— 同一现象，
这里是从零复现了一遍。★ 第一版 D06 曾写下 `assert not converged`，
把一个环境依赖的结果、而且是**坏的那一侧**焊进门里；**是那条断言自己红了**
才把问题揪出来，不是复核。）

★★ **近零升力必须用绝对量**。M0.803/α−0.1 的参考 cl 是 **−0.0288**，
所以"−47 %"的分母近零 —— 绝对差只有 **0.0135**。这与 α = 0 两例是同一条：
**相对还是绝对，取决于基线，不存在一个总是对的选择。**

★ 一条被测量否掉的预期：我按教科书预期"全位势的激波应落在 Euler **下游**"。
NACA0012 上**实测是上游**（4/5），而 RAE2822 上是**下游**（D06）——
**没有一致方向，所以没有可登记的符号。**

★★★ **M0.75 的激波位置对这份参考不可分辨。** pyFP3D 与 CFL3D 差 **0.0047**，
而 CFL3D 自身的 L2→L3 差是 **0.018026**（且**非单调**：0.285 → 0.278 → 0.296），
差比参照物自身的不确定度小 **3.8 倍** ⇒ 与 D07 对实验那次同一条判据。
★ 对照 D06：RAE2822 的参考紧到 0.00088，同样的比较在那里**可分辨**（35×）——
**参照物的分辨率是逐例的性质，不是数据集的性质。**

## ★ 这份参考自己的短板，如实记下

2-D 参考表只有 `delta_L2_L3`，**没有** 3-D 那套"相邻差是否在缩小"的渐近检验
（D07 上第四档把 cl 从"收敛得极好"翻成"不渐近"）。所以这里的 `delta_L2_L3`
是**最后一步的差**，不是经检验的误差棒。M0.50 的 `cd` 相对差 **−48 %** 就是例子。
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
from tests.D._cfl3d_cp import cp_rms, read_2d_cp
from tests._gate_evidence import assert_matches_committed, fmt
from tests.conftest import REPO_ROOT, gate_figures_enabled

REF_DIR = os.path.join(str(REPO_ROOT), "cases", "reference_data", "cfl3d",
                       "euler_naca0012")
#: 生产配方，与 `bench/run_m1_gate.py` / D11 逐字一致
RECIPE = dict(upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80,
              n_picard_seed=5)

#: (名字, case-id, M, α, 这条腿是设门还是只记录)
CASES = (
    ("M0.50_a2.00", "n0012_m0500_a2.00", 0.50, 2.00, "gate"),
    ("M0.72_a0.00", "n0012_m0720_a0.00", 0.72, 0.00, "gate"),
    ("M0.75_a0.00", "n0012_m0750_a0.00", 0.75, 0.00, "gate"),
    ("M0.80_a1.25", "n0012_m0800_a1.25", 0.80, 1.25, "gate"),
    ("M0.778_a2.03", "n0012_m0778_a2.03", 0.778, 2.03, "record"),
    ("M0.803_a-0.10", "n0012_m0803_am0.10", 0.803, -0.10, "record"),
)
LEVELS = ("coarse", "medium")

# —— 判据（标定，不是实测值；实测基线见 docstring 的表）——
SUBCRIT_CL_MAX = 0.010      # 亚临界 |Δcl|/cl：实测 medium 0.0004
SUBCRIT_CONTRACT = 2.0      # coarse -> medium 必须收缩：实测 0.0235/0.0004 = 57x
ZERO_LIFT_CL_MAX = 0.005    # α=0 的 |cl| 绝对带：实测 0.0012 / 0.0013
TRANSONIC_CL_MAX = 0.05     # M0.80/α1.25（M1 靶子）|Δcl|/cl：实测 0.0198
#: ★ 近零升力用**绝对**带（参考 cl = −0.0288，相对判据的分母近零）
NEAR_ZERO_ABS_MAX = 0.05    # M0.803/α−0.1 |Δcl| 绝对：实测 0.0135
#: ★ 与参考自身不确定度比较时的"可分辨"下限，形式与 D07 的 RESOLVE_FRAC 相同：
#: 差必须超过参照物自身 L2→L3 差的这个倍数才谈得上比较。
RESOLVABLE_FRAC = 1.0
#: —— W1.3 / H14 的三条（标定；依据见 TestDAlembertOnTheProductionPath 的表）——
#: ★ 网格已提交 ⇒ 无网格实现噪声，带可以比 C06 紧。
DALEMBERT_CD_MAX = 3.0e-3      # 无激波腿 medium |cd|：实测 1.564e-3 / 1.266e-3
DALEMBERT_CONTRACTION = 2.0    # coarse -> medium：实测 3.18x / 3.23x
WAVE_DRAG_MIN = 5.0e-3         # 带激波腿 medium |cd| 下界：实测 7.6e-3 … 3.8e-2


def _read_reference(path=None):
    """从**已提交的 CSV** 读 CFL3D 参考，**L4（最细）**档。

    ★★★ 这是 F06 那次的教训：门里写死一个字面量、再拿它和另一个字面量断言，
    改参考文件什么都不会红。本函数**真的解析文件**，并且
    `test_gate_actually_reads_the_reference` 拿一份**被扰动的副本**喂给它、
    要求读数跟着走 —— 行为验证，不是"值等于它自己的来源"。
    """
    import csv
    p = path or os.path.join(REF_DIR, "forces.csv")
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            #: ★ 读**最细**档。2026-09-05 起参考有 L4；继续读 L3 会把
            #:   pyFP3D 与一个被参考自己取代了的档位相比。
            if r["level"] == "L4":
                out[r["case"]] = dict(cl=float(r["cl"]), cd=float(r["cd"]),
                                      cm=float(r["cm_quarter_chord"]))
    if not out:
        raise RuntimeError(f"{p}: no L4 rows -- the reference layout changed")
    return out


def _read_reference_delta(path=None):
    """参考自身的 L2→L3 差 **和它的渐近判定**。

    ★★★ 2026-09-05 起参考数据集自己带 `ratio` / `asymptotic` / `error_bar`。
    本函数一并读出来，因为**"最后一步的差"只有在差在缩小时才是误差棒**，
    而这道门用到的那一个恰恰不是（见下）。
    """
    import csv
    p = path or os.path.join(REF_DIR, "grid_convergence.csv")
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            vals = [float(r[lv]) for lv in ("L1", "L2", "L3", "L4")
                    if r.get(lv)]
            out[(r["case"], r["quantity"])] = dict(
                #: ★★★ **散布**（各档极差），不是单一的最后差。
                #: 对一个判定不稳的量，"最后一步的差"可能恰好很小而该量仍在
                #: 各档间游走；**它游走的范围才是参照物的不确定度**，
                #: 而且无论收敛与否都有定义。
                scatter=(max(vals) - min(vals)) if len(vals) >= 2 else float("nan"),
                delta=abs(float(r["delta_L2_L3"])) if r["delta_L2_L3"] else float("nan"),
                asymptotic=r["asymptotic"].split(" (")[0])
    return out


def _one(level, m_inf, alpha):
    mc, wc = cut_wake(read_mesh(os.path.join(
        str(REPO_ROOT), "cases", "meshes", "naca0012_2.5d", f"{level}.msh")))
    dz = float(np.ptp(mc.nodes[:, 2]))
    r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=alpha, **RECIPE)
    phi = np.asarray(r["phi"])
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], phi,
                                alpha_deg=alpha, u_inf=1.0, s_ref=dz,
                                m_inf=m_inf)
    cur = section_cp_curve(mc, phi, z=float(np.mean(mc.nodes[:, 2])),
                           smooth_passes=1, m_inf=m_inf)
    up = shock_report(cur, m_inf)["upper"]
    #: ★ `sigma_freeze_report` 是已知缺陷的**报告**通道；把它带进证据里，
    #:   因为跨声速腿全部命中它，而那正是这些腿只能 RECORDED 的理由之一。
    frz = r.get("sigma_freeze_report") or {}
    return dict(
        # ★ 键名是 `cd_pressure`，不是 `cd` —— 对无粘的 CFL3D 侧
        #   cd == cd_pressure（cd_friction 恒 0），所以这是同名量比较。
        cl=float(f["cl"]), cd=float(f["cd_pressure"]),
        x_shock=(float(up["x_shock"]) if up.get("has_shock") else float("nan")),
        converged=bool(r.get("converged")),
        residual=float(np.asarray(r["residual_history"], float)[-1]),
        n_limited=int(r.get("n_limited", 0)),
        n_floored=int(r.get("n_floored", 0)),
        # ★ 真实键名（读出来的，不是回忆的）：frozen_in_transient /
        #   selection_churn / tau_settled / last_sigma_delta / tail_n_shock
        sigma_frozen=bool(frz.get("frozen_in_transient", False)),
        #: ★ Cp 曲线带出来，因为图必须是 Cp 分布而不是力系数柱状图
        #:   （使用者裁决 2026-09-05）：一个 cl 是一个数，Cp 显示差在哪里。
        curve=cur,
        sigma_churn=bool(frz.get("selection_churn", False)),
    )


def _with_cp_rms(d, cid):
    """★★ Cp RMS 是**主要**比较量（使用者裁决 2026-09-05），所以它必须进证据
    CSV 并被回归锁覆盖 —— 只画在 PNG 里的量，改一行代码没有任何东西会红。"""
    rc = read_2d_cp(REF_DIR, cid, "L4", "none")
    for side in ("upper", "lower"):
        if side in rc:
            d[f"cp_rms_{side}"] = cp_rms(
                rc[side][0], rc[side][1],
                d["curve"][f"x_{side}"], d["curve"][f"cp_{side}"])[0]
        else:
            d[f"cp_rms_{side}"] = float("nan")
    return d


@pytest.fixture(scope="module")
def runs():
    return {(nm, lv): _with_cp_rms(_one(lv, m, a), cid)
            for nm, cid, m, a, _k in CASES for lv in LEVELS}


@pytest.fixture(scope="module")
def ref():
    return _read_reference()


class TestWhatIsGateable:
    r"""★ 只有**无粘全位势确实负责**的那一格：亚临界升力，和 α=0 的对称性。"""

    def test_subcritical_lift_matches_cfl3d_euler(self, runs, ref):
        """M0.50/α2.0：全位势与 Euler 在亚临界下应当一致，因为两者的模型差
        （熵、涡量）在无激波时消失。实测 medium 0.04 %。"""
        got = runs[("M0.50_a2.00", "medium")]["cl"]
        want = ref["n0012_m0500_a2.00"]["cl"]
        rel = abs(got - want) / abs(want)
        assert runs[("M0.50_a2.00", "medium")]["converged"]
        assert rel <= SUBCRIT_CL_MAX, (
            f"subcritical cl {got:.6f} vs CFL3D Euler {want:.6f} "
            f"= {rel*100:.3f} % > {SUBCRIT_CL_MAX*100:.1f} %")

    def test_subcritical_lift_contracts_with_refinement(self, runs, ref):
        """★ 一个不收缩的一致只是巧合。coarse 2.35 % -> medium 0.04 %。"""
        want = ref["n0012_m0500_a2.00"]["cl"]
        e = {lv: abs(runs[("M0.50_a2.00", lv)]["cl"] - want) / abs(want)
             for lv in LEVELS}
        assert e["coarse"] / e["medium"] >= SUBCRIT_CONTRACT, (
            f"subcritical cl error must contract: coarse {e['coarse']*100:.3f} % "
            f"-> medium {e['medium']*100:.3f} % is only "
            f"{e['coarse']/max(e['medium'],1e-15):.2f}x")

    @pytest.mark.parametrize("name", ["M0.72_a0.00", "M0.75_a0.00"])
    def test_zero_incidence_carries_no_lift(self, runs, name):
        """★★ **绝对**带，不是相对：两侧 cl 都 ≈ 0，相对判据要除以零。"""
        got = runs[(name, "medium")]["cl"]
        assert runs[(name, "medium")]["converged"]
        assert abs(got) <= ZERO_LIFT_CL_MAX, (
            f"{name}: symmetric section at alpha 0 gives cl {got:+.6f}, "
            f"outside the absolute band {ZERO_LIFT_CL_MAX}")


    def test_m1_target_condition_lift(self, runs, ref):
        """★★ M0.80/α1.25 = **M1 靶子本身**。@16 线程收敛，Δcl −1.98 %。

        ★ 本条**不断言 `converged`** —— 那个旗标在 8/12 线程是 False、16 线程
        是 True，而 cl 只差 0.36 %。断言旗标就是把线程数焊进门里。"""
        got = runs[("M0.80_a1.25", "medium")]["cl"]
        want = ref["n0012_m0800_a1.25"]["cl"]
        rel = abs(got - want) / abs(want)
        assert rel <= TRANSONIC_CL_MAX, (
            f"M1-target cl {got:.6f} vs CFL3D Euler {want:.6f} = "
            f"{rel*100:.2f} % > {TRANSONIC_CL_MAX*100:.0f} %")


class TestDAlembertOnTheProductionPath:
    r"""★★★ W1.3 / H14（2026-09-06）：**杂散阻力 —— 生产 Newton 路径上的免费误差范数。**

    审计实测：`cd` / `cd_pressure` 在 D05 D06 D07 D08 D10 **五个门里被读取**、
    写进证据 CSV，**而没有一条 `assert` 碰过它**。本条是第一条。

    判据不需要任何参考数据：**无激波的无粘二维解，压差阻力精确为 0**
    （d'Alembert；二维无诱导阻力，无激波则无波阻）。⇒ `|cd|` 就是误差本身。

    ★★★ **判据自己划定作用域，靠的是实测量而不是假设**：用 `x_shock` 是否为
    NaN 把「无激波腿」挑出来 —— 不是按 M∞ 猜。实测结果正好说明为什么必须这样：

    | 腿 | x_shock | cd coarse → medium | 比值 |
    |---|---|---|---|
    | M0.50/α2.0 | — | 4.970e-03 → 1.564e-03 | **3.18x** |
    | M0.72/α0 | — | 4.082e-03 → 1.266e-03 | **3.23x** |
    | M0.75/α0 | 0.292 / 0.300 | 4.115e-03 → 1.311e-03 | 3.14x |
    | M0.80/α1.25 | 0.614 / 0.599 | 1.985e-02 → 1.818e-02 | 1.09x |
    | M0.778/α2.03 | 0.596 / 0.702 | 2.592e-02 → 3.799e-02 | **0.68x** |
    | M0.803/α−0.1 | 0.462 / 0.490 | 7.829e-03 → 7.625e-03 | 1.03x |

    ⇒ **两族的行为截然不同**，而分界线是激波而不是马赫数：跨声速腿的 `cd`
    **不收缩**（甚至增长），因为那是**真实的波阻**，不是误差。
    ★ M0.75/α0 检测到一道弱激波但 `cd` 仍按无激波的方式收缩 ⇒ 那道激波的波阻
    **低于离散误差**。它**被排除在设门之外**（它的 x_shock 非 NaN），
    只作为记录 —— 一个量的作用域要按判别式划，不按「看起来像」划。

    ★ 这里的网格是**已提交**的（`naca0012_2.5d`），所以没有 C06 那种
    网格实现噪声的问题，带可以比 C06 紧。
    """

    #: 无激波腿由 x_shock 为 NaN **判别**出来，不是按名字写死
    def _shock_free(self, runs):
        return [nm for nm, _c, _m, _a, _k in CASES
                if all(not np.isfinite(runs[(nm, lv)]["x_shock"])
                       for lv in LEVELS)]

    def test_the_shock_free_set_is_what_we_think_it_is(self, runs):
        """前提断言：判别式必须挑出恰好那两条腿。"""
        got = self._shock_free(runs)
        assert got == ["M0.50_a2.00", "M0.72_a0.00"], (
            f"shock-free legs moved: {got} -- the criterion's DOMAIN changed, "
            "re-read the table in this class's docstring before touching bands")

    def test_shock_free_pressure_drag_approaches_zero(self, runs):
        for nm in self._shock_free(runs):
            got = abs(runs[(nm, "medium")]["cd"])
            assert got <= DALEMBERT_CD_MAX, (
                f"{nm}: spurious cd {got:.4e} > {DALEMBERT_CD_MAX:.1e} "
                f"-- exact is 0 (no shock, 2-D, inviscid)")

    def test_shock_free_pressure_drag_contracts(self, runs):
        for nm in self._shock_free(runs):
            c = abs(runs[(nm, "coarse")]["cd"])
            m = abs(runs[(nm, "medium")]["cd"])
            assert c / m >= DALEMBERT_CONTRACTION, (
                f"{nm}: spurious cd contracts only {c / m:.2f}x "
                f"(< {DALEMBERT_CONTRACTION}) -- {c:.3e} -> {m:.3e}")

    def test_the_transonic_legs_are_outside_this_criterion(self, runs):
        """★★ G-DOMAIN **被断言出来，不是被声明出来**：带激波的腿 `cd` 既大
        又不收缩 ⇒ 那是真实波阻，把 d'Alembert 判据套上去就是**把物理当缺陷**。
        这一条红了意味着作用域变了，要重读上表。"""
        for nm, _c, _m, _a, _k in CASES:
            if nm in self._shock_free(runs):
                continue
            if nm == "M0.75_a0.00":
                continue          # 弱激波，波阻低于离散误差 —— 见 docstring
            c = abs(runs[(nm, "coarse")]["cd"])
            m = abs(runs[(nm, "medium")]["cd"])
            assert m >= WAVE_DRAG_MIN, (
                f"{nm}: cd {m:.3e} < {WAVE_DRAG_MIN:.0e} -- a shocked leg with "
                "no wave drag; the shock-free discriminator may be wrong")
            assert c / m < DALEMBERT_CONTRACTION, (
                f"{nm}: shocked-leg cd contracts {c / m:.2f}x like a pure error "
                "norm -- if that is real, the wave drag has vanished")


class TestWhatIsOnlyRecorded:
    r"""★★ 跨声速举力腿：两条不收敛、五条全部熵修正冻结。**记录，不设门。**"""

    @pytest.mark.parametrize("name", ["M0.778_a2.03", "M0.803_a-0.10"])
    def test_transonic_legs_are_recorded_with_their_solver_state(self, runs,
                                                                 ref, name):
        """断言的不是"值对"，而是**状态被如实带出来了** —— 一条腿的
        收敛旗标、限幅计数和 sigma 冻结状态必须存在且可读。

        ★ 这道断言防的是"悄悄把一条没收敛的腿当成结果发表"，
        它正是本轮在 CFL3D 侧踩过的那个缺陷（SUMMARY 存在 ≠ 跑完了）。"""
        r = runs[(name, "medium")]
        #: ★ 断言的是**状态可读**，不是状态的值 —— `converged` 随线程数翻转。
        assert isinstance(r["converged"], bool)
        assert np.isfinite(r["residual"])
        assert isinstance(r["sigma_frozen"], bool)
        assert isinstance(r["sigma_churn"], bool)

    def test_near_zero_lift_uses_an_absolute_band(self, runs, ref):
        """★★ M0.803/α−0.1 的"−47 %"分母近零（参考 cl −0.0288）。
        绝对差 0.0135 —— 大，但只有相对数字的三分之一那么骇人。
        **RECORDED**：它超出 M0.50 的绝对误差 100 倍，是真实分歧。"""
        got = runs[("M0.803_a-0.10", "medium")]["cl"]
        want = ref["n0012_m0803_am0.10"]["cl"]
        assert abs(got - want) <= NEAR_ZERO_ABS_MAX, (
            f"near-zero-lift |dcl| {abs(got-want):.6f} exceeds the absolute "
            f"band {NEAR_ZERO_ABS_MAX}")

    def test_m075_shock_is_unresolvable_against_this_reference(self, runs):
        """★★★ pyFP3D 与 CFL3D 差 0.0047，而 CFL3D 自身的 L2→L3 差是 0.018。

        差比参照物自身的散布小 3.8 倍 ⇒ **不可分辨**。

        ★★ **两次更正 2026-09-05**。(1) 这条原来把单一的 `delta_L2_L3`
        （0.018）叫作参考的"不确定度"。(2) 加了第四档后该量的判定是
        **unstable** —— 三档三元组给 ratio 2.483、四档三元组给 0.217，
        **多一档就翻**，所以这条阶梯**决定不了**它是否收敛。

        ⇒ 改用**散布**（各档极差）：参考自身四档读
        **0.28502 / 0.27776 / 0.29578 / 0.29970**，极差 **0.021941**，
        而且非单调。pyFP3D 与它差 0.0047 = 散布的 **0.21 倍** ⇒ **不可分辨**，
        而且理由比原来更硬：参照物**本身就在 0.022 弦长的范围里游走**。

        ★ 本条断言那个"定不了"的状态仍然如此 —— 若参考哪天在这条阶梯上
        收敛了，才谈得上写一条真正的激波位置判据。"""
        d = _read_reference_delta()
        ref = d[("n0012_m0750_a0.00", "x_shock_upper")]
        assert ref["asymptotic"] in ("no", "unstable"), (
            f"the reference's own M0.75 shock is now {ref['asymptotic']!r} -- "
            f"its convergence has become decidable on this ladder, so the "
            f"'doubly unresolvable' reading no longer applies and this row "
            f"should be re-specified")
        unc = ref["scatter"]
        got = runs[("M0.75_a0.00", "medium")]["x_shock"]
        import csv
        with open(os.path.join(REF_DIR, "shock.csv")) as fh:
            want = next(float(r["x_shock"]) for r in csv.DictReader(fh)
                        if r["case"] == "n0012_m0750_a0.00"
                        and r["level"] == "L4" and r["surface"] == "upper")
        diff = abs(got - want)
        assert diff < RESOLVABLE_FRAC * unc, (
            f"the pyFP3D-vs-CFL3D shock difference {diff:.6f} now EXCEEDS the "
            f"reference's own across-rung SCATTER {unc:.6f} -- it has become "
            f"resolvable, so this RECORDED row should be re-specified as a "
            f"real criterion")


class TestReferenceIsLoadBearing:
    r"""★★★ F06 的教训：门必须**真的读**参考文件，并且**行为验证**这一点。"""

    def test_gate_actually_reads_the_reference(self, tmp_path):
        """把一份**被扰动**的副本喂给读取器，要求读数跟着走。

        ★ 只断言"读到的值 == 文件里的值"是按构造恒真的（F06 实测：把参考改成
        0.64，全套依然绿）。必须让扰动**传播**。"""
        import csv
        src = os.path.join(REF_DIR, "forces.csv")
        dst = tmp_path / "forces.csv"
        with open(src) as fh:
            rows = list(csv.DictReader(fh))
            cols = rows[0].keys()
        for r in rows:
            if r["case"] == "n0012_m0500_a2.00" and r["level"] == "L4":
                r["cl"] = f"{float(r['cl']) + 0.1234:.6f}"
        with open(dst, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(cols))
            w.writeheader()
            w.writerows(rows)
        base = _read_reference()["n0012_m0500_a2.00"]["cl"]
        pert = _read_reference(str(dst))["n0012_m0500_a2.00"]["cl"]
        assert abs((pert - base) - 0.1234) < 1e-9, (
            "the reader did not follow a perturbed reference -- the gate is "
            "not actually reading the committed file")


class TestCommittedEvidenceIsLoadBearing:
    r"""★★★ 新鲜计算 vs 已提交 `summary.csv`。设计见 `tests/_gate_evidence.py`。

    ★★★ **这条锁要在 16 线程下跑**（2026-09-06 实测，W1.3 顺带查出）。
    已提交 CSV 的 `n_threads` 列写着 **16**，而本门 docstring 第 1 节记着
    M0.80 那条腿**在 8 线程上 80 步封顶、16 线程上收敛到 |R| 2.9e-13，cl 差
    0.36 %**。⇒ 8 线程跑这条锁**必红，而那是环境红不是回归红**：
    实测同一份代码 `@8t` FAIL / `@16t` PASS。

    ★ **登记未修**（H29）：本门刻意「不对 `converged` 设任何断言」，因为旗标不稳；
    但**证据锁本身仍然继承了那份不稳定**。修法（属于 W2 而不是 W1）是让锁按
    `n_threads` 分档、或对那条腿放宽到跨线程分歧之上 —— 两者都动既有判据，
    要走预注册。在此之前：**跑 D 类锁请用 16 线程**（纪律 #1 的上限，也是本 CSV 的口径）。
    
    ★★★ **`residual` 退出数值锁 —— 2026-09-06 修掉的一处判据缺陷**（W1 收口实测）。

    残差是**解算过程的诊断**，不是**答案**：在机器零附近它的末位就是舍入本身。
    逐列实测（8 线程 vs 已提交的 16 线程证据，**同一份代码**）：
    D05 的 `residual` 最坏相对差 **3.47e+06**、D06 **2.97e-01**、D07 **1.55e-02**
    —— 而 D07 的**十个物理列全部 0.00e+00，逐位相同**。
    D07 那一行最说明问题：**7.285e-15 对 7.174e-15 被算成「相对差 1.55e-02」，
    绝对差却只有 1e-16**。用 `rel_tol = 1e-6` 锁这种量，锁的是舍入对舍入。

    ⇒ 该列**留在证据 CSV 里**（它是有用的记录），但**不进数值锁**；有意义的断言是
    **量级**（「这条腿收敛了」），那由本门自己的收敛/状态断言负责。
    ★ 这**不是放宽容差**：`tests/_gate_evidence.py` 明写「若容差需要放宽，那本身就是
    一个要查的信号」—— 查的结果是**这一列从一开始就不该进值锁**。
    """

    #: ★ cp_rms_* 在前：它是主要比较量
    MEASURED = ("cp_rms_upper", "cp_rms_lower", "cl", "cd", "x_shock")

    def test_matches_committed_summary(self, runs, gate_evidence_dir):
        fresh = {f"{nm}|{lv}": {k: fmt(runs[(nm, lv)][k]) for k in self.MEASURED}
                 for nm, _c, _m, _a, _k in CASES for lv in LEVELS}
        assert_matches_committed(
            gate_evidence_dir, fresh, self.MEASURED,
            key_of=lambda r: f"{r['case']}|{r['level']}",
            refresh_hint="PYFP3D_GATE_FIGURES=1 pytest "
                         "tests/D/test_D05_euler_naca0012.py")


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图/CSV 证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_evidence(runs, ref, gate_evidence_dir):
    """写 `summary.csv` + 一张对照图。

    ★ 刷新是**两遍**（`tests/_gate_evidence.py` 第 3 条）：先带
    `PYFP3D_GATE_FIGURES=1` 跑一遍刷新，再**不带**标志跑一遍验证 ——
    同一次运行里回归锁在本腿写盘之前就跑了，必然红一次。
    """
    import csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        # ★★ 线程数进证据：收敛旗标是"算例 x 线程数"的联合性质，
        #    不带线程数的旗标是无意义的（项目里 gated count 那条同理）。
        nthr = os.environ.get("NUMBA_NUM_THREADS", "unset")
        w.writerow(["case", "level", "kind", "n_threads", "mach", "alpha_deg",
                    "cp_rms_upper", "cp_rms_lower",
                    "cl", "cd", "x_shock", "residual", "converged",
                    "n_limited", "n_floored", "sigma_frozen", "sigma_churn",
                    "cl_ref_cfl3d", "d_cl_rel"])
        for nm, cid, m, a, kind in CASES:
            for lv in LEVELS:
                d = runs[(nm, lv)]
                cref = ref[cid]["cl"]
                dcl = ("" if abs(cref) < 1e-12 else
                       f"{(d['cl'] - cref) / abs(cref):.6e}")
                w.writerow([nm, lv, kind, nthr, m, a,
                            fmt(d["cp_rms_upper"]), fmt(d["cp_rms_lower"]),
                            fmt(d["cl"]), fmt(d["cd"]), fmt(d["x_shock"]),
                            fmt(d["residual"]), int(d["converged"]),
                            d["n_limited"], d["n_floored"],
                            int(d["sigma_frozen"]), int(d["sigma_churn"]),
                            f"{cref:.6f}", dcl])

    # ★★★ Cp 分布是主图。柱状图只留一格作附注 —— 一个 cl 是一个数，
    #     两条曲线差 2 % 可能来自激波位置、前缘峰或尾缘卸载中的任何一处，
    #     甚至来自互相抵消的两处；**Cp 显示差在哪里**。
    fig, axes = plt.subplots(2, 4, figsize=(19, 8.4))
    for k, (nm, cid, m, a, kind) in enumerate(CASES):
        ax = axes.ravel()[k]
        rc = read_2d_cp(REF_DIR, cid, "L4", "none")
        for side, ls in (("upper", "-"), ("lower", "--")):
            if side in rc:
                ax.plot(rc[side][0], rc[side][1], ls, color="#c0392b", lw=1.6,
                        label="CFL3D Euler L4" if side == "upper" else None)
        for lv, sty in (("coarse", dict(lw=0.9, ls=":", color="#9aa0a6")),
                        ("medium", dict(lw=1.4, color="#1a56db"))):
            cur = runs[(nm, lv)]["curve"]
            for side in ("upper", "lower"):
                ax.plot(cur[f"x_{side}"], cur[f"cp_{side}"],
                        label=f"pyFP3D {lv}" if side == "upper" else None,
                        **sty)
        rms = {}
        for side in ("upper", "lower"):
            if side in rc:
                cur = runs[(nm, "medium")]["curve"]
                rms[side] = cp_rms(rc[side][0], rc[side][1],
                                   cur[f"x_{side}"], cur[f"cp_{side}"])[0]
        lab = "GATED" if kind == "gate" else "RECORDED"
        ax.set_title(f'{nm}  [{lab}]\nCp RMS upper {rms.get("upper", float("nan")):.4f}'
                     f'  lower {rms.get("lower", float("nan")):.4f}',
                     fontsize=9, color=("black" if kind == "gate" else "#c0392b"))
        ax.invert_yaxis(), ax.grid(alpha=.3), ax.set_xlabel("x/c")
        if k == 0:
            ax.set_ylabel("$C_p$"), ax.legend(fontsize=7)
    ax = axes.ravel()[6]
    xs = [nm for nm, *_ in CASES]
    ax.barh(xs, [abs(runs[(nm, "medium")]["cl"] - ref[cid]["cl"])
                 for nm, cid, *_ in CASES], color="#5b8def")
    ax.set_xlabel("|d cl| (absolute)"), ax.grid(alpha=.3, axis="x")
    ax.set_title("footnote: absolute cl gap\n(relative is undefined at alpha 0)",
                 fontsize=8)
    axes.ravel()[7].axis("off")
    axes.ravel()[7].text(0.02, 0.98,
        "D05  pyFP3D inviscid vs CFL3D Euler\n\n"
        "Cp is the primary comparison: a cl is one\n"
        "number, and two curves 2 % apart may differ\n"
        "at the shock, the LE peak or the TE unloading\n"
        "-- or at two places that cancel.\n\n"
        "GATED (5 of 7 legs): subcritical cl +0.04 %,\n"
        "alpha = 0 on an ABSOLUTE band, and the M1\n"
        "target M0.80/alpha1.25 at -1.98 %.\n\n"
        "RECORDED: M0.778 (+43.9 %) and M0.803, whose\n"
        "'-47 %' is a near-zero denominator -- its\n"
        "absolute gap is 0.0135.\n\n"
        "*** The convergence flag is NOT gated: it\n"
        "flips in BOTH directions with thread count\n"
        "while cl moves by 0.15-0.36 %.  n_threads is\n"
        "recorded in summary.csv.",
        va="top", ha="left", fontsize=8.2, family="monospace")
    fig.suptitle("D05  pyFP3D inviscid vs CFL3D Euler, NACA0012 -- "
                 "Cp distributions at every condition", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "d05_vs_cfl3d_euler.png"),
                dpi=130)
    plt.close(fig)
