r"""D06 — pyFP3D 无粘 vs CFL3D Euler — RAE2822（2D-5..2D-6）。

Case 7 M0.725/α2.55 · Case 9/10 M0.730/α3.19（`euler_rae2822/`，2 工况 × 3 档）。
α 一律**实验值不修正**（裁决二）。本门读参考的**最细档 L4**
（2026-09-05 参考补第四档；★ 力系数与 Cp 必须读**同一档**）。

---

## ★★★ 实测（先测后登记），以及第一次测量为什么是错的

第一次在 **6 线程**上测，得到"M0.730 medium 不收敛"，我据此写了一条
`assert not converged`。**那条断言自己红了**（测试跑在 8 线程），才把问题揪出来。
项目纪律 #1 是**上限 16，含 BLAS/OMP**。@16 线程重测：

★ 2026-10-10：证据按真 BLAS=1 重刷（登记见
docs/dev_phase_six/20261010-1300-evidence-relock-blas1-prereg.md §3）。M0.725 medium 是
**近简并敏感腿**：跨 BLAS 档落在邻近吸引子上，cl 1.016752 → 1.016375、dx +0.0319 →
+0.0318、Δcl +3.58 % → +3.55 %。
★ 同日：本表 **CFL3D 列与全部派生量已改读 L4**（M0.730 行整行同步：cl 1.113852、
Δcl +0.27 % → +0.06 %、dx +0.0051 → +0.0052；第一行参考格 0.982077 → 0.981571）——
原 L3 档口径的遗留自此核销，判据带不受影响。

| 工况（medium） | pyFP3D cl | CFL3D cl（L4） | Δcl | dx | @6t | @8t | **@16t** |
|---|---|---|---|---|---|---|---|
| M0.725/α2.55 | 1.016375 | 0.981571 | **+3.55 %** | +0.0318 | ✓ | ✓ | ✓ |
| M0.730/α3.19 | 1.113852 | 1.113160 | **+0.06 %** | +0.0052 | ✗ | ✓ | **✓** |

M0.730 medium 逐线程：4 线程 ✗（|R| 5.5e-06，80 步封顶）· 6 线程 ✗（2.9e-06）·
8 线程 ✓（**2.96e-13**）· 12 线程 ✓（3.1e-13），而 **cl 在四者间只差 0.15 %**。

★★★ **收敛旗标是算例与线程数的联合性质，不是算例的性质**，而**答案是稳定的**。
⇒ 本门**不对 `converged` 设任何断言**，只把它连同线程数写进证据；
**设门依据是分歧的大小**。

★ **误差跨零翻号**：两例都是 coarse 负、medium 正（−6.79 → +3.55、
−33.2 → +0.06）。跨零的"一致"更像**穿越**而非收缩 —— 与 D07 上 cl 三档
ratio 0.153、第四档翻成 3.161 是同一个签名。所以只在**幅度**上设带，
并明写它不单调收缩。

★ **方向不设门，因为它在两个翼型上相反**：NACA0012 上 pyFP3D 的激波在 CFL3D
**上游**（D05，4/5），RAE2822 上在**下游**（+0.032 / +0.005）。

## ★ 与 D05 相反的一点：这份参考紧到能分辨

RAE2822 参考自身的激波位置**散布**（L1–L4 极差）是 **0.0032 / 0.0024**，pyFP3D 与它差
0.0318 / 0.0052 ⇒ **10.0× / 2.2×，可分辨**。D05 的 M0.75 是 0.21×，不可分辨。
**同一条判据，两个相反的结论 —— 参照物的分辨率是逐例的性质。**
（口径：0.00088 / 0.00100 与 "35×" 是旧的 L2→L3 最后差，2026-09-05 已更正为**散布** ——
见 `test_shock_difference_is_resolvable_against_this_reference` 的更正记录；数值为 2026-10-10 现行 CSV。）

★ 2026-10-10（W2.1a）：证据加 **xcoarse** 第三档 —— 只进证据 / 锁 / 图 /
`grid_convergence_ours.csv`（RECORDED，用参照侧 `implied_order` 仪器），
**不进任何判据**（本门判据全部点名层级，不经过 `LEVELS`）。

★ 2026-10-11（W2.1b）：|S−D| 判据改写为 **min(标定带, u_val)**（u_val = √(u_num²
+ u_D²)，u_num 来自三级梯子、u_D 来自参照散布，带永不变宽）。逐腿体制：
M0.725 **floor**（u_val ≈ 0.10 > 带 0.06，带不变）、M0.730 **no-bar**
（梯子 coarse 腿 clamping ⇒ u_num 无定义，带 = 标定带 0.03）。
逐腿读数见 `grid_convergence_ours.csv` 的 regime 列。
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
from tests.D._grid_convergence import leg_uncertainty, mesh_h
from tests._gate_evidence import assert_matches_committed, blas_thread_count, fmt
from tests.conftest import REPO_ROOT, gate_figures_enabled

REF_DIR = os.path.join(str(REPO_ROOT), "cases", "reference_data", "cfl3d",
                       "euler_rae2822")
RECIPE = dict(upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=400,
              n_picard_seed=5)
CASES = (
    ("M0.725_a2.55", "rae2822_m0725_a2.55", 0.725, 2.55, "gate"),
    ("M0.730_a3.19", "rae2822_m0730_a3.19", 0.730, 3.19, "gate"),
)
#: 证据层级（runs / summary.csv / 证据锁 / 图 / grid_convergence_ours.csv）。
#: xcoarse 只进证据：本门判据全部用 "coarse" / "medium" 字面量点名层级、
#: 不经过 LEVELS，所以扩档碰不到判据域 —— W2.1a 预注册 C2。
#: （D05 的同结构判据域由 GATE_LEVELS 显式隔开。）
LEVELS = ("xcoarse", "coarse", "medium")

# —— 判据（标定，不是实测值）——
#: ★ 2026-10-11（W2.1b）：|S−D| 腿的带改写为 **min(标定带, u_val)**，
#:   u_val = √(u_num² + u_D²)（u_num = 我方梯子最细差，u_D = 参照散布）——
#:   **带永不变宽**；逐腿体制写在注释尾巴与 `grid_convergence_ours.csv` 的 regime 列。
CL_BAND = 0.06              # M0.725 medium |Δcl|：实测 0.0355；W2.1b：u_val 0.10 > 带 ⇒ floor（带不变）
CL_BAND_M0730 = 0.03        # M0.730 medium |Δcl|：实测 0.0006；W2.1b：梯子 coarse 腿是 clamping ⇒ u_num 无定义 ⇒ no-bar（带 = 标定带）
RESOLVABLE_FRAC = 1.0       # 与参考自身散布比较的可分辨下限（W2.1b 不改写：它是前提锁，不是 |S−D| 带）


def _read_reference(path=None):
    """读已提交的 CFL3D 参考（**L4**，最细档）。★ 行为验证见 `TestReferenceIsLoadBearing`。"""
    import csv
    p = path or os.path.join(REF_DIR, "forces.csv")
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            #: ★ 读**最细**档（2026-09-05 起参考有 L4）
            if r["level"] == "L4":
                out[r["case"]] = dict(cl=float(r["cl"]), cd=float(r["cd"]))
    if not out:
        raise RuntimeError(f"{p}: no L4 rows -- the reference layout changed")
    return out


def _reference_uncertainty():
    """参考自身的**散布**（各档极差）。

    ★★ 不是单一的 `delta_L2_L3`：对一个判定不稳的量，最后一步的差可能恰好
    很小而该量仍在各档间游走 —— **游走的范围才是不确定度**，且无论收敛与否
    都有定义。
    """
    import csv
    out = {}
    with open(os.path.join(REF_DIR, "grid_convergence.csv")) as fh:
        for r in csv.DictReader(fh):
            vals = [float(r[lv]) for lv in ("L1", "L2", "L3", "L4")
                    if r.get(lv)]
            out[(r["case"], r["quantity"])] = (max(vals) - min(vals)
                                               if len(vals) >= 2
                                               else float("nan"))
    return out


def _reference_shock(case):
    import csv
    with open(os.path.join(REF_DIR, "shock.csv")) as fh:
        for r in csv.DictReader(fh):
            if (r["case"] == case and r["level"] == "L4"
                    and r["surface"] == "upper" and r.get("x_shock")):
                return float(r["x_shock"])
    return float("nan")


MESH_DIR = os.path.join(str(REPO_ROOT), "cases", "meshes", "rae2822_2.5d")


def _one(level, m_inf, alpha):
    mc, wc = cut_wake(read_mesh(os.path.join(
        str(REPO_ROOT), "cases", "meshes", "rae2822_2.5d", f"{level}.msh")))
    dz = float(np.ptp(mc.nodes[:, 2]))
    r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=alpha, **RECIPE)
    phi = np.asarray(r["phi"])
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], phi,
                                alpha_deg=alpha, u_inf=1.0, s_ref=dz,
                                m_inf=m_inf)
    #: ★★ RAE2822 是**有弯度**的：上下面必须用局部外法向判，不能用 y > 0。
    #:   那正是 D12 修掉的库缺陷（37 个点被错分）。`section_cp_curve` 已修，
    #:   这里记下来是因为任何新写的截面处理都会再踩一次。
    cur = section_cp_curve(mc, phi, z=float(np.mean(mc.nodes[:, 2])),
                           smooth_passes=1, m_inf=m_inf)
    up = shock_report(cur, m_inf)["upper"]
    frz = r.get("sigma_freeze_report") or {}
    return dict(
        cl=float(f["cl"]), cd=float(f["cd_pressure"]),
        x_shock=(float(up["x_shock"]) if up.get("has_shock") else float("nan")),
        converged=bool(r.get("converged")),
        residual=float(np.asarray(r["residual_history"], float)[-1]),
        n_limited=int(r.get("n_limited", 0)),
        n_floored=int(r.get("n_floored", 0)),
        sigma_frozen=bool(frz.get("frozen_in_transient", False)),
        #: ★ Cp 曲线带出来，因为图必须是 Cp 分布而不是力系数柱状图
        #:   （使用者裁决 2026-09-05）：一个 cl 是一个数，Cp 显示差在哪里。
        curve=cur,
    )


def _with_cp_rms(d, cid):
    """★★ Cp RMS 是主要比较量 ⇒ 必须进证据 CSV，不能只活在 PNG 里。"""
    rc = read_2d_cp(REF_DIR, cid, "L4", "none")
    for side in ("upper", "lower"):
        d[f"cp_rms_{side}"] = cp_rms(
            rc[side][0], rc[side][1],
            d["curve"][f"x_{side}"], d["curve"][f"cp_{side}"])[0]
    return d


@pytest.fixture(scope="module")
def runs():
    return {(nm, lv): _with_cp_rms(_one(lv, m, a), cid)
            for nm, cid, m, a, _k in CASES for lv in LEVELS}


@pytest.fixture(scope="module")
def ref():
    return _read_reference()


class TestWhatIsGateable:
    def test_m0725_lift_magnitude(self, runs, ref):
        """幅度带，**不是**收缩判据：误差跨零（coarse −6.79 % → medium +3.55 %），
        所以不能要求单调收缩，只能要求 medium 的幅度在带内。

        ★ 2026-10-11（W2.1b）：带 = min(0.06, u_val)；该腿梯子差大
        （u_val ≈ 0.10 相对）⇒ **floor**，带不变。"""
        got = runs[("M0.725_a2.55", "medium")]["cl"]
        want = ref["rae2822_m0725_a2.55"]["cl"]
        rel = abs(got - want) / abs(want)
        u = leg_uncertainty([runs[("M0.725_a2.55", lv)]["cl"] / abs(want)
                             for lv in LEVELS], mesh_h(MESH_DIR),
                            _reference_uncertainty()[("rae2822_m0725_a2.55", "cl")]
                            / abs(want), CL_BAND, sd=rel)
        assert rel <= u["band"], (
            f"M0.725 cl {got:.6f} vs CFL3D Euler {want:.6f} = {rel*100:.2f} % "
            f"> band {u['band']*100:.2f} % = min({CL_BAND*100:.0f} %, u_val);"
            f"regime {u['regime']}, u_num {u['u_num']}, u_D {u['u_D']}")

    def test_m0730_lift_magnitude(self, runs, ref):
        """@16 线程 Δcl = +0.06 %（L4 口径）。★ 同样**不断言 `converged`**。

        ★ 2026-10-11（W2.1b）：梯子 coarse 腿是 clamping ⇒ u_num 无定义 ⇒
        **no-bar**，带 = 标定带 0.03（本条现在是没有误差棒的锁，写明）。"""
        got = runs[("M0.730_a3.19", "medium")]["cl"]
        want = ref["rae2822_m0730_a3.19"]["cl"]
        rel = abs(got - want) / abs(want)
        u = leg_uncertainty([runs[("M0.730_a3.19", lv)]["cl"] / abs(want)
                             for lv in LEVELS], mesh_h(MESH_DIR),
                            _reference_uncertainty()[("rae2822_m0730_a3.19", "cl")]
                            / abs(want), CL_BAND_M0730, sd=rel)
        assert rel <= u["band"], (
            f"M0.730 cl {got:.6f} vs CFL3D Euler {want:.6f} = {rel*100:.2f} % "
            f"> band {u['band']*100:.2f} % = min({CL_BAND_M0730*100:.0f} %, u_val);"
            f"regime {u['regime']}, u_num {u['u_num']}, u_D {u['u_D']}")

    def test_shock_difference_is_resolvable_against_this_reference(self, runs):
        """★★ 与 D05 的 M0.75 相反：这里的差远大于参考自身的散布，**可分辨**。

        ★ **数字更正 2026-09-05**：原来按单一的 `delta_L2_L3` = 0.00088 说
        "35×"。改用**散布**（四档 0.63302 → 0.62986，极差 **0.003169**）后是
        **9.8×**（2026-10-10 现行口径 10.0×）—— **结论不变，倍数缩水四倍**。单一的最后差把参照物说得
        比它实际更精确。

        锁住的是"可分辨"这个前提本身。若参考日后放松到不可分辨，本条会红，
        那时上面那条激波比较就该改成 RECORDED —— 与 D05 的镜像断言配对。"""
        unc = _reference_uncertainty()[("rae2822_m0725_a2.55", "x_shock_upper")]
        diff = abs(runs[("M0.725_a2.55", "medium")]["x_shock"]
                   - _reference_shock("rae2822_m0725_a2.55"))
        assert diff > RESOLVABLE_FRAC * unc, (
            f"the shock difference {diff:.6f} has fallen to or below the "
            f"reference's own uncertainty {unc:.6f} -- it is no longer "
            f"resolvable and must be re-specified as RECORDED")


class TestWhatIsOnlyRecorded:
    def test_convergence_flag_is_recorded_not_gated(self, runs):
        """★★★ 这里曾经有一条 `assert not converged`。

        它把一个**环境依赖的结果、而且是坏的那一侧**焊进了门里：M0.730 medium
        在 4/6 线程不收敛、8/12/16 线程收敛，而 cl 在四者间只差 0.15 %。
        **是那条断言自己红了**才暴露问题。⇒ 现在只断言旗标**可读**，
        它的值连同线程数写进证据。"""
        for nm, _c, _m, _a, _k in CASES:
            r = runs[(nm, "medium")]
            assert isinstance(r["converged"], bool)
            assert np.isfinite(r["residual"])

    def test_direction_of_the_shock_offset_is_not_gated(self, runs):
        """★ 记录方向，**不设门**：NACA0012 上 pyFP3D 的激波在上游，
        RAE2822 上在下游。没有一致方向 ⇒ 没有可登记的符号。
        本条只断言两个符号都被算出来了、可读。"""
        for nm, _c, _m, _a, _k in CASES:
            d = runs[(nm, "medium")]["x_shock"]
            assert np.isfinite(d), f"{nm}: no shock detected at medium"


class TestReferenceIsLoadBearing:
    def test_gate_actually_reads_the_reference(self, tmp_path):
        """扰动副本必须传播 —— F06 的教训，行为验证。"""
        import csv
        src = os.path.join(REF_DIR, "forces.csv")
        with open(src) as fh:
            rows = list(csv.DictReader(fh))
            cols = list(rows[0].keys())
        for r in rows:
            if r["case"] == "rae2822_m0725_a2.55" and r["level"] == "L4":
                r["cl"] = f"{float(r['cl']) + 0.4321:.6f}"
        dst = tmp_path / "forces.csv"
        with open(dst, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        base = _read_reference()["rae2822_m0725_a2.55"]["cl"]
        pert = _read_reference(str(dst))["rae2822_m0725_a2.55"]["cl"]
        assert abs((pert - base) - 0.4321) < 1e-9, (
            "the reader did not follow a perturbed reference")


class TestCommittedEvidenceIsLoadBearing:
    MEASURED = ("cp_rms_upper", "cp_rms_lower", "cl", "cd", "x_shock",
                "blas_threads")

    def test_matches_committed_summary(self, runs, gate_evidence_dir):
        #: ★ blas_threads 不在 runs 里（它是库实际在用的 BLAS 线程数，不是算出来
        #:   的量 —— H30 实测答案只随它动），从库读、跟着 MEASURED 进锁。
        fresh = {f"{nm}|{lv}": {**{k: fmt(runs[(nm, lv)][k]) for k in self.MEASURED
                                   if k != "blas_threads"},
                                "blas_threads": fmt(blas_thread_count())}
                 for nm, _c, _m, _a, _k in CASES for lv in LEVELS}
        assert_matches_committed(
            gate_evidence_dir, fresh, self.MEASURED,
            key_of=lambda r: f"{r['case']}|{r['level']}",
            refresh_hint="PYFP3D_GATE_FIGURES=1 pytest "
                         "tests/D/test_D06_euler_rae2822.py")


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图/CSV 证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_evidence(runs, ref, gate_evidence_dir):
    """★ 刷新是两遍：先带标志刷新，再不带标志验证。
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
    import csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    unc = _reference_uncertainty()
    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        # ★★ 线程数进证据：收敛旗标是"算例 x 线程数"的联合性质，
        #    不带线程数的旗标是无意义的（项目里 gated count 那条同理）。
        # ★ n_threads 是 numba 的旋钮；blas_threads 是库实际在用的 BLAS
        #    线程数（H30 实测答案只随它动）—— 两个都记。
        nthr = os.environ.get("NUMBA_NUM_THREADS", "unset")
        w.writerow(["case", "level", "kind", "n_threads", "blas_threads",
                    "mach", "alpha_deg", "cp_rms_upper", "cp_rms_lower", "cl",
                    "cd", "x_shock", "residual", "converged", "n_limited",
                    "n_floored", "sigma_frozen", "cl_ref_cfl3d", "d_cl_rel",
                    "x_shock_ref", "d_x_shock", "ref_shock_uncertainty",
                    "dx_over_ref_uncertainty"])
        for nm, cid, m, a, kind in CASES:
            xr = _reference_shock(cid)
            u = unc.get((cid, "x_shock_upper"), float("nan"))
            for lv in LEVELS:
                d = runs[(nm, lv)]
                dx = d["x_shock"] - xr
                w.writerow([nm, lv, kind, nthr, blas_thread_count(), m, a,
                            fmt(d["cp_rms_upper"]), fmt(d["cp_rms_lower"]),
                            fmt(d["cl"]), fmt(d["cd"]),
                            fmt(d["x_shock"]), fmt(d["residual"]),
                            int(d["converged"]), d["n_limited"],
                            d["n_floored"], int(d["sigma_frozen"]),
                            f"{ref[cid]['cl']:.6f}",
                            f"{(d['cl'] - ref[cid]['cl']) / abs(ref[cid]['cl']):.6e}",
                            f"{xr:.6f}", f"{dx:+.6f}", f"{u:.6f}",
                            f"{abs(dx) / u:.2f}"])

    #: ★ W2.1a/W2.1b：我方三级网格收敛读数（RECORDED，不进值锁）—— 估计器与
    #:   判定词汇都是参照侧的（tests/D/_grid_convergence.py）；u_D = 参照散布，
    #:   sd / band = 该腿判据的 |S−D| 与标定带（原生单位）；cp_rms_* 参照侧
    #:   没有同名阶梯量 ⇒ u_D 缺 ⇒ no-bar。
    from tests.D._grid_convergence import write_grid_convergence_ours
    cids = {nm: cid for nm, cid, *_ in CASES}

    def _u_d_of(case, q):
        key = {"cl": "cl", "cd": "cd", "x_shock": "x_shock_upper"}.get(q)
        return None if key is None else unc.get((cids[case], key))

    def _sd_of(case, q):
        d, cid = runs[(case, "medium")], cids[case]
        if q in ("cl", "cd"):
            return abs(d[q] - ref[cid][q])
        if q == "x_shock":
            xr = _reference_shock(cid)
            return (None if not (np.isfinite(d["x_shock"]) and np.isfinite(xr))
                    else abs(d["x_shock"] - xr))
        return d[q]  # cp_rms_* 本身就是 |S−D|

    def _band_of(case, q):
        if q != "cl":
            return None  # cd / x_shock / cp_rms 在本门没有标定带
        return ((CL_BAND if case == "M0.725_a2.55" else CL_BAND_M0730)
                * abs(ref[cids[case]]["cl"]))

    write_grid_convergence_ours(
        os.path.join(str(gate_evidence_dir), "grid_convergence_ours.csv"),
        runs, [nm for nm, *_ in CASES], LEVELS,
        MESH_DIR, u_d_of=_u_d_of, sd_of=_sd_of, band_of=_band_of)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for k, (nm, cid, m, a, kind) in enumerate(CASES):
        ax = axes[k]
        rc = read_2d_cp(REF_DIR, cid, "L4", "none")
        for side, ls in (("upper", "-"), ("lower", "--")):
            ax.plot(rc[side][0], rc[side][1], ls, color="#c0392b", lw=1.6,
                    label="CFL3D Euler L4" if side == "upper" else None)
        for lv, sty in (("xcoarse", dict(lw=0.7, ls="--", color="#c7ccd4")),
                        ("coarse", dict(lw=0.9, ls=":", color="#9aa0a6")),
                        ("medium", dict(lw=1.4, color="#1a56db"))):
            cur = runs[(nm, lv)]["curve"]
            for side in ("upper", "lower"):
                ax.plot(cur[f"x_{side}"], cur[f"cp_{side}"],
                        label=f"pyFP3D {lv}" if side == "upper" else None,
                        **sty)
        cur = runs[(nm, "medium")]["curve"]
        ru = cp_rms(rc["upper"][0], rc["upper"][1],
                    cur["x_upper"], cur["cp_upper"])[0]
        rl = cp_rms(rc["lower"][0], rc["lower"][1],
                    cur["x_lower"], cur["cp_lower"])[0]
        ax.axvline(_reference_shock(cid), color="#c0392b", ls="-.", lw=0.9)
        ax.axvline(runs[(nm, "medium")]["x_shock"], color="#1a56db",
                   ls="-.", lw=0.9)
        ax.set_title(f'{nm}  [GATED]\nCp RMS upper {ru:.4f}  lower {rl:.4f}',
                     fontsize=9)
        ax.invert_yaxis(), ax.grid(alpha=.3), ax.set_xlabel("x/c")
        if k == 0:
            ax.set_ylabel("$C_p$"), ax.legend(fontsize=7)
    axes[2].axis("off")
    axes[2].text(0.02, 0.98,
        "D06  pyFP3D inviscid vs CFL3D Euler, RAE2822\n\n"
        "Cp is the primary comparison; the dash-dot\n"
        "lines are the detected shocks (red CFL3D,\n"
        "blue pyFP3D).\n\n"
        "BOTH legs GATED: dcl +3.55 % and +0.06 %.\n\n"
        "*** The shock difference IS resolvable here\n"
        "(10x the reference's own L1-L4 spread of\n"
        "0.0032), where on NACA0012 M0.75 the same\n"
        "comparison is 0.21x and UNRESOLVABLE.\n"
        "A reference's resolution is a PER-CASE\n"
        "property, not a dataset property.\n\n"
        "*** The shock sits DOWNSTREAM of CFL3D here\n"
        "and UPSTREAM on NACA0012 -- three geometries,\n"
        "two directions, so no registered sign.",
        va="top", ha="left", fontsize=8.2, family="monospace")
    fig.suptitle("D06  pyFP3D inviscid vs CFL3D Euler, RAE2822 -- "
                 "Cp distributions", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "d06_vs_cfl3d_euler.png"),
                dpi=130)
    plt.close(fig)
