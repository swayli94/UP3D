r"""D07 — pyFP3D 无粘 vs CFL3D Euler — ONERA M6（3D-1 + 3D-2）。

工况：**M 0.8395 / α 3.06**，AGARD AR-138 TEST 2308 逐字（α 实验值不修正）。
参考 `euler_onera_m6/` 四档 + REF，本门读 **L4**（最细，11.278 M 点）。

**3D-2 = M 0.50 / α 3.06** 于 2026-09-05 补齐，参考在
`euler_onera_m6_m050/`（四档，**4/4 力系数量带误差棒**）。数据请求点名它的理由是
"项目的 M6 / 翼身**亚声速**锚点（B9、B32 的 M0.5）**目前完全没有外部参照**"。

★★★ **我对这条的预期是错的，而项目自己的笔记里就写着为什么。** 我预期
"亚临界是全速势最占便宜的地方"（D05 二维 M0.50 实测 +0.04 %），因此这条会给出
一道**紧**门。实测 −1.54 %（coarse）→ **−6.89 %**（medium），**误差随加密变大**，
与跨声速那条同型。原因是 **M∞ = 0.50 时 M_max 已到 0.9964** —— 这个工况
**根本不是亚临界**。CLAUDE.md 的操作笔记里有这条：
*"M6 at α 3.06 is NOT subcritical at M0.70 — measured M_max 1.5358"*。
**该查它，而不是从"M0.50 听起来是亚声速"推断。**

★★ **但它给出了比一道紧门更有价值的东西。** 两个工况的 Cp 展向结构几乎相同：

| | M 0.50 | M 0.8395 |
|---|---|---|
| 上表面 Cp RMS（y/b 0.20 → 0.99）| 0.0799 → 0.1491 | 0.127 → 0.244 |
| **外/内比** | **1.87** | **1.92** |
| M_max（medium）| 0.996 | 1.995 |

⇒ **越往翼尖越重的吸力亏损不是激波现象** —— 它在没有激波的 M0.50 上形状几乎
一致，只是幅度约 60 %。病因被缩小到"**沿展向作用、外翼更重、与马赫数无关**"。
★ 一个已知候选：生产配置的 `tip_taper`（`vanish_smooth`, 0.05·b_semi）在两条上
都开着，CLAUDE.md 记它代价 **−1.3 % cl_p** 且作用在翼尖 —— 但那只占 −6.89 %
的五分之一，**至多是部分原因，不归因**。

★ **3D-2 走的是另一条求解路径**，这是它的意义而不是麻烦：`_m6_case` 驱动
`solve_newton_transonic`，后者直接拒绝 M 0.50（"upward ramp only: m_inf 0.5 <
m_start 0.7"）—— Mach 阶梯是为跨声速存在的。这条工况要锚定的那些亚声速锚点
（B9、B32 的 M0.5）走的正是 `solve_newton_lifting`，所以本门也走它。

---

## ★★★ 实测：激波位置随加密**收缩**，而 cl 反而**发散**

| 量 | coarse | medium | CFL3D L4 |
|---|---|---|---|
| cl | 0.277536（−6.06 %）| 0.268366（**−9.17 %**）| 0.295448 |
| 激波 η 0.44 | +0.0344 | **+0.0036** | 0.5912 |
| 激波 η 0.65 | +0.0615 | **+0.0279** | 0.5115 |
| 激波 η 0.90 | +0.1300 | **+0.0718** | 0.3157 |
| M_max | 1.6109 | 1.9952 | — |
| 收敛 | ✓ \|R\| 3.3e-11 | ✓ \|R\| 7.2e-15 | — |

★★ **两个量朝相反方向走**：三个站位的激波位置都单调收缩（9.5× / 2.2× / 1.8×），
而 cl 的误差从 6.06 % 涨到 9.17 %。⇒ **不能拿"整体一致"下结论** —— 设门的是
收缩的那一半，发散的那一半只能记录。

★ 激波全部在 CFL3D **下游**（+0.004…+0.072）。这与 RAE2822 同向、与 NACA0012
**相反**（D05 上 4/5 在上游）。三个几何、两个方向 ⇒ **没有可登记的符号**，
这已是第三次确认。

## ★★ 两个口径，都是查出来的不是假定的

1. **马赫数**：参考是 **M 0.8395**（TEST 2308 逐字），而项目已有的 M6 锚点
   （G8.2 / P14）是 **M 0.84**。本门在 **0.8395** 上重跑 —— 复用 0.84 的锚点
   会是跨工况比较。
2. **参考面积**：pyFP3D 的 M6 网格是**米制**、CFL3D 网格归一到**根弦 = 1**，
   两边各自除以自身单位下的面积 ⇒ **cl 可比** ✓。但 pyFP3D 用的是**离散**平面
   面积 0.760177，参考用**解析**值 0.752951，差 **+0.96 %**，而且
   **不随加密收缩**（coarse +0.960 % / medium +0.965 %）⇒ 是**定义差异**
   （应为圆翼尖端帽的投影面积），不是离散误差。**这是一条 0.96 % 的系统偏置，
   方向使 pyFP3D 的 cl 偏低**，即它解释了 −9.17 % 里的约十分之一。
"""
import os

import numpy as np
import pytest

from tests.conftest import REPO_ROOT, gate_figures_enabled
from tests.D._cfl3d_cp import cp_rms, read_3d_cp
from tests._gate_evidence import assert_matches_committed, fmt

REF_DIR = os.path.join(str(REPO_ROOT), "cases", "reference_data", "cfl3d",
                       "euler_onera_m6")
M_INF, ALPHA = 0.8395, 3.06
LEVELS = ("coarse", "medium")
ETAS = (0.44, 0.65, 0.90)
#: ★ 参考数据集的**七个实验站位** —— Cp 图画全部七个，
#:   而设门只用 ETAS 那三个（其余站位的参考激波前提未逐一核过）。
CP_ETAS = (0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99)

# —— 判据（标定，不是实测值）——
SHOCK_CONTRACT_MIN = 1.5     # 每个站位 coarse->medium 的收缩：实测 9.5 / 2.2 / 1.8
SHOCK_MAX_MEDIUM = 0.10      # medium 的最大激波偏差：实测 0.0036 / 0.0279 / 0.0718
CL_RECORD_MAX = 0.15         # cl 只记录，带子宽到只拦住量级性变化：实测 0.0917


def _read_reference(path=None):
    import csv
    p = path or os.path.join(REF_DIR, "forces.csv")
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            out[r["level"]] = dict(cl=float(r["cl"]), cd=float(r["cd"]))
    if "L4" not in out:
        raise RuntimeError(f"{p}: no L4 row -- the reference layout changed")
    return out


def _reference_shocks():
    import csv
    out = {}
    with open(os.path.join(REF_DIR, "shock.csv")) as fh:
        for r in csv.DictReader(fh):
            if (r["level"] == "L4" and r["surface"] == "upper"
                    and r["x_shock"]):
                out[float(r["eta_requested"])] = dict(
                    x=float(r["x_shock"]), premise=r["detector_premise"])
    return out


def _one(level):
    #: ★ 复用 tests/E 自己的 M6 案例函数，不另起一条流水线。
    #:   它返回 FOUR 个值 —— 返回元数是要核对的，不是回忆的。
    from tests.E.test_E01_p8_newton_anchors import _m6_case
    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.mesh.wake_cut import cut_wake
    from pyfp3d.post.surface import planform_area
    from pyfp3d.meshgen.wing3d import B_SEMI
    from pyfp3d.post.section_cut import section_cp_curve
    r, forces, shocks, wall = _m6_case(level, m_inf=M_INF, alpha=ALPHA)
    mc, _wc = cut_wake(read_mesh(os.path.join(
        str(REPO_ROOT), "cases", "meshes", "onera_m6", f"{level}.msh")))
    #: ★★ Cp 分布是主图（使用者裁决 2026-09-05）：三维上尤其如此 ——
    #:   一个 cl 把七个站位、上下表面、激波位置全压成一个数。
    curves = {e: section_cp_curve(mc, r["phi"], eta=e, b_semi=B_SEMI,
                                  m_inf=M_INF)
              for e in CP_ETAS}
    return dict(
        cl=float(forces["cl"]), shocks={k: float(v) for k, v in shocks.items()},
        converged=bool(r.get("converged")),
        residual=float(np.asarray(r["residual_history"], float)[-1]),
        mach_max=float(np.sqrt(r["mach2_max"])),
        s_ref_discrete=float(planform_area(mc.nodes,
                                           mc.boundary_faces["wall"])),
        curves=curves)


@pytest.fixture(scope="module")
def runs():
    """★★ 每站位的 Cp RMS 进证据 —— Cp 是主要比较量，只画在 PNG 里的量
    改一行代码不会有任何东西变红。"""
    out = {}
    for lv in LEVELS:
        d = _one(lv)
        for e in CP_ETAS:
            rc = read_3d_cp(REF_DIR, "L4", e, "none")
            for side in ("upper", "lower"):
                d[f"cp_rms_{side}_eta{e:.2f}"] = cp_rms(
                    rc[side][0], rc[side][1], d["curves"][e][f"x_{side}"],
                    d["curves"][e][f"cp_{side}"])[0]
        out[lv] = d
    return out


@pytest.fixture(scope="module")
def ref():
    return _read_reference()


class TestWhatIsGateable:
    r"""★ 只有**收缩**的那一半。"""

    @pytest.mark.parametrize("eta", ETAS)
    def test_shock_offset_contracts_with_refinement(self, runs, eta):
        """三个站位都必须收缩。实测 9.5× / 2.2× / 1.8×。"""
        xr = _reference_shocks()[eta]["x"]
        d = {lv: abs(runs[lv]["shocks"][eta] - xr) for lv in LEVELS}
        assert d["coarse"] / d["medium"] >= SHOCK_CONTRACT_MIN, (
            f"eta {eta}: shock offset must contract, got coarse {d['coarse']:.4f} "
            f"-> medium {d['medium']:.4f} = "
            f"{d['coarse']/max(d['medium'],1e-12):.2f}x")

    @pytest.mark.parametrize("eta", ETAS)
    def test_shock_offset_magnitude_at_medium(self, runs, eta):
        xr = _reference_shocks()[eta]["x"]
        d = abs(runs["medium"]["shocks"][eta] - xr)
        assert d <= SHOCK_MAX_MEDIUM, (
            f"eta {eta}: medium shock offset {d:.4f} > {SHOCK_MAX_MEDIUM}")

    def test_the_reference_shocks_used_here_pass_their_own_premise(self):
        """★★ 参考侧的 η 0.99 因 Cp\\*-掠过被撤回（见 `euler_onera_m6/`）。
        本门用的三个站位必须是**前提通过**的那些 —— 拿一个被撤回的参考值设门
        就是把参考的缺陷搬进门里。"""
        sh = _reference_shocks()
        for eta in ETAS:
            assert sh[eta]["premise"] == "ok", (
                f"eta {eta}: the CFL3D reference shock now FAILS its own "
                f"detector premise ({sh[eta]['premise']}) -- it can no longer "
                f"anchor a gate")


class TestWhatIsOnlyRecorded:
    def test_lift_error_grows_with_refinement(self, runs, ref):
        """★★★ cl 的误差**变大**：−6.06 % → −9.17 %。**记录，不设门。**

        红了（误差超过 15 %）说明量级变了，需要重新诊断；
        而如果它开始**收缩**，本条不会红 —— 那需要人去把这条腿升级为设门，
        所以下面单独有一条记录收缩比的断言。"""
        want = ref["L4"]["cl"]
        e = {lv: abs(runs[lv]["cl"] - want) / abs(want) for lv in LEVELS}
        assert e["medium"] <= CL_RECORD_MAX, (
            f"cl error at medium {e['medium']*100:.2f} % exceeds the recording "
            f"band {CL_RECORD_MAX*100:.0f} % -- the disagreement changed scale")
        assert e["medium"] > e["coarse"], (
            f"cl error now CONTRACTS ({e['coarse']*100:.2f} % -> "
            f"{e['medium']*100:.2f} %) -- an improvement.  Re-specify this leg "
            f"as a gated comparison.")

    def test_reference_area_caliber_offset_is_recorded(self, runs):
        """★★ pyFP3D 用离散平面面积、参考用解析值，差 +0.96 %，**不随加密收缩**
        ⇒ 定义差异（圆翼尖端帽），不是离散误差。它使 pyFP3D 的 cl 系统偏低，
        约占 −9.17 % 的十分之一。锁住"它仍然是定义差异"这一点。"""
        ana = 0.5 * (0.8059 + 0.4529) * 1.1963
        off = {lv: runs[lv]["s_ref_discrete"] / ana - 1.0 for lv in LEVELS}
        assert all(0.005 < o < 0.02 for o in off.values()), off
        assert abs(off["medium"] - off["coarse"]) < 0.002, (
            f"the reference-area offset now CHANGES with refinement "
            f"({off['coarse']*100:.3f} % -> {off['medium']*100:.3f} %) -- it "
            f"was a definition difference, now it looks like discretisation")


    def test_the_suction_deficit_grows_outboard(self, runs):
        """★★★ 这是 Cp 图揭示、而 cl 藏住的东西（使用者裁决 2026-09-05：
        Cp 分布才是验证计算的东西）。

        实测 medium 上表面 Cp RMS 沿展向单调增大：
        0.127 / 0.162 / 0.161 / 0.184 / 0.222 / 0.223 / 0.244 ——
        pyFP3D 的前缘吸力峰与超声速平台在每个站位都偏浅，**越往翼尖越重**。
        ⇒ −9.17 % 的 cl 亏损**不是均匀偏置**，参考面积口径那 0.96 % 只是零头。
        **记录，不设门**：它是一个待解释的结构，不是一个要收紧的目标。"""
        r = [runs["medium"][f"cp_rms_upper_eta{e:.2f}"] for e in CP_ETAS]
        assert all(np.isfinite(v) for v in r), r
        assert max(r) < 0.5, f"upper-surface Cp RMS blew up: {r}"
        assert r[-1] > r[0], (
            f"the outboard Cp RMS is no longer the largest ({r[0]:.4f} at "
            f"y/b 0.20 vs {r[-1]:.4f} at 0.99) -- the recorded spanwise "
            f"structure changed and needs re-reading")


class TestReferenceIsLoadBearing:
    def test_gate_actually_reads_the_reference(self, tmp_path):
        import csv
        src = os.path.join(REF_DIR, "forces.csv")
        with open(src) as fh:
            rows = list(csv.DictReader(fh))
            cols = list(rows[0].keys())
        for r in rows:
            if r["level"] == "L4":
                r["cl"] = f"{float(r['cl']) + 0.1357:.6f}"
        dst = tmp_path / "forces.csv"
        with open(dst, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        base = _read_reference()["L4"]["cl"]
        pert = _read_reference(str(dst))["L4"]["cl"]
        assert abs((pert - base) - 0.1357) < 1e-9


class TestCommittedEvidenceIsLoadBearing:
    MEASURED = (tuple(f"cp_rms_upper_eta{e:.2f}" for e in CP_ETAS)
                + ("cl", "mach_max", "s_ref_discrete"))

    def test_matches_committed_summary(self, runs, gate_evidence_dir):
        fresh = {lv: {k: fmt(runs[lv][k]) for k in self.MEASURED}
                 for lv in LEVELS}
        assert_matches_committed(
            gate_evidence_dir, fresh, self.MEASURED,
            key_of=lambda r: r["level"],
            refresh_hint="PYFP3D_GATE_FIGURES=1 pytest "
                         "tests/D/test_D07_euler_onera_m6.py")


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图/CSV 证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_evidence(runs, ref, gate_evidence_dir):
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sh = _reference_shocks()
    nthr = os.environ.get("NUMBA_NUM_THREADS", "unset")
    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["level", "n_threads", "mach", "alpha_deg", "cl",
                    "cl_ref_cfl3d_L4", "d_cl_rel", "residual", "converged",
                    "mach_max", "s_ref_discrete", "s_ref_analytic"]
                   + [f"cp_rms_upper_eta{e:.2f}" for e in CP_ETAS]
                   + [f"cp_rms_lower_eta{e:.2f}" for e in CP_ETAS]
                   + [f"x_shock_eta{e}" for e in ETAS]
                   + [f"d_x_shock_eta{e}" for e in ETAS])
        ana = 0.5 * (0.8059 + 0.4529) * 1.1963
        for lv in LEVELS:
            d = runs[lv]
            w.writerow([lv, nthr, M_INF, ALPHA, fmt(d["cl"]),
                        f'{ref["L4"]["cl"]:.6f}',
                        f'{(d["cl"] - ref["L4"]["cl"]) / ref["L4"]["cl"]:.6e}',
                        fmt(d["residual"]), int(d["converged"]),
                        fmt(d["mach_max"]), fmt(d["s_ref_discrete"]),
                        f"{ana:.6f}"]
                       + [fmt(d[f"cp_rms_upper_eta{e:.2f}"]) for e in CP_ETAS]
                       + [fmt(d[f"cp_rms_lower_eta{e:.2f}"]) for e in CP_ETAS]
                       + [fmt(d["shocks"][e]) for e in ETAS]
                       + [f'{d["shocks"][e] - sh[e]["x"]:+.6f}' for e in ETAS])

    fig, axes = plt.subplots(2, 4, figsize=(19, 8.4))
    for k, e in enumerate(CP_ETAS):
        ax = axes.ravel()[k]
        rc = read_3d_cp(REF_DIR, "L4", e, "none")
        for side, ls in (("upper", "-"), ("lower", "--")):
            ax.plot(rc[side][0], rc[side][1], ls, color="#c0392b", lw=1.5,
                    label="CFL3D Euler L4" if side == "upper" else None)
        for lv, sty in (("coarse", dict(lw=0.9, ls=":", color="#9aa0a6")),
                        ("medium", dict(lw=1.4, color="#1a56db"))):
            cur = runs[lv]["curves"][e]
            for side in ("upper", "lower"):
                ax.plot(cur[f"x_{side}"], cur[f"cp_{side}"],
                        label=f"pyFP3D {lv}" if side == "upper" else None,
                        **sty)
        cur = runs["medium"]["curves"][e]
        ru = cp_rms(rc["upper"][0], rc["upper"][1],
                    cur["x_upper"], cur["cp_upper"])[0]
        gated = e in ETAS
        if gated:
            ax.axvline(sh[e]["x"], color="#c0392b", ls="-.", lw=0.9)
            ax.axvline(runs["medium"]["shocks"][e], color="#1a56db",
                       ls="-.", lw=0.9)
        ax.set_title(f'y/b = {e:.2f}   Cp RMS upper {ru:.4f}'
                     + ("   [shock GATED]" if gated else ""), fontsize=9,
                     color=("black" if gated else "#555555"))
        ax.invert_yaxis(), ax.grid(alpha=.3), ax.set_xlabel("x/c")
        if k == 0:
            ax.set_ylabel("$C_p$"), ax.legend(fontsize=7)
    axes.ravel()[7].axis("off")
    axes.ravel()[7].text(0.02, 0.98,
        "D07  pyFP3D inviscid vs CFL3D Euler\n"
        "ONERA M6, M 0.8395 / alpha 3.06 (TEST 2308)\n\n"
        "Cp at all SEVEN measured stations; the three\n"
        "with dash-dot shock markers are the GATED ones\n"
        "(their reference shock passes the dataset's own\n"
        "Cp*-grazing premise; y/b = 0.99 is WITHDRAWN on\n"
        "the reference side, so it is drawn but not\n"
        "gated).\n\n"
        "*** THE TWO HALVES GO OPPOSITE WAYS.\n"
        "Shock offsets CONTRACT under refinement --\n"
        "9.5x / 2.2x / 1.8x at y/b 0.44 / 0.65 / 0.90 --\n"
        "while cl DIVERGES, -6.06 % to -9.17 %.\n"
        "So the shocks are gated and the lift recorded;\n"
        "neither half may speak for the whole.\n\n"
        "*** Reference-area caliber: pyFP3D uses the\n"
        "DISCRETE planform 0.760177, the reference the\n"
        "ANALYTIC 0.752951 -- +0.96 %, and it does NOT\n"
        "shrink with refinement, so it is a definition\n"
        "difference (the rounded tip cap) worth about a\n"
        "tenth of the -9.17 %.",
        va="top", ha="left", fontsize=8.2, family="monospace")
    fig.suptitle("D07  pyFP3D inviscid vs CFL3D Euler, ONERA M6 -- "
                 "Cp at the seven measured stations", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "d07_vs_cfl3d_euler.png"),
                dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
#  3D-2 — M 0.50 / alpha 3.06（另一条求解路径：solve_newton_lifting）
# ---------------------------------------------------------------------------

REF_DIR_M050 = os.path.join(str(REPO_ROOT), "cases", "reference_data", "cfl3d",
                            "euler_onera_m6_m050")
M_INF_SUB = 0.50
#: 记录带（标定）：实测 medium −6.89 %
SUBSONIC_CL_RECORD_MAX = 0.15
#: 展向结构：外/内 Cp RMS 之比，实测 1.87（跨声速 1.92）
SPANWISE_RATIO_MIN = 1.2


def _read_reference_m050(path=None):
    import csv
    p = path or os.path.join(REF_DIR_M050, "forces.csv")
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            out[r["level"]] = dict(cl=float(r["cl"]), cd=float(r["cd"]))
    if "L4" not in out:
        raise RuntimeError(f"{p}: no L4 row -- the 3D-2 reference is missing")
    return out


def _one_m050(level):
    """★ `solve_newton_lifting`，不是 `solve_newton_transonic` —— 见 docstring。
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
    from bench.recipes import NEWTON_M6_RECIPE
    from pyfp3d.constraints.wake import tip_taper_factors
    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.mesh.wake_cut import cut_wake
    from pyfp3d.meshgen.wing3d import B_SEMI
    from pyfp3d.post.section_cut import section_cp_curve
    from pyfp3d.post.surface import planform_area, wall_force_coefficients
    from pyfp3d.solve.newton import solve_newton_lifting
    #: ★ 网格缺失必须 **skip 而不是 error**（W0.1 / H1, 2026-09-06）。主腿经
    #: `_m6_case` 继承 `tests/E/test_E01…:159-160` 的守卫，**这条 m050 腿绕过了它**
    #: ⇒ `onera_m6/*.msh` 被 gitignore，干净 clone 上这里是硬 FileNotFoundError。
    p = os.path.join(str(REPO_ROOT), "cases", "meshes", "onera_m6",
                     f"{level}.msh")
    if not os.path.exists(p):
        pytest.skip(f"onera_m6/{level}.msh not generated; run "
                    "cases/meshes/onera_m6/generate_onera_m6.py")
    mc, wc = cut_wake(read_mesh(p))
    taper = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth",
                              0.05 * B_SEMI)
    kw = dict(NEWTON_M6_RECIPE["newton_kw"], tip_taper=taper)
    r = solve_newton_lifting(mc, wc, m_inf=M_INF_SUB, alpha_deg=ALPHA, **kw)
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], r["phi"],
                                alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF_SUB)
    out = dict(cl=float(f["cl"]), converged=bool(r.get("converged")),
               residual=float(np.asarray(r["residual_history"], float)[-1]),
               mach_max=float(np.sqrt(r["mach2_max"])))
    if level == "medium":
        for e in CP_ETAS:
            cur = section_cp_curve(mc, r["phi"], eta=e, b_semi=B_SEMI,
                                   m_inf=M_INF_SUB)
            rc = read_3d_cp(REF_DIR_M050, "L4", e, "none")
            out[f"cp_rms_upper_eta{e:.2f}"] = cp_rms(
                rc["upper"][0], rc["upper"][1],
                cur["x_upper"], cur["cp_upper"])[0]
    return out


@pytest.fixture(scope="module")
def runs_m050():
    return {lv: _one_m050(lv) for lv in LEVELS}


class TestSubsonicCondition3D2:
    r"""★ 3D-2：本工况填的是一个**完全空白**（M6 亚声速无任何外部参照）。"""

    def test_both_levels_converge(self, runs_m050):
        """★ 与跨声速那条不同，这里两档都收敛得很干净（|R| ~8e-15）——
        所以下面的读数不带"收敛与否"的混淆项。"""
        for lv in LEVELS:
            assert runs_m050[lv]["converged"], f"{lv} did not converge"
            assert runs_m050[lv]["residual"] < 1e-10

    def test_this_condition_is_not_actually_subcritical(self, runs_m050):
        """★★★ 锁住那个把我的预期推翻的事实：**M∞ = 0.50 而 M_max ≈ 1.0**。

        实测 0.9592（coarse）/ **0.9964**（medium）。⇒ "亚声速"是个误称，
        这条工况处在临界边缘，所以不能指望它像二维 M0.50 那样给 +0.04 %。
        若哪天 M_max 明显掉下来，本条会红 —— 那说明几何或配置变了，
        而上面那些读数的前提也就变了。"""
        assert 0.85 < runs_m050["medium"]["mach_max"] < 1.15, (
            f'M_max {runs_m050["medium"]["mach_max"]:.4f} is no longer near '
            f'critical -- the "not actually subcritical" reading that explains '
            f'this condition\'s behaviour needs re-checking')

    def test_lift_deficit_is_recorded_and_grows(self, runs_m050):
        """★★ 与 3D-1 同型：误差**随加密变大**（−1.54 % → −6.89 %）。
        **记录，不设门。** 若它开始收缩，本条会红 —— 那是好消息，
        届时应把这条腿升级为设门。"""
        want = _read_reference_m050()["L4"]["cl"]
        e = {lv: abs(runs_m050[lv]["cl"] - want) / abs(want) for lv in LEVELS}
        assert e["medium"] <= SUBSONIC_CL_RECORD_MAX, (
            f'cl error {e["medium"]*100:.2f} % exceeds the recording band')
        assert e["medium"] > e["coarse"], (
            f'the 3D-2 cl error now CONTRACTS ({e["coarse"]*100:.2f} % -> '
            f'{e["medium"]*100:.2f} %) -- an improvement; re-specify as gated')

    def test_the_spanwise_deficit_is_not_a_shock_phenomenon(self, runs_m050):
        """★★★ 本门最有价值的一条，而它需要**两个工况**才能说。

        上表面 Cp RMS 沿展向增大，两个工况的外/内比几乎相同：
        **M 0.50 给 1.87，M 0.8395 给 1.92** —— 而 M 0.50 上**没有激波**
        （M_max 0.996）。⇒ 那个越往翼尖越重的吸力亏损**与激波无关**，
        病因被缩小到"沿展向作用、外翼更重、与马赫数无关"。

        ★ 已知候选：`tip_taper`（−1.3 % cl_p，作用在翼尖，两条都开着）——
        但它只占 −6.89 % 的五分之一，**至多是部分原因，不归因**。"""
        r = [runs_m050["medium"][f"cp_rms_upper_eta{e:.2f}"] for e in CP_ETAS]
        assert all(np.isfinite(v) for v in r), r
        assert r[-1] / r[0] >= SPANWISE_RATIO_MIN, (
            f"the spanwise growth is gone at M 0.50 ({r[0]:.4f} at y/b 0.20 vs "
            f"{r[-1]:.4f} at 0.99, ratio {r[-1]/r[0]:.2f}) -- the deficit may "
            f"be shock-related after all; re-read against 3D-1's 1.92")

    def test_reference_reader_follows_a_perturbed_copy(self, tmp_path):
        import csv
        src = os.path.join(REF_DIR_M050, "forces.csv")
        with open(src) as fh:
            rows = list(csv.DictReader(fh))
            cols = list(rows[0].keys())
        for r in rows:
            if r["level"] == "L4":
                r["cl"] = f"{float(r['cl']) + 0.2222:.6f}"
        dst = tmp_path / "forces.csv"
        with open(dst, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        base = _read_reference_m050()["L4"]["cl"]
        pert = _read_reference_m050(str(dst))["L4"]["cl"]
        assert abs((pert - base) - 0.2222) < 1e-9
