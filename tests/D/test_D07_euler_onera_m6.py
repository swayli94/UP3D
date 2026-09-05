r"""D07 — pyFP3D 无粘 vs CFL3D Euler — ONERA M6（3D-1）。

工况：**M 0.8395 / α 3.06**，AGARD AR-138 TEST 2308 逐字（α 实验值不修正）。
参考 `euler_onera_m6/` 四档 + REF，本门读 **L4**（最细，11.278 M 点）。

★ 占位期的 docstring 还列了 M0.50/α3.06 —— **参考数据集里没有那个工况**，
所以本门不含它；要补就得先补参考。

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
from tests._gate_evidence import assert_matches_committed, fmt

REF_DIR = os.path.join(str(REPO_ROOT), "cases", "reference_data", "cfl3d",
                       "euler_onera_m6")
M_INF, ALPHA = 0.8395, 3.06
LEVELS = ("coarse", "medium")
ETAS = (0.44, 0.65, 0.90)

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
    r, forces, shocks, wall = _m6_case(level, m_inf=M_INF, alpha=ALPHA)
    mc, _wc = cut_wake(read_mesh(os.path.join(
        str(REPO_ROOT), "cases", "meshes", "onera_m6", f"{level}.msh")))
    return dict(
        cl=float(forces["cl"]), shocks={k: float(v) for k, v in shocks.items()},
        converged=bool(r.get("converged")),
        residual=float(np.asarray(r["residual_history"], float)[-1]),
        mach_max=float(np.sqrt(r["mach2_max"])),
        s_ref_discrete=float(planform_area(mc.nodes,
                                           mc.boundary_faces["wall"])))


@pytest.fixture(scope="module")
def runs():
    return {lv: _one(lv) for lv in LEVELS}


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
    MEASURED = ("cl", "residual", "mach_max", "s_ref_discrete")

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
                       + [fmt(d["shocks"][e]) for e in ETAS]
                       + [f'{d["shocks"][e] - sh[e]["x"]:+.6f}' for e in ETAS])

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    x = np.arange(len(ETAS))
    for k, lv in enumerate(LEVELS):
        ax[0].bar(x + k * 0.35 - 0.18,
                  [abs(runs[lv]["shocks"][e] - sh[e]["x"]) for e in ETAS],
                  0.33, label=lv)
    ax[0].set_xticks(x), ax[0].set_xticklabels([f"eta {e}" for e in ETAS])
    ax[0].set_ylabel("|shock offset| vs CFL3D L4")
    ax[0].set_title("GATED: contracts 9.5x / 2.2x / 1.8x", fontsize=9)
    ax[0].legend(fontsize=8), ax[0].grid(alpha=.3, axis="y")
    ax[1].bar(LEVELS, [abs(runs[lv]["cl"] - ref["L4"]["cl"])
                       / ref["L4"]["cl"] * 100 for lv in LEVELS],
              color="#c0392b")
    ax[1].set_ylabel("|d cl| %")
    ax[1].set_title("RECORDED: cl error GROWS, 6.06 % -> 9.17 %", fontsize=9,
                    color="#c0392b")
    ax[1].grid(alpha=.3, axis="y")
    fig.suptitle("D07  pyFP3D inviscid vs CFL3D Euler, ONERA M6 TEST 2308 -- "
                 "the shock converges while the lift does not", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "d07_vs_cfl3d_euler.png"),
                dpi=130)
    plt.close(fig)
