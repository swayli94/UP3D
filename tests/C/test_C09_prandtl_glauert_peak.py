"""C09 -- M0.3 球的 Cp 峰对 Prandtl-Glauert（C 类）。

★★★ 这是本项目**唯一已存在的可压缩检验**，而我在 phase-6 归类时**漏掉过它**
（正则没匹配 pg），还把 PG 列成「缺、可选」。
★ 它**不是精确解**，是小扰动**渐近**关系 ⇒ 只能判趋势/量级，**不能当真值**。
"""
import csv
import numpy as np
import pytest
from tests._gate_evidence import assert_matches_committed, fmt
from tests.conftest import REPO_ROOT, gate_figures_enabled
import matplotlib
import matplotlib.pyplot as plt
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.physics.isentropic import pressure_coefficient
from pyfp3d.post.surface import wall_tangential_gradient_quadratic
from pyfp3d.solve.picard import (
    solve_laplace,
    solve_laplace_lifting,
    solve_subsonic,
    solve_subsonic_lifting,
)
from tests._sphere_case import M_INF_SPHERE, run_sphere_pair, sphere_medium


class TestG31SphereCompressible:
    r"""★★★ 判据**常跑**（W0.4 / H19，2026-09-06）。

    改这一条之前，整个类挂在 `skipif(not gate_figures_enabled())` 下，于是
    **PG 的 2 % 判据在平时跑里根本不执行** —— 与其余 C 门的形状正好相反
    （它们判据常跑、出图 opt-in）。而计算早就在共享 fixture `sphere_medium`
    里，所以把断言移出来是**零额外计算**。
    ★ 图与 CSV 仍只在 `PYFP3D_GATE_FIGURES=1` 时写盘（见文件末的导出腿），
    并且**用同一个 module 级 fixture** ⇒ 图与断言同源，由构造保证。
    """

    def test_g31_cp_peak_vs_pg(self, sphere_medium):
        """G3.1: |Cp_peak(FP, M=0.3) - Cp_peak(incompressible)/beta| < 2%."""
        case = sphere_medium
        beta = np.sqrt(1.0 - M_INF_SPHERE**2)
        cp_peak_pg = float(case["cp_inc"].min()) / beta
        cp_peak_c = float(case["cp_c"].min())
        rel = abs(cp_peak_c - cp_peak_pg) / abs(cp_peak_pg)

        assert case["result_c"]["converged"]
        assert rel < 0.02, (
            f"Cp peak {cp_peak_c:.5f} vs PG-corrected {cp_peak_pg:.5f} "
            f"({100 * rel:.2f}% >= 2%)"
        )


class TestCommittedEvidenceIsLoadBearing:
    r"""★★★ 新鲜计算 vs 已提交 `G3.1/summary.csv`。设计与四个坑见 `tests/_gate_evidence.py`。

    ★★ **本门此前在平时跑里一条活断言都没有** —— 实测 `pytest tests/C/test_C09...`
    只报 `1 skipped`：唯一的测试类整个挂在 `skipif(not gate_figures_enabled())` 下。
    ⇒ **本条是它的第一条活断言**，也是本轮转换里价值最高的一处。
    ★ 零额外计算：用共享 fixture `sphere_medium`。
    """

    def test_fresh_run_reproduces_the_committed_summary(self, sphere_medium):
        import os

        case = sphere_medium
        beta = np.sqrt(1.0 - M_INF_SPHERE ** 2)
        cp_inc = float(case["cp_inc"].min())
        fresh = {(): {"cp_peak_incompressible": cp_inc,
                      "cp_peak_pg_corrected": cp_inc / beta,
                      "cp_peak_compressible": float(case["cp_c"].min()),
                      "mach2_max": float(case["result_c"]["mach2_max"])}}
        n = assert_matches_committed(
            os.path.join(str(REPO_ROOT), "cases", "gates",
                         "C09_prandtl_glauert_peak", "G3.1"),
            fresh, ("cp_peak_incompressible", "cp_peak_pg_corrected",
                    "cp_peak_compressible", "mach2_max"), key_of=lambda r: (),
            refresh_hint="PYFP3D_GATE_FIGURES=1 pytest "
                         "tests/C/test_C09_prandtl_glauert_peak.py")
        assert n == 4, f"比了 {n} 个数，应为 4"


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="C/D 证据图只在 PYFP3D_GATE_FIGURES=1 时写盘 —— "
                           "平时不写，避免每次 pytest 都脏工作树（2026-08-24）")
def test_export_c09_evidence(sphere_medium, gate_evidence_dir):
    """★ 只导出，不判定 —— 判据在 `TestG31SphereCompressible`（常跑）。

    ★★ 同一个 module 级 fixture ⇒ 图 / CSV / 断言来自**同一次计算**，构造保证。
    """
    case = sphere_medium
    beta = np.sqrt(1.0 - M_INF_SPHERE**2)
    cp_peak_inc = float(case["cp_inc"].min())
    cp_peak_pg = cp_peak_inc / beta
    cp_peak_c = float(case["cp_c"].min())
    rel = abs(cp_peak_c - cp_peak_pg) / abs(cp_peak_pg)
    rc = case["result_c"]

    gate_dir = gate_evidence_dir / "G3.1"
    gate_dir.mkdir(parents=True, exist_ok=True)
    with open(gate_dir / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cp_peak_incompressible", "cp_peak_pg_corrected",
                    "cp_peak_compressible", "rel_diff_pct", "n_picard",
                    "picard_converged", "mach2_max"])
        w.writerow([fmt(cp_peak_inc), fmt(cp_peak_pg),
                    fmt(cp_peak_c), fmt(100 * rel),
                    rc["n_picard"], rc["converged"], fmt(rc['mach2_max'])])

    # V3.1: Cp(theta) line cut -- incompressible, PG-corrected, compressible;
    # the amplification must be symmetric fore/aft.
    theta = np.degrees(np.arccos(np.clip(case["cos_theta"], -1, 1)))
    order = np.argsort(theta)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(theta[order], case["cp_inc"][order], ".", ms=2, color="0.6",
            label="incompressible (same mesh)")
    ax.plot(theta[order], case["cp_inc"][order] / beta, "-", lw=1.0,
            color="tab:green", label="PG-corrected incompressible")
    ax.plot(theta[order], case["cp_c"][order], ".", ms=2,
            color="tab:red", label=f"full potential M={M_INF_SPHERE}")
    ax.invert_yaxis()
    ax.set_xlabel("theta (deg, from +x stagnation)")
    ax.set_ylabel("Cp")
    ax.set_title(f"V3.1 sphere Cp at M={M_INF_SPHERE} (medium): "
                 f"peak diff {100 * rel:.2f}% vs PG")
    ax.legend()
    fig.savefig(gate_dir / "v3_1_sphere_cp.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
