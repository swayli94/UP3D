"""C09 -- M0.3 球的 Cp 峰对 Prandtl-Glauert（C 类）。

★★★ 这是本项目**唯一已存在的可压缩检验**，而我在 phase-6 归类时**漏掉过它**
（正则没匹配 pg），还把 PG 列成「缺、可选」。
★ 它**不是精确解**，是小扰动**渐近**关系 ⇒ 只能判趋势/量级，**不能当真值**。
"""
import csv
import numpy as np
import pytest
from tests.conftest import gate_figures_enabled
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
    @pytest.mark.skipif(not gate_figures_enabled(),
                        reason="C/D 证据图只在 PYFP3D_GATE_FIGURES=1 时写盘 —— 平时不写，避免每次 pytest 都脏工作树（2026-08-24）")
    def test_g31_cp_peak_vs_pg(self, sphere_medium, gate_evidence_dir):
        """G3.1: |Cp_peak(FP, M=0.3) - Cp_peak(incompressible)/beta| < 2%."""
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
            w.writerow([f"{cp_peak_inc:.6f}", f"{cp_peak_pg:.6f}",
                        f"{cp_peak_c:.6f}", f"{100 * rel:.4f}",
                        rc["n_picard"], rc["converged"],
                        f"{rc['mach2_max']:.4f}"])

        # V3.1: Cp(theta) line cut -- incompressible, PG-corrected,
        # compressible; the amplification must be symmetric fore/aft.
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
        fig.savefig(gate_dir / "v3_1_sphere_cp.png", dpi=150,
                    bbox_inches="tight")
        plt.close(fig)

        assert rc["converged"]
        assert rel < 0.02, (
            f"Cp peak {cp_peak_c:.5f} vs PG-corrected {cp_peak_pg:.5f} "
            f"({100 * rel:.2f}% >= 2%)"
        )
