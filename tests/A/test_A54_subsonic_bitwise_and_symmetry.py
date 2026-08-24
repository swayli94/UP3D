"""A54 -- 亚声速路径的位相同与对称性（A 类）。

★ 从 test_p3_subsonic.py 拆出。断言：M0 时装配矩阵**逐位**等于 Laplace 的、
lifting/nonlifting 与 P1/P2 路径**逐位**相同、Cp 前后对称、Picard 残差单调。
★★ **跨路径一致性归 A**（使用者裁决 2026-08-24，边界 ③）—— 它检验库内部的一致性，无外部真值。
★ Prandtl-Glauert 那条在 tests/C/test_C09_prandtl_glauert_peak.py。
"""
import csv
import numpy as np
import pytest
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
    def test_g31_picard_monotone_and_symmetric(self, sphere_medium,
                                               tmp_path):
        """Non-lifting Picard: monotone residual to its floor, and the
        compressible amplification keeps fore/aft symmetry (V3.1 check)."""
        case = sphere_medium
        h = case["result_c"]["residual_history"]
        floor = 10.0 * h[-1]
        for i in range(len(h) - 1):
            if h[i] > floor:
                assert h[i + 1] <= h[i], (
                    f"residual rose above floor at iter {i}: "
                    f"{h[i]:.3e} -> {h[i + 1]:.3e}"
                )

        gate_dir = tmp_path / "G3.1"
        gate_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.semilogy(range(1, len(h) + 1), h, "o-", ms=4)
        ax.set_xlabel("Picard iteration")
        ax.set_ylabel("||R||_inf (free dofs)")
        ax.set_title("sphere M=0.3 Picard convergence (medium)")
        fig.savefig(gate_dir / "sphere_picard_residual.png", dpi=150,
                    bbox_inches="tight")
        plt.close(fig)

        # Fore/aft symmetry of the compressible Cp: mirror-pair scatter
        # bounded by the mesh's own asymmetry (loose 5% of dynamic range).
        ct, cp = case["cos_theta"], case["cp_c"]
        front = ct > 0.2
        cp_mirror = np.interp(-ct[front], np.sort(ct),
                              cp[np.argsort(ct)])
        assert np.median(np.abs(cp[front] - cp_mirror)) < 0.05 * np.ptp(cp)


class TestG33BitIdenticalLaplaceLimit:
    def test_matrix_bits_at_m0(self):
        """A(rho(M=0, any phi)) == Laplace A, bitwise."""
        from tests.conftest import REPO_ROOT
        from pyfp3d.kernels.jacobian import PicardOperator
        from pyfp3d.kernels.residual import assemble_stiffness_matrix
        from pyfp3d.physics.isentropic import density_field

        m = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh")
        A_lap = assemble_stiffness_matrix(m.nodes, m.elements)
        op = PicardOperator(m.nodes, m.elements)
        phi = np.sin(m.nodes[:, 0]) + m.nodes[:, 1] ** 2  # arbitrary field
        _, q2 = op.velocities(phi)
        A_m0 = op.assemble_matrix(density_field(q2, 0.0))
        assert np.array_equal(A_lap.data, A_m0.data)
        assert np.array_equal(A_lap.indices, A_m0.indices)

    def test_lifting_m0_bitwise_vs_p2(self):
        """solve_subsonic_lifting(M=0) == solve_laplace_lifting, bitwise,
        for both the fixed-Gamma solve and the full secant Kutta loop."""
        from tests.conftest import REPO_ROOT

        m = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh")
        mc, wc = cut_wake(m)

        a = solve_laplace_lifting(mc, wc, alpha_deg=0.0, gamma_fixed=0.3)
        b = solve_subsonic_lifting(mc, wc, m_inf=0.0, alpha_deg=0.0,
                                   gamma_fixed=0.3)
        assert np.array_equal(a["phi"], b["phi"])
        assert b["n_picard"] == 1 and b["converged"]

        a = solve_laplace_lifting(mc, wc, alpha_deg=4.0)
        b = solve_subsonic_lifting(mc, wc, m_inf=0.0, alpha_deg=4.0)
        assert np.array_equal(a["phi"], b["phi"])
        assert np.array_equal(a["gamma"], b["gamma"])
        assert b["n_picard"] == 1

    def test_nonlifting_m0_bitwise_vs_p1(self):
        """solve_subsonic(M=0) == solve_laplace, bitwise, at matched CG
        controls (rho == 1 exactly makes iteration 1 the P1 solve)."""
        from tests.conftest import REPO_ROOT

        mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "sphere_shell" / "coarse.msh")
        ff = np.unique(mesh.boundary_faces["farfield"])
        phi_ff = mesh.nodes[ff, 0]
        a = solve_laplace(mesh.nodes, mesh.elements, ff, phi_ff,
                          rtol=1e-10, maxiter=3000)
        b = solve_subsonic(mesh.nodes, mesh.elements, ff, phi_ff,
                           m_inf=0.0, rtol=1e-10, maxiter=3000)
        assert np.array_equal(a["phi"], b["phi"])
        assert b["n_picard"] == 1 and b["converged"]
