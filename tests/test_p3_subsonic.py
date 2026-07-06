"""
P3 gates G3.1 and G3.3 (roadmap "P3 -- Subsonic compressible").

G3.1: sphere at M_inf = 0.3 -- compressible Cp suction peak vs the
Prandtl-Glauert-corrected incompressible peak computed on the SAME mesh
with the SAME surface recovery, < 2%. Comparing solver-vs-solver on one
mesh cancels the common discretization/recovery error (the known G1.6
flat-facet wall bias hits both solves identically), so this isolates the
compressibility machinery, which is what P3 adds. Both solves use plain
freestream Dirichlet phi = x at the far field (r_out = 20 R; the doublet
correction ~1/r^2 is common-mode).

G3.3: the nu == 0 / M_inf -> 0 path is bit-identical to the P1/P2
Laplace drivers -- same assembled matrix bits, same solve bits, same
Kutta trajectory (structural: rho == 1.0 bitwise from the density law,
beta == 1.0 bitwise in the vortex far field, seeded AMG setup in
solve/linear.py::build_amg_preconditioner). The rest of G3.3 ("all P1/P2
gates still green") is the full suite itself.

Artifacts: artifacts/G3.1/ (V3.1 Cp line cut + summary.csv + V3.3-style
residual plot for the sphere).
"""

import csv

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")
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

M_INF_SPHERE = 0.3


def run_sphere_pair(mesh_path):
    """Incompressible + compressible solves on one sphere mesh, wall Cp
    from the quadratic tangential recovery for both."""
    mesh = read_mesh(mesh_path)
    nodes, elements = mesh.nodes, mesh.elements
    wall = mesh.boundary_faces["wall"]
    wn = np.unique(wall)
    ff = np.unique(mesh.boundary_faces["farfield"])
    phi_ff = nodes[ff, 0]  # freestream Dirichlet, same for both solves

    r_inc = solve_laplace(nodes, elements, ff, phi_ff, rtol=1e-11, maxiter=3000)
    g = wall_tangential_gradient_quadratic(nodes, wall, r_inc["phi"])
    cp_inc = 1.0 - np.sum(g[wn] ** 2, axis=1)

    r_c = solve_subsonic(
        nodes, elements, ff, phi_ff, m_inf=M_INF_SPHERE,
        phi_init=nodes[:, 0].copy(), rtol=1e-12,
    )
    g = wall_tangential_gradient_quadratic(nodes, wall, r_c["phi"])
    q2 = np.sum(g[wn] ** 2, axis=1)
    cp_c = np.array([pressure_coefficient(q, M_INF_SPHERE) for q in q2])

    r = np.linalg.norm(nodes[wn], axis=1)
    cos_theta = nodes[wn, 0] / r
    return {
        "cos_theta": cos_theta, "cp_inc": cp_inc, "cp_c": cp_c,
        "result_c": r_c,
    }


@pytest.fixture(scope="module")
def sphere_medium(request):
    from .conftest import REPO_ROOT
    return run_sphere_pair(
        REPO_ROOT / "cases" / "meshes" / "sphere_shell" / "medium.msh"
    )


class TestG31SphereCompressible:
    def test_g31_cp_peak_vs_pg(self, sphere_medium, artifacts_dir):
        """G3.1: |Cp_peak(FP, M=0.3) - Cp_peak(incompressible)/beta| < 2%."""
        case = sphere_medium
        beta = np.sqrt(1.0 - M_INF_SPHERE**2)
        cp_peak_inc = float(case["cp_inc"].min())
        cp_peak_pg = cp_peak_inc / beta
        cp_peak_c = float(case["cp_c"].min())
        rel = abs(cp_peak_c - cp_peak_pg) / abs(cp_peak_pg)

        rc = case["result_c"]
        gate_dir = artifacts_dir / "G3.1"
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

    def test_g31_picard_monotone_and_symmetric(self, sphere_medium,
                                               artifacts_dir):
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

        gate_dir = artifacts_dir / "G3.1"
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
        from .conftest import REPO_ROOT
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
        from .conftest import REPO_ROOT

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
        from .conftest import REPO_ROOT

        mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "sphere_shell" / "coarse.msh")
        ff = np.unique(mesh.boundary_faces["farfield"])
        phi_ff = mesh.nodes[ff, 0]
        a = solve_laplace(mesh.nodes, mesh.elements, ff, phi_ff,
                          rtol=1e-10, maxiter=3000)
        b = solve_subsonic(mesh.nodes, mesh.elements, ff, phi_ff,
                           m_inf=0.0, rtol=1e-10, maxiter=3000)
        assert np.array_equal(a["phi"], b["phi"])
        assert b["n_picard"] == 1 and b["converged"]
