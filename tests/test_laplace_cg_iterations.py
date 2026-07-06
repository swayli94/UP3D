"""
Gate G1.2 (formerly G1.3): AMG-preconditioned CG iteration count is roughly mesh-independent.

Runs the same smooth (non-polynomial) MMS-style problem from
test_laplace_mms.py on 3 structured-cube levels and checks that the CG+AMG
iteration count does not blow up with problem size the way an
unpreconditioned Krylov method would (condition number of the plain Laplace
stiffness matrix scales like O(h^-2), i.e. iteration count like O(1/h) for
CG without a multigrid preconditioner). AMG should instead give iteration
counts that grow only slowly (empirically ~log(1/h)) with mesh size.
"""

import numpy as np
import pytest

from .mesh_utils import cube_boundary_mask, generate_structured_cube_mesh
from .test_laplace_mms import phi_exact_fn
from pyfp3d.solve.picard import solve_laplace

LEVELS = (8, 16, 32)


def run_cg_case(n: int) -> dict:
    nodes, elements = generate_structured_cube_mesh(n=n, L=1.0)
    boundary = cube_boundary_mask(nodes, L=1.0)
    dirichlet_nodes = np.where(boundary)[0]
    phi_exact = phi_exact_fn(nodes[:, 0], nodes[:, 1], nodes[:, 2])

    result = solve_laplace(
        nodes, elements, dirichlet_nodes, phi_exact[dirichlet_nodes], rtol=1e-10, maxiter=2000,
    )
    return {
        "n": n,
        "n_nodes": len(nodes),
        "n_cg_iterations": result["n_cg_iterations"],
        "residual_norm": result["residual_norm"],
    }


class TestCGMeshIndependence:
    """Gate G1.2: CG+AMG iteration count roughly mesh-independent."""

    def test_iteration_count_bounded_growth(self):
        levels = [run_cg_case(n) for n in LEVELS]

        for level in levels:
            assert level["residual_norm"] < 1e-7, (
                f"n={level['n']}: not converged, residual_norm={level['residual_norm']:.3e}"
            )

        iters = np.array([lvl["n_cg_iterations"] for lvl in levels])
        n_nodes = np.array([lvl["n_nodes"] for lvl in levels])

        # An 8x/level increase in node count (each level doubles resolution
        # in 3D) should not come with anywhere near a proportional increase
        # in iteration count if AMG is doing its job -- cap the whole-range
        # growth ratio well below what an unpreconditioned method would show.
        growth_ratio = iters.max() / iters.min()
        assert growth_ratio < 2.0, (
            f"CG+AMG iteration count grew {growth_ratio:.2f}x over a "
            f"{n_nodes.max() / n_nodes.min():.0f}x increase in node count "
            f"(iters={iters.tolist()}, n_nodes={n_nodes.tolist()}) -- looks "
            f"mesh-dependent, not mesh-independent."
        )


class TestCGIterationArtifacts:
    """Generate visual artifacts for G1.2."""

    def test_export_cg_convergence_overlay(self, gate_artifacts_dir):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        import pyamg
        import scipy.sparse.linalg as spla

        from pyfp3d.kernels.residual import assemble_stiffness_matrix
        from pyfp3d.solve.linear import apply_dirichlet

        fig, ax = plt.subplots(figsize=(7, 6))
        summary_rows = []

        for n in LEVELS:
            nodes, elements = generate_structured_cube_mesh(n=n, L=1.0)
            boundary = cube_boundary_mask(nodes, L=1.0)
            dirichlet_nodes = np.where(boundary)[0]
            phi_exact = phi_exact_fn(nodes[:, 0], nodes[:, 1], nodes[:, 2])

            A = assemble_stiffness_matrix(nodes, elements)
            b = np.zeros(len(nodes))
            A_free, b_free, free, _ = apply_dirichlet(A, b, dirichlet_nodes, phi_exact[dirichlet_nodes])

            ml = pyamg.smoothed_aggregation_solver(A_free)
            M = ml.aspreconditioner()
            residuals = []
            spla.cg(A_free, b_free, M=M, rtol=1e-10, maxiter=2000,
                    callback=lambda xk: residuals.append(np.linalg.norm(b_free - A_free @ xk)))

            residuals = np.array(residuals)
            residuals = residuals / residuals[0] if len(residuals) else residuals
            ax.semilogy(np.arange(len(residuals)), np.maximum(residuals, 1e-16),
                        "o-", label=f"n={n} ({len(nodes)} nodes)")
            summary_rows.append((n, len(nodes), len(residuals)))

        ax.set_xlabel("CG iteration")
        ax.set_ylabel("relative residual")
        ax.set_title("G1.2: CG+AMG convergence history (mesh-independence)")
        ax.legend()
        ax.grid(True, which="both", alpha=0.3)

        output_file = gate_artifacts_dir / "cg_convergence_overlay.png"
        fig.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close(fig)
        assert output_file.exists()

        csv_file = gate_artifacts_dir / "summary.csv"
        with open(csv_file, "w") as f:
            f.write("n,n_nodes,n_cg_iterations\n")
            for n, n_nodes, n_iter in summary_rows:
                f.write(f"{n},{n_nodes},{n_iter}\n")
        assert csv_file.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
