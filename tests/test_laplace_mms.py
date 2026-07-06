"""
Gate G1.1: Method of Manufactured Solutions (MMS) convergence.

Manufactured solution phi_exact = sin(pi x) sin(pi y) sin(pi z) on the unit
cube, full Dirichlet BC from the exact solution, consistent FEM load vector
b_i = integral(f * N_i dV) with f = -Laplacian(phi_exact) = 3 pi^2 phi_exact,
via a standard 4-point degree-2-exact tet quadrature rule.

Note on the choice of exact solution: a harmonic *polynomial* (e.g.
x^2 - 0.5y^2 - 0.5z^2, body_source_rhs=None) was tried first and rejected --
on the structured Kuhn-triangulated cube from tests/mesh_utils.py, the P1
Galerkin stiffness matrix reduces to a stencil that reproduces harmonic
quadratics *exactly* (to machine precision) at every mesh level, the same
way a central finite difference is exact for a quadratic. That gives zero
signal about convergence order. sin*cos products (design.md Sec 10, V1) are
not polynomials, so no such superconvergence applies, and a real O(h^2) L2 /
O(h) H1 trend appears.
"""

import numpy as np
import pytest

from pyfp3d.mesh.metrics import compute_tet_volumes
from pyfp3d.solve.picard import solve_laplace

from .mesh_utils import cube_boundary_mask, generate_structured_cube_mesh

# 4-point, degree-2-exact symmetric tet quadrature rule (barycentric coords).
_A = (5.0 - np.sqrt(5.0)) / 20.0
_B = (5.0 + 3.0 * np.sqrt(5.0)) / 20.0
_QUAD_BARY = np.array(
    [
        [_B, _A, _A, _A],
        [_A, _B, _A, _A],
        [_A, _A, _B, _A],
        [_A, _A, _A, _B],
    ]
)
_QUAD_WEIGHT = np.full(4, 0.25)  # weights sum to 1; scaled by V_e per element


def phi_exact_fn(x, y, z):
    return np.sin(np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z)


def _source_fn(x, y, z):
    """f = -Laplacian(phi_exact) = 3 pi^2 sin(pi x) sin(pi y) sin(pi z)."""
    return 3.0 * np.pi**2 * phi_exact_fn(x, y, z)


def _consistent_load_vector(nodes: np.ndarray, elements: np.ndarray, volumes: np.ndarray) -> np.ndarray:
    """Assemble b_i = sum_e integral_e(f * N_i dV) via the 4-point tet rule."""
    n_nodes = len(nodes)
    b = np.zeros(n_nodes, dtype=np.float64)
    elem_nodes = nodes[elements]  # (n_tets, 4, 3)

    for q in range(4):
        bary = _QUAD_BARY[q]
        x_q = np.einsum("i,eij->ej", bary, elem_nodes)
        f_q = _source_fn(x_q[:, 0], x_q[:, 1], x_q[:, 2])
        w = _QUAD_WEIGHT[q] * volumes
        for i in range(4):
            np.add.at(b, elements[:, i], w * f_q * bary[i])

    return b


def _lumped_nodal_volumes(elements: np.ndarray, volumes: np.ndarray, n_nodes: int) -> np.ndarray:
    lumped = np.zeros(n_nodes, dtype=np.float64)
    for local in range(4):
        np.add.at(lumped, elements[:, local], volumes / 4.0)
    return lumped


def run_mms_case(n: int) -> dict:
    """Solve the MMS problem on an n x n x n structured cube; return L2 error, h, cg iters."""
    nodes, elements = generate_structured_cube_mesh(n=n, L=1.0)
    boundary = cube_boundary_mask(nodes, L=1.0)
    dirichlet_nodes = np.where(boundary)[0]

    phi_exact = phi_exact_fn(nodes[:, 0], nodes[:, 1], nodes[:, 2])
    volumes = compute_tet_volumes(nodes, elements)
    b = _consistent_load_vector(nodes, elements, volumes)

    result = solve_laplace(
        nodes, elements, dirichlet_nodes, phi_exact[dirichlet_nodes],
        body_source_rhs=b, rtol=1e-12, maxiter=3000,
    )
    phi = result["phi"]

    lumped = _lumped_nodal_volumes(elements, volumes, len(nodes))
    l2_error = float(np.sqrt(np.sum(lumped * (phi - phi_exact) ** 2) / np.sum(lumped)))

    return {
        "n": n,
        "h": 1.0 / n,
        "l2_error": l2_error,
        "n_cg_iterations": result["n_cg_iterations"],
        "residual_norm": result["residual_norm"],
    }


class TestLaplaceMMS:
    """Gate G1.1: MMS L2 convergence slope >= 1.9 over 3 mesh levels."""

    def test_mms_convergence_slope(self):
        levels = [run_mms_case(n) for n in (4, 8, 16)]

        for level in levels:
            assert level["residual_norm"] < 1e-8, (
                f"n={level['n']}: linear solve not converged, "
                f"residual_norm={level['residual_norm']:.3e}"
            )

        h = np.array([lvl["h"] for lvl in levels])
        err = np.array([lvl["l2_error"] for lvl in levels])
        slope, _ = np.polyfit(np.log(h), np.log(err), 1)

        assert slope >= 1.9, f"MMS L2 convergence slope {slope:.3f} < 1.9 (errors: {err})"


class TestLaplaceMMSArtifacts:
    """Generate visual artifacts for G1.1."""

    def test_export_mms_convergence(self, gate_artifacts_dir):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        levels = [run_mms_case(n) for n in (4, 8, 16)]
        h = np.array([lvl["h"] for lvl in levels])
        err = np.array([lvl["l2_error"] for lvl in levels])
        slope, intercept = np.polyfit(np.log(h), np.log(err), 1)

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.loglog(h, err, "o-", label="measured L2 error")
        fit = np.exp(intercept) * h**slope
        ax.loglog(h, fit, "k--", label=f"fit slope = {slope:.2f}")
        ax.set_xlabel("h (cube edge / n)")
        ax.set_ylabel("L2 error")
        ax.set_title("G1.1: Laplace MMS convergence")
        ax.legend()
        ax.grid(True, which="both", alpha=0.3)

        output_file = gate_artifacts_dir / "mms_convergence_loglog.png"
        fig.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close(fig)
        assert output_file.exists()

        csv_file = gate_artifacts_dir / "summary.csv"
        with open(csv_file, "w") as f:
            f.write("n,h,l2_error,cg_iterations\n")
            for lvl in levels:
                f.write(f"{lvl['n']},{lvl['h']},{lvl['l2_error']:.6e},{lvl['n_cg_iterations']}\n")
            f.write(f"fit_slope,{slope:.4f},,\n")
        assert csv_file.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
