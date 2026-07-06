"""
P1 driver: Laplace (rho == 1) is linear, so the Picard loop degenerates to a
single assemble + solve -- no outer iteration, no relaxation.

Reference: docs/roadmap.md P1 ("Picard-degenerate driver, single linear
solve"); design.md Sec 8 (full Picard loop lands in P3+).
"""

from typing import Dict, Optional

import numpy as np

from pyfp3d.kernels.residual import assemble_residual, assemble_stiffness_matrix
from pyfp3d.solve.linear import apply_dirichlet, solve_cg_amg


def solve_laplace(
    nodes: np.ndarray,
    elements: np.ndarray,
    dirichlet_nodes: np.ndarray,
    dirichlet_values: np.ndarray,
    body_source_rhs: Optional[np.ndarray] = None,
    rtol: float = 1e-10,
    maxiter: int = 500,
) -> Dict[str, object]:
    """
    Solve the Laplace problem A phi = b with Dirichlet BCs (far field) and
    natural (do-nothing) BCs elsewhere (solid walls, symmetry planes).

    Args:
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) tetrahedral connectivity
        dirichlet_nodes: node indices with prescribed phi (e.g. far field)
        dirichlet_values: prescribed phi values at dirichlet_nodes
        body_source_rhs: optional (n_nodes,) assembled RHS vector -- an MMS
            load vector, or a wall-flux boundary correction
            (solve/wall_correction.py); the physical full-potential equation
            itself has no volumetric source term
        rtol, maxiter: CG convergence controls (see solve.linear.solve_cg_amg)

    Returns:
        dict with keys: phi, n_cg_iterations, residual_norm
    """
    n_nodes = len(nodes)
    A = assemble_stiffness_matrix(nodes, elements)
    b = np.zeros(n_nodes, dtype=np.float64) if body_source_rhs is None else np.asarray(
        body_source_rhs, dtype=np.float64
    )

    A_free, b_free, free, phi = apply_dirichlet(A, b, dirichlet_nodes, dirichlet_values)
    x_free, n_iter = solve_cg_amg(A_free, b_free, rtol=rtol, maxiter=maxiter)
    phi[free] = x_free

    # Residual of the equilibrium equation b - A@phi == 0, restricted to the
    # free dofs. Dirichlet (far-field) rows are excluded: assemble_residual()
    # returns the raw operator A@phi with no boundary condition applied, so
    # at a Dirichlet node it is the natural-BC flux imbalance -- an O(1)
    # quantity that never vanishes and isn't part of the system actually
    # being solved (that row is overridden by the prescribed value instead).
    residual_full = b - assemble_residual(nodes, elements, phi)
    residual_norm = float(np.max(np.abs(residual_full[free])))

    return {
        "phi": phi,
        "n_cg_iterations": n_iter,
        "residual_norm": residual_norm,
    }
