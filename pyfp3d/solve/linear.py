"""
Linear algebra for SPD Laplace/Picard systems: Dirichlet elimination plus
AMG-preconditioned CG.

Reference: design.md Sec 7 point 6 (SciPy/PyAMG land), Sec 8 (Nonlinear
solution strategy -- linear solve step).
"""

from typing import Tuple

import numpy as np
import pyamg
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def apply_dirichlet(
    A: sp.spmatrix,
    b: np.ndarray,
    dirichlet_nodes: np.ndarray,
    dirichlet_values: np.ndarray,
) -> Tuple[sp.csr_matrix, np.ndarray, np.ndarray, np.ndarray]:
    """
    Eliminate Dirichlet dofs via a principal submatrix (free-dof) reduction.

    A principal submatrix of an SPD matrix is SPD, so the reduced system
    stays solvable with CG+AMG -- simpler and more robust than in-place
    row/column zeroing on a sparse CSR matrix.

    Args:
        A: (n, n) SPD sparse matrix
        b: (n,) right-hand side
        dirichlet_nodes: indices with prescribed values
        dirichlet_values: prescribed values at dirichlet_nodes

    Returns:
        (A_free, b_free, free_idx, phi_dirichlet):
          A_free, b_free: reduced SPD system on the free dofs
          free_idx: indices of the free dofs (into the full-size arrays)
          phi_dirichlet: (n,) array with dirichlet_values scattered in
              (zero at free dofs -- fill those in after solving)
    """
    n = A.shape[0]
    is_dirichlet = np.zeros(n, dtype=bool)
    is_dirichlet[dirichlet_nodes] = True
    free = np.where(~is_dirichlet)[0]

    phi_dirichlet = np.zeros(n, dtype=np.float64)
    phi_dirichlet[dirichlet_nodes] = dirichlet_values

    A_csr = A.tocsr()
    b_free = b[free] - A_csr[free][:, dirichlet_nodes] @ dirichlet_values
    A_free = A_csr[free][:, free].tocsr()

    return A_free, b_free, free, phi_dirichlet


def solve_cg_amg(
    A: sp.spmatrix,
    b: np.ndarray,
    rtol: float = 1e-10,
    maxiter: int = 500,
) -> Tuple[np.ndarray, int]:
    """
    Solve A x = b with CG, preconditioned by PyAMG smoothed aggregation.

    Args:
        A: (n, n) SPD sparse matrix
        b: (n,) right-hand side
        rtol: relative residual tolerance (scipy.sparse.linalg.cg convention)
        maxiter: maximum CG iterations

    Returns:
        (x, n_iterations)

    Raises:
        RuntimeError if CG fails to converge within maxiter.
    """
    ml = pyamg.smoothed_aggregation_solver(A)
    M = ml.aspreconditioner()

    n_iter = [0]

    def _count(_xk):
        n_iter[0] += 1

    x, info = spla.cg(A, b, M=M, rtol=rtol, maxiter=maxiter, callback=_count)

    if info != 0:
        raise RuntimeError(
            f"CG did not converge (info={info}) after {n_iter[0]} iterations"
        )

    return x, n_iter[0]
