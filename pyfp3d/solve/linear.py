"""
Linear algebra: Dirichlet elimination, AMG-preconditioned CG for the SPD
Laplace/Picard systems, and (P8) preconditioned GMRES for the exact Newton
Jacobian -- which is NONSYMMETRIC and indefinite inside supersonic zones
(design.md Sec 6.3: the Term-3 upstream coupling has no transpose partner),
so CG does not apply there.

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


def build_amg_preconditioner(A: sp.spmatrix):
    """Smoothed-aggregation AMG preconditioner with a PINNED RNG seed.

    pyamg's prolongation smoother estimates the spectral radius from an
    UNSEEDED np.random starting vector, so two setups on a bit-identical
    matrix produce (slightly) different hierarchies -- measured 2e-11 phi
    scatter between repeated identical Laplace solves. Pinning the seed
    (and restoring the caller's RNG state) makes every solve in the code
    base bit-reproducible run-to-run, which the G3.3 "M_inf -> 0 is
    bit-identical to Laplace" gate relies on.

    Returns:
        (ml, M): the pyamg hierarchy and its preconditioner LinearOperator
    """
    state = np.random.get_state()
    np.random.seed(0)
    try:
        ml = pyamg.smoothed_aggregation_solver(A)
    finally:
        np.random.set_state(state)
    return ml, ml.aspreconditioner()


def build_ilu_preconditioner(
    A: sp.spmatrix,
    drop_tol: float = 1e-4,
    fill_factor: float = 10.0,
) -> spla.LinearOperator:
    """Incomplete-LU preconditioner for the Newton Jacobian (fallback when
    the Term-1 AMG hierarchy stalls on a strong-shock indefinite block).
    Note it factors the MATRIX it is given; when the GMRES operator is the
    low-rank-corrected J_ff + B K, preconditioning with ILU(J_ff) alone is
    fine -- GMRES absorbs the rank-n_stations perturbation.

    Robustness ladder (B11): the default drop can leave the incomplete factor
    exactly singular on the level-set fused matrix (its aux weld / wake-LS /
    TE-Kutta rows carry zero or negative diagonals -- measured on the M6-family
    2.5D MEDIUM and 3D MEDIUM lifting matrices). On a `spilu` RuntimeError,
    retry with a gentler drop, then MODIFIED ILU (`ILU_MILU="SMILU_2"`, which
    compensates dropped mass onto the diagonal), then a small DIAGONAL-SHIFTED
    MILU (`spilu(A + eps·I)` -- a Manteuffel shift can never produce a zero
    pivot; GMRES corrects the O(eps) perturbation). The default (first) attempt
    is unchanged, so a matrix that factors plainly (conforming Newton Jacobian,
    LS coarse) pays nothing, and the fallbacks stay at modest fill (not the
    huge factor a large fill_factor would build)."""
    Acsc = A.tocsc()
    n = A.shape[0]
    milu = {"ILU_MILU": "SMILU_2"}
    # (shift, spilu_kwargs); the shift for the last attempt is built lazily
    # only if the earlier attempts all fail (so a plainly-factorable matrix
    # pays nothing).
    attempts = (
        (0.0, {"drop_tol": drop_tol, "fill_factor": fill_factor}),
        (0.0, {"drop_tol": drop_tol * 1e-2, "fill_factor": fill_factor * 2.0}),
        (0.0, {"drop_tol": 0.0, "fill_factor": fill_factor * 2.0, "options": milu}),
        ("shift", {"drop_tol": 0.0, "fill_factor": fill_factor * 2.0,
                   "options": milu}),
    )
    last = None
    for shift, kw in attempts:
        try:
            if shift == "shift":  # Manteuffel shift: A + eps*I near-diagonal
                eps = 1e-3 * float(np.abs(Acsc.diagonal()).mean() or 1.0)
                mat = (Acsc + eps * sp.eye(n, format="csc")).tocsc()
            else:
                mat = Acsc
            ilu = spla.spilu(mat, **kw)
            return spla.LinearOperator((n, n), matvec=ilu.solve)
        except RuntimeError as exc:  # "Factor is exactly singular"
            last = exc
    raise last


def solve_gmres(
    A,
    b: np.ndarray,
    M=None,
    rtol: float = 1e-6,
    atol: float = 0.0,
    x0: np.ndarray = None,
    restart: int = 60,
    maxiter: int = 50,
    on_fail: str = "retry",
) -> Tuple[np.ndarray, int, int]:
    """
    Solve A x = b with restarted, preconditioned GMRES (the P8 Newton
    linear solve; design.md Sec 8.1). A may be any scipy sparse matrix or
    LinearOperator (the coupled Newton step uses the eliminated low-rank
    operator x -> J_ff x + B (K x)).

    Args:
        A: (n, n) matrix or LinearOperator, need not be symmetric
        b: (n,) right-hand side
        M: preconditioner LinearOperator (approximate inverse of A), e.g.
            build_amg_preconditioner(A_picard_free)[1] on the SPD Term-1
            block, or build_ilu_preconditioner
        rtol: relative residual tolerance -- the Newton driver passes the
            Eisenstat-Walker forcing eta_k here
        atol: absolute tolerance floor (scipy convention)
        x0: initial guess
        restart: Krylov subspace size between restarts
        maxiter: maximum number of restart cycles
        on_fail: "retry" (default) -- one automatic retry at doubled
            restart, then raise RuntimeError; "return" -- hand back the
            best iterate immediately (an INEXACT Newton step is still a
            descent-quality step; the P8 driver prefers taking it over
            burning thousands of stagnating Krylov iterations -- measured
            2.4 h on the first medium G4.1 run, dominated by retries)

    Returns:
        (x, n_inner_iterations, info) -- info 0 = converged (scipy
        convention); nonzero only reachable with on_fail="return".

    Raises:
        RuntimeError if GMRES fails to converge under on_fail="retry".
    """
    n_iter = [0]

    def _count(_pr_norm):
        n_iter[0] += 1

    x, info = spla.gmres(
        A, b, x0=x0, rtol=rtol, atol=atol, restart=restart, maxiter=maxiter,
        M=M, callback=_count, callback_type="pr_norm",
    )
    if info != 0 and on_fail == "retry":
        n_first = n_iter[0]
        n_iter[0] = 0
        x, info = spla.gmres(
            A, b, x0=x, rtol=rtol, atol=atol, restart=2 * restart,
            maxiter=maxiter, M=M, callback=_count,
            callback_type="pr_norm",
        )
        if info != 0:
            raise RuntimeError(
                f"GMRES did not converge (info={info}) after "
                f"{n_first} + {n_iter[0]} inner iterations "
                f"(restart {restart} then {2 * restart}, rtol={rtol:g})"
            )
        n_iter[0] += n_first
    return x, n_iter[0], info


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
    _, M = build_amg_preconditioner(A)

    n_iter = [0]

    def _count(_xk):
        n_iter[0] += 1

    x, info = spla.cg(A, b, M=M, rtol=rtol, maxiter=maxiter, callback=_count)

    if info != 0:
        raise RuntimeError(
            f"CG did not converge (info={info}) after {n_iter[0]} iterations"
        )

    return x, n_iter[0]
