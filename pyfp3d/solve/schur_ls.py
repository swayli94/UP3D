"""B14: Schur-eliminated aux block + AMG(SPD Picard main block) for the fused
level-set system (design_track_b.md §5.3; roadmap track_b.md §B14).

The fused wake_ls matrix is structurally nonsymmetric only because of the aux
rows (wake-LS g1+g2 + nonlinear TE-Kutta, negative diagonals) -- pyamg cannot
precondition it directly, and the B11 SPD spring surrogate STALLS on the
lifting operator because its jump==0 prior kills the global circulation mode
(gamma 0.0033 vs 0.139). B14 removes the mismatch STRUCTURALLY: eliminate the
SMALL aux block exactly,

    K x_m = J_mm x_m - J_ma . J_aa^{-1} . (J_am x_m)        (matrix-free)
    r     = b_m - J_ma . J_aa^{-1} b_a
    x_a   = J_aa^{-1} (b_a - J_am x_m)                      (back-substitution)

with `lu_aa = splu(J_aa)` a thin-strip factorization (n_ext ~8k at M6 medium,
milliseconds; the TE-Kutta rows re-linearize every Newton step / Picard outer,
so the strip refactors per construction). GMRES on K then faces "elliptic +
cut-strip-localized correction" -- the operator shape the conforming path
already preconditions with AMG on the SPD Picard Term-1 block
(solve/newton.py) -- so the preconditioner here is the exact analogue, built
on `op.assemble_matrix(rho_own)` restricted to the main-free set. NO springs:
no aux DOF survives into the preconditioned system.

Direction is load-bearing (roadmap B14 corrections): eliminate the SMALL aux
block. Eliminating main instead would need A_mm^{-1} = an AMG inner solve per
operator application.

Honesty identity (locked by tests/test_b14_schur_ls.py): the returned
x = [x_m; x_a] satisfies the aux rows EXACTLY (to LU roundoff) by
back-substitution, and the full-system main-row residual equals the reduced
GMRES residual r - K x_m -- so the reduced rtol IS the full-system relative
residual, and stall reporting needs no extra bookkeeping.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from pyfp3d.solve.linear import build_amg_preconditioner, solve_gmres


class SchurReducedSystem:
    """Exact elimination of the aux tail block of a free-reduced fused matrix.

    Constructed FRESH per Newton step / Picard outer (the TE-Kutta rows change
    per linearization). Relies on the B3 load-bearing fact that aux DOFs
    (global ids >= n_main) are NEVER Dirichlet, and that `free` is sorted
    ascending -- so within A_free = A[free][:, free] the aux DOFs are the
    contiguous tail and no permutation is needed.

    Args:
        A_free: the free-reduced fused matrix (main-free + aux, square).
        free: sorted free-DOF indices into the n_total numbering.
        n_main: mvop.n_main (aux DOFs are the ids >= n_main).
        n_aux_expected: mvop.n_ext; a mismatch means an aux DOF became
            Dirichlet somewhere upstream -- fail loudly, never silently
            mis-split.
    """

    def __init__(self, A_free: sp.spmatrix, free: np.ndarray, n_main: int,
                 n_aux_expected: Optional[int] = None) -> None:
        free = np.asarray(free, dtype=np.int64)
        n_mf = int(np.searchsorted(free, n_main))
        n_aux = int(free.size - n_mf)
        if n_aux_expected is not None and n_aux != int(n_aux_expected):
            raise ValueError(
                f"Schur split expected {n_aux_expected} aux DOFs in the free "
                f"set but found {n_aux}: an aux DOF appears in the Dirichlet "
                "set (aux is never Dirichlet -- B3) or `free` is unsorted.")
        A = A_free.tocsr()
        self.n_mf = n_mf
        self.n_aux = n_aux
        self.main_free = free[:n_mf]          # global node ids, all < n_main
        self.J_mm = A[:n_mf, :n_mf].tocsr()
        self.J_ma = A[:n_mf, n_mf:].tocsr()
        self.J_am = A[n_mf:, :n_mf].tocsr()
        self.J_aa = A[n_mf:, n_mf:].tocsc()
        # The thin-strip factorization; raises on a singular strip (GB14.1
        # measures conditioning -- "measure, don't assume").
        self.lu_aa = spla.splu(self.J_aa)

    def operator(self) -> spla.LinearOperator:
        """The reduced operator K, matrix-free (one thin solve per matvec)."""

        def matvec(x: np.ndarray) -> np.ndarray:
            return self.J_mm @ x - self.J_ma @ self.lu_aa.solve(self.J_am @ x)

        return spla.LinearOperator((self.n_mf, self.n_mf), matvec=matvec,
                                   dtype=np.float64)

    def reduce_rhs(self, b_free: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Split b_free and form the reduced RHS r = b_m - J_ma J_aa^-1 b_a."""
        b_m = b_free[:self.n_mf]
        b_a = b_free[self.n_mf:]
        return b_m - self.J_ma @ self.lu_aa.solve(b_a), b_a

    def back_substitute(self, x_m: np.ndarray, b_a: np.ndarray) -> np.ndarray:
        """Aux solution from the main one: exact aux-row satisfaction."""
        return self.lu_aa.solve(b_a - self.J_am @ x_m)

    def solve(self, b_free: np.ndarray, M: Optional[spla.LinearOperator],
              rtol: float, x0_main: Optional[np.ndarray] = None,
              restart: int = 60, maxiter: int = 50,
              ) -> Tuple[np.ndarray, int, int]:
        """Solve A_free x = b_free via the reduced system.

        Returns (x_full, n_inner_iterations, info) with the solve_gmres
        conventions (info != 0 = stalled; the best iterate is still returned,
        with the aux part back-substituted so the aux rows hold exactly).
        """
        r, b_a = self.reduce_rhs(b_free)
        x_m, n_it, info = solve_gmres(
            self.operator(), r, M=M, rtol=rtol, x0=x0_main,
            restart=restart, maxiter=maxiter, on_fail="return")
        x_m = np.atleast_1d(x_m)
        x_a = self.back_substitute(x_m, b_a)
        return np.concatenate([x_m, x_a]), n_it, info


def main_block_preconditioner(mvop, rho_pair, main_free: np.ndarray
                              ) -> spla.LinearOperator:
    """AMG on the SPD single-valued Picard Term-1 block, main-free restricted.

    The exact conforming analogue (solve/newton.py builds AMG on
    op.assemble_matrix(rho_t) reduced to free DOFs). NO springs -- the whole
    point of B14 vs the B11 `_amg_surrogate_preconditioner`.

    Args:
        mvop: the MultivaluedOperator (supplies op + own_side_field).
        rho_pair: (rho_up, rho_lo) per-element densities, or None (Laplace
            limit / subcritical rho == 1).
        main_free: global node ids of the free MAIN DOFs (all < n_main),
            i.e. SchurReducedSystem.main_free.
    """
    rho_own = None if rho_pair is None else mvop.own_side_field(*rho_pair)
    A_pic = mvop.op.assemble_matrix(rho_own).tocsr()
    return build_amg_preconditioner(A_pic[main_free][:, main_free].tocsr())[1]


def jaa_diagnostic(schur: SchurReducedSystem) -> Dict[str, float]:
    """GB14.1 pre-registered diagnostic: measured J_aa conditioning.

    1-norm condition estimate cond1 = onenormest(J_aa) * onenormest(J_aa^-1),
    the inverse applied through the already-computed lu_aa. Diagnostic only
    (demo/tests) -- never in the production solve loop.
    """
    n = schur.n_aux
    norm1 = float(spla.onenormest(schur.J_aa))
    inv_op = spla.LinearOperator(
        (n, n),
        matvec=schur.lu_aa.solve,
        rmatvec=lambda b: schur.lu_aa.solve(b, trans="T"),
        dtype=np.float64)
    inv_norm1 = float(spla.onenormest(inv_op))
    return {"n_aux": float(n), "norm1": norm1, "inv_norm1": inv_norm1,
            "cond1": norm1 * inv_norm1}
