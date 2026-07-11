"""
Level-set (Track B) solve drivers -- parallel to solve/picard.py, never
imported by the conforming path.

B2 ships the NON-LIFTING driver only: a single extended assemble + direct
solve on the multivalued operator, enough for the B2 consistency gates
(V0 freestream, V1 MMS convergence, Laplace a=0 -> cl ~ 0). The extended
matrix is nonsymmetric (the aux-row weld couples aux -> main one-way), so
it is solved with a sparse direct LU (scipy `spsolve`), not CG/AMG -- for
the coarse/medium B2 meshes the direct factor is cheap and isolates
"is the assembly correct" from any preconditioner-convergence question.
GMRES + AMG(main block) (solve/linear.py) is the scaling path B3+ adopts
(design_track_b.md section 5.3).

The lifting driver with implicit Kutta (the g1+g2 wake LS closure, no Gamma
outer loop) is B3.
"""

from typing import Dict, Optional

import numpy as np
import scipy.sparse.linalg as spla

from pyfp3d.solve.linear import apply_dirichlet
from pyfp3d.wake.multivalued import MultivaluedOperator


def solve_multivalued_laplace(
    mvop: MultivaluedOperator,
    dirichlet_nodes: np.ndarray,
    dirichlet_values: np.ndarray,
    body_source_rhs: Optional[np.ndarray] = None,
) -> Dict[str, object]:
    """Non-lifting Laplace solve on a level-set cut mesh (B2).

    Assembles the extended multivalued matrix (rho == 1) with the B2
    continuity closure and solves A x = b with Dirichlet BCs on the main
    DOFs (far field). With the weld closure the solution is single valued,
    so this reproduces the ordinary single-valued Laplace/MMS solution --
    the B2 consistency check.

    Args:
        mvop: MultivaluedOperator for the mesh + wake level set
        dirichlet_nodes: node ids with prescribed phi (< n_main); their aux
            DOFs are pinned to the same value through the weld rows
        dirichlet_values: prescribed phi at dirichlet_nodes
        body_source_rhs: optional (n_main,) consistent load vector on the
            main DOFs (an MMS source); aux rows get 0. The physical FP
            equation has no volume source.

    Returns:
        dict: phi_ext (n_total), phi (n_main main potential), te_jump
        (Gamma at TE nodes), residual_norm (inf-norm of the reduced free-DOF
        residual), n_total, n_ext.
    """
    A = mvop.assemble_matrix(rho_tilde=None, closure="continuity")
    b = np.zeros(mvop.n_total, dtype=np.float64)
    if body_source_rhs is not None:
        b[: mvop.n_main] = np.asarray(body_source_rhs, dtype=np.float64)

    dirichlet_nodes = np.asarray(dirichlet_nodes, dtype=np.int64)
    dirichlet_values = np.asarray(dirichlet_values, dtype=np.float64)
    A_free, b_free, free, phi_ext = apply_dirichlet(
        A, b, dirichlet_nodes, dirichlet_values
    )

    x = spla.spsolve(A_free.tocsc(), b_free)
    x = np.atleast_1d(x)
    phi_ext[free] = x
    residual = float(np.max(np.abs(b_free - A_free @ x))) if len(x) else 0.0

    return {
        "phi_ext": phi_ext,
        "phi": mvop.main_potential(phi_ext),
        "te_jump": mvop.te_jump(phi_ext),
        "residual_norm": residual,
        "n_total": mvop.n_total,
        "n_ext": mvop.n_ext,
    }
