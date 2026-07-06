"""
Laplace element assembly (P1: rho == 1, no upwinding, no wake).

Galerkin weak form (design.md (6.1)-(6.2)) with rho_tilde == 1:

    R_i = sum_e (grad(phi)_e . grad(N_i)) V_e
    A_ij = sum_e (grad(N_i) . grad(N_j)) V_e

Reference: design.md Sec 6 (Spatial discretization), Sec 7 (Numba kernel
architecture).
"""

import os

import numba
import numpy as np
import scipy.sparse as sp

from pyfp3d.mesh.metrics import compute_tet_volumes, element_gradients


def _njit(*args, **kwargs):
    if os.environ.get("PYFP3D_NOJIT", "0") == "1":
        def _identity(func):
            return func
        return _identity
    return numba.njit(*args, **kwargs)


@_njit(cache=True, fastmath=True)
def _assemble_stiffness_triplets(nodes, elements, volumes):
    n_tets = len(elements)
    rows = np.empty(n_tets * 16, dtype=np.int64)
    cols = np.empty(n_tets * 16, dtype=np.int64)
    vals = np.empty(n_tets * 16, dtype=np.float64)

    idx = 0
    for e in range(n_tets):
        tet = elements[e]
        grads = element_gradients(nodes, elements, e)
        V = volumes[e]
        for i in range(4):
            for j in range(4):
                rows[idx] = tet[i]
                cols[idx] = tet[j]
                vals[idx] = V * (
                    grads[i, 0] * grads[j, 0]
                    + grads[i, 1] * grads[j, 1]
                    + grads[i, 2] * grads[j, 2]
                )
                idx += 1

    return rows, cols, vals


def assemble_stiffness_matrix(nodes: np.ndarray, elements: np.ndarray) -> sp.csr_matrix:
    """
    Assemble the SPD Laplace stiffness matrix A_ij (design.md (6.2), rho == 1).

    Args:
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) tetrahedral connectivity

    Returns:
        A: (n_nodes, n_nodes) sparse CSR matrix
    """
    n_nodes = len(nodes)
    volumes = compute_tet_volumes(nodes, elements)
    rows, cols, vals = _assemble_stiffness_triplets(nodes, elements, volumes)
    A = sp.coo_matrix((vals, (rows, cols)), shape=(n_nodes, n_nodes)).tocsr()
    return A


@_njit(cache=True, fastmath=True)
def assemble_residual(nodes: np.ndarray, elements: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """
    Assemble the Laplace residual R_i (design.md (6.1), rho_tilde == 1).

    For phi == x (uniform freestream), R must be machine-zero on any mesh
    (design.md Sec 3, "freestream preservation") -- this is the first
    regression check after any kernel change (agent-rules.md hard rule #1).

    Args:
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) tetrahedral connectivity
        phi: (n_nodes,) nodal potential

    Returns:
        R: (n_nodes,) assembled residual
    """
    n_nodes = len(nodes)
    n_tets = len(elements)
    R = np.zeros(n_nodes, dtype=np.float64)
    volumes = compute_tet_volumes(nodes, elements)

    for e in range(n_tets):
        tet = elements[e]
        grads = element_gradients(nodes, elements, e)
        V = volumes[e]

        grad_phi = np.zeros(3, dtype=np.float64)
        for i in range(4):
            for d in range(3):
                grad_phi[d] += phi[tet[i]] * grads[i, d]

        for i in range(4):
            dotp = (
                grad_phi[0] * grads[i, 0]
                + grad_phi[1] * grads[i, 1]
                + grad_phi[2] * grads[i, 2]
            )
            R[tet[i]] += V * dotp

    return R
