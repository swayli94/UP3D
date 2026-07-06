"""
Galerkin P1 residual assembly (design.md (6.1)) and the stiffness-matrix
entry point.

Since P3 the hot paths are the colored-`prange` kernels driven by
kernels/jacobian.py::PicardOperator (precomputed B_e/V_e, no hot-loop
allocation -- design.md Sec 7 rules 2/4, agent-rules hard rule #3). The
original serial P1 kernels are kept below as the *reference
implementation*: they recompute geometry per call and stay deliberately
simple, and the P3 regression suite asserts the fast path reproduces them
(tests/test_p3_assembly.py).

Reference: design.md Sec 6 (Spatial discretization), Sec 7 (Numba kernel
architecture).
"""

import os

import numba
import numpy as np
import scipy.sparse as sp

from pyfp3d.mesh.metrics import compute_tet_volumes, element_gradients

if os.environ.get("PYFP3D_NOJIT", "0") == "1":
    prange = range

    def _njit(*args, **kwargs):
        def _identity(func):
            return func
        return _identity
else:
    from numba import prange

    def _njit(*args, **kwargs):
        return numba.njit(*args, **kwargs)


# ---------------------------------------------------------------------------
# Fast colored kernels (P3+): consume precomputed B_e/V_e and a coloring.
# ---------------------------------------------------------------------------

@_njit(cache=True, fastmath=True, parallel=True)
def assemble_residual_colored(
    color_offsets: np.ndarray,
    color_elems: np.ndarray,
    elements: np.ndarray,
    B: np.ndarray,
    V: np.ndarray,
    rho_tilde: np.ndarray,
    grad_elem: np.ndarray,
    R: np.ndarray,
) -> None:
    """
    R_i = sum_e rho_tilde_e (grad(phi)_e . grad(N_i)) V_e  (design.md (6.1)).

    Serial outer loop over colors, prange within a color: no two same-color
    elements share a node, so the nodal scatter-add is race-free AND the
    per-node accumulation order is fixed by the color sequence (bit-
    deterministic across runs/thread counts). Zero allocation: R is
    overwritten in place, grad_elem comes from the gradient sweep
    (kernels/gradient.py).
    """
    R[:] = 0.0
    n_colors = len(color_offsets) - 1
    for c in range(n_colors):
        for k in prange(color_offsets[c], color_offsets[c + 1]):
            e = color_elems[k]
            w = rho_tilde[e] * V[e]
            gx = grad_elem[e, 0]
            gy = grad_elem[e, 1]
            gz = grad_elem[e, 2]
            for i in range(4):
                R[elements[e, i]] += w * (
                    gx * B[e, i, 0] + gy * B[e, i, 1] + gz * B[e, i, 2]
                )


def assemble_stiffness_matrix(nodes: np.ndarray, elements: np.ndarray) -> sp.csr_matrix:
    """
    Assemble the SPD Laplace stiffness matrix A_ij (design.md (6.2), rho == 1).

    Since P3 this delegates to the colored fast path (a one-shot
    PicardOperator workspace), so every consumer -- including the P1/P2
    Laplace drivers -- runs the SAME assembly code as the compressible
    Picard loop: the G3.3 "nu == 0 / M_inf -> 0 is bit-identical to
    Laplace" guarantee is structural, not tested-in.

    Args:
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) tetrahedral connectivity

    Returns:
        A: (n_nodes, n_nodes) sparse CSR matrix
    """
    from pyfp3d.kernels.jacobian import PicardOperator

    return PicardOperator(nodes, elements).assemble_matrix()


# ---------------------------------------------------------------------------
# Reference implementations (P1): serial, geometry recomputed per call.
# Kept for the old-vs-new assembly regression and for NOJIT debugging.
# ---------------------------------------------------------------------------

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


def assemble_stiffness_matrix_reference(
    nodes: np.ndarray, elements: np.ndarray
) -> sp.csr_matrix:
    """The original P1 serial triplet assembly (rho == 1), kept as the
    independent reference the P3 fast path is regression-tested against."""
    n_nodes = len(nodes)
    volumes = compute_tet_volumes(nodes, elements)
    rows, cols, vals = _assemble_stiffness_triplets(nodes, elements, volumes)
    A = sp.coo_matrix((vals, (rows, cols)), shape=(n_nodes, n_nodes)).tocsr()
    return A


@_njit(cache=True, fastmath=True)
def assemble_residual(nodes: np.ndarray, elements: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """
    Assemble the Laplace residual R_i (design.md (6.1), rho_tilde == 1) --
    serial reference implementation, geometry recomputed per call.

    For phi == x (uniform freestream), R must be machine-zero on any mesh
    (design.md Sec 3, "freestream preservation") -- this is the first
    regression check after any kernel change (agent-rules.md hard rule #1).

    Hot loops (the Picard iteration) use assemble_residual_colored via
    PicardOperator instead; this one-shot form remains for the P1 drivers
    and the freestream regression.

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
