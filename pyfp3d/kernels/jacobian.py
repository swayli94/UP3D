"""
Picard (frozen-density) matrix assembly, A_ij = sum_e rho_tilde_e
(grad N_i . grad N_j) V_e (design.md (6.2)), on a precomputed symbolic
sparsity pattern with colored-`prange` writes.

Retires the P1 assembly tech debt (roadmap P3): per-element geometry
(B_e, V_e) is precomputed once, the CSR sparsity and the elem_to_csr
scatter map are built once, and the per-iteration hot kernel allocates
nothing (design.md Sec 7 rules 2 and 4; agent-rules hard rule #3). The
coloring guarantees no two same-color elements share a node, hence no two
same-color elements touch the same CSR entry -- prange within a color is
race-free, and the per-nnz accumulation order is fixed by the color
sequence alone, so assembly is bit-deterministic across runs and thread
counts.

The exact Newton Jacobian (design.md (6.3)) lands here in P6.
"""

import os
from typing import Optional

import numba
import numpy as np
import scipy.sparse as sp

from pyfp3d.mesh.coloring import color_partition_csr
from pyfp3d.mesh.metrics import precompute_element_geometry

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


def build_csr_pattern(elements: np.ndarray, n_nodes: int) -> sp.csr_matrix:
    """Symbolic CSR sparsity of the P1 Galerkin operator (node-to-node
    coupling through shared elements), with canonical (sorted, deduped)
    indices and zeroed data. Built once per mesh."""
    el = np.asarray(elements, dtype=np.int64)
    rows = np.repeat(el, 4, axis=1).reshape(-1)
    cols = np.tile(el, (1, 4)).reshape(-1)
    pattern = sp.coo_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, cols)),
        shape=(n_nodes, n_nodes),
    ).tocsr()
    pattern.sum_duplicates()
    pattern.sort_indices()
    pattern.data[:] = 0.0
    return pattern


@_njit(cache=True)
def build_elem_to_csr(
    elements: np.ndarray, indptr: np.ndarray, indices: np.ndarray
) -> np.ndarray:
    """elem_to_csr[e, i, j] -> nnz index of (row=elements[e,i],
    col=elements[e,j]) in the canonical CSR pattern (design.md Sec 7 rule 2:
    precomputed once, assembly then writes straight into A.data)."""
    n_tets = len(elements)
    emap = np.empty((n_tets, 4, 4), dtype=np.int64)
    for e in range(n_tets):
        for i in range(4):
            row = elements[e, i]
            lo = indptr[row]
            hi = indptr[row + 1]
            for j in range(4):
                col = elements[e, j]
                a = lo
                b = hi
                while a < b:
                    mid = (a + b) // 2
                    if indices[mid] < col:
                        a = mid + 1
                    else:
                        b = mid
                emap[e, i, j] = a
    return emap


@_njit(cache=True, fastmath=True, parallel=True)
def assemble_matrix_data_colored(
    color_offsets: np.ndarray,
    color_elems: np.ndarray,
    elem_to_csr: np.ndarray,
    B: np.ndarray,
    V: np.ndarray,
    rho_tilde: np.ndarray,
    data: np.ndarray,
) -> None:
    """Scatter rho_tilde_e V_e (B_e B_e^T) into CSR data (must be zeroed by
    the caller). Serial outer loop over colors, prange within a color."""
    n_colors = len(color_offsets) - 1
    for c in range(n_colors):
        for k in prange(color_offsets[c], color_offsets[c + 1]):
            e = color_elems[k]
            w = rho_tilde[e] * V[e]
            for i in range(4):
                bi0 = B[e, i, 0]
                bi1 = B[e, i, 1]
                bi2 = B[e, i, 2]
                for j in range(4):
                    data[elem_to_csr[e, i, j]] += w * (
                        bi0 * B[e, j, 0] + bi1 * B[e, j, 1] + bi2 * B[e, j, 2]
                    )


class PicardOperator:
    """Per-mesh assembly workspace: everything that design.md Sec 7/Sec 8
    says to compute ONCE (B_e, V_e, element coloring, CSR sparsity,
    elem_to_csr scatter map, preallocated sweep buffers), plus the
    per-iteration assemble/residual entry points that allocate nothing but
    the returned csr_matrix wrapper.

    Usage (Picard outer loop, design.md Sec 8):
        op = PicardOperator(nodes, elements)
        grad, q2 = op.velocities(phi)          # element sweep
        rho = density_field(q2 / u_inf**2, m_inf)   # physics/isentropic.py
        A = op.assemble_matrix(rho)            # colored assembly, (6.2)
        R = op.assemble_residual(phi, rho)     # colored assembly, (6.1)
    """

    def __init__(self, nodes: np.ndarray, elements: np.ndarray):
        self.n_nodes = len(nodes)
        self.n_tets = len(elements)
        self.elements = np.ascontiguousarray(elements)
        self.B, self.V = precompute_element_geometry(nodes, self.elements)
        self.color_offsets, self.color_elems, self.n_colors = color_partition_csr(
            self.elements
        )
        self._pattern = build_csr_pattern(self.elements, self.n_nodes)
        self.elem_to_csr = build_elem_to_csr(
            self.elements, self._pattern.indptr, self._pattern.indices
        )
        # Preallocated hot-loop buffers (design.md Sec 7 rule 4).
        self._grad = np.empty((self.n_tets, 3), dtype=np.float64)
        self._q2 = np.empty(self.n_tets, dtype=np.float64)
        self._R = np.empty(self.n_nodes, dtype=np.float64)
        self.rho_ones = np.ones(self.n_tets, dtype=np.float64)

    def velocities(self, phi: np.ndarray):
        """Element-constant grad(phi)_e and q^2_e (views into the workspace
        buffers -- consume before the next call)."""
        from pyfp3d.kernels.gradient import element_velocity_q2

        element_velocity_q2(self.elements, self.B, phi, self._grad, self._q2)
        return self._grad, self._q2

    def assemble_matrix(self, rho_tilde: Optional[np.ndarray] = None) -> sp.csr_matrix:
        """Weighted stiffness matrix (6.2); rho_tilde=None means the Laplace
        limit rho == 1. Returns a fresh csr_matrix sharing the pattern's
        index arrays but owning its data."""
        if rho_tilde is None:
            rho_tilde = self.rho_ones
        data = np.zeros(self._pattern.nnz, dtype=np.float64)
        assemble_matrix_data_colored(
            self.color_offsets, self.color_elems, self.elem_to_csr,
            self.B, self.V, rho_tilde, data,
        )
        return sp.csr_matrix(
            (data, self._pattern.indices, self._pattern.indptr),
            shape=self._pattern.shape,
        )

    def assemble_residual(
        self, phi: np.ndarray, rho_tilde: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """Nonlinear residual (6.1) at (phi, rho_tilde); returns a view into
        the workspace buffer -- consume before the next call."""
        from pyfp3d.kernels.residual import assemble_residual_colored

        if rho_tilde is None:
            rho_tilde = self.rho_ones
        grad, _ = self.velocities(phi)
        assemble_residual_colored(
            self.color_offsets, self.color_elems, self.elements,
            self.B, self.V, rho_tilde, grad, self._R,
        )
        return self._R
