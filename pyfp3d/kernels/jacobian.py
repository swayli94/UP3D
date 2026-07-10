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

P8 adds the exact Newton Jacobian (design.md (6.3)) on top:

    dR_i/dphi_k = sum_e V_e [ rho_tilde_e (grad N_i . grad N_k)      Term 1
                  + (drho_tilde_e/dphi_k) (grad phi_e . grad N_i) ]  Terms 2+3

with the frozen-selection walk-flux sensitivities (s_e, s_u) =
(drho_tilde_e/dq2_e, drho_tilde_e/dq2_u(e)) from
kernels/upwind.py::rho_tilde_sensitivities (P7) and the DOF chain
dq2/dphi_k = (2/u_inf^2) grad phi . grad N_k supplied here. Term 2
(rows = cols = nodes(e)) has exactly the Term-1 footprint, so it shares
the CSR pattern, coloring and elem_to_csr scatter (one fused colored
kernel). Term 3 (rows = nodes(e), cols = nodes(u(e)), graph-distance <= 4)
falls OUTSIDE the element pattern; it is assembled as active-set COO
(s_u == 0 on subsonic/floored/self-upstream elements by construction, so
the extra entries are 16 per SUPERSONIC element only) and added sparsely
-- rebuilt from scratch every Newton step, so upstream-selection churn
between steps can never corrupt a reused pattern (design.md Sec 6.3:
measure before building any wider colored-CSR machinery).
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


@_njit(cache=True, fastmath=True, parallel=True)
def assemble_jacobian_data_colored(
    color_offsets: np.ndarray,
    color_elems: np.ndarray,
    elem_to_csr: np.ndarray,
    B: np.ndarray,
    V: np.ndarray,
    rho_tilde: np.ndarray,
    grad: np.ndarray,
    s_e: np.ndarray,
    inv_u2: float,
    data: np.ndarray,
) -> None:
    """Newton Terms 1+2 (design.md (6.3), local block) scattered into the
    Picard CSR pattern (data must be zeroed by the caller):

        data[(i,k)] += V_e [ rho_tilde_e (B_i . B_k)
                       + (grad phi_e . B_i) s_e (2 inv_u2) (grad phi_e . B_k) ]

    Same coloring/race-freedom/bit-determinism argument as
    assemble_matrix_data_colored: Term 2's footprint is exactly Term 1's."""
    n_colors = len(color_offsets) - 1
    for c in range(n_colors):
        for k in prange(color_offsets[c], color_offsets[c + 1]):
            e = color_elems[k]
            w1 = rho_tilde[e] * V[e]
            g0 = grad[e, 0]
            g1 = grad[e, 1]
            g2 = grad[e, 2]
            w2 = 2.0 * inv_u2 * s_e[e] * V[e]
            for i in range(4):
                bi0 = B[e, i, 0]
                bi1 = B[e, i, 1]
                bi2 = B[e, i, 2]
                gi = g0 * bi0 + g1 * bi1 + g2 * bi2
                for j in range(4):
                    bj0 = B[e, j, 0]
                    bj1 = B[e, j, 1]
                    bj2 = B[e, j, 2]
                    # Term 1 kept as its own expression (same form as
                    # assemble_matrix_data_colored) and Term 2 added only
                    # when active: a fused w1*(..) + w2*(..) lets fastmath
                    # contract to an FMA whose rounding differs from the
                    # Picard kernel, so J would NOT reduce to A bitwise at
                    # s_e == 0 (masked/limited elements).
                    val = w1 * (bi0 * bj0 + bi1 * bj1 + bi2 * bj2)
                    if w2 != 0.0:
                        val += w2 * gi * (g0 * bj0 + g1 * bj1 + g2 * bj2)
                    data[elem_to_csr[e, i, j]] += val


@_njit(cache=True, fastmath=True, parallel=True)
def fill_term3_coo(
    active: np.ndarray,
    elements: np.ndarray,
    B: np.ndarray,
    V: np.ndarray,
    grad: np.ndarray,
    upstream: np.ndarray,
    s_u: np.ndarray,
    inv_u2: float,
    rows: np.ndarray,
    cols: np.ndarray,
    vals: np.ndarray,
) -> None:
    """Newton Term 3 (design.md (6.3) upstream coupling, Lopez B.4) as COO
    triplets: for each active element e (s_u != 0), 16 entries coupling
    rows = nodes(e) to cols = nodes(u(e)):

        val[(i,k)] = V_e (grad phi_e . B_e,i) s_u (2 inv_u2) (grad phi_u . B_u,k)

    Each active index a owns the fixed slot slice [16a, 16a+16) of the
    preallocated triplet arrays -- prange-safe with no atomics; duplicate
    (row, col) pairs across elements are summed by scipy's COO->CSR."""
    n_active = len(active)
    for a in prange(n_active):
        e = active[a]
        u = upstream[e]
        ge0 = grad[e, 0]
        ge1 = grad[e, 1]
        ge2 = grad[e, 2]
        gu0 = grad[u, 0]
        gu1 = grad[u, 1]
        gu2 = grad[u, 2]
        w = 2.0 * inv_u2 * s_u[e] * V[e]
        base = 16 * a
        for i in range(4):
            gi = ge0 * B[e, i, 0] + ge1 * B[e, i, 1] + ge2 * B[e, i, 2]
            ri = elements[e, i]
            wi = w * gi
            for k in range(4):
                idx = base + 4 * i + k
                rows[idx] = ri
                cols[idx] = elements[u, k]
                vals[idx] = wi * (
                    gu0 * B[u, k, 0] + gu1 * B[u, k, 1] + gu2 * B[u, k, 2]
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

    def assemble_newton_jacobian(
        self,
        phi: np.ndarray,
        rho_tilde: np.ndarray,
        s_e: np.ndarray,
        s_u: np.ndarray,
        upstream: np.ndarray,
        u_inf: float = 1.0,
        lim_mask: Optional[np.ndarray] = None,
    ) -> sp.csr_matrix:
        """Exact Newton Jacobian dR/dphi (design.md (6.3)) at frozen
        upstream selection: Term 1 (Picard) + Term 2 (local density
        sensitivity) fused into the shared colored CSR pattern, plus the
        active-set Term 3 (upstream coupling) as COO (module docstring).

        Args:
            phi: (n_nodes,) potential on the cut mesh (grad recomputed here
                so the Jacobian is exactly consistent with the state)
            rho_tilde: (n_tets,) upwinded densities at phi (forward flux)
            s_e, s_u, upstream: frozen-selection sensitivities from
                UpwindOperator.rho_tilde_sensitivities at the SAME state --
                s = drho_tilde/dq2 with q2 NORMALIZED by u_inf^2, so the
                DOF chain applied here is (2/u_inf^2) grad phi . grad N_k
            u_inf: freestream speed (the q2-normalization scale)
            lim_mask: optional (n_tets,) bool, True where the speed limiter
                is INACTIVE (q2l == q2). The limiter is a flat clamp, so
                d(q2_limited)/d(q2) is 0 on limited elements: s_e is masked
                by lim_mask[e] and s_u by lim_mask[upstream[e]]. None means
                no limiting (exact when n_limited == 0 -- all converged
                states). Sensitivity inputs are not mutated.

        Records self.newton_nnz and self.n_term3_active (the N2
        measurement: Jacobian density vs the Picard matrix). Does not touch
        any Picard buffer -- the forward assemble_matrix/assemble_residual
        paths are byte-identical before and after a call."""
        if lim_mask is not None:
            s_e = s_e * lim_mask
            s_u = s_u * lim_mask[upstream]
        inv_u2 = 1.0 / (u_inf * u_inf)

        from pyfp3d.kernels.gradient import element_velocity_q2

        grad = np.empty((self.n_tets, 3), dtype=np.float64)
        q2_scratch = np.empty(self.n_tets, dtype=np.float64)
        element_velocity_q2(self.elements, self.B, phi, grad, q2_scratch)

        data = np.zeros(self._pattern.nnz, dtype=np.float64)
        assemble_jacobian_data_colored(
            self.color_offsets, self.color_elems, self.elem_to_csr,
            self.B, self.V, rho_tilde, grad, s_e, inv_u2, data,
        )
        J = sp.csr_matrix(
            (data, self._pattern.indices, self._pattern.indptr),
            shape=self._pattern.shape,
        )

        active = np.nonzero(s_u)[0].astype(np.int64)
        self.n_term3_active = int(len(active))
        if len(active) > 0:
            n_trip = 16 * len(active)
            rows = np.empty(n_trip, dtype=np.int64)
            cols = np.empty(n_trip, dtype=np.int64)
            vals = np.empty(n_trip, dtype=np.float64)
            fill_term3_coo(
                active, self.elements, self.B, self.V, grad, upstream,
                s_u, inv_u2, rows, cols, vals,
            )
            J = (J + sp.coo_matrix(
                (vals, (rows, cols)), shape=self._pattern.shape,
            ).tocsr()).tocsr()
        self.newton_nnz = int(J.nnz)
        return J

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
