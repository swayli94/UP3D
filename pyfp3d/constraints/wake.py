"""
Wake jump constraint: master-slave elimination of the "+"-side copies and
the Kutta update target (design.md Sec 4, eqs (4.3)-(4.4)).

The cut mesh (mesh/wake_cut.py) appends the "+"-side (slave) copies AFTER
the original nodes, so the reduced dof space is exactly the ORIGINAL node
set: phi_full = T phi_red + g(Gamma), with

    T[i, i]         = 1   for i < n_orig          (identity block)
    T[slave_k, m_k] = 1                            (slave rows)
    g[slave_k]      = Gamma(station of k), else 0

Reduced system: (T^T A T) phi_red = T^T (b - A g). A and T are fixed;
Gamma enters only through g, so a Gamma update is RHS-only -- per station
j the vector h_j = T^T A g_j (g_j = indicator of station-j slaves) is
precomputed once and b_red(Gamma) = T^T b - sum_j Gamma_j h_j.

T has full column rank, so T^T A T stays SPD and CG+AMG applies unchanged.
Folding the slave rows into their masters also sums the sheet's two
one-sided flux integrals into one interior-like equation, which is exactly
the weak mass-flux continuity [rho dphi/dn] = 0 of (4.2).
"""

from typing import Tuple

import numpy as np
import scipy.sparse as sp

from pyfp3d.mesh.wake_cut import WakeCut


class WakeConstraint:
    """Precomputed elimination operator for one cut mesh."""

    def __init__(self, A: sp.spmatrix, wc: WakeCut):
        """
        Args:
            A: (n_cut, n_cut) stiffness matrix assembled on the CUT mesh
            wc: WakeCut from mesh/wake_cut.py
        """
        n_cut = A.shape[0]
        n_red = wc.n_nodes_orig
        n_slaves = len(wc.slave_nodes)
        if n_cut != n_red + n_slaves:
            raise ValueError(
                f"A is {n_cut}x{n_cut} but WakeCut implies "
                f"{n_red} + {n_slaves} nodes"
            )

        rows = np.concatenate([np.arange(n_red, dtype=np.int64), wc.slave_nodes])
        cols = np.concatenate([np.arange(n_red, dtype=np.int64), wc.master_nodes])
        vals = np.ones(n_cut, dtype=np.float64)
        self.T = sp.coo_matrix((vals, (rows, cols)), shape=(n_cut, n_red)).tocsr()

        self.wc = wc
        self.n_reduced = n_red
        self.update_matrix(A)

    def update_matrix(self, A: sp.spmatrix) -> None:
        """Recompute A_reduced = T^T A T and the per-station RHS vectors
        h_j = T^T A g_j for a NEW A on the same cut mesh. T is purely
        topological (mesh-only), so the P3 Picard loop calls this once per
        density update instead of rebuilding the whole constraint."""
        self.A_reduced = (self.T.T @ (A @ self.T)).tocsr()

        n_cut = self.T.shape[0]
        wc = self.wc
        A_csr = A.tocsr()
        self._h = np.empty((wc.n_stations, self.n_reduced), dtype=np.float64)
        for j in range(wc.n_stations):
            g_j = np.zeros(n_cut, dtype=np.float64)
            g_j[wc.slave_nodes[wc.node_station == j]] = 1.0
            self._h[j] = self.T.T @ (A_csr @ g_j)

    def reduced_rhs(self, b: np.ndarray, gamma: np.ndarray) -> np.ndarray:
        """T^T b - sum_j Gamma_j h_j  (Gamma is RHS-only by construction)."""
        gamma = np.atleast_1d(np.asarray(gamma, dtype=np.float64))
        return self.T.T @ b - gamma @ self._h

    def expand(self, phi_red: np.ndarray, gamma: np.ndarray) -> np.ndarray:
        """phi on the cut mesh: masters copied to slaves plus the jump."""
        gamma = np.atleast_1d(np.asarray(gamma, dtype=np.float64))
        phi = np.empty(self.n_reduced + len(self.wc.slave_nodes), dtype=np.float64)
        phi[: self.n_reduced] = phi_red
        phi[self.wc.slave_nodes] = (
            phi_red[self.wc.master_nodes] + self.wc.gamma_at_nodes(gamma)
        )
        return phi

    def to_reduced_dirichlet(
        self, dirichlet_nodes: np.ndarray, dirichlet_values: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Map cut-mesh Dirichlet data to the reduced space.

        Slave entries are dropped: an eliminated "+"-side far-field node
        takes the value master + Gamma automatically, which matches the
        vortex correction's branch-cut jump across the wake (see
        constraints/dirichlet.py -- the cut lies ON the wake sheet).
        """
        nodes = np.asarray(dirichlet_nodes, dtype=np.int64)
        keep = nodes < self.n_reduced
        red_nodes, idx = np.unique(nodes[keep], return_index=True)
        return red_nodes, np.asarray(dirichlet_values, dtype=np.float64)[keep][idx]


def kutta_targets(phi: np.ndarray, wc: WakeCut) -> np.ndarray:
    """Per-station Kutta target (design.md (4.4)), one node off the TE:

        Gamma_j = mean over the station's TE nodes of
                  (phi[upper wall probe] - phi[lower wall probe])

    Probes are wall nodes (free dofs), NOT wake nodes -- at wake nodes the
    jump equals the prescribed Gamma identically by (4.3), so measuring
    there would be circular. On a quasi-2D extrusion one station spans
    both z-planes, so the mean also filters the O(h) spanwise noise of the
    3-tet prism split (roadmap G2.5 evidence note).
    """
    jumps = phi[wc.kutta_upper] - phi[wc.kutta_lower]
    n_st = wc.n_stations
    sums = np.bincount(wc.te_station, weights=jumps, minlength=n_st)
    counts = np.bincount(wc.te_station, minlength=n_st)
    return sums / counts
