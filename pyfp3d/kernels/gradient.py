"""
Element-constant P1 velocity sweep: grad(phi)_e and q^2_e for all elements
from the precomputed shape-gradient matrices B_e.

This is the first sweep of every Picard iteration (design.md Sec 8) and is
embarrassingly parallel (per-element writes, no scatter), so it uses a plain
`prange` with no coloring. Zero allocation: callers pass preallocated output
arrays (design.md Sec 7 rule 4).

Reference: design.md Sec 6 ("one 4x3 gemv for grad(phi)_e"), Sec 7.
"""

import os

import numba
import numpy as np

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


@_njit(cache=True, fastmath=True, parallel=True)
def element_velocity_q2(
    elements: np.ndarray,
    B: np.ndarray,
    phi: np.ndarray,
    grad_out: np.ndarray,
    q2_out: np.ndarray,
) -> None:
    """
    grad(phi)_e = sum_i phi_i B_e[i, :] and q^2_e = |grad(phi)_e|^2 per tet.

    Args:
        elements: (n_tets, 4) connectivity
        B: (n_tets, 4, 3) precomputed basis gradients
            (mesh.metrics.precompute_element_geometry)
        phi: (n_nodes,) nodal potential
        grad_out: (n_tets, 3) preallocated output, overwritten
        q2_out: (n_tets,) preallocated output, overwritten
    """
    n_tets = len(elements)
    for e in prange(n_tets):
        gx = 0.0
        gy = 0.0
        gz = 0.0
        for i in range(4):
            p = phi[elements[e, i]]
            gx += p * B[e, i, 0]
            gy += p * B[e, i, 1]
            gz += p * B[e, i, 2]
        grad_out[e, 0] = gx
        grad_out[e, 1] = gy
        grad_out[e, 2] = gz
        q2_out[e] = gx * gx + gy * gy + gz * gz
