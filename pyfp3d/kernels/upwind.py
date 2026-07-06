"""
Artificial-density upwinding for supersonic zones (design.md Sec 3,
roadmap P4): per element,

    rho_tilde_e = rho_e - nu_e (rho_e - rho_u(e))                   (3.1')
    nu_e = C * max(0, 1 - M_c^2 / M_e^2)                            (3.2)

with u(e) the upstream element found by a MULTI-HOP walk through face
neighbors (an evidence-based extension of the design.md Sec 3 step-1
single-neighbor rule): starting at e, repeatedly step to the face
neighbor with the most negative centroid displacement along V_e, until
the ACCUMULATED streamwise displacement reaches the element's own
streamwise extent (capped at `max_hops`). Rationale, measured on the M0
prism-split quasi-2D meshes: a single face-neighbor hop reaches only
0.25-0.39 of the element's streamwise extent in the supersonic pocket
(stacked sliver tets), so the one-hop rho_e - rho_u bias delivers ~1/3
of the intended Delta-l * d(rho)/dl and the effective nu falls BELOW the
(M^2-1)/M^2 stability threshold for M >~ 1.2 -- the observed consequence
was a violent instability growing out of the pocket peak that no
omega/C tuning could tame (blow-up at M_inf = 0.75 for every omega in
0.3-0.9, C in 1.5-2.0), while the marginally-stable M_inf = 0.70 case
crawled to convergence. The walk restores a genuine one-cell upstream
reach. Face adjacency and centroids are precomputed once; the
per-iteration sweeps are prange over elements with zero allocation
(Sec 7 rules 3/4).

Exactness of the subcritical no-op (gate G4.2): where M_e <= M_c the
switch gives nu_e = 0.0 EXACTLY (max against 0), and
rho_e - 0.0 * (rho_e - rho_u) == rho_e bitwise, so a subcritical solve
with the upwind machinery in the loop is bit-identical to the P3 path.
Freestream preservation (design.md Sec 3 properties): phi = x gives
q^2 = 1, subcritical for all M_inf < 1, hence nu == 0 and the residual
stays machine-zero -- covered by the V0-with-rho regression.

Boundary/degenerate handling: the upstream search only considers real
face-neighbors; if no neighbor face has inflow (all dots >= 0, e.g. a
stagnation-like or wall-corner element), u(e) = e and the upwind term
vanishes identically.
"""

import os

import numba
import numpy as np

from pyfp3d.mesh.metrics import build_face_adjacency
from pyfp3d.physics.isentropic import GAMMA, mach_number_squared

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
def upstream_elements(
    face_neighbors: np.ndarray,
    centroids: np.ndarray,
    nodes: np.ndarray,
    elements: np.ndarray,
    grad: np.ndarray,
    out: np.ndarray,
    max_hops: int = 4,
) -> None:
    """u(e) = multi-hop upstream walk (see module docstring): step to the
    face neighbor with the most negative centroid displacement along
    V_e, until the accumulated displacement reaches ~the element's own
    streamwise extent (both projected on the UNnormalized grad[e], so no
    normalization is needed). u(e) = e when no neighbor lies upstream.
    Recomputed every nonlinear iteration (design.md Sec 3)."""
    n_tets = len(face_neighbors)
    for e in prange(n_tets):
        gx = grad[e, 0]
        gy = grad[e, 1]
        gz = grad[e, 2]
        # Streamwise extent of e, projected on grad[e].
        pmin = 1e300
        pmax = -1e300
        for k in range(4):
            nd = elements[e, k]
            p = nodes[nd, 0] * gx + nodes[nd, 1] * gy + nodes[nd, 2] * gz
            if p < pmin:
                pmin = p
            if p > pmax:
                pmax = p
        target = 0.8 * (pmax - pmin)

        cur = np.int64(e)
        reach = 0.0
        hops = 0
        walking = True
        while walking and hops < max_hops and reach < target:
            best = np.int64(-1)
            best_disp = 0.0
            for f in range(4):
                nb = face_neighbors[cur, f]
                if nb >= 0:
                    d = (
                        (centroids[nb, 0] - centroids[cur, 0]) * gx
                        + (centroids[nb, 1] - centroids[cur, 1]) * gy
                        + (centroids[nb, 2] - centroids[cur, 2]) * gz
                    )
                    if d < best_disp:
                        best_disp = d
                        best = np.int64(nb)
            if best < 0:
                walking = False
            else:
                cur = best
                reach += -best_disp
                hops += 1
        out[e] = cur


@_njit(cache=True, fastmath=True, parallel=True)
def rho_tilde_sweep(
    q2: np.ndarray,
    rho: np.ndarray,
    upstream: np.ndarray,
    m_inf: float,
    m_crit: float,
    C: float,
    gamma: float,
    rho_floor: float,
    nu_out: np.ndarray,
    rho_tilde_out: np.ndarray,
) -> None:
    """nu switch (3.2) + upwinded density (3.1') per element.

    Same guarded form as physics.isentropic.upwind_factor:
    nu = C * max(0, 1 - M_c^2 / max(M^2, M_c^2)) -- exactly 0.0 in the
    subcritical range, no division hazard at stagnation.

    rho_floor: positivity safeguard. With C > 1 a transient iterate can
    drive rho_tilde <= 0 at a nascent shock (measured: the M 0.70 -> 0.75
    continuation step made the Picard matrix indefinite and CG failed
    outright); flooring keeps the weighted stiffness SPD. Inactive at any
    subcritical state (rho_tilde == rho ~ 1) and, in practice, at
    converged transonic states -- the caller monitors n_floored."""
    n = len(q2)
    mc2 = m_crit * m_crit
    for e in prange(n):
        u = upstream[e]
        m2 = mach_number_squared(q2[e], m_inf, gamma)
        m2u = mach_number_squared(q2[u], m_inf, gamma)
        nu_e = C * max(0.0, 1.0 - mc2 / max(m2, mc2))
        nu_u = C * max(0.0, 1.0 - mc2 / max(m2u, mc2))
        # Shock-point operator: the first subsonic element downstream of
        # a shock has nu_e = 0 (purely central) while sitting on the
        # largest density jump in the field -- measured to flip-flop and
        # pin the transonic residual at a ~4e-4 limit-cycle floor. Taking
        # nu = max(nu_e, nu_upstream) extends the upstream-side
        # dissipation one cell through the shock (the unstructured analog
        # of Murman's shock-point operator). Exactly 0 in subcritical
        # flow (both switches vanish) -- gate G4.2 unaffected.
        nu = max(nu_e, nu_u)
        nu_out[e] = nu
        rt = rho[e] - nu * (rho[e] - rho[u])
        rho_tilde_out[e] = rt if rt > rho_floor else rho_floor


class UpwindOperator:
    """Per-mesh upwinding workspace: face adjacency + outward face normals
    precomputed once, upstream/nu/rho_tilde buffers preallocated; the
    per-iteration entry point allocates nothing.

    Usage (inside the Picard loop, design.md Sec 8):
        upw = UpwindOperator(nodes, elements)
        rho_tilde = upw.rho_tilde(grad, q2, rho, m_inf, C=1.5, m_crit=0.95)
        # monitors: upw.nu (per element), upw.nu_max, upw.n_supersonic
    """

    def __init__(self, nodes: np.ndarray, elements: np.ndarray):
        self.face_neighbors, _ = build_face_adjacency(np.ascontiguousarray(elements))
        self.centroids = nodes[elements].mean(axis=1)
        self._nodes = np.ascontiguousarray(nodes, dtype=np.float64)
        self._elements = np.ascontiguousarray(elements)
        n_tets = len(elements)
        self._upstream = np.empty(n_tets, dtype=np.int64)
        self.nu = np.empty(n_tets, dtype=np.float64)
        self._rho_tilde = np.empty(n_tets, dtype=np.float64)
        self.nu_max = 0.0
        self.n_supersonic = 0
        self.n_floored = 0

    def rho_tilde(
        self,
        grad: np.ndarray,
        q2: np.ndarray,
        rho: np.ndarray,
        m_inf: float,
        C: float,
        m_crit: float,
        gamma: float = GAMMA,
        rho_floor: float = 0.05,
    ) -> np.ndarray:
        """Upwinded element densities (view into the workspace buffer --
        consume before the next call). Also refreshes the nu/floor
        monitors (nu_max, n_supersonic, n_floored)."""
        upstream_elements(self.face_neighbors, self.centroids,
                          self._nodes, self._elements, grad,
                          self._upstream)
        rho_tilde_sweep(q2, rho, self._upstream, m_inf, m_crit, C, gamma,
                        rho_floor, self.nu, self._rho_tilde)
        self.nu_max = float(self.nu.max())
        self.n_supersonic = int(np.count_nonzero(self.nu > 0.0))
        self.n_floored = int(np.count_nonzero(self._rho_tilde == rho_floor))
        return self._rho_tilde
