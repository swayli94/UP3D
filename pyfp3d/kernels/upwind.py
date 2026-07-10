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
from pyfp3d.physics.isentropic import (
    GAMMA,
    density_derivative_wrt_q_sq,
    mach_number_squared,
    mach_squared_derivative_wrt_q_sq,
)

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
def compute_face_normals(
    nodes: np.ndarray,
    elements: np.ndarray,
    out: np.ndarray,
) -> None:
    """Outward unit normal of each tet face, indexed to match
    `build_face_adjacency` (face f is opposite local node f; face nodes are
    metrics.face_defs[f]). `out` is (n_tets, 4, 3), preallocated. Precomputed
    once per mesh (geometry only); the P6 weighted upwind (3.4) uses these to
    weight face neighbours by inflow alignment max(0, -V.n_hat)."""
    n_tets = len(elements)
    # face_defs[f] = the 3 face-node local indices for face f (opposite node f)
    fd = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]], dtype=np.int64)
    for e in prange(n_tets):
        for f in range(4):
            a = elements[e, fd[f, 0]]
            b = elements[e, fd[f, 1]]
            c = elements[e, fd[f, 2]]
            apex = elements[e, f]
            e1x = nodes[b, 0] - nodes[a, 0]
            e1y = nodes[b, 1] - nodes[a, 1]
            e1z = nodes[b, 2] - nodes[a, 2]
            e2x = nodes[c, 0] - nodes[a, 0]
            e2y = nodes[c, 1] - nodes[a, 1]
            e2z = nodes[c, 2] - nodes[a, 2]
            nx = e1y * e2z - e1z * e2y
            ny = e1z * e2x - e1x * e2z
            nz = e1x * e2y - e1y * e2x
            nn = np.sqrt(nx * nx + ny * ny + nz * nz)
            if nn > 0.0:
                nx /= nn
                ny /= nn
                nz /= nn
            # Orient outward: face centroid minus apex points out of the tet.
            ox = (nodes[a, 0] + nodes[b, 0] + nodes[c, 0]) / 3.0 - nodes[apex, 0]
            oy = (nodes[a, 1] + nodes[b, 1] + nodes[c, 1]) / 3.0 - nodes[apex, 1]
            oz = (nodes[a, 2] + nodes[b, 2] + nodes[c, 2]) / 3.0 - nodes[apex, 2]
            if nx * ox + ny * oy + nz * oz < 0.0:
                nx = -nx
                ny = -ny
                nz = -nz
            out[e, f, 0] = nx
            out[e, f, 1] = ny
            out[e, f, 2] = nz


@_njit(cache=True, fastmath=True, parallel=True)
def rho_tilde_weighted_sweep(
    q2: np.ndarray,
    rho: np.ndarray,
    grad: np.ndarray,
    nodes: np.ndarray,
    elements: np.ndarray,
    centroids: np.ndarray,
    face_neighbors: np.ndarray,
    m_inf: float,
    m_crit: float,
    C: float,
    gamma: float,
    rho_floor: float,
    p_weight: float,
    reach_frac: float,
    reach_gain_max: float,
    nu_out: np.ndarray,
    rho_tilde_out: np.ndarray,
) -> None:
    """P6 differentiable artificial density (design.md §3.2, Eq. 3.4).

    Replaces the P4 integer walk `u(e)` + per-element `rho[u]` by a *smooth*
    upstream-weighted blend over the face neighbours, so rho_tilde varies
    continuously in space (no cell-to-cell selection flip -> no sawtooth) and
    is differentiable in phi at fixed neighbour set (Newton-ready, P7).

    Weighting is by the neighbour's CENTROID DISPLACEMENT along the reversed
    streamline, NOT the face normal: on prism-split sliver tets a face can
    face upstream (`-V.n_hat > 0`) while its neighbour centroid is not actually
    upstream, giving a net-DOWNSTREAM blend (`reach < 0`) that makes the term
    anti-dissipative and blows the transonic solve up (measured). Using the
    centroid displacement mirrors the validated walk's criterion and keeps
    `reach >= 0` by construction:

        a_f    = -V_hat . (c_nb - c_e) / |c_nb - c_e|   (cos to reversed flow)
        w_f    = max(0, a_f)^p                          (upstream neighbours)
        rho_up = sum_f w_f rho[nb_f] / sum_f w_f
        m2_up  = sum_f w_f M^2[nb_f] / sum_f w_f
        reach  = sum_f w_f (-V_hat . (c_nb - c_e)) / sum_f w_f   (>= 0)

    Reach compensation restores the multi-hop walk's effective dissipation
    without a discrete walk (the walk existed only to reach ~one streamwise
    extent on sliver tets, design.md §3 / P4 hardening trail):

        G = clip(reach_frac * extent_e / max(reach, eps), 1, reach_gain_max)
        rho_tilde = rho_e - max(nu_e, nu_up) * G * (rho_e - rho_up)

    extent_e is the element's streamwise span projected on V_hat. The
    shock-point operator stays a HARD max(nu_e, nu_up): it is exactly 0 when
    both are subcritical (bit-exact no-op, gate G4.2) and matches Lopez's
    mu_s = mu_c max(0, mu, mu_up); the sawtooth came from the flipping rho_up,
    not from this max, so smoothing it is unnecessary (and would break the
    no-op, since max_eps(0,0) = eps != 0).
    """
    n = len(q2)
    mc2 = m_crit * m_crit
    for e in prange(n):
        gx = grad[e, 0]
        gy = grad[e, 1]
        gz = grad[e, 2]
        gmag = np.sqrt(gx * gx + gy * gy + gz * gz)
        m2 = mach_number_squared(q2[e], m_inf, gamma)
        nu_e = C * max(0.0, 1.0 - mc2 / max(m2, mc2))
        # Stagnation / no gradient: no streamwise direction, upwind term
        # vanishes. Also the subcritical fast path is handled after nu_up.
        if gmag <= 0.0:
            nu_out[e] = nu_e
            rho_tilde_out[e] = rho[e]
            continue
        vhx = gx / gmag
        vhy = gy / gmag
        vhz = gz / gmag

        w_sum = 0.0
        rho_up_acc = 0.0
        m2_up_acc = 0.0
        reach_acc = 0.0
        for f in range(4):
            nb = face_neighbors[e, f]
            if nb < 0:
                continue
            dcx = centroids[nb, 0] - centroids[e, 0]
            dcy = centroids[nb, 1] - centroids[e, 1]
            dcz = centroids[nb, 2] - centroids[e, 2]
            dlen = np.sqrt(dcx * dcx + dcy * dcy + dcz * dcz)
            if dlen <= 0.0:
                continue
            # Streamwise upstream distance and its cosine (both use the same
            # reversed-flow projection; only genuinely upstream neighbours,
            # a_f > 0, contribute -> reach >= 0).
            up_dist = -(vhx * dcx + vhy * dcy + vhz * dcz)
            a = up_dist / dlen
            if a <= 0.0:
                continue
            w = a ** p_weight
            w_sum += w
            rho_up_acc += w * rho[nb]
            m2_up_acc += w * mach_number_squared(q2[nb], m_inf, gamma)
            reach_acc += w * up_dist

        if w_sum <= 0.0:
            # No upstream neighbour (wall corner / all boundary): no upwind term.
            nu_out[e] = nu_e
            rho_tilde_out[e] = rho[e]
            continue

        rho_up = rho_up_acc / w_sum
        m2_up = m2_up_acc / w_sum
        reach = reach_acc / w_sum
        nu_u = C * max(0.0, 1.0 - mc2 / max(m2_up, mc2))
        nu = nu_e if nu_e > nu_u else nu_u
        nu_out[e] = nu
        if nu <= 0.0:
            rho_tilde_out[e] = rho[e]
            continue

        # Streamwise extent of e projected on V_hat, for reach compensation.
        pmin = 1e300
        pmax = -1e300
        for k in range(4):
            nd = elements[e, k]
            pr = (nodes[nd, 0] * vhx + nodes[nd, 1] * vhy + nodes[nd, 2] * vhz)
            if pr < pmin:
                pmin = pr
            if pr > pmax:
                pmax = pr
        reach_target = reach_frac * (pmax - pmin)
        g = reach_target / max(reach, 1e-30)
        if g < 1.0:
            g = 1.0
        elif g > reach_gain_max:
            g = reach_gain_max

        rt = rho[e] - nu * g * (rho[e] - rho_up)
        rho_tilde_out[e] = rt if rt > rho_floor else rho_floor


def build_upstream_neighborhoods(face_neighbors: np.ndarray, depth: int = 3):
    """BFS face-neighbour neighbourhood of each element out to `depth` hops,
    as a CSR (offsets, indices) pair (excludes e itself). One-time geometric
    precompute for the streamline-kernel upstream sample: on sliver prism-split
    tets the genuinely-upstream cell (~0.8 streamwise extent away) is several
    hops off, so a single face ring cannot reach it -- the multi-hop walk
    existed for exactly this reason (design.md §3 / P4 hardening trail)."""
    n = len(face_neighbors)
    offsets = np.zeros(n + 1, dtype=np.int64)
    idx_lists = []
    for e in range(n):
        seen = {e}
        frontier = [e]
        for _ in range(depth):
            nxt = []
            for c in frontier:
                for f in range(4):
                    nb = int(face_neighbors[c, f])
                    if nb >= 0 and nb not in seen:
                        seen.add(nb)
                        nxt.append(nb)
            frontier = nxt
        seen.discard(e)
        lst = sorted(seen)
        idx_lists.append(lst)
        offsets[e + 1] = offsets[e] + len(lst)
    indices = (np.concatenate([np.array(l, dtype=np.int64) for l in idx_lists])
               if idx_lists else np.empty(0, dtype=np.int64))
    return offsets, indices


@_njit(cache=True, fastmath=True, parallel=True)
def rho_tilde_kernel_sweep(
    q2: np.ndarray,
    rho: np.ndarray,
    grad: np.ndarray,
    nodes: np.ndarray,
    elements: np.ndarray,
    centroids: np.ndarray,
    nb_off: np.ndarray,
    nb_idx: np.ndarray,
    m_inf: float,
    m_crit: float,
    C: float,
    gamma: float,
    rho_floor: float,
    reach_frac: float,
    sigma_s_frac: float,
    sigma_p_frac: float,
    nu_out: np.ndarray,
    rho_tilde_out: np.ndarray,
) -> None:
    """P6 differentiable artificial density (design.md §3.2): streamline-kernel
    upstream sample.

    rho_up is a Gaussian-kernel average of rho over the multi-ring neighbourhood
    (`nb_off`/`nb_idx`), centred on the point one streamwise extent upstream:

        target = reach_frac * extent_e            (streamwise reach ~ the walk)
        up_c   = -V_hat . (c_c - c_e)              (upstream distance of cand. c)
        perp2  = |c_c - c_e|^2 - up_c^2            (off-streamline distance^2)
        w_c    = exp(-0.5[(up_c-target)^2/sig_s^2 + perp2/sig_p^2]),  up_c > 0
        rho_up = sum_c w_c rho[c] / sum_c w_c

    This samples GENUINELY upstream cells (not the near ring), so it keeps the
    walk's stabilising reach, while the smooth kernel removes the walk's
    cell-to-cell selection flip (the sawtooth) and is C-inf in V_hat=grad/|grad|
    (Newton-ready at fixed neighbourhood, P7). No reach-amplification factor is
    needed -- the reach is in the sampling, not a multiplier (the near-ring blend
    + amplification was measured to diverge; a genuine-reach sample does not).

    Shock-point operator: HARD max(nu_e, nu_up) with nu_up from the same kernel
    blend of neighbour M^2 -- exactly 0 subcritically (bit no-op, gate G4.2).
    """
    n = len(q2)
    mc2 = m_crit * m_crit
    for e in prange(n):
        gx = grad[e, 0]
        gy = grad[e, 1]
        gz = grad[e, 2]
        gmag = np.sqrt(gx * gx + gy * gy + gz * gz)
        m2 = mach_number_squared(q2[e], m_inf, gamma)
        nu_e = C * max(0.0, 1.0 - mc2 / max(m2, mc2))
        if gmag <= 0.0:
            nu_out[e] = nu_e
            rho_tilde_out[e] = rho[e]
            continue
        vhx = gx / gmag
        vhy = gy / gmag
        vhz = gz / gmag

        # Streamwise extent of e projected on V_hat -> kernel length scales.
        pmin = 1e300
        pmax = -1e300
        for k in range(4):
            nd = elements[e, k]
            pr = nodes[nd, 0] * vhx + nodes[nd, 1] * vhy + nodes[nd, 2] * vhz
            if pr < pmin:
                pmin = pr
            if pr > pmax:
                pmax = pr
        extent = pmax - pmin
        target = reach_frac * extent
        inv2ss = 1.0 / (2.0 * (sigma_s_frac * extent) ** 2 + 1e-300)
        inv2sp = 1.0 / (2.0 * (sigma_p_frac * extent) ** 2 + 1e-300)

        w_sum = 0.0
        rho_up_acc = 0.0
        m2_up_acc = 0.0
        for j in range(nb_off[e], nb_off[e + 1]):
            c = nb_idx[j]
            dcx = centroids[c, 0] - centroids[e, 0]
            dcy = centroids[c, 1] - centroids[e, 1]
            dcz = centroids[c, 2] - centroids[e, 2]
            up = -(vhx * dcx + vhy * dcy + vhz * dcz)
            if up <= 0.0:
                continue
            d2 = dcx * dcx + dcy * dcy + dcz * dcz
            perp2 = d2 - up * up
            if perp2 < 0.0:
                perp2 = 0.0
            w = np.exp(-((up - target) * (up - target) * inv2ss + perp2 * inv2sp))
            w_sum += w
            rho_up_acc += w * rho[c]
            m2_up_acc += w * mach_number_squared(q2[c], m_inf, gamma)

        if w_sum <= 0.0:
            nu_out[e] = nu_e
            rho_tilde_out[e] = rho[e]
            continue

        rho_up = rho_up_acc / w_sum
        m2_up = m2_up_acc / w_sum
        nu_u = C * max(0.0, 1.0 - mc2 / max(m2_up, mc2))
        nu = nu_e if nu_e > nu_u else nu_u
        nu_out[e] = nu
        if nu <= 0.0:
            rho_tilde_out[e] = rho[e]
            continue
        rt = rho[e] - nu * (rho[e] - rho_up)
        rho_tilde_out[e] = rt if rt > rho_floor else rho_floor


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


@_njit(cache=True, fastmath=True, parallel=True)
def rho_tilde_sensitivities_sweep(
    q2: np.ndarray,
    rho: np.ndarray,
    upstream: np.ndarray,
    m_inf: float,
    m_crit: float,
    C: float,
    gamma: float,
    rho_floor: float,
    se_out: np.ndarray,
    su_out: np.ndarray,
) -> None:
    """P7: exact per-element sensitivities of the walk flux at FROZEN
    upstream selection (design.md §6.3 / López B.3–B.6),

        s_e = ∂rho_tilde_e/∂q²_e ,   s_u = ∂rho_tilde_e/∂q²_u(e) ,

    so that the Newton Term-2/Term-3 chain is
    ∂rho_tilde_e/∂φ_k = s_e·(2 ∇φ_e·∇N_k|_e) + s_u·(2 ∇φ_u·∇N_k|_u) —
    the geometric factor lives with the caller (P8 assembly / FD tests).

    Branch structure of rho_tilde_sweep, differentiated branch-wise with
    u(e) frozen (López differentiates through ρ_up/ν/ρ only, never the
    selection; the C⁰ kinks at ν_e = ν_u and at the switch threshold are
    measure-zero and freeze near the solution):

      floored (rt <= rho_floor, flat clamp)        : s_e = s_u = 0
      u == e (upwind term identically 0)           : s_e = ρ'_e, s_u = 0
      subsonic (ν_e = ν_u = 0)                     : s_e = ρ'_e, s_u = 0
      accelerating (ν = ν_e >= ν_u, ν_e > 0)       : s_e = ρ'_e(1−ν) − (ρ_e−ρ_u)·ν'_e
                                                     s_u = ν·ρ'_u
      shock-point (ν = ν_u > ν_e)                  : s_e = ρ'_e(1−ν)
                                                     s_u = ν·ρ'_u − (ρ_e−ρ_u)·ν'_u

    with ρ' = dρ/dq² and, on the ACTIVE switch branch only (M² > M_c²),
    ν' = C·(M_c²/M⁴)·(dM²/dq²) > 0 (else 0 — the max(0,·) outer clamp).
    Must mirror rho_tilde_sweep's guarded ν form exactly so the FD check
    against the shipped flux is meaningful."""
    n = len(q2)
    mc2 = m_crit * m_crit
    for e in prange(n):
        u = upstream[e]
        rho_e = rho[e]
        drho_e = density_derivative_wrt_q_sq(q2[e], m_inf, gamma)
        if u == e:
            # rho_tilde == rho_e (upwind term identically 0), but the
            # floor still applies to rho_e itself in rho_tilde_sweep.
            if rho_e <= rho_floor:
                se_out[e] = 0.0
                su_out[e] = 0.0
            else:
                se_out[e] = drho_e
                su_out[e] = 0.0
            continue
        m2 = mach_number_squared(q2[e], m_inf, gamma)
        m2u = mach_number_squared(q2[u], m_inf, gamma)
        nu_e = C * max(0.0, 1.0 - mc2 / max(m2, mc2))
        nu_u = C * max(0.0, 1.0 - mc2 / max(m2u, mc2))
        nu = max(nu_e, nu_u)
        rho_u = rho[u]
        rt = rho_e - nu * (rho_e - rho_u)
        if rt <= rho_floor:
            se_out[e] = 0.0
            su_out[e] = 0.0
            continue
        if nu == 0.0:
            se_out[e] = drho_e
            su_out[e] = 0.0
            continue
        drho_u = density_derivative_wrt_q_sq(q2[u], m_inf, gamma)
        djump = rho_e - rho_u
        if nu_e >= nu_u:
            # accelerating: nu tracks the element's own switch
            dnu_e = C * (mc2 / (m2 * m2)) * mach_squared_derivative_wrt_q_sq(
                q2[e], m_inf, gamma) if m2 > mc2 else 0.0
            se_out[e] = drho_e * (1.0 - nu) - djump * dnu_e
            su_out[e] = nu * drho_u
        else:
            # shock-point: nu tracks the upstream switch (nu_u > nu_e >= 0
            # implies M²_u > M_c², so the upstream branch is active)
            dnu_u = C * (mc2 / (m2u * m2u)) * mach_squared_derivative_wrt_q_sq(
                q2[u], m_inf, gamma)
            se_out[e] = drho_e * (1.0 - nu)
            su_out[e] = nu * drho_u - djump * dnu_u


class UpwindOperator:
    """Per-mesh upwinding workspace: face adjacency + outward face normals
    precomputed once, upstream/nu/rho_tilde buffers preallocated; the
    per-iteration entry point allocates nothing.

    Usage (inside the Picard loop, design.md Sec 8):
        upw = UpwindOperator(nodes, elements)
        rho_tilde = upw.rho_tilde(grad, q2, rho, m_inf, C=1.5, m_crit=0.95)
        # monitors: upw.nu (per element), upw.nu_max, upw.n_supersonic
    """

    def __init__(self, nodes: np.ndarray, elements: np.ndarray,
                 weighted: bool = False, mode: str = "kernel",
                 p_weight: float = 3.0, reach_frac: float = 1.0,
                 reach_gain_max: float = 4.0, nbr_depth: int = 3,
                 sigma_s_frac: float = 0.35, sigma_p_frac: float = 0.35):
        self.face_neighbors, _ = build_face_adjacency(np.ascontiguousarray(elements))
        self.centroids = nodes[elements].mean(axis=1)
        self._nodes = np.ascontiguousarray(nodes, dtype=np.float64)
        self._elements = np.ascontiguousarray(elements)
        n_tets = len(elements)
        # P6 differentiable flux. weighted=False restores the P4 integer walk.
        # mode selects the weighted operator: "kernel" (streamline-Gaussian
        # sample over a multi-ring neighbourhood -- the shipped P6 path, real
        # reach + smooth) or "blend" (single-ring centroid blend -- kept for
        # the record; measured transiently unstable on sliver tets).
        self.weighted = weighted
        self.mode = mode
        self.p_weight = p_weight
        self.reach_frac = reach_frac
        self.reach_gain_max = reach_gain_max
        self.sigma_s_frac = sigma_s_frac
        self.sigma_p_frac = sigma_p_frac
        self._nb_off = None
        self._nb_idx = None
        if weighted and mode == "kernel":
            self._nb_off, self._nb_idx = build_upstream_neighborhoods(
                self.face_neighbors, depth=nbr_depth)
        self._upstream = np.empty(n_tets, dtype=np.int64)
        self.nu = np.empty(n_tets, dtype=np.float64)
        self._rho_tilde = np.empty(n_tets, dtype=np.float64)
        # P7 frozen-selection sensitivity buffers (rho_tilde_sensitivities)
        self._se = np.empty(n_tets, dtype=np.float64)
        self._su = np.empty(n_tets, dtype=np.float64)
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
        monitors (nu_max, n_supersonic, n_floored). Uses the P6 streamline
        kernel by default; `weighted=False` restores the P4 walk."""
        if self.weighted and self.mode == "kernel":
            rho_tilde_kernel_sweep(
                q2, rho, grad, self._nodes, self._elements, self.centroids,
                self._nb_off, self._nb_idx, m_inf, m_crit, C, gamma, rho_floor,
                self.reach_frac, self.sigma_s_frac, self.sigma_p_frac,
                self.nu, self._rho_tilde)
        elif self.weighted:
            rho_tilde_weighted_sweep(
                q2, rho, grad, self._nodes, self._elements, self.centroids,
                self.face_neighbors, m_inf, m_crit, C, gamma, rho_floor,
                self.p_weight, self.reach_frac, self.reach_gain_max,
                self.nu, self._rho_tilde)
        else:
            upstream_elements(self.face_neighbors, self.centroids,
                              self._nodes, self._elements, grad,
                              self._upstream)
            rho_tilde_sweep(q2, rho, self._upstream, m_inf, m_crit, C, gamma,
                            rho_floor, self.nu, self._rho_tilde)
        self.nu_max = float(self.nu.max())
        self.n_supersonic = int(np.count_nonzero(self.nu > 0.0))
        self.n_floored = int(np.count_nonzero(self._rho_tilde == rho_floor))
        return self._rho_tilde

    def rho_tilde_sensitivities(
        self,
        grad: np.ndarray,
        q2: np.ndarray,
        rho: np.ndarray,
        m_inf: float,
        C: float,
        m_crit: float,
        gamma: float = GAMMA,
        rho_floor: float = 0.05,
    ):
        """P7 (gate G7.3): exact ∂rho_tilde/∂q² sensitivities of the WALK
        flux at frozen upstream selection — the P8 Newton Term-2/Term-3
        physics factor (design.md §6.3). Freezes u(e) with the same
        upstream walk as `rho_tilde` on the same grad, then sweeps the
        branch-wise derivative (see rho_tilde_sensitivities_sweep).

        Returns (s_e, s_u, upstream): views into workspace buffers —
        consume before the next call. s_e = ∂rho_tilde_e/∂q²_e,
        s_u = ∂rho_tilde_e/∂q²_{u(e)} (0 whenever u(e) == e).

        Walk mode only: the kernel-mode (weighted) flux has a dense
        neighbourhood dependence — its Jacobian is a P8 measurement item,
        not a P7 deliverable (design.md §3.2)."""
        if self.weighted:
            raise NotImplementedError(
                "rho_tilde_sensitivities is defined for the walk flux "
                "(weighted=False); kernel-mode sensitivities are a P8 "
                "candidate, not implemented.")
        upstream_elements(self.face_neighbors, self.centroids,
                          self._nodes, self._elements, grad,
                          self._upstream)
        rho_tilde_sensitivities_sweep(q2, rho, self._upstream, m_inf,
                                      m_crit, C, gamma, rho_floor,
                                      self._se, self._su)
        return self._se, self._su, self._upstream
