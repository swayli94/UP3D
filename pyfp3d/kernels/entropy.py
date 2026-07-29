"""Entropy-corrected (non-isentropic) full-potential density factor.

Phase-two gate GS1b.3. Pre-registration:
docs/dev_phase_two/20260729-0700-s1b-entropy-implementation.md.

WHY THIS EXISTS (measured, not assumed). At a fixed transonic condition the
computed lift does not converge under refinement (M0.7875/alpha1.25:
cl 0.3725 / 0.5234 / 0.5686 over the three levels), and the only shock-weakening
mechanism in the code -- the artificial density -- is ~O(h): at the finest level
upwind_c = 1.5 and 3.0 give nearly the same answer (GS1b.2 Q5). A mechanism that
vanishes with h cannot cure an error that does not, so the weakening has to come
from the physics. Isentropic full potential over-predicts the density jump across
a shock by 5.1 % / 6.9 % / 8.9 % at M1 = 1.30 / 1.35 / 1.40, exactly the strengths
the airfoil cases run at (GS1b.2 Q1).

THE RELATION. For a perfect gas whose total enthalpy is constant along
streamlines, rho0 = p0/(R*T0) and T0 is preserved across a shock, so

    rho02/rho01 = p02/p01   EXACTLY,   and   rho_s = sigma * rho_isen(q^2)

with sigma = p02/p01 the Rankine-Hugoniot total-pressure ratio at the pre-shock
Mach number. Applying that factor makes the full-potential mass-conservation jump
reproduce Rankine-Hugoniot identically (verified to 0.0000 % for M1 = 1.15..1.60
in bench/s1_duct/run_entropy_premise.py). It is NOT sigma^(1/(gamma-1)) -- that
form (inherited from some phase-one notes) applies the exponent twice; the
numerical check settles it: at M1 = 1.35 sigma^1 gives the R-H ratio 1.60278
while sigma^(1/(gamma-1)) gives 1.38162.

HOW sigma IS BUILT (both steps reuse the existing upwinding data -- no new
geometry or topology machinery):

  detection   Each element already has a unique upstream donor u(e) from
              `upwind.upstream_elements`. A shock is a supersonic-to-subsonic
              transition along a streamline, so: donor supersonic AND element
              subsonic => `e` is a post-shock element. M1 is the PEAK Mach over
              the supersonic run upstream of it (NOT the donor's own Mach -- the
              artificial density smears the shock, and the pointwise value reads
              the smearing; see shock_factor_sweep). s_e = p02/p01(M1), else 1.
  transport   Entropy is constant along a streamline, so sigma_e = sigma_u(e)*s_e,
              i.e. sigma is the product of the local factors along the donor
              chain. Computed by pointer doubling (see transport_sigma) in
              O(n log depth), with donor-cycle detection.

sigma is meant to be FROZEN over a Newton step by the caller (residual, Jacobian
and line search sharing one sigma, refreshed between steps): unlike the upstream
selection, which is piecewise constant in phi, sigma depends CONTINUOUSLY on phi
through p02/p01(M1), so a live sigma with no d(sigma)/d(phi) term would make the
Jacobian genuinely inexact rather than exact-almost-everywhere.
"""

import os

import numba
import numpy as np

from pyfp3d.physics.isentropic import GAMMA, mach_number_squared

if os.environ.get("PYFP3D_NOJIT", "0") == "1":
    prange = range

    def _njit(*args, **kwargs):
        def deco(fn):
            return fn
        return deco
else:
    from numba import prange

    def _njit(*args, **kwargs):
        return numba.njit(*args, **kwargs)


@_njit(cache=True, fastmath=True)
def total_pressure_ratio(m1: float, gamma: float = GAMMA) -> float:
    """Rankine-Hugoniot total-pressure ratio p02/p01 for a normal shock at
    pre-shock Mach m1. Returns 1.0 at or below m1 = 1 (no shock, no entropy)."""
    if m1 <= 1.0:
        return 1.0
    m2 = m1 * m1
    a = ((gamma + 1.0) * m2) / ((gamma - 1.0) * m2 + 2.0)
    b = (gamma + 1.0) / (2.0 * gamma * m2 - (gamma - 1.0))
    return a ** (gamma / (gamma - 1.0)) * b ** (1.0 / (gamma - 1.0))


@_njit(cache=True, fastmath=True, parallel=True)
def shock_factor_sweep(
    q2: np.ndarray,
    upstream: np.ndarray,
    m_inf: float,
    gamma: float,
    max_walk: int,
    s_out: np.ndarray,
    m1_out: np.ndarray,
) -> None:
    """Per-element local entropy factor s_e (see module docstring, `detection`).

    An element is POST-SHOCK when its donor is supersonic and it is itself
    subsonic. M1 is then the PEAK Mach over the supersonic run immediately
    upstream along the donor chain -- not the donor's own Mach.

    ★ Why the peak, measured: the first implementation used the donor's Mach and
    produced sigma_min = 0.9968 (a 0.32 % density cut, i.e. M1 ~ 1.14) on a state
    whose M_max was 1.366. The artificial density smears the shock over two or
    three cells, so along a streamline the Mach runs 1.37 -> 1.25 -> 1.10 -> 0.95
    and the pointwise test fires at the LAST supersonic cell -- inside the shock
    structure, where the Mach has already dropped. It read the smearing, not the
    shock. Walking upstream through the supersonic cells and taking the peak
    recovers the physical pre-shock Mach: on an airfoil the pocket accelerates
    from the sonic line to the shock, so the peak sits at the shock foot.

    The walk stops at the first subsonic cell (so it never leaves the pocket it
    belongs to, and a second, upstream shock cannot contaminate it), at a chain
    root, or after `max_walk` hops (cycle protection).

    m1_out records the detected pre-shock Mach (0.0 where no shock was detected)
    for diagnostics.
    """
    n = len(q2)
    for e in prange(n):
        s_out[e] = 1.0
        m1_out[e] = 0.0
        u = upstream[e]
        if u == e:
            continue
        m2e = mach_number_squared(q2[e], m_inf, gamma)
        if m2e >= 1.0:
            continue                      # still supersonic: not post-shock
        m2u = mach_number_squared(q2[u], m_inf, gamma)
        if m2u <= 1.0:
            continue                      # donor subsonic: no shock crossed
        m2_peak = m2u
        c = u
        for _ in range(max_walk):
            nxt = upstream[c]
            if nxt == c:
                break                     # chain root
            m2n = mach_number_squared(q2[nxt], m_inf, gamma)
            if m2n <= 1.0:
                break                     # front of the supersonic pocket
            if m2n > m2_peak:
                m2_peak = m2n
            c = nxt
        m1 = np.sqrt(m2_peak)
        s_out[e] = total_pressure_ratio(m1, gamma)
        m1_out[e] = m1


@_njit(cache=True, parallel=True)
def transport_sigma(
    s: np.ndarray,
    upstream: np.ndarray,
    n_round: int,
    sigma_out: np.ndarray,
    anc_a: np.ndarray,
    anc_b: np.ndarray,
    prod_b: np.ndarray,
) -> int:
    """sigma_e = product of s over the donor chain from e up to its root, by
    POINTER DOUBLING.

    The obvious implementation -- sweep `sigma_e = sigma_u(e) * s_e` until
    nothing moves -- was rejected for two measured-in-advance reasons: its cost is
    O(n * depth) (the depth is the whole downstream extent, hundreds of cells, so
    ~2e8 serial operations per Newton step on the fine mesh), and a donor CYCLE
    (possible on tie-degenerate meshes) makes it multiply sigma down toward zero
    once per sweep -- a silent catastrophe in the density.

    Pointer doubling fixes both. With A[e] the current ancestor and P[e] the
    product of s over the path from e up to (not including) A[e]:

        P'[e] = P[e] * P[A[e]],    A'[e] = A[A[e]]

    doubles the covered depth each round, so ceil(log2(depth)) rounds suffice.
    A root r has upstream[r] == r and s[r] == 1 by construction (the detector
    skips u == e), so P'[r] = P[r] * P[r] = 1 and A'[r] = r: the recursion is
    idempotent once every chain has reached its root.

    ★ The convergence test is "every ancestor is a GENUINE root"
    (upstream[A[e]] == A[e]), NOT "the ancestor pointers stopped moving". The
    weaker test is wrong and the first implementation shipped it: in a 2-cycle
    u(0) = 1, u(1) = 0 the pointers settle at A[0] = 0, A[1] = 1 after one round
    -- each element becomes its OWN ancestor and looks exactly like a root --
    while the product keeps SQUARING every round (P[0] -> (s0*s1)^2 -> ^4 ...).
    A unit test with that donor map reported "converged" with a corrupted sigma
    until the test was written. With the genuine-root test a cycle can never
    settle, the round cap binds, and this returns n_round -- which the caller must
    treat as "not converged" rather than use the numbers (the GS1.4
    clamp-not-silent contract).

    Deterministic under threading: each round reads only the previous round's
    buffers, so the result does not depend on scheduling (phase-one discipline
    #12 -- a non-reproducible kernel breaks every A/B).

    Returns the number of rounds taken (< n_round means converged).
    """
    n = len(s)
    for e in prange(n):
        anc_a[e] = upstream[e]
        sigma_out[e] = s[e]
    for it in range(n_round):
        for e in prange(n):
            a = anc_a[e]
            prod_b[e] = sigma_out[e] * sigma_out[a]
            anc_b[e] = anc_a[a]
        n_unsettled = 0
        for e in prange(n):
            a = anc_b[e]
            if upstream[a] != a:
                n_unsettled += 1      # a += reduction: the form numba supports
            anc_a[e] = a
            sigma_out[e] = prod_b[e]
        if n_unsettled == 0:
            return it + 1
    return n_round


class EntropyOperator:
    """Per-mesh workspace for the entropy factor (allocates once).

    Usage inside a driver, once per Newton step (sigma frozen over the step):

        ent = EntropyOperator(n_elements)
        sigma = ent.sigma(q2, upstream, m_inf, gamma)   # view into the buffer
        rho_s = sigma * rho_isentropic                  # then feed rho_tilde

    Monitors after each call: `n_shock` (elements with a detected shock),
    `m1_max` (strongest detected pre-shock Mach), `sigma_min`, `n_rounds`, and
    `converged` -- False means a donor cycle defeated the transport, and the
    caller must refuse to report convergence rather than use the numbers.
    """

    #: pointer-doubling rounds. Each round doubles the covered chain depth, so
    #: 24 rounds cover 16.7 M hops -- far beyond any mesh here; the loop exits as
    #: soon as every chain has reached its root, so the cap only binds on a donor
    #: cycle, and then `converged` is False and the caller must not use sigma.
    N_ROUND_DEFAULT = 24
    #: hops the pre-shock walk may take through a supersonic pocket. The pocket
    #: is at most a chord long, i.e. a few hundred cells at the finest level;
    #: the walk stops on its own at the sonic line.
    MAX_WALK_DEFAULT = 512

    def __init__(self, n_elements: int, n_round: int = N_ROUND_DEFAULT,
                 max_walk: int = MAX_WALK_DEFAULT):
        self.n_round = int(n_round)
        self.max_walk = int(max_walk)
        self._s = np.ones(n_elements, dtype=np.float64)
        self._m1 = np.zeros(n_elements, dtype=np.float64)
        self._sigma = np.ones(n_elements, dtype=np.float64)
        self._anc_a = np.zeros(n_elements, dtype=np.int64)
        self._anc_b = np.zeros(n_elements, dtype=np.int64)
        self._prod_b = np.ones(n_elements, dtype=np.float64)
        self.n_shock = 0
        self.m1_max = 0.0
        self.sigma_min = 1.0
        self.n_rounds = 0
        self.converged = True

    def sigma(self, q2: np.ndarray, upstream: np.ndarray, m_inf: float,
              gamma: float = GAMMA) -> np.ndarray:
        """Entropy factor per element (view into the workspace buffer -- copy it
        if it must outlive the next call, which the frozen-sigma callers do)."""
        shock_factor_sweep(q2, upstream, m_inf, gamma, self.max_walk,
                           self._s, self._m1)
        self.n_rounds = transport_sigma(
            self._s, np.ascontiguousarray(upstream, dtype=np.int64),
            self.n_round, self._sigma, self._anc_a, self._anc_b, self._prod_b)
        self.converged = self.n_rounds < self.n_round
        self.n_shock = int(np.count_nonzero(self._m1 > 0.0))
        self.m1_max = float(self._m1.max()) if len(self._m1) else 0.0
        self.sigma_min = float(self._sigma.min()) if len(self._sigma) else 1.0
        return self._sigma
