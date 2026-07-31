"""Entropy-corrected (non-isentropic) full-potential density factor.

★★ GS1b.10 (2026-07-31) built a CHAIN-FREE replacement for everything below --
a flux-weighted FV upwind transport plus an additive entropy production density
-- and it is a MEASURED NEGATIVE, so this donor-chain version is what stands.
What the chain-free version fixed: no donor churn, sigma self-consistency
exactly 0.0 (the GS1b.9 P2 that failed here), four knobs removed, Newton in 8
steps instead of 23. What it could not do: LOCALISE the production. Additivity
itself was verified to 0.004-0.49 % (Q1), but charging every supersonic
compression over-produces (1-sigma 6.22 % where M_max allows 3.57 %) because the
airfoil pocket's 0.3-chord ISENTROPIC deceleration is not a shock, and all three
discriminators failed: the artificial-density activity does not separate the
populations (a continuum, p50 -0.003 to max 0.0615), the streamwise second
difference passes on synthetic streamlines but fails on a real tet mesh (noise
keeps d2/d1 away from 1), and an absolute per-cell threshold is strongly
NON-MONOTONE in the threshold (1-sigma 5.96 / 4.67 / 12.45 / 1.26 / 0.94 % over
lo = 0.004..0.030, with the shock direction flipping sign) -- a fitting
parameter by that round's own Q5. Reading: on this discretisation the
"is this cell part of a shock" information is not cleanly present in the local
field, because the shock is smeared over 2-3 IRREGULAR tet cells and the
per-cell compression noise is the same size as the discrimination needed.
The chain-free code is in git history at 774ef96; the round file is
docs/dev_phase_two/20260731-0000-s1b-chainfree-sigma.md.

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
from typing import Optional

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
    knee_frac: float,
    lim: np.ndarray,
    s_out: np.ndarray,
    m1_out: np.ndarray,
) -> None:
    """Per-element local entropy factor s_e (see module docstring, `detection`).

    An element is POST-SHOCK when its donor is supersonic and it is itself
    subsonic. M1 must then be the Mach at the UPSTREAM EDGE OF THE NUMERICAL
    SHOCK STRUCTURE, and getting that right took two measured corrections -- the
    two obvious choices are wrong in opposite directions:

      donor's own Mach     TOO WEAK. The artificial density smears the shock over
                           two or three cells, so along a streamline the Mach runs
                           1.37 -> 1.25 -> 1.10 -> 0.95 and the pointwise test
                           fires at the LAST supersonic cell, INSIDE the
                           structure. Measured: sigma_min 0.9968 (a 0.32 % cut)
                           on a state whose M_max was 1.366.
      pocket-wide peak     TOO STRONG. On an airfoil the pocket Mach is NOT
                           monotone: at M0.80/alpha1.25 it rises to 1.408 at
                           x/c 0.33 and then FALLS along the pocket to ~1.15-1.29
                           at the shock (x/c 0.63). Taking the peak read 1.3897
                           where the pre-shock value is ~1.29 -- sigma 0.9598
                           instead of 0.9812, i.e. the correction ~2x too strong,
                           which put the G4.1 shock at 0.548 against an Euler
                           anchor of 0.60-0.63.
                           ★ The Laval-nozzle bench CANNOT catch this: there the
                           flow accelerates monotonically to the shock, so
                           peak == pre-shock by construction. It took the airfoil
                           gate to expose it.

    What separates the two is the RATE. Walking upstream from the post-shock cell,
    the Mach climbs steeply while inside the shock structure and then flattens on
    leaving it (measured per-hop rises: nozzle +0.13, +0.01, negative; airfoil
    +0.09, +0.05, +0.02, then ~0.01). So the walk stops at that KNEE: it continues
    while the current hop's rise is at least `knee_frac` of the largest rise seen
    so far in this walk, and stops otherwise.

    The walk also stops at the first subsonic cell (so it never leaves the pocket
    it belongs to, and a second, upstream shock cannot contaminate it), at a chain
    root, or after `max_walk` hops (cycle protection).

    `knee_frac` is a knob and is treated as one: its sensitivity is measured
    (EntropyOperator.KNEE_FRAC_DEFAULT documents the band), and its deletion
    condition is "replace with a shock width derived from nu, or drop it if the
    result is insensitive across [0.2, 0.5]".

    `lim[e]` is True where the element's q2 was NOT touched by the m_cap limiter
    (the `lim = q2l == q2n` convention of NewtonWorkspace.eval_residual). Where it
    is False the system has ALREADY declared that cell's speed non-physical, so no
    entropy production may be computed from it: the walk stops there, and a
    post-shock cell whose own or whose donor's speed was limited gets s = 1.0.

    ★ That guard is not a safety belt, it is the fix for a measured default-path
    defect (pre-registered 20260731-2000, evidence bench/gate_results/m3_budget.csv).
    refresh_sigma passes the LIMITED field q2l, so at a limited cell the recovered
    "pre-shock Mach" was the cap itself -- m1_max read exactly 2.9999999999999996 at
    M6 medium -- giving s = sigma_RH(3.0) = 0.32834, a 67 % density cut invented out
    of a limiter. Two measured routes then carry that to a collapse (both locked in
    tests/test_s1b_entropy.py::test_mcap_the_mechanism_itself_is_locked): separate
    capped shocks along an ACYCLIC chain multiply one factor each (twelve give
    1.570e-06, transport still converged), while a capped shock inside a DONOR CYCLE
    goes to exactly 0.0 with converged=False, because pointer doubling squares the
    accumulated product every round. M6 medium read sigma_min exactly 0.0, so it was
    the cycle route -- hence "the G8.2 donor-cycle signature" -- with 57 floored
    cells and |R| stalled at 2.49e-06. M6 coarse is healthy only because its Newton
    states have no limited cells at all (its PICARD seed does: 595 of them).

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
        if not lim[e]:
            continue                      # own speed non-physical: no correction
        m2e = mach_number_squared(q2[e], m_inf, gamma)
        if m2e >= 1.0:
            continue                      # still supersonic: not post-shock
        if not lim[u]:
            continue                      # donor speed non-physical (m_cap)
        m2u = mach_number_squared(q2[u], m_inf, gamma)
        if m2u <= 1.0:
            continue                      # donor subsonic: no shock crossed
        m_cur = np.sqrt(m2u)
        rise_max = 0.0
        c = u
        for _ in range(max_walk):
            nxt = upstream[c]
            if nxt == c:
                break                     # chain root
            if not lim[nxt]:
                break                     # limited: stop before reading the cap
            m2n = mach_number_squared(q2[nxt], m_inf, gamma)
            if m2n <= 1.0:
                break                     # front of the supersonic pocket
            m_n = np.sqrt(m2n)
            rise = m_n - m_cur
            if rise <= 0.0:
                break                     # past the local maximum
            if rise < knee_frac * rise_max:
                break                     # the KNEE: out of the shock structure
            if rise > rise_max:
                rise_max = rise
            m_cur = m_n
            c = nxt
        m1 = m_cur
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
    #: knee fraction for the pre-shock walk (see shock_factor_sweep). 0.3 is the
    #: measured middle of an insensitive band; sensitivity is reported by
    #: bench/s1_duct/run_entropy_knee.py.
    KNEE_FRAC_DEFAULT = 0.3

    def __init__(self, n_elements: int, n_round: int = N_ROUND_DEFAULT,
                 max_walk: int = MAX_WALK_DEFAULT,
                 knee_frac: float = KNEE_FRAC_DEFAULT):
        self.n_round = int(n_round)
        self.max_walk = int(max_walk)
        self.knee_frac = float(knee_frac)
        self._s = np.ones(n_elements, dtype=np.float64)
        self._m1 = np.zeros(n_elements, dtype=np.float64)
        self._sigma = np.ones(n_elements, dtype=np.float64)
        self._anc_a = np.zeros(n_elements, dtype=np.int64)
        self._anc_b = np.zeros(n_elements, dtype=np.int64)
        self._prod_b = np.ones(n_elements, dtype=np.float64)
        #: the inert all-true limiter mask for `sigma(lim=None)`. Same reasoning
        #: as UpwindOperator._sigma_ones: numba specialises on the argument TYPE,
        #: so passing None would force a second compilation of the hot sweep.
        self._lim_ones = np.ones(n_elements, dtype=np.bool_)
        self.n_shock = 0
        self.m1_max = 0.0
        self.sigma_min = 1.0
        self.n_rounds = 0
        self.converged = True

    def sigma(self, q2: np.ndarray, upstream: np.ndarray, m_inf: float,
              gamma: float = GAMMA, lim: Optional[np.ndarray] = None
              ) -> np.ndarray:
        """Entropy factor per element (view into the workspace buffer -- copy it
        if it must outlive the next call, which the frozen-sigma callers do).

        `lim` is the m_cap limiter mask (True = this element's q2 was NOT limited),
        i.e. NewtonWorkspace.eval_residual's `lim = q2l == q2n`. Passing it is what
        keeps the correction off cells whose speed the system has already declared
        non-physical -- see shock_factor_sweep. `lim=None` means "nothing was
        limited" and reproduces the pre-20260731-2000 behaviour bit for bit; it is
        the control leg of that round's A/B, and callers with a real limiter mask
        should pass it.
        """
        up = np.ascontiguousarray(upstream, dtype=np.int64)
        n = len(self._s)
        if len(up) != n:
            raise ValueError(
                f"upstream map has {len(up)} entries, expected {n}")
        # ★ Guard, not paranoia: both kernels index arrays with these values, so
        # an out-of-range entry is a SEGMENTATION FAULT, not an exception. That
        # happened (GS1b.3, 2026-07-29): the Picard driver passed
        # UpwindOperator._upstream before any walk had filled it -- `np.empty`
        # int64 garbage. Checking two extrema per call is free next to the sweeps.
        if n and (up.min() < 0 or up.max() >= n):
            raise ValueError(
                f"upstream map out of range [{up.min()}, {up.max()}] for "
                f"{n} elements -- it was probably read before a walk filled it "
                f"(use UpwindOperator.upstream_map(grad))")
        if lim is None:
            lm = self._lim_ones
        else:
            lm = np.ascontiguousarray(lim, dtype=np.bool_)
            if len(lm) != n:
                raise ValueError(
                    f"lim mask has {len(lm)} entries, expected {n}")
        shock_factor_sweep(q2, up, m_inf, gamma, self.max_walk,
                           self.knee_frac, lm, self._s, self._m1)
        self.n_rounds = transport_sigma(
            self._s, up, self.n_round, self._sigma, self._anc_a, self._anc_b,
            self._prod_b)
        self.converged = self.n_rounds < self.n_round
        self.n_shock = int(np.count_nonzero(self._m1 > 0.0))
        self.m1_max = float(self._m1.max()) if len(self._m1) else 0.0
        self.sigma_min = float(self._sigma.min()) if len(self._sigma) else 1.0
        return self._sigma
