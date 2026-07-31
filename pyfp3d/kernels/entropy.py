"""Entropy-corrected (non-isentropic) full-potential density factor.

Phase-two gates GS1b.3 (first delivery) and GS1b.10 (this construction).
Pre-registrations: docs/dev_phase_two/20260729-0700-s1b-entropy-implementation.md
and 20260731-0000-s1b-chainfree-sigma.md.

WHY THIS EXISTS (measured, not assumed). At a fixed transonic condition the computed
lift does not converge under refinement, and the only shock-weakening mechanism in the
code -- the artificial density -- is ~O(h): at the finest level upwind_c = 1.5 and 3.0
give nearly the same answer (GS1b.2 Q5). A mechanism that vanishes with h cannot cure
an error that does not, so the weakening has to come from the physics. Isentropic full
potential over-predicts the density jump across a shock by 5.1 / 6.9 / 8.9 % at
M1 = 1.30 / 1.35 / 1.40, exactly the strengths the airfoil cases run at (GS1b.2 Q1).

THE RELATION. For a perfect gas whose total enthalpy is constant along streamlines,
rho0 = p0/(R*T0) and T0 is preserved across a shock, so

    rho02/rho01 = p02/p01   EXACTLY,   and   rho_s = sigma * rho_isen(q^2)

with sigma = p02/p01 the Rankine-Hugoniot total-pressure ratio. Applying that factor
makes the full-potential mass-conservation jump reproduce Rankine-Hugoniot identically
(verified to 0.0000 % for M1 = 1.15..1.60 in bench/s1_duct/run_entropy_premise.py). It
is NOT sigma^(1/(gamma-1)) -- that form applies the exponent twice; at M1 = 1.35
sigma^1 gives the R-H ratio 1.60278 while sigma^(1/(gamma-1)) gives 1.38162.

HOW sigma IS BUILT, and the three designs it took. The first two hung sigma on the
UPWIND DONOR MAP -- one most-upwind face neighbour per element -- and that map is a
DISCRETE selection: between adjacent Mach steps it changes for 0.06 % (coarse) /
0.34 % (medium) of elements, and since sigma was a PRODUCT along the donor chain, one
flip near the shock re-routed an entire downstream chain. Measured consequences: the
refresh limit-cycled (max|dsigma| pinned at 2.9e-2 with the residual stalled at
~5e-6), sigma_min swung 0.932-0.986 between neighbouring conditions, and the converged
ON answer became RECIPE-DEPENDENT by 0.118 c of shock position at medium, where the
isentropic answer is recipe-independent to four decimals. Smoothing the author's own
hard tests did not help (GS1b.4, S1 FAIL) because the churn is one layer below them,
and driving sigma to a self-consistent fixed point did not converge either (GS1b.9,
P2 FAIL -- the map sigma -> phi -> sigma is not a contraction).

This construction removes the donor chain from BOTH places it was used:

  transport   a flux-weighted first-order FINITE-VOLUME UPWIND (face_inflow_weights
              plus sigma_transport_sweep). The weights max(0, -v.n)A vary CONTINUOUSLY
              as the velocity rotates, so there is no selection to flip. One
              Gauss-Seidel sweep is EXACT because elements are visited in order of
              their own phi, a topological order of the upwind operator (dphi/ds > 0
              along a streamline).
  production  no pre-shock Mach is identified at all. Shock entropy is ADDITIVE along
              the path, so each cell contributes the production of its own compression
              via g(M) = d(-ln sigma_RH)/dM, and the product across a smeared shock
              telescopes to sigma_RH(M1) exactly (GS1b.10 Q1: N = 1..16 cells reproduce
              it to +0.49..0.004 %, improving with smearing). That also retires the two
              opposite-direction bugs the walk produced -- reading inside the smeared
              structure (too weak) and reading the pocket peak (2x too strong) -- and
              all four of its knobs (knee_frac, eps_m, f_lo, f_hi), with no new
              softening width needed because g(1) = 0 to machine precision.

sigma is FROZEN over a Newton step by the caller: unlike the flux's upwind selection,
which is piecewise constant in phi, sigma depends continuously on phi, so a live sigma
with no d(sigma)/d(phi) term would make the Jacobian genuinely inexact rather than
exact-almost-everywhere.
"""

import os

import numba
import numpy as np

from pyfp3d.mesh.metrics import precompute_face_normals
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
    """Rankine-Hugoniot total-pressure ratio p02/p01 for a normal shock at pre-shock
    Mach m1. Returns 1.0 at or below m1 = 1 (no shock, no entropy)."""
    if m1 <= 1.0:
        return 1.0
    m2 = m1 * m1
    a = ((gamma + 1.0) * m2) / ((gamma - 1.0) * m2 + 2.0)
    b = (gamma + 1.0) / (2.0 * gamma * m2 - (gamma - 1.0))
    return a ** (gamma / (gamma - 1.0)) * b ** (1.0 / (gamma - 1.0))


@_njit(cache=True, fastmath=True)
def entropy_production_density(m: float, gamma: float = GAMMA) -> float:
    """g(M) = d(-ln sigma_RH)/dM -- entropy production per unit Mach of compression,
    analytic from the R-H total-pressure ratio (agrees with a central difference of
    log(total_pressure_ratio) to 1e-9).

    This is what lets the construction avoid identifying a pre-shock Mach: the integral
    of g from 1 to M1 IS -ln sigma_RH(M1), so distributing the production cell by cell
    telescopes to the correct total however the shock is smeared. And g(1) = 0 to
    machine precision (the R-H entropy jump vanishes like (M-1)^3), so the sonic line
    needs no softening width.
    """
    if m <= 1.0:
        return 0.0
    u = m * m
    dln_a = 1.0 / u - (gamma - 1.0) / ((gamma - 1.0) * u + 2.0)
    dln_b = -2.0 * gamma / (2.0 * gamma * u - (gamma - 1.0))
    return -2.0 * m * (gamma / (gamma - 1.0) * dln_a + dln_b / (gamma - 1.0))


@_njit(cache=True, fastmath=True, parallel=True)
def face_inflow_weights(
    grad: np.ndarray,
    face_normals: np.ndarray,
    face_areas: np.ndarray,
    w_out: np.ndarray,
) -> None:
    """w[e, f] = max(0, -v_e . n_ef) * A_ef, the inflow weight of face f of element e.

    Continuous in the velocity direction by construction -- that is the point. The
    previous designs picked the single most-upwind neighbour, a discrete choice whose
    flips were the measured root of every sigma instability (module docstring).
    """
    n = len(grad)
    for e in prange(n):
        for f in range(4):
            d = -(grad[e, 0] * face_normals[e, f, 0]
                  + grad[e, 1] * face_normals[e, f, 1]
                  + grad[e, 2] * face_normals[e, f, 2])
            w_out[e, f] = d * face_areas[e, f] if d > 0.0 else 0.0


@_njit(cache=True)
def sigma_transport_sweep(
    q2: np.ndarray,
    weights: np.ndarray,
    face_neighbors: np.ndarray,
    order: np.ndarray,
    m_inf: float,
    gamma: float,
    sigma_out: np.ndarray,
    m_up_out: np.ndarray,
    d1_out: np.ndarray,
    chi_out: np.ndarray,
) -> None:
    """One Gauss-Seidel sweep of the FV upwind entropy transport, in flow order:

        sigma_e = (sum_in w sigma_nbr / sum_in w) * exp(-g(Mbar) * dM_supersonic)

    ONE sweep is exact: `order` sorts elements by their own phi, and dphi/ds > 0 along
    a streamline, so every inflow neighbour is visited before the element receiving
    from it. Serial and in a fixed order, hence bit-reproducible -- which the project's
    bit-identity A/Bs require (phase-one discipline #12).

    Domain-boundary inflow (face_neighbors < 0) carries sigma = 1 and the element's own
    Mach, i.e. contributes no entropy: correct for farfield inflow, and harmless at a
    wall where v is tangential so the inflow weight is ~0 anyway.
    """
    n = len(q2)
    for k in range(n):
        e = order[k]
        m_e = np.sqrt(mach_number_squared(q2[e], m_inf, gamma))
        w_sum = 0.0
        sig_in = 0.0
        m_in = 0.0
        d1_in = 0.0
        for f in range(4):
            w = weights[e, f]
            if w <= 0.0:
                continue
            nb = face_neighbors[e, f]
            w_sum += w
            if nb < 0:
                sig_in += w
                m_in += w * m_e
                # boundary inflow: no upstream compression to compare against
            else:
                sig_in += w * sigma_out[nb]
                m_in += w * np.sqrt(mach_number_squared(q2[nb], m_inf, gamma))
                d1_in += w * d1_out[nb]
        if w_sum <= 0.0:
            sigma_out[e] = 1.0
            m_up_out[e] = m_e
            d1_out[e] = 0.0
            chi_out[e] = 0.0
            continue
        m_up = m_in / w_sum
        m_up_out[e] = m_up
        d1 = m_up - m_e                      # this cell's compression
        d1_out[e] = d1
        d2 = d1_in / w_sum                   # the compression one cell upstream
        # ★ SHOCK DISCRIMINATOR (GS1b.10, candidate (i)): entropy is produced by a
        # SHOCK, not by a smooth supersonic compression -- and the two are NOT
        # distinguishable locally, only by RATE. Without this gate the construction
        # charges the airfoil pocket's 0.3-chord isentropic deceleration (M 1.41 ->
        # 1.29) and over-produces: measured 1-sigma 6.22 % where the strongest shock
        # the field can hold allows 3.57 %.
        #
        # chi = max(0, 1 - d2/d1) compares this cell's compression with the previous
        # cell's: equal rates (a smooth compression) give chi = 0, a rate that jumps
        # (a shock front) gives chi -> 1. Zero-knob and dimensionless, reusing the FV
        # weights one extra time. Verified on synthetic streamlines (smooth 1.41 ->
        # 1.29 then a shock to 0.95, over 10/30/60 smooth and 1/2/3 shock cells):
        # the total reproduces sigma_RH(1.29) to -0.01..+0.50 %, against -2.35..-2.56 %
        # with no gate. ★ My own reasoning had predicted this would FAIL inside the
        # shock structure (where consecutive drops are comparable, so chi -> 0); the
        # synthetic measurement refuted that prediction, which is why it was measured
        # rather than argued.
        chi = 1.0 - d2 / d1 if d1 > 0.0 else 0.0
        if chi < 0.0:
            chi = 0.0
        elif chi > 1.0:
            chi = 1.0
        chi_out[e] = chi
        m_lo = m_e if m_e > 1.0 else 1.0
        d_m = m_up - m_lo
        if d_m > 0.0 and chi > 0.0:
            m_bar = 0.5 * (m_up + m_lo)
            s_e = np.exp(-chi * entropy_production_density(m_bar, gamma) * d_m)
        else:
            s_e = 1.0
        sigma_out[e] = (sig_in / w_sum) * s_e


def _face_areas(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """(n_tets, 4) face areas, face f opposite local node f -- the same local ordering
    build_face_adjacency and precompute_face_normals use."""
    p = nodes[elements]                                  # (n, 4, 3)
    out = np.empty((len(elements), 4), dtype=np.float64)
    for f in range(4):
        idx = [i for i in range(4) if i != f]
        a = p[:, idx[1]] - p[:, idx[0]]
        b = p[:, idx[2]] - p[:, idx[0]]
        out[:, f] = 0.5 * np.linalg.norm(np.cross(a, b), axis=1)
    return out


class EntropyOperator:
    """Per-mesh workspace for the entropy factor (allocates once, geometry once).

    Usage inside a driver, once per Newton step (sigma frozen over the step):

        ent = EntropyOperator(nodes, elements, face_neighbors)
        sigma = ent.sigma(q2, grad, phi_cut, m_inf, gamma)   # view into the buffer
        rho_s = sigma * rho_isentropic                       # then feed rho_tilde

    Monitors after each call: `n_shock` (cells carrying a non-negligible factor),
    `m1_max` (the strongest inflow-weighted upwind Mach where entropy was produced),
    `sigma_min`, and `converged` -- always True here, because the FV transport has no
    donor cycle to fail on, unlike the chain product it replaces. The flag is kept so
    the drivers' clamp-not-silent wiring stays unchanged.
    """

    def __init__(self, nodes: np.ndarray, elements: np.ndarray,
                 face_neighbors: np.ndarray):
        self._fn = np.ascontiguousarray(face_neighbors)
        self._normals = precompute_face_normals(
            np.ascontiguousarray(nodes, dtype=np.float64),
            np.ascontiguousarray(elements))
        self._areas = _face_areas(np.asarray(nodes, dtype=np.float64),
                                  np.asarray(elements))
        n = len(elements)
        self._w = np.zeros((n, 4), dtype=np.float64)
        self._sigma = np.ones(n, dtype=np.float64)
        self._m_up = np.zeros(n, dtype=np.float64)
        self._d1 = np.zeros(n, dtype=np.float64)
        self._chi = np.zeros(n, dtype=np.float64)
        self._elements = np.ascontiguousarray(elements)
        self.n_shock = 0
        self.m1_max = 0.0
        self.sigma_min = 1.0
        self.converged = True

    def sigma(self, q2: np.ndarray, grad: np.ndarray, phi_cut: np.ndarray,
              m_inf: float, gamma: float = GAMMA) -> np.ndarray:
        """Entropy factor per element (view into the workspace buffer -- copy it if it
        must outlive the next call, which the frozen-sigma callers do).

        `phi_cut` supplies the flow ordering: elements are visited in order of their
        own mean phi, a topological order of the upwind operator because
        dphi/ds = |grad phi| > 0 along a streamline. That is what makes one sweep exact
        instead of needing hundreds of Jacobi passes.
        """
        face_inflow_weights(np.ascontiguousarray(grad, dtype=np.float64),
                            self._normals, self._areas, self._w)
        phi_e = np.asarray(phi_cut, dtype=np.float64)[self._elements].mean(axis=1)
        order = np.argsort(phi_e, kind="stable").astype(np.int64)
        self._sigma[:] = 1.0
        self._d1[:] = 0.0
        self._chi[:] = 0.0
        sigma_transport_sweep(np.ascontiguousarray(q2, dtype=np.float64),
                              self._w, self._fn, order, m_inf, gamma,
                              self._sigma, self._m_up, self._d1, self._chi)
        prod = self._sigma < 1.0 - 1e-9
        self.n_shock = int(np.count_nonzero(prod))
        self.m1_max = float(self._m_up[prod].max()) if prod.any() else 0.0
        self.sigma_min = float(self._sigma.min())
        self.converged = True
        return self._sigma
