"""Two-dimensional integral boundary-layer closure by **fitted correlations**
(Drela-Giles), for the chordwise strip core only.

GS4.1 round 3, pre-registration
`docs/dev_phase_four/20260819-0500-gs41-a2-correlation-closure-prereg.md`.

★★ AUTHORITY, fixed before this file was written (pre-registration section 2):

  - `closures.py` -- profile family, 6-state, carries crossflow -- is the SOLE
    authority for the three-dimensional IBL (`ibl3.py`). Frozen research branch.
  - **this module** -- correlations, 2-state, no crossflow -- is the SOLE
    authority for the two-dimensional strip core.

They are NOT two implementations of one model: different state spaces, different
dimensionality, different provenance. Never compare them as if they were. And
neither is derived from the other -- **this module must not import `closures.py`
and vice versa** (asserted in the round's guard G-AUTHORITY). What keeps them
from silently drifting apart is that BOTH are checked against the same external
oracle, the Falkner-Skan ODE, which lives inside neither of them.

★ What this closure is, in one sentence: the correlations are fits to the
Falkner-Skan similarity family, so Blasius is a fixed point of the resulting
system BY CONSTRUCTION -- measured, the two source terms cancel to 2.5e-4. That
is the exact mirror of round 1's finding that the profile family cannot hold
Blasius and drifts to H = 2.7083, and it is the whole content of route (a2).

★★ Consequence, and the reason a pass here proves less than round 1's failure:
verifying these correlations against Blasius or Falkner-Skan checks a fit
against its own training data. Treat such a check as a TRANSCRIPTION test, not
as validation.

State: `(theta, H)`. Marched form (incompressible edge, u_e prescribed):

    dtheta/dxi = c_f/2 - (H + 2) (theta/u_e) du_e/dxi
    dH/dxi     = [2 c_D - H* c_f/2 + H*(H-1)(theta/u_e) du_e/dxi]
                 / (theta * dH*/dH)

★★ The direct (prescribed-u_e) form is SINGULAR at laminar separation, because
`dH*/dH` vanishes at H = 4 (measured -0.00038 at H = 3.99). That is a known
property of the direct method and the motivation for the quasi-simultaneous
coupling of GS4.2 -- it is not a defect to be patched here. `march_state`
refuses to step past `H_SEPARATION_GUARD` rather than propping the equation up
with an artificial limiter.

LAMINAR ONLY. The turbulent correlations and the shear-stress lag equation are a
second model and a separate round; the strip core's turbulent branch continues
to use the profile closure. Do not read this module as covering turbulent flow.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Correlation constants -- the laminar closure of
#
#   Drela & Giles, "Viscous-Inviscid Analysis of Transonic and Low Reynolds
#   Number Airfoils", AIAA Journal 25(10), Eqs. (10), (11), (12).
#
# ★ CITATION ADDED 2026-08-19 (GS4.1 round 5 dry run). Round 3 wrote these from
# memory and verified them against an independent Falkner-Skan ODE rather than
# against a source, because the source was not in the tree; the verdict recorded
# their provenance as unattributed. With the paper available they check out
# VERBATIM against Eqs. (10)-(12), so the code was right and only the citation
# was missing.
#
# ★★ They are the PUBLISHED forms, not the ones in XFOIL 6.99. XFOIL's
# xblsys.f carries later code revisions of the same closures -- HSL is built on
# (Hk - 4.35), CFL on 0.0727*(5.5-Hk)^3/(Hk+1), and HST is marked "new
# correlation 29 Nov 91" with the form used here sitting commented out above it.
# Measured against the Falkner-Skan ODE at Blasius, the published forms here are
# the more accurate ones: c_f*Re_theta is -0.08 % against XFOIL CFL's -2.89 %,
# and in adverse gradient (m = -0.05) +0.53 % against -7.35 %. So this is not a
# stale copy -- do not "update" it to the XFOIL code forms without measuring.
#
# DIL is the one that agrees between both: the 0.00205 / 5.5 / 0.207 group here
# is identical in XFOIL's DIL.
# ---------------------------------------------------------------------------

H_KINK = 4.0            # branch split of all three correlations (= separation)
HS_BASE = 1.515
HS_A_LO, HS_A_HI = 0.076, 0.040
CF_BASE, CF_A_LO = -0.067, 0.01977
CF_H_SPLIT, CF_A_HI = 7.4, 0.022
CD_BASE, CD_A_LO, CD_P_LO = 0.207, 0.00205, 5.5
CD_A_HI, CD_B_HI = 0.0016, 0.02

#: Refuse to march past this. Below H = 4 the direct form is well posed; at 4 it
#: is singular. 3.9 leaves dH*/dH = -0.004, already four decades below its
#: Blasius value, so a leg reaching here is reported rather than continued.
H_SEPARATION_GUARD = 3.9
#: The (H - 1) denominator in the c_f correlation.
H_MIN = 1.05


class ClosureRangeError(ValueError):
    """The march left the correlation's valid range -- report, do not clamp."""


def h_star(H):
    """Kinetic-energy shape parameter `H* = theta*/theta`."""
    if H < H_KINK:
        return HS_BASE + HS_A_LO * (H_KINK - H) ** 2 / H
    return HS_BASE + HS_A_HI * (H - H_KINK) ** 2 / H


def dh_star_dh(H):
    """`dH*/dH`. Vanishes at H = 4 -- see the module docstring."""
    if H < H_KINK:
        return -HS_A_LO * (H_KINK - H) * (H + H_KINK) / H ** 2
    return HS_A_HI * (H - H_KINK) * (H + H_KINK) / H ** 2


def re_theta_cf_half(H):
    """`Re_theta * c_f/2`."""
    if H < CF_H_SPLIT:
        return CF_BASE + CF_A_LO * (CF_H_SPLIT - H) ** 2 / (H - 1.0)
    return CF_BASE + CF_A_HI * (1.0 - 1.4 / (H - 6.0)) ** 2


def re_theta_2cd_over_hstar(H):
    """`Re_theta * 2 c_D / H*`."""
    if H < H_KINK:
        return CD_BASE + CD_A_LO * (H_KINK - H) ** CD_P_LO
    d = H - H_KINK
    return CD_BASE - CD_A_HI * d ** 2 / (1.0 + CD_B_HI * d ** 2)


def packet(theta, H, ue, rho=1.0, mu=1.0e-5):
    """Closure quantities at one station.

    Returns a dict rather than a 30-slot array: this closure has five outputs,
    and naming them costs nothing at the call rate the strip core uses.
    """
    if not np.isfinite(H) or H <= H_MIN:
        raise ClosureRangeError(
            f"H = {H!r} at or below the correlation's floor {H_MIN} -- the "
            "(H-1) denominator diverges; reporting rather than clamping")
    re_theta = max(rho * ue * theta / mu, 1.0e-12)
    hs = h_star(H)
    return {
        "re_theta": re_theta,
        "H_star": hs,
        "dH_star_dH": dh_star_dh(H),
        "cf": 2.0 * re_theta_cf_half(H) / re_theta,
        "cD": re_theta_2cd_over_hstar(H) * hs / (2.0 * re_theta),
    }


def rhs(theta, H, ue, due, rho=1.0, mu=1.0e-5):
    """`(dtheta/dxi, dH/dxi)` -- the two-equation system, explicit.

    No global system, no quadrature, no state Jacobian: this is why the route
    exists. Compare `closures.py`'s per-call cost, measured at 7.82 us.
    """
    p = packet(theta, H, ue, rho=rho, mu=mu)
    grad = theta / ue * due
    dtheta = 0.5 * p["cf"] - (H + 2.0) * grad
    num = (2.0 * p["cD"] - p["H_star"] * 0.5 * p["cf"]
           + p["H_star"] * (H - 1.0) * grad)
    den = theta * p["dH_star_dH"]
    if den == 0.0:
        raise ClosureRangeError(
            f"dH*/dH vanished at H = {H!r}: the DIRECT two-equation form is "
            "singular at separation by construction (GS4.2's motivation), not "
            "a defect to patch here")
    return dtheta, num / den


def zpg_fixed_point():
    """The zero-pressure-gradient self-similar `H` of THIS closure family.

    Under the similarity ansatz the two source terms must cancel, and both
    carry one factor `1/Re_theta`, so the condition collapses to the
    Re-independent algebraic statement

        `Re_theta 2c_D/H*  ==  Re_theta c_f/2`

    -- no march, no discretization, no Reynolds number. It is the direct
    counterpart of `strip2d.similarity_fixed_point` for the profile family, so
    each closure's ZPG prediction can be compared against the SAME external
    truth (Blasius). ★ That is comparing each family to the oracle, which is the
    sanctioned pattern; it is NOT treating the two as interchangeable
    implementations of one model, which the round forbids.
    """
    from scipy.optimize import brentq

    def resid(H):
        return re_theta_2cd_over_hstar(H) - re_theta_cf_half(H)

    return brentq(resid, 1.5, 3.9, xtol=1e-14)


def blasius_state(x, ue=1.0, rho=1.0, mu=1.0e-5, H=2.591100):
    """`(theta, H)` at station `x` on a Blasius plate.

    `theta sqrt(Re_x)/x = 0.664115` is the similarity constant; the caller may
    override `H` with a value it computed itself (the round's script passes its
    own ODE solution, so no reference value is taken on trust here).
    """
    re_x = rho * ue * x / mu
    return 0.664115 * x / np.sqrt(re_x), H
