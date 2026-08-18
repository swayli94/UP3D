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

# ---------------------------------------------------------------------------
# TURBULENT closure (GS4.1 round 5). Every constant carries its source, because
# the previous attempt at these was written from memory and drove a flat plate
# to an unphysical H = 0.60.
#
# Sources:
#   XFOIL 6.99, Xfoil699src/src/xblsys.f and xbl.f  (reference/XFOIL6.99.zip)
#   Drela & Giles, AIAA J 25(10)                    (reference/drela-giles-*.pdf)
#
# ★★ LOCAL EQUILIBRIUM ONLY: Ctau is set to CtauEQ everywhere. The shear-stress
# LAG equation is deliberately absent -- a zero-pressure-gradient plate cannot
# test it (equilibrium IS Ctau = CtauEQ, so the lag term vanishes identically),
# and shipping code this round cannot verify is what the scope boundary exists
# to prevent. The gate's "two-equation + lag" is met by rounds 5 and 6 together.
# ---------------------------------------------------------------------------

GACON = 6.70                       # xbl.f:1559
GBCON = 0.75                       # xbl.f:1560
GCCON = 18.0                       # xbl.f:1561; USED on the turbulent WALL
CTCON = 0.5 / (GACON ** 2 * GBCON)  # xbl.f:1569  -> 0.01485...
HSMIN = 1.500                      # xblsys.f:2394 DATA HSMIN
DHSINF = 0.015                     # xblsys.f:2394 DATA DHSINF

#: Physical band a turbulent flat plate must stay inside. Leaving it is reported,
#: never clamped -- the memory attempt left it and that is how the bug surfaced.
H_TURB_LO, H_TURB_HI = 1.05, 4.0


def h_kinematic(H, mach=0.0):
    """Whitfield kinematic shape parameter. `xblsys.f:2278` HKIN / D-G eq (9).

    At M = 0 this is the identity, which is the only regime this module runs in;
    it is written out so the incompressible assumption is visible rather than
    implied.
    """
    msq = mach * mach
    return (H - 0.29 * msq) / (1.0 + 0.113 * msq)


def h_star_turb(H, re_theta, mach=0.0, arm="new"):
    """Turbulent `H*`. `xblsys.f:2388` HST.

    `arm="new"` is the correlation XFOIL 6.99 actually runs, marked there
    "new correlation 29 Nov 91". `arm="old"` is the form commented out directly
    above it, which is the same generation as this module's laminar set. Round 4
    measured the published laminar forms MORE accurate than XFOIL's later code
    forms, so which generation is better on the turbulent side cannot be settled
    by analogy -- both arms are provided and the round measures them.
    """
    hk = h_kinematic(H, mach)
    ho = 3.0 + 400.0 / re_theta if re_theta > 400.0 else 4.0
    rtz = re_theta if re_theta > 200.0 else 200.0
    if arm == "old":
        hs = 1.505 + 4.0 / re_theta
        if hk < ho:
            hs += (0.165 - 1.6 / np.sqrt(re_theta)) * (ho - hk) ** 1.6 / hk
        return hs
    if hk < ho:                                   # attached branch
        hr = (ho - hk) / (ho - 1.0)
        hs = (2.0 - HSMIN - 4.0 / rtz) * hr ** 2 * 1.5 / (hk + 0.5) \
            + HSMIN + 4.0 / rtz
    else:                                         # separated branch
        grt = np.log(rtz)
        hdif = hk - ho
        rtmp = hdif + 4.0 / grt
        hs = hdif ** 2 * (0.007 * grt / rtmp ** 2 + DHSINF / hk) \
            + HSMIN + 4.0 / rtz
    fm = 1.0 + 0.014 * mach * mach                # Whitfield compressibility
    return (hs + 0.028 * mach * mach) / fm


def cf_turb(H, re_theta, mach=0.0):
    """Turbulent `c_f` (Coles). `xblsys.f:2483` CFT / D-G eq (13) (Swafford)."""
    hk = h_kinematic(H, mach)
    fc = np.sqrt(1.0 + 0.5 * 0.4 * mach * mach)
    grt = max(np.log(re_theta / fc), 3.0)
    gex = -1.74 - 0.31 * hk
    arg = max(-20.0, -1.33 * hk)
    cfo = 0.3 * np.exp(arg) * (grt / 2.3026) ** gex
    return (cfo + 1.1e-4 * (np.tanh(4.0 - hk / 0.875) - 1.0)) / fc


def slip_velocity(hs, H, mach=0.0):
    """Normalized slip velocity `Us`. `xblsys.f:825`, `GBCON` from `xbl.f:1560`."""
    hk = h_kinematic(H, mach)
    return 0.5 * hs * (1.0 - (hk - 1.0) / (GBCON * H))


def ctau_eq(hs, H, us, re_theta, mach=0.0):
    """Equilibrium shear coefficient `CtauEQ`. `xblsys.f:856-877` (BLVAR).

    XFOIL stores its square root (`CQ2 = SQRT(...)`); this returns `CtauEQ`
    itself, which is what the dissipation relation consumes.

        HKB = Hk - 1
        HKC = Hk - 1 - GCCON/Re_theta        (floored at 0.01)
        CtauEQ = CTCON * H* * HKB * HKC^2 / ((1 - Us) * H * Hk^2)

    ★★ HKB and HKC are TWO DIFFERENT quantities. Round 5 wrote (Hk-1)**3,
    collapsing them and dropping the -GCCON/Re_theta correction, because the
    subtraction sits inside XFOIL's `IF(ITYP.EQ.2)` branch and I read ITYP=2 as
    the wake. The source says two lines above (`xblsys.f:794-795`) that
    ITYP = 1 is laminar and ITYP = 2 is TURBULENT, and `xbl.f:810` calls
    BLVAR(2) for every station past transition. So the correction applies on
    the turbulent wall -- exactly where this function is used.

    The omission made CtauEQ high by 14 % at Re_theta 578, 4.7 % at 2000 and
    1.1 % at 1e4, and it feeds c_D. Fixed in round 8; see
    docs/dev_phase_four/20260820-0500-gcc-fix-prereg.md.
    """
    hk = h_kinematic(H, mach)
    hkb = hk - 1.0
    hkc = max(hk - 1.0 - GCCON / re_theta, 0.01)
    return CTCON * hs * hkb * hkc * hkc / ((1.0 - us) * H * hk * hk)


def cd_turb(cf, us, ctau, hs):
    """Turbulent `c_D`. `xblsys.f:2375` DIT.

    ★★ DIT returns `DI = 2 c_D / H*`, NOT `c_D`. Reading it as `c_D` puts a
    spurious factor of `H*/2` into the dissipation and drove the previous
    attempt's flat plate to H = 0.60. `dissipation_identity` below is the
    machine check for exactly that (guard G-DI).
    """
    return 0.5 * cf * us + ctau * (1.0 - us)


def dissipation_identity(cd, hs):
    """`DI = 2 c_D / H*` -- returned so a caller can assert the relation that
    the memory attempt got wrong."""
    return 2.0 * cd / hs


def packet_turb(theta, H, ue, rho=1.0, mu=1.0e-5, mach=0.0, arm="new"):
    """Turbulent closure quantities at one station, at local equilibrium."""
    if not np.isfinite(H) or not (H_TURB_LO <= H <= H_TURB_HI):
        raise ClosureRangeError(
            f"turbulent H = {H!r} outside the physical band "
            f"[{H_TURB_LO}, {H_TURB_HI}] -- reporting rather than clamping "
            "(the memory attempt reached H = 0.60 here)")
    re_theta = max(rho * ue * theta / mu, 1.0e-12)
    hs = h_star_turb(H, re_theta, mach, arm)
    cf = cf_turb(H, re_theta, mach)
    us = slip_velocity(hs, H, mach)
    ct = ctau_eq(hs, H, us, re_theta, mach)
    return {"re_theta": re_theta, "H_star": hs, "cf": cf, "Us": us,
            "Ctau_eq": ct, "cD": cd_turb(cf, us, ct, hs)}


def rhs_turb(theta, H, ue, due, rho=1.0, mu=1.0e-5, mach=0.0, arm="new",
             rel_step=1.0e-6):
    """`(dtheta/dxi, dH/dxi)` on the turbulent branch.

    ★ One structural difference from the laminar branch: `H*` depends on
    `Re_theta` as well as `H` (`HSL` sets `HS_RT = 0`, `HST` does not), so the
    kinetic-energy equation carries an extra term,

        theta (dH*/dH H' + dH*/dRe_theta Re_theta') = 2 c_D - H* c_f/2
                                                      + H*(H-1) theta u_e'/u_e

    with `Re_theta' = (rho/mu)(u_e theta' + theta u_e')`. Dropping it would be a
    silent modelling change, not a simplification.

    The two `H*` partials are taken by central difference. That is a deliberate
    choice: they enter a marching right-hand side rather than a Newton Jacobian,
    so difference accuracy is ample, and transcribing HST's analytic
    derivatives would add a large surface of hand-copied algebra to a round
    whose whole point is that hand-copying is what failed last time.
    """
    p = packet_turb(theta, H, ue, rho=rho, mu=mu, mach=mach, arm=arm)
    hs, ret = p["H_star"], p["re_theta"]
    grad = theta / ue * due

    dtheta = 0.5 * p["cf"] - (H + 2.0) * grad

    dh = rel_step * max(abs(H), 1.0)
    dhs_dH = (h_star_turb(H + dh, ret, mach, arm)
              - h_star_turb(H - dh, ret, mach, arm)) / (2.0 * dh)
    dr = rel_step * ret
    dhs_dre = (h_star_turb(H, ret + dr, mach, arm)
               - h_star_turb(H, ret - dr, mach, arm)) / (2.0 * dr)

    dre = rho / mu * (ue * dtheta + theta * due)
    num = (2.0 * p["cD"] - hs * 0.5 * p["cf"] + hs * (H - 1.0) * grad
           - theta * dhs_dre * dre)
    den = theta * dhs_dH
    if den == 0.0:
        raise ClosureRangeError(
            f"dH*/dH vanished at H = {H!r} on the turbulent branch")
    return dtheta, num / den


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
