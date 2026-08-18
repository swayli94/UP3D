"""XFOIL's two-point implicit BL discretisation, turbulent intervals only.

Binding text: docs/dev_phase_four/20260821-0100-two-point-scheme-prereg.md and
addendum #1 (20260821-0130, turbulent-only + consistency on the real u_e).

★★ This is an INSTRUMENT, not a product path. It exists to turn round 10's loose
bound on the discretisation into a reading, and its registered expiry is round
11's verdict. It lives in bench/ and NOT in pyfp3d/ precisely so that it cannot
quietly become a second marcher nobody chose.

It shares `closures_2d` with `strip2d.march_correlation` and shares no
discretisation code with it -- that is what makes the consistency check an
independent oracle rather than a tautology.

Source: xblsys.f BLDIF (1554-1979) + BLMID (1128-1180), at M = 0 with no wall
blowing, so MA = HCA = HWA = 0.
"""

import numpy as np
from scipy.optimize import root

from pyfp3d.viscous import closures_2d as C2

HUPWT = 1.0          # xblsys.f:1590
HDCON_WALL = 5.0     # xblsys.f:1591  HDCON = 5*HUPWT/HK2**2
HLSQ_CAP = 15.0      # xblsys.f:1598


def pack(T, D, S, u, rho, mu):
    """Every per-station quantity BLDIF reads, from OUR closure."""
    H = D / T
    HK = C2.h_kinematic(H, 0.0)                      # identity at M = 0
    RT = rho * u * T / mu
    HS = C2.h_star_turb(H, RT)
    US = C2.slip_velocity(HS, H)
    CF = C2.cf_turb(H, RT)                           # max(CFT, CFL) -- momentum
    CFw = C2.cf_turb_wall(H, RT)                     # raw CFT -- dissipation
    cD = C2.cd_turb(CFw, US, S * S, HS, C2.dfac_low_hk(H, RT), RT)
    return {"T": T, "D": D, "S": S, "U": u, "H": H, "HK": HK, "RT": RT,
            "HS": HS, "US": US, "CF": CF,
            "DI": C2.dissipation_identity(cD, HS),
            "CQ": np.sqrt(C2.ctau_eq(HS, H, US, RT)),
            "DE": C2.bl_thickness(T, H)}


def residuals(p1, p2, X1, X2):
    """REZT / REZH / REZC for one turbulent interval. `xblsys.f` line numbers in
    the comments; `M = 0` and no blowing, so MA = HCA = HWA = 0."""
    HK1, HK2 = p1["HK"], p2["HK"]

    # --- upwinding weight, xblsys.f:1590-1600 --------------------------------
    HDCON = HDCON_WALL * HUPWT / HK2 ** 2
    ARG = abs((HK2 - 1.0) / (HK1 - 1.0))
    HL = np.log(max(ARG, 1.0e-300))
    HLSQ = min(HL * HL, HLSQ_CAP)
    UPW = 1.0 - 0.5 * np.exp(-HLSQ * HDCON)

    XLOG = np.log(X2 / X1)
    ULOG = np.log(p2["U"] / p1["U"])
    TLOG = np.log(p2["T"] / p1["T"])
    HLOG = np.log(p2["HS"] / p1["HS"])

    HA = 0.5 * (p1["H"] + p2["H"])
    XA = 0.5 * (X1 + X2)
    TA = 0.5 * (p1["T"] + p2["T"])
    XOT1, XOT2 = X1 / p1["T"], X2 / p2["T"]

    # --- momentum, xblsys.f:1876-1886 ---------------------------------------
    # ★ CFM is BLMID's midpoint c_f at the ARITHMETIC mean state (xblsys.f:1155)
    #   -- NOT the upwind-weighted HKA used by the lag block below. Addendum #1
    #   named this in advance because BLDIF defines HKA twice.
    HKA_mid = 0.5 * (HK1 + HK2)
    RTA_mid = 0.5 * (p1["RT"] + p2["RT"])
    CFM = C2.cf_turb(HKA_mid, RTA_mid)               # BLMID takes max(CFT, CFL)
    CFX_T = 0.5 * CFM * XA / TA + 0.25 * (p1["CF"] * XOT1 + p2["CF"] * XOT2)
    REZT = TLOG + (HA + 2.0) * ULOG - XLOG * 0.5 * CFX_T

    # --- kinetic energy, xblsys.f:1930-1943 ---------------------------------
    HSA = 0.5 * (p1["HS"] + p2["HS"])
    DIX = (1.0 - UPW) * p1["DI"] * XOT1 + UPW * p2["DI"] * XOT2
    CFX_H = (1.0 - UPW) * p1["CF"] * XOT1 + UPW * p2["CF"] * XOT2
    REZH = HLOG + (1.0 - HA) * ULOG + XLOG * (0.5 * CFX_H - DIX)

    # --- lag, xblsys.f:1690-1771 --------------------------------------------
    SA = (1.0 - UPW) * p1["S"] + UPW * p2["S"]
    CQA = (1.0 - UPW) * p1["CQ"] + UPW * p2["CQ"]
    CFA = (1.0 - UPW) * p1["CF"] + UPW * p2["CF"]
    HKA = (1.0 - UPW) * HK1 + UPW * HK2              # ★ upwind-weighted, :1693
    USA = 0.5 * (p1["US"] + p2["US"])
    RTA = 0.5 * (p1["RT"] + p2["RT"])
    DEA = 0.5 * (p1["DE"] + p2["DE"])
    DA = 0.5 * (p1["D"] + p2["D"])
    ALD = 1.0                                        # wall layer, :1705
    HKC = max(HKA - 1.0 - C2.GCCON / RTA, 0.01)
    HR = HKC / (C2.GACON * ALD * HKA)
    UQ = (0.5 * CFA - HR * HR) / (C2.GBCON * DA)
    SCC = C2.SCCON * 1.333 / (1.0 + USA)
    SLOG = np.log(p2["S"] / p1["S"])
    DXI = X2 - X1
    REZC = (SCC * (CQA - SA * ALD) * DXI
            - DEA * 2.0 * SLOG
            + DEA * 2.0 * (UQ * DXI - ULOG) * C2.DUXCON)
    return np.array([REZT, REZH, REZC]), UPW


def step(y1, X1, X2, u1, u2, rho, mu, guess=None):
    """Solve the three residuals for (T2, D2, S2). Newton via scipy's hybrid
    method on a numerical Jacobian -- BLDIF's analytic derivatives are not
    transcribed, and are not needed for a 3x3."""
    p1 = pack(*y1, u1, rho, mu)
    g = np.array(guess if guess is not None else y1, dtype=float)

    def F(z):
        T2, D2, S2 = np.exp(z)                       # positivity by construction
        return residuals(p1, pack(T2, D2, S2, u2, rho, mu), X1, X2)[0]

    sol = root(F, np.log(g), method="hybr", tol=1e-13)
    r = float(np.max(np.abs(sol.fun)))
    # ★ Addendum #3: accept on the RESIDUAL, not on sol.success. scipy's hybr
    # reports failure when its xtol is unattainable at this scale even though
    # |F| is at machine zero -- one leg threw with |F| = 9.6e-16. Throwing away
    # a converged solve in the reporting layer is a family this project has paid
    # for before.
    return np.exp(sol.x), r <= 1.0e-12, r


def march(stations, y0, x0, ue_of, rho, mu):
    """March station to station. `stations` are the X (arc length) values and
    `ue_of(x)` returns u_e; both come from XFOIL, never from an interpolant of
    ours in X-ANSWER."""
    y = np.array(y0, dtype=float)
    X = float(x0)
    out = {k: [] for k in ("x", "theta", "H", "cf", "ctau", "re_theta", "upw",
                           "resid")}
    for X2 in np.asarray(stations, dtype=float):
        y2, ok, r = step(y, X, X2, ue_of(X), ue_of(X2), rho, mu, guess=y)
        if not ok:
            raise RuntimeError(f"two-point step failed at X = {X2:.6f}, "
                               f"|F| = {r:.3e}")
        p = pack(*y2, ue_of(X2), rho, mu)
        _, upw = residuals(pack(*y, ue_of(X), rho, mu), p, X, X2)
        out["x"].append(X2); out["theta"].append(y2[0]); out["H"].append(p["H"])
        out["cf"].append(p["CF"]); out["ctau"].append(y2[2] ** 2)
        out["re_theta"].append(p["RT"]); out["upw"].append(upw)
        out["resid"].append(r)
        y, X = y2, X2
    return {k: np.asarray(v) for k, v in out.items()}
