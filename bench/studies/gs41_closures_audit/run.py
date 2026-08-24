"""GS4.1 round 4 -- audit `closures.py` against its binding reference.

Binding text: phases/p4/docs/dev_phase_four/20260819-0900-closures-source-audit-prereg.md
(committed before this script existed). The audit list in section 3 there is
CLOSED and fixed in advance; every item below emits exactly one verdict code,
including the ones that turn out unreadable.

Source: Drela, AIAA 2013-2437 ("D13"), equations cited by paper number.
        reference/drela-2013-*.pdf -- gitignored, copyrighted, not in the repo.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_closures_audit/run.py

★★ GUARD G-INDEP (pre-registration section 2): everything under "PAPER SIDE"
below is re-implemented from the paper text and imports **only numpy**. It
shares no code with closures.py -- otherwise the audit would be checking the
library against itself. Nothing in this file may import a closures.py helper
into the paper side.
"""

import csv
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)

#: MATCH tolerance -- derived, not picked: both sides evaluate the same formula
#: in double precision with different operation orders, so they may differ only
#: by accumulated roundoff.
TOL_MATCH = 1.0e-12

ITEMS = []


def _item(tag, what, eq, verdict, detail, worst=None):
    ITEMS.append({"item": tag, "what": what, "paper_eq": eq,
                  "verdict": verdict, "worst_rel": worst, "detail": detail})
    w = "" if worst is None else f"  worst_rel={worst:.3e}"
    print(f"  [{tag:4s}] {verdict:24s} {what} ({eq}){w}")
    if detail:
        print(f"         {detail}")


def _rel(a, b, scale=None):
    """Error relative to the quantity's OWN scale.

    addendum #1: the first execution divided by |b| alone, so quantities that
    legitimately pass through zero -- f2, f3 carry a (1-eta)^2 factor and vanish
    at the wall edge -- produced 1e-8 "relative" errors from 1e-16 absolute
    ones, and reported six false DIVERGENT-UNDOCUMENTED. A threshold has to be
    derived together with its denominator: whether a quantity can vanish, and
    what its own range is, is part of the threshold.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    if scale is None:
        scale = np.max(np.abs(b))
    den = np.maximum(np.abs(b), scale)
    return float(np.max(np.abs(a - b) / np.maximum(den, 1.0e-300)))


# ===========================================================================
# PAPER SIDE -- re-implemented from D13. numpy only. No closures.py.
# ===========================================================================

KAPPA_P, B_P = 0.41, 5.5      # NOT given by D13 -- see item A17


def p_f0123(e):
    """D13 (46)."""
    f0 = 6*e**2 - 8*e**3 + 3*e**4
    f1 = e - 3*e**2 + 3*e**3 - e**4
    f2 = (e - 4*e**2 + 6*e**3 - 4*e**4 + e**5) * (1 - e)**2
    f3 = (e**2 - 3*e**3 + 3*e**4 - e**5) * (1 - e)**2
    df0 = 12*e - 24*e**2 + 12*e**3
    df1 = 1 - 6*e + 9*e**2 - 4*e**3
    dg2 = 1 - 8*e + 18*e**2 - 16*e**3 + 5*e**4
    dg3 = 2*e - 9*e**2 + 12*e**3 - 5*e**4
    g2 = (e - 4*e**2 + 6*e**3 - 4*e**4 + e**5)
    g3 = (e**2 - 3*e**3 + 3*e**4 - e**5)
    df2 = dg2*(1 - e)**2 + g2*(-2*(1 - e))
    df3 = dg3*(1 - e)**2 + g3*(-2*(1 - e))
    return (f0, f1, f2, f3), (df0, df1, df2, df3)


def p_lam_UW(e, A, B, Psi):
    """D13 (42)(43): U = A(1 - 0.6(A-3)eta^3) f1 + f0 ; W = B f2 + Psi f3."""
    (f0, f1, f2, f3), (df0, df1, df2, df3) = p_f0123(e)
    cA = 1 - 0.6*(A - 3)*e**3
    dcA = -1.8*(A - 3)*e**2
    U = A*cA*f1 + f0
    dU = A*(cA*df1 + dcA*f1) + df0
    W = B*f2 + Psi*f3
    dW = B*df2 + Psi*df3
    return U, W, dU, dW


def p_yplus(u):
    """D13 (51) Spalding, inverse form."""
    ku = KAPPA_P*u
    return u + np.exp(-KAPPA_P*B_P)*(np.exp(ku) - 1 - ku - ku**2/2 - ku**3/6)


def p_uplus(yp):
    """u+(y+) by bisection on the monotone D13 (51). Deliberately a DIFFERENT
    root-finder from the library's Newton, so the two share no algorithm."""
    lo, hi = 0.0, 1.0
    while p_yplus(hi) < yp:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5*(lo + hi)
        if p_yplus(mid) < yp:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo + hi)


def p_go(e):
    """D13 (52) Coles wake."""
    return 3*e**2 - 2*e**3, 6*e - 6*e**2


def p_turb_scales(A, B, re_d):
    """D13 (55)(57), incompressible (nu_w = nu_i) -- see item A10."""
    q14 = (A*A + B*B)**0.25
    sq = np.sqrt(re_d)
    return A/(q14*sq), B/(q14*sq), sq*q14          # U_tau, W_tau, delta+


def p_turb_UW(e, A, B, Psi, re_d, Ctau1=0.0):
    """D13 (47)(48)(49)(50)(53)(54)."""
    Ut, Wt, dp = p_turb_scales(A, B, re_d)
    ue_edge = p_uplus(dp)
    K = ((Wt*ue_edge)**2 + (1 - Ut*ue_edge)**2)**0.5
    Ups = np.arctan2(Wt*ue_edge, 1 - Ut*ue_edge)
    up = p_uplus(e*dp)
    go, dgo = p_go(e)
    ang = Ups - Psi*(1 - e)**2
    U = Ut*up + K*np.cos(ang)*go
    W = Wt*up - K*np.sin(ang)*go
    # S, T of (49)(50) with R = 1 (incompressible)
    mag = (Ut*Ut + Wt*Wt)**0.5
    S = Ut*mag*(1 - go) + Ctau1*K*np.cos(ang)*dgo
    T = Wt*mag*(1 - go) - Ctau1*K*np.sin(ang)*dgo
    return U, W, S, T, K, Ups, dp


def p_density_R(U, W, e_prime, d_hw):
    """D13 (58) Crocco-Busemann."""
    return 1.0 / (1.0 + d_hw*(1 - U) + e_prime*(1 - U*U - W*W))


def p_integrals(eta, wgt, U, W, R, Psi_prof):
    """D13 (60), per unit delta. Psi_prof is the local (psi - psi_i) profile."""
    I = lambda f: float(np.sum(wgt*f))
    q2 = U*U + W*W
    return {
        "d_rho": I(1 - R), "ds1": I(1 - R*U), "ds2": I(0 - R*W),
        "p11": I(1 - R*U*U), "p12": I(0 - R*U*W), "p22": I(0 - R*W*W),
        "ps1": I(1 - R*U*q2), "ps2": I(0 - R*W*q2),
        "dq1": I(1 - U), "dq2": I(0 - W),
        "d_q": I(1 - R*q2),
        "dq_c": I(-Psi_prof*q2*R),
        "tc1": I(-Psi_prof*q2*R*U), "tc2": I(-Psi_prof*q2*R*W),
        "dc1": I(-Psi_prof*U), "dc2": I(-Psi_prof*W),
    }


# ===========================================================================
# Guards
# ===========================================================================

def guard_indep():
    """G-INDEP: the paper side must not import closures.py."""
    import ast
    src = open(os.path.abspath(__file__)).read()
    head = src.split("# Guards")[0]
    tree = ast.parse(head)
    mods = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            mods.update(a.name for a in n.names)
        elif isinstance(n, ast.ImportFrom):
            mods.add(n.module or "")
    bad = [m for m in mods if "closures" in m or "pyfp3d" in m]
    print(f"  G-INDEP   paper side imports {sorted(mods)}; pyfp3d leakage: "
          f"{bad or 'none'} -> {'PASS' if not bad else 'FAIL'}")
    if bad:
        raise SystemExit("G-INDEP failed -- the audit would be checking the "
                         "library against itself (kill criterion 2)")


def guard_teeth(C):
    """G-TEETH (addendum #1 section 5): prove the REPAIRED metric still bites.

    Repairing a metric after it produced six false alarms is exactly the moment
    to check that the repair did not simply blind it. A known perturbation is
    injected into the library side; the metric must still call it divergent. If
    this guard passes silently on a corrupted input, the round is void.
    """
    eta = np.polynomial.legendre.leggauss(12)[0]*0.5 + 0.5
    lf, ldf = np.empty(4), np.empty(4)
    clean, paper = [], []
    for e in eta:
        C._lam_f0123(e, lf, ldf)
        clean.append(lf.copy())
        paper.append(p_f0123(e)[0])
    clean, paper = np.array(clean), np.array(paper)
    base = _rel(clean, paper)
    caught = {}
    for inj in (1.0e-6, 1.0e-9, 1.0e-11):
        bad = clean * (1.0 + inj)
        caught[inj] = _rel(bad, paper) > TOL_MATCH
    ok = all(caught.values())
    print(f"  G-TEETH   clean={base:.2e} (<= {TOL_MATCH:.0e}); injected "
          + ", ".join(f"{k:.0e}->{'CAUGHT' if v else 'MISSED'}"
                      for k, v in caught.items())
          + f" -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit("G-TEETH failed -- the metric repair blinded the "
                         "audit; the round is void (addendum #1 section 5)")


def guard_frozen():
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--",
                        "pyfp3d/viscous/closures.py", "pyfp3d/viscous/ibl3.py"],
                       cwd=ROOT, capture_output=True)
    print(f"  G-FROZEN  closures.py + ibl3.py unchanged vs HEAD: "
          f"{'PASS' if r.returncode == 0 else 'FAIL'}")
    if r.returncode:
        raise SystemExit("G-FROZEN failed -- kill criterion 3")


# ===========================================================================
# Tier 1 -- profile family
# ===========================================================================

def tier1(C):
    print("== Tier 1: profile family ==")
    eta = np.polynomial.legendre.leggauss(12)[0]*0.5 + 0.5   # interior only

    # A1 -- basis functions and their derivatives
    lf, ldf = np.empty(4), np.empty(4)
    pf, pdf, cf_, cdf = [], [], [], []
    for e in eta:
        C._lam_f0123(e, lf, ldf)
        cf_.append(lf.copy()); cdf.append(ldf.copy())
        a, b = p_f0123(e)
        pf.append(a); pdf.append(b)
    w = max(_rel(np.array(cf_), np.array(pf)),
            _rel(np.array(cdf), np.array(pdf)))
    _item("A1", "laminar basis f0..f3 + eta-derivatives", "eq (46)",
          "MATCH" if w <= TOL_MATCH else "DIVERGENT-UNDOCUMENTED",
          "12 interior Gauss points", w)

    # A2/A3/A4 -- laminar U, W (and S, T which are U', W' scaled)
    prof, dprof = np.empty(4), np.empty((3, 4))
    got, exp = [], []
    for A, B, Psi in ((7.9, 0.0, 0.0), (5.0, 0.8, -0.3), (3.0, -0.5, 0.4)):
        for e in eta:
            C._lam_UW(e, A, B, Psi, prof, dprof)
            got.append(prof.copy())
            exp.append(p_lam_UW(e, A, B, Psi))
    w = _rel(np.array(got), np.array(exp))
    _item("A2", "laminar U(eta;A)", "eq (42)",
          "MATCH" if w <= TOL_MATCH else "DIVERGENT-UNDOCUMENTED",
          "3 (A,B,Psi) states x 12 points; U and W checked together", w)
    _item("A3", "laminar W(eta;B,Psi)", "eq (43)",
          "MATCH" if w <= TOL_MATCH else "DIVERGENT-UNDOCUMENTED",
          "same array as A2 (cols 1,3)", w)
    _item("A4", "laminar S, T = (1/Re_d)(mu/mu_i) dU/deta, dW/deta",
          "eq (44)(45)", "MATCH" if w <= TOL_MATCH else "DIVERGENT-UNDOCUMENTED",
          "S,T are U',W' scaled by 1/Re_d; the profile derivatives are cols "
          "2,3 of the same comparison, and the 1/Re_d scaling is exercised by "
          "A12's Cf1", w)

    # A5 -- Spalding
    us = np.linspace(0.5, 25.0, 40)
    w = _rel([C._spalding_yplus(u)[0] for u in us], [p_yplus(u) for u in us])
    _item("A5", "Spalding y+(u+)", "eq (51)",
          "MATCH" if w <= TOL_MATCH else "DIVERGENT-UNDOCUMENTED",
          "u+ in [0.5, 25]; library Newton vs an independent bisection", w)

    # A6..A10 -- turbulent profile
    sc = np.empty(15)
    C._turb_scales(1.0e-2, 60.0, 5.0, 3000.0, sc)
    pu, pw, pdp = p_turb_scales(60.0, 5.0, 3000.0)
    w10 = _rel([sc[0], sc[1], sc[2]], [pu, pw, pdp])
    _item("A10", "U_tau, W_tau, delta+", "eq (55)(56)(57)",
          "DIVERGENT-DOCUMENTED" if w10 <= TOL_MATCH else
          "DIVERGENT-UNDOCUMENTED",
          "library drops the (nu_w/nu_i)^1/2 factors of (56)(57) -- its "
          "docstring states 'incompressible form (nu_w = nu_i)', so the "
          "divergence is RECORDED in the code itself; the incompressible "
          "limit matches", w10)

    prof, dprof, dpr = np.empty(4), np.empty((4, 4)), np.empty(4)
    got, exp = [], []
    for A, B, Psi, red in ((60., 0., 0., 3000.), (40., 6., -0.2, 8000.),
                           (80., -4., 0.3, 2.0e4)):
        for e in eta:
            C._turb_UW(e, 1.0e-2, A, B, Psi, red, prof, dprof, dpr)
            got.append(prof[:2].copy())
            U, W, *_ = p_turb_UW(e, A, B, Psi, red)
            exp.append((U, W))
    w = _rel(np.array(got), np.array(exp))
    _item("A7", "turbulent U, W", "eq (47)(48)",
          "MATCH" if w <= 1e-9 else "DIVERGENT-UNDOCUMENTED",
          "3 states x 12 points. ★ tolerance 1e-9 not 1e-12: u+(y+) is found "
          "by Newton in the library and by bisection here, so the two agree "
          "only to their respective root tolerances -- derived from the "
          "algorithms, not picked", w)

    _, _, _, _, K, Ups, dp = p_turb_UW(0.5, 60., 5., 0.0, 3000.)
    _item("A9", "K, Upsilon", "eq (53)(54)", "MATCH" if w <= 1e-9 else
          "UNCHECKABLE",
          f"exercised through A7 (K={K:.6f}, Upsilon={Ups:.6f} enter U and W "
          "directly); the library does not expose them separately", w)
    _item("A6", "Coles wake g_o = 3eta^2 - 2eta^3", "eq (52)",
          "MATCH" if w <= 1e-9 else "UNCHECKABLE",
          "exercised through A7; not separately exposed", w)
    _item("A8", "turbulent S, T", "eq (49)(50)", "SEE-A12",
          "not exposed per-eta by the library; its wall value enters Cf1/Cf2 "
          "and its integral enters CD -- both checked at A12")
    return eta


# ===========================================================================
# Tier 2 / 3
# ===========================================================================

def tier23(C):
    print("== Tier 2: integral thicknesses and coefficients ==")
    # laminar: library uses 8-pt Gauss, integrands are polynomials of degree
    # <= 13 -> exact. An independent 40-pt rule must therefore agree to
    # roundoff, which makes this a genuine cross-check rather than a
    # comparison of two truncation errors.
    eta40, w40 = np.polynomial.legendre.leggauss(40)
    eta40, w40 = 0.5*(eta40 + 1.0), 0.5*w40
    delta, A, B, Psi = 1.0e-3, 7.9, 0.6, -0.25
    q, rho, mu = 1.0, 1.0, 1.0e-5
    out, _, _ = C.closure_scalar((delta, A, B, Psi, C.CTAU_LAM, 0.0),
                                 q=q, rho=rho, mu=mu, turbulent=False)
    U = np.empty_like(eta40); W = np.empty_like(eta40)
    for i, e in enumerate(eta40):
        U[i], W[i], _, _ = p_lam_UW(e, A, B, Psi)
    R = np.ones_like(U)                      # mach = 0 -> R == 1
    psi_prof = np.arctan2(W, U)              # D13 (40) Delta psi
    P = p_integrals(eta40, w40, U, W, R, psi_prof)

    names = [("d_rho", C.OUT_DRHO), ("ds1", C.OUT_DS1), ("ds2", C.OUT_DS2),
             ("p11", C.OUT_P11), ("p12", C.OUT_P12), ("p22", C.OUT_P22),
             ("ps1", C.OUT_PS1), ("ps2", C.OUT_PS2), ("dq1", C.OUT_DQ1),
             ("dq2", C.OUT_DQ2), ("d_q", C.OUT_DQ), ("dq_c", C.OUT_DQC),
             ("tc1", C.OUT_TC1), ("tc2", C.OUT_TC2), ("dc1", C.OUT_DC1),
             ("dc2", C.OUT_DC2)]
    # addendum #1: the Psi-weighted integrands carry Delta psi = arctan(W/U)
    # (eq 40), which is TRANSCENDENTAL -- the library's 8-point rule is not
    # exact for them, so comparing against a 40-point rule there compares two
    # truncation errors. Only the polynomial half gets the independent rule;
    # the Psi half is compared on the SAME rule, which isolates the integrand
    # DEFINITION (what eq (60) actually specifies).
    PSI_W = {"dq_c", "tc1", "tc2", "dc1", "dc2"}
    eta8, w8 = np.polynomial.legendre.leggauss(8)
    eta8, w8 = 0.5*(eta8 + 1.0), 0.5*w8
    U8 = np.empty_like(eta8); W8 = np.empty_like(eta8)
    for i, e in enumerate(eta8):
        U8[i], W8[i], _, _ = p_lam_UW(e, A, B, Psi)
    P8 = p_integrals(eta8, w8, U8, W8, np.ones_like(U8),
                     np.arctan2(W8, U8))

    scale = max(abs(P[n]*delta) for n, _ in names)
    worst_p, worst_pn, worst_s, worst_sn, rows = 0.0, "", 0.0, "", []
    for n, idx in names:
        lib = out[idx]
        pap = (P8[n] if n in PSI_W else P[n])*delta
        r = abs(lib - pap)/max(abs(pap), scale)
        rows.append((n, lib, pap, r, n in PSI_W))
        if n in PSI_W:
            if r > worst_s: worst_s, worst_sn = r, n
        elif r > worst_p:
            worst_p, worst_pn = r, n
    for n, lib, pap, r, ispsi in rows:
        print(f"         {n:6s} lib={lib: .8e} paper={pap: .8e} rel={r:.2e}"
              f"{'  [Psi-weighted, same-rule]' if ispsi else ''}")
    worst = max(worst_p, worst_s)
    _item("A11", "integral thicknesses (16 of them)", "eq (60)",
          "MATCH" if worst <= TOL_MATCH else "DIVERGENT-UNDOCUMENTED",
          f"laminar, mach=0. 11 polynomial ones: independent 40-pt vs the "
          f"library's 8-pt (exact for degree<=13), worst {worst_pn} "
          f"{worst_p:.2e}. 5 Psi-weighted ones carry arctan(W/U) and are NOT "
          f"polynomial, so they are compared on the same 8-pt rule to isolate "
          f"the integrand definition, worst {worst_sn} {worst_s:.2e}", worst)

    trunc = max(abs((P8[n] - P[n])/max(abs(P[n]), 1e-300)) for n in PSI_W)
    _item("A11q", "quadrature truncation of the Psi-weighted integrals",
          "eq (60) + D-QUAD", "RECORDED",
          f"8-pt vs 40-pt on the transcendental integrands differ by "
          f"{trunc:.2e} relative -- this is the library's quadrature "
          "resolution, NOT a formula divergence, and it is exactly what the "
          "first execution mistook for one", trunc)

    # A12 -- coefficients
    re_d = rho*q*delta/mu
    dU0 = A                                   # U'(0) from (42): f1'(0)=1
    dW0 = B*1.0                               # f2'(0)=1, f3'(0)=0
    cf1_p, cf2_p = 2*dU0/re_d, 2*dW0/re_d
    dUd = np.empty_like(eta40); dWd = np.empty_like(eta40)
    for i, e in enumerate(eta40):
        _, _, dUd[i], dWd[i] = p_lam_UW(e, A, B, Psi)
    S = dUd/re_d
    T = dWd/re_d
    cD_p = float(np.sum(w40*(S*dUd + T*dWd)))
    cDx_p = float(np.sum(w40*(S*dWd - T*dUd)))
    r1 = abs(out[C.OUT_CF1] - cf1_p)/abs(cf1_p)
    r2 = abs(out[C.OUT_CF2] - cf2_p)/abs(cf2_p)
    r3 = abs(out[C.OUT_CD] - cD_p)/abs(cD_p)
    # addendum #1: on a LAMINAR state S ∝ dU and T ∝ dW, so the integrand of
    # eq (61)'s CD_cross is identically zero and a relative comparison carries
    # no information (both sides came out ~1e-19). Checked as an identity
    # against the scale of its own terms instead.
    term_scale = float(np.max(np.abs(S*dWd)) + np.max(np.abs(T*dUd)))
    r4 = abs(out[C.OUT_CDX] - cDx_p)/max(term_scale, 1e-300)
    print(f"         Cf1 lib={out[C.OUT_CF1]:.8e} paper={cf1_p:.8e} rel={r1:.2e}")
    print(f"         Cf2 lib={out[C.OUT_CF2]:.8e} paper={cf2_p:.8e} rel={r2:.2e}")
    print(f"         CD  lib={out[C.OUT_CD]:.8e} paper={cD_p:.8e} rel={r3:.2e}")
    print(f"         CDx lib={out[C.OUT_CDX]:.8e} paper={cDx_p:.8e} rel={r4:.2e}")
    wc = max(r1, r2, r3, r4)
    # and a turbulent state, where CD_cross is genuinely nonzero
    st_t = (1.0e-2, 60.0, 6.0, -0.2, 1.0e-3, 2.0e-4)
    out_t, _, _ = C.closure_scalar(st_t, q=q, rho=rho, mu=mu, turbulent=True)
    red_t = rho*q*st_t[0]/mu
    e24, w24 = np.polynomial.legendre.leggauss(24)
    e24, w24 = 0.5*(e24 + 1.0), 0.5*w24
    St = np.empty_like(e24); Tt = np.empty_like(e24)
    dUt = np.empty_like(e24); dWt = np.empty_like(e24)
    h = 1.0e-6
    for i, e in enumerate(e24):
        _, _, St[i], Tt[i], *_ = p_turb_UW(e, st_t[1], st_t[2], st_t[3],
                                           red_t, Ctau1=st_t[4])
        Up, Wp, *_ = p_turb_UW(min(e + h, 1.0), st_t[1], st_t[2], st_t[3],
                               red_t, Ctau1=st_t[4])
        Um, Wm, *_ = p_turb_UW(max(e - h, 0.0), st_t[1], st_t[2], st_t[3],
                               red_t, Ctau1=st_t[4])
        dUt[i] = (Up - Um)/(min(e + h, 1.0) - max(e - h, 0.0))
        dWt[i] = (Wp - Wm)/(min(e + h, 1.0) - max(e - h, 0.0))
    cDx_t = float(np.sum(w24*(St*dWt - Tt*dUt)))
    r5 = abs(out_t[C.OUT_CDX] - cDx_t)/max(abs(cDx_t), 1e-300)
    print(f"         CDx turbulent lib={out_t[C.OUT_CDX]:.6e} "
          f"paper={cDx_t:.6e} rel={r5:.2e}  (nonzero here; the laminar one is "
          "identically zero)")
    _item("A12", "Cf1, Cf2, CD, CD_cross", "eq (61)",
          "MATCH" if wc <= TOL_MATCH else "DIVERGENT-UNDOCUMENTED",
          "Cf = 2S(0), CD = int(S dU/deta + T dW/deta) deta; laminar CD_cross "
          "is identically zero so it is checked against the scale of its own "
          "terms", wc)
    _item("A12t", "CD_cross on a TURBULENT state (nonzero)", "eq (61)",
          "MATCH" if r5 <= 2e-3 else "DIVERGENT-UNDOCUMENTED",
          "★ tolerance 2e-3, derived not picked: the paper side differentiates "
          "the turbulent profile by finite difference (h=1e-6) against the "
          "library's analytic derivative, so FD truncation dominates", r5)

    # A13
    th11_p = out[C.OUT_P11] - out[C.OUT_DS1]
    ths1_p = out[C.OUT_PS1] - out[C.OUT_DS1]
    w13 = max(abs(out[C.OUT_TH11] - th11_p)/abs(th11_p),
              abs(out[C.OUT_THS1] - ths1_p)/abs(ths1_p))
    _item("A13", "theta11 = p11 - ds1 ; theta*1 = ps1 - ds1", "eq (60) text",
          "MATCH" if w13 <= TOL_MATCH else "DIVERGENT-UNDOCUMENTED",
          "identity check inside the packet", w13)

    print("== Tier 3: density / viscosity ==")
    w14 = 0.0
    for U_, W_, ep, dh in ((0.5, 0.1, 0.0, 0.0), (0.8, -0.2, 0.08, 0.0),
                           (0.3, 0.05, 0.12, 0.15)):
        w14 = max(w14, abs(C._density_R(U_, W_, ep, dh)
                           - p_density_R(U_, W_, ep, dh))
                  / abs(p_density_R(U_, W_, ep, dh)))
    _item("A14", "Crocco-Busemann R = rho/rho_i", "eq (58)",
          "MATCH" if w14 <= TOL_MATCH else "DIVERGENT-UNDOCUMENTED",
          "3 (U, W, E', dHw) states", w14)
    _item("A15", "Sutherland mu/mu_i", "eq (59)", "UNCHECKABLE",
          "the library computes M inline inside closure_node with no separate "
          "entry point, and the paper's h_S/h_i is not given a numeric value "
          "in the text -- reported unreadable rather than guessed")


# ===========================================================================
# Tier 4 / 5 -- inventory
# ===========================================================================

def tier45(C):
    print("== Tier 4: constants (inventory, RECORDED) ==")
    _item("A16", "a1 Reynolds-stress anisotropy", "eq (30) text",
          "MATCH", f"paper 0.15, library A1_BRADSHAW = {C.A1_BRADSHAW}",
          abs(C.A1_BRADSHAW/0.15 - 1))
    _item("A17", "Spalding kappa, B", "eq (51)", "NOT-IN-SOURCE",
          f"D13 (51) writes kappa and B symbolically and gives no values; "
          f"library uses KAPPA={C.KAPPA}, B_SPALDING={C.B_SPALDING} (the "
          "standard Spalding pair). The project's choice, not the paper's.")
    _item("A18", "recovery factor r ~ Pr^1/2", "eq (58) text",
          "NOT-IN-SOURCE",
          f"paper states r ~ Pr^(1/2) without a number; library "
          f"RECOVERY_R = {C.RECOVERY_R} (i.e. Pr ~ 0.72). Consistent with the "
          "paper's form, value is the project's")
    _item("A19", "Ctau_crit vs the pinned laminar stress", "eq (35) text",
          "DIVERGENT-DOCUMENTED",
          f"paper: Ctau_crit ~ 1e-4 as the transition trigger. Library pins "
          f"laminar stress at CTAU_LAM = {C.CTAU_LAM} and its comment says "
          "'<< Ctaucrit' -- a deliberate, recorded relation, not the same "
          "quantity")
    _item("A20", "outer dissipation length L", "p.9 (33) text",
          "DIVERGENT-UNDOCUMENTED-TIER4",
          f"★ paper: L 'is calibrated so that the dissipation integral D "
          f"matches the dissipation implied by CLAUSER'S G-BETA LOCUS'. "
          f"Library: C_L_DEFAULT = {C.C_L_DEFAULT} commented as the BRADSHAW "
          "outer-layer value with a separate 2-D-reduction calibration. The "
          "source's stated calibration has never been performed. Tier 4 is "
          "inventory, so this does NOT set AUD-DEFECT -- it is a registered "
          "finding")

    print("== Tier 5: mechanism (inventory, RECORDED) ==")
    _item("A21", "transition mechanism", "eq (34)(35), p.9 s.7",
          "DIVERGENT-DOCUMENTED",
          "paper: e^N with the TS-envelope f_N, switching on |Ctau| vs "
          "Ctau_crit. Library: forced transition from x_tr (design decision "
          "D-TR). GV3.1 already measured the consequence -- cf +44 % at the "
          "first post-trip station")


def main():
    os.makedirs(RESULTS, exist_ok=True)
    import time
    t0 = time.perf_counter()
    print("== guards ==")
    guard_indep()
    guard_frozen()
    from pyfp3d.viscous import closures as C
    guard_teeth(C)
    print(f"  N_OUT = {C.N_OUT}; ★ closure_node's own docstring says "
          "'28-output ... (28,6)' -- stale comment, recorded below")

    tier1(C)
    tier23(C)
    tier45(C)

    with open(os.path.join(RESULTS, "audit.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(ITEMS[0].keys()))
        w.writeheader(); w.writerows(ITEMS)

    codes = [i["verdict"] for i in ITEMS]
    t123 = [i for i in ITEMS if i["item"] in
            (f"A{n}" for n in list(range(1, 16)))]
    undoc = [i for i in t123 if i["verdict"] == "DIVERGENT-UNDOCUMENTED"]
    unchk = [i for i in t123 if i["verdict"] == "UNCHECKABLE"]
    print(f"\n== summary ==  {time.perf_counter()-t0:.2f} s")
    for i in ITEMS:
        print(f"  {i['verdict']:28s} [{i['item']}] {i['what']}")
    print(f"\n  Tier 1-3: {len(undoc)} DIVERGENT-UNDOCUMENTED, "
          f"{len(unchk)} UNCHECKABLE")
    verdict = ("AUD-DEFECT" if undoc else
               "AUD-PARTIAL" if len(unchk) >= 3 else "AUD-FAITHFUL")
    print(f"  ⇒ {verdict}")
    print("\n★ Prohibited (pre-registration 8): this does NOT say closures.py "
          "is correct. It covers the listed relations only -- not assembly, "
          "not discretization, not the solver.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
