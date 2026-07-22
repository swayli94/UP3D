"""Drela IBL3 closure relations: velocity/stress profile families, integral
thicknesses, skin-friction and dissipation coefficients, and the shear-stress
transport closure, all with analytic derivatives w.r.t. the six primary
unknowns (delta, A, B, Psi, Ctau1, Ctau2).

Binding reference: Drela, AIAA 2013-2437 ("D13"), equations cited by their
paper numbers. Implementation decisions for the points D13 leaves
under-specified are recorded in docs/design_track_v.md §3–§4 (D-CT, D-QUAD,
D-TR, D-STRESS-2); the binding text is the design doc, this module implements
it.

Structure per design.md §7: pure scalar physics as njit'd functions
(unit-testable in isolation, PYFP3D_NOJIT=1 debuggable); the per-node closure
packet is evaluated with a Gauss-Legendre eta-quadrature engine shared by
values and derivatives (analytic differentiation under the integral sign), so
FD checks converge to quadrature accuracy, not to noise.

Nomenclature (D13 §II/§III): state v = (delta, A, B, Psi, Ctau1, Ctau2);
delta = thickness scale, A/B = wall-slope direction (U'(0), W'(0) in viscous
units), Psi = profile twist, Ctau1/Ctau2 = outer stress scales. External
parameters per node: edge speed q, density rho, viscosity mu, Mach (for the
Crocco-Busemann density profile (58); M=0 recovers incompressible R==1).

The 30-output packet (N_OUT=30), indices fixed by OUT_* constants below:
  0..15  integral thicknesses (D13 (60)): d_rho, ds1, ds2, p11, p12, p22,
         ps1, ps2, dq1, dq2, d_q, dq_c, tc1, tc2, dc1, dc2
  16..20 coefficients (D13 (61)):       cf1, cf2, cD, cDx, cDc
  21..24 stress-transport factors (D-CT): a_k, s_p1, s_p2, s_d
  25..27 derived:                       theta11, theta_star1, H1
  28..29 stress-flux factors (D-CT):    ku1, ku2 (K_tau advection integrals)
Derivatives: (30, 6) analytic w.r.t. v (u_e parameters held fixed).

Laminar vs turbulent branch is chosen by the `turbulent` flag (forced
transition in V1: per-node regime from x_tr, design doc §2.4 / D-TR); the
stress-transport factors are evaluated the same way in both regimes (in the
laminar branch the stress equations are PINNED by the solver, so the factors
are computed but unused upstream of x_tr).
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


# ---------------------------------------------------------------------------
# Constants (recorded per evidence discipline; sources: D13 + standard values)
# ---------------------------------------------------------------------------

KAPPA = 0.41          # von Karman constant (Spalding law (51))
B_SPALDING = 5.5      # Spalding log-law intercept
A1_BRADSHAW = 0.15    # Reynolds stress anisotropy ratio a1 (D13 (30))
C_L_DEFAULT = 0.09    # outer dissipation length L = C_L * delta (D-CT-2;
                      # Bradshaw outer-layer value; 2-D-reduction calibration
                      # recorded in the GV1.1 VERDICT)
CTAU_LAM = 1.0e-8     # pinned laminar stress level (D-TR; << Ctaucrit)
GAMMA_AIR = 1.4
RECOVERY_R = 0.85     # r ~ Pr^1/2 (D13 (58))

# Safety floors for Newton/FD probe states that overshoot out of the physical
# branch (never active at a valid solution; derivatives are masked to zero in
# the floored direction so FD-vs-analytic stays consistent there too).
DELTA_MIN = 1.0e-8    # thickness scale floor
RE_D_MIN = 1.0e-3     # Re_delta floor (keeps sqrt/log-free divisions finite)

N_OUT = 30
# packet indices
OUT_DRHO = 0
OUT_DS1 = 1
OUT_DS2 = 2
OUT_P11 = 3
OUT_P12 = 4
OUT_P22 = 5
OUT_PS1 = 6
OUT_PS2 = 7
OUT_DQ1 = 8
OUT_DQ2 = 9
OUT_DQ = 10
OUT_DQC = 11
OUT_TC1 = 12
OUT_TC2 = 13
OUT_DC1 = 14
OUT_DC2 = 15
OUT_CF1 = 16
OUT_CF2 = 17
OUT_CD = 18
OUT_CDX = 19
OUT_CDC = 20
OUT_AK = 21
OUT_SP1 = 22
OUT_SP2 = 23
OUT_SD = 24
OUT_TH11 = 25
OUT_THS1 = 26
OUT_H1 = 27
OUT_KU1 = 28
OUT_KU2 = 29


# ---------------------------------------------------------------------------
# Gauss-Legendre eta-quadrature tables on [0, 1] (D-QUAD)
# ---------------------------------------------------------------------------

def _gauss_table(n):
    x, w = np.polynomial.legendre.leggauss(n)
    eta = 0.5 * (x + 1.0)
    wgt = 0.5 * w
    return np.ascontiguousarray(eta), np.ascontiguousarray(wgt)


# Laminar: integrands are polynomials of degree <= 13 in eta (products of the
# degree-4/5 Bernstein-type basis functions) -> 8 points exact to degree 15.
ETA_LAM, W_LAM = _gauss_table(8)
# Turbulent: Spalding profile is non-polynomial; 24 points resolve the
# near-wall log-region variation for delta+ up to O(1e5) (unit-tested).
ETA_TURB, W_TURB = _gauss_table(24)


# ---------------------------------------------------------------------------
# Laminar profile family (D13 (42)-(46))
# ---------------------------------------------------------------------------

@_njit(cache=True, fastmath=True)
def _lam_f0123(eta, out_f, out_df):
    """Basis functions f0..f3 and their eta-derivatives at one station."""
    e = eta
    e2 = e * e
    e3 = e2 * e
    e4 = e3 * e
    e5 = e4 * e
    f0 = 6.0 * e2 - 8.0 * e3 + 3.0 * e4
    f1 = e - 3.0 * e2 + 3.0 * e3 - e4
    g2 = e - 4.0 * e2 + 6.0 * e3 - 4.0 * e4 + e5
    g3 = e2 - 3.0 * e3 + 3.0 * e4 - e5
    h = (1.0 - e) * (1.0 - e)
    f2 = g2 * h
    f3 = g3 * h
    df0 = 12.0 * e - 24.0 * e2 + 12.0 * e3
    df1 = 1.0 - 6.0 * e + 9.0 * e2 - 4.0 * e3
    dg2 = 1.0 - 8.0 * e + 18.0 * e2 - 16.0 * e3 + 5.0 * e4
    dg3 = 2.0 * e - 9.0 * e2 + 12.0 * e3 - 5.0 * e4
    dh = -2.0 * (1.0 - e)
    df2 = dg2 * h + g2 * dh
    df3 = dg3 * h + g3 * dh
    out_f[0] = f0
    out_f[1] = f1
    out_f[2] = f2
    out_f[3] = f3
    out_df[0] = df0
    out_df[1] = df1
    out_df[2] = df2
    out_df[3] = df3


@_njit(cache=True, fastmath=True)
def _lam_UW(eta, A, B, Psi, prof, dprof):
    """Laminar (U, W, dU/deta, dW/deta) and derivatives w.r.t. (A, B, Psi).

    prof = (U, W, U', W'); dprof rows: d/dA, d/dB, d/dPsi of each entry.
    U = A*(1 - 0.6*(A-3)*eta^3)*f1 + f0 ; W = B*f2 + Psi*f3  (D13 (42)(43)).
    """
    f = np.empty(4)
    df = np.empty(4)
    _lam_f0123(eta, f, df)
    f0, f1, f2, f3 = f[0], f[1], f[2], f[3]
    df0, df1, df2, df3 = df[0], df[1], df[2], df[3]
    e3 = eta * eta * eta
    e2 = eta * eta
    cA = 1.0 - 0.6 * (A - 3.0) * e3          # multiplier of f1 in U
    dcA_dA = -0.6 * e3
    # U and its eta derivative
    U = A * cA * f1 + f0
    Up = A * (cA * df1 + (-1.8 * (A - 3.0) * e2) * f1) + df0
    # dU/dA: d[A*cA]/dA = cA + A*dcA_dA = 1 - 0.6*(2A-3)*e3
    dU_dA = (cA + A * dcA_dA) * f1
    # dU'/dA
    dUp_dA = (cA * df1 - 1.8 * (A - 3.0) * e2 * f1) + A * (
        dcA_dA * df1 - 1.8 * e2 * f1
    )
    W = B * f2 + Psi * f3
    Wp = B * df2 + Psi * df3
    prof[0] = U
    prof[1] = W
    prof[2] = Up
    prof[3] = Wp
    # rows: 0=d/dA, 1=d/dB, 2=d/dPsi
    dprof[0, 0] = dU_dA
    dprof[0, 1] = 0.0
    dprof[0, 2] = dUp_dA
    dprof[0, 3] = 0.0
    dprof[1, 0] = 0.0
    dprof[1, 1] = f2
    dprof[1, 2] = 0.0
    dprof[1, 3] = df2
    dprof[2, 0] = 0.0
    dprof[2, 1] = f3
    dprof[2, 2] = 0.0
    dprof[2, 3] = df3


# ---------------------------------------------------------------------------
# Turbulent profile family (D13 (47)-(57)): Spalding wall law + Coles wake
# ---------------------------------------------------------------------------

@_njit(cache=True, fastmath=True)
def _spalding_yplus(u):
    """y+_S(u+) and dy+/du+ (D13 (51)); vectorized-scalar, monotone."""
    ek = np.exp(-KAPPA * B_SPALDING)
    ku = KAPPA * u
    if ku > 60.0:
        # overflow guard: exponential dominates
        e_ = np.exp(ku)
        y = u + ek * e_
        dy = 1.0 + KAPPA * ek * e_
        return y, dy
    e_ = np.exp(ku)
    poly = 1.0 + ku + 0.5 * ku * ku + (ku * ku * ku) / 6.0
    dpoly = KAPPA * (1.0 + ku + 0.5 * ku * ku)
    y = u + ek * (e_ - poly)
    dy = 1.0 + ek * (KAPPA * e_ - dpoly)
    return y, dy


@_njit(cache=True, fastmath=True)
def _spalding_uplus(yplus):
    """u+(y+) by Newton on the monotone y+_S; returns (u+, du+/dy+)."""
    if yplus <= 0.0:
        return 0.0, 1.0
    # initial guess: viscous below the buffer layer, log law above
    if yplus < 11.0:
        u = yplus
    else:
        u = np.log(yplus) / KAPPA + B_SPALDING
    for _ in range(60):
        y, dy = _spalding_yplus(u)
        res = y - yplus
        step = res / dy
        u -= step
        if u < 0.0:
            u = 0.0
        if abs(step) < 1.0e-14 * (1.0 + abs(u)):
            break
    y, dy = _spalding_yplus(u)
    return u, 1.0 / dy


@_njit(cache=True, fastmath=True)
def _turb_scales(delta, A, B, re_d, out):
    """U_tau, W_tau, delta+ and derivatives w.r.t. (delta, A, B).

    D13 (55)(57) incompressible form (nu_w = nu_i):
      U_tau = A / (A^2+B^2)^{1/4} / sqrt(Re_delta)
      W_tau = B / (A^2+B^2)^{1/4} / sqrt(Re_delta)
      delta+ = sqrt(Re_delta) * (A^2+B^2)^{1/4}
    out = (U_tau, W_tau, delta+, dU_tau/d.., dW_tau/d.., ddelta+/d..)
    rows 6..14 hold d/d(delta, A, B) triples packed flat.
    """
    a2b2 = A * A + B * B
    a2b2 = max(a2b2, 1.0e-12)
    q14 = a2b2 ** 0.25
    sqr = np.sqrt(max(re_d, RE_D_MIN))
    u_t = A / (q14 * sqr)
    w_t = B / (q14 * sqr)
    dp = sqr * q14
    # derivatives of q14: dq14/dA = q14 * A / (2 a2b2)
    dq_dA = q14 * A / (2.0 * a2b2)
    dq_dB = q14 * B / (2.0 * a2b2)
    dsqr_dD = 0.5 * sqr / delta
    # U_tau = A q14^{-1} sqr^{-1}
    du_dD = -u_t / (2.0 * delta)
    du_dA = 1.0 / (q14 * sqr) - A * dq_dA / (q14 * q14 * sqr)
    du_dB = -A * dq_dB / (q14 * q14 * sqr)
    dw_dD = -w_t / (2.0 * delta)
    dw_dA = -B * dq_dA / (q14 * q14 * sqr)
    dw_dB = 1.0 / (q14 * sqr) - B * dq_dB / (q14 * q14 * sqr)
    ddp_dD = dsqr_dD * q14
    ddp_dA = sqr * dq_dA
    ddp_dB = sqr * dq_dB
    out[0] = u_t
    out[1] = w_t
    out[2] = dp
    out[3] = du_dD
    out[4] = du_dA
    out[5] = du_dB
    out[6] = dw_dD
    out[7] = dw_dA
    out[8] = dw_dB
    out[9] = ddp_dD
    out[10] = ddp_dA
    out[11] = ddp_dB


@_njit(cache=True, fastmath=True)
def _turb_UW(eta, delta, A, B, Psi, re_d, prof, dprof):
    """Turbulent (U, W, U', W') + derivatives w.r.t. (delta, A, B, Psi).

    D13 (47)(48): U = U_t u+(eta*d+) + K cos(Ups - Psi(1-eta)^2) g_o
                  W = W_t u+(eta*d+) - K sin(Ups - Psi(1-eta)^2) g_o
    dprof rows: d/dDelta, d/dA, d/dB, d/dPsi of each entry (4x4).
    """
    sc = np.empty(12)
    _turb_scales(delta, A, B, re_d, sc)
    u_t, w_t, dp = sc[0], sc[1], sc[2]
    # u+ at eta*d+ and at the edge d+
    up_e, dup_e = _spalding_uplus(dp)
    # K, Upsilon (D13 (53)(54))
    wx = w_t * up_e
    ux = 1.0 - u_t * up_e
    K = np.sqrt(wx * wx + ux * ux)
    K = max(K, 1.0e-30)
    Ups = np.arctan2(wx, ux)
    # d(K, Ups)/ds via chain rule; helpers: d(wx), d(ux)
    # derivatives of u_t, w_t, dp w.r.t. (delta, A, B)
    du_t = sc[3:6]
    dw_t = sc[6:9]
    ddp = sc[9:12]
    # d(up_e)/ds = dup_e * d(dp)/ds
    dwx = np.empty(3)
    dux = np.empty(3)
    for i in range(3):
        dwx[i] = dw_t[i] * up_e + w_t * dup_e * ddp[i]
        dux[i] = -(du_t[i] * up_e + u_t * dup_e * ddp[i])
    dK = np.empty(3)
    dUps = np.empty(3)
    den = max(ux * ux + wx * wx, 1.0e-30)
    for i in range(3):
        dK[i] = (wx * dwx[i] + ux * dux[i]) / K
        dUps[i] = (ux * dwx[i] - wx * dux[i]) / den
    # local u+ and its eta derivative
    y = eta * dp
    up, dup = _spalding_uplus(y)
    # shape functions
    g_o = 3.0 * eta * eta - 2.0 * eta * eta * eta
    dg_o = 6.0 * eta * (1.0 - eta)
    t = (1.0 - eta) * (1.0 - eta)
    ang = Ups - Psi * t
    ca = np.cos(ang)
    sa = np.sin(ang)
    U = u_t * up + K * ca * g_o
    W = w_t * up - K * sa * g_o
    # eta derivatives
    d_up_deta = dup * dp
    d_ang_deta = 2.0 * Psi * (1.0 - eta)
    Up_ = u_t * d_up_deta + K * (-sa * d_ang_deta * g_o + ca * dg_o)
    Wp_ = w_t * d_up_deta - K * (ca * d_ang_deta * g_o + sa * dg_o)
    prof[0] = U
    prof[1] = W
    prof[2] = Up_
    prof[3] = Wp_
    # state derivatives: rows d/dDelta, d/dA, d/dB; then d/dPsi
    for i in range(3):
        # d(u_t*up): up depends on state via y = eta*dp
        d_up = dup * eta * ddp[i]
        dU = du_t[i] * up + u_t * d_up + dK[i] * ca * g_o + K * (-sa) * dUps[i] * g_o
        dW = dw_t[i] * up + w_t * d_up - dK[i] * sa * g_o - K * ca * dUps[i] * g_o
        # eta-derivative of dU/dstate: Up_ = u_t*dup*dp + K(-sa*ang'*g_o + ca*dg_o)
        # (d ang'/ds = 0 for s in (delta, A, B); dup differentiates via
        # _spalding_d2up(y) * eta * ddp)
        dUp_s = (
            du_t[i] * d_up_deta
            + u_t * _spalding_d2up(y) * eta * ddp[i] * dp
            + u_t * dup * ddp[i]
            + dK[i] * (-sa * d_ang_deta * g_o + ca * dg_o)
            + K * (-ca * dUps[i] * d_ang_deta * g_o - sa * dUps[i] * dg_o)
        )
        dWp_s = (
            dw_t[i] * d_up_deta
            + w_t * _spalding_d2up(y) * eta * ddp[i] * dp
            + w_t * dup * ddp[i]
            - dK[i] * (ca * d_ang_deta * g_o + sa * dg_o)
            - K * (-sa * dUps[i] * d_ang_deta * g_o + ca * dUps[i] * dg_o)
        )
        dprof[i, 0] = dU
        dprof[i, 1] = dW
        dprof[i, 2] = dUp_s
        dprof[i, 3] = dWp_s
    # d/dPsi: d ang/dPsi = -t ; d ang'/dPsi = 2(1-eta) ; d sa = ca*(-t), d ca = sa*t
    dU_psi = K * sa * t * g_o
    dW_psi = K * ca * t * g_o
    dUp_psi = K * (ca * t * d_ang_deta * g_o - 2.0 * sa * (1.0 - eta) * g_o + sa * t * dg_o)
    dWp_psi = -K * (sa * t * d_ang_deta * g_o + 2.0 * ca * (1.0 - eta) * g_o - ca * t * dg_o)
    dprof[3, 0] = dU_psi
    dprof[3, 1] = dW_psi
    dprof[3, 2] = dUp_psi
    dprof[3, 3] = dWp_psi


@_njit(cache=True, fastmath=True)
def _spalding_d2up(yplus):
    """d^2u+/dy+^2 of the Spalding law (for eta-derivative chains)."""
    u, dup = _spalding_uplus(yplus)
    ek = np.exp(-KAPPA * B_SPALDING)
    ku = KAPPA * u
    if ku > 60.0:
        d2y = KAPPA * KAPPA * ek * np.exp(ku)
    else:
        e_ = np.exp(ku)
        d2poly = KAPPA * KAPPA * (1.0 + ku)
        d2y = ek * (KAPPA * KAPPA * e_ - d2poly)
    # d2u/dy2 = - (d2y/du2) / (dy/du)^3
    y, dy = _spalding_yplus(u)
    return -d2y / (dy * dy * dy)


# ---------------------------------------------------------------------------
# Density / viscosity profiles (D13 (58)(59)); M=0 -> R == 1 exactly
# ---------------------------------------------------------------------------

@_njit(cache=True, fastmath=True)
def _density_R(U, W, e_prime, d_hw):
    """R = rho/rho_i (Crocco-Busemann (58)); e_prime = r(g-1)/2 M_i^2."""
    inv = 1.0 + d_hw * (1.0 - U) + e_prime * (1.0 - U * U - W * W)
    return 1.0 / inv


# ---------------------------------------------------------------------------
# The nodal closure packet
# ---------------------------------------------------------------------------

@_njit(cache=True, fastmath=True)
def closure_node(state, q, rho, mu, mach, turbulent, c_l, out, dout):
    """Evaluate the 28-output closure packet and its (28,6) state derivatives.

    state = (delta, A, B, Psi, Ctau1, Ctau2); q/rho/mu/mach = edge parameters
    (held fixed; derivatives w.r.t. them are NOT produced in V1).
    turbulent: 0 = laminar family, 1 = turbulent family (forced-transition
    regime flag, D-TR). c_l: dissipation length constant (D-CT-2).

    The eta-quadrature engine evaluates every integrand of D13 (60)(61) and
    its analytic state derivative at the same Gauss points; outputs are the
    weighted sums. All thicknesses scale linearly in delta (integrals are
    per-unit-delta); the delta-derivative includes that scaling explicitly.
    """
    delta = state[0]
    dmult = 1.0
    if delta < DELTA_MIN:
        # probe state overshot to nonpositive thickness: evaluate at the
        # floor and mask d/d(delta) to zero (consistent FD-vs-analytic).
        delta = DELTA_MIN
        dmult = 0.0
    A = state[1]
    B = state[2]
    Psi = state[3]
    ct1 = state[4]
    ct2 = state[5]
    re_d = rho * q * delta / mu
    re_d = max(re_d, RE_D_MIN)
    e_prime = RECOVERY_R * 0.5 * (GAMMA_AIR - 1.0) * mach * mach
    d_hw = 0.0  # adiabatic wall (V1 scope; overheat ratio is a follow-up)
    if turbulent:
        eta_g = ETA_TURB
        w_g = W_TURB
        n_g = ETA_TURB.shape[0]
    else:
        eta_g = ETA_LAM
        w_g = W_LAM
        n_g = ETA_LAM.shape[0]

    # accumulators for the 28 outputs and their derivatives w.r.t.
    # (delta, A, B, Psi, ct1, ct2)
    acc = np.zeros(N_OUT)
    dacc = np.zeros((N_OUT, 6))

    # stress-sector outer weight shape w(eta) = 4 eta (1-eta) (D-CT-1)
    prof = np.empty(4)
    dprof = np.empty((4, 4))

    for k in range(n_g):
        eta = eta_g[k]
        wk = w_g[k]
        if turbulent:
            _turb_UW(eta, delta, A, B, Psi, re_d, prof, dprof)
        else:
            _lam_UW(eta, A, B, Psi, prof, dprof)
        U = prof[0]
        Wv = prof[1]
        Up = prof[2]
        Wp = prof[3]
        # derivative arrays: turbulent rows (delta, A, B, Psi); laminar rows
        # (A, B, Psi) -> remap; delta enters only via the overall delta factor.
        dU = np.zeros(6)
        dW = np.zeros(6)
        dUp = np.zeros(6)
        dWp = np.zeros(6)
        if turbulent:
            dU[0] = dprof[0, 0]
            dU[1] = dprof[1, 0]
            dU[2] = dprof[2, 0]
            dU[3] = dprof[3, 0]
            dW[0] = dprof[0, 1]
            dW[1] = dprof[1, 1]
            dW[2] = dprof[2, 1]
            dW[3] = dprof[3, 1]
            dUp[0] = dprof[0, 2]
            dUp[1] = dprof[1, 2]
            dUp[2] = dprof[2, 2]
            dUp[3] = dprof[3, 2]
            dWp[0] = dprof[0, 3]
            dWp[1] = dprof[1, 3]
            dWp[2] = dprof[2, 3]
            dWp[3] = dprof[3, 3]
        else:
            dU[1] = dprof[0, 0]
            dU[2] = dprof[1, 0]
            dU[3] = dprof[2, 0]
            dW[2] = dprof[1, 1]
            dW[3] = dprof[2, 1]
            dUp[1] = dprof[0, 2]
            dUp[2] = dprof[1, 2]
            dUp[3] = dprof[2, 2]
            dWp[2] = dprof[1, 3]
            dWp[3] = dprof[2, 3]

        R = _density_R(U, Wv, e_prime, d_hw)
        dR = np.zeros(6)
        # R = 1/(1 + dHw(1-U) + E'(1-U^2-W^2)); dR/ds = -R^2 * d(inv)/ds
        if e_prime != 0.0 or d_hw != 0.0:
            cof = -R * R
            for s in range(6):
                dinv = (-d_hw - 2.0 * e_prime * U) * dU[s] - 2.0 * e_prime * Wv * dW[s]
                dR[s] = cof * dinv

        # angle-deviation profile dpsi = atan2(W, U) (D13 (40); the "-Psi"
        # symbols inside D13 (60) are this pointwise deviation)
        u2w2 = U * U + Wv * Wv
        u2w2 = max(u2w2, 1.0e-30)
        dpsi = np.arctan2(Wv, U)
        d_dpsi = np.zeros(6)
        ddpsi_deta = (U * Wp - Wv * Up) / u2w2
        d_ddpsi = np.zeros(6)
        for s in range(6):
            d_dpsi[s] = (U * dW[s] - Wv * dU[s]) / u2w2
            # d(ddpsi_deta)/ds (product rule through the quotient)
            num = dU[s] * Wp + U * dWp[s] - dW[s] * Up - Wv * dUp[s]
            dnum = 2.0 * (U * dU[s] + Wv * dW[s])
            d_ddpsi[s] = (num * u2w2 - (U * Wp - Wv * Up) * dnum) / (u2w2 * u2w2)

        # ---- integrands of the 16 thicknesses (D13 (60)) ----
        # g[0..15] in packet order 0..15; each scaled by delta afterwards.
        g = np.empty(16)
        dg = np.zeros((16, 6))
        g[0] = 1.0 - R
        g[1] = 1.0 - R * U
        g[2] = -R * Wv
        g[3] = 1.0 - R * U * U
        g[4] = -R * U * Wv
        g[5] = -R * Wv * Wv
        g[6] = 1.0 - R * U * u2w2
        g[7] = -R * Wv * u2w2
        g[8] = 1.0 - U
        g[9] = -Wv
        g[10] = 1.0 - R * u2w2
        g[11] = -dpsi * u2w2 * R
        g[12] = -dpsi * u2w2 * R * U
        g[13] = -dpsi * u2w2 * R * Wv
        g[14] = -dpsi * U
        g[15] = -dpsi * Wv
        for s in range(6):
            dg[0, s] = -dR[s]
            dg[1, s] = -dR[s] * U - R * dU[s]
            dg[2, s] = -dR[s] * Wv - R * dW[s]
            dg[3, s] = -dR[s] * U * U - 2.0 * R * U * dU[s]
            dg[4, s] = -dR[s] * U * Wv - R * (dU[s] * Wv + U * dW[s])
            dg[5, s] = -dR[s] * Wv * Wv - 2.0 * R * Wv * dW[s]
            dg[6, s] = -dR[s] * U * u2w2 - R * (
                dU[s] * u2w2 + U * 2.0 * (U * dU[s] + Wv * dW[s])
            )
            dg[7, s] = -dR[s] * Wv * u2w2 - R * (
                dW[s] * u2w2 + Wv * 2.0 * (U * dU[s] + Wv * dW[s])
            )
            dg[8, s] = -dU[s]
            dg[9, s] = -dW[s]
            dg[10, s] = -dR[s] * u2w2 - R * 2.0 * (U * dU[s] + Wv * dW[s])
            dg[11, s] = -d_dpsi[s] * u2w2 * R - dpsi * (
                2.0 * (U * dU[s] + Wv * dW[s]) * R + u2w2 * dR[s]
            )
            dg[12, s] = -d_dpsi[s] * u2w2 * R * U - dpsi * (
                (2.0 * (U * dU[s] + Wv * dW[s]) * R + u2w2 * dR[s]) * U
                + u2w2 * R * dU[s]
            )
            dg[13, s] = -d_dpsi[s] * u2w2 * R * Wv - dpsi * (
                (2.0 * (U * dU[s] + Wv * dW[s]) * R + u2w2 * dR[s]) * Wv
                + u2w2 * R * dW[s]
            )
            dg[14, s] = -d_dpsi[s] * U - dpsi * dU[s]
            dg[15, s] = -d_dpsi[s] * Wv - dpsi * dW[s]

        # ---- shear profiles S, T (D13 (44)(45) laminar / (49)(50) turbulent)
        if turbulent:
            sc = np.empty(12)
            _turb_scales(delta, A, B, re_d, sc)
            u_t, w_t = sc[0], sc[1]
            dp = sc[2]
            du_t = sc[3:6]
            dw_t = sc[6:9]
            ddp = sc[9:12]
            up_e, dup_e = _spalding_uplus(dp)
            wx = w_t * up_e
            ux = 1.0 - u_t * up_e
            K = max(np.sqrt(wx * wx + ux * ux), 1.0e-30)
            Ups = np.arctan2(wx, ux)
            t1 = (1.0 - eta) * (1.0 - eta)
            ang = Ups - Psi * t1
            g_o = 3.0 * eta * eta - 2.0 * eta * eta * eta
            dg_o = 6.0 * eta * (1.0 - eta)
            mag = np.sqrt(u_t * u_t + w_t * w_t)
            mag = max(mag, 1.0e-30)
            ca = np.cos(ang)
            sa = np.sin(ang)
            S = R * (u_t * mag * (1.0 - g_o) + ct1 * K * ca * dg_o)
            T = R * (w_t * mag * (1.0 - g_o) - ct1 * K * sa * dg_o)
            # derivatives
            dS = np.zeros(6)
            dT = np.zeros(6)
            dmag = np.zeros(6)
            dmag[0] = (u_t * du_t[0] + w_t * dw_t[0]) / mag
            dmag[1] = (u_t * du_t[1] + w_t * dw_t[1]) / mag
            dmag[2] = (u_t * du_t[2] + w_t * dw_t[2]) / mag
            dK = np.zeros(6)
            dUp2 = np.zeros(6)
            den = max(ux * ux + wx * wx, 1.0e-30)
            dwx = np.empty(3)
            dux = np.empty(3)
            for i in range(3):
                dwx[i] = dw_t[i] * up_e + w_t * dup_e * ddp[i]
                dux[i] = -(du_t[i] * up_e + u_t * dup_e * ddp[i])
            for i in range(3):
                dK[i] = (wx * dwx[i] + ux * dux[i]) / K
                dUp2[i] = (ux * dwx[i] - wx * dux[i]) / den
            dang = np.zeros(6)
            for i in range(3):
                dang[i] = dUp2[i]
            dang[3] = -t1
            for s in range(6):
                dvisc1 = du_t[s] * mag + u_t * dmag[s] if s < 3 else 0.0
                dvisc2 = dw_t[s] * mag + w_t * dmag[s] if s < 3 else 0.0
                base1 = dvisc1 * (1.0 - g_o)
                base2 = dvisc2 * (1.0 - g_o)
                dKterm1 = dK[s] * ca - K * sa * dang[s]
                dKterm2 = -dK[s] * sa - K * ca * dang[s]
                ds1 = base1 + (1.0 if s == 4 else 0.0) * K * ca * dg_o + ct1 * dKterm1 * dg_o
                ds2 = base2 - (1.0 if s == 4 else 0.0) * K * sa * dg_o + ct1 * dKterm2 * dg_o
                dS[s] = dR[s] * (u_t * mag * (1.0 - g_o) + ct1 * K * ca * dg_o) + R * ds1
                dT[s] = dR[s] * (w_t * mag * (1.0 - g_o) - ct1 * K * sa * dg_o) + R * ds2
        else:
            # laminar (44)(45): S = (1/Re_d) dU/deta, T = (1/Re_d) dW/deta
            inv_re = 1.0 / re_d
            S = inv_re * Up
            T = inv_re * Wp
            dS = np.zeros(6)
            dT = np.zeros(6)
            dre_d = np.zeros(6)
            dre_d[0] = re_d / delta
            for s in range(6):
                dS[s] = -dre_d[s] * inv_re * inv_re * Up + inv_re * dUp[s]
                dT[s] = -dre_d[s] * inv_re * inv_re * Wp + inv_re * dWp[s]

        # ---- dissipation coefficients (D13 (61)) ----
        # cD = int(S U' + T W'), cDx = int(S W' - T U'),
        # cDc = int(S (dpsi U)' + T (dpsi W)')
        ddpsiU_deta = ddpsi_deta * U + dpsi * Up
        ddpsiW_deta = ddpsi_deta * Wv + dpsi * Wp
        # derivatives of ddpsiU_deta, ddpsiW_deta
        d_ddpsiU = np.zeros(6)
        d_ddpsiW = np.zeros(6)
        for s in range(6):
            d_ddpsiU[s] = d_ddpsi[s] * U + ddpsi_deta * dU[s] + d_dpsi[s] * Up + dpsi * dUp[s]
            d_ddpsiW[s] = d_ddpsi[s] * Wv + ddpsi_deta * dW[s] + d_dpsi[s] * Wp + dpsi * dWp[s]

        cD_g = S * Up + T * Wp
        cDx_g = S * Wp - T * Up
        cDc_g = S * ddpsiU_deta + T * ddpsiW_deta
        d_cD = np.zeros(6)
        d_cDx = np.zeros(6)
        d_cDc = np.zeros(6)
        for s in range(6):
            d_cD[s] = dS[s] * Up + S * dUp[s] + dT[s] * Wp + T * dWp[s]
            d_cDx[s] = dS[s] * Wp + S * dWp[s] - dT[s] * Up - T * dUp[s]
            d_cDc[s] = dS[s] * ddpsiU_deta + S * d_ddpsiU[s] + dT[s] * ddpsiW_deta + T * d_ddpsiW[s]

        # ---- stress-transport factors (D-CT-1) ----
        wsig = 4.0 * eta * (1.0 - eta)
        ak_g = R * wsig
        sp1_g = R * wsig * Up
        sp2_g = R * wsig * Wp
        sd_g = (R * wsig) ** 1.5
        ku1_g = R * wsig * U
        ku2_g = R * wsig * Wv
        d_ak = np.zeros(6)
        d_sp1 = np.zeros(6)
        d_sp2 = np.zeros(6)
        d_sd = np.zeros(6)
        d_ku1 = np.zeros(6)
        d_ku2 = np.zeros(6)
        for s in range(6):
            d_ak[s] = dR[s] * wsig
            d_sp1[s] = dR[s] * wsig * Up + R * wsig * dUp[s]
            d_sp2[s] = dR[s] * wsig * Wp + R * wsig * dWp[s]
            d_sd[s] = 1.5 * np.sqrt(R * wsig) * wsig * dR[s]
            d_ku1[s] = dR[s] * wsig * U + R * wsig * dU[s]
            d_ku2[s] = dR[s] * wsig * Wv + R * wsig * dW[s]

        # ---- accumulate (weighted) ----
        for j in range(16):
            acc[j] += wk * g[j]
            for s in range(6):
                dacc[j, s] += wk * dg[j, s]
        acc[OUT_CD] += wk * cD_g
        acc[OUT_CDX] += wk * cDx_g
        acc[OUT_CDC] += wk * cDc_g
        acc[OUT_AK] += wk * ak_g
        acc[OUT_SP1] += wk * sp1_g
        acc[OUT_SP2] += wk * sp2_g
        acc[OUT_SD] += wk * sd_g
        acc[OUT_KU1] += wk * ku1_g
        acc[OUT_KU2] += wk * ku2_g
        for s in range(6):
            dacc[OUT_CD, s] += wk * d_cD[s]
            dacc[OUT_CDX, s] += wk * d_cDx[s]
            dacc[OUT_CDC, s] += wk * d_cDc[s]
            dacc[OUT_AK, s] += wk * d_ak[s]
            dacc[OUT_SP1, s] += wk * d_sp1[s]
            dacc[OUT_SP2, s] += wk * d_sp2[s]
            dacc[OUT_SD, s] += wk * d_sd[s]
            dacc[OUT_KU1, s] += wk * d_ku1[s]
            dacc[OUT_KU2, s] += wk * d_ku2[s]

        # Cf at the wall is taken from the eta=0 endpoint of S, T (evaluated
        # once outside the loop below for exactness).

    # mask d/d(delta) everywhere when the DELTA_MIN floor is active
    if dmult == 0.0:
        for j in range(N_OUT):
            dacc[j, 0] = 0.0

    # ---- wall values (D13 (61)): cf1 = 2 S(0), cf2 = 2 T(0) ----
    if turbulent:
        sc = np.empty(12)
        _turb_scales(delta, A, B, re_d, sc)
        u_t, w_t = sc[0], sc[1]
        mag = np.sqrt(u_t * u_t + w_t * w_t)
        mag = max(mag, 1.0e-30)
        R0 = _density_R(0.0, 0.0, e_prime, d_hw)
        S0 = R0 * u_t * mag
        T0 = R0 * w_t * mag
        cf1 = 2.0 * S0
        cf2 = 2.0 * T0
        dcf1 = np.zeros(6)
        dcf2 = np.zeros(6)
        du_t = sc[3:6]
        dw_t = sc[6:9]
        dmag = np.zeros(6)
        for i in range(3):
            dmag[i] = (u_t * du_t[i] + w_t * dw_t[i]) / mag
        for s in range(3):
            dcf1[s] = 2.0 * R0 * (du_t[s] * mag + u_t * dmag[s])
            dcf2[s] = 2.0 * R0 * (dw_t[s] * mag + w_t * dmag[s])
        # R0 depends on (A,B,Psi) via dHw*(1-U)+E'(1-U^2-W^2) at U=W=0 -> const
    else:
        inv_re = 1.0 / re_d
        # U'(0) = A, W'(0) = B exactly in the laminar family
        cf1 = 2.0 * inv_re * A
        cf2 = 2.0 * inv_re * B
        dcf1 = np.zeros(6)
        dcf2 = np.zeros(6)
        dcf1[0] = -cf1 / delta
        dcf1[1] = 2.0 * inv_re
        dcf2[0] = -cf2 / delta
        dcf2[2] = 2.0 * inv_re
    dcf1[0] *= dmult
    dcf2[0] *= dmult

    # ---- scale thicknesses by delta (integrals were per-unit-delta) ----
    for j in range(16):
        v = acc[j]
        acc[j] = delta * v
        for s in range(6):
            dacc[j, s] = (dmult if s == 0 else 0.0) * v + delta * dacc[j, s]

    out[0:16] = acc[0:16]
    out[OUT_CF1] = cf1
    out[OUT_CF2] = cf2
    out[OUT_CD] = acc[OUT_CD]
    out[OUT_CDX] = acc[OUT_CDX]
    out[OUT_CDC] = acc[OUT_CDC]
    out[OUT_AK] = acc[OUT_AK]
    out[OUT_SP1] = acc[OUT_SP1]
    out[OUT_SP2] = acc[OUT_SP2]
    out[OUT_SD] = acc[OUT_SD]
    out[OUT_KU1] = acc[OUT_KU1]
    out[OUT_KU2] = acc[OUT_KU2]
    for s in range(6):
        dout[0:16, s] = dacc[0:16, s]
        dout[OUT_CF1, s] = dcf1[s]
        dout[OUT_CF2, s] = dcf2[s]
        dout[OUT_CD, s] = dacc[OUT_CD, s]
        dout[OUT_CDX, s] = dacc[OUT_CDX, s]
        dout[OUT_CDC, s] = dacc[OUT_CDC, s]
        dout[OUT_AK, s] = dacc[OUT_AK, s]
        dout[OUT_SP1, s] = dacc[OUT_SP1, s]
        dout[OUT_SP2, s] = dacc[OUT_SP2, s]
        dout[OUT_SD, s] = dacc[OUT_SD, s]
        dout[OUT_KU1, s] = dacc[OUT_KU1, s]
        dout[OUT_KU2, s] = dacc[OUT_KU2, s]

    # ---- derived: theta11 = p11 - ds1, theta*1 = ps1 - ds1, H = ds1/theta11
    th = acc[OUT_P11] - acc[OUT_DS1]
    ths = acc[OUT_PS1] - acc[OUT_DS1]
    out[OUT_TH11] = th
    out[OUT_THS1] = ths
    th_safe = th if abs(th) > 1.0e-30 else 1.0e-30
    out[OUT_H1] = acc[OUT_DS1] / th_safe
    for s in range(6):
        dout[OUT_TH11, s] = dacc[OUT_P11, s] - dacc[OUT_DS1, s]
        dout[OUT_THS1, s] = dacc[OUT_PS1, s] - dacc[OUT_DS1, s]
        dout[OUT_H1, s] = (dacc[OUT_DS1, s] * th - acc[OUT_DS1] * (dacc[OUT_P11, s] - dacc[OUT_DS1, s])) / (th_safe * th_safe)


# ---------------------------------------------------------------------------
# Stress-equation source (D13 (32) + D-CT closure, design doc §3.2)
# ---------------------------------------------------------------------------

@_njit(cache=True, fastmath=True)
def stress_source(state, q, rho, c_l, sp, sd, comp, out, dout):
    """S_tau = 2 a1 (P - D) for stress component `comp` (0 or 1), and its
    state derivatives. P_c = rho q^3 Ctau_c sp_c ; D_c = rho q^3 (delta/L)
    |Ctau_c|^{1/2} Ctau_c sd with L = c_l * delta  =>  delta/L = 1/c_l
    exactly, so D carries no explicit delta factor (the state dependence of
    sp_c, sd is the caller's, via the product rule in ibl3.py).
    out = (P, D, S); dout = d/d(state), 3x6.
    """
    ct = state[4 + comp]
    P = rho * q ** 3 * ct * sp
    ctm = max(ct, 1.0e-30)
    D = rho * q ** 3 / c_l * np.sqrt(ctm) * ct * sd
    S = 2.0 * A1_BRADSHAW * (P - D)
    out[0] = P
    out[1] = D
    out[2] = S
    for s in range(6):
        dP = 0.0
        dD = 0.0
        if s == 4 + comp:
            dP = rho * q ** 3 * sp
            dD = rho * q ** 3 / c_l * sd * 1.5 * np.sqrt(ctm)
        dout[0, s] = dP
        dout[1, s] = dD
        dout[2, s] = 2.0 * A1_BRADSHAW * (dP - dD)


# ---------------------------------------------------------------------------
# Vectorized driver: evaluate packets for many nodes
# ---------------------------------------------------------------------------

@_njit(cache=True, fastmath=True, parallel=True)
def closure_all(states, q, rho, mu, mach, turbulent_flags, c_l, outs, douts):
    """closure_node over all surface nodes; outs (n,30), douts (n,30,6)
    preallocated by the caller and written in place (design.md §7 rule 4).
    Per-node embarrassingly parallel loop (no scatter) -> plain prange.
    """
    n = states.shape[0]
    for i in prange(n):
        closure_node(
            states[i], q[i], rho[i], mu[i], mach[i], turbulent_flags[i], c_l,
            outs[i], douts[i],
        )


# ---------------------------------------------------------------------------
# Convenience: single-state Python wrapper (tests / seeding)
# ---------------------------------------------------------------------------

def closure_scalar(state, q=1.0, rho=1.0, mu=1.0e-5, mach=0.0, turbulent=False,
                   c_l=C_L_DEFAULT):
    """Python-level single-node evaluation; returns (out (28,), dout (28,6))."""
    st = np.ascontiguousarray(state, dtype=np.float64)
    out = np.empty(N_OUT)
    dout = np.empty((N_OUT, 6))
    closure_node(st, q, rho, mu, mach, 1 if turbulent else 0, c_l, out, dout)
    return out, dout


# ---------------------------------------------------------------------------
# Blasius seeding map (GV1.1 inflow Dirichlet; design doc §7(a))
# ---------------------------------------------------------------------------

def blasius_A(target_H=2.5906, mu=1.0e-5, q=1.0, rho=1.0):
    """Solve H_lam(A) = target_H for A (1-D damped Newton on the family)."""
    A = 8.0
    for _ in range(60):
        out, dout = closure_scalar((1.0e-3, A, 0.0, 0.0, CTAU_LAM, 0.0),
                                   q=q, rho=rho, mu=mu, turbulent=False)
        H = out[OUT_H1]
        dH_dA = dout[OUT_H1, 1]
        step = (H - target_H) / dH_dA
        step = max(-2.0, min(2.0, step))
        A -= step
        if abs(step) < 1.0e-13:
            break
    return A


def blasius_seed(x, q=1.0, rho=1.0, mu=1.0e-5):
    """(delta, A, 0, 0, CTAU_LAM, 0) matching Blasius theta and H at station x.

    theta_Blasius = 0.664 x / sqrt(Re_x); A from H=2.5906; delta from theta.
    The profile family is not exactly Blasius: H and theta are matched
    exactly, cf is matched approximately (residual recorded in the gate).
    """
    re_x = rho * q * x / mu
    theta_target = 0.664 * x / np.sqrt(re_x)
    A = blasius_A(mu=mu, q=q, rho=rho)
    out, _ = closure_scalar((1.0e-3, A, 0.0, 0.0, CTAU_LAM, 0.0),
                            q=q, rho=rho, mu=mu, turbulent=False)
    theta_hat = out[OUT_TH11] / 1.0e-3  # per-unit-delta theta
    delta = theta_target / theta_hat
    return np.array([delta, A, 0.0, 0.0, CTAU_LAM, 0.0])
