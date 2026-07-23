"""Track V1 Stage 2: closure packet unit tests (binding: docs/design_track_v.md,
Drela AIAA 2013-2437 equations as cited in pyfp3d/viscous/closures.py).

Coverage:
- analytic (28,6) state derivatives vs central-step FD, laminar/turbulent x
  2-D/3-D states (noise-aware: exact-zero entries skipped below 1e-8);
- structural crossflow zeros in pure 2-D states (B = Psi = Ctau2 = 0);
- laminar thickness integrals vs independent trapezoidal quadrature;
- Spalding inner-law round trip y+ -> u+ -> y+;
- turbulent profile edge conditions U(1) = 1, W(1) = 0;
- blasius_seed / blasius_A self-consistency (H matched, theta matched);
- stress_source FD + equilibrium Ctau_eq = (c_l sp/sd)^2 (D-CT);
- safety floors: nonpositive delta stays finite with masked d/d(delta);
- closure_all parity with per-node closure_scalar, zero-length safety.

Runs in both lanes: default JIT and PYFP3D_NOJIT=1.
"""

import numpy as np
import pytest

from pyfp3d.viscous import closures as C

Q0 = 1.0
RHO0 = 1.0
MU0 = 1.0e-5

LAM_2D = np.array([0.05, 8.03, 0.0, 0.0, C.CTAU_LAM, C.CTAU_LAM])
LAM_3D = np.array([0.05, 8.03, 1.2, 0.35, C.CTAU_LAM, C.CTAU_LAM])
TURB_2D = np.array([1.0, 16.0, 0.0, 0.0, 0.03, 0.0])       # Re_delta = 1e5
TURB_3D = np.array([1.0, 16.0, 2.0, 0.4, 0.03, 0.01])

# outputs that must vanish identically in a pure 2-D state (integrals odd in
# the crossflow profile W, or products thereof)
CROSSFLOW_ZERO_OUTPUTS = (
    C.OUT_DS2, C.OUT_P12, C.OUT_PS2, C.OUT_DQ2, C.OUT_TC2, C.OUT_DC2,
    C.OUT_CF2, C.OUT_CDX, C.OUT_CDC, C.OUT_SP2,
)


def _fd_worst_rel(state, turbulent, mach=0.0, eps=1.0e-6):
    """Worst relative |fd - analytic| over the (28,6) block (central FD);
    entries with both |fd| and |analytic| below 1e-8 are exact zeros (noise
    floor) and skipped."""
    out0, d0, _ = C.closure_scalar(state, q=Q0, rho=RHO0, mu=MU0, mach=mach,
                                   turbulent=turbulent)
    assert np.all(np.isfinite(out0))
    assert np.all(np.isfinite(d0))
    worst = 0.0
    for s in range(6):
        st = np.array(state, dtype=float)
        st[s] += eps
        out1, _, _ = C.closure_scalar(st, q=Q0, rho=RHO0, mu=MU0, mach=mach,
                                      turbulent=turbulent)
        st[s] -= 2.0 * eps
        out2, _, _ = C.closure_scalar(st, q=Q0, rho=RHO0, mu=MU0, mach=mach,
                                      turbulent=turbulent)
        fd = (out1 - out2) / (2.0 * eps)
        for j in range(C.N_OUT):
            if abs(fd[j]) < 1.0e-8 and abs(d0[j, s]) < 1.0e-8:
                continue
            denom = max(abs(fd[j]), abs(d0[j, s]), 1.0e-12)
            worst = max(worst, abs(fd[j] - d0[j, s]) / denom)
    return worst


# ---------------------------------------------------------------------------
# FD consistency of the analytic derivative block
# ---------------------------------------------------------------------------

def test_fd_laminar_2d():
    assert _fd_worst_rel(LAM_2D, False) < 5.0e-5


def test_fd_laminar_3d():
    assert _fd_worst_rel(LAM_3D, False) < 5.0e-5


def test_fd_turbulent_2d():
    assert _fd_worst_rel(TURB_2D, True) < 5.0e-5


def test_fd_turbulent_3d():
    assert _fd_worst_rel(TURB_3D, True) < 5.0e-5


def test_fd_compressible_mach():
    # Crocco-Busemann density path (mach > 0) stays FD-consistent
    assert _fd_worst_rel(LAM_3D, False, mach=0.3) < 5.0e-5
    assert _fd_worst_rel(TURB_3D, True, mach=0.3) < 5.0e-5


# ---------------------------------------------------------------------------
# Structural zeros in pure 2-D states
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("state,turb", [(LAM_2D, False), (TURB_2D, True)])
def test_crossflow_structural_zeros(state, turb):
    out, _, _ = C.closure_scalar(state, q=Q0, rho=RHO0, mu=MU0, turbulent=turb)
    for j in CROSSFLOW_ZERO_OUTPUTS:
        assert abs(out[j]) < 1.0e-14, f"output {j} = {out[j]} not ~0 in 2-D"


# ---------------------------------------------------------------------------
# Laminar thicknesses vs independent trapezoidal quadrature
# ---------------------------------------------------------------------------

def test_laminar_thickness_independent_quadrature():
    n = 20001
    eta = np.linspace(0.0, 1.0, n)
    A, B, Psi = LAM_3D[1], LAM_3D[2], LAM_3D[3]
    U = np.empty(n)
    W = np.empty(n)
    prof = np.empty(4)
    dprof = np.empty((4, 4))
    for i in range(n):
        C._lam_UW(eta[i], A, B, Psi, prof, dprof)
        U[i] = prof[0]
        W[i] = prof[1]
    ds1_trap = np.trapezoid(1.0 - U, eta)
    th11_trap = np.trapezoid(U * (1.0 - U), eta)
    out, _, _ = C.closure_scalar(LAM_3D, q=Q0, rho=RHO0, mu=MU0, turbulent=False)
    delta = LAM_3D[0]
    assert abs(out[C.OUT_DS1] / delta - ds1_trap) < 1.0e-6
    assert abs(out[C.OUT_TH11] / delta - th11_trap) < 1.0e-6
    # shape factor from the independent integration must match too
    H_trap = ds1_trap / th11_trap
    assert abs(out[C.OUT_H1] - H_trap) < 1.0e-5


# ---------------------------------------------------------------------------
# Spalding law round trip
# ---------------------------------------------------------------------------

def test_spalding_roundtrip():
    yplus = np.logspace(-2, 6, 41)
    for y in yplus:
        u, _ = C._spalding_uplus(y)
        y_back, _ = C._spalding_yplus(u)
        assert abs(y_back - y) / y < 1.0e-10, f"y+={y}: {y_back}"


def test_turbulent_edge_conditions():
    # U(eta=1) = 1, W(eta=1) = 0 identically in the turbulent family
    delta, A, B, Psi = TURB_3D[0], TURB_3D[1], TURB_3D[2], TURB_3D[3]
    re_d = RHO0 * Q0 * delta / MU0
    prof = np.empty(4)
    dprof = np.empty((4, 4))
    dprof_re = np.empty(4)
    C._turb_UW(1.0, delta, A, B, Psi, re_d, prof, dprof, dprof_re)
    assert abs(prof[0] - 1.0) < 1.0e-12
    assert abs(prof[1]) < 1.0e-12


# ---------------------------------------------------------------------------
# Blasius seeding map
# ---------------------------------------------------------------------------

def test_blasius_seed_consistency():
    x = 1.0
    state = C.blasius_seed(x, q=Q0, rho=RHO0, mu=MU0)
    assert state[0] > 0.0
    assert state[2] == 0.0 and state[3] == 0.0
    out, _, _ = C.closure_scalar(state, q=Q0, rho=RHO0, mu=MU0, turbulent=False)
    # H matched to the Blasius target by construction of blasius_A
    assert abs(out[C.OUT_H1] - 2.5906) < 1.0e-6
    # theta matched to 0.664 x / sqrt(Re_x)
    theta_target = 0.664 * x / np.sqrt(RHO0 * Q0 * x / MU0)
    assert abs(out[C.OUT_TH11] - theta_target) / theta_target < 1.0e-8


# ---------------------------------------------------------------------------
# Stress-equation source (D-CT)
# ---------------------------------------------------------------------------

def test_stress_source_fd():
    state = TURB_3D
    out_c, _, _ = C.closure_scalar(state, q=Q0, rho=RHO0, mu=MU0, turbulent=True)
    for comp, sp, sd in ((0, out_c[C.OUT_SP1], out_c[C.OUT_SD]),
                         (1, out_c[C.OUT_SP2], out_c[C.OUT_SD])):
        o = np.empty(3)
        d = np.empty((3, 6))
        C.stress_source(state, Q0, RHO0, C.C_L_DEFAULT, sp, sd, comp, o, d)
        assert np.all(np.isfinite(o)) and np.all(np.isfinite(d))
        eps = 1.0e-6
        for s in range(6):
            st = np.array(state)
            st[s] += eps
            o1 = np.empty(3)
            d1 = np.empty((3, 6))
            C.stress_source(st, Q0, RHO0, C.C_L_DEFAULT, sp, sd, comp, o1, d1)
            st[s] -= 2.0 * eps
            o2 = np.empty(3)
            d2 = np.empty((3, 6))
            C.stress_source(st, Q0, RHO0, C.C_L_DEFAULT, sp, sd, comp, o2, d2)
            fd = (o1 - o2) / (2.0 * eps)
            for j in range(3):
                if abs(fd[j]) < 1.0e-8 and abs(d[j, s]) < 1.0e-8:
                    continue
                denom = max(abs(fd[j]), abs(d[j, s]), 1.0e-12)
                assert abs(fd[j] - d[j, s]) / denom < 5.0e-5


def test_stress_source_equilibrium():
    # S = 2 a1 (P - D) = 0 at Ctau_eq = (c_l sp / sd)^2 (D-CT equilibrium)
    out_c, _, _ = C.closure_scalar(TURB_2D, q=Q0, rho=RHO0, mu=MU0, turbulent=True)
    sp, sd = out_c[C.OUT_SP1], out_c[C.OUT_SD]
    ct_eq = (C.C_L_DEFAULT * sp / sd) ** 2
    state = TURB_2D.copy()
    state[4] = ct_eq
    o = np.empty(3)
    d = np.empty((3, 6))
    C.stress_source(state, Q0, RHO0, C.C_L_DEFAULT, sp, sd, 0, o, d)
    assert abs(o[2]) / max(abs(o[0]), 1.0e-30) < 1.0e-10
    assert ct_eq > 0.0


# ---------------------------------------------------------------------------
# Safety floors
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("turb", [False, True])
def test_nonpositive_delta_floored(turb):
    for d0 in (0.0, -1.0e-3, 1.0e-10):
        state = np.array([d0, 8.0, 0.5, 0.2, 1.0e-8, 1.0e-8])
        out, dout, _ = C.closure_scalar(state, q=Q0, rho=RHO0, mu=MU0,
                                        turbulent=turb)
        assert np.all(np.isfinite(out))
        assert np.all(np.isfinite(dout))
        assert np.all(dout[:, 0] == 0.0)  # d/d(delta) masked below the floor
    # two sub-floor probes must give identical outputs (piecewise constant)
    oa, _, _ = C.closure_scalar(np.array([1.0e-10, 8.0, 0.5, 0.2, 1e-8, 1e-8]),
                                q=Q0, rho=RHO0, mu=MU0, turbulent=turb)
    ob, _, _ = C.closure_scalar(np.array([2.0e-10, 8.0, 0.5, 0.2, 1e-8, 1e-8]),
                                q=Q0, rho=RHO0, mu=MU0, turbulent=turb)
    assert np.allclose(oa, ob, rtol=0.0, atol=0.0)


# ---------------------------------------------------------------------------
# closure_all parity / robustness
# ---------------------------------------------------------------------------

def test_closure_all_parity():
    rng = np.random.default_rng(0)
    n = 50
    states = np.zeros((n, 6))
    states[:, 0] = rng.uniform(0.01, 2.0, n)
    states[:, 1] = rng.uniform(2.0, 25.0, n)
    states[:, 2] = rng.uniform(-3.0, 3.0, n)
    states[:, 3] = rng.uniform(-0.5, 0.5, n)
    states[:, 4] = rng.uniform(1.0e-8, 0.05, n)
    states[:, 5] = rng.uniform(-0.02, 0.02, n)
    q = np.full(n, Q0)
    rho = np.full(n, RHO0)
    mu = np.full(n, MU0)
    mach = np.zeros(n)
    flags = (rng.random(n) > 0.5).astype(np.int64)
    outs = np.empty((n, C.N_OUT))
    douts = np.empty((n, C.N_OUT, 6))
    douts_e = np.empty((n, C.N_OUT, 2))
    C.closure_all(states, q, rho, mu, mach, flags, C.C_L_DEFAULT, outs, douts,
                  douts_e)
    for i in range(n):
        o, d, de = C.closure_scalar(states[i], q=Q0, rho=RHO0, mu=MU0,
                                    turbulent=bool(flags[i]))
        assert np.allclose(outs[i], o, rtol=0.0, atol=0.0)
        assert np.allclose(douts[i], d, rtol=0.0, atol=0.0)
        assert np.allclose(douts_e[i], de, rtol=0.0, atol=0.0)


def test_closure_all_zero_length():
    empty6 = np.empty((0, 6))
    empty1 = np.empty(0)
    flags = np.empty(0, dtype=np.int64)
    outs = np.empty((0, C.N_OUT))
    douts = np.empty((0, C.N_OUT, 6))
    douts_e = np.empty((0, C.N_OUT, 2))
    C.closure_all(empty6, empty1, empty1, empty1, empty1, flags,
                  C.C_L_DEFAULT, outs, douts, douts_e)  # must not raise
