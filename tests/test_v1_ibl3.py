"""Track V1 Stage 3: IBL3 solver tests (binding: docs/design_track_v.md §5,
Drela AIAA 2013-2437).

Coverage:
- analytic Jacobian vs central FD (Jacobian-vector products on random
  directions; includes Dirichlet inflow rows and laminar pin rows);
- bit-determinism of the colored assembly;
- Newton convergence on the laminar flat plate (Blasius-seeded inflow),
  family-consistent H and delta* growth;
- quasi-2-D structural lock: B, Psi, Ctau2 stay at machine zero from a 2-D
  seed (design doc GV1.1(d) mechanism);
- forced-turbulent plate smoke: converges from a power-law seed, physical
  Cf/H ballpark (the +-5% reference comparison is the GV1.1(b) gate, not
  this unit test);
- inflow Dirichlet values retained; laminar Ctau pinned at CTAU_LAM.

Runs in both lanes: default JIT and PYFP3D_NOJIT=1.
"""

import numpy as np
import pytest

from pyfp3d.viscous import closures as C
from pyfp3d.viscous.ibl3 import IBL3Solver
from pyfp3d.viscous.surface_mesh import (
    SurfaceMesh,
    structured_rectangle_surface,
)

Q0 = 1.0
RHO0 = 1.0
MU0 = 1.0e-5


def _plate_mesh(nx=10, nz=5, x0=0.2, x1=1.2, z0=-0.3, z1=0.3):
    xyz, tris = structured_rectangle_surface(x0, x1, z0, z1, nx, nz)
    return SurfaceMesh.from_wall_faces(xyz, tris)


def _lam_seed_field(xyz):
    n = xyz.shape[0]
    U = np.zeros((n, 6))
    for i in range(n):
        U[i] = C.blasius_seed(max(xyz[i, 0], 1.0e-3), q=Q0, rho=RHO0, mu=MU0)
    return U


def _turb_seed_state(x, q=Q0, rho=RHO0, mu=MU0):
    """Power-law flat-plate turbulent seed: delta/x = 0.37 Re_x^-0.2,
    cf = 0.0576 Re_x^-0.2; A from U_tau^2 = cf/2 with U_tau = sqrt(A/Re_delta)
    (B=0 limit of D13 (57)); Ctau from the D-CT equilibrium of the closure.
    """
    re_x = rho * q * x / mu
    delta = 0.37 * x * re_x ** -0.2
    cf = 0.0576 * re_x ** -0.2
    re_d = rho * q * delta / mu
    # U_tau^2 = cf/2 and U_tau = sqrt(A/Re_delta) (B=0 limit of D13 (57))
    A = 0.5 * cf * re_d
    state = np.array([delta, A, 0.0, 0.0, 1.0e-3, 0.0])
    out, _ = C.closure_scalar(state, q=q, rho=rho, mu=mu, turbulent=True)
    ct_eq = (C.C_L_DEFAULT * out[C.OUT_SP1] / out[C.OUT_SD]) ** 2
    state[4] = max(ct_eq, 1.0e-6)
    return state


def _make_solver(sm, flags, x0=0.2):
    n = sm.xyz.shape[0]
    u_e = np.zeros((n, 3))
    u_e[:, 0] = Q0
    inflow = np.abs(sm.xyz[:, 0] - x0) < 1.0e-12
    if flags[0] == 1:
        st = _turb_seed_state(x0)
    else:
        st = C.blasius_seed(x0, q=Q0, rho=RHO0, mu=MU0)
    return IBL3Solver(
        sm, u_e, RHO0, MU0, 0.0, flags, inflow, st,
    )


# ---------------------------------------------------------------------------
# Analytic Jacobian vs FD
# ---------------------------------------------------------------------------

def test_jacobian_fd_random_directions():
    sm = _plate_mesh(nx=5, nz=3)
    n = sm.xyz.shape[0]
    rng = np.random.default_rng(1)
    flags = (rng.random(n) > 0.5).astype(np.int64)
    solver = _make_solver(sm, flags)
    # nontrivial but physical state field
    U = _lam_seed_field(sm.xyz)
    U[:, 1] += 0.5 * np.sin(3.0 * sm.xyz[:, 0])
    U[:, 2] = 0.3 * np.cos(5.0 * sm.xyz[:, 2])
    U[:, 3] = 0.1 * np.sin(4.0 * sm.xyz[:, 0])
    U[:, 4] = np.where(flags == 1, 0.02, C.CTAU_LAM)
    U[:, 5] = 0.005 * np.where(flags == 1, 1.0, 0.0)
    R0, J = solver.residual_jacobian(U)
    assert np.all(np.isfinite(R0))
    assert np.all(np.isfinite(J.data))
    eps = 1.0e-6
    for _ in range(4):
        v = rng.standard_normal((n, 6))
        v /= np.max(np.abs(v))
        Rp = solver.residual(U + eps * v)
        Rm = solver.residual(U - eps * v)
        fd = (Rp - Rm).ravel() / (2.0 * eps)
        an = J @ v.ravel()
        err = np.max(np.abs(fd - an)) / max(np.max(np.abs(an)), 1.0e-12)
        assert err < 1.0e-5, f"Jacobian FD mismatch: {err:.3e}"


def test_jacobian_fd_with_pseudotime():
    """Same FD check but through the pseudo-time residual path (D-PT terms
    and their nodal block-diagonal Jacobian included)."""
    sm = _plate_mesh(nx=5, nz=3)
    n = sm.xyz.shape[0]
    rng = np.random.default_rng(2)
    flags = (rng.random(n) > 0.5).astype(np.int64)
    solver = _make_solver(sm, flags)
    U = _lam_seed_field(sm.xyz)
    U[:, 2] = 0.2 * np.cos(5.0 * sm.xyz[:, 2])
    U[:, 3] = 0.1 * np.sin(4.0 * sm.xyz[:, 0])
    U[:, 4] = np.where(flags == 1, 0.02, C.CTAU_LAM)
    U_old = _lam_seed_field(sm.xyz) * 0.9
    h = np.sqrt(sm.node_area)
    dt_inv = solver._q / (3.0 * h)
    R0, J = solver.residual_jacobian(U, dt_inv=dt_inv, U_old=U_old)
    eps = 1.0e-6
    for _ in range(4):
        v = rng.standard_normal((n, 6))
        v /= np.max(np.abs(v))
        Rp = solver.residual(U + eps * v, dt_inv=dt_inv, U_old=U_old)
        Rm = solver.residual(U - eps * v, dt_inv=dt_inv, U_old=U_old)
        fd = (Rp - Rm).ravel() / (2.0 * eps)
        an = J @ v.ravel()
        err = np.max(np.abs(fd - an)) / max(np.max(np.abs(an)), 1.0e-12)
        assert err < 1.0e-5, f"pseudo-time Jacobian FD mismatch: {err:.3e}"


def test_assembly_bit_deterministic():
    sm = _plate_mesh(nx=5, nz=3)
    n = sm.xyz.shape[0]
    flags = np.zeros(n, dtype=np.int64)
    solver = _make_solver(sm, flags)
    U = _lam_seed_field(sm.xyz)
    R1, J1 = solver.residual_jacobian(U)
    R2, J2 = solver.residual_jacobian(U)
    assert np.array_equal(R1, R2)
    assert np.array_equal(J1.data, J2.data)


# ---------------------------------------------------------------------------
# Laminar flat plate: Newton convergence + physics smoke
# ---------------------------------------------------------------------------

def test_laminar_plate_newton():
    sm = _plate_mesh(nx=12, nz=5)
    n = sm.xyz.shape[0]
    flags = np.zeros(n, dtype=np.int64)
    solver = _make_solver(sm, flags)
    U0 = _lam_seed_field(sm.xyz)
    U, info = solver.solve(U0, tol=1.0e-9, max_iter=60)
    assert info["converged"]
    assert info["final_residual"] < 1.0e-9 * max(
        info["residual_history"][0], 1.0
    )
    # inflow Dirichlet retained
    inflow = np.abs(sm.xyz[:, 0] - 0.2) < 1.0e-12
    st_bc = C.blasius_seed(0.2, q=Q0, rho=RHO0, mu=MU0)
    assert np.allclose(U[inflow], st_bc, atol=1.0e-12)
    # laminar pin: Ctau rows pinned at CTAU_LAM everywhere outside inflow
    assert np.allclose(U[~inflow, 4], C.CTAU_LAM, atol=1.0e-12)
    assert np.allclose(U[~inflow, 5], C.CTAU_LAM, atol=1.0e-12)
    # quasi-2-D lock: crossflow unknowns remain at machine zero
    interior = ~inflow
    assert np.max(np.abs(U[interior, 2])) < 1.0e-12
    assert np.max(np.abs(U[interior, 3])) < 1.0e-12
    # family-consistent H in the interior (loose smoke; GV1.1(a) is the gate)
    outs = np.empty((n, C.N_OUT))
    douts = np.empty((n, C.N_OUT, 6))
    q = np.full(n, Q0)
    rho = np.full(n, RHO0)
    mu = np.full(n, MU0)
    mach = np.zeros(n)
    C.closure_all(U, q, rho, mu, mach, flags, C.C_L_DEFAULT, outs, douts)
    H = outs[interior, C.OUT_H1]
    assert np.all(H > 2.55) and np.all(H < 2.75)
    # delta* grows downstream (sqrt-x-like)
    ds1 = outs[:, C.OUT_DS1]
    x_line = np.unique(sm.xyz[:, 0])
    ds_line = [ds1[np.abs(sm.xyz[:, 0] - xx) < 1.0e-12].mean() for xx in x_line]
    assert np.all(np.diff(ds_line) > 0.0)


# ---------------------------------------------------------------------------
# Turbulent flat plate: convergence smoke
# ---------------------------------------------------------------------------

def test_turbulent_plate_newton():
    sm = _plate_mesh(nx=12, nz=5)
    n = sm.xyz.shape[0]
    flags = np.ones(n, dtype=np.int64)
    solver = _make_solver(sm, flags)
    U0 = np.zeros((n, 6))
    for i in range(n):
        U0[i] = _turb_seed_state(max(sm.xyz[i, 0], 1.0e-3))
    U, info = solver.solve(U0, tol=1.0e-9, max_iter=80)
    assert info["converged"]
    outs = np.empty((n, C.N_OUT))
    douts = np.empty((n, C.N_OUT, 6))
    q = np.full(n, Q0)
    rho = np.full(n, RHO0)
    mu = np.full(n, MU0)
    mach = np.zeros(n)
    C.closure_all(U, q, rho, mu, mach, flags, C.C_L_DEFAULT, outs, douts)
    cf = outs[:, C.OUT_CF1]
    H = outs[:, C.OUT_H1]
    assert np.all(cf > 0.0)
    assert np.all(H > 1.2) and np.all(H < 2.0)
    # quasi-2-D lock also in the turbulent branch
    assert np.max(np.abs(U[:, 2])) < 1.0e-10
    assert np.max(np.abs(U[:, 3])) < 1.0e-10
    assert np.max(np.abs(U[:, 5])) < 1.0e-10
