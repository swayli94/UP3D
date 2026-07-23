"""Track V V5 Stage 2 -- the J_BL,phi block (binding: docs/roadmap/
track_v.md GV5.1 + the 2026-07-22 pre-registered FD note; design:
cases/analysis/v5_tight_coupling/PRE_REGISTRATION.md; modules under test:
pyfp3d/viscous/ibl3.py (residual_edge_jacobian) and pyfp3d/viscous/
tight.py (edge_data_jacobian, assemble_j_bl_phi)).

The phi -> IBL-residual chain at the frozen k=1 state U:

    phi --G--> u_e --D_ue--> (q, rho, mu, M, u_hat) --J_e--> R_BL

with J_e the IBL3 edge-data Jacobian (state fixed, veps/veps_s FROZEN by
decision 5 -- every FD gate below measures the omission by running the
reference twice, veps frozen / naturally recomputed), D_ue the per-node
isentropic packet Jacobian and G the linear recovery operator.

Covers:

  - the closure packet's dout_e (d/d re_d, d/d e_prime) through the J2
    chain of _nodal_fluxes (re_d = rho q delta/mu, e_prime = r(g-1)/2 M^2)
    vs central FD of closure_scalar in (q, rho, mu, M): laminar,
    turbulent, M = 0 and RE_D-floored states (the floor pins re_d, so the
    masked column is FD-consistent);
  - update_edge_data: bit-identity with a fresh construction;
  - J_e D_ue vs central FD of solver.residual through update_edge_data
    with the full isentropic chain (q -> rho, M recomputed), 4 random
    max-normalized directions, epsilon ladder 1e-5/1e-6/1e-7;
  - the full J_BL,phi = J_e D_ue G vs central FD through
    edge_velocity_per_zone + the isentropic chain, same protocol.

FD protocol mirrors the Stage-1 gate (tests/test_v5_tight_jacobian.py):
per-direction errors normalized by max|analytic|, sweet spot = min over
the ladder, tolerance 1e-5 (packet test: 1e-6). Nodes with q <= 1e-12
(the _frames fallback; D_ue is degenerate there by construction) have
their rows -- and their element neighbours' rows, which the fallback
column feeds -- masked out of the comparison (count printed, bounded).

Runs in both lanes: default JIT and PYFP3D_NOJIT=1. Under NOJIT the three
k=1-fixture tests (bit-identity, J_e.D_ue and J_BL,phi FD gates) are
skipped -- the pure-Python k=1 fixture (FP + IBL solves) is JIT-lane
only, the tests/test_v3_coupling.py:204 precedent; the closure-packet
FD tests (no k=1 state) still run.
"""

import os

import numpy as np
import pytest

from pyfp3d.physics.isentropic import density_field, mach_squared_field
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.coupling import _turb_seed
from pyfp3d.viscous.ibl3 import IBL3Solver
from pyfp3d.viscous.tight import (
    assemble_j_bl_phi,
    edge_data_jacobian,
    edge_velocity_operator,
)
from pyfp3d.viscous.transpiration import edge_velocity_per_zone
from tests.v5_state import build_k1_state, build_naca_case

NOJIT = os.environ.get("PYFP3D_NOJIT", "0") == "1"
K1_JIT_ONLY = pytest.mark.skipif(
    NOJIT, reason="the k=1 fixture (FP + IBL solves) is JIT-lane only"
)


@pytest.fixture(scope="module")
def naca_case():
    return build_naca_case()


@pytest.fixture(scope="module")
def k1_state(naca_case):
    return build_k1_state(naca_case)


# ---------------------------------------------------------------------------
# (a) the closure packet's edge derivatives (dout_e through the J2 chain)
# ---------------------------------------------------------------------------


def _j2_chain(delta, q, rho, mu, mach):
    """d(re_d, e_prime)/d(q, rho, mu, M) mirroring ibl3._nodal_fluxes
    (ibl3.py:264-284; the DELTA_MIN / RE_D_MIN floor logic of
    closure_node, closures.py:490-508)."""
    de_f = max(delta, C.DELTA_MIN)
    dre = np.zeros(3)  # d re_d / d(q, rho, mu); re_d pinned at the floor
    if rho * q * de_f / mu >= C.RE_D_MIN:
        dre[:] = (
            rho * de_f / mu,
            q * de_f / mu,
            -rho * q * de_f / (mu * mu),
        )
    dep = C.RECOVERY_R * (C.GAMMA_AIR - 1.0) * mach  # d e_prime / dM
    return dre, dep


def _packet_fd(state, q, rho, mu, mach, turbulent, h):
    """Central FD of closure_scalar outputs w.r.t. (q, rho, mu, M); the
    steps are scaled per parameter (dmu = h*mu in particular, so the
    ladder never drives mu through zero)."""
    steps = (h * max(q, 0.5), h * max(rho, 0.5), h * mu,
             h * max(mach, 0.5))
    cols = []
    for p in range(4):
        d = [0.0, 0.0, 0.0, 0.0]
        d[p] = steps[p]
        op, _, _ = C.closure_scalar(
            state, q=q + d[0], rho=rho + d[1], mu=mu + d[2],
            mach=mach + d[3], turbulent=turbulent,
        )
        om, _, _ = C.closure_scalar(
            state, q=q - d[0], rho=rho - d[1], mu=mu - d[2],
            mach=mach - d[3], turbulent=turbulent,
        )
        cols.append((op - om) / (2.0 * steps[p]))
    return np.stack(cols)  # (4, N_OUT)


@pytest.mark.parametrize(
    "kind,state,q,rho,mu,mach,turbulent",
    [
        # laminar Blasius-ish, compressible
        ("lam", C.blasius_seed(0.5, q=1.0, rho=1.0, mu=1.0e-5),
         0.9, 1.0, 1.0e-5, 0.3, False),
        # turbulent, compressible
        ("turb", _turb_seed(0.5, 1.0, 1.0, 1.0e-5),
         1.1, 1.0, 1.0e-5, 0.5, True),
        # M = 0: e_prime = 0 and d e_prime/dM = 0 (the symmetric FD is
        # exactly zero, the analytic M column likewise)
        ("M0", C.blasius_seed(0.5, q=1.0, rho=1.0, mu=1.0e-5),
         0.9, 1.0, 1.0e-5, 0.0, False),
        # RE_D-floored: re_d = 1 * 1e-3 * 1e-4 / 1.0 = 1e-7 << RE_D_MIN,
        # so the re_d column is pinned (zero) and FD-consistent
        ("floored", np.array([1.0e-4, 8.0, 0.0, 0.0, C.CTAU_LAM, 0.0]),
         1.0e-3, 1.0, 1.0, 0.3, False),
    ],
)
def test_closure_packet_edge_fd(kind, state, q, rho, mu, mach, turbulent):
    """dout_e through the J2 chain vs central FD of closure_scalar in
    (q, rho, mu, M), ladder 1e-5/1e-6/1e-7, tol 1e-6 (scaled)."""
    out0, _, dout_e = C.closure_scalar(
        state, q=q, rho=rho, mu=mu, mach=mach, turbulent=turbulent
    )
    dre, dep = _j2_chain(state[0], q, rho, mu, mach)
    an = np.stack(
        [dout_e[:, 0] * dre[0], dout_e[:, 0] * dre[1],
         dout_e[:, 0] * dre[2], dout_e[:, 1] * dep]
    )  # (4, N_OUT): columns (q, rho, mu, M)
    scale = max(float(np.max(np.abs(an))), 1.0)
    worst = 0.0
    errs = []
    for h in (1e-5, 1e-6, 1e-7):
        fd = _packet_fd(state, q, rho, mu, mach, turbulent, h)
        err = float(np.max(np.abs(fd - an))) / scale
        errs.append(err)
        worst = max(worst, err)
    print(f"packet edge FD [{kind}]: ladder errs = "
          f"{errs[0]:.3e} / {errs[1]:.3e} / {errs[2]:.3e}")
    assert min(errs) < 1e-6, f"packet edge FD [{kind}]: {min(errs):.3e}"


# ---------------------------------------------------------------------------
# (b) update_edge_data bit-identity with fresh construction
# ---------------------------------------------------------------------------


@K1_JIT_ONLY
def test_update_edge_data_bit_identical(k1_state):
    """solver.update_edge_data(u2, ...) -> residual_edge_jacobian must be
    bit-for-bit identical to a fresh IBL3Solver built on u2 (frames,
    veps/veps_s and the J_e tables all rebuilt from the new data)."""
    st = k1_state
    solver, sm, case, cfg = st["solver"], st["sm"], st["case"], st["cfg"]
    U = st["U"]
    mu_arr = np.full(sm.n_node, st["mu"])
    rng = np.random.default_rng(31)
    ue2 = st["ue_surf"] + 0.05 * rng.standard_normal(st["ue_surf"].shape)
    q2 = np.sum(ue2 ** 2, axis=1)
    rho2 = density_field(q2, cfg.m_inf, cfg.gamma_air)
    mach2 = np.sqrt(mach_squared_field(q2, cfg.m_inf, cfg.gamma_air))
    try:
        solver.update_edge_data(ue2, rho2, mu_arr, mach2)
        R1, J1 = solver.residual_edge_jacobian(U)
    finally:
        solver.update_edge_data(
            st["ue_surf"], st["rho_e"], mu_arr, st["mach_e"])
    fresh = IBL3Solver(
        sm, ue2, rho2, mu_arr, mach2, case.turbulent_flags,
        st["inflow_mask"], st["inflow_state"],
        eps_diff=cfg.eps_diff, eps_diff_s=cfg.eps_diff_s,
    )
    R2, J2 = fresh.residual_edge_jacobian(U)
    assert np.array_equal(R1, R2)
    assert np.array_equal(J1.indptr, J2.indptr)
    assert np.array_equal(J1.indices, J2.indices)
    assert np.array_equal(J1.data, J2.data)


# ---------------------------------------------------------------------------
# shared helpers for the (c)/(d) FD gates
# ---------------------------------------------------------------------------


def _fallback_row_mask(st):
    """True on rows kept in the FD comparison: drop the 6 rows of every
    q <= 1e-12 fallback node (the _frames fallback; D_ue degenerate) and
    of every node sharing an element with one."""
    sm, q = st["sm"], st["q"]
    bad = np.where(q <= 1.0e-12)[0]
    rows = np.ones(6 * sm.n_node, dtype=bool)
    if bad.size:
        drop = set(int(i) for i in bad)
        for e in range(len(sm.triangles)):
            if np.isin(sm.triangles[e], bad).any():
                drop.update(int(i) for i in sm.triangles[e])
        for i in drop:
            rows[6 * i: 6 * i + 6] = False
    return rows, bad


def _resid_at_ue(st, ue, freeze_veps):
    """R_BL(U) at perturbed edge data through the production path:
    isentropic rho/M recomputed from the perturbed q (coupling.py:723-726),
    update_edge_data in place, base state restored before returning."""
    solver, cfg = st["solver"], st["cfg"]
    mu_arr = np.full(st["sm"].n_node, st["mu"])
    veps0, veps_s0 = st["_veps0"]
    q2 = np.sum(ue ** 2, axis=1)
    solver.update_edge_data(
        ue,
        density_field(q2, cfg.m_inf, cfg.gamma_air),
        mu_arr,
        np.sqrt(mach_squared_field(q2, cfg.m_inf, cfg.gamma_air)),
    )
    if freeze_veps:
        solver._veps, solver._veps_s = veps0, veps_s0
    R = solver.residual(st["U"]).ravel()
    solver.update_edge_data(
        st["ue_surf"], st["rho_e"], mu_arr, st["mach_e"])
    solver._veps, solver._veps_s = veps0, veps_s0
    return R


@pytest.fixture(scope="module")
def j_e_state(k1_state):
    """Base J_e / D_ue / row mask / frozen veps for the (c)/(d) gates."""
    st = k1_state
    st["_veps0"] = (st["solver"]._veps, st["solver"]._veps_s)
    _, J_e = st["solver"].residual_edge_jacobian(st["U"])
    D_ue = edge_data_jacobian(st["ue_surf"], st["cfg"].m_inf,
                              st["cfg"].gamma_air)
    row_mask, bad = _fallback_row_mask(st)
    return st, J_e, D_ue, row_mask, bad


# ---------------------------------------------------------------------------
# (c) J_e D_ue vs FD through update_edge_data
# ---------------------------------------------------------------------------


@K1_JIT_ONLY
def test_j_e_d_ue_fd_random_directions(j_e_state):
    """J_e D_ue (the u_e -> R_BL block) vs central FD of solver.residual
    through update_edge_data with the full isentropic chain; 4 random
    max-normalized directions, ladder 1e-5/1e-6/1e-7, tol 1e-5. The FD
    reference runs with veps/veps_s frozen at the base values (J_e
    freezes them, decision 5); the unfrozen contribution is printed."""
    st, J_e, D_ue, row_mask, bad = j_e_state
    n_s = st["sm"].n_node
    print(f"fallback nodes (q <= 1e-12): {len(bad)}; "
          f"masked rows: {int((~row_mask).sum())} / {6 * n_s}")
    assert (~row_mask).sum() <= 0.02 * 6 * n_s
    J = (J_e @ D_ue).tocsr()
    rng = np.random.default_rng(71)
    worst_sweet = 0.0
    for d in range(4):
        v = rng.standard_normal(3 * n_s)
        v /= np.max(np.abs(v))
        an = (J @ v)[row_mask]
        scale = max(float(np.max(np.abs(an))), 1e-12)
        errs = []
        for eps in (1e-5, 1e-6, 1e-7):
            du = (eps * v).reshape(n_s, 3)
            fd = (_resid_at_ue(st, st["ue_surf"] + du, True)
                  - _resid_at_ue(st, st["ue_surf"] - du, True)) / (2.0 * eps)
            errs.append(float(np.max(np.abs(fd[row_mask] - an))) / scale)
        # the frozen-veps omission (decision 5), measured at eps = 1e-6
        du = (1e-6 * v).reshape(n_s, 3)
        fd_unf = (_resid_at_ue(st, st["ue_surf"] + du, False)
                  - _resid_at_ue(st, st["ue_surf"] - du, False)) / 2.0e-6
        fd_fr = (_resid_at_ue(st, st["ue_surf"] + du, True)
                 - _resid_at_ue(st, st["ue_surf"] - du, True)) / 2.0e-6
        veps_eff = float(np.max(np.abs(fd_unf[row_mask] - fd_fr[row_mask])))
        print(f"J_e.D_ue FD direction {d}: eps=1e-5/1e-6/1e-7 errs = "
              f"{errs[0]:.3e} / {errs[1]:.3e} / {errs[2]:.3e}   "
              f"veps-freeze contribution {veps_eff:.3e} (scaled "
              f"{veps_eff / scale:.3e})")
        worst_sweet = max(worst_sweet, min(errs))
    print(f"J_e.D_ue FD worst sweet-spot: {worst_sweet:.3e}")
    assert worst_sweet < 1e-5, f"J_e.D_ue FD mismatch: {worst_sweet:.3e}"


# ---------------------------------------------------------------------------
# (d) the full J_BL,phi = J_e D_ue G vs FD through edge_velocity_per_zone
# ---------------------------------------------------------------------------


@K1_JIT_ONLY
def test_j_bl_phi_fd_random_directions(j_e_state):
    """J_BL,phi = J_e D_ue G vs central FD of R_BL(U) through
    edge_velocity_per_zone + the isentropic chain (the production phi ->
    u_e -> edge-data -> residual path), same FD protocol as (c)."""
    st, J_e, D_ue, row_mask, bad = j_e_state
    mc, case, sm, cfg = st["mc"], st["case"], st["sm"], st["cfg"]
    n_s = sm.n_node
    G = edge_velocity_operator(
        mc.nodes, case.wall_faces, sm.volume_node_of,
        elements=case.elements, le_band_mask=st["le_mask_vol"],
        n_smooth_passes=cfg.n_smooth_passes,
    )
    J = assemble_j_bl_phi(J_e, D_ue, G)
    assert J.shape == (6 * n_s, st["n_cut"])

    def ue_at_phi(phi):
        ue_vol = edge_velocity_per_zone(
            mc.nodes, case.wall_faces, phi,
            elements=case.elements, le_band_mask=st["le_mask_vol"],
            n_smooth_passes=cfg.n_smooth_passes,
        )
        return ue_vol[sm.volume_node_of]

    rng = np.random.default_rng(83)
    worst_sweet = 0.0
    for d in range(4):
        v = rng.standard_normal(st["n_cut"])
        v /= np.max(np.abs(v))
        an = (J @ v)[row_mask]
        scale = max(float(np.max(np.abs(an))), 1e-12)
        errs = []
        for eps in (1e-5, 1e-6, 1e-7):
            fd = (_resid_at_ue(st, ue_at_phi(st["phi"] + eps * v), True)
                  - _resid_at_ue(st, ue_at_phi(st["phi"] - eps * v), True)
                  ) / (2.0 * eps)
            errs.append(float(np.max(np.abs(fd[row_mask] - an))) / scale)
        du_p = ue_at_phi(st["phi"] + 1e-6 * v)
        du_m = ue_at_phi(st["phi"] - 1e-6 * v)
        fd_unf = (_resid_at_ue(st, du_p, False)
                  - _resid_at_ue(st, du_m, False)) / 2.0e-6
        fd_fr = (_resid_at_ue(st, du_p, True)
                 - _resid_at_ue(st, du_m, True)) / 2.0e-6
        veps_eff = float(np.max(np.abs(fd_unf[row_mask] - fd_fr[row_mask])))
        print(f"J_BL,phi FD direction {d}: eps=1e-5/1e-6/1e-7 errs = "
              f"{errs[0]:.3e} / {errs[1]:.3e} / {errs[2]:.3e}   "
              f"veps-freeze contribution {veps_eff:.3e} (scaled "
              f"{veps_eff / scale:.3e})")
        worst_sweet = max(worst_sweet, min(errs))
    print(f"J_BL,phi FD worst sweet-spot: {worst_sweet:.3e}")
    assert worst_sweet < 1e-5, f"J_BL,phi FD mismatch: {worst_sweet:.3e}"
