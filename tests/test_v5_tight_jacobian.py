"""Track V V5 Stage 1 -- augmented-Newton fixed operators and the J_phi,BL
FD gate (binding: docs/roadmap/track_v.md GV5.1 + the 2026-07-22
pre-registered FD note; design: cases/analysis/v5_tight_coupling/
PRE_REGISTRATION.md; module under test: pyfp3d/viscous/tight.py).

Covers, on the 2.5-D NACA0012 coarse wake-cut strip (M0.5, alpha 2, Re
3.0e6, x_tr 0.05/0.05 -- the GV3.1 configuration, coarse recorded):

  - W vs assemble_transpiration_rhs, L vs transpiration_from_delta_star,
    S/P scatter/pin structure, D vs the closure packet's OUT_DS1 row --
    all to machine precision (the chain factors are EXACT linear pieces,
    not approximations);
  - G (the phi -> u_e recovery operator) vs edge_velocity_per_zone on a
    random phi: quadratic rows and the two-zone (LE-band) path both to
    machine precision;
  - the J_phi,BL FD gate (the B19/B31/V1 pattern, tests/test_v1_ibl3.py:
    84-109, tests/test_v2_newton_rhs_channel.py:141-180): state point =
    the loose-loop k=1 state (inviscid-converged phi + ONE IBL solve,
    mirroring coupling.py:689-806); reference = central FD of
    NewtonWorkspace.eval_residual with the external_rhs attribute set
    post-construction (the delta* = 0 bit-identity with None checked
    first); epsilon ladder 1e-5/1e-6/1e-7, 4 random max-normalized
    directions, tolerance < 1e-5 (the PRE_REGISTRATION FD protocol;
    unmasked rows expected at 1e-6-1e-8 per the P8/P14 precedent).

Runs in both lanes: default JIT and PYFP3D_NOJIT=1. Under NOJIT the two
k=1-fixture tests (the external_rhs roundtrip and the J_phi,BL FD gate)
are skipped -- the pure-Python k=1 fixture (FP + IBL solves) is JIT-lane
only, the tests/test_v3_coupling.py:204 precedent; the cheap operator
unit tests (W/S/P/L/D/G, no k=1 state) still run.
"""

import os

import numpy as np
import pytest

from pyfp3d.viscous import closures as C
from pyfp3d.viscous.tight import (
    assemble_j_phi_bl,
    closure_ds1_jacobian,
    edge_velocity_operator,
    pin_mask_matrix,
    surface_divergence_delta_operator,
    surface_scatter_matrix,
    wall_load_matrix,
)
from pyfp3d.viscous.transpiration import (
    assemble_transpiration_rhs,
    edge_velocity_per_zone,
    transpiration_from_delta_star,
)
from tests.v5_state import (
    M_CRIT,
    M_CAP,
    RHO_FLOOR,
    UPWIND_C,
    build_k1_state,
    build_naca_case,
)

NOJIT = os.environ.get("PYFP3D_NOJIT", "0") == "1"


@pytest.fixture(scope="module")
def naca_case():
    """The 2.5-D coarse strip + IBL case wiring (tests/v5_state.py)."""
    return build_naca_case()


@pytest.fixture(scope="module")
def k1_state(naca_case):
    """The loose-loop k=1 state point (tests/v5_state.py; shared with the
    Stage-2 gate so both probe the same pre-registered state)."""
    return build_k1_state(naca_case)


def _rel_err(got, ref, floor=1e-30):
    return float(np.max(np.abs(got - ref)) / max(np.max(np.abs(ref)), floor))


# ---------------------------------------------------------------------------
# W: the wall-load matrix
# ---------------------------------------------------------------------------


def test_wall_load_matrix_matches_vector_assembly(naca_case):
    """W @ m_vol == assemble_transpiration_rhs(nodes, wall_faces, m_vol)
    to machine precision (transpiration.py:66-114; the quadrature reads
    diagonal a/6, off-diagonal a/12 per vertex pair)."""
    mc, _, _, case = naca_case
    W = wall_load_matrix(mc.nodes, case.wall_faces)
    rng = np.random.default_rng(7)
    for _ in range(3):
        m_vol = rng.standard_normal(len(mc.nodes))
        ref = assemble_transpiration_rhs(mc.nodes, case.wall_faces, m_vol)
        err = _rel_err(W @ m_vol, ref)
        assert err < 1e-13, f"W mismatch: {err:.3e}"
    # the GV2.1(b) discipline: zero in, exact zero out
    z = np.zeros(len(mc.nodes))
    assert np.array_equal(W @ z, assemble_transpiration_rhs(
        mc.nodes, case.wall_faces, z))


# ---------------------------------------------------------------------------
# S / P: scatter and pin mask
# ---------------------------------------------------------------------------


def test_scatter_and_pin_mask(naca_case):
    """S @ m_surf reproduces the coupling.py:832-833 scatter; P is the
    coupling.py:822-831 diagonal 0/1 mask (identity at None)."""
    _, _, _, case = naca_case
    sm = case.sm
    n_cut = len(case.nodes)
    S = surface_scatter_matrix(sm, n_cut)
    rng = np.random.default_rng(8)
    m_surf = rng.standard_normal(sm.n_node)
    ref = np.zeros(n_cut)
    ref[sm.volume_node_of] = m_surf
    assert np.array_equal(S @ m_surf, ref)

    P = pin_mask_matrix(sm.n_node, None)
    assert np.array_equal(P @ m_surf, m_surf)
    pin = np.zeros(sm.n_node, dtype=bool)
    pin[::7] = True
    Pp = pin_mask_matrix(sm.n_node, pin)
    masked = Pp @ m_surf
    assert np.all(masked[pin] == 0.0)
    assert np.array_equal(masked[~pin], m_surf[~pin])


# ---------------------------------------------------------------------------
# L: the surface-divergence operator at fixed (rho_e, u_e)
# ---------------------------------------------------------------------------


def test_surface_divergence_operator_matches(naca_case):
    """L @ ds == transpiration_from_delta_star(sm, rho_e, ue, ds) to
    machine precision (transpiration.py:117-174)."""
    _, _, _, case = naca_case
    sm = case.sm
    rng = np.random.default_rng(9)
    rho_e = 0.5 + rng.random(sm.n_node)
    ue = rng.standard_normal((sm.n_node, 3))
    L = surface_divergence_delta_operator(sm, rho_e, ue)
    for _ in range(3):
        ds = rng.standard_normal(sm.n_node)
        ref = transpiration_from_delta_star(sm, rho_e, ue, ds)
        err = _rel_err(L @ ds, ref)
        assert err < 1e-13, f"L mismatch: {err:.3e}"


# ---------------------------------------------------------------------------
# D: the closure delta* derivative row
# ---------------------------------------------------------------------------


def test_closure_ds1_jacobian_rows(naca_case):
    """D's row i is douts[i, OUT_DS1, :] against the flat 6-per-node
    state layout (closures.py:870)."""
    _, _, _, case = naca_case
    sm = case.sm
    rng = np.random.default_rng(10)
    douts = rng.standard_normal((sm.n_node, C.N_OUT, 6))
    D = closure_ds1_jacobian(douts)
    assert D.shape == (sm.n_node, 6 * sm.n_node)
    v = rng.standard_normal(6 * sm.n_node)
    ref = np.einsum(
        "ij,ij->i", douts[:, C.OUT_DS1, :], v.reshape(sm.n_node, 6)
    )
    err = _rel_err(D @ v, ref)
    assert err < 1e-14, f"D mismatch: {err:.3e}"


# ---------------------------------------------------------------------------
# G: the phi -> u_e recovery operator (both zone paths)
# ---------------------------------------------------------------------------


def test_edge_velocity_operator_quadratic_rows(naca_case):
    """G with no zone mask == wall_tangential_gradient_quadratic /
    edge_velocity_per_zone (all-quadratic default) on a random phi."""
    mc, _, _, case = naca_case
    sm = case.sm
    G = edge_velocity_operator(mc.nodes, case.wall_faces, sm.volume_node_of)
    assert G.shape == (3 * sm.n_node, len(mc.nodes))
    rng = np.random.default_rng(13)
    for _ in range(2):
        phi = rng.standard_normal(len(mc.nodes))
        got = (G @ phi).reshape(sm.n_node, 3)
        ref = edge_velocity_per_zone(
            mc.nodes, case.wall_faces, phi
        )[sm.volume_node_of]
        err = _rel_err(got, ref)
        assert err < 1e-11, f"G quadratic-row mismatch: {err:.3e}"


def test_edge_velocity_operator_two_zone_rows(naca_case):
    """G with the A4 LE-band mask == edge_velocity_per_zone's two-zone
    path (quadratic elsewhere, P1 + crease-gated smoothing + nodal
    average in the band) on a random phi."""
    mc, _, cfg, case = naca_case
    sm = case.sm
    n_cut = len(mc.nodes)
    le_mask_vol = np.zeros(n_cut, dtype=bool)
    le_mask_vol[sm.volume_node_of[case.le_band_surf]] = True
    assert case.le_band_surf.any() and (~case.le_band_surf).any()
    G = edge_velocity_operator(
        mc.nodes,
        case.wall_faces,
        sm.volume_node_of,
        elements=case.elements,
        le_band_mask=le_mask_vol,
        n_smooth_passes=cfg.n_smooth_passes,
    )
    rng = np.random.default_rng(14)
    for _ in range(2):
        phi = rng.standard_normal(len(mc.nodes))
        got = (G @ phi).reshape(sm.n_node, 3)
        ref = edge_velocity_per_zone(
            mc.nodes,
            case.wall_faces,
            phi,
            elements=case.elements,
            le_band_mask=le_mask_vol,
            n_smooth_passes=cfg.n_smooth_passes,
        )[sm.volume_node_of]
        err = _rel_err(got, ref)
        assert err < 1e-11, f"G two-zone mismatch: {err:.3e}"


# ---------------------------------------------------------------------------
# The J_phi,BL FD gate
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    NOJIT, reason="the k=1 fixture (FP + IBL solves) is JIT-lane only"
)
def test_external_rhs_attribute_roundtrip(k1_state):
    """Precondition of the FD protocol: NewtonWorkspace.external_rhs is a
    plain attribute read per eval_residual call (solve/newton.py:197-203,
    372-377), so post-construction assignment works; the exact-zero
    transpiration load (delta* = 0) is then bit-identical to the
    channel-absent None (the GV2.1(b) discipline through the setter)."""
    st = k1_state
    ws, mc, case = st["ws"], st["mc"], st["case"]
    args = (UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
    ws.external_rhs = assemble_transpiration_rhs(
        mc.nodes, case.wall_faces, np.zeros(st["n_cut"])
    )
    R_z, F_z, _ = ws.eval_residual(st["phi_free"], st["gamma"], *args)
    ws.external_rhs = None
    R_n, F_n, _ = ws.eval_residual(st["phi_free"], st["gamma"], *args)
    assert np.array_equal(R_z, R_n)
    assert np.array_equal(F_z, F_n)


@pytest.mark.skipif(
    NOJIT, reason="the k=1 fixture (FP + IBL solves) is JIT-lane only"
)
def test_j_phi_bl_fd_random_directions(k1_state):
    """GV5.1 Stage-1 gate: J_phi,BL = -(T^T W S P L D)[free, :] vs central
    FD of eval_residual through the external_rhs channel, epsilon ladder
    1e-5/1e-6/1e-7, 4 max-normalized random directions, tolerance < 1e-5
    (the PRE_REGISTRATION FD protocol). The FD reference perturbs U only;
    (phi, gamma, q, rho_e, mach_e) stay frozen -- the loose-loop
    under-relaxation and the max(ds, 0) floor are DROPPED from the tight
    chain (PRE_REGISTRATION decision 2), the closure's DELTA_MIN floor
    masks its own derivatives consistently."""
    st = k1_state
    mc, cfg, case, sm = st["mc"], st["cfg"], st["case"], st["sm"]
    ws = st["ws"]
    n_cut, n_s = st["n_cut"], sm.n_node
    rho_e, ue_surf, q, mach_e = (
        st["rho_e"], st["ue_surf"], st["q"], st["mach_e"])
    mu_arr = np.full(n_s, st["mu"])
    args = (UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)

    W = wall_load_matrix(mc.nodes, case.wall_faces)
    S = surface_scatter_matrix(sm, n_cut)
    P = pin_mask_matrix(n_s, case.outflow_pin_surf)  # airfoil: identity
    L = surface_divergence_delta_operator(sm, rho_e, ue_surf)
    D = closure_ds1_jacobian(st["douts"])
    J = assemble_j_phi_bl(ws, W, S, P, L, D)
    assert J.shape == (ws.n_free, 6 * n_s)

    def residual_at(U_pert):
        """R_free(phi, gamma; U) through the production residual path:
        closure delta* -> m_surf -> scatter -> wall RHS -> eval_residual
        with external_rhs set."""
        outs = np.empty((n_s, C.N_OUT), dtype=np.float64)
        douts = np.empty((n_s, C.N_OUT, 6), dtype=np.float64)
        douts_e = np.empty((n_s, C.N_OUT, 2), dtype=np.float64)
        C.closure_all(
            np.ascontiguousarray(U_pert), q, rho_e, mu_arr, mach_e,
            case.turbulent_flags, C.C_L_DEFAULT, outs, douts, douts_e,
        )
        m_surf = transpiration_from_delta_star(
            sm, rho_e, ue_surf, outs[:, C.OUT_DS1]
        )
        ws.external_rhs = assemble_transpiration_rhs(
            mc.nodes, case.wall_faces, S @ m_surf
        )
        R_free, _, _ = ws.eval_residual(st["phi_free"], st["gamma"], *args)
        return R_free

    rng = np.random.default_rng(51)
    worst_sweet = 0.0
    try:
        for d in range(4):
            v = rng.standard_normal((n_s, 6))
            v /= np.max(np.abs(v))
            an = J @ v.ravel()
            scale = max(float(np.max(np.abs(an))), 1e-12)
            errs = []
            for eps in (1e-5, 1e-6, 1e-7):
                fd = (residual_at(st["U"] + eps * v)
                      - residual_at(st["U"] - eps * v)) / (2.0 * eps)
                errs.append(float(np.max(np.abs(fd - an))) / scale)
            print(f"J_phi,BL FD direction {d}: eps=1e-5/1e-6/1e-7 errs = "
                  f"{errs[0]:.3e} / {errs[1]:.3e} / {errs[2]:.3e}")
            worst_sweet = max(worst_sweet, min(errs))
    finally:
        ws.external_rhs = None
    print(f"J_phi,BL FD worst sweet-spot: {worst_sweet:.3e} "
          f"(ibl converged={st['ibl_info'].get('converged')}, "
          f"final_residual={st['ibl_info'].get('final_residual'):.3e})")
    assert worst_sweet < 1e-5, f"J_phi,BL FD mismatch: {worst_sweet:.3e}"
