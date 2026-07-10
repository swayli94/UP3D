"""P8/N2: exact Newton Jacobian assembly (design.md (6.3)) --
finite-difference verification of PicardOperator.assemble_newton_jacobian
against the SHIPPED walk flux at frozen upstream selection.

What is verified is the full assembled dR/dphi:

    J @ delta  vs  central difference of  R(phi) = assemble_residual(
        phi, rho_tilde_sweep(density chain at FROZEN u(e)))

so Term 1 (Picard block), Term 2 (local density sensitivity s_e) and
Term 3 (upstream coupling s_u, the long-range COO block) are all exercised
through the same JVP. The FD oracle freezes the selection exactly like
P7's _frozen_rho_tilde (tests/test_p7_diff_flux.py): differencing the raw
driver residual would re-run the walk and turn selection flips into fake
derivative errors.

Kink protocol inherited from P7 (G7.3): rho_tilde is C0 but not C1 at the
max(nu_e, nu_u) tie and at the switch threshold M^2 = M_c^2; FD probes
straddling that measure-zero locus read branch AVERAGES. Fields are
generic (noise-broken), kink-neighbourhood ELEMENT masks are lifted to
residual ROWS (rho_tilde_e enters R_i only for i in nodes(e)), and the
excluded fraction is asserted small.

Runs under PYFP3D_NOJIT=1 as well. The converged-pocket check on the real
G4.1 field (the G8.1 FD clause) is gated behind PYFP3D_TRANSONIC_GATES=1.
"""

import os

import numpy as np
import pytest

from pyfp3d.kernels.jacobian import PicardOperator
from pyfp3d.kernels.upwind import UpwindOperator, rho_tilde_sweep
from pyfp3d.physics.isentropic import (
    GAMMA,
    density_field,
    mach_number_squared,
)

from .mesh_utils import generate_structured_cube_mesh

M_CRIT = 0.95
UPWIND_C = 1.5
RHO_FLOOR = 0.05
REL_TOL = 1e-6          # G8.1 FD acceptance (same figure as G7.3)

run_gates = pytest.mark.skipif(
    os.environ.get("PYFP3D_TRANSONIC_GATES", "0") != "1",
    reason="converged-pocket Jacobian FD needs the coarse transonic solve; "
           "set PYFP3D_TRANSONIC_GATES=1 for the gate-closure run",
)


def _cube():
    return generate_structured_cube_mesh(n=6, L=1.0)


def _degeneracy_breaker(n_nodes, amp=0.005, seed=42):
    """See test_p7_diff_flux: lattice-symmetric phi parks whole element
    slabs exactly on the max(nu_e, nu_u) tie kink; generic nodal noise
    empties the tie set."""
    rng = np.random.default_rng(seed)
    return amp * rng.standard_normal(n_nodes)


def _frozen_residual(op, phi, m_inf, upstream, rho_floor=RHO_FLOOR):
    """R(phi) through the SHIPPED density chain and flux at FROZEN
    selection: velocities -> density_field -> rho_tilde_sweep (verbatim,
    u(e) held fixed) -> assemble_residual. u_inf = 1 and no speed limiter
    (cube fields stay far below M_cap), matching the Jacobian call."""
    _, q2 = op.velocities(phi)
    q2 = q2.copy()
    rho = density_field(q2, m_inf)
    nu = np.empty_like(q2)
    rt = np.empty_like(q2)
    rho_tilde_sweep(q2, rho, upstream, m_inf, M_CRIT, UPWIND_C, GAMMA,
                    rho_floor, nu, rt)
    return op.assemble_residual(phi, rt).copy()


def _kink_row_mask(elements, q2, upstream, m_inf, n_nodes,
                   kink_guard=3e-5):
    """Rows of R that an FD probe could corrupt: nodes of elements within
    the kink guard of the max(nu_e, nu_u) tie or the M^2 = M_c^2 switch
    threshold (rho_tilde_e reaches R_i only for i in nodes(e)). Same
    element-level criterion as P7's _jvp_check."""
    mc2 = M_CRIT ** 2
    m2 = mach_number_squared(q2, m_inf, GAMMA)
    m2u = m2[upstream]
    nu_e = UPWIND_C * np.maximum(0.0, 1.0 - mc2 / np.maximum(m2, mc2))
    nu_u = UPWIND_C * np.maximum(0.0, 1.0 - mc2 / np.maximum(m2u, mc2))
    active = np.maximum(nu_e, nu_u) > 0.0
    self_up = upstream == np.arange(len(q2))
    near_tie = active & ~self_up & (np.abs(nu_e - nu_u) < kink_guard)
    near_thresh = ~self_up & ((np.abs(m2 - mc2) < kink_guard)
                              | (np.abs(m2u - mc2) < kink_guard))
    bad_elem = near_tie | near_thresh
    row_bad = np.zeros(n_nodes, dtype=bool)
    if bad_elem.any():
        row_bad[np.asarray(elements)[bad_elem].ravel()] = True
    return ~row_bad, int(bad_elem.sum())


def _jacobian_fd_check(phi, m_inf, seed=0, n_dirs=3, eps=1e-6,
                       max_kink_frac=0.02):
    """Directional (JVP) check of the assembled Newton Jacobian against a
    central difference of the frozen-selection residual. Returns
    (max_rel_err, op, upw, s_u, upstream, q2)."""
    nodes, elements = _cube()
    op = PicardOperator(nodes, elements)
    upw = UpwindOperator(nodes, elements, weighted=False)

    grad, q2 = op.velocities(phi)
    grad, q2 = grad.copy(), q2.copy()
    rho = density_field(q2, m_inf)
    rho_t = upw.rho_tilde(grad, q2, rho, m_inf, UPWIND_C, M_CRIT,
                          rho_floor=RHO_FLOOR).copy()
    s_e, s_u, upstream = upw.rho_tilde_sensitivities(
        grad, q2, rho, m_inf, UPWIND_C, M_CRIT, rho_floor=RHO_FLOOR)
    s_e, s_u, upstream = s_e.copy(), s_u.copy(), upstream.copy()

    J = op.assemble_newton_jacobian(phi, rho_t, s_e, s_u, upstream)

    # Same P7 acceptance: the kink ELEMENT set must be a tiny fraction of
    # the mesh (generic fields keep it near-empty); each kink element then
    # masks its 4 residual rows -- on a small cube a handful of elements
    # is a visible ROW fraction, which is expected and harmless.
    valid, n_bad_elems = _kink_row_mask(elements, q2, upstream, m_inf,
                                        op.n_nodes)
    assert n_bad_elems <= max_kink_frac * op.n_tets, (
        f"{n_bad_elems} elements near a max-kink -- field too "
        "symmetry-degenerate for an FD check (see module docstring)")

    rng = np.random.default_rng(seed)
    max_rel = 0.0
    for _ in range(n_dirs):
        delta = rng.standard_normal(op.n_nodes)
        delta /= np.abs(delta).max()
        jvp = J @ delta
        r_p = _frozen_residual(op, phi + eps * delta, m_inf, upstream)
        r_m = _frozen_residual(op, phi - eps * delta, m_inf, upstream)
        fd = (r_p - r_m) / (2.0 * eps)
        scale = np.abs(fd[valid]).max()
        assert scale > 0.0
        rel = np.abs(jvp[valid] - fd[valid]).max() / scale
        max_rel = max(max_rel, rel)
    return max_rel, op, upw, s_u, upstream, q2


def _regime_counts(q2, upstream, m_inf):
    """(subsonic, accelerating, shock-point) element counts at the frozen
    selection (same classification as test_p7_diff_flux._regimes)."""
    m2 = mach_number_squared(q2, m_inf, GAMMA)
    m2u = m2[upstream]
    mc2 = M_CRIT ** 2
    nu_e = UPWIND_C * np.maximum(0.0, 1.0 - mc2 / np.maximum(m2, mc2))
    nu_u = UPWIND_C * np.maximum(0.0, 1.0 - mc2 / np.maximum(m2u, mc2))
    self_up = upstream == np.arange(len(q2))
    subsonic = (~self_up) & (np.maximum(nu_e, nu_u) == 0.0)
    accel = (~self_up) & (nu_e >= nu_u) & (nu_e > 0.0)
    shockpt = (~self_up) & (nu_u > nu_e)
    return int(subsonic.sum()), int(accel.sum()), int(shockpt.sum())


# ------------------------------------------------------------- JVP ladder


def test_full_jacobian_jvp_subsonic():
    """nu == 0 everywhere: J must still carry Term 2 (drho/dq2 != 0) and
    the Term-3 block must be EMPTY (s_u == 0 by construction)."""
    nodes, _ = _cube()
    phi = nodes[:, 0].copy()                     # q2 = 1, subcritical at 0.5
    rel, op, upw, s_u, upstream, q2 = _jacobian_fd_check(phi, m_inf=0.5)
    assert rel < REL_TOL
    assert np.all(s_u == 0.0)
    assert op.n_term3_active == 0


def test_full_jacobian_jvp_accelerating():
    """Supersonic accelerating pocket: nu = nu_e branch, Term 3 active."""
    nodes, _ = _cube()
    x = nodes[:, 0]
    phi = x + 0.3 * x ** 2 + _degeneracy_breaker(len(nodes))
    rel, op, upw, s_u, upstream, q2 = _jacobian_fd_check(phi, m_inf=0.8)
    assert rel < REL_TOL
    _, n_accel, _ = _regime_counts(q2, upstream, 0.8)
    assert n_accel > 0
    assert op.n_term3_active > 0                 # long-range block exercised


def test_full_jacobian_jvp_shockpoint():
    """Decelerating field: nu = nu_u shock-point branch in the Jacobian."""
    nodes, _ = _cube()
    x = nodes[:, 0]
    phi = 1.6 * x - 0.3 * x ** 2 + _degeneracy_breaker(len(nodes))
    rel, op, upw, s_u, upstream, q2 = _jacobian_fd_check(phi, m_inf=0.8)
    assert rel < REL_TOL
    _, _, n_shockpt = _regime_counts(q2, upstream, 0.8)
    assert n_shockpt > 0


def test_full_jacobian_jvp_mixed():
    """All three regimes in one field, 5 random directions."""
    nodes, _ = _cube()
    x = nodes[:, 0]
    phi = (1.25 * x + 0.12 * np.sin(2.0 * np.pi * x)
           + _degeneracy_breaker(len(nodes)))
    rel, op, upw, s_u, upstream, q2 = _jacobian_fd_check(phi, m_inf=0.8,
                                                         n_dirs=5)
    assert rel < REL_TOL
    n_sub, n_accel, n_shockpt = _regime_counts(q2, upstream, 0.8)
    assert n_sub > 0 and n_accel > 0 and n_shockpt > 0


# --------------------------------------------------------- block structure


def test_subsonic_jacobian_shares_picard_pattern():
    """Subsonic: Term 3 empty, so J lives exactly on the Picard CSR
    pattern, and its symmetric Term-1 part equals assemble_matrix."""
    nodes, elements = _cube()
    op = PicardOperator(nodes, elements)
    upw = UpwindOperator(nodes, elements, weighted=False)
    phi = nodes[:, 0].copy()
    grad, q2 = op.velocities(phi)
    grad, q2 = grad.copy(), q2.copy()
    rho = density_field(q2, 0.5)
    rho_t = upw.rho_tilde(grad, q2, rho, 0.5, UPWIND_C, M_CRIT).copy()
    s_e, s_u, upstream = upw.rho_tilde_sensitivities(
        grad, q2, rho, 0.5, UPWIND_C, M_CRIT)
    s_e, s_u, upstream = s_e.copy(), s_u.copy(), upstream.copy()

    J = op.assemble_newton_jacobian(phi, rho_t, s_e, s_u, upstream)
    A = op.assemble_matrix(rho_t)
    assert op.n_term3_active == 0
    assert np.array_equal(J.indptr, A.indptr)
    assert np.array_equal(J.indices, A.indices)
    # J - A is exactly the Term-2 rank-one-per-element block: with s_e all
    # equal to rho'(q2) (q2 uniform here), it is nonzero and NOT symmetric
    # in general, but J must reduce to A when the sensitivities vanish.
    J0 = op.assemble_newton_jacobian(phi, rho_t, np.zeros_like(s_e),
                                     np.zeros_like(s_u), upstream)
    assert np.array_equal(J0.data, A.data)


def test_limiter_mask_gates_the_chain():
    """lim_mask=False on an element zeroes its Term-2 row block (flat
    clamp => zero derivative through q2_limited), and s inputs are not
    mutated."""
    nodes, elements = _cube()
    op = PicardOperator(nodes, elements)
    upw = UpwindOperator(nodes, elements, weighted=False)
    x = nodes[:, 0]
    phi = x + 0.3 * x ** 2 + _degeneracy_breaker(len(nodes))
    grad, q2 = op.velocities(phi)
    grad, q2 = grad.copy(), q2.copy()
    rho = density_field(q2, 0.8)
    rho_t = upw.rho_tilde(grad, q2, rho, 0.8, UPWIND_C, M_CRIT).copy()
    s_e, s_u, upstream = upw.rho_tilde_sensitivities(
        grad, q2, rho, 0.8, UPWIND_C, M_CRIT)
    s_e, s_u, upstream = s_e.copy(), s_u.copy(), upstream.copy()
    s_e_ref, s_u_ref = s_e.copy(), s_u.copy()

    all_limited = np.zeros(op.n_tets, dtype=bool)
    J_lim = op.assemble_newton_jacobian(phi, rho_t, s_e, s_u, upstream,
                                        lim_mask=all_limited)
    A = op.assemble_matrix(rho_t)
    assert op.n_term3_active == 0                # s_u fully masked
    assert np.array_equal(J_lim.data, A.data)    # pure Term 1 remains
    assert np.array_equal(s_e, s_e_ref) and np.array_equal(s_u, s_u_ref)


# ------------------------------------------- frozen-assignment machinery


def test_frozen_sweep_matches_live_at_freeze_state():
    """rho_tilde_frozen at the state the assignment was taken from equals
    the live walk flux BITWISE (the freeze forces the branch max() would
    pick anyway) -- the property that lets the N5 driver freeze without
    perturbing the accepted iterate."""
    from pyfp3d.kernels.upwind import UpwindOperator

    nodes, elements = _cube()
    op = PicardOperator(nodes, elements)
    upw = UpwindOperator(nodes, elements, weighted=False)
    x = nodes[:, 0]
    phi = (1.25 * x + 0.12 * np.sin(2.0 * np.pi * x)
           + _degeneracy_breaker(len(nodes)))
    grad, q2 = op.velocities(phi)
    grad, q2 = grad.copy(), q2.copy()
    rho = density_field(q2, 0.8)
    live = upw.rho_tilde(grad, q2, rho, 0.8, UPWIND_C, M_CRIT).copy()
    up, br = upw.freeze_upwind_state(grad, q2, rho, 0.8, UPWIND_C, M_CRIT)
    frozen = upw.rho_tilde_frozen(q2, rho, up, br, 0.8, UPWIND_C,
                                  M_CRIT).copy()
    assert np.array_equal(live, frozen)
    assert set(np.unique(br)).issubset({0, 1, 2, 3})
    assert (br == 1).sum() > 0 and (br == 2).sum() > 0   # regimes present


def test_frozen_jacobian_jvp():
    """Full-Jacobian JVP against FD of the FROZEN flux residual: with the
    assignment frozen there is no max-tie kink at all (only the sonic
    threshold nu' guard), so the check needs no tie exclusion -- exactly
    why the N5 finish phase converges quadratically."""
    from pyfp3d.kernels.upwind import UpwindOperator, rho_tilde_frozen_sweep

    nodes, elements = _cube()
    op = PicardOperator(nodes, elements)
    upw = UpwindOperator(nodes, elements, weighted=False)
    x = nodes[:, 0]
    phi = (1.25 * x + 0.12 * np.sin(2.0 * np.pi * x)
           + _degeneracy_breaker(len(nodes)))
    grad, q2 = op.velocities(phi)
    grad, q2 = grad.copy(), q2.copy()
    rho = density_field(q2, 0.8)
    rho_t = upw.rho_tilde(grad, q2, rho, 0.8, UPWIND_C, M_CRIT).copy()
    up, br = upw.freeze_upwind_state(grad, q2, rho, 0.8, UPWIND_C, M_CRIT)
    s_e, s_u = upw.rho_tilde_frozen_sensitivities(q2, rho, up, br, 0.8,
                                                  UPWIND_C, M_CRIT)
    J = op.assemble_newton_jacobian(phi, rho_t, s_e.copy(), s_u.copy(), up)

    def _frozen_res(p):
        _, q2p = op.velocities(p)
        q2p = q2p.copy()
        rhop = density_field(q2p, 0.8)
        nu = np.empty_like(q2p)
        rt = np.empty_like(q2p)
        rho_tilde_frozen_sweep(q2p, rhop, up, br, 0.8, M_CRIT, UPWIND_C,
                               GAMMA, RHO_FLOOR, nu, rt)
        return op.assemble_residual(p, rt).copy()

    # only the sonic-threshold guard (nu' -> 0 at M^2 = M_c^2) remains a
    # kink; exclude rows of elements whose active-side M^2 straddles it
    m2 = mach_number_squared(q2, 0.8, GAMMA)
    mc2 = M_CRIT ** 2
    m2_src = np.where(br == 2, m2[up], m2)
    bad = (br != 0) & (br != 3) & (np.abs(m2_src - mc2) < 3e-5)
    row_bad = np.zeros(op.n_nodes, dtype=bool)
    if bad.any():
        row_bad[np.asarray(elements)[bad].ravel()] = True
    valid = ~row_bad

    rng = np.random.default_rng(1)
    eps = 1e-6
    for _ in range(3):
        delta = rng.standard_normal(op.n_nodes)
        delta /= np.abs(delta).max()
        jvp = J @ delta
        fd = (_frozen_res(phi + eps * delta)
              - _frozen_res(phi - eps * delta)) / (2.0 * eps)
        scale = np.abs(fd[valid]).max()
        rel = np.abs(jvp[valid] - fd[valid]).max() / scale
        assert rel < REL_TOL, f"frozen JVP rel err {rel:.3e}"


# ------------------------------------------------------------- guardrails


def test_forward_paths_untouched_by_newton_assembly():
    """assemble_matrix/assemble_residual before and after a Newton
    Jacobian assembly are byte-identical (G4.2 bit-identity is structural:
    the Newton path shares no forward buffer)."""
    nodes, elements = _cube()
    op = PicardOperator(nodes, elements)
    upw = UpwindOperator(nodes, elements, weighted=False)
    x = nodes[:, 0]
    phi = 1.25 * x + _degeneracy_breaker(len(nodes))
    grad, q2 = op.velocities(phi)
    grad, q2 = grad.copy(), q2.copy()
    rho = density_field(q2, 0.82)
    rho_t = upw.rho_tilde(grad, q2, rho, 0.82, UPWIND_C, M_CRIT).copy()

    A_ref = op.assemble_matrix(rho_t)
    R_ref = op.assemble_residual(phi, rho_t).copy()

    s_e, s_u, upstream = upw.rho_tilde_sensitivities(
        grad, q2, rho, 0.82, UPWIND_C, M_CRIT)
    op.assemble_newton_jacobian(phi, rho_t, s_e.copy(), s_u.copy(),
                                upstream.copy())

    rho_t_again = upw.rho_tilde(grad, q2, rho, 0.82, UPWIND_C,
                                M_CRIT).copy()
    assert np.array_equal(rho_t, rho_t_again)
    A_again = op.assemble_matrix(rho_t_again)
    R_again = op.assemble_residual(phi, rho_t_again).copy()
    assert np.array_equal(A_ref.data, A_again.data)
    assert np.array_equal(R_ref, R_again)


# ------------------------------------------- gated: real converged pocket


@run_gates
def test_jacobian_fd_converged_pocket(mesh_dir):
    """G8.1 FD clause on the REAL converged coarse M0.80/alpha1.25 field
    (supersonic pocket with both accelerating and shock-point branches):
    assembled-Jacobian JVP vs frozen-selection FD, rel < 1e-6. The field
    is the NEWTON solution (residual < 1e-9 -- the actual discrete
    solution; the Picard state is a stall artifact, see the roadmap P4
    erratum). Converged-field kink exclusion is tiny."""
    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.mesh.wake_cut import cut_wake
    from pyfp3d.solve.newton import solve_newton_transonic

    mesh = read_mesh(mesh_dir / "naca0012_2.5d" / "coarse.msh")
    mc, wc = cut_wake(mesh)
    r = solve_newton_transonic(
        mc, wc, m_inf=0.80, alpha_deg=1.25, dm=0.025, dm_min=0.003,
        freeze_tol=1e-6,
        newton_kw=dict(freeze_refresh_max=8, precond="direct",
                       n_newton_max=60))
    assert r["converged"]
    phi = np.asarray(r["phi"], dtype=np.float64)

    op = PicardOperator(mc.nodes, mc.elements)
    upw = UpwindOperator(mc.nodes, mc.elements, weighted=False)
    grad, q2 = op.velocities(phi)
    grad, q2 = grad.copy(), q2.copy()
    rho = density_field(q2, 0.80)
    rho_t = upw.rho_tilde(grad, q2, rho, 0.80, UPWIND_C, M_CRIT).copy()
    s_e, s_u, upstream = upw.rho_tilde_sensitivities(
        grad, q2, rho, 0.80, UPWIND_C, M_CRIT)
    s_e, s_u, upstream = s_e.copy(), s_u.copy(), upstream.copy()
    J = op.assemble_newton_jacobian(phi, rho_t, s_e, s_u, upstream)
    assert op.n_term3_active > 0

    valid, _ = _kink_row_mask(mc.elements, q2, upstream, 0.80, op.n_nodes)
    assert (~valid).sum() <= 0.02 * op.n_nodes

    rng = np.random.default_rng(0)
    eps = 1e-6
    for _ in range(3):
        delta = rng.standard_normal(op.n_nodes)
        delta /= np.abs(delta).max()
        jvp = J @ delta
        r_p = _frozen_residual(op, phi + eps * delta, 0.80, upstream)
        r_m = _frozen_residual(op, phi - eps * delta, 0.80, upstream)
        fd = (r_p - r_m) / (2.0 * eps)
        scale = np.abs(fd[valid]).max()
        rel = np.abs(jvp[valid] - fd[valid]).max() / scale
        assert rel < REL_TOL, f"converged-pocket JVP rel err {rel:.3e}"
