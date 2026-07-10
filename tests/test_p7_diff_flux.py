"""P7 (gate G7.3): frozen-selection differentiability of the shipped walk
flux — finite-difference verification of ∂ρ̃/∂φ (design.md §3.1/§6.3,
López B.3–B.6).

The Newton prerequisite is the sensitivity of the P4 walk artificial
density at FROZEN upstream selection u(e):

    ∂ρ̃_e/∂φ_k = s_e·(2 ∇φ_e·∇N_k|_e) + s_u·(2 ∇φ_u·∇N_k|_u)

with (s_e, s_u) from `rho_tilde_sensitivities_sweep`. Verified here as a
directional (JVP) central difference against the SHIPPED `rho_tilde_sweep`
with the selection held frozen — the forward flux is reused verbatim, so
what is tested is exactly the derivative P8's Term-2/Term-3 assembly needs.

Unit-level checks on a structured cube — no transonic solve (the walk's
transonic behaviour is locked by G4.1/G4.3; its subcritical bit-no-op by
G4.2). Runs under PYFP3D_NOJIT=1 as well (pure-Python kernels).
"""

import numpy as np

from pyfp3d.kernels.jacobian import PicardOperator
from pyfp3d.kernels.upwind import (
    UpwindOperator,
    rho_tilde_sweep,
)
from pyfp3d.physics.isentropic import (
    GAMMA,
    density_field,
    mach_number_squared,
    mach_squared_derivative_wrt_q_sq,
)

from .mesh_utils import generate_structured_cube_mesh

M_CRIT = 0.95
UPWIND_C = 1.5
RHO_FLOOR = 0.05
REL_TOL = 1e-6          # G7.3 acceptance


def _cube():
    return generate_structured_cube_mesh(n=6, L=1.0)


def _degeneracy_breaker(n_nodes, amp=0.005, seed=42):
    """Small deterministic nodal perturbation. On the structured cube any
    lattice-symmetric phi (separable, or with cell-symmetric cross terms)
    gives whole families of element pairs with IDENTICAL gradients, parking
    e and its frozen u(e) exactly on the max(nu_e, nu_u) tie kink where a
    central difference is a branch average, not a derivative (measured:
    rel err ~1e-5 from ties alone). Generic nodal noise removes every such
    tie; the JVP identity being verified holds for arbitrary phi."""
    rng = np.random.default_rng(seed)
    return amp * rng.standard_normal(n_nodes)


def _frozen_rho_tilde(q2, m_inf, upstream, rho_floor=RHO_FLOOR):
    """The SHIPPED flux at a frozen selection: density law + rho_tilde_sweep
    reused verbatim (byte-identical physics), only u(e) held fixed."""
    rho = density_field(q2, m_inf)
    nu = np.empty_like(q2)
    rt = np.empty_like(q2)
    rho_tilde_sweep(q2, rho, upstream, m_inf, M_CRIT, UPWIND_C, GAMMA,
                    rho_floor, nu, rt)
    return rt


def _jvp_check(phi, m_inf, seed=0, n_dirs=3, eps=5e-7, rho_floor=RHO_FLOOR,
               max_kink_frac=0.02):
    """Directional-derivative (JVP) check of (s_e, s_u) against a central
    difference of the shipped flux at frozen selection. Returns
    (max_rel_err, s_e, s_u, upstream, q2) — rel err is field-scale
    normalized: ||J·δ − FD||_inf / ||FD||_inf over FD-valid elements.

    Kink guard: ρ̃ is C⁰ but not C¹ exactly AT the max(ν_e, ν_u) tie and at
    the switch threshold M² = M_c² (the measure-zero locus of design.md
    §3.1 — López's Newton converges because the active set freezes there).
    A central difference straddling such a kink returns a branch AVERAGE,
    not a derivative, so elements within an ε-neighbourhood of a kink are
    excluded from the comparison; their fraction is asserted tiny (and is
    0 for the generic fields used below — only symmetry-degenerate fields
    populate the locus)."""
    nodes, elements = _cube()
    op = PicardOperator(nodes, elements)
    upw = UpwindOperator(nodes, elements, weighted=False)

    grad, q2 = op.velocities(phi)
    grad, q2 = grad.copy(), q2.copy()
    rho = density_field(q2, m_inf)
    s_e, s_u, upstream = upw.rho_tilde_sensitivities(
        grad, q2, rho, m_inf, UPWIND_C, M_CRIT, rho_floor=rho_floor)
    s_e, s_u, upstream = s_e.copy(), s_u.copy(), upstream.copy()

    # kink-neighbourhood mask (branch could flip inside the FD stencil):
    # the flip zone is |Δν| ≲ eps·|d(Δν)/dt| ~ eps·ν'·(|dq2_e|+|dq2_u|)
    # ≈ 4e-6 at eps=5e-7 on these fields; 3e-5 keeps ≥ 5× margin while
    # excluding only a handful of elements on generic (non-degenerate) fields.
    kink_guard = 3e-5
    mc2 = M_CRIT ** 2
    m2 = mach_number_squared(q2, m_inf, GAMMA)
    m2u = m2[upstream]
    nu_e = UPWIND_C * np.maximum(0.0, 1.0 - mc2 / np.maximum(m2, mc2))
    nu_u = UPWIND_C * np.maximum(0.0, 1.0 - mc2 / np.maximum(m2u, mc2))
    active = np.maximum(nu_e, nu_u) > 0.0
    self_up = upstream == np.arange(len(q2))
    # self-upstream elements have a zero upwind jump: rho_tilde == rho_e for
    # ANY nu, so the max tie is irrelevant there (no kink) — not excluded.
    near_tie = active & ~self_up & (np.abs(nu_e - nu_u) < kink_guard)
    near_thresh = ~self_up & ((np.abs(m2 - mc2) < kink_guard)
                              | (np.abs(m2u - mc2) < kink_guard))
    valid = ~(near_tie | near_thresh)
    n_excl = int((~valid).sum())
    assert n_excl <= max_kink_frac * len(valid), (
        f"{n_excl} elements near a max-kink — field too symmetry-degenerate "
        "for an FD check (see docstring)")

    rng = np.random.default_rng(seed)
    max_rel = 0.0
    for _ in range(n_dirs):
        delta = rng.standard_normal(len(phi))
        delta /= np.abs(delta).max()
        gradd, _ = op.velocities(delta)
        gradd = gradd.copy()

        # analytic JVP: dq2_e = 2 grad_e . gradd_e (chain through u frozen)
        dq2 = 2.0 * np.einsum("ij,ij->i", grad, gradd)
        jvp = s_e * dq2 + s_u * dq2[upstream]

        _, q2p = op.velocities(phi + eps * delta)
        rt_p = _frozen_rho_tilde(q2p.copy(), m_inf, upstream, rho_floor)
        _, q2m = op.velocities(phi - eps * delta)
        rt_m = _frozen_rho_tilde(q2m.copy(), m_inf, upstream, rho_floor)
        fd = (rt_p - rt_m) / (2.0 * eps)

        scale = np.abs(fd[valid]).max()
        assert scale > 0.0
        rel = np.abs(jvp[valid] - fd[valid]).max() / scale
        max_rel = max(max_rel, rel)
    return max_rel, s_e, s_u, upstream, q2


def _regimes(q2, upstream, m_inf):
    """Classify elements by the frozen-selection branch actually taken.
    Also returns the active-switch TIE set (nu_e == nu_u > 0): the genuine
    C⁰ kink of max(nu_e, nu_u), where a central difference straddles two
    branches and is not a derivative — the FD fields below are designed to
    keep this set empty (the measure-zero locus of design.md §3.1)."""
    m2 = mach_number_squared(q2, m_inf, GAMMA)
    m2u = m2[upstream]
    mc2 = M_CRIT ** 2
    nu_e = UPWIND_C * np.maximum(0.0, 1.0 - mc2 / np.maximum(m2, mc2))
    nu_u = UPWIND_C * np.maximum(0.0, 1.0 - mc2 / np.maximum(m2u, mc2))
    self_up = upstream == np.arange(len(q2))
    subsonic = (~self_up) & (np.maximum(nu_e, nu_u) == 0.0)
    accel = (~self_up) & (nu_e >= nu_u) & (nu_e > 0.0)
    shockpt = (~self_up) & (nu_u > nu_e)
    ties = (~self_up) & (nu_e > 0.0) & (nu_e == nu_u)
    return subsonic, accel, shockpt, ties


# ---------------------------------------------------------------- physics


def test_mach_squared_derivative_matches_central_difference():
    """New physics scalar dM²/dq² vs central diff of mach_number_squared."""
    eps = 1e-7
    for m_inf in (0.3, 0.5, 0.8, 0.84):
        for q2 in (0.2, 0.7, 1.0, 1.5, 2.2):
            fd = (mach_number_squared(q2 + eps, m_inf, GAMMA)
                  - mach_number_squared(q2 - eps, m_inf, GAMMA)) / (2 * eps)
            an = mach_squared_derivative_wrt_q_sq(q2, m_inf, GAMMA)
            assert an > 0.0                      # M² monotone in q²
            assert abs(an - fd) / abs(fd) < 1e-8


# ---------------------------------------------------------- regime cases


def test_subsonic_field_sensitivities_are_pure_density_derivative():
    """(a) ν ≡ 0: s_u ≡ 0 and s_e ≡ ρ'(q²) everywhere; FD confirms."""
    nodes, _ = _cube()
    phi = nodes[:, 0].copy()                     # q² = 1, subcritical at 0.5
    rel, s_e, s_u, upstream, q2 = _jvp_check(phi, m_inf=0.5)
    assert rel < REL_TOL
    assert np.all(s_u == 0.0)
    # same closed form as the scalar density_derivative_wrt_q_sq:
    # rho' = -(M_inf^2/2) rho^(2-gamma)
    expected = -0.5 * 0.5 ** 2 * density_field(q2, 0.5) ** (2.0 - GAMMA)
    assert np.allclose(s_e, expected, rtol=1e-14, atol=0)


def test_accelerating_supersonic_field():
    """(b) q² grows along the flow → ν = ν_e branch (B.3 + B.4).

    Note the y/z terms: a purely 1D field on the structured cube gives
    whole slabs of elements with IDENTICAL gradients, so e and its frozen
    u(e) tie exactly at the max(ν_e, ν_u) kink and the central difference
    straddles two branches (measured: rel err ~1e-5, branch-average, not a
    derivative bug). Genuinely 3D fields keep the tie set empty."""
    nodes, _ = _cube()
    x = nodes[:, 0]
    phi = x + 0.3 * x ** 2 + _degeneracy_breaker(len(nodes))
    rel, s_e, s_u, upstream, q2 = _jvp_check(phi, m_inf=0.8)
    assert rel < REL_TOL
    subsonic, accel, shockpt, ties = _regimes(q2, upstream, 0.8)
    assert ties.sum() == 0                        # FD branch-unambiguous
    assert accel.sum() > 0                        # branch exercised
    assert np.any(s_u[accel] != 0.0)              # upstream coupling nonzero


def test_decelerating_field_populates_shock_point_branch():
    """(c) q² decays along the flow → ν = ν_u shock-point branch."""
    nodes, _ = _cube()
    x = nodes[:, 0]
    phi = 1.6 * x - 0.3 * x ** 2 + _degeneracy_breaker(len(nodes))
    rel, s_e, s_u, upstream, q2 = _jvp_check(phi, m_inf=0.8)
    assert rel < REL_TOL
    subsonic, accel, shockpt, ties = _regimes(q2, upstream, 0.8)
    assert ties.sum() == 0
    assert shockpt.sum() > 0                      # branch exercised


def test_mixed_field_all_regimes_nonempty():
    """(d) oscillating gradient exercises all three regimes in one field."""
    nodes, _ = _cube()
    x = nodes[:, 0]
    phi = (1.25 * x + 0.12 * np.sin(2.0 * np.pi * x)
           + _degeneracy_breaker(len(nodes)))
    rel, s_e, s_u, upstream, q2 = _jvp_check(phi, m_inf=0.8, n_dirs=5)
    assert rel < REL_TOL
    subsonic, accel, shockpt, ties = _regimes(q2, upstream, 0.8)
    assert ties.sum() == 0
    assert subsonic.sum() > 0
    assert accel.sum() > 0
    assert shockpt.sum() > 0


def test_floored_elements_have_zero_sensitivities():
    """(e) the ρ̃ floor is a flat clamp: floored elements must report
    s_e = s_u = 0, and the FD (which reuses the floored flux) agrees."""
    nodes, _ = _cube()
    x = nodes[:, 0]
    phi = x + 0.3 * x ** 2
    # rho_floor above the whole density range (ρ ≤ ρ_stag ≈ 1.35 at M0.8)
    # floors every element -> all sensitivities exactly zero.
    nodes_, elements = _cube()
    op = PicardOperator(nodes_, elements)
    upw = UpwindOperator(nodes_, elements, weighted=False)
    grad, q2 = op.velocities(phi)
    rho = density_field(q2, 0.8)
    s_e, s_u, upstream = upw.rho_tilde_sensitivities(
        grad, q2, rho, 0.8, UPWIND_C, M_CRIT, rho_floor=2.0)
    assert np.all(s_e == 0.0)
    assert np.all(s_u == 0.0)
    # and the frozen-selection FD of the floored flux is identically zero
    rt = _frozen_rho_tilde(q2.copy(), 0.8, upstream.copy(), rho_floor=2.0)
    assert np.all(rt == 2.0)


# ------------------------------------------------------------- guardrails


def test_forward_walk_flux_unchanged_by_sensitivity_call():
    """Calling rho_tilde_sensitivities must not perturb the forward path
    (separate buffers; same walk selection)."""
    nodes, elements = _cube()
    op = PicardOperator(nodes, elements)
    upw = UpwindOperator(nodes, elements, weighted=False)
    phi = 1.25 * nodes[:, 0].copy()
    grad, q2 = op.velocities(phi)
    rho = density_field(q2, 0.82)
    ref = upw.rho_tilde(grad, q2, rho, 0.82, UPWIND_C, M_CRIT).copy()
    upw.rho_tilde_sensitivities(grad, q2, rho, 0.82, UPWIND_C, M_CRIT)
    again = upw.rho_tilde(grad, q2, rho, 0.82, UPWIND_C, M_CRIT).copy()
    assert np.array_equal(ref, again)


def test_kernel_mode_sensitivities_not_implemented():
    """Walk mode only: the dense kernel-flux Jacobian is a P8 measurement
    item (design.md §3.2), not a P7 deliverable."""
    import pytest

    nodes, elements = _cube()
    upw = UpwindOperator(nodes, elements, weighted=True, mode="kernel")
    op = PicardOperator(nodes, elements)
    phi = nodes[:, 0].copy()
    grad, q2 = op.velocities(phi)
    rho = density_field(q2, 0.5)
    with pytest.raises(NotImplementedError):
        upw.rho_tilde_sensitivities(grad, q2, rho, 0.5, UPWIND_C, M_CRIT)
