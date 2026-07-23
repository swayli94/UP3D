"""Track V V5 Stage 3 -- the augmented Newton system: full-system FD gate +
coarse smoke Newton (binding: docs/roadmap/track_v.md GV5.1 + the
2026-07-22 pre-registered FD note; design: cases/analysis/v5_tight_coupling/
PRE_REGISTRATION.md; modules under test: pyfp3d/viscous/tight.py Stage-3
operators and pyfp3d/viscous/tight_driver.py).

The pre-registered full state/residual/Jacobian (tight_driver.py module
docstring): x = (phi_free, Gamma, U), F = (F_phi, F_Gamma, F_BL) with the
transpiration delta* = delta*(U) at FROZEN base closure edge inputs (the
pre-registered "no closure calculus" semantics) and veps/veps_s frozen at
the base values (decision 5; the natural-recompute omission is measured,
not gated).

Covers, on the 2.5-D NACA0012 coarse wake-cut strip (M0.5, alpha 2, Re
3.0e6 -- the GV3.1 configuration, coarse recorded):

  - Div vs transpiration_from_delta_star and Drhou vs FD of the
    isentropic rho_e(q^2) u_e(phi) chain -- machine precision / FD sweet
    spot (the Stage-1 operator discipline: verify each factor against its
    reference BEFORE composing);
  - the targeted per-block check (PRE_REGISTRATION FD protocol): the
    J_phiphi augmentation A_free/A_gam vs FD of F_phi w.r.t. (phi_free,
    Gamma) at fixed U, isolated by subtracting the bare-inviscid
    eval_residual;
  - the FULL-SYSTEM FD gate: 4 random max-normalized full-state
    directions v over (phi_free, Gamma, U_flat), central FD ladder
    1e-5/1e-6/1e-7 of augmented_residual vs augmented_jacobian @ v,
    scaled max-norm, tol < 1e-5 (expected 1e-8-1e-10 sweet spot);
    q <= 1e-12 fallback rows masked per the Stage-2 convention
    (tests/test_v5_tight_edge.py:205-219; zero at the k=1 state), masked
    count printed and bounded at 2 %; the veps natural-recompute
    omission measured at eps = 1e-6 and printed (decision 5);
  - the coarse smoke Newton: from the k=1 seed, <= 3 iterations, finite
    iterates + strictly decreasing merit; per-block norms/lam/ds_change
    printed (quadratic tail NOT asserted -- the convergence pass bands
    are the binding medium run's, PRE_REGISTRATION convergence protocol).

Runs in both lanes: default JIT and PYFP3D_NOJIT=1. Under NOJIT the
k=1-fixture tests (targeted augmentation, full-system FD gate, smoke
Newton) are skipped -- the pure-Python k=1 fixture (FP + IBL solves) is
JIT-lane only, the tests/test_v3_coupling.py:204 precedent; the cheap
operator unit tests (Div/Drhou, no k=1 state) still run.
"""

import os

import numpy as np
import pytest
import scipy.sparse as sp

from pyfp3d.viscous import tight_driver as td
from pyfp3d.viscous.tight import (
    edge_velocity_operator,
    rhou_jacobian,
    surface_divergence_vector_operator,
)
from pyfp3d.physics.isentropic import (
    density_derivative_wrt_q_sq_field,
    density_field,
)
from pyfp3d.viscous.transpiration import transpiration_from_delta_star
from tests.v5_state import (
    M_CRIT,
    M_CAP,
    RHO_FLOOR,
    UPWIND_C,
    build_k1_state,
    build_naca_case,
)

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


@pytest.fixture(scope="module")
def pack(k1_state):
    """The Stage-3 pack at the pre-registered k=1 seed."""
    return td.build_tight_pack(k1_state, UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)


def _rel_err(got, ref, floor=1e-30):
    return float(np.max(np.abs(got - ref)) / max(np.max(np.abs(ref)), floor))


# ---------------------------------------------------------------------------
# Div: the frozen-delta* surface-divergence operator
# ---------------------------------------------------------------------------


def test_surface_divergence_vector_operator_matches(naca_case):
    """Div(ds) @ (rho_e u_e).ravel() == transpiration_from_delta_star(
    sm, rho_e, ue, ds) to machine precision (transpiration.py:137-174 --
    the same lumped-divergence idiom as the Stage-1 L test)."""
    _, _, _, case = naca_case
    sm = case.sm
    rng = np.random.default_rng(17)
    rho_e = 0.5 + rng.random(sm.n_node)
    ue = rng.standard_normal((sm.n_node, 3))
    ds = rng.standard_normal(sm.n_node)
    Div = surface_divergence_vector_operator(sm, ds)
    assert Div.shape == (sm.n_node, 3 * sm.n_node)
    for _ in range(3):
        got = Div @ (rho_e[:, None] * ue).ravel()
        ref = transpiration_from_delta_star(sm, rho_e, ue, ds)
        err = _rel_err(got, ref)
        assert err < 1e-13, f"Div mismatch: {err:.3e}"


def test_rhou_jacobian_fd(naca_case):
    """Drhou = (rho I + 2 rho' u u^T) G vs central FD of the isentropic
    rho_e(q^2(phi)) u_e(phi) chain through G, random phi directions,
    ladder 1e-6/1e-7 (no division by q: the map is smooth everywhere, so
    no masking)."""
    mc, _, cfg, case = naca_case
    sm = case.sm
    G = edge_velocity_operator(mc.nodes, case.wall_faces, sm.volume_node_of)
    rng = np.random.default_rng(19)
    phi = rng.standard_normal(len(mc.nodes))

    def rhou(phi_):
        u = (G @ phi_).reshape(sm.n_node, 3)
        r = density_field(
            np.einsum("ij,ij->i", u, u), cfg.m_inf, cfg.gamma_air
        )
        return (r[:, None] * u).ravel()

    ue = (G @ phi).reshape(sm.n_node, 3)
    q2 = np.einsum("ij,ij->i", ue, ue)
    drho = density_derivative_wrt_q_sq_field(q2, cfg.m_inf, cfg.gamma_air)
    Dr = rhou_jacobian(ue, drho, G, cfg.m_inf, cfg.gamma_air)
    assert Dr.shape == (3 * sm.n_node, len(mc.nodes))
    for d in range(3):
        v = rng.standard_normal(len(mc.nodes))
        v /= np.max(np.abs(v))
        an = Dr @ v
        errs = [
            _rel_err((rhou(phi + eps * v) - rhou(phi - eps * v)) / (2 * eps), an)
            for eps in (1e-6, 1e-7)
        ]
        print(f"Drhou FD direction {d}: {errs[0]:.3e} / {errs[1]:.3e}")
        assert min(errs) < 1e-8, f"Drhou FD mismatch: {min(errs):.3e}"


# ---------------------------------------------------------------------------
# shared FD helpers for the k=1 gates
# ---------------------------------------------------------------------------


def _fallback_row_mask(st):
    """The Stage-2 mask (tests/test_v5_tight_edge.py:205-219): drop the 6
    rows of every q <= 1e-12 fallback node and of every node sharing an
    element with one; F_phi/F_Gamma rows are never masked (Div/Drhou
    carry no division by q)."""
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


def _f_phi_transp_only(pack, x):
    """F_phi^tight(x) - F_phi^bare(phi, gamma): the transpiration part
    -(T^T b)[free] in isolation (the bare eval_residual is read with
    external_rhs None -- the driver restores it, so the attribute is
    None here)."""
    phi_free, gamma, _ = pack.split_x(x)
    F_tight = td.block_residuals(pack, x)[0]
    R_bare, _, _ = pack.ws.eval_residual(
        phi_free, gamma, pack.upwind_c, pack.m_crit, pack.m_cap,
        pack.rho_floor,
    )
    return F_tight - R_bare


# ---------------------------------------------------------------------------
# the targeted per-block check: the J_phiphi augmentation
# ---------------------------------------------------------------------------


@K1_JIT_ONLY
def test_j_phi_phi_augmentation_fd(pack):
    """A_free = A_cut T_free and A_gam = A_cut T_gam vs FD of the
    transpiration-only F_phi at fixed U (PRE_REGISTRATION FD protocol,
    targeted check; ladder 1e-5/1e-6/1e-7, tol 1e-5). The Gamma column
    is checked too -- its map T_gam = T[:, dir] V_red + G_jump is the
    one structural surprise of the kutta/Gamma assembly (the slave-jump
    G_jump part is invisible to phi_free-only probes)."""
    bl = td.augmented_jacobian_blocks(pack, pack.x_base())
    rng = np.random.default_rng(41)
    x0 = pack.x_base()
    n_f, n_g = pack.n_free, pack.n_st

    def fd(v, which, eps):
        def at(s):
            pf, g, U = pack.split_x(x0)
            if which == "phi":
                pf = pf + s * eps * v
            else:
                g = g + s * eps * v
            return _f_phi_transp_only(
                pack, np.concatenate([pf, g, U.ravel()])
            )

        return (at(1) - at(-1)) / (2 * eps)

    for which, A, n in (("phi", bl["A_free"], n_f), ("gam", bl["A_gam"], n_g)):
        assert A.shape == (n_f, n)
        worst = 0.0
        for d in range(3):
            v = rng.standard_normal(n)
            v /= np.max(np.abs(v))
            an = A @ v
            errs = [
                _rel_err(fd(v, which, eps), an)
                for eps in (1e-5, 1e-6, 1e-7)
            ]
            print(f"J_phiphi^aug [{which}] direction {d}: "
                  f"{errs[0]:.3e} / {errs[1]:.3e} / {errs[2]:.3e}")
            worst = max(worst, min(errs))
        print(f"J_phiphi^aug [{which}] worst sweet-spot: {worst:.3e}")
        assert worst < 1e-5, f"J_phiphi^aug [{which}] FD mismatch: {worst:.3e}"


# ---------------------------------------------------------------------------
# the full-system FD gate
# ---------------------------------------------------------------------------


@K1_JIT_ONLY
def test_full_system_fd_random_directions(pack, k1_state):
    """GV5.1 Stage-3 gate: the assembled augmented Jacobian vs central FD
    of augmented_residual over FULL-STATE random directions (phi_free,
    Gamma, U_flat), max-normalized, ladder 1e-5/1e-6/1e-7, scaled
    max-norm, tol < 1e-5 (the PRE_REGISTRATION FD protocol). The FD
    reference runs with veps/veps_s frozen (decision 5); the
    natural-recompute omission is measured at eps = 1e-6 and printed."""
    x0 = pack.x_base()
    J = td.augmented_jacobian(pack, x0)
    n = len(x0)
    assert J.shape == (n, n)
    row_mask, bad = _fallback_row_mask(k1_state)
    n_bl = 6 * pack.n_s
    # full-vector keep mask: F_phi/F_Gamma rows always kept
    keep = np.ones(n, dtype=bool)
    keep[n - n_bl:] = row_mask
    print(f"fallback nodes (q <= 1e-12): {len(bad)}; "
          f"masked rows: {int((~keep).sum())} / {n}")
    assert (~keep).sum() <= 0.02 * n_bl

    rng = np.random.default_rng(97)
    worst_sweet = 0.0
    for d in range(4):
        v = rng.standard_normal(n)
        v /= np.max(np.abs(v))
        an = (J @ v)[keep]
        scale = max(float(np.max(np.abs(an))), 1e-12)
        errs = []
        for eps in (1e-5, 1e-6, 1e-7):
            fd = (td.augmented_residual(pack, x0 + eps * v)
                  - td.augmented_residual(pack, x0 - eps * v)) / (2 * eps)
            errs.append(float(np.max(np.abs(fd[keep] - an))) / scale)
        # the frozen-veps omission (decision 5), measured at eps = 1e-6
        fd_unf = (td.augmented_residual(pack, x0 + 1e-6 * v, veps_frozen=False)
                  - td.augmented_residual(pack, x0 - 1e-6 * v, veps_frozen=False)) / 2.0e-6
        fd_fr = (td.augmented_residual(pack, x0 + 1e-6 * v)
                 - td.augmented_residual(pack, x0 - 1e-6 * v)) / 2.0e-6
        veps_eff = float(np.max(np.abs(fd_unf[keep] - fd_fr[keep])))
        print(f"full-system FD direction {d}: eps=1e-5/1e-6/1e-7 errs = "
              f"{errs[0]:.3e} / {errs[1]:.3e} / {errs[2]:.3e}   "
              f"veps-freeze contribution {veps_eff:.3e} (scaled "
              f"{veps_eff / scale:.3e})")
        worst_sweet = max(worst_sweet, min(errs))
    print(f"full-system FD worst sweet-spot: {worst_sweet:.3e}")
    assert worst_sweet < 1e-5, f"full-system FD mismatch: {worst_sweet:.3e}"


# ---------------------------------------------------------------------------
# the coarse smoke Newton
# ---------------------------------------------------------------------------


@K1_JIT_ONLY
def test_newton_tight_smoke(pack):
    """Coarse smoke: <= 3 augmented Newton iterations from the k=1 seed;
    finite iterates and a strictly decreasing merit (the quadratic tail
    and the N_aug <= 2 pass band are the binding medium run's, not this
    smoke's). Per-block norms, step lengths and delta*-changes printed
    for the record."""
    res = td.newton_tight(pack, max_iter=3, verbose=True)
    assert np.all(np.isfinite(res["x"]))
    hist = res["history"]
    merits = [h["merit"] for h in hist]
    for i in range(1, len(merits)):
        assert merits[i] < merits[i - 1], (
            f"merit not decreasing at iteration {i}: "
            f"{merits[i - 1]:.3e} -> {merits[i]:.3e}"
        )
    for h in hist:
        bm = h["block_max"]
        print(f"smoke iter {h['iter']}: |F_phi|={bm[0]:.3e} "
              f"|F_Gam|={bm[1]:.3e} |F_BL|={bm[2]:.3e} "
              f"merit={h['merit']:.3e} lam={h['lam']} "
              f"ds_change={h['ds_change']:.3e}")
    print(f"smoke: converged={res['converged']} n_iter={res['n_iter']} "
          f"ds_change_last={res['ds_change_last']:.3e}")
