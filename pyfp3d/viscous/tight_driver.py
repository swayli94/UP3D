"""Track V V5 Stage 3 -- the augmented (tight) Newton driver (binding:
docs/roadmap/track_v.md GV5.1 + the 2026-07-22 pre-registered FD note;
design: cases/analysis/v5_tight_coupling/PRE_REGISTRATION.md).

The pre-registered full state, residual and Jacobian (PRE_REGISTRATION
"System", "Jacobian blocks"):

    x = (phi_free, Gamma, U)          n_free + n_st + 6 n_s unknowns
    F = (F_phi, F_Gamma, F_BL)

    F_phi   = R_free(phi, Gamma) - [T^T b(m(phi, U))]_free   (inviscid
              reduced residual with the V2 transpiration channel LIVE,
              solve/newton.py:372-379; b = assemble_transpiration_rhs,
              m = P . div_Gamma(rho_e(phi) u_e(phi) delta*(U)) -- the
              pre-registered m: the FLUX rho_e u_e is a function of phi
              (the augmentation chain), while delta* = delta*(U) carries
              NO phi-dependence: the closure is evaluated with its edge
              inputs (q, rho_e, mach_e) frozen at the pack base state,
              so the augmentation -T^T W d m/d phi is the explicit
              rho_e(q^2(phi)) u_e(phi) algebra with NO closure calculus
              (PRE_REGISTRATION lines 20-34), and the dout_e/re_d/e_prime
              edge chain of delta* is absent from BOTH F and J by
              construction -- consistent and FD-clean)
    F_Gamma = the Kutta residual (probe estimator: taper * kutta_targets
              - gamma, solve/newton.py:418-421)
    F_BL    = IBL3Solver.residual(U) at the edge data recovered from phi
              (coupling.py:710-726), veps/veps_s FROZEN at the pack base
              values within a Newton step (decision 5)

assembled as the UN-eliminated block system (solve_newton_lifting
eliminates the Kutta block row via kutta_blocks -- delta_gamma =
K_tilde dphi + F_tilde, solve/newton.py:459-501; here the row is kept
explicit, with D = dF_Gamma/dGamma = -I exactly for the probe
estimator):

    [ J_ff + A_free   B + A_gam    J_phi,BL ]
    [     K              -I           0     ]
    [ J_BL,phi T_free  J_BL,phi T_gam  J_BL,BL ]

with

    J_ff, B   = ws.assemble_coupled (solve/newton.py:430-457; B already
                carries the far-field vortex columns J_red[:, dir] V_red
                and the wake-jump block H_J)
    A_cut     = assemble_j_phi_phi_aug = -(T^T W S P Div Drhou)[free, :],
                the phi-derivative of the transpiration wall RHS at
                frozen (U, delta*) (tight.py, Stage-3 operators)
    J_phi,BL  = assemble_j_phi_bl (Stage 1), J_BL,phi = J_e D_ue G
                (Stage 2), J_BL,BL = solver.residual_jacobian(U)
    K         = ws.K (solve/newton.py:256-280, the probe row dF/dphi_free)

Column maps (the SAME d phi_cut/d unknown for every phi_cut-derivative
block, mirroring B's construction at solve/newton.py:453-456 with
phi_cut = T phi_red + G_jump Gamma, constraints/wake.py:50-53/104-112):

    T_free = d phi_cut/d phi_free = T[:, ws.free]
    T_gam  = d phi_cut/d Gamma    = T[:, ws.dir_red] @ V_red + G_jump

where V_red is the affine far-field basis (set_mach,
solve/newton.py:295-324) and G_jump the per-station slave-indicator
matrix of reduce_operator (constraints/wake.py:76-80). B's H_J term is
exactly (T^T J) G_jump, so the kutta/far-field phi-dependence seen by
the inviscid blocks and by the coupling blocks is one and the same map.

Frozen/lagged pieces (the pre-registered decisions): the inflow
Dirichlet band and its seed states are frozen at the pack build state
(the coupling.py:728-763 k=1 freeze); the loose-loop under-relaxation
and the max(ds, 0) floor are dropped (decision 2); the zone partition of
G is state-independent; veps/veps_s are frozen within a Newton step and
recomputed between steps by the caller (decision 5 -- the omission is
measured in the full-system FD gate by running augmented_residual with
veps_frozen=False). Only the PROBE kutta estimator is wired (the
pressure estimator's un-eliminated rows are [Kp_free, D] of
kutta_blocks, not [K, -I]).

Convergence criterion: per-block max-norm |F_j|_inf <= tol_abs +
tol * |F_j(x0)|_inf (an absolute-plus-relative test; a purely relative
test is unmeaningful for blocks the seed already converged -- F_Gamma
of an inviscid-converged seed is O(1e-11)). The merit for the line
search is the PRE-REGISTERED unnormalized sum
merit = |F_phi|_2^2 + |F_Gamma|_2^2 + |F_BL|_2^2 (PRE_REGISTRATION
convergence protocol; the P8/P14 shared-merit idiom with F_Gamma in
Gamma-like units, solve/newton.py:104-116/750/1068) -- a
seed-normalized merit sum_j (|F_j|/|F_j(x0)|)^2 was tried first and is
numerically BROKEN at this seed: F_Gamma(x0) ~ 5.6e-17 (inviscid
converged), so its normalized term explodes by ~1e13 on the first
Gamma-moving step and the line search vetoes every full step (measured
2026-07-23 on the coarse k=1 seed: lam collapsed to 0.12/0.06 with the
block max-norms stuck). With the P8/P14 merit the blocks contribute at
their natural scales. The backtracking is the P8/P14 safety-only idiom
(solve/newton.py:1053-1087): accept the full step unless the merit
grows, else halve and take the least-bad finite try. The halving budget
is 30 (not the FP driver's 5, solve/newton.py:1062): at the stalled k=1
IBL floor the steady BL block is massively ill-conditioned (measured
2026-07-23: cond(J_BL,BL) ~ 4e10, 501/1236 singular values below 1e-6 of
max -- the pre-registered basin risk), so the raw Newton Delta_U is
O(5e2-6e3) and needs ~12 halvings to re-enter the finite basin; the
accept/least-bad/raise logic itself is unchanged. The delta*-change of
the last accepted step is recorded for the tol_ds = 1e-3 cross-check
(PRE_REGISTRATION convergence protocol).
"""

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

from pyfp3d.physics.isentropic import (
    density_derivative_wrt_q_sq_field,
    density_field,
    mach_squared_field,
)
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.tight import (
    assemble_j_bl_phi,
    assemble_j_phi_bl,
    assemble_j_phi_phi_aug,
    closure_ds1_jacobian,
    edge_data_jacobian,
    edge_velocity_operator,
    pin_mask_matrix,
    rhou_jacobian,
    surface_divergence_delta_operator,
    surface_divergence_vector_operator,
    surface_scatter_matrix,
    wall_load_matrix,
)
from pyfp3d.viscous.transpiration import (
    assemble_transpiration_rhs,
    transpiration_from_delta_star,
)

__all__ = [
    "TightPack",
    "build_tight_pack",
    "augmented_residual",
    "block_residuals",
    "augmented_jacobian_blocks",
    "augmented_jacobian",
    "newton_tight",
]


@dataclass
class TightPack:
    """Everything the augmented residual/Jacobian needs, fixed at the
    pack build state (the pre-registered k=1 seed: inviscid-converged
    (phi, Gamma) + ONE IBL solve, tests/v5_state.py::build_k1_state).

    The pack owns NO solver state: augmented_residual/jacobian mutate
    solver edge data and ws.external_rhs per evaluation and restore them
    before returning, so both are pure functions of x.
    """

    ws: object  # solve/newton.py NewtonWorkspace (probe estimator)
    solver: object  # viscous/ibl3.py IBL3Solver at the base edge data
    case: object  # viscous/coupling.py airfoil case wiring
    cfg: object  # viscous/coupling.py CouplingConfig
    sm: object  # viscous/surface_mesh.py SurfaceMesh
    nodes: np.ndarray  # (n_cut, 3) cut-mesh nodes
    wall_faces: np.ndarray  # (F, 3) wall connectivity (cut-mesh ids)
    # fixed sparse operators (tight.py)
    W: sp.csr_matrix
    S: sp.csr_matrix
    P: sp.csr_matrix
    G: sp.csr_matrix
    T_free: sp.csr_matrix  # d phi_cut/d phi_free = T[:, free]
    T_gam: sp.csr_matrix  # d phi_cut/d Gamma = T[:, dir] V_red + G_jump
    # base state (the k=1 seed)
    phi_free0: np.ndarray
    gamma0: np.ndarray
    U0: np.ndarray  # (n_s, 6)
    ue0: np.ndarray
    q0: np.ndarray
    rho_e0: np.ndarray
    mach_e0: np.ndarray
    mu: float  # 1/re_chord (coupling.py:681)
    veps0: Tuple[float, float]  # (solver._veps, solver._veps_s) at base
    # dims and the GV3.1 Newton-leg constants (tests/v5_state.py:32-33)
    n_cut: int
    n_s: int
    n_free: int
    n_st: int
    upwind_c: float
    m_crit: float
    m_cap: float
    rho_floor: float

    # -- state (un)packing ----------------------------------------------------
    def x_base(self) -> np.ndarray:
        """The pre-registered seed x0 = (phi_free0, gamma0, U0.ravel())."""
        return np.concatenate(
            [self.phi_free0, self.gamma0, self.U0.ravel()]
        )

    def split_x(self, x: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """x -> (phi_free, gamma, U) with U the (n_s, 6) IBL state."""
        x = np.asarray(x, dtype=np.float64)
        n_f, n_g = self.n_free, self.n_st
        if x.shape != (n_f + n_g + 6 * self.n_s,):
            raise ValueError(
                f"x must be ({n_f + n_g + 6 * self.n_s},), got {x.shape}"
            )
        return (
            x[:n_f],
            x[n_f: n_f + n_g],
            x[n_f + n_g:].reshape(self.n_s, 6),
        )

    def split_F(self, F: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """F -> (F_phi, F_Gamma, F_BL flat) block views."""
        n_f, n_g = self.n_free, self.n_st
        return F[:n_f], F[n_f: n_f + n_g], F[n_f + n_g:]


def build_tight_pack(
    state: Dict,
    upwind_c: float = 1.5,
    m_crit: float = 0.95,
    m_cap: float = 3.0,
    rho_floor: float = 0.05,
) -> TightPack:
    """Assemble the pack from a k=1 state mapping (the tests/v5_state.py
    ::build_k1_state keys: mc, cfg, case, sm, ws, solver, phi_free,
    gamma, U, ue_surf, q, rho_e, mach_e, mu, n_cut, le_mask_vol).

    The Newton-leg constants default to the GV3.1 settings
    (cases/analysis/v3_loose_coupling/run.py:67). Restrictions: the PROBE
    kutta estimator only (the un-eliminated row [K, -I] is hard-coded),
    set_mach already called, ws.external_rhs None (the driver owns the
    attribute during evaluations).
    """
    mc, cfg, case, sm = state["mc"], state["cfg"], state["case"], state["sm"]
    ws, solver = state["ws"], state["solver"]
    if ws.kutta_estimator != "probe":
        raise ValueError(
            "tight_driver wires the un-eliminated PROBE kutta row "
            "[K, -I] only (kutta_blocks, solve/newton.py:486-487); got "
            f"{ws.kutta_estimator!r}"
        )
    if ws.m_inf is None:
        raise ValueError("ws.set_mach must be called before packing")
    if ws.external_rhs is not None:
        raise ValueError("ws.external_rhs must be None at pack time")
    n_cut, n_s = state["n_cut"], sm.n_node
    gamma = np.asarray(state["gamma"], dtype=np.float64)
    if gamma.shape != (ws.n_st,):
        raise ValueError(f"gamma must be ({ws.n_st},), got {gamma.shape}")

    W = wall_load_matrix(mc.nodes, case.wall_faces)
    S = surface_scatter_matrix(sm, n_cut)
    P = pin_mask_matrix(n_s, case.outflow_pin_surf)
    G = edge_velocity_operator(
        mc.nodes,
        case.wall_faces,
        sm.volume_node_of,
        elements=case.elements,
        le_band_mask=state["le_mask_vol"],
        n_smooth_passes=cfg.n_smooth_passes,
    )
    # the phi_cut column maps (module docstring): T_free = T[:, free],
    # T_gam = T[:, dir_red] @ V_red + G_jump with G_jump the
    # slave-indicator matrix of constraints/wake.py:76-80.
    T = ws.con.T
    T_free = T[:, ws.free].tocsr()
    wc = ws.wc
    G_jump = sp.coo_matrix(
        (
            np.ones(len(wc.slave_nodes), dtype=np.float64),
            (wc.slave_nodes, wc.node_station),
        ),
        shape=(n_cut, ws.n_st),
    ).tocsr()
    T_gam = (
        T[:, ws.dir_red] @ sp.csr_matrix(ws.V_red) + G_jump
    ).tocsr()
    return TightPack(
        ws=ws,
        solver=solver,
        case=case,
        cfg=cfg,
        sm=sm,
        nodes=np.asarray(mc.nodes, dtype=np.float64),
        wall_faces=case.wall_faces,
        W=W,
        S=S,
        P=P,
        G=G,
        T_free=T_free,
        T_gam=T_gam,
        phi_free0=np.asarray(state["phi_free"], dtype=np.float64).copy(),
        gamma0=gamma.copy(),
        U0=np.ascontiguousarray(state["U"], dtype=np.float64).copy(),
        ue0=np.asarray(state["ue_surf"], dtype=np.float64).copy(),
        q0=np.asarray(state["q"], dtype=np.float64).copy(),
        rho_e0=np.asarray(state["rho_e"], dtype=np.float64).copy(),
        mach_e0=np.asarray(state["mach_e"], dtype=np.float64).copy(),
        mu=float(state["mu"]),
        veps0=(float(solver._veps), float(solver._veps_s)),
        n_cut=n_cut,
        n_s=n_s,
        n_free=ws.n_free,
        n_st=ws.n_st,
        upwind_c=float(upwind_c),
        m_crit=float(m_crit),
        m_cap=float(m_cap),
        rho_floor=float(rho_floor),
    )


# ---------------------------------------------------------------------------
# internal chains (mirrors of the coupling.py:710-835 loose-loop pieces)
# ---------------------------------------------------------------------------


def _phi_cut_of(pack: TightPack, phi_free, gamma) -> np.ndarray:
    """phi_cut = con.expand(phi_red, gamma) with the reduced vector
    assembled from the unknowns (solve/newton.py:339-343)."""
    ws = pack.ws
    phi_red = np.empty(ws.n_red, dtype=np.float64)
    phi_red[ws.free] = phi_free
    phi_red[ws.dir_red] = ws.vals0_red + ws.V_red @ gamma
    return ws.con.expand(phi_red, gamma)


def _edge_data_of(pack: TightPack, phi_cut):
    """phi_cut -> (u_e, q, rho_e, mach_e) through G + the isentropic
    packet (coupling.py:723-726)."""
    ue = (pack.G @ phi_cut).reshape(pack.n_s, 3)
    q2 = np.einsum("ij,ij->i", ue, ue)
    rho = density_field(q2, pack.cfg.m_inf, pack.cfg.gamma_air)
    mach = np.sqrt(mach_squared_field(q2, pack.cfg.m_inf, pack.cfg.gamma_air))
    return ue, np.sqrt(q2), rho, mach


def _closure_packet(pack: TightPack, U):
    """outs/douts of closures.closure_all at U with the closure's edge
    inputs (q, rho_e, mach_e) FROZEN at the pack base state -- the
    pre-registered delta*(U) semantics (module docstring): the
    transpiration delta* carries no phi-dependence, so the phi-side
    augmentation is the explicit rho_e u_e algebra with no closure
    calculus. (The U-derivative douts is unaffected by the freeze: it is
    d closure/d U at fixed edge inputs, which is what J_phi,BL's D
    consumes.)"""
    n_s = pack.n_s
    outs = np.empty((n_s, C.N_OUT), dtype=np.float64)
    douts = np.empty((n_s, C.N_OUT, 6), dtype=np.float64)
    douts_e = np.empty((n_s, C.N_OUT, 2), dtype=np.float64)
    C.closure_all(
        np.ascontiguousarray(U, dtype=np.float64),
        pack.q0,
        pack.rho_e0,
        np.full(n_s, pack.mu, dtype=np.float64),
        pack.mach_e0,
        pack.case.turbulent_flags,
        C.C_L_DEFAULT,
        outs,
        douts,
        douts_e,
    )
    return outs, douts


def _transpiration_rhs(pack: TightPack, rho, ue, ds) -> np.ndarray:
    """m_surf -> pin -> scatter -> wall RHS (the production reference
    functions, transpiration.py:66-114/165-174; the operators W/S/P are
    their machine-precision matrix forms, Stage-1 unit tests)."""
    m_surf = transpiration_from_delta_star(pack.sm, rho, ue, ds)
    return assemble_transpiration_rhs(
        pack.nodes, pack.wall_faces, pack.S @ (pack.P @ m_surf)
    )


def _restore_solver(pack: TightPack) -> None:
    """Base edge data + frozen veps back onto the solver."""
    pack.solver.update_edge_data(
        pack.ue0,
        pack.rho_e0,
        np.full(pack.n_s, pack.mu, dtype=np.float64),
        pack.mach_e0,
    )
    pack.solver._veps, pack.solver._veps_s = pack.veps0


def block_residuals(
    pack: TightPack, x: np.ndarray, veps_frozen: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(F_phi, F_Gamma, F_BL) at x (module docstring). Pure in x: the
    solver edge data and ws.external_rhs are restored before returning.

    veps_frozen=True (decision 5, the pre-registered system): after each
    update_edge_data the solver's _veps/_veps_s are reset to the pack
    base values, so the F_BL evaluations match the frozen-veps analytic
    blocks (J_e / residual_jacobian do not differentiate veps). False
    lets veps/veps_s recompute naturally -- the FD gate measures the
    omission with it. The transpiration delta* is closure(U) at the
    FROZEN base edge inputs (the pre-registered delta*(U); the m-flux
    rho_e u_e stays current).
    """
    ws, solver = pack.ws, pack.solver
    phi_free, gamma, U = pack.split_x(x)
    phi_cut = _phi_cut_of(pack, phi_free, gamma)
    ue, q, rho, mach = _edge_data_of(pack, phi_cut)
    outs, _ = _closure_packet(pack, U)
    rhs = _transpiration_rhs(pack, rho, ue, outs[:, C.OUT_DS1])
    mu_arr = np.full(pack.n_s, pack.mu, dtype=np.float64)

    assert ws.external_rhs is None, "no nested augmented evaluations"
    try:
        ws.external_rhs = rhs
        R_free, F_gam, _ = ws.eval_residual(
            phi_free, gamma, pack.upwind_c, pack.m_crit, pack.m_cap,
            pack.rho_floor,
        )
    finally:
        ws.external_rhs = None
    try:
        solver.update_edge_data(ue, rho, mu_arr, mach)
        if veps_frozen:
            solver._veps, solver._veps_s = pack.veps0
        F_bl = solver.residual(U).ravel()
    finally:
        _restore_solver(pack)
    return R_free, F_gam, F_bl


def augmented_residual(
    pack: TightPack, x: np.ndarray, veps_frozen: bool = True
) -> np.ndarray:
    """F(x) = (F_phi, F_Gamma, F_BL) concatenated (block_residuals)."""
    return np.concatenate(block_residuals(pack, x, veps_frozen=veps_frozen))


def augmented_jacobian_blocks(pack: TightPack, x: np.ndarray) -> Dict:
    """The named blocks of the augmented Jacobian at x (module
    docstring). The inviscid diagonal/off-diagonal (J_ff, B) come from
    ws.assemble_coupled at the current state (the external_rhs is
    additive and lagged inside eval_residual, so it does not touch the
    blocks -- the GV2.1(c) FD-exactness clause); the kutta row is the
    un-eliminated probe row [K, -I]; the coupling blocks are the
    Stage-1/2/3 assembled products with the shared column maps
    T_free/T_gam.
    """
    ws, solver = pack.ws, pack.solver
    phi_free, gamma, U = pack.split_x(x)
    phi_cut = _phi_cut_of(pack, phi_free, gamma)
    ue, q, rho, mach = _edge_data_of(pack, phi_cut)
    outs, douts = _closure_packet(pack, U)
    q2 = q * q
    m_inf, gamma_air = pack.cfg.m_inf, pack.cfg.gamma_air
    mu_arr = np.full(pack.n_s, pack.mu, dtype=np.float64)

    _, _, state = ws.eval_residual(
        phi_free, gamma, pack.upwind_c, pack.m_crit, pack.m_cap,
        pack.rho_floor,
    )
    J_ff, B = ws.assemble_coupled(
        state, pack.upwind_c, pack.m_crit, pack.rho_floor
    )

    # the phi-side transpiration augmentation at frozen (U, delta*)
    Div = surface_divergence_vector_operator(pack.sm, outs[:, C.OUT_DS1])
    drho = density_derivative_wrt_q_sq_field(q2, m_inf, gamma_air)
    Drhou = rhou_jacobian(ue, drho, pack.G, m_inf, gamma_air)
    A_cut = assemble_j_phi_phi_aug(ws, pack.W, pack.S, pack.P, Div, Drhou)
    A_free = (A_cut @ pack.T_free).tocsr()
    A_gam = (A_cut @ pack.T_gam).tocsr()

    # the BL rows at the current edge data, veps frozen (decision 5)
    try:
        solver.update_edge_data(ue, rho, mu_arr, mach)
        solver._veps, solver._veps_s = pack.veps0
        _, J_e = solver.residual_edge_jacobian(U)
        _, J_blbl = solver.residual_jacobian(U)
    finally:
        _restore_solver(pack)
    D_ue = edge_data_jacobian(ue, m_inf, gamma_air)
    J_blphi_cut = assemble_j_bl_phi(J_e, D_ue, pack.G)
    J_blphi_free = (J_blphi_cut @ pack.T_free).tocsr()
    J_blphi_gam = (J_blphi_cut @ pack.T_gam).tocsr()

    # the phi <- BL block at frozen edge data
    L = surface_divergence_delta_operator(pack.sm, rho, ue)
    D = closure_ds1_jacobian(douts)
    J_phibl = assemble_j_phi_bl(ws, pack.W, pack.S, pack.P, L, D)

    return {
        "J_ff": J_ff,
        "B": B.tocsr(),
        "K": ws.K,
        "minus_I": -sp.eye(pack.n_st, format="csr"),
        "A_cut": A_cut,
        "A_free": A_free,
        "A_gam": A_gam,
        "J_phibl": J_phibl,
        "J_blphi_cut": J_blphi_cut,
        "J_blphi_free": J_blphi_free,
        "J_blphi_gam": J_blphi_gam,
        "J_blbl": J_blbl,
    }


def augmented_jacobian(pack: TightPack, x: np.ndarray) -> sp.csr_matrix:
    """The assembled (n_free + n_st + 6 n_s)^2 augmented Jacobian at x
    (augmented_jacobian_blocks through the module-docstring layout)."""
    bl = augmented_jacobian_blocks(pack, x)
    return sp.bmat(
        [
            [bl["J_ff"] + bl["A_free"], bl["B"] + bl["A_gam"], bl["J_phibl"]],
            [bl["K"], bl["minus_I"], None],
            [bl["J_blphi_free"], bl["J_blphi_gam"], bl["J_blbl"]],
        ],
        format="csr",
    )


def _merit(pack: TightPack, F: np.ndarray) -> float:
    """The pre-registered unnormalized merit |F_phi|_2^2 + |F_Gamma|_2^2
    + |F_BL|_2^2 (the P8/P14 shared-merit idiom; module docstring)."""
    return float(F @ F)


def _block_max(pack: TightPack, F: np.ndarray) -> np.ndarray:
    return np.array(
        [float(np.max(np.abs(F_j))) if F_j.size else 0.0
         for F_j in pack.split_F(F)]
    )


# ---------------------------------------------------------------------------
# GV5.1b scaled + damped linear step (cases/analysis/v5_1b_scaled_newton/
# PRE_REGISTRATION.md): solver-internal only -- the assembled F/J are
# bit-identical to GV5.1 (the FD verdicts stand); these helpers only
# transform the assembled operators. Design inputs: the committed IBL-floor
# diagnosis (cases/analysis/v5_ibl_floor/results/findings.md).
# ---------------------------------------------------------------------------

MU0 = 1.0e-6       # Levenberg-style damping, initial value
MU_MIN = 1.0e-12   # schedule lower bound (accept: mu <- max(mu/3, MU_MIN))
MU_MAX = 1.0e2     # schedule upper bound (reject: mu <- min(10*mu, MU_MAX))
FLOOR_REL_TOL = 1.0e-4  # floor-reached: merit relative decrease below this
FLOOR_CONSEC = 3        # ... over this many consecutive accepted steps


def mu_on_accept(mu: float) -> float:
    """The pre-registered damping schedule, accepted step."""
    return max(mu / 3.0, MU_MIN)


def mu_on_reject(mu: float) -> float:
    """The pre-registered damping schedule, line-search rejection."""
    return min(mu * 10.0, MU_MAX)


class FloorStop:
    """The pre-registered floor-reached stop: the merit relative decrease
    stays below FLOOR_REL_TOL over FLOOR_CONSEC consecutive accepted
    steps (replaces GV5.1's lambda-collapse crawl to the iteration
    cap)."""

    def __init__(self) -> None:
        self._hits = 0

    def update(self, rel_decrease: float) -> bool:
        self._hits = (self._hits + 1
                      if rel_decrease < FLOOR_REL_TOL else 0)
        return self._hits >= FLOOR_CONSEC


def equilibrate_rc(J: sp.csr_matrix):
    """One-pass row/column 2-norm equilibration of the assembled sparse J
    (zero-safe: a zero row/column gets scale 1 -- the pre-registered
    recipe, the diagnosis Q3 method). Returns (rn, cn, Jsc) with
    Jsc = diag(1/rn) @ J @ diag(1/cn)."""
    J = J.tocsr()
    rn = np.sqrt(np.asarray(J.multiply(J).sum(axis=1)).ravel())
    rn = np.where(rn > 0.0, rn, 1.0)
    J1 = sp.diags(1.0 / rn) @ J
    cn = np.sqrt(np.asarray(J1.multiply(J1).sum(axis=0)).ravel())
    cn = np.where(cn > 0.0, cn, 1.0)
    return rn, cn, (J1 @ sp.diags(1.0 / cn)).tocsr()


def scaled_damped_step(J: sp.csr_matrix, F: np.ndarray, mu: float,
                       scaling: Optional[str] = None) -> np.ndarray:
    """(J~ + mu I) dy = -F~ via splu (sparsity preserved, no normal
    equations), unscaled dx = C dy. scaling="rowcol": J~ = R J C and
    F~ = R F (equilibrate_rc); scaling=None: J~ = J, F~ = F, C = I, so
    mu = 0 recovers the undamped splu step exactly."""
    if scaling == "rowcol":
        rn, cn, Js = equilibrate_rc(J)
        Fs = F / rn
    elif scaling is None:
        cn, Js, Fs = None, J.tocsr(), F
    else:
        raise ValueError(f"unknown scaling {scaling!r}")
    n = Js.shape[0]
    d_y = -sla.splu(
        (Js.tocsc() + mu * sp.eye(n, format="csc")).tocsc()
    ).solve(Fs)
    return d_y if cn is None else d_y / cn


def newton_tight(
    pack: TightPack,
    x0: Optional[np.ndarray] = None,
    tol: float = 1.0e-8,
    tol_abs: float = 1.0e-10,
    max_iter: int = 10,
    line_search: bool = True,
    max_backtracks: int = 30,
    verbose: bool = False,
    scaling: Optional[str] = None,
    lm_damping: bool = False,
    floor_stop: bool = False,
) -> Dict:
    """The augmented Newton loop on the full state (PRE_REGISTRATION
    convergence protocol): seed = the pack base x0 (inviscid-converged
    (phi, Gamma) + ONE IBL solve), splu on the assembled CSR, safety-only
    backtracking on the unnormalized augmented merit
    |F_phi|^2 + |F_Gamma|^2 + |F_BL|^2 (the P8/P14 accept-or-least-bad
    idiom, solve/newton.py:1053-1087, with the halving budget raised
    from 5 to max_backtracks=30 for the stalled-seed BL step scale --
    module docstring). Backtracking PROBES are guarded the IBL3Solver
    way (ibl3.py:1497-1513): a probe that throws (a nonphysical x_try
    can divide by zero inside the closure quadrature -- observed on the
    GV5.1 medium polish, iter 6) or that returns a non-finite residual
    or merit counts as merit=+inf and simply keeps halving; only probe
    rejection is affected, never an accepted step.

    Convergence (module docstring): every block's max-norm <= tol_abs +
    tol * |F_j(x0)|_inf. Per-block norms are tracked separately; the
    delta*-change of the LAST accepted step is recorded against the
    tol_ds = 1e-3 cross-check. Returns a dict(x, converged, n_iter,
    history, block_max0, ds_change_last, merit0, termination); history
    rows carry (iter, block_max, merit, lam, ds_change, wall_s, mu) with
    wall_s the cumulative seconds since the loop start. termination is
    "converged" or "cap" on the legacy path, plus "floor_reached" on the
    GV5.1b path.

    GV5.1b (cases/analysis/v5_1b_scaled_newton/PRE_REGISTRATION.md):
    scaling="rowcol" + lm_damping=True + floor_stop=True switch the
    linear step to the scaled + damped path -- per iteration R, C from
    the current J (equilibrate_rc), (R J C + mu I) dy = -R F by splu,
    dx = C dy; mu follows the deterministic schedule (MU0; x10 on a
    line-search rejection with a retry, /3 on an accepted step, bounded
    [MU_MIN, MU_MAX]); the P8/P14 backtracking + probe guard applies
    unchanged on each damped trial step; floor_stop adds the
    pre-registered floor-reached termination (FloorStop). The DEFAULTS
    (scaling=None, lm_damping=False, floor_stop=False) keep the legacy
    path bit-for-bit (the committed GV5.1 runner reproduces).
    """
    x = pack.x_base() if x0 is None else np.asarray(x0, dtype=np.float64).copy()
    t0 = time.perf_counter()
    F = augmented_residual(pack, x)
    block_max0 = _block_max(pack, F)
    merit = _merit(pack, F)
    merit0 = merit

    def ds_at(x_):
        # the system's delta* (frozen-edge closure; the tol_ds cross-check
        # tracks the state movement, not the closure's edge chain)
        _, _, U = pack.split_x(x_)
        outs, _ = _closure_packet(pack, U)
        return outs[:, C.OUT_DS1].copy()

    ds_prev = ds_at(x)
    history = []
    converged = False
    ds_change_last = 0.0
    n_it_done = 0
    termination = None

    def converged_at(F_):
        bm = _block_max(pack, F_)
        return bool(np.all(bm <= tol_abs + tol * block_max0)), bm

    converged, bm = converged_at(F)
    history.append(
        {"iter": 0, "block_max": bm, "merit": merit, "lam": None,
         "ds_change": 0.0, "wall_s": time.perf_counter() - t0, "mu": None}
    )
    if verbose:
        print(
            f"newton_tight iter 0: |F_phi|={bm[0]:.3e} |F_Gam|={bm[1]:.3e} "
            f"|F_BL|={bm[2]:.3e} merit={merit:.3e}"
        )

    if scaling is None and not lm_damping and not floor_stop:
        # ---- the legacy GV5.1 path (bit-for-bit; the committed runner
        # reproduces) -------------------------------------------------
        for it in range(1, max_iter + 1):
            if converged:
                termination = "converged"
                break
            n_it_done = it
            J = augmented_jacobian(pack, x)
            try:
                d = -sla.splu(J.tocsc()).solve(F)
            except RuntimeError as exc:
                raise RuntimeError(
                    f"tight Newton factorization failed at iteration {it}: {exc}"
                ) from exc
            if not np.all(np.isfinite(d)):
                raise RuntimeError(
                    f"tight Newton step non-finite at iteration {it}"
                )
            # P8/P14 safety-only backtracking (module docstring). A far
            # probe can be nonphysical enough to THROW inside the closure
            # quadrature -- the ZeroDivisionError reaches the jit boundary,
            # where numba's CPUDispatcher re-raises it as SystemError
            # ("returned a result with an exception set"), so BOTH are
            # caught -- or to return non-finite values: count it as
            # merit=+inf and keep halving -- the IBL3Solver.solve
            # halving-on-nonfinite idiom (ibl3.py:1497-1513, where a NaN
            # probe simply fails the decrease test) and the P8/P14
            # globalization precedent. Only probe REJECTION is guarded;
            # accepted steps are untouched.
            lam = 1.0
            best = None
            accepted = False
            for _ in range(max_backtracks):
                x_try = x + lam * d
                try:
                    F_try = augmented_residual(pack, x_try)
                    merit_try = _merit(pack, F_try)
                except (ArithmeticError, SystemError):
                    F_try, merit_try = None, np.inf
                if F_try is not None and not np.all(np.isfinite(F_try)):
                    F_try, merit_try = None, np.inf
                if np.isfinite(merit_try) and (
                    not line_search or merit_try <= merit
                ):
                    accepted = True
                    break
                if np.isfinite(merit_try) and (
                    best is None or merit_try < best[0]
                ):
                    best = (merit_try, lam, x_try, F_try)
                lam *= 0.5
            if not accepted:
                if best is None:
                    raise RuntimeError(
                        f"tight Newton line search failed at iteration {it} "
                        f"(merit {merit:.3e})"
                    )
                merit_try, lam, x_try, F_try = best
            x, F, merit = x_try, F_try, merit_try
            ds_new = ds_at(x)
            ds_change_last = float(np.max(np.abs(ds_new - ds_prev)))
            ds_prev = ds_new
            converged, bm = converged_at(F)
            history.append(
                {"iter": it, "block_max": bm, "merit": merit,
                 "lam": float(lam), "ds_change": ds_change_last,
                 "wall_s": time.perf_counter() - t0, "mu": None}
            )
            if verbose:
                print(
                    f"newton_tight iter {it}: |F_phi|={bm[0]:.3e} "
                    f"|F_Gam|={bm[1]:.3e} |F_BL|={bm[2]:.3e} "
                    f"merit={merit:.3e} lam={lam:.2f} "
                    f"ds_change={ds_change_last:.3e}"
                )
    else:
        # ---- the GV5.1b scaled + damped path (pre-registered) ---------
        mu = MU0
        floor = FloorStop() if floor_stop else None
        for it in range(1, max_iter + 1):
            if converged:
                termination = "converged"
                break
            n_it_done = it
            J = augmented_jacobian(pack, x)
            merit_prev = merit
            mu_retries = 0
            while True:
                d = scaled_damped_step(
                    J, F, mu if lm_damping else 0.0, scaling=scaling)
                if not np.all(np.isfinite(d)):
                    raise RuntimeError(
                        f"tight Newton step non-finite at iteration {it}"
                    )
                # the P8/P14 backtracking + probe guard, unchanged, on
                # each damped trial step (the legacy idiom above)
                lam = 1.0
                best = None
                accepted = False
                for _ in range(max_backtracks):
                    x_try = x + lam * d
                    try:
                        F_try = augmented_residual(pack, x_try)
                        merit_try = _merit(pack, F_try)
                    except (ArithmeticError, SystemError):
                        F_try, merit_try = None, np.inf
                    if F_try is not None and not np.all(np.isfinite(F_try)):
                        F_try, merit_try = None, np.inf
                    if np.isfinite(merit_try) and (
                        not line_search or merit_try <= merit
                    ):
                        accepted = True
                        break
                    if np.isfinite(merit_try) and (
                        best is None or merit_try < best[0]
                    ):
                        best = (merit_try, lam, x_try, F_try)
                    lam *= 0.5
                if accepted or not lm_damping:
                    break
                # a line-search rejection in the damped path: raise mu and
                # retry the damped solve (the pre-registered schedule); at
                # the cap, fall back to the legacy accept-or-least-bad
                if mu >= MU_MAX:
                    if best is None:
                        raise RuntimeError(
                            f"tight Newton line search failed at iteration "
                            f"{it} (merit {merit:.3e}, mu at cap {MU_MAX:.0e})"
                        )
                    break
                mu = mu_on_reject(mu)
                mu_retries += 1
            mu_used = mu
            if not accepted:
                if best is None:
                    raise RuntimeError(
                        f"tight Newton line search failed at iteration {it} "
                        f"(merit {merit:.3e})"
                    )
                # the least-bad fallback: a REJECTED step taken per the
                # legacy idiom -- mu stays (no accept-decrease)
                merit_try, lam, x_try, F_try = best
            elif lm_damping:
                mu = mu_on_accept(mu)
            x, F, merit = x_try, F_try, merit_try
            ds_new = ds_at(x)
            ds_change_last = float(np.max(np.abs(ds_new - ds_prev)))
            ds_prev = ds_new
            converged, bm = converged_at(F)
            history.append(
                {"iter": it, "block_max": bm, "merit": merit,
                 "lam": float(lam), "ds_change": ds_change_last,
                 "wall_s": time.perf_counter() - t0, "mu": float(mu_used),
                 "mu_retries": mu_retries,
                 "accepted": bool(accepted)}
            )
            if verbose:
                print(
                    f"newton_tight iter {it}: |F_phi|={bm[0]:.3e} "
                    f"|F_Gam|={bm[1]:.3e} |F_BL|={bm[2]:.3e} "
                    f"merit={merit:.3e} lam={lam:.2f} mu={mu_used:.1e} "
                    f"retries={mu_retries} ds_change={ds_change_last:.3e}"
                )
            rel_dec = ((merit_prev - merit)
                       / max(abs(merit_prev), 1.0e-300))
            if not converged and floor is not None and floor.update(rel_dec):
                termination = "floor_reached"
                break
    if termination is None:
        termination = "converged" if converged else "cap"
    return {
        "x": x,
        "converged": converged,
        "n_iter": n_it_done,
        "history": history,
        "block_max0": block_max0,
        "ds_change_last": ds_change_last,
        "merit0": merit0,
        "termination": termination,
    }
