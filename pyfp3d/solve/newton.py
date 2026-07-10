"""
Fully-coupled (phi_red, Gamma) Newton driver (design.md Sec 8.1, roadmap
P8): the exact Jacobian (6.3) at frozen upstream selection, with the
circulation solved as an unknown alongside the potential -- replacing the
P5-fragile Gamma-secant (whose secant-density coupling was the P5 medium
root cause) with one coupled system.

Unknowns x = (phi_free, Gamma). ONE state-reconstruction + residual path
(eval_residual) is shared by the driver, the Jacobian assembly inputs and
the FD tests -- residual/Jacobian consistency is what buys the quadratic
convergence:

    vals_red(Gamma) = vals0_red + V_red @ Gamma       (far field, LINEAR)
    phi_red[free] = phi_free; phi_red[dir] = vals_red(Gamma)
    phi_cut = T phi_red + g(Gamma)                    (wake jump)
    density chain exactly as solve_subsonic_lifting (velocities ->
        limit_q2_field -> density_field -> walk rho_tilde)
    R_free = (T^T assemble_residual)[free]
    F      = kutta_targets(phi_cut) - Gamma

Newton block system and its EXACT elimination (the (2,2) block is -I):

    [ J_ff   B  ] [dphi]   [-R_free]        dGamma = K dphi + F
    [ K     -I  ] [dGam] = [-F     ]  ==>   (J_ff + B K) dphi = -R_free - B F

with  J_ff = (T^T J T)[free, free]
      B    = (T^T J T)[free, dir] @ V_red  +  (T^T J G)[free, :]
             ^far-field vortex column          ^wake-jump column H_J
      K    = dF/dphi_free (sparse +-1/n_j at the TE Kutta probes).

The far-field column is the roadmap-flagged easy-to-miss term: in the
Picard loop it is folded silently into the RHS through A_coupling; the
Newton Gamma-column must carry it explicitly (with the FULL J, not A).

Linear solve: GMRES on the eliminated low-rank operator, right-hand
preconditioned by AMG built on the SPD Picard (Term-1) block -- the exact
Jacobian is nonsymmetric/indefinite in supersonic zones, so CG does not
apply. Eisenstat-Walker (choice 2) forcing keeps early Newton steps cheap.

Globalization: plain full Newton steps (Lopez) with a SAFETY-ONLY
backtracking line search on |R|^2 + |F|^2; no damping_theta anywhere (that
is a Picard stabilizer -- it would destroy quadratic convergence); an
optional consistent pseudo-transient term (ptc_dtau) is the fallback for
hard transients. Transonic runs wrap this in the upward-only Mach
continuation of solve_newton_transonic (M_crit and upwind_c held FIXED
within the ramp, Lopez Tables 4.7/4.8/4.13).
"""

import time
from typing import Dict, Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from pyfp3d.constraints.dirichlet import farfield_dirichlet, freestream_phi
from pyfp3d.constraints.wake import WakeConstraint, kutta_targets
from pyfp3d.kernels.jacobian import PicardOperator
from pyfp3d.kernels.upwind import UpwindOperator
from pyfp3d.physics.isentropic import (
    density_field,
    limit_q2_field,
    mach_squared_field,
)
from pyfp3d.solve.linear import (
    build_amg_preconditioner,
    build_ilu_preconditioner,
    solve_gmres,
)


class NewtonWorkspace:
    """Per-case precomputation for the coupled Newton solve: operators,
    the free/Dirichlet split, the Kutta row K, and the per-Mach-level
    far-field basis (vals0_red, V_red). Reused across the Mach levels of
    solve_newton_transonic (set_mach re-derives only the beta-dependent
    far-field basis)."""

    def __init__(self, mesh_cut, wc, alpha_deg: float = 0.0,
                 u_inf: float = 1.0, gamma_air: float = 1.4,
                 vortex_center=(0.25, 0.0),
                 farfield_spanwise_gamma: bool = False):
        self.mesh_cut = mesh_cut
        self.wc = wc
        self.alpha_deg = float(alpha_deg)
        self.u_inf = float(u_inf)
        self.gamma_air = float(gamma_air)
        self.vortex_center = vortex_center
        self.spanwise = bool(farfield_spanwise_gamma)

        self.op = PicardOperator(mesh_cut.nodes, mesh_cut.elements)
        self.upw = UpwindOperator(mesh_cut.nodes, mesh_cut.elements,
                                  weighted=False)
        self.con = WakeConstraint(self.op.assemble_matrix(), wc)
        self.n_red = self.con.n_reduced
        self.n_st = wc.n_stations

        # Dirichlet pattern is Gamma-, rho- and beta-independent: fix the
        # split once (identical construction to solve_subsonic_lifting).
        dir_nodes, _ = farfield_dirichlet(
            mesh_cut, wc, alpha_deg, np.zeros(self.n_st), u_inf,
            vortex_center, beta=1.0, spanwise_gamma=self.spanwise,
        )
        dir_red, _ = self.con.to_reduced_dirichlet(
            dir_nodes, np.zeros(len(dir_nodes)))
        is_dir = np.zeros(self.n_red, dtype=bool)
        is_dir[dir_red] = True
        self.dir_red = dir_red
        self.free = np.where(~is_dir)[0]
        self.n_free = len(self.free)
        # cut-node -> reduced-value mapping for far-field value arrays
        # (mirrors to_reduced_dirichlet: drop slaves, unique-sort); the
        # ff_nodes ordering from farfield_dirichlet is deterministic, so
        # one keep/idx pair reduces every value array consistently with
        # dir_red's ordering.
        keep = dir_nodes < self.n_red
        _, idx = np.unique(dir_nodes[keep], return_index=True)
        self._ff_keep = keep
        self._ff_idx = idx

        # K = dF/dphi_free: station row j holds +1/n_j at that station's
        # upper TE probes and -1/n_j at the lower (kutta_targets is the
        # per-station MEAN of probe jumps, so shared probes -- adjacent
        # stations reusing a node, the P5-recorded robustness item --
        # simply appear in both rows, which IS the derivative).
        counts = np.bincount(wc.te_station, minlength=self.n_st).astype(
            np.float64)
        probe_nodes = np.concatenate([wc.kutta_upper, wc.kutta_lower])
        assert np.all(probe_nodes < self.n_red), (
            "Kutta probes must be original (non-slave) nodes")
        assert not is_dir[probe_nodes].any(), (
            "Kutta probes must be free dofs")
        cols_free = np.searchsorted(self.free, probe_nodes)
        assert np.array_equal(self.free[cols_free], probe_nodes)
        rows = np.concatenate([wc.te_station, wc.te_station])
        w = 1.0 / counts[wc.te_station]
        vals = np.concatenate([w, -w])
        self.K = sp.coo_matrix(
            (vals, (rows, cols_free)), shape=(self.n_st, self.n_free),
        ).tocsr()

        self.m_inf = None
        self.beta = None
        self.vals0_red = None
        self.V_red = None
        # lumped element volumes on the free reduced dofs (ptc option)
        m_lumped = np.zeros(self.op.n_nodes, dtype=np.float64)
        np.add.at(m_lumped, np.asarray(self.op.elements).reshape(-1),
                  np.repeat(self.op.V / 4.0, 4))
        self.m_lumped_free = (self.con.T.T @ m_lumped)[self.free]

    def _reduce_ff_values(self, values: np.ndarray) -> np.ndarray:
        return np.asarray(values, dtype=np.float64)[self._ff_keep][self._ff_idx]

    def set_mach(self, m_inf: float) -> None:
        """Fix the Mach level: Prandtl-Glauert beta and the affine
        far-field basis vals_red(Gamma) = vals0_red + V_red @ Gamma.
        farfield_dirichlet is exactly linear in the station Gammas (the
        vortex potential is linear in its strength and the spanwise
        interpolant is linear in the ordinates), so unit-Gamma probing
        recovers the exact columns; test_farfield_values_linear_in_gamma
        guards the linearity assumption."""
        if not 0.0 <= m_inf < 1.0:
            raise ValueError(f"Newton driver needs 0 <= M_inf < 1, got {m_inf}")
        self.m_inf = float(m_inf)
        self.beta = float(np.sqrt(1.0 - m_inf ** 2))
        _, vals0 = farfield_dirichlet(
            self.mesh_cut, self.wc, self.alpha_deg, np.zeros(self.n_st),
            self.u_inf, self.vortex_center, beta=self.beta,
            spanwise_gamma=self.spanwise,
        )
        self.vals0_red = self._reduce_ff_values(vals0)
        V = np.empty((len(self.dir_red), self.n_st), dtype=np.float64)
        for j in range(self.n_st):
            e_j = np.zeros(self.n_st)
            e_j[j] = 1.0
            _, vals_j = farfield_dirichlet(
                self.mesh_cut, self.wc, self.alpha_deg, e_j, self.u_inf,
                self.vortex_center, beta=self.beta,
                spanwise_gamma=self.spanwise,
            )
            V[:, j] = self._reduce_ff_values(vals_j) - self.vals0_red
        self.V_red = V
        self._V_red_sp = sp.csr_matrix(V)

    def eval_residual(self, phi_free: np.ndarray, gamma: np.ndarray,
                      upwind_c: float, m_crit: float, m_cap: float,
                      rho_floor: float):
        """The single state-reconstruction + residual path (module
        docstring). Returns (R_free, F, state) with state caching
        everything the Jacobian assembly at this point needs."""
        gamma = np.asarray(gamma, dtype=np.float64)
        phi_red = np.empty(self.n_red, dtype=np.float64)
        phi_red[self.free] = phi_free
        phi_red[self.dir_red] = self.vals0_red + self.V_red @ gamma
        phi_cut = self.con.expand(phi_red, gamma)

        grad, q2 = self.op.velocities(phi_cut)
        grad = grad.copy()
        q2n = q2 / self.u_inf ** 2
        q2l = limit_q2_field(q2n, self.m_inf, m_cap, self.gamma_air)
        lim = q2l == q2n
        n_limited = int(np.count_nonzero(~lim))
        rho = density_field(q2l, self.m_inf, self.gamma_air)
        if self.m_inf > 0.0:
            rho_t = self.upw.rho_tilde(
                grad, q2l, rho, self.m_inf, upwind_c, m_crit,
                self.gamma_air, rho_floor).copy()
            n_floored = self.upw.n_floored
            nu_max = self.upw.nu_max
            n_nu_active = self.upw.n_supersonic
        else:
            rho_t = rho
            n_floored = 0
            nu_max = 0.0
            n_nu_active = 0

        R_red = self.con.T.T @ self.op.assemble_residual(phi_cut, rho_t)
        R_free = R_red[self.free]
        F = kutta_targets(phi_cut, self.wc) - gamma
        state = {
            "phi_red": phi_red, "phi_cut": phi_cut, "grad": grad,
            "q2l": q2l, "lim": lim, "rho": rho, "rho_t": rho_t,
            "n_limited": n_limited, "n_floored": n_floored,
            "nu_max": nu_max, "n_nu_active": n_nu_active,
        }
        return R_free, F, state

    def assemble_coupled(self, state, upwind_c: float, m_crit: float,
                         rho_floor: float):
        """Exact Jacobian blocks at the state returned by eval_residual:
        (J_ff, B) with B = J_red[free, dir] @ V_red + H_J[free, :]."""
        if self.m_inf > 0.0:
            s_e, s_u, upstream = self.upw.rho_tilde_sensitivities(
                state["grad"], state["q2l"], state["rho"], self.m_inf,
                upwind_c, m_crit, self.gamma_air, rho_floor)
        else:
            n_tets = self.op.n_tets
            s_e = np.zeros(n_tets)
            s_u = np.zeros(n_tets)
            upstream = np.arange(n_tets, dtype=np.int64)
        lim = state["lim"]
        J = self.op.assemble_newton_jacobian(
            state["phi_cut"], state["rho_t"], s_e, s_u, upstream,
            u_inf=self.u_inf, lim_mask=None if lim.all() else lim)
        J_red, H_J = self.con.reduce_operator(J)
        J_ff = J_red[self.free][:, self.free].tocsr()
        J_fd = J_red[self.free][:, self.dir_red].tocsr()
        B = (J_fd @ self._V_red_sp + H_J.tocsr()[self.free, :]).tocsc()
        return J_ff, B


def _ew_forcing(r_norm, r_norm_prev, eta_prev, eta0=1e-2, gamma_ew=0.9,
                eta_max=1e-2, eta_floor=1e-10):
    """Eisenstat-Walker choice 2 with the standard safeguard."""
    if r_norm_prev is None:
        return eta0
    eta = gamma_ew * (r_norm / r_norm_prev) ** 2
    safeguard = gamma_ew * eta_prev ** 2
    if safeguard > 0.1:
        eta = max(eta, safeguard)
    return float(min(eta_max, max(eta, eta_floor)))


def solve_newton_lifting(
    mesh_cut,
    wc,
    m_inf: float,
    alpha_deg: float = 0.0,
    u_inf: float = 1.0,
    gamma_air: float = 1.4,
    vortex_center=(0.25, 0.0),
    farfield_spanwise_gamma: bool = False,
    upwind_c: float = 1.5,
    m_crit: float = 0.95,
    m_cap: float = 3.0,
    rho_floor: float = 0.05,
    phi_init: Optional[np.ndarray] = None,
    gamma_init: Optional[np.ndarray] = None,
    n_picard_seed: int = 5,
    n_newton_max: int = 30,
    tol_residual: float = 1e-10,
    tol_gamma: float = 1e-8,
    ew_eta0: float = 1e-2,
    ew_gamma: float = 0.9,
    ew_eta_max: float = 1e-2,
    amg_rebuild_every: int = 2,
    precond: str = "amg",
    gmres_restart: int = 60,
    gmres_maxiter: int = 50,
    line_search: bool = True,
    ptc_dtau: Optional[float] = None,
    rtol_seed: float = 1e-7,
    workspace: Optional[NewtonWorkspace] = None,
    verbose: bool = False,
) -> Dict[str, object]:
    """
    Fully-coupled (phi_red, Gamma) Newton solve at ONE Mach level (module
    docstring). Subsonic: converges from a short Picard warm start (or
    freestream) in a handful of quadratic steps. Transonic levels are
    driven through solve_newton_transonic (warm-started, upward ramp).

    The upwind machinery is always active for m_inf > 0 (subcritically it
    is exactly rho_tilde == rho and s_u == 0, the G4.2 bit-identity
    regime, so there is one code path); upwind_c must be > 0 then. At
    m_inf = 0 the problem is linear and one Newton step reproduces the P2
    lifting Laplace solution.

    `n_picard_seed` > 0 runs that many Picard outer iterations (loose
    tol) purely as an initial guess; `phi_init`/`gamma_init` (cut-mesh
    phi) take precedence and skip the seed. `ptc_dtau` adds the CONSISTENT
    pseudo-transient diag(m_lumped_free/dtau) to J_ff only -- it
    multiplies dphi and vanishes at convergence, so the converged state is
    exactly the Newton state (globalization fallback; default off; note it
    breaks terminal quadratic convergence while active). `precond` "amg"
    builds the hierarchy on the SPD Picard block (rebuilt every
    `amg_rebuild_every` Newton steps), "ilu" factors J_ff itself.

    Convergence: ||R_free||_inf < tol_residual AND ||F||_inf < tol_gamma,
    refused while any element is speed-limited or density-floored (a
    clamped state is not a converged flow -- mirrors the P5 gate
    semantics).

    Returns the Picard-compatible result keys (phi on the cut mesh, gamma,
    converged, residual_history, mach2_max, nu_max, n_nu_active,
    n_limited, n_floored) plus F_history, newton_orders (per-step observed
    order p_k), eta_history, n_gmres_total, n_newton, jacobian_nnz,
    n_term3_active, timings (per-stage wall-clock seconds).
    """
    if m_inf > 0.0 and upwind_c <= 0.0:
        raise ValueError(
            "the Newton driver runs with the walk upwind machinery active "
            "for m_inf > 0 (bit-inert subcritically); pass upwind_c > 0")

    ws = workspace
    if ws is None:
        ws = NewtonWorkspace(mesh_cut, wc, alpha_deg, u_inf, gamma_air,
                             vortex_center, farfield_spanwise_gamma)
    ws.set_mach(m_inf)

    timings = {"seed": 0.0, "residual": 0.0, "jacobian": 0.0,
               "amg_setup": 0.0, "gmres": 0.0}

    # ---- initial state -------------------------------------------------
    t0 = time.perf_counter()
    if phi_init is not None:
        phi_red0 = np.asarray(phi_init, dtype=np.float64)[:ws.n_red].copy()
        gamma = (np.zeros(ws.n_st) if gamma_init is None
                 else np.atleast_1d(np.asarray(gamma_init,
                                               dtype=np.float64)).copy())
    elif n_picard_seed > 0:
        from pyfp3d.solve.picard import solve_subsonic_lifting

        seed = solve_subsonic_lifting(
            mesh_cut, wc, m_inf=m_inf, alpha_deg=alpha_deg, u_inf=u_inf,
            gamma_air=gamma_air, vortex_center=vortex_center,
            upwind_c=upwind_c, m_crit=m_crit, m_cap=m_cap,
            n_picard_max=n_picard_seed, tol_rho=1e-3, rtol=rtol_seed,
            farfield_spanwise_gamma=farfield_spanwise_gamma,
        )
        phi_red0 = np.asarray(seed["phi"], dtype=np.float64)[:ws.n_red].copy()
        gamma = np.atleast_1d(np.asarray(seed["gamma"],
                                         dtype=np.float64)).copy()
    else:
        phi_red0 = freestream_phi(mesh_cut.nodes[:ws.n_red], alpha_deg,
                                  u_inf)
        gamma = (np.zeros(ws.n_st) if gamma_init is None
                 else np.atleast_1d(np.asarray(gamma_init,
                                               dtype=np.float64)).copy())
    phi_free = phi_red0[ws.free].copy()
    timings["seed"] = time.perf_counter() - t0

    # ---- Newton loop ---------------------------------------------------
    residual_history = []
    F_history = []
    eta_history = []
    newton_orders = []
    n_gmres_total = 0
    converged = False
    M_pre = None
    eta_prev = None

    t0 = time.perf_counter()
    R_free, F, state = ws.eval_residual(phi_free, gamma, upwind_c, m_crit,
                                        m_cap, rho_floor)
    timings["residual"] += time.perf_counter() - t0
    merit = float(R_free @ R_free + F @ F)

    for it in range(n_newton_max):
        r_norm = float(np.max(np.abs(R_free)))
        f_norm = float(np.max(np.abs(F))) if ws.n_st > 0 else 0.0
        residual_history.append(r_norm)
        F_history.append(f_norm)
        if len(residual_history) >= 3:
            r2, r1, r0 = (residual_history[-1], residual_history[-2],
                          residual_history[-3])
            if r2 > 0.0 and r1 > 0.0 and r0 > 0.0 and r1 != r0:
                newton_orders.append(
                    float(np.log(r2 / r1) / np.log(r1 / r0)))
        if verbose:
            print(f"  newton {it:2d}: |R|={r_norm:.3e} |F|={f_norm:.3e} "
                  f"lim={state['n_limited']} flr={state['n_floored']}")
        if (r_norm < tol_residual and f_norm < tol_gamma
                and state["n_limited"] == 0 and state["n_floored"] == 0):
            converged = True
            break

        t0 = time.perf_counter()
        J_ff, B = ws.assemble_coupled(state, upwind_c, m_crit, rho_floor)
        if ptc_dtau is not None:
            J_ff = (J_ff + sp.diags(ws.m_lumped_free / ptc_dtau)).tocsr()
        timings["jacobian"] += time.perf_counter() - t0

        t0 = time.perf_counter()
        if precond == "amg":
            if M_pre is None or it % amg_rebuild_every == 0:
                A_pic = ws.op.assemble_matrix(state["rho_t"])
                A_red = (ws.con.T.T @ (A_pic @ ws.con.T)).tocsr()
                A_ff = A_red[ws.free][:, ws.free].tocsr()
                if ptc_dtau is not None:
                    A_ff = (A_ff
                            + sp.diags(ws.m_lumped_free / ptc_dtau)).tocsr()
                _, M_pre = build_amg_preconditioner(A_ff)
        elif precond == "ilu":
            M_pre = build_ilu_preconditioner(J_ff)
        else:
            raise ValueError(f"unknown precond {precond!r}")
        timings["amg_setup"] += time.perf_counter() - t0

        K = ws.K
        A_op = spla.LinearOperator(
            (ws.n_free, ws.n_free),
            matvec=lambda x: J_ff @ x + B @ (K @ x),
        )
        rhs = -R_free - B @ F
        eta = _ew_forcing(r_norm, residual_history[-2] if it > 0 else None,
                          eta_prev, ew_eta0, ew_gamma, ew_eta_max)
        eta_prev = eta
        eta_history.append(eta)
        t0 = time.perf_counter()
        dphi, n_it = solve_gmres(A_op, rhs, M=M_pre, rtol=eta,
                                 restart=gmres_restart,
                                 maxiter=gmres_maxiter)
        timings["gmres"] += time.perf_counter() - t0
        n_gmres_total += n_it
        dgamma = K @ dphi + F

        # safety-only backtracking: accept the full step unless the merit
        # blows up or goes non-finite (Lopez runs plain Newton; the load
        # stepping is the real globalization).
        lam = 1.0
        accepted = False
        for _ in range(5):
            t0 = time.perf_counter()
            R_try, F_try, state_try = ws.eval_residual(
                phi_free + lam * dphi, gamma + lam * dgamma,
                upwind_c, m_crit, m_cap, rho_floor)
            timings["residual"] += time.perf_counter() - t0
            merit_try = float(R_try @ R_try + F_try @ F_try)
            if np.isfinite(merit_try) and (not line_search
                                           or merit_try <= merit
                                           or lam <= 0.0625):
                accepted = True
                break
            lam *= 0.5
        if not accepted:                     # pragma: no cover - safety net
            raise RuntimeError(
                f"Newton line search failed at iteration {it} "
                f"(merit {merit:.3e} -> {merit_try:.3e})")
        phi_free = phi_free + lam * dphi
        gamma = gamma + lam * dgamma
        R_free, F, state = R_try, F_try, state_try
        merit = merit_try

    q2n = state["q2l"]
    mach2_max = float(np.max(mach_squared_field(q2n, m_inf, gamma_air)))
    return {
        "phi": state["phi_cut"],
        "gamma": gamma,
        "converged": converged,
        "n_newton": len(residual_history) - (1 if converged else 0),
        "residual_history": residual_history,
        "F_history": F_history,
        "newton_orders": newton_orders,
        "eta_history": eta_history,
        "n_gmres_total": n_gmres_total,
        "mach2_max": mach2_max,
        "nu_max": state["nu_max"],
        "n_nu_active": state["n_nu_active"],
        "n_limited": state["n_limited"],
        "n_floored": state["n_floored"],
        "jacobian_nnz": getattr(ws.op, "newton_nnz", None),
        "n_term3_active": getattr(ws.op, "n_term3_active", None),
        "timings": timings,
        "workspace": ws,
    }


def solve_newton_transonic(
    mesh_cut,
    wc,
    m_inf: float,
    alpha_deg: float,
    m_start: float = 0.70,
    dm: float = 0.05,
    dm_min: float = 0.01,
    upwind_c: float = 1.5,
    m_crit: float = 0.95,
    upwind_c_post=None,
    newton_kw: Optional[dict] = None,
    verbose: bool = False,
) -> Dict[str, object]:
    """
    Transonic Newton solve by upward Mach continuation (design.md Sec 8.1;
    Lopez Sec 4.4/4.8). MECHANICAL SKELETON delivered with the P8
    subsonic milestone -- the transonic tuning (level schedules, EW/PTC
    parameters against G8.1) is the N5 work package; interfaces are final.

    Level 0 (subcritical m_start) seeds with the standard Picard warm
    start; every later level warm-starts (phi, Gamma) from the previous
    level -- ALWAYS ramping M_inf upward from a subcritical converged
    state (design.md Sec 12 risk 2: continuation direction selects the FP
    solution branch at M ~ 0.82-0.85). Within the ramp M_crit and
    upwind_c are held FIXED (Lopez Tables 4.7/4.8); `upwind_c_post` is an
    optional list of DECREASING upwind_c values applied at the final Mach
    only (Lopez Table 4.13's 2.0 -> 1.6 dissipation sharpening). A level
    that fails to converge is retried with the Mach step halved (down to
    dm_min).

    Returns the final-level result dict plus level_history (per-level
    (m, n_newton, |R|_final)).
    """
    from pyfp3d.solve.continuation import mach_schedule

    if m_inf < m_start:
        raise ValueError(
            f"upward ramp only: m_inf {m_inf} < m_start {m_start} "
            "(design.md Sec 12 risk 2 -- never continue downward)")
    kw = dict(newton_kw or {})
    kw.setdefault("upwind_c", upwind_c)
    kw.setdefault("m_crit", m_crit)

    ws = None
    result = None
    last_good = None
    level_history = []
    levels = mach_schedule(m_inf, m_start, dm)
    i = 0
    while i < len(levels):
        m = levels[i]
        lvl_kw = dict(kw)
        if last_good is not None:
            # warm-start from the last CONVERGED level (a failed level's
            # state is not a valid continuation seed)
            lvl_kw.update(phi_init=last_good["phi"],
                          gamma_init=last_good["gamma"], n_picard_seed=0)
        if verbose:
            print(f"newton continuation level M = {m:.4f}")
        result = solve_newton_lifting(mesh_cut, wc, m_inf=m,
                                      alpha_deg=alpha_deg, workspace=ws,
                                      verbose=verbose, **lvl_kw)
        ws = result["workspace"]
        level_history.append((m, result["n_newton"],
                              result["residual_history"][-1]))
        if not result["converged"]:
            m_prev = m_start if i == 0 else levels[i - 1]
            dm_new = 0.5 * (m - m_prev)
            if i == 0 or dm_new < dm_min:
                break
            levels.insert(i, m_prev + dm_new)
            continue
        last_good = result
        i += 1

    if result["converged"] and upwind_c_post:
        for c_post in upwind_c_post:
            lvl_kw = dict(kw)
            lvl_kw.update(upwind_c=c_post, phi_init=result["phi"],
                          gamma_init=result["gamma"], n_picard_seed=0)
            result = solve_newton_lifting(mesh_cut, wc, m_inf=m_inf,
                                          alpha_deg=alpha_deg, workspace=ws,
                                          verbose=verbose, **lvl_kw)
            level_history.append((m_inf, result["n_newton"],
                                  result["residual_history"][-1]))
            if not result["converged"]:
                break

    result = dict(result)
    result["level_history"] = level_history
    result.pop("workspace", None)
    return result
