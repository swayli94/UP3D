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
import scipy.linalg as sla
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
from pyfp3d.solve.timing import (
    finalize,
    new_timings,
    snapshot,
    step_delta,
    sum_timings,
)


class _EliminatedKuttaRow:
    """K_tilde = -D^{-1} K_p as an operator supporting `@` on vectors AND
    dense matrices -- the pressure-estimator counterpart of the probe
    path's constant sparse K (whose D = -I makes K_tilde = K exactly).
    Used in exactly the three driver expressions the probe K is used in:
    the coupled matvec J_ff x + B (K x), the Woodbury S = I + K @ JB, and
    the step map dgamma = K @ dphi + F_tilde."""

    def __init__(self, Kp_free: sp.csr_matrix, D_lu):
        self._Kp = Kp_free
        self._lu = D_lu
        self.shape = Kp_free.shape

    def __matmul__(self, x):
        return -sla.lu_solve(self._lu, self._Kp @ x)


class NewtonWorkspace:
    """Per-case precomputation for the coupled Newton solve: operators,
    the free/Dirichlet split, the Kutta row K, and the per-Mach-level
    far-field basis (vals0_red, V_red). Reused across the Mach levels of
    solve_newton_transonic (set_mach re-derives only the beta-dependent
    far-field basis).

    `kutta_estimator` selects the Kutta closure residual F (P14):
      "probe" (default, BIT-IDENTICAL): F = tip_taper * kutta_targets - Gamma
        with the constant sparse row K and dF/dGamma = -I exactly.
      "pressure": F = sigma * station-mean(|q_u|^2 - |q_l|^2) on the
        wall-adjacent TE control volumes (constraints/te_pressure.py).
        sigma_j = 1/max(|D_jj|, 0.1 median) is frozen at the FIRST residual
        evaluation of this workspace (the seed state) -- once per workspace,
        NOT per Mach level (the transonic ramp reuses the workspace, and
        re-freezing would add a second ramp-history dependence on top of
        warm starts) -- so F stays in Gamma-like units for the shared merit
        |R|^2 + |F|^2 and the f_norm < tol_gamma check. sigma cancels in
        the eliminated step (kutta_blocks), so the Newton ITERATES are
        sigma-independent; only the merit weighting sees it.
        With a non-unit tip_taper (B31; the P13/G13.2 tip unload ported to
        the homogeneous row): scaling the pressure row alone is a no-op
        (P_j = 0 is homogeneous), so the row is BLENDED with a Gamma pin
        (the B8 row-blend form -- but welding an EXPLICIT Newton unknown,
        so unlike the LS span_blend there is no global re-leveling
        channel):

          F_j = taper_j * sigma_j * F_raw_j + (1 - taper_j) * s_j * Gamma_j

        taper_j = 1 is exactly the production row; taper_j -> 0 pins
        Gamma_j to 0. The pin slope s_j = sign(diag D0)_j (frozen WITH
        sigma at the first residual evaluation, recorded as
        kutta_weld_sign) is the measured orientation of the row's own
        Gamma-sensitivity: the naive unsigned weld -(1-t)*Gamma is only
        correct under the dF/dGamma ~ -I orientation, while the measured
        conforming meshes carry diag D > 0 (NACA 2.5D coarse: +80; M6
        wing-body conforming coarse: median +566, uniform sign, 0 flips),
        where it AMPLIFIES mid-taper loading (measured kappa = 1.08 at
        t = 0.7) and the blended Gamma block sigma(tD - (1-t)) crosses
        zero at t = 1/(1+D). With the row's own sign the blend unloads
        monotonically toward the pinned endpoint (at the fixed-phi
        linearization Gamma_b = t * Gamma*, the probe taper's Gamma_eff
        = taper * Gamma_Kutta semantics; the re-equilibrated R = 0
        manifold unloads HARDER -- the pressure row's Gamma-slope along
        the manifold is the Schur-complement slope, much smaller than
        the frozen direct slope, so mid-blend the pin dominates; measured
        kappa(0.7) ~ 0.14 on NACA 2.5D coarse) and the raw Gamma block
        stays ~ diag(D0)-conditioned at every t. In raw units (F/sigma)
        the weld slope is s_j/sigma_j = sign(dj)*max(|dj|, 0.1 median) =
        the floored frozen dF_raw/dGamma diagonal, so the blend is
        dimensionally consistent (the pin enters with the row's own
        slope, the B8 weld analog); because that ratio is physical, the
        untapered path's post-hoc sigma-independence of the iterates does
        NOT carry over -- sigma must be frozen, which it is.

    `external_rhs` (Track V V2 transpiration channel; roadmap
    track_v.md "V2 -- Transpiration channel through all three drivers"):
    an optional (n_cut,) external volume-RHS vector on the CUT mesh --
    the transpiration wall load of viscous/transpiration.py. It is
    subtracted inside eval_residual (R_red = T^T (R - b_ext)), so the
    converged state solves T^T R(phi, Gamma) = T^T b_ext. LAGGED: b_ext
    is frozen per solve and independent of (phi, Gamma), so the coupled
    Jacobian is UNCHANGED (the GV2.1(c) FD-exactness clause) -- tight
    coupling (d b_ext/d BL terms) is V5's augmentation, not this
    channel. None (default) is bit-identical to the pre-V2 workspace.
"""

    def __init__(self, mesh_cut, wc, alpha_deg: float = 0.0,
                 u_inf: float = 1.0, gamma_air: float = 1.4,
                 vortex_center=(0.25, 0.0),
                 farfield_spanwise_gamma: bool = False,
                 tip_taper: Optional[np.ndarray] = None,
                 kutta_estimator: str = "probe",
                 external_rhs: Optional[np.ndarray] = None):
        self.mesh_cut = mesh_cut
        self.wc = wc
        self.alpha_deg = float(alpha_deg)
        self.u_inf = float(u_inf)
        self.gamma_air = float(gamma_air)
        self.vortex_center = vortex_center
        self.spanwise = bool(farfield_spanwise_gamma)
        # P13/G13.2 tip-edge desingularization: per-station loading taper
        # F_j (constraints/wake.py::tip_taper_factors). The accepted
        # circulation becomes Gamma_eff = F_j * Gamma_Kutta, so BOTH the
        # residual F and its derivative K = dF/dphi scale by F_j. None (the
        # default) = all ones = bit-identical to the untapered path.
        self.tip_taper = (
            np.ones(wc.n_stations, dtype=np.float64) if tip_taper is None
            else np.asarray(tip_taper, dtype=np.float64).copy()
        )
        if self.tip_taper.shape != (wc.n_stations,):
            raise ValueError(
                f"tip_taper must be ({wc.n_stations},), got "
                f"{self.tip_taper.shape}")
        if kutta_estimator not in ("probe", "pressure"):
            raise ValueError(f"unknown kutta_estimator {kutta_estimator!r}")
        self.kutta_estimator = kutta_estimator
        # V2 transpiration channel (class docstring): lagged external RHS
        # on the cut mesh, subtracted in eval_residual. None = channel
        # absent = bit-identical to the pre-V2 path.
        if external_rhs is not None:
            external_rhs = np.asarray(external_rhs, dtype=np.float64)
            if external_rhs.shape != (len(mesh_cut.nodes),):
                raise ValueError(
                    f"external_rhs must be ({len(mesh_cut.nodes)},) on the "
                    f"cut mesh, got {external_rhs.shape}")
        self.external_rhs = external_rhs
        # B31: a non-unit taper activates the probe-path K row scaling
        # below resp. the pressure-path row blend (eval_residual /
        # kutta_blocks, class docstring). All-ones (incl. the None
        # default) keeps BOTH estimators bit-identical to untapered.
        self._taper_active = not np.all(self.tip_taper == 1.0)
        self.kutta_sigma = None        # frozen at the first eval_residual
        self.kutta_sigma_sign_flips = 0
        # B31 blend pin slope (class docstring): sign(diag D0)_j, frozen
        # alongside sigma; None until the first pressure-path evaluation.
        self.kutta_weld_sign = None

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

        # P14 pressure estimator: build the TE control volumes ONLY when
        # requested (the probe path's construction work is untouched). The
        # Dirichlet-contact assert inside guarantees dF/dGamma carries only
        # the slave-jump chain -- no V_red far-field term to add.
        self.cvs = None
        if kutta_estimator == "pressure":
            from pyfp3d.constraints.te_pressure import TEControlVolumes

            self.cvs = TEControlVolumes(mesh_cut, wc,
                                        dirichlet_nodes=dir_nodes)

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
        # F_j = taper_j * mean(probe jumps)_j - Gamma_j, so dF/dphi picks up
        # the SAME constant factor (taper is phi-independent). Scaling the
        # rows keeps the Jacobian exact; taper == 1 leaves K untouched.
        if self._taper_active:
            self.K = sp.diags(self.tip_taper) @ self.K

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
                      rho_floor: float, frozen=None):
        """The single state-reconstruction + residual path (module
        docstring). Returns (R_free, F, state) with state caching
        everything the Jacobian assembly at this point needs.

        `frozen` = (upstream, branch) from
        UpwindOperator.freeze_upwind_state switches the flux to the
        frozen-assignment sweep (the N5 Newton finish phase: within a
        fixed assignment the residual is smooth, so Newton converges
        quadratically where the live walk/branch churn limit-cycles on
        tie-degenerate meshes)."""
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
        if self.m_inf > 0.0 and frozen is not None:
            rho_t = self.upw.rho_tilde_frozen(
                q2l, rho, frozen[0], frozen[1], self.m_inf, upwind_c,
                m_crit, self.gamma_air, rho_floor).copy()
            n_floored = self.upw.n_floored
            nu_max = self.upw.nu_max
            n_nu_active = self.upw.n_supersonic
        elif self.m_inf > 0.0:
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

        R_vol = self.op.assemble_residual(phi_cut, rho_t)
        if self.external_rhs is not None:
            # V2 transpiration channel (class docstring): lagged, so the
            # Jacobian below is untouched; the converged state solves
            # T^T R = T^T b_ext.
            R_vol = R_vol - self.external_rhs
        R_red = self.con.T.T @ R_vol
        R_free = R_red[self.free]
        if self.kutta_estimator == "pressure":
            # P14: sigma-scaled pressure-equality residual; sigma frozen at
            # the first evaluation of this workspace (class docstring).
            F_raw = self.cvs.residual_stations(phi_cut)
            if self.kutta_sigma is None:
                # sigma is a MERIT WEIGHTING only (it cancels in the
                # eliminated step -- kutta_blocks), so a rough seed state
                # is good enough: |D_jj| sets the scale, the floor bounds
                # stations sitting near a transient sign change (a 5-step
                # Picard seed CAN carry flipped-sign stations near the tip
                # -- measured on M6 medium M0.5 -- without the estimator
                # being degenerate; the per-step exact D carries the true
                # signs). Flip count recorded for diagnostics.
                D0 = self.cvs.gamma_jacobian(phi_cut, mode="exact")
                dj = np.diag(D0)
                sign_ref = np.sign(np.median(dj))
                self.kutta_sigma_sign_flips = int(
                    np.count_nonzero(np.sign(dj) != sign_ref))
                adj = np.abs(dj)
                self.kutta_sigma = 1.0 / np.maximum(
                    adj, 0.1 * np.median(adj))
                # B31: the blend's pin slope is the SIGN of the row's own
                # frozen Gamma-sensitivity (class docstring); a measure-
                # zero dj == 0 takes the mesh-wide median sign.
                s = np.sign(dj)
                self.kutta_weld_sign = np.where(s == 0.0, sign_ref, s)
            if not self._taper_active:
                F = self.kutta_sigma * F_raw
            else:
                # B31 row blend (class docstring): F_j = taper_j *
                # sigma_j * F_raw_j + (1 - taper_j) * s_j * Gamma_j --
                # the production row at taper_j = 1, a unit-slope Gamma
                # pin (with the row's own orientation) at taper_j = 0.
                # Gamma enters EXPLICITLY here (an argument of
                # eval_residual), unlike the probe path.
                F = (self.kutta_sigma * (self.tip_taper * F_raw)
                     + (1.0 - self.tip_taper)
                     * (self.kutta_weld_sign * gamma))
        else:
            # Gamma_eff = taper * Gamma_Kutta (P13/G13.2; taper == 1 by
            # default, so this is the untapered Kutta residual bit-for-bit).
            F = self.tip_taper * kutta_targets(phi_cut, self.wc) - gamma
        state = {
            "phi_red": phi_red, "phi_cut": phi_cut, "grad": grad,
            "q2l": q2l, "lim": lim, "rho": rho, "rho_t": rho_t,
            "n_limited": n_limited, "n_floored": n_floored,
            "nu_max": nu_max, "n_nu_active": n_nu_active,
        }
        return R_free, F, state

    def assemble_coupled(self, state, upwind_c: float, m_crit: float,
                         rho_floor: float, frozen=None):
        """Exact Jacobian blocks at the state returned by eval_residual:
        (J_ff, B) with B = J_red[free, dir] @ V_red + H_J[free, :].
        Pass the SAME `frozen` given to eval_residual (consistency)."""
        if self.m_inf > 0.0 and frozen is not None:
            upstream = frozen[0]
            s_e, s_u = self.upw.rho_tilde_frozen_sensitivities(
                state["q2l"], state["rho"], frozen[0], frozen[1],
                self.m_inf, upwind_c, m_crit, self.gamma_air, rho_floor)
        elif self.m_inf > 0.0:
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

    def kutta_blocks(self, state, F):
        """(K_tilde, F_tilde) for the eliminated coupled step

            delta_gamma = K_tilde @ dphi + F_tilde,
            (J_ff + B K_tilde) dphi = -R - B F_tilde,

        exact for BOTH estimators (elimination of the second block row
        K dphi + D dgamma = -F with K = dF/dphi_free, D = dF/dGamma).

        Probe: D = -I exactly, so this returns (self.K, F) verbatim and
        the pre-P14 driver expressions execute unchanged (G14.4 by
        structure). Pressure: K_tilde = -D^{-1} K_p, F_tilde = -D^{-1}
        F_raw, with (K_p, D) the EXACT state-dependent rows rebuilt every
        Newton step (unlike the constant probe K); sigma cancels
        (-(sigma D)^{-1} (sigma F_raw) = -D^{-1} F_raw), so the Newton
        iterates are sigma-independent. D carries only the slave-jump
        chain: the CV construction asserted no control-volume dof is a
        far-field Dirichlet node, so restricting K_p to free columns
        drops nothing. B31 blend (class docstring): the code-level row is
        F_j = taper_j * sigma_j * F_raw_j + (1 - taper_j) * s_j * Gamma_j,
        so F / sigma recovers the RAW blend G_j = taper_j * F_raw_j +
        (1 - taper_j) * (s_j/sigma_j) * Gamma_j and the raw eliminated
        blocks are dG/dphi_free = diag(taper) K_p and dG/dGamma =
        diag(taper) D + diag((1 - taper) * s/sigma) -- the pin's
        derivative is the frozen-slope diagonal on the Gamma block
        (Gamma_j is an explicit Newton unknown; the slave-jump chain sits
        entirely in the pressure part)."""
        if self.kutta_estimator == "probe":
            return self.K, F
        Kp_cut, D = self.cvs.newton_rows(state["phi_cut"])
        Kp_free = (Kp_cut @ self.con.T).tocsr()[:, self.free]
        if not self._taper_active:
            D_lu = sla.lu_factor(D)
            F_tilde = -sla.lu_solve(D_lu, F / self.kutta_sigma)
            return _EliminatedKuttaRow(Kp_free, D_lu), F_tilde
        t = self.tip_taper
        Kt_free = (sp.diags(t) @ Kp_free).tocsr()
        Db = t[:, None] * D
        Db[np.diag_indices_from(Db)] += ((1.0 - t) * self.kutta_weld_sign
                                         / self.kutta_sigma)
        D_lu = sla.lu_factor(Db)
        F_tilde = -sla.lu_solve(D_lu, F / self.kutta_sigma)
        return _EliminatedKuttaRow(Kt_free, D_lu), F_tilde


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
    direct_refactor_every: int = 1,
    direct_reuse_rtol: float = 1e-8,
    gmres_restart: int = 60,
    gmres_maxiter: int = 10,
    line_search: bool = True,
    ptc_dtau: Optional[float] = None,
    rtol_seed: float = 1e-7,
    tol_residual_loose: Optional[float] = None,
    tol_residual_rel: Optional[float] = None,
    accept_on_stall: bool = False,
    freeze_tol: Optional[float] = None,
    freeze_refresh_max: int = 2,
    freeze_max_reverts: int = 3,
    workspace: Optional[NewtonWorkspace] = None,
    tip_taper: Optional[np.ndarray] = None,
    kutta_estimator: str = "probe",
    external_rhs: Optional[np.ndarray] = None,
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

    `tol_residual_loose` / `tol_residual_rel` / `accept_on_stall`
    (G10.2, all default off = bit-identical): LOOSE acceptance criteria
    for levels whose only job is seeding the next Mach-continuation
    level. With the same 0-limited/0-floored and ||F|| < tol_gamma
    guards, the level is also accepted once ||R||_inf < tol_residual_
    loose, or the residual has dropped by 1/tol_residual_rel from the
    level-entry value, or the live-stall detector fires (the same
    plateau test that triggers the freeze) -- i.e. intermediate-level
    stall becomes "advance" instead of "freeze and polish". All three
    require AT LEAST ONE Newton step at this level (measured A/B trap:
    warm-started fold-zone levels ENTER below any absolute threshold,
    and zero-step acceptance degenerates the continuation into a level
    skip -- the M0.7875 medium ramp then has no tracking seed and the
    final level diverges) and never apply inside a frozen phase. The
    result records which criterion fired as `accept_reason`
    ("tol"/"loose_tol"/"rel_drop"/"stall").

    `direct_refactor_every` (the N6 3D-cost fix): with precond="direct",
    refactor the LU only every k-th Newton step and drive the steps in
    between with GMRES on the FRESH coupled operator preconditioned by
    the stale LU, converged to `direct_reuse_rtol` (tight -- orders below
    the eta that stalls the shock-position soft mode, affordable only
    because the stale LU is a near-exact preconditioner). A reuse step
    whose GMRES fails falls back to refactor + exact Woodbury in the same
    iteration, so robustness never degrades below the every-step-direct
    path. Default 1 = factor every step, bit-identical to the N5 G8.1
    behavior. Motivation: on true-3D meshes the LU fill makes each
    factorization ~100x more expensive than on the thin 2.5D family
    (measured 18.6 s vs ~0.2 s at the same ~6e4 dofs -- the M6 medium
    all-levels run was 1606 s with 97% in splu).

    `freeze_tol` (the N5 churn fix): once ||R||_inf drops below it (with
    0 limited/floored), the upwind SELECTION AND BRANCH are frozen
    (UpwindOperator.freeze_upwind_state) and Newton finishes on the
    now-smooth frozen system. Rationale: on tie-degenerate meshes (the
    2.5D prism-split family) ~1e3 elements sit in the near-tie band of
    max(nu_e, nu_u) and the live sweep's branch/selection churn
    limit-cycles the residual at ~1e-8 (measured on the medium G4.1
    ramp); within a fixed assignment the residual is smooth and the
    terminal quadratic rate is restored. This is Lopez's frozen-selection
    architecture made persistent. On reaching tol the LIVE residual is
    re-evaluated and reported as `residual_unfrozen` with
    `n_assignment_stale` flip counts; if it misses tol the freeze is
    refreshed from the fresh state (up to `freeze_refresh_max` times --
    active-set iteration) and the best result reported. None = off
    (bit-identical to the pre-freeze driver).

    `freeze_max_reverts` (A3, backported from newton_ls.py's B15 net): a
    freeze that keeps diverging is worse than no freeze, so after this many
    reverts the freeze is DISARMED for the rest of the level and the live
    path finishes it. The freeze may only ever help; it can never cost
    convergence. Inert on every committed conforming run (no level has ever
    reverted more than twice).

    `kutta_estimator` (P14): "probe" (default, bit-identical) or
    "pressure" -- the wall-adjacent-CV pressure-equality Kutta residual
    (NewtonWorkspace docstring for the sigma scaling; kutta_blocks for the
    exact eliminated step). Forwarded through solve_newton_transonic via
    newton_kw; a reused workspace must have been built with the same flag.
    With a non-unit `tip_taper` the pressure row runs the B31 Gamma-pin
    blend (NewtonWorkspace docstring); the probe path's per-station
    scaling is unchanged.

    `external_rhs` (Track V V2 transpiration channel): optional lagged
    external volume RHS on the cut mesh (NewtonWorkspace docstring).
    None (default) = channel absent = bit-identical. Forwarded through
    solve_newton_transonic via newton_kw like the other workspace flags;
    a reused workspace must carry the vector itself.

    Returns the Picard-compatible result keys (phi on the cut mesh, gamma,
    converged, residual_history, mach2_max, nu_max, n_nu_active,
    n_limited, n_floored) plus F_history, newton_orders (per-step observed
    order p_k), eta_history, n_gmres_total, n_newton, jacobian_nnz,
    n_term3_active, timings (per-stage wall-clock seconds), and
    kutta_estimator/kutta_sigma (the frozen pressure scaling, None on the
    probe path) and kutta_weld_sign (the frozen B31 blend pin slope, None
    on the probe path and recorded under the pressure estimator).
    """
    if m_inf > 0.0 and upwind_c <= 0.0:
        raise ValueError(
            "the Newton driver runs with the walk upwind machinery active "
            "for m_inf > 0 (bit-inert subcritically); pass upwind_c > 0")

    ws = workspace
    if ws is None:
        ws = NewtonWorkspace(mesh_cut, wc, alpha_deg, u_inf, gamma_air,
                             vortex_center, farfield_spanwise_gamma,
                             tip_taper=tip_taper,
                             kutta_estimator=kutta_estimator,
                             external_rhs=external_rhs)
    elif ws.kutta_estimator != kutta_estimator:
        raise ValueError(
            f"workspace was built with kutta_estimator="
            f"{ws.kutta_estimator!r} but {kutta_estimator!r} was requested "
            "-- pass the flag consistently (transonic ramps forward it via "
            "newton_kw)")
    elif external_rhs is not None and ws.external_rhs is not external_rhs:
        raise ValueError(
            "the workspace carries its own external_rhs (or none): build "
            "the NewtonWorkspace with the vector and pass workspace=... "
            "instead of forwarding external_rhs alongside it")
    ws.set_mach(m_inf)

    # Canonical Track-A schema (solve/timing.py) PLUS the three legacy keys
    # `jacobian`/`amg_setup`/`gmres` this driver has always reported --
    # cases/demo/p8_newton reads them to draw a committed (gated, expensive)
    # figure, so they stay, bit for bit. The canonical aliases are the ones
    # to read: `assembly` = jacobian, `precond` = amg_setup (which has always
    # included the splu factorization, despite the name), `linsolve` = gmres.
    t_wall0 = time.perf_counter()
    timings = new_timings()
    timings.update({"jacobian": 0.0, "amg_setup": 0.0, "gmres": 0.0})
    step_records = []
    n_refactor = 0

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
    gamma_history = []
    clamp_history = []             # (n_limited, n_floored) per iteration
    eta_history = []
    newton_orders = []
    n_gmres_total = 0
    n_gmres_stalled = 0
    converged = False
    accept_reason = None
    M_pre = None
    eta_prev = None
    lu_direct = None
    lu_age = 0                     # Newton steps driven by lu_direct

    t0 = time.perf_counter()
    R_free, F, state = ws.eval_residual(phi_free, gamma, upwind_c, m_crit,
                                        m_cap, rho_floor)
    timings["residual"] += time.perf_counter() - t0
    merit = float(R_free @ R_free + F @ F)

    frozen = None
    n_freeze_refresh = 0
    residual_unfrozen = None
    n_assignment_stale = 0
    assignment_cycle = False
    prev_live_residual = None
    freeze_point = None            # (phi_free, gamma, r_norm) at freeze
    freeze_cooldown = 0
    n_freeze_reverts = 0
    r_level_best = np.inf
    freeze_armed = freeze_tol is not None

    # A1: each record carries the state ON ENTRY to a Newton step plus the
    # cost OF that step, so it is closed once the step has been taken (the
    # freeze revert/refresh paths `continue` without solving -- their records
    # legitimately show n_lin_solves = 0).
    prev_snap = snapshot(timings)
    prev_counts = [0, 0]           # n_gmres_total, n_refactor at last close

    def _close_step(rec):
        nonlocal prev_snap
        rec.update(step_delta(timings, prev_snap))
        rec["n_lin_iters"] = n_gmres_total - prev_counts[0]
        rec["n_refactor"] = n_refactor - prev_counts[1]
        rec["wall_cum_s"] = time.perf_counter() - t_wall0
        prev_snap = snapshot(timings)
        prev_counts[0] = n_gmres_total
        prev_counts[1] = n_refactor

    for it in range(n_newton_max):
        r_norm = float(np.max(np.abs(R_free)))
        f_norm = float(np.max(np.abs(F))) if ws.n_st > 0 else 0.0
        residual_history.append(r_norm)
        F_history.append(f_norm)
        gamma_history.append(gamma.copy())
        clamp_history.append((state["n_limited"], state["n_floored"]))
        if step_records:
            _close_step(step_records[-1])
        rec = {
            "i": it,
            "residual": r_norm,
            "F": f_norm,
            "gamma_mean": float(np.mean(gamma)) if ws.n_st > 0 else 0.0,
            "gamma_root": float(gamma[0]) if ws.n_st > 0 else 0.0,
            "n_limited": int(state["n_limited"]),
            "n_floored": int(state["n_floored"]),
            "frozen": frozen is not None,
            "n_lin_solves": 0,
            "eta": None,
            "lam": None,
        }
        step_records.append(rec)
        if len(residual_history) >= 3:
            r2, r1, r0 = (residual_history[-1], residual_history[-2],
                          residual_history[-3])
            if r2 > 0.0 and r1 > 0.0 and r0 > 0.0 and r1 != r0:
                newton_orders.append(
                    float(np.log(r2 / r1) / np.log(r1 / r0)))
        if verbose:
            print(f"  newton {it:2d}: |R|={r_norm:.3e} |F|={f_norm:.3e} "
                  f"lim={state['n_limited']} flr={state['n_floored']}"
                  + (" [frozen]" if frozen is not None else ""))
        r_level_best = min(r_level_best, r_norm)

        # frozen-phase safety net: the frozen sweep is floor-free by
        # design, so a bad assignment can send its Newton path into
        # clamp land (measured: refresh with 928 stale assignments at
        # M0.775 medium -> M_max 3.0, 3800 limited). Revert to the
        # freeze point, unfreeze, and let live Newton (+ level fail-fast
        # below) take over after a cooldown.
        if frozen is not None and freeze_point is not None:
            frozen_diverged = (
                r_norm > 10.0 * freeze_point[2]
                or state["n_limited"] > 0 or state["n_floored"] > 0)
            if frozen_diverged:
                phi_free = freeze_point[0].copy()
                gamma = freeze_point[1].copy()
                frozen = None
                freeze_point = None
                freeze_cooldown = 5
                n_freeze_reverts += 1
                # A3 (C2, backported from newton_ls.py): r_level_best is
                # scoped to the CURRENT SELECTION EPOCH. Frozen-phase and
                # live-phase residuals are not comparable -- the frozen
                # system descends to ~1e-11, and reverting/refreshing
                # legitimately returns the residual to the live scale. Carry
                # the frozen best across that boundary and the fail-fast
                # below reads a spurious 1e8x "blow-up" and kills a healthy
                # freeze cycle (measured on the LS path, M6 medium M0.65).
                r_level_best = np.inf
                # A3 (C3): a freeze that keeps diverging is worse than no
                # freeze -- disarm permanently and let the live path finish.
                if n_freeze_reverts >= freeze_max_reverts:
                    freeze_armed = False
                if verbose:
                    print("  [freeze reverted: frozen path diverged]")
                t0 = time.perf_counter()
                R_free, F, state = ws.eval_residual(
                    phi_free, gamma, upwind_c, m_crit, m_cap, rho_floor)
                timings["residual"] += time.perf_counter() - t0
                merit = float(R_free @ R_free + F @ F)
                continue
        if freeze_cooldown > 0:
            freeze_cooldown -= 1

        # level fail-fast: a warm start this bad will not recover within
        # the budget -- fail early so the continuation halves the Mach
        # step instead of reporting a wandered state
        if r_norm > 100.0 * r_level_best and r_norm > 1e-3:
            if verbose:
                print("  [level fail-fast: residual 100x above best]")
            break
        # Freeze trigger: either the residual is already small
        # (< freeze_tol), or the LIVE Newton has stalled -- no 4x progress
        # over the last 6 steps below the settledness cap. The stall
        # trigger adapts to the level's churn floor (measured 1e-6 at
        # M0.70 but 1.4e-5 at M0.80 on the medium mesh -- no absolute
        # threshold fits all levels), while the cap keeps a far-from-
        # solution transient from locking in a garbage assignment
        # (measured: freezing at 1e-4 mid-churn sent M0.75 to 1.2e-2).
        # plateau, not growth: freezing a diverging transient locks a
        # garbage assignment (measured M0.80 medium -> M_max 3.0 blowup)
        live_stalled = (
            len(residual_history) >= 7
            and r_norm < 1e-3
            and r_norm > 0.25 * min(residual_history[-7:-1])
            and r_norm < 2.0 * min(residual_history[-7:-1]))
        if (frozen is None and freeze_armed and freeze_tol is not None
                and freeze_cooldown == 0
                and m_inf > 0.0 and (r_norm < freeze_tol or live_stalled)
                and state["n_limited"] == 0 and state["n_floored"] == 0):
            frozen = ws.upw.freeze_upwind_state(
                state["grad"], state["q2l"], state["rho"], m_inf,
                upwind_c, m_crit, ws.gamma_air, rho_floor)
            freeze_point = (phi_free.copy(), gamma.copy(), r_norm)
            r_level_best = np.inf      # new selection epoch (A3/C2)
            # the frozen flux equals the live flux AT the freeze state, so
            # (R_free, F, state) stay valid for this iteration
        accept = None
        if (f_norm < tol_gamma and state["n_limited"] == 0
                and state["n_floored"] == 0):
            if r_norm < tol_residual:
                accept = "tol"
            elif frozen is None and it > 0:
                if (tol_residual_loose is not None
                        and r_norm < tol_residual_loose):
                    accept = "loose_tol"
                elif (tol_residual_rel is not None
                      and r_norm <= tol_residual_rel
                      * residual_history[0]):
                    accept = "rel_drop"
                elif accept_on_stall and live_stalled:
                    accept = "stall"
        if accept is not None:
            if frozen is None:
                converged = True
                accept_reason = accept
                break
            # honesty check: the frozen system converged -- how far is the
            # LIVE (re-selected, re-branched) residual?
            R_live, F_live, state_live = ws.eval_residual(
                phi_free, gamma, upwind_c, m_crit, m_cap, rho_floor)
            residual_unfrozen = float(np.max(np.abs(R_live)))
            up_new, br_new = ws.upw.freeze_upwind_state(
                state_live["grad"], state_live["q2l"], state_live["rho"],
                m_inf, upwind_c, m_crit, ws.gamma_air, rho_floor)
            n_assignment_stale = int(
                np.count_nonzero(up_new != frozen[0])
                + np.count_nonzero((br_new != frozen[1])
                                   & (up_new == frozen[0])))
            if verbose:
                print(f"  frozen system converged; live |R|="
                      f"{residual_unfrozen:.3e}, stale assignments="
                      f"{n_assignment_stale}")
            if residual_unfrozen < tol_residual:
                converged = True
                break
            if (prev_live_residual is not None
                    and residual_unfrozen > 0.5 * prev_live_residual):
                # active-set two-cycle: the remaining stale elements flip
                # assignment at each other's converged state, so the live
                # residual has hit the discretization's intrinsic
                # assignment-discontinuity floor -- accept the frozen
                # solution and report the floor honestly
                assignment_cycle = True
                converged = True
                break
            if n_freeze_refresh >= freeze_refresh_max:
                converged = True
                break
            # refresh the freeze from the fresh state and keep going
            prev_live_residual = residual_unfrozen
            n_freeze_refresh += 1
            frozen = (up_new, br_new)
            freeze_point = (phi_free.copy(), gamma.copy(),
                            residual_unfrozen)
            r_level_best = np.inf      # new selection epoch (A3/C2)
            R_free, F, state = R_live, F_live, state_live
            merit = float(R_free @ R_free + F @ F)
            continue

        t0 = time.perf_counter()
        J_ff, B = ws.assemble_coupled(state, upwind_c, m_crit, rho_floor,
                                      frozen=frozen)
        if ptc_dtau is not None:
            J_ff = (J_ff + sp.diags(ws.m_lumped_free / ptc_dtau)).tocsr()
        timings["assembly"] += time.perf_counter() - t0

        # (K, F_elim) = the eliminated Kutta blocks. Probe: literally
        # (ws.K, F) -- the historical expressions below run unchanged.
        # Pressure: state-dependent, rebuilt EVERY Newton step (the probe
        # K's loop-invariance does not carry over -- P14 trap A).
        K, F_elim = ws.kutta_blocks(state, F)
        rhs = -R_free - B @ F_elim
        if precond == "direct":
            # exact Newton step: sparse LU of J_ff + Woodbury for the
            # rank-n_st coupling B K. The transonic Jacobian carries a
            # near-singular shock-position mode that grows under mesh
            # refinement -- an eta-accurate Krylov step leaves an O(1)
            # error in that mode and Newton stalls (measured on the
            # medium G4.1 ramp: frozen-system residual flat at 3.7e-6
            # with GMRES converging to eta); the exact step restores the
            # quadratic rate. Cost: one splu + (1 + n_st) back-solves per
            # step -- seconds at 6e4 dofs, the Lopez-scale 2D/medium-3D
            # regime. GMRES+AMG remains the large-mesh path.
            dphi = None
            if lu_direct is not None and lu_age < direct_refactor_every:
                # lagged-LU reuse step (N6): fresh coupled operator,
                # stale-LU preconditioner, tight rtol; fall back to a
                # refactor + exact step if GMRES does not converge
                A_op = spla.LinearOperator(
                    (ws.n_free, ws.n_free),
                    matvec=lambda x: J_ff @ x + B @ (K @ x),
                )
                M_lag = spla.LinearOperator(
                    (ws.n_free, ws.n_free), matvec=lu_direct.solve)
                t0 = time.perf_counter()
                dphi, n_it, gmres_info = solve_gmres(
                    A_op, rhs, M=M_lag, rtol=direct_reuse_rtol,
                    restart=gmres_restart, maxiter=2, on_fail="return")
                timings["linsolve"] += time.perf_counter() - t0
                n_gmres_total += n_it
                if gmres_info != 0:
                    dphi = None            # stale LU exhausted: refactor
            if dphi is None:
                t0 = time.perf_counter()
                lu_direct = spla.splu(J_ff.tocsc())
                lu_age = 0
                timings["precond"] += time.perf_counter() - t0
                n_refactor += 1
                t0 = time.perf_counter()
                z = lu_direct.solve(rhs)
                JB = lu_direct.solve(B.toarray())
                S = np.eye(ws.n_st) + K @ JB
                dphi = z - JB @ np.linalg.solve(S, K @ z)
                timings["linsolve"] += time.perf_counter() - t0
            lu_age += 1
            eta_history.append(0.0)
            rec["eta"] = 0.0
            rec["n_lin_solves"] = 1
        else:
            t0 = time.perf_counter()
            if precond == "amg":
                if M_pre is None or it % amg_rebuild_every == 0:
                    A_pic = ws.op.assemble_matrix(state["rho_t"])
                    A_red = (ws.con.T.T @ (A_pic @ ws.con.T)).tocsr()
                    A_ff = A_red[ws.free][:, ws.free].tocsr()
                    if ptc_dtau is not None:
                        A_ff = (A_ff + sp.diags(
                            ws.m_lumped_free / ptc_dtau)).tocsr()
                    _, M_pre = build_amg_preconditioner(A_ff)
            elif precond == "ilu":
                M_pre = build_ilu_preconditioner(J_ff)
            else:
                raise ValueError(f"unknown precond {precond!r}")
            timings["precond"] += time.perf_counter() - t0

            A_op = spla.LinearOperator(
                (ws.n_free, ws.n_free),
                matvec=lambda x: J_ff @ x + B @ (K @ x),
            )
            eta = _ew_forcing(r_norm,
                              residual_history[-2] if it > 0 else None,
                              eta_prev, ew_eta0, ew_gamma, ew_eta_max)
            eta_prev = eta
            eta_history.append(eta)
            rec["eta"] = eta
            rec["n_lin_solves"] = 1
            t0 = time.perf_counter()
            dphi, n_it, gmres_info = solve_gmres(
                A_op, rhs, M=M_pre, rtol=eta, restart=gmres_restart,
                maxiter=gmres_maxiter, on_fail="return")
            timings["linsolve"] += time.perf_counter() - t0
            n_gmres_total += n_it
            if gmres_info != 0:
                # take the best iterate as an INEXACT Newton step (still
                # a useful direction; the line search guards it) instead
                # of burning stagnating Krylov cycles
                n_gmres_stalled += 1
        dgamma = K @ dphi + F_elim

        # safety-only backtracking: accept the full step unless the merit
        # blows up or goes non-finite (Lopez runs plain Newton; the load
        # stepping is the real globalization). If NO tried step improves,
        # take the LEAST-BAD one -- unconditionally accepting the
        # smallest lam allowed sustained residual growth (measured on the
        # medium ramp's churn phases).
        lam = 1.0
        best = None
        accepted = False
        for _ in range(5):
            t0 = time.perf_counter()
            R_try, F_try, state_try = ws.eval_residual(
                phi_free + lam * dphi, gamma + lam * dgamma,
                upwind_c, m_crit, m_cap, rho_floor, frozen=frozen)
            timings["residual"] += time.perf_counter() - t0
            merit_try = float(R_try @ R_try + F_try @ F_try)
            if np.isfinite(merit_try) and (not line_search
                                           or merit_try <= merit):
                accepted = True
                break
            if np.isfinite(merit_try) and (best is None
                                           or merit_try < best[0]):
                best = (merit_try, lam, R_try, F_try, state_try)
            lam *= 0.5
        if not accepted:
            if best is None:                 # pragma: no cover - safety net
                raise RuntimeError(
                    f"Newton line search failed at iteration {it} "
                    f"(merit {merit:.3e} -> {merit_try:.3e})")
            merit_try, lam, R_try, F_try, state_try = best
        rec["lam"] = float(lam)
        phi_free = phi_free + lam * dphi
        gamma = gamma + lam * dgamma
        R_free, F, state = R_try, F_try, state_try
        merit = merit_try

    if step_records:
        _close_step(step_records[-1])
    if converged and accept_reason is None:
        accept_reason = "tol"          # frozen-phase acceptance paths
    q2n = state["q2l"]
    mach2_max = float(np.max(mach_squared_field(q2n, m_inf, gamma_air)))
    # Legacy aliases over the canonical buckets, so cases/demo/p8_newton keeps
    # reading the keys it has always read (same seconds, same figure -- that
    # demo is gated and expensive, it must not need a re-run). `kutta` stays
    # 0.0: the Kutta residual F is evaluated inside ws.eval_residual, so its
    # cost already sits in `residual` and cannot be split out without double
    # counting.
    timings["jacobian"] = timings["assembly"]
    timings["amg_setup"] = timings["precond"]
    timings["gmres"] = timings["linsolve"]
    timings["kutta"] = 0.0
    finalize(timings, time.perf_counter() - t_wall0)
    return {
        "phi": state["phi_cut"],
        "gamma": gamma,
        "converged": converged,
        "accept_reason": accept_reason,
        "n_newton": len(residual_history) - (1 if converged else 0),
        "residual_history": residual_history,
        "F_history": F_history,
        "gamma_history": gamma_history,
        "clamp_history": clamp_history,
        "newton_orders": newton_orders,
        "eta_history": eta_history,
        "n_gmres_total": n_gmres_total,
        "n_gmres_stalled": n_gmres_stalled,
        "froze": frozen is not None,
        "n_freeze_refresh": n_freeze_refresh,
        "residual_unfrozen": residual_unfrozen,
        "n_assignment_stale": n_assignment_stale,
        "assignment_cycle": assignment_cycle,
        "n_freeze_reverts": n_freeze_reverts,
        "mach2_max": mach2_max,
        "nu_max": state["nu_max"],
        "n_nu_active": state["n_nu_active"],
        "n_limited": state["n_limited"],
        "n_floored": state["n_floored"],
        "jacobian_nnz": getattr(ws.op, "newton_nnz", None),
        "n_term3_active": getattr(ws.op, "n_term3_active", None),
        "n_refactor": n_refactor,
        "step_records": step_records,
        "timings": timings,
        "kutta_estimator": ws.kutta_estimator,
        "kutta_sigma": ws.kutta_sigma,
        "kutta_sigma_sign_flips": ws.kutta_sigma_sign_flips,
        "kutta_weld_sign": ws.kutta_weld_sign,
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
    freeze_tol: Optional[float] = 1e-6,
    intermediate_tol: Optional[float] = None,
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

    `intermediate_tol` (G10.2, opt-in; None = current behaviour
    BIT-IDENTICAL): level-adaptive tolerance for the INTERMEDIATE ramp
    levels, whose only role is to warm-start the next level. When set,
    the ORIGINAL-SCHEDULE levels except the final Mach level run with
    loose acceptance (tol_residual_loose=intermediate_tol after >= 1
    Newton step, a 1e3 relative-drop acceptance, stall-accept instead
    of freeze-and-polish; freeze_tol=None) -- the 0-limited/0-floored
    and ||F|| guards stay. The final level (and any upwind_c_post
    stages) keeps the full tol-1e-10 + freeze/honesty machinery, so
    the CONVERGED solution semantics are unchanged; and every level
    INSERTED by dm-halving runs STRICT too -- the halving cascade is
    the robustness fallback, and the A/B measured that loose retry
    levels (1 step each) cannot repair a fold-zone seed: the loose
    ramp reaches the final level with an untracked state, and only a
    strictly re-converged inserted level restores the default path's
    seed quality (the P8 "warm-start only from CONVERGED levels"
    trap, measured on the G8.1 NACA-medium A/B).

    Returns the final-level result dict plus level_history (per-level
    (m, n_newton, |R|_final)) and level_results (per-level dicts with
    the full residual/F/Gamma histories, wall_s and the level's timings
    -- the returned top-level `timings` covers only the FINAL level).
    """
    from pyfp3d.solve.continuation import mach_schedule

    if m_inf < m_start:
        raise ValueError(
            f"upward ramp only: m_inf {m_inf} < m_start {m_start} "
            "(design.md Sec 12 risk 2 -- never continue downward)")
    kw = dict(newton_kw or {})
    kw.setdefault("upwind_c", upwind_c)
    kw.setdefault("m_crit", m_crit)
    kw.setdefault("freeze_tol", freeze_tol)

    ws = None
    result = None
    last_good = None
    level_history = []
    level_results = []
    levels = mach_schedule(m_inf, m_start, dm)
    # G10.2 loose-level flags: original-schedule intermediates only; the
    # final level is always levels[-1] (dm-halving inserts BELOW the
    # failed level, never past m_inf) and inserted retry levels run
    # strict (docstring)
    loose = [True] * (len(levels) - 1) + [False]
    i = 0
    while i < len(levels):
        m = levels[i]
        lvl_kw = dict(kw)
        if intermediate_tol is not None and loose[i]:
            lvl_kw.update(tol_residual_loose=intermediate_tol,
                          tol_residual_rel=1e-3, accept_on_stall=True,
                          freeze_tol=None)
        if last_good is not None:
            # warm-start from the last CONVERGED level (a failed level's
            # state is not a valid continuation seed)
            lvl_kw.update(phi_init=last_good["phi"],
                          gamma_init=last_good["gamma"], n_picard_seed=0)
        if verbose:
            print(f"newton continuation level M = {m:.4f}")
        t_lvl = time.perf_counter()
        result = solve_newton_lifting(mesh_cut, wc, m_inf=m,
                                      alpha_deg=alpha_deg, workspace=ws,
                                      verbose=verbose, **lvl_kw)
        wall_lvl = time.perf_counter() - t_lvl
        ws = result["workspace"]
        level_history.append((m, result["n_newton"],
                              result["residual_history"][-1]))
        level_results.append({
            "m": m, "upwind_c": lvl_kw["upwind_c"],
            "converged": result["converged"],
            "accept_reason": result["accept_reason"],
            "n_newton": result["n_newton"],
            "residual_history": result["residual_history"],
            "F_history": result["F_history"],
            "gamma_history": result["gamma_history"],
            "step_records": result["step_records"],   # A1
            "n_lin_iters": result["n_gmres_total"],
            "n_refactor": result["n_refactor"],
            "wall_s": wall_lvl,
            "timings": result["timings"],
        })
        if not result["converged"]:
            m_prev = m_start if i == 0 else levels[i - 1]
            dm_new = 0.5 * (m - m_prev)
            if i == 0 or dm_new < dm_min:
                break
            levels.insert(i, m_prev + dm_new)
            loose.insert(i, False)      # retry levels run strict (G10.2)
            continue
        last_good = result
        i += 1

    if result["converged"] and upwind_c_post:
        for c_post in upwind_c_post:
            lvl_kw = dict(kw)
            lvl_kw.update(upwind_c=c_post, phi_init=result["phi"],
                          gamma_init=result["gamma"], n_picard_seed=0)
            t_lvl = time.perf_counter()
            result = solve_newton_lifting(mesh_cut, wc, m_inf=m_inf,
                                          alpha_deg=alpha_deg, workspace=ws,
                                          verbose=verbose, **lvl_kw)
            wall_lvl = time.perf_counter() - t_lvl
            level_history.append((m_inf, result["n_newton"],
                                  result["residual_history"][-1]))
            level_results.append({
                "m": m_inf, "upwind_c": lvl_kw["upwind_c"],
                "converged": result["converged"],
                "accept_reason": result["accept_reason"],
                "n_newton": result["n_newton"],
                "residual_history": result["residual_history"],
                "F_history": result["F_history"],
                "gamma_history": result["gamma_history"],
                "step_records": result["step_records"],   # A1
                "n_lin_iters": result["n_gmres_total"],
                "n_refactor": result["n_refactor"],
                "wall_s": wall_lvl,
                "timings": result["timings"],
            })
            if not result["converged"]:
                break

    result = dict(result)
    result["level_history"] = level_history
    result["level_results"] = level_results
    # A1: ramp total across all levels; the top-level `timings` (inherited
    # from dict(result)) is the FINAL level only -- the documented footgun.
    result["timings_total"] = sum_timings(
        [lr["timings"] for lr in level_results])
    result.pop("workspace", None)
    return result
