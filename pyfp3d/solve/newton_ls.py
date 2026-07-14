"""
Level-set (Track B) Newton solver -- the post-B6 re-derivation
(design_track_b.md section 5.5). Closes the FP-fold cases (e.g. medium
M0.7875) that the B6 Picard reaches only as a bounded stall: a Picard method,
conforming or level-set, does not reach the isolated fold solution (P4-erratum
/ P8 record), so the medium quantitative gate needs Newton.

Structure (simpler than the P8 conforming Newton, section 5.5): NO Gamma DOF
(the implicit Kutta makes the TE jump a solution mode -- no delta-Gamma
elimination, no Woodbury, no far-field vortex column when the far field is the
B6 neumann outlet); the wake-LS rows are LINEAR in phi (constant Jacobian
block); Terms 1-3 (kernels/jacobian.py physics) are reused per side through the
DOF indirection.

The residual is the multivalued fixed-point residual R = A(phi) phi - b, where
A is `MultivaluedOperator.assemble_matrix` at the current phi -- exact for all
rows including the TE Kutta row (its factorization makes row.x = |q_u|^2 -
|q_l|^2 exactly). The Jacobian adds, on top of that Picard matrix:
  - Terms 2+3 on the mass rows (per side, P7 frozen-selection sensitivities);
  - the EXACT quadratic derivative of the TE Kutta row (replacing the
    frozen-mean linearized row -- the frozen row is exact as a value but not
    as a slope, which is what quadratic convergence needs).
The wake-LS rows need no correction (already exact-linear).

Nonsymmetric (wake-LS aux rows + supersonic upwind) -> sparse-direct LU
(scipy spsolve) by default, consistent with B2/B3's choice. B11 (2026-07-14)
adds the `precond` kwarg (None = the bit-identical direct default; "ilu"/"amg"
run preconditioned GMRES on the fused Jacobian): the iterative escape for when
the LU dominates at 3D sizes (P8/N6 measured true-3D LU fill ~100x the 2.5D
cost). **Use "ilu"** -- it factors the real fused Jacobian and converges; "amg"
(the SPD surrogate `picard_ls._amg_surrogate_preconditioner`) STALLS on the
wake_ls-closure operator (measured; convection-like aux rows), so it is not the
Newton escape (design_track_b.md §5.3). NOTE the lagged-LU direct path
(`direct_refactor_every`, solve/newton.py) is a SEPARATE follow-up (roadmap
"LS Newton on M6 = DEFERRED"): B11 ships the iterative escape, not the
direct-reuse one.
"""

from typing import Dict, Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from pyfp3d.constraints.dirichlet import freestream_phi, vortex_phi_2d
from pyfp3d.kernels.cut_assembly import (
    newton_terms23_side_coo,
    te_kutta_jacobian_coo,
    te_kutta_residual,
    te_weld_coo,
)
from pyfp3d.solve.linear import build_ilu_preconditioner, solve_gmres
from pyfp3d.solve.picard_ls import (
    _amg_surrogate_preconditioner,
    _farfield_main,
    _farfield_split,
    _neumann_outlet_rhs,
    solve_multivalued_lifting,
)


def _seed_from_picard(mvop, mesh, m_inf, alpha_deg, u_inf, gamma_air,
                      upwind_c, m_crit, farfield, phi_init, gamma_init,
                      n_seed, damping_theta, omega_rho, precond=None):
    """Warm start: a short B6 Picard-LS solve to land near the solution before
    Newton (the transonic ramp / a lower-Mach state is the caller's job)."""
    r = solve_multivalued_lifting(
        mvop, mesh, m_inf, alpha_deg=alpha_deg, u_inf=u_inf,
        gamma_air=gamma_air, upwind_c=upwind_c, m_crit=m_crit,
        farfield=farfield, damping_theta=damping_theta, omega_rho=omega_rho,
        n_outer_max=n_seed, gamma_init=gamma_init, phi_init=phi_init,
        tol_residual=None, precond=precond,
    )
    return r["phi_ext"], r["gamma"]


def solve_multivalued_newton(
    mvop,
    mesh,
    m_inf: float,
    alpha_deg: float = 1.25,
    u_inf: float = 1.0,
    gamma_air: float = 1.4,
    upwind_c: float = 1.5,
    m_crit: float = 0.95,
    m_cap: float = 3.0,
    rho_floor: float = 0.05,
    farfield: str = "neumann",
    vortex_center=(0.25, 0.0),
    phi_init: Optional[np.ndarray] = None,
    gamma_init: float = 0.0,
    n_seed: int = 40,
    seed_damping_theta: Optional[float] = 0.2,
    seed_omega_rho: float = 0.5,
    n_newton_max: int = 40,
    tol_residual: float = 1e-10,
    lam_min: float = 0.05,
    freeze: bool = True,
    verbose: bool = False,
    tip_taper: Optional[np.ndarray] = None,
    precond: Optional[str] = None,
    seed_precond: Optional[str] = None,
    linear_rtol: float = 1e-10,
    gmres_restart: int = 60,
    gmres_maxiter: int = 50,
    amg_rebuild_every: int = 2,
) -> Dict[str, object]:
    """Newton solve on the level-set path (B6-Newton; design_track_b.md 5.5).

    Args:
        phi_init: warm-start extended state (e.g. the previous Mach level's).
            None -> seed with a short B6 Picard-LS solve.
        n_seed: Picard-LS seed iterations when phi_init is None (or to settle
            a warm start onto the current Mach before Newton).
        farfield: "neumann" (B6 transonic default -- no Gamma feedback, the
            far-field column vanishes), "vortex" (option a; adds the low-rank
            Gamma(z) far-field column -- handled by refreshing the RHS each
            Newton step, frozen within the step) or "freestream".
        freeze: freeze the per-side upwind selection between Newton steps is
            NOT yet needed at 2.5D (the seed lands close and live re-selection
            is stable); the flag is reserved. Live selection each step.
        tol_residual: converged when max|R_free| < tol.
        precond (B11): inner linear solver -- None (default) = the bit-identical
            sparse-direct `spsolve`; **use "ilu"** for the iterative escape from
            the true-3D splu wall (factors the real fused Jacobian, converges).
            "amg" is wired but STALLS on this wake_ls operator (measured; §5.3),
            not the Newton escape. Inexact steps (on_fail="return") are absorbed
            by the existing backtracking line search. `seed_precond` (default
            None) is the same knob for the Picard warm-start solve.
            `linear_rtol=1e-10` keeps the Newton terminus unperturbed at the
            tol_residual gate.

    Returns: dict with phi_ext, phi, gamma, cl_kj, te_jump, converged,
    residual_history, n_newton, n_limited, n_floored, mach2_max, nu_max,
    precond, n_gmres_total, n_gmres_stalled (B11 monitors).
    """
    if farfield not in ("neumann", "vortex", "freestream"):
        raise ValueError(f"farfield={farfield!r} unknown")
    if precond not in (None, "ilu", "amg"):
        raise ValueError(f"precond={precond!r} unknown (None|'ilu'|'amg')")
    beta = float(np.sqrt(1.0 - m_inf**2))

    # --- far-field split (fixed across the solve) --------------------------
    if farfield == "neumann":
        ff_nodes = _farfield_split(mesh, alpha_deg, u_inf)[3]      # inflow
        ff_vals = freestream_phi(mesh.nodes[ff_nodes], alpha_deg, u_inf)
        b_base = _neumann_outlet_rhs(mesh, alpha_deg, u_inf, mvop.n_total)
    else:
        ff_nodes = np.unique(mesh.boundary_faces["farfield"])
        b_base = np.zeros(mvop.n_total)
        ff_vals = None  # (re)built per step for the vortex option

    # --- warm start --------------------------------------------------------
    phi_ext, gamma = _seed_from_picard(
        mvop, mesh, m_inf, alpha_deg, u_inf, gamma_air, upwind_c, m_crit,
        farfield, phi_init, gamma_init, n_seed, seed_damping_theta,
        seed_omega_rho, precond=seed_precond)

    te_aux = mvop.cm.ext_dof_of_node[mvop.cm.te_nodes].astype(np.int64)
    te_main = np.asarray(mvop.cm.te_nodes, dtype=np.int64)
    n_total = mvop.n_total

    # B8 tip taper (roadmap Track B B8; P13/G13.2 finding (8)). The blended
    # TE Kutta residual is F_i*(|q_u|^2-|q_l|^2) + (1-F_i)*(phi_aux-phi_main)
    # and its Jacobian is F*J_kutta + (1-F)*weld. tip_taper=None (or all ones)
    # -> pure pressure Kutta, bit-identical. Because D_keep zeros the TE aux
    # rows of the Picard matrix and R[te_aux] is overwritten below, the blend
    # is applied ONLY here -- assemble_matrix keeps tip_taper=None.
    tt = None
    if tip_taper is not None:
        tt = np.asarray(tip_taper, dtype=np.float64)
        if tt.shape != te_main.shape:
            raise ValueError(
                f"tip_taper must be ({te_main.size},), got {tt.shape}")
        if np.all(tt == 1.0):
            tt = None
    f_of_row = np.ones(n_total, dtype=np.float64)   # F_i keyed by TE aux dof
    if tt is not None:
        f_of_row[te_aux] = tt

    def _te_residual(phi):
        r = te_kutta_residual(mvop, phi)
        if tt is None:
            return r
        return tt * r + (1.0 - tt) * (phi[te_aux] - phi[te_main])
    is_dir = np.zeros(n_total, dtype=bool)
    is_dir[ff_nodes] = True
    free = np.flatnonzero(~is_dir)

    # Terms 2/3 are the Jacobian of the MASS-conservation residual, so their
    # rows must follow the SAME row mapping assemble_matrix applies to
    # mass_conservation_coo (design_track_b.md 5.5): the non-TE aux rows are
    # DROPPED (replaced by the phi-independent wake LS) and the TE aux rows are
    # REROUTED to the TE main row (their mass balance moved there when the TE
    # aux row became the Kutta condition). Without this, Terms 2/3 leak the
    # dropped mass Jacobian onto the linear wake-LS rows (FD-caught).
    nonte_aux = np.asarray(mvop._nonte_rows, dtype=np.int64)
    row_map = np.arange(n_total, dtype=np.int64)
    row_map[nonte_aux] = -1                        # drop
    row_map[te_aux] = mvop.cm.te_nodes             # reroute TE aux -> TE main

    def _remap(r, c, d):
        if not len(r):
            return r, c, d
        rr = row_map[r]
        keep = rr >= 0
        return rr[keep], c[keep], d[keep]
    # a diagonal that zeros the TE aux rows (to be replaced by the exact Kutta
    # Jacobian) while keeping every other row of the Picard matrix
    keep_row = np.ones(n_total, dtype=np.float64)
    keep_row[te_aux] = 0.0
    D_keep = sp.diags(keep_row)

    residual_history = []
    converged = False
    M_pre = None
    n_gmres_total = 0
    n_gmres_stalled = 0

    for it in range(n_newton_max):
        # far-field Dirichlet values (vortex option refreshes with gamma)
        if farfield == "vortex":
            ff_nodes_v, ff_vals = _farfield_main(mesh, alpha_deg, gamma, u_inf,
                                                 vortex_center, beta)
            phi_ext[ff_nodes_v] = ff_vals
            dir_nodes = ff_nodes_v
        else:
            phi_ext[ff_nodes] = ff_vals
            dir_nodes = ff_nodes

        # per-side density + frozen-selection sensitivities at the current phi
        up, lo = mvop.newton_side_data(phi_ext, m_inf, upwind_c, m_crit,
                                       gamma_air, u_inf, m_cap, rho_floor)
        rho_up, rho_lo = up["rho_tilde"], lo["rho_tilde"]

        # Picard matrix at phi (mass + wake-LS + frozen-mean Kutta row)
        A = mvop.assemble_matrix(rho_tilde=(rho_up, rho_lo), closure="wake_ls",
                                 te_kutta="pressure", phi_ext=phi_ext)

        # residual R = A phi - b (exact on every row incl. the Kutta row);
        # overwrite the Kutta rows with the true nonlinear value for safety
        R = A @ phi_ext - b_base
        R[te_aux] = _te_residual(phi_ext)
        res = float(np.max(np.abs(R[free])))
        residual_history.append(res)
        if verbose:
            print(f"    newton {it}: |R|={res:.3e} gamma={gamma:.5f} "
                  f"lim/flr={mvop.n_limited}/{mvop.n_floored}")
        if res < tol_residual and mvop.n_limited == 0 and mvop.n_floored == 0:
            converged = True
            break

        # Jacobian: Picard matrix, minus the frozen Kutta rows, plus the exact
        # Kutta rows, plus Terms 2+3 on the mass rows (per side).
        J = (D_keep @ A).tocoo()
        parts_r = [J.row]
        parts_c = [J.col]
        parts_d = [J.data]
        # blended TE Kutta Jacobian: F*J_kutta (+ (1-F)*weld appended below)
        kr, kc, kd = te_kutta_jacobian_coo(mvop, phi_ext)
        if tt is not None and len(kr):
            kd = kd * f_of_row[kr]
        jac_parts = [(kr, kc, kd)]
        if tt is not None:
            jac_parts.append(te_weld_coo(mvop.cm, tt))     # (1-F)*[aux-main]
        jac_parts += [
            _remap(*newton_terms23_side_coo(mvop.op, up, u_inf)),
            _remap(*newton_terms23_side_coo(mvop.op, lo, u_inf)),
        ]
        for r_, c_, d_ in jac_parts:
            if len(r_):
                parts_r.append(r_)
                parts_c.append(c_)
                parts_d.append(d_)
        J = sp.coo_matrix((np.concatenate(parts_d),
                           (np.concatenate(parts_r), np.concatenate(parts_c))),
                          shape=(n_total, n_total)).tocsr()

        # reduce Dirichlet and solve J_free d = -R_free (nonsymmetric)
        if precond is None:
            J_free = J[free][:, free].tocsc()
            d_free = np.atleast_1d(spla.spsolve(J_free, -R[free]))
        else:
            # B11: preconditioned GMRES on the fused Jacobian. Inexact steps
            # (on_fail="return") are caught by the backtracking line search.
            J_free = J[free][:, free].tocsr()
            if precond == "ilu":
                M_pre = build_ilu_preconditioner(J_free)
            elif it % amg_rebuild_every == 0 or M_pre is None:
                M_pre = _amg_surrogate_preconditioner(
                    mvop, J, (rho_up, rho_lo), free)
            d_free, n_it, info = solve_gmres(
                J_free, -R[free], M=M_pre, rtol=linear_rtol,
                restart=gmres_restart, maxiter=gmres_maxiter, on_fail="return")
            d_free = np.atleast_1d(d_free)
            n_gmres_total += n_it
            n_gmres_stalled += int(info != 0)

        # backtracking line search on the free-DOF residual (safety only)
        lam = 1.0
        base = res
        phi_try = phi_ext.copy()
        while lam >= lam_min:
            phi_try[free] = phi_ext[free] + lam * d_free
            phi_try[dir_nodes] = phi_ext[dir_nodes]
            up_t, lo_t = mvop.newton_side_data(phi_try, m_inf, upwind_c, m_crit,
                                               gamma_air, u_inf, m_cap,
                                               rho_floor)
            A_t = mvop.assemble_matrix(rho_tilde=(up_t["rho_tilde"],
                                                  lo_t["rho_tilde"]),
                                       closure="wake_ls", te_kutta="pressure",
                                       phi_ext=phi_try)
            R_t = A_t @ phi_try - b_base
            R_t[te_aux] = _te_residual(phi_try)
            if np.max(np.abs(R_t[free])) < base or lam <= lam_min:
                break
            lam *= 0.5
        phi_ext[free] = phi_ext[free] + lam * d_free
        gamma = float(np.mean(mvop.te_jump(phi_ext)))

    mach2_max = float(np.max(mvop.element_mach2(phi_ext, m_inf, gamma_air,
                                               u_inf)))
    return {
        "phi_ext": phi_ext,
        "phi": mvop.main_potential(phi_ext),
        "gamma": gamma,
        "cl_kj": 2.0 * gamma / u_inf,
        "te_jump": mvop.te_jump(phi_ext),
        "converged": converged,
        "residual_history": residual_history,
        "n_newton": len(residual_history),
        "n_limited": mvop.n_limited,
        "n_floored": mvop.n_floored,
        "nu_max": mvop.nu_max,
        "mach2_max": mach2_max,
        "farfield": farfield,
        "precond": precond,
        "n_gmres_total": n_gmres_total,
        "n_gmres_stalled": n_gmres_stalled,
    }
