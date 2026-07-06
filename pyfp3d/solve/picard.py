"""
P1 driver: Laplace (rho == 1) is linear, so the Picard loop degenerates to a
single assemble + solve -- no outer iteration, no relaxation.

Reference: docs/roadmap.md P1 ("Picard-degenerate driver, single linear
solve"); design.md Sec 8 (full Picard loop lands in P3+).
"""

from typing import Dict, Optional

import numpy as np

from pyfp3d.kernels.residual import assemble_residual, assemble_stiffness_matrix
from pyfp3d.solve.linear import apply_dirichlet, solve_cg_amg


def solve_laplace(
    nodes: np.ndarray,
    elements: np.ndarray,
    dirichlet_nodes: np.ndarray,
    dirichlet_values: np.ndarray,
    body_source_rhs: Optional[np.ndarray] = None,
    rtol: float = 1e-10,
    maxiter: int = 500,
) -> Dict[str, object]:
    """
    Solve the Laplace problem A phi = b with Dirichlet BCs (far field) and
    natural (do-nothing) BCs elsewhere (solid walls, symmetry planes).

    Args:
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) tetrahedral connectivity
        dirichlet_nodes: node indices with prescribed phi (e.g. far field)
        dirichlet_values: prescribed phi values at dirichlet_nodes
        body_source_rhs: optional (n_nodes,) assembled RHS vector -- an MMS
            load vector, or a wall-flux boundary correction
            (solve/wall_correction.py); the physical full-potential equation
            itself has no volumetric source term
        rtol, maxiter: CG convergence controls (see solve.linear.solve_cg_amg)

    Returns:
        dict with keys: phi, n_cg_iterations, residual_norm
    """
    n_nodes = len(nodes)
    A = assemble_stiffness_matrix(nodes, elements)
    b = np.zeros(n_nodes, dtype=np.float64) if body_source_rhs is None else np.asarray(
        body_source_rhs, dtype=np.float64
    )

    A_free, b_free, free, phi = apply_dirichlet(A, b, dirichlet_nodes, dirichlet_values)
    x_free, n_iter = solve_cg_amg(A_free, b_free, rtol=rtol, maxiter=maxiter)
    phi[free] = x_free

    # Residual of the equilibrium equation b - A@phi == 0, restricted to the
    # free dofs. Dirichlet (far-field) rows are excluded: assemble_residual()
    # returns the raw operator A@phi with no boundary condition applied, so
    # at a Dirichlet node it is the natural-BC flux imbalance -- an O(1)
    # quantity that never vanishes and isn't part of the system actually
    # being solved (that row is overridden by the prescribed value instead).
    residual_full = b - assemble_residual(nodes, elements, phi)
    residual_norm = float(np.max(np.abs(residual_full[free])))

    return {
        "phi": phi,
        "n_cg_iterations": n_iter,
        "residual_norm": residual_norm,
    }


def solve_laplace_lifting(
    mesh_cut,
    wc,
    alpha_deg: float = 0.0,
    u_inf: float = 1.0,
    vortex_center=(0.25, 0.0),
    gamma_fixed: Optional[np.ndarray] = None,
    omega_gamma: float = 0.9,
    tol_gamma: float = 1e-8,
    max_kutta_updates: int = 30,
    rtol: float = 1e-10,
    maxiter: int = 2000,
) -> Dict[str, object]:
    """
    Incompressible lifting Laplace solve on a wake-cut mesh: outer Kutta
    loop around a linear solve whose matrix is assembled ONCE (Gamma is
    RHS-only, design.md Sec 4/Sec 8; roadmap P2).

    Per Kutta update (skipped if gamma_fixed is given):
        solve  ->  Gamma_target_j = phi[TE upper probe] - phi[TE lower probe]
    then update Gamma. The first update is plain under-relaxation
    (omega_gamma); from the second on, a per-station SECANT (Aitken)
    step jumps to the affine map's fixed point directly. Measured on the
    M0 NACA0012 meshes the map's slope is b ~ 0.93 (the probes sit one
    node off the TE, so the smooth-flow jump there is nearly Gamma
    itself), which would need O(100) plain relaxed updates for 1e-8 --
    the secant reaches it in 3-4 solves because the Laplace problem is
    exactly affine in Gamma. The under-relaxed form stays as the fallback
    (and the P3+ nonlinear outer loop can re-tune omega_gamma).
    Convergence: ||Gamma_target - Gamma||_inf < tol_gamma (the actual
    Kutta mismatch); the G2.3 gate budget is < 20 updates.

    Args:
        mesh_cut, wc: output of mesh/wake_cut.cut_wake()
        alpha_deg, u_inf: freestream incidence (deg) and speed
        vortex_center: 2D point-vortex location for the far-field correction
        gamma_fixed: prescribe Gamma per station and SKIP the Kutta loop
            (G2.1/G2.2 use this); scalar or (n_stations,)
        omega_gamma: Kutta under-relaxation (design.md Sec 4: 0.7-1.0)
        tol_gamma: convergence tolerance on ||Delta Gamma||_inf
        max_kutta_updates: hard cap on outer updates
        rtol, maxiter: inner CG controls

    Returns:
        dict: phi (cut-mesh nodal potential), gamma (per station),
        n_kutta_updates, gamma_history, residual_norm (inf-norm of the
        reduced free-dof residual at the final state), n_cg_total,
        kutta_converged
    """
    from pyfp3d.constraints.dirichlet import farfield_dirichlet
    from pyfp3d.constraints.wake import WakeConstraint, kutta_targets

    A = assemble_stiffness_matrix(mesh_cut.nodes, mesh_cut.elements)
    con = WakeConstraint(A, wc)
    n_red = con.n_reduced
    b_zero = np.zeros(A.shape[0], dtype=np.float64)

    # Dirichlet pattern is Gamma-independent: fix the free/constrained
    # split and the AMG hierarchy once, update only values per iteration.
    dir_nodes, _ = farfield_dirichlet(
        mesh_cut, wc, alpha_deg, np.zeros(wc.n_stations), u_inf, vortex_center
    )
    dir_red, _ = con.to_reduced_dirichlet(dir_nodes, np.zeros(len(dir_nodes)))
    is_dir = np.zeros(n_red, dtype=bool)
    is_dir[dir_red] = True
    free = np.where(~is_dir)[0]
    A_csr = con.A_reduced
    A_free = A_csr[free][:, free].tocsr()
    A_coupling = A_csr[free][:, dir_red].tocsr()

    import pyamg
    import scipy.sparse.linalg as spla

    ml = pyamg.smoothed_aggregation_solver(A_free)
    M = ml.aspreconditioner()

    def _solve_for(gamma: np.ndarray):
        nodes_d, vals_d = farfield_dirichlet(
            mesh_cut, wc, alpha_deg, gamma, u_inf, vortex_center
        )
        _, vals_red = con.to_reduced_dirichlet(nodes_d, vals_d)
        b_red = con.reduced_rhs(b_zero, gamma)
        b_free = b_red[free] - A_coupling @ vals_red
        n_it = [0]
        x, info = spla.cg(
            A_free, b_free, M=M, rtol=rtol, maxiter=maxiter,
            callback=lambda _x: n_it.__setitem__(0, n_it[0] + 1),
        )
        if info != 0:
            raise RuntimeError(f"CG did not converge (info={info})")
        phi_red = np.empty(n_red, dtype=np.float64)
        phi_red[free] = x
        phi_red[dir_red] = vals_red
        res = float(np.max(np.abs(b_free - A_free @ x))) if len(x) else 0.0
        return phi_red, res, n_it[0]

    if gamma_fixed is not None:
        gamma = np.broadcast_to(
            np.atleast_1d(np.asarray(gamma_fixed, dtype=np.float64)),
            (wc.n_stations,),
        ).copy()
        phi_red, res, n_cg = _solve_for(gamma)
        return {
            "phi": con.expand(phi_red, gamma),
            "gamma": gamma,
            "n_kutta_updates": 0,
            "gamma_history": [gamma.copy()],
            "residual_norm": res,
            "n_cg_total": n_cg,
            "kutta_converged": True,
        }

    gamma = np.zeros(wc.n_stations, dtype=np.float64)
    history = [gamma.copy()]
    n_cg_total = 0
    converged = False
    prev_gamma = prev_target = None
    phi_red, res = None, np.inf
    for _ in range(max_kutta_updates):
        phi_red, res, n_cg = _solve_for(gamma)
        n_cg_total += n_cg
        target = kutta_targets(con.expand(phi_red, gamma), wc)
        if float(np.max(np.abs(target - gamma))) < tol_gamma:
            converged = True
            break
        if prev_gamma is None:
            new_gamma = gamma + omega_gamma * (target - gamma)
        else:
            # Secant on the affine map target(Gamma) = a + b Gamma:
            # fixed point Gamma* = (target - b Gamma) / (1 - b).
            dg = gamma - prev_gamma
            b = np.where(
                np.abs(dg) > 1e-14, (target - prev_target) / np.where(dg == 0, 1, dg), 0.0
            )
            denom = 1.0 - b
            relaxed = gamma + omega_gamma * (target - gamma)
            new_gamma = np.where(
                np.abs(denom) > 1e-3, (target - b * gamma) / np.where(denom == 0, 1, denom),
                relaxed,
            )
        prev_gamma, prev_target = gamma, target
        gamma = new_gamma
        history.append(gamma.copy())

    if not converged and phi_red is not None:
        # Loop exhausted with gamma updated after the last solve: resolve
        # once so the returned phi is consistent with the returned gamma.
        phi_red, res, n_cg = _solve_for(gamma)
        n_cg_total += n_cg

    return {
        "phi": con.expand(phi_red, gamma),
        "gamma": gamma,
        "n_kutta_updates": len(history) - 1,
        "gamma_history": history,
        "residual_norm": res,
        "n_cg_total": n_cg_total,
        "kutta_converged": converged,
    }
