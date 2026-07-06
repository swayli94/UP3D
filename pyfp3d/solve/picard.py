"""
Picard drivers.

P1: Laplace (rho == 1) is linear, so the loop degenerates to a single
assemble + solve. P2: lifting Laplace, matrix assembled once, Kutta outer
loop with secant acceleration. P3: the density Picard loop of design.md
Sec 8 (subsonic: nu == 0, rho_tilde == rho) -- `solve_subsonic`
(non-lifting) and `solve_subsonic_lifting` (NESTED: outer density update,
inner P2 secant Kutta at frozen rho; see its docstring for why nesting
replaced the interleaved form).

Reference: docs/roadmap.md P1-P3; design.md Sec 8 (Nonlinear solution
strategy).
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

    import scipy.sparse.linalg as spla

    from pyfp3d.solve.linear import build_amg_preconditioner

    _, M = build_amg_preconditioner(A_free)

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


def solve_subsonic(
    nodes: np.ndarray,
    elements: np.ndarray,
    dirichlet_nodes: np.ndarray,
    dirichlet_values: np.ndarray,
    m_inf: float,
    gamma_air: float = 1.4,
    u_inf: float = 1.0,
    phi_init: Optional[np.ndarray] = None,
    omega: float = 1.0,
    n_picard_max: int = 40,
    tol_rho: float = 1e-10,
    rtol: float = 1e-12,
    maxiter: int = 3000,
) -> Dict[str, object]:
    """
    Non-lifting subsonic full-potential solve: relaxed-Picard density outer
    loop (design.md Sec 8, nu == 0 so rho_tilde == rho and every Picard
    matrix is SPD). Per iteration: element velocity sweep -> isentropic
    density -> colored assembly of A(rho) -> CG+AMG solve -> relax.

    Convergence is on the density lag ||rho(phi_new) - rho_matrix||_inf <
    tol_rho -- the exact quantity Picard freezes -- with the nonlinear
    residual ||R||_inf on free dofs recorded per iteration (gate G3.2
    monotonicity evidence). At m_inf = 0 the density is exactly 1.0
    everywhere (physics/isentropic.py::density_field), so iteration 1
    reproduces the P1 Laplace solve bit-identically and the loop exits
    after it (gate G3.3).

    Args:
        nodes, elements: mesh arrays
        dirichlet_nodes, dirichlet_values: far-field Dirichlet data
        m_inf: freestream Mach number
        gamma_air: specific heat ratio
        u_inf: freestream speed (q^2 is nondimensionalized by u_inf^2
            before entering the density law)
        phi_init: initial potential; None = freestream-like start with
            rho == 1 on the first matrix (design.md Sec 8 initial guess)
        omega: Picard under-relaxation (1.0 subsonic; omega == 1.0 replaces
            phi exactly instead of computing phi + 1.0*(phi_new - phi))
        n_picard_max: outer-iteration cap
        tol_rho: density-lag convergence tolerance
        rtol, maxiter: inner CG controls

    Returns:
        dict: phi, n_picard, converged, residual_history (per-iteration
        nonlinear ||R||_inf on free dofs), drho_history, rho (element
        densities), mach2_max, n_cg_total
    """
    import scipy.sparse.linalg as spla

    from pyfp3d.kernels.jacobian import PicardOperator
    from pyfp3d.physics.isentropic import density_field, mach_squared_field
    from pyfp3d.solve.linear import build_amg_preconditioner

    op = PicardOperator(nodes, elements)
    n_nodes = op.n_nodes

    is_dir = np.zeros(n_nodes, dtype=bool)
    is_dir[dirichlet_nodes] = True
    free = np.where(~is_dir)[0]
    vals_d = np.asarray(dirichlet_values, dtype=np.float64)

    if phi_init is not None:
        phi = np.asarray(phi_init, dtype=np.float64).copy()
        phi[dirichlet_nodes] = vals_d
        _, q2 = op.velocities(phi)
        rho = density_field(q2 / u_inf**2, m_inf, gamma_air).copy()
    else:
        phi = np.zeros(n_nodes, dtype=np.float64)
        phi[dirichlet_nodes] = vals_d
        rho = np.ones(op.n_tets, dtype=np.float64)

    residual_history = []
    drho_history = []
    n_cg_total = 0
    converged = False

    for _ in range(n_picard_max):
        A = op.assemble_matrix(rho)
        A_csr = A.tocsr()
        A_free = A_csr[free][:, free].tocsr()
        b_free = -A_csr[free][:, dirichlet_nodes] @ vals_d

        _, M = build_amg_preconditioner(A_free)
        n_it = [0]
        x, info = spla.cg(
            A_free, b_free, M=M, rtol=rtol, maxiter=maxiter,
            callback=lambda _x: n_it.__setitem__(0, n_it[0] + 1),
        )
        if info != 0:
            raise RuntimeError(f"CG did not converge (info={info})")
        n_cg_total += n_it[0]

        phi_new = np.empty(n_nodes, dtype=np.float64)
        phi_new[free] = x
        phi_new[dirichlet_nodes] = vals_d
        if omega == 1.0:
            phi = phi_new
        else:
            phi = phi + omega * (phi_new - phi)

        _, q2 = op.velocities(phi)
        rho_new = density_field(q2 / u_inf**2, m_inf, gamma_air).copy()
        R = op.assemble_residual(phi, rho_new)
        residual_history.append(float(np.max(np.abs(R[free]))))
        drho = float(np.max(np.abs(rho_new - rho)))
        drho_history.append(drho)
        rho = rho_new

        if drho < tol_rho:
            converged = True
            break

    _, q2 = op.velocities(phi)
    mach2_max = float(np.max(mach_squared_field(q2 / u_inf**2, m_inf, gamma_air)))

    return {
        "phi": phi,
        "n_picard": len(residual_history),
        "converged": converged,
        "residual_history": residual_history,
        "drho_history": drho_history,
        "rho": rho,
        "mach2_max": mach2_max,
        "n_cg_total": n_cg_total,
    }


def solve_subsonic_lifting(
    mesh_cut,
    wc,
    m_inf: float,
    alpha_deg: float = 0.0,
    u_inf: float = 1.0,
    gamma_air: float = 1.4,
    vortex_center=(0.25, 0.0),
    gamma_fixed: Optional[np.ndarray] = None,
    omega: float = 1.0,
    omega_gamma: float = 0.9,
    tol_gamma: float = 1e-8,
    tol_rho: float = 1e-10,
    n_picard_max: int = 40,
    rtol: float = 1e-10,
    maxiter: int = 3000,
    amg_rebuild_every: int = 4,
    forcing: float = 0.0,
    upwind_c: float = 1.5,
    m_crit: float = 0.95,
    phi_init: Optional[np.ndarray] = None,
    gamma_init: Optional[np.ndarray] = None,
    kutta_per_outer: Optional[int] = None,
    m_cap: float = 3.0,
    omega_rho: float = 1.0,
    pseudo_dt: Optional[float] = None,
    pseudo_dt_growth: float = 1.1,
    pseudo_dt_max_ratio: float = 100.0,
    tol_residual: Optional[float] = None,
) -> Dict[str, object]:
    """
    Lifting subsonic full-potential solve on a wake-cut mesh: NESTED
    Picard (design.md Sec 8) with the Prandtl-Glauert-scaled vortex far
    field (beta = sqrt(1 - M_inf^2), constraints/dirichlet.py).

    Outer loop = density update: velocity sweep -> rho -> colored assembly
    of A(rho) -> wake reduction A_red = T^T A T (T built once, matrix
    updated in place). Inner loop = the P2 secant Kutta iteration at
    FROZEN rho (where the map target(Gamma) is exactly affine, so the
    secant reasoning of solve_laplace_lifting applies verbatim); matrix
    and AMG hierarchy are per-outer, reused across the inner solves.
    Nesting (rather than interleaving one Gamma step per density step) is
    what makes the G3.2 "monotone residual" criterion meaningful: each
    recorded residual sits at a Kutta-converged state, so the history
    tracks the density-lag contraction alone -- interleaved Gamma jumps
    were measured to inject 10x residual spikes all the way down.

    Accelerations (design.md Sec 8 "in order of implementation value"):
    the AMG hierarchy is rebuilt only every `amg_rebuild_every` outer
    iterations (rho drifts little between rebuilds), and for m_inf > 0
    each inner CG is warm-started from the previous solution. With
    `forcing` = eta > 0 the inner solves also use an inexact-Picard
    forcing term, atol = eta*||b - A x0|| (cut your own starting residual
    by 1/eta; final accuracy still governed by rtol) -- measured ~2x
    faster at eta = 0.05 on the medium NACA case, at the price of
    bounded eta-noise wiggles in the residual tail. The default eta = 0
    keeps the inner solves exact so the G3.2 "monotone residual"
    criterion holds as written (measured strictly monotone to the 5e-11
    floor on the medium mesh). NOTE a loose RELATIVE inner tolerance is
    NOT a valid alternative: warm-started CG then exits at x0 without
    ever computing the density correction (measured false convergence).
    All accelerations are exact no-ops at m_inf = 0 (atol = 0, x0 = None,
    matrix never changes so the reused hierarchy equals a rebuilt one),
    preserving gate G3.3: at m_inf = 0, rho == 1.0 bitwise and beta ==
    1.0 bitwise, so outer iteration 1 IS solve_laplace_lifting
    bit-for-bit and the loop exits after it.

    P4 transonic support (design.md Sec 3): when `upwind_c` > 0 and
    m_inf > 0, the matrix/residual weight is the UPWINDED density
    rho_tilde_e = rho_e - nu_e (rho_e - rho_u(e)) with the (3.2) switch
    nu = upwind_c * max(0, 1 - m_crit^2/M^2) (kernels/upwind.py). In the
    subcritical range nu == 0.0 exactly and rho_tilde == rho bitwise, so
    a subcritical solve with the machinery active is bit-identical to the
    P3 path (gate G4.2); upwind_c = 0 (or m_inf = 0) bypasses the sweep
    entirely. `phi_init`/`gamma_init` seed the loop for Mach continuation
    (solve/continuation.py); phi_init also warm-starts the first inner CG.

    `kutta_per_outer`: None (default) converges the inner secant Kutta
    loop per outer iteration -- correct at subcritical, where the frozen-
    rho map is benign. At transonic this is UNSTABLE: solving the Kutta
    condition exactly against a stale rho_tilde overshoots Gamma, the
    shock strengthens in response, and the composite outer map runs away
    (measured on the coarse NACA mesh, M 0.70 -> 0.75 continuation step:
    Gamma escalated 0.115 -> 4.99). Transonic runs pass kutta_per_outer=1:
    ONE plain under-relaxed Gamma step per density update (no secant --
    its slope estimate is polluted when the matrix changes between
    solves), which co-converges Gamma with the shock. Transonic tuning
    ladder (design.md Sec 12.4): raise upwind_c -> lower omega (and
    omega_gamma) -> Mach continuation.

    Converged when the density lag ||rho_tilde(phi_new) -
    rho_tilde_matrix||_inf < tol_rho AND the inner Kutta loop met
    tol_gamma.

    Returns:
        dict: phi (cut-mesh), gamma, n_picard (outer iterations),
        converged, kutta_converged, residual_history (per-outer nonlinear
        ||T^T R||_inf on free reduced dofs), drho_history, gamma_history,
        rho (UPWINDED element densities actually used as the matrix
        weight), mach2_max, nu_max, n_nu_active (elements with nu > 0),
        n_cg_total, n_solves_total
    """
    import scipy.sparse.linalg as spla

    from pyfp3d.constraints.dirichlet import farfield_dirichlet, freestream_phi
    from pyfp3d.constraints.wake import WakeConstraint, kutta_targets
    from pyfp3d.kernels.jacobian import PicardOperator
    from pyfp3d.physics.isentropic import (
        density_field,
        limit_q2_field,
        mach_squared_field,
    )
    from pyfp3d.solve.linear import build_amg_preconditioner

    if not 0.0 <= m_inf < 1.0:
        raise ValueError(f"solve_subsonic_lifting needs 0 <= M_inf < 1, got {m_inf}")
    beta = float(np.sqrt(1.0 - m_inf**2))

    op = PicardOperator(mesh_cut.nodes, mesh_cut.elements)
    use_upwind = upwind_c > 0.0 and m_inf > 0.0
    if use_upwind:
        from pyfp3d.kernels.upwind import UpwindOperator
        upw = UpwindOperator(mesh_cut.nodes, mesh_cut.elements)
    con = WakeConstraint(op.assemble_matrix(), wc)
    n_red = con.n_reduced
    b_zero = np.zeros(op.n_nodes, dtype=np.float64)

    # Dirichlet pattern is Gamma- and rho-independent: fix the split once.
    dir_nodes, _ = farfield_dirichlet(
        mesh_cut, wc, alpha_deg, np.zeros(wc.n_stations), u_inf, vortex_center,
        beta=beta,
    )
    dir_red, _ = con.to_reduced_dirichlet(dir_nodes, np.zeros(len(dir_nodes)))
    is_dir = np.zeros(n_red, dtype=bool)
    is_dir[dir_red] = True
    free = np.where(~is_dir)[0]

    # Pseudo-transient term (design.md Sec 8 acceleration 4, pulled into
    # the P4 Picard loop as the transonic stabilizer of last resort):
    # A' = A + diag(m_lumped/dtau), b' = b + diag(...) phi_k -- an implicit
    # pseudo-time step that bounds the per-iteration update and breaks the
    # shock-position limit cycle; consistent at steady state (the added
    # term vanishes when phi stops changing). None = off (P3 bit-path).
    d_tau_free = None
    m_red_free = None
    dt_cur = pseudo_dt
    if pseudo_dt is not None:
        m_lumped = np.zeros(op.n_nodes, dtype=np.float64)
        v_quarter = np.repeat(op.V / 4.0, 4)
        np.add.at(m_lumped, np.asarray(op.elements).reshape(-1), v_quarter)
        m_red_free = (con.T.T @ m_lumped)[free]
        d_tau_free = m_red_free / dt_cur

    def _solve_for(gamma: np.ndarray, A_free, A_coupling, M, x0):
        nodes_d, vals_d = farfield_dirichlet(
            mesh_cut, wc, alpha_deg, gamma, u_inf, vortex_center, beta=beta,
        )
        _, vals_red = con.to_reduced_dirichlet(nodes_d, vals_d)
        b_red = con.reduced_rhs(b_zero, gamma)
        b_free = b_red[free] - A_coupling @ vals_red
        if d_tau_free is not None and x0 is not None:
            b_free = b_free + d_tau_free * x0
        # Inexact-Picard forcing term (design.md Sec 8 acceleration 2):
        # each inner solve must cut ITS OWN starting residual by 1/eta;
        # the final accuracy is still governed by rtol (scipy cg stops at
        # max(rtol*||b||, atol)). Exact no-op at m_inf = 0 or eta = 0
        # (atol = 0.0, x0 = None): bit-identical to the P2 call -- G3.3.
        if m_inf > 0.0 and forcing > 0.0:
            r0 = b_free - A_free @ x0 if x0 is not None else b_free
            atol_k = forcing * float(np.linalg.norm(r0))
        else:
            atol_k = 0.0
        n_it = [0]
        x, info = spla.cg(
            A_free, b_free, M=M, rtol=rtol, atol=atol_k, maxiter=maxiter,
            x0=x0, callback=lambda _x: n_it.__setitem__(0, n_it[0] + 1),
        )
        if info != 0:
            raise RuntimeError(f"CG did not converge (info={info})")
        phi_red = np.empty(n_red, dtype=np.float64)
        phi_red[free] = x
        phi_red[dir_red] = vals_red
        return phi_red, n_it[0]

    if gamma_fixed is not None:
        gamma = np.broadcast_to(
            np.atleast_1d(np.asarray(gamma_fixed, dtype=np.float64)),
            (wc.n_stations,),
        ).copy()
    elif gamma_init is not None:
        gamma = np.atleast_1d(np.asarray(gamma_init, dtype=np.float64)).copy()
    else:
        gamma = np.zeros(wc.n_stations, dtype=np.float64)

    # Freestream start (design.md Sec 8), or the continuation seed.
    if phi_init is not None:
        phi_cut = np.asarray(phi_init, dtype=np.float64).copy()
    else:
        phi_cut = freestream_phi(mesh_cut.nodes, alpha_deg, u_inf)
    grad, q2 = op.velocities(phi_cut)
    q2l = limit_q2_field(q2 / u_inf**2, m_inf, m_cap, gamma_air)
    rho = density_field(q2l, m_inf, gamma_air)
    if use_upwind:
        rho_t = upw.rho_tilde(grad, q2l, rho, m_inf, upwind_c, m_crit,
                              gamma_air).copy()
    else:
        rho_t = rho

    residual_history = []
    drho_history = []
    n_limited = 0
    gamma_history = [gamma.copy()]
    n_cg_total = 0
    n_solves_total = 0
    # A continuation seed also warm-starts the first inner CG solve.
    phi_red = phi_cut[:n_red].copy() if phi_init is not None else None
    M = None
    converged = False
    kutta_converged = gamma_fixed is not None
    max_kutta_updates = 30

    for outer in range(n_picard_max):
        A = op.assemble_matrix(rho_t)
        con.update_matrix(A)
        A_csr = con.A_reduced
        A_free = A_csr[free][:, free].tocsr()
        A_coupling = A_csr[free][:, dir_red].tocsr()
        if d_tau_free is not None:
            import scipy.sparse as _sp
            A_free = (A_free + _sp.diags(d_tau_free)).tocsr()
        if M is None or outer % amg_rebuild_every == 0:
            _, M = build_amg_preconditioner(A_free)
        x0 = phi_red[free] if (m_inf > 0.0 and phi_red is not None) else None

        # Inner Kutta loop at frozen rho: exactly the P2 secant iteration
        # (the map target(Gamma) is affine when the matrix is fixed).
        # kutta_per_outer caps it (transonic: 1 relaxed step, no secant).
        prev_gamma = prev_target = None
        n_inner = (max_kutta_updates if kutta_per_outer is None
                   else kutta_per_outer)
        if gamma_fixed is not None:
            phi_red_new, n_cg = _solve_for(gamma, A_free, A_coupling, M, x0)
            n_cg_total += n_cg
            n_solves_total += 1
        else:
            for _ in range(n_inner):
                phi_red_new, n_cg = _solve_for(gamma, A_free, A_coupling, M, x0)
                n_cg_total += n_cg
                n_solves_total += 1
                if m_inf > 0.0:
                    x0 = phi_red_new[free]
                target = kutta_targets(con.expand(phi_red_new, gamma), wc)
                kutta_converged = (
                    float(np.max(np.abs(target - gamma))) < tol_gamma
                )
                if kutta_converged:
                    break
                if prev_gamma is None or kutta_per_outer is not None:
                    new_gamma = gamma + omega_gamma * (target - gamma)
                else:
                    dg = gamma - prev_gamma
                    b = np.where(
                        np.abs(dg) > 1e-14,
                        (target - prev_target) / np.where(dg == 0, 1, dg),
                        0.0,
                    )
                    denom = 1.0 - b
                    relaxed = gamma + omega_gamma * (target - gamma)
                    new_gamma = np.where(
                        np.abs(denom) > 1e-3,
                        (target - b * gamma) / np.where(denom == 0, 1, denom),
                        relaxed,
                    )
                prev_gamma, prev_target = gamma, target
                gamma = new_gamma
                gamma_history.append(gamma.copy())

        if omega == 1.0 or phi_red is None:
            phi_red = phi_red_new
        else:
            phi_red = phi_red + omega * (phi_red_new - phi_red)
        phi_cut = con.expand(phi_red, gamma)

        grad, q2 = op.velocities(phi_cut)
        q2l = limit_q2_field(q2 / u_inf**2, m_inf, m_cap, gamma_air)
        n_limited = int(np.count_nonzero(q2l != q2 / u_inf**2))
        rho_new = density_field(q2l, m_inf, gamma_air)
        if use_upwind:
            rho_t_new = upw.rho_tilde(grad, q2l, rho_new, m_inf, upwind_c,
                                      m_crit, gamma_air).copy()
        else:
            rho_t_new = rho_new
        R_red = con.T.T @ op.assemble_residual(phi_cut, rho_t_new)
        residual_history.append(float(np.max(np.abs(R_red[free]))))
        drho = float(np.max(np.abs(rho_t_new - rho_t)))
        drho_history.append(drho)
        # Density under-relaxation (transonic stabilizer; omega_rho = 1.0
        # adopts rho_t_new exactly -- bit-identical subcritical path).
        if omega_rho == 1.0:
            rho, rho_t = rho_new, rho_t_new
        else:
            rho = rho_new
            rho_t = rho_t + omega_rho * (rho_t_new - rho_t)

        # SER ramp of the pseudo-time step: grow while the density lag
        # falls (approaching pure Picard, which is locally fine near the
        # converged shock), back off on a rebound.
        if pseudo_dt is not None and len(drho_history) >= 2:
            if drho < drho_history[-2]:
                dt_cur = min(dt_cur * pseudo_dt_growth,
                             pseudo_dt * pseudo_dt_max_ratio)
            elif drho > 2.0 * drho_history[-2]:
                dt_cur = max(dt_cur / 2.0, pseudo_dt)
            d_tau_free = m_red_free / dt_cur

        # With pseudo-time the per-iteration change is bounded by dtau, so
        # a small drho alone can be damping, not convergence -- require
        # the nonlinear residual too when tol_residual is given.
        res_ok = (tol_residual is None
                  or residual_history[-1] < tol_residual)
        if kutta_converged and drho < tol_rho and res_ok:
            converged = True
            break

    _, q2 = op.velocities(phi_cut)
    mach2_max = float(np.max(mach_squared_field(q2 / u_inf**2, m_inf, gamma_air)))

    return {
        "phi": phi_cut,
        "gamma": gamma,
        "n_picard": len(residual_history),
        "converged": converged,
        "kutta_converged": kutta_converged,
        "residual_history": residual_history,
        "drho_history": drho_history,
        "gamma_history": gamma_history,
        "rho": rho_t,
        "mach2_max": mach2_max,
        "nu_max": upw.nu_max if use_upwind else 0.0,
        "n_nu_active": upw.n_supersonic if use_upwind else 0,
        "n_limited": n_limited,
        "n_floored": upw.n_floored if use_upwind else 0,
        "n_cg_total": n_cg_total,
        "n_solves_total": n_solves_total,
    }
