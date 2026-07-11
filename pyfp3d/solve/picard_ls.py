"""
Level-set (Track B) solve drivers -- parallel to solve/picard.py, never
imported by the conforming path.

B2 ships the NON-LIFTING driver only: a single extended assemble + direct
solve on the multivalued operator, enough for the B2 consistency gates
(V0 freestream, V1 MMS convergence, Laplace a=0 -> cl ~ 0). The extended
matrix is nonsymmetric (the aux-row weld couples aux -> main one-way), so
it is solved with a sparse direct LU (scipy `spsolve`), not CG/AMG -- for
the coarse/medium B2 meshes the direct factor is cheap and isolates
"is the assembly correct" from any preconditioner-convergence question.
GMRES + AMG(main block) (solve/linear.py) is the scaling path B3+ adopts
(design_track_b.md section 5.3).

The lifting driver with implicit Kutta (the g1+g2 wake LS closure, no Gamma
outer loop) is B3.
"""

from typing import Dict, Optional

import numpy as np
import scipy.sparse.linalg as spla

from pyfp3d.constraints.dirichlet import freestream_phi, vortex_phi_2d
from pyfp3d.solve.linear import apply_dirichlet
from pyfp3d.wake.multivalued import MultivaluedOperator


def solve_multivalued_laplace(
    mvop: MultivaluedOperator,
    dirichlet_nodes: np.ndarray,
    dirichlet_values: np.ndarray,
    body_source_rhs: Optional[np.ndarray] = None,
) -> Dict[str, object]:
    """Non-lifting Laplace solve on a level-set cut mesh (B2).

    Assembles the extended multivalued matrix (rho == 1) with the B2
    continuity closure and solves A x = b with Dirichlet BCs on the main
    DOFs (far field). With the weld closure the solution is single valued,
    so this reproduces the ordinary single-valued Laplace/MMS solution --
    the B2 consistency check.

    Args:
        mvop: MultivaluedOperator for the mesh + wake level set
        dirichlet_nodes: node ids with prescribed phi (< n_main); their aux
            DOFs are pinned to the same value through the weld rows
        dirichlet_values: prescribed phi at dirichlet_nodes
        body_source_rhs: optional (n_main,) consistent load vector on the
            main DOFs (an MMS source); aux rows get 0. The physical FP
            equation has no volume source.

    Returns:
        dict: phi_ext (n_total), phi (n_main main potential), te_jump
        (Gamma at TE nodes), residual_norm (inf-norm of the reduced free-DOF
        residual), n_total, n_ext.
    """
    A = mvop.assemble_matrix(rho_tilde=None, closure="continuity")
    b = np.zeros(mvop.n_total, dtype=np.float64)
    if body_source_rhs is not None:
        b[: mvop.n_main] = np.asarray(body_source_rhs, dtype=np.float64)

    dirichlet_nodes = np.asarray(dirichlet_nodes, dtype=np.int64)
    dirichlet_values = np.asarray(dirichlet_values, dtype=np.float64)
    A_free, b_free, free, phi_ext = apply_dirichlet(
        A, b, dirichlet_nodes, dirichlet_values
    )

    x = spla.spsolve(A_free.tocsc(), b_free)
    x = np.atleast_1d(x)
    phi_ext[free] = x
    residual = float(np.max(np.abs(b_free - A_free @ x))) if len(x) else 0.0

    return {
        "phi_ext": phi_ext,
        "phi": mvop.main_potential(phi_ext),
        "te_jump": mvop.te_jump(phi_ext),
        "residual_norm": residual,
        "n_total": mvop.n_total,
        "n_ext": mvop.n_ext,
    }


def _farfield_main(mesh, alpha_deg, gamma, u_inf, vortex_center, beta):
    """B-path far field = freestream + PG vortex on the far-field MAIN DOFs
    ONLY (design_track_b.md section 5.4 option a). The wake-jump aux DOFs are
    left FREE (natural/Neumann): pinning them to the vortex lower branch was
    measured to drain the circulation (the jump decayed downstream instead of
    staying constant), whereas a free outlet lets the wake-LS carry a
    constant jump and the implicit Kutta stiffen -- coarse->medium Gamma
    0.018->0.109 vs conforming 0.1175 on the NACA case."""
    ff = np.unique(mesh.boundary_faces["farfield"])
    xy = mesh.nodes[ff]
    vals = freestream_phi(xy, alpha_deg, u_inf) + vortex_phi_2d(
        xy, gamma, vortex_center, beta=beta
    )
    return ff, vals


def solve_multivalued_lifting(
    mvop: MultivaluedOperator,
    mesh,
    m_inf: float,
    alpha_deg: float = 2.0,
    u_inf: float = 1.0,
    gamma_air: float = 1.4,
    vortex_center=(0.25, 0.0),
    omega_gamma: float = 0.5,
    omega_rho: float = 0.5,
    tol_gamma: float = 1e-6,
    tol_rho: float = 1e-6,
    n_outer_max: int = 80,
    gamma_init: float = 0.0,
) -> Dict[str, object]:
    """Subsonic lifting solve on a level-set cut mesh with IMPLICIT Kutta
    (Track B, B3; design_track_b.md D2). NO Gamma secant and no master-slave
    Gamma constraint: the TE jump is carried by the multivalued aux DOFs, the
    g1+g2 wake LS holds it constant downstream, and its value emerges from the
    doubled-TE mass conservation. The only outer coupling is (a) the density
    lag (subcritical: nu == 0, rho_tilde == rho, per-side on cut elements) and
    (b) an RHS-only Gamma(z) refresh of the far-field vortex from the extracted
    TE jump -- both ride the same relaxed fixed-point loop (design_track_b.md
    section 5.4, "not a new outer loop").

    The extended system is structurally nonsymmetric (the wake-LS aux rows),
    so it is solved by sparse-direct LU. Implicit Kutta is O(h)-soft on coarse
    meshes and converges to the conforming circulation under refinement (the
    B3 gate is convergence-based, not same-mesh <1%; design_track_b.md section
    7 / roadmap B3).

    ★ Compressibility is carried by the BULK density, NOT the far-field vortex:
    measured on the medium NACA case, PG-scaling the far-field vortex (beta <
    1) leaves gamma unchanged (0.1086 -> 0.1086 -- the soft Kutta does not
    propagate the outer-vortex stretch), while the isentropic bulk density
    alone raises it 0.1086 -> 0.1256 (the physical compressible lift rise, ~93%
    of the conforming M0.5 gamma, the SAME convergence ratio as incompressible).
    The cut-strip per-side density limit-cycles, so it is under-relaxed
    (omega_rho); at omega_rho = 1 the loop diverges after ~80 iterations.

    Returns:
        dict: phi_ext, phi (main), gamma (extracted TE circulation),
        cl_kj (= 2 gamma / (u_inf c), c = 1), te_jump, n_outer, converged,
        gamma_history, drho_history, residual_norm, mach2_max.
    """
    if not 0.0 <= m_inf < 1.0:
        raise ValueError(f"needs 0 <= M_inf < 1, got {m_inf}")
    beta = float(np.sqrt(1.0 - m_inf**2))
    op = mvop.op

    phi_ext = np.zeros(mvop.n_total, dtype=np.float64)
    phi_ext[: mvop.n_main] = freestream_phi(mesh.nodes, alpha_deg, u_inf)
    gamma = float(gamma_init)
    gamma_history = [gamma]
    drho_history = []
    rho_up = rho_lo = None
    residual = np.inf
    converged = False

    for outer in range(n_outer_max):
        A = mvop.assemble_matrix(
            rho_tilde=(None if rho_up is None else (rho_up, rho_lo)),
            closure="wake_ls",
        )
        ff, vals = _farfield_main(mesh, alpha_deg, gamma, u_inf,
                                  vortex_center, beta)
        A_free, b_free, free, phi_ext = apply_dirichlet(
            A, np.zeros(mvop.n_total), ff, vals
        )
        x = np.atleast_1d(spla.spsolve(A_free.tocsc(), b_free))
        phi_ext[free] = x
        residual = float(np.max(np.abs(b_free - A_free @ x))) if len(x) else 0.0

        # Gamma refresh (implicit Kutta -> extracted TE jump), relaxed.
        gamma_new = float(np.mean(mvop.te_jump(phi_ext)))
        dgamma = abs(gamma_new - gamma)
        gamma = (1.0 - omega_gamma) * gamma + omega_gamma * gamma_new
        gamma_history.append(gamma)

        # Density lag (subcritical per-side), UNDER-RELAXED. The per-side
        # density on the thin cut strip limit-cycles; adopting it fully
        # (omega_rho = 1) holds for ~80 outer iterations then diverges (the
        # medium NACA M0.5 case collapsed gamma 0.126 -> 0.010, M_max
        # blowing up), so it is relaxed like the conforming transonic path.
        # drho is measured on the OWN-SIDE (junk-free) density -- the raw
        # rho_upper/rho_lower carry other-side garbage that never enters the
        # matrix but would pollute the lag.
        rho_up_new, rho_lo_new = mvop.element_densities(
            phi_ext, m_inf, gamma_air, u_inf
        )
        if rho_up is None:
            drho = 0.0
            rho_up, rho_lo = rho_up_new, rho_lo_new
        else:
            used_new = mvop.own_side_field(rho_up_new, rho_lo_new)
            used_old = mvop.own_side_field(rho_up, rho_lo)
            drho = float(np.max(np.abs(used_new - used_old)))
            rho_up = rho_up + omega_rho * (rho_up_new - rho_up)
            rho_lo = rho_lo + omega_rho * (rho_lo_new - rho_lo)
        drho_history.append(drho)

        if outer > 0 and dgamma < tol_gamma and drho < tol_rho:
            converged = True
            break

    mach2_max = float(np.max(mvop.element_mach2(phi_ext, m_inf, gamma_air, u_inf)))

    return {
        "phi_ext": phi_ext,
        "phi": mvop.main_potential(phi_ext),
        "gamma": gamma,
        "cl_kj": 2.0 * gamma / u_inf,
        "te_jump": mvop.te_jump(phi_ext),
        "n_outer": len(drho_history),
        "converged": converged,
        "gamma_history": gamma_history,
        "drho_history": drho_history,
        "residual_norm": residual,
        "mach2_max": mach2_max,
    }
