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


def _farfield_split(mesh, alpha_deg, u_inf):
    """Split the far-field boundary faces into inflow / outflow by the sign of
    the outward flux u_inf . n_hat (design_track_b.md section 5.4 option b,
    the Lopez rectangular-domain inlet-Dirichlet / outlet-Neumann form adapted
    to pyFP3D's spherical/circular far field).

    Outward normals are oriented away from the far-field boundary's own
    geometric centre (a sphere/circle, so its centroid is its centre), which
    is robust for both the quasi-2D circle and the 3D half-sphere without
    needing the domain-interior point.

    Returns:
        faces:        (n_f, 3) far-field triangles (global node ids)
        area:         (n_f,)  triangle areas
        u_dot_n:      (n_f,)  u_inf . n_hat_out per face
        inflow_nodes: unique node ids on faces with u_dot_n < 0
    """
    faces = np.asarray(mesh.boundary_faces["farfield"], dtype=np.int64)
    p = mesh.nodes[faces]                       # (n_f, 3, 3)
    v0, v1, v2 = p[:, 0], p[:, 1], p[:, 2]
    cross = np.cross(v1 - v0, v2 - v0)
    area = 0.5 * np.linalg.norm(cross, axis=1)
    n_hat = cross / (np.linalg.norm(cross, axis=1, keepdims=True) + 1e-300)
    centroid = p.mean(axis=1)                   # (n_f, 3)
    ff_centre = np.unique(faces)
    ctr = mesh.nodes[ff_centre].mean(axis=0)
    outward = np.sign(np.einsum("fi,fi->f", n_hat, centroid - ctr))
    outward[outward == 0.0] = 1.0
    n_hat = n_hat * outward[:, None]

    a = np.deg2rad(alpha_deg)
    u_vec = np.array([u_inf * np.cos(a), u_inf * np.sin(a), 0.0])
    u_dot_n = n_hat @ u_vec
    inflow_nodes = np.unique(faces[u_dot_n < 0.0].ravel())
    return faces, area, u_dot_n, inflow_nodes


def _neumann_outlet_rhs(mesh, alpha_deg, u_inf, n_total, rho_inf=1.0):
    """Consistent P1 load vector for the Lopez outlet Neumann flux
    q = rho_inf (u_inf . n_hat) on the OUTFLOW far-field faces (option b).

    In the Galerkin weak form the boundary term is oint rho (grad phi . n) v dS;
    prescribing the freestream flux rho_inf (u_inf . n) on the outflow adds
    sum_faces rho_inf (u.n) * area/3 to each of the face's three main-DOF rows.
    (rho_inf = 1.0 in the stagnation-normalised isentropic density used here --
    density_isentropic(u_inf^2, M_inf) == 1.) The freestream potential already
    satisfies this flux exactly, so it imposes "the perturbation has zero net
    outflow flux" without a vortex correction. Inflow faces are Dirichlet
    (phi_inf), so their rows are overwritten by apply_dirichlet and any RHS
    entry there is discarded."""
    faces, area, u_dot_n, _ = _farfield_split(mesh, alpha_deg, u_inf)
    out = u_dot_n > 0.0
    b = np.zeros(n_total, dtype=np.float64)
    contrib = rho_inf * u_dot_n[out] * area[out] / 3.0
    np.add.at(b, faces[out].ravel(), np.repeat(contrib, 3))
    return b


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
    te_kutta: str = "pressure",
    farfield: str = "vortex",
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

    Far-field options (design_track_b.md section 5.4, arbitrated as the B4.5
    A/B stage). farfield="vortex" (default, option a): spherical Dirichlet
    freestream + PG vortex on the far-field MAIN DOFs, with the extracted
    Gamma(z) refreshed into the vortex each outer iteration (RHS-only).
    farfield="neumann" (option b, Lopez): inflow far-field is Dirichlet
    freestream (NO vortex), outflow is a Neumann outlet carrying the freestream
    flux rho_inf(u.n); no Gamma-into-far-field feedback (the attractive
    workflow property -- the alpha-sweep loop needs no vortex refresh).
    farfield="freestream": Dirichlet freestream on the WHOLE far field, no
    vortex -- the crudest truncation, kept for the B4.5 domain-size study
    (it is the upper bound on the truncation bias at a given domain radius).

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
    if farfield not in ("vortex", "neumann", "freestream"):
        raise ValueError(f"farfield={farfield!r} unknown")
    beta = float(np.sqrt(1.0 - m_inf**2))
    op = mvop.op

    # Option b (Neumann outlet): the outflow flux RHS is a FIXED vector
    # (freestream flux, solution-independent), assembled once.
    if farfield == "neumann":
        ff_split_nodes = _farfield_split(mesh, alpha_deg, u_inf)[3]
        ff_split_vals = freestream_phi(mesh.nodes[ff_split_nodes], alpha_deg,
                                       u_inf)
        neumann_rhs = _neumann_outlet_rhs(mesh, alpha_deg, u_inf,
                                          mvop.n_total)

    phi_ext = np.zeros(mvop.n_total, dtype=np.float64)
    phi_ext[: mvop.n_main] = freestream_phi(mesh.nodes, alpha_deg, u_inf)
    # Seed the aux DOFs with a ZERO jump (aux = main). Leaving them at 0 would
    # manufacture a huge fake jump, and the first TE-Kutta linearization reads
    # the seed: with a zero jump q_u = q_l = u_inf, so s = q_u + q_l = 2 u_inf
    # and the row starts as the classical linearized Kutta.
    cut_nodes = np.flatnonzero(mvop.cm.ext_dof_of_node >= 0)
    phi_ext[mvop.cm.ext_dof_of_node[cut_nodes]] = phi_ext[cut_nodes]
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
            te_kutta=te_kutta,
            phi_ext=phi_ext,          # re-linearizes the TE Kutta row (B4)
        )
        if farfield == "vortex":
            ff, vals = _farfield_main(mesh, alpha_deg, gamma, u_inf,
                                      vortex_center, beta)
            b = np.zeros(mvop.n_total)
        elif farfield == "neumann":
            ff, vals, b = ff_split_nodes, ff_split_vals, neumann_rhs
        else:  # "freestream": Dirichlet freestream on the whole far field
            ff = np.unique(mesh.boundary_faces["farfield"])
            vals = freestream_phi(mesh.nodes[ff], alpha_deg, u_inf)
            b = np.zeros(mvop.n_total)
        A_free, b_free, free, phi_ext = apply_dirichlet(A, b, ff, vals)
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
        "farfield": farfield,
    }
