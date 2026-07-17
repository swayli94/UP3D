"""
Level-set (Track B) solve drivers -- parallel to solve/picard.py, never
imported by the conforming path.

B2 ships the NON-LIFTING driver only: a single extended assemble + direct
solve on the multivalued operator, enough for the B2 consistency gates
(V0 freestream, V1 MMS convergence, Laplace a=0 -> cl ~ 0). The extended
matrix is nonsymmetric (the aux-row weld couples aux -> main one-way), so
by default it is solved with a sparse direct LU (scipy `spsolve`), not
CG/AMG -- for the coarse/medium meshes the direct factor is cheap and
isolates "is the assembly correct" from any preconditioner-convergence
question.

B11 (2026-07-14) lands the deferred scaling path (design_track_b.md §5.3):
every driver here takes a `precond` kwarg (None = the bit-identical direct
default; "ilu"/"amg" run preconditioned GMRES via solve/linear.solve_gmres).
**★ Measured boundary: the iterative escapes are COARSE-only** -- ILU converges
on the coarse fused matrix (434 iters, exact) but DIVERGES at 2.5D medium
lifting and factor-fails at M6 medium; AMG is built on an SPD surrogate
(`_amg_surrogate_preconditioner`) and converges only on the SPD
`continuity`/Laplace system, stalling on the `wake_ls` lifting operator
(convection-like wake-LS + TE-Kutta rows). So at medium/M6 sizes sparse-direct
is the only converging tool and the cost driver is the NUMBER of
factorizations. B13 (2026-07-14) attacks exactly that: `direct_refactor_every`
on `solve_multivalued_lifting` (the B12 lagged-LU mechanism applied to the
Picard outer loop -- refactor every k-th outer, stale-exact-LU-preconditioned
GMRES in between; default 1 = bit-identical per-outer spsolve). The transonic
driver inherits it through **kwargs.

B14 (2026-07-17) adds `precond="schur"` to `solve_multivalued_lifting` (and,
via kwargs, the transonic driver): exact Schur elimination of the aux
thin-strip block + AMG on the SPD Picard main block, NO springs
(solve/schur_ls.py) -- the structural fix for the surrogate's stall.
`solve_multivalued_laplace` keeps None|'ilu'|'amg' only: its weld-closure aux
block is the trivially invertible identity redirection and the SPD surrogate
already works there.

The lifting driver with implicit Kutta (the g1+g2 wake LS closure, no Gamma
outer loop) is B3.
"""

import time
from typing import Dict, Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from pyfp3d.constraints.dirichlet import freestream_phi, vortex_phi_2d
from pyfp3d.solve.linear import (
    apply_dirichlet,
    build_amg_preconditioner,
    build_ilu_preconditioner,
    solve_gmres,
)
from pyfp3d.solve.schur_ls import SchurReducedSystem, main_block_preconditioner
from pyfp3d.solve.timing import (
    finalize,
    new_timings,
    phase,
    snapshot,
    step_delta,
    sum_timings,
)
from pyfp3d.wake.multivalued import MultivaluedOperator


def _amg_surrogate_preconditioner(mvop, A_fused, rho_pair, free):
    """SA-AMG preconditioner for the FUSED extended system, built on an SPD
    SURROGATE (B11; design_track_b.md §5.3).

    pyamg smoothed aggregation assumes a near-SPD operator, but the fused
    level-set matrix is structurally nonsymmetric (the aux weld / wake-LS /
    TE-Kutta rows, with negative diagonals). So instead of preconditioning it
    directly we build AMG on an SPD SURROGATE

        S = op.assemble_matrix(rho_own)   [the SPD single-valued Picard block]
          + Σ_a  k_a · (e_h - e_a)(e_h - e_a)ᵀ   [aux<->host springs]

    where each aux dof a is tied to its host node h (`cm.ext_dof_of_node`) by
    an SPD spring of strength k_a = |A_fused[a,a]| (a plain block-diagonal aux
    singleton makes SA isolate each aux in its own aggregate; the spring makes
    SA AGGREGATE aux WITH its coincident host, §5.3 "把 N_ext 个辅助 DOF 当普通
    节点处理即可").

    ★ MEASURED BOUNDARY (2026-07-14): this works as a GMRES preconditioner for
    the `continuity`-closure (Laplace) system, but NOT for the `wake_ls`-closure
    LIFTING/transonic/Newton operator -- there the aux rows are the g1+g2
    wake-LS + nonlinear TE-Kutta rows (convection-like, not SPD springs), the
    surrogate cannot model them, and GMRES stalls at the restart cap (coarse
    M0.5 lifting: gamma 0.0033 vs 0.139, all outers stalled, 455 s).
    **On the lifting path use `precond="ilu"`** -- it factors the real fused
    matrix and converges (434 iters, exact). AMG stays wired for the SPD
    Laplace case and as the recorded §5.3 knob; the Núñez symmetric row
    assignment (which would restore genuine AMG applicability) is the recorded
    not-prebuilt fallback.

    Args:
        mvop: the MultivaluedOperator (supplies op + cm + n_main/n_total).
        A_fused: the fused extended matrix (its aux diagonal sets k_a).
        rho_pair: (rho_up, rho_lo) per-element densities, or None (Laplace
            limit / subcritical rho == 1).
        free: free-dof indices (into n_total), the apply_dirichlet split.

    Returns:
        M: the pyamg preconditioner LinearOperator (restricted to free dofs).
    """
    rho_own = None if rho_pair is None else mvop.own_side_field(*rho_pair)
    a = mvop.op.assemble_matrix(rho_own).tocoo()
    n_main, n_total = mvop.n_main, mvop.n_total
    host = np.flatnonzero(mvop.cm.ext_dof_of_node >= 0)
    aux = mvop.cm.ext_dof_of_node[host].astype(np.int64)
    k = np.maximum(np.abs(A_fused.diagonal()[aux]), 1e-30)
    # SPD aux<->host spring [[k,-k],[-k,k]] per pair.
    rows = np.concatenate([a.row, host, aux, host, aux])
    cols = np.concatenate([a.col, host, aux, aux, host])
    data = np.concatenate([a.data, k, k, -k, -k])
    S = sp.coo_matrix((data, (rows, cols)),
                      shape=(n_total, n_total)).tocsr()
    return build_amg_preconditioner(S[free][:, free].tocsr())[1]


def solve_multivalued_laplace(
    mvop: MultivaluedOperator,
    dirichlet_nodes: np.ndarray,
    dirichlet_values: np.ndarray,
    body_source_rhs: Optional[np.ndarray] = None,
    precond: Optional[str] = None,
    linear_rtol: float = 1e-10,
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

    if precond is None:
        x = np.atleast_1d(spla.spsolve(A_free.tocsc(), b_free))
    else:
        # B11: iterative escape from the splu wall. One-shot solve => "retry"
        # (a Laplace GMRES that stalls is a bug, not a step to accept).
        if precond == "ilu":
            M = build_ilu_preconditioner(A_free)
        elif precond == "amg":
            M = _amg_surrogate_preconditioner(mvop, A, None, free)
        else:
            raise ValueError(f"precond={precond!r} unknown (None|'ilu'|'amg')")
        x, _, _ = solve_gmres(A_free, b_free, M=M, rtol=linear_rtol,
                              on_fail="retry")
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


def farfield_aux_dofs(mesh, cm):
    """Far-field boundary nodes that carry a wake-jump aux DOF, and those aux
    DOF ids (B16).

    A wake level set has no outflow clip (cut_elements.py: the cut test only
    requires the crossing to be downstream of the TE and within the span), so
    the sheet reaches the far-field boundary and the outer nodes it crosses
    each get an aux DOF. Those aux are governed only by a near-singular wake-LS
    row on a giant outer tet; at a converged freestream Picard state they hold
    garbage (measured B16/GB16.1 on the B9 wing-body: |jump| 22-53 vs the
    physical Gamma ~ 0.06). The Picard fixed point tolerates it -- it solves
    those rows to zero garbage-and-all -- but the Newton residual reads it as
    an O(1) local inconsistency (8 neighbouring far-field MAIN rows |R| ~ 84),
    which is exactly why the committed LS Newton recipes churn on the wing-body
    (B9 recorded follow-up).

    Under farfield in ("freestream", "vortex") the far-field MAIN DOFs are
    Dirichlet; the matching move (B16 farfield_aux="pin", solve/newton_ls.py) is
    to pin these aux to the branch value their host carries. Independent of
    Gamma / Mach / geometry -- a pure structural constraint. Returns empty
    arrays when the sheet does not reach the boundary.

    Returns:
        hosts: (k,) far-field node ids carrying an aux DOF (all < cm.n_main)
        aux:   (k,) their aux DOF ids (all >= cm.n_main)
    """
    ff = np.unique(np.asarray(mesh.boundary_faces["farfield"], dtype=np.int64))
    aux = cm.ext_dof_of_node[ff]
    has = aux >= 0
    return ff[has], aux[has]


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
    gamma_farfield_fixed: Optional[float] = None,
    te_kutta: str = "pressure",
    farfield: str = "vortex",
    farfield_aux: str = "pin_gamma",
    upwind_c: float = 0.0,
    m_crit: float = 0.95,
    m_cap: float = 3.0,
    rho_floor: float = 0.05,
    damping_theta: Optional[float] = None,
    damping_scope: str = "supersonic",
    omega: float = 1.0,
    tol_residual: Optional[float] = None,
    phi_init: Optional[np.ndarray] = None,
    tip_taper: Optional[np.ndarray] = None,
    precond: Optional[str] = None,
    linear_rtol: float = 1e-10,
    gmres_restart: int = 60,
    gmres_maxiter: int = 50,
    amg_rebuild_every: int = 5,
    direct_refactor_every: int = 1,
    direct_reuse_rtol: float = 1e-10,
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

    Far-field options (design_track_b.md section 5.4, arbitrated as the B5
    A/B stage). farfield="vortex" (default, option a): spherical Dirichlet
    freestream + PG vortex on the far-field MAIN DOFs, with the extracted
    Gamma(z) refreshed into the vortex each outer iteration (RHS-only).
    farfield="neumann" (option b, Lopez): inflow far-field is Dirichlet
    freestream (NO vortex), outflow is a Neumann outlet carrying the freestream
    flux rho_inf(u.n); no Gamma-into-far-field feedback (the attractive
    workflow property -- the alpha-sweep loop needs no vortex refresh).
    farfield="freestream": Dirichlet freestream on the WHOLE far field, no
    vortex -- the crudest truncation, kept for the B5 domain-size study
    (it is the upper bound on the truncation bias at a given domain radius).

    `farfield_aux` (B17): how the far-field-BOUNDARY aux DOFs are treated when
    the wake sheet reaches the outflow (no LS outflow clip; farfield_aux_dofs).
    Acts ONLY on farfield="freestream"; inert on vortex/neumann (so the default
    leaves every committed 2.5D vortex/neumann run bit-identical). "legacy"
    leaves them FREE -- the pre-B17 Picard behaviour (the B9/B16 freestream demos
    pass it explicitly to keep their committed numbers). The Picard fixed point
    absorbs the near-singular outer wake-LS rows (it solves them to zero
    garbage-and-all, |jump| 22-53 at the converged state) which is why Picard
    "works" where the Newton residual reads the same rows as an O(1)
    inconsistency (B16/GB16.1).
    "pin" adds the far-field aux to the Dirichlet set at the host node's own
    single-valued phi_inf (jump -> 0 on the outflow ring), the SAME pin B16
    applies on the Newton path (solve/newton_ls.py) -- so a Picard-vs-Newton A/B
    is on ONE far-field discretization. ★ B17 MEASURED that jump -> 0 REMOVES the
    wake circulation the outflow physically carries: the medium wing-body cl_p
    drops from 0.2165 (legacy/conforming) to 0.1691, a 22% resolution-dependent
    error (the coarse "match" was a coincidence -- jump=0 there cancelled the
    coarse legacy's outer-tet garbage). "pin" is kept as the B16 reproduction.
    "pin_gamma" (default) is the FIX: aux = host phi_inf - side*gamma, i.e.
    jump -> gamma, carrying the circulation OUT (refreshed with the live gamma
    each outer). It cures the same near-singular outer wake-LS rows (identical
    Dirichlet conditioning) while keeping cl_p monotone-convergent to conforming
    (coarse 0.2087, medium 0.2117). All three act on farfield="freestream" only;
    on vortex/neumann they are inert (neumann's outer aux ARE
    wake-LS-constrained, and farfield="vortex" already pins aux to jump=gamma on
    the Newton path -- see solve/newton_ls.py).

    ★ Compressibility is carried by the BULK density, NOT the far-field vortex:
    measured on the medium NACA case, PG-scaling the far-field vortex (beta <
    1) leaves gamma unchanged (0.1086 -> 0.1086 -- the soft Kutta does not
    propagate the outer-vortex stretch), while the isentropic bulk density
    alone raises it 0.1086 -> 0.1256 (the physical compressible lift rise, ~93%
    of the conforming M0.5 gamma, the SAME convergence ratio as incompressible).
    The cut-strip per-side density limit-cycles, so it is under-relaxed
    (omega_rho); at omega_rho = 1 the loop diverges after ~80 iterations.

    Transonic (B6, upwind_c > 0): the artificial density runs PER SIDE on the
    cut elements with a same-side-restricted upstream walk
    (`MultivaluedOperator.element_rho_tilde`, design_track_b.md section
    5.2/D10), the isentropic q^2 limiter (m_cap) is applied per side, and
    `damping_theta` adds theta*diag damping LOCALIZED to the nu > 0 rows
    (damping_scope="supersonic"; scope "fluid" -- the conforming P4 form on
    every fluid row -- is kept only as the measured NEGATIVE result: it
    throttles the implicit-Kutta circulation mode, see the inline note).
    `gamma_farfield_fixed` freezes the far-field VORTEX strength while the TE
    jump stays free -- the fold-zone stabilizer used by
    `solve_multivalued_transonic(farfield_lag="level")`. `tol_residual` adds a
    nonlinear-residual convergence route (None = the B3 lag-based criterion
    only, bit-identical). `phi_init` warm-starts the FULL extended state
    (continuation). All defaults keep the B3/B4 subsonic behavior bitwise.

    `precond` (B11; design_track_b.md §5.3) selects the inner linear solver.
    None (default) is the pre-B11 sparse-direct `spsolve`, bit-identical.
    **Use "ilu"** for the iterative escape from the M6-fine splu wall (roadmap
    "no precond option" caveat) -- preconditioned GMRES (solve/linear.solve_gmres)
    on the fused nonsymmetric matrix, converging (434 iters coarse). "amg"
    (the SPD surrogate `_amg_surrogate_preconditioner`) is wired but STALLS on
    this `wake_ls` lifting operator (measured; its aux rows are convection-like,
    not SPD) -- it is for the SPD Laplace case only. Each outer warm-starts
    GMRES from the previous free-dof iterate (AMG surrogate rebuilt every
    `amg_rebuild_every` outers). `linear_rtol=1e-10` is tight enough that the
    fixed point is unperturbed at gate precision (the lag criteria live at
    tol_gamma/tol_rho=1e-6, tol_residual~1e-7); inexact GMRES steps
    (`on_fail="return"`) only perturb the Picard path and are counted as
    stalls, never raised.

    "schur" (B14; design_track_b.md §5.3, roadmap §B14) is the STRUCTURAL
    iterative escape: eliminate the aux thin-strip block exactly per outer
    (`lu_aa = splu(J_aa)`, n_ext-sized, milliseconds) and run GMRES on the
    reduced main-free operator preconditioned by AMG on the SPD single-valued
    Picard block -- NO springs, so the B11 surrogate's jump==0 bias is
    structurally absent and the circulation mode survives. A stalled reduced
    GMRES falls back to a full fused spsolve in the same outer
    (`n_schur_fallback`), so robustness never drops below per-outer-direct.

    `direct_refactor_every` (B13): lagged-LU direct-reuse on the `precond=None`
    path -- the B12 mechanism applied to the Picard OUTER loop, which after B12
    is the M6-medium cost driver (one 17.5 s spsolve per outer; B11 measured
    ILU diverging at 2.5D medium and factor-failing at M6 medium, so at these
    sizes sparse-direct is the only converging tool and the cost is the NUMBER
    of factorizations). 1 (default) = the bit-identical per-outer `spsolve`.
    > 1 = refactor the LU every k-th outer and drive the outers in between with
    GMRES on the FRESH matrix preconditioned by the stale (exact) LU, converged
    to `direct_reuse_rtol` and warm-started from the previous iterate; a reuse
    step whose GMRES fails falls back to refactor + exact solve in the same
    outer (robustness never below per-outer-direct; early outers with large
    density moves are expected to trigger it). Ignored when precond is
    "ilu"/"amg"/"schur". `direct_reuse_rtol` defaults to 1e-10, NOT B12's 1e-8: a
    Picard fixed point is pinned only by its lag tolerances (1e-6), so an
    inexact reuse step SHIFTS the stopping point (measured |dgamma| 8e-8 at
    rtol 1e-8), whereas Newton's terminus is pinned by tol_residual regardless
    -- 1e-10 restores <1e-8 gamma agreement for ~1-2 extra Krylov iters on a
    near-exact preconditioner.

    Returns:
        dict: phi_ext, phi (main), gamma (extracted TE circulation),
        cl_kj (= 2 gamma / (u_inf c), c = 1), te_jump, n_outer, converged,
        gamma_history, drho_history, residual_history, residual_norm,
        mach2_max, nu_max, n_nu_active, n_limited, n_floored, rho_tilde,
        precond, n_gmres_total, n_gmres_stalled (B11 monitors),
        n_refactor (B13), n_schur_fallback (B14).
    """
    if precond not in (None, "ilu", "amg", "schur"):
        raise ValueError(
            f"precond={precond!r} unknown (None|'ilu'|'amg'|'schur')")
    if not 0.0 <= m_inf < 1.0:
        raise ValueError(f"needs 0 <= M_inf < 1, got {m_inf}")
    if farfield not in ("vortex", "neumann", "freestream"):
        raise ValueError(f"farfield={farfield!r} unknown")
    if farfield_aux not in ("legacy", "pin", "pin_gamma"):
        raise ValueError(f"farfield_aux={farfield_aux!r} unknown")
    # B17: the pin acts only on farfield="freestream" (the B9/B16 wing-body BC).
    # It is INERT on vortex/neumann -- exactly like B16's neumann inertness --
    # so the default "pin_gamma" leaves EVERY committed vortex/neumann Picard run
    # (all the 2.5D NACA cases) bit-identical. A freestream Picard that wants the
    # pre-B17 free-aux behaviour must pass farfield_aux="legacy" explicitly (the
    # B9/B16 demos do, to keep their committed legacy numbers reproducible).
    pin_aux = farfield_aux in ("pin", "pin_gamma") and farfield == "freestream"
    pin_gamma = farfield_aux == "pin_gamma"
    if len(mvop.cm.te_nodes) == 0:
        # A wake level-set that matches NO trailing-edge wall node carries no
        # Kutta condition at all: Gamma is unpinned, the TE aux rows are empty,
        # and the solve silently produces NaN (measured on ONERA M6 medium: a
        # hand-rolled TE polyline off by ~2e-4 in x -> 0 TE nodes -> 340k
        # limited cells, gamma = mean([]) = NaN). Fail loudly instead. Build the
        # polyline from the AUTHORITATIVE geometry (e.g. meshgen.wing3d.x_te),
        # whose endpoints are exact wall nodes.
        raise ValueError(
            "the wake level-set matched 0 TE nodes on this mesh: the TE "
            "polyline does not lie on the wall. Without TE nodes there is no "
            "Kutta condition and the solution is meaningless. Check that the "
            "polyline comes from the mesh's own geometry (e.g. "
            "pyfp3d.meshgen.wing3d.x_te / B_SEMI) and that wall_nodes was "
            "passed to CutElementMap.")
    t_wall0 = time.perf_counter()
    timings = new_timings()
    beta = float(np.sqrt(1.0 - m_inf**2))
    op = mvop.op
    # B6: artificial-density upwinding, per side on the cut elements
    # (design_track_b.md §5.2/D10). upwind_c = 0 keeps the B3/B4 subcritical
    # path bit-identical (no upwind operator is even built).
    use_upwind = upwind_c > 0.0 and m_inf > 0.0

    # Option b (Neumann outlet): the outflow flux RHS is a FIXED vector
    # (freestream flux, solution-independent), assembled once.
    if farfield == "neumann":
        ff_split_nodes = _farfield_split(mesh, alpha_deg, u_inf)[3]
        ff_split_vals = freestream_phi(mesh.nodes[ff_split_nodes], alpha_deg,
                                       u_inf)
        neumann_rhs = _neumann_outlet_rhs(mesh, alpha_deg, u_inf,
                                          mvop.n_total)

    phi_ext = np.zeros(mvop.n_total, dtype=np.float64)
    if phi_init is not None:
        # Continuation warm start (B6): the previous Mach level's full
        # extended state, aux DOFs (i.e. the wake jump) included.
        phi_ext[:] = np.asarray(phi_init, dtype=np.float64)
    else:
        phi_ext[: mvop.n_main] = freestream_phi(mesh.nodes, alpha_deg, u_inf)
        # Seed the aux DOFs with a ZERO jump (aux = main). Leaving them at 0
        # would manufacture a huge fake jump, and the first TE-Kutta
        # linearization reads the seed: with a zero jump q_u = q_l = u_inf, so
        # s = q_u + q_l = 2 u_inf and the row starts as the classical
        # linearized Kutta.
        cut_nodes = np.flatnonzero(mvop.cm.ext_dof_of_node >= 0)
        phi_ext[mvop.cm.ext_dof_of_node[cut_nodes]] = phi_ext[cut_nodes]
    gamma = float(gamma_init)
    gamma_history = [gamma]
    drho_history = []
    residual_history = []
    step_records = []
    rho_up = rho_lo = None
    residual = np.inf
    dgamma = np.inf
    # B11 iterative-solver state (inert when precond is None).
    M_pre = None
    n_gmres_total = 0
    n_gmres_stalled = 0
    # B13 lagged-LU state (inert when direct_refactor_every <= 1).
    lu_direct = None
    lu_age = 0
    n_refactor = 0
    # B14 Schur state (inert unless precond == "schur").
    n_schur_fallback = 0
    converged = False
    # ---- B6 damping: LOCALIZED to the supersonic zone --------------------
    #
    # ★ The conforming P4 stabilizer (theta*diag(A_free) on the whole reduced
    # system, solve/picard.py) does NOT transfer to the level-set path, and the
    # reason is structural, not a tuning issue (measured 2026-07-12; see the
    # `damping_scope="fluid"` branch, kept so the negative result stays
    # reproducible). The error of a damped Picard step contracts as
    # (A+D)^{-1}D, whose eigenvalue for a mode of stiffness lambda is
    # theta*d/(lambda + theta*d): near 0 for the stiff LOCAL modes (the shock)
    # but near 1 for the SMOOTH GLOBAL ones. Diagonal damping is a Jacobi
    # smoother -- it barely touches smooth modes.
    #
    # On the CONFORMING path that is harmless, because the smooth global mode
    # that matters -- the circulation -- is NOT in the damped matrix at all:
    # Gamma is an OUTER secant unknown. On the B path the implicit Kutta makes
    # Gamma a SOLUTION MODE, so the same damping throttles precisely the thing
    # the whole track exists to compute: measured on NACA coarse M0.60, Gamma
    # crawls 0.0005 -> 0.017 in 160 outers (converged value 0.0968) and the
    # residual never falls, while the undamped run converges in 35 outers.
    #
    # The fix follows from the mechanism: damp ONLY where the instability lives.
    # The rows carrying nodes of supersonic (nu > 0) elements get theta*diag;
    # everything else is untouched, so the global circulation mode converges at
    # its undamped rate while the shock's local modes are smoothed. Measured:
    # the M > 2 blow-up cells sit at x 0.19-0.87, y 0.02-0.23 -- the pocket ABOVE
    # the airfoil, with ZERO cells on the wake sheet -- so the instability is the
    # ordinary transonic one, and localizing the cure to it costs nothing.
    #
    # Aux rows are never damped: the wake-LS/TE-Kutta rows are CONSTRAINT rows
    # whose diagonal is negative by construction (the LS row is +k on the phi_u
    # columns and -k on phi_l, and a "+"-side node's own aux dof IS its phi_l
    # column), so theta*diag there would be ANTI-damping.
    if damping_scope not in ("supersonic", "fluid"):
        raise ValueError(f"damping_scope={damping_scope!r} unknown")
    is_main = np.zeros(mvop.n_total, dtype=bool)
    is_main[: mvop.n_main] = True
    elements = np.asarray(op.elements, dtype=np.int64)
    dampable = is_main.copy()      # refreshed each outer when scope=supersonic

    # B17: far-field aux pin (freestream only). Precompute the aux DOF ids and
    # their host phi_inf once -- structural, Gamma/Mach-independent -- and assert
    # they are disjoint from the TE (Kutta) aux so the pin never overwrites a
    # Kutta row. Same guard as solve/newton_ls.py.
    if pin_aux:
        ff_aux_hosts, ff_aux_dofs = farfield_aux_dofs(mesh, mvop.cm)
        ff_aux_vals = freestream_phi(mesh.nodes[ff_aux_hosts], alpha_deg, u_inf)
        te_aux = mvop.cm.ext_dof_of_node[mvop.cm.te_nodes]
        te_aux = te_aux[te_aux >= 0]
        if np.intersect1d(ff_aux_dofs, te_aux).size:
            raise RuntimeError(
                "a TE aux DOF lies on the far-field boundary; pinning it would "
                "overwrite a Kutta row -- the wake reaches the outflow at a TE "
                "station, which the wing-body geometry should preclude")

    for outer in range(n_outer_max):
        t_phase0 = snapshot(timings)
        n_gmres_step = n_gmres_total
        n_refactor_step = n_refactor
        t_assembly0 = time.perf_counter()
        A = mvop.assemble_matrix(
            rho_tilde=(None if rho_up is None else (rho_up, rho_lo)),
            closure="wake_ls",
            te_kutta=te_kutta,
            phi_ext=phi_ext,          # re-linearizes the TE Kutta row (B4)
            tip_taper=tip_taper,      # B8 tip-edge desingularization (None=no-op)
        )
        if farfield == "vortex":
            # B6: near the FP fold the live Gamma->far-field-vortex feedback
            # has loop gain > 1 (measured runaway; see solve_multivalued_
            # transonic). gamma_farfield_fixed freezes the VORTEX strength
            # (the O(Gamma/R) outer-boundary consistency only) while the
            # implicit Kutta still sets the TE jump freely -- the transonic
            # driver refreshes it between Mach levels / polish rounds.
            g_ff = (gamma if gamma_farfield_fixed is None
                    else gamma_farfield_fixed)
            ff, vals = _farfield_main(mesh, alpha_deg, g_ff, u_inf,
                                      vortex_center, beta)
            b = np.zeros(mvop.n_total)
        elif farfield == "neumann":
            ff, vals, b = ff_split_nodes, ff_split_vals, neumann_rhs
        else:  # "freestream": Dirichlet freestream on the whole far field
            ff = np.unique(mesh.boundary_faces["farfield"])
            vals = freestream_phi(mesh.nodes[ff], alpha_deg, u_inf)
            if pin_aux:
                # B17: also pin the far-field-boundary aux. ff are all main
                # nodes and ff_aux_dofs are all >= n_main, so the concatenated
                # set has no duplicate DOF. Two branches:
                #   "pin"       -> aux = host phi_inf  (jump -> 0, the B16 rule)
                #   "pin_gamma" -> aux = host phi_inf - side*gamma (jump -> gamma,
                #                  carrying the wake circulation OUT through the
                #                  outflow -- refreshed with the live gamma each
                #                  outer, the physically correct ring value).
                if pin_gamma:
                    side = mvop.cm.node_side[ff_aux_hosts]
                    aux_vals = ff_aux_vals - side * gamma
                else:
                    aux_vals = ff_aux_vals
                ff = np.concatenate([ff, ff_aux_dofs])
                vals = np.concatenate([vals, aux_vals])
            b = np.zeros(mvop.n_total)
        phi_prev = phi_ext
        A_free, b_free, free, phi_ext = apply_dirichlet(A, b, ff, vals)
        x_prev = phi_prev[free]

        # Nonlinear (fixed-point) residual of the PREVIOUS iterate against the
        # operator freshly assembled at that iterate's density + TE
        # linearization -- the honest transonic convergence monitor. The
        # post-solve ||b - A x|| is only the direct solver's ~1e-13 linear
        # residual and says nothing about the nonlinearity.
        residual = float(np.max(np.abs(b_free - A_free @ x_prev)))
        residual_history.append(residual)

        A_solve, b_solve = A_free, b_free
        if damping_theta is not None:
            import scipy.sparse as sp_

            d = damping_theta * A_free.diagonal()
            d[~dampable[free]] = 0.0
            np.clip(d, 0.0, None, out=d)
            A_solve = (A_free + sp_.diags(d)).tocsr()
            b_solve = b_free + d * x_prev
        # The residual above is one sparse matvec on an already-assembled
        # operator -- it rides along with the assembly it is measured against.
        timings["assembly"] += time.perf_counter() - t_assembly0

        t_lin0 = time.perf_counter()
        if precond is None:
            if direct_refactor_every <= 1:
                x = np.atleast_1d(spla.spsolve(A_solve.tocsc(), b_solve))
            else:
                # B13 lagged-LU (the B12 mechanism on the Picard outer loop):
                # refactor the LU only every k-th outer; in between, GMRES on
                # the FRESH matrix preconditioned by the stale (exact) LU,
                # warm-started from the previous iterate. The under-relaxed
                # density (omega_rho) keeps the matrix drift small, so the
                # stale LU stays near-exact; a failed reuse step refactors in
                # the same outer (never below per-outer-direct robustness).
                x = None
                if lu_direct is not None and lu_age < direct_refactor_every:
                    n = free.size
                    A_op = spla.LinearOperator((n, n), matvec=A_solve.dot)
                    M_lag = spla.LinearOperator((n, n), matvec=lu_direct.solve)
                    x_try, n_it, info = solve_gmres(
                        A_op, b_solve, M=M_lag, rtol=direct_reuse_rtol,
                        x0=x_prev, restart=gmres_restart, maxiter=2,
                        on_fail="return")
                    n_gmres_total += n_it
                    if info == 0:
                        x = np.atleast_1d(x_try)
                    else:
                        n_gmres_stalled += 1   # stale LU exhausted: refactor
                if x is None:
                    t_lu = time.perf_counter()
                    lu_direct = spla.splu(A_solve.tocsc())
                    dt_lu = time.perf_counter() - t_lu
                    timings["precond"] += dt_lu
                    t_lin0 += dt_lu        # the factorization is not a solve
                    lu_age = 0
                    n_refactor += 1
                    x = np.atleast_1d(lu_direct.solve(b_solve))
                lu_age += 1
        else:
            # B11: preconditioned GMRES on the fused matrix, warm-started from
            # the previous outer's free-dof iterate; an inexact step
            # (on_fail="return") just perturbs the fixed-point path.
            # B14 "schur": exact aux-block elimination per outer + AMG on the
            # SPD Picard main block; the Schur split is built on A_solve --
            # the theta*diag damping is main-rows-only (dampable is a main
            # mask), so J_aa is undamped, which is correct (aux diagonals are
            # negative; damping there would be anti-damping).
            t_pre = time.perf_counter()
            schur = None
            if precond == "ilu":
                M_pre = build_ilu_preconditioner(A_solve)
            elif precond == "schur":
                schur = SchurReducedSystem(A_solve, free, mvop.n_main,
                                           n_aux_expected=mvop.n_ext)
                if outer % amg_rebuild_every == 0 or M_pre is None:
                    M_pre = main_block_preconditioner(
                        mvop,
                        (None if rho_up is None else (rho_up, rho_lo)),
                        schur.main_free)
            elif outer % amg_rebuild_every == 0 or M_pre is None:
                M_pre = _amg_surrogate_preconditioner(
                    mvop, A, (None if rho_up is None else (rho_up, rho_lo)),
                    free)
            dt_pre = time.perf_counter() - t_pre
            timings["precond"] += dt_pre
            t_lin0 += dt_pre
            if precond == "schur":
                x, n_it, info = schur.solve(
                    b_solve, M=M_pre, rtol=linear_rtol,
                    x0_main=x_prev[:schur.n_mf],
                    restart=gmres_restart, maxiter=gmres_maxiter)
                n_gmres_total += n_it
                if info != 0:
                    # Safety net (the lagged-LU pattern): a stalled reduced
                    # GMRES never leaves an inexact Picard step behind -- the
                    # fixed point MOVES under inexact steps (B13, |dgamma|
                    # 8e-8 at rtol 1e-8) -- refactor-and-solve the full fused
                    # system in the same outer instead.
                    n_gmres_stalled += 1
                    n_schur_fallback += 1
                    x = np.atleast_1d(spla.spsolve(A_solve.tocsc(), b_solve))
                x = np.atleast_1d(x)
            else:
                x, n_it, info = solve_gmres(
                    A_solve, b_solve, M=M_pre, rtol=linear_rtol, x0=x_prev,
                    restart=gmres_restart, maxiter=gmres_maxiter,
                    on_fail="return")
                x = np.atleast_1d(x)
                n_gmres_total += n_it
                n_gmres_stalled += int(info != 0)
        timings["linsolve"] += time.perf_counter() - t_lin0
        if omega != 1.0:
            x = x_prev + omega * (x - x_prev)
        phi_ext[free] = x

        # Gamma refresh (implicit Kutta -> extracted TE jump), relaxed.
        with phase(timings, "kutta"):
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
        with phase(timings, "residual"):
            if use_upwind:
                rho_up_new, rho_lo_new = mvop.element_rho_tilde(
                    phi_ext, m_inf, upwind_c, m_crit, gamma_air, u_inf,
                    m_cap, rho_floor
                )
                if damping_theta is not None and damping_scope == "supersonic":
                    # Rows to damp = main DOFs of the nodes of the nu > 0
                    # elements (the supersonic pocket + the shock), refreshed
                    # every outer so the damped set tracks the moving shock.
                    hot = np.zeros(mvop.n_total, dtype=bool)
                    hot[np.unique(elements[mvop.nu_active_elements()])] = True
                    dampable = hot & is_main
            else:
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

        step_records.append({
            "i": outer,
            "residual": residual,
            "drho": drho,
            "dgamma": dgamma,
            "gamma_mean": gamma,
            "gamma_root": gamma,
            "n_lin_solves": 1,
            "n_lin_iters": n_gmres_total - n_gmres_step,
            "n_refactor": n_refactor - n_refactor_step,
            "n_limited": int(mvop.n_limited) if use_upwind else 0,
            "n_floored": int(mvop.n_floored) if use_upwind else 0,
            "wall_cum_s": time.perf_counter() - t_wall0,
            **step_delta(timings, t_phase0),
        })

        # Two convergence routes: (a) the B3 fixed-point lags (dgamma, drho)
        # both settle; (b) the honest NONLINEAR residual of the previous
        # iterate is already below tol_residual -- the criterion that matters
        # under B6's localized damping, where the nu > 0 damped set churns by
        # a few elements per outer and keeps drho in a ~1e-6 micro limit cycle
        # long after the residual has fallen to 1e-8 (measured M0.70 coarse).
        # tol_residual=None (the B3/B4 default) keeps route (a) alone,
        # bit-identical.
        res_converged = tol_residual is not None and residual < tol_residual
        if outer > 0 and ((dgamma < tol_gamma and drho < tol_rho)
                          or res_converged):
            converged = True
            break

    with phase(timings, "residual"):
        mach2_max = float(
            np.max(mvop.element_mach2(phi_ext, m_inf, gamma_air, u_inf)))
    finalize(timings, time.perf_counter() - t_wall0)

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
        "residual_history": residual_history,
        "residual_norm": residual,
        "step_records": step_records,
        "timings": timings,
        "mach2_max": mach2_max,
        "farfield": farfield,
        "nu_max": mvop.nu_max if use_upwind else 0.0,
        "n_nu_active": mvop.n_nu_active if use_upwind else 0,
        "n_limited": mvop.n_limited if use_upwind else 0,
        "n_floored": mvop.n_floored if use_upwind else 0,
        "rho_tilde": (rho_up, rho_lo),
        "precond": precond,
        "n_gmres_total": n_gmres_total,
        "n_gmres_stalled": n_gmres_stalled,
        "n_refactor": n_refactor,
        "n_schur_fallback": n_schur_fallback,
    }


# B6 transonic recipe: the conforming TRANSONIC_DEFAULTS (solve/continuation.py)
# carried over where they mean the same thing -- the artificial density
# (C = 1.5, M_crit = 0.95) and theta = 0.2, but the damping is LOCALIZED to
# the supersonic zone (damping_scope="supersonic" -- the whole-field P4 form
# throttles the implicit-Kutta circulation mode, design_track_b.md section
# 10.2). omega_rho = 0.5 is B3's: the per-side cut-strip density limit-cycles
# and must be under-relaxed (already needed subsonically).
#
# ★ Strong-shock cases (fold neighborhood, e.g. NACA coarse M >= 0.775) must
# ALSO pass farfield="neumann": the live option-a Gamma -> vortex feedback has
# loop gain > 1 there and runs away monotonically through the solution, and
# the per-level lagged variant's outer map has no fixed point (section 10.3).
# The vortex default is kept for the subsonic/weak-transonic range (B5's
# arbitrated verdict; more accurate on the compact 15c domain).
B_TRANSONIC_DEFAULTS = dict(
    upwind_c=1.5,
    m_crit=0.95,
    damping_theta=0.2,
    omega_rho=0.5,
    m_start=0.60,
    dm=0.05,
)


def solve_multivalued_transonic(
    mvop: MultivaluedOperator,
    mesh,
    m_target: float,
    alpha_deg: float = 1.25,
    u_inf: float = 1.0,
    gamma_air: float = 1.4,
    m_start: float = B_TRANSONIC_DEFAULTS["m_start"],
    dm: float = B_TRANSONIC_DEFAULTS["dm"],
    upwind_c: float = B_TRANSONIC_DEFAULTS["upwind_c"],
    m_crit: float = B_TRANSONIC_DEFAULTS["m_crit"],
    damping_theta: float = B_TRANSONIC_DEFAULTS["damping_theta"],
    omega_rho: float = B_TRANSONIC_DEFAULTS["omega_rho"],
    n_outer_seed: int = 120,
    n_outer_level: int = 400,
    tol_rho: float = 1e-6,
    tol_gamma: float = 1e-6,
    tol_residual: Optional[float] = 1e-7,
    farfield_lag: str = "live",
    n_ff_polish: int = 0,
    n_outer_polish: int = 200,
    **kwargs,
) -> Dict[str, object]:
    """Mach-continuation transonic solve on the level-set path (Track B, B6).

    ★ The continuation is MUCH simpler than the conforming one
    (`solve/continuation.py::solve_transonic_lifting`) and that is the Track B
    payoff made concrete: with the implicit Kutta there is no Gamma secant and
    no per-station Kutta-closure loop to keep alive across Mach levels, so a
    level is just a warm-started Picard solve. The P5 failure mode that cost a
    whole investigation (a SINGLE station's secant not closing at the top Mach
    level, st133 32% under-circulated) cannot occur here -- Gamma is not an
    unknown of an outer loop, it is read off the converged TE jump.

    The ramp itself is kept (roadmap G10.3 verdict: KEEP the Mach ramp -- a
    cold transonic solve transits clamped states and, in the fold zone,
    diverges outright). Each level warm-starts from the previous level's FULL
    extended state (aux DOFs = the wake jump included) and its Gamma.

    Args:
        m_target: the target M_inf
        m_start, dm: the ascending Mach schedule (m_start, ..., m_target)
        upwind_c, m_crit: artificial density (design.md §3)
        damping_theta: LOCAL theta*diag(A) damping on the fluid rows (P4's fix)
        n_outer_seed: outer budget at the first (subcritical) level
        n_outer_level: outer budget at each subsequent level
        farfield_lag: "live" keeps the per-outer Gamma -> far-field-vortex
            refresh (the B3/B5 loop); "level" freezes the far-field vortex at
            the PREVIOUS level's emergent Gamma for the whole level (the
            implicit Kutta still sets the TE jump freely) -- ★ near the FP
            fold the live loop has gain > 1 and runs away monotonically
            THROUGH the solution (measured, NACA coarse M0.80: Gamma climbs
            0.14 -> 0.37 past both the conforming-Picard 0.18 and the Newton
            0.23 at flat residual ~5e-5, then blows up; under-relaxation
            cannot fix a monotone gain > 1 loop -- 1 + omega(lambda-1) > 1
            for every omega > 0). The stale-vortex bias is only the
            O(dGamma/R) outer-boundary term (B5), removed by the polish
            below. Ignored for farfield="neumann"/"freestream" (no vortex).
        n_ff_polish: with farfield_lag="level", rounds of far-field polish at
            the TARGET Mach -- refresh the vortex to the newest emergent
            Gamma, re-converge (the P5 fixed-Gamma polish pattern: each round
            is one step of the outer map at the O(1/R) coupling strength,
            contractive even where the live per-outer loop is not).
        n_outer_polish: outer budget per polish round.

    The B11 inner-solver knobs (`precond`, `linear_rtol`, `gmres_restart`,
    `gmres_maxiter`, `amg_rebuild_every`) ride `**kwargs` into every level's
    `solve_multivalued_lifting`, so `precond="amg"` here escapes the splu wall
    over the whole Mach ramp; the GMRES warm start compounds with the
    continuation `phi_init` warm start (each level starts near-converged).

    Returns:
        The final level's `solve_multivalued_lifting` dict, plus `levels`
        (per-level m_inf / gamma / mach2_max / n_outer / converged / monitors,
        polish rounds appended with m_inf=m_target) and `mach_schedule`.
    """
    from pyfp3d.solve.continuation import mach_schedule

    if farfield_lag not in ("live", "level"):
        raise ValueError(f"farfield_lag={farfield_lag!r} unknown")
    lagged = (farfield_lag == "level"
              and kwargs.get("farfield", "vortex") == "vortex")
    schedule = mach_schedule(m_target, m_start=m_start, dm=dm)
    phi_ext = None
    gamma = 0.0
    levels = []
    res = None
    t_wall0 = time.perf_counter()

    def _run(m, budget, g_ff):
        t0 = time.perf_counter()
        r = solve_multivalued_lifting(
            mvop, mesh, m, alpha_deg=alpha_deg, u_inf=u_inf,
            gamma_air=gamma_air,
            upwind_c=upwind_c, m_crit=m_crit,
            damping_theta=damping_theta, omega_rho=omega_rho,
            n_outer_max=budget,
            tol_rho=tol_rho, tol_gamma=tol_gamma, tol_residual=tol_residual,
            gamma_init=gamma, phi_init=phi_ext,
            gamma_farfield_fixed=g_ff,
            **kwargs,
        )
        r["_wall_s"] = time.perf_counter() - t0
        return r

    def _record(m, r):
        # A1: wall_s, timings, residual_history and step_records were all
        # dropped here before -- the ramp reported only a residual_norm.
        levels.append({
            "m_inf": float(m),
            "gamma": float(r["gamma"]),
            "cl_kj": float(r["cl_kj"]),
            "mach_max": float(np.sqrt(r["mach2_max"])),
            "n_outer": int(r["n_outer"]),
            "converged": bool(r["converged"]),
            "n_limited": int(r["n_limited"]),
            "n_floored": int(r["n_floored"]),
            "residual_norm": float(r["residual_norm"]),
            "wall_s": float(r["_wall_s"]),
            "timings": r["timings"],
            "residual_history": list(r["residual_history"]),
            "step_records": r["step_records"],
            "n_lin_iters": int(r["n_gmres_total"]),
            "n_refactor": int(r["n_refactor"]),
        })

    for k, m in enumerate(schedule):
        # Level 0 always runs the LIVE vortex even in lagged mode: it is
        # subcritical (the B3-proven stable regime) and a frozen Gamma_ff = 0
        # first level briefly transits a limited state (measured on medium).
        res = _run(m, n_outer_seed if k == 0 else n_outer_level,
                   gamma if (lagged and k > 0) else None)
        phi_ext = res["phi_ext"]
        gamma = res["gamma"]
        _record(m, res)

    if lagged:
        for _ in range(n_ff_polish):
            res = _run(m_target, n_outer_polish, gamma)
            phi_ext = res["phi_ext"]
            gamma = res["gamma"]
            _record(m_target, res)

    res = dict(res)
    res.pop("_wall_s", None)
    res["levels"] = levels
    res["mach_schedule"] = [float(m) for m in schedule]
    res["timings_total"] = sum_timings([lv["timings"] for lv in levels])
    res["wall_s"] = time.perf_counter() - t_wall0
    return res
