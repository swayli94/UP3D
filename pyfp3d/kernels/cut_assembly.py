"""
Multivalued (CutFEM-style) FE assembly for the level-set embedded wake
(Track B, B2). design_track_b.md sections 2.1/2.5 (Lopez dissertation
eqs. 3.33-3.34, section 3.5.5) / D6.

The key simplification the dissertation proves (section 3.5.5, "full
integration"): a cut element needs NO geometric sub-cell integration. The
same P1 element matrix K_e = V_e rho_e (B_e B_e^T) is assembled TWICE, once
for the upper copy (DOF vector `dofs_upper`) and once for the lower copy
(`dofs_lower`); each mass-conservation row is taken from the copy on that
node's OWN side (Lopez section 3.5, "main DOF rows = mass conservation from
the node's-side copy"). This module expresses that as a sparse CORRECTION
to the ordinary single-valued matrix, which is what makes it cheap and
obviously consistent:

    multivalued row a (side sigma_a) uses columns dofs_side(sigma_a):
        col b -> main(b)  if sigma_b == sigma_a
                 aux(b)   if sigma_b != sigma_a
    single-valued row a always uses col b -> main(b).

So the ONLY difference from the ordinary assembly is that, on a cut
element, the entries K[a,b] whose two nodes lie on OPPOSITE sides move their
column from main(b) to aux(b). Everything else -- non-cut elements, and
same-side entries of cut elements -- is byte-identical to the single-valued
operator. `multivalued_redirection_coo` returns exactly that move as COO
triplets (subtract at main(b), add at aux(b)); `wake/multivalued.py` adds
them to `PicardOperator.assemble_matrix()`'s output and appends the aux-row
closure.

Consistency (the B2 gate rationale): weld every aux DOF to its main DOF
(aux_k = main_j, the B2 continuity closure) and the moved entries fold
straight back to their main column -- the extended system reduces EXACTLY
to the single-valued one. So freestream / MMS / a=0 must reproduce the
single-valued solution to machine precision; any bug in the side marking,
the aux numbering or the row/col mapping breaks that immediately. The
physical wake jump ([phi] != 0) enters only in B3, when the aux rows carry
the g1+g2 wake least-squares condition instead of the weld.

Vectorized numpy over the cut elements only (O(10^2)-O(10^4), not a hot
kernel -- same policy as wake/cut_elements.py and mesh/wake_cut.py). Nothing
here is imported by the conforming solver paths.
"""

import numpy as np


def cut_element_stiffness(
    B_cut: np.ndarray, V_cut: np.ndarray, rho_cut: np.ndarray
) -> np.ndarray:
    """(n_cut, 4, 4) element stiffness K[e,a,b] = rho_e V_e (grad N_a . grad
    N_b) for the cut elements, from PicardOperator geometry (B (n,4,3),
    V (n,)) and the per-element weight rho (Laplace: ones)."""
    return np.einsum("e,ead,ebd->eab", V_cut * rho_cut, B_cut, B_cut)


def multivalued_redirection_coo(op, cm, rho_tilde=None):
    """COO triplets converting the single-valued cut-element contribution
    (already present in ``op.assemble_matrix(rho_tilde)``) into the
    multivalued one.

    For every cut-element entry (a, b) whose nodes are on OPPOSITE sides,
    subtract K[a,b] from (main(a), main(b)) and add it at (main(a),
    aux(b)) -- i.e. redirect the column from the node's main DOF to its aux
    DOF. Same-side entries and non-cut elements are untouched (their columns
    already point at the correct main DOF).

    Args:
        op: PicardOperator (supplies elements, B, V)
        cm: CutElementMap (supplies cut_elems, node_side, ext_dof_of_node)
        rho_tilde: (n_tets,) element weight, or None for the Laplace limit
            rho == 1 (must match the weight passed to assemble_matrix)

    Returns:
        (rows, cols, data) int64/int64/float64 COO triplets into the
        extended (n_total x n_total) matrix. Empty arrays if no cut element.
    """
    cut = np.asarray(cm.cut_elems, dtype=np.int64)
    if len(cut) == 0:
        e = np.empty(0, dtype=np.int64)
        return e, e, np.empty(0, dtype=np.float64)

    conn = np.asarray(op.elements, dtype=np.int64)[cut]          # (n_cut, 4)
    rho_cut = (
        np.ones(len(cut), dtype=np.float64)
        if rho_tilde is None
        else np.asarray(rho_tilde, dtype=np.float64)[cut]
    )
    K = cut_element_stiffness(op.B[cut], op.V[cut], rho_cut)     # (n_cut,4,4)

    side = cm.node_side[conn]                                    # (n_cut, 4)
    aux = cm.ext_dof_of_node[conn]                               # (n_cut, 4)
    opp = side[:, :, None] != side[:, None, :]                   # (n_cut,4,4)
    e_i, a_i, b_i = np.nonzero(opp)

    row = conn[e_i, a_i]                                         # main(a)
    col_main = conn[e_i, b_i]                                    # main(b)
    col_aux = aux[e_i, b_i]                                      # aux(b)
    kval = K[e_i, a_i, b_i]

    rows = np.concatenate([row, row])
    cols = np.concatenate([col_main, col_aux])
    data = np.concatenate([-kval, kval])
    return rows, cols, data


def nonte_aux_rows(cm):
    """Extended-DOF ids of the aux rows that carry the WAKE LS condition:
    the aux DOFs of cut nodes that are NOT trailing-edge nodes. TE aux DOFs
    are excluded -- their row stays lower-side mass conservation so the TE
    jump (= Gamma) is free (implicit Kutta; design_track_b.md sections
    2.1/2.3)."""
    is_te = np.zeros(cm.n_main, dtype=bool)
    is_te[cm.te_nodes] = True
    cut_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
    nonte = cut_nodes[~is_te[cut_nodes]]
    return cm.ext_dof_of_node[nonte].astype(np.int64)


def mass_conservation_coo(op, cm, rho_upper=None, rho_lower=None):
    """Extended (n_total x n_total) multivalued MASS-CONSERVATION matrix as
    COO triplets (Lopez's "replacement" assembly, design_track_b.md section
    2.1, BEFORE the non-TE aux rows are overwritten by the wake LS):

      - cut elements -> two copies, dofs_upper and dofs_lower (the same
        element matrix scattered twice; section 2.5 full integration), each
        weighted by ITS OWN side density (design_track_b.md section 5.2/D10:
        a cut element has two velocity states, so rho differs per side);
      - te_lower ("below-TE fan") elements -> the TE node reference uses its
        aux DOF (the lower value; Lopez fig. 3.6c, D11), lower-side density;
      - every other element -> ordinary single-valued (main DOFs), that
        element's own-side density.

    The resulting rows are: main rows = the node's own-side copy; TE aux rows
    = lower-side mass conservation (kept); non-TE aux rows = the other-side
    copy (discarded and replaced by wake_ls_coo).

    Args:
        rho_upper, rho_lower: (n_tets,) densities of the upper/lower fields
            (from MultivaluedOperator.element_densities). Both None -> the
            Laplace limit rho == 1.
    """
    el = np.asarray(op.elements, dtype=np.int64)
    n_tets = len(el)
    if rho_upper is None and rho_lower is None:
        rho_upper = rho_lower = np.ones(n_tets, dtype=np.float64)
    rho_upper = np.asarray(rho_upper, dtype=np.float64)
    rho_lower = np.asarray(rho_lower, dtype=np.float64)

    is_special = np.zeros(n_tets, dtype=bool)
    is_special[cm.cut_elems] = True
    is_special[cm.te_lower_elems] = True
    plain = np.flatnonzero(~is_special)
    # A plain element is single-sided; pick its density from that side.
    plain_side_plus = cm.node_side[el[plain]].max(axis=1) == 1
    rho_plain = np.where(plain_side_plus, rho_upper[plain], rho_lower[plain])

    is_te = np.zeros(cm.n_main, dtype=bool)
    is_te[cm.te_nodes] = True
    tel = np.asarray(cm.te_lower_elems, dtype=np.int64)
    tel_conn = el[tel]
    tel_dof = np.where(is_te[tel_conn], cm.ext_dof_of_node[tel_conn], tel_conn)

    r_parts, c_parts, d_parts = [], [], []

    def _scatter(elem_ids, dofvec, rho_e):
        if len(elem_ids) == 0:
            return
        K = cut_element_stiffness(op.B[elem_ids], op.V[elem_ids], rho_e)
        r_parts.append(np.repeat(dofvec, 4, axis=1).reshape(-1))
        c_parts.append(np.tile(dofvec, (1, 4)).reshape(-1))
        d_parts.append(K.reshape(-1))

    _scatter(plain, el[plain], rho_plain)
    _scatter(cm.cut_elems, cm.dofs_upper, rho_upper[cm.cut_elems])
    _scatter(cm.cut_elems, cm.dofs_lower, rho_lower[cm.cut_elems])
    _scatter(tel, tel_dof, rho_lower[tel])

    return (
        np.concatenate(r_parts),
        np.concatenate(c_parts),
        np.concatenate(d_parts),
    )


def wake_ls_coo(op, cm, wake_normal: np.ndarray, u_hat: np.ndarray):
    """Wake least-squares COO (Lopez eqs. 3.43-3.51), ADDED to the mass
    conservation on BOTH the upper and lower rows of every non-TE wake node.

    Per cut element the LS functional (Lopez eq. 3.45) is

        Pi_e = 1/2 V_e [ (n_hat . grad[phi])^2 + (u_hat . grad[phi])^2 ]

    with grad[phi] = sum_c B_c delta_c the element-constant jump gradient and
    delta_c = phi_u^c - phi_l^c = x[dofs_upper[c]] - x[dofs_lower[c]] the
    nodal jump. The dissertation takes the stationarity w.r.t. BOTH DOF sets
    (eqs. 3.46 AND 3.47), giving all FOUR Jacobian blocks (3.48-3.51):

        R_u^a = +V_e sum_c P^e_ac delta_c   -> row dofs_upper[a]
        R_l^a = -V_e sum_c P^e_ac delta_c   -> row dofs_lower[a]
        P^e_ac = (n_hat . B_a)(n_hat . B_c) + (u_hat . B_a)(u_hat . B_c)

    (the sign of R_l is fixed by J_lu = -(..) / J_ll = +(..) in 3.50/3.51).

    ROW ASSIGNMENT (Lopez section 3.5.4, REPLACEMENT -- the aux DOF rows carry
    the wake BC INSTEAD of mass conservation): the aux DOF of a node ABOVE the
    wake IS phi_l (eq. 3.34), so its row is R_l; the aux DOF of a node BELOW
    the wake IS phi_u (eq. 3.33), so its row is R_u. That is why all FOUR
    Jacobian blocks (3.48-3.51) appear even though each aux row carries only
    ONE residual: each row couples to BOTH the phi_u and phi_l columns. The
    main rows keep pure mass conservation. (Measured: superposing the LS onto
    the mass-conservation rows instead corrupts the fluid equations and
    over-circulates 3x -- Gamma 0.34 vs 0.12 on NACA0012 a=2.)

    A replacement row may be scaled by any nonzero constant, so the -side_a
    factor from d(delta_a)/d(phi_aux^a) is dropped: the row is written as
    +k on the phi_u^c columns and -k on the phi_l^c columns.

    TE nodes are excluded (Lopez section 3.5.4): they belong to the body, so
    the slip BC applies and their aux DOF carries mass conservation only --
    that exception is what lets the TE jump (= Gamma) be free (implicit Kutta).

    Args:
        op: PicardOperator
        cm: CutElementMap
        wake_normal: (n_cut, 3) "+"-side wake normal per cut element
            (WakeLevelSet.surface_normals at the element centroids)
        u_hat: (3,) freestream unit direction

    Returns:
        (rows, cols, data) into the extended matrix, to be ADDED to the
        mass-conservation triplets.
    """
    cut = np.asarray(cm.cut_elems, dtype=np.int64)
    if len(cut) == 0:
        e = np.empty(0, dtype=np.int64)
        return e, e, np.empty(0, dtype=np.float64)

    el = np.asarray(op.elements, dtype=np.int64)
    conn = el[cut]                                              # (nc, 4)
    B = op.B[cut]                                               # (nc, 4, 3)
    V = op.V[cut]
    n = np.asarray(wake_normal, dtype=np.float64)               # (nc, 3)
    u = np.asarray(u_hat, dtype=np.float64)
    u = u / np.linalg.norm(u)

    Bn = np.einsum("ead,ed->ea", B, n)                         # B_a . n_hat
    Bu = np.einsum("ead,d->ea", B, u)                          # B_a . u_hat
    P = Bn[:, :, None] * Bn[:, None, :] + Bu[:, :, None] * Bu[:, None, :]
    VP = V[:, None, None] * P                                   # (nc, 4, 4)

    up = np.asarray(cm.dofs_upper, dtype=np.int64)              # (nc, 4)
    lo = np.asarray(cm.dofs_lower, dtype=np.int64)
    is_te = np.zeros(cm.n_main, dtype=bool)
    is_te[cm.te_nodes] = True
    row_ok = ~is_te[conn]                                       # (nc, 4)

    aux = cm.ext_dof_of_node[conn]                              # (nc, 4)
    e_i, a_i, c_i = np.nonzero(row_ok[:, :, None] & np.ones((1, 1, 4), bool))
    k = VP[e_i, a_i, c_i]
    # aux row of node a;  delta_c = x[up[c]] - x[lo[c]]
    row = aux[e_i, a_i]
    rows = np.concatenate([row, row])
    cols = np.concatenate([up[e_i, c_i], lo[e_i, c_i]])
    data = np.concatenate([k, -k])
    return rows, cols, data


def te_kutta_coo(mvop, phi_ext, weights=None):
    """TE Kutta row: the NONLINEAR pressure-equality (Bernoulli) condition,
    linearized about the current iterate. One row per TE node, placed on that
    node's AUX dof (Track B, B4; design_track_b.md section 9.5 route (b)).

    The physical Kutta condition at a trailing edge is equal pressure on the
    two sides, i.e. equal speed:

        |q_u|^2 - |q_l|^2 = 0      (isentropic p is a monotone function of q^2)

    which factorizes EXACTLY as

        (q_u + q_l) . (q_u - q_l) = 0

    with q_u, q_l the recovered nodal velocities on the TE's UPPER and LOWER
    control volumes (MultivaluedOperator.te_velocities). Freezing the mean
    s = q_u + q_l at the previous iterate leaves a row that is LINEAR in phi:

        s . (q_u - q_l) = 0,

    so it drops straight into the Picard loop (re-linearized each outer
    iteration, exactly like the density lag) and converges to the exact
    nonlinear condition. At a freestream start s ~ 2 u_inf and the row reduces
    to the classical linearized Kutta.

    ★ Why this and not the wake LS: q_u and q_l come from DIFFERENT element
    sets, so q_u - q_l is NOT the gradient of a jump field and does NOT vanish
    for a constant jump. The wake LS contracts grad[phi] on a SINGLE cut
    element, which IS identically zero for a constant jump (partition of unity,
    measured 1.9e-16), so it cannot pin Gamma -- the B3 blocker
    (design_track_b.md section 9.2). This row is non-degenerate in Gamma and is
    what replaces the old TE aux row (lower-side mass conservation), whose
    control volume was up/down asymmetric on a symmetric airfoil.

    Note the mesh is NOT symmetric at the TE and the eps shift biases the
    upper/lower split, so the condition is deliberately a POINTWISE PHYSICAL
    statement (equal pressure) that needs no symmetry -- symmetrizing the
    control volume is not an available route (user-arbitrated 2026-07-12).

    Returns:
        (rows, cols, data) into the extended matrix; rows are the TE aux dofs.
        The RHS is zero.
    """
    cm = mvop.cm
    op = mvop.op
    qu, ql = mvop.te_velocities(phi_ext)
    s = qu + ql                                          # (n_te, 3), frozen

    # Track B B8 tip-taper: an optional per-TE-node factor F_i in [0, 1]
    # scaling this (homogeneous) pressure-equality row. Used ONLY in the
    # row BLEND F*K + (1-F)*W (see te_weld_coo): scaling K alone is a no-op
    # because the RHS is zero, so this must be paired with the weld. weights
    # is None on every existing path -> bit-identical.
    if weights is not None:
        weights = np.asarray(weights, dtype=np.float64)
        if weights.shape != (len(cm.te_nodes),):
            raise ValueError(
                f"te_kutta_coo weights must be ({len(cm.te_nodes)},), "
                f"got {weights.shape}"
            )

    rows, cols, data = [], [], []
    for i, t in enumerate(cm.te_nodes):
        cv = mvop._te_cv[i]
        row = cm.ext_dof_of_node[t]
        f_i = 1.0 if weights is None else float(weights[i])
        for key, sign in (("upper", 1.0), ("lower", -1.0)):
            e = cv[f"{key}_elems"]
            d = cv[f"{key}_dofs"]
            if len(e) == 0:
                continue
            w = op.V[e] / op.V[e].sum()                  # volume weights
            # d(s . q_side)/d(x[d[e,a]]) = w_e * (s . B[e,a,:])
            coef = f_i * sign * w[:, None] * np.einsum("ead,d->ea", op.B[e], s[i])
            rows.append(np.full(d.size, row, dtype=np.int64))
            cols.append(d.reshape(-1))
            data.append(coef.reshape(-1))

    if not rows:
        e = np.empty(0, dtype=np.int64)
        return e, e, np.empty(0, dtype=np.float64)
    return (np.concatenate(rows), np.concatenate(cols),
            np.concatenate(data))


def te_weld_coo(cm, weights):
    """TE-restricted continuity-weld triplets scaled by (1 - F_i), the
    second term of the Track B B8 tip-taper row blend (roadmap Track B B8;
    P13/G13.2 finding (8)).

    Per TE node i the blended Kutta row on its aux DOF is

        F_i * [ s . (q_u - q_l) ]  +  (1 - F_i) * [ phi_aux - phi_main ]  = 0

    The first term is `te_kutta_coo(mvop, phi_ext, weights=F)`; this function
    supplies the second, i.e. the same weld `aux_k - main_j = 0` that
    `continuity_closure_coo` applies to every cut node, but here restricted to
    the TE nodes and weighted by (1 - F_i). F_i = 1 (inboard) -> no weld,
    pure pressure Kutta (bit-identical to the untapered path); F_i = 0 (tip)
    -> pure weld -> the jump is pinned to 0 at that node -> the tip is
    unloaded (the level-set analogue of the conforming Gamma_eff -> 0).

    Because the weld row is NOT proportional to the (homogeneous) pressure
    row, the blend genuinely changes the solution, unlike scaling the
    pressure row alone. The RHS stays 0 (both terms are homogeneous).

    Args:
        cm: the CutElementMap (needs `te_nodes`, `ext_dof_of_node`).
        weights: (n_te,) tip-taper factors F_i, ordered like `cm.te_nodes`.

    Returns:
        (rows, cols, data) into the extended matrix; empty if all F_i == 1.
    """
    te = np.asarray(cm.te_nodes, dtype=np.int64)
    weights = np.asarray(weights, dtype=np.float64)
    if weights.shape != te.shape:
        raise ValueError(
            f"te_weld_coo weights must be ({te.size},), got {weights.shape}"
        )
    one_minus_f = 1.0 - weights
    sel = one_minus_f != 0.0                 # only blended (near-tip) nodes
    if not np.any(sel):
        e = np.empty(0, dtype=np.int64)
        return e, e, np.empty(0, dtype=np.float64)
    nodes = te[sel]
    aux_k = cm.ext_dof_of_node[nodes]
    coef = one_minus_f[sel]
    rows = np.concatenate([aux_k, aux_k])
    cols = np.concatenate([aux_k, nodes])
    data = np.concatenate([coef, -coef])
    return rows, cols, data


def newton_terms23_side_coo(op, side, u_inf=1.0):
    """LS-Newton density-sensitivity blocks (Terms 2 + 3, design.md (6.3) /
    Lopez B.4) for ONE side of the multivalued assembly, as COO triplets on
    the EXTENDED DOFs (B6-Newton; design_track_b.md section 5.5).

    Term 2 (local): dR_a/dphi_b += 2 inv_u2 s_e V_e (gradphi_e.B_a)(gradphi_e.B_b)
        rows = dofvec[e], cols = readvec[e].
    Term 3 (upstream): dR_a/dphi_k += 2 inv_u2 s_u V_e (gradphi_e.B_a)
        (gradphi_u.B_u,k), rows = dofvec[e], cols = readvec[u(e)] (SAME side --
        the walk graph is side-restricted, so u(e) is in this side's set).

    ROWS AND COLUMNS COME FROM DIFFERENT MAPS (B19 Leg A, from the A3/C1
    finding). `dofvec` is where the element's residual LANDS (the
    `mass_conservation_coo` scatter); `readvec` is which DOF the element's side
    FIELD READS (`side_potentials`' per-node rule). They coincide on cut
    elements -- which is why using one for both survived every 2.5-D gate --
    but differ on **mixed-side plain** elements, which exist only in 3-D. There
    the residual depends on cut nodes' aux DOFs through the side density while
    the old code scattered that sensitivity onto a MAIN column, so J was not
    dR/dphi (rel err 1.146e-01 on the affected directions vs 6.33e-10
    elsewhere; `cases/analysis/c1_ls_jacobian_fd/`).

    `side` is a `MultivaluedOperator.newton_side_data` dict (s_e/s_u/upstream/
    grad/lim_mask/keep/dofvec/readvec). Sensitivities are masked by lim_mask
    exactly as P8 (flat clamp -> zero derivative on limited elements). Only
    elements in this side's set (keep) with a nonzero sensitivity emit.

    This is the per-side, DOF-indirected analogue of
    PicardOperator.assemble_newton_jacobian's Term 2/3 -- vectorized over the
    active elements (cut + supersonic, O(1e2-1e3)), not a hot kernel.
    """
    B, V = op.B, op.V
    keep = side["keep"]
    dof = side["dofvec"]                         # rows: where R lands
    read = side.get("readvec")                   # cols: what the field reads
    if read is None:                             # pre-B19 dicts
        read = dof
    grad = side["grad"]
    lim = side["lim_mask"]
    inv_u2 = 1.0 / (u_inf * u_inf)
    # Two gradient factors, matching the two DOF maps (B19 Leg A). The element
    # residual is rho_tilde(grad_READ) * V * (grad_SCATTER . B_a), so the ROW
    # factor uses the scatter field's gradient and the COLUMN factor uses the
    # read (side) field's. Identical on cut elements, different on mixed-side
    # plain ones.
    grad_row = side.get("grad_row")
    if grad_row is None:                         # pre-B19 dicts
        grad_row = grad
    gB_row = np.einsum("ed,ead->ea", grad_row, B)     # rows (index a)
    gB = np.einsum("ed,ead->ea", grad, B)             # cols (index b / k)

    s_e = side["s_e"] * lim                      # limiter flat clamp -> 0
    up = side["upstream"]
    s_u = side["s_u"] * lim[up]

    rows_l, cols_l, data_l = [], [], []

    # -- Term 2 -------------------------------------------------------------
    e2 = np.flatnonzero(keep & (s_e != 0.0))
    if len(e2):
        w2 = (2.0 * inv_u2) * s_e[e2] * V[e2]           # (m,)
        blk = w2[:, None, None] * (gB_row[e2][:, :, None] * gB[e2][:, None, :])
        rows_l.append(np.repeat(dof[e2], 4, axis=1).reshape(-1))
        cols_l.append(np.tile(read[e2], (1, 4)).reshape(-1))
        data_l.append(blk.reshape(-1))

    # -- Term 3 -------------------------------------------------------------
    e3 = np.flatnonzero(keep & (s_u != 0.0))
    if len(e3):
        u3 = up[e3]
        w3 = (2.0 * inv_u2) * s_u[e3] * V[e3]            # (m,)
        # (gradphi_e.B_e,a) outer (gradphi_u.B_u,k)
        blk = w3[:, None, None] * (gB_row[e3][:, :, None] * gB[u3][:, None, :])
        de = dof[e3]                                     # rows: nodes of e
        du = read[u3]                                    # cols: nodes of u(e)
        rows_l.append(np.repeat(de, 4, axis=1).reshape(-1))
        cols_l.append(np.tile(du, (1, 4)).reshape(-1))
        data_l.append(blk.reshape(-1))

    if not rows_l:
        e = np.empty(0, dtype=np.int64)
        return e, e, np.empty(0, dtype=np.float64)
    return (np.concatenate(rows_l).astype(np.int64),
            np.concatenate(cols_l).astype(np.int64),
            np.concatenate(data_l))


def te_kutta_jacobian_coo(mvop, phi_ext):
    """EXACT Jacobian of the nonlinear TE pressure-equality Kutta residual
    R = |q_u|^2 - |q_l|^2 (B6-Newton), replacing the frozen-mean linearized
    row of `te_kutta_coo` on the TE aux rows.

    dR/dphi = 2 q_u . dq_u/dphi - 2 q_l . dq_l/dphi, i.e. the same control-
    volume structure as te_kutta_coo but with the UPPER part contracted
    against s = 2 q_u and the LOWER part against s = 2 q_l (each side its own
    velocity), instead of both against the frozen mean q_u + q_l. Quadratic in
    phi, so this exact derivative is what gives Newton its quadratic rate; the
    frozen-mean row is exact only as a value (row . x = R), not as a slope.

    Returns (rows, cols, data) on the TE aux rows (RHS handled by the driver:
    the Newton residual on these rows is R itself, computed separately).
    """
    cm = mvop.cm
    op = mvop.op
    qu, ql = mvop.te_velocities(phi_ext)

    rows, cols, data = [], [], []
    for i, t in enumerate(cm.te_nodes):
        cv = mvop._te_cv[i]
        row = cm.ext_dof_of_node[t]
        for key, q in (("upper", 2.0 * qu[i]), ("lower", -2.0 * ql[i])):
            e = cv[f"{key}_elems"]
            d = cv[f"{key}_dofs"]
            if len(e) == 0:
                continue
            w = op.V[e] / op.V[e].sum()
            coef = w[:, None] * np.einsum("ead,d->ea", op.B[e], q)
            rows.append(np.full(d.size, row, dtype=np.int64))
            cols.append(d.reshape(-1))
            data.append(coef.reshape(-1))

    if not rows:
        e = np.empty(0, dtype=np.int64)
        return e, e, np.empty(0, dtype=np.float64)
    return (np.concatenate(rows), np.concatenate(cols),
            np.concatenate(data))


def te_kutta_residual(mvop, phi_ext):
    """R = |q_u|^2 - |q_l|^2 per TE node (the nonlinear pressure-equality
    Kutta residual; B6-Newton). Ordered like cm.te_nodes / the TE aux rows."""
    qu, ql = mvop.te_velocities(phi_ext)
    return np.einsum("ij,ij->i", qu, qu) - np.einsum("ij,ij->i", ql, ql)


def continuity_closure_coo(cm):
    """Aux-row COO triplets for the B2 continuity ("weld") closure:
    aux_k - main_j = 0 for every aux DOF k belonging to cut node j. This
    ties the two copies together, so a non-lifting solve stays single
    valued. B3 replaces these rows with the g1+g2 wake least-squares
    condition (design_track_b.md D2); the RHS of every weld row is 0.

    Returns (rows, cols, data) into the extended matrix.
    """
    cut_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0).astype(np.int64)
    aux_k = cm.ext_dof_of_node[cut_nodes]
    rows = np.concatenate([aux_k, aux_k])
    cols = np.concatenate([aux_k, cut_nodes])
    data = np.concatenate([
        np.ones(len(aux_k), dtype=np.float64),
        -np.ones(len(aux_k), dtype=np.float64),
    ])
    return rows, cols, data
