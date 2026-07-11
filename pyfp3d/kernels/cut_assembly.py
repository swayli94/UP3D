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
