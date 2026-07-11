"""
Multivalued FE operator on a level-set cut mesh (Track B, B2).

Wraps a `PicardOperator` (the ordinary single-valued P1 assembly, untouched)
and a `CutElementMap` (B1) into the EXTENDED (n_total = n_main + n_ext) DOF
space, where each cut-element node carries an auxiliary DOF holding the
value on the OTHER side of the wake (Lopez dissertation eqs. 3.33-3.34;
design_track_b.md sections 2.1/2.5). The extended matrix is

    [ mass conservation (multivalued, side-selected) | ... ]   main rows
    [ aux-row closure                                | ... ]   aux rows

built as a sparse correction to `op.assemble_matrix(rho)` plus an aux-row
block (see kernels/cut_assembly.py for the derivation). At B2 the aux rows
are the continuity ("weld") closure aux_k = main_j, which keeps a
non-lifting solve single valued and lets freestream / MMS / a=0 reproduce
the single-valued solution to machine precision. B3 swaps that block for
the g1+g2 wake least-squares condition (implicit Kutta, design_track_b.md
D2), at which point [phi] becomes nonzero and the mesh carries lift.

The extended matrix is structurally NONSYMMETRIC (the aux rows use a
different operator than the main columns feeding them -- design_track_b.md
section 5.3), so it is solved with a direct/GMRES path, never CG. This
module owns only assembly + DOF bookkeeping; drivers live in
solve/picard_ls.py. Nothing here is imported by the conforming solver paths.
"""

import numpy as np
import scipy.sparse as sp

from pyfp3d.kernels.cut_assembly import (
    continuity_closure_coo,
    multivalued_redirection_coo,
)
from pyfp3d.kernels.jacobian import PicardOperator


class MultivaluedOperator:
    """Extended-DOF assembly on a cut mesh.

    Args:
        nodes: (n_nodes, 3)
        elements: (n_tets, 4)
        cm: CutElementMap for this mesh + wake level set (B1)

    Attributes:
        op: the underlying single-valued PicardOperator
        cm: the CutElementMap
        n_main, n_ext, n_total: DOF counts
    """

    def __init__(self, nodes: np.ndarray, elements: np.ndarray, cm):
        self.op = PicardOperator(nodes, elements)
        self.cm = cm
        self.n_main = cm.n_main
        self.n_ext = cm.n_ext_dofs
        self.n_total = cm.n_total_dofs
        if self.op.n_nodes != self.n_main:
            raise ValueError(
                "CutElementMap.n_main != mesh node count -- the map was built "
                "on a different mesh"
            )
        # Aux-row closure block is state-independent (weld); cache it.
        self._closure_coo = continuity_closure_coo(cm)

    def assemble_matrix(self, rho_tilde=None, closure: str = "continuity"):
        """Extended (n_total x n_total) multivalued matrix.

        Args:
            rho_tilde: (n_tets,) element weight (Laplace: None -> rho == 1)
            closure: aux-row block. Only "continuity" (the B2 weld) is
                implemented; B3 adds "wake_ls".

        Returns:
            scipy CSR of shape (n_total, n_total).
        """
        if closure != "continuity":
            raise NotImplementedError(
                f"closure={closure!r} not available until B3 (wake LS)"
            )
        a_main = self.op.assemble_matrix(rho_tilde).tocoo()
        r_row, r_col, r_data = multivalued_redirection_coo(
            self.op, self.cm, rho_tilde
        )
        c_row, c_col, c_data = self._closure_coo

        rows = np.concatenate([a_main.row, r_row, c_row])
        cols = np.concatenate([a_main.col, r_col, c_col])
        data = np.concatenate([a_main.data, r_data, c_data])
        A = sp.coo_matrix(
            (data, (rows, cols)), shape=(self.n_total, self.n_total)
        ).tocsr()
        A.sum_duplicates()
        return A

    # -- DOF bookkeeping / post-processing -----------------------------------

    def main_potential(self, phi_ext: np.ndarray) -> np.ndarray:
        """The single-valued (own-side) potential on the n_main mesh nodes."""
        return np.asarray(phi_ext)[: self.n_main]

    def te_jump(self, phi_ext: np.ndarray) -> np.ndarray:
        """[phi] = phi_upper - phi_lower at the TE nodes.

        The main DOF holds a node's OWN-side value, the aux DOF the other
        side, so [phi] = side * (main - aux). TE nodes are shifted "+"
        (design_track_b.md section 2.3), so this is main - aux there; the
        sign factor keeps it correct for any queried node. B2 non-lifting
        solves give ~0 (the weld forbids a jump); the real Gamma emerges in
        B3.
        """
        return self.node_jump(phi_ext, self.cm.te_nodes)

    def node_jump(self, phi_ext: np.ndarray, node_ids: np.ndarray) -> np.ndarray:
        """side * (main - aux) at the given cut nodes (0 where a node has no
        aux DOF -- i.e. is not a cut-element node)."""
        phi_ext = np.asarray(phi_ext, dtype=np.float64)
        node_ids = np.asarray(node_ids, dtype=np.int64)
        aux = self.cm.ext_dof_of_node[node_ids]
        has_aux = aux >= 0
        jump = np.zeros(len(node_ids), dtype=np.float64)
        j = node_ids[has_aux]
        jump[has_aux] = self.cm.node_side[j] * (phi_ext[j] - phi_ext[aux[has_aux]])
        return jump
