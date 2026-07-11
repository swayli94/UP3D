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
    mass_conservation_coo,
    multivalued_redirection_coo,
    nonte_aux_rows,
    te_kutta_coo,
    wake_ls_coo,
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

    def __init__(self, nodes: np.ndarray, elements: np.ndarray, cm,
                 levelset=None):
        self.op = PicardOperator(nodes, elements)
        self.cm = cm
        self.levelset = levelset
        self.n_main = cm.n_main
        self.n_ext = cm.n_ext_dofs
        self.n_total = cm.n_total_dofs
        if self.op.n_nodes != self.n_main:
            raise ValueError(
                "CutElementMap.n_main != mesh node count -- the map was built "
                "on a different mesh"
            )
        # Continuity (weld) aux-row block is state-independent; cache it.
        self._closure_coo = continuity_closure_coo(cm)
        # Wake-LS geometry (per cut element): centroid normals + freestream
        # direction. Constant in phi, so the LS block is cached too.
        self._nonte_rows = nonte_aux_rows(cm)
        self._elem_plus = None  # lazily built own-side element mask
        self._ls_coo = None
        if levelset is not None and len(cm.cut_elems) > 0:
            el = np.asarray(elements, dtype=np.int64)
            centroids = nodes[el[cm.cut_elems]].mean(axis=1)
            n_hat = levelset.surface_normals(centroids)
            self._ls_coo = wake_ls_coo(self.op, cm, n_hat, levelset.direction)
        self._te_cv = self._build_te_control_volumes()

    # -- TE control volumes (B4) ---------------------------------------------

    def _build_te_control_volumes(self):
        """Per TE node, the two control volumes the B4 TE-Kutta velocities are
        recovered on:

          UPPER = the WALL-ADJACENT elements that reference the TE through its
                  MAIN dof  (i.e. the UPPER body surface at the TE),
          LOWER = the WALL-ADJACENT elements that reference it through its AUX
                  dof        (i.e. the LOWER body surface at the TE).

        Each entry carries the element id and the DOF vector of THAT copy, so a
        per-side velocity can be recovered exactly as the copy's assembly sees
        it.

        ★ WALL-ADJACENT (an element with a wall FACE, >= 3 wall nodes) is the
        whole point, and it is what makes the condition accurate: the Kutta
        condition is a statement about the UPPER- and LOWER-SURFACE velocities
        at the TE, not about a volume average over the whole element fan.
        Measured on NACA0012 a=2 (incompressible, vs the conforming 0.1175 /
        0.1200): recovering over the FULL fan gives Gamma 0.1407 / 0.1355
        (~15% high, because interior and wake elements pollute the average),
        while the wall-adjacent recovery gives **0.1177 / 0.1191** -- 0.2% /
        0.8%.

        The two volumes are NOT mirror images (the mesh is not symmetric at the
        TE and the eps shift biases the split), which is exactly why the TE
        condition must be a POINTWISE PHYSICAL statement (equal pressure)
        rather than anything that assumes symmetry -- symmetrising the control
        volume is not an available route (user-arbitrated 2026-07-12).
        """
        cm = self.cm
        el = np.asarray(self.op.elements, dtype=np.int64)
        n_tets = len(el)
        cut_index = np.full(n_tets, -1, dtype=np.int64)
        cut_index[cm.cut_elems] = np.arange(len(cm.cut_elems))
        is_tel = np.zeros(n_tets, dtype=bool)
        is_tel[cm.te_lower_elems] = True
        is_te = np.zeros(cm.n_main, dtype=bool)
        is_te[cm.te_nodes] = True
        tel_dof = np.where(is_te[el], cm.ext_dof_of_node[el], el)  # aux at TE

        # Wall-adjacent = the element carries a wall FACE (>= 3 wall nodes).
        if cm.wall_nodes is None:
            on_wall = np.ones(len(el), dtype=bool)   # degenerate (no body)
        else:
            is_w = np.zeros(cm.n_main, dtype=bool)
            is_w[cm.wall_nodes] = True
            on_wall = is_w[el].sum(axis=1) >= 3

        cvs = []
        for t in cm.te_nodes:
            fan = np.flatnonzero((el == t).any(axis=1) & on_wall)
            eu, du, elo, dlo = [], [], [], []
            for e in fan:
                k = cut_index[e]
                if k >= 0:                       # cut: it has BOTH copies
                    eu.append(e); du.append(cm.dofs_upper[k])
                    elo.append(e); dlo.append(cm.dofs_lower[k])
                elif is_tel[e]:                  # below-TE fan -> lower only
                    elo.append(e); dlo.append(tel_dof[e])
                else:                            # plain "+" fan -> upper only
                    eu.append(e); du.append(el[e])
            cvs.append({
                "upper_elems": np.array(eu, dtype=np.int64),
                "upper_dofs": np.array(du, dtype=np.int64).reshape(-1, 4),
                "lower_elems": np.array(elo, dtype=np.int64),
                "lower_dofs": np.array(dlo, dtype=np.int64).reshape(-1, 4),
            })
        return cvs

    def te_velocities(self, phi_ext):
        """(q_upper, q_lower), each (n_te, 3): the volume-weighted recovered
        nodal velocity at every TE node, taken separately over the TE's UPPER
        and LOWER control volumes.

        ★ These are NOT two evaluations of one element's gradient -- they come
        from DIFFERENT element sets (the flow over the upper surface vs the
        lower surface at the TE). That is exactly why the pressure-equality
        Kutta built on them is NON-degenerate in Gamma, whereas the wake LS
        (which contracts grad[phi] on a SINGLE cut element) is identically zero
        for a constant jump and cannot pin Gamma at all (design_track_b.md
        section 9.2).
        """
        x = np.asarray(phi_ext, dtype=np.float64)
        B, V = self.op.B, self.op.V
        qu = np.zeros((len(self._te_cv), 3))
        ql = np.zeros((len(self._te_cv), 3))
        for i, cv in enumerate(self._te_cv):
            for key, out in (("upper", qu), ("lower", ql)):
                e = cv[f"{key}_elems"]
                d = cv[f"{key}_dofs"]
                if len(e) == 0:
                    continue
                grad = np.einsum("ead,ea->ed", B[e], x[d])   # (n, 3) per copy
                out[i] = (V[e][:, None] * grad).sum(axis=0) / V[e].sum()
        return qu, ql

    def assemble_matrix(self, rho_tilde=None, closure: str = "continuity",
                        te_kutta: str = "pressure", phi_ext=None):
        """Extended (n_total x n_total) multivalued matrix.

        Args:
            rho_tilde: (n_tets,) element weight (Laplace: None -> rho == 1)
            closure: aux-row block. Only "continuity" (the B2 weld) is
                implemented; B3 adds "wake_ls".

        Returns:
            scipy CSR of shape (n_total, n_total).
        """
        if closure == "continuity":
            a_main = self.op.assemble_matrix(rho_tilde).tocoo()
            r_row, r_col, r_data = multivalued_redirection_coo(
                self.op, self.cm, rho_tilde
            )
            c_row, c_col, c_data = self._closure_coo
            rows = np.concatenate([a_main.row, r_row, c_row])
            cols = np.concatenate([a_main.col, r_col, c_col])
            data = np.concatenate([a_main.data, r_data, c_data])
        elif closure == "wake_ls":
            if self._ls_coo is None:
                raise ValueError(
                    "wake_ls closure needs a WakeLevelSet (pass levelset= to "
                    "MultivaluedOperator) and a non-empty cut set"
                )
            # Extended mass conservation on EVERY row (main and aux), with the
            # g1+g2 wake LS ADDED on the non-TE wake rows (Lopez: "assembled
            # together with Eqs. (3.9) and (3.11)" -- a superposition, not a
            # replacement; see kernels/cut_assembly.py::wake_ls_coo). TE nodes
            # get mass conservation only, which is what leaves the TE jump
            # (= Gamma) free -- the implicit Kutta. rho_tilde is a
            # (rho_upper, rho_lower) pair, or None for the Laplace limit.
            if rho_tilde is None:
                rho_up = rho_lo = None
            else:
                rho_up, rho_lo = rho_tilde
            m_row, m_col, m_data = mass_conservation_coo(
                self.op, self.cm, rho_up, rho_lo
            )
            drop = np.zeros(self.n_total, dtype=bool)
            drop[self._nonte_rows] = True
            keep = ~drop[m_row]
            m_row, m_col, m_data = m_row[keep], m_col[keep], m_data[keep]
            ls_row, ls_col, ls_data = self._ls_coo
            extra = [(ls_row, ls_col, ls_data)]

            if te_kutta == "pressure":
                # B4: the TE aux row carries the NONLINEAR pressure-equality
                # Kutta instead of lower-side mass conservation (which is what
                # the wake LS cannot supply -- it is blind to a constant jump,
                # design_track_b.md section 9). To keep mass conserved at the
                # TE, the lower-side mass-conservation entries are re-routed
                # onto the TE MAIN row, which then carries the TOTAL
                # (upper + lower) balance over the whole TE fan. This also
                # avoids arbitrarily choosing WHICH side loses its equation.
                if phi_ext is None:
                    raise ValueError(
                        "te_kutta='pressure' needs phi_ext (the current "
                        "iterate) to linearize |q_u|^2 = |q_l|^2"
                    )
                te_aux = self.cm.ext_dof_of_node[self.cm.te_nodes]
                reroute = np.zeros(self.n_total, dtype=np.int64) - 1
                reroute[te_aux] = self.cm.te_nodes
                hit = reroute[m_row] >= 0
                m_row = m_row.copy()
                m_row[hit] = reroute[m_row[hit]]
                extra.append(te_kutta_coo(self, phi_ext))
            elif te_kutta != "mass":
                raise NotImplementedError(f"te_kutta={te_kutta!r} unknown")

            rows = np.concatenate([m_row] + [e[0] for e in extra])
            cols = np.concatenate([m_col] + [e[1] for e in extra])
            data = np.concatenate([m_data] + [e[2] for e in extra])
        else:
            raise NotImplementedError(f"closure={closure!r} unknown")

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

    def side_potentials(self, phi_ext: np.ndarray):
        """(phi_upper, phi_lower), two (n_main,) nodal fields: the "+" and
        "-" side potential at every mesh node. Away from the cut they
        coincide (single valued); at a cut node the own side is the main DOF
        and the other side is the aux DOF."""
        phi_ext = np.asarray(phi_ext, dtype=np.float64)
        phi_up = phi_ext[: self.n_main].copy()
        phi_lo = phi_ext[: self.n_main].copy()
        cut_nodes = np.flatnonzero(self.cm.ext_dof_of_node >= 0)
        aux = self.cm.ext_dof_of_node[cut_nodes]
        side = self.cm.node_side[cut_nodes]
        plus = side == 1
        # "+" node: main = upper, aux = lower;  "-" node: main = lower, aux = upper
        phi_lo[cut_nodes[plus]] = phi_ext[aux[plus]]
        phi_up[cut_nodes[~plus]] = phi_ext[aux[~plus]]
        return phi_up, phi_lo

    def element_densities(self, phi_ext, m_inf, gamma_air=1.4, u_inf=1.0):
        """(rho_upper, rho_lower) per element from the two side fields
        (subcritical isentropic density; nu == 0 so rho_tilde == rho). Feed
        both to assemble_matrix(rho_tilde=(rho_up, rho_lo), closure='wake_ls')
        so each cut-element copy is weighted by its own side.

        NOTE the two arrays are only VALID on their own side's elements
        (rho_upper on "+"/cut elements, rho_lower on "-"/cut/te_lower); each
        carries junk on the other side, where the side field mixes an aux
        (other-side) value into a bulk element and manufactures a spurious
        gradient. mass_conservation_coo consumes each only where it is valid,
        so the junk never enters the matrix -- but use own_side_field() for
        any diagnostic (drho, M_max)."""
        from pyfp3d.physics.isentropic import density_field

        phi_up, phi_lo = self.side_potentials(phi_ext)
        _, q2u = self.op.velocities(phi_up)
        rho_up = density_field(q2u / u_inf**2, m_inf, gamma_air).copy()
        _, q2l = self.op.velocities(phi_lo)
        rho_lo = density_field(q2l / u_inf**2, m_inf, gamma_air).copy()
        return rho_up, rho_lo

    def own_side_field(self, field_up, field_lo):
        """Pick each element's own-side value from an (upper, lower) pair,
        discarding the invalid other-side junk (see element_densities): "+"
        and cut elements take the upper array, "-" elements the lower. For a
        clean per-element diagnostic (drho lag, M_max) whose value is the one
        actually assembled."""
        if self._elem_plus is None:
            side = self.cm.node_side[np.asarray(self.op.elements)]
            is_cut = np.zeros(self.op.n_tets, dtype=bool)
            is_cut[self.cm.cut_elems] = True
            self._elem_plus = (side.max(axis=1) == 1) | is_cut
        return np.where(self._elem_plus, field_up, field_lo)

    def element_mach2(self, phi_ext, m_inf, gamma_air=1.4, u_inf=1.0):
        """Own-side element M^2 (max of the two sides on cut elements),
        junk-free -- for reporting M_max on the cut mesh."""
        from pyfp3d.physics.isentropic import mach_squared_field

        phi_up, phi_lo = self.side_potentials(phi_ext)
        _, q2u = self.op.velocities(phi_up)
        m2u = mach_squared_field(q2u / u_inf**2, m_inf, gamma_air).copy()
        _, q2l = self.op.velocities(phi_lo)
        m2l = mach_squared_field(q2l / u_inf**2, m_inf, gamma_air).copy()
        m2 = self.own_side_field(m2u, m2l)
        m2[self.cm.cut_elems] = np.maximum(
            m2u[self.cm.cut_elems], m2l[self.cm.cut_elems]
        )
        return m2

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
