"""
Cut-element identification and extended-DOF numbering (Track B, B1).

Given a mesh (which needs NO embedded wake sheet and NO "wake" tag) and a
WakeLevelSet, classify:

  - node sides: +1 (upper) / -1 (lower), never 0 -- nodes lying on (or
    within eps of) the wake surface are shifted to the "+" side (Lopez
    dissertation section 3.5.3 / design_track_b.md D4). eps is RELATIVE to
    the node's local edge length (eps_rel * h_node), so the rule is mesh-
    scale free and deterministic across an alpha sweep.
  - TE nodes: nodes on the trailing-edge line itself (|s| and |d| both
    below te_tol_rel * h_node, optionally restricted to wall nodes). They
    are shifted "+" like other on-surface nodes but flagged: per Lopez
    section 3.5.4 their aux DOF carries mass conservation, NOT the wake
    boundary conditions (consumed by B2 assembly).
  - cut (wake) elements: elements whose nodes change side AND whose s = 0
    crossing lies downstream of the TE (any crossing with d > 0). The
    downstream test is what excludes the ahead-of-leading-edge region,
    where the domain is connected across the wake plane's extension.
  - te_lower_elems: NOT cut, but contain a TE node with all their other
    nodes on the "-" side -- the below-TE fan whose TE reference must use
    the aux DOF in B2 (Lopez fig. 3.6c); recorded here so assembly never
    re-derives topology.
  - extended DOFs: one aux DOF per unique cut-element node (per-node
    granularity, shared across neighboring cut elements -- design_track_b.md
    D3), numbered n_main + k.

DOF-assignment convention (Lopez eqs. 3.33-3.34): for a cut element,
dofs_upper uses the main DOF at "+" nodes and the aux DOF at "-" nodes;
dofs_lower is the complement. The main DOF of a node always holds the
value on the node's OWN side.

All preprocessing is vectorized numpy (one-time, not a hot kernel; same
policy as mesh/wake_cut.py). Nothing here is imported by the conforming
solver paths.
"""

from typing import Dict, Optional

import numpy as np

from pyfp3d.wake.levelset import WakeLevelSet

# The 6 edges of a tet as local vertex index pairs.
_TET_EDGES = np.array(
    [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)], dtype=np.int64
)


def _node_min_edge_length(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """(n_nodes,) minimum incident edge length per node (inf if isolated)."""
    el = np.asarray(elements, dtype=np.int64)
    i = el[:, _TET_EDGES[:, 0]].ravel()
    j = el[:, _TET_EDGES[:, 1]].ravel()
    lengths = np.linalg.norm(nodes[i] - nodes[j], axis=1)
    h = np.full(len(nodes), np.inf)
    np.minimum.at(h, i, lengths)
    np.minimum.at(h, j, lengths)
    return h


class CutElementMap:
    """Cut-element census, side marking, and extended-DOF numbering.

    Args:
        nodes: (n_nodes, 3)
        elements: (n_tets, 4)
        levelset: WakeLevelSet
        wall_nodes: optional 1D array of wall node ids; when given, TE-node
            detection is restricted to them (recommended: the geometric
            |s|,|d| test alone is already sharp, this makes it airtight).
        eps_rel: on-surface shift threshold relative to local edge length.
        te_tol_rel: TE-node detection threshold relative to local edge length.

    Attributes (all set in __init__):
        s_raw, d: (n_nodes,) level-set offset / downstream coordinate
        s: (n_nodes,) shifted offset actually used for classification
        node_side: (n_nodes,) int8, +1 / -1, never 0
        shifted_nodes: ids whose |s_raw| < eps (includes the TE nodes)
        te_nodes: sorted TE node ids
        cut_elems: sorted cut (wake) element ids
        te_lower_elems: sorted below-TE fan element ids (disjoint from
            cut_elems by construction)
        n_main, n_ext_dofs, n_total_dofs
        ext_dof_of_node: (n_nodes,) aux DOF id or -1
        dofs_upper, dofs_lower: (n_cut, 4) DOF assignment per cut element
    """

    def __init__(
        self,
        nodes: np.ndarray,
        elements: np.ndarray,
        levelset: WakeLevelSet,
        wall_nodes: Optional[np.ndarray] = None,
        eps_rel: float = 1e-6,
        te_tol_rel: float = 1e-3,
    ):
        nodes = np.asarray(nodes, dtype=np.float64)
        el = np.asarray(elements, dtype=np.int64)
        n_nodes = len(nodes)
        self.n_main = n_nodes

        s_raw, d = levelset.evaluate(nodes)
        h_node = _node_min_edge_length(nodes, el)
        self.s_raw = s_raw
        self.d = d

        # TE nodes: on the TE line itself (s ~ 0 AND d ~ 0).
        te_tol = te_tol_rel * h_node
        te_mask = (np.abs(s_raw) < te_tol) & (np.abs(d) < te_tol)
        if wall_nodes is not None:
            on_wall = np.zeros(n_nodes, dtype=bool)
            on_wall[np.asarray(wall_nodes, dtype=np.int64)] = True
            te_mask &= on_wall
        self.te_nodes = np.flatnonzero(te_mask)

        # eps side-shift: on-surface nodes go "+" (deterministic).
        eps = eps_rel * h_node
        shift_mask = np.abs(s_raw) < eps
        s = s_raw.copy()
        s[shift_mask] = eps[shift_mask]
        self.s = s
        self.shifted_nodes = np.flatnonzero(shift_mask)
        self.node_side = np.where(s > 0.0, 1, -1).astype(np.int8)

        # Sign-change candidates, then the downstream-crossing test.
        side_e = self.node_side[el]                      # (n_tets, 4)
        cand = np.flatnonzero(side_e.min(axis=1) != side_e.max(axis=1))
        s_e = s[el[cand]]                                # (n_cand, 4)
        d_e = d[el[cand]]
        si = s_e[:, _TET_EDGES[:, 0]]
        sj = s_e[:, _TET_EDGES[:, 1]]
        di = d_e[:, _TET_EDGES[:, 0]]
        dj = d_e[:, _TET_EDGES[:, 1]]
        crossing = (si * sj) < 0.0                       # (n_cand, 6)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(crossing, si / (si - sj), 0.0)
        d_cross = np.where(crossing, di + t * (dj - di), -np.inf)
        self.cut_elems = cand[d_cross.max(axis=1) > 0.0]

        # Below-TE fan: not cut, contains a TE node, all non-TE nodes "-".
        is_te = np.zeros(n_nodes, dtype=bool)
        is_te[self.te_nodes] = True
        te_e = is_te[el]                                 # (n_tets, 4)
        has_te = te_e.any(axis=1)
        all_nonte_lower = np.all(te_e | (side_e == -1), axis=1)
        fan_mask = has_te & all_nonte_lower
        fan_mask[self.cut_elems] = False
        self.te_lower_elems = np.flatnonzero(fan_mask)

        # Per-node aux DOF over the unique cut-element nodes.
        cut_nodes = np.unique(el[self.cut_elems])
        self.n_ext_dofs = len(cut_nodes)
        self.n_total_dofs = self.n_main + self.n_ext_dofs
        self.ext_dof_of_node = np.full(n_nodes, -1, dtype=np.int64)
        self.ext_dof_of_node[cut_nodes] = self.n_main + np.arange(
            self.n_ext_dofs, dtype=np.int64
        )

        # Lopez eqs. 3.33-3.34: upper copy = main at "+", aux at "-".
        cut_conn = el[self.cut_elems]
        ext = self.ext_dof_of_node[cut_conn]
        upper_side = self.node_side[cut_conn] == 1
        self.dofs_upper = np.where(upper_side, cut_conn, ext)
        self.dofs_lower = np.where(upper_side, ext, cut_conn)

    def summary(self) -> Dict[str, int]:
        """Census dict for gate artifacts / asserts."""
        return {
            "n_main": int(self.n_main),
            "n_ext_dofs": int(self.n_ext_dofs),
            "n_cut_elems": int(len(self.cut_elems)),
            "n_te_nodes": int(len(self.te_nodes)),
            "n_te_lower_elems": int(len(self.te_lower_elems)),
            "n_shifted_nodes": int(len(self.shifted_nodes)),
            "n_side_plus": int(np.sum(self.node_side == 1)),
            "n_side_minus": int(np.sum(self.node_side == -1)),
        }
