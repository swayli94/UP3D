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
    crossing lies ON the (bounded) sheet -- downstream of the TE (d > 0)
    AND within the sheet's spanwise extent (0 <= q <= span_length). The
    downstream test excludes the ahead-of-leading-edge region, where the
    domain is connected across the wake plane's extension; the spanwise
    test excludes the region outboard of the wing TIP, where the sheet has
    ended (Gamma(tip) = 0 -- the conforming path's free-edge rule; without
    it a swept-wing level set cuts the whole plane extension beyond the
    tip, i.e. P5's far-field branch-ray artifact re-created). Both tests
    are inert on the quasi-2D meshes. The INBOARD end of the extent test
    (q >= 0) can be replaced by a caller-supplied `inboard_clip` (B25:
    the conforming sheet's fragment topology -- the sheet runs inboard
    until it hits the fuselage surface / the symmetry plane, so no free
    edge remains inside the fluid).
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
        inboard_clip: optional callable replacing the INBOARD spanwise clip
            (q >= 0) in the on-sheet test; None (default) keeps the legacy
            q >= 0 test BIT-IDENTICAL. When given, it is called with the
            s = 0 crossing points, shape (n_candidates, 6, 3) (one row per
            sign-changing tet edge), and must return a boolean
            (n_candidates, 6) mask: True = the crossing is on the sheet.
            The downstream (d > 0) and tip (q <= span_length) tests are
            unaffected. (B25: meshgen.fuselage.make_inboard_clip trims the
            sheet against the fuselage surface / symmetry plane.)
        outboard_fringe: optional OUTBOARD extension of the spanwise clip.
            ⊘ MEASURED NEGATIVE (GB31.3, 2026-07-22): both the faded-fringe
            (C1) and the no-clip sentinel (C3) forms were closed negative as
            an LS wing-tip C-class cure -- C1 gives an inboard backflow
            -19.5%, C3 diverges on the coarse mesh (the sheet reaches q~14).
            RETAINED default 0.0 = bit-identical (evidence machinery; see
            cases/analysis/b31_tip_termination/VERDICT.md); do NOT re-propose
            as a tip cure without new evidence. Mechanics below:
            (B31 candidate C1, GB31.3): the tip test becomes
            q_cross <= span_length + outboard_fringe, so the sheet
            continues a fringe width w past the tip and -- paired with
            MultivaluedOperator's span_blend fade -- ends there at a
            pure-weld edge instead of as a Heaviside at the tip. The aux
            rows are upwind-convective along +x, so welding OUTBOARD of
            the tip cannot re-level the inboard lifting solution (the
            mechanical distinction from the B8 inboard span_blend dead
            end, which welded the lifting sheet and cost 20% of the
            lift). 0.0 (default) keeps the legacy clip BIT-IDENTICAL;
            beyond_tip_elems then reports the elements rejected by the
            EXTENDED clip (the still-further-outboard plane extension).
            SENTINEL outboard_fringe = np.inf (B31 candidate C3): the tip
            test is DROPPED entirely (the downstream d > 0 and inboard
            tests stay), so the sheet never ends in the fluid -- the
            B25 "no free edge inside the domain" analogue at the tip.
            Geometry note: the B28 production sheet is the FLAT y = 0
            plane ruled along +x, so removing the clip extends the SAME
            plane outboard to the far field (the q parameterization
            slants along the TE sweep; the surface is y = 0 regardless).
            This is NOT the conforming model, which deliberately keeps a
            single-valued free edge at the tip (mesh/wake_cut.py); it is
            the unrolled planar tip-vortex continuation, and it
            re-exposes the far-field branch-ray question P5 fixed on the
            conforming path (the GB31.3 regression gate).

    Attributes (all set in __init__):
        s_raw, d, q: (n_nodes,) level-set offset / downstream coordinate /
            spanwise arclength (unclamped)
        s: (n_nodes,) shifted offset actually used for classification
        node_side: (n_nodes,) int8, +1 / -1, never 0
        shifted_nodes: ids whose |s_raw| < eps (includes the TE nodes)
        te_nodes: sorted TE node ids
        cut_elems: sorted cut (wake) element ids
        beyond_tip_elems: sorted ids of elements the wake PLANE crosses
            downstream of the TE but OUTBOARD of the sheet's spanwise end
            (span_length + outboard_fringe; empty on quasi-2D meshes; on a
            swept wing these are the beyond-tip elements the spanwise clip
            rejects; with the np.inf C3 sentinel the outboard rejection is
            gone, so only inboard-clip rejections can remain here)
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
        inboard_clip=None,
        outboard_fringe: float = 0.0,
    ):
        nodes = np.asarray(nodes, dtype=np.float64)
        el = np.asarray(elements, dtype=np.int64)
        n_nodes = len(nodes)
        self.n_main = n_nodes

        # Kept for the B4 TE Kutta: its control volumes are the WALL-adjacent
        # elements at the TE (the upper/lower body surface), not the whole fan.
        self.wall_nodes = (None if wall_nodes is None
                           else np.asarray(wall_nodes, dtype=np.int64))

        s_raw, d, q = levelset.evaluate(nodes)
        h_node = _node_min_edge_length(nodes, el)
        self.s_raw = s_raw
        self.d = d
        self.q = q
        self.span_length = levelset.span_length
        if outboard_fringe < 0.0:
            raise ValueError("outboard_fringe must be >= 0")
        self.outboard_fringe = float(outboard_fringe)

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
        # A3 (C7b): TE detection uses te_tol_rel*h (1e-3) but the side shift
        # uses eps_rel*h (1e-6), so a TE node landing in the gap between them
        # on the "-" side would be classed lower, and ext_dof_of_node[t] --
        # the UPPER value by construction -- would then be used as the LOWER
        # TE reference, inverting Lopez fig. 3.6c / D11. TE nodes are on the
        # sheet by definition, so shift them deterministically to "+" too.
        # No-op on every committed mesh (their TE nodes have s_raw exactly 0).
        shift_mask |= te_mask
        s = s_raw.copy()
        s[shift_mask] = eps[shift_mask]
        self.s = s
        self.shifted_nodes = np.flatnonzero(shift_mask)
        self.node_side = np.where(s > 0.0, 1, -1).astype(np.int8)

        # Sign-change candidates, then the on-sheet crossing test
        # (downstream of the TE AND within the sheet's spanwise extent).
        side_e = self.node_side[el]                      # (n_tets, 4)
        cand = np.flatnonzero(side_e.min(axis=1) != side_e.max(axis=1))
        s_e = s[el[cand]]                                # (n_cand, 4)
        d_e = d[el[cand]]
        q_e = q[el[cand]]
        si = s_e[:, _TET_EDGES[:, 0]]
        sj = s_e[:, _TET_EDGES[:, 1]]
        di = d_e[:, _TET_EDGES[:, 0]]
        dj = d_e[:, _TET_EDGES[:, 1]]
        qi = q_e[:, _TET_EDGES[:, 0]]
        qj = q_e[:, _TET_EDGES[:, 1]]
        crossing = (si * sj) < 0.0                       # (n_cand, 6)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(crossing, si / (si - sj), 0.0)
        d_cross = di + t * (dj - di)
        q_cross = qi + t * (qj - qi)
        if np.isinf(self.outboard_fringe):
            # C3 sentinel: the tip clip is REMOVED -- the sheet runs
            # outboard to the far field (the downstream and inboard tests
            # below are unaffected).
            on_sheet = crossing & (d_cross > 0.0)
        else:
            on_sheet = (
                crossing
                & (d_cross > 0.0)
                & (q_cross <= self.span_length + self.outboard_fringe)
            )
        if inboard_clip is None:
            # Legacy inboard clip: the sheet starts at the q = 0 TE end.
            on_sheet = on_sheet & (q_cross >= 0.0)
        else:
            # Caller-supplied inboard clip (B25 fragment topology),
            # evaluated at the s = 0 crossing POINTS.
            xyz_e = nodes[el[cand]]                     # (n_cand, 4, 3)
            pi = xyz_e[:, _TET_EDGES[:, 0], :]
            pj = xyz_e[:, _TET_EDGES[:, 1], :]
            p_cross = pi + t[:, :, None] * (pj - pi)
            inb = np.asarray(inboard_clip(p_cross), dtype=bool)
            assert inb.shape == crossing.shape, (
                f"inboard_clip returned {inb.shape}, expected "
                f"{crossing.shape}")
            on_sheet = on_sheet & inb
        self.cut_elems = cand[on_sheet.any(axis=1)]
        # Candidates rejected purely by the spanwise clip (the beyond-tip
        # wake-plane extension): reported, not cut. With an inboard_clip
        # this set also collects the inboard-rejected elements; with an
        # outboard_fringe (B31/C1) "the clip" is the EXTENDED one, so only
        # the still-further-outboard extension is reported here.
        downstream = crossing & (d_cross > 0.0)
        self.beyond_tip_elems = cand[
            downstream.any(axis=1) & ~on_sheet.any(axis=1)
        ]

        # Below-TE fan: contains a TE node, all non-TE nodes "-". These sit
        # entirely at/below the wake plane and only TOUCH the TE vertex -- the
        # sheet does not pass through their interior, so they are LOWER-side
        # elements (their TE reference is the aux DOF; Lopez fig. 3.6c), NOT
        # cut elements.
        #
        # ★ TE eps-shift exception (Lopez section 3.5.4, p.57: "Moving these
        # nodes upwards would result in cutting all wake elements under the
        # wake that are in touch with the trailing edge, which yields
        # inaccurate results"). The TE node is on the sheet, so the eps shift
        # sends it "+", which manufactures a sign change against the fan's "-"
        # nodes and a marginally-downstream crossing (d_cross -> 0+). Left
        # alone, those elements are classified CUT and handed a spurious UPPER
        # copy BELOW the wake, right at the TE -- which corrupts the doubled-TE
        # mass conservation that the implicit Kutta depends on. Measured on
        # NACA0012 a=2: 3 of the 6 below-TE fan elements were being cut, and
        # the emergent circulation OVER-shot by ~45% (Gamma 0.186 vs the
        # conforming/thin-airfoil 0.120), mesh-converged. So the fan is
        # subtracted from the cut set here, BEFORE the aux DOFs are numbered.
        is_te = np.zeros(n_nodes, dtype=bool)
        is_te[self.te_nodes] = True
        te_e = is_te[el]                                 # (n_tets, 4)
        has_te = te_e.any(axis=1)
        all_nonte_lower = np.all(te_e | (side_e == -1), axis=1)
        fan_mask = has_te & all_nonte_lower
        self.te_lower_elems = np.flatnonzero(fan_mask)
        self.cut_elems = self.cut_elems[~fan_mask[self.cut_elems]]

        # Per-node aux DOF over the unique cut-element nodes.
        cut_nodes = np.unique(el[self.cut_elems])
        self.n_ext_dofs = len(cut_nodes)
        self.n_total_dofs = self.n_main + self.n_ext_dofs
        self.ext_dof_of_node = np.full(n_nodes, -1, dtype=np.int64)
        self.ext_dof_of_node[cut_nodes] = self.n_main + np.arange(
            self.n_ext_dofs, dtype=np.int64
        )

        # A3 (C7): every TE node MUST own an aux DOF. If one is not a node of
        # any cut element, `MultivaluedOperator`'s pressure-Kutta reroute
        # writes reroute[-1] -- silently redirecting the LAST aux DOF's
        # mass-conservation row onto a TE main row and leaving that aux DOF
        # with no equation at all, and te_kutta_coo emits row = -1 triplets.
        # B1's gate locks this invariant on the committed meshes; assert it
        # here so a new mesh fails loudly at construction instead of
        # producing a quietly wrong system.
        # Scoped to meshes that actually carry a wake: with no cut elements
        # there are no aux DOFs at all and no multivalued system is built
        # (the degenerate classification-only case, e.g. a lone TE-fan tet).
        if len(self.te_nodes) > 0 and self.n_ext_dofs > 0:
            assert (self.ext_dof_of_node[self.te_nodes] >= 0).all(), (
                "TE node without an aux DOF: some TE node is not a node of "
                "any cut element (wake sheet does not reach the TE?)")

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
            "n_beyond_tip_elems": int(len(self.beyond_tip_elems)),
            "n_shifted_nodes": int(len(self.shifted_nodes)),
            "n_side_plus": int(np.sum(self.node_side == 1)),
            "n_side_minus": int(np.sum(self.node_side == -1)),
        }
