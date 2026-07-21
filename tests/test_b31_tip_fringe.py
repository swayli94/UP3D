"""
Track B / B31 candidate C1 (GB31.3): outboard faded fringe of the LS
wake-sheet termination.

docs/roadmap.md Track B B31; cases/analysis/b31_tip_termination/
PRE_REGISTRATION.md section 3 (GB31.3).

The embedded sheet ends at the spanwise clip as a HEAVISIDE (F1): the last
cut ring carries a finite jump |delta| ~ 0.026 that the single-valued
beyond-tip elements drop over one cell -- the tip-edge singularity
(honest p = +0.62). C1 extends the clip OUTBOARD by a fringe width w
(`CutElementMap(outboard_fringe=w)`) and fades the aux coupling 1 -> 0
across the fringe by re-aiming the B8 span_blend taper at the NEW sheet
end (z_tip = span_length + w, r_blend = w), so the sheet terminates at a
pure-weld edge. Mechanism note (pre-reg U2): the aux rows are
upwind-convective along +x, so an OUTBOARD weld cannot re-level the
inboard lifting solution -- the mechanical distinction from the B8
inboard span_blend dead end (welded the lifting sheet, -20% global cl).

These fast tests pin, on a dependency-free structured-cube mesh whose
sheet tip sits MID-DOMAIN (tests/mesh_utils.py):

  (a) default bit-identity: outboard_fringe omitted == 0.0, classification
      bitwise; span_blend=None untouched;
  (b) census response: fringe = k * h grows cut_elems monotonically
      (superset), shrinks beyond_tip_elems, and the newly cut elements
      reach beyond the (physical) tip;
  (c) fade algebra: w == 1 inboard (q <= span_length), strictly < 1 inside
      the fringe, 0 beyond the extended clip; r_blend > fringe violates
      the inboard invariant and must raise; the legacy (no-fringe) blend
      call is pinned bitwise;

plus, on the gitignored M6 coarse mesh (skip-guarded, the real swept-TE
production geometry), the census response and the A3 TE/aux invariant
(cut_elements.py: every TE node owns an aux DOF) under the fringe.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.constraints.wake import tip_taper_factors
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

from .mesh_utils import generate_structured_cube_mesh

REPO_ROOT = Path(__file__).parent.parent
M6_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6"
FORM = "vanish_smooth"


# ---------------------------------------------------------------------------
# Synthetic cube with a mid-domain sheet tip (no mesh files, never skips)
# ---------------------------------------------------------------------------

def _cube_tip_case(n=6, L=2.0, z0=0.5, z1=1.5):
    """(nodes, elements, wls, h): structured n^3 cube with a wake sheet
    ending MID-DOMAIN.

    TE = straight segment at x = 0.9, y = 1.0, z in [z0, z1]; the sheet is
    the y = 1 half-plane x > 0.9 with q = z - z0. The cube spans z in
    [0, 2], so with the default z0 = 0.5 BOTH plane extensions (inboard
    q < 0 and beyond-tip q > span) live inside the domain, while z0 = 0.0
    puts the inboard sheet end ON the domain boundary (the q >= 0 clip
    inert -- the C3 inf-clip test mesh). x = 0.9 is off-grid, so no node
    lies on the TE line (te_nodes empty here; the TE-side invariants are
    exercised on M6 below). Nodes ON y = 1 read s = 0 and get the D4 eps
    shift, so the s = 0 crossings land at the y = 1 grid planes.
    """
    nodes, elements = generate_structured_cube_mesh(n, L)
    te = np.array([[0.9, 1.0, z0], [0.9, 1.0, z1]])
    wls = WakeLevelSet(te, direction=(1.0, 0.0, 0.0))
    return nodes, elements, wls, L / n


def _nonte_nodes(cm):
    cut_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
    is_te = np.zeros(cm.n_main, dtype=bool)
    is_te[cm.te_nodes] = True
    return cut_nodes[~is_te[cut_nodes]]


class TestFringeDefaultBitIdentical:
    """(a) outboard_fringe=0.0 (or omitted) reproduces the legacy
    classification bit-for-bit; span_blend=None stays untouched."""

    def test_classification_bitwise(self):
        nodes, elements, wls, _ = _cube_tip_case()
        cm_def = CutElementMap(nodes, elements, wls)
        cm_zero = CutElementMap(nodes, elements, wls, outboard_fringe=0.0)
        assert cm_def.outboard_fringe == 0.0
        assert cm_zero.outboard_fringe == 0.0
        assert cm_def.summary() == cm_zero.summary()
        for attr in ("cut_elems", "beyond_tip_elems", "te_lower_elems",
                     "node_side", "ext_dof_of_node", "dofs_upper",
                     "dofs_lower"):
            assert np.array_equal(getattr(cm_def, attr),
                                  getattr(cm_zero, attr)), attr

    def test_span_blend_none_untouched(self):
        """A fringe map assembled WITHOUT span_blend keeps the raw wake-LS
        block (no blend, n_span_blended == 0) -- the default-off path."""
        nodes, elements, wls, h = _cube_tip_case()
        cm = CutElementMap(nodes, elements, wls, outboard_fringe=h)
        m_def = MultivaluedOperator(nodes, elements, cm, levelset=wls)
        m_none = MultivaluedOperator(nodes, elements, cm, levelset=wls,
                                     span_blend=None)
        assert m_def.n_span_blended == 0
        assert all(np.array_equal(a, b)
                   for a, b in zip(m_def._ls_coo, m_none._ls_coo))

    def test_negative_fringe_rejected(self):
        nodes, elements, wls, h = _cube_tip_case()
        with pytest.raises(ValueError, match="outboard_fringe"):
            CutElementMap(nodes, elements, wls, outboard_fringe=-h)


class TestFringeCensus:
    """(b) census response to the fringe width."""

    def _maps_with_h(self):
        """(cm0, cm1, cm2) at fringes 0, h, 2h on the default cube."""
        nodes, elements, wls, h = _cube_tip_case()
        return [CutElementMap(nodes, elements, wls,
                              outboard_fringe=k * h)
                for k in (0.0, 1.0, 2.0)]

    def test_cut_grows_beyond_shrinks_monotonically(self):
        cm0, cm1, cm2 = self._maps_with_h()
        n0, n1, n2 = (len(cm.cut_elems) for cm in (cm0, cm1, cm2))
        assert n0 < n1 < n2, "fringe must strictly grow the cut set"
        b0, b1, b2 = (len(cm.beyond_tip_elems) for cm in (cm0, cm1, cm2))
        assert b0 > b1 > b2, "fringe must strictly shrink the rejected set"
        assert b2 > 0, "inboard-rejected elements keep the census non-blind"

    def test_fringe_cut_set_is_superset(self):
        cm0, cm1, cm2 = self._maps_with_h()
        assert set(cm0.cut_elems) <= set(cm1.cut_elems)
        assert set(cm1.cut_elems) <= set(cm2.cut_elems)

    def test_new_cut_elements_reach_beyond_tip(self):
        nodes, elements, wls, h = _cube_tip_case()
        el = np.asarray(elements, dtype=np.int64)
        cm0 = CutElementMap(nodes, elements, wls)
        cm1 = CutElementMap(nodes, elements, wls, outboard_fringe=h)
        new = np.array(sorted(set(cm1.cut_elems) - set(cm0.cut_elems)),
                       dtype=np.int64)
        assert len(new) > 0
        # q is linear along an edge, so an element whose every on-sheet
        # crossing was rejected by the LEGACY clip must reach past the tip
        # with at least one node (the b1 beyond-tip invariant).
        q_max = cm1.q[el[new]].max(axis=1)
        assert np.all(q_max > cm1.span_length)
        # and each new element contains an on-sheet crossing inside the
        # EXTENDED clip, so its nearest node cannot lie past it either
        # (q_cross is a convex combination of nodal q along the edge).
        q_min = cm1.q[el[new]].min(axis=1)
        assert np.all(q_min <= cm1.span_length + h)


class TestFringeFade:
    """(c) fade-factor algebra of the fringe-aimed span_blend."""

    def _build(self, fringe_k=1.0, r_blend_k=1.0):
        nodes, elements, wls, h = _cube_tip_case()
        cm = CutElementMap(nodes, elements, wls,
                           outboard_fringe=fringe_k * h)
        base = MultivaluedOperator(nodes, elements, cm, levelset=wls)
        blend = MultivaluedOperator(nodes, elements, cm, levelset=wls,
                                    span_blend=(FORM, r_blend_k * h))
        return cm, h, base, blend

    def test_fade_factor_profile(self):
        cm, h, _, blend = self._build()
        nonte = _nonte_nodes(cm)
        w = tip_taper_factors(cm.q[nonte], cm.span_length + h, FORM, h)
        inb = cm.q[nonte] <= cm.span_length
        assert np.any(inb) and np.all(w[inb] == 1.0), \
            "inboard (q <= span_length) rows must be untouched"
        inside = (cm.q[nonte] > cm.span_length) \
            & (cm.q[nonte] < cm.span_length + h)
        assert np.any(inside), "cube must place cut nodes inside the fringe"
        assert np.all((w[inside] > 0.0) & (w[inside] < 1.0))
        out = cm.q[nonte] > cm.span_length + h
        if np.any(out):
            assert np.all(w[out] == 0.0), \
                "nodes past the extended clip are pure-welded"
        assert blend.n_span_blended == int(np.count_nonzero(w < 1.0))
        assert blend.n_span_blended > 0

    def test_ls_rows_scaled_and_weld_tail(self):
        """Same algebra as the legacy blend (test_b8_span_blend): LS rows
        scaled by w_j, weld tail (1 - w_j) * s_j with s_j the row's
        UNBLENDED LS magnitude -- only z_tip moved."""
        cm, h, base, blend = self._build()
        nonte = _nonte_nodes(cm)
        w = tip_taper_factors(cm.q[nonte], cm.span_length + h, FORM, h)
        w_of_row = np.ones(base.n_total)
        w_of_row[cm.ext_dof_of_node[nonte]] = w
        r0, c0, d0 = base._ls_coo
        r1, c1, d1 = blend._ls_coo
        n0 = len(r0)
        assert np.array_equal(r0, r1[:n0]) and np.array_equal(c0, c1[:n0])
        assert np.array_equal(d1[:n0], d0 * w_of_row[r0])
        # weld tail blocks (+coef at aux, -coef at main)
        s_of_row = np.zeros(base.n_total)
        np.add.at(s_of_row, r0, np.abs(d0))
        s_of_row *= 0.5
        sel = w < 1.0
        n_w = int(sel.sum())
        aux_b = cm.ext_dof_of_node[nonte[sel]]
        assert np.array_equal(r1[-2 * n_w:-n_w], aux_b)
        assert np.array_equal(r1[-n_w:], aux_b)
        assert np.array_equal(c1[-n_w:], nonte[sel])
        assert np.allclose(d1[-2 * n_w:-n_w],
                           (1.0 - w[sel]) * s_of_row[aux_b], rtol=1e-14)
        assert np.allclose(d1[-n_w:], -(1.0 - w[sel]) * s_of_row[aux_b],
                           rtol=0)

    def test_inboard_rows_bitwise_in_assembled_matrix(self):
        """w = 1 (inboard) rows of the assembled system must be bitwise
        the unblended rows -- locality by construction."""
        cm, h, base, blend = self._build()
        nonte = _nonte_nodes(cm)
        w = tip_taper_factors(cm.q[nonte], cm.span_length + h, FORM, h)
        aux_keep = cm.ext_dof_of_node[nonte[w == 1.0]]
        rng = np.random.default_rng(0)
        phi = rng.standard_normal(base.n_total)
        A0 = base.assemble_matrix(closure="wake_ls", phi_ext=phi).tocsr()
        A1 = blend.assemble_matrix(closure="wake_ls", phi_ext=phi).tocsr()
        for r in aux_keep[:: max(1, len(aux_keep) // 50)]:
            s0, e0 = A0.indptr[r], A0.indptr[r + 1]
            s1, e1 = A1.indptr[r], A1.indptr[r + 1]
            assert np.array_equal(A0.indices[s0:e0], A1.indices[s1:e1])
            assert np.array_equal(A0.data[s0:e0], A1.data[s1:e1])

    def test_legacy_blend_call_bitwise_without_fringe(self):
        """No-fringe map: the blend must use the LEGACY z_tip = span_length
        call, bit-identical (the B8 semantics, incl. beyond-tip straddler
        nodes welded at w = 0)."""
        nodes, elements, wls, h = _cube_tip_case()
        cm = CutElementMap(nodes, elements, wls)
        base = MultivaluedOperator(nodes, elements, cm, levelset=wls)
        blend = MultivaluedOperator(nodes, elements, cm, levelset=wls,
                                    span_blend=(FORM, h))
        nonte = _nonte_nodes(cm)
        w = tip_taper_factors(cm.q[nonte], cm.span_length, FORM, h)
        w_of_row = np.ones(base.n_total)
        w_of_row[cm.ext_dof_of_node[nonte]] = w
        r0, c0, d0 = base._ls_coo
        r1, c1, d1 = blend._ls_coo
        n0 = len(r0)
        assert np.array_equal(r0, r1[:n0]) and np.array_equal(c0, c1[:n0])
        assert np.array_equal(d1[:n0], d0 * w_of_row[r0])
        assert blend.n_span_blended == int(np.count_nonzero(w < 1.0))

    def test_rblend_larger_than_fringe_raises(self):
        """r_blend > fringe would fade INBOARD rows -- the C1 inboard
        invariant must reject it loudly at construction."""
        nodes, elements, wls, h = _cube_tip_case()
        cm = CutElementMap(nodes, elements, wls, outboard_fringe=h)
        with pytest.raises(AssertionError, match="inboard"):
            MultivaluedOperator(nodes, elements, cm, levelset=wls,
                                span_blend=(FORM, 2.0 * h))


# ---------------------------------------------------------------------------
# C3 sentinel: outboard_fringe = np.inf drops the tip clip entirely (the
# sheet runs outboard to the far field; downstream/inboard tests stay).
# ---------------------------------------------------------------------------

class TestNoClipInf:
    def _maps(self, z0=0.0):
        nodes, elements, wls, h = _cube_tip_case(z0=z0)
        cms = [CutElementMap(nodes, elements, wls, outboard_fringe=k * h)
               for k in (0.0, 1.0, 2.0)]
        cm_inf = CutElementMap(nodes, elements, wls,
                               outboard_fringe=np.inf)
        return cms, cm_inf, np.asarray(elements, dtype=np.int64)

    def test_beyond_tip_empty_when_plane_reaches_boundary(self):
        """z0 = 0 puts the inboard end on the domain boundary, so no
        inboard rejection can fire either: with the tip clip dropped,
        NOTHING downstream is rejected -- beyond_tip_elems is empty."""
        _, cm_inf, _ = self._maps()
        assert len(cm_inf.cut_elems) > 0
        assert len(cm_inf.beyond_tip_elems) == 0

    def test_cut_set_strict_superset_of_every_fringe(self):
        (cm0, cm1, cm2), cm_inf, el = self._maps()
        for cm in (cm0, cm1, cm2):
            assert set(cm.cut_elems) <= set(cm_inf.cut_elems)
        assert set(cm0.cut_elems) < set(cm_inf.cut_elems)
        # the extension reaches genuinely outboard: elements cut by the
        # inf sentinel but NOT by the h-fringe must have an on-sheet
        # crossing past span + h, hence (q linear along the edge) a node
        # reaching there -- the far-field branch P5 audits, deliberately
        # admitted.
        new = np.array(sorted(set(cm_inf.cut_elems) - set(cm1.cut_elems)),
                       dtype=np.int64)
        assert len(new) > 0
        assert np.all(cm_inf.q[el[new]].max(axis=1)
                      > cm_inf.span_length + 1.0 / 6.0)
        # and the corridor is continuous: cut-element nodal q covers the
        # whole outboard range without gaps
        cen_q = cm_inf.q[el[cm_inf.cut_elems]].mean(axis=1)
        hist, _ = np.histogram(cen_q, bins=8,
                               range=(0.0, cen_q.max()))
        assert np.all(hist > 0), "no-clip corridor has spanwise gaps"

    def test_inf_span_blend_is_noop(self):
        """Documented interaction: with no free end there is nothing to
        fade -- z_tip = inf gives w == 1 everywhere, the blend is inert."""
        nodes, elements, wls, h = _cube_tip_case(z0=0.0)
        cm = CutElementMap(nodes, elements, wls, outboard_fringe=np.inf)
        base = MultivaluedOperator(nodes, elements, cm, levelset=wls)
        blend = MultivaluedOperator(nodes, elements, cm, levelset=wls,
                                    span_blend=(FORM, h))
        assert blend.n_span_blended == 0
        assert all(np.array_equal(a, b)
                   for a, b in zip(base._ls_coo, blend._ls_coo))

    def test_inboard_clip_still_applies_with_inf(self):
        """The default-z0 cube HAS an inboard (q < 0) plane extension:
        the legacy q >= 0 clip must still reject it under the inf
        sentinel (only the TIP test is dropped)."""
        (cm0, _, _), cm_inf, el = self._maps(z0=0.5)
        assert len(cm_inf.beyond_tip_elems) > 0, \
            "inboard-rejected elements keep the census non-blind"
        # every rejected element's on-plane crossings are all inboard
        q_max = cm_inf.q[el[cm_inf.beyond_tip_elems]].max(axis=1)
        assert np.all(q_max < 0.5)   # all nodes below the z0 = 0.5 TE end


# ---------------------------------------------------------------------------
# M6 coarse (gitignored): the real swept-TE production geometry -- census
# response + the A3 TE/aux invariant under the fringe.
# ---------------------------------------------------------------------------

def _m6_levelset(alpha_deg: float = 3.06) -> WakeLevelSet:
    from pyfp3d.meshgen.wing3d import B_SEMI, x_te

    a = np.radians(alpha_deg)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    return WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))


@pytest.fixture(scope="module")
def m6():
    path = M6_DIR / "coarse.msh"
    if not path.exists():
        pytest.skip(f"{path} missing (run cases/meshes/onera_m6/"
                    "generate_onera_m6.py)")
    return read_mesh(path)


class TestM6Fringe:
    W_FRAC = 0.05   # the C1 probe width: w = 0.05 * B_SEMI

    def _maps(self, mesh):
        from pyfp3d.meshgen.wing3d import B_SEMI

        wls = _m6_levelset()
        wall = np.unique(mesh.boundary_faces["wall"])
        cm0 = CutElementMap(mesh.nodes, mesh.elements, wls, wall_nodes=wall)
        cm1 = CutElementMap(mesh.nodes, mesh.elements, wls, wall_nodes=wall,
                            outboard_fringe=self.W_FRAC * B_SEMI)
        return wls, cm0, cm1

    def test_census_response(self, m6):
        from pyfp3d.meshgen.wing3d import B_SEMI

        _, cm0, cm1 = self._maps(m6)
        el = np.asarray(m6.elements, dtype=np.int64)
        assert len(cm1.cut_elems) > len(cm0.cut_elems)
        assert set(cm0.cut_elems) <= set(cm1.cut_elems)
        assert len(cm1.beyond_tip_elems) < len(cm0.beyond_tip_elems)
        new = np.array(sorted(set(cm1.cut_elems) - set(cm0.cut_elems)),
                       dtype=np.int64)
        assert np.all(cm1.q[el[new]].max(axis=1) > cm1.span_length)
        # nothing new may reach past the EXTENDED clip
        assert np.all(cm1.q[el[cm1.cut_elems]].min(axis=1)
                      <= cm1.span_length + self.W_FRAC * B_SEMI)

    def test_te_aux_invariant_holds_with_fringe(self, m6):
        """(d) cut_elements.py:262-265 under the fringe: every TE node
        still owns an aux DOF (the constructor asserts it; pin it here)."""
        _, _, cm1 = self._maps(m6)
        assert len(cm1.te_nodes) > 0
        assert np.all(cm1.ext_dof_of_node[cm1.te_nodes] >= 0)

    def test_noclip_census(self, m6):
        """C3 inf sentinel on the real geometry: the clip rejection
        vanishes entirely (beyond_tip_elems empty -- M6 has no inboard
        q < 0 extension), the cut set is a strict superset reaching the
        far-field boundary, and the A3 TE/aux invariant survives."""
        wls, cm0, _ = self._maps(m6)
        wall = np.unique(m6.boundary_faces["wall"])
        ci = CutElementMap(m6.nodes, m6.elements, wls, wall_nodes=wall,
                           outboard_fringe=np.inf)
        el = np.asarray(m6.elements, dtype=np.int64)
        assert set(cm0.cut_elems) < set(ci.cut_elems)
        assert len(ci.beyond_tip_elems) == 0
        ff = np.unique(m6.boundary_faces["farfield"])
        assert np.count_nonzero(np.isin(el[ci.cut_elems], ff).any(axis=1)) \
            > 0, "the no-clip extension must reach the far field"
        assert np.all(ci.ext_dof_of_node[ci.te_nodes] >= 0)

    def test_fade_on_real_geometry(self, m6):
        from pyfp3d.meshgen.wing3d import B_SEMI

        wls, _, cm1 = self._maps(m6)
        w = self.W_FRAC * B_SEMI
        blend = MultivaluedOperator(m6.nodes, m6.elements, cm1, levelset=wls,
                                    span_blend=(FORM, w))
        assert blend.n_span_blended > 0
        nonte = _nonte_nodes(cm1)
        w_exp = tip_taper_factors(cm1.q[nonte], cm1.span_length + w,
                                  FORM, w)
        inb = cm1.q[nonte] <= cm1.span_length
        assert np.all(w_exp[inb] == 1.0)
        inside = (cm1.q[nonte] > cm1.span_length) \
            & (cm1.q[nonte] < cm1.span_length + w)
        assert np.any(inside), "M6 coarse must place nodes in the fringe"
        assert np.all(w_exp[inside] < 1.0)
