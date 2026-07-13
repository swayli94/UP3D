"""
Track B / B8 re-spec: spanwise termination blend of the wake LS rows.

docs/roadmap.md Track B B8 (re-spec); run_b8_termination_diagnosis.py.

The B8 row blend (TE rows) was measured NOT to bound the LS tip edge. The
termination diagnosis found why, and what actually diverges:

  * the committed p = +1.34 was 5x inflated by a metric artifact --
    `element_mach2` reads mixed-side PLAIN (beyond-tip) elements from the
    aux-substituted side field, but their assembly uses MAIN dofs (the
    honest exponent is +0.62, the same object as the conforming +0.52);
  * the REAL residual singularity is the sheet's spanwise TERMINATION: the
    last cut ring carries a FINITE jump (|delta| ~ 0.026, h- and
    taper-independent) that drops to zero across one single-valued element
    -- and the TE tip_taper cannot reach it (it welds TE nodes; this jump
    lives on the downstream termination ring).

The span blend replaces each near-tip wake-LS aux row by

    w_j * [wake LS row]  +  (1 - w_j) * s_j * [phi_aux - phi_main]

with w_j = tip_taper_factors(q_j, span_length, form, r_blend) and s_j the
row's own LS magnitude (the LS entries are O(h), a bare weld O(1): without
s_j the blend would drift toward the weld under refinement).

These fast tests pin the ALGEBRA (default-inert bitwise, convex-limit
behavior, weld scaling, beyond-tip nodes always welded) plus the honest
element_mach2 reading. The physical mechanism probe (honest off-body
tip-peak exponent under refinement, delta_ring decay, locality) lives in
the demo `cases/demo/b8_tip_taper_ls/run_b8_span_blend.py`.

All tests need the gitignored M6 coarse mesh (the spanwise termination does
not exist on the quasi-2D families) and skip when it is absent.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.constraints.wake import tip_taper_factors
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M6_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6"
ALPHA, M = 3.06, 0.5
FORM, R_BLEND_FRAC = "vanish_smooth", 0.05


@pytest.fixture(scope="module")
def m6():
    """(mesh, wls, cm, mvop_base, mvop_blend) on M6 coarse; module-cached."""
    path = M6_DIR / "coarse.msh"
    if not path.exists():
        pytest.skip(f"{path} missing (run cases/meshes/onera_m6/"
                    "generate_onera_m6.py)")
    from pyfp3d.meshgen.wing3d import B_SEMI, x_te

    mesh = read_mesh(path)
    a = np.radians(ALPHA)
    wls = WakeLevelSet(
        np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
        direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    base = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    blend = MultivaluedOperator(
        mesh.nodes, mesh.elements, cm, levelset=wls,
        span_blend=(FORM, R_BLEND_FRAC * B_SEMI))
    return mesh, wls, cm, base, blend


def _phi_seed(mvop):
    rng = np.random.default_rng(0)
    return rng.standard_normal(mvop.n_total)


def _nonte_nodes_and_w(cm, r_blend):
    cut_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
    is_te = np.zeros(cm.n_main, dtype=bool)
    is_te[cm.te_nodes] = True
    nonte = cut_nodes[~is_te[cut_nodes]]
    return nonte, tip_taper_factors(cm.q[nonte], cm.span_length, FORM, r_blend)


class TestDefaultInert:
    """span_blend=None must be bit-identical everywhere."""

    def test_ls_coo_bitwise(self, m6):
        mesh, wls, cm, base, _ = m6
        m_none = MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                     levelset=wls, span_blend=None)
        assert all(np.array_equal(a, b)
                   for a, b in zip(base._ls_coo, m_none._ls_coo))

    def test_assembled_matrix_bitwise(self, m6):
        mesh, wls, cm, base, _ = m6
        m_none = MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                     levelset=wls, span_blend=None)
        phi = _phi_seed(base)
        A0 = base.assemble_matrix(closure="wake_ls", phi_ext=phi)
        A1 = m_none.assemble_matrix(closure="wake_ls", phi_ext=phi)
        A0.sort_indices()
        A1.sort_indices()
        assert np.array_equal(A0.data, A1.data)
        assert np.array_equal(A0.indices, A1.indices)
        assert np.array_equal(A0.indptr, A1.indptr)


class TestBlendAlgebra:
    def test_blended_count_matches_weights(self, m6):
        _, _, cm, _, blend = m6
        from pyfp3d.meshgen.wing3d import B_SEMI

        _, w = _nonte_nodes_and_w(cm, R_BLEND_FRAC * B_SEMI)
        assert blend.n_span_blended == int(np.count_nonzero(w < 1.0))
        assert blend.n_span_blended > 0

    def test_beyond_tip_nodes_always_welded(self, m6):
        """Cut nodes OUTBOARD of the tip (q > span_length -- B1's straddler
        cells) get w = 0 for ANY r_blend: the sheet must not carry a jump
        beyond its own end. This is the semantic core of the re-spec."""
        _, _, cm, _, _ = m6
        from pyfp3d.meshgen.wing3d import B_SEMI

        nonte, w = _nonte_nodes_and_w(cm, R_BLEND_FRAC * B_SEMI)
        beyond = cm.q[nonte] > cm.span_length
        assert np.any(beyond), "M6 coarse must have straddler cut nodes"
        assert np.all(w[beyond] == 0.0)
        _, w_tiny = _nonte_nodes_and_w(cm, 1e-12)
        assert np.all(w_tiny[beyond] == 0.0)

    def test_inboard_rows_untouched(self, m6):
        """w = 1 rows must be bitwise the unblended LS rows (locality by
        construction)."""
        _, _, cm, base, blend = m6
        from pyfp3d.meshgen.wing3d import B_SEMI

        nonte, w = _nonte_nodes_and_w(cm, R_BLEND_FRAC * B_SEMI)
        aux_keep = cm.ext_dof_of_node[nonte[w == 1.0]]
        phi = _phi_seed(base)
        A0 = base.assemble_matrix(closure="wake_ls", phi_ext=phi).tocsr()
        A1 = blend.assemble_matrix(closure="wake_ls", phi_ext=phi).tocsr()
        for r in aux_keep[:: max(1, len(aux_keep) // 50)]:
            s0, e0 = A0.indptr[r], A0.indptr[r + 1]
            s1, e1 = A1.indptr[r], A1.indptr[r + 1]
            assert np.array_equal(A0.indices[s0:e0], A1.indices[s1:e1])
            assert np.array_equal(A0.data[s0:e0], A1.data[s1:e1])

    def test_weld_scaled_by_row_ls_magnitude(self, m6):
        """The weld coefficient must be (1 - w_j) * s_j with s_j the row's
        UNBLENDED LS magnitude -- the h-invariance normalization."""
        _, _, cm, base, blend = m6
        from pyfp3d.meshgen.wing3d import B_SEMI

        ls_row0, _, ls_data0 = base._ls_coo
        s_of_row = np.zeros(base.n_total)
        np.add.at(s_of_row, ls_row0, np.abs(ls_data0))
        s_of_row *= 0.5

        nonte, w = _nonte_nodes_and_w(cm, R_BLEND_FRAC * B_SEMI)
        sel = w < 1.0
        aux_b = cm.ext_dof_of_node[nonte[sel]]
        # the blended _ls_coo ends with the two weld blocks (+coef, -coef)
        n_w = int(sel.sum())
        rows, cols, data = blend._ls_coo
        wp, wm = data[-2 * n_w:-n_w], data[-n_w:]
        assert np.allclose(wp, (1.0 - w[sel]) * s_of_row[aux_b], rtol=1e-14)
        assert np.allclose(wm, -wp, rtol=0)
        assert np.array_equal(rows[-2 * n_w:-n_w], aux_b)
        assert np.array_equal(cols[-n_w:], nonte[sel])

    def test_ls_rows_scaled_by_w(self, m6):
        """The LS part of a blended row must be exactly w_j x the unblended
        row (checked on the raw COO, before duplicate summing)."""
        _, _, cm, base, blend = m6
        from pyfp3d.meshgen.wing3d import B_SEMI

        nonte, w = _nonte_nodes_and_w(cm, R_BLEND_FRAC * B_SEMI)
        w_of_row = np.ones(base.n_total)
        w_of_row[cm.ext_dof_of_node[nonte]] = w
        r0, c0, d0 = base._ls_coo
        r1, c1, d1 = blend._ls_coo
        n0 = len(r0)
        assert np.array_equal(r0, r1[:n0]) and np.array_equal(c0, c1[:n0])
        assert np.array_equal(d1[:n0], d0 * w_of_row[r0])


class TestHonestMach:
    """element_mach2(mixed_plain=...) -- the B8 metric fix (opt-in)."""

    def test_default_bitwise(self, m6):
        _, _, _, base, _ = m6
        phi = _phi_seed(base)
        a = base.element_mach2(phi, M)
        b = base.element_mach2(phi, M, mixed_plain="side")
        assert np.array_equal(a, b, equal_nan=True)

    def test_main_changes_only_mixed_plain(self, m6):
        mesh, _, cm, base, _ = m6
        phi = _phi_seed(base)
        a = base.element_mach2(phi, M)
        b = base.element_mach2(phi, M, mixed_plain="main")
        el = np.asarray(mesh.elements, dtype=np.int64)
        special = np.zeros(len(el), dtype=bool)
        special[cm.cut_elems] = True
        special[cm.te_lower_elems] = True
        side_e = cm.node_side[el]
        fix = ~special & (side_e.min(axis=1) != side_e.max(axis=1))
        diff = ~np.isclose(np.nan_to_num(a, nan=-1.0),
                           np.nan_to_num(b, nan=-1.0))
        assert np.any(diff), "the honest reading must change something"
        assert np.all(fix[diff])

    def test_main_equals_main_field_on_fixed_elements(self, m6):
        mesh, _, cm, base, _ = m6
        from pyfp3d.physics.isentropic import mach_squared_field

        phi = _phi_seed(base)
        b = base.element_mach2(phi, M, mixed_plain="main")
        el = np.asarray(mesh.elements, dtype=np.int64)
        special = np.zeros(len(el), dtype=bool)
        special[cm.cut_elems] = True
        special[cm.te_lower_elems] = True
        side_e = cm.node_side[el]
        fix = ~special & (side_e.min(axis=1) != side_e.max(axis=1))
        _, q2m = base.op.velocities(base.main_potential(phi))
        m2m = mach_squared_field(q2m, M)
        assert np.array_equal(np.nan_to_num(b[fix], nan=-1.0),
                              np.nan_to_num(m2m[fix], nan=-1.0))

    def test_bad_arg_raises(self, m6):
        _, _, _, base, _ = m6
        with pytest.raises(ValueError):
            base.element_mach2(_phi_seed(base), M, mixed_plain="junk")
