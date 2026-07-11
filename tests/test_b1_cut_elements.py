"""
Track B / B1 gate: level-set wake + cut-element identification
(docs/roadmap.md Track B B1; docs/design_track_b.md sections 3, 5.7).

Dual-mesh rule (roadmap Track B working rules):

  (a) wake-embedded M0 mesh -- every wake-plane node lies exactly on the
      level set, so the eps side-shift (D4) is stress-tested at scale, and
      the census can be cross-validated EXACTLY against the conforming
      preprocessor: the union of cut elements and the below-TE fan must
      equal the sheet's minus-side element star from cut_wake().

  (b) wake-free M3 mesh (cases/meshes/naca0012_wakefree_2.5d) -- generic
      cuts through generic elements, the actual Track B workflow form.
      No conforming counterpart exists: acceptance = census/corridor
      consistency + side classification against raw geometry + an
      alpha-sweep re-aim (update_direction) that never touches the mesh.

Synthetic single-tet and two-segment-TE cases pin the classification
rules (downstream test, TE fan, per-segment frames) independently of any
mesh file.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.wake import CutElementMap, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
M3_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_wakefree_2.5d"


def _wall_nodes(mesh):
    return np.unique(mesh.boundary_faces["wall"])


def _levelset_for(mesh, alpha_deg: float = 0.0) -> WakeLevelSet:
    """Straight wake from the quasi-2D TE segment at incidence alpha."""
    z = mesh.nodes[:, 2]
    a = np.radians(alpha_deg)
    te = np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]])
    return WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))


# ---------------------------------------------------------------------------
# Synthetic pins (no mesh files)
# ---------------------------------------------------------------------------

class TestSyntheticClassification:
    def _strip(self, x0):
        """Two-tet strip around the wake plane y=0 at x offset x0 from the
        TE point (1,0,0): nodes above and below, wake along +x."""
        nodes = np.array([
            [1.0 + x0, -0.5, 0.0],
            [1.0 + x0 + 1.0, -0.5, 0.0],
            [1.0 + x0, 0.5, 0.0],
            [1.0 + x0, 0.0, 1.0],
            [1.0 + x0 + 1.0, 0.5, 1.0],
        ])
        elements = np.array([[0, 1, 2, 3], [1, 2, 3, 4]])
        return nodes, elements

    def test_downstream_element_is_cut(self):
        nodes, elements = self._strip(x0=0.5)
        wls = WakeLevelSet([1.0, 0.0, 0.0], direction=(1.0, 0.0, 0.0))
        cm = CutElementMap(nodes, elements, wls)
        assert set(cm.cut_elems) == {0, 1}
        assert cm.n_ext_dofs == 5
        # side marking straight from geometry (no on-plane nodes here
        # except node 3, which the shift sends "+")
        assert cm.node_side[0] == -1 and cm.node_side[2] == 1
        assert cm.node_side[3] == 1  # on-plane -> "+" (D4)

    def test_upstream_element_not_cut(self):
        """Same strip moved AHEAD of the TE: sign change but every s=0
        crossing has d < 0 -- ahead-of-LE region, must NOT be cut."""
        nodes, elements = self._strip(x0=-5.0)
        wls = WakeLevelSet([1.0, 0.0, 0.0], direction=(1.0, 0.0, 0.0))
        cm = CutElementMap(nodes, elements, wls)
        assert len(cm.cut_elems) == 0
        assert cm.n_ext_dofs == 0

    def test_te_fan_split(self):
        """A tet with the TE node and all other nodes below goes to
        te_lower_elems, not cut_elems (Lopez fig. 3.6c)."""
        nodes = np.array([
            [1.0, 0.0, 0.0],       # the TE point itself
            [0.5, -0.3, 0.0],
            [0.9, -0.4, 1.0],
            [0.4, -0.2, 0.5],
        ])
        elements = np.array([[0, 1, 2, 3]])
        wls = WakeLevelSet([1.0, 0.0, 0.0], direction=(1.0, 0.0, 0.0))
        cm = CutElementMap(nodes, elements, wls)
        assert len(cm.cut_elems) == 0
        assert list(cm.te_lower_elems) == [0]
        assert list(cm.te_nodes) == [0]

    def test_dof_assignment_convention(self):
        """dofs_upper: main at '+', aux at '-'; dofs_lower complement."""
        nodes, elements = self._strip(x0=0.5)
        wls = WakeLevelSet([1.0, 0.0, 0.0], direction=(1.0, 0.0, 0.0))
        cm = CutElementMap(nodes, elements, wls)
        for k, e in enumerate(cm.cut_elems):
            conn = elements[e]
            for loc, n in enumerate(conn):
                if cm.node_side[n] == 1:
                    assert cm.dofs_upper[k, loc] == n
                    assert cm.dofs_lower[k, loc] == cm.ext_dof_of_node[n]
                else:
                    assert cm.dofs_upper[k, loc] == cm.ext_dof_of_node[n]
                    assert cm.dofs_lower[k, loc] == n

    def test_two_segment_te_frames(self):
        """Swept two-segment TE (D9): each panel classifies in its own
        frame; a point above panel 2's tilted surface must read s > 0
        even where panel 1's frame would disagree."""
        te = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 1.0], [1.2, 0.2, 2.0]])
        wls = WakeLevelSet(te, direction=(1.0, 0.0, 0.0))
        # Above panel 1 (z=0.5): plain +y offset.
        s1, d1 = wls.evaluate(np.array([[2.0, 0.1, 0.5]]))
        assert s1[0] > 0 and d1[0] > 0
        # Panel 2 tilts: its ruled surface through (1.1, 0.1, 1.5) --
        # a point slightly BELOW that local surface must read s < 0
        # even though its absolute y is positive.
        s2, _ = wls.evaluate(np.array([[2.0, 0.05, 1.5]]))
        assert s2[0] < 0

    def test_update_direction_reaims_without_mesh(self):
        nodes, elements = self._strip(x0=0.5)
        wls = WakeLevelSet([1.0, 0.0, 0.0], direction=(1.0, 0.0, 0.0))
        s0, _ = wls.evaluate(nodes)
        wls.update_direction((np.cos(0.3), np.sin(0.3), 0.0))
        s1, _ = wls.evaluate(nodes)
        assert not np.allclose(s0, s1)


# ---------------------------------------------------------------------------
# (a) wake-embedded M0 mesh: eps stress test + exact conforming cross-check
# ---------------------------------------------------------------------------

# B1 closure needs the medium gate green too (roadmap working rule 0);
# both levels are committed .msh for both families.
@pytest.fixture(scope="module", params=["coarse", "medium"])
def m0(request):
    mesh = read_mesh(M0_DIR / f"{request.param}.msh")
    mesh_cut, wc = cut_wake(mesh)
    # cut_wake returns a new mesh; re-read the pristine one for the map
    mesh_orig = read_mesh(M0_DIR / f"{request.param}.msh")
    wls = _levelset_for(mesh_orig)
    cm = CutElementMap(mesh_orig.nodes, mesh_orig.elements, wls,
                       wall_nodes=_wall_nodes(mesh_orig))
    return mesh_orig, mesh_cut, wc, cm


# M3 medium/fine .msh are gitignored (M6 pattern; ~40 s to regenerate:
# python cases/meshes/naca0012_wakefree_2.5d/generate_naca0012_wakefree.py);
# coarse is committed so the always-on B1 checks never skip.
@pytest.fixture(scope="module", params=["coarse", "medium"])
def m3(request):
    path = M3_DIR / f"{request.param}.msh"
    if not path.exists():
        pytest.skip(f"{path.name} not generated (gitignored; see fixture note)")
    return read_mesh(path)


class TestM0EmbeddedMesh:
    def test_eps_shift_covers_all_sheet_nodes(self, m0):
        """Every conforming sheet node lies ON the level set and must have
        been eps-shifted to '+' -- the D4 stress test at scale."""
        _, _, wc, cm = m0
        sheet = wc.master_nodes  # original ids of every duplicated node
        assert len(sheet) > 0
        shifted = set(cm.shifted_nodes)
        assert set(sheet) <= shifted
        assert np.all(cm.node_side[sheet] == 1)

    def test_census_equals_conforming_minus_star(self, m0):
        """cut_elems | te_lower_elems == elements that kept a master
        (minus-side) sheet-node reference through cut_wake -- an exact,
        element-by-element cross-validation against the shipped
        preprocessor's side classification."""
        mesh_orig, mesh_cut, wc, cm = m0
        el_orig = mesh_orig.elements.astype(np.int64)
        el_cut = mesh_cut.elements.astype(np.int64)
        is_master = np.zeros(len(mesh_orig.nodes), dtype=bool)
        is_master[wc.master_nodes] = True
        kept = is_master[el_orig] & (el_cut == el_orig)
        expected = set(np.flatnonzero(kept.any(axis=1)))
        computed = set(cm.cut_elems) | set(cm.te_lower_elems)
        assert computed == expected

    def test_te_nodes_match_conforming(self, m0):
        _, _, wc, cm = m0
        assert set(cm.te_nodes) == set(wc.te_nodes)

    def test_side_marking_matches_geometry(self, m0):
        """Off-surface nodes: side == sign(y) at alpha = 0."""
        mesh_orig, _, _, cm = m0
        y = mesh_orig.nodes[:, 1]
        off = np.abs(cm.s_raw) >= 1e-6  # exclude the shifted band
        assert np.all(cm.node_side[off] == np.sign(y[off]).astype(np.int8))

    def test_ext_dof_count_and_range(self, m0):
        mesh_orig, _, _, cm = m0
        cut_nodes = np.unique(
            mesh_orig.elements.astype(np.int64)[cm.cut_elems])
        assert cm.n_ext_dofs == len(cut_nodes)
        ext = cm.ext_dof_of_node[cut_nodes]
        assert ext.min() == cm.n_main
        assert ext.max() == cm.n_main + cm.n_ext_dofs - 1
        assert len(np.unique(ext)) == cm.n_ext_dofs


# ---------------------------------------------------------------------------
# (b) wake-free M3 mesh: generic cuts, corridor coverage, alpha re-aim
# ---------------------------------------------------------------------------

class TestM3WakeFreeMesh:
    def test_m3_ingest_without_wake_tag(self, m3):
        """M3 gate clause: the preprocessor path must not require a wake
        tag."""
        assert "wake" not in m3.boundary_faces
        assert {"wall", "farfield", "symmetry"} <= set(m3.boundary_faces)
        assert len(m3.elements) > 0

    def test_m3_corridor_sizing_evidence(self):
        """M3 gate clause: corridor sizing vs M0's sheet-adjacent h_wake
        (committed stats CSV is the evidence)."""
        import csv
        stats = {}
        with open(M3_DIR / "coarse_stats.csv") as f:
            for row in csv.reader(f):
                stats[row[0]] = row[1]
        h_wake = float(stats["h_wake_target"])
        median = float(stats["corridor_median_edge"])
        assert median <= 1.25 * h_wake
        assert int(stats["corridor_n_triangles"]) > 50

    @pytest.mark.parametrize("alpha", [0.0, 4.0])
    def test_census_and_corridor_coverage(self, m3, alpha):
        """Generic cuts: a connected cut corridor must run from the TE to
        the far field, at both alpha = 0 and a re-aimed alpha = 4 deg on
        the SAME mesh (the workflow payoff)."""
        wls = _levelset_for(m3, alpha_deg=alpha)
        cm = CutElementMap(m3.nodes, m3.elements, wls,
                           wall_nodes=_wall_nodes(m3))
        assert len(cm.cut_elems) > 0
        assert np.all(cm.node_side != 0)
        # ext DOFs = unique cut nodes
        cut_nodes = np.unique(m3.elements.astype(np.int64)[cm.cut_elems])
        assert cm.n_ext_dofs == len(cut_nodes)
        # corridor coverage: bin cut-element centroids by downstream
        # distance; every bin from TE to the far field must be populated
        cen = m3.nodes[m3.elements[cm.cut_elems]].mean(axis=1)
        a = np.radians(alpha)
        d_cen = (cen - np.array([1.0, 0.0, 0.0])) @ np.array(
            [np.cos(a), np.sin(a), 0.0])
        assert d_cen.min() < 0.1          # reaches the TE
        assert d_cen.max() > 13.0         # reaches the far field (~14.5)
        hist, _ = np.histogram(d_cen, bins=40, range=(0.0, 13.5))
        assert np.all(hist > 0), "cut corridor has gaps"
        # TE exists on the wall and the below-TE fan is present
        assert len(cm.te_nodes) >= 2      # both z-planes of the layer
        assert len(cm.te_lower_elems) > 0

    def test_alpha_reaim_changes_cut_set(self, m3):
        wls = _levelset_for(m3, alpha_deg=0.0)
        cm0 = CutElementMap(m3.nodes, m3.elements, wls,
                            wall_nodes=_wall_nodes(m3))
        wls.update_direction((np.cos(np.radians(4.0)),
                              np.sin(np.radians(4.0)), 0.0))
        cm4 = CutElementMap(m3.nodes, m3.elements, wls,
                            wall_nodes=_wall_nodes(m3))
        set0, set4 = set(cm0.cut_elems), set(cm4.cut_elems)
        assert set0 != set4
        # both corridors share the TE neighborhood but diverge downstream
        assert len(set0 & set4) > 0
