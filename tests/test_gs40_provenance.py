"""
GS4.0 locks: the honesty fields, the mesh manifest, and the node-list check.

These three cover the defects the 2026-08-16 independent audit found
(docs/inspection/20260816-2200-independent-audit-zh.md §7.1, §7.3, §7.6).
Every one of them is an INSTRUMENT defect -- code whose job is to tell the
truth about a measurement -- and all three had the same shape: the instrument
kept reporting after the thing it measured went away, and nothing went red.

★ Why these are ungated and cheap by construction: the failure modes worth
locking are an empty list, a ramp that died at level 0, a ramp that died
mid-way, and two meshes that differ only at the tip. Every one of them is
cheaper to BUILD than to provoke, and this season two reporting-layer defects
each destroyed an expensive solve precisely because they were only ever
exercised live. So none of this file solves anything.
"""

import json

import numpy as np
import pytest

from pyfp3d.mesh.manifest import (
    describe_difference,
    manifest_path,
    mesh_fingerprint,
    read_manifest,
    write_manifest,
)
from pyfp3d.mesh.reader import read_mesh, write_mesh
from pyfp3d.solve.newton import _ramp_honesty_fields


# --------------------------------------------------------------------------
# F1 -- which Mach does the returned state actually live at?
# --------------------------------------------------------------------------
def _lv(*pairs):
    return [{"m": m, "converged": c} for m, c in pairs]


class TestRampHonestyFields:
    """`solve_newton_transonic` returns the FAILED level's state when a ramp
    dies early, at a LOWER Mach than requested. Until GS4.0 the conforming
    driver never said so -- the fields existed only on the deleted level-set
    driver (`newton_ls.py:1186-1193` at d224223), so three bench consumers read
    them through `.get(..., <default>)` and silently got the REQUESTED Mach.
    """

    def test_empty_level_results_gives_none_not_an_exception(self):
        """A reporting layer that RAISES is how a 40-minute solve was lost.

        Unreachable today (the ramp loop always runs at least once), which is
        exactly why it is asserted rather than argued -- the next refactor does
        not have to rediscover the intent.
        """
        out = _ramp_honesty_fields([], 0.84, False)
        assert out["m_final"] is None
        assert out["m_last_converged"] is None
        assert out["target_reached"] is False

    def test_died_at_first_level_reports_no_converged_mach(self):
        out = _ramp_honesty_fields(_lv((0.70, False)), 0.84, False)
        assert out["m_final"] == pytest.approx(0.70)
        assert out["m_last_converged"] is None, (
            "nothing converged, so there is no converged Mach to report; "
            "reporting 0.70 here would be the original defect in a new place")
        assert out["target_reached"] is False

    def test_died_mid_ramp_separates_the_two_machs(self):
        out = _ramp_honesty_fields(
            _lv((0.70, True), (0.75, True), (0.80, False)), 0.84, False)
        assert out["m_final"] == pytest.approx(0.80), (
            "the RETURNED state is the failed level's, so the row's cl/M_max "
            "belong to 0.80")
        assert out["m_last_converged"] == pytest.approx(0.75)
        assert out["target_reached"] is False

    def test_inserted_retry_levels_do_not_break_the_derivation(self):
        """The recorded difference from the level-set original.

        That driver used `levels[i - 1]`. This ramp INSERTS retry levels into
        `levels` on failure (newton.py, the dm-halving branch), so the positional
        index stops meaning "the previous level". Derived from level_results
        instead; the sequence of CONVERGED levels stays monotone increasing, so
        "last converged" and "highest converged" remain the same level.
        """
        out = _ramp_honesty_fields(
            _lv((0.70, True), (0.80, False), (0.75, True), (0.80, False)),
            0.84, False)
        assert out["m_last_converged"] == pytest.approx(0.75), (
            "positional levels[i-1] would have picked the 0.80 failure here")
        assert out["m_final"] == pytest.approx(0.80)

    def test_target_reached_only_when_both_conditions_hold(self):
        reached = _ramp_honesty_fields(
            _lv((0.70, True), (0.84, True)), 0.84, True)
        assert reached["target_reached"] is True
        assert reached["m_final"] == pytest.approx(0.84)
        assert reached["m_last_converged"] == pytest.approx(0.84)

        #: an upwind_c_post leg runs AT m_inf and can fail there -- m_final is
        #: then the target while the state is not usable. Both clauses matter.
        post_failed = _ramp_honesty_fields(
            _lv((0.70, True), (0.84, True), (0.84, False)), 0.84, False)
        assert post_failed["m_final"] == pytest.approx(0.84)
        assert post_failed["target_reached"] is False, (
            "m_final == m_inf is not sufficient: the returned state did not "
            "converge")

    def test_types_are_what_a_reporting_layer_assumes(self):
        """Return ARITY and argument TYPES, not just the signature.

        A dry-check that read only argument NAMES passed lists where arrays
        were needed and unpacked 2 returns from a 4-return function, and it
        raised AFTER a completed M6 medium solve and BEFORE the row was
        appended. Same class, so it gets an assertion.
        """
        out = _ramp_honesty_fields(_lv((0.7, True), (0.8, False)), 0.84, False)
        assert set(out) == {"m_final", "m_last_converged", "target_reached"}
        assert isinstance(out["target_reached"], bool)
        for k in ("m_final", "m_last_converged"):
            assert out[k] is None or isinstance(out[k], float)


# --------------------------------------------------------------------------
# F6 -- mesh provenance
# --------------------------------------------------------------------------
@pytest.fixture
def small_mesh_path():
    from pathlib import Path
    p = (Path(__file__).parent.parent / "cases" / "meshes"
         / "sphere_shell" / "coarse.msh")
    if not p.exists():
        pytest.skip(f"{p} not present (regenerate with its generate_*.py)")
    return p


class TestMeshFingerprint:
    def test_fingerprint_is_deterministic(self, small_mesh_path):
        a = mesh_fingerprint(small_mesh_path)
        b = mesh_fingerprint(small_mesh_path)
        assert a == b, "same file, two calls -- nothing here may vary with time"
        assert a["n_tets"] > 0 and a["n_nodes"] > 0
        assert set(a["n_tris"]) == {"wall", "farfield"}

    def test_nothing_time_or_host_dependent_leaks_in(self, small_mesh_path):
        """Two byte-identical meshes written at different times must agree.

        A fingerprint containing mtime/hostname/git-SHA would compare unequal
        for identical meshes, i.e. it would answer a DIFFERENT question than
        the one asked ("is this the same mesh?").
        """
        fp = mesh_fingerprint(small_mesh_path)
        assert set(fp) == {"schema", "sha256", "bytes", "n_nodes", "n_tets",
                           "n_tris"}

    def test_one_changed_byte_changes_the_hash(self, small_mesh_path, tmp_path):
        raw = small_mesh_path.read_bytes()
        clone = tmp_path / "clone.msh"
        clone.write_bytes(raw)
        assert (mesh_fingerprint(clone)["sha256"]
                == mesh_fingerprint(small_mesh_path)["sha256"])

        #: flip a digit inside a coordinate; still a valid file, different mesh
        tampered = tmp_path / "tampered.msh"
        idx = raw.rindex(b"\n", 0, len(raw) // 2)
        tampered.write_bytes(raw[:idx] + b" " + raw[idx:])
        assert (mesh_fingerprint(tampered)["sha256"]
                != mesh_fingerprint(small_mesh_path)["sha256"])

    def test_missing_mesh_is_loud(self, tmp_path):
        """Quietly returning "unknown" would put this module in the same class
        as the defect it exists to prevent."""
        with pytest.raises(FileNotFoundError):
            mesh_fingerprint(tmp_path / "does_not_exist.msh")

    def test_manifest_path_keeps_the_full_filename(self):
        """`coarse.msh` and `coarse_flat.msh` are exactly the pair whose
        confusion motivated this module -- their sidecars must not collide."""
        a = manifest_path("cases/meshes/onera_m6/coarse.msh")
        b = manifest_path("cases/meshes/onera_m6/coarse_flat.msh")
        assert a != b
        assert a.name == "coarse.msh.manifest.json"


class TestWriteMeshEmitsManifest:
    def test_write_mesh_writes_a_matching_sidecar(self, small_mesh_path,
                                                  tmp_path):
        mesh = read_mesh(small_mesh_path)
        out = tmp_path / "written.msh"
        write_mesh(mesh, out)

        side = manifest_path(out)
        assert side.exists(), (
            "write_mesh must emit provenance unconditionally -- the 3-D meshes "
            "are gitignored, so opt-in provenance is no provenance")
        loaded = json.loads(side.read_text())
        assert loaded == mesh_fingerprint(out)
        assert loaded["n_tets"] == len(mesh.elements)
        assert loaded["n_nodes"] == len(mesh.nodes)

    def test_sidecar_does_not_disturb_the_mesh_bytes(self, small_mesh_path,
                                                     tmp_path):
        mesh = read_mesh(small_mesh_path)
        a, b = tmp_path / "a.msh", tmp_path / "b.msh"
        write_mesh(mesh, a)
        write_mesh(mesh, b)
        assert a.read_bytes() == b.read_bytes()
        assert read_mesh(a).elements.shape == mesh.elements.shape

    def test_extra_provenance_is_quarantined_from_the_hash(self, small_mesh_path,
                                                           tmp_path):
        mesh = read_mesh(small_mesh_path)
        out = tmp_path / "x.msh"
        write_mesh(mesh, out)
        plain = json.loads(manifest_path(out).read_text())

        write_manifest(out, generator="test", level="coarse",
                       counts=np.array([1, 2, 3]))
        stamped = json.loads(manifest_path(out).read_text())
        assert stamped["sha256"] == plain["sha256"], (
            "non-reproducible provenance must not enter the identity")
        assert stamped["extra"]["level"] == "coarse"
        assert stamped["extra"]["counts"] == [1, 2, 3], "ndarray must serialise"

    def test_read_manifest_is_none_when_absent(self, tmp_path):
        """The normal state for a mesh generated before GS4.0."""
        assert read_manifest(tmp_path / "nothing.msh") is None


class TestDescribeDifference:
    def test_same_hash_reads_as_same_mesh(self):
        fp = {"sha256": "abc", "n_tets": 10, "n_tris": {"wall": 3}}
        assert describe_difference(fp, dict(fp)).startswith("SAME MESH")

    def test_difference_names_the_field_that_moved(self):
        """The flat/round pair, in miniature: the tip is where they differ."""
        a = {"sha256": "a" * 64, "n_tets": 68624, "bytes": 1,
             "n_tris": {"wall": 6890}}
        b = {"sha256": "b" * 64, "n_tets": 55531, "bytes": 2,
             "n_tris": {"wall": 5152}}
        msg = describe_difference(a, b)
        assert msg.startswith("DIFFERENT MESH")
        assert "n_tets 68624 vs 55531" in msg
        assert "n_tris[wall] 6890 vs 5152" in msg


# --------------------------------------------------------------------------
# F3 -- the fast tier's exclusion list must not rot
# --------------------------------------------------------------------------
class TestCapabilityLockNodePaths:
    """The fast tier prints "what this tier does NOT cover" and the file's own
    header calls that list part of the OUTPUT, not a comment. It was never
    validated, and the audit found an entry pointing at a file archived in
    phase 3 -- untagged, while every sibling carried [ARCHIVED].
    """

    @staticmethod
    def _module():
        import importlib.util
        from pathlib import Path
        p = Path(__file__).parent.parent / "bench" / "run_capability_locks.py"
        spec = importlib.util.spec_from_file_location("_caplocks", p)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_every_live_node_path_exists(self):
        mod = self._module()
        assert mod._check_node_paths() >= len(mod.LOCKS)

    def test_a_rotted_entry_is_caught(self, monkeypatch):
        """Verified to go RED before the real entry was fixed; kept so it stays
        able to. An assertion that has never failed is not a lock."""
        mod = self._module()
        monkeypatch.setattr(
            mod, "NOT_COVERED",
            mod.NOT_COVERED + (("tests/test_gone_forever.py", "rotted"),))
        with pytest.raises(SystemExit, match="rotted"):
            mod._check_node_paths()

    def test_cannot_pass_by_being_empty(self, monkeypatch):
        """The 8th question: what does this criterion give on an empty sample?

        Without the count assertion an empty list passes vacuously -- the same
        shape as `f_bias == 1.000` when n = 1, and as a `d2h == 0` check whose
        baseline was already 0.
        """
        mod = self._module()
        monkeypatch.setattr(mod, "LOCKS", ())
        monkeypatch.setattr(mod, "NOT_COVERED", ())
        with pytest.raises(AssertionError):
            mod._check_node_paths()
