"""
Regression tests for the gate-artifact fixtures in conftest.py.

The `artifacts_dir` fixture used to hand out a `tempfile.TemporaryDirectory`,
so every PNG/CSV a gate test produced was deleted at teardown and the repo's
`artifacts/` directory stayed permanently empty. These tests lock in the fix:
output lands in the persistent (gitignored) `<repo>/artifacts/` by default, and
PYFP3D_ARTIFACTS_DIR redirects it (e.g. to a CI upload directory).

*** What this directory IS and IS NOT (corrected 2026-08-24) ***
It is SCRATCH -- the mechanism that forces headless output so a visual gate can
never decay into a GUI-only check. It is NOT evidence: it is gitignored, so
nothing written here reaches HEAD, and discipline 3 says a claim without a
committed artifact is not evidence. The workflow rule used to name this
directory as the evidence form for visual gates, which made the two rules
contradict each other and left 11 documents citing `artifacts/G1.3/` and
`artifacts/G2.{1..5}/` -- paths that exist in no fresh clone. C/D-class figures
now go to a TRACKED `results/` directory instead.
"""

import os
from pathlib import Path

import pytest

REPO_ARTIFACTS = Path(__file__).parent.parent / "artifacts"


@pytest.mark.skipif(
    "PYFP3D_ARTIFACTS_DIR" in os.environ,
    reason="artifacts explicitly redirected by PYFP3D_ARTIFACTS_DIR",
)
def test_gate_artifacts_land_in_persistent_repo_dir(gate_artifacts_dir):
    assert REPO_ARTIFACTS.resolve() in gate_artifacts_dir.resolve().parents, (
        f"gate artifacts dir {gate_artifacts_dir} is not under {REPO_ARTIFACTS} -- "
        "artifacts written there would not survive the test run"
    )
    probe = gate_artifacts_dir / "persistence_probe.txt"
    probe.write_text("gate artifacts must outlive the test run\n")
    assert probe.exists()


def test_artifacts_dir_env_override(tmp_path, monkeypatch, request):
    redirected = tmp_path / "ci_upload"
    monkeypatch.setenv("PYFP3D_ARTIFACTS_DIR", str(redirected))
    artifacts_dir = request.getfixturevalue("artifacts_dir")
    assert artifacts_dir == redirected
    assert artifacts_dir.is_dir()


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
