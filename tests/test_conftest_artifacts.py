"""
Regression tests for the gate-artifact fixtures in conftest.py.

The `artifacts_dir` fixture used to hand out a `tempfile.TemporaryDirectory`,
so every PNG/CSV a gate test produced was deleted at teardown and the repo's
`artifacts/` directory stayed permanently empty -- violating the
CLAUDE.md/roadmap workflow rule that every visual gate leaves inspectable
headless artifacts behind. These tests lock in the fix: artifacts land in
the persistent (gitignored) `<repo>/artifacts/` by default, and
PYFP3D_ARTIFACTS_DIR redirects them (e.g. to a CI upload directory).
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
