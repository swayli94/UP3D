import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = "tests/fixtures/single_tetrahedron.json"


def invoke(path: str) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env["PYTHONNOUSERSITE"] = "1"
    return subprocess.run(
        [sys.executable, "-m", "pyfp3d", path],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def test_documented_cli_reads_real_file_and_writes_json() -> None:
    completed = invoke(FIXTURE)

    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    result = json.loads(completed.stdout)
    assert set(result) == {"velocities", "free_stream_preservation_error"}
    assert result["velocities"] == [[0.75, -1.25, 2.5]]
    assert result["free_stream_preservation_error"] == 0.0


def test_missing_file_fails_without_success_json(tmp_path: Path) -> None:
    completed = invoke(str(tmp_path / "missing.json"))

    assert completed.returncode != 0
    assert completed.stdout == ""
    assert "pyfp3d:" in completed.stderr
    assert "missing.json" in completed.stderr


def test_invalid_document_fails_without_success_json(tmp_path: Path) -> None:
    path = tmp_path / "invalid.json"
    path.write_text("{}", encoding="utf-8")

    completed = invoke(str(path))

    assert completed.returncode != 0
    assert completed.stdout == ""
    assert "missing required field" in completed.stderr


def test_cli_preserves_large_integer_free_stream(tmp_path: Path) -> None:
    source = json.loads((ROOT / FIXTURE).read_text(encoding="utf-8"))
    large_integer = 9007199254740993
    source["free_stream_velocity"] = [large_integer, -1.25, 2.5]
    path = tmp_path / "large-integer.json"
    path.write_text(json.dumps(source), encoding="utf-8")

    completed = invoke(str(path))

    assert completed.returncode == 0, completed.stderr
    result = json.loads(completed.stdout)
    assert result["velocities"] == [[large_integer, -1.25, 2.5]]
    assert isinstance(result["velocities"][0][0], int)
    assert result["free_stream_preservation_error"] == 0.0
