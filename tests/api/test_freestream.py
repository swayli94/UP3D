from fractions import Fraction
import json
import math
from pathlib import Path

import pytest

from pyfp3d import run_freestream


FIXTURE = Path(__file__).parents[1] / "fixtures" / "single_tetrahedron.json"


def document() -> dict:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def write_document(tmp_path: Path, value: object) -> Path:
    path = tmp_path / "input.json"
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def test_single_tetrahedron_output_contract() -> None:
    source = document()

    result = run_freestream(FIXTURE)

    assert set(result) == {"velocities", "free_stream_preservation_error"}
    assert result["velocities"] == [source["free_stream_velocity"]]
    assert result["free_stream_preservation_error"] == 0.0
    assert all(math.isfinite(component) for component in result["velocities"][0])
    recomputed = max(
        abs(component - expected)
        for velocity in result["velocities"]
        for component, expected in zip(velocity, source["free_stream_velocity"])
    )
    assert result["free_stream_preservation_error"] == recomputed


def test_multiple_cells(tmp_path: Path) -> None:
    source = document()
    source["tetrahedra"] = [[0, 1, 2, 3], [3, 2, 1, 0]]

    result = run_freestream(write_document(tmp_path, source))

    assert result["velocities"] == [source["free_stream_velocity"]] * 2
    assert result["free_stream_preservation_error"] == 0.0


def test_malformed_json_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "broken.json"
    path.write_text("{", encoding="utf-8")

    with pytest.raises(json.JSONDecodeError):
        run_freestream(path)


@pytest.mark.parametrize("root", [[], "not an object", 1])
def test_non_object_documents_are_rejected(tmp_path: Path, root: object) -> None:
    with pytest.raises(ValueError, match="JSON object"):
        run_freestream(write_document(tmp_path, root))


def test_missing_field_is_diagnosed(tmp_path: Path) -> None:
    source = document()
    del source["vertices"]

    with pytest.raises(ValueError, match="vertices"):
        run_freestream(write_document(tmp_path, source))


@pytest.mark.parametrize("version", [True, 1.0, 2])
def test_format_version_must_be_integer_one(tmp_path: Path, version: object) -> None:
    source = document()
    source["format_version"] = version

    with pytest.raises(ValueError, match="format_version"):
        run_freestream(write_document(tmp_path, source))


@pytest.mark.parametrize("field", ["vertices", "free_stream_velocity"])
def test_nonfinite_numbers_are_rejected(tmp_path: Path, field: str) -> None:
    source = document()
    if field == "vertices":
        source[field][0][0] = float("nan")
    else:
        source[field][1] = float("inf")

    with pytest.raises(ValueError, match=field):
        run_freestream(write_document(tmp_path, source))


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("vertices", [[0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        ("vertices", [[0.0, 0.0, 0.0]] * 3),
        ("free_stream_velocity", [1.0, 2.0]),
        ("tetrahedra", []),
        ("tetrahedra", [[0, 1, 2]]),
    ],
)
def test_invalid_shapes_and_empty_cells_are_rejected(
    tmp_path: Path, field: str, value: object,
) -> None:
    source = document()
    source[field] = value

    with pytest.raises(ValueError, match=field):
        run_freestream(write_document(tmp_path, source))


@pytest.mark.parametrize("index", [True, 1.0, -1, 4])
def test_invalid_indices_are_rejected(tmp_path: Path, index: object) -> None:
    source = document()
    source["tetrahedra"] = [[0, 1, 2, index]]

    with pytest.raises(ValueError, match=r"tetrahedra\[0\]"):
        run_freestream(write_document(tmp_path, source))


def test_repeated_indices_are_rejected(tmp_path: Path) -> None:
    source = document()
    source["tetrahedra"] = [[0, 1, 2, 2]]

    with pytest.raises(ValueError, match="distinct"):
        run_freestream(write_document(tmp_path, source))


def test_boolean_coordinates_are_rejected(tmp_path: Path) -> None:
    source = document()
    source["vertices"][0][0] = True

    with pytest.raises(ValueError, match="vertices"):
        run_freestream(write_document(tmp_path, source))


def test_degenerate_cell_is_rejected(tmp_path: Path) -> None:
    source = document()
    source["vertices"][3] = [1.0, 1.0, 0.0]

    with pytest.raises(ValueError, match="degenerate"):
        run_freestream(write_document(tmp_path, source))


def test_missing_file_raises_os_error(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        run_freestream(tmp_path / "missing.json")


def test_large_integer_free_stream_preserves_original_numeric_semantics(
    tmp_path: Path,
) -> None:
    source = document()
    large_integer = 9007199254740993
    source["free_stream_velocity"] = [large_integer, -1.25, 2.5]

    result = run_freestream(write_document(tmp_path, source))

    assert result["velocities"] == [[large_integer, -1.25, 2.5]]
    assert isinstance(result["velocities"][0][0], int)
    recomputed = max(
        abs(Fraction(component) - Fraction(expected))
        for velocity in result["velocities"]
        for component, expected in zip(velocity, source["free_stream_velocity"])
    )
    assert result["free_stream_preservation_error"] == float(recomputed) == 0.0


def test_exact_coplanarity_rejects_large_integer_cell(tmp_path: Path) -> None:
    source = document()
    source["vertices"] = [
        [54, 42, -42],
        [92837491, 77263052, 13974157],
        [-40060229, 70427411, 71680763],
        [52777208, 147690421, 85654962],
    ]

    with pytest.raises(ValueError, match="degenerate"):
        run_freestream(write_document(tmp_path, source))


def test_tiny_nondegenerate_cell_is_accepted(tmp_path: Path) -> None:
    source = document()
    source["vertices"] = [
        [0.0, 0.0, 0.0],
        [1e-200, 0.0, 0.0],
        [0.0, 1e-200, 0.0],
        [0.0, 0.0, 1e-200],
    ]

    result = run_freestream(write_document(tmp_path, source))

    assert result["velocities"] == [source["free_stream_velocity"]]
    assert result["free_stream_preservation_error"] == 0.0


def test_boolean_free_stream_is_rejected(tmp_path: Path) -> None:
    source = document()
    source["free_stream_velocity"][0] = True

    with pytest.raises(ValueError, match="free_stream_velocity"):
        run_freestream(write_document(tmp_path, source))


@pytest.mark.parametrize("field", ["vertices", "tetrahedra", "free_stream_velocity"])
def test_non_list_containers_are_rejected(tmp_path: Path, field: str) -> None:
    source = document()
    source[field] = {"not": "a list"}

    with pytest.raises(ValueError, match=field):
        run_freestream(write_document(tmp_path, source))


@pytest.mark.parametrize("field", ["format_version", "vertices", "free_stream_velocity"])
def test_non_numeric_strings_are_rejected(tmp_path: Path, field: str) -> None:
    source = document()
    if field == "format_version":
        source[field] = "1"
    elif field == "vertices":
        source[field][0][0] = "0"
    else:
        source[field][0] = "fast"

    with pytest.raises(ValueError, match=field):
        run_freestream(write_document(tmp_path, source))
