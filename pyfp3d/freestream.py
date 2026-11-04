from __future__ import annotations

from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any


Number = int | float


def _finite_number(value: Any, field: str) -> Number:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field} must contain only JSON numbers")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{field} must contain only finite numbers")
    return value


def _vector(value: Any, field: str) -> list[Number]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{field} must be an array of three finite numbers")
    return [_finite_number(component, field) for component in value]


def _scaled_volume_six(vertices: list[list[Number]], cell: list[int]) -> Fraction:
    v0, v1, v2, v3 = [
        [Fraction(component) for component in vertices[index]] for index in cell
    ]
    a = [v1[axis] - v0[axis] for axis in range(3)]
    b = [v2[axis] - v0[axis] for axis in range(3)]
    c = [v3[axis] - v0[axis] for axis in range(3)]
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _validate(document: Any) -> tuple[list[list[int]], list[Number]]:
    if not isinstance(document, dict):
        raise ValueError("input must be a JSON object")

    required = {"format_version", "vertices", "tetrahedra", "free_stream_velocity"}
    missing = sorted(required - document.keys())
    if missing:
        raise ValueError(f"missing required field(s): {', '.join(missing)}")

    version = document["format_version"]
    if isinstance(version, bool) or not isinstance(version, int) or version != 1:
        raise ValueError("format_version must be the integer 1")

    raw_vertices = document["vertices"]
    if not isinstance(raw_vertices, list) or len(raw_vertices) < 4:
        raise ValueError("vertices must contain at least four three-dimensional points")
    vertices = [_vector(point, "vertices") for point in raw_vertices]

    raw_cells = document["tetrahedra"]
    if not isinstance(raw_cells, list) or not raw_cells:
        raise ValueError("tetrahedra must contain at least one cell")
    cells: list[list[int]] = []
    for cell_index, cell in enumerate(raw_cells):
        field = f"tetrahedra[{cell_index}]"
        if not isinstance(cell, list) or len(cell) != 4:
            raise ValueError(f"{field} must contain exactly four vertex indices")
        if any(isinstance(index, bool) or not isinstance(index, int) for index in cell):
            raise ValueError(f"{field} indices must be integers, not booleans")
        if any(index < 0 or index >= len(vertices) for index in cell):
            raise ValueError(f"{field} contains an out-of-range vertex index")
        if len(set(cell)) != 4:
            raise ValueError(f"{field} indices must be distinct")
        cells.append(cell)

    free_stream = _vector(document["free_stream_velocity"], "free_stream_velocity")

    for cell_index, cell in enumerate(cells):
        if _scaled_volume_six(vertices, cell) == 0:
            raise ValueError(f"tetrahedra[{cell_index}] is degenerate")

    return cells, free_stream


def run_freestream(path: str | Path) -> dict[str, object]:
    input_path = Path(path)
    with input_path.open("r", encoding="utf-8") as stream:
        document = json.load(stream)

    cells, free_stream = _validate(document)
    velocities = [list(free_stream) for _ in cells]
    error = max(
        (
            abs(Fraction(component) - Fraction(expected))
            for velocity in velocities
            for component, expected in zip(velocity, free_stream)
        ),
        default=Fraction(0),
    )
    return {
        "velocities": velocities,
        "free_stream_preservation_error": float(error),
    }
