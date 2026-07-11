"""
M4 deliverable: wake-FREE ONERA M6 half-wing tet mesh family (roadmap
Track M M4; design_track_b.md section 5.7 dual-mesh rule -- the 3D half).

Identical geometry, sizing policy and refinement ladder to the M1 family
(cases/meshes/onera_m6) except `embed_wake=False`: the chord-plane sheet is
built ONLY as the source of the wake-corridor Distance size field and is
neither fragmented into the fluid BRep nor embedded, so the tet mesh does
not conform to it and no "wake" tag exists. This is Track B's deliverable
form in 3D -- generic cuts through generic elements, with the swept TE and
the wing tip that the 2.5D M3 family cannot exercise (design_track_b.md D9
TE-polyline ruled surface; the spanwise clip that keeps the level set from
cutting beyond the tip).

Because the corridor size field is the SAME field the M1 family uses
(same h_wake, same Distance thresholds, sheet at the same chord-plane
location), element counts land close to M1's at equal h_wall -- which is
what makes the B4.5 A/B against the P5/P8 baselines a controlled
comparison. A wide alpha-sweep corridor (the 2.5D M3 wedge) is deliberately
NOT used here: in 3D the wedge volume scales with span and the element cost
is prohibitive (recorded in the roadmap M4 row) -- 3D alpha sweeps re-aim
the level set within the near-nominal band the corridor already resolves.

Tags: wall, farfield, symmetry. NO wake tag.

No .msh committed (gitignored like the M1 family; coarse+medium regenerate
in ~30 s). The per-level stats CSVs are the committed evidence.

Usage:
    python generate_onera_m6_wakefree.py --level coarse
    python generate_onera_m6_wakefree.py --all          # includes fine
(default: coarse + medium)
"""

import argparse
import csv
import time
from pathlib import Path

import numpy as np

from pyfp3d.mesh.reader import write_mesh, mesh_stats
from pyfp3d.mesh.metrics import compute_aspect_ratios, compute_min_dihedral_angles
from pyfp3d.meshgen.wing3d import B_SEMI, MAC, onera_m6_wing_mesh, x_te

R_FAR = 15.0 * MAC
QUALITY_BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}


def _level_params(h_wall: float) -> dict:
    """Same policy as the M1 family (one parameter per level)."""
    return dict(
        h_wall=h_wall,
        h_wake=3.0 * h_wall,
        h_edge=0.5 * h_wall,
        h_far=min(2.5, 120.0 * h_wall),
    )


LEVELS = {
    "coarse": _level_params(0.030),
    "medium": _level_params(0.015),
    "fine": _level_params(0.0075),
}


def _corridor_stats(mesh, h_wake: float) -> dict:
    """Edge lengths of elements whose centroid lies in the wake corridor
    (behind the TE, within the sheet's span, |y| < 0.15 m) -- the M4 gate's
    sizing evidence against the M1 family's sheet-adjacent h_wake."""
    cen = mesh.nodes[mesh.elements].mean(axis=1)
    z = cen[:, 2]
    in_span = (z > 0.0) & (z < B_SEMI)
    behind = cen[:, 0] > np.array([x_te(zi) for zi in z])
    near_plane = np.abs(cen[:, 1]) < 0.15
    near_wing = cen[:, 0] < x_te(B_SEMI) + 2.0
    sel = in_span & behind & near_plane & near_wing
    el = mesh.elements[sel].astype(np.int64)
    if len(el) == 0:
        return {"corridor_n_elements": 0}
    p = mesh.nodes
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    e = np.concatenate([
        np.linalg.norm(p[el[:, i]] - p[el[:, j]], axis=1) for i, j in pairs
    ])
    return {
        "corridor_n_elements": int(len(el)),
        "corridor_median_edge": float(np.median(e)),
        "corridor_p90_edge": float(np.percentile(e, 90)),
        "h_wake_target": h_wake,
    }


def generate_level(out_dir: Path, level: str) -> Path:
    p = LEVELS[level]
    t0 = time.perf_counter()
    mesh = onera_m6_wing_mesh(
        h_wall=p["h_wall"], h_far=p["h_far"], h_wake=p["h_wake"],
        h_edge=p["h_edge"], r_far=R_FAR,
        name=f"onera_m6_wakefree_{level}", embed_wake=False,
    )
    gen_seconds = time.perf_counter() - t0
    assert "wake" not in mesh.boundary_faces, "M4 mesh must carry no wake tag"

    out_path = out_dir / f"{level}.msh"
    write_mesh(mesh, out_path)

    stats = mesh_stats(mesh)
    dihedral = compute_min_dihedral_angles(mesh.nodes, mesh.elements)
    aspect = compute_aspect_ratios(mesh.nodes, mesh.elements)
    stats["min_dihedral_deg"] = float(dihedral.min())
    stats["p01_dihedral_deg"] = float(np.percentile(dihedral, 1))
    stats["max_aspect_ratio"] = float(aspect.max())
    stats["p99_aspect_ratio"] = float(np.percentile(aspect, 99))
    for tag, faces in mesh.boundary_faces.items():
        stats[f"n_{tag}_faces"] = len(faces)
    stats["h_wall"] = p["h_wall"]
    stats["gen_seconds"] = round(gen_seconds, 1)
    stats.update(_corridor_stats(mesh, p["h_wake"]))

    ok = (stats["min_dihedral_deg"] >= QUALITY_BOUNDS["min_dihedral_deg"]
          and stats["max_aspect_ratio"] <= QUALITY_BOUNDS["max_aspect_ratio"])
    stats["quality_within_bounds"] = ok

    with open(out_dir / f"{level}_stats.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for k, v in stats.items():
            writer.writerow([k, v])

    if not ok:
        raise AssertionError(
            f"{level}: element quality outside bounds "
            f"(min dihedral {stats['min_dihedral_deg']:.2f} deg, "
            f"max aspect {stats['max_aspect_ratio']:.1f}; "
            f"bounds {QUALITY_BOUNDS})"
        )
    return out_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", action="append", choices=sorted(LEVELS))
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    levels = sorted(LEVELS) if args.all else (args.level or ["coarse", "medium"])

    out_dir = Path(__file__).parent
    for level in levels:
        out_path = generate_level(out_dir, level)
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
