"""
GV5.2 deliverable: single-layer extruded RAE2822 quasi-2D mesh family
(docs/roadmap/track_v.md V5 GV5.2; pre-registration
bench/studies/v5_2_rae2822/PRE_REGISTRATION.md).

Geometry: the Cook, McDonald & Firmin (AGARD-AR-138) Table 6.1 measured
ordinates, committed verbatim as rae2822.dat beside this script
(source: https://www.grc.nasa.gov/www/wind/valid/raetaf/geom.txt ,
fetched 2026-07-24, sha256
32a2803a211223c2899ff036b32223367fa5566c4bfb999a2e556813e1be43c5).
The point-set path
(pyfp3d/meshgen/planar.py::pointset_airfoil_coordinates) PCHIP-resamples
each side onto the cosine-clustered grid; the TE closes sharp at (1, 0).

Everything else mirrors the NACA0012 M0 family (embedded wake line
TE -> farfield, circular far field R_FAR centered at mid-chord, the same
Distance+Threshold size recipe, single-layer z extrusion, prism -> 3 tets
min-global-index rule, tags wall / farfield / symmetry / wake). Levels:
coarse + medium ONLY (fine registered, not built).

Usage:
    python generate_rae2822.py --level coarse --level medium
    python generate_rae2822.py --all
"""

import argparse
import csv
import os
import sys
from pathlib import Path

# resolve pyfp3d from THIS worktree (the site-packages editable install may
# point at a sibling worktree)
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..")))

import numpy as np

from pyfp3d.mesh.reader import write_mesh, mesh_stats
from pyfp3d.meshgen.extrude import extrude_single_layer
from pyfp3d.meshgen.planar import (
    airfoil_wake_2d,
    load_airfoil_ordinates,
    pointset_airfoil_coordinates,
    te_wedge_angle_deg,
)

R_FAR = 15.0  # far-field radius in chords, centered at (0.5, 0)
DAT = Path(__file__).parent / "rae2822.dat"


def _level_params(h_wall: float) -> dict:
    """Everything derived from the single wall-size parameter h."""
    return dict(
        h_wall=h_wall,
        h_wake=3.0 * h_wall,
        h_far=min(3.0, 150.0 * h_wall),
        dist_min=0.1,
        dist_max=6.0,
        wake_dist_max=1.5,
        dz=2.0 * h_wall,
        n_half=max(80, int(round(2.0 / h_wall))),
    )


LEVELS = {
    "coarse": _level_params(0.020),
    "medium": _level_params(0.010),
}


def generate_level(out_dir: Path, level: str, inspect: bool = True) -> Path:
    p = LEVELS[level]
    x, z_lo, z_up = load_airfoil_ordinates(DAT)
    # Cook Table 6.1 convention: the lower ordinate is tabulated
    # POSITIVE-DOWN (distance below the chord line) -- negate for the
    # physical z (verified: thickness 12.1% @ 37.9%, camber 1.3% @ 75.7%).
    z_lo = -z_lo
    coords = pointset_airfoil_coordinates(x, z_lo, z_up, n_half=p["n_half"])
    points2d, triangles, edge_groups, interior_groups = airfoil_wake_2d(
        coords, model_name=f"rae2822_2d_{level}",
        r_far=R_FAR,
        h_wall=p["h_wall"], h_far=p["h_far"], h_wake=p["h_wake"],
        dist_min=p["dist_min"], dist_max=p["dist_max"],
        wake_dist_max=p["wake_dist_max"],
    )
    mesh = extrude_single_layer(
        points2d, triangles, edge_groups,
        interior_edge_groups=interior_groups,
        dz=p["dz"], z0=0.0, name=f"rae2822_2.5d_{level}",
    )

    out_path = out_dir / f"{level}.msh"
    write_mesh(mesh, out_path)

    stats = mesh_stats(mesh)
    stats["n_2d_triangles"] = len(triangles)
    stats["n_wake_faces"] = len(mesh.boundary_faces["wake"])
    stats["dz"] = p["dz"]
    stats["h_wall"] = p["h_wall"]
    stats["te_wedge_deg_ordinates"] = te_wedge_angle_deg(x, z_lo, z_up)
    with open(out_dir / f"{level}_stats.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for k, v in stats.items():
            writer.writerow([k, v])

    if inspect:
        _write_inspection_png(out_dir / f"{level}_layer.png", points2d,
                              triangles, edge_groups, interior_groups, level)
    return out_path


def _write_inspection_png(path, points2d, triangles, edge_groups,
                          interior_groups, level):
    """Headless artifact (roadmap Sec 0.1): 2D layer with tagged boundaries and
    the embedded wake line, full domain + airfoil zoom."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    styles = {"wall": ("tab:red", 1.4), "farfield": ("tab:blue", 1.2),
              "wake": ("tab:green", 1.2)}
    all_groups = dict(edge_groups)
    all_groups.update(interior_groups)
    for ax, (x0, x1, y0, y1) in zip(
        axes,
        [(-R_FAR + 0.5, R_FAR + 0.5, -R_FAR, R_FAR), (-0.3, 2.0, -0.6, 0.6)],
    ):
        ax.triplot(points2d[:, 0], points2d[:, 1], triangles,
                   linewidth=0.2, color="0.65")
        for tag, edges in all_groups.items():
            color, lw = styles.get(tag, ("k", 1.0))
            seg = points2d[np.asarray(edges)]
            for s in seg:
                ax.plot(s[:, 0], s[:, 1], color=color, linewidth=lw)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
    axes[0].set_title(f"rae2822_2.5d {level}: 2D layer + tags")
    axes[1].set_title("airfoil zoom (wall=red, wake=green, farfield=blue)")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", action="append", choices=sorted(LEVELS))
    parser.add_argument("--all", action="store_true", help="generate all levels")
    args = parser.parse_args()
    levels = sorted(LEVELS) if args.all else (args.level or ["coarse", "medium"])

    out_dir = Path(__file__).parent
    for level in levels:
        out_path = generate_level(out_dir, level)
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
