"""
M0 deliverable: single-layer extruded NACA0012 quasi-2D mesh family
(docs/roadmap.md Track M, gate link G2.5).

The 2D airfoil domain (closed sharp TE, circular far field of radius R_FAR
centered at mid-chord) is meshed with vanilla Gmsh, with the wake line
TE -> farfield *embedded* in the surface (gmsh.model.mesh.embed) so triangle
edges conform to it. Then exactly ONE cell layer is extruded in z and every
prism is subdivided into 3 tets with the globally consistent
min-global-index diagonal rule (pyfp3d/meshgen/extrude.py); the wake edges
become the tagged interior face sheet "wake" (nodes NOT duplicated -- that
stays in the solver preprocessor, mesh/wake_cut.py, per hard rule 8).

Tags: wall (airfoil), farfield (outer circle), symmetry (both z-planes),
wake (interior sheet TE -> farfield).

One parameter per level (h = h_wall); everything else scales with it.
Single-layer targets ~15k / 60k / 240k tets (~5k / 20k / 80k 2D triangles).

Usage:
    python generate_naca0012.py --level coarse --level medium
    python generate_naca0012.py --all
"""

import argparse
import csv
from pathlib import Path

import numpy as np

from pyfp3d.mesh.reader import write_mesh, mesh_stats
from pyfp3d.meshgen.extrude import extrude_single_layer
from pyfp3d.meshgen.planar import naca0012_wake_2d

R_FAR = 15.0  # far-field radius in chords, centered at (0.5, 0)


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
    "fine": _level_params(0.005),
}


def generate_level(out_dir: Path, level: str, inspect: bool = True) -> Path:
    p = LEVELS[level]
    points2d, triangles, edge_groups, interior_groups = naca0012_wake_2d(
        r_far=R_FAR,
        h_wall=p["h_wall"], h_far=p["h_far"], h_wake=p["h_wake"],
        dist_min=p["dist_min"], dist_max=p["dist_max"],
        wake_dist_max=p["wake_dist_max"], n_half=p["n_half"],
    )
    mesh = extrude_single_layer(
        points2d, triangles, edge_groups,
        interior_edge_groups=interior_groups,
        dz=p["dz"], z0=0.0, name=f"naca0012_2.5d_{level}",
    )

    out_path = out_dir / f"{level}.msh"
    write_mesh(mesh, out_path)

    stats = mesh_stats(mesh)
    stats["n_2d_triangles"] = len(triangles)
    stats["n_wake_faces"] = len(mesh.boundary_faces["wake"])
    stats["dz"] = p["dz"]
    stats["h_wall"] = p["h_wall"]
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
    axes[0].set_title(f"naca0012_2.5d {level}: 2D layer + tags")
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
