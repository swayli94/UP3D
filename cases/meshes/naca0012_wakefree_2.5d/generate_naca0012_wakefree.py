"""
M3 deliverable: wake-free ("O-mesh style") single-layer extruded NACA0012
quasi-2D mesh family (docs/roadmap.md Track M M3; design_track_b.md §5.7).

Identical pipeline to the M0 family (cases/meshes/naca0012_2.5d) except the
wake line is NOT embedded: the mesh topology knows nothing about the wake --
this is Track B's level-set deliverable form (generic cuts through generic
elements). Because there is then no conforming sheet to attract refinement,
a fan of size-field-only lines from the TE spanning the alpha-sweep envelope
CORRIDOR_ALPHA_DEG keeps element size ~h_wake in the wedge the level-set
wake will sweep through (design_track_b.md D8). The fan lines are never
embedded, so triangle edges do NOT conform to them.

Tags: wall (airfoil), farfield (outer circle), symmetry (both z-planes).
No "wake" tag exists -- ingest must work without it (M3 gate).

The stats CSV records corridor_median_edge / corridor_p90_edge (edge lengths
of triangles whose centroid lies in the corridor wedge behind the TE) as the
M3 gate's sizing evidence against M0's sheet-adjacent h_wake = 3*h_wall.

Usage:
    python generate_naca0012_wakefree.py --level coarse --level medium
    python generate_naca0012_wakefree.py --all
"""

import argparse
import csv
from pathlib import Path

import numpy as np

from pyfp3d.mesh.reader import write_mesh, mesh_stats
from pyfp3d.meshgen.extrude import extrude_single_layer
from pyfp3d.meshgen.planar import naca0012_wake_2d

R_FAR = 15.0  # far-field radius in chords, centered at (0.5, 0) -- same as M0
CORRIDOR_ALPHA_DEG = (-6.0, 6.0)  # alpha-sweep envelope covered by the corridor
CORRIDOR_N_LINES = 5


def _level_params(h_wall: float) -> dict:
    """Everything derived from the single wall-size parameter h (M0 policy)."""
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


def _corridor_edge_stats(points2d: np.ndarray, triangles: np.ndarray,
                         wake_dist_max: float) -> dict:
    """Edge-length stats of triangles whose centroid lies in the corridor
    wedge behind the TE (the M3 sizing-evidence clause)."""
    te = np.array([1.0, 0.0])
    cen = points2d[triangles].mean(axis=1)
    rel = cen - te
    r = np.linalg.norm(rel, axis=1)
    ang = np.degrees(np.arctan2(rel[:, 1], rel[:, 0]))
    in_wedge = (
        (r > 1e-9) & (r <= wake_dist_max)
        & (ang >= CORRIDOR_ALPHA_DEG[0]) & (ang <= CORRIDOR_ALPHA_DEG[1])
    )
    tri = triangles[in_wedge]
    if len(tri) == 0:
        return {"corridor_n_triangles": 0}
    p = points2d
    e = np.concatenate([
        np.linalg.norm(p[tri[:, 0]] - p[tri[:, 1]], axis=1),
        np.linalg.norm(p[tri[:, 1]] - p[tri[:, 2]], axis=1),
        np.linalg.norm(p[tri[:, 2]] - p[tri[:, 0]], axis=1),
    ])
    return {
        "corridor_n_triangles": int(len(tri)),
        "corridor_median_edge": float(np.median(e)),
        "corridor_p90_edge": float(np.percentile(e, 90)),
    }


def generate_level(out_dir: Path, level: str, inspect: bool = True) -> Path:
    p = LEVELS[level]
    points2d, triangles, edge_groups, interior_groups = naca0012_wake_2d(
        r_far=R_FAR,
        h_wall=p["h_wall"], h_far=p["h_far"], h_wake=p["h_wake"],
        dist_min=p["dist_min"], dist_max=p["dist_max"],
        wake_dist_max=p["wake_dist_max"], n_half=p["n_half"],
        embed_wake=False,
        corridor_alpha_deg=CORRIDOR_ALPHA_DEG,
        corridor_n_lines=CORRIDOR_N_LINES,
    )
    assert interior_groups == {}, "wake-free build must not return a wake group"
    mesh = extrude_single_layer(
        points2d, triangles, edge_groups,
        interior_edge_groups=None,
        dz=p["dz"], z0=0.0, name=f"naca0012_wakefree_2.5d_{level}",
    )
    assert "wake" not in mesh.boundary_faces

    out_path = out_dir / f"{level}.msh"
    write_mesh(mesh, out_path)

    stats = mesh_stats(mesh)
    stats["n_2d_triangles"] = len(triangles)
    stats["dz"] = p["dz"]
    stats["h_wall"] = p["h_wall"]
    stats["h_wake_target"] = p["h_wake"]
    stats["corridor_alpha_deg"] = f"{CORRIDOR_ALPHA_DEG[0]}..{CORRIDOR_ALPHA_DEG[1]}"
    stats.update(_corridor_edge_stats(points2d, triangles, p["wake_dist_max"]))
    with open(out_dir / f"{level}_stats.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for k, v in stats.items():
            writer.writerow([k, v])

    if inspect:
        _write_inspection_png(out_dir / f"{level}_layer.png", points2d,
                              triangles, edge_groups, level)
    return out_path


def _write_inspection_png(path, points2d, triangles, edge_groups, level):
    """Headless artifact (roadmap Sec 0.1): 2D layer + tags; the corridor
    wedge is drawn dashed for reference (it is NOT a mesh entity)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    styles = {"wall": ("tab:red", 1.4), "farfield": ("tab:blue", 1.2)}
    for ax, (x0, x1, y0, y1) in zip(
        axes,
        [(-R_FAR + 0.5, R_FAR + 0.5, -R_FAR, R_FAR), (-0.3, 2.0, -0.6, 0.6)],
    ):
        ax.triplot(points2d[:, 0], points2d[:, 1], triangles,
                   linewidth=0.2, color="0.65")
        for tag, edges in edge_groups.items():
            color, lw = styles.get(tag, ("k", 1.0))
            seg = points2d[np.asarray(edges)]
            for s in seg:
                ax.plot(s[:, 0], s[:, 1], color=color, linewidth=lw)
        for ang in np.radians(CORRIDOR_ALPHA_DEG):
            L = R_FAR - 0.5
            ax.plot([1.0, 1.0 + L * np.cos(ang)], [0.0, L * np.sin(ang)],
                    "g--", linewidth=0.9)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
    axes[0].set_title(f"naca0012_wakefree_2.5d {level}: no wake in topology")
    axes[1].set_title("airfoil zoom (corridor wedge dashed, size-field only)")
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
