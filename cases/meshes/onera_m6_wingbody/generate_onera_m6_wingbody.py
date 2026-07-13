"""
Track M / M2: ONERA M6 wing + simplified axisymmetric fuselage, fused into a
wing-body half model, meshed WAKE-FREE for the level-set (Track B) solver path.

Run:  python cases/meshes/onera_m6_wingbody/generate_onera_m6_wingbody.py
      [--levels coarse medium fine]

Artifacts per level:
  <level>.msh                 volume mesh (GITIGNORED -- M1/M5 policy)
  <level>_stats.csv           counts, quality, junction + fuselage metrics
  <level>_wingbody.png        overview + junction zoom + cross sections

Why wake-free: roadmap M2's own note says the wake-fuselage junction is the
case a pre-embedded conforming sheet handles worst, and to schedule M2 with
Track B's embedded (level-set) wake, which removes the need to embed a sheet
at all. So nothing here fragments or embeds a wake surface; the wake lives
only as a corridor size field, and the solver builds it from the analytic TE
polyline (wingbody.te_polyline) via WakeLevelSet. The mesh has NO "wake" group.

Ladder: self-similar by construction -- every length scale is a fixed multiple
of h_wall and h_far is NOT clamped (the M1b defect, which put `coarse` off the
refinement ray and invalidated every past M6 three-point Richardson).

The fuselage is one splined body of revolution, not four fused primitives, so
its skin carries no seam creases: the FUSELAGE CREASE metric below is the
evidence (it must DECAY under refinement -- O(h) faceting of a smooth surface --
exactly the M5 tip-cap argument, applied to the body).
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from pyfp3d.mesh.metrics import compute_aspect_ratios, compute_min_dihedral_angles
from pyfp3d.mesh.reader import write_mesh
from pyfp3d.meshgen.fuselage import FuselageParams, radius_at
from pyfp3d.meshgen.wing3d import B_SEMI, TIP_CAP_RADIUS, x_te
from pyfp3d.meshgen.wingbody import (
    junction_z,
    onera_m6_wingbody_mesh,
    te_polyline,
)
from pyfp3d.post.surface import wall_crease_angles

OUT_DIR = Path(__file__).resolve().parent

# Same wall sizes as the M1 / M5 wing ladders, so the wing-body is a
# controlled A/B against the wing-alone families at equal h_wall.
LEVELS = {"coarse": 0.030, "medium": 0.015, "fine": 0.0075}
RICHARDSON_LADDER = ("coarse", "medium", "fine")

QUALITY_BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}

#: The fuselage skin is a SMOOTH surface of revolution, so its crease angle is
#: O(h * curvature) and must fall with refinement. A body-of-revolution seam
#: (what fusing 4 primitives would leave) would park at a fixed angle instead.
FUSELAGE_CREASE_MAX_DEG = {"coarse": 25.0, "medium": 15.0, "fine": 10.0}

FUSELAGE = FuselageParams()


def _level_params(h_wall: float) -> dict:
    """Every scale a fixed multiple of h_wall -- NO h_far clamp (M1b)."""
    return {
        "h_wall": h_wall,
        "h_far": 120.0 * h_wall,
        "h_wake": 3.0 * h_wall,
        "h_edge": 0.5 * h_wall,
        "h_tip": 0.25 * h_wall,
        "h_junction": 0.5 * h_wall,
    }


def _fuselage_crease(mesh) -> tuple:
    """Crease angle across the fuselage skin triangulation (p99, max)."""
    ang, _ = wall_crease_angles(mesh.nodes, mesh.elements,
                                mesh.boundary_faces["fuselage"])
    return float(np.percentile(ang, 99)), float(ang.max()), float(np.median(ang))


def _junction_stats(mesh, p: FuselageParams) -> dict:
    """How well the wing-fuselage junction and the wake polyline's inboard end
    are resolved."""
    wall = np.unique(mesh.boundary_faces["wall"])
    z_junc = junction_z(p)
    te_in = np.array([x_te(z_junc), 0.0, z_junc])
    d_te_in = float(np.linalg.norm(mesh.nodes[wall] - te_in, axis=1).min())
    corner = np.array([x_te(B_SEMI), 0.0, B_SEMI])
    d_tip_te = float(np.linalg.norm(mesh.nodes[wall] - corner, axis=1).min())
    # fuselage skin fidelity: distance of fuselage nodes from R(x)
    fus = np.unique(mesh.boundary_faces["fuselage"])
    rr = np.hypot(mesh.nodes[fus, 1], mesh.nodes[fus, 2])
    R = np.array([radius_at(p, float(x)) for x in mesh.nodes[fus, 0]])
    return {
        "junction_z": z_junc,
        "d_wake_inboard_end": d_te_in,
        "d_tip_te_corner": d_tip_te,
        "n_wall_nodes_near_junction": int((mesh.nodes[wall, 2] < z_junc + 0.05).sum()),
        "fuselage_radius_err_max": float(np.abs(rr - R).max()),
        "wall_z_max": float(mesh.nodes[wall, 2].max()),
    }


def generate_level(level: str, render: bool = True) -> dict:
    prm = _level_params(LEVELS[level])
    print(f"[{level}] h_wall={prm['h_wall']:.4f} ...", flush=True)
    mesh = onera_m6_wingbody_mesh(fuselage=FUSELAGE, tip_cap="round", **prm)

    dihedral = compute_min_dihedral_angles(mesh.nodes, mesh.elements)
    aspect = compute_aspect_ratios(mesh.nodes, mesh.elements)
    cre_p99, cre_max, cre_med = _fuselage_crease(mesh)

    stats = {
        "level": level,
        **prm,
        "n_nodes": len(mesh.nodes),
        "n_tets": len(mesh.elements),
        "min_dihedral_deg": float(dihedral.min()),
        "max_aspect_ratio": float(aspect.max()),
        "fuselage_crease_p99_deg": cre_p99,
        "fuselage_crease_max_deg": cre_max,
        "fuselage_crease_median_deg": cre_med,
        **{k: len(v) for k, v in
           (("n_tris_" + g, f) for g, f in mesh.boundary_faces.items())},
        **_junction_stats(mesh, FUSELAGE),
    }

    # --- gates ---------------------------------------------------------
    assert "wake" not in mesh.boundary_faces, "must be wake-free (LS path)"
    assert stats["min_dihedral_deg"] >= QUALITY_BOUNDS["min_dihedral_deg"], (
        f"min dihedral {stats['min_dihedral_deg']:.2f} deg below bound "
        f"{QUALITY_BOUNDS['min_dihedral_deg']}"
    )
    assert stats["max_aspect_ratio"] <= QUALITY_BOUNDS["max_aspect_ratio"], (
        f"max aspect {stats['max_aspect_ratio']:.1f} above bound "
        f"{QUALITY_BOUNDS['max_aspect_ratio']}"
    )
    assert cre_p99 <= FUSELAGE_CREASE_MAX_DEG[level], (
        f"fuselage crease p99 {cre_p99:.1f} deg above bound "
        f"{FUSELAGE_CREASE_MAX_DEG[level]} -- the skin has a seam, not a "
        "smooth surface of revolution"
    )

    write_mesh(mesh, str(OUT_DIR / f"{level}.msh"))
    with open(OUT_DIR / f"{level}_stats.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["key", "value"])
        for k, v in stats.items():
            w.writerow([k, v])
    if render:
        _write_png(OUT_DIR / f"{level}_wingbody.png", mesh, level)

    print(f"[{level}] {stats['n_tets']} tets | dihedral "
          f"{stats['min_dihedral_deg']:.2f} deg | aspect "
          f"{stats['max_aspect_ratio']:.1f} | fuselage crease p99 "
          f"{cre_p99:.1f} deg", flush=True)
    return stats


def _write_png(path: Path, mesh, level: str) -> None:
    """Overview + junction zoom + cross sections. The junction panel is the
    M2 visual check: the wing must emerge from the body with no crack."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    def edges(tris):
        t = np.asarray(tris, dtype=np.int64)
        return np.unique(np.sort(np.concatenate(
            [t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]), axis=1), axis=0)

    z_junc = junction_z(FUSELAGE)
    fig = plt.figure(figsize=(20, 6.6))

    ax = fig.add_subplot(1, 3, 1, projection="3d")
    for g, col in (("fuselage", "0.62"), ("wall", "tab:blue")):
        e = edges(mesh.boundary_faces[g])
        ax.add_collection3d(Line3DCollection(mesh.nodes[e], colors=col,
                                             linewidths=0.2))
    ax.set_xlim(FUSELAGE.x_nose_tip - 0.1, FUSELAGE.x_tail_tip + 0.1)
    ax.set_ylim(-0.6, 0.6)
    ax.set_zlim(0.0, 1.3)
    ax.view_init(elev=25, azim=-60)
    ax.set_title(f"{level}: wing (blue) + fuselage (grey)")
    for a, lab in ((ax.set_xlabel, "x"), (ax.set_ylabel, "y"),
                   (ax.set_zlabel, "z")):
        a(lab)

    ax = fig.add_subplot(1, 3, 2, projection="3d")
    for g, col in (("fuselage", "0.62"), ("wall", "tab:blue")):
        e = edges(mesh.boundary_faces[g])
        near = ((mesh.nodes[e, 2].max(axis=1) < 3 * z_junc)
                & (mesh.nodes[e, 0].max(axis=1) < x_te(z_junc) + 0.2)
                & (mesh.nodes[e, 0].min(axis=1) > -0.05))
        ax.add_collection3d(Line3DCollection(mesh.nodes[e[near]], colors=col,
                                             linewidths=0.35))
    ax.set_xlim(-0.05, x_te(z_junc) + 0.2)
    ax.set_ylim(-0.25, 0.25)
    ax.set_zlim(0.0, 3 * z_junc)
    ax.view_init(elev=28, azim=-70)
    ax.set_title(f"JUNCTION zoom (z < {3 * z_junc:.2f})\n"
                 "the wing must emerge with no crack")
    for a, lab in ((ax.set_xlabel, "x"), (ax.set_ylabel, "y"),
                   (ax.set_zlabel, "z")):
        a(lab)

    ax = fig.add_subplot(1, 3, 3)
    allw = np.vstack([mesh.boundary_faces["wall"],
                      mesh.boundary_faces["fuselage"]]).astype(np.int64)
    for x0, c in ((0.4, "tab:red"), (0.6, "tab:green"), (0.75, "tab:purple")):
        segs = _slice_x(mesh.nodes, allw, x0)
        if len(segs):
            ax.add_collection(LineCollection(segs, colors=c, linewidths=0.8,
                                             label=f"x = {x0:.2f}"))
    ax.axvline(z_junc, color="k", ls="--", lw=0.8,
               label=f"junction z = {z_junc:.2f}")
    ax.set_xlim(-0.02, 0.55)
    ax.set_ylim(-0.2, 0.2)
    ax.set_aspect("equal")
    ax.legend(fontsize=8)
    ax.set_xlabel("z")
    ax.set_ylabel("y")
    ax.set_title("cross sections: fuselage arc + wing")

    fig.suptitle(f"onera_m6_wingbody {level} (wake-free, round tip cap)")
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def _slice_x(nodes, tris, x0):
    """Intersection segments of the given triangles with x = x0, in (z, y)."""
    segs = []
    for t in tris:
        pts = nodes[t]
        d = pts[:, 0] - x0
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if d[a] * d[b] < 0.0:
                w = d[a] / (d[a] - d[b])
                p = pts[a] + w * (pts[b] - pts[a])
                segs.append(p[[2, 1]])
    if not segs:
        return np.empty((0, 2, 2))
    return np.asarray(segs).reshape(-1, 2, 2)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=["coarse"],
                    choices=list(LEVELS))
    ap.add_argument("--no-render", action="store_true")
    args = ap.parse_args()

    print(f"fuselage: r_f={FUSELAGE.r_f} length={FUSELAGE.length:.3f} "
          f"junction_z={junction_z(FUSELAGE):.3f}")
    print(f"TE polyline (level-set wake): "
          f"{np.array2string(te_polyline(FUSELAGE), precision=4)}")
    rows = [generate_level(lv, render=not args.no_render) for lv in args.levels]

    with open(OUT_DIR / "ladder_stats.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT_DIR}/ladder_stats.csv")


if __name__ == "__main__":
    main()
