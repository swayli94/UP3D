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
refinement ray and invalidated every past M6 three-point Richardson). h_body
and h_body_tip are multiples of h_wall for the same reason -- see _level_params.

The fuselage is one splined body of revolution, not four fused primitives, so
its skin carries no seam creases: the FUSELAGE CREASE metric below is the
evidence (it must DECAY under refinement -- O(h) faceting of a smooth surface --
exactly the M5 tip-cap argument, applied to the body).

BODY RE-SPEC 2026-07-16 (user-directed), on the body delivered 2026-07-13:
  * it was 2.13 long with its nose tip 0.30 (= 0.37 C_ROOT) ahead of the wing
    root LE, so the nose's own displacement flow sat on the wing. The body is
    now 5 root chords long with the wing root chord CENTERED on it (2 root
    chords of body fore and aft) and a slender 2-diameter ellipsoid nose --
    all three as RULES in fuselage.py, not station numbers.
  * the body skin no longer runs at h_wall everywhere: it relaxes to
    h_body = 2 h_wall away from the wing, where the flow is a smooth
    displacement field, and is driven by the local body RADIUS at the two tips
    (wingbody._fuselage_field). Net: the body nearly doubled in length for
    ~the same fuselage triangle count.
  * FAR FIELD (same-day follow-up): r_far was 15 MAC -- the wing-ALONE
    convention -- which left the 5-chord body spanning ~42% of the far-field
    diameter. Now R_FAR = 25 MAC = 3.51 body lengths clear of each tip, with
    h_far and every fixed refinement distance scaled with it so the 2.78×
    domain costs ~nothing and no near-body gradient moves (wingbody.R_FAR /
    H_FAR_IN_H_WALL / H_FAR_REF).
All of it is measured in <level>_stats.csv (body_*, fus_edge_*, r_far_* and
farfield_clearance_* keys) and drawn in the BODY SIZING panel of
<level>_wingbody.png -- the re-spec is a committed artifact, not prose.
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
from pyfp3d.meshgen.wing3d import B_SEMI, C_ROOT, MAC, TIP_CAP_RADIUS, x_te
from pyfp3d.meshgen.wingbody import (
    H_FAR_IN_H_WALL,
    R_FAR,
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
    """Every scale a fixed multiple of h_wall -- NO h_far clamp (M1b).

    h_body is a MULTIPLE of h_wall, never a clamp, for the same reason: the
    gate metric below is the fuselage crease angle, and it is only evidence of
    "faceting of a smooth surface, not a seam" if it DECAYS with the level. A
    clamped body size would park it at a constant and forge a seam.
    """
    return {
        "h_wall": h_wall,
        "h_far": H_FAR_IN_H_WALL * h_wall,
        "h_wake": 3.0 * h_wall,
        "h_edge": 0.5 * h_wall,
        "h_tip": 0.25 * h_wall,
        "h_junction": 0.5 * h_wall,
        "h_body": 2.0 * h_wall,
        "h_body_tip": 0.25 * h_wall,
    }


def _fuselage_crease(mesh) -> tuple:
    """Crease angle across the fuselage skin triangulation (p99, max)."""
    ang, _ = wall_crease_angles(mesh.nodes, mesh.elements,
                                mesh.boundary_faces["fuselage"])
    return float(np.percentile(ang, 99)), float(ang.max()), float(np.median(ang))


def _fuselage_sizing(mesh, p: FuselageParams) -> dict:
    """Median fuselage-skin edge length, bucketed to measure the 2026-07-16
    graded body sizing. Without this the grading would be a prose claim.

    Buckets are all taken on the CONSTANT-RADIUS section, so they compare
    like with like: the nose and afterbody carry their own radius-driven
    refinement (wingbody._fuselage_field), which would otherwise read as
    "the far body is finer than the near body" and measure nothing.
    """
    tris = np.asarray(mesh.boundary_faces["fuselage"], dtype=np.int64)
    pts = mesh.nodes[tris]
    edge = np.linalg.norm(pts[:, [1, 2, 0]] - pts, axis=2).mean(axis=1)
    xc = pts[:, :, 0].mean(axis=1)

    def med(sel):
        return float(np.median(edge[sel])) if sel.any() else float("nan")

    # Chordwise distance from the wing root chord (0 under the wing).
    dx = np.maximum(0.0, np.maximum(xc - C_ROOT, -xc))
    cyl = (xc >= p.x_nose_end) & (xc <= p.x_body_end)
    return {
        "fus_edge_med_under_wing": med(cyl & (dx == 0.0)),
        "fus_edge_med_cyl_far": med(cyl & (dx > 0.5)),
        "fus_edge_med_nose_tip": med(xc < p.x_nose_tip + 0.1),
        "fus_edge_med_tail_tip": med(xc > p.x_tail_tip - 0.1),
        "fus_edge_med_all": med(np.ones_like(xc, dtype=bool)),
    }


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


def _body_geometry(p: FuselageParams) -> dict:
    """The 2026-07-16 proportion + far-field re-spec, as numbers in the
    record. The far-field keys are what make "the boundary is far enough from
    the BODY, not just from the wing" checkable rather than asserted."""
    xc = 0.5 * (p.x_nose_tip + max(x_te(B_SEMI), p.x_tail_tip))
    return {
        "body_length": p.length,
        "body_length_in_root_chords": p.length / C_ROOT,
        "x_nose_tip": p.x_nose_tip,
        "x_tail_tip": p.x_tail_tip,
        "nose_tip_to_root_le": 0.0 - p.x_nose_tip,     # x_le(0) = 0
        "root_te_to_tail_tip": p.x_tail_tip - C_ROOT,
        "l_nose_in_diameters": p.l_nose / (2.0 * p.r_f),
        "nose_tip_curv_radius": p.r_f ** 2 / p.l_nose,
        "r_far": R_FAR,
        "r_far_in_mac": R_FAR / MAC,
        "farfield_clearance_nose": p.x_nose_tip - (xc - R_FAR),
        "farfield_clearance_tail": (xc + R_FAR) - p.x_tail_tip,
        "farfield_clearance_in_body_lengths":
            ((xc + R_FAR) - p.x_tail_tip) / p.length,
        "body_span_frac_of_farfield_diameter": p.length / (2.0 * R_FAR),
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
        **_fuselage_sizing(mesh, FUSELAGE),
        **_body_geometry(FUSELAGE),
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
    # Graded body sizing (2026-07-16 re-spec), both directions: the far body
    # must actually be coarser, and the body under the wing must NOT be.
    ratio = stats["fus_edge_med_cyl_far"] / stats["fus_edge_med_under_wing"]
    assert ratio >= 1.5, (
        f"far body only {ratio:.2f}x coarser than the body under the wing -- "
        "the graded body sizing is not taking effect"
    )
    assert stats["fus_edge_med_under_wing"] <= 1.2 * prm["h_wall"], (
        f"body under the wing at {stats['fus_edge_med_under_wing']:.4f} > "
        f"h_wall {prm['h_wall']:.4f} -- h_body has leaked into the junction "
        "region, where the wing's own field must win"
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
          f"{cre_p99:.1f} deg | body edge under-wing "
          f"{stats['fus_edge_med_under_wing']:.4f} vs far-cylinder "
          f"{stats['fus_edge_med_cyl_far']:.4f}", flush=True)
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

    def true_aspect(ax):
        """Draw to scale. The default cube would stretch the 5-chord body's
        y/z by ~4x and hide exactly what these panels are for."""
        ax.set_box_aspect([hi - lo for lo, hi in
                           (ax.get_xlim(), ax.get_ylim(), ax.get_zlim())])

    z_junc = junction_z(FUSELAGE)
    fig = plt.figure(figsize=(26, 6.6))

    ax = fig.add_subplot(1, 4, 1, projection="3d")
    for g, col in (("fuselage", "0.62"), ("wall", "tab:blue")):
        e = edges(mesh.boundary_faces[g])
        ax.add_collection3d(Line3DCollection(mesh.nodes[e], colors=col,
                                             linewidths=0.2))
    ax.set_xlim(FUSELAGE.x_nose_tip - 0.1, FUSELAGE.x_tail_tip + 0.1)
    ax.set_ylim(-0.6, 0.6)
    ax.set_zlim(0.0, 1.3)
    ax.view_init(elev=25, azim=-60)
    ax.set_title(f"{level}: wing (blue) + fuselage (grey)\n"
                 f"body {FUSELAGE.length / C_ROOT:.1f} root chords, "
                 "wing centered; body skin coarsens away from the wing")
    for a, lab in ((ax.set_xlabel, "x"), (ax.set_ylabel, "y"),
                   (ax.set_zlabel, "z")):
        a(lab)
    true_aspect(ax)

    ax = fig.add_subplot(1, 4, 2, projection="3d")
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
    true_aspect(ax)

    ax = fig.add_subplot(1, 4, 3)
    _plot_body_sizing(ax, mesh, level)

    ax = fig.add_subplot(1, 4, 4)
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
    fig.subplots_adjust(wspace=0.32)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def _plot_body_sizing(ax, mesh, level: str) -> None:
    """Fuselage skin element size along the body, against the body radius.

    The evidence panel for the 2026-07-16 re-spec: the skin must sit at h_wall
    under the wing, relax to h_body away from it, and dive back down at the
    two tips -- where it is not the wing that sets the size but R(x), since
    the faceting angle of a surface of revolution is h / R.
    """
    p = FUSELAGE
    prm = _level_params(LEVELS[level])
    tris = np.asarray(mesh.boundary_faces["fuselage"], dtype=np.int64)
    pts = mesh.nodes[tris]
    edge = np.linalg.norm(pts[:, [1, 2, 0]] - pts, axis=2).mean(axis=1)
    xc = pts[:, :, 0].mean(axis=1)

    bins = np.linspace(p.x_nose_tip, p.x_tail_tip, 60)
    idx = np.clip(np.digitize(xc, bins) - 1, 0, len(bins) - 2)
    xm = 0.5 * (bins[:-1] + bins[1:])
    med = np.array([np.median(edge[idx == i]) if (idx == i).any() else np.nan
                    for i in range(len(xm))])

    ax.plot(xm, med, "o-", ms=3, lw=1.2, color="tab:blue",
            label="median skin edge")
    for h, c, lab in ((prm["h_wall"], "tab:green", "h_wall"),
                      (prm["h_body"], "tab:red", "h_body = 2 h_wall"),
                      (prm["h_body_tip"], "tab:purple", "h_body_tip")):
        ax.axhline(h, color=c, ls="--", lw=0.9, label=lab)
    ax.axvspan(0.0, C_ROOT, color="tab:green", alpha=0.12,
               label="wing root chord")
    ax.set_xlim(p.x_nose_tip, p.x_tail_tip)
    ax.set_ylim(0.0, 1.7 * prm["h_body"])
    ax.set_xlabel("x")
    ax.set_ylabel("element size")
    ax.legend(fontsize=7, loc="upper center", ncol=2)
    ax.set_title("BODY SIZING: fine under the wing, h_body away,\n"
                 "R-driven at the tips (facet angle = h / R)")

    axr = ax.twinx()
    rr = [radius_at(p, float(x)) for x in xm]
    axr.plot(xm, rr, color="0.55", lw=1.0)
    axr.fill_between(xm, rr, color="0.55", alpha=0.15)
    axr.set_ylim(0.0, 6.0 * p.r_f)
    axr.set_yticks([])
    axr.text(p.x_center, 0.02, "body radius R(x)", color="0.4", fontsize=7,
             ha="center")


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
