"""Track V / GV3.3: standalone fuselage-alone body-of-revolution mesh family.

Run:  python cases/meshes/fuselage_bor/generate_fuselage_bor.py
      [--levels coarse medium]

Artifacts per level:
  <level>.msh                 volume mesh (GITIGNORED -- M1/M5 policy)
  <level>_stats.csv           counts, quality, radius-fidelity metrics
  <level>_fuselage_bor.png    meridian radius fidelity + sizing panels

Why this family exists (roadmap track_v.md V3, GV3.3): the loose FP+IBL
coupling needs ONE genuinely-3-D closed-surface exercise -- a body of
revolution at alpha = 0, non-lifting, wake-free, with the wall tagged as a
single "wall" group so the coupling's surface machinery applies unchanged.
No wing, no junction, no symmetry plane: the sphere is the FULL 2pi domain
(the wing-body family's z >= 0 half-model is NOT reused -- axisymmetry must
be exercised, not imposed).

Geometry: FuselageParams() defaults (meshgen/fuselage.py), the SAME body as
the ONERA M6 wing-body (5 root chords, 2-diameter ellipsoid nose, cone
afterbody, tail sphere cap), built by `add_fuselage_solid` (full-2pi
revolve; the split variant exists only for the conforming wake sheet, which
this family does not have). Far field: sphere R_FAR = 25 MAC centered on
the body mid-length, h_far = 200 h_body -- the wingbody.py far-field
policy, so h_far/r_far and every size gradient match that family.

Sizing: the body is the SUBJECT here (no wing to prioritize), so the whole
skin runs at h_body = the level size, except the two tips, where the size
follows the local body radius (wingbody._fuselage_field's tip-ball law,
replicated: a constant h over the collapsing R(x) facets the afterbody
into a polygonal needle).

Ladder: self-similar -- h_far and every fixed distance scale with h_body
(no h_far clamp, the M1b defect). coarse/medium only; GV3.3 is a smoke,
not a Richardson study.
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from pyfp3d.mesh.metrics import (
    compute_aspect_ratios,
    compute_min_dihedral_angles,
    compute_tet_volumes,
)
from pyfp3d.mesh.reader import write_mesh
from pyfp3d.meshgen.fuselage import (
    FuselageParams,
    add_fuselage_solid,
    radius_at,
)
from pyfp3d.meshgen.wing3d import MAC, _collect_3d
from pyfp3d.meshgen.wingbody import H_FAR_IN_H_WALL, R_FAR
from pyfp3d.post.surface import wall_crease_angles

OUT_DIR = Path(__file__).resolve().parent

#: Skin element size per level (h_body == h_wall: the body is the subject).
LEVELS = {"coarse": 0.030, "medium": 0.015}

QUALITY_BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}

#: Smooth-revolution-skin crease bound (decays with refinement; a seam
#: would park -- the M5 tip-cap argument, as in the wing-body family).
FUSELAGE_CREASE_MAX_DEG = {"coarse": 25.0, "medium": 15.0}

FUSELAGE = FuselageParams()


def fuselage_bor_mesh(h_body: float, fuselage: FuselageParams = FUSELAGE,
                      r_far: float = R_FAR, name: str = "fuselage_bor",
                      verbose: bool = False):
    """Full-2pi body of revolution in a full sphere; groups wall/farfield.

    h_body: skin size away from the tips. h_body_tip = 0.25 h_body at the
    nose/tail tips (the wingbody tip-ball law); h_far = 200 h_body.
    Sphere center: the body mid-length station.
    """
    import gmsh

    p = fuselage
    h_far = H_FAR_IN_H_WALL * h_body
    h_body_tip = 0.25 * h_body
    xc = p.x_center

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.model.add(name)
        occ = gmsh.model.occ

        body = add_fuselage_solid(occ, p)
        ball = occ.addSphere(xc, 0.0, 0.0, r_far)
        fluid, _ = occ.cut([(3, ball)], body)
        occ.synchronize()

        vols = gmsh.model.getEntities(3)
        assert len(vols) == 1, f"expected one fluid volume, got {vols}"

        # --- classify boundary surfaces: sphere vs body skin -------------
        groups = {"wall": [], "farfield": []}
        for dim, tag in gmsh.model.getEntities(2):
            bb = gmsh.model.getBoundingBox(dim, tag)
            extent = max(bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2])
            if extent > 1.2 * r_far:
                groups["farfield"].append(tag)
            else:
                groups["wall"].append(tag)
        for g, tags in groups.items():
            assert tags, f"surface classification found no '{g}' faces"

        # --- size field: radius-driven tips + distance ramp (Max-composed,
        #     the wingbody._fuselage_field law; see its docstring) ---------
        field = gmsh.model.mesh.field
        f_dist = field.add("Distance")
        field.setNumbers(f_dist, "SurfacesList", groups["wall"])
        field.setNumber(f_dist, "Sampling", 200)
        ramp = field.add("Threshold")
        field.setNumber(ramp, "InField", f_dist)
        field.setNumber(ramp, "SizeMin", h_body_tip)
        field.setNumber(ramp, "SizeMax", h_far)
        field.setNumber(ramp, "DistMin", 0.05)
        field.setNumber(ramp, "DistMax", 0.55 * r_far)

        def _tip_ball(x_tip: float, ramp_len: float) -> int:
            f = field.add("Ball")
            field.setNumber(f, "XCenter", x_tip)
            field.setNumber(f, "YCenter", 0.0)
            field.setNumber(f, "ZCenter", 0.0)
            field.setNumber(f, "Radius", 0.2 * p.r_f)
            field.setNumber(f, "Thickness", ramp_len)
            field.setNumber(f, "VIn", h_body_tip)
            field.setNumber(f, "VOut", h_body)
            return f

        tips = field.add("Min")
        field.setNumbers(tips, "FieldsList",
                         [_tip_ball(p.x_nose_tip, 0.5 * p.l_nose),
                          _tip_ball(p.x_tail_tip, p.l_tail + p.r_tail)])
        f_max = field.add("Max")
        field.setNumbers(f_max, "FieldsList", [ramp, tips])
        field.setAsBackgroundMesh(f_max)

        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.Optimize", 1)
        # wake-free closed body: Netgen on (the wingbody wake-free policy)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

        gmsh.model.mesh.generate(3)
        mesh = _collect_3d(groups, name=name)
        _asserts(mesh, p, r_far=r_far, xc=xc)
        return mesh
    finally:
        gmsh.finalize()


def _asserts(mesh, p: FuselageParams, r_far: float, xc: float) -> None:
    """Generation-time invariants (wingbody._wingbody_asserts subset)."""
    vols = compute_tet_volumes(mesh.nodes, mesh.elements)
    assert vols.min() > 0.0, "non-positive tet volume"
    assert set(mesh.boundary_faces) == {"wall", "farfield"}, \
        f"groups must be exactly wall+farfield: {set(mesh.boundary_faces)}"

    far_nodes = np.unique(mesh.boundary_faces["farfield"])
    r = np.linalg.norm(mesh.nodes[far_nodes] - np.array([xc, 0.0, 0.0]),
                       axis=1)
    assert np.all(np.abs(r - r_far) < 1e-3 * r_far), \
        "far-field nodes not on the sphere"

    # Wall skin genuinely lies on the revolution surface (GV3.3's geometry
    # anchor: the axisymmetry assertions downstream are meaningless unless
    # the wall IS the analytic body).
    wall_nodes = np.unique(mesh.boundary_faces["wall"])
    fx = mesh.nodes[wall_nodes, 0]
    frr = np.hypot(mesh.nodes[wall_nodes, 1], mesh.nodes[wall_nodes, 2])
    fR = np.array([radius_at(p, float(x)) for x in fx])
    assert np.abs(frr - fR).max() < 0.01 * p.r_f + 3e-3, \
        "wall group contains faces off the revolution surface"


def _stats(mesh, p: FuselageParams, level: str, h_body: float) -> dict:
    dihedral = compute_min_dihedral_angles(mesh.nodes, mesh.elements)
    aspect = compute_aspect_ratios(mesh.nodes, mesh.elements)
    ang, _ = wall_crease_angles(mesh.nodes, mesh.elements,
                                mesh.boundary_faces["wall"])
    wall = np.unique(mesh.boundary_faces["wall"])
    rr = np.hypot(mesh.nodes[wall, 1], mesh.nodes[wall, 2])
    R = np.array([radius_at(p, float(x)) for x in mesh.nodes[wall, 0]])

    tris = np.asarray(mesh.boundary_faces["wall"], dtype=np.int64)
    pts = mesh.nodes[tris]
    edge = np.linalg.norm(pts[:, [1, 2, 0]] - pts, axis=2).mean(axis=1)
    xtri = pts[:, :, 0].mean(axis=1)

    def med(sel):
        return float(np.median(edge[sel])) if sel.any() else float("nan")

    cyl = (xtri >= p.x_nose_end) & (xtri <= p.x_body_end)
    return {
        "level": level,
        "h_body": h_body,
        "h_far": H_FAR_IN_H_WALL * h_body,
        "r_far": R_FAR,
        "r_far_in_mac": R_FAR / MAC,
        "n_nodes": len(mesh.nodes),
        "n_tets": len(mesh.elements),
        "n_tris_wall": len(mesh.boundary_faces["wall"]),
        "n_tris_farfield": len(mesh.boundary_faces["farfield"]),
        "min_dihedral_deg": float(dihedral.min()),
        "max_aspect_ratio": float(aspect.max()),
        "wall_crease_p99_deg": float(np.percentile(ang, 99)),
        "wall_crease_max_deg": float(ang.max()),
        "wall_radius_err_max": float(np.abs(rr - R).max()),
        "wall_edge_med_cylinder": med(cyl),
        "wall_edge_med_nose_tip": med(xtri < p.x_nose_tip + 0.1),
        "wall_edge_med_tail_tip": med(xtri > p.x_tail_tip - 0.1),
        "body_length": p.length,
        "x_nose_tip": p.x_nose_tip,
        "x_tail_tip": p.x_tail_tip,
        "farfield_clearance_in_body_lengths":
            ((p.x_center + R_FAR) - p.x_tail_tip) / p.length,
    }


def _panel(mesh, p: FuselageParams, level: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    wall = np.unique(mesh.boundary_faces["wall"])
    x = mesh.nodes[wall, 0]
    rr = np.hypot(mesh.nodes[wall, 1], mesh.nodes[wall, 2])
    xs = np.linspace(p.x_nose_tip, p.x_tail_tip, 400)
    Rs = np.array([radius_at(p, float(v)) for v in xs])

    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    ax = axes[0]
    ax.plot(x, rr, ".", ms=1.5, label="wall nodes")
    ax.plot(xs, Rs, "r-", lw=1.4, label="analytic R(x)")
    ax.set_xlabel("x")
    ax.set_ylabel("sqrt(y^2+z^2)")
    ax.set_title(f"GV3.3 BoR wall radius fidelity ({level})")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_aspect("equal", adjustable="datalim")

    tris = np.asarray(mesh.boundary_faces["wall"], dtype=np.int64)
    pts = mesh.nodes[tris]
    edge = np.linalg.norm(pts[:, [1, 2, 0]] - pts, axis=2).mean(axis=1)
    xtri = pts[:, :, 0].mean(axis=1)
    bins = np.linspace(p.x_nose_tip, p.x_tail_tip, 60)
    idx = np.digitize(xtri, bins)
    meds = [np.median(edge[idx == i]) if np.any(idx == i) else np.nan
            for i in range(1, len(bins))]
    ax = axes[1]
    ax.plot(0.5 * (bins[1:] + bins[:-1]), meds, "-", lw=1.4)
    ax.set_xlabel("x")
    ax.set_ylabel("median wall edge length")
    ax.set_title("sizing: radius-driven tips, h_body on the cylinder")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    path = OUT_DIR / f"{level}_fuselage_bor.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}", flush=True)


def generate_level(level: str) -> dict:
    h_body = LEVELS[level]
    print(f"[{level}] h_body={h_body:.4f} r_far={R_FAR:.2f} ...", flush=True)
    mesh = fuselage_bor_mesh(h_body)
    stats = _stats(mesh, FUSELAGE, level, h_body)

    assert stats["min_dihedral_deg"] >= QUALITY_BOUNDS["min_dihedral_deg"], \
        (f"min dihedral {stats['min_dihedral_deg']:.2f} deg below bound "
         f"{QUALITY_BOUNDS['min_dihedral_deg']}")
    assert stats["max_aspect_ratio"] <= QUALITY_BOUNDS["max_aspect_ratio"], \
        (f"max aspect {stats['max_aspect_ratio']:.1f} above bound "
         f"{QUALITY_BOUNDS['max_aspect_ratio']}")
    assert stats["wall_crease_p99_deg"] <= FUSELAGE_CREASE_MAX_DEG[level], \
        (f"wall crease p99 {stats['wall_crease_p99_deg']:.1f} deg above "
         f"bound {FUSELAGE_CREASE_MAX_DEG[level]} -- the skin has a seam, "
         "not a smooth surface of revolution")

    write_mesh(mesh, str(OUT_DIR / f"{level}.msh"))
    with open(OUT_DIR / f"{level}_stats.csv", "w") as f:
        f.write("metric,value\n")
        for k, v in stats.items():
            f.write(f"{k},{v}\n")
    print(f"  wrote {OUT_DIR / (level + '_stats.csv')}", flush=True)
    _panel(mesh, FUSELAGE, level)
    print(f"  n_nodes={stats['n_nodes']} n_tets={stats['n_tets']} "
          f"wall_tris={stats['n_tris_wall']} "
          f"dihedral={stats['min_dihedral_deg']:.2f} "
          f"crease_p99={stats['wall_crease_p99_deg']:.2f}", flush=True)
    return stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=list(LEVELS),
                    choices=list(LEVELS))
    args = ap.parse_args()
    for level in args.levels:
        generate_level(level)


if __name__ == "__main__":
    main()
