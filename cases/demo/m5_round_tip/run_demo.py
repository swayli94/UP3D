"""Track M / M5 -- the ONERA M6 tip cap is now ROUND, and the mesh proves it.

WHY. P13/G13.3 ran the three-region box study that finally localized what was
still blocking 3D grid convergence, and it was not in the wake at all:

    region                            peak-Mach exponent p    verdict
    wake free edge (G13.2's fix)              +0.045          bounded
    wing interior (control)                   -0.014          converged
    tip-cap edge (WALL)                       +0.321          DIVERGES

`meshgen/wing3d.py` closed the wing with a FLAT tip cap -- a documented
deliberate simplification, standard for FP validation meshes. A flat cap meets
the upper and lower surfaces at a sharp convex edge, and in potential flow a
sharp convex edge is a singularity: refinement RESOLVES it instead of removing
it, the lift sequence never enters the asymptotic range, and no Richardson may
be extrapolated. It is not a P11 problem -- isoparametric ELEMENTS cannot
regularize a genuinely sharp geometric EDGE. The geometry was wrong.

THE FIX (`tip_cap="round"`): close the wing with the half body of revolution
swept by the tip section about its own chord line, so the cap radius at chord
station x IS the section half thickness t(x). Then

  * the cap meets the wing tangentially (dy/dz = 0 on both sides at z = B_SEMI)
    -- there is no edge; and
  * t -> 0 at the LE and the TE, so the cap degenerates to a POINT at each,
    which means the TE line, the wake sheet, the tip TE corner, the Kutta
    stations and B_SEMI are all exactly what they were. The A/B against the M1
    family is therefore controlled: the tip WALL moved, and nothing else did.

WHAT THIS DEMO ASSERTS. That the cap is round in the only geometry the solver
ever sees -- the TRIANGULATION. The discriminator is the crease angle across
the tip-section seam (the locus that IS the sharp edge in M1): facets
approximating a smooth surface crease by O(h * curvature) and HALVE when h
halves, while a real edge creases by its own turning angle and does not move.
Both families are measured by the same code here, because the M1 number is what
the M5 number has to be better than.

This demo is mesh-only (no solve, ~1 min). The FLOW consequence -- the tip-cap
box exponent, and whether the lift sequence finally becomes asymptotic -- is
P13/G13.3's, and lives in
cases/demo/p13_tip_edge_singularity/run_g133_roundtip.py.

  python cases/demo/m5_round_tip/run_demo.py
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from cases.demo._common import (  # noqa: E402
    CRITICAL, CheckList, INK_2, S1_BLUE, S2_AQUA,
    apply_style, finish, plt, write_csv,
)
from pyfp3d.mesh.metrics import (  # noqa: E402
    compute_aspect_ratios, compute_min_dihedral_angles,
)
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.meshgen.wing3d import (  # noqa: E402
    B_SEMI, C_TIP, ONERA_D_UPPER, TIP_CAP_RADIUS, x_le, x_te,
)
from pyfp3d.post.surface import wall_crease_angles  # noqa: E402

OUT = Path(__file__).parent / "results"
MESHES = REPO / "cases" / "meshes"

#: Levels are keyed by h_wall (0.030 / 0.015), NOT by file name -- because the
#: two families spell their self-similar coarse end differently. M1's shipped
#: `coarse` still carries the M1b h_far CLAMP (h_far 2.5 instead of 3.6), so it
#: is NOT on its own refinement ray and is the wrong thing to compare against;
#: `coarse_ss` is the level M1b re-cut without the clamp, and is the one that
#: sits on the RICHARDSON_LADDER. M5 needed no such repair -- its h_far was
#: never clamped -- so its coarse IS its ladder coarse.
LEVELS = ("coarse", "medium", "fine")
FAMILIES = {
    "flat (M1)": ("onera_m6", {"coarse": "coarse_ss", "medium": "medium",
                               "fine": "fine"}),
    "round (M5)": ("onera_m6_roundtip", {"coarse": "coarse", "medium": "medium",
                                         "fine": "fine"}),
}


def mesh_path(fam, level):
    d, names = FAMILIES[fam]
    return MESHES / d / f"{names[level]}.msh"

SEAM_XI = (0.05, 0.95)
QUALITY_BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}

#: Chordwise station of the tip section's max thickness -- where the cap has its
#: largest radius, and so the clearest place to draw it.
XI_TMAX = float(ONERA_D_UPPER[np.argmax(ONERA_D_UPPER[:, 1]), 0])


def seam_crease(mesh):
    """Crease angle across the tip-section seam at z = B_SEMI, away from the LE
    and TE (both sharp BY DESIGN in either family -- the TE carries the Kutta
    condition, so it is not the cap's business)."""
    ang, mid = wall_crease_angles(mesh.nodes, mesh.elements,
                                  mesh.boundary_faces["wall"])
    xi = (mid[:, 0] - x_le(B_SEMI)) / C_TIP
    on_seam = ((np.abs(mid[:, 2] - B_SEMI) < 1e-9)
               & (xi > SEAM_XI[0]) & (xi < SEAM_XI[1]))
    return ang[on_seam]


def measure(mesh):
    seam = seam_crease(mesh)
    wall = np.unique(mesh.boundary_faces["wall"])
    return dict(
        ntet=len(mesh.elements),
        seam_max=float(seam.max()),
        seam_p99=float(np.percentile(seam, 99)),
        n_seam=len(seam),
        z_wall_max=float(mesh.nodes[wall, 2].max()),
        min_dihedral=float(compute_min_dihedral_angles(
            mesh.nodes, mesh.elements).min()),
        max_aspect=float(compute_aspect_ratios(mesh.nodes, mesh.elements).max()),
    )


def tip_section_outline(mesh, half_width=0.02):
    """Wall nodes within a thin chordwise slab at the tip section's max-thickness
    station, projected to the (z, y) plane -- i.e. the shape the cap actually
    has, read off the triangulation rather than off the CAD."""
    x0 = x_le(B_SEMI) + XI_TMAX * C_TIP
    wall = np.unique(mesh.boundary_faces["wall"])
    p = mesh.nodes[wall]
    near = (np.abs(p[:, 0] - x0) < half_width * C_TIP) & (p[:, 2] > B_SEMI - 0.075)
    return p[near]


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    cl = CheckList("Track M / M5: the ONERA M6 tip cap is round")

    missing = [str(mesh_path(fam, lv)) for fam in FAMILIES for lv in LEVELS
               if not mesh_path(fam, lv).exists()]
    if missing:
        print("missing meshes: " + ", ".join(missing) + "\n"
              "  python cases/meshes/onera_m6/generate_onera_m6.py "
              "--level coarse_ss --level medium\n"
              "  python cases/meshes/onera_m6_roundtip/"
              "generate_onera_m6_roundtip.py")
        return 1

    M = {(fam, lv): read_mesh(mesh_path(fam, lv))
         for fam in FAMILIES for lv in LEVELS}
    D = {k: measure(m) for k, m in M.items()}

    def seq(fam, key):
        return [D[(fam, lv)][key] for lv in LEVELS]

    def exponent(fam, key):
        """crease ~ h^q, fitted against the ACTUAL element counts (so a mesh
        family knocked off its refinement ray would show up here, not hide)."""
        x = np.log(np.array(seq(fam, "ntet"), float) ** (1 / 3))
        return float(np.polyfit(x, np.log(seq(fam, key)), 1)[0])

    flat, rnd = seq("flat (M1)", "seam_max"), seq("round (M5)", "seam_max")
    flat99, rnd99 = seq("flat (M1)", "seam_p99"), seq("round (M5)", "seam_p99")
    q_flat, q_rnd = exponent("flat (M1)", "seam_p99"), exponent("round (M5)", "seam_p99")

    # ---- the defect, restated on this metric so the fix has a baseline -------
    cl.add("M5", "★ the FLAT cap is a real edge (it does not refine away)",
           "seam crease " + " -> ".join(f"{v:.1f}" for v in flat)
           + f" deg over the ladder; exponent q = {q_flat:+.2f}",
           "> 80 deg at EVERY level, q ~ 0 (no decay whatsoever)",
           min(flat) > 80.0 and abs(q_flat) < 0.05,
           note="a sharp convex edge creases by its own turning angle; halving "
                "h RESOLVES it and removes nothing. This is the object "
                "P13/G13.3 measured DIVERGING in the flow (p = +0.321). Note "
                "the MEDIAN seam crease is 0.00 deg -- the rest of the seam is "
                "already smooth, so this metric is reading the edge and nothing "
                "else")

    # ---- the fix -------------------------------------------------------------
    cl.add("M5", "★ the ROUND cap creases like O(h) => no edge in the limit",
           f"seam crease p99 " + " -> ".join(f"{v:.1f}" for v in rnd99)
           + f" deg; exponent q = {q_rnd:+.2f} (max: "
           + " -> ".join(f"{v:.1f}" for v in rnd) + ")",
           "q <= -0.7 -- the crease DECAYS with h, i.e. there is no edge in "
           "the limit", q_rnd <= -0.7,
           note="facets approximating a SMOOTH surface crease by "
                "O(h * curvature). Fitted on the p99, not the max: the max is a "
                "single worst facet at the thin end of the seam window (where "
                "the cap radius is smallest and the facet size is set by h_edge, "
                "not h_tip), so it decays more slowly (46.8 -> 25.0 -> 18.1) "
                "while the p99 halves cleanly. Both DECAY -- which is the "
                "claim; the flat cap's does not move at all")
    cl.add("M5", "the cap is resolved to its apex",
           "; ".join(f"{lv}: wall z_max {D[('round (M5)', lv)]['z_wall_max']:.6f}"
                     for lv in LEVELS)
           + f" (apex {B_SEMI + TIP_CAP_RADIUS:.6f})",
           "reaches the apex to within 1% of the cap radius",
           all(D[("round (M5)", lv)]["z_wall_max"]
               > B_SEMI + 0.99 * TIP_CAP_RADIUS for lv in LEVELS),
           note=f"cap radius = the tip section's max half thickness, "
                f"{TIP_CAP_RADIUS * 1000:.1f} mm = {100 * TIP_CAP_RADIUS / B_SEMI:.1f}% "
                "of the semi-span; h_tip = 0.25 h_wall exists so the coarse cap "
                "is not one element wide (which would discretize it back flat)")

    # ---- nothing else moved --------------------------------------------------
    for lv in LEVELS:
        mesh = M[("round (M5)", lv)]
        wall = np.unique(mesh.boundary_faces["wall"])
        corner = np.array([x_te(B_SEMI), 0.0, B_SEMI])
        d_corner = float(np.linalg.norm(mesh.nodes[wall] - corner, axis=1).min())
        wake = np.unique(mesh.boundary_faces["wake"])
        cl.add("M5", f"NOTHING ELSE MOVED ({lv}): TE corner + wake sheet",
               f"tip TE corner offset {d_corner:.2e} m; wake |y|max "
               f"{np.abs(mesh.nodes[wake, 1]).max():.2e}, z_max "
               f"{mesh.nodes[wake, 2].max():.6f}",
               "corner exact, wake in the chord plane, ending at B_SEMI",
               d_corner < 1e-9
               and np.abs(mesh.nodes[wake, 1]).max() < 1e-9
               and mesh.nodes[wake, 2].max() <= B_SEMI + 1e-9,
               note="the cap radius vanishes at the TE, so it degenerates to a "
                    "point exactly at the corner the wake attaches to. That is "
                    "what makes the M1 A/B controlled")

    _, wc = cut_wake(M[("round (M5)", "coarse")])
    cl.add("M5", "cut_wake keeps the M1 semantics",
           f"{len(wc.station_z)} TE stations, outermost z = "
           f"{wc.station_z.max():.4f} < B_SEMI = {B_SEMI:.4f}",
           "tip TE corner stays a FREE edge node (Gamma(tip) = 0 discretely)",
           len(wc.station_z) > 50 and wc.station_z.max() < B_SEMI)

    cl.add("M5", "element quality inside the M1 bounds",
           "; ".join(f"{lv}: min dihedral "
                     f"{D[('round (M5)', lv)]['min_dihedral']:.2f} deg, max aspect "
                     f"{D[('round (M5)', lv)]['max_aspect']:.1f}" for lv in LEVELS),
           f"min dihedral >= {QUALITY_BOUNDS['min_dihedral_deg']}, "
           f"max aspect <= {QUALITY_BOUNDS['max_aspect_ratio']}",
           all(D[("round (M5)", lv)]["min_dihedral"]
               >= QUALITY_BOUNDS["min_dihedral_deg"]
               and D[("round (M5)", lv)]["max_aspect"]
               <= QUALITY_BOUNDS["max_aspect_ratio"] for lv in LEVELS),
           note="a small-radius feature does cost some quality; the gate is "
                "that it stays inside the bounds the solver has always run on")

    cost = [D[("round (M5)", lv)]["ntet"] / D[("flat (M1)", lv)]["ntet"]
            for lv in LEVELS]
    cl.add("M5", "COST: the cap refinement, and it is level-independent",
           "; ".join(f"{lv}: {D[('round (M5)', lv)]['ntet']:,} vs "
                     f"{D[('flat (M1)', lv)]['ntet']:,} tets (x{c:.2f})"
                     for lv, c in zip(LEVELS, cost)),
           "recorded; the RATIO must be level-independent or the ladder is not "
           "self-similar", max(cost) - min(cost) < 0.05,
           note="h_tip scales with h_wall, so the cap costs the same FRACTION at "
                "every level -- the refinement ray is preserved. Measured "
                "against M1's SELF-SIMILAR coarse (coarse_ss); against its "
                "shipped clamped `coarse` the ratio would read a spurious 1.07, "
                "because that mesh has an over-refined far field of its own")

    rows = [(fam, lv, D[(fam, lv)]["ntet"], f"{D[(fam, lv)]['seam_max']:.2f}",
             f"{D[(fam, lv)]['seam_p99']:.2f}", D[(fam, lv)]["n_seam"],
             f"{D[(fam, lv)]['z_wall_max']:.6f}",
             f"{D[(fam, lv)]['min_dihedral']:.2f}",
             f"{D[(fam, lv)]['max_aspect']:.1f}")
            for fam in FAMILIES for lv in LEVELS]
    write_csv(OUT, "m5_tip_cap.csv",
              "family,level,n_tets,seam_crease_max_deg,seam_crease_p99_deg,"
              "n_seam_edges,wall_z_max,min_dihedral_deg,max_aspect_ratio", rows)

    # ---- figure --------------------------------------------------------------
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.4))
    h = np.array([D[("flat (M1)", lv)]["ntet"] for lv in LEVELS]) ** (1 / 3)

    ax.plot(h, flat99, "x--", color=CRITICAL, lw=1.8, ms=8,
            label=f"flat cap (M1): a real edge (q = {q_flat:+.2f})")
    ax.plot(h, rnd99, "o-", color=S2_AQUA, lw=1.9, ms=6,
            label=f"round cap (M5): O(h) faceting (q = {q_rnd:+.2f})")
    ax.plot(h, rnd99[0] * h[0] / h, ":", color=INK_2, lw=1.1,
            label="O(h) reference")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"mesh density $n_{\rm tet}^{1/3}\propto 1/h$")
    ax.set_ylabel("wall crease on the tip seam, p99 [deg]")
    ax.set_title("A sharp edge does not refine away; a smooth surface does",
                 loc="left")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, which="both", alpha=0.25)

    for fam, color, mk in (("flat (M1)", CRITICAL, "x"),
                           ("round (M5)", S1_BLUE, "o")):
        p = tip_section_outline(M[(fam, "medium")])
        ax2.plot(p[:, 2], p[:, 1], mk, color=color, ms=3.0, alpha=0.75,
                 label=f"{fam} wall nodes")
    ax2.axvline(B_SEMI, color=INK_2, lw=1.0, ls=":")
    ax2.text(B_SEMI, 0.031, " z = B_SEMI", fontsize=7.5, color=INK_2)
    th = np.linspace(0, 2 * np.pi, 200)
    ax2.plot(B_SEMI + TIP_CAP_RADIUS * np.sin(th),
             TIP_CAP_RADIUS * np.cos(th), "-", color=S2_AQUA, lw=1.0,
             label="cap circle, radius $t(x)$")
    ax2.set_aspect("equal")
    ax2.set_xlabel("z  (span)")
    ax2.set_ylabel("y  (thickness)")
    ax2.set_title("The tip, at the max-thickness station (medium mesh)",
                  loc="left")
    ax2.legend(fontsize=8, frameon=False, loc="lower left")
    ax2.grid(True, alpha=0.25)
    finish(fig, OUT, "m5_tip_cap.png")

    return cl.report(OUT, "checks_m5.csv")


if __name__ == "__main__":
    raise SystemExit(main())
