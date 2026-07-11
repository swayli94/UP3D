"""
Track B / B3 + B4 demo -- the level-set embedded wake, visualized.

The mesh topology knows NOTHING about the wake. A level set cuts through it,
the cut elements are enriched with auxiliary DOFs, and lift emerges. This demo
makes every part of that visible and self-checks it.

  1. flowfield_lift_vs_nolift.png (M0 embedded) and
     flowfield_lift_vs_nolift_m3.png (M3 wake-free) -- NACA0012 at alpha = 0 (no
     lift) and alpha = 4 (lift), on the SAME mesh with the SAME level set, on
     BOTH mesh families. Speed field (smooth flow off the TE = Kutta) and the
     perturbation potential rendered DISCONTINUOUSLY, i.e. exactly as the
     multivalued DOFs store it: a crisp branch cut carrying [phi] = Gamma at
     alpha = 4, and NO jump at all at alpha = 0. The M3 panel shows the coarser
     wake-free triangulation the level set cuts through generically.

  3. levelset_region.png -- WHERE the level set acts. The unstructured mesh
     near the TE and along the wake corridor, with each element coloured by
     its role: cut (enriched, assembled twice), below-TE fan (its TE reference
     is the aux DOF), and the TE's wall-adjacent upper/lower control volumes
     (where the B4 Kutta condition lives). Nothing else in the mesh is touched.

  4. wake_jump.png -- HOW the jump survives from the TE to the far field. The
     nodal jump [phi] at every cut node vs downstream distance: the g1+g2 wake
     LS convects it unchanged (constant = Gamma) all the way out, and the
     far-field aux DOFs are left FREE so it exits rather than being drained.
     Also shows the multivalued STORAGE: main vs aux DOF values at the sheet.

  5. wall_cp.png -- surface Cp at both incidences on BOTH mesh families
     (solid = M0 embedded, dotted = M3 wake-free, colour = surface, grey dashed
     = conforming reference), inverted axis. Uses the D11 per-side DOF mapping
     (lower-surface TE triangles must read the TE's AUX value).

  6. dual_mesh_embedded_vs_free.png -- the SAME level-set path on the
     wake-EMBEDDED M0 mesh (which HAS a `wake` tag, wake nodes lying exactly on
     the sheet) and on the wake-FREE M3 mesh (NO `wake` tag anywhere, generic
     cuts through generic elements -- the actual Track B workflow target). Same
     lift to ~2%. This is the payoff: no pre-embedded wake surface is needed.

Standalone + self-checking:  python cases/demo/b3_levelset_lifting/run_demo.py
Outputs: cases/demo/b3_levelset_lifting/results/{*.png, summary.csv, checks.csv}
Exit code 0 iff every acceptance check passes. Runtime ~2 min.

Dual-mesh (roadmap Track B working rule): the flow field (fig 1-2) and Cp
(fig 5) run on BOTH the wake-embedded M0 family (which also permits the strict
same-mesh A/B against the conforming solver) and the wake-free M3 family (where
no conforming counterpart exists at all); the enrichment map and wake-jump
figures (3-4) are shown on M0. Runtime ~2 min.
"""

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, DIV_BLUE_RED, INK, INK_2, MESH_DIR, MUTED,
    REPO_ROOT, S1_BLUE, S2_AQUA, S3_YELLOW, S5_VIOLET, SEQ_BLUE,
    CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from pyfp3d.constraints.dirichlet import freestream_phi  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.physics.isentropic import pressure_coefficient_incompressible  # noqa: E402
from pyfp3d.post.section_cut import section_cut  # noqa: E402
from pyfp3d.post.surface import (  # noqa: E402
    triangle_tangential_gradients, wall_outward_normals,
)
from pyfp3d.solve.picard import solve_laplace_lifting  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
NACA_DIR = MESH_DIR / "naca0012_2.5d"                 # wake-EMBEDDED ("C-grid")
NACA_FREE_DIR = MESH_DIR / "naca0012_wakefree_2.5d"   # wake-FREE ("O-grid")
LEVEL = "medium"
ALPHAS = (0.0, 4.0)


# ---------------------------------------------------------------------------
# setup
# ---------------------------------------------------------------------------
def build(mesh):
    """The level set is a HORIZONTAL ruled sheet from the TE line (chord
    plane, design.md Sec 4). The mesh is never modified."""
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return wls, cm, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                        levelset=wls)


def solve_all(mesh):
    cases = {}
    for a in ALPHAS:
        wls, cm, mvop = build(mesh)
        t0 = time.time()
        r = solve_multivalued_lifting(mvop, mesh, 0.0, alpha_deg=a)
        cases[a] = {"wls": wls, "cm": cm, "mvop": mvop, "res": r,
                    "wall_s": time.time() - t0}
        print(f"    alpha={a:.0f}: Gamma = {r['gamma']:+.5f}  "
              f"({r['n_outer']} outers, {cases[a]['wall_s']:.0f}s)")
    return cases


# ---------------------------------------------------------------------------
# 1. flow field -- lift vs no lift, with the jump rendered as it is STORED
# ---------------------------------------------------------------------------
def _section(mesh, cm, mvop, phi_ext, alpha):
    """2D symmetry-plane section + the two side potentials, plus the ORIGINAL
    node ids (smuggled through as a field so we can index the DOF tables).

    The SPEED is taken per node from its OWN side. phi_up is only valid on the
    "+" side (below the wake it mixes an upper aux value into lower elements
    and manufactures a fake gradient), so |grad phi_up| would paint a spurious
    streak along the wake. Note |grad phi| itself is CONTINUOUS across the
    sheet -- only phi jumps -- so the own-side speed is a genuine smooth
    nodal field."""
    from pyfp3d.post.surface import nodal_gradient_recovery

    phi_up, phi_lo = mvop.side_potentials(phi_ext)
    pert_u = phi_up - freestream_phi(mesh.nodes, alpha)
    pert_l = phi_lo - freestream_phi(mesh.nodes, alpha)
    qu = np.linalg.norm(
        nodal_gradient_recovery(mesh.nodes, mesh.elements, phi_up), axis=1)
    ql = np.linalg.norm(
        nodal_gradient_recovery(mesh.nodes, mesh.elements, phi_lo), axis=1)
    q = np.where(cm.node_side == 1, qu, ql)          # own side
    ids = np.arange(len(mesh.nodes), dtype=np.float64)
    return section_cut(mesh, {"pu": pert_u, "pl": pert_l, "q": q, "id": ids},
                       z=float(mesh.nodes[:, 2].min()))


def _explode(sec, wls):
    """Render the section as a DISCONTINUOUS (per-element) field -- which is
    exactly how the multivalued method stores it. Each triangle takes the side
    field of the side it lies on, so the wake jump appears as the crisp
    discontinuity it physically is, instead of being smeared by nodal
    averaging."""
    p = sec.points2d
    tri = sec.triangles
    cen = p[tri].mean(axis=1)
    cen3 = np.column_stack([cen, np.full(len(cen), sec.z)])
    s, _, _ = wls.evaluate(cen3)
    upper = s > 0.0
    vals = np.where(upper[:, None], sec.fields["pu"][tri],
                    sec.fields["pl"][tri])          # (n_tri, 3)
    verts = p[tri]                                   # (n_tri, 3, 2)
    return verts, vals


def demo_flowfield(mesh, cases, checks, *, tag, mesh_label, out_name, step):
    print(f"[{step}] flow field ({mesh_label}): lift vs no lift "
          f"(same mesh, same level set)")
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.4))
    # ONE colour scale for both rows, so "no lift" reads as a genuinely flat
    # field instead of an auto-scaled picture of round-off.
    lim = 0.75 * abs(cases[max(ALPHAS)]["res"]["gamma"])

    for row, a in enumerate(ALPHAS):
        c = cases[a]
        sec = _section(mesh, c["cm"], c["mvop"], c["res"]["phi_ext"], a)
        gamma = c["res"]["gamma"]

        ax = axes[row, 0]
        t = ax.tricontourf(sec.points2d[:, 0], sec.points2d[:, 1],
                           sec.triangles, sec.fields["q"],
                           levels=np.linspace(0.0, 1.6, 25), cmap=SEQ_BLUE,
                           extend="max")
        ax.set_xlim(-0.35, 1.65)
        ax.set_ylim(-0.7, 0.7)
        ax.set_aspect("equal")
        ax.set_title(f"alpha = {a:.0f} deg   speed |grad phi| / U_inf "
                     f"(own side)")
        ax.set_ylabel("y/c")
        ax.set_xlabel("x/c" if row == 1 else "")
        ax.tick_params(labelbottom=(row == 1))
        ax.grid(False)
        cb = fig.colorbar(t, ax=ax, shrink=0.88, pad=0.02)
        cb.ax.tick_params(labelsize=8, colors=MUTED)
        cb.outline.set_edgecolor(BASELINE)

        ax = axes[row, 1]
        verts, vals = _explode(sec, c["wls"])
        pc = PolyCollection(verts, array=vals.mean(axis=1), cmap=DIV_BLUE_RED,
                            edgecolors="none")
        pc.set_clim(-lim, lim)
        ax.add_collection(pc)
        ax.plot([1.0, 5.0], [0.0, 0.0], color=INK, lw=1.1, ls="--", zorder=3)
        ax.set_xlim(-1.5, 5.0)
        ax.set_ylim(-2.2, 2.2)
        ax.set_aspect("equal")
        jump = "NO jump ([phi] = 0)" if abs(gamma) < 1e-3 \
            else f"branch cut carries [phi] = Gamma = {gamma:.3f}"
        ax.set_title(f"alpha = {a:.0f} deg   phi - phi_inf   ({jump})")
        ax.set_xlabel("x/c" if row == 1 else "")
        ax.tick_params(labelbottom=(row == 1))
        ax.grid(False)
        cb = fig.colorbar(pc, ax=ax, shrink=0.88, pad=0.02)
        cb.ax.tick_params(labelsize=8, colors=MUTED)
        cb.outline.set_edgecolor(BASELINE)

    fig.suptitle(f"Level-set embedded wake on the {mesh_label} mesh: the SAME "
                 "mesh and the SAME level set, with and without lift\n"
                 "(right column drawn per-element = exactly how the "
                 "multivalued DOFs store the field)",
                 fontsize=12.5, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.92), h_pad=2.4)
    finish(fig, OUT, out_name)

    g0, g4 = cases[0.0]["res"]["gamma"], cases[4.0]["res"]["gamma"]
    checks.add(f"no-lift[{tag}]", f"Gamma at alpha = 0 ({mesh_label})",
               f"{g0:.2e}", "|Gamma| < 1e-3 (symmetric => no circulation)",
               bool(abs(g0) < 1e-3))
    checks.add(f"lift[{tag}]", f"Gamma at alpha = 4 ({mesh_label})",
               f"{g4:.4f}", "> 0.2 (the level-set path develops circulation)",
               bool(g4 > 0.2))


# ---------------------------------------------------------------------------
# 2. where the level set acts
# ---------------------------------------------------------------------------
def _tri_class(mesh, cm, mvop, sec):
    """Classify each SECTION triangle by the role of the tet it belongs to,
    using the real 3D census (not a re-derivation)."""
    el = np.asarray(mesh.elements, dtype=np.int64)
    ids = sec.fields["id"].astype(np.int64)
    tri_nodes = ids[sec.triangles]                       # (n_tri, 3) global

    def faces_of(elems):
        f = set()
        for e in elems:
            n = el[e]
            for drop in range(4):
                f.add(tuple(sorted(np.delete(n, drop))))
        return f

    cut_f = faces_of(cm.cut_elems)
    tel_f = faces_of(cm.te_lower_elems)

    kind = np.zeros(len(tri_nodes), dtype=np.int64)      # 0 = untouched
    for i, t in enumerate(tri_nodes):
        key = tuple(sorted(t))
        if key in cut_f:
            kind[i] = 1
        elif key in tel_f:
            kind[i] = 2
    return kind


def demo_levelset_region(mesh, cases, checks):
    print("[3/6] where the level set acts (mesh + enriched elements)")
    c = cases[4.0]
    cm, mvop = c["cm"], c["mvop"]
    sec = _section(mesh, cm, mvop, c["res"]["phi_ext"], 4.0)
    kind = _tri_class(mesh, cm, mvop, sec)
    verts = sec.points2d[sec.triangles]

    styles = {
        0: (None, "untouched mesh (single-valued)"),
        1: (S2_AQUA, "CUT: enriched, assembled twice (upper + lower copy)"),
        2: (S3_YELLOW, "below-TE fan: TE reference = AUX dof"),
    }
    # The B4 TE-Kutta control volumes are WALL-adjacent tets; they generally
    # have no face on the symmetry plane, so mark them by centroid instead of
    # trying to colour a section triangle.
    el3 = np.asarray(mesh.elements, dtype=np.int64)
    cen3 = mesh.nodes[el3].mean(axis=1)
    cv_up = np.concatenate([cv["upper_elems"] for cv in mvop._te_cv])
    cv_lo = np.concatenate([cv["lower_elems"] for cv in mvop._te_cv])

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4))
    for ax, (xlim, ylim, ttl) in zip(axes, [
        ((0.93, 1.10), (-0.085, 0.085), "zoom: the trailing edge"),
        ((0.6, 6.0), (-0.9, 0.9), "the whole wake corridor, TE -> far field"),
    ]):
        pc = PolyCollection(verts, facecolors="none", edgecolors=BASELINE,
                            linewidths=0.35, zorder=1)
        ax.add_collection(pc)
        for k, (col, _lab) in styles.items():
            if col is None:
                continue
            m = kind == k
            if not m.any():
                continue
            ax.add_collection(PolyCollection(
                verts[m], facecolors=col, edgecolors=INK_2, linewidths=0.3,
                alpha=0.85, zorder=2))
        for e, col, mk in ((cv_up, CRITICAL, "^"), (cv_lo, S1_BLUE, "v")):
            ax.plot(cen3[e, 0], cen3[e, 1], mk, ms=9, mfc=col, mec=INK_2,
                    mew=0.6, ls="", zorder=6)
        ax.plot([1.0, 8.0], [0.0, 0.0], color=INK, lw=1.2, ls="--", zorder=3)
        te = np.array([1.0, 0.0])
        ax.plot(*te, "k*", ms=14, zorder=5)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal")
        ax.set_title(ttl)
        ax.set_xlabel("x/c")
        ax.grid(False)
    axes[0].set_ylabel("y/c")

    handles = [Line2D([], [], marker="s", ls="", ms=9, mfc=col, mec=INK_2,
                      label=lab) for k, (col, lab) in styles.items()
               if col is not None]
    handles += [
        Line2D([], [], marker="^", ls="", ms=9, mfc=CRITICAL, mec=INK_2,
               label="B4 TE Kutta control volume: UPPER (wall-adjacent)"),
        Line2D([], [], marker="v", ls="", ms=9, mfc=S1_BLUE, mec=INK_2,
               label="B4 TE Kutta control volume: LOWER (wall-adjacent)"),
        Line2D([], [], marker="s", ls="", ms=9, mfc="none", mec=BASELINE,
               label=styles[0][1]),
        Line2D([], [], color=INK, ls="--", label="level set (wake sheet)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, -0.02))

    s = cm.summary()
    fig.suptitle(
        "The level set touches ONE layer of elements. The mesh is never "
        "modified.\n"
        "(the cut layer sits just BELOW the sheet: the eps side-shift sends "
        "on-sheet nodes '+', i.e. the sheet effectively lies at y = -eps)\n"
        f"{s['n_cut_elems']} cut tets of {len(mesh.elements)} "
        f"({100*s['n_cut_elems']/len(mesh.elements):.1f}%), "
        f"{s['n_ext_dofs']} aux DOFs of {s['n_main']} main "
        f"({100*s['n_ext_dofs']/s['n_main']:.1f}%)",
        fontsize=12.5, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0.09, 1, 0.90))
    finish(fig, OUT, "levelset_region.png")

    checks.add("enrichment-cost", "aux DOFs / main DOFs",
               f"{100*s['n_ext_dofs']/s['n_main']:.1f}%",
               "< 15% (the enrichment is a thin strip)",
               bool(s["n_ext_dofs"] / s["n_main"] < 0.15))
    checks.add("te-kutta-cv", "TE Kutta control volumes are wall-adjacent",
               f"{len(mvop._te_cv[0]['upper_elems'])}u / "
               f"{len(mvop._te_cv[0]['lower_elems'])}l",
               "both non-empty", bool(
                   len(mvop._te_cv[0]["upper_elems"]) > 0
                   and len(mvop._te_cv[0]["lower_elems"]) > 0))


# ---------------------------------------------------------------------------
# 3. how the jump survives to the far field + how it is STORED
# ---------------------------------------------------------------------------
def demo_wake_jump(mesh, cases, checks):
    print("[4/6] the wake jump: convected TE -> far field, and how it is stored")
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.0))

    ax = axes[0]
    for a, col in ((0.0, S5_VIOLET), (4.0, S1_BLUE)):
        c = cases[a]
        cm, mvop = c["cm"], c["mvop"]
        cut_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
        jump = mvop.node_jump(c["res"]["phi_ext"], cut_nodes)
        d = cm.d[cut_nodes]
        m = d > 0
        ax.plot(d[m], jump[m], ".", ms=3.5, color=col, alpha=0.55,
                label=f"alpha = {a:.0f} deg   nodal [phi]")
        g = c["res"]["gamma"]
        ax.axhline(g, color=col, ls="--", lw=1.6,
                   label=f"            Gamma = {g:+.4f} (from the TE Kutta)")
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set_xlabel("downstream distance from the TE,  d / c")
    ax.set_ylabel("nodal jump  [phi] = phi_u - phi_l")
    ax.set_title("g1 + g2 CONVECT the jump: it is constant from the TE\n"
                 "to the far field (and the far-field aux DOFs stay FREE)")
    ax.legend(fontsize=8.5, loc="center right")

    # storage: main vs aux at the sheet nodes
    ax = axes[1]
    c = cases[4.0]
    cm, mvop = c["cm"], c["mvop"]
    x = c["res"]["phi_ext"]
    sheet = np.flatnonzero((cm.ext_dof_of_node >= 0) & (np.abs(cm.s_raw) < 1e-9)
                           & (cm.d > 0))
    order = np.argsort(cm.d[sheet])
    sheet = sheet[order]
    aux = cm.ext_dof_of_node[sheet]
    # PERTURBATION potential -- phi itself grows like U*x (~15 here) and would
    # hide the jump; subtracting the freestream makes the Gamma gap visible.
    pert = freestream_phi(mesh.nodes, 4.0)
    ax.plot(cm.d[sheet], x[sheet] - pert[sheet], "-", color=CRITICAL, lw=2.2,
            label="MAIN dof  (the node's OWN side = upper)")
    ax.plot(cm.d[sheet], x[aux] - pert[sheet], "-", color=S1_BLUE, lw=2.2,
            label="AUX dof   (the OTHER side = lower)")
    ax.fill_between(cm.d[sheet], x[aux] - pert[sheet], x[sheet] - pert[sheet],
                    color=S2_AQUA, alpha=0.45,
                    label=f"[phi] = Gamma = {c['res']['gamma']:.4f}"
                          "   <- what the aux DOFs buy")
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set_xlabel("downstream distance from the TE,  d / c")
    ax.set_ylabel("phi - phi_inf   (perturbation potential)")
    ax.set_title("Multivalued STORAGE (Lopez eq. 3.33-3.34): ONE mesh, one\n"
                 "extra dof per cut node -- NOT two meshes")
    ax.legend(fontsize=8.5, loc="upper right")

    fig.suptitle("The wake jump: carried by the aux DOFs, convected by the "
                 "wake LS, valued by the TE Kutta",
                 fontsize=12.5, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    finish(fig, OUT, "wake_jump.png")

    c = cases[4.0]
    cm, mvop = c["cm"], c["mvop"]
    cut_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
    jump = mvop.node_jump(c["res"]["phi_ext"], cut_nodes)
    d = cm.d[cut_nodes]
    near = jump[(d > 0.05) & (d < 1.0)].mean()
    far = jump[d > 5.0].mean()
    drift = abs(far - near) / abs(near)
    checks.add("jump-convected", "|[phi](far) - [phi](near)| / [phi]",
               f"{drift*100:.1f}%",
               "< 10% (no drain: the far-field aux DOFs are free)",
               bool(drift < 0.10))
    checks.add("jump-is-gamma", "[phi] near the TE vs the reported Gamma",
               f"{abs(near - c['res']['gamma'])/abs(c['res']['gamma'])*100:.1f}%",
               "< 10%",
               bool(abs(near - c["res"]["gamma"]) / abs(c["res"]["gamma"]) < 0.10))


# ---------------------------------------------------------------------------
# 4. wall Cp (with the D11 per-side mapping) vs the conforming path
# ---------------------------------------------------------------------------
def _cp_b(mesh, cm, mvop, phi_ext):
    """Wall Cp on the level-set path. D11: the lower-surface TE triangles must
    read the TE node's AUX dof, not its main dof -- using phi_main alone makes
    cl_pressure junk (measured -3.35)."""
    wall = mesh.boundary_faces["wall"]
    phi_up, phi_lo = mvop.side_potentials(phi_ext)
    n_out = wall_outward_normals(mesh.nodes, mesh.elements, wall)
    gu, area, _ = triangle_tangential_gradients(mesh.nodes, wall, phi_up)
    gl, _, _ = triangle_tangential_gradients(mesh.nodes, wall, phi_lo)
    grad = np.where((n_out[:, 1] > 0)[:, None], gu, gl)
    q2 = np.sum(grad * grad, axis=1)
    cp = np.array([pressure_coefficient_incompressible(v) for v in q2])
    xc = mesh.nodes[wall].mean(axis=1)[:, 0]
    return xc, cp, n_out[:, 1] > 0, cp, area, n_out


def _cl_p(mesh, cp, area, n_out, a):
    lift = np.array([-np.sin(np.radians(a)), np.cos(np.radians(a)), 0.0])
    dz = float(np.ptp(mesh.nodes[:, 2]))
    return float((-(cp * area) @ n_out / dz) @ lift)


def demo_wall_cp(mesh_m0, cases_m0, mesh_m3, cases_m3, checks, step):
    print(f"[{step}] wall Cp (D11 per-side mapping): BOTH meshes vs conforming")
    mesh_cut, wc = cut_wake(mesh_m0)

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2), sharey=True)
    rows = []
    for ax, a in zip(axes, ALPHAS):
        # --- M0 (wake-EMBEDDED): level set (solid) + conforming (grey dashed)
        c0 = cases_m0[a]
        xc0, cp0, up0, _, area0, n0 = _cp_b(mesh_m0, c0["cm"], c0["mvop"],
                                            c0["res"]["phi_ext"])
        rc = solve_laplace_lifting(mesh_cut, wc, alpha_deg=a)
        wall_c = mesh_cut.boundary_faces["wall"]
        n_c = wall_outward_normals(mesh_cut.nodes, mesh_cut.elements, wall_c)
        g_c, area_c, _ = triangle_tangential_gradients(
            mesh_cut.nodes, wall_c, rc["phi"])
        cp_c = np.array([pressure_coefficient_incompressible(v)
                         for v in np.sum(g_c * g_c, axis=1)])
        xc_c = mesh_cut.nodes[wall_c].mean(axis=1)[:, 0]
        up_c = n_c[:, 1] > 0

        # --- M3 (wake-FREE): level set only (no conforming counterpart)
        c3 = cases_m3[a]
        xc3, cp3, up3, _, area3, n3 = _cp_b(mesh_m3, c3["cm"], c3["mvop"],
                                            c3["res"]["phi_ext"])

        # conforming = grey dashed reference underneath
        for m in (up_c, ~up_c):
            o = np.argsort(xc_c[m])
            ax.plot(xc_c[m][o], cp_c[m][o], "--", color=MUTED, lw=1.2,
                    zorder=1)
        # colour = surface (red upper / blue lower); linestyle = mesh
        # (solid = M0 embedded, dotted = M3 wake-free)
        for xc, cp, up, ls, mk in ((xc0, cp0, up0, "-", None),
                                   (xc3, cp3, up3, ":", None)):
            for m, col in ((up, CRITICAL), (~up, S1_BLUE)):
                o = np.argsort(xc[m])
                ax.plot(xc[m][o], cp[m][o], ls, color=col, lw=2.0,
                        zorder=3 if ls == "-" else 4)
        ax.invert_yaxis()
        ax.set_xlabel("x/c")
        ax.set_title(f"alpha = {a:.0f} deg      "
                     f"Gamma:  M0 {c0['res']['gamma']:+.4f}   "
                     f"M3 {c3['res']['gamma']:+.4f}   "
                     f"(conforming {float(np.mean(rc['gamma'])):+.4f})",
                     fontsize=9.5)

        cl0 = _cl_p(mesh_m0, cp0, area0, n0, a)
        cl3 = _cl_p(mesh_m3, cp3, area3, n3, a)
        clc = _cl_p(mesh_cut, cp_c, area_c, n_c, a)
        rows.append((f"{a:.0f}", f"{c0['res']['gamma']:.6f}",
                     f"{c3['res']['gamma']:.6f}",
                     f"{float(np.mean(rc['gamma'])):.6f}",
                     f"{cl0:.6f}", f"{cl3:.6f}", f"{clc:.6f}"))
        if a > 0:
            checks.add("cp-vs-conforming",
                       f"cl_p alpha={a:.0f}: M0 level set vs conforming",
                       f"{cl0:.4f} vs {clc:.4f}",
                       "within 3%", bool(abs(cl0 - clc) / abs(clc) < 0.03))
            checks.add("D11", f"cl_p vs cl_KJ = 2*Gamma (M0, alpha={a:.0f})",
                       f"{cl0:.4f} vs {2*c0['res']['gamma']:.4f}",
                       "within 5% (per-side wall mapping is correct)",
                       bool(abs(cl0 - 2*c0["res"]["gamma"])
                            / abs(2*c0["res"]["gamma"]) < 0.05))
            checks.add("cp-wakefree",
                       f"cl_p alpha={a:.0f}: M3 wake-free vs conforming",
                       f"{cl3:.4f} vs {clc:.4f}",
                       "within 5%", bool(abs(cl3 - clc) / abs(clc) < 0.05))

    axes[0].set_ylabel("Cp   (inverted axis)")
    handles = [
        Line2D([], [], color=CRITICAL, lw=2, label="upper surface"),
        Line2D([], [], color=S1_BLUE, lw=2, label="lower surface"),
        Line2D([], [], color=INK, lw=2, ls="-",
               label="M0 wake-embedded (level set)"),
        Line2D([], [], color=INK, lw=2, ls=":",
               label="M3 wake-free (level set)"),
        Line2D([], [], color=MUTED, lw=1.4, ls="--",
               label="conforming solver (M0, reference)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Surface Cp on BOTH mesh families: the wake-embedded (M0) and "
                 "wake-free (M3) level-set paths\n"
                 "reproduce the conforming solver (alpha = 0 collapses "
                 "upper/lower; alpha = 4 shows the loading)",
                 fontsize=12.2, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0.07, 1, 0.91))
    finish(fig, OUT, "wall_cp.png")

    write_csv(OUT, "summary.csv",
              "alpha_deg,gamma_m0,gamma_m3,gamma_conforming,"
              "cl_p_m0,cl_p_m3,cl_p_conforming", rows)


# ---------------------------------------------------------------------------
# 5. dual-mesh: the SAME level set on a wake-EMBEDDED and a wake-FREE mesh
# ---------------------------------------------------------------------------
def demo_dual_mesh(mesh_m0, cases_m0, mesh_m3, cases_m3, checks, step):
    print(f"[{step}] dual-mesh: embedded (M0) vs wake-free (M3) cut structure")
    r3 = cases_m3[4.0]["res"]

    panels = [
        ("M0  wake-EMBEDDED  (the mesh HAS a `wake` tag;\n"
         "wake nodes lie exactly ON the sheet)", mesh_m0,
         cases_m0[4.0]["cm"], cases_m0[4.0]["mvop"],
         cases_m0[4.0]["res"], "wake" in mesh_m0.boundary_faces),
        ("M3  wake-FREE  (NO `wake` tag anywhere;\n"
         "the level set makes GENERIC cuts)", mesh_m3,
         cases_m3[4.0]["cm"], cases_m3[4.0]["mvop"], r3, False),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.6))
    for ax, (ttl, mesh, cm, mvop, res, embedded) in zip(axes, panels):
        sec = _section(mesh, cm, mvop, res["phi_ext"], 4.0)
        kind = _tri_class(mesh, cm, mvop, sec)
        verts = sec.points2d[sec.triangles]
        ax.add_collection(PolyCollection(verts, facecolors="none",
                                         edgecolors=BASELINE, linewidths=0.4,
                                         zorder=1))
        for k, col in ((1, S2_AQUA), (2, S3_YELLOW)):
            m = kind == k
            if m.any():
                ax.add_collection(PolyCollection(
                    verts[m], facecolors=col, edgecolors=INK_2,
                    linewidths=0.3, alpha=0.85, zorder=2))
        # the wake nodes: on M0 they lie ON the sheet, on M3 they are wherever
        # the generic cut falls.
        cut_nodes = sec.fields["id"].astype(np.int64)[
            np.unique(sec.triangles)]
        onsheet = cut_nodes[cm.ext_dof_of_node[cut_nodes] >= 0]
        p = mesh.nodes[onsheet]
        ax.plot(p[:, 0], p[:, 1], ".", ms=3, color=CRITICAL, zorder=4,
                label="cut-element nodes (carry an aux DOF)")
        ax.plot([1.0, 1.25], [0.0, 0.0], color=INK, lw=1.4, ls="--", zorder=3)
        ax.plot(1.0, 0.0, "k*", ms=13, zorder=5)
        ax.set_xlim(0.95, 1.12)
        ax.set_ylim(-0.05, 0.05)
        ax.set_aspect("equal")
        g = float(res["gamma"])
        ax.set_title(f"{ttl}\nGamma = {g:+.4f}", fontsize=10)
        ax.set_xlabel("x/c")
        ax.legend(fontsize=8, loc="lower left")
    axes[0].set_ylabel("y/c")

    g0 = cases_m0[4.0]["res"]["gamma"]
    fig.suptitle("Same level-set path, two mesh kinds -- the wake-free mesh is "
                 "the Track B workflow target\n"
                 f"embedded Gamma {g0:.4f}  vs  wake-free Gamma {r3['gamma']:.4f}"
                 f"  ({abs(r3['gamma']-g0)/g0*100:.1f}% apart -- "
                 "no `wake` tag, generic cuts, same lift)",
                 fontsize=12.2, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    finish(fig, OUT, "dual_mesh_embedded_vs_free.png")

    checks.add("wake-free", "M3 mesh has NO `wake` tag",
               "wake" in mesh_m3.boundary_faces,
               "False (topology knows nothing about the wake)",
               "wake" not in mesh_m3.boundary_faces)
    rel = abs(r3["gamma"] - g0) / abs(g0)
    checks.add("dual-mesh", "wake-free Gamma vs embedded Gamma",
               f"{r3['gamma']:.4f} vs {g0:.4f} ({rel*100:.1f}%)",
               "within 5% (generic cuts reproduce the embedded result)",
               bool(rel < 0.05))


# ---------------------------------------------------------------------------
def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    checks = CheckList("Track B / B3 + B4 -- level-set embedded wake")

    path0 = NACA_DIR / f"{LEVEL}.msh"
    path3 = NACA_FREE_DIR / "coarse.msh"    # committed; M3 medium is gitignored
    if not path0.exists() or not path3.exists():
        print(f"missing mesh: {path0 if not path0.exists() else path3}")
        return 1

    print(f"[0/6] wake-EMBEDDED {LEVEL} M0 mesh (has a `wake` tag)")
    mesh_m0 = read_mesh(path0)
    cases_m0 = solve_all(mesh_m0)
    print("[0/6] wake-FREE coarse M3 mesh (NO `wake` tag -- the workflow form)")
    mesh_m3 = read_mesh(path3)
    assert "wake" not in mesh_m3.boundary_faces
    cases_m3 = solve_all(mesh_m3)

    demo_flowfield(mesh_m0, cases_m0, checks, tag="M0",
                   mesh_label="M0 (embedded)",
                   out_name="flowfield_lift_vs_nolift.png", step="1/6")
    demo_flowfield(mesh_m3, cases_m3, checks, tag="M3",
                   mesh_label="M3 (wake-free)",
                   out_name="flowfield_lift_vs_nolift_m3.png", step="2/6")
    demo_levelset_region(mesh_m0, cases_m0, checks)
    demo_wake_jump(mesh_m0, cases_m0, checks)
    demo_wall_cp(mesh_m0, cases_m0, mesh_m3, cases_m3, checks, step="5/6")
    demo_dual_mesh(mesh_m0, cases_m0, mesh_m3, cases_m3, checks, step="6/6")

    return checks.report(OUT)


if __name__ == "__main__":
    sys.exit(main())
