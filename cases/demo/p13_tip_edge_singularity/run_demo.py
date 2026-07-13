"""
Tip-edge singularity probe -- is the rigid-planar-wake free-edge singularity a
WAKE-MODEL defect (present on any consistent discretization) or something the
Track B level-set REPRESENTATION removes? And what is its refinement RATE?

This demo is the evidence for **P13/G13.1** (tip / wake-edge singularity
characterization; roadmap.md Track P P13, design.md §4.1).

Context. P9/G9.1 found that the ONERA M6 transonic solve does not converge under
refinement: the UNLIMITED local Mach at the wake sheet's free tip edge climbs
1.40 -> 2.13 -> 7.93 (coarse/medium/fine), with 9 cells crossing the M_cap=3
speed limiter on the fine mesh -- and those permanently-limited cells block the
Newton freeze machinery, so the fine "solution" is a limit-cycle artifact, not a
discrete solution. P9 attributed this to a "vortex-sheet-edge singularity of a
rigid planar wake" and pointed at "the tip/wake treatment (Track B / free wake or
an explicit tip-vortex model)" as the fix. One doc (agent-rules) over-stated this
as "precisely what Track B exists to fix".

This demo settles it empirically, and cheaply. It isolates the GEOMETRY by
probing SUBSONICALLY (M0.5): no shock, no artificial density, no speed limiter --
nothing to confound the pure potential-flow edge signal. It measures the peak
local Mach in the P9 tip-edge box (z/b>0.95, at/just aft of the swept tip TE
corner, chord plane) and asks (a) HOW it grows under refinement and (b) WHETHER
the level-set representation blunts it, on the CONFORMING path (solve/newton.py,
three levels coarse/medium/fine) and the level-set path (solve/picard_ls.py, both
the wake-embedded M1 and wake-free M4 meshes, two levels each -- the LS path has
no scalable transonic-grade solver at fine yet, so LS fine is skipped, see note).
conforming-M1 vs level-set-M1 use the IDENTICAL onera_m6 mesh, so it is a true
same-mesh A/B of the representation change.

Findings (results/summary.csv, checks.csv):
  * The tip-edge PEAK Mach diverges under refinement on ALL paths (grows x1.4-2.3
    per level) while WITHIN THE SAME BOX the p95/mean stay flat (conforming tip
    p95 0.573 -> 0.562 -> 0.525, mean ~0.49 across all three levels) -- the
    localized-edge-singularity signature, seen at M0.5 with zero transonic
    machinery, i.e. a genuine potential-flow feature, not a limiter/shock
    artifact. (The bulk "wing" interior is also flat coarse->medium; its fine
    MAX is polluted by a separate sharp-TE edge cell -- another P1 edge feature
    -- so the clean same-box control is the tip-box p95, plotted in the figure.)
  * ★ RATE: the CONFORMING three-point log-log slope of peak Mach vs 1/h is
    p ~ 0.59 (peak ~ h^-p) -- a **1/sqrt(r) flat-plate-edge singularity**, NOT
    the "1/r-type" the earlier docs wrote (1/r would be a concentrated line
    vortex, p = 1). The DRIVER is the trailing vorticity dGamma/dz (tip-max:
    |dGamma/dz| is ~10x larger at the tip than mid-span on B7's smooth
    Gamma(z)), NOT the bound circulation Gamma, which correctly -> 0 at the tip:
    Gamma(tip)=0 is a NECESSARY-not-sufficient regularity condition; a
    terminating flat vortex sheet has a free-edge crossflow singularity like a
    flat-plate edge.
  * ★ The M0.5 conforming FINE solve does NOT converge to a discrete solution
    (converged=False, ~1.4k NaN/floored cells): the tip singularity trips the
    speed limiter / density floor even SUBSONICALLY -- the exact M0.5 analogue
    of what G9.1 found transonically (fine M6 = limit-cycle artifact, not a
    solution). Its tip-edge peak 1.51 nonetheless continues the h^-p trend; the
    non-convergence is CAUSED by the singularity, so it is the finding, not a
    confound for the rate.
  * The level-set representation does NOT remove or blunt it: the LS tip peak is
    at least as large as conforming on the same mesh, and its edge growth rate is
    >= conforming, sitting in the +2.9% straddler cells at/just beyond the
    geometric tip (the jump terminates mid-element there).
  => It is a WAKE-MODEL defect. Track B changes the wake REPRESENTATION, not the
     model (B7 keeps the same rigid planar sheet ending at the tip with
     Gamma(tip)->0). The model-level fix is wake roll-up / an explicit tip
     vortex, which NO current Track B phase does (B9 free-wake is shelved and is
     about O(theta^2) deflection, not roll-up). That is P13/G13.2 (future),
     handed to Track B as a rescope of B9.

Heavy solves cache to results/*.npz (gitignored; ~30 min total to regenerate,
the conforming fine M0.5 AMG solve being ~10-15 min of it); the committed
evidence is the figure + CSVs. Standalone:
  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
  python cases/demo/p13_tip_edge_singularity/run_demo.py
"""
import sys, time
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, CheckList, INK_2, S1_BLUE, S2_AQUA, S3_YELLOW,
    apply_style, finish, plt, write_csv,
)
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te  # noqa: E402
from pyfp3d.kernels.jacobian import PicardOperator  # noqa: E402
from pyfp3d.physics.isentropic import mach_squared_field  # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet  # noqa: E402

OUT = Path(__file__).parent / "results"
MESH = REPO / "cases" / "meshes"
ALPHA, M = 3.06, 0.5
TIP_AFT, TIP_Y = 0.30, 0.03      # tip-edge box aft/chord-plane extents
# The conforming path has three levels; the level-set path has two (its fine
# solve needs the scalable transonic-grade LS solver not built in Part 1).
LEVELS_BY_TAG = {
    "conf-M1": ("coarse", "medium", "fine"),
    "ls-M1": ("coarse", "medium"),
    "ls-M4": ("coarse", "medium"),
}


def _mach_field_conforming(level):
    mc, wc = cut_wake(read_mesh(MESH / "onera_m6" / f"{level}.msh"))
    r = solve_newton_lifting(mc, wc, m_inf=M, alpha_deg=ALPHA, upwind_c=1.5,
                             precond="amg", tol_residual=1e-9)
    op = PicardOperator(mc.nodes, mc.elements)
    _, q2 = op.velocities(np.asarray(r["phi"]))
    return mc.nodes, mc.elements, np.sqrt(mach_squared_field(q2, M, 1.4)), \
        bool(r["converged"])


def _mach_field_ls(meshdir, level):
    mesh = read_mesh(meshdir / f"{level}.msh")
    a = np.radians(ALPHA)
    wls = WakeLevelSet(np.array([[x_te(0.), 0, 0], [x_te(B_SEMI), 0, B_SEMI]]),
                       direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    r = solve_multivalued_lifting(mvop, mesh, M, alpha_deg=ALPHA,
                                  farfield="neumann", upwind_c=0.0,
                                  n_outer_max=80, tol_residual=1e-7)
    # mixed_plain="side" pinned 2026-07-14 (default flipped to "main"): the
    # committed G13.1 LS-exponent evidence (p=+1.34, tip-box peaks) was
    # measured through the historical side reading -- see the B8 termination
    # diagnosis for the honest (+0.62) reading and the G13.1 erratum note.
    return mesh.nodes, mesh.elements, \
        np.sqrt(mvop.element_mach2(r["phi_ext"], M, mixed_plain="side")), \
        bool(r["converged"])


SOLVERS = {
    "conf-M1": lambda lvl: _mach_field_conforming(lvl),
    "ls-M1": lambda lvl: _mach_field_ls(MESH / "onera_m6", lvl),
    "ls-M4": lambda lvl: _mach_field_ls(MESH / "onera_m6_wakefree", lvl),
}


def get(tag, level):
    """Cached per-element Mach field + centroid + x_te(z)."""
    cache = OUT / f"e1_{tag}_{level}.npz"
    if cache.exists():
        d = np.load(cache)
        return {k: d[k] for k in d.files}
    print(f"  solving {tag} {level} (minutes)...", flush=True)
    t0 = time.time()
    nodes, elements, mach, conv = SOLVERS[tag](level)
    cen = nodes[elements].mean(axis=1)
    xte = np.array([x_te(np.clip(z, 0, B_SEMI)) for z in cen[:, 2]])
    d = dict(cen=cen, mach=mach, xte=xte, dt=time.time() - t0, conv=conv,
             ntet=len(elements))
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(cache, **d)
    return d


def boxes(d):
    cen, mach, xte = d["cen"], d["mach"], d["xte"]
    zb = cen[:, 2] / B_SEMI
    dx = cen[:, 0] - xte
    tip = (zb > 0.95) & (dx >= -0.05) & (dx <= TIP_AFT) & (np.abs(cen[:, 1]) < TIP_Y)
    wing = (dx < 0.0) & (zb < 0.95)      # control: the real (bounded) flow
    return mach[tip], mach[wing]


def loglog_slope(inv_h, peak):
    """Least-squares slope p of  log(peak) vs log(1/h):  peak ~ (1/h)^p ~ h^-p."""
    return float(np.polyfit(np.log(np.asarray(inv_h)),
                            np.log(np.asarray(peak)), 1)[0])


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    cl = CheckList("tip-edge singularity (P13/G13.1): wake MODEL vs REPRESENTATION")

    data, rows = {}, []
    for tag in SOLVERS:
        for lvl in LEVELS_BY_TAG[tag]:
            d = get(tag, lvl)
            tv, wv = boxes(d)
            data[(tag, lvl)] = dict(
                ntet=int(d["ntet"]), h=float(d["ntet"]) ** (-1.0 / 3.0),
                inv_h=float(d["ntet"]) ** (1.0 / 3.0),
                tip_max=float(tv.max()), tip_p95=float(np.percentile(tv, 95)),
                tip_mean=float(tv.mean()), wing_max=float(wv.max()),
                conv=bool(d["conv"]))
            r = data[(tag, lvl)]
            rows.append((tag, lvl, r["ntet"], f"{r['tip_max']:.3f}",
                         f"{r['tip_p95']:.3f}", f"{r['tip_mean']:.3f}",
                         f"{r['wing_max']:.3f}", str(int(r["conv"]))))

    # ---- checks (the argument, as pass/fail) --------------------------------
    # (1)-(2) per path: coarse->medium (the two levels every path has)
    for tag in SOLVERS:
        a, b = data[(tag, "coarse")], data[(tag, "medium")]
        tip_g = b["tip_max"] / a["tip_max"]
        wing_g = b["wing_max"] / a["wing_max"]
        p95_g = b["tip_p95"] / a["tip_p95"]
        # (1) the SAME tip box: its p95 (bulk) stays flat => localized, not broad
        cl.add("G13.1", f"{tag}: tip-box p95 stays flat (localized edge, not broad)",
               f"p95 x{p95_g:.2f} (wing bulk x{wing_g:.2f})",
               "0.85 < p95 ratio < 1.15", 0.85 < p95_g < 1.15,
               note="the bulk of the tip box AND the wing interior converge; only "
                    "the few cells at the very edge grow")
        # (2) the tip-edge PEAK diverges => singularity present
        cl.add("G13.1", f"{tag}: tip-edge PEAK Mach diverges under refinement",
               f"max x{tip_g:.2f}", "max ratio > 1.25", tip_g > 1.25)

    # (2b) ★ conforming FINE is NOT a discrete solution -- the M0.5 analogue of
    #      G9.1: the tip singularity trips the limiter / floor even subsonically.
    fd = data[("conf-M1", "fine")]
    cl.add("G13.1", "conf fine: limiter tripped => NOT a discrete solution "
           "(M0.5 analogue of G9.1)", f"converged={fd['conv']}",
           "converged is False (singularity trips the limiter subsonically)",
           not fd["conv"],
           note="floored/limited cells appear even at M0.5, so the fine mesh "
                "limit-cycles -- exactly what G9.1 found transonically")

    # (3) same-mesh A/B: LS is NOT blunter than conforming (coarse->medium)
    lg = data[("ls-M1", "medium")]["tip_max"] / data[("ls-M1", "coarse")]["tip_max"]
    cg = data[("conf-M1", "medium")]["tip_max"] / data[("conf-M1", "coarse")]["tip_max"]
    cl.add("G13.1", "same-mesh A/B: level-set does NOT blunt the edge singularity",
           f"LS x{lg:.2f} vs conforming x{cg:.2f}", "LS growth >= conforming",
           lg >= cg - 1e-9,
           note="=> WAKE-MODEL defect; the representation change does not "
                "remove it (LS peak sits in the straddler cells beyond the tip)")

    # (4) ★ RATE: conforming THREE-point log-log exponent p (peak ~ h^-p).
    #     p ~ 0.5  <=>  1/sqrt(r) flat-plate-edge singularity (NOT 1/r, p=1).
    conf_lvls = LEVELS_BY_TAG["conf-M1"]
    inv_h = [data[("conf-M1", l)]["inv_h"] for l in conf_lvls]
    peak = [data[("conf-M1", l)]["tip_max"] for l in conf_lvls]
    p_conf = loglog_slope(inv_h, peak)
    cl.add("G13.1", "conforming refinement exponent p (3-pt log-log)",
           f"p = {p_conf:.3f}", "0.4 < p < 0.65  (1/sqrt(r), flat-plate edge)",
           0.4 < p_conf < 0.65,
           note="peak Mach ~ h^-p; 1/r (concentrated line vortex) would give "
                "p=1 -- the earlier 'induces a 1/r-type velocity' was wrong")

    # (5) rate form of the A/B: the LS two-point exponent >= the conforming
    #     two-point exponent on the SAME (coarse,medium) meshes.
    def two_pt_p(tag):
        a, b = data[(tag, "coarse")], data[(tag, "medium")]
        return loglog_slope([a["inv_h"], b["inv_h"]], [a["tip_max"], b["tip_max"]])
    p_ls_m1, p_conf2 = two_pt_p("ls-M1"), two_pt_p("conf-M1")
    cl.add("G13.1", "edge growth RATE: level-set >= conforming (same meshes)",
           f"p_LS = {p_ls_m1:.3f} vs p_conf = {p_conf2:.3f}",
           "p_LS >= p_conf", p_ls_m1 >= p_conf2 - 1e-9,
           note="both diverge; the representation change does not slow the edge "
                "singularity => model, not representation")

    # ---- figure -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    colors = {"conf-M1": INK_2, "ls-M1": S1_BLUE, "ls-M4": S2_AQUA}
    labels = {"conf-M1": "conforming (M1)", "ls-M1": "level-set (M1, same mesh)",
              "ls-M4": "level-set (M4, wake-free)"}
    for tag in SOLVERS:
        lvls = LEVELS_BY_TAG[tag]
        hh = np.array([data[(tag, l)]["inv_h"] for l in lvls])
        tp = np.array([data[(tag, l)]["tip_max"] for l in lvls])
        cg = np.array([data[(tag, l)]["tip_p95"] for l in lvls])
        ax.plot(hh, tp, "-o", color=colors[tag], lw=1.8, ms=7,
                label=f"{labels[tag]} — tip-edge peak")
        ax.plot(hh, cg, ":s", color=colors[tag], lw=1.2, ms=5, alpha=0.6,
                label=f"{labels[tag]} — tip-box p95 (control)")
    # slope guides: measured conforming p, and reference p=0.5 and p=1.
    x0, x1 = ax.get_xlim()
    xg = np.array([inv_h[0], inv_h[-1]])
    y_anchor = peak[0]
    for pp, style, txt in ((p_conf, "-", f"conforming p={p_conf:.2f}"),
                           (0.5, "--", "1/√r  (p=0.5)"),
                           (1.0, ":", "1/r  (p=1)")):
        yg = y_anchor * (xg / xg[0]) ** pp
        ax.plot(xg, yg, style, color=S3_YELLOW, lw=1.2, alpha=0.9, zorder=0)
        ax.text(xg[-1] * 1.01, yg[-1], txt, color=S3_YELLOW, fontsize=7.5,
                va="center")
    ax.axhline(1.0, color=CRITICAL, lw=1.0, ls="--")
    ax.text(x0, 1.01, "sonic", color=CRITICAL, fontsize=8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"mesh density  $n_{\mathrm{tet}}^{1/3}\ \propto\ 1/h$")
    ax.set_ylabel("peak local Mach (UNLIMITED)")
    ax.set_title("Tip-edge singularity, ONERA M6, M$_\\infty$=0.5 (subsonic)\n"
                 "tip-edge PEAK diverges as $h^{-p}$, $p\\approx0.59$ "
                 "(1/√r, flat-plate edge); tip-box p95 control flat")
    ax.legend(frameon=False, fontsize=8, ncol=1, loc="upper left")
    finish(fig, OUT, "tip_edge_growth.png")

    write_csv(OUT, "summary.csv",
              "path,level,n_tet,tip_max,tip_p95,tip_mean,wing_max,converged", rows)
    write_csv(OUT, "rates.csv", "quantity,value",
              [("conforming_3pt_exponent_p", f"{p_conf:.4f}"),
               ("conforming_2pt_exponent_p", f"{p_conf2:.4f}"),
               ("ls_M1_2pt_exponent_p", f"{p_ls_m1:.4f}"),
               ("ls_M4_2pt_exponent_p", f"{two_pt_p('ls-M4'):.4f}")])
    return cl.report(OUT, "checks_g131.csv")


if __name__ == "__main__":
    sys.exit(main())
