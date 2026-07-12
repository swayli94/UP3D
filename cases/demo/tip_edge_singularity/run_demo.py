"""
Tip-edge singularity probe -- is the rigid-planar-wake free-edge singularity a
WAKE-MODEL defect (present on any consistent discretization) or something the
Track B level-set REPRESENTATION removes?

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
nothing to confound the pure potential-flow 1/r edge signal. It measures the
peak local Mach in the P9 tip-edge box (z/b>0.95, at/just aft of the swept tip TE
corner, chord plane) and asks how it grows coarse->medium, on the CONFORMING path
(solve/newton.py) and the level-set path (solve/picard_ls.py, both the
wake-embedded M1 and wake-free M4 meshes -- conforming cannot run wake-free).
conforming-M1 vs level-set-M1 use the IDENTICAL onera_m6 mesh, so it is a true
same-mesh A/B of the representation change.

Finding (results/summary.csv):
  * The tip-edge PEAK Mach diverges under refinement on ALL THREE paths (max
    grows x1.4-2.3 coarse->medium) while the surrounding field and the wing
    control stay flat (p95/mean/wing ratio ~1) -- the localized-edge-singularity
    signature, now seen at M0.5 with zero transonic machinery, i.e. it is a
    genuine potential-flow feature, not a limiter/shock artifact.
  * The level-set representation does NOT remove or blunt it: the LS tip peak is
    at least as large as conforming on the same mesh, and sits in the +2.9%
    straddler cells at/just beyond the geometric tip (the jump terminates
    mid-element there).
  => It is a WAKE-MODEL defect. Track B changes the wake REPRESENTATION, not the
     model (B7 keeps the same rigid planar sheet ending at the tip with
     Gamma(tip)->0). The model-level fix is wake roll-up / an explicit tip
     vortex, which NO current Track B phase does (B9 free-wake is shelved and is
     about O(theta^2) deflection, not roll-up).

Heavy solves cache to results/*.npz (gitignored; ~20 min total to regenerate);
the committed evidence is the figure + CSVs. Standalone:
  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
  python cases/demo/tip_edge_singularity/run_demo.py
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
LEVELS = ("coarse", "medium")
TIP_AFT, TIP_Y = 0.30, 0.03      # tip-edge box aft/chord-plane extents


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
    return mesh.nodes, mesh.elements, np.sqrt(mvop.element_mach2(r["phi_ext"], M)), \
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


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    cl = CheckList("tip-edge singularity: wake MODEL vs REPRESENTATION")

    data, rows = {}, []
    for tag in SOLVERS:
        for lvl in LEVELS:
            d = get(tag, lvl)
            tv, wv = boxes(d)
            data[(tag, lvl)] = dict(
                ntet=int(d["ntet"]), h=float(d["ntet"]) ** (-1.0 / 3.0),
                tip_max=float(tv.max()), tip_p95=float(np.percentile(tv, 95)),
                tip_mean=float(tv.mean()), wing_max=float(wv.max()))
            r = data[(tag, lvl)]
            rows.append((tag, lvl, r["ntet"], f"{r['tip_max']:.3f}",
                         f"{r['tip_p95']:.3f}", f"{r['tip_mean']:.3f}",
                         f"{r['wing_max']:.3f}"))

    # ---- checks (the argument, as pass/fail) --------------------------------
    for tag in SOLVERS:
        a, b = data[(tag, "coarse")], data[(tag, "medium")]
        tip_g = b["tip_max"] / a["tip_max"]
        wing_g = b["wing_max"] / a["wing_max"]
        p95_g = b["tip_p95"] / a["tip_p95"]
        # (1) the wing control is mesh-converged (real flow is bounded)
        cl.add("TIP", f"{tag}: wing-control peak converges (real flow bounded)",
               f"x{wing_g:.2f}", "0.9 < ratio < 1.3", 0.9 < wing_g < 1.3)
        # (2) the tip-edge PEAK diverges -- singularity present
        cl.add("TIP", f"{tag}: tip-edge PEAK Mach diverges under refinement",
               f"max x{tip_g:.2f} (p95 x{p95_g:.2f})", "max ratio > 1.25",
               tip_g > 1.25,
               note="p95/mean flat => a localized edge singularity, not a "
                    "broad region")
    # (3) the decisive same-mesh A/B: LS is NOT blunter than conforming
    lg = data[("ls-M1", "medium")]["tip_max"] / data[("ls-M1", "coarse")]["tip_max"]
    cg = data[("conf-M1", "medium")]["tip_max"] / data[("conf-M1", "coarse")]["tip_max"]
    cl.add("TIP", "same-mesh A/B: level-set does NOT blunt the edge singularity",
           f"LS x{lg:.2f} vs conforming x{cg:.2f}", "LS growth >= conforming",
           lg >= cg - 1e-9,
           note="=> WAKE-MODEL defect; the representation change does not "
                "remove it (LS peak sits in the straddler cells beyond the tip)")

    # ---- figure -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    colors = {"conf-M1": INK_2, "ls-M1": S1_BLUE, "ls-M4": S2_AQUA}
    labels = {"conf-M1": "conforming (M1)", "ls-M1": "level-set (M1, same mesh)",
              "ls-M4": "level-set (M4, wake-free)"}
    for tag in SOLVERS:
        hh = np.array([1.0 / data[(tag, l)]["h"] for l in LEVELS])
        tp = np.array([data[(tag, l)]["tip_max"] for l in LEVELS])
        wg = np.array([data[(tag, l)]["wing_max"] for l in LEVELS])
        ax.plot(hh, tp, "-o", color=colors[tag], lw=1.8, ms=7,
                label=f"{labels[tag]} — tip edge")
        ax.plot(hh, wg, ":s", color=colors[tag], lw=1.2, ms=5, alpha=0.6,
                label=f"{labels[tag]} — wing (control)")
    ax.axhline(1.0, color=CRITICAL, lw=1.0, ls="--")
    ax.text(ax.get_xlim()[0], 1.01, "sonic", color=CRITICAL, fontsize=8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"mesh density  $n_{\mathrm{tet}}^{1/3}\ \propto\ 1/h$")
    ax.set_ylabel("peak local Mach (UNLIMITED)")
    ax.set_title("Tip-edge singularity, ONERA M6, M$_\\infty$=0.5 (subsonic, no "
                 "limiter)\ntip-edge PEAK diverges on BOTH paths; wing control "
                 "stays flat")
    ax.legend(frameon=False, fontsize=8, ncol=1)
    finish(fig, OUT, "tip_edge_growth.png")

    write_csv(OUT, "summary.csv",
              "path,level,n_tet,tip_max,tip_p95,tip_mean,wing_max", rows)
    return cl.report(OUT)


if __name__ == "__main__":
    sys.exit(main())
