"""
Track B / B17: the far-field aux pin must carry jump = GAMMA, not 0.

B16 cured the LS Newton wing-body churn by pinning the far-field-boundary aux
(the outer nodes a wake level set crosses on its way out -- no outflow clip) to
remove the near-singular outer wake-LS rows. On farfield="freestream" it pinned
them to jump=0. ★ B17 shows that jump=0 REMOVES the wake circulation the outflow
physically carries -> a resolution-dependent lift error. The GB16.4 "medium
Newton-pin stalls 22% low" that B16 flagged as an UNRESOLVED non-convergence is
actually this BC-modelling error: BOTH the Newton pin AND an independent Picard
pin converge to the SAME wrong 0.169 (so it is not a solver stall), and the
coarse "match to conforming" was a coincidence (jump=0 there cancelled the coarse
legacy's outer-tet garbage). The fix is farfield_aux="pin_gamma" (the new
default): aux = host phi_inf - side*gamma (jump=gamma), the same Dirichlet
conditioning cure with the physically correct ring value.

The measured triangle (cl_p wing, ONERA M6 wing-body, M0.5, alpha 3.06):

                    coarse    medium
  conforming        0.2089    0.2173     (P14 Newton, the reference)
  legacy free-aux   0.1853    0.2165     coarse polluted by |jump|=53 outer tet
  pin  jump=0 (B16) 0.2086    0.1690     NON-monotone: kills the outflow circ
  pin_gamma (B17)   0.2087    0.2117/15  monotone to conforming, both solvers

Parts:
  1 (coarse, ungated)  GB17.1 the coarse 4-way triangle + the far-field jump
                       collapse (legacy 53 -> pin 0 -> pin_gamma gamma).
  2 (coarse, ungated)  GB17.2 post-processing audit: cl_p (surface integral) and
                       cl_KJ (circulation integral) move TOGETHER, so the pin
                       gap is a flow-state change, NOT a post artifact; and the
                       plotted sectional cl(z) is the Gamma-based 2*Gamma/(u*c).
  3 (gated, medium)    GB17.3 the medium triangle: pin jump=0 lands 0.169 on
                       BOTH solvers (the BC error, not a stall); GB17.4 pin_gamma
                       closes it (Newton 0.2115 ~ Picard 0.2117 ~ conforming,
                       monotone); GB17.5 the spanwise Gamma(z) uniform-offset is
                       removed.

Cost: ~2-3 min coarse; medium under PYFP3D_TRANSONIC_GATES=1 (four solves incl.
one Newton, lagged-LU), cached to gitignored results/*.npz. Meshes are
gitignored -> a level whose mesh is absent is skipped.

Run:  python cases/demo/b17_farfield_pin_gamma/run_demo.py
      PYFP3D_TRANSONIC_GATES=1 NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 \
        OPENBLAS_NUM_THREADS=16 python cases/demo/b17_farfield_pin_gamma/run_demo.py
"""

import os
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cases.demo._common import (CheckList, CRITICAL, MUTED, S1_BLUE, S2_AQUA,
                                S3_YELLOW, S4_ROSE, apply_style, finish,
                                write_csv)
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.meshgen.fuselage import FuselageParams
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at
from pyfp3d.meshgen.wingbody import te_polyline
from pyfp3d.post.surface import planform_area, sectional_cl_from_gamma
from pyfp3d.post.surface_ls import section_cp_curve_levelset
from pyfp3d.post.unified import wall_forces
from pyfp3d.solve.newton_ls import solve_multivalued_newton
from pyfp3d.solve.picard_ls import farfield_aux_dofs, solve_multivalued_lifting
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
GATED = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"

ALPHA, M = 3.06, 0.5
FUS = FuselageParams()
LS_DIR = REPO_ROOT / "cases/meshes/onera_m6_wingbody"
B9_DIR = REPO_ROOT / "cases/demo/b9_wingbody/results"
# committed conforming reference (P14 Newton pressure-Kutta; b9_summary.csv)
CONF = {"coarse": 0.2089, "medium": 0.2173}
ETAS = (0.20, 0.44, 0.65, 0.90)

checks = CheckList("B17 far-field aux pin carries jump=gamma")


def cl_kj_exposed(z, gamma, s_ref):
    """Exposed-span cl_KJ (the same reducer B9/B16 use)."""
    o = np.argsort(z)
    zz = np.concatenate([z[o], [B_SEMI]])
    gg = np.concatenate([np.asarray(gamma)[o], [0.0]])
    return 2.0 * float(np.trapezoid(gg, zz)) / s_ref


def build(level):
    mp = LS_DIR / f"{level}.msh"
    if not mp.exists():
        return None, None
    mesh = read_mesh(str(mp))
    a = np.radians(ALPHA)
    wls = WakeLevelSet(te_polyline(FUS), direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return mesh, mvop


def _clp(mesh, mvop, phi, s_ref):
    return float(wall_forces(mesh, mvop=mvop, phi_ext=phi, alpha_deg=ALPHA,
                             s_ref=s_ref, m_inf=M, wall_tag="wall")["cl"])


def _ff_jump(mesh, mvop, phi):
    hosts, aux = farfield_aux_dofs(mesh, mvop.cm)
    if not aux.size:
        return 0.0
    side = mvop.cm.node_side[hosts]
    return float(np.max(np.abs(side * (phi[hosts] - phi[aux]))))


# --------------------------------------------------------------------------
# solves (each cached to a gitignored npz)
# --------------------------------------------------------------------------
def picard(level, mesh, mvop, aux_mode, dre):
    """A freestream Picard solve with the given far-field aux mode. Reuses the
    B9 committed legacy cache for aux_mode='legacy'."""
    if aux_mode == "legacy":
        b9 = B9_DIR / f"ls_{level}.npz"
        if b9.exists():
            return np.load(b9)["phi_ext"]
    cache = OUT / f"picard_{aux_mode}_{level}.npz"
    if cache.exists():
        return np.load(cache)["phi_ext"]
    r = solve_multivalued_lifting(mvop, mesh, M, alpha_deg=ALPHA,
                                  farfield="freestream", farfield_aux=aux_mode,
                                  n_outer_max=80, tol_residual=1e-7,
                                  direct_refactor_every=dre)
    np.savez(cache, phi_ext=r["phi_ext"])
    print(f"  [picard {aux_mode} {level}] conv={r['converged']} "
          f"n_outer={r['n_outer']} res={r['residual_norm']:.1e}", flush=True)
    return r["phi_ext"]


def newton_pin_gamma(level, mesh, mvop, seed_phi, dre):
    """Newton with farfield_aux='pin_gamma' (the new default), warm-started from
    the Picard pin_gamma state so it settles the junction quickly. Records the
    gamma trajectory + residual history."""
    cache = OUT / f"newton_pin_gamma_{level}.npz"
    if cache.exists():
        d = np.load(cache)
        return d["phi_ext"], list(d["gtraj"]), list(d["res"])
    gseed = float(np.mean(mvop.te_jump(seed_phi)))
    r = solve_multivalued_newton(mvop, mesh, M, alpha_deg=ALPHA,
                                 farfield="freestream", farfield_aux="pin_gamma",
                                 phi_init=seed_phi.copy(), gamma_init=gseed,
                                 n_seed=3, n_newton_max=40, tol_residual=1e-8,
                                 freeze_tol=1e-5, direct_refactor_every=dre)
    gtraj = [sr["gamma_mean"] for sr in r["step_records"]]
    np.savez(cache, phi_ext=r["phi_ext"], gtraj=np.asarray(gtraj),
             res=np.asarray(r["residual_history"]))
    print(f"  [newton pin_gamma {level}] conv={r['converged']} "
          f"res={r['residual_history'][-1]:.1e} "
          f"nlim={r['step_records'][-1]['n_limited']} "
          f"nflr={r['step_records'][-1]['n_floored']}", flush=True)
    return r["phi_ext"], gtraj, list(r["residual_history"])


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------
def fig_jump_collapse(level, mesh, mvop, states):
    """|far-field ring jump| for legacy / pin / pin_gamma: 53 -> 0 -> gamma."""
    hosts, aux = farfield_aux_dofs(mesh, mvop.cm)
    x = mesh.nodes[hosts, 0]
    side = mvop.cm.node_side[hosts]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    cols = {"legacy": CRITICAL, "pin (jump=0)": MUTED, "pin_gamma": S2_AQUA}
    for lab, phi in states.items():
        j = np.abs(side * (phi[hosts] - phi[aux]))
        ax.semilogy(x, np.maximum(j, 1e-16), "o", ms=6, label=lab,
                    color=cols.get(lab, S1_BLUE))
    ax.set_xlabel("x (far-field host node)"); ax.set_ylabel("|ring jump| [phi]")
    ax.set_title(f"B17 far-field ring jump, {level}: legacy garbage -> pin 0 -> "
                 "pin_gamma=gamma")
    ax.legend(fontsize=8)
    finish(fig, OUT, f"b17_jump_collapse_{level}.png")


def fig_triangle(states_clp):
    """The 4-way cl_p triangle at coarse+medium: pin(jump0) non-monotone down,
    pin_gamma monotone to conforming."""
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    xs = [0, 1]
    order = [("conforming", S1_BLUE, "o-"), ("legacy", S4_ROSE, "s--"),
             ("pin (jump=0)", MUTED, "^--"), ("pin_gamma", S2_AQUA, "D-")]
    for lab, col, sty in order:
        ys = [states_clp["coarse"].get(lab), states_clp["medium"].get(lab)]
        if all(y is not None for y in ys):
            ax.plot(xs, ys, sty, color=col, lw=1.6, ms=7, label=lab)
    ax.set_xticks(xs); ax.set_xticklabels(["coarse", "medium"])
    ax.set_ylabel("cl_p (wing)")
    ax.set_title("B17 cl_p triangle: pin jump=0 kills lift; pin_gamma is monotone")
    ax.legend(fontsize=8)
    finish(fig, OUT, "b17_triangle.png")


def fig_spanwise(level, mesh, mvop, phi_pin0, phi_ping, phi_leg):
    """Gamma(z) and sectional cl(z): the pin jump=0 uniform offset vs pin_gamma
    (restored, tracks legacy/conforming)."""
    z = mesh.nodes[mvop.cm.te_nodes, 2]
    o = np.argsort(z)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    chord = np.array([chord_at(zz) for zz in z[o]])
    for phi, lab, col in ((phi_leg, "legacy (~conforming)", S4_ROSE),
                          (phi_ping, "pin_gamma (fix)", S2_AQUA),
                          (phi_pin0, "pin jump=0 (B16)", MUTED)):
        g = mvop.te_jump(phi)[o]
        ax[0].plot(z[o], g, "-", lw=1.4, color=col, label=lab)
        ax[1].plot(z[o], sectional_cl_from_gamma(g, chord=chord), "-", lw=1.4,
                   color=col)
    ax[0].set_xlabel("z"); ax[0].set_ylabel("Gamma(z)"); ax[0].legend(fontsize=8)
    ax[1].set_xlabel("z"); ax[1].set_ylabel("sectional cl = 2*Gamma/(u*c)")
    fig.suptitle(f"B17 spanwise, {level}: pin jump=0 uniform offset vs pin_gamma")
    finish(fig, OUT, f"b17_spanwise_{level}.png")


# --------------------------------------------------------------------------
# post-processing audit (GB17.2): cl_p and cl_KJ move together (the pin gap is a
# flow-state change, not a post artifact), and the plotted sectional cl is
# Gamma-based.
# --------------------------------------------------------------------------
def postproc_audit(level, mesh, mvop, phi_pin0, phi_ref, s_ref):
    """Compare the two INDEPENDENT lift reductions on the pin-jump0 vs a
    reference (legacy/pin_gamma) state: cl_p (surface pressure integral) and
    cl_KJ (circulation integral). If they move by ~the same %, the state genuinely
    differs (not a post artifact)."""
    z = mesh.nodes[mvop.cm.te_nodes, 2]
    clp0 = _clp(mesh, mvop, phi_pin0, s_ref)
    clpR = _clp(mesh, mvop, phi_ref, s_ref)
    kj0 = cl_kj_exposed(z, mvop.te_jump(phi_pin0), s_ref)
    kjR = cl_kj_exposed(z, mvop.te_jump(phi_ref), s_ref)
    dclp = abs(clpR - clp0) / abs(clpR) * 100
    dkj = abs(kjR - kj0) / abs(kjR) * 100
    write_csv(OUT, f"postproc_audit_{level}.csv",
              "reduction,pin_jump0,reference,gap_pct",
              [("cl_p_surface_integral", f"{clp0:.4f}", f"{clpR:.4f}", f"{dclp:.1f}"),
               ("cl_kj_circulation", f"{kj0:.4f}", f"{kjR:.4f}", f"{dkj:.1f}")])
    return dclp, dkj


# --------------------------------------------------------------------------
def run_coarse():
    level = "coarse"
    mesh, mvop = build(level)
    if mesh is None:
        print(f"[skip] {level}: wing-body mesh absent (regenerate)")
        return False
    print(f"=== {level} (Parts 1-2) ===", flush=True)
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    phi_leg = picard(level, mesh, mvop, "legacy", dre=1)
    phi_pin0 = picard(level, mesh, mvop, "pin", dre=1)
    phi_ping = picard(level, mesh, mvop, "pin_gamma", dre=1)

    clp = {"conforming": CONF[level],
           "legacy": _clp(mesh, mvop, phi_leg, s_ref),
           "pin (jump=0)": _clp(mesh, mvop, phi_pin0, s_ref),
           "pin_gamma": _clp(mesh, mvop, phi_ping, s_ref)}
    jl, jp, jg = (_ff_jump(mesh, mvop, p) for p in (phi_leg, phi_pin0, phi_ping))
    print(f"  cl_p: conf {clp['conforming']:.4f} | legacy {clp['legacy']:.4f} "
          f"(ff|jump| {jl:.1f}) | pin0 {clp['pin (jump=0)']:.4f} (ff {jp:.1e}) | "
          f"pin_gamma {clp['pin_gamma']:.4f} (ff {jg:.3f})", flush=True)
    write_csv(OUT, f"triangle_{level}.csv", "path,cl_p_wing,ff_max_jump",
              [("conforming", f"{clp['conforming']:.4f}", ""),
               ("legacy", f"{clp['legacy']:.4f}", f"{jl:.3f}"),
               ("pin_jump0", f"{clp['pin (jump=0)']:.4f}", f"{jp:.2e}"),
               ("pin_gamma", f"{clp['pin_gamma']:.4f}", f"{jg:.4f}")])
    fig_jump_collapse(level, mesh, mvop,
                      {"legacy": phi_leg, "pin (jump=0)": phi_pin0,
                       "pin_gamma": phi_ping})

    d_ping = abs(clp["pin_gamma"] - CONF[level]) / CONF[level] * 100
    checks.add("GB17.1", f"{level}_pin_gamma_vs_conforming",
               f"pin_gamma {clp['pin_gamma']:.4f} vs conforming {CONF[level]} "
               f"({d_ping:.1f}%)",
               "pin_gamma matches conforming at coarse (<1%); the pin ring jump "
               f"collapses legacy {jl:.0f} -> pin_gamma {jg:.3f}=gamma", d_ping < 1.0)
    checks.add("GB17.1", f"{level}_legacy_contaminated",
               f"legacy {clp['legacy']:.4f} (ff|jump| {jl:.0f}) vs pin_gamma "
               f"{clp['pin_gamma']:.4f}",
               "the coarse legacy far-field aux carry O(50) garbage -> a 12% "
               "lift deficit (B9's 'coarse resolution' was CONTAMINATION)",
               jl > 10.0)

    # GB17.2 post-processing audit (pin jump=0 vs pin_gamma)
    dclp, dkj = postproc_audit(level, mesh, mvop, phi_pin0, phi_ping, s_ref)
    checks.add("GB17.2", f"{level}_clp_kj_move_together",
               f"cl_p gap {dclp:.1f}% vs cl_KJ gap {dkj:.1f}%",
               "the surface-pressure and circulation reductions move together "
               "=> the pin gap is a flow-state change, NOT a post artifact",
               abs(dclp - dkj) < 3.0)
    return clp


def run_medium(clp_coarse):
    level = "medium"
    mesh, mvop = build(level)
    if mesh is None:
        print(f"[skip] {level}: wing-body mesh absent (regenerate)")
        return
    print(f"=== {level} (Part 3, gated) ===", flush=True)
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    phi_leg = picard(level, mesh, mvop, "legacy", dre=1000)
    phi_pin0 = picard(level, mesh, mvop, "pin", dre=1000)
    phi_ping = picard(level, mesh, mvop, "pin_gamma", dre=1000)
    phi_nwt, gtraj, res = newton_pin_gamma(level, mesh, mvop, phi_ping, dre=1000)

    clp = {"conforming": CONF[level],
           "legacy": _clp(mesh, mvop, phi_leg, s_ref),
           "pin (jump=0)": _clp(mesh, mvop, phi_pin0, s_ref),
           "pin_gamma": _clp(mesh, mvop, phi_ping, s_ref)}
    clp_nwt = _clp(mesh, mvop, phi_nwt, s_ref)
    z = mesh.nodes[mvop.cm.te_nodes, 2]
    print(f"  cl_p: conf {clp['conforming']:.4f} | legacy {clp['legacy']:.4f} | "
          f"pin0 {clp['pin (jump=0)']:.4f} | pin_gamma(Picard) {clp['pin_gamma']:.4f} "
          f"| pin_gamma(Newton) {clp_nwt:.4f}", flush=True)
    write_csv(OUT, f"triangle_{level}.csv",
              "path,cl_p_wing,cl_kj_exposed,solver",
              [("conforming", f"{clp['conforming']:.4f}", "", "Newton(P14)"),
               ("legacy", f"{clp['legacy']:.4f}",
                f"{cl_kj_exposed(z, mvop.te_jump(phi_leg), s_ref):.4f}", "Picard"),
               ("pin_jump0", f"{clp['pin (jump=0)']:.4f}",
                f"{cl_kj_exposed(z, mvop.te_jump(phi_pin0), s_ref):.4f}", "Picard"),
               ("pin_gamma", f"{clp['pin_gamma']:.4f}",
                f"{cl_kj_exposed(z, mvop.te_jump(phi_ping), s_ref):.4f}", "Picard"),
               ("pin_gamma", f"{clp_nwt:.4f}",
                f"{cl_kj_exposed(z, mvop.te_jump(phi_nwt), s_ref):.4f}", "Newton")])
    clp_all = {"coarse": clp_coarse,
               "medium": {**clp, "pin_gamma": clp["pin_gamma"]}}
    fig_triangle(clp_all)
    fig_spanwise(level, mesh, mvop, phi_pin0, phi_ping, phi_leg)
    write_csv(OUT, "newton_pin_gamma_trajectory_medium.csv",
              "step,gamma_mean,residual",
              [(i, f"{g:.6f}", f"{r:.3e}") for i, (g, r) in
               enumerate(zip(gtraj, res))])

    # GB17.3: pin jump=0 lands 0.169 on BOTH solvers (the BC error, not a stall).
    # The Newton-pin-jump0 number is B16's committed 0.1690; here we show the
    # Picard-pin-jump0 independently lands there too.
    d_pin0_solvers = abs(clp["pin (jump=0)"] - 0.1690) / 0.1690 * 100
    checks.add("GB17.3", "pin_jump0_both_solvers",
               f"Picard pin0 {clp['pin (jump=0)']:.4f} vs B16 Newton pin0 0.1690 "
               f"({d_pin0_solvers:.1f}%)",
               "pin jump=0 lands ~0.169 on BOTH Picard and Newton => GB16.4 is a "
               "BC-modelling error (killed outflow circulation), NOT a Newton stall",
               d_pin0_solvers < 2.0)

    # GB17.4: pin_gamma closes the triangle (Newton ~ Picard ~ conforming, and
    # monotone coarse->medium).
    d_solvers = abs(clp_nwt - clp["pin_gamma"]) / clp["pin_gamma"] * 100
    d_conf = abs(clp["pin_gamma"] - CONF[level]) / CONF[level] * 100
    monotone = clp_coarse["pin_gamma"] < clp["pin_gamma"] < CONF[level]
    checks.add("GB17.4", "pin_gamma_closes_triangle",
               f"Newton {clp_nwt:.4f} ~ Picard {clp['pin_gamma']:.4f} "
               f"({d_solvers:.1f}%); vs conforming {d_conf:.1f}%; monotone={monotone}",
               "pin_gamma: both solvers agree (<1%) and cl_p is monotone "
               "coarse->medium toward conforming (residual gap is truncation, "
               "not the jump=0 error)",
               d_solvers < 1.0 and monotone)

    # GB17.5: the spanwise Gamma(z) uniform offset (jump=0) is removed by pin_gamma
    g0 = mvop.te_jump(phi_pin0)
    gg = mvop.te_jump(phi_ping)
    off = float(np.mean((gg - g0) / np.maximum(np.abs(gg), 1e-9))) * 100
    checks.add("GB17.5", "spanwise_offset_removed",
               f"pin jump=0 Gamma(z) is ~{off:.0f}% below pin_gamma (uniform)",
               "RECORDED: the pin jump=0 spanwise Gamma is a uniform multiplicative "
               "deficit (the flow lost circulation everywhere); pin_gamma restores it",
               True)

    # GB17.2 medium postproc audit
    dclp, dkj = postproc_audit(level, mesh, mvop, phi_pin0, phi_ping, s_ref)
    checks.add("GB17.2", "medium_clp_kj_move_together",
               f"cl_p gap {dclp:.1f}% vs cl_KJ gap {dkj:.1f}%",
               "medium: surface-pressure and circulation reductions move together "
               "(~22% both) => a real flow-state change, not a post artifact",
               abs(dclp - dkj) < 3.0)

    # GB17.6 RECORDED: does farfield="vortex" (the physically-consistent lifting
    # far field, freestream + PG vortex on the MAIN dofs) close the pin_gamma
    # residual gap to conforming? Measured: NO -- it BRACKETS conforming from the
    # other side (medium +2.5% vs pin_gamma's -2.6%), and its free far-field aux
    # churn at coarse (needs its own pin). The 2-3% is far-field truncation, not
    # a bug -- freestream pin_gamma stays the recommended (both-resolutions-clean)
    # BC. (User-requested B17 side-evaluation.)
    cache = OUT / f"picard_vortex_{level}.npz"
    if cache.exists():
        phi_vtx = np.load(cache)["phi_ext"]
    else:
        rv = solve_multivalued_lifting(mvop, mesh, M, alpha_deg=ALPHA,
                                       farfield="vortex", farfield_aux="legacy",
                                       n_outer_max=80, tol_residual=1e-7,
                                       direct_refactor_every=1000)
        phi_vtx = rv["phi_ext"]
        np.savez(cache, phi_ext=phi_vtx)
    clp_vtx = _clp(mesh, mvop, phi_vtx, s_ref)
    d_vtx = (clp_vtx - CONF[level]) / CONF[level] * 100
    d_ping_signed = (clp["pin_gamma"] - CONF[level]) / CONF[level] * 100
    print(f"  [vortex {level}] cl_p {clp_vtx:.4f} vs conforming {CONF[level]} "
          f"({d_vtx:+.1f}%); pin_gamma {clp['pin_gamma']:.4f} ({d_ping_signed:+.1f}%)",
          flush=True)
    write_csv(OUT, "farfield_bc_bracket_medium.csv",
              "far_field,cl_p_wing,vs_conforming_pct,note",
              [("conforming", f"{CONF[level]:.4f}", "0.0", "reference (P14 Newton)"),
               ("freestream_pin_gamma", f"{clp['pin_gamma']:.4f}",
                f"{d_ping_signed:+.1f}", "recommended; clean at both resolutions"),
               ("vortex_legacy", f"{clp_vtx:.4f}", f"{d_vtx:+.1f}",
                "brackets from above; churns at coarse (free aux)")])
    checks.add("GB17.6", "vortex_brackets_conforming",
               f"vortex {clp_vtx:.4f} ({d_vtx:+.1f}%) vs pin_gamma "
               f"{clp['pin_gamma']:.4f} ({d_ping_signed:+.1f}%)",
               "RECORDED: vortex does NOT close the gap -- it brackets conforming "
               "from the OTHER side; the 2-3% is far-field truncation. freestream "
               "pin_gamma stays recommended", True)


# --------------------------------------------------------------------------
def main():
    apply_style()
    clp_coarse = run_coarse()
    if not clp_coarse:
        print("no levels ran (mesh absent)")
        sys.exit(0)
    if GATED:
        run_medium(clp_coarse)
    else:
        print("\n[medium gated: set PYFP3D_TRANSONIC_GATES=1 for Part 3]")
    sys.exit(checks.report(OUT, "checks.csv"))


if __name__ == "__main__":
    main()
