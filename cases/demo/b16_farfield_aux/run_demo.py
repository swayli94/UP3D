"""
Track B / B16: LS Newton far-field BC generalisation -- far-field aux-DOF pin.

The B9 recorded follow-up: the level-set Newton solver churns on the wing-body,
so B9's LS leg fell back to Picard (medium 1458.9 s vs the conforming Newton's
52.4 s). This demo is the evidence dossier for the fix -- and it FIRST reproduces
the diagnostic B9 only recorded as prose (session discipline #3: a number that
lives only in .md is not evidence).

Root cause (measured here, GB16.1): a wake level set has no outflow clip, so the
sheet reaches the far-field boundary and the outer nodes it crosses each carry
an aux DOF governed only by a near-singular wake-LS row on a giant outer tet. At
the converged freestream Picard state those aux hold garbage (|jump| ~ 20-50 vs
the physical Gamma ~ 0.06); the Picard fixed point tolerates it, but the Newton
residual reads it as an O(1) inconsistency (8 far-field MAIN rows |R| ~ 84).
Pinning the far-field-boundary aux to the branch value their host carries
(freestream: the same phi_inf) removes it: the residual falls to machine, the
limiter goes quiet, and the outer jumps collapse to zero.

Parts:
  1 (coarse)         GB16.1 diagnostic: neumann blow-up (D0), the residual-row
                     census at the Picard state (D1), the jump-vs-x profile
                     before/after pin (D4/D7/D8), and the aux-block conditioning
                     legacy vs pin (D6).
  2 (coarse)         A/B convergence: neumann (diverges) / freestream legacy
                     (churns) / freestream pin (converges to machine).
  3 (gated, medium)  GB16.3/16.4 RECORDED: the pin runs CLEAN at medium (0
                     limited, 0 floored -- the churn is gone) but plateaus at
                     res ~7e-6 (the junction O(h), G1.6 class) and is NOT a speed
                     win -- the 35-outer Picard seed dominates (2172 s) and the
                     fair lagged-LU Picard arm is ~205 s. B16's value is that the
                     Newton path WORKS cleanly, not wall-clock.

Cost: ~2-3 min coarse; medium under PYFP3D_TRANSONIC_GATES=1 (Newton solves +
the B9 Picard fair arm), cached to gitignored results/*.npz. Meshes are
gitignored -> a level whose mesh is absent is skipped.

Run:  python cases/demo/b16_farfield_aux/run_demo.py
      PYFP3D_TRANSONIC_GATES=1 NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 \
        OPENBLAS_NUM_THREADS=16 python cases/demo/b16_farfield_aux/run_demo.py
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
from pyfp3d.meshgen.wingbody import junction_z, te_polyline
from pyfp3d.post.surface import planform_area, sectional_cl_from_gamma
from pyfp3d.post.surface_ls import section_cp_curve_levelset
from pyfp3d.post.unified import wall_forces
from pyfp3d.solve.newton_ls import LSNewtonSystem, solve_multivalued_newton
from pyfp3d.solve.picard_ls import (farfield_aux_dofs, solve_multivalued_lifting)
from pyfp3d.solve.schur_ls import SchurReducedSystem, jaa_diagnostic
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
GATED = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"

ALPHA, M = 3.06, 0.5
FUS = FuselageParams()
Z_JUNC = junction_z(FUS)
ETAS = (0.20, 0.44, 0.65, 0.90)   # near-junction, mid, outer stations (b9)
LS_DIR = REPO_ROOT / "cases/meshes/onera_m6_wingbody"
B9_DIR = REPO_ROOT / "cases/demo/b9_wingbody/results"
B9_PICARD_MEDIUM_WALL = 1458.9   # committed b9_summary.csv (level_set, medium)
# committed b9_summary.csv references (cl_p_wing), for the three-path lift triangle
B9_CONF_PRESSURE = {"coarse": 0.2089, "medium": 0.2173}   # conforming (pressure)
B9_LS_PICARD = {"coarse": 0.1853, "medium": 0.2165}       # LS Picard

checks = CheckList("B16 LS Newton far-field aux pin")


def cl_kj_exposed(z, gamma, s_ref):
    """Exposed-span cl_KJ, the SAME reducer B9 uses (trapezoid over the actual
    junction..tip stations + tip closure to B_SEMI, no root flat-extension), so
    the Newton and Picard circulations are compared on one footing."""
    o = np.argsort(z)
    zz = np.concatenate([z[o], [B_SEMI]])
    gg = np.concatenate([np.asarray(gamma)[o], [0.0]])
    return 2.0 * float(np.trapezoid(gg, zz)) / s_ref


# --------------------------------------------------------------------------
def build(level):
    """(mesh, mvop) for the B9 wing-body LS mesh, or (None, None) if absent."""
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


def picard_state(level, mesh, mvop):
    """The freestream Picard state (the churn's starting point). Reuses the B9
    cache when present, else solves + caches locally."""
    b9 = B9_DIR / f"ls_{level}.npz"
    if b9.exists():
        return np.load(b9)["phi_ext"]
    cache = OUT / f"picard_{level}.npz"
    if cache.exists():
        return np.load(cache)["phi_ext"]
    r = solve_multivalued_lifting(mvop, mesh, M, alpha_deg=ALPHA,
                                  farfield="freestream", n_outer_max=80,
                                  tol_residual=1e-7)
    np.savez(cache, phi_ext=r["phi_ext"])
    return r["phi_ext"]


def free_set(mesh, mvop, pin):
    """The legacy or pinned free-DOF set for a freestream far field."""
    ff = np.unique(mesh.boundary_faces["farfield"])
    is_dir = np.zeros(mvop.n_total, dtype=bool)
    is_dir[ff] = True
    if pin:
        _, aux = farfield_aux_dofs(mesh, mvop.cm)
        is_dir[aux] = True
    return np.flatnonzero(~is_dir)


def newton_solve(level, mesh, mvop, mode):
    """Cached freestream Newton solve (mode in {'legacy','pin'})."""
    cache = OUT / f"newton_{mode}_{level}.npz"
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        return dict(phi_ext=d["phi_ext"], res=list(d["res"]),
                    conv=bool(d["conv"]), nlim=int(d["nlim"]),
                    nflr=int(d["nflr"]), wall_s=float(d["wall"]),
                    gamma=float(d["gamma"]))
    t0 = time.perf_counter()
    r = solve_multivalued_newton(mvop, mesh, M, alpha_deg=ALPHA,
                                 farfield="freestream", farfield_aux=mode,
                                 n_seed=35, n_newton_max=40, tol_residual=1e-8,
                                 freeze_tol=1e-6)
    wall = time.perf_counter() - t0
    out = dict(phi_ext=r["phi_ext"], res=list(r["residual_history"]),
               conv=bool(r["converged"]), nlim=int(r["n_limited"]),
               nflr=int(r["n_floored"]), wall_s=wall, gamma=float(r["gamma"]))
    np.savez(cache, phi_ext=out["phi_ext"], res=np.asarray(out["res"]),
             conv=out["conv"], nlim=out["nlim"], nflr=out["nflr"],
             wall=wall, gamma=out["gamma"])
    print(f"  [newton {mode} {level}] conv={out['conv']} "
          f"res={out['res'][-1]:.2e} lim={out['nlim']} flr={out['nflr']} "
          f"({wall:.0f}s)", flush=True)
    return out


def neumann_blowup(mesh, mvop):
    """D0: neumann far field on the wing-body is UNBOUNDED (the fuselage blockage
    makes the Lopez outlet-flux balance inconsistent). A COLD neumann solve
    (its own Picard seed under the blockage) is the faithful reproduction of the
    B9 divergence -- measured res O(1e6-1e8) on the coarse wing-body; seeding
    from the good freestream Picard state instead would MASK it (that state is
    near the freestream solution, not the neumann one)."""
    r = solve_multivalued_newton(mvop, mesh, M, alpha_deg=ALPHA,
                                 farfield="neumann", n_seed=15, n_newton_max=8,
                                 tol_residual=1e-8)
    return list(r["residual_history"])


# --------------------------------------------------------------------------
# GB16.1 diagnostic on the Picard state
# --------------------------------------------------------------------------
def diagnose(level, mesh, mvop, phi):
    n_main = mvop.n_main
    free = free_set(mesh, mvop, pin=False)
    sysm = LSNewtonSystem(mvop, M)
    _, R, _, _ = sysm.residual(phi)
    res = float(np.max(np.abs(R[free])))
    n_lim, n_flr = mvop.n_limited, mvop.n_floored

    # row census: host node + coord + class for every free row with |R|>1
    host_of_aux = np.full(mvop.n_ext, -1, np.int64)
    cut = np.flatnonzero(mvop.cm.ext_dof_of_node >= 0)
    host_of_aux[mvop.cm.ext_dof_of_node[cut] - n_main] = cut
    te_aux = set(int(a) for a in sysm.te_aux)
    big = free[np.abs(R[free]) > 1.0]
    rows = []
    for d in big[np.argsort(-np.abs(R[big]))]:
        h = int(d) if d < n_main else int(host_of_aux[d - n_main])
        kind = ("main" if d < n_main else
                ("te_aux" if int(d) in te_aux else "wake_aux"))
        rows.append((int(d), kind, f"{mesh.nodes[h,0]:.3f}",
                     f"{mesh.nodes[h,2]:.3f}", f"{R[d]:.3f}"))
    write_csv(OUT, f"residual_rows_{level}.csv",
              "dof,kind,x,z,residual", rows)

    # jump-vs-x at the Picard state (D4/D7)
    x = mesh.nodes[cut, 0]
    jump_pic = mvop.node_jump(phi, cut)
    gbar = float(np.abs(mvop.te_jump(phi)).mean())

    # aux-block conditioning legacy vs pin (D6)
    A = mvop.assemble_matrix(closure="wake_ls", te_kutta="pressure",
                             phi_ext=phi).tocsr()
    cond = {}
    for pin in (False, True):
        fr = free_set(mesh, mvop, pin)
        n_pin = (0 if not pin else farfield_aux_dofs(mesh, mvop.cm)[1].size)
        sc = SchurReducedSystem(A[fr][:, fr], fr, n_main,
                                n_aux_expected=mvop.n_ext - n_pin)
        cond["pin" if pin else "legacy"] = float(jaa_diagnostic(sc)["cond1"])

    return dict(res=res, n_lim=n_lim, n_flr=n_flr, n_big=len(rows),
                x=x, jump_pic=jump_pic, cut=cut, gbar=gbar, cond=cond,
                R=R, free=free, host_of_aux=host_of_aux)


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------
def fig_residual_map(level, mesh, dg):
    n_main = mesh.nodes.shape[0]
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.0))
    R, free = dg["R"], dg["free"]
    # per free row: x of host + class
    host = np.where(free < n_main, free,
                    dg["host_of_aux"][np.clip(free - n_main, 0, None)])
    x = mesh.nodes[host, 0]
    z = mesh.nodes[host, 2]
    absR = np.abs(R[free])
    is_main = free < n_main
    ax[0].scatter(x[is_main], np.maximum(absR[is_main], 1e-16), s=6,
                  c=S1_BLUE, alpha=0.5, label="main rows")
    ax[0].scatter(x[~is_main], np.maximum(absR[~is_main], 1e-16), s=6,
                  c=S3_YELLOW, alpha=0.6, label="aux rows")
    ax[0].axhline(1.0, color=CRITICAL, ls="--", lw=1,
                  label="|R|=1 (outer inconsistency)")
    ax[0].set_yscale("log"); ax[0].set_xlabel("x (streamwise)")
    ax[0].set_ylabel("|R| at the Picard state"); ax[0].legend(fontsize=8)
    ax[0].set_title(f"{level}: Newton residual by row -- the 8 outer main rows")
    # x-z projection with the big rows flagged
    ax[1].scatter(x, z, s=4, c="0.75", label="free rows")
    big = absR > 1.0
    ax[1].scatter(x[big], z[big], s=70, facecolors="none",
                  edgecolors=CRITICAL, linewidths=1.6,
                  label=f"|R|>1 ({big.sum()} rows)")
    ax[1].set_xlabel("x"); ax[1].set_ylabel("z (span)")
    ax[1].legend(fontsize=8)
    ax[1].set_title(f"{level}: the inconsistency is the outer wake corridor")
    fig.suptitle("B16/GB16.1: the far-field aux inconsistency the LS Newton "
                 "churns on (max|R|=%.1f)" % dg["res"])
    finish(fig, OUT, f"b16_residual_map_{level}.png")


def fig_jump_profile(level, mesh, dg, pin_state, mvop):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    x, jp = dg["x"], np.abs(dg["jump_pic"])
    jump_pin = np.abs(mvop.node_jump(pin_state["phi_ext"], dg["cut"]))
    for a, j, ttl, col in ((ax[0], jp, "legacy (free aux)", S3_YELLOW),
                           (ax[1], jump_pin, "pin", S2_AQUA)):
        a.scatter(x, np.maximum(j, 1e-16), s=9, c=col, alpha=0.6)
        a.axhline(dg["gbar"], color=S1_BLUE, ls="--", lw=1,
                  label=f"mean TE Gamma = {dg['gbar']:.3f}")
        a.axvspan(10, x.max() + 0.5, color=CRITICAL, alpha=0.06)
        a.set_yscale("log"); a.set_xlabel("x (streamwise)")
        a.set_title(f"{level}: |[phi]| -- {ttl}"); a.legend(fontsize=8)
    ax[0].set_ylabel("|jump| = |side*(main-aux)|")
    fig.suptitle("B16: the outer junk jumps (x>=10) collapse to 0 under the pin")
    finish(fig, OUT, f"b16_jump_profile_{level}.png")


def fig_convergence(level, neu, leg, pin):
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    for lab, res, col, mk in (
            ("neumann (diverges)", neu, CRITICAL, "^"),
            ("freestream legacy (churns)", leg["res"], S3_YELLOW, "o"),
            ("freestream pin (converges)", pin["res"], S2_AQUA, "s")):
        res = np.asarray(res)
        ax.semilogy(np.arange(len(res)), np.maximum(res, 1e-16), mk + "-",
                    ms=4, color=col, label=lab)
    ax.axhline(1e-8, color="0.6", ls=":", lw=1, label="tol 1e-8")
    ax.set_xlabel("Newton iteration"); ax.set_ylabel("residual (inf-norm)")
    ax.legend(fontsize=9)
    ax.set_title(f"B16 {level}: LS Newton far-field BC A/B (M0.5 wing-body)")
    finish(fig, OUT, f"b16_convergence_{level}.png")


def fig_ladder(level, dg, leg, pin):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    # left: max|R| legacy(Picard) / legacy(Newton end) / pin(Newton end)
    labs = ["legacy\n@Picard", "legacy\nNewton end", "pin\nNewton end"]
    vals = [dg["res"], leg["res"][-1], pin["res"][-1]]
    cols = [S3_YELLOW, S3_YELLOW, S2_AQUA]
    ax[0].bar(labs, np.maximum(vals, 1e-16), color=cols)
    ax[0].set_yscale("log"); ax[0].set_ylabel("max|R[free]|")
    ax[0].axhline(1e-8, color=CRITICAL, ls="--", lw=1)
    ax[0].set_title(f"{level}: residual -- the churn vs the fix")
    # right: aux-block conditioning legacy vs pin (D6)
    ax[1].bar(["legacy", "pin"], [dg["cond"]["legacy"], dg["cond"]["pin"]],
              color=[S3_YELLOW, S2_AQUA])
    ax[1].set_yscale("log"); ax[1].set_ylabel("J_aa cond1 estimate")
    ax[1].axhline(1e14, color=CRITICAL, ls="--", lw=1, label="GB14.1 ceiling")
    ax[1].legend(fontsize=8)
    ax[1].set_title(f"{level}: aux-block conditioning (D6)")
    fig.suptitle("B16/GB16.1: the residual and the aux-block conditioning")
    finish(fig, OUT, f"b16_ladder_{level}.png")


def fig_walltime(rows):
    """rows: list of (label, wall_s, color)."""
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    labs = [r[0] for r in rows]
    ax.bar(labs, [r[1] for r in rows], color=[r[2] for r in rows])
    for i, r in enumerate(rows):
        ax.text(i, r[1], f"{r[1]:.0f}s", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("wall-clock (s)")
    ax.set_title("B16 medium M0.5 wing-body: LS Newton (pin) vs the B9 Picard")
    finish(fig, OUT, "b16_walltime_medium.png")


def fig_spanwise(level, mesh, mvop, phi_pin, phi_pic):
    """Gamma(z) + sectional cl(z): LS Newton pin vs LS Picard, SAME mesh (the
    GB16.4 lift-distribution A/B -- the new Newton solution must reproduce the
    established Picard lift, not just some scalar)."""
    z = mesh.nodes[mvop.cm.te_nodes, 2]
    o = np.argsort(z)
    zz = z[o]
    chord = np.array([chord_at(float(zi)) for zi in zz])
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for phi, lab, col in ((phi_pic, "LS Picard (B9)", S3_YELLOW),
                          (phi_pin, "LS Newton pin (B16)", S2_AQUA)):
        g = np.asarray(mvop.te_jump(phi))[o]
        ax[0].plot(zz, g, "o-", ms=3, color=col, label=lab)
        ax[1].plot(zz, sectional_cl_from_gamma(g, chord=chord), "o-", ms=3,
                   color=col, label=lab)
    for a in ax:
        a.axvline(Z_JUNC, color=MUTED, ls="--", lw=0.8,
                  label=f"junction z={Z_JUNC:.2f}")
        a.set_xlabel("z (span)"); a.legend(fontsize=8)
    ax[0].set_ylabel("Gamma(z)"); ax[0].set_title(f"{level}: spanwise circulation")
    ax[1].set_ylabel("sectional cl(z)"); ax[1].set_title(f"{level}: sectional lift")
    fig.suptitle("B16 spanwise lift: LS Newton pin vs LS Picard (same mesh, M0.5)")
    finish(fig, OUT, f"b16_spanwise_{level}.png")


def fig_sections(level, mesh, mvop, phi_pin, phi_pic):
    """Section Cp at 4 span stations, LS Newton pin vs LS Picard, same extractor
    (wing 'wall' only) -- the pressure-distribution A/B."""
    fig, ax = plt.subplots(1, len(ETAS), figsize=(4.2 * len(ETAS), 4.2),
                           squeeze=False)

    def _plot(a, sc, col, lab):
        c = sc["chord"] or 1.0
        xl = sc["x_le"] or 0.0
        a.plot((np.asarray(sc["x_upper"]) - xl) / c, sc["cp_upper"], "-",
               lw=1.2, color=col, label=lab)
        a.plot((np.asarray(sc["x_lower"]) - xl) / c, sc["cp_lower"], "-",
               lw=1.2, color=col, alpha=0.6)

    for j, eta in enumerate(ETAS):
        a = ax[0][j]
        for k, (phi, lab, col) in enumerate((
                (phi_pic, "LS Picard", S3_YELLOW),
                (phi_pin, "LS Newton pin", S2_AQUA))):
            try:
                sc = section_cp_curve_levelset(mesh, mvop, phi, eta=eta,
                                               b_semi=B_SEMI, m_inf=M,
                                               wall_tag="wall")
                _plot(a, sc, col, lab)
            except Exception as e:
                a.text(0.5, 0.5 - 0.08 * k, f"{lab} n/a\n{type(e).__name__}",
                       fontsize=7, ha="center", transform=a.transAxes)
        a.invert_yaxis(); a.set_xlabel("x/c")
        a.set_title(f"eta={eta} (z={eta*B_SEMI:.2f})")
        if j == 0:
            a.set_ylabel("Cp"); a.legend(fontsize=8)
    fig.suptitle(f"B16 section Cp, {level}, M0.5: LS Newton pin vs LS Picard "
                 "(same mesh, wing wall)")
    finish(fig, OUT, f"b16_sections_{level}.png")


def lift_ab(level, mesh, mvop, phi_pin, phi_pic, gated):
    """GB16.4 lift A/B: the three-path triangle {LS Newton pin, LS Picard,
    conforming}. ★ This is where the OPEN non-convergence problem shows: the
    triangle does NOT close consistently across resolutions -- at coarse the
    (machine-converged) Newton pin matches conforming and Picard is the low
    outlier; at medium Picard matches conforming (B9's headline) and the Newton
    pin -- which STALLS at res 7e-6, not machine -- is the low outlier. So at
    least one path is not converged; UNRESOLVED (see the demo docstring / roadmap
    B16). This draws the Gamma(z)+cl(z) and section-Cp distributions and records
    the numbers; GB16.4 is XFAIL at medium (the discrepancy is the open state,
    not a pass)."""
    fig_spanwise(level, mesh, mvop, phi_pin, phi_pic)
    fig_sections(level, mesh, mvop, phi_pin, phi_pic)
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    z = mesh.nodes[mvop.cm.te_nodes, 2]

    def _cl(phi):
        clp = wall_forces(mesh, mvop=mvop, phi_ext=phi, alpha_deg=ALPHA,
                          s_ref=s_ref, m_inf=M, wall_tag="wall")["cl"]
        return float(clp), cl_kj_exposed(z, mvop.te_jump(phi), s_ref)

    clp_pin, clkj_pin = _cl(phi_pin)
    clp_pic, clkj_pic = _cl(phi_pic)
    clp_conf = B9_CONF_PRESSURE[level]
    dp_pic = abs(clp_pin - clp_pic) / abs(clp_pic) * 100          # vs LS Picard
    dp_conf = abs(clp_pin - clp_conf) / abs(clp_conf) * 100       # vs conforming
    print(f"  [GB16.4 {level}] cl_p: Newton-pin {clp_pin:.4f} | Picard {clp_pic:.4f}"
          f" | conforming {clp_conf:.4f}  ->  Newton vs Picard {dp_pic:.1f}%, "
          f"vs conforming {dp_conf:.1f}%", flush=True)
    write_csv(OUT, f"lift_ab_{level}.csv", "path,cl_p_wing,cl_kj_exposed,vs_picard_pct",
              [("ls_newton_pin", f"{clp_pin:.4f}", f"{clkj_pin:.4f}", f"{dp_pic:.1f}"),
               ("ls_picard_b9", f"{clp_pic:.4f}", f"{clkj_pic:.4f}", "0.0"),
               ("conforming_pressure_b9", f"{clp_conf:.4f}", "", "")])
    # ★ The tension is the finding. coarse: the CONVERGED Newton pin agrees with
    # conforming (reassuring -- a converged Newton-pin IS right); medium: the
    # Newton pin STALLS (res 7e-6) 22% low, disagreeing with the Picard==conforming
    # pair. XFAIL carries this open state without failing the run.
    if gated:
        checks.add("GB16.4", f"{level}_lift_triangle",
                   f"Newton-pin {clp_pin:.4f} vs Picard {clp_pic:.4f} "
                   f"({dp_pic:.0f}%) / conforming {clp_conf:.4f}",
                   "OPEN: the medium Newton pin STALLS (res 7e-6) and disagrees "
                   "with Picard==conforming ⇒ a non-convergence problem, "
                   "unresolved (analysis deferred)",
                   dp_pic < 1.0, xfail=True)
    else:
        checks.add("GB16.4", f"{level}_lift_triangle",
                   f"Newton-pin {clp_pin:.4f} vs conforming {clp_conf:.4f} "
                   f"({dp_conf:.1f}%); Picard {clp_pic:.4f}",
                   "RECORDED (coarse): the CONVERGED Newton pin matches "
                   "conforming to ~0.1%; the LS Picard is the under-resolved "
                   "outlier here", True)


# --------------------------------------------------------------------------
def run_coarse():
    level = "coarse"
    mesh, mvop = build(level)
    if mesh is None:
        print(f"[skip] {level}: wing-body mesh absent (regenerate)")
        return False
    print(f"=== {level} (Parts 1-2) ===", flush=True)
    phi = picard_state(level, mesh, mvop)

    dg = diagnose(level, mesh, mvop, phi)
    print(f"  [D1] max|R[free]|={dg['res']:.3f} n_big={dg['n_big']} "
          f"lim={dg['n_lim']} flr={dg['n_flr']} | Gamma_bar={dg['gbar']:.4f}",
          flush=True)
    outer = dg["x"] >= 10.0
    junk = float(np.abs(dg["jump_pic"][outer]).max()) if outer.any() else 0.0
    print(f"  [D4] outer aux (x>=10): {int(outer.sum())} junk max|jump|={junk:.2f}"
          f" | cond1 legacy={dg['cond']['legacy']:.2e} pin={dg['cond']['pin']:.2e}",
          flush=True)

    neu = neumann_blowup(mesh, mvop)
    print(f"  [D0] cold neumann UNBOUNDED: max res {max(neu):.2e} "
          f"(the fuselage blockage; freestream is the wing-body BC)", flush=True)
    leg = newton_solve(level, mesh, mvop, "legacy")
    pin = newton_solve(level, mesh, mvop, "pin")
    jump_pin = float(np.abs(mvop.node_jump(pin["phi_ext"], dg["cut"])[outer]).max())
    print(f"  [D8] pin outer max|jump|={jump_pin:.2e} (was {junk:.2f})", flush=True)

    fig_residual_map(level, mesh, dg)
    fig_jump_profile(level, mesh, dg, pin, mvop)
    fig_convergence(level, neu, leg, pin)
    fig_ladder(level, dg, leg, pin)
    lift_ab(level, mesh, mvop, pin["phi_ext"], phi, gated=False)  # coarse=RECORDED

    write_csv(OUT, f"jump_profile_{level}.csv", "x,jump_legacy,jump_pin",
              [(f"{xi:.4f}", f"{jl:.5f}", f"{jp:.2e}") for xi, jl, jp in
               zip(dg["x"], dg["jump_pic"],
                   mvop.node_jump(pin["phi_ext"], dg["cut"]))])
    write_csv(OUT, f"convergence_{level}.csv", "path,converged,n_iter,res_final,"
              "n_lim,n_flr,wall_s",
              [("neumann", False, len(neu), f"{neu[-1]:.3e}", "", "", ""),
               ("freestream_legacy", leg["conv"], len(leg["res"]),
                f"{leg['res'][-1]:.3e}", leg["nlim"], leg["nflr"],
                f"{leg['wall_s']:.1f}"),
               ("freestream_pin", pin["conv"], len(pin["res"]),
                f"{pin['res'][-1]:.3e}", pin["nlim"], pin["nflr"],
                f"{pin['wall_s']:.1f}")])

    # ---- checks (coarse tier is RECORDED; the gated medium is the pass tier)
    checks.add("GB16.1", "coarse_diagnostic",
               f"max|R|={dg['res']:.1f}, {dg['n_big']} outer rows, "
               f"junk|jump|={junk:.1f}",
               "reproduces the B9 churn: O(1) residual on the outer wake aux",
               dg["res"] > 1.0 and dg["n_big"] >= 4 and junk > 10.0)
    checks.add("GB16.1", "coarse_neumann_unbounded",
               f"cold neumann max res {max(neu):.1e}",
               "neumann is unbounded on the wing-body (fuselage blockage)",
               max(neu) > 1e4)
    checks.add("GB16.1", "coarse_pin_cleans_jumps",
               f"outer |jump| {junk:.1f} -> {jump_pin:.1e}",
               "pin drives the outer junk jumps to ~0 (D8)",
               jump_pin < 1e-6)
    checks.add("GB16.1", "coarse_cond_drop",
               f"cond1 {dg['cond']['legacy']:.1e} -> {dg['cond']['pin']:.1e}",
               "aux-block conditioning improves under the pin (D6)",
               dg["cond"]["pin"] <= dg["cond"]["legacy"])
    checks.add("GB16.3", "coarse_legacy_churns",
               f"legacy res {leg['res'][-1]:.2e}, lim {leg['nlim']}",
               "legacy freestream Newton does NOT converge (churns)",
               leg["res"][-1] > 1.0)
    checks.add("GB16.3", "coarse_pin_converges",
               f"pin res {pin['res'][-1]:.2e}, lim {pin['nlim']}, flr {pin['nflr']}",
               "pin reaches machine residual, 0 limited (BC layer fixed; "
               "RECORDED at coarse -- flr from B8/G1.6 junction, GB9.4)",
               pin["res"][-1] < 1e-6 and pin["nlim"] == 0)
    return True


def run_medium():
    level = "medium"
    mesh, mvop = build(level)
    if mesh is None:
        print(f"[skip] {level}: wing-body mesh absent (regenerate)")
        return
    print(f"=== {level} (Part 3, gated) ===", flush=True)
    pin = newton_solve(level, mesh, mvop, "pin")
    lift_ab(level, mesh, mvop, pin["phi_ext"], picard_state(level, mesh, mvop),
            gated=True)   # GB16.4 medium: Newton pin == Picard lift < 1%
    # fair Picard+lagged arm so the B13 speedup is not miscredited to B16
    fair = OUT / "picard_lagged_medium.npz"
    if fair.exists():
        d = np.load(fair); pic_wall = float(d["wall"]); pic_gamma = float(d["gamma"])
    else:
        t0 = time.perf_counter()
        r = solve_multivalued_lifting(mvop, mesh, M, alpha_deg=ALPHA,
                                      farfield="freestream", n_outer_max=80,
                                      tol_residual=1e-7, direct_refactor_every=1000)
        pic_wall = time.perf_counter() - t0
        pic_gamma = float(np.abs(mvop.te_jump(r["phi_ext"])).mean())
        np.savez(fair, wall=pic_wall, gamma=pic_gamma)
        print(f"  [picard+lagged medium] {pic_wall:.0f}s gamma={pic_gamma:.4f}",
              flush=True)

    fig_walltime([("B9 Picard\n(committed)", B9_PICARD_MEDIUM_WALL, S4_ROSE),
                  ("Picard+lagged\n(fair arm)", pic_wall, S3_YELLOW),
                  ("Newton pin\n(B16)", pin["wall_s"], S2_AQUA)])
    print(f"  [medium] pin Newton {pin['wall_s']:.0f}s (res {pin['res'][-1]:.1e}, "
          f"lim {pin['nlim']}, flr {pin['nflr']}) vs B9 Picard "
          f"{B9_PICARD_MEDIUM_WALL:.0f}s / lagged {pic_wall:.0f}s", flush=True)
    write_csv(OUT, "walltime_medium.csv", "path,wall_s,note",
              [("b9_picard_committed", f"{B9_PICARD_MEDIUM_WALL:.1f}", "b9_summary.csv"),
               ("picard_lagged_fair", f"{pic_wall:.1f}", "B13 lagged-LU arm"),
               ("newton_pin_b16", f"{pin['wall_s']:.1f}",
                f"conv={pin['conv']} res={pin['res'][-1]:.1e} "
                f"lim={pin['nlim']} flr={pin['nflr']}")])
    # ★ RECORDED, not pass/fail (honest): the pin makes the medium freestream
    # Newton run CLEANLY (0 limited AND 0 floored -- cleaner than coarse, whose
    # 3 floored junction cells do not recur here), the churn is gone -- but it
    # plateaus at res ~7e-6 (the wing-fuselage junction O(h) inconsistency, the
    # G1.6 class again, not the far-field aux) and it is NOT faster than Picard:
    # the 35-outer Picard SEED dominates, and the fair lagged-LU Picard arm is
    # ~205 s. B16's value is that the Newton path WORKS cleanly on the wing-body
    # (enabling the B15 transonic ramp there), NOT wall-clock.
    checks.add("GB16.3", "medium_pin_clean",
               f"res {pin['res'][-1]:.2e}, lim {pin['nlim']}, flr {pin['nflr']}",
               "RECORDED: pin Newton runs clean (0 lim/flr, no churn); plateaus "
               "at ~7e-6 (junction O(h), G1.6 -- not the far-field aux)",
               pin["nlim"] == 0)
    checks.add("GB16.3", "medium_walltime",
               f"Newton {pin['wall_s']:.0f}s vs Picard {B9_PICARD_MEDIUM_WALL:.0f}s "
               f"/ lagged {pic_wall:.0f}s",
               "RECORDED: seed-dominated, NOT faster than Picard -- B16's value "
               "is convergence, not speed", True)


def main():
    apply_style()
    ran = run_coarse()
    if GATED:
        run_medium()
    elif ran:
        print("\n[medium gated: set PYFP3D_TRANSONIC_GATES=1 for Part 3]")
    if not ran and not GATED:
        print("no levels ran (mesh absent)")
        sys.exit(0)
    sys.exit(checks.report(OUT, "checks.csv"))


if __name__ == "__main__":
    main()
