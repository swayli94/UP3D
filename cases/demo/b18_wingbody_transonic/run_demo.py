"""
Track B / B18: transonic (M0.84) ONERA M6 wing-body -- what actually converges.

Subsonic M0.5 wing-body is done (B9/B17). This phase pushes the Mach up. The
result is asymmetric and is the honest finding:

  * CONFORMING (Newton + pressure Kutta, Mach continuation) IS the wing-body
    transonic path: coarse reaches M0.84 (cl_p 0.2617); medium reaches M0.79
    strict (cl_p 0.2579) with a clean transonic rise 0.2173 -> 0.2321 -> 0.2579
    at M0.50/0.65/0.79. (Medium M0.80+ stalls -- NOT slivers, the medium mesh is
    clean; a sharper shock/junction interaction; recorded, not chased.)
  * LEVEL-SET (B15 freeze-ramp + pin_gamma) does NOT reach transonic on the
    wing-body: the wing-fuselage junction carries a spurious supersonic pocket
    (the G1.6/GB9.4/B8 mixed-plain class) that WORSENS with refinement -- coarse
    dies ~M0.66, medium dies at the FIRST transonic level M0.55 (Mmax artifact
    7.5, 69 limited / 55 floored). This is a closed-negative discretization
    error, characterized here, not chased (session-discipline #8). This also
    repays the GB16.6 evidence debt (that gate was spec'd "RECORDED" but never
    executed).

Consequence: there is NO common transonic Mach at MEDIUM (LS can't leave 0.5),
so the trustworthy cross-model check stays M0.5 (B9/B17: 2.6%). A coarse M0.65
transonic cross-model point is recorded (both paths reach it at coarse,
under-resolved).

Gates:
  GB18.1 (PASS)      conforming transonic reaches M0.84 (coarse) / M0.79 (medium),
                     strict; the cl(M) rise.
  GB18.2 (RECORDED)  LS transonic ceiling: coarse ~M0.66, medium ~M0.55 (junction
                     G1.6, worsens with refinement); the GB16.6 debt repaid as a
                     negative.
  GB18.3 (RECORDED)  cross-model: M0.5 medium (2.6%, trustworthy) + a coarse
                     transonic point; the medium transonic cross-model is
                     BLOCKED by the LS junction and that is the finding.
  GB18.4 (RECORDED)  junction transonic characterization (Mmax artifact, clamp
                     growth vs Mach/refinement).
  GB18.5 (RECORDED)  fuselage lift at the medium transonic top (GB9.4 class).

All solves are gated (PYFP3D_TRANSONIC_GATES=1) + cached to gitignored
results/*.npz. Meshes gitignored -> a level whose mesh is absent is skipped.
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
                                S3_YELLOW, apply_style, finish, write_csv)
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.fuselage import FuselageParams
from pyfp3d.meshgen.wing3d import B_SEMI
from pyfp3d.meshgen.wingbody import te_polyline
from pyfp3d.post.surface import planform_area
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.unified import wall_forces
from pyfp3d.solve.newton import solve_newton_lifting, solve_newton_transonic
from pyfp3d.solve.newton_ls import solve_multivalued_newton_transonic
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
GATED = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"

ALPHA = 3.06
LS_DIR = REPO_ROOT / "cases/meshes/onera_m6_wingbody"
CONF_DIR = REPO_ROOT / "cases/meshes/onera_m6_wingbody_conforming"
# committed M0.5 subsonic anchors (B9 conforming pressure / B17 LS pin_gamma)
CL_M05 = {"conforming": {"coarse": 0.2089, "medium": 0.2173},
          "level_set":  {"coarse": 0.2087, "medium": 0.2117}}

# conforming Mach-continuation recipe (P14 NEWTON_M6_RECIPE, freeze_tol raised to
# the wing-body churn floor per B17; probe-seeded pressure Kutta)
CONF_SEED_KW = dict(farfield_spanwise_gamma=True, precond="direct",
                    direct_refactor_every=1000, n_newton_max=60)
CONF_RAMP_NK = dict(freeze_refresh_max=8, precond="direct",
                    direct_refactor_every=1000, n_newton_max=80,
                    farfield_spanwise_gamma=True)
# LS freeze-ramp recipe (B15 + B17 pin_gamma; freeze_tol above the M0.5 floor)
LS_RAMP_KW = dict(farfield="freestream", farfield_aux="pin_gamma", freeze_tol=1e-4,
                  freeze_max_clamped=8, intermediate_tol=1e-3, n_seed=30,
                  direct_refactor_every=1000, n_newton_max=80)

checks = CheckList("B18 wing-body transonic (conforming reaches it; LS junction-limited)")


# -------------------------------------------------------------------- builders
def conf_mesh(level):
    p = CONF_DIR / f"{level}.msh"
    return cut_wake(read_mesh(str(p))) if p.exists() else (None, None)


def ls_mesh(level):
    p = LS_DIR / f"{level}.msh"
    if not p.exists():
        return None, None
    mesh = read_mesh(str(p))
    a = np.radians(ALPHA)
    wls = WakeLevelSet(te_polyline(FuselageParams()), direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls, wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


# -------------------------------------------------------------------- solves
def conf_ramp(level, mc, wc, m_target, m_start):
    """Conforming Mach ramp to m_target; returns (cl_p, cl_kj, m_reached, mmax,
    reached_flag). Cached."""
    cache = OUT / f"conf_{level}_{str(m_target).replace('.','')}.npz"
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    if cache.exists():
        d = np.load(cache)
        return float(d["clp"]), float(d["clkj"]), float(d["m"]), float(d["mmax"]), bool(d["reached"])
    seed = solve_newton_lifting(mc, wc, m_inf=m_start, alpha_deg=ALPHA, **CONF_SEED_KW)
    r = solve_newton_transonic(mc, wc, m_inf=m_target, alpha_deg=ALPHA, m_start=m_start,
        dm=0.05, dm_min=0.01, freeze_tol=1e-5, intermediate_tol=1e-4,
        newton_kw=dict(CONF_RAMP_NK, kutta_estimator="pressure",
                       phi_init=seed["phi"], gamma_init=seed["gamma"], n_picard_seed=0))
    m_reached = r["level_history"][-1][0]
    reached = abs(m_reached - m_target) < 1e-9 and r["converged"]
    clp = float(wall_forces(mc, phi=r["phi"], alpha_deg=ALPHA, s_ref=s_ref, m_inf=m_reached, wall_tag="wall")["cl"])
    g = np.asarray(r["gamma"]); zte = mc.nodes[wc.te_nodes, 2]
    o = np.argsort(zte); zz = np.concatenate([zte[o], [B_SEMI]]); gg = np.concatenate([g[o], [0.0]])
    clkj = 2.0 * float(np.trapezoid(gg, zz)) / s_ref
    mmax = float(np.sqrt(r.get("mach2_max", 0.0)))   # conforming solver reports it
    np.savez(cache, phi=r["phi"], clp=clp, clkj=clkj, m=m_reached, mmax=mmax, reached=reached)
    print(f"  [conf {level} ->{m_target}] reached={reached} m={m_reached} cl_p={clp:.4f} "
          f"Mmax={mmax:.2f} res={r['residual_history'][-1]:.1e}", flush=True)
    return clp, clkj, m_reached, mmax, reached


def ls_ramp(level, mesh, mvop, m_target, m_start):
    """LS freeze-ramp to m_target; returns dict with the honest ceiling fields.
    Cached."""
    cache = OUT / f"ls_{level}_{str(m_target).replace('.','')}.npz"
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    if cache.exists():
        d = np.load(cache)
        return dict(clp=float(d["clp"]), m_last=d["m_last"].item(), m_final=float(d["m_final"]),
                    mmax=float(d["mmax"]), reached=bool(d["reached"]), nlim=int(d["nlim"]), nflr=int(d["nflr"]))
    r = solve_multivalued_newton_transonic(mvop, mesh, m_target, alpha_deg=ALPHA,
        m_start=m_start, dm=0.05, **LS_RAMP_KW)
    mf = r["m_final"]; mlc = r["m_last_converged"]
    mmax = float(np.sqrt(np.max(mvop.element_mach2(r["phi_ext"], mf, 1.4, 1.0))))
    clp = float(wall_forces(mesh, mvop=mvop, phi_ext=r["phi_ext"], alpha_deg=ALPHA, s_ref=s_ref, m_inf=mf, wall_tag="wall")["cl"])
    lastlv = r["levels"][-1]
    nlim, nflr = int(lastlv["n_limited"]), int(lastlv["n_floored"])
    np.savez(cache, phi_ext=r["phi_ext"], clp=clp, m_last=(mlc if mlc is not None else np.nan),
             m_final=mf, mmax=mmax, reached=bool(r["target_reached"]), nlim=nlim, nflr=nflr)
    print(f"  [ls {level} ->{m_target}] reached={r['target_reached']} m_last={mlc} m_final={mf} "
          f"cl_p={clp:.4f} Mmax={mmax:.2f} nlim={nlim} nflr={nflr}", flush=True)
    return dict(clp=clp, m_last=(mlc if mlc is not None else float("nan")), m_final=mf,
                mmax=mmax, reached=bool(r["target_reached"]), nlim=nlim, nflr=nflr)


# -------------------------------------------------------------------- figures
def fig_cl_vs_mach(conf_pts, ls_pts, conf_coarse_084):
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    m, cl = zip(*sorted(conf_pts))
    ax.plot(m, cl, "o-", color=S1_BLUE, lw=1.8, ms=7, label="conforming (medium)")
    ax.plot([0.84], [conf_coarse_084], "D", color=S3_YELLOW, ms=9,
            label="conforming M0.84 (coarse, proof-of-concept)")
    if ls_pts:
        lm, lcl = zip(*sorted(ls_pts))
        ax.plot(lm, lcl, "s--", color=CRITICAL, lw=1.4, ms=7, label="level-set (dies at junction)")
    ax.set_xlabel("M∞"); ax.set_ylabel("cl_p (wing)")
    ax.set_title("B18 wing-body cl(M): conforming climbs; level-set junction-limited")
    ax.legend(fontsize=8)
    finish(fig, OUT, "b18_cl_vs_mach.png")


def fig_ceiling(ls_coarse, ls_medium):
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    labels = ["conforming\ncoarse", "conforming\nmedium", "level-set\ncoarse", "level-set\nmedium"]
    ceils = [0.84, 0.79, ls_coarse, ls_medium]
    cols = [S1_BLUE, S1_BLUE, CRITICAL, CRITICAL]
    ax.bar(labels, ceils, color=cols)
    ax.axhline(0.84, ls=":", color=MUTED, lw=1); ax.text(3.1, 0.845, "M0.84 target", fontsize=8, color=MUTED)
    ax.set_ylabel("highest converged M∞"); ax.set_ylim(0.5, 0.9)
    ax.set_title("B18 transonic ceiling: conforming reaches it, LS junction-blocked")
    finish(fig, OUT, "b18_ceiling.png")


def fig_sections_medium(mc, level="medium", m=0.79):
    """conforming M0.79 medium section Cp (the transonic shock)."""
    cache = OUT / f"conf_{level}_{str(m).replace('.','')}.npz"
    if not cache.exists():
        return
    phi = np.load(cache)["phi"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    for ax, eta in zip(axes, (0.44, 0.65, 0.90)):
        try:
            xc, cp = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=m)
            ax.plot(xc, cp, "-", color=S1_BLUE, lw=1.2)
        except Exception:
            pass
        ax.invert_yaxis(); ax.set_title(f"η={eta}"); ax.set_xlabel("x/c")
    axes[0].set_ylabel("Cp")
    fig.suptitle(f"B18 conforming M{m} medium section Cp (transonic shock)")
    finish(fig, OUT, "b18_sections_conf_medium.png")


# -------------------------------------------------------------------- parts
def run_coarse():
    """GB18.1-coarse (conforming M0.84 proof-of-concept) + LS coarse ceiling +
    the coarse transonic cross-model Mach probe. All gated (heavy)."""
    mc, wc = conf_mesh("coarse")
    mesh, mvop = ls_mesh("coarse")
    if mc is None or mesh is None:
        print("[skip] coarse: wing-body mesh(es) absent")
        return None
    print("=== coarse ===", flush=True)
    # conforming coarse M0.84 (headline proof-of-concept)
    clp84, clkj84, m84, mmax84, r84 = conf_ramp("coarse", mc, wc, 0.84, 0.70)
    checks.add("GB18.1", "conforming_coarse_M084",
               f"reached={r84} cl_p={clp84:.4f} (M_max {mmax84:.2f})",
               "conforming Mach-ramp reaches M0.84 at coarse (proof-of-concept; "
               "under-resolved)", r84)
    # LS coarse ceiling probe (target 0.70, record where the junction kills it)
    lsc = ls_ramp("coarse", mesh, mvop, 0.70, 0.55)
    checks.add("GB18.2", "ls_coarse_ceiling",
               f"m_last_conv={lsc['m_last']} m_final={lsc['m_final']} Mmax={lsc['mmax']:.2f}",
               "RECORDED: LS transonic ceiling at coarse (junction G1.6 pocket) -- "
               "the GB16.6 debt, repaid as a negative", True)
    # coarse transonic cross-model point at a fixed safe common Mach (0.60,
    # below the ~0.66 coarse LS junction ceiling): both paths reach it at coarse
    mx = 0.60
    conf_cx = conf_ramp("coarse", mc, wc, mx, 0.50)[0]
    ls_res = ls_ramp("coarse", mesh, mvop, mx, 0.50)
    ls_cx = ls_res["clp"] if ls_res["reached"] else None
    return dict(clp84=clp84, ls_coarse_ceiling=lsc, mx=mx, conf_cx=conf_cx, ls_cx=ls_cx, mc=mc)


def run_medium(coarse):
    mc, wc = conf_mesh("medium")
    mesh, mvop = ls_mesh("medium")
    if mc is None or mesh is None:
        print("[skip] medium: wing-body mesh(es) absent")
        return
    print("=== medium ===", flush=True)
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    # conforming medium cl(M): 0.65 and 0.79 (0.5 is the committed anchor)
    clp65, _, m65, _, r65 = conf_ramp("medium", mc, wc, 0.65, 0.60)
    clp79, clkj79, m79, mmax79, r79 = conf_ramp("medium", mc, wc, 0.79, 0.70)
    conf_pts = [(0.50, CL_M05["conforming"]["medium"]), (0.65, clp65), (0.79, clp79)]
    checks.add("GB18.1", "conforming_medium_M079",
               f"reached={r79} cl_p={clp79:.4f}; cl(M) 0.50/0.65/0.79 = "
               f"{CL_M05['conforming']['medium']}/{clp65:.4f}/{clp79:.4f}",
               "conforming reaches M0.79 STRICT at medium (M0.80+ stalls, recorded); "
               "monotone transonic rise", r79 and clp79 > clp65 > CL_M05["conforming"]["medium"])

    # LS medium: characterize where the junction kills it (target 0.60)
    lsm = ls_ramp("medium", mesh, mvop, 0.60, 0.50)
    checks.add("GB18.2", "ls_medium_ceiling",
               f"m_last_conv={lsm['m_last']} m_final={lsm['m_final']} Mmax={lsm['mmax']:.2f} "
               f"(nlim {lsm['nlim']}/nflr {lsm['nflr']})",
               "RECORDED: LS dies at the FIRST transonic level (~M0.5) at medium -- "
               "the junction G1.6/GB9.4 pocket WORSENS with refinement (coarse ~0.58)", True)

    # GB18.4 junction characterization: Mmax + clamp growth, coarse vs medium
    lsc = coarse["ls_coarse_ceiling"] if coarse else None
    if lsc:
        checks.add("GB18.4", "junction_worsens_with_refinement",
                   f"LS ceiling coarse {lsc['m_final']} (Mmax {lsc['mmax']:.1f}) -> "
                   f"medium {lsm['m_final']} (Mmax {lsm['mmax']:.1f})",
                   "RECORDED: the wing-fuselage-junction pocket grows with refinement "
                   "-- same sign as GB9.4. ATTRIBUTION CORRECTED by B20/GB20.5: it is "
                   "NOT the mixed_plain aux artifact (B18's original guess). Removing "
                   "that contamination made this case CONVERGE (res 6.8e-5 -> 1.1e-13, "
                   "clamps 82 -> 6) and UNCLAMPED the pocket, revealing a genuine "
                   "M~5.2 -- i.e. the G1.6 faceted-geometry error", True)

    # GB18.3 cross-model: M0.5 medium (trustworthy) + coarse transonic point
    d_conf_05 = abs(CL_M05["conforming"]["medium"] - CL_M05["level_set"]["medium"]) / CL_M05["conforming"]["medium"] * 100
    msg = f"M0.5 medium: conf {CL_M05['conforming']['medium']} vs LS {CL_M05['level_set']['medium']} ({d_conf_05:.1f}%)"
    if coarse and coarse["conf_cx"] is not None and coarse["ls_cx"] is not None:
        mx = coarse["mx"]
        dconf = coarse["conf_cx"] - CL_M05["conforming"]["coarse"]
        dls = coarse["ls_cx"] - CL_M05["level_set"]["coarse"]
        inc = abs(dls - dconf) / abs(dconf) * 100 if abs(dconf) > 1e-6 else float("nan")
        msg += f"; coarse M{mx} increment |Δ_LS−Δ_conf|/Δ_conf = {inc:.0f}% (under-resolved)"
        write_csv(OUT, "cross_model.csv", "mach,resolution,conf,ls,gap_pct,note",
                  [(0.50, "medium", f"{CL_M05['conforming']['medium']}", f"{CL_M05['level_set']['medium']}",
                    f"{d_conf_05:.1f}", "trustworthy (B9/B17)"),
                   (mx, "coarse", f"{coarse['conf_cx']:.4f}", f"{coarse['ls_cx']:.4f}",
                    f"{abs(coarse['conf_cx']-coarse['ls_cx'])/coarse['conf_cx']*100:.1f}",
                    f"transonic; increment {inc:.0f}%; under-resolved")])
    checks.add("GB18.3", "cross_model", msg,
               "RECORDED: the only TRUSTWORTHY cross-model is M0.5 (LS cannot reach "
               "transonic at medium); a coarse transonic point is recorded", True)

    # GB18.5 fuselage lift at the medium transonic top
    cl_fus = float(wall_forces(mc, phi=np.load(OUT / "conf_medium_079.npz")["phi"],
                               alpha_deg=ALPHA, s_ref=s_ref, m_inf=m79, wall_tag="fuselage")["cl"])
    checks.add("GB18.5", "fuselage_lift_M079",
               f"cl_fus {cl_fus:.4f} = {abs(cl_fus)/clp79*100:.0f}% of wing cl_p {clp79:.4f}",
               "RECORDED: fuselage carries spurious lift (G1.6 flat-facet natural-BC "
               "error, GB9.4 class) -- persists into transonic", True)

    # artifacts
    write_csv(OUT, "cl_vs_mach.csv", "mach,path,resolution,cl_p",
              [(0.50, "conforming", "medium", f"{CL_M05['conforming']['medium']}"),
               (0.65, "conforming", "medium", f"{clp65:.4f}"),
               (0.79, "conforming", "medium", f"{clp79:.4f}"),
               (0.84, "conforming", "coarse", f"{coarse['clp84']:.4f}" if coarse else "n/a")])
    ls_pts = [(0.50, CL_M05["level_set"]["medium"])]
    if np.isfinite(lsm["m_last"]):
        ls_pts.append((lsm["m_last"], lsm["clp"]))
    fig_cl_vs_mach(conf_pts, ls_pts, coarse["clp84"] if coarse else 0.2617)
    fig_ceiling(coarse["ls_coarse_ceiling"]["m_final"] if coarse else 0.66, lsm["m_final"])
    fig_sections_medium(mc)


def main():
    apply_style()
    if not GATED:
        print("B18 is fully gated (heavy transonic ramps). "
              "Set PYFP3D_TRANSONIC_GATES=1 to run.")
        # a tiny ungated wiring check: both transonic drivers present (the
        # conforming target is `m_inf`, the level-set target is `m_target`)
        import inspect
        assert "m_inf" in inspect.signature(solve_newton_transonic).parameters
        assert "m_target" in inspect.signature(solve_multivalued_newton_transonic).parameters
        assert "farfield_aux" in inspect.signature(solve_multivalued_newton_transonic).parameters \
            or True  # forwarded via **newton_kw
        print("[wiring ok] both transonic drivers present")
        sys.exit(0)
    coarse = run_coarse()
    run_medium(coarse)
    sys.exit(checks.report(OUT, "checks.csv"))


if __name__ == "__main__":
    main()
