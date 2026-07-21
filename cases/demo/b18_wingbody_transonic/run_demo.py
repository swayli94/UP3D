"""
Track B / B18 -- transonic (M0.84) ONERA M6 wing-body, REFRESHED by B29
(2026-07-20): the level-set PRODUCTION side is now the B28 flat fragment
(inboard clip + sheet geometry dragged flat at y=0). REFRESHED again by
B32 (2026-07-22): the conforming PRODUCTION recipe now carries the
P13/G13.2 tip taper (vanish_smooth, r_c = 0.05*b_semi) through the B31
pressure-row blend -- user-adjudicated from the B31 verdict (the taper
cures the 0.83 dying level on the analysis grid; adoption cost cl_p ~ -3%
at the top legs, reported per leg).

Historical layers: the B18 baseline said "conforming reaches it, LS is
junction-limited" (a spurious wing-fuselage junction supersonic pocket,
closed-negative). B25/B26/B27 changed the verdict: the pocket is the B23
inboard free-edge singularity of the wake sheet, and clipping the sheet to
the conforming fragment topology (``CutElementMap(inboard_clip=
make_inboard_clip(FUS))``) heals it, co-locating the LS ceiling with the
conforming one. B28 then proved the residual cross-model fuselage-lift gap
is wake-sheet POSITION sensitivity (not an error): dragging the LS sheet
GEOMETRY flat at y=0 (``WakeLevelSet(sheet_direction=(1,0,0))``) while the
physics convects with the flow makes the out-band cl_fus agree with the
conforming oracle to 7.25%. B29 adopts that flat fragment as the production
LS configuration here (user-adjudicated, B28 VERDICT section 6) and
re-solves the C side; the conforming legs and the side-A historical
comparison are untouched (committed caches bit-reproduce). B32 adopts the
tip taper as the production CONFORMING configuration (user-adjudicated,
B31 VERDICT section 8) and re-solves the conforming legs (NEW conf_taper_*
caches; the pre-adoption conf_* caches are stale); the LS legs and the
side-A comparison are untouched (their caches bit-reproduce).

  * CONFORMING (Newton + pressure Kutta + tip taper, Mach continuation) --
    the reference path: coarse reaches M0.84; medium reaches M0.79 strict
    with a clean transonic rise at M0.50/0.65/0.75/0.79 plus the B32
    ceiling-climb record leg toward 0.84 (the B31 0.83 dying level is
    cured by the production taper on the analysis grid). The M0.50
    conforming anchor is re-solved LIVE with the taper (the B9 no-taper
    value is historical).
  * LEVEL-SET without the clip (side A = the B18 recipe re-run on the current
    code): STILL pocket-limited -- medium dies at M0.5125, class (a) (the
    pocket erupts at 0.55: Mmax 13.1, beyond freeze_max_clamped=8); coarse
    dies at M0.84, class (b) (the junction-strip pocket stalls Newton). The
    A side climbs past the historical B18 committed anchors (died 0.50/0.55):
    that difference is the B21/B22 freeze-capture fix (the B26 T1 finding),
    not physics drift -- the pocket's true kill line on the A medium side is
    0.55 (Mmax 13.1 > freeze_max_clamped=8).
  * LEVEL-SET production (side C = B25 inboard clip + B28 FLAT sheet): the
    sheet geometry is dragged along +x at y=0 while the physics convects
    with the flow -- the cross-model-consistent configuration (B28-F1).
    The M0.5 subsonic anchors move to the re-run B9 demo values
    (0.2115 / 0.2184), tightening the M0.5 cross-model gap 2.6% -> 0.5%.
  * CROSS-MODEL: M0.5 (B32 taper conf / B28 flat LS) + M0.65 (gated <= 5%)
    + M0.75 (recorded, shock-sensitive) at medium, plus the coarse M0.60
    point (under-resolved).

Gates:
  GB18.1 (PASS)      conforming transonic reaches M0.84 (coarse) / M0.79
                     (medium), strict; the monotone cl(M) rise; the B32
                     climb leg records the post-adoption medium ceiling
                     (>= 0.80 expected per the B31 analysis-grid cure).
  GB18.2 (PASS)      LS ceiling A vs C: C coarse reaches 0.84 while A dies
                     below it; C medium climbs past 0.70 while A dies
                     at/below 0.55. The pocket was the limiter.
  GB18.3 (PASS/REC)  cross-model: M0.65 medium |cl_p gap| <= 5% is a REAL
                     gate; M0.5 anchor, M0.75 medium and the coarse
                     M0.60 point are recorded.
  GB18.4 (RECORDED)  pocket attribution: the A side cites the committed B26
                     evidence (cases/analysis/b26_ls_transonic_ceiling/
                     results/g1_summary.csv + g1_peaks.csv -- the A-side
                     dying peak sits in the junction strip x ~ 2.1-2.25,
                     q ~ 0, next to the fuselage); the C-side dying peak is
                     measured LIVE on the flat sheet (expected: the wing TIP,
                     P13 class). The flat-side junction corridor at M0.5 was
                     measured clean in B28 (corrM 0.64, n_sup 0).
  GB18.5 (RECORDED)  fuselage lift: the conforming cl_fus at the medium
                     transonic top (live) + the LIVE flat C-side ceiling
                     decomposition (pocket band / out-band / polar caps, the
                     B23 W2 definitions). The B26 tilted-sheet "x2 out-band"
                     reading (P11 watch item) is RETIRED by B28: sheet
                     position sensitivity, not a lesion.

All solves are gated (PYFP3D_TRANSONIC_GATES=1) + cached to gitignored
results/*.npz (the flat C side uses NEW ls_flat_* caches; the B26 tilted
ls_C_* caches are stale and unused; the B32 conforming taper legs use NEW
conf_taper_* caches, the pre-adoption conf_* ones are stale). Meshes
gitignored -> a level whose mesh is absent is skipped.
"""

import json
import os
import sys
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
from pyfp3d.constraints.wake import tip_taper_factors
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.fuselage import FuselageParams, make_inboard_clip
from pyfp3d.meshgen.wing3d import B_SEMI
from pyfp3d.meshgen.wingbody import junction_z, te_polyline
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.surface import _cp_from_q2, _pressure_force, planform_area
from pyfp3d.post.surface_ls import _d11_wall_state
from pyfp3d.post.unified import wall_forces
from pyfp3d.solve.continuation import mach_schedule
from pyfp3d.solve.newton import solve_newton_lifting, solve_newton_transonic
from pyfp3d.solve.newton_ls import solve_multivalued_newton_transonic
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
GATED = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"

ALPHA = 3.06
FUS = FuselageParams()
Z_JUNC = junction_z(FUS)
LS_DIR = REPO_ROOT / "cases/meshes/onera_m6_wingbody"
CONF_DIR = REPO_ROOT / "cases/meshes/onera_m6_wingbody_conforming"
# committed M0.5 subsonic anchors. LS = B28 flat-fragment (the re-run B9
# demo); conforming = HISTORICAL pre-adoption record (B9 pressure, no
# taper) -- superseded by the B32 taper adoption: the demo re-solves the
# conforming M0.5 leg LIVE with the production taper below and uses that
# value everywhere (this dict's conforming entries are kept for the
# coarse increment reference + the paper trail).
CL_M05 = {"conforming": {"coarse": 0.2089, "medium": 0.2173},
          "level_set":  {"coarse": 0.2115, "medium": 0.2184}}

# conforming Mach-continuation recipe (P14 NEWTON_M6_RECIPE, freeze_tol raised to
# the wing-body churn floor per B17; probe-seeded pressure Kutta)
CONF_SEED_KW = dict(farfield_spanwise_gamma=True, precond="direct",
                    direct_refactor_every=1000, n_newton_max=60)
CONF_RAMP_NK = dict(freeze_refresh_max=8, precond="direct",
                    direct_refactor_every=1000, n_newton_max=80,
                    farfield_spanwise_gamma=True)
# B32 tip-taper adoption (user-adjudicated from the B31 verdict): the
# conforming production recipe unloads the tip Gamma with the P13/G13.2
# vanish_smooth taper (r_c = 0.05*b_semi) through the B31 pressure-row
# blend (solve/newton.py NewtonWorkspace docstring; the weld sign is
# per-step refreshed since B32). Applied uniformly over M. NEW
# conf_taper_* caches -- the pre-adoption conf_* caches are stale (the
# ls_flat_* precedent).
CONF_TAPER = ("vanish_smooth", 0.05)
# LS freeze-ramp recipe (B15 + B17 pin_gamma; freeze_tol above the M0.5 floor)
# -- frozen verbatim since B18; B26 re-measured the ceiling with it.
LS_RAMP_KW = dict(farfield="freestream", farfield_aux="pin_gamma", freeze_tol=1e-4,
                  freeze_max_clamped=8, intermediate_tol=1e-3, n_seed=30,
                  direct_refactor_every=1000, n_newton_max=80)
DM = 0.05

checks = CheckList("B18 wing-body transonic, B29 refresh: LS production = "
                   "flat-fragment (inboard clip + flat sheet); B32 refresh: "
                   "conforming production = pressure Kutta + tip taper")


# -------------------------------------------------------------------- builders
def conf_mesh(level):
    p = CONF_DIR / f"{level}.msh"
    return cut_wake(read_mesh(str(p))) if p.exists() else (None, None)


def ls_mesh(level, clip=False):
    """LS mesh + multivalued operator. clip=False = the B18 default (side A,
    q>=0 inboard clip, the junction pocket; tilted sheet -- the historical
    comparison side). clip=True = the PRODUCTION side C (B29): the B25
    inboard fragment clip + the B28 flat sheet (geometry dragged along +x
    at y=0, physics convecting with the flow)."""
    p = LS_DIR / f"{level}.msh"
    if not p.exists():
        return None, None
    mesh = read_mesh(str(p))
    a = np.radians(ALPHA)
    wls = WakeLevelSet(te_polyline(FUS), direction=(np.cos(a), np.sin(a), 0.0),
                       sheet_direction=(1.0, 0.0, 0.0) if clip else None)
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]),
                       inboard_clip=make_inboard_clip(FUS) if clip else None)
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


# ------------------------------------------------------ fuselage decomposition
# (GB18.5, B29): the B23 W2 definitions -- pocket band |z - Z_JUNC| < 0.06 &
# x > 1.0 / outside / polar caps (0.10); same core as wall_forces.
X_BAND0, BW, POLE_MARGIN = 1.0, 0.06, 0.10


def fuselage_parts(mesh, mvop, phi_ext, s_ref, m_inf):
    """cl_fus split into pocket-band / out-band / polar-cap contributions on
    the LS side (the `_d11_wall_state` + `_cp_from_q2` + `_pressure_force`
    core, B23 W2; the mvop selects the sheet sides)."""
    wall = np.asarray(mesh.boundary_faces["fuselage"], dtype=np.int64)
    q2, _, area, n_out = _d11_wall_state(mesh, mvop, phi_ext, wall, 1.0)
    cp = _cp_from_q2(q2, m_inf, 1.4)
    cents = mesh.nodes[wall].mean(axis=1)

    def cl_of(mask):
        _, cl, _ = _pressure_force(cp[mask], area[mask], n_out[mask],
                                   s_ref, ALPHA)
        return float(cl)

    band = (np.abs(cents[:, 2] - Z_JUNC) < BW) & (cents[:, 0] > X_BAND0)
    x_nose = FUS.x_center - FUS.length / 2.0
    x_tail = FUS.x_center + FUS.length / 2.0
    pole = (cents[:, 0] < x_nose + POLE_MARGIN) | \
           (cents[:, 0] > x_tail - POLE_MARGIN)
    return dict(band=cl_of(band), out=cl_of(~band), poles=cl_of(pole))


# -------------------------------------------------------------------- solves
def conf_ramp(level, mc, wc, m_target, m_start):
    """Conforming Mach ramp to m_target (B32: the production recipe carries
    the tip taper -- CONF_TAPER above); returns (cl_p, cl_kj, m_reached,
    mmax, reached_flag, n_limited, n_floored). Cached as conf_taper_*."""
    cache = OUT / f"conf_taper_{level}_{str(m_target).replace('.','')}.npz"
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    if cache.exists():
        d = np.load(cache)
        return (float(d["clp"]), float(d["clkj"]), float(d["m"]),
                float(d["mmax"]), bool(d["reached"]),
                int(d["nlim"]), int(d["nflr"]))
    seed = solve_newton_lifting(mc, wc, m_inf=m_start, alpha_deg=ALPHA, **CONF_SEED_KW)
    taper = tip_taper_factors(wc.station_z, B_SEMI, CONF_TAPER[0],
                              CONF_TAPER[1] * B_SEMI)
    r = solve_newton_transonic(mc, wc, m_inf=m_target, alpha_deg=ALPHA, m_start=m_start,
        dm=DM, dm_min=0.01, freeze_tol=1e-5, intermediate_tol=1e-4,
        newton_kw=dict(CONF_RAMP_NK, kutta_estimator="pressure",
                       tip_taper=taper,
                       phi_init=seed["phi"], gamma_init=seed["gamma"], n_picard_seed=0))
    m_reached = r["level_history"][-1][0]
    reached = abs(m_reached - m_target) < 1e-9 and r["converged"]
    clp = float(wall_forces(mc, phi=r["phi"], alpha_deg=ALPHA, s_ref=s_ref, m_inf=m_reached, wall_tag="wall")["cl"])
    g = np.asarray(r["gamma"]); zte = mc.nodes[wc.te_nodes, 2]
    o = np.argsort(zte); zz = np.concatenate([zte[o], [B_SEMI]]); gg = np.concatenate([g[o], [0.0]])
    clkj = 2.0 * float(np.trapezoid(gg, zz)) / s_ref
    mmax = float(np.sqrt(r.get("mach2_max", 0.0)))   # conforming solver reports it
    nlim, nflr = int(r["n_limited"]), int(r["n_floored"])
    np.savez(cache, phi=r["phi"], gamma=r["gamma"], clp=clp, clkj=clkj,
             m=m_reached, mmax=mmax, reached=reached, nlim=nlim, nflr=nflr)
    print(f"  [conf {level} ->{m_target}] reached={reached} m={m_reached} cl_p={clp:.4f} "
          f"Mmax={mmax:.2f} clamps={nlim}+{nflr} res={r['residual_history'][-1]:.1e}", flush=True)
    return clp, clkj, m_reached, mmax, reached, nlim, nflr


_LV_KEYS = ("m_inf", "tag", "gamma", "cl_kj", "mach_max", "n_newton",
            "converged", "accept_reason", "residual_norm", "n_limited",
            "n_floored", "froze", "n_freeze_refresh", "n_freeze_reverts",
            "n_refactor", "n_schur_fallback", "wall_s", "n_lin_iters",
            "n_lin_solves")


def _classify(reached, last_level, n_halvings):
    """B26 failure taxonomy: (a) strict-gate rejection (res < 1e-6 but
    clamps > 0); (b) Newton non-convergence; suffix "+dm" after the halving
    cascade."""
    if reached:
        return "reached"
    tight = last_level["residual_norm"] < 1e-6
    clamps = last_level["n_limited"] + last_level["n_floored"]
    cls = "a" if (tight and clamps > 0) else "b"
    return cls + ("+dm" if n_halvings > 0 else "")


def ls_ramp(side, level, mesh, mvop, m_target, m_start):
    """LS freeze-ramp to m_target (the B18/B26 recipe, frozen verbatim).
    side = "A" (default clip, pocket, tilted sheet -- historical) or "C"
    (PRODUCTION: inboard fragment clip + flat sheet, B29). Returns dict with
    the honest ceiling fields + the slimmed per-level records (for the
    consistency diff vs the committed g1_levels). Cached -- the flat C side
    uses NEW ls_flat_* caches (the B26 tilted ls_C_* ones are stale)."""
    tag = "flat" if side == "C" else side
    cache = OUT / f"ls_{tag}_{level}_{str(m_target).replace('.','')}.npz"
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        return dict(clp=float(d["clp"]), m_last=d["m_last"].item(), m_final=float(d["m_final"]),
                    mmax=float(d["mmax"]), reached=bool(d["reached"]), nlim=int(d["nlim"]),
                    nflr=int(d["nflr"]), res=float(d["res"]), cls=str(d["cls"]),
                    levels=json.loads(str(d["levels_json"])),
                    mach_schedule=[float(x) for x in d["mach_schedule"]])
    r = solve_multivalued_newton_transonic(mvop, mesh, m_target, alpha_deg=ALPHA,
        m_start=m_start, dm=DM, **LS_RAMP_KW)
    mf = r["m_final"]; mlc = r["m_last_converged"]
    mmax = float(np.sqrt(np.max(mvop.element_mach2(r["phi_ext"], mf, 1.4, 1.0))))
    clp = float(wall_forces(mesh, mvop=mvop, phi_ext=r["phi_ext"], alpha_deg=ALPHA, s_ref=s_ref, m_inf=mf, wall_tag="wall")["cl"])
    lastlv = r["levels"][-1]
    nlim, nflr = int(lastlv["n_limited"]), int(lastlv["n_floored"])
    res = float(lastlv["residual_norm"])
    levels = [{k: lv[k] for k in _LV_KEYS} for lv in r["levels"]]
    sched = [float(x) for x in r["mach_schedule"]]
    n_halvings = len(sched) - len(mach_schedule(m_target, m_start=m_start, dm=DM))
    cls = _classify(bool(r["target_reached"]), lastlv, n_halvings)
    np.savez(cache, phi_ext=r["phi_ext"], clp=clp, m_last=(mlc if mlc is not None else np.nan),
             m_final=mf, mmax=mmax, reached=bool(r["target_reached"]), nlim=nlim, nflr=nflr,
             res=res, cls=cls, levels_json=json.dumps(levels), mach_schedule=np.asarray(sched))
    print(f"  [ls {side} {level} ->{m_target}] reached={r['target_reached']} m_last={mlc} "
          f"m_final={mf} cl_p={clp:.4f} Mmax={mmax:.2f} cls={cls} "
          f"(nlim {nlim}/nflr {nflr})", flush=True)
    return dict(clp=clp, m_last=(mlc if mlc is not None else float("nan")), m_final=mf,
                mmax=mmax, reached=bool(r["target_reached"]), nlim=nlim, nflr=nflr,
                res=res, cls=cls, levels=levels, mach_schedule=sched)


# -------------------------------------------------------------------- figures
def fig_cl_vs_mach(conf_pts, ls_pts, conf_coarse_084, ls_coarse_084, ls_ceiling_med, a_deaths):
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    m, cl = zip(*sorted(conf_pts))
    ax.plot(m, cl, "o-", color=S1_BLUE, lw=1.8, ms=7, label="conforming (medium)")
    lm, lcl = zip(*sorted(ls_pts))
    ax.plot(lm, lcl, "s-", color=S2_AQUA, lw=1.8, ms=7, label="level-set flat-fragment (medium)")
    ax.plot([0.84], [conf_coarse_084], "D", color=S3_YELLOW, ms=9,
            label="conforming M0.84 (coarse)")
    ax.plot([0.84], [ls_coarse_084], "D", color=S2_AQUA, ms=9, mfc="none", mew=2,
            label="level-set flat-frag M0.84 (coarse)")
    mx, mcl = ls_ceiling_med
    ax.plot([mx], [mcl], "s", color=S2_AQUA, ms=10, mfc="none", mew=2.2,
            label=f"LS flat-frag medium ceiling {mx:.4g} (cl_p at the dying state)")
    for mx, lbl in a_deaths:
        ax.plot([mx], [0.18], "x", color=CRITICAL, ms=10, mew=2.2)
        ax.annotate(lbl, (mx, 0.18), textcoords="offset points", xytext=(4, 8),
                    fontsize=8, color=CRITICAL)
    ax.set_ylim(0.16, 0.30)
    ax.set_xlabel("M∞"); ax.set_ylabel("cl_p (wing)")
    ax.set_title("B18/B29 wing-body cl(M): pocket-healed flat-frag LS climbs to the conforming ceiling")
    ax.legend(fontsize=8)
    finish(fig, OUT, "b18_cl_vs_mach.png")


def fig_ceiling(a_coarse, c_coarse, a_medium, c_medium, conf_medium=0.79):
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    labels = ["conforming", "level-set\nA (no clip)", "level-set\nC (flat frag)"]
    xs = np.arange(3)
    w = 0.36
    vals_coarse = [0.84, a_coarse, c_coarse]
    vals_medium = [conf_medium, a_medium, c_medium]
    ax.bar(xs - w / 2, vals_coarse, width=w, color=S3_YELLOW, label="coarse")
    ax.bar(xs + w / 2, vals_medium, width=w, color=S1_BLUE, label="medium")
    for x, v in zip(xs - w / 2, vals_coarse):
        ax.text(x, v + 0.008, f"{v:.3g}", ha="center", fontsize=8)
    for x, v in zip(xs + w / 2, vals_medium):
        ax.text(x, v + 0.008, f"{v:.3g}", ha="center", fontsize=8)
    ax.axhline(0.84, ls=":", color=MUTED, lw=1)
    ax.set_xticks(xs); ax.set_xticklabels(labels)
    ax.set_ylabel("highest converged M∞"); ax.set_ylim(0.45, 0.9)
    ax.set_title("B18/B29 transonic ceiling: the pocket was the limiter -- "
                 "flat-frag LS reaches the conforming site")
    ax.legend(fontsize=8)
    finish(fig, OUT, "b18_ceiling.png")


def fig_sections_medium(mc, level="medium", m=0.79):
    """conforming M0.79 medium section Cp (the transonic shock)."""
    cache = OUT / f"conf_taper_{level}_{str(m).replace('.','')}.npz"
    if not cache.exists():
        return
    phi = np.load(cache)["phi"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    for ax, eta in zip(axes, (0.44, 0.65, 0.90)):
        try:
            d = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=m)
            ax.plot(d["x_upper"], d["cp_upper"], "-", color=S1_BLUE, lw=1.2,
                    label="upper")
            ax.plot(d["x_lower"], d["cp_lower"], "-", color=S2_AQUA, lw=1.2,
                    label="lower")
        except Exception:
            pass
        ax.invert_yaxis(); ax.set_title(f"η={eta}"); ax.set_xlabel("x/c")
    axes[0].set_ylabel("Cp")
    axes[0].legend(fontsize=8)
    fig.suptitle(f"B18 conforming M{m} medium section Cp (transonic shock)")
    finish(fig, OUT, "b18_sections_conf_medium.png")


# -------------------------------------------------------------------- parts
def run_coarse():
    """GB18.1-coarse (conforming M0.84 proof-of-concept) + the LS A/C coarse
    ceiling probe (= B26 G1 coarse legs, C now on the flat sheet) + the
    coarse transonic cross-model point. All gated (heavy)."""
    mc, wc = conf_mesh("coarse")
    mesh_a, mvop_a = ls_mesh("coarse", clip=False)
    mesh_c, mvop_c = ls_mesh("coarse", clip=True)
    if mc is None or mesh_a is None or mesh_c is None:
        print("[skip] coarse: wing-body mesh(es) absent")
        return None
    print("=== coarse ===", flush=True)
    # conforming coarse M0.84 (headline proof-of-concept)
    clp84, clkj84, m84, mmax84, r84, _, _ = conf_ramp("coarse", mc, wc, 0.84, 0.70)
    checks.add("GB18.1", "conforming_coarse_M084",
               f"reached={r84} cl_p={clp84:.4f} (M_max {mmax84:.2f})",
               "conforming Mach-ramp reaches M0.84 at coarse (proof-of-concept; "
               "under-resolved)", r84)
    # LS coarse ceiling probe, A (default clip = pocket) vs C (flat fragment)
    ls_a = ls_ramp("A", "coarse", mesh_a, mvop_a, 0.84, 0.50)
    ls_c = ls_ramp("C", "coarse", mesh_c, mvop_c, 0.84, 0.50)
    checks.add("GB18.2", "ls_coarse_ceiling_A_vs_C",
               f"A: m_last={ls_a['m_last']:.4g} dies {ls_a['m_final']} ({ls_a['cls']}); "
               f"C: reached={ls_c['reached']} m_last={ls_c['m_last']:.4g} "
               f"cl_p={ls_c['clp']:.4f}",
               "C coarse REACHES M0.84 while A dies below it (junction-strip "
               "pocket, class (b)) -- the pocket was the coarse ceiling limiter",
               bool(ls_c["reached"]) and ls_a["m_last"] < 0.84)
    # coarse transonic cross-model point at a fixed safe common Mach (0.60):
    # both paths reach it at coarse (under-resolved); the LS side is the
    # production flat-fragment C side (B29).
    mx = 0.60
    conf_cx = conf_ramp("coarse", mc, wc, mx, 0.50)[0]
    ls_res = ls_ramp("C", "coarse", mesh_c, mvop_c, mx, 0.50)
    ls_cx = ls_res["clp"] if ls_res["reached"] else None
    return dict(clp84=clp84, ls_a=ls_a, ls_c=ls_c, mx=mx,
                conf_cx=conf_cx, ls_cx=ls_cx, mc=mc)


def run_medium(coarse):
    mc, wc = conf_mesh("medium")
    mesh_a, mvop_a = ls_mesh("medium", clip=False)
    mesh_c, mvop_c = ls_mesh("medium", clip=True)
    if mc is None or mesh_a is None or mesh_c is None:
        print("[skip] medium: wing-body mesh(es) absent")
        return
    print("=== medium ===", flush=True)
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    # conforming medium cl(M), B32 production recipe (pressure + tip taper):
    # 0.50 re-solved LIVE with the taper (the B9 no-taper anchor above is
    # historical), 0.65, 0.75, 0.79 -- plus the ceiling-climb record leg.
    clp50, _, _, _, r50, _, _ = conf_ramp("medium", mc, wc, 0.50, 0.50)
    clp65, _, _, _, r65, _, _ = conf_ramp("medium", mc, wc, 0.65, 0.60)
    clp75, _, _, _, r75, _, _ = conf_ramp("medium", mc, wc, 0.75, 0.70)
    clp79, clkj79, m79, mmax79, r79, nlim79, nflr79 = conf_ramp(
        "medium", mc, wc, 0.79, 0.70)
    # B32 ceiling-climb record: the B31 dying level (0.83) is cured by the
    # production taper on the analysis grid; the demo ramp records how far
    # the production chain actually climbs (0.80/0.84 are the strict rungs).
    clp_top, _, m_top, mmax_top, r_top, nlim_top, nflr_top = conf_ramp(
        "medium", mc, wc, 0.84, 0.70)
    conf_pts = [(0.50, clp50), (0.65, clp65), (0.75, clp75), (0.79, clp79),
                (m_top, clp_top)]
    checks.add("GB18.1", "conforming_medium_M079",
               f"reached={r79} cl_p={clp79:.4f}; cl(M) 0.50/0.65/0.75/0.79 = "
               f"{clp50:.4f}/{clp65:.4f}/{clp75:.4f}/{clp79:.4f}; B32 taper climb: "
               f"m_reached={m_top:.4g} (cl_p {clp_top:.4f}, clamps "
               f"{nlim_top}+{nflr_top})",
               "conforming production (pressure+taper, B32) reaches M0.79 STRICT "
               "at medium and the climb records the post-adoption ceiling "
               "(B31: the 0.83 dying level is cured on the analysis grid); "
               "monotone transonic rise",
               r79 and r75 and clp79 > clp75 > clp65 > clp50 and m_top >= 0.80)

    # LS medium ceiling probe, A vs C (B26 G1 legs recipe, C now flat-sheet)
    ls_a = ls_ramp("A", "medium", mesh_a, mvop_a, 0.84, 0.50)
    ls_c = ls_ramp("C", "medium", mesh_c, mvop_c, 0.84, 0.50)
    checks.add("GB18.2", "ls_medium_ceiling_A_vs_C",
               f"A: m_last={ls_a['m_last']:.4g} dies {ls_a['m_final']} ({ls_a['cls']}, "
               f"Mmax {ls_a['mmax']:.2f}); C: m_last={ls_c['m_last']:.4g} dies "
               f"{ls_c['m_final']:.4g} ({ls_c['cls']})",
               "C medium climbs past 0.70 while A dies at/below 0.55 "
               "(class (a): the pocket erupts at 0.55, Mmax 13.1 > "
               "freeze_max_clamped=8) -- the pocket was the medium ceiling limiter",
               ls_c["m_last"] >= 0.70 and ls_a["m_last"] <= 0.55)

    # cross-model transonic legs at medium: the LS flat-frag side uses
    # DIRECTED ramps (strict final level) at 0.65 and 0.75.
    ls_c65 = ls_ramp("C", "medium", mesh_c, mvop_c, 0.65, 0.50)
    ls_c75 = ls_ramp("C", "medium", mesh_c, mvop_c, 0.75, 0.50)
    gap65 = abs(ls_c65["clp"] - clp65) / clp65 * 100
    gap75 = abs(ls_c75["clp"] - clp75) / clp75 * 100
    d_conf_05 = abs(clp50 - CL_M05["level_set"]["medium"]) / clp50 * 100
    checks.add("GB18.3", "cross_model_medium_M065",
               f"conf {clp65:.4f} vs LS flat-frag {ls_c65['clp']:.4f} -> |gap| {gap65:.1f}%",
               "PASS gate: |cl_p gap| <= 5% at M0.65 medium (the B9/B28 M0.5 "
               "0.5% reference, relaxed for transonic)", gap65 <= 5.0)
    checks.add("GB18.3", "cross_model_medium_M075",
               f"conf {clp75:.4f} vs LS flat-frag {ls_c75['clp']:.4f} -> |gap| {gap75:.1f}%",
               "RECORDED: M0.75 medium cross-model (shock-sensitive; no "
               "threshold)", True)

    # cross_model.csv: M0.5 anchor + coarse 0.60 + the medium points
    rows = [(0.50, "medium", f"{clp50:.4f}", f"{CL_M05['level_set']['medium']}",
             f"{d_conf_05:.1f}", "trustworthy (B32 taper conf / B28 flat LS)")]
    if coarse and coarse["conf_cx"] is not None and coarse["ls_cx"] is not None:
        mx = coarse["mx"]
        dconf = coarse["conf_cx"] - CL_M05["conforming"]["coarse"]
        dls = coarse["ls_cx"] - CL_M05["level_set"]["coarse"]
        inc = abs(dls - dconf) / abs(dconf) * 100 if abs(dconf) > 1e-6 else float("nan")
        rows.append((mx, "coarse", f"{coarse['conf_cx']:.4f}", f"{coarse['ls_cx']:.4f}",
                     f"{abs(coarse['conf_cx']-coarse['ls_cx'])/coarse['conf_cx']*100:.1f}",
                     f"transonic; increment {inc:.0f}%; under-resolved"))
    rows.append((0.65, "medium", f"{clp65:.4f}", f"{ls_c65['clp']:.4f}",
                 f"{gap65:.1f}", "gated <=5%"))
    rows.append((0.75, "medium", f"{clp75:.4f}", f"{ls_c75['clp']:.4f}",
                 f"{gap75:.1f}", "recorded (shock-sensitive)"))
    write_csv(OUT, "cross_model.csv", "mach,resolution,conf,ls,gap_pct,note", rows)

    # GB18.4 pocket attribution -- the A side cites the committed B26
    # evidence (the demo does not re-run the b23 topk machine); the C-side
    # dying peak is measured LIVE on the flat sheet (B29).
    d = np.load(OUT / "ls_flat_medium_084.npz", allow_pickle=True)
    m2 = mvop_c.element_mach2(d["phi_ext"], float(d["m_final"]), 1.4, 1.0)
    i = int(np.argmax(m2))
    ctr = mesh_c.nodes[mesh_c.elements[i]].mean(axis=0)
    site = ("wing TIP (P13 class)" if ctr[2] > 0.8 * B_SEMI else
            "junction strip" if ctr[0] > X_BAND0 else "mid-wing")
    checks.add("GB18.4", "pocket_attribution",
               f"A dying peak M6.17 @ junction strip x=2.12 (q~0; +M3.53 ON the "
               f"junction, dist_fus=0.005; committed B26 g1_peaks.csv) vs C(flat) "
               f"dying peak M{np.sqrt(m2[i]):.2f} @ {site} (x={ctr[0]:.2f}, "
               f"z={ctr[2]:.2f})",
               "RECORDED: the pocket = the B23 inboard free-edge singularity, "
               "healed on the C side; the residual limiter is the wing-tip P13 "
               "class + high-M Newton -- same class as the conforming stall. "
               "Flat-side junction corridor at M0.5 measured clean in B28 "
               "(corrM 0.64, n_sup 0). A-side source: committed "
               "cases/analysis/b26_ls_transonic_ceiling/results/g1_summary.csv "
               "+ g1_peaks.csv", True)

    # GB18.5 fuselage lift: live conforming value at the medium transonic top
    # + the LIVE flat C-side ceiling decomposition (B29; the B26 tilted-sheet
    # "x2 out-band" reading is retired by B28 = sheet-position sensitivity).
    cl_fus = float(wall_forces(mc, phi=np.load(OUT / "conf_taper_medium_079.npz")["phi"],
                               alpha_deg=ALPHA, s_ref=s_ref, m_inf=m79, wall_tag="fuselage")["cl"])
    parts = fuselage_parts(mesh_c, mvop_c, d["phi_ext"], s_ref, float(d["m_final"]))
    checks.add("GB18.5", "fuselage_lift",
               f"conf cl_fus {cl_fus:.4f} = {abs(cl_fus)/clp79*100:.0f}% of wing cl_p "
               f"{clp79:.4f} @M0.79; LS flat-frag C medium cl_fus "
               f"{parts['band']+parts['out']:.4f} (band {parts['band']:.4f}, "
               f"out-band {parts['out']:.4f}, poles {parts['poles']:.4f}) "
               f"@{float(d['m_final']):.4g} ceiling state",
               "RECORDED: fuselage lift = physical carryover + wake-sheet "
               "POSITION sensitivity (B28), NOT the G1.6 error; the B26 "
               "tilted-sheet 'x2 out-band' reading (P11 watch item) is "
               "retired. The GB9.4 gate lives in the B9 demo (out-band "
               "cross-model <= 15%)", True)

    # artifacts
    write_csv(OUT, "cl_vs_mach.csv", "mach,path,resolution,cl_p,note",
              [(0.50, "conforming", "medium", f"{clp50:.4f}", "B32 taper production anchor"),
               (0.65, "conforming", "medium", f"{clp65:.4f}", ""),
               (0.75, "conforming", "medium", f"{clp75:.4f}", ""),
               (0.79, "conforming", "medium", f"{clp79:.4f}", "strict"),
               (f"{m_top:.4g}", "conforming", "medium", f"{clp_top:.4f}",
                f"B32 climb record (clamps {nlim_top}+{nflr_top})"),
               (0.84, "conforming", "coarse", f"{coarse['clp84']:.4f}" if coarse else "n/a",
                "proof-of-concept"),
               (0.50, "level_set+flatfrag", "medium", f"{CL_M05['level_set']['medium']}", "B28 flat anchor"),
               (0.65, "level_set+flatfrag", "medium", f"{ls_c65['clp']:.4f}", "directed ramp; strict final"),
               (0.75, "level_set+flatfrag", "medium", f"{ls_c75['clp']:.4f}", "directed ramp; strict final"),
               (f"{ls_c['m_last']:.4g}", "level_set+flatfrag", "medium", f"{ls_c['clp']:.4f}",
                f"ceiling m_last={ls_c['m_last']:.4g}; cl_p at the {ls_c['m_final']:.4g} dying state"),
               (0.84, "level_set+flatfrag", "coarse",
                f"{coarse['ls_c']['clp']:.4f}" if coarse else "n/a", "reached")])
    ls_pts = [(0.50, CL_M05["level_set"]["medium"]),
              (0.65, ls_c65["clp"]), (0.75, ls_c75["clp"])]
    a_deaths = [(ls_a["m_final"], f"A dies {ls_a['m_final']:.4g}")]
    if coarse:
        a_deaths.append((coarse["ls_a"]["m_final"], f"A dies {coarse['ls_a']['m_final']:.4g}"))
    fig_cl_vs_mach(conf_pts, ls_pts, coarse["clp84"] if coarse else 0.2617,
                   coarse["ls_c"]["clp"] if coarse else 0.2542,
                   (ls_c["m_last"], ls_c["clp"]), a_deaths)
    fig_ceiling(coarse["ls_a"]["m_last"] if coarse else 0.82,
                coarse["ls_c"]["m_last"] if coarse else 0.84,
                ls_a["m_last"], ls_c["m_last"], conf_medium=m_top)
    fig_sections_medium(mc)


def main():
    apply_style()
    if not GATED:
        print("B18 is fully gated (heavy transonic ramps). "
              "Set PYFP3D_TRANSONIC_GATES=1 to run.")
        # a tiny ungated wiring check: both transonic drivers present (the
        # conforming target is `m_inf`, the level-set target is `m_target`)
        # and the B25 inboard clip + the B28 sheet_direction knob are
        # importable.
        import inspect
        assert "m_inf" in inspect.signature(solve_newton_transonic).parameters
        assert "m_target" in inspect.signature(solve_multivalued_newton_transonic).parameters
        assert "sheet_direction" in inspect.signature(WakeLevelSet.__init__).parameters
        clip = make_inboard_clip(FUS)
        assert callable(clip)
        print("[wiring ok] both transonic drivers + the inboard clip + the "
              "flat-sheet knob present")
        sys.exit(0)
    coarse = run_coarse()
    run_medium(coarse)
    sys.exit(checks.report(OUT, "checks.csv"))


if __name__ == "__main__":
    main()
