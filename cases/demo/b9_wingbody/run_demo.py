"""
Track B / B9 (re-spec 2026-07-17): wing-body CROSS-MODEL comparison at M0.5.

Runs BOTH wake models on the ONERA M6 wing-body (coarse + medium, alpha 3.06,
the committed M6 subsonic convention) and compares the four things the phase
was re-spec'd around:

  1. spanwise lift            Gamma(z) and sectional cl(z)
  2. section Cp               same-extractor sweep, both models (A2/V14.6 rule)
  3. convergence histories    residual vs iteration and vs wall-clock
  4. timing breakdown         the A1 phase schema (seed/assembly/precond/
                              linsolve/residual/kutta)

Solvers -- the best-known recipe per path (measured 2026-07-17):
  * CONFORMING: coupled Newton, probe estimator -> probe-seeded PRESSURE
    estimator (P14). Mesh = the NEW wake-embedded variant
    onera_m6_wingbody_conforming (cut_wake + P14 pressure Kutta, all unchanged).
  * LEVEL-SET: Picard (solve_multivalued_lifting), farfield="freestream".
    Two measured findings drive this: (a) farfield="neumann" (the Lopez
    inlet-Dirichlet/outlet-Neumann outlet, fine for a thin wing) DIVERGES on
    the wing-body -- the fuselage blockage makes the outlet flux balance
    unbounded (res ~ 1e43); the 25-MAC domain makes full-freestream Dirichlet
    accurate. (b) the subsonic LS NEWTON oscillates on the wing-body (limiter
    fires at M0.5, drifts off the converged Picard state), so the LS path uses
    its proven subsonic solver, Picard (B7); the LS Newton is the transonic
    ramp tool (B15). This is itself a "best practice from all tracks" outcome.

Cross-model discipline (A2/V14.6): cl_KJ uses a demo-local EXPOSED-SPAN reducer
(trapezoid over the actual junction..tip stations + tip closure, NO root
flat-extension) applied identically to BOTH paths; section Cp uses the shared
section extractors, wing 'wall' only. Fuselage lift is reported separately.

The heavy solves cache to gitignored results/*.npz (the committed evidence is
the figures + CSVs). Meshes are gitignored -> the demo skips a level whose mesh
is absent. Cost: ~2 min coarse, ~10-20 min medium (cold).

Run:  python cases/demo/b9_wingbody/run_demo.py [--levels coarse medium]
"""

import argparse
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

from cases.demo._common import CheckList, apply_style, finish, write_csv
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.fuselage import FuselageParams
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at, x_te
from pyfp3d.meshgen.wingbody import junction_z, te_polyline
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.surface import planform_area, sectional_cl_from_gamma
from pyfp3d.post.surface_ls import section_cp_curve_levelset
from pyfp3d.post.unified import wall_forces
from pyfp3d.solve.newton import solve_newton_lifting
from pyfp3d.solve.picard_ls import solve_multivalued_lifting
from pyfp3d.solve.timing import PHASES
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA, M = 3.06, 0.5
FUS = FuselageParams()
Z_JUNC = junction_z(FUS)
CONF_DIR = REPO_ROOT / "cases/meshes/onera_m6_wingbody_conforming"
LS_DIR = REPO_ROOT / "cases/meshes/onera_m6_wingbody"
# near-junction, mid, outer stations (eta_junction = 0.15/1.1963 = 0.125)
ETAS = (0.20, 0.44, 0.65, 0.90)
CONF_KW = dict(farfield_spanwise_gamma=True, precond="direct",
               direct_refactor_every=1000, n_newton_max=60)

checks = CheckList("B9 wing-body cross-model (LS + conforming), M0.5")


# --------------------------------------------------------------------------
# exposed-span cl_KJ reducer -- SAME for both paths (no root flat-extension)
# --------------------------------------------------------------------------
def cl_kj_exposed(z, gamma, s_ref):
    """CL_KJ = 2/S * integral Gamma(z) dz over the ACTUAL stations
    (junction..tip) + tip closure to B_SEMI (Gamma -> 0). No z=0 root
    extension: the wing-body carries no bound circulation across the body,
    so cl_kj_3d's flat carry-through strip [0, r_f] is a spurious add."""
    o = np.argsort(z)
    zz, gg = z[o], gamma[o]
    zz = np.concatenate([zz, [B_SEMI]])
    gg = np.concatenate([gg, [0.0]])
    return 2.0 * float(np.trapezoid(gg, zz)) / s_ref


# --------------------------------------------------------------------------
# solves (cached)
# --------------------------------------------------------------------------
def solve_conforming(level):
    """probe Newton -> probe-seeded pressure Newton (P14). Returns a dict of
    both estimators' states + metrics + histories + timings."""
    cache = OUT / f"conf_{level}.npz"
    mp = CONF_DIR / f"{level}.msh"
    if not mp.exists():
        return None
    mc, wc = cut_wake(read_mesh(str(mp)))
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    z = mc.nodes[wc.te_nodes, 2]
    rec = {"s_ref": s_ref, "z": z}
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        for est in ("probe", "pressure"):
            rec[est] = dict(phi=d[f"{est}_phi"], gamma=d[f"{est}_gamma"],
                            res=list(d[f"{est}_res"]), wall_s=float(d[f"{est}_wall"]),
                            timings={p: float(d[f"{est}_t_{p}"]) for p in PHASES},
                            n=int(d[f"{est}_n"]), conv=bool(d[f"{est}_conv"]),
                            nlim=int(d[f"{est}_nlim"]), nflr=int(d[f"{est}_nflr"]))
        _finish_conf(rec, mc, wc)
        return rec
    save = {}
    for est in ("probe", "pressure"):
        seed = ({} if est == "probe" else
                dict(phi_init=rec["probe"]["phi"], gamma_init=rec["probe"]["gamma"],
                     n_picard_seed=0))
        t0 = time.perf_counter()
        r = solve_newton_lifting(mc, wc, m_inf=M, alpha_deg=ALPHA,
                                 kutta_estimator=est, **seed, **CONF_KW)
        wall_s = time.perf_counter() - t0
        rec[est] = dict(phi=r["phi"], gamma=np.asarray(r["gamma"]),
                        res=list(r["residual_history"]), wall_s=wall_s,
                        timings=r["timings"], n=r["n_newton"],
                        conv=bool(r["converged"]), nlim=int(r["n_limited"]),
                        nflr=int(r["n_floored"]))
        save[f"{est}_phi"] = r["phi"]; save[f"{est}_gamma"] = np.asarray(r["gamma"])
        save[f"{est}_res"] = np.asarray(r["residual_history"])
        save[f"{est}_wall"] = wall_s; save[f"{est}_n"] = r["n_newton"]
        save[f"{est}_conv"] = r["converged"]; save[f"{est}_nlim"] = r["n_limited"]
        save[f"{est}_nflr"] = r["n_floored"]
        for p in PHASES:
            save[f"{est}_t_{p}"] = r["timings"].get(p, 0.0)
        print(f"  [conf {level} {est}] conv={r['converged']} n={r['n_newton']} "
              f"lim={r['n_limited']} flr={r['n_floored']} ({wall_s:.0f}s)", flush=True)
    np.savez(cache, **save)
    _finish_conf(rec, mc, wc)
    return rec


def _finish_conf(rec, mc, wc):
    for est in ("probe", "pressure"):
        r = rec[est]
        clp = wall_forces(mc, phi=r["phi"], alpha_deg=ALPHA, s_ref=rec["s_ref"],
                          m_inf=M, wall_tag="wall")["cl"]
        clf = wall_forces(mc, phi=r["phi"], alpha_deg=ALPHA, s_ref=rec["s_ref"],
                          m_inf=M, wall_tag="fuselage")["cl"]
        clkj = cl_kj_exposed(rec["z"], r["gamma"], rec["s_ref"])
        r.update(cl_p=float(clp), cl_fus=float(clf), cl_kj=float(clkj))
    rec["mc"], rec["wc"] = mc, wc


def solve_ls(level):
    """LS Picard, farfield=freestream (the measured wing-body recipe)."""
    cache = OUT / f"ls_{level}.npz"
    mp = LS_DIR / f"{level}.msh"
    if not mp.exists():
        return None
    mesh = read_mesh(str(mp))
    a = np.radians(ALPHA)
    wls = WakeLevelSet(te_polyline(FUS), direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    z = mesh.nodes[cm.te_nodes, 2]
    rec = {"s_ref": s_ref, "z": z, "mesh": mesh, "mvop": mvop, "cm": cm}
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        rec.update(phi_ext=d["phi_ext"], res=list(d["res"]),
                   wall_s=float(d["wall"]), n=int(d["n"]), conv=bool(d["conv"]),
                   timings={p: float(d[f"t_{p}"]) for p in PHASES})
    else:
        t0 = time.perf_counter()
        r = solve_multivalued_lifting(mvop, mesh, M, alpha_deg=ALPHA,
                                      farfield="freestream", farfield_aux="legacy",
                                      n_outer_max=80, tol_residual=1e-7)
        # B17 erratum: farfield_aux="legacy" (free aux) pins this to B9's
        # COMMITTED behaviour. B17 found the far-field aux pin matters here --
        # legacy leaves them free (the wake-LS carries the jump: at medium
        # jump~0.105~Gamma -> cl_p 0.2165~conforming; at COARSE the giant outer
        # tet is near-singular -> garbage |jump|=53 -> cl_p 0.1853, so B9's
        # "coarse 12.8% = resolution" was largely far-field CONTAMINATION, not
        # resolution). farfield_aux="pin_gamma" (the new default) carries jump
        # =Gamma cleanly at both (coarse 0.2087, medium 0.2117). See B17.
        wall_s = time.perf_counter() - t0
        rec.update(phi_ext=r["phi_ext"], res=list(r["residual_history"]),
                   wall_s=wall_s, n=r["n_outer"], conv=bool(r["converged"]),
                   timings=r["timings"])
        save = dict(phi_ext=r["phi_ext"], res=np.asarray(r["residual_history"]),
                    wall=wall_s, n=r["n_outer"], conv=r["converged"])
        for p in PHASES:
            save[f"t_{p}"] = r["timings"].get(p, 0.0)
        np.savez(cache, **save)
        print(f"  [ls {level}] conv={r['converged']} n_outer={r['n_outer']} "
              f"({wall_s:.0f}s)", flush=True)
    g = mvop.te_jump(rec["phi_ext"])
    rec["gamma"] = g
    rec["cl_p"] = float(wall_forces(mesh, mvop=mvop, phi_ext=rec["phi_ext"],
                        alpha_deg=ALPHA, s_ref=s_ref, m_inf=M, wall_tag="wall")["cl"])
    rec["cl_fus"] = float(wall_forces(mesh, mvop=mvop, phi_ext=rec["phi_ext"],
                          alpha_deg=ALPHA, s_ref=s_ref, m_inf=M, wall_tag="fuselage")["cl"])
    rec["cl_kj"] = cl_kj_exposed(z, g, s_ref)
    return rec


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------
def fig_spanwise(conf, ls, level):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for rec, lab, col, gk in ((conf["pressure"], "conforming (pressure)", "tab:blue", "gamma"),
                              (ls, "level-set (Picard)", "tab:red", "gamma")):
        z = conf["z"] if lab.startswith("conf") else ls["z"]
        g = rec[gk] if isinstance(rec, dict) else rec
        o = np.argsort(z)
        ax[0].plot(z[o], np.asarray(g)[o], "o-", ms=3, color=col, label=lab)
        c = np.array([chord_at(float(zz)) for zz in z[o]])
        ax[1].plot(z[o], sectional_cl_from_gamma(np.asarray(g)[o], chord=c),
                   "o-", ms=3, color=col, label=lab)
    for a in ax:
        a.axvline(Z_JUNC, color="0.5", ls="--", lw=0.8, label=f"junction z={Z_JUNC}")
        a.set_xlabel("z (span)"); a.grid(alpha=0.3); a.legend(fontsize=8)
    ax[0].set_ylabel("Gamma(z)"); ax[0].set_title(f"{level}: spanwise circulation")
    ax[1].set_ylabel("sectional cl(z)"); ax[1].set_title(f"{level}: sectional lift")
    fig.suptitle("B9 spanwise lift: conforming-pressure vs level-set (M0.5)")
    finish(fig, OUT, f"b9_spanwise_{level}.png")


def fig_sections(conf, ls, level):
    fig, ax = plt.subplots(1, len(ETAS), figsize=(4.2 * len(ETAS), 4.2),
                           squeeze=False)
    mc, wc = conf["mc"], conf["wc"]
    phi = conf["pressure"]["phi"]
    def _plot(a, sc, col, lab):
        c = sc["chord"] or 1.0
        xl = (sc["x_le"] or 0.0)
        xu = (np.asarray(sc["x_upper"]) - xl) / c
        xlo = (np.asarray(sc["x_lower"]) - xl) / c
        a.plot(xu, sc["cp_upper"], "-", lw=1.0, color=col, label=lab)
        a.plot(xlo, sc["cp_lower"], "-", lw=1.0, color=col, alpha=0.6)

    for j, eta in enumerate(ETAS):
        a = ax[0][j]
        try:
            sc = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M,
                                  wall_tag="wall")
            _plot(a, sc, "tab:blue", "conforming")
        except Exception as e:
            a.text(0.5, 0.5, f"conf n/a\n{type(e).__name__}", fontsize=7,
                   ha="center", transform=a.transAxes)
        try:
            sl = section_cp_curve_levelset(ls["mesh"], ls["mvop"], ls["phi_ext"],
                                           eta=eta, b_semi=B_SEMI, m_inf=M,
                                           wall_tag="wall")
            _plot(a, sl, "tab:red", "level-set")
        except Exception as e:
            a.text(0.5, 0.4, f"LS n/a\n{type(e).__name__}", fontsize=7,
                   ha="center", transform=a.transAxes)
        a.invert_yaxis(); a.set_xlabel("x/c"); a.grid(alpha=0.3)
        a.set_title(f"eta={eta} (z={eta*B_SEMI:.2f})")
        if j == 0:
            a.set_ylabel("Cp"); a.legend(fontsize=8)
    fig.suptitle(f"B9 section Cp, {level}, M0.5 (same extractor, wing wall only)")
    finish(fig, OUT, f"b9_sections_{level}.png")


def fig_convergence(conf, ls, level):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    series = [("conf probe", conf["probe"]["res"], conf["probe"]["wall_s"], "tab:cyan"),
              ("conf pressure", conf["pressure"]["res"], conf["pressure"]["wall_s"], "tab:blue"),
              ("level-set", ls["res"], ls["wall_s"], "tab:red")]
    for lab, res, wall, col in series:
        res = np.asarray(res)
        it = np.arange(len(res))
        ax[0].semilogy(it, res, "o-", ms=3, color=col, label=lab)
        tt = np.linspace(0, wall, len(res)) if len(res) > 1 else [wall]
        ax[1].semilogy(tt, res, "o-", ms=3, color=col, label=lab)
    ax[0].set_xlabel("iteration"); ax[1].set_xlabel("wall-clock (s)")
    for a in ax:
        a.set_ylabel("residual (inf-norm)"); a.grid(alpha=0.3, which="both")
        a.legend(fontsize=8)
    ax[0].set_title(f"{level}: residual vs iteration")
    ax[1].set_title(f"{level}: residual vs wall-clock")
    fig.suptitle("B9 convergence: conforming Newton vs level-set Picard (M0.5)")
    finish(fig, OUT, f"b9_convergence_{level}.png")


def fig_timing(rows, level):
    fig, ax = plt.subplots(figsize=(9, 4.8))
    labels = [r[0] for r in rows]
    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    cmap = plt.get_cmap("tab10")
    for i, p in enumerate(PHASES + ("other",)):
        vals = np.array([r[1].get(p, 0.0) for r in rows])
        ax.bar(x, vals, bottom=bottom, label=p, color=cmap(i))
        bottom += vals
    for xi, r in zip(x, rows):
        ax.text(xi, bottom[xi], f"{bottom[xi]:.0f}s", ha="center", va="bottom",
                fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("wall-clock (s)"); ax.legend(fontsize=8, ncol=2)
    ax.set_title(f"B9 timing breakdown ({level}, A1 phase schema)")
    finish(fig, OUT, f"b9_time_breakdown_{level}.png")


# --------------------------------------------------------------------------
def run_level(level):
    print(f"=== {level} ===", flush=True)
    conf = solve_conforming(level)
    ls = solve_ls(level)
    if conf is None or ls is None:
        print(f"[skip] {level}: "
              f"{'conforming ' if conf is None else ''}"
              f"{'level-set ' if ls is None else ''}mesh absent")
        return None

    fig_spanwise(conf, ls, level)
    fig_sections(conf, ls, level)
    fig_convergence(conf, ls, level)
    trows = [("conf probe", conf["probe"]["timings"]),
             ("conf pressure", conf["pressure"]["timings"]),
             ("level-set", ls["timings"])]
    fig_timing(trows, level)

    cp, lp = conf["pressure"], ls
    xm_p = abs(cp["cl_p"] - lp["cl_p"]) / abs(lp["cl_p"]) * 100
    xm_kj = abs(cp["cl_kj"] - lp["cl_kj"]) / abs(lp["cl_kj"]) * 100
    print(f"  [{level}] conf cl_p(wing)={cp['cl_p']:.4f} cl_kj={cp['cl_kj']:.4f} "
          f"cl_fus/wing={abs(cp['cl_fus']/cp['cl_p']):.3f} | "
          f"LS cl_p={lp['cl_p']:.4f} cl_kj={lp['cl_kj']:.4f} "
          f"cl_fus/wing={abs(lp['cl_fus']/lp['cl_p']):.3f} | "
          f"cross-model cl_p {xm_p:.1f}% cl_kj {xm_kj:.1f}%", flush=True)

    # gates
    gated = (level == "medium")
    checks.add("GB9.1", f"{level}_conforming_converged",
               f"probe {conf['probe']['n']} / pressure {cp['n']} steps",
               "both Newton solves converge, 0 lim/flr",
               conf["probe"]["conv"] and cp["conv"]
               and cp["nlim"] == 0 and cp["nflr"] == 0)
    checks.add("GB9.2", f"{level}_ls_converged",
               f"Picard {lp['n']} outer, res in cache",
               "LS Picard converges (freestream)", lp["conv"])
    # GB9.4 is an XFAIL (documented-open, like G1.6): the pre-registered
    # "fuselage carries no lift <= 5% cl_wing" band FAILS -- the fuselage
    # pressure-lift is 16% (conf) / 20% (LS) of the wing's, and the LS value
    # GROWS with refinement (0.164 -> 0.205) while conforming stays flat.
    # That resolution- and model-sensitivity is the signature of the G1.6
    # smooth-wall flat-facet natural-BC error on the fuselage (the GB9.6
    # guardrail's whole subject), not clean physics. The band is NOT moved
    # after the fact (house rule); it is recorded as the open state a
    # wing-body body-surface claim carries until P11/Option C.
    gb94_pass = (abs(cp["cl_fus"]) <= 0.05 * abs(cp["cl_p"])
                 and abs(lp["cl_fus"]) <= 0.05 * abs(lp["cl_p"]))
    checks.add("GB9.4", f"{level}_fuselage_lift",
               f"conf {abs(cp['cl_fus']/cp['cl_p']):.3f} / "
               f"LS {abs(lp['cl_fus']/lp['cl_p']):.3f}",
               "|cl_fus| <= 0.05 cl_p_wing" if gated else "RECORDED (coarse)",
               gb94_pass if gated else True, xfail=gated,
               note=("fuselage lift resolution/model-sensitive -> G1.6 "
                     "fuselage-Cp error (GB9.6); band NOT moved") if gated else "")
    checks.add("GB9.5", f"{level}_cross_model",
               f"cl_p {xm_p:.1f}% cl_kj {xm_kj:.1f}%",
               "conf-pressure vs LS < 1% (gated at medium)" if gated
               else "RECORDED (coarse)",
               (xm_p < 1.0 and xm_kj < 1.0) if gated else True)

    return dict(level=level, conf=conf, ls=ls, xm_p=xm_p, xm_kj=xm_kj)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=["coarse", "medium"])
    args = ap.parse_args()
    apply_style()

    results = [r for r in (run_level(lv) for lv in args.levels) if r]
    if not results:
        print("no levels ran (meshes absent); nothing to check")
        sys.exit(0)

    # summary + cross-model CSVs
    srows = []
    for r in results:
        for path, rec, est in (("conf_probe", r["conf"]["probe"], None),
                               ("conf_pressure", r["conf"]["pressure"], None),
                               ("level_set", r["ls"], None)):
            srows.append((
                r["level"], path,
                f"{rec.get('conv', rec.get('conv', '')):}",
                rec.get("n", ""),
                f"{rec['cl_p']:.4f}", f"{rec['cl_fus']:.4f}", f"{rec['cl_kj']:.4f}",
                rec.get("nlim", ""), rec.get("nflr", ""),
                f"{rec['wall_s']:.1f}",
                *[f"{rec['timings'].get(p, 0.0):.1f}" for p in PHASES]))
    write_csv(OUT, "b9_summary.csv",
              "level,path,converged,n_iter,cl_p_wing,cl_p_fus,cl_kj_exposed,"
              "n_limited,n_floored,wall_s," + ",".join(f"t_{p}" for p in PHASES),
              srows)

    cm_rows = []
    for r in results:
        for path, rec in (("conforming_pressure", r["conf"]["pressure"]),
                          ("level_set", r["ls"])):
            cm_rows.append((r["level"], path,
                            f"{rec['cl_p']:.4f}", f"{rec['cl_kj']:.4f}",
                            f"{rec['cl_fus']:.4f}"))
    write_csv(OUT, "cross_model_m05.csv",
              "level,path,cl_p_wing,cl_kj_exposed,cl_p_fus", cm_rows)

    sys.exit(checks.report(OUT, "checks.csv"))


if __name__ == "__main__":
    main()
