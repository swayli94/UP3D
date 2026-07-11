"""
P9 evidence demo: grid-convergence & accuracy-gap discrimination
(roadmap P9, gates G9.1/G9.2/G9.3 -- decision bands PRE-REGISTERED in
roadmap.md before any fine-mesh number existed; this script only
evaluates them).

  G9.1  ONERA M6 coarse/medium/fine Newton at M0.84/alpha3.06
        (NEWTON_M6_RECIPE + the G10.2 intermediate tolerance) +
        Richardson extrapolation of cl_KJ and cl_p with the measured
        per-level h-ratio from tet counts. Verdict bands on cl_KJ_inf:
        >= 0.283 resolution-dominated / <= 0.278 floor confirmed /
        else inconclusive.
  G9.2  NACA0012 2.5D coarse/medium/fine subsonic (M0.5/alpha2) Newton
        lift oracle vs the corrected-panel PG-KT midpoint: PASS iff
        |error| decreases monotonically AND fine within +-1%.
  G9.3  attribution split of the 0.019 external gap (cl_KJ 0.2692 vs
        Tranair/KRATOS 0.288) into resolution / floor / unattributed
        shares -- written to verdict.csv; the user arbitrates the P11
        go/no-go.

Heavy-run protocol (roadmap P9 non-goals): fine meshes stay LOCAL
(gitignored; regenerate with the mesh scripts' --level fine); the M6
fine solve (~450k dofs, the phase's technical unknown) runs ONCE and is
npz-cached in results/ exactly like the P5 medium cache -- the demo
re-solves only when the cache is absent or PYFP3D_P9_RESOLVE=1. Committed
evidence = the CSVs/PNGs. Suite budget untouched (no tests here).

Timing protocol: NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
OPENBLAS_NUM_THREADS=16.
"""

import csv
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import (  # noqa: E402
    BASELINE, CRITICAL, CheckList, INK_2, MESH_DIR, REFERENCE_DIR, S1_BLUE,
    S2_AQUA, S3_YELLOW, apply_style, finish, write_csv,
)
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.wing3d import B_SEMI
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.shock import shock_report
from pyfp3d.post.surface import cl_kj_3d, planform_area, wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting, solve_newton_transonic

OUT = Path(__file__).parent / "results"
OUT.mkdir(exist_ok=True)

M_M6, ALPHA_M6 = 0.84, 3.06
M_SUB, ALPHA_SUB = 0.5, 2.0
CL_EXTERNAL_REF = 0.288          # Tranair/KRATOS (Lopez Table 4.15)
BAND_HI, BAND_LO = 0.283, 0.278  # pre-registered G9.1 verdict bands
LEVELS = ("coarse", "medium", "fine")

# NEWTON_M6_RECIPE (tests/test_p8_newton.py) + the G10.2 level-adaptive
# intermediate tolerance (promoted by the committed A/B in
# cases/demo/p10_newton_usability/ -- see roadmap P10/G10.2)
M6_RECIPE = dict(
    dm=0.05, dm_min=0.01, freeze_tol=1e-6, intermediate_tol=1e-5,
    newton_kw=dict(freeze_refresh_max=8, precond="direct",
                   direct_refactor_every=1000, n_newton_max=60,
                   farfield_spanwise_gamma=True),
)


def require_mesh(geom, level):
    p = MESH_DIR / geom / f"{level}.msh"
    if not p.exists():
        raise SystemExit(
            f"{p} missing -- generate it first (fine meshes are local by "
            f"design): python cases/meshes/{geom}/generate_"
            f"{'onera_m6' if geom == 'onera_m6' else 'naca0012'}.py "
            f"--level {level}")
    return p


# ------------------------------------------------------------------ G9.1

def run_m6_level(level):
    """One M6 Newton continuation run, npz-cached (phi/gamma + walls)."""
    cache = OUT / f"g91_m6_{level}.npz"
    path = require_mesh("onera_m6", level)
    t0 = time.perf_counter()
    mc, wc = cut_wake(read_mesh(path))
    wall_mesh = time.perf_counter() - t0
    n_tets = len(mc.elements)

    if cache.exists() and os.environ.get("PYFP3D_P9_RESOLVE", "0") != "1":
        z = np.load(cache, allow_pickle=True)
        r = dict(z["result_scalars"].item())
        phi, gamma = z["phi"], z["gamma"]
        wall_solve = float(z["wall_solve"])
        print(f"  [m6 {level}] loaded cache {cache.name}")
    else:
        t0 = time.perf_counter()
        res = solve_newton_transonic(mc, wc, m_inf=M_M6, alpha_deg=ALPHA_M6,
                                     verbose=True, **M6_RECIPE)
        wall_solve = time.perf_counter() - t0
        phi, gamma = np.asarray(res["phi"]), np.asarray(res["gamma"])
        r = {
            "converged": bool(res["converged"]),
            "n_limited": int(res["n_limited"]),
            "n_floored": int(res["n_floored"]),
            "residual_final": float(res["residual_history"][-1]),
            "F_final": float(res["F_history"][-1]),
            "mach2_max": float(res["mach2_max"]),
            "residual_unfrozen": (None if res["residual_unfrozen"] is None
                                  else float(res["residual_unfrozen"])),
            "n_levels": len(res["level_results"]),
            "level_walls": [float(lr["wall_s"])
                            for lr in res["level_results"]],
            "level_ms": [float(lr["m"]) for lr in res["level_results"]],
            "level_steps": [int(lr["n_newton"])
                            for lr in res["level_results"]],
            "level_accepts": [str(lr["accept_reason"])
                              for lr in res["level_results"]],
            "quad_drops": [float(d) for d in np.divide(
                res["residual_history"][1:],
                res["residual_history"][:-1])] if len(
                    res["residual_history"]) > 1 else [],
        }
        np.savez_compressed(cache, phi=phi, gamma=gamma,
                            wall_solve=wall_solve,
                            result_scalars=np.array(r, dtype=object))
        print(f"  [m6 {level}] solved in {wall_solve:.0f}s, cached")

    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
        alpha_deg=ALPHA_M6, s_ref=s_ref, m_inf=M_M6)
    shocks = {}
    for eta in (0.44, 0.65, 0.90):
        c = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_M6)
        shocks[eta] = shock_report(c, M_M6)["upper"]["x_shock"]
    return dict(
        level=level, n_tets=n_tets, n_nodes=len(mc.nodes),
        h=float(n_tets) ** (-1.0 / 3.0),
        cl_p=float(forces["cl"]),
        cl_kj=float(cl_kj_3d(gamma, wc.station_z, s_ref=s_ref,
                             b_semi=B_SEMI)),
        m_max=float(np.sqrt(r["mach2_max"])), shocks=shocks,
        wall_mesh=wall_mesh, wall_solve=wall_solve, meta=r,
    )


def solve_order_and_extrapolant(h, f):
    """Fit f_i = f_inf + C h_i^p through three points with UNEQUAL
    refinement ratio (h from tet counts): solve the observed order p by
    bisection on g(p) = (f1-f2)/(f2-f3) - (h1^p-h2^p)/(h2^p-h3^p),
    then f_inf = f3 - (f2-f3) h3^p / (h2^p - h3^p).
    Requires a monotone sequence; returns (p, f_inf) or (None, None)."""
    f1, f2, f3 = f
    h1, h2, h3 = h
    d12, d23 = f1 - f2, f2 - f3
    if d12 == 0.0 or d23 == 0.0 or np.sign(d12) != np.sign(d23):
        return None, None
    target = d12 / d23

    def g(p):
        return (h1 ** p - h2 ** p) / (h2 ** p - h3 ** p) - target

    lo, hi = 1e-3, 10.0
    if g(lo) * g(hi) > 0:
        return None, None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(lo) * g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    p = 0.5 * (lo + hi)
    f_inf = f3 - d23 * h3 ** p / (h2 ** p - h3 ** p)
    return float(p), float(f_inf)


# ------------------------------------------------------------------ G9.2

def load_cl_bracket(alpha_deg):
    with open(REFERENCE_DIR / "naca0012_m05" / "cl_reference.csv") as f:
        for row in csv.DictReader(f):
            if abs(float(row["alpha_deg"]) - alpha_deg) < 1e-9:
                return float(row["cl_pg"]), float(row["cl_kt"])
    raise ValueError(f"alpha {alpha_deg} not in cl_reference.csv")


def run_naca_level(level):
    path = require_mesh("naca0012_2.5d", level)
    mc, wc = cut_wake(read_mesh(path))
    t0 = time.perf_counter()
    r = solve_newton_lifting(mc, wc, m_inf=M_SUB, alpha_deg=ALPHA_SUB,
                             precond="direct")
    wall_solve = time.perf_counter() - t0
    dz = float(np.ptp(mc.nodes[:, 2]))
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
        alpha_deg=ALPHA_SUB, s_ref=dz, m_inf=M_SUB)
    return dict(
        level=level, n_tets=len(mc.elements), n_nodes=len(mc.nodes),
        converged=bool(r["converged"]),
        residual_final=float(r["residual_history"][-1]),
        cl_p=float(forces["cl"]),
        cl_kj=2.0 * float(np.mean(r["gamma"])),
        wall_solve=wall_solve,
    )


# --------------------------------------------------------------- figures

def fig_g91(recs, rich):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    h = np.asarray([r["h"] for r in recs])
    for ax, key, title in ((axes[0], "cl_kj", "M6 cl_KJ vs h"),
                           (axes[1], "cl_p", "M6 cl_p (pressure) vs h")):
        f = np.asarray([r[key] for r in recs])
        ax.plot(h, f, "o-", color=S1_BLUE, label="Newton solutions")
        p, f_inf = rich[key]
        if f_inf is not None:
            ax.axhline(f_inf, color=S2_AQUA, ls="--", lw=1.4,
                       label=f"Richardson $f_\\infty$={f_inf:.4f} "
                             f"(p={p:.2f})")
        if key == "cl_kj":
            ax.axhline(CL_EXTERNAL_REF, color=INK_2, ls=":", lw=1.2,
                       label="Tranair/KRATOS 0.288")
            ax.axhspan(BAND_LO, BAND_HI, color=S3_YELLOW, alpha=0.15,
                       label="pre-registered inconclusive band")
        ax.set_xlabel(r"$h \propto N_{tet}^{-1/3}$")
        ax.set_ylabel(key)
        ax.set_xlim(left=0.0)
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=8)
        for r_, x, y in zip(recs, h, f):
            ax.annotate(r_["level"], (x, y), textcoords="offset points",
                        xytext=(4, -10), fontsize=8, color=INK_2)
    finish(fig, OUT, "g91_m6_grid_convergence.png")


def fig_g92(recs, cl_mid, cl_pg, cl_kt):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    h = np.asarray([float(r["n_tets"]) ** (-1.0 / 3.0) for r in recs])
    f = np.asarray([r["cl_p"] for r in recs])
    axes[0].plot(h, f, "o-", color=S1_BLUE, label="Newton cl_p")
    axes[0].axhline(cl_mid, color=INK_2, ls=":", lw=1.2, label="PG-KT midpoint")
    axes[0].axhspan(cl_pg, cl_kt, color=BASELINE, alpha=0.4,
                    label="PG-KT bracket")
    axes[0].set_xlabel(r"$h \propto N_{tet}^{-1/3}$")
    axes[0].set_ylabel("cl")
    axes[0].set_xlim(left=0.0)
    axes[0].set_title("NACA0012 2.5D M0.5/alpha2 lift vs h")
    axes[0].legend(frameon=False, fontsize=8)
    err = np.abs(f / cl_mid - 1.0) * 100.0
    axes[1].bar([r["level"] for r in recs], err, color=S1_BLUE, width=0.55)
    axes[1].axhline(1.0, color=CRITICAL, ls="--", lw=1.2,
                    label="G9.2 fine bound (1%)")
    axes[1].set_ylabel("|cl error| vs midpoint  [%]")
    axes[1].set_title("sharp-TE lift oracle: monotone error decrease")
    axes[1].legend(frameon=False, fontsize=8)
    finish(fig, OUT, "g92_naca_lift_oracle.png")


# ------------------------------------------------------------------ main

def main():
    apply_style()
    cl = CheckList("P9 grid-convergence & accuracy-gap discrimination")

    # ---- G9.1 -----------------------------------------------------------
    print("G9.1: ONERA M6 three-point grid study")
    m6 = [run_m6_level(lv) for lv in LEVELS]
    h = [r["h"] for r in m6]
    rich = {}
    for key in ("cl_kj", "cl_p"):
        f = [r[key] for r in m6]
        mono = (f[0] < f[1] < f[2]) or (f[0] > f[1] > f[2])
        cl.add("G9.1", f"{key} three-point monotone",
               " -> ".join(f"{v:.4f}" for v in f), "monotone", mono)
        rich[key] = (solve_order_and_extrapolant(h, f) if mono
                     else (None, None))
    for r in m6:
        meta = r["meta"]
        cl.add("G9.1", f"{r['level']} converged Newton solution",
               f"|R|={meta['residual_final']:.1e} lim/flr="
               f"{meta['n_limited']}/{meta['n_floored']}",
               "converged, 0/0", meta["converged"]
               and meta["n_limited"] == 0 and meta["n_floored"] == 0)
    cl.add("G9.1", "fine solve within budget",
           f"{m6[2]['wall_solve']:.0f}s", "<= 7200 s (run-once)",
           m6[2]["wall_solve"] <= 7200.0)

    p_kj, cl_kj_inf = rich["cl_kj"]
    if cl_kj_inf is None:
        verdict = "invalid (non-monotone sequence)"
    elif cl_kj_inf >= BAND_HI:
        verdict = "resolution-dominated"
    elif cl_kj_inf <= BAND_LO:
        verdict = "floor confirmed"
    else:
        verdict = "inconclusive"
    cl.add("G9.1", "pre-registered verdict on cl_KJ_inf",
           f"{cl_kj_inf:.4f}" if cl_kj_inf is not None else "n/a",
           f">= {BAND_HI} resolution / <= {BAND_LO} floor / else "
           "inconclusive", True, note=verdict)
    fig_g91(m6, rich)

    # ---- G9.2 -----------------------------------------------------------
    print("G9.2: NACA0012 2.5D subsonic lift oracle")
    cl_pg, cl_kt = load_cl_bracket(ALPHA_SUB)
    cl_mid = 0.5 * (cl_pg + cl_kt)
    naca = [run_naca_level(lv) for lv in LEVELS]
    errs = [abs(r["cl_p"] / cl_mid - 1.0) for r in naca]
    for r in naca:
        cl.add("G9.2", f"{r['level']} converged",
               f"|R|={r['residual_final']:.1e}", "converged",
               r["converged"])
    cl.add("G9.2", "|error| decreases monotonically",
           " -> ".join(f"{100*e:.2f}%" for e in errs),
           "monotone decrease", errs[0] > errs[1] > errs[2])
    cl.add("G9.2", "fine within +-1% of PG-KT midpoint",
           f"{100*errs[2]:.2f}%", "<= 1%", errs[2] <= 0.01)
    fig_g92(naca, cl_mid, cl_pg, cl_kt)

    # ---- G9.3 -----------------------------------------------------------
    gap_total = CL_EXTERNAL_REF - m6[1]["cl_kj"]     # vs the P8 medium anchor
    if cl_kj_inf is not None:
        resolution_share = cl_kj_inf - m6[1]["cl_kj"]
        floor_share = CL_EXTERNAL_REF - cl_kj_inf
    else:
        resolution_share = float("nan")
        floor_share = float("nan")
    write_csv(OUT, "verdict.csv",
              "quantity,value,note",
              [("cl_kj_coarse", f"{m6[0]['cl_kj']:.4f}", ""),
               ("cl_kj_medium", f"{m6[1]['cl_kj']:.4f}",
                "P8 anchor 0.2692"),
               ("cl_kj_fine", f"{m6[2]['cl_kj']:.4f}", ""),
               ("observed_order_p_clkj",
                "n/a" if p_kj is None else f"{p_kj:.3f}", ""),
               ("cl_kj_inf",
                "n/a" if cl_kj_inf is None else f"{cl_kj_inf:.4f}",
                "Richardson"),
               ("cl_p_inf",
                "n/a" if rich["cl_p"][1] is None
                else f"{rich['cl_p'][1]:.4f}", "Richardson"),
               ("external_ref", f"{CL_EXTERNAL_REF:.3f}", "Tranair/KRATOS"),
               ("gap_vs_medium", f"{gap_total:.4f}", ""),
               ("resolution_share", f"{resolution_share:.4f}",
                "cl_kj_inf - medium"),
               ("floor_share", f"{floor_share:.4f}",
                "external_ref - cl_kj_inf"),
               ("g91_verdict", verdict, "pre-registered bands"),
               ("g92_fine_err_pct", f"{100*errs[2]:.3f}",
                "sharp-TE 2D lift floor test")])

    # ---- committed summary tables ---------------------------------------
    write_csv(OUT, "g91_m6_levels.csv",
              "level,n_tets,n_nodes,h,cl_kj,cl_p,m_max,shock_044,shock_065,"
              "shock_090,wall_mesh_s,wall_solve_s,n_cont_levels,"
              "level_accepts",
              [(r["level"], r["n_tets"], r["n_nodes"], f"{r['h']:.6f}",
                f"{r['cl_kj']:.5f}", f"{r['cl_p']:.5f}", f"{r['m_max']:.3f}",
                f"{r['shocks'][0.44]:.3f}", f"{r['shocks'][0.65]:.3f}",
                f"{r['shocks'][0.90]:.3f}", f"{r['wall_mesh']:.1f}",
                f"{r['wall_solve']:.1f}", r["meta"]["n_levels"],
                ";".join(r["meta"]["level_accepts"])) for r in m6])
    write_csv(OUT, "g92_naca_levels.csv",
              "level,n_tets,n_nodes,cl_p,cl_kj,err_vs_midpoint_pct,"
              "wall_solve_s",
              [(r["level"], r["n_tets"], r["n_nodes"], f"{r['cl_p']:.5f}",
                f"{r['cl_kj']:.5f}", f"{100*e:.3f}",
                f"{r['wall_solve']:.1f}") for r, e in zip(naca, errs)])
    print(f"\nG9.3 verdict: {verdict}  (cl_KJ_inf="
          f"{'n/a' if cl_kj_inf is None else f'{cl_kj_inf:.4f}'}, "
          f"resolution share {resolution_share:+.4f}, floor share "
          f"{floor_share:+.4f} of the 0.019 gap)")
    return cl.report(OUT)


if __name__ == "__main__":
    sys.exit(main())
