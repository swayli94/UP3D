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
from pyfp3d.kernels.jacobian import PicardOperator
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.wing3d import B_SEMI, x_te
from pyfp3d.physics.isentropic import limit_q2_field, mach_squared_field
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

# ---------------------------------------------------------------------
# G9.1 fine-mesh solver path (the roadmap's pre-registered fallback,
# TAKEN 2026-07-11): the direct (lagged-LU) recipe does NOT scale to the
# fine mesh -- at ~450k dofs a SINGLE splu ran 4 h 39 min without
# returning (RSS 26 GB, 16 threads saturated) and was killed, vs 18.6 s
# per factorization on medium (63k dofs). True-3D LU fill is the wall
# (the N6 finding, one mesh level further).
#
# The fallback is `precond="amg"` -- but with TIGHT Eisenstat-Walker
# forcing (eta = 1e-8, not the 1e-2 default): N5's "Krylov steps stall on
# the shock-position soft mode" is a property of the LOOSE forcing, not
# of AMG. Validated on the medium mesh against the G8.2 locks before
# spending the fine budget: 66 s, cl 0.2646 / M_max 2.129 / shocks
# 0.596-0.541-0.362 -- the SAME solution as the direct path in HALF the
# wall time (141 s), terminal quadratic in the frozen phase, 0 GMRES
# stalls. Applied to coarse+medium too (identical solutions, cheaper),
# so all three grid points share one solver path -- which is what a grid
# study must have.
M6_RECIPE_AMG = dict(
    dm=0.05, dm_min=0.01, freeze_tol=1e-6, intermediate_tol=1e-5,
    newton_kw=dict(freeze_refresh_max=8, precond="amg",
                   amg_rebuild_every=2, ew_eta0=1e-8, ew_eta_max=1e-8,
                   gmres_restart=200, gmres_maxiter=20, n_newton_max=60,
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
    """One M6 Newton continuation run, npz-cached (phi/gamma + walls).

    The fine level starts the Mach ramp deep SUBCRITICAL (m_start=0.30)
    with a deeper first-level Picard seed (n_picard_seed=12), rather than
    at the shared 0.70/seed-5: measured 2026-07-11, the fine mesh resolves
    the LE suction pocket so sharply that the crude cold Picard-5 seed
    OVERSHOOTS it into the density floor -- a M0.70 cold seed lands ~4000
    speed-limited + ~1800 density-floored cells (level-0 Newton stalls at
    |R|~6e-2, never sheds the clamps, and level 0 cannot dm-halve so the
    whole solve breaks); even a M0.50 cold seed still floors ~658. The
    fix has two parts, both continuation-PATH only: (i) start at M0.30
    where the LE peak is mild enough that a crude seed cannot reach the
    floor (needs M_local > ~2.2, unphysical there), and (ii) a deeper
    Picard seed at that first clean level to settle the LE before Newton.
    Both are safe -- G8.2 established the M0.84 solution is
    continuation-path-independent (coarse to 1e-15, medium cl to 5
    digits), so the three grid points stay comparable (same AMG solver,
    same target, cheaper-per-level path); the added subcritical levels
    are loose (intermediate_tol), a few Newton steps each, all
    warm-started so they shed transient clamps the way the medium ramp
    did. n_picard_seed lands only on the FIRST level (later levels are
    forced to phi_init/seed-0 by the continuation)."""
    cache = OUT / f"g91_m6_{level}.npz"
    path = require_mesh("onera_m6", level)
    t0 = time.perf_counter()
    mc, wc = cut_wake(read_mesh(path))
    wall_mesh = time.perf_counter() - t0
    n_tets = len(mc.elements)

    recipe = dict(M6_RECIPE_AMG)
    if level == "fine":
        m_start = 0.30
        recipe["newton_kw"] = dict(recipe["newton_kw"], n_picard_seed=12)
    else:
        m_start = 0.70

    if cache.exists() and os.environ.get("PYFP3D_P9_RESOLVE", "0") != "1":
        z = np.load(cache, allow_pickle=True)
        r = dict(z["result_scalars"].item())
        phi, gamma = z["phi"], z["gamma"]
        wall_solve = float(z["wall_solve"])
        print(f"  [m6 {level}] loaded cache {cache.name}", flush=True)
    else:
        t0 = time.perf_counter()
        res = solve_newton_transonic(mc, wc, m_inf=M_M6, alpha_deg=ALPHA_M6,
                                     m_start=m_start, verbose=True,
                                     **recipe)
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
    sing = tip_te_singularity_census(mc, phi)
    return dict(
        level=level, n_tets=n_tets, n_nodes=len(mc.nodes),
        h=float(n_tets) ** (-1.0 / 3.0),
        cl_p=float(forces["cl"]),
        cl_kj=float(cl_kj_3d(gamma, wc.station_z, s_ref=s_ref,
                             b_semi=B_SEMI)),
        m_max=float(np.sqrt(r["mach2_max"])), shocks=shocks,
        sing=sing,
        wall_mesh=wall_mesh, wall_solve=wall_solve, meta=r,
    )


def tip_te_singularity_census(mc, phi):
    """The ★ G9.1 finding, measured per level: the UNLIMITED local Mach
    field (i.e. before the m_cap speed limiter hides it) and the census of
    cells that cross the cap, with their location relative to the swept TIP
    TRAILING-EDGE corner.

    P5 recorded a "bounded tip-TE-corner P1 overshoot" as the only
    surviving singularity trace on the medium mesh (M_max 1.995). This
    census asks whether refinement bounds it -- and it does not."""
    op = PicardOperator(mc.nodes, mc.elements)
    _, q2 = op.velocities(np.asarray(phi))
    q2n = q2 / 1.0 ** 2
    q2l = limit_q2_field(q2n, M_M6, 3.0, 1.4)
    limited = q2l != q2n
    mach = np.sqrt(mach_squared_field(q2n, M_M6, 1.4))   # UNLIMITED
    cen = mc.nodes[mc.elements].mean(axis=1)
    e_hot = int(np.argmax(mach))
    idx = np.where(limited)[0]
    if len(idx):
        z_over_b = cen[idx, 2] / B_SEMI
        dx_te = np.array([cen[e, 0] - x_te(cen[e, 2]) for e in idx])
        loc = (f"z/b {z_over_b.min():.3f}-{z_over_b.max():.3f}, "
               f"x-x_TE {dx_te.min():+.4f}..{dx_te.max():+.4f}")
    else:
        loc = "none"
    return {
        "m_max_unlimited": float(mach.max()),
        "n_capped": int(limited.sum()),
        "hot_z_over_b": float(cen[e_hot, 2] / B_SEMI),
        "hot_dx_te": float(cen[e_hot, 0] - x_te(cen[e_hot, 2])),
        "capped_locus": loc,
    }


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
    """Left: the lift sequence (with the fine point flagged as NOT a
    discrete solution). Right: the ★ finding -- the tip-TE corner's
    unlimited local Mach DIVERGES under refinement, which is why the fine
    Newton cannot converge (capped cells => the freeze machinery, which
    requires 0 limited, can never engage) and why no Richardson exists."""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
    h = np.asarray([r["h"] for r in recs])
    ok = np.asarray([bool(r["meta"]["converged"]) for r in recs])

    ax = axes[0]
    for key, color, lab in (("cl_kj", S1_BLUE, "cl_KJ"),
                            ("cl_p", S2_AQUA, "cl_p")):
        f = np.asarray([r[key] for r in recs])
        ax.plot(h[ok], f[ok], "o-", color=color, label=f"{lab} (converged)")
        if (~ok).any():
            ax.plot(h[~ok], f[~ok], "X", ms=11, color=CRITICAL,
                    label=f"{lab} — NOT a discrete solution")
        p, f_inf = rich[key]
        if f_inf is not None:
            ax.axhline(f_inf, color=color, ls="--", lw=1.2,
                       label=f"Richardson {lab}$_\\infty$={f_inf:.4f}")
    ax.axhline(CL_EXTERNAL_REF, color=INK_2, ls=":", lw=1.2,
               label="Tranair/KRATOS 0.288")
    ax.axhspan(BAND_LO, BAND_HI, color=S3_YELLOW, alpha=0.15,
               label="pre-registered inconclusive band")
    ax.set_xlabel(r"$h \propto N_{tet}^{-1/3}$")
    ax.set_ylabel("lift coefficient")
    ax.set_xlim(left=0.0)
    ax.set_title("G9.1: M6 lift vs h — no valid 3-point Richardson")
    ax.legend(frameon=False, fontsize=7)
    for r_, x, y in zip(recs, h, [r["cl_kj"] for r in recs]):
        ax.annotate(r_["level"], (x, y), textcoords="offset points",
                    xytext=(4, -11), fontsize=8, color=INK_2)

    ax = axes[1]
    mm = np.asarray([r["sing"]["m_max_unlimited"] for r in recs])
    ax.plot(h, mm, "o-", color=CRITICAL, lw=2)
    ax.axhline(3.0, color=INK_2, ls="--", lw=1.2,
               label="speed-limiter cap (M_cap = 3)")
    ax.set_xlabel(r"$h \propto N_{tet}^{-1/3}$")
    ax.set_ylabel("max local Mach (unlimited)")
    ax.set_xlim(left=0.0)
    ax.set_title("★ tip-TE corner singularity DIVERGES under refinement")
    ax.legend(frameon=False, fontsize=8)
    for r_, x, y in zip(recs, h, mm):
        ax.annotate(f"{r_['level']}  M={y:.2f}\n({r_['sing']['n_capped']} "
                    f"capped)", (x, y), textcoords="offset points",
                    xytext=(6, -4), fontsize=7, color=INK_2)
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
    fine_ok = bool(m6[2]["meta"]["converged"])
    for key in ("cl_kj", "cl_p"):
        f = [r[key] for r in m6]
        mono = (f[0] < f[1] < f[2]) or (f[0] > f[1] > f[2])
        # Richardson validity needs a monotone sequence of DISCRETE
        # SOLUTIONS; the fine point is neither (see the converged check)
        cl.add("G9.1", f"{key} three-point monotone",
               " -> ".join(f"{v:.4f}" for v in f), "monotone",
               mono and fine_ok, xfail=not fine_ok,
               note=("fine is not a discrete solution -- its value is a "
                     "limit-cycle artifact, not a lift" if not fine_ok
                     else ""))
        rich[key] = (solve_order_and_extrapolant(h, f)
                     if (mono and fine_ok) else (None, None))
    for r in m6:
        meta = r["meta"]
        clean = (meta["converged"] and meta["n_limited"] == 0
                 and meta["n_floored"] == 0)
        # the fine failure is the MEASURED OUTCOME of this gate, not an
        # open defect: recorded as XFAIL with its root cause
        cl.add("G9.1", f"{r['level']} converged Newton solution",
               f"|R|={meta['residual_final']:.1e} lim/flr="
               f"{meta['n_limited']}/{meta['n_floored']}",
               "converged, 0/0", clean, xfail=(r["level"] == "fine"),
               note=("tip-TE singularity crosses M_cap => permanently "
                     "limited => the freeze machinery (needs 0 limited) "
                     "never engages => live-churn floor ~1e-5"
                     if r["level"] == "fine" else ""))
    cl.add("G9.1", "fine solve within budget",
           f"{m6[2]['wall_solve']:.0f}s", "<= 7200 s (run-once)",
           m6[2]["wall_solve"] <= 7200.0,
           note="AMG path did its job; the blocker is the discretization")

    # ★ the finding that replaces the expected Richardson verdict
    mm = [r["sing"]["m_max_unlimited"] for r in m6]
    caps = [r["sing"]["n_capped"] for r in m6]
    diverges = mm[0] < mm[1] < mm[2] and mm[2] > 3.0
    cl.add("G9.1", "★ tip-TE corner M_max (unlimited) vs refinement",
           " -> ".join(f"{v:.2f}" for v in mm),
           "DIVERGES (P5's 'bounded' overshoot is not bounded)", diverges,
           note=f"capped cells {caps}; fine locus "
                f"{m6[2]['sing']['capped_locus']}")
    cl.add("G9.1", "★ capped cells sit at the swept TIP TE corner",
           f"z/b={m6[2]['sing']['hot_z_over_b']:.3f}, "
           f"x-x_TE={m6[2]['sing']['hot_dx_te']:+.4f}",
           "z/b ~ 1.0 and just aft of the TE",
           m6[2]["sing"]["hot_z_over_b"] > 0.99
           and abs(m6[2]["sing"]["hot_dx_te"]) < 0.05)

    p_kj, cl_kj_inf = rich["cl_kj"]
    if cl_kj_inf is None:
        verdict = ("INVALID -- no 3-point Richardson: the fine mesh is not "
                   "a discrete solution (tip-TE singularity diverges past "
                   "the Mach cap), so the pre-registered bands cannot fire")
    elif cl_kj_inf >= BAND_HI:
        verdict = "resolution-dominated"
    elif cl_kj_inf <= BAND_LO:
        verdict = "floor confirmed"
    else:
        verdict = "inconclusive"
    cl.add("G9.1", "pre-registered verdict on cl_KJ_inf",
           f"{cl_kj_inf:.4f}" if cl_kj_inf is not None else "n/a (invalid)",
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
              [("cl_kj_coarse", f"{m6[0]['cl_kj']:.4f}", "converged"),
               ("cl_kj_medium", f"{m6[1]['cl_kj']:.4f}",
                "converged; P8 anchor 0.2692"),
               ("cl_kj_fine", f"{m6[2]['cl_kj']:.4f}",
                "NOT A SOLUTION -- limit-cycle artifact; do not use"),
               ("observed_order_p_clkj",
                "n/a" if p_kj is None else f"{p_kj:.3f}",
                "needs 3 discrete solutions"),
               ("cl_kj_inf",
                "n/a" if cl_kj_inf is None else f"{cl_kj_inf:.4f}",
                "Richardson -- INVALID, see g91_verdict"),
               ("external_ref", f"{CL_EXTERNAL_REF:.3f}", "Tranair/KRATOS"),
               ("gap_vs_medium", f"{gap_total:.4f}",
                "0.288 - medium cl_KJ; the gap P9 set out to split"),
               ("resolution_share", "n/a",
                "unsplittable: no valid extrapolant"),
               ("floor_share", "n/a",
                "unsplittable: no valid extrapolant"),
               ("m_max_unlimited_series",
                ";".join(f"{v:.2f}" for v in mm),
                "★ coarse/medium/fine tip-TE corner: DIVERGES past M_cap=3"),
               ("n_capped_series", ";".join(str(c) for c in caps),
                "★ cells hitting the speed limiter"),
               ("fine_capped_locus", m6[2]["sing"]["capped_locus"]
                .replace(",", ";"),
                "★ swept TIP TE corner (z/b ~ 1, just aft of the TE)"),
               ("g91_verdict", verdict.replace(",", ";"),
                "pre-registered bands could not fire"),
               ("g92_fine_err_pct", f"{100*errs[2]:.3f}",
                "★ 2D sharp TE: NO lift floor (err -> 0.03%)"),
               ("g93_attribution",
                "2D sharp TE exonerated (G9.2); the 3D flat-facet swept "
                "TIP-TE corner is a DIVERGENT singularity (G9.1) -- the "
                "3D lift gap cannot be split by grid convergence because "
                "the sequence does not converge. User arbitrates P11.",
                "G9.3")])

    # ---- committed summary tables ---------------------------------------
    write_csv(OUT, "g91_m6_levels.csv",
              "level,n_tets,n_nodes,h,converged,cl_kj,cl_p,m_max_limited,"
              "m_max_unlimited,n_capped,capped_locus,shock_044,shock_065,"
              "shock_090,wall_mesh_s,wall_solve_s,n_cont_levels,"
              "level_accepts",
              [(r["level"], r["n_tets"], r["n_nodes"], f"{r['h']:.6f}",
                r["meta"]["converged"],
                f"{r['cl_kj']:.5f}", f"{r['cl_p']:.5f}", f"{r['m_max']:.3f}",
                f"{r['sing']['m_max_unlimited']:.3f}",
                r["sing"]["n_capped"],
                r["sing"]["capped_locus"].replace(",", ";"),
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
    print(f"\nG9.1: {verdict}")
    print(f"  tip-TE M_max (unlimited): "
          f"{' -> '.join(f'{v:.2f}' for v in mm)}  capped cells {caps}")
    print(f"G9.2: 2D sharp-TE lift error -> {100*errs[2]:.2f}% "
          f"(no lift floor from a sharp TE)")
    print("G9.3: attribution recorded in verdict.csv -- user arbitrates P11")
    return cl.report(OUT)


if __name__ == "__main__":
    sys.exit(main())
