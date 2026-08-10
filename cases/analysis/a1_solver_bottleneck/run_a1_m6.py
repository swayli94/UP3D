"""
Track A / A1 -- solver bottleneck study, 3-D ONERA M6 leg (GATED, ~1 h once).

The interesting bottleneck story lives in 3-D transonic: the B15 measurement was
that the level-set Picard ramp parks on the shock-position residual plateau and
burns its budget there (M6 medium M0.84 ≈ 38 min), while LS Newton reaches a
strict solution in ≈ 11 min. A1 now MEASURES where those minutes go, per phase.

Runs (M6 medium, M0.84 / alpha 3.06), each on ITS OWN committed recipe -- see
the ANCHORS and LS_*_M6_KW blocks below; that pairing is what makes GA1.5 a
reproduction rather than a coincidence:
  * conforming Newton   on onera_m6/medium.msh        (NEWTON_M6_RECIPE)
  * level-set Picard     on onera_m6_wakefree/medium.msh (m6_medium_ls_workflow)
  * level-set Newton     on onera_m6_wakefree/medium.msh (B15 B_NEWTON_M6_DEFAULTS)
All three have committed wall-clock anchors (~145 s conforming Newton -- the
POST-P10-promotion G8.2, not the superseded 250 s; 657.4 s LS Newton and
2304.7 s LS Picard per B15), so GA1.5 checks the harness reproduces them to
±25%. Conforming Picard in 3-D transonic is the P5 45-75 min beast with NO
committed M0.84 anchor and a fold risk -- it runs ONLY under the extra opt-in
`PYFP3D_A1_CONF_PICARD_3D=1`.

MUST run at the 16-thread cap: the anchors were measured there, so a different
thread count measures SMT rather than harness fidelity.

*** COST CAUTION ***  ~1 h of solves on first run; each method is cached to a
gitignored results/*.npz (PNG/CSV are the committed evidence). Re-solve with
PYFP3D_A1_RESOLVE=1. The .msh are gitignored -- regenerate (~30 s) with
  python cases/meshes/onera_m6/generate_onera_m6.py
  python cases/meshes/onera_m6/generate_onera_m6_wakefree.py   (if present)

Run (gated, 16 threads):
    NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
        PYFP3D_TRANSONIC_GATES=1 \
        python cases/analysis/a1_solver_bottleneck/run_a1_m6.py
"""

import csv
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(HERE))

from cases.demo._common import CheckList, apply_style               # noqa: E402
import _bench as B                                                  # noqa: E402
import _figs as F                                                   # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te                      # noqa: E402
from pyfp3d.post.surface import cl_kj_3d, planform_area             # noqa: E402
from pyfp3d.solve.newton_ls import (                                # noqa: E402
    B_NEWTON_M6_DEFAULTS,
    solve_multivalued_newton_transonic,
)
from pyfp3d.solve.picard_ls import solve_multivalued_transonic      # noqa: E402
from pyfp3d.solve.timing import PHASES                              # noqa: E402
from pyfp3d.wake import (                                           # noqa: E402
    CutElementMap,
    MultivaluedOperator,
    WakeLevelSet,
)

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

GATED = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"
RESOLVE = os.environ.get("PYFP3D_A1_RESOLVE", "0") == "1"
CONF_PICARD_3D = os.environ.get("PYFP3D_A1_CONF_PICARD_3D", "0") == "1"

M_INF = 0.84
ALPHA = 3.06
PREFIX = "a1_m6"        # figure/CSV namespace; the 2.5-D leg owns plain "a1"

# Committed wall-clock anchors (@16 threads) -- the GA1.5 reproduction band.
# Each MUST be the wall clock of the run that the recipe below reproduces:
#   conf_newton  ~145 s  demo_report/track_p.md: G8.2 closed at 249.2 s
#                        end-to-end (2026-07-11), but P10 then PROMOTED
#                        `intermediate_tol=1e-5` into NEWTON_M6_RECIPE --
#                        239.5 s -> 140.3 s solve, "gated G8.2 now ~145 s",
#                        final level identical to 4 digits. The recipe below
#                        IS the post-promotion one, so 250 (the pre-promotion
#                        number) was the wrong partner for it and made GA1.5
#                        fail a faithful run at -40%.
#   ls_newton    657.4 s  b15_ls_newton_ramp/results/summary.csv (GB15.4)
#   ls_picard   2304.7 s  m6_medium_ls_workflow/results/summary.csv, quoted
#                         as the committed Picard reference by B15
ANCHORS = {"conf_newton": 145.0, "ls_newton": 657.4, "ls_picard": 2304.7}

# The committed conforming M6 Newton recipe (tests/test_p8_newton.py
# NEWTON_M6_RECIPE, G8.2) -- inlined so this script does not import the test
# package. `precond="direct"` + direct_refactor_every=1000 is the medium
# lagged-LU path (NOT the fine amg trap -- that is a different mesh scale).
NEWTON_M6_RECIPE = dict(
    dm=0.05, dm_min=0.01, freeze_tol=1e-6, intermediate_tol=1e-5,
    newton_kw=dict(freeze_refresh_max=8, precond="direct",
                   direct_refactor_every=1000, n_newton_max=60,
                   farfield_spanwise_gamma=True),
)

# The committed M6-medium LEVEL-SET recipes, verbatim from the runs the anchors
# above were measured on. `_bench.run_ls` must NOT be used here: it carries the
# 2.5-D NACA transonic defaults, which are a materially different recipe on this
# mesh -- freeze_tol 1e-6 sits BELOW the Mach-rising churn floor (2.7e-4 by
# M0.70) so the freeze never arms (the P9/G9.1 wall); no direct_refactor_every
# means no lagged LU (the true-3D splu 18.6 s/step trap, B12/B13); and for a
# budget-bound Picard that parks on the plateau, dm + n_outer_level ARE the wall
# clock. Measuring those against these anchors would not be a reproduction.
LS_NEWTON_M6_KW = dict(n_seed=40, n_newton_max=80, tol_residual=1e-10,
                       **B_NEWTON_M6_DEFAULTS)          # B15, GB15.4
LS_PICARD_M6_KW = dict(m_start=0.60, dm=0.04, n_outer_seed=120,
                       n_outer_level=200, tol_residual=1e-5,
                       direct_refactor_every=1000)      # m6_medium_ls_workflow


def _summarise(rec):
    """Row for a1_m6_runs.csv: phase split + %s + the cost ratios."""
    row = B.csv_row(rec)
    row["wake_model_headline"] = True
    return row


def _build_ls_m6(mesh_path):
    """(mesh, MultivaluedOperator) for the M6 half wing.

    NOT `_bench.build_ls`: that one puts a straight TE polyline at x=1 spanning
    the MESH z-extent, which is right for the 2.5-D NACA slab (chord-1 TE, wing
    spans the whole domain) and wrong here on both counts -- the M6 TE is swept
    (x_te(z)) and the mesh z-extent is the far-field box, not the wing. It
    matched 0 TE nodes and raised. Polyline + direction are the committed M6
    convention (b7_onera_m6, m6_medium_ls_workflow).
    """
    mesh = read_mesh(str(mesh_path))
    a = np.radians(ALPHA)
    wls = WakeLevelSet(
        np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
        direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                     levelset=wls)


def _wrap_ls(r, mvop, mesh, mesh_path, method):
    """Normalise an LS ramp result into the `_bench` record schema.

    The M6 leg calls the LS solvers directly (committed M6 recipes, see
    LS_*_M6_KW) rather than through `_bench.run_ls`, so it reuses the
    normalisation the same way `_wrap_conforming` does.
    """
    timings, levels = r["timings_total"], r["levels"]
    n_outer = sum(B._lvl_iters(lv) for lv in levels)
    rec = B._base_record(method, "level-set", mesh_path, "transonic",
                         M_INF, ALPHA)
    rec.update(B._finish(r, timings, levels, timings["wall"], n_outer,
                         path="ls"))
    B._attach_trajectory(rec, r, levels)
    rec["phi_ext"] = r["phi_ext"]
    rec["gamma_scalar"] = float(r["gamma"])
    rec["_mvop"] = mvop
    rec["_te_jump"] = np.atleast_1d(r["te_jump"])
    return rec


def _cl_kj_3d(rec, s_ref):
    """Replace the 2.5-D sectional cl_KJ that `_bench.add_forces` writes.

    `add_forces` is shared with the 2.5-D leg, where cl_KJ = 2*Gamma is the
    right sectional formula. On a 3-D half wing it is meaningless: the lift
    is the spanwise integral of Gamma(z) over S_ref (post/surface.cl_kj_3d,
    the convention every committed M6 demo uses). The V6 gap is derived from
    cl_KJ, so it is recomputed with it.
    """
    z, g = rec.get("_span_z"), rec.get("_span_gamma")
    if z is None or len(np.atleast_1d(z)) == 0:
        return
    rec["cl_kj"] = float(cl_kj_3d(np.atleast_1d(g), np.atleast_1d(z),
                                  s_ref=s_ref, b_semi=B_SEMI))
    denom = abs(rec["cl_kj"]) if abs(rec["cl_kj"]) > 1e-9 else 1.0
    rec["v6_consistency"] = abs(rec["cl_p"] - rec["cl_kj"]) / denom


def main():
    if not GATED:
        print("A1 3-D leg is gated. Set PYFP3D_TRANSONIC_GATES=1 to run "
              "(~1 h of solves). Nothing done.")
        return 0
    apply_style()
    cl = CheckList("Track A / A1 solver bottleneck -- 3-D ONERA M6 (GA1.5)")
    recs, run_rows, level_rows = {}, [], []

    conf_mesh = REPO_ROOT / "cases/meshes/onera_m6/medium.msh"
    free_mesh = REPO_ROOT / "cases/meshes/onera_m6_wakefree/medium.msh"

    plan = [("conf_newton", conf_mesh, "conforming"),
            ("ls_picard", free_mesh, "level-set"),
            ("ls_newton", free_mesh, "level-set")]
    if CONF_PICARD_3D:
        plan.insert(0, ("conf_picard", conf_mesh, "conforming"))

    for method, mesh_path, _ in plan:
        if not mesh_path.exists():
            print(f"[skip] {mesh_path} not generated (gitignored) -- {method}")
            continue
        rec = _solve_cached(method, mesh_path)
        if rec is None:
            continue
        recs[method] = rec
        run_rows.append(_summarise(rec))
        for lr in rec.get("_level_rows") or []:
            level_rows.append({"method": method, "m": lr["m"],
                               "n_iter": lr["n_iter"],
                               "wall_s": round(lr["wall_s"], 3),
                               "residual": lr["residual"]})

    _write_csv("a1_m6_runs.csv", run_rows,
               B.CSV_COLUMNS + ["wake_model_headline"])
    if level_rows:
        _write_csv("a1_m6_levels.csv", level_rows, list(level_rows[0].keys()))

    # Figures: shared builders, but written under the "a1_m6" prefix -- both
    # legs share results/, so the default "a1" names would overwrite the 2.5-D
    # leg's committed figures with 3-D content (they did, once).
    F.fig_time_breakdown(recs, OUT, "Where the wall clock goes (3-D)",
                         f"M6 medium M{M_INF}", prefix=PREFIX)
    F.fig_convergence(recs, OUT, f"ONERA M6 medium M{M_INF} transonic",
                      prefix=PREFIX)
    F.fig_ramp(recs, OUT, f"Mach-ramp anatomy (3-D, M{M_INF})", prefix=PREFIX)
    F.fig_linear_solver(recs, OUT, f"Linear-algebra anatomy (3-D, M{M_INF})",
                        prefix=PREFIX)
    if any("_span_z" in r for r in recs.values()):
        F.fig_spanwise(recs, OUT,
                       f"Spanwise loading Γ(z) — four methods (M{M_INF})",
                       prefix=PREFIX)

    # GA1.5 -- reproduce the committed anchors within +-25%
    for method, anchor in ANCHORS.items():
        if method not in recs:
            continue
        wall = recs[method]["wall_s"]
        cl.add("GA1.5", f"{method} wall within 25% of committed anchor",
               f"{wall:.0f}s vs {anchor:.0f}s",
               "|wall-anchor|/anchor < 0.25",
               abs(wall - anchor) / anchor < 0.25)
    # GA1.1 in 3-D: unaccounted time still < 5%
    if run_rows:
        worst = max(run_rows, key=lambda r: r["pct_other"])
        cl.add("GA1.1", "unaccounted time < 5% of wall (3-D)",
               f"{worst['pct_other']:.1f}%", "max pct_other < 5%",
               worst["pct_other"] < 5.0, note=worst["method"])
    return cl.report(OUT, fname="checks_m6.csv")


def _solve_cached(method, mesh_path):
    cache = OUT / f"a1_m6_{method}.npz"
    if cache.exists() and not RESOLVE:
        print(f"  [{method}] reuse cache {cache.name} "
              "(PYFP3D_A1_RESOLVE=1 to re-solve)")
        d = np.load(cache, allow_pickle=True)
        rec = {k: (d[k].item() if d[k].shape == () else d[k])
               for k in d.files}
        # The trajectory fields are lists of dicts on the fresh path but come
        # back from npz as object ARRAYS, and both this script and _figs test
        # them with `x or []` / `if x:` -- ambiguous on an ndarray. Restore
        # list-ness at the cache boundary so the cached and fresh paths are
        # the same object (_figs is shared with the 2.5-D leg; don't make it
        # carry the cache's type quirk).
        for k in ("_level_rows", "_steps"):
            if k in rec:
                rec[k] = list(rec[k])
        return rec

    print(f"  [{method}] solving M6 medium M{M_INF} (minutes)...", flush=True)
    t0 = time.perf_counter()
    if method == "conf_picard":
        mc, wc = B.build_conforming(mesh_path)
        s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        rec = B.run_conforming(mc, wc, mesh_path, method, "transonic",
                               M_INF, ALPHA, m_start=0.60)
        rec["_mesh_for_post"] = mc
        _span_conforming(rec, wc)
    elif method == "conf_newton":
        mc, wc = B.build_conforming(mesh_path)
        s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        from pyfp3d.solve.newton import solve_newton_transonic
        r = solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                   **NEWTON_M6_RECIPE)
        rec = _wrap_conforming(r, mc, wc, mesh_path, method, s_ref)
    else:
        mesh, mvop = _build_ls_m6(mesh_path)
        s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
        kw = LS_PICARD_M6_KW if method == "ls_picard" else LS_NEWTON_M6_KW
        solve = (solve_multivalued_transonic if method == "ls_picard"
                 else solve_multivalued_newton_transonic)
        r = solve(mvop, mesh, m_target=M_INF, alpha_deg=ALPHA,
                  farfield="neumann", **kw)
        rec = _wrap_ls(r, mvop, mesh, mesh_path, method)
        _span_ls(rec, mvop, mesh)
    rec["wall_s"] = time.perf_counter() - t0
    B.add_forces(rec, rec["_mesh_for_post"], M_INF, ALPHA, s_ref)
    _cl_kj_3d(rec, s_ref)
    _save_cache(cache, rec)
    print(f"  [{method}] {rec['wall_s'] / 60:.1f} min "
          f"gamma={rec['gamma_scalar']:.4f} conv={rec['converged']}", flush=True)
    return rec


def _wrap_conforming(r, mc, wc, mesh_path, method, s_ref):
    timings = r.get("timings_total", r["timings"])
    levels = r.get("level_results")
    n_outer = sum(lv["n_newton"] for lv in levels) if levels else r["n_newton"]
    wall = timings["wall"]
    rec = B._base_record(method, "conforming", mesh_path, "transonic",
                         M_INF, ALPHA)
    rec.update(B._finish(r, timings, levels, wall, n_outer, path="conforming"))
    B._attach_trajectory(rec, r, levels)
    rec["phi"] = r["phi"]
    rec["gamma_scalar"] = float(np.mean(r["gamma"]))
    rec["_mesh_for_post"] = mc
    _span_conforming(rec, wc, r["gamma"])
    return rec


def _span_conforming(rec, wc, gamma=None):
    g = np.atleast_1d(rec.get("_gamma_stations", gamma))
    rec["_span_z"] = np.asarray(wc.station_z, dtype=float)
    rec["_span_gamma"] = g


def _span_ls(rec, mvop, mesh):
    z = mesh.nodes[mvop.cm.te_nodes, 2]
    o = np.argsort(z)
    rec["_span_z"] = z[o]
    rec["_span_gamma"] = np.atleast_1d(rec["_te_jump"])[o]
    rec["_mesh_for_post"] = mesh


def _save_cache(cache, rec):
    """Persist the harvestable fields (drop live solver objects/arrays)."""
    keep = {k: v for k, v in rec.items()
            if not k.startswith("_") and k not in ("phi", "phi_ext")}
    keep["_span_z"] = rec.get("_span_z", np.array([]))
    keep["_span_gamma"] = rec.get("_span_gamma", np.array([]))
    keep["_level_rows"] = np.array(rec.get("_level_rows") or [], dtype=object)
    keep["_steps"] = np.array(rec.get("_steps") or [], dtype=object)
    # store phi/phi_ext + the mvop-free bits so add_forces still works on reload
    if "phi" in rec:
        keep["phi"] = rec["phi"]
    if "phi_ext" in rec:
        keep["phi_ext"] = rec["phi_ext"]
    np.savez(cache, **{k: np.array(v, dtype=object) if isinstance(v, (list, dict))
                       else v for k, v in keep.items()})


def _write_csv(name, rows, cols):
    if not rows:
        return
    with open(OUT / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {name} ({len(rows)} rows)")


if __name__ == "__main__":
    sys.exit(main())
