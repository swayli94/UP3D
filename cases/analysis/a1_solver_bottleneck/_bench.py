"""
Track A / A1 -- shared benchmark machinery for the conforming-vs-level-set,
Picard-vs-Newton cost comparison.

One `run_case` per (wake model, method, regime), returning a NORMALISED record
so the four drivers -- whose result dicts differ in almost every key -- can be
rolled into one CSV and one set of figures. The timing schema is the canonical
`pyfp3d.solve.timing` one every driver now reports (A1 instrumentation refactor).

Nothing here re-tunes a solver: the transonic legs use the committed recipes
(`continuation.TRANSONIC_DEFAULTS`, `picard_ls.B_TRANSONIC_DEFAULTS`, the P8
Newton recipe) verbatim.

SCOPE: these are the 2.5-D NACA recipes. The 3-D M6 leg has its OWN committed
recipes -- `B_NEWTON_M6_DEFAULTS` (B15) and the m6_medium_ls_workflow Picard
settings -- which differ materially on that mesh (freeze_tol, lagged LU, dm and
the per-level budget). `run_a1_m6.py` therefore calls the LS solvers directly
and must NOT route through `run_ls`; see the LS_*_M6_KW block there. GA1.5's
anchor comparison is only a reproduction if each leg runs ITS OWN committed
recipe.
"""

import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                          # noqa: E402
from pyfp3d.post.shock import shock_report                        # noqa: E402
from pyfp3d.post.unified import section_cp, wall_forces           # noqa: E402
from pyfp3d.solve.continuation import (                           # noqa: E402
    TRANSONIC_DEFAULTS,
    solve_transonic_lifting,
)
from pyfp3d.solve.newton import (                                 # noqa: E402
    solve_newton_lifting,
    solve_newton_transonic,
)
from pyfp3d.solve.newton_ls import (                              # noqa: E402
    solve_multivalued_newton,
    solve_multivalued_newton_transonic,
)
from pyfp3d.solve.picard import solve_subsonic_lifting            # noqa: E402
from pyfp3d.solve.picard_ls import (                              # noqa: E402
    B_TRANSONIC_DEFAULTS,
    solve_multivalued_lifting,
    solve_multivalued_transonic,
)
from pyfp3d.solve.timing import PHASES                            # noqa: E402
from pyfp3d.wake import (                                         # noqa: E402
    CutElementMap,
    MultivaluedOperator,
    WakeLevelSet,
)

# Committed conforming-Newton transonic recipe (tests/test_p8_newton.py
# NEWTON_TRANSONIC_RECIPE, G8.1): exact direct steps because the shock-position
# soft mode leaves Krylov-eta steps stalled; fine Mach steps near the top.
# Inlined so the harness does not import the test package.
NEWTON_TRANSONIC_RECIPE = dict(
    dm=0.025, dm_min=0.003, freeze_tol=1e-6,
    newton_kw=dict(freeze_refresh_max=8, precond="direct", n_newton_max=60),
)

# The four methods, in a fixed order (bars/curves read the same everywhere).
METHODS = ["conf_picard", "conf_newton", "ls_picard", "ls_newton"]
METHOD_LABEL = {
    "conf_picard": "conforming Picard",
    "conf_newton": "conforming Newton",
    "ls_picard": "level-set Picard",
    "ls_newton": "level-set Newton",
}


# --------------------------------------------------------------------------
# mesh builders
# --------------------------------------------------------------------------

def build_conforming(mesh_path):
    """(mesh_cut, wake_constraint) for the conforming wake-cut path."""
    mesh = read_mesh(str(mesh_path))
    mc, wc = cut_wake(mesh)
    return mc, wc


def build_ls(mesh_path):
    """(mesh, MultivaluedOperator) for the level-set path. The TE polyline
    endpoints ride the mesh's own span extent so they land on wall nodes."""
    mesh = read_mesh(str(mesh_path))
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def span_extent(mesh):
    return float(np.ptp(mesh.nodes[:, 2]))


# --------------------------------------------------------------------------
# case runners -- each returns a normalised record
# --------------------------------------------------------------------------

def _base_record(method, wake_model, mesh_path, regime, m_inf, alpha):
    return dict(method=method, wake_model=wake_model,
                mesh=Path(mesh_path).parent.name + "/" + Path(mesh_path).stem,
                regime=regime, m_inf=float(m_inf), alpha_deg=float(alpha))


def run_conforming(mc, wc, mesh_path, method, regime, m_inf, alpha,
                   m_start=0.60, u_inf=1.0):
    """conforming Picard/Newton, subsonic single level or transonic ramp."""
    t0 = time.perf_counter()
    if regime == "subsonic":
        if method == "conf_picard":
            r = solve_subsonic_lifting(mc, wc, m_inf=m_inf, alpha_deg=alpha,
                                       n_picard_max=60)
        else:
            r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=alpha)
        timings, levels = r["timings"], None
        n_outer = r.get("n_picard", r.get("n_newton"))
    else:
        if method == "conf_picard":
            r = solve_transonic_lifting(
                mc, wc, m_inf=m_inf, alpha_deg=alpha, m_start=m_start,
                dm=0.05, n_picard_seed=TRANSONIC_DEFAULTS["n_picard_seed"],
                n_picard_eval=TRANSONIC_DEFAULTS["n_picard_eval"],
                max_gamma_evals=TRANSONIC_DEFAULTS["max_gamma_evals"])
        else:
            # committed converging recipe (G8.1); dm comes from the recipe
            r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=alpha,
                                       m_start=m_start,
                                       **NEWTON_TRANSONIC_RECIPE)
        timings = r["timings_total"]
        levels = r["level_results"]
        n_outer = sum(_lvl_iters(lv) for lv in levels)
    wall = time.perf_counter() - t0
    rec = _base_record(method, "conforming", mesh_path, regime, m_inf, alpha)
    rec.update(_finish(r, timings, levels, wall, n_outer, path="conforming"))
    _attach_trajectory(rec, r, levels)
    rec["phi"] = r["phi"]
    rec["gamma_scalar"] = float(np.mean(r["gamma"]))
    rec["_wc"] = wc
    rec["_gamma_stations"] = np.atleast_1d(r["gamma"])
    return rec


def run_ls(mesh, mvop, mesh_path, method, regime, m_inf, alpha,
           m_start=0.60, u_inf=1.0, farfield="vortex"):
    """level-set Picard/Newton, subsonic single level or transonic ramp."""
    t0 = time.perf_counter()
    if regime == "subsonic":
        if method == "ls_picard":
            r = solve_multivalued_lifting(
                mvop, mesh, m_inf=m_inf, alpha_deg=alpha, n_outer_max=120,
                tol_residual=1e-8, farfield=farfield)
        else:
            r = solve_multivalued_newton(
                mvop, mesh, m_inf=m_inf, alpha_deg=alpha, n_seed=40,
                farfield=farfield)
        timings, levels = r["timings"], None
        n_outer = r["n_outer"] if method == "ls_picard" else r["n_newton"]
    else:
        if method == "ls_picard":
            # B_TRANSONIC_DEFAULTS carries its own m_start/dm; drop them so the
            # shared ramp start (m_start arg) wins and the four ramps line up.
            ls_kw = {k: v for k, v in B_TRANSONIC_DEFAULTS.items()
                     if k not in ("m_start", "dm")}
            r = solve_multivalued_transonic(
                mvop, mesh, m_target=m_inf, alpha_deg=alpha, m_start=m_start,
                dm=0.05, tol_residual=1e-5, farfield=farfield, **ls_kw)
        else:
            r = solve_multivalued_newton_transonic(
                mvop, mesh, m_target=m_inf, alpha_deg=alpha, m_start=m_start,
                dm=0.05, freeze_tol=1e-6, farfield=farfield)
        timings = r["timings_total"]
        levels = r["levels"]
        n_outer = sum(_lvl_iters(lv) for lv in levels)
    wall = time.perf_counter() - t0
    rec = _base_record(method, "level-set", mesh_path, regime, m_inf, alpha)
    rec.update(_finish(r, timings, levels, wall, n_outer, path="ls"))
    _attach_trajectory(rec, r, levels)
    rec["phi_ext"] = r["phi_ext"]
    rec["gamma_scalar"] = float(r["gamma"])
    rec["_mvop"] = mvop
    rec["_te_jump"] = np.atleast_1d(r["te_jump"])
    return rec


def _lvl_iters(lv):
    return int(lv.get("n_outer", lv.get("n_newton", 0)))


def _attach_trajectory(rec, r, levels):
    """Per-iteration step records (flattened over levels) + per-level rows,
    for the convergence and ramp figures."""
    if levels is None:
        rec["_steps"] = r.get("step_records", [])
        rec["_levels_steps"] = None
        rec["_level_rows"] = None
    else:
        flat, rows, cum = [], [], 0.0
        for lv in levels:
            # each level is a separate driver call whose step wall_cum_s starts
            # from ITS OWN clock; offset by the prior levels' wall so the
            # flattened trace is monotonic across the whole ramp (otherwise the
            # "convergence per wall-second" axis resets every Mach level).
            steps = lv.get("step_records", [])
            for s in steps:
                s = dict(s)
                s["wall_cum_s"] = cum + s.get("wall_cum_s", 0.0)
                flat.append(s)
            rows.append({
                "m": float(lv["m_inf"] if "m_inf" in lv else lv["m"]),
                "n_iter": _lvl_iters(lv),
                "wall_s": float(lv.get("wall_s", 0.0)),
                "residual": float(lv.get("residual_norm",
                                         lv.get("residual_history", [np.nan])[-1]
                                         if lv.get("residual_history")
                                         else np.nan)),
            })
            cum += float(lv.get("wall_s", 0.0))
        rec["_steps"] = flat
        rec["_levels_steps"] = flat
        rec["_level_rows"] = rows


def _finish(r, timings, levels, wall, n_outer, path):
    """Common normalised fields: convergence, iteration + linear-algebra
    counts, cl_p, and the per-phase timing split."""
    out = {
        "converged": bool(r["converged"]),
        "wall_s": float(wall),
        "n_outer_total": int(n_outer),
        "n_levels": 0 if levels is None else len(levels),
    }
    for p in PHASES + ("other",):
        out["t_" + p] = float(timings[p])
    out["t_wall_accounted"] = float(timings["wall"])
    # linear-algebra tallies. For a RAMP the top-level counters (n_gmres_total,
    # n_solves_total, n_refactor) describe only the FINAL level -- the same
    # final-level-only footgun as `timings` -- so aggregate across levels when
    # they are present; only fall back to the single-solve counters otherwise.
    if levels:
        out["n_lin_iters"] = int(sum(lv.get("n_lin_iters", 0) for lv in levels))
        out["n_lin_solves"] = int(sum(
            lv.get("n_lin_solves",
                   sum(s.get("n_lin_solves", 0)
                       for s in lv.get("step_records", [])))
            for lv in levels))
        out["n_refactor"] = int(sum(lv.get("n_refactor", 0) for lv in levels))
    else:
        out["n_lin_iters"] = int(r.get("n_gmres_total", r.get("n_cg_total", 0)))
        out["n_lin_solves"] = int(r.get("n_solves_total", 0)) or _count_solves(
            r, None)
        out["n_refactor"] = int(r.get("n_refactor", 0))
    out["n_stalled"] = int(r.get("n_gmres_stalled", 0))
    out["mach_max"] = float(np.sqrt(r["mach2_max"])) if "mach2_max" in r \
        else float("nan")
    out["n_limited"] = int(r.get("n_limited", 0))
    out["n_floored"] = int(r.get("n_floored", 0))
    return out


def _count_solves(r, levels):
    if "step_records" in r:
        return sum(s.get("n_lin_solves", 0) for s in r["step_records"])
    if levels:
        return sum(s.get("n_lin_solves", 0)
                   for lv in levels for s in lv.get("step_records", []))
    return 0


# --------------------------------------------------------------------------
# derived metrics + post-processing for the figures
# --------------------------------------------------------------------------

def add_forces(rec, mesh, m_inf, alpha, s_ref):
    """cl_p from the wall pressure integral, cl_kj from the circulation, and
    the V6 consistency gap between them -- via the unified post dispatch."""
    if rec["wake_model"] == "conforming":
        f = wall_forces(mesh, phi=rec["phi"], alpha_deg=alpha, s_ref=s_ref,
                        m_inf=m_inf)
    else:
        f = wall_forces(mesh, mvop=rec["_mvop"], phi_ext=rec["phi_ext"],
                        alpha_deg=alpha, s_ref=s_ref, m_inf=m_inf)
    rec["cl_p"] = float(f["cl"])
    # sectional cl_KJ (2.5-D): 2 gamma / (U c), chord = 1
    rec["cl_kj"] = 2.0 * rec["gamma_scalar"]
    denom = abs(rec["cl_kj"]) if abs(rec["cl_kj"]) > 1e-9 else 1.0
    rec["v6_consistency"] = abs(rec["cl_p"] - rec["cl_kj"]) / denom
    return rec


def section_curve(rec, mesh, m_inf, z, b_semi=None):
    """Cp(x/c) on the upper+lower wall at a spanwise station, either path."""
    if rec["wake_model"] == "conforming":
        return section_cp(mesh, phi=rec["phi"], z=z, b_semi=b_semi,
                          m_inf=m_inf)
    return section_cp(mesh, mvop=rec["_mvop"], phi_ext=rec["phi_ext"], z=z,
                      b_semi=b_semi, m_inf=m_inf)


def shock_x(rec, mesh, m_inf, z):
    """Upper-surface shock x/c (nan if subcritical / no shock found)."""
    try:
        curve = section_curve(rec, mesh, m_inf, z)
        rep = shock_report(curve, m_inf)
        return float(rep["upper"].get("x_shock", float("nan")))
    except Exception:
        return float("nan")


CSV_COLUMNS = [
    "regime", "m_inf", "wake_model", "method", "mesh", "converged",
    "wall_s", "n_levels", "n_outer_total",
    "t_seed", "t_assembly", "t_precond", "t_linsolve", "t_residual",
    "t_kutta", "t_other",
    "pct_seed", "pct_assembly", "pct_precond", "pct_linsolve",
    "pct_residual", "pct_kutta", "pct_other",
    "dominant_phase", "dominant_pct",
    "n_lin_solves", "n_lin_iters", "lin_iters_per_solve",
    "lin_solves_per_outer", "n_refactor", "n_stalled",
    "s_per_outer", "mach_max", "n_limited", "n_floored",
    "gamma", "cl_p", "cl_kj", "v6_consistency",
]


def csv_row(rec):
    """Flatten a record into the a1_runs.csv schema (derived %/ratios here)."""
    wall = rec["wall_s"]
    phases = {p: rec["t_" + p] for p in list(PHASES) + ["other"]}
    pct = {p: (100.0 * v / wall if wall > 0 else 0.0)
           for p, v in phases.items()}
    dom = max(PHASES, key=lambda p: phases[p])
    n_outer = max(rec["n_outer_total"], 1)
    n_solves = max(rec["n_lin_solves"], 1)
    row = {
        "regime": rec["regime"], "m_inf": rec["m_inf"],
        "wake_model": rec["wake_model"], "method": rec["method"],
        "mesh": rec["mesh"], "converged": rec["converged"],
        "wall_s": round(wall, 3), "n_levels": rec["n_levels"],
        "n_outer_total": rec["n_outer_total"],
        "dominant_phase": dom, "dominant_pct": round(pct[dom], 1),
        "n_lin_solves": rec["n_lin_solves"], "n_lin_iters": rec["n_lin_iters"],
        "lin_iters_per_solve": round(rec["n_lin_iters"] / n_solves, 1),
        "lin_solves_per_outer": round(rec["n_lin_solves"] / n_outer, 2),
        "n_refactor": rec["n_refactor"], "n_stalled": rec["n_stalled"],
        "s_per_outer": round(wall / n_outer, 4),
        "mach_max": round(rec["mach_max"], 4),
        "n_limited": rec["n_limited"], "n_floored": rec["n_floored"],
        "gamma": round(rec["gamma_scalar"], 6),
        "cl_p": round(rec.get("cl_p", float("nan")), 5),
        "cl_kj": round(rec.get("cl_kj", float("nan")), 5),
        "v6_consistency": round(rec.get("v6_consistency", float("nan")), 5),
    }
    for p in list(PHASES) + ["other"]:
        row["t_" + p] = round(phases[p], 3)
        row["pct_" + p] = round(pct[p], 1)
    return row
