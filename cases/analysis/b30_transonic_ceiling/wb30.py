"""B30 shared harness: (b)-class wing-body transonic ceiling attribution
(PRE_REGISTRATION.md, 2026-07-21).

Two warm-started single-level strict solves at the production recipes'
dying levels, zero pyfp3d changes:

  LS    = the B29 production side (B25 inboard clip + B28 flat sheet),
          m = 0.7875 from the committed ls_flat_medium_084.npz seed
          (phi_ext @ 0.775); recipe = LS_RAMP_KW minus the ramp-only keys.
  CONF  = conforming, m = 0.80 from conf_medium_079.npz (phi @ 0.79; Gamma
          re-solved from zero under the warm phi -- pre-reg risk T6);
          recipe = CONF_RAMP_NK + freeze_tol=1e-5 + kutta pressure.

plus the census machinery replicating the solvers' OWN clamp accounting
(limit_q2_field caps q2 at q2(M=3); rho_tilde floor at 0.05) -- the LS
per-side loop of multivalued.py:472-486 and the conforming chain of
newton.py:286-311 -- validated against the recorded n_limited/n_floored
(pre-reg GB30.2 self-check). Per-iteration counts come from the returned
step_records (LS) / clamp_history (CONF); element-level per-iteration
sets are NOT available without a library change (pre-registered
limitation, GB30.2 item 3).
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
B23 = REPO_ROOT / "cases/analysis/b23_junction_discriminator"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(B23))

from wb_common import FUS, LS_MESH_DIR, load_mesh  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.meshgen.fuselage import make_inboard_clip  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI  # noqa: E402
from pyfp3d.meshgen.wingbody import te_polyline  # noqa: E402
from pyfp3d.physics.isentropic import (density_field, limit_q2_field,  # noqa: E402
                                      mach_squared_field)
from pyfp3d.post.surface import planform_area  # noqa: E402
from pyfp3d.post.unified import wall_forces  # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting  # noqa: E402
from pyfp3d.solve.newton_ls import solve_multivalued_newton  # noqa: E402
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet  # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

CONF_DIR = REPO_ROOT / "cases/meshes/onera_m6_wingbody_conforming"
DEMO_OUT = REPO_ROOT / "cases/demo/b18_wingbody_transonic/results"

ALPHA = 3.06
# production recipes, verbatim minus the ramp-only keys (pre-reg section 2)
LS_LEVEL_KW = dict(farfield="freestream", farfield_aux="pin_gamma",
                   freeze_tol=1e-4, freeze_max_clamped=8,
                   direct_refactor_every=1000, n_newton_max=80)
CONF_LEVEL_KW = dict(freeze_refresh_max=8, precond="direct",
                     direct_refactor_every=1000, n_newton_max=80,
                     farfield_spanwise_gamma=True, freeze_tol=1e-5,
                     kutta_estimator="pressure")
UPWIND_DEFAULT = dict(upwind_c=1.5, m_crit=0.95, m_cap=3.0, rho_floor=0.05)

# committed B29 anchors (checks.csv / B29 ledger) for GB30.1
ANCHORS = {
    "conf_medium_079": dict(clp=0.2579, reached=True),
    "conf_coarse_084": dict(clp=0.2617, reached=True),
    "ls_flat_coarse_084": dict(clp=0.2551, m_last=0.84, reached=True),
    "ls_flat_medium_084": dict(m_last=0.775, m_final=0.7875, reached=False,
                               cls="a+dm"),
}


# ------------------------------------------------------------------ builders
def build_ls_flat(mesh, alpha_deg: float = ALPHA):
    """The B29 PRODUCTION side: B25 inboard fragment clip + B28 flat sheet
    (geometry dragged along +x at y=0, physics convecting with the flow)."""
    a = np.radians(alpha_deg)
    wls = WakeLevelSet(te_polyline(FUS), direction=(np.cos(a), np.sin(a), 0.0),
                       sheet_direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]),
                       inboard_clip=make_inboard_clip(FUS))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return wls, cm, mvop


def build_conf(level: str):
    return cut_wake(read_mesh(str(CONF_DIR / f"{level}.msh")))


# ------------------------------------------------------------------ runners
def _cache_load(cache):
    d = np.load(cache, allow_pickle=True)
    rec = dict(json.loads(str(d["summary"])))
    rec["phi"] = d["phi"]
    rec["gamma"] = d["gamma"]
    rec["cached"] = True
    print(f"  [{cache.stem}] CACHED (m={rec['m_inf']}, "
          f"converged={rec['converged']})", flush=True)
    return rec


def _cache_save(cache, rec, phi, gamma):
    summary = {k: v for k, v in rec.items()
               if k not in ("phi", "gamma", "workspace", "step_records",
                            "timings", "residual_history", "F_history",
                            "gamma_history", "clamp_history")}
    np.savez(cache, phi=phi, gamma=np.atleast_1d(gamma),
             summary=json.dumps(summary, default=str))


def run_ls_level(mesh, mvop, m_inf, phi_init, upwind_c=1.5, tag="g2"):
    """One STRICT LS level (production recipe) warm-started from phi_init.
    Cached to results/{tag}_ls_m{...}.npz."""
    cache = OUT / f"{tag}_ls_m{str(round(m_inf, 4)).replace('.', '_')}.npz"
    if cache.exists():
        return _cache_load(cache)
    t0 = time.perf_counter()
    r = solve_multivalued_newton(
        mvop, mesh, m_inf, alpha_deg=ALPHA, phi_init=phi_init,
        gamma_init=0.0, n_seed=0, upwind_c=upwind_c, verbose=True,
        **LS_LEVEL_KW)
    wall_s = time.perf_counter() - t0
    rec = dict(m_inf=float(m_inf), upwind_c=float(upwind_c),
               converged=bool(r["converged"]),
               accept_reason=str(r["accept_reason"]),
               residual_norm=float(r["residual_history"][-1]),
               n_newton=int(r["n_newton"]),
               n_limited=int(r["n_limited"]), n_floored=int(r["n_floored"]),
               mach_max=float(np.sqrt(r["mach2_max"])),
               gamma=float(r["gamma"]), cl_kj=float(r["cl_kj"]),
               froze=bool(r["froze"]),
               n_freeze_refresh=int(r["n_freeze_refresh"]),
               n_freeze_reverts=int(r["n_freeze_reverts"]),
               step_nlim=[s["n_limited"] for s in r["step_records"]],
               step_nflr=[s["n_floored"] for s in r["step_records"]],
               step_res=[float(x) for x in r["residual_history"]],
               wall_s=wall_s)
    _cache_save(cache, rec, r["phi_ext"], r["gamma"])
    rec["phi"], rec["gamma"] = r["phi_ext"], r["gamma"]
    rec["cached"] = False
    print(f"  [ls m={m_inf}] converged={rec['converged']} "
          f"res={rec['residual_norm']:.2e} nlim={rec['n_limited']} "
          f"nflr={rec['n_floored']} Mmax={rec['mach_max']:.2f} "
          f"({wall_s:.0f}s)", flush=True)
    return rec


def run_conf_level(mc, wc, m_inf, phi_init, gamma_init=None, upwind_c=1.5,
                   tag="g2"):
    """One STRICT conforming level (production recipe) warm-started from
    (phi_init, gamma_init). NB (T6, measured 2026-07-21): gamma_init=None
    -> zeros is a MEASURED BAD seed -- with the wake jump and the vortex
    far-field both killed the level diverges (lim 22k at it 0, |R| 1e10 at
    it 1). Reconstruct Gamma from the cached cut phi via
    `kutta_targets(phi, wc)` instead.
    Cached to results/{tag}_conf_m{...}.npz."""
    cache = OUT / f"{tag}_conf_m{str(round(m_inf, 4)).replace('.', '_')}.npz"
    if cache.exists():
        return _cache_load(cache)
    t0 = time.perf_counter()
    r = solve_newton_lifting(
        mc, wc, m_inf=m_inf, alpha_deg=ALPHA, phi_init=phi_init,
        gamma_init=gamma_init, n_picard_seed=0, upwind_c=upwind_c,
        verbose=True, **CONF_LEVEL_KW)
    wall_s = time.perf_counter() - t0
    clamps = r["clamp_history"]
    rec = dict(m_inf=float(m_inf), upwind_c=float(upwind_c),
               converged=bool(r["converged"]),
               accept_reason=str(r["accept_reason"]),
               residual_norm=float(r["residual_history"][-1]),
               n_newton=int(r["n_newton"]),
               n_limited=int(clamps[-1][0]) if clamps else 0,
               n_floored=int(clamps[-1][1]) if clamps else 0,
               mach_max=float(np.sqrt(r.get("mach2_max", 0.0))),
               gamma=float(np.mean(r["gamma"])),
               froze=bool(r.get("froze", False)),
               step_nlim=[int(c[0]) for c in clamps],
               step_nflr=[int(c[1]) for c in clamps],
               step_res=[float(x) for x in r["residual_history"]],
               wall_s=wall_s)
    _cache_save(cache, rec, r["phi"], r["gamma"])
    rec["phi"], rec["gamma"] = r["phi"], r["gamma"]
    rec["cached"] = False
    print(f"  [conf m={m_inf}] converged={rec['converged']} "
          f"res={rec['residual_norm']:.2e} nlim={rec['n_limited']} "
          f"nflr={rec['n_floored']} Mmax={rec['mach_max']:.2f} "
          f"({wall_s:.0f}s)", flush=True)
    return rec


# ------------------------------------------------------- census (GB30.2)
def clamp_masks_ls(mesh, mvop, phi_ext, m_inf, upwind_c=1.5, m_crit=0.95,
                   m_cap=3.0, rho_floor=0.05, gamma_air=1.4, u_inf=1.0):
    """(limited, floored) element masks replicating multivalued.py:472-486
    (the solver's own accounting -- the pre-reg self-check validates the
    counts against the recorded n_limited/n_floored)."""
    in_upper, in_lower = mvop._side_element_sets()
    upw_u, upw_l = mvop._side_upwind()
    phi_up, phi_lo = mvop.side_potentials(phi_ext)
    n_tets = mesh.elements.shape[0]
    limited = np.zeros(n_tets, dtype=bool)
    floored = np.zeros(n_tets, dtype=bool)
    for phi_s, upw, keep in ((phi_up, upw_u, in_upper),
                             (phi_lo, upw_l, in_lower)):
        grad, q2 = mvop.op.velocities(phi_s)
        grad, q2 = mvop._apply_main_density(phi_ext, grad, q2)
        q2n = q2 / u_inf ** 2
        q2l = limit_q2_field(q2n, m_inf, m_cap, gamma_air)
        limited |= (q2l != q2n) & keep
        rho = density_field(q2l, m_inf, gamma_air)
        rt = upw.rho_tilde(grad, q2l, rho, m_inf, upwind_c, m_crit,
                           gamma_air, rho_floor)
        floored |= (rt == rho_floor) & keep
    return limited, floored


def clamp_masks_conf(ws, phi, gamma, m_inf, upwind_c=1.5, m_crit=0.95,
                     m_cap=3.0, rho_floor=0.05, u_inf=1.0):
    """(limited, floored) element masks replicating newton.py:286-311."""
    phi_red = np.asarray(phi, dtype=np.float64)[:ws.n_red]
    phi_cut = ws.con.expand(phi_red, np.atleast_1d(np.asarray(gamma,
                                                             dtype=np.float64)))
    grad, q2 = ws.op.velocities(phi_cut)
    q2n = q2 / u_inf ** 2
    q2l = limit_q2_field(q2n, m_inf, m_cap, ws.gamma_air)
    limited = q2l != q2n
    rho = density_field(q2l, m_inf, ws.gamma_air)
    rho_t = ws.upw.rho_tilde(grad, q2l, rho, m_inf, upwind_c, m_crit,
                             ws.gamma_air, rho_floor)
    return limited, rho_t == rho_floor


def _radius_at_vec(p, x):
    """Vectorized fuselage.radius_at (same piecewise meridian profile:
    ellipsoid nose / cylinder / cone / sphere cap, 0 outside)."""
    x = np.asarray(x, dtype=np.float64)
    R = np.zeros_like(x)
    nose = (x > p.x_nose_tip) & (x < p.x_nose_end)
    t = (p.x_nose_end - x[nose]) / p.l_nose
    R[nose] = p.r_f * np.sqrt(np.maximum(0.0, 1.0 - t * t))
    body = (x >= p.x_nose_end) & (x <= p.x_body_end)
    R[body] = p.r_f
    tail = (x > p.x_body_end) & (x <= p.x_tail_start)
    t = (x[tail] - p.x_body_end) / p.l_tail
    R[tail] = p.r_f * (1.0 - t) + p.r_tail * t
    cap = (x > p.x_tail_start) & (x < p.x_tail_tip)
    dx = x[cap] - p.x_tail_start
    R[cap] = np.sqrt(np.maximum(0.0, p.r_tail ** 2 - dx ** 2))
    return R


def region_codes(cents, b_semi: float = B_SEMI):
    """Mutually exclusive region per element centroid (priority: tip_box >
    corridor > near_fus > field), definitions per the B26/B23 machines."""
    code = np.full(len(cents), "field", dtype=object)
    rr = np.hypot(cents[:, 1], cents[:, 2])
    dist_fus = np.abs(rr - _radius_at_vec(FUS, cents[:, 0]))
    code[dist_fus < 0.10] = "near_fus"
    code[(cents[:, 2] < 0.5) & (cents[:, 0] > 0.8)] = "corridor"
    code[cents[:, 2] > 0.95 * b_semi] = "tip_box"
    return code, dist_fus


def census_rows(leg, cls_name, mask, m2, cents, code):
    """One row per (clamp class x region): count + peak M + argmax coords."""
    rows = []
    for region in ("tip_box", "corridor", "near_fus", "field"):
        sel = mask & (code == region)
        row = dict(leg=leg, clamp=cls_name, region=region,
                   count=int(np.count_nonzero(sel)))
        if row["count"]:
            i = int(np.argmax(np.where(sel, m2, -1.0)))
            row.update(peak_mach=float(np.sqrt(m2[i])),
                       pk_x=float(cents[i, 0]), pk_y=float(cents[i, 1]),
                       pk_z=float(cents[i, 2]))
        rows.append(row)
    return rows


def oscillation_stats(rec, tail=20):
    """Count-level oscillation summary of the last `tail` Newton steps."""
    out = {}
    for key in ("step_nlim", "step_nflr"):
        s = np.asarray(rec.get(key, []), dtype=float)
        s = s[-tail:] if s.size > tail else s
        out[key + "_min"] = float(s.min()) if s.size else 0.0
        out[key + "_max"] = float(s.max()) if s.size else 0.0
        out[key + "_last"] = float(s[-1]) if s.size else 0.0
        out[key + "_range"] = out[key + "_max"] - out[key + "_min"]
    return out


def mach2_ls(mvop, phi_ext, m_inf):
    return np.asarray(mvop.element_mach2(phi_ext, m_inf, 1.4, 1.0))


def mach2_conf(ws, phi, gamma, m_inf, u_inf=1.0):
    phi_red = np.asarray(phi, dtype=np.float64)[:ws.n_red]
    phi_cut = ws.con.expand(phi_red, np.atleast_1d(np.asarray(gamma,
                                                             dtype=np.float64)))
    _grad, q2 = ws.op.velocities(phi_cut)
    return mach_squared_field(q2 / u_inf ** 2, m_inf, ws.gamma_air)


def cl_p_ls(mesh, mvop, phi_ext, m_inf):
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    return float(wall_forces(mesh, mvop=mvop, phi_ext=phi_ext,
                             alpha_deg=ALPHA, s_ref=s_ref, m_inf=m_inf,
                             wall_tag="wall")["cl"])


def cl_p_conf(mc, phi, m_inf):
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    return float(wall_forces(mc, phi=phi, alpha_deg=ALPHA, s_ref=s_ref,
                             m_inf=m_inf, wall_tag="wall")["cl"])


def corridor_mmax(m2, cents):
    """Peak local Mach in the junction corridor (z<0.5, x>0.8; B26 def)."""
    corr = (cents[:, 2] < 0.5) & (cents[:, 0] > 0.8)
    if not np.count_nonzero(corr):
        return 0.0, (0.0, 0.0, 0.0)
    idx = np.where(corr)[0]
    i = idx[int(np.argmax(m2[idx]))]
    return float(np.sqrt(m2[i])), tuple(float(v) for v in cents[i])
