#!/usr/bin/env python
"""GV5.2 — RAE2822 transonic VII vs committed experiment.

Pre-registered (cases/analysis/v5_2_rae2822/PRE_REGISTRATION.md, committed
BEFORE the first code change). Question: does the loose VII loop (the
committed GV3.1 recipe with the GV3.2 transonic-point Newton-driver
protocol, NEWTON_ARGS IMPORTED from the GV3 runner) reproduce the
committed RAE2822 experimental Cp?

Bands (pre-registered):
  (a) A4 TE-wedge pre-check (RECORDED-with-fallback): crease-angle wedge
      + linear/quadratic recovery availability on both meshes + the
      ordinate-direct measure; quadratic unavailable -> the TE-band u_e
      falls back to linear+smoothed, recorded.
  (b) shock location (PASS/FAIL, medium binding, coarse recorded):
      x_shock = argmax |dCp/dx| on the computed upper wall Cp vs the
      committed brackets widened +/-0.03c (P1 [0.495, 0.580], P2
      [0.520, 0.605]).
  (c) Cp RMS (RECORDED): per point per side, with the A4 medium input
      band (~2.5 % peak-rel u_e) annotated.
  (d) convergence/guards (RECORDED): loose <= 10 outer, final IBL
      floors, cl/cd_p histories, pre-shock peak Mach vs the M_shock
      <= 1.3 validity envelope.

Conditions = the committed datasets' labeled conditions verbatim:
P1 = (M 0.725, alpha 2.55, Re 6.5e6), P2 = (M 0.73, alpha 3.19,
Re 6.5e6); forced transition x_tr/c = 0.03 both sides.

Threading default 16 (agent-rules); this session ran with 8 (temporary).

Run:  python cases/analysis/v5_2_rae2822/run.py
"""
import importlib.util
import os
import sys
import time
from pathlib import Path

for _var in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_var, "16")

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.meshgen.planar import (  # noqa: E402
    load_airfoil_ordinates,
    te_wedge_angle_deg,
)
from pyfp3d.post.surface import (  # noqa: E402
    triangle_tangential_gradients,
    wall_crease_angles,
    wall_force_coefficients,
    wall_outward_normals,
    wall_tangential_gradient,
    wall_tangential_gradient_quadratic,
)
from pyfp3d.viscous.coupling import (  # noqa: E402
    CouplingConfig,
    build_airfoil_case,
    run_loose_coupling,
)
from pyfp3d.solve.newton import (  # noqa: E402
    solve_newton_lifting,
    solve_newton_transonic,
)


def _load(mod_name, rel):
    spec = importlib.util.spec_from_file_location(mod_name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gv3 = _load("gv3_run", "cases/analysis/v3_loose_coupling/run.py")
NEWTON_ARGS = gv3.NEWTON_ARGS  # IMPORTED, not re-invented (pre-registration)

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
RAE_DIR = ROOT / "cases" / "meshes" / "rae2822_2.5d"
EXP_DIR = ROOT / "cases" / "reference_data" / "rae2822_experiment"

GAMMA = 1.4
RE = 6.5e6
X_TR = 0.03
LEVELS = ("coarse", "medium")
POINTS = {
    "P1": {"m": 0.725, "alpha": 2.55,
           "exp": EXP_DIR / "ExpCase7_RAE2822_M0.725_AoA2.55_Rec6.5e6.dat",
           "bracket": (0.525, 0.55), "band": (0.495, 0.580)},
    "P2": {"m": 0.73, "alpha": 3.19,
           "exp": EXP_DIR / "Expe_RAE2822_M0.73_AoA3.19_Rec6.5e6.dat",
           "bracket": (0.55, 0.575), "band": (0.520, 0.605)},
}

# the A4 medium input band (committed a4_ue_error_band), annotated on (c)
A4_BAND_NOTE = "A4 medium input band: ~2.5% peak-rel u_e; LE band 4-7%"


# ---------------------------------------------------------------------------
# experiment parsing + comparison helpers (unit-tested)
# ---------------------------------------------------------------------------

def load_experiment_cp(path):
    """Parse a committed Tecplot ASCII Cp file -> ((x_up, cp_up), (x_lo,
    cp_lo)), each sorted by ascending x. Handles both committed layouts:
    named "Exp Upper"/"Exp Lower" zones (P1) and a single zone wrapping
    around the LE (P2, split at the x-turning point)."""
    zones = []
    cur = None
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            up = s.upper()
            if up.startswith("ZONE"):
                if cur is not None:
                    zones.append(cur)
                cur = ([], [], s)
                continue
            if not up[0].isdigit() and not up.startswith(("-", "+", ".")):
                continue  # VARIABLES / AUXDATA / TITLE / ... lines
            vals = s.replace(",", " ").split()
            try:
                xv, cv = float(vals[0]), float(vals[1])
            except (ValueError, IndexError):
                continue
            if cur is None:
                cur = ([], [], "")
            cur[0].append(xv)
            cur[1].append(cv)
    if cur is not None:
        zones.append(cur)
    if not zones:
        raise ValueError(f"{path}: no data rows parsed")

    xs_u: list = []
    cp_u: list = []
    xs_l: list = []
    cp_l: list = []
    for xs, cps, title in zones:
        x = np.asarray(xs, dtype=np.float64)
        c = np.asarray(cps, dtype=np.float64)
        tl = title.lower()
        if "upper" in tl:
            pair_u, pair_l = (x, c), None
        elif "lower" in tl:
            pair_u, pair_l = None, (x, c)
        else:
            i = int(np.argmin(x))
            pair_u = (x[: i + 1], c[: i + 1])
            pair_l = (x[i + 1:], c[i + 1:])
        if pair_u is not None:
            xs_u.append(pair_u[0])
            cp_u.append(pair_u[1])
        if pair_l is not None:
            xs_l.append(pair_l[0])
            cp_l.append(pair_l[1])

    def _join(xs, cps):
        x = np.concatenate(xs)
        c = np.concatenate(cps)
        o = np.argsort(x)
        return x[o], c[o]

    return _join(xs_u, cp_u), _join(xs_l, cp_l)


def cp_rms(x_cfd, cp_cfd, x_exp, cp_exp):
    """RMS(Cp_cfd - Cp_exp) at the experimental stations (computed Cp
    linearly interpolated; stations outside the computed range dropped)."""
    x_cfd = np.asarray(x_cfd, dtype=np.float64)
    cp_cfd = np.asarray(cp_cfd, dtype=np.float64)
    m = (x_exp >= x_cfd.min()) & (x_exp <= x_cfd.max())
    d = np.interp(x_exp[m], x_cfd, cp_cfd) - cp_exp[m]
    return float(np.sqrt(np.mean(d * d)))


def shock_x(x, cp, x_lo=0.2, x_hi=0.9):
    """x/c of the shock = argmax of the COMPRESSION gradient (positive
    dCp/dx) inside the mid-chord window [x_lo, x_hi] (pre-registered
    addendum 2026-07-24: the unwindowed max |dCp/dx| picks the LE
    suction spike, not the shock)."""
    x = np.asarray(x, dtype=np.float64)
    cp = np.asarray(cp, dtype=np.float64)
    g = np.gradient(cp, x)
    m = (x >= x_lo) & (x <= x_hi)
    gw = np.where(m, g, -np.inf)
    return float(x[int(np.argmax(gw))])


# ---------------------------------------------------------------------------
# wall Cp curves + peak Mach
# ---------------------------------------------------------------------------

def _wall_cp_sides(mc, phi, m_inf, alpha, s_ref):
    """Per-side binned wall Cp curves ((x, cp) upper, (x, cp) lower),
    averaged per (rounded-x, side) bin from the per-triangle Cp.

    Side split by the outward-normal y sign (the D11 idiom, surface_ls.py):
    robust on reflex-cambered sections whose aft lower surface sits above
    the chord line (RAE2822), where a centroid-y split would mislabel the
    aft lower triangles as upper."""
    wf = mc.boundary_faces["wall"]
    f = wall_force_coefficients(mc.nodes, mc.elements, wf, phi,
                                alpha_deg=alpha, s_ref=s_ref, m_inf=m_inf)
    cp = f["cp_tri"]
    cen = mc.nodes[wf].mean(axis=1)
    side_up = wall_outward_normals(mc.nodes, mc.elements, wf)[:, 1] > 0.0
    out = []
    for mask in (side_up, ~side_up):
        xr = np.round(cen[mask, 0], 4)
        keys, inv = np.unique(xr, return_inverse=True)
        acc = np.zeros(keys.size)
        cnt = np.zeros(keys.size)
        np.add.at(acc, inv, cp[mask])
        np.add.at(cnt, inv, 1.0)
        out.append((keys, acc / np.maximum(cnt, 1.0)))
    return out[0], out[1]


def _peak_mach(mc, phi, m_inf, x_max):
    """Pre-shock peak Mach on the upper wall (total-enthalpy relation,
    u_inf = 1 -> a_inf = 1/M_inf), over triangles with centroid
    x <= x_max; returns (M_peak, x_at_peak). Upper-side selection by the
    outward-normal y sign (reflex-camber robust, see _wall_cp_sides)."""
    wf = mc.boundary_faces["wall"]
    grad_tri, _, _ = triangle_tangential_gradients(mc.nodes, wf, phi)
    q2 = np.sum(grad_tri * grad_tri, axis=1)
    cen = mc.nodes[wf].mean(axis=1)
    up = wall_outward_normals(mc.nodes, mc.elements, wf)[:, 1] > 0.0
    mask = up & (cen[:, 0] <= x_max)
    a0_sq = 1.0 / m_inf**2 + 0.5 * (GAMMA - 1.0)
    a2 = a0_sq - 0.5 * (GAMMA - 1.0) * q2[mask]
    a2 = np.maximum(a2, 1e-12)
    M = np.sqrt(q2[mask] / a2)
    i = int(np.argmax(M))
    return float(M[i]), float(cen[mask, 0][i])


# ---------------------------------------------------------------------------
# band (a): the A4 TE-wedge pre-check
# ---------------------------------------------------------------------------

def _te_precheck(level, mesh):
    """A4-method TE structural pre-check, on the UNCUT mesh (the A4
    runner measures crease/gradient availability pre-cut; on the cut mesh
    the TE strips no longer share the TE edge and the crease measure is
    meaningless)."""
    wf = mesh.boundary_faces["wall"]
    ang, _ = wall_crease_angles(mesh.nodes, mesh.elements, wf)
    wedge_mesh = 180.0 - float(np.nanmax(ang))
    phi_dummy = mesh.nodes[:, 0].copy()
    lin_ok, quad_ok = True, True
    lin_note = quad_note = "ok"
    try:
        wall_tangential_gradient(mesh.nodes, wf, phi_dummy)
    except ValueError as e:
        lin_ok, lin_note = False, str(e)[:60]
    try:
        wall_tangential_gradient_quadratic(mesh.nodes, wf, phi_dummy)
    except ValueError as e:
        quad_ok, quad_note = False, str(e)[:60]
    x, z_lo, z_up = load_airfoil_ordinates(RAE_DIR / "rae2822.dat")
    wedge_ord = te_wedge_angle_deg(x, z_lo, z_up)
    return {
        "level": level,
        "te_wedge_deg_mesh": wedge_mesh,
        "te_wedge_deg_ordinates": wedge_ord,
        "linear_available": lin_ok,
        "quadratic_available": quad_ok,
        "quadratic_note": quad_note,
        "fallback_linear_smoothed": (not quad_ok),
    }


# ---------------------------------------------------------------------------
# one (level, point) VII leg
# ---------------------------------------------------------------------------

def _make_transonic_driver(mc, wc, m_inf, alpha):
    """GV3.2 Newton driver + the committed transonic rescue chain
    (pre-registration addenda 2026-07-24 #2 and #3). Per FP call, ordered
    cheap -> deep, FIRST success wins:

    warm-started (k >= 1): single-shot strict -> single-shot
    stall-accept -> cold continuation strict -> cold continuation
    stall-accept -> the GV3.3 loud raise.
    cold start (k = 0):    continuation strict -> continuation
    stall-accept -> raise.

    Continuation = `solve_newton_transonic` upward Mach ramp from
    m_start = 0.70 (the library's designated transonic path, strict
    final level at NEWTON_ARGS' tol_residual). Stall-accept =
    `accept_on_stall=True` -- the library's honesty-guarded plateau
    acceptance (accept_reason "stall" only under f_norm < tol_gamma,
    zero upwind-limiter/floor activity, and the live_stalled plateau
    rule), for the |R| ~ 1e-9 shock-cell plateaus the strict 1e-10
    cannot reach; every attempt's (path, accept_reason, converged) is
    logged for the summary."""
    log = []

    def _single(rhs, seed, stall):
        kw = dict(NEWTON_ARGS)
        if stall:
            kw["accept_on_stall"] = True
        r = solve_newton_lifting(
            mc, wc, m_inf=m_inf, alpha_deg=alpha, external_rhs=rhs,
            phi_init=seed.phi, gamma_init=seed.gamma, **kw)
        log.append(("single_stall_accept" if stall else "single_strict",
                    str(r.get("accept_reason")), bool(r["converged"])))
        return r

    def _cont(rhs, stall):
        kw = dict(external_rhs=rhs, **NEWTON_ARGS)
        if stall:
            kw["accept_on_stall"] = True
        r = solve_newton_transonic(
            mc, wc, m_inf=m_inf, alpha_deg=alpha, newton_kw=kw)
        log.append(("continuation_stall_accept" if stall
                    else "continuation_strict",
                    str(r.get("accept_reason")), bool(r["converged"])))
        return r

    def solve(rhs, seed):
        attempts = (
            [lambda: _single(rhs, seed, False),
             lambda: _single(rhs, seed, True),
             lambda: _cont(rhs, False),
             lambda: _cont(rhs, True)] if seed is not None else
            [lambda: _cont(rhs, False), lambda: _cont(rhs, True)])
        r = attempts[0]()
        for attempt in attempts[1:]:
            if r["converged"]:
                break
            r = attempt()
        return r["phi"], r["gamma"], r

    return solve, log


def _run_point(level, mc, wc, pname, pt):
    cfg = CouplingConfig(re_chord=RE, m_inf=pt["m"], alpha_deg=pt["alpha"],
                         x_tr_upper=X_TR, x_tr_lower=X_TR)
    case = build_airfoil_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg)
    dz = float(np.ptp(mc.nodes[:, 2]))
    s_ref = 1.0 * dz

    def probe(phi, gamma, k):
        f = wall_force_coefficients(
            mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
            alpha_deg=pt["alpha"], s_ref=s_ref, m_inf=pt["m"])
        return {"cl": f["cl"], "cd_p": f["cd_pressure"]}

    driver, path_log = _make_transonic_driver(mc, wc, pt["m"], pt["alpha"])
    t0 = time.perf_counter()
    res = run_loose_coupling(driver, case, cfg, probe=probe)
    wall = time.perf_counter() - t0
    n_stall = sum(1 for _, a, c in path_log if a == "stall" and c)
    n_cont = sum(1 for p, _, _ in path_log if p.startswith("continuation"))
    print(f"    fp calls: {len(path_log)} ({n_cont} continuation, "
          f"{n_stall} stall-accepted)", flush=True)

    up, lo = _wall_cp_sides(mc, res.phi, pt["m"], pt["alpha"], s_ref)
    xs = shock_x(up[0], up[1])
    lo_b, hi_b = pt["band"]
    in_band = bool(lo_b <= xs <= hi_b)
    m_peak, x_peak = _peak_mach(mc, res.phi, pt["m"], xs)

    exp_up, exp_lo = load_experiment_cp(pt["exp"])
    rms_up = cp_rms(up[0], up[1], exp_up[0], exp_up[1])
    rms_lo = cp_rms(lo[0], lo[1], exp_lo[0], exp_lo[1])

    # artifacts
    hist = res.history
    pd.DataFrame(hist).to_csv(
        RESULTS / f"convergence_{pname}_{level}.csv", index=False)
    rows = []
    for side, (xe, ce), (xc, cc) in (
            ("upper", exp_up, up), ("lower", exp_lo, lo)):
        m = (xe >= xc.min()) & (xe <= xc.max())
        rows.append(pd.DataFrame({
            "side": side, "x_c": xe[m], "cp_exp": ce[m],
            "cp_cfd": np.interp(xe[m], xc, cc),
            "diff": np.interp(xe[m], xc, cc) - ce[m],
        }))
    pd.concat(rows).to_csv(
        RESULTS / f"cp_compare_{pname}_{level}.csv", index=False)

    return {
        "level": level, "point": pname,
        "m_inf": pt["m"], "alpha": pt["alpha"],
        "converged": bool(res.converged), "n_outer": int(res.n_outer),
        "wall_s": wall,
        "fp_calls": int(len(path_log)),
        "fp_continuation": int(n_cont),
        "fp_stall_accepted": int(n_stall),
        "cl_final": float(hist[-1]["cl"]) if hist else float("nan"),
        "ibl_final_residual": float(hist[-1]["ibl_final_residual"])
        if hist and "ibl_final_residual" in hist[-1] else float("nan"),
        "x_shock": xs, "band_lo": lo_b, "band_hi": hi_b,
        "shock_in_band": in_band,
        "m_peak_preshock": m_peak, "x_at_peak": x_peak,
        "outside_envelope": bool(m_peak > 1.3),
        "cp_rms_upper": rms_up, "cp_rms_lower": rms_lo,
        "failure": "",
    }


def _failure_row(level, pname, pt, exc):
    """Pre-registered §6 recipe-limit clause: an FP-driver failure at a
    point reads as that point RECORDED; the remaining legs still read.
    The summary row carries the failure note; metrics stay NaN."""
    return {
        "level": level, "point": pname,
        "m_inf": pt["m"], "alpha": pt["alpha"],
        "converged": False, "n_outer": -1,
        "wall_s": float("nan"), "fp_calls": -1,
        "fp_continuation": -1, "fp_stall_accepted": -1,
        "cl_final": float("nan"), "ibl_final_residual": float("nan"),
        "x_shock": float("nan"),
        "band_lo": pt["band"][0], "band_hi": pt["band"][1],
        "shock_in_band": False,
        "m_peak_preshock": float("nan"), "x_at_peak": float("nan"),
        "outside_envelope": False,
        "cp_rms_upper": float("nan"), "cp_rms_lower": float("nan"),
        "failure": str(exc)[:100],
    }


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    print(f"GV5.2 NEWTON_ARGS imported from the GV3 runner: {NEWTON_ARGS}")
    te_rows = []
    rows = []
    for level in LEVELS:
        print(f"\n==================== GV5.2 level={level} ====================")
        mesh = read_mesh(RAE_DIR / f"{level}.msh")
        te = _te_precheck(level, mesh)
        te_rows.append(te)
        print(f"  band (a): wedge mesh {te['te_wedge_deg_mesh']:.2f} deg / "
              f"ordinates {te['te_wedge_deg_ordinates']:.2f} deg; "
              f"quadratic_available={te['quadratic_available']} "
              f"({te['quadratic_note']})")
        mc, wc = cut_wake(mesh)
        del mesh
        for pname, pt in POINTS.items():
            print(f"  --- {level} {pname}: M={pt['m']} alpha={pt['alpha']} "
                  f"Re={RE:.2e} x_tr={X_TR} ---", flush=True)
            try:
                r = _run_point(level, mc, wc, pname, pt)
            except RuntimeError as e:
                r = _failure_row(level, pname, pt, e)
                print(f"    POINT RECORDED (§6 recipe-limit clause): {e}",
                      flush=True)
            rows.append(r)
            print(f"    converged={r['converged']} outer={r['n_outer']} "
                  f"wall={r['wall_s']:.0f}s x_shock={r['x_shock']:.4f} "
                  f"band=[{r['band_lo']},{r['band_hi']}] "
                  f"in_band={r['shock_in_band']} "
                  f"M_peak={r['m_peak_preshock']:.3f} "
                  f"rms_up={r['cp_rms_upper']:.4f} "
                  f"rms_lo={r['cp_rms_lower']:.4f}", flush=True)
        del mc, wc
        import gc
        gc.collect()

    pd.DataFrame(te_rows).to_csv(RESULTS / "te_wedge.csv", index=False)
    pd.DataFrame(rows).to_csv(RESULTS / "summary.csv", index=False)
    shock_rows = [{k: r[k] for k in ("level", "point", "x_shock", "band_lo",
                                     "band_hi", "shock_in_band",
                                     "m_peak_preshock", "outside_envelope")}
                  for r in rows]
    for level in LEVELS:
        pd.DataFrame([r for r in shock_rows if r["level"] == level]).to_csv(
            RESULTS / f"shock_{level}.csv", index=False)

    # verdict (pre-registered): band (b) PASS = both points in band, medium
    med = [r for r in rows if r["level"] == "medium"]
    b_pass = all(r["shock_in_band"] for r in med) and len(med) == 2
    print("\n--- summary ---")
    for r in rows:
        print(r)
    print(f"\nband (b) medium binding: "
          f"{'PASS' if b_pass else 'FAIL (see clauses)'}")


if __name__ == "__main__":
    main()
