"""GV3.1/GV3.2 loose FP+IBL coupling gates (Track V3).

Binding text: docs/roadmap/track_v.md GV3.1/GV3.2 (2026-07-22 re-spec);
pre-registered bands: cases/analysis/v3_loose_coupling/PRE_REGISTRATION.md
(committed before the first execution). Regenerates every CSV/PNG in
results/ and exits 0 iff all pre-registered assertions hold (honest FAIL
otherwise).

  GV3.1: per-side station delta*/c and C_f vs the committed XFOIL reference
         (banded window x/c in (0.05, 0.95], |err| <= 25 % / 15 % at EVERY
         station, medium level binding, coarse RECORDED) + viscous Delta-cl
         (direction > 0 AND magnitude within a factor 2 of XFOIL's own
         inviscid->viscous decrement);
  GV3.2: loose loop converges (raw successive delta* max rel change < 1e-3)
         in <= 10 outer iterations at the recorded omega (executed at
         omega = 1.0) + one RECORDED transonic-attached point (M 0.72,
         alpha 2, coarse, conforming-Newton driver).

Run:  python cases/analysis/v3_loose_coupling/run.py
"""

import os
import sys

# resolve pyfp3d from THIS worktree (the site-packages editable install may
# point at a sibling worktree)
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..")))

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.physics.isentropic import density_field
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    build_airfoil_case,
    make_newton_lifting_driver,
    make_picard_lifting_driver,
    run_loose_coupling,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
NACA_DIR = os.path.join(REPO, "cases", "meshes", "naca0012_2.5d")
XFOIL_DIR = os.path.join(REPO, "cases", "reference_data",
                         "naca0012_viscous_xfoil")

M_INF, ALPHA, RE = 0.5, 2.0, 3.0e6
M_TRANSONIC = 0.72
# Newton leg settings (pre-registered, identical to GV2.1's CASE_ARGS)
NEWTON_ARGS = dict(upwind_c=1.5, m_crit=0.95, m_cap=3.0, rho_floor=0.05,
                   tol_residual=1e-10)

# windows (pre-registered)
BAND_LO, BAND_HI = 0.05, 0.95
TOL_DS, TOL_CF = 0.25, 0.15      # per-station |rel err| bands (medium)
DCL_FACTOR = 2.0                  # Delta-cl magnitude band vs XFOIL's own

SUMMARY = []  # (gate, metric, band, measured, verdict)


def _record(gate, metric, band, measured, ok):
    verdict = "RECORDED" if ok is None else ("PASS" if ok else "FAIL")
    SUMMARY.append((gate, metric, band, measured, verdict))
    print(f"  [{verdict:8s}] {gate} {metric}: measured={measured} "
          f"(band: {band})", flush=True)


def _write_csv(name, header, rows):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write(",".join(str(x) for x in row) + "\n")
    print(f"  wrote {path}", flush=True)


# ---------------------------------------------------------------------------
# XFOIL reference
# ---------------------------------------------------------------------------

def load_xfoil():
    """Committed reference: per-surface delta*/c and C_f profiles + the
    polar/inviscid summaries. Returns dict with sorted per-surface arrays
    and the scalar decrements."""
    prof = {"upper": {"x": [], "ds": [], "cf": []},
            "lower": {"x": [], "ds": [], "cf": []}}
    with open(os.path.join(
            XFOIL_DIR, "delta_star_cf_alpha2_m05_xtr005.csv")) as f:
        for row in csv.DictReader(f):
            s = prof[row["surface"]]
            s["x"].append(float(row["x_c"]))
            s["ds"].append(float(row["dstar_over_c"]))
            s["cf"].append(float(row["cf"]))
    for s in prof.values():
        idx = np.argsort(s["x"])
        for k in s:
            s[k] = np.asarray(s[k])[idx]

    with open(os.path.join(XFOIL_DIR, "polar_summary.csv")) as f:
        polar = {r["case"]: r for r in csv.DictReader(f)}
    with open(os.path.join(XFOIL_DIR, "inviscid_summary.csv")) as f:
        inv = next(csv.DictReader(f))
    cl_inv = float(inv["cl"])
    cl_visc = float(polar["xtr005"]["cl"])
    return {
        "prof": prof,
        "cl_inviscid": cl_inv,
        "cl_viscous": cl_visc,
        "dcl": cl_inv - cl_visc,
        "cl_viscous_xtr030": float(polar["xtr030"]["cl"]),
    }


# ---------------------------------------------------------------------------
# one loose-coupling run (one mesh level / Mach / driver kind)
# ---------------------------------------------------------------------------

def run_case(level, m_inf, driver_kind):
    print(f"--- case: {level} M={m_inf} alpha={ALPHA} driver={driver_kind} "
          f"---", flush=True)
    mc, wc = cut_wake(read_mesh(os.path.join(NACA_DIR, f"{level}.msh")))
    cfg = CouplingConfig(re_chord=RE, m_inf=m_inf, alpha_deg=ALPHA)
    case = build_airfoil_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg)
    dz = float(np.ptp(mc.nodes[:, 2]))
    s_ref = 1.0 * dz  # chord x span thickness (quasi-2D sectional loads)

    def probe(phi, gamma, k):
        f = wall_force_coefficients(
            mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
            alpha_deg=ALPHA, s_ref=s_ref, m_inf=m_inf)
        return {"cl": f["cl"], "cd_p": f["cd_pressure"]}

    if driver_kind == "picard":
        driver = make_picard_lifting_driver(mc, wc, m_inf, ALPHA)
    else:
        driver = make_newton_lifting_driver(mc, wc, m_inf, ALPHA,
                                            **NEWTON_ARGS)
    res = run_loose_coupling(driver, case, cfg, probe=probe)
    print(f"    converged={res.converged} n_outer={res.n_outer} "
          f"omega={cfg.omega}", flush=True)
    return {"level": level, "m_inf": m_inf, "driver": driver_kind,
            "case": case, "cfg": cfg, "res": res, "s_ref": s_ref,
            "dz": dz}


# ---------------------------------------------------------------------------
# GV3.1 profile comparison
# ---------------------------------------------------------------------------

def side_station_profiles(case, res, m_inf):
    """Per-side station rows: (xc, side_name, ds, cf_fs, cf_loc) with the
    station's span copies of that side averaged (TE copies split by
    side_node).

    cf is compared in the FREESTREAM reference frame (XFOIL's DUMP
    convention, xoper.f CF = TAU/(0.5*QINF**2)): our closure's OUT_CF1 is
    LOCAL-normalized (D13 (61), cf = 2 tau/(rho_e u_e^2)), so it is
    converted per node BEFORE averaging -- cf_fs = cf_loc * rho_e * |u_e|^2
    (u_inf = rho_inf = 1). See the PRE_REGISTRATION addendum 2026-07-22.
    """
    st = case.stations
    ds = res.outs[:, C.OUT_DS1]
    cf_loc = res.outs[:, C.OUT_CF1]
    q2 = np.sum(res.ue_surf ** 2, axis=1)
    rho_e = density_field(q2, m_inf, 1.4)
    cf_fs = cf_loc * rho_e * q2
    rows = []
    for side_val, side_name in ((1, "upper"), (-1, "lower")):
        for r in range(len(st.xc)):
            mask = (st.station_of == r) & (st.side_node == side_val)
            if not np.any(mask):
                continue
            rows.append((float(st.xc[r]), side_name,
                         float(np.mean(ds[mask])),
                         float(np.mean(cf_fs[mask])),
                         float(np.mean(cf_loc[mask]))))
    return rows


def window_of(xc, pin_x):
    if xc <= pin_x:
        return "pinned"   # Dirichlet inflow band: boundary data, not solution
    if xc <= BAND_LO:
        return "le"
    if xc > BAND_HI:
        return "te"
    return "banded"


def compare_profiles(run, xf, gate, binding):
    """Interpolate XFOIL onto our station x/c per side; per-station relative
    errors; write the CSV; record the windowed metrics.

    Binding metric: cf in the freestream frame (see side_station_profiles).
    The as-registered RAW local-frame comparison is kept as RECORDED rows
    (the pre-registration addendum documents the normalization fix).
    Stations inside the Dirichlet inflow band (x/c <= inflow_band_x) are
    labeled "pinned" and excluded from every error statistic -- they are
    prescribed boundary data (per-node Blasius states), not solution."""
    st_rows = side_station_profiles(run["case"], run["res"], run["m_inf"])
    pin_x = run["cfg"].inflow_band_x
    csv_rows, worst, worst_loc = [], {}, {}
    for xc, side, ds, cf, cf_l in st_rows:
        p = xf["prof"][side]
        ds_x = float(np.interp(xc, p["x"], p["ds"]))
        cf_x = float(np.interp(xc, p["x"], p["cf"]))
        e_ds = (ds - ds_x) / ds_x if ds_x > 0.0 else np.nan
        e_cf = (cf - cf_x) / cf_x if cf_x > 0.0 else np.nan
        e_cl = (cf_l - cf_x) / cf_x if cf_x > 0.0 else np.nan
        win = window_of(xc, pin_x)
        csv_rows.append((f"{xc:.6f}", side, f"{ds:.6e}", f"{ds_x:.6e}",
                         f"{e_ds:+.4f}", f"{cf:.6e}", f"{cf_x:.6e}",
                         f"{e_cf:+.4f}", f"{cf_l:.6e}", f"{e_cl:+.4f}", win))
        if win == "banded":
            for metric, e, table in (("ds", e_ds, worst),
                                     ("cf", e_cf, worst),
                                     ("cf_local", e_cl, worst_loc)):
                key = (metric, side)
                if np.isfinite(e) and (
                        key not in table or abs(e) > abs(table[key][0])):
                    table[key] = (e, xc)
    _write_csv(f"gv3_1_profiles_{run['level']}.csv",
               "x_c,side,ds_ours,ds_xfoil,err_ds,cf_ours,cf_xfoil,err_cf,"
               "cf_ours_localframe,err_cf_localframe,window", csv_rows)

    for (metric, side), (e, xc) in sorted(worst.items()):
        band = TOL_DS if metric == "ds" else TOL_CF
        name = f"{'delta*' if metric == 'ds' else 'cf'} {side} worst " \
               f"|err| (x/c={xc:.3f})"
        if binding:
            _record(gate, name, f"<= {band:.0%} every banded station",
                    f"{abs(e):.3f}", abs(e) <= band)
        else:
            _record(gate, f"{run['level']} {name}", "recorded",
                    f"{abs(e):.3f}", None)
    # as-registered raw local-frame cf (superseded; kept for transparency)
    for (metric, side), (e, xc) in sorted(worst_loc.items()):
        _record(gate, f"{run['level']} cf LOCAL-frame {side} worst |err| "
                f"(x/c={xc:.3f}) [as-registered, superseded]",
                "recorded", f"{abs(e):.3f}", None)

    # RECORDED: LE / near-TE window max |err| (context, not gated)
    for win in ("le", "te"):
        es = [abs(float(r[4])) for r in csv_rows
              if r[10] == win and np.isfinite(float(r[4]))]
        if es:
            _record(gate, f"{run['level']} delta* max|err| {win}-zone",
                    "recorded", f"{max(es):.3f}", None)
    return csv_rows


def gate_delta_cl(run, xf, binding):
    """Delta-cl = cl_inviscid(k=0 probe) - cl_coupled(final). Direction > 0
    and magnitude within a factor 2 of XFOIL's own decrement."""
    hist = run["res"].history
    cl_inv = float(hist[0]["cl"])
    cl_fin = float(hist[-1]["cl"])
    dcl = cl_inv - cl_fin
    ratio = dcl / xf["dcl"] if xf["dcl"] > 0.0 else np.nan
    tag = "" if binding else f"{run['level']} "
    _record("GV3.1", f"{tag}cl_inviscid (k=0)", "recorded",
            f"{cl_inv:.4f}", None)
    _record("GV3.1", f"{tag}cl_coupled (final)", "recorded",
            f"{cl_fin:.4f}", None)
    _record("GV3.1", f"{tag}Delta-cl vs XFOIL own "
            f"({xf['dcl']:.4f})", "recorded", f"{dcl:.4f}", None)
    ok = bool(dcl > 0.0 and 1.0 / DCL_FACTOR <= ratio <= DCL_FACTOR)
    if binding:
        _record("GV3.1", "Delta-cl direction + factor-2 band",
                f"0 < dcl, 0.5 <= ratio <= {DCL_FACTOR}",
                f"ratio={ratio:.3f}", ok)
    else:
        _record("GV3.1", f"{run['level']} Delta-cl ratio", "recorded",
                f"{ratio:.3f}", None)
    return cl_inv, cl_fin, dcl


def write_history(run, name):
    head = ("k,ds_max,ds_change_rel,ds_neg_floored,mdot_max,ibl_n_iter,"
            "ibl_converged,ibl_final_residual,inflow_n_pinned,cl,cd_p")
    rows = []
    for h in run["res"].history:
        rows.append(tuple(h.get(k, "") for k in head.split(",")))
    _write_csv(name, head, rows)


# ---------------------------------------------------------------------------
# panels
# ---------------------------------------------------------------------------

def _panel_profiles(runs, xf):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    style = {"medium": "-", "coarse": "--"}
    for run in runs:
        pin_x = run["cfg"].inflow_band_x
        rows = [r for r in side_station_profiles(
            run["case"], run["res"], run["m_inf"]) if r[0] > pin_x]
        for side in ("upper", "lower"):
            sel = [(x, d, c) for x, s, d, c, _ in rows if s == side]
            x = [s[0] for s in sel]
            col = 0 if side == "upper" else 1
            axes[0, col].plot(x, [s[1] for s in sel], style[run["level"]],
                              label=f"ours {run['level']}")
            axes[1, col].plot(x, [s[2] for s in sel], style[run["level"]],
                              label=f"ours {run['level']}")
    for col, side in enumerate(("upper", "lower")):
        p = xf["prof"][side]
        axes[0, col].plot(p["x"], p["ds"], "k-", lw=1.6, label="XFOIL")
        axes[1, col].plot(p["x"], p["cf"], "k-", lw=1.6, label="XFOIL")
        axes[0, col].set_ylabel("delta*/c")
        axes[1, col].set_ylabel("C_f")
        axes[1, col].set_xlabel("x/c")
        for ax in (axes[0, col], axes[1, col]):
            ax.axvspan(BAND_LO, BAND_HI, color="g", alpha=0.06)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        axes[0, col].set_title(f"GV3.1 delta*/c {side} (M0.5 a2 Re3e6 "
                               f"xtr0.05; banded window shaded)")
    fig.tight_layout()
    path = os.path.join(RESULTS, "gv3_1_profiles.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}", flush=True)


def _panel_cl(runs, xf):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for run in runs:
        h = run["res"].history
        ks = [r["k"] for r in h]
        cls = [r.get("cl", np.nan) for r in h]
        ax.plot(ks, cls, "o-", label=f"ours {run['level']}")
    ax.axhline(xf["cl_inviscid"], color="k", ls=":", label="XFOIL inviscid")
    ax.axhline(xf["cl_viscous"], color="k", ls="--",
               label="XFOIL viscous xtr005")
    ax.set_xlabel("outer iteration k")
    ax.set_ylabel("c_l (pressure integral)")
    ax.set_title("GV3.1/GV3.2 cl history of the loose loop")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(RESULTS, "gv3_1_cl_history.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}", flush=True)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    xf = load_xfoil()
    print(f"XFOIL reference: cl_inviscid={xf['cl_inviscid']:.4f} "
          f"cl_viscous(xtr005)={xf['cl_viscous']:.4f} "
          f"Delta-cl={xf['dcl']:.4f}", flush=True)

    # -- GV3.1 medium (binding) ---------------------------------------------
    run_m = run_case("medium", M_INF, "picard")
    write_history(run_m, "gv3_1_history_medium.csv")
    compare_profiles(run_m, xf, "GV3.1", binding=True)
    gate_delta_cl(run_m, xf, binding=True)

    # GV3.2 loop convergence (from the same medium run)
    res_m = run_m["res"]
    _record("GV3.2", "loose loop converged <= 10 outer iters",
            "<= 10 (omega recorded)", f"n_outer={res_m.n_outer} "
            f"omega={run_m['cfg'].omega} converged={res_m.converged}",
            bool(res_m.converged and res_m.n_outer <= 10))
    floor = max(h.get("ibl_final_residual", 0.0) for h in res_m.history[1:])
    _record("GV3.2", "IBL residual floor over outer iters", "recorded",
            f"{floor:.3e}", None)

    # -- GV3.1 coarse cross-check (RECORDED) --------------------------------
    run_c = run_case("coarse", M_INF, "picard")
    write_history(run_c, "gv3_1_history_coarse.csv")
    compare_profiles(run_c, xf, "GV3.1", binding=False)
    gate_delta_cl(run_c, xf, binding=False)

    # -- GV3.2 transonic recorded point (coarse, Newton driver) -------------
    run_t = run_case("coarse", M_TRANSONIC, "newton")
    write_history(run_t, "gv3_2_transonic_history.csv")
    res_t = run_t["res"]
    _record("GV3.2", "transonic M0.72 coarse Newton: outer count",
            "recorded", f"n_outer={res_t.n_outer} "
            f"converged={res_t.converged} omega={run_t['cfg'].omega}",
            None)
    floor_t = max(h.get("ibl_final_residual", 0.0)
                  for h in res_t.history[1:])
    _record("GV3.2", "transonic IBL residual floor", "recorded",
            f"{floor_t:.3e}", None)
    _record("GV3.2", "transonic cl final", "recorded",
            f"{res_t.history[-1].get('cl', float('nan')):.4f}", None)

    # -- panels + summary ----------------------------------------------------
    _panel_profiles([run_m, run_c], xf)
    _panel_cl([run_m, run_c], xf)
    _write_csv("summary.csv", "gate,metric,band,measured,verdict", SUMMARY)
    n_fail = sum(1 for *_, v in SUMMARY if v == "FAIL")
    n_pass = sum(1 for *_, v in SUMMARY if v == "PASS")
    n_rec = sum(1 for *_, v in SUMMARY if v == "RECORDED")
    print(f"\nGV3.1/GV3.2: {n_pass} PASS / {n_fail} FAIL / {n_rec} RECORDED",
          flush=True)
    if n_fail:
        print("HONEST FAIL -- see summary.csv", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
