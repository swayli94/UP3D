"""GV5.0 M6 subsonic loose-coupling bridge (Track V5 entry check).

Binding text: docs/roadmap/track_v.md GV5.0 (2026-07-22, user-directed);
pre-registered: cases/analysis/v5_m6_bridge/PRE_REGISTRATION.md (committed
before the first execution). ALL metrics RECORDED (no pass/fail; no
delta*(z) truth data exists for the M6). Regenerates every CSV/PNG in
results/ and exits 0 unless a hard guardrail fires (FP non-convergence /
non-finite state -- those raise, GV3.3 discipline).

  ONERA M6 (coarse + medium), conforming path, M0.5 / alpha 3.06, forced
  transition x_tr/c = 0.05 both sides, Re_MAC = 11.72e6. Outputs: delta*(z)
  spanwise distributions at fixed x/c stations (feeds GV5.3's band
  pre-registration), the crossflow field (first live 3-D exercise of the
  Psi, B equations), Delta-CL viscous-inviscid (direction recorded), and
  the 3-D loose-loop iteration count vs the GV3.2 2.5-D count (4-5). The
  wing-tip band (z > 0.95*b_semi, the production tip_taper radius) is
  pinned + masked.

Run:  python cases/analysis/v5_m6_bridge/run.py [--levels coarse medium]
"""

import argparse
import os
import sys
import time

# resolve pyfp3d from THIS worktree (the site-packages editable install may
# point at a sibling worktree)
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..")))

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pyfp3d.constraints.wake import tip_taper_factors
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at, x_le
from pyfp3d.post.surface import (
    cl_kj_3d,
    planform_area,
    wall_force_coefficients,
)
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    build_wing_case,
    run_loose_coupling,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
M6_DIR = os.path.join(REPO, "cases", "meshes", "onera_m6")

M_INF, ALPHA = 0.5, 3.06
MAC = 0.64607
RE_MAC = 11.72e6
RE_CHORD = RE_MAC / MAC            # per meter (meshes are in NASA meters)
X_TR = 0.05                        # forced transition, both sides
TIP_FRAC = 0.05                    # tip mask = production tip_taper radius
XC_STATIONS = (0.2, 0.4, 0.6, 0.8)
XC_WIN = 0.03                      # station half-window in local x/c
Z_BIN = 0.05                       # spanwise bin width, fraction of b_semi
GV32_25D_OUTER = "4 (coarse) / 5 (medium)"   # the 2.5-D reference count

# production M6 subsonic Newton recipe (P14 demo, tier 1)
M6_NEWTON_KW = dict(farfield_spanwise_gamma=True, precond="direct",
                    direct_refactor_every=1000, n_newton_max=60)

SUMMARY = []  # (gate, metric, band, measured, verdict)


def _record(gate, metric, band, measured, ok=None):
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
# FP driver: production M6 Newton + pressure Kutta + tip taper; the k=0 cold
# start chains probe -> pressure seeding (P14 recipe), outer iterations
# warm-start the pressure solve.
# ---------------------------------------------------------------------------

def make_m6_driver(mc, wc, taper):
    from pyfp3d.solve.newton import solve_newton_lifting

    def solve(rhs, seed):
        kw = dict(M6_NEWTON_KW, tip_taper=taper, external_rhs=rhs)
        if seed is None or seed.phi is None:
            r0 = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                      **kw)
            if not r0["converged"]:
                return r0["phi"], r0.get("gamma"), r0
            r = solve_newton_lifting(
                mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                kutta_estimator="pressure",
                phi_init=r0["phi"], gamma_init=r0["gamma"],
                n_picard_seed=0, **kw)
        else:
            r = solve_newton_lifting(
                mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                kutta_estimator="pressure",
                phi_init=seed.phi, gamma_init=seed.gamma,
                n_picard_seed=0, **kw)
        return r["phi"], r["gamma"], r

    return solve


# ---------------------------------------------------------------------------
# one level
# ---------------------------------------------------------------------------

def run_level(level):
    print(f"--- GV5.0 {level}: ONERA M6 M={M_INF} alpha={ALPHA} "
          f"Re_MAC={RE_MAC:.3e} ---", flush=True)
    mc, wc = cut_wake(read_mesh(os.path.join(M6_DIR, f"{level}.msh")))
    wall = mc.boundary_faces["wall"]
    taper = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth",
                              TIP_FRAC * B_SEMI)
    cfg = CouplingConfig(re_chord=RE_CHORD, m_inf=M_INF, alpha_deg=ALPHA,
                         x_tr_upper=X_TR, x_tr_lower=X_TR)
    case = build_wing_case(mc.nodes, mc.elements, wall, cfg,
                           x_le=x_le, chord_at=chord_at,
                           tip_mask_frac=TIP_FRAC)
    sm = case.sm
    print(f"    IBL surface: {sm.n_node} nodes / {sm.n_tri} tris; "
          f"tip-masked {int(case.outflow_pin_surf.sum())}; "
          f"LE-band {int(case.inflow_candidates.sum())}", flush=True)

    s_ref = planform_area(mc.nodes, wall)
    o = np.argsort(wc.station_z)

    def probe(phi, gamma, k):
        f = wall_force_coefficients(mc.nodes, mc.elements, wall, phi,
                                    alpha_deg=ALPHA, s_ref=s_ref,
                                    m_inf=M_INF)
        cl_kj = cl_kj_3d(np.asarray(gamma)[o], wc.station_z[o], s_ref,
                         B_SEMI)
        return {"cl_p": f["cl"], "cl_kj": float(cl_kj),
                "cd_p": f["cd_pressure"]}

    driver = make_m6_driver(mc, wc, taper)
    t0 = time.perf_counter()
    res = run_loose_coupling(driver, case, cfg, probe=probe)
    wall_s = time.perf_counter() - t0
    print(f"    converged={res.converged} n_outer={res.n_outer} "
          f"omega={cfg.omega} wall={wall_s:.0f}s", flush=True)

    # -- per-node surface tables ---------------------------------------------
    x, y, z = sm.xyz[:, 0], sm.xyz[:, 1], sm.xyz[:, 2]
    zc = np.clip(z, 0.0, B_SEMI)
    chord_n = chord_at(zc)
    xc_n = (x - x_le(zc)) / chord_n
    zb = z / B_SEMI
    cent_y = sm.xyz[sm.triangles].mean(axis=1)[:, 1]
    ysum = np.zeros(sm.n_node)
    ycnt = np.zeros(sm.n_node)
    np.add.at(ysum, sm.triangles.reshape(-1), np.repeat(cent_y, 3))
    np.add.at(ycnt, sm.triangles.reshape(-1), 1.0)
    side = np.where(ysum / np.maximum(ycnt, 1.0) >= 0.0, 1, -1)
    tip = case.outflow_pin_surf
    # inflow pin set = candidates | tip pin (run_loose_coupling closed-body
    # branch); flag it so statistics can drop boundary data
    pinned = case.inflow_candidates | tip

    U = res.U
    ds = res.outs[:, C.OUT_DS1]
    h1 = res.outs[:, C.OUT_H1]
    cf1 = res.outs[:, C.OUT_CF1]
    dsc = ds / chord_n

    _write_csv(
        f"gv5_0_surface_{level}.csv",
        "x,y,z,xc_local,z_over_b,side,ds,ds_over_c,H,B,Psi,Ctau1,Ctau2,"
        "cf1,pinned,tip_masked",
        [(f"{x[i]:.6f}", f"{y[i]:.6f}", f"{z[i]:.6f}", f"{xc_n[i]:.4f}",
          f"{zb[i]:.4f}", int(side[i]), f"{ds[i]:.6e}", f"{dsc[i]:.6e}",
          f"{h1[i]:.4f}", f"{U[i, 2]:.6e}", f"{U[i, 3]:.6e}",
          f"{U[i, 4]:.6e}", f"{U[i, 5]:.6e}", f"{cf1[i]:.6e}",
          int(pinned[i]), int(tip[i])) for i in range(sm.n_node)])

    head = ("k,ds_max,ds_change_rel,ds_neg_floored,mdot_max,ibl_n_iter,"
            "ibl_converged,ibl_final_residual,inflow_n_pinned,cl_p,cl_kj,"
            "cd_p")
    _write_csv(f"gv5_0_history_{level}.csv", head,
               [tuple(h.get(k, "") for k in head.split(","))
                for h in res.history])

    # -- (1) delta*(z) spanwise at fixed x/c stations -------------------------
    edges = np.arange(0.0, 1.0 + 0.5 * Z_BIN, Z_BIN)
    zcen = 0.5 * (edges[:-1] + edges[1:])
    iz = np.clip(np.digitize(zb, edges) - 1, 0, len(zcen) - 1)
    rows = []
    for s in XC_STATIONS:
        for side_val, side_name in ((1, "upper"), (-1, "lower")):
            sel = (np.abs(xc_n - s) <= XC_WIN) & (side == side_val)
            for b in range(len(zcen)):
                m = sel & (iz == b)
                if m.sum() < 1:
                    continue
                rows.append((f"{s:.2f}", side_name, f"{zcen[b]:.3f}",
                             int(m.sum()), f"{float(np.mean(dsc[m])):.6e}",
                             f"{float(np.std(dsc[m])):.6e}",
                             int(np.any(tip[m]))))
    _write_csv(f"gv5_0_dstar_spanwise_{level}.csv",
               "xc_station,side,z_over_b,n_nodes,ds_over_c_mean,"
               "ds_over_c_std,tip_masked", rows)

    # -- (2) crossflow --------------------------------------------------------
    live = ~tip
    a_mag = float(np.max(np.abs(U[live, 1]))) if np.any(live) else 0.0
    b_mag = float(np.max(np.abs(U[live, 2]))) if np.any(live) else 0.0
    i_b = int(np.argmax(np.abs(U[live, 2]))) if np.any(live) else 0
    c1_mag = float(np.max(np.abs(U[live, 4]))) if np.any(live) else 0.0
    c2_mag = float(np.max(np.abs(U[live, 5]))) if np.any(live) else 0.0
    psi_mag = float(np.max(np.abs(U[live, 3]))) if np.any(live) else 0.0
    zb_b = float(zb[live][i_b]) if np.any(live) else float("nan")
    zone = ("tip" if zb_b > 0.9 else "root" if zb_b < 0.1 else "mid")
    _record("GV5.0", f"{level} max|B|/max|A| (argmax z/b={zb_b:.3f}, "
            f"{zone})", "recorded",
            f"{b_mag:.4e} / {a_mag:.4e} = {b_mag / max(a_mag, 1e-30):.4f}")
    _record("GV5.0", f"{level} max|Ct2|/max|Ct1|, max|Psi|", "recorded",
            f"{c2_mag / max(c1_mag, 1e-30):.4f}, {psi_mag:.4e}")

    # -- (3) Delta-CL ----------------------------------------------------------
    hist = res.history
    cl0_p, cl1_p = float(hist[0]["cl_p"]), float(hist[-1]["cl_p"])
    cl0_kj, cl1_kj = float(hist[0]["cl_kj"]), float(hist[-1]["cl_kj"])
    dcl_p, dcl_kj = cl1_p - cl0_p, cl1_kj - cl0_kj
    for name, v0, dv in (("cl_p", cl0_p, dcl_p), ("cl_kj", cl0_kj, dcl_kj)):
        rel = abs(dv) / max(abs(v0), 1e-30)
        flag = ("DOWN" if dv < 0 else "UP") + (
            " (input-limited: |dcl| < A4 2.5% floor)"
            if rel < 0.025 else "")
        _record("GV5.0", f"{level} Delta-{name} viscous-inviscid "
                f"({v0:.4f} -> {v0 + dv:.4f})", "recorded, expect DOWN",
                f"{dv:+.4f} ({100 * rel:.2f}%) {flag}")

    # -- (4) loop cost ---------------------------------------------------------
    _record("GV5.0", f"{level} loose loop n_outer / converged (GV3.2 2.5-D "
            f"count: {GV32_25D_OUTER})", "recorded",
            f"{res.n_outer} / {res.converged}")
    floor = max((h.get("ibl_final_residual", 0.0) for h in hist[1:]),
                default=0.0)
    mdot = [float(h.get("mdot_max", 0.0)) for h in hist]
    _record("GV5.0", f"{level} IBL residual floor / mdot_max first->last",
            "recorded",
            f"{floor:.3e} / {mdot[0]:.3e} -> {mdot[-1]:.3e}"
            + (f" (growth x{mdot[-1] / mdot[1]:.2f} over k)"
               if len(mdot) > 2 and mdot[1] > 0 else ""))
    _record("GV5.0", f"{level} re_chord (Re_MAC/MAC) / wall time",
            "recorded", f"{RE_CHORD:.4e} per m / {wall_s:.0f}s")
    d = res.driver_info
    _record("GV5.0", f"{level} final FP solve: converged / n_newton",
            "recorded",
            f"{d.get('converged')} / {d.get('n_newton')}")

    _panels(level, xc_n, zb, side, dsc, U, cf1, tip, res, zcen)
    return {"level": level, "res": res, "wall_s": wall_s}


# ---------------------------------------------------------------------------
# panels
# ---------------------------------------------------------------------------

def _panels(level, xc_n, zb, side, dsc, U, cf1, tip, res, zcen):
    live = ~tip
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    for j, s in enumerate(XC_STATIONS):
        ax = axes[j // 2, j % 2]
        for side_val, side_name, mk in ((1, "upper", "o"),
                                        (-1, "lower", "s")):
            sel = (np.abs(xc_n - s) <= XC_WIN) & (side == side_val) & live
            ax.plot(zb[sel], dsc[sel], mk, ms=2.6, alpha=0.55,
                    label=side_name)
        ax.axvspan(1.0 - TIP_FRAC, 1.0, color="r", alpha=0.07)
        ax.set_xlabel("z / b_semi")
        ax.set_ylabel("delta* / c_local")
        ax.set_title(f"x/c = {s:.2f} +/- {XC_WIN}")
        ax.grid(alpha=0.3)
        if j == 0:
            ax.legend(fontsize=8)
    fig.suptitle(f"GV5.0 delta*(z) spanwise, M6 {level} M0.5 a3.06 "
                 f"(tip band shaded; pinned/masked nodes excluded)")
    fig.tight_layout()
    path = os.path.join(RESULTS, f"gv5_0_dstar_spanwise_{level}.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}", flush=True)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    sc = axes[0].scatter(xc_n[live], zb[live], c=dsc[live], s=2,
                         cmap="viridis")
    axes[0].set_title("delta* / c_local")
    fig.colorbar(sc, ax=axes[0])
    blog = np.log10(np.maximum(np.abs(U[live, 2]), 1e-12))
    sc = axes[1].scatter(xc_n[live], zb[live], c=blog, s=2, cmap="magma")
    axes[1].set_title("log10 |B| (crossflow profile)")
    fig.colorbar(sc, ax=axes[1])
    c1 = np.maximum(np.abs(U[live, 4]), 1e-12)
    sc = axes[2].scatter(xc_n[live], zb[live],
                         c=np.abs(U[live, 5]) / c1, s=2, cmap="magma")
    axes[2].set_title("|C_t2| / |C_t1|")
    fig.colorbar(sc, ax=axes[2])
    for ax in axes:
        ax.set_xlabel("x / c_local")
        ax.set_ylabel("z / b_semi")
    fig.suptitle(f"GV5.0 surface fields, M6 {level} (tip band masked)")
    fig.tight_layout()
    path = os.path.join(RESULTS, f"gv5_0_crossflow_{level}.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}", flush=True)


def _panel_cl(runs):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
    for run, ls in zip(runs, ("-", "--")):
        h = run["res"].history
        ks = [r["k"] for r in h]
        axes[0].plot(ks, [r.get("cl_p", np.nan) for r in h], "o-", ls=ls,
                     label=f"{run['level']}")
        axes[1].plot(ks, [r.get("cl_kj", np.nan) for r in h], "o-", ls=ls,
                     label=f"{run['level']}")
    axes[0].set_ylabel("cl_p (pressure integral)")
    axes[1].set_ylabel("cl_KJ (Gamma integration)")
    for ax in axes:
        ax.set_xlabel("outer iteration k (0 = inviscid baseline)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("GV5.0 cl history of the loose loop (M6 M0.5 a3.06)")
    fig.tight_layout()
    path = os.path.join(RESULTS, "gv5_0_cl_history.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}", flush=True)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=None,
                    choices=["coarse", "medium"])
    args = ap.parse_args()
    levels = args.levels or [
        lv for lv in ("coarse", "medium")
        if os.path.exists(os.path.join(M6_DIR, f"{lv}.msh"))]
    print(f"GV5.0 levels: {levels}", flush=True)
    runs = [run_level(lv) for lv in levels]
    if runs:
        _panel_cl(runs)
    _write_csv("summary.csv", "gate,metric,band,measured,verdict", SUMMARY)
    n_fail = sum(1 for *_, v in SUMMARY if v == "FAIL")
    n_rec = sum(1 for *_, v in SUMMARY if v == "RECORDED")
    print(f"\nGV5.0: {n_rec} RECORDED / {n_fail} FAIL", flush=True)
    if n_fail:
        print("HONEST FAIL -- see summary.csv", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
