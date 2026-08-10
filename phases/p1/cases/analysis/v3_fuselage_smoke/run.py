"""GV3.3 fuselage body-of-revolution smoke (Track V3).

Binding text: docs/roadmap/track_v.md GV3.3 (2026-07-22, user-directed);
pre-registered bands: cases/analysis/v3_fuselage_smoke/PRE_REGISTRATION.md
(committed before the first execution). Regenerates every CSV/PNG in
results/ and exits 0 iff all pre-registered assertions hold.

  (a) azimuthal delta* collapse at fixed x-stations (sigma/mu <= 0.15 at
      every window station) -- the surface-FE scatter measure on an
      unstructured triangulation, the genuinely-3-D GV1.1(d) analogue;
  (b) crossflow unknowns ~ 0 (max|B| <= 0.05 max|A|, max|C_t2| <=
      0.05 max|C_t1|) -- axisymmetric flow on a 3-D surface;
  (c)/(d) RECORDED: stagnation-band seeding; transpiration on/off Cp
      delta + tail-cone H rise (indicated tail separation recorded and
      masked, not chased).

Run:  python cases/analysis/v3_fuselage_smoke/run.py [--levels coarse medium]
"""

import argparse
import os
import sys

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

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    build_closed_body_case,
    make_picard_nonlifting_driver,
    run_loose_coupling,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BOR_DIR = os.path.join(REPO, "cases", "meshes", "fuselage_bor")

M_INF, ALPHA = 0.3, 0.0
RE_BODY = 3.0e6           # per body length (pre-registered)
X_TR, STAG_BAND = 0.05, 0.05
CV_BAND = 0.15            # (a) azimuthal sigma/mu of delta*
CROSS_BAND = 0.05         # (b) crossflow/streamwise unknown maxima
WIN_LO, WIN_HI = 0.20, 0.95   # (a) window, fraction of body length

SUMMARY = []


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


def run_level(level):
    print(f"--- GV3.3 {level}: fuselage BoR M={M_INF} alpha={ALPHA} ---",
          flush=True)
    mesh = read_mesh(os.path.join(BOR_DIR, f"{level}.msh"))
    nodes, elements = mesh.nodes, mesh.elements
    wall = mesh.boundary_faces["wall"]
    ff = np.unique(mesh.boundary_faces["farfield"])
    phi_ff = nodes[ff, 0].copy()          # u_inf along +x, alpha = 0

    wall_x = nodes[np.unique(wall), 0]
    x0, x1 = float(wall_x.min()), float(wall_x.max())
    body_len = x1 - x0
    re_chord = RE_BODY / body_len
    cfg = CouplingConfig(re_chord=re_chord, m_inf=M_INF, alpha_deg=ALPHA)
    case = build_closed_body_case(nodes, elements, wall, cfg,
                                  x_tr_frac=X_TR, stag_band_frac=STAG_BAND)
    driver = make_picard_nonlifting_driver(nodes, elements, ff, phi_ff,
                                           m_inf=M_INF)

    cp_snap = {}

    def probe(phi, gamma, k):
        f = wall_force_coefficients(nodes, elements, wall, phi,
                                    alpha_deg=ALPHA, m_inf=M_INF)
        cp_snap[k] = f["cp_tri"]
        return {"cp_mean": float(np.mean(f["cp_tri"])),
                "cd_p": f["cd_pressure"]}

    res = run_loose_coupling(driver, case, cfg, probe=probe)
    print(f"    converged={res.converged} n_outer={res.n_outer} "
          f"omega={cfg.omega}", flush=True)

    sm = case.sm
    x_s = sm.xyz[:, 0]
    ds = res.outs[:, C.OUT_DS1]
    h1 = res.outs[:, C.OUT_H1]
    U = res.U

    # -- (a) azimuthal collapse ---------------------------------------------
    tris = np.asarray(wall, dtype=np.int64)
    pts = nodes[tris]
    edge = np.linalg.norm(pts[:, [1, 2, 0]] - pts, axis=2).mean()
    bw = 2.0 * float(edge)
    edges = np.arange(x0, x1 + 0.5 * bw, bw)
    ibin = np.digitize(x_s, edges)
    rows, worst_cv, worst_cv_x = [], 0.0, None
    for b in range(1, len(edges)):
        sel = ibin == b
        if sel.sum() < 6:
            continue
        xc = 0.5 * (edges[b - 1] + edges[b])
        frac = (xc - x0) / body_len
        v = ds[sel]
        mu = float(np.mean(v))
        if mu <= 0.0:
            continue
        sd = float(np.std(v))
        cv = sd / mu
        mm = float((v.max() - v.min()) / mu)
        win = WIN_LO <= frac <= WIN_HI
        rows.append((f"{xc:.5f}", f"{frac:.4f}", int(sel.sum()),
                     f"{mu:.6e}", f"{sd:.6e}", f"{cv:.4f}", f"{mm:.4f}",
                     "window" if win else "outer"))
        if win and cv > worst_cv:
            worst_cv, worst_cv_x = cv, frac
    _write_csv(f"gv3_3_azimuthal_{level}.csv",
               "x_station,x_frac,n_nodes,ds_mean,ds_std,ds_cv,"
               "ds_maxmin_over_mean,window", rows)
    tag = "" if level == "coarse" else f"{level} "
    ok = worst_cv <= CV_BAND
    _record("GV3.3(a)", f"{tag}azimuthal sigma/mu(delta*) worst "
            f"(x/L={worst_cv_x:.3f})" if worst_cv_x else
            f"{tag}azimuthal sigma/mu(delta*)",
            f"<= {CV_BAND} every window station",
            f"{worst_cv:.4f}", ok if level == "coarse" else None)

    # -- (b) crossflow -------------------------------------------------------
    a_mag = float(np.max(np.abs(U[:, 1])))
    b_mag = float(np.max(np.abs(U[:, 2])))
    psi_mag = float(np.max(np.abs(U[:, 3])))
    c1_mag = float(np.max(np.abs(U[:, 4])))
    c2_mag = float(np.max(np.abs(U[:, 5])))
    r_b = b_mag / max(a_mag, 1e-30)
    r_c = c2_mag / max(c1_mag, 1e-30)
    ok_b = bool(r_b <= CROSS_BAND and r_c <= CROSS_BAND)
    _record("GV3.3(b)", f"{tag}max|B|/max|A|, max|Ct2|/max|Ct1|",
            f"<= {CROSS_BAND}", f"{r_b:.4f}, {r_c:.4f}",
            ok_b if level == "coarse" else None)
    _record("GV3.3(b)", f"{tag}max|Psi|, max|DS2|/|DS1|, max|CF2|/|CF1|",
            "recorded",
            f"{psi_mag:.3e}, "
            f"{float(np.max(np.abs(res.outs[:, C.OUT_DS2]))):.3e}/"
            f"{float(np.max(np.abs(ds))):.3e}, "
            f"{float(np.max(np.abs(res.outs[:, C.OUT_CF2]))):.3e}/"
            f"{float(np.max(np.abs(res.outs[:, C.OUT_CF1]))):.3e}", None)

    # -- (c) stagnation-band seeding (RECORDED) ------------------------------
    h1rec = res.history[1] if len(res.history) > 1 else {}
    _record("GV3.3(c)", f"{tag}inflow pinned nodes / LE-band nodes",
            "recorded",
            f"{h1rec.get('inflow_n_pinned', '?')} / "
            f"{int(case.le_band_surf.sum())}", None)

    # -- (d) transpiration on/off + tail H rise (RECORDED) -------------------
    ks = sorted(cp_snap)
    dcp = np.abs(cp_snap[ks[-1]] - cp_snap[ks[0]])
    _record("GV3.3(d)", f"{tag}wall |dCp| on/off: max / mean",
            "recorded", f"{float(dcp.max()):.4f} / "
            f"{float(np.mean(dcp)):.4f}", None)
    n_neg = int(np.count_nonzero(ds < 0.0))
    _record("GV3.3(d)", f"{tag}nodes with delta* < 0 (separation "
            "indicator)", "recorded", n_neg, None)
    h_rows = []
    for r in rows:
        sel = ibin == np.digitize(float(r[0]), edges)
        h_rows.append((float(r[0]), float(np.mean(h1[sel]))))
    if h_rows:
        hx = np.array([r[0] for r in h_rows])
        hv = np.array([r[1] for r in h_rows])
        i_tail = int(np.argmax(hv))
        _record("GV3.3(d)", f"{tag}H(x) max {hv.max():.2f} at "
                f"x/L={(hx[i_tail]-x0)/body_len:.3f} (cyl "
                f"{float(np.mean(hv[(hx > x0 + 0.3 * body_len) & (hx < x0 + 0.6 * body_len)])):.2f})",
                "recorded", "tail-cone APG rise" if hx[i_tail] > x0 + 0.6 * body_len else "no tail rise", None)

    # -- history + loop record -----------------------------------------------
    head = "k,ds_max,ds_change_rel,mdot_max,ibl_n_iter,ibl_final_residual,cp_mean,cd_p"
    _write_csv(f"gv3_3_history_{level}.csv", head,
               [tuple(h.get(k, "") for k in head.split(","))
                for h in res.history])
    _record("GV3.3", f"{tag}loose loop n_outer / omega / converged",
            "recorded",
            f"{res.n_outer} / {cfg.omega} / {res.converged}", None)
    _record("GV3.3", f"{tag}Re per body length (re_chord)",
            "recorded", f"{RE_BODY:.3e} ({re_chord:.4e})", None)

    tri_xf = (nodes[tris].mean(axis=1)[:, 0] - x0) / body_len
    _panel(level, x_s, ds, h1, U, x0, body_len, rows, cp_snap, ks,
           tri_xf)
    return res


def _panel(level, x_s, ds, h1, U, x0, body_len, rows, cp_snap, ks,
           tri_xf):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    xf = (x_s - x0) / body_len

    ax = axes[0, 0]
    ax.plot(xf, ds, ".", ms=1.0, alpha=0.3, label="wall nodes")
    wx = [float(r[0]) for r in rows]
    wm = [float(r[3]) for r in rows]
    ws = [float(r[4]) for r in rows]
    ax.errorbar([(v - x0) / body_len for v in wx], wm, yerr=ws, fmt="k-",
                lw=1.4, capsize=2, label="station mean +/- sigma")
    ax.axvspan(WIN_LO, WIN_HI, color="g", alpha=0.06)
    ax.set_xlabel("(x - x_nose)/L")
    ax.set_ylabel("delta*")
    ax.set_title(f"GV3.3(a) azimuthal collapse ({level})")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    ax.plot(xf, h1, ".", ms=1.0, alpha=0.3)
    ax.set_xlabel("(x - x_nose)/L")
    ax.set_ylabel("H = delta*/theta")
    ax.set_title("GV3.3(d) shape factor (tail-cone APG rise)")
    ax.set_ylim(-5.0, 15.0)  # clip stagnation/pin-edge theta->0 artifacts
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.plot(tri_xf, cp_snap[ks[0]], ".", ms=1.5, alpha=0.35,
            label=f"off (k={ks[0]})")
    ax.plot(tri_xf, cp_snap[ks[-1]], ".", ms=1.5, alpha=0.35,
            label=f"on (k={ks[-1]})")
    ax.set_title("GV3.3(d) Cp on/off transpiration")
    ax.set_xlabel("(x - x_nose)/L")
    ax.set_ylabel("Cp")
    ax.legend(fontsize=8, markerscale=3)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    ax.plot(xf, np.abs(U[:, 2]), ".", ms=1.0, alpha=0.3, label="|B|")
    ax.plot(xf, np.abs(U[:, 1]), ".", ms=1.0, alpha=0.3, label="|A|")
    ax.plot(xf, np.abs(U[:, 5]), ".", ms=1.0, alpha=0.3, label="|C_t2|")
    ax.set_yscale("log")
    ax.set_xlabel("(x - x_nose)/L")
    ax.set_title("GV3.3(b) crossflow vs streamwise unknowns")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    path = os.path.join(RESULTS, f"gv3_3_panels_{level}.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=None,
                    choices=["coarse", "medium"])
    args = ap.parse_args()
    levels = args.levels or [
        lv for lv in ("coarse", "medium")
        if os.path.exists(os.path.join(BOR_DIR, f"{lv}.msh"))]
    print(f"GV3.3 levels: {levels}", flush=True)
    for level in levels:
        run_level(level)
    _write_csv("gv3_3_summary.csv", "gate,metric,band,measured,verdict",
               SUMMARY)
    n_fail = sum(1 for *_, v in SUMMARY if v == "FAIL")
    n_pass = sum(1 for *_, v in SUMMARY if v == "PASS")
    n_rec = sum(1 for *_, v in SUMMARY if v == "RECORDED")
    print(f"\nGV3.3: {n_pass} PASS / {n_fail} FAIL / {n_rec} RECORDED",
          flush=True)
    if n_fail:
        print("HONEST FAIL -- see gv3_3_summary.csv", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
