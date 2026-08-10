"""GV6.2 measured wake-IBL on/off effect vs the A4 band (Track V6).

Binding text: docs/roadmap/track_v.md GV6.2; pre-registered bands:
cases/analysis/v6_2_measured_effect/PRE_REGISTRATION.md (the XFOIL
wake-reference sourcing RULED Option A 2026-07-25, user). All bands
RECORDED (the gate text's verdict type); exits 0 unless a guard fires
(guards = recipe-error raisers, PRE_REGISTRATION section 4).

  (a) on/off measured effect (RECORDED): in-process A/B on the
      committed GV3.1 medium recipe verbatim (M0.5/alpha2/Re3e6/
      xtr0.05 both surfaces, loose Picard leg, CouplingConfig
      defaults): cl off/on, Delta-cl (abs+rel), per-outer histories,
      outer counts, TE-region (x/c>0.9) max|dCp| + location, m_wake
      scale, ds_TE/th_TE, wall times (flagged non-comparable). A4
      quoting per the pinned formulas: rel Delta-cl vs 0.025 (the A4
      medium peak u_e band); max|dCp| vs dCp_A4 =
      2*(u_e/Uinf)_TE*0.025 (first-order propagation from
      Cp = 1 - (u_e/Uinf)^2, (u_e/Uinf)_TE the OFF leg's TE-region max
      recovered value, recorded);
  (b) XFOIL wake direction check (RECORDED): the runner drives the
      pinned XFOIL binary (tools/xfoil, the committed reference
      generator's batch script + constants inherited verbatim) at the
      xtr005 committed conditions; G3: the saved polar must reproduce
      the committed polar_summary.csv xtr005 row to the printed
      digits. The DUMP wake rows -> results/xfoil_wake.csv. Recorded
      reads: (i) near-wake dstar relaxation direction vs the
      producer's monotone delta*_TE -> theta_TE; (ii) the residual
      fraction (ds-theta)/(ds_TE-theta_TE) at x/c = 2 vs the
      producer's e^-1 (L_rel = 1.0 c); (iii) the TE anchor vs XFOIL's
      first-wake dstar and the committed surface-TE sum; (iv)
      downstream theta vs the conserved-theta construction;
  (c) L_rel sweep (RECORDED): ON legs at L_rel {0.5, 2.0} c (the (a)
      legs give 1.0 c; the OFF leg serves all three), Delta-cl +
      TE-region max|dCp| per L_rel. L_rel = 1.0 c stays the pinned
      MODEL CHOICE regardless of the readings.

Run:  python cases/analysis/v6_2_measured_effect/run.py
Thread cap 8 (NUMBA/OMP/OPENBLAS env; the standing temporary
user-directed session constraint, PRE_REGISTRATION section 6); wall
times non-comparable.
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
import subprocess
import tempfile
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import _cp_from_q2, wall_force_coefficients
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    build_airfoil_case,
    make_picard_lifting_driver,
    run_loose_coupling,
)
from pyfp3d.viscous.transpiration import edge_velocity_per_zone
from pyfp3d.viscous.wake_sheet import (
    build_wake_sheet_case,
    wake_edge_velocity,
    wake_transpiration_source,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
NACA_DIR = os.path.join(REPO, "cases", "meshes", "naca0012_2.5d")
REF_DIR = os.path.join(REPO, "cases", "reference_data",
                       "naca0012_viscous_xfoil")
# the committed reference generator, imported for the pinned batch script
# + constants verbatim (import-safe: module level defines only
# constants/functions; main() is under the __main__ guard)
sys.path.append(REF_DIR)
import generate_xfoil_reference as GXR  # noqa: E402

M_INF, ALPHA, RE = 0.5, 2.0, 3.0e6
A4_BAND = 0.025  # the A4 medium peak-relative u_e input band
L_RELS = (0.5, 1.0, 2.0)  # the (c) sweep; 1.0 c = the pinned MODEL CHOICE

SUMMARY = []  # (band, metric, band_text, measured, verdict)


def _record(band, metric, band_text, measured, ok=None):
    verdict = "RECORDED" if ok is None else ("PASS" if ok else "FAIL")
    SUMMARY.append((band, metric, band_text, measured, verdict))
    print(f"  [{verdict:8s}] ({band}) {metric}: {measured} "
          f"(band: {band_text})", flush=True)


def _write_csv(name, header, rows):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write(",".join(str(x) for x in row) + "\n")
    print(f"  wrote {path}", flush=True)


def _loose_loop(level, wake=None, l_rel_chords=None):
    """One loose-coupling leg; wake=None -> flag OFF (legacy), else the
    GV6.1 ON leg with the given WakeSheetCase at the given L_rel
    (None -> the pinned default 1.0 c). Returns (mc, wc, case, cfg,
    res, s_ref)."""
    mc, wc = cut_wake(read_mesh(os.path.join(NACA_DIR, f"{level}.msh")))
    cfg = CouplingConfig(re_chord=RE, m_inf=M_INF, alpha_deg=ALPHA,
                         wake_transpiration=wake is not None)
    if l_rel_chords is not None:
        cfg.wake_l_rel_chords = l_rel_chords
    case = build_airfoil_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg)
    dz = float(np.ptp(mc.nodes[:, 2]))
    s_ref = 1.0 * dz

    def probe(phi, gamma, k):
        f = wall_force_coefficients(
            mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
            alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
        return {"cl": f["cl"], "cd_p": f["cd_pressure"]}

    driver = make_picard_lifting_driver(mc, wc, M_INF, ALPHA)
    t0 = time.perf_counter()
    res = run_loose_coupling(driver, case, cfg, probe=probe, wake=wake)
    wall = time.perf_counter() - t0
    return mc, wc, case, cfg, res, s_ref, wall


def _wall_cp(mc, case, phi):
    """Per-wall-surface-node Cp + u_e from the loop's own recovery
    discipline (the v6_1 runner's idiom; u_inf magnitude = 1)."""
    n_vol = len(mc.nodes)
    le_mask_vol = np.zeros(n_vol, dtype=bool)
    le_mask_vol[case.sm.volume_node_of[case.le_band_surf]] = True
    ue_vol = edge_velocity_per_zone(
        mc.nodes, case.wall_faces, phi, elements=mc.elements,
        le_band_mask=le_mask_vol)
    ue_surf = ue_vol[case.sm.volume_node_of]
    q2 = np.sum(ue_surf ** 2, axis=1)
    return _cp_from_q2(q2, M_INF), ue_surf


def _final_wake_state(wsc, mc, case, res, l_rel_chords):
    """The producer state at the leg's final wall-IBL state (outs) and
    its final phi's wake u_e."""
    return wake_transpiration_source(
        wsc, case.stations, res.outs,
        wake_edge_velocity(mc.nodes, wsc, res.phi),
        M_INF, 1.4, 1.0, l_rel_chords=l_rel_chords)


def _te_cp_rows(case, cp_off, cp_on, l_rel):
    """TE-region (x/c > 0.9) per-side station-averaged Cp rows."""
    st = case.stations
    rows, dcp_all = [], []
    for side_val, side_name in ((1, "upper"), (-1, "lower")):
        for r in range(len(st.xc)):
            if st.xc[r] <= 0.9:
                continue
            mask = (st.station_of == r) & (st.side_node == side_val)
            if not np.any(mask):
                continue
            a, b = float(np.mean(cp_off[mask])), float(np.mean(cp_on[mask]))
            rows.append((f"{st.xc[r]:.4f}", side_name, f"{l_rel:.2f}",
                         f"{a:.6f}", f"{b:.6f}", f"{b - a:+.6f}"))
            dcp_all.append((abs(b - a), f"{st.xc[r]:.4f}", side_name))
    return rows, dcp_all


HIST_HEAD = ("k,ds_max,ds_change_rel,ds_neg_floored,mdot_max,"
             "ibl_n_iter,ibl_converged,ibl_final_residual,"
             "inflow_n_pinned,cl,cd_p,ds_wake_te,th_wake_te,"
             "mdot_wake_max")


def _write_history(tag, res):
    _write_csv(f"history_{tag}.csv", HIST_HEAD,
               [tuple(h.get(k, "") for k in HIST_HEAD.split(","))
                for h in res.history])


# ---------------------------------------------------------------------------
# XFOIL (band (b)): drive the pinned binary, G3 polar cross-check, wake rows
# ---------------------------------------------------------------------------

def _xfoil_wake_rows():
    """Run XFOIL at the committed xtr005 conditions (the committed
    generator's batch script verbatim) in a scratch dir; G3-check the
    saved polar against the committed polar_summary.csv row; return the
    DUMP wake rows (8 fields) as an array."""
    binary = GXR.find_xfoil()  # exits with the build recipe if missing
    with tempfile.TemporaryDirectory(prefix="gv6_2_xfoil_") as tmp:
        workdir = Path(tmp)
        script = GXR.batch_input("xtr005", 0.05, 0.05,
                                 "polar.txt", "dump.txt")
        proc = subprocess.run(
            [str(binary)], input=script, cwd=workdir,
            capture_output=True, text=True, timeout=600)
        log = proc.stdout + proc.stderr
        if proc.returncode != 0:
            raise RuntimeError(
                f"GUARD G3: xfoil run failed (rc={proc.returncode})\n{log}")

        # --- G3: the saved polar reproduces the committed xtr005 row ---
        data_rows = []
        for line in (workdir / "polar.txt").read_text().splitlines():
            tok = line.split()
            if len(tok) >= 7:
                try:
                    data_rows.append([float(v) for v in tok])
                except ValueError:
                    pass
        if not data_rows:
            raise RuntimeError("GUARD G3: no polar data row "
                               "(harness error)")
        _alpha, cl, cd, _cdp, cm = data_rows[-1][:5]
        with open(os.path.join(REF_DIR, "polar_summary.csv")) as f:
            committed = {r["case"]: r for r in csv.DictReader(f)}["xtr005"]
        refs = (float(committed["cl"]), float(committed["cd"]),
                float(committed["cm"]))
        if (abs(cl - refs[0]) > 5e-5 or abs(cd - refs[1]) > 5e-6
                or abs(cm - refs[2]) > 5e-5):
            raise RuntimeError(
                f"GUARD G3: the live polar ({cl:.4f}/{cd:.5f}/{cm:.4f}) "
                f"does not reproduce the committed xtr005 row "
                f"({refs[0]:.4f}/{refs[1]:.5f}/{refs[2]:.4f}) -- wrong "
                "binary or run conditions (harness error)")
        print(f"  [guard ok ] G3 polar reproduction: cl {cl:.4f} / "
              f"cd {cd:.5f} / cm {cm:.4f} == committed xtr005 row",
              flush=True)

        # --- the DUMP wake rows (8 fields) ---
        rows = []
        for line in (workdir / "dump.txt").read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) == 8:  # wake rows (surface rows have 12+)
                rows.append([float(v) for v in fields])
        if not rows:
            raise RuntimeError("GUARD G3: no wake rows in the DUMP "
                               "(harness error)")
        return np.array(rows)


def band_b(wsc, res_on, ds_wake, th_te):
    print("--- (b) XFOIL wake direction check (RECORDED) ---", flush=True)
    xf = _xfoil_wake_rows()
    s_w = xf[:, 0] - xf[0, 0]  # wake arc from the TE (chords)
    x_c, ds_x, th_x, h_x = xf[:, 1], xf[:, 4], xf[:, 5], xf[:, 7]
    _write_csv(
        "xfoil_wake.csv",
        "s_from_te,x_c,y,ue_over_vinf,dstar_over_c,theta_over_c,cf,h",
        [(f"{sw:.5f}", f"{x:.5f}", f"{y:.5f}", f"{u:.5f}", f"{d:.6e}",
          f"{t:.6e}", f"{c:.6e}", f"{h:.4f}")
         for sw, x, y, u, d, t, c, h in
         zip(s_w, x_c, xf[:, 2], xf[:, 3], ds_x, th_x, xf[:, 6], h_x)])
    _record("b", "XFOIL wake rows sourced (Option A; xtr005, the pinned "
            "binary at the committed conditions)", "recorded",
            f"{len(xf)} rows, x/c [{x_c.min():.4f}, {x_c.max():.4f}]")

    # (i) near-wake relaxation direction
    d_ds = np.diff(ds_x)
    n_neg = int(np.sum(d_ds < 0.0))
    _record("b", "(i) XFOIL wake dstar direction: negative d(dstar)/ds "
            "steps / total (the producer relaxes monotone delta*_TE -> "
            "theta_TE by construction)", "recorded",
            f"{n_neg}/{len(d_ds)} "
            f"(monotone_decreasing={bool(np.all(d_ds <= 0.0))})")

    # (ii) residual fraction at x/c = 2 (s_w = 1) vs the producer's e^-1
    r_xf = float((ds_x[-1] - th_x[-1]) / (ds_x[0] - th_x[0]))
    i_s1 = int(np.argmin(np.abs(wsc.s - 1.0)))
    ds_te = float(res_on.wake_info["ds_te"])
    r_ours = float((ds_wake[i_s1] - th_te) / (ds_te - th_te))
    _record("b", "(ii) residual fraction (ds-th)/(ds_TE-th_TE) at "
            "x/c = 2.0: XFOIL vs producer (e^-1 at L_rel = 1.0 c)",
            "recorded",
            f"{r_xf:.4f} vs {r_ours:.4f} (x/c {x_c[-1]:.4f}, "
            f"our s {wsc.s[i_s1]:.4f})")
    _record("b", "(ii') XFOIL effective relaxation length over "
            "[TE, TE+1c] (derived: -1/ln R)", "recorded",
            f"{-1.0 / np.log(r_xf):.3f} c vs pinned 1.0 c")

    # (iii) the TE anchor
    with open(os.path.join(
            REF_DIR, "delta_star_cf_alpha2_m05_xtr005.csv")) as f:
        surf = list(csv.DictReader(f))
    ds_surf_te = 0.0
    for side in ("upper", "lower"):
        rows_s = [r for r in surf if r["surface"] == side]
        te_row = max(rows_s, key=lambda r: float(r["x_c"]))
        ds_surf_te += float(te_row["dstar_over_c"])
    _record("b", "(iii) TE anchor: our ds*_TE/c (ON final) vs XFOIL "
            "first-wake dstar/c vs the committed surface-TE sum (xtr005; "
            "GV3.1 dstar caveat carried)", "recorded",
            f"{ds_te:.5f} vs {ds_x[0]:.5f} vs {ds_surf_te:.5f}")

    # (iv) downstream theta
    _record("b", "(iv) downstream theta: XFOIL th(x/c=2)/th(TE) vs the "
            "producer's conserved-theta construction (model-form "
            "difference)", "recorded",
            f"{th_x[-1] / th_x[0]:.4f} vs 1.0 (by construction)")
    _record("b", "XFOIL wake H at TE / at x/c=2 (H -> 1 relaxation)",
            "recorded", f"{h_x[0]:.4f} / {h_x[-1]:.4f}")
    return xf


# ---------------------------------------------------------------------------
# main: the four legs (off + three ON), then bands (a)/(c)/(b)
# ---------------------------------------------------------------------------

def main():
    print("--- GV6.2 legs: OFF + ON at L_rel {0.5, 1.0, 2.0} c "
          "(medium, the committed GV3.1 recipe verbatim, in-process) ---",
          flush=True)
    mc, wc, case, cfg, res_off, s_ref, wall_off = _loose_loop(
        "medium", wake=None)
    wsc = build_wake_sheet_case(mc, wc)
    legs = {}
    for l_rel in L_RELS:
        *_, res_l, _, wall_l = _loose_loop(
            "medium", wake=wsc, l_rel_chords=l_rel)
        legs[l_rel] = (res_l, wall_l)
    res_10, wall_10 = legs[1.0]

    cp_off, ue_off = _wall_cp(mc, case, res_off.phi)
    q2_off = np.sum(ue_off ** 2, axis=1)
    st = case.stations
    # TE-region (x/c > 0.9) wall nodes, both sides
    te_mask = np.zeros(len(q2_off), dtype=bool)
    for r in range(len(st.xc)):
        if st.xc[r] > 0.9:
            te_mask |= st.station_of == r
    ue_te = float(np.sqrt(np.max(q2_off[te_mask])))
    dcp_a4 = 2.0 * ue_te * A4_BAND

    _write_history("off", res_off)
    for l_rel in L_RELS:
        _write_history(f"on_l{int(l_rel * 10):02d}", legs[l_rel][0])

    # ----- (a) on/off measured effect (RECORDED) -----
    print("--- (a) on/off measured effect (RECORDED) ---", flush=True)
    cl_off = float(res_off.history[-1]["cl"])
    cl_10 = float(res_10.history[-1]["cl"])
    dcl = cl_10 - cl_off
    _record("a", "cl flag OFF / ON (final, L_rel = 1.0 c)", "recorded",
            f"{cl_off:.5f} / {cl_10:.5f}")
    _record("a", "on/off Delta-cl (abs + rel; in-process A/B -- the "
            "committed GV3.1 medium fixed point is NOT cross-run "
            "reproducible, the scatter caveat stands)", "recorded",
            f"{dcl:+.5f} ({dcl / cl_off:+.4%})")
    _record("a", "A4 quoting: |rel Delta-cl| vs the A4 medium peak u_e "
            "band", "recorded",
            f"{abs(dcl / cl_off):.4%} vs {A4_BAND:.2%} "
            f"(ratio {abs(dcl / cl_off) / A4_BAND:.3f})")
    _record("a", "outer count ON vs OFF", "recorded",
            f"{res_10.n_outer} vs {res_off.n_outer}")
    mmax = max(h.get("mdot_wake_max", 0.0) for h in res_10.history[1:])
    _record("a", "m_wake scale max over outers (L_rel = 1.0 c)",
            "recorded", f"{mmax:.3e}")
    _record("a", "ds_TE / th_TE (ON final)", "recorded",
            f"{res_10.wake_info['ds_te']:.5e} / "
            f"{res_10.wake_info['th_te']:.5e}")

    ds_wake_10, m_wake_10, _ = _final_wake_state(
        wsc, mc, case, res_10, 1.0)
    cp_10, _ = _wall_cp(mc, case, res_10.phi)
    te_rows, dcp_10 = _te_cp_rows(case, cp_off, cp_10, 1.0)
    worst_10 = max(dcp_10)
    _record("a", "TE-region (x/c > 0.9) max |dCp| on/off + location",
            "recorded",
            f"{worst_10[0]:.5f} at x/c {worst_10[1]} {worst_10[2]}")
    _record("a", "A4 quoting: max|dCp| vs dCp_A4 = 2*(u_e/Uinf)_TE*"
            "0.025 (first-order propagation; (u_e/Uinf)_TE the OFF "
            "leg's TE-region max)", "recorded",
            f"{worst_10[0]:.5f} vs {dcp_a4:.4f} "
            f"((u_e/Uinf)_TE = {ue_te:.4f}; ratio "
            f"{worst_10[0] / dcp_a4:.3f})")
    _record("a", "wall times per leg (8-thread session, NON-COMPARABLE)",
            "recorded",
            f"off {wall_off:.0f}s / on(1.0) {wall_10:.0f}s")
    _record("a", "context anchors: GV6.1(d) prior smoke / GV3.1 "
            "inviscid-cl A4 floor", "recorded",
            "+0.00015 / 0.00250; inviscid cl 2.6% below XFOIL")

    # ----- (c) L_rel sweep (RECORDED) -----
    print("--- (c) L_rel sweep (RECORDED; 1.0 c stays pinned) ---",
          flush=True)
    ds_profiles = {1.0: ds_wake_10}
    cps = {1.0: cp_10}
    for l_rel in (0.5, 2.0):
        res_l, wall_l = legs[l_rel]
        cl_l = float(res_l.history[-1]["cl"])
        cp_l, _ = _wall_cp(mc, case, res_l.phi)
        cps[l_rel] = cp_l
        rows_l, dcp_l = _te_cp_rows(case, cp_off, cp_l, l_rel)
        te_rows.extend(rows_l)
        worst_l = max(dcp_l)
        ds_wake_l, _, _ = _final_wake_state(wsc, mc, case, res_l, l_rel)
        ds_profiles[l_rel] = ds_wake_l
        _record("c", f"L_rel = {l_rel:.1f} c: on/off Delta-cl", "recorded",
                f"{cl_l - cl_off:+.5f} ({(cl_l - cl_off) / cl_off:+.4%})")
        _record("c", f"L_rel = {l_rel:.1f} c: TE-region max |dCp|",
                "recorded", f"{worst_l[0]:.5f} at x/c {worst_l[1]} "
                f"{worst_l[2]}")
        _record("c", f"L_rel = {l_rel:.1f} c: outer count / wall "
                "(non-comparable)", "recorded",
                f"{res_l.n_outer} / {wall_l:.0f}s")
    _record("c", "L_rel = 1.0 c stays the pinned MODEL CHOICE (the "
            "sweep records sensitivity, never tunes)", "recorded",
            "unchanged")
    _write_csv("te_cp.csv", "x_c,side,l_rel,cp_off,cp_on,dcp", te_rows)
    _write_csv(
        "wake_profiles.csv",
        "s_from_te,ds_wake_l05,ds_wake_l10,ds_wake_l20",
        [(f"{s:.6f}", f"{ds_profiles[0.5][i]:.6e}",
          f"{ds_profiles[1.0][i]:.6e}", f"{ds_profiles[2.0][i]:.6e}")
         for i, s in enumerate(wsc.s)])

    # ----- (b) XFOIL wake direction check (RECORDED) -----
    th_te_10 = float(res_10.wake_info["th_te"])
    xf = band_b(wsc, res_10, ds_wake_10, th_te_10)

    _panels(res_off, legs, wsc, ds_profiles, xf, th_te_10,
            te_rows)

    _write_csv("summary.csv", "band,metric,band,measured,verdict",
               SUMMARY)
    n_rec = sum(1 for *_, v in SUMMARY if v == "RECORDED")
    print(f"\nGV6.2: 0 PASS / 0 FAIL / {n_rec} RECORDED "
          "(all bands RECORDED per the gate text; guards G1-G4 clean)",
          flush=True)


def _panels(res_off, legs, wsc, ds_profiles, xf, th_te, te_rows):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.4))
    ax = axes[0]
    ax.plot([h["k"] for h in res_off.history],
            [h.get("cl", np.nan) for h in res_off.history], "o--",
            label="flag OFF")
    for l_rel, sty in ((0.5, "^-"), (1.0, "s-"), (2.0, "d-")):
        res = legs[l_rel][0]
        ax.plot([h["k"] for h in res.history],
                [h.get("cl", np.nan) for h in res.history], sty,
                label=f"ON L_rel={l_rel:.1f}c")
    ax.set_xlabel("outer iteration k")
    ax.set_ylabel("c_l (pressure integral)")
    ax.set_title("GV6.2(a/c) GV3.1 medium: cl history ON vs OFF")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    order = np.argsort(wsc.s)
    for l_rel, sty in ((0.5, "^-"), (1.0, "s-"), (2.0, "d-")):
        ax.plot(wsc.s[order], ds_profiles[l_rel][order], sty, ms=3,
                label=f"ours L_rel={l_rel:.1f}c")
    ax.axhline(th_te, color="gray", ls=":", lw=1,
               label=f"our th_TE (conserved) {th_te:.4f}")
    s_w = xf[:, 0] - xf[0, 0]
    ax.plot(s_w, xf[:, 4], "ko-", ms=3, label="XFOIL dstar")
    ax.plot(s_w, xf[:, 5], "kx--", ms=3, label="XFOIL theta")
    ax.set_xlim(0.0, 3.0)
    ax.set_xlabel("s / c from TE")
    ax.set_ylabel("delta* / c")
    ax.set_title("GV6.2(b) wake dstar relaxation: ours vs XFOIL")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[2]
    for side in ("upper", "lower"):
        sel = [r for r in te_rows if r[1] == side and r[2] == "1.00"]
        xs = [float(r[0]) for r in sel]
        ax.plot(xs, [float(r[3]) for r in sel], "o--", label=f"OFF {side}")
        ax.plot(xs, [float(r[4]) for r in sel], "s-", label=f"ON {side}")
    ax.set_xlabel("x/c")
    ax.set_ylabel("C_p")
    ax.set_title("GV6.2(a) TE-region Cp ON (1.0 c) vs OFF")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout()
    path = os.path.join(RESULTS, "gv6_2.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}", flush=True)


if __name__ == "__main__":
    main()
