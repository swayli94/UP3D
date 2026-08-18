"""GS4.1 round 14 (G13) -- is the post-transition relaxation one defect or two?

Binding text: docs/dev_phase_four/20260821-0800-g13-prereg.md (committed before
any code). Re-scoped there, with the reason: the plate's two external correlations
cannot adjudicate an attribution question, because round 9 leg A measured that a
MORE faithful closure agreed with them WORSE.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_relaxation/run.py
"""

import csv
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "bench", "studies", "gs41_lag_xfoil"))
import run as R9                                                    # noqa: E402

RHO, MU, BAND = R9.RHO, R9.MU, 0.02
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    from pyfp3d.viscous import closures_2d as C2
    from pyfp3d.viscous import strip2d as S

    print("== guards ==")
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", "pyfp3d/"],
                       cwd=ROOT, capture_output=True)
    print(f"  G-FROZEN-LIB  pyfp3d/ unchanged: "
          f"{'PASS' if r.returncode == 0 else 'FAIL -- and that is a RESULT'}")
    print("  G-PROV        reference = the locally rebuilt XFOIL 6.99, used whole")
    print("  G-SAMEAXIS    airfoil axis = XFOIL's OWN Re_theta = rho*ue*theta_X/mu "
          "(theta from its dump, NOT ours -- Re_theta contains theta, which is "
          "under comparison: pre-registration section 2)")

    with tempfile.TemporaryDirectory() as wd:
        surf, _ = R9.run_xfoil(wd)
    sa, sb, _ = R9.split_at_stagnation(surf)

    rows = []
    for side in (sa, sb):
        i_tr = R9.transition_index(side)
        if not R9.guard_xtr(side, i_tr):
            raise SystemExit("G-XTR failed")
        ue_fn, _ = R9.ue_interp(side)
        st = side.s[i_tr:]
        y0 = (float(side.theta[i_tr]), float(side.H[i_tr]),
              float(np.sqrt(max(side.ctau[i_tr], 1e-12))))
        lag = S.march_correlation(st[1:], y0, float(st[0]), ue_fn, rho=RHO, mu=MU,
                                 n_substep=8000, x_tr=float(st[0]) * 0.999,
                                 lag=True)
        idx = np.arange(i_tr + 1, len(side.s))
        ret_x = RHO * side.ue[idx] * side.theta[idx] / MU     # XFOIL's own
        for j, k in enumerate(idx):
            # CtauEQ at OUR state, for the relaxation measure
            p = C2.packet_turb(float(lag.theta[j]), float(lag.H[j]),
                               float(side.ue[k]), rho=RHO, mu=MU)
            rows.append({"side": side.name, "x": side.x[k],
                         "re_theta_xfoil": ret_x[j],
                         "theta_x": side.theta[k], "theta_o": lag.theta[j],
                         "ctau_x": side.ctau[k], "ctau_o": lag.ctau[j],
                         "ctau_eq_o": p["Ctau_eq"],
                         "delta_o": C2.bl_thickness(float(lag.theta[j]),
                                                    float(lag.H[j])),
                         "s": side.s[k]})

    # ---- P-AXIS -----------------------------------------------------------
    dth = np.array([abs(r["theta_o"] / r["theta_x"] - 1.0) for r in rows])
    axis_ok = float(np.median(dth)) <= 0.02
    _record("P-AXIS", "our theta vs XFOIL's, per station",
            "median <= 2 %, else the Re_theta axis is NOT common",
            f"median {100*np.median(dth):.2f} %, worst {100*dth.max():.2f} %",
            "P-AXIS PASS" if axis_ok else "P-AXIS FAIL -> P-SAME UNDEFINED")

    # ---- P-SEED -----------------------------------------------------------
    for side in ("side_a", "side_b"):
        z = [r for r in rows if r["side"] == side]
        r0 = z[0]
        _record("P-SEED", f"{side}: Ctau at the first station after transition",
                "ratio ours/XFOIL ~ 1",
                f"ours {r0['ctau_o']:.6e} vs XFOIL {r0['ctau_x']:.6e} = "
                f"{r0['ctau_o']/r0['ctau_x']:.4f}",
                "P-SEED PASS" if abs(r0["ctau_o"] / r0["ctau_x"] - 1.0) <= BAND
                else "P-SEED FAIL -- the transition seed is a suspect")

    # ---- P-RATE -----------------------------------------------------------
    # approach length of Ctau -> CtauEQ, normalised by Delta (a sourced length)
    for side in ("side_a", "side_b"):
        z = [r for r in rows if r["side"] == side]
        s = np.array([r["s"] for r in z])
        eq = np.array([r["ctau_eq_o"] for r in z])
        L = {}
        for who, key in (("ours", "ctau_o"), ("XFOIL", "ctau_x")):
            g = np.abs(np.array([r[key] for r in z]) / eq - 1.0)
            m = (g > 1e-4) & (g < 0.9) & (np.arange(g.size) < g.size // 2)
            if m.sum() < 4:
                L[who] = np.nan
                continue
            sl = np.polyfit(s[m], np.log(g[m]), 1)[0]
            L[who] = -1.0 / sl if sl < 0 else np.nan
        d = float(np.median([r["delta_o"] for r in z[:10]]))
        ok = np.isfinite(L["ours"]) and np.isfinite(L["XFOIL"]) and \
            abs(L["ours"] / L["XFOIL"] - 1.0) <= 0.25
        _record("P-RATE", f"{side}: approach length of Ctau -> CtauEQ, in Delta",
                "ratio ours/XFOIL within 25 %",
                f"ours {L['ours']/d if np.isfinite(L['ours']) else float('nan'):.2f} Delta, "
                f"XFOIL {L['XFOIL']/d if np.isfinite(L['XFOIL']) else float('nan'):.2f} Delta "
                f"(Delta = {d:.3e}); ratio "
                f"{L['ours']/L['XFOIL'] if np.isfinite(L['ours']) and np.isfinite(L['XFOIL']) else float('nan'):.3f}",
                "P-RATE PASS" if ok else "P-RATE FAIL -- Delta is a suspect")

    # ---- P-SAME -----------------------------------------------------------
    plate = list(csv.DictReader(open(os.path.join(
        ROOT, "bench", "studies", "gs41_five_terms", "results",
        "window_fiveterms.csv"))))
    out = [float(p["re_theta"]) for p in plate if p["inside"] == "0"]
    air = [r["re_theta_xfoil"] for r in rows
           if abs(r["ctau_o"] / r["ctau_x"] - 1.0) > BAND
           or abs(r["theta_o"] / r["theta_x"] - 1.0) > BAND]
    if not axis_ok:
        _record("P-SAME", "the two Re_theta ranges", "overlap?",
                "axis not common", "P-SAME UNDEFINED")
    elif not air:
        _record("P-SAME", "the two Re_theta ranges", "overlap?",
                f"plate {min(out):.0f}-{max(out):.0f}; airfoil: NO station "
                f"outside the band", "P-SAME UNDEFINED (airfoil set empty)")
    else:
        lo, hi = max(min(out), min(air)), min(max(out), max(air))
        ov = hi > lo
        _record("P-SAME", "plate's out-of-band Re_theta vs the airfoil's",
                "overlap = CONSISTENT WITH one defect (never proof, section 2); "
                "no overlap = two different defects",
                f"plate {min(out):.0f}-{max(out):.0f} ({len(out)} st); "
                f"airfoil {min(air):.0f}-{max(air):.0f} ({len(air)} st); "
                f"overlap {'YES ' + f'{lo:.0f}-{hi:.0f}' if ov else 'NO'}",
                "P-SAME: consistent with one defect" if ov
                else "P-SAME: TWO DIFFERENT DEFECTS")

    _record("P-CORR", "the plate against Coles-Fernholz and the 1/4 power law",
            "RECORDED ONLY -- round 9 leg A measured that a MORE faithful "
            "closure agreed with these WORSE, so their applicability to this "
            "family is not established",
            f"{sum(1 for p in plate if p['inside']=='1')}/{len(plate)} inside",
            "RECORDED")

    with open(os.path.join(RESULTS, "relaxation.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n== summary ==  {time.perf_counter()-t0:.2f} s")
    for tag, metric, band, measured, verdict in SUMMARY:
        print(f"  {verdict:34s} [{tag}] {measured}")
    print("\n★ The plate is ZPG and the airfoil is not, so they are NOT the same "
          "physical problem: numerical overlap counts as CONSISTENT WITH one "
          "defect, never as proof (pre-registration section 2).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
