"""GS4.1 round 8 -- fix the omitted GCC/Re_theta in CtauEQ, and record what moves.

Binding text: docs/dev_phase_four/20260820-0500-gcc-fix-prereg.md (committed
before any code). No new functionality: no lag equation, no XFOIL comparison.

★★ The fix is justified by the SOURCE, not by the outcome. It stands even if
E-CF gets worse; no criterion here is phrased as "did E-CF improve"
(pre-registration section 2).
"""

import csv
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
R6 = os.path.join(ROOT, "bench", "studies", "gs41_repaired_criteria")
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)

RHO, MU, U_INF = 1.0, 1.0e-5, 1.0
X0, X_TR, X1 = 0.05, 5.0, 400.0
RE_HI, ATTRACT_TOL, SEED_PERTURB = 1.0e4, 0.01, 1.15
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


# --- independent transcription of XFOIL's CQ2 block; imports NO closures_2d ---
def xfoil_ctau_eq(hs, H, us, re_theta, hk, ctcon, gccon):
    """xblsys.f:856-877, BLVAR, ITYP=2 (turbulent wall). Written from the source
    for this check; shares no code with the library (round 4's discipline)."""
    hkc = hk - 1.0 - gccon / re_theta
    if hkc < 0.01:
        hkc = 0.01
    hkb = hk - 1.0
    usb = 1.0 - us
    return ctcon * hs * hkb * hkc ** 2 / (usb * H * hk ** 2)


def cf_coles_fernholz(ret):
    return 2.0 * ((1.0 / 0.384) * np.log(ret) + 4.127) ** -2


def cf_power_law(ret):
    return 0.024 * ret ** -0.25


def band_at(ret):
    a, b = cf_coles_fernholz(ret), cf_power_law(ret)
    s = abs(a - b) / (0.5 * (a + b))
    return max(0.03, 2.0 * s), s


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    from pyfp3d.viscous import closures_2d as C2
    from pyfp3d.viscous import strip2d as S

    print("== guards ==")
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--",
                        "pyfp3d/viscous/closures.py", "pyfp3d/viscous/ibl3.py"],
                       cwd=ROOT, capture_output=True)
    print(f"  G-FROZEN  ibl3.py + closures.py unchanged: "
          f"{'PASS' if r.returncode == 0 else 'FAIL'}")
    if r.returncode:
        raise SystemExit("G-FROZEN failed")

    # F-SOURCE
    worst = 0.0
    for H, ret in ((1.30, 600.0), (1.40, 2000.0), (1.35, 5000.0),
                   (1.50, 578.0), (2.20, 800.0), (1.33, 1.0e4)):
        hs = C2.h_star_turb(H, ret)
        us = C2.slip_velocity(hs, H)
        lib = C2.ctau_eq(hs, H, us, ret)
        ref = xfoil_ctau_eq(hs, H, us, ret, C2.h_kinematic(H), C2.CTCON,
                            C2.GCCON)
        worst = max(worst, abs(lib / ref - 1.0))
    _record("F-SOURCE", "ctau_eq vs an independent transcription of CQ2",
            "<= 1e-14 (same formula, two orderings)", f"{worst:.2e}",
            "PASS" if worst <= 1e-14 else "F-FAIL -> kill 1")
    if worst > 1e-14:
        raise SystemExit("F-SOURCE failed")

    # F-USED -- the gap round 5's G-SOURCE left: printed the constant, never
    # checked it was used.
    hs = C2.h_star_turb(1.4, 2000.0)
    us = C2.slip_velocity(hs, 1.4)
    with_g = C2.ctau_eq(hs, 1.4, us, 2000.0)
    saved, C2.GCCON = C2.GCCON, 0.0
    without_g = C2.ctau_eq(hs, 1.4, us, 2000.0)
    C2.GCCON = saved
    rel = abs(with_g / without_g - 1.0)
    _record("F-USED", "GCCON actually reaches ctau_eq's return value",
            "setting GCCON=0 must change the result",
            f"{100*rel:.2f} % change at Re_theta 2000",
            "PASS" if rel > 1e-6 else "F-FAIL")

    # F-LEGACY -- GCC touches only the turbulent branch
    h = C2.zpg_fixed_point()
    y0 = C2.blasius_state(0.01, H=2.591100)
    st = S.march_correlation(np.array([1.0, 100.0]), y0, 0.01,
                             S.flat_plate_ue(U_INF), n_substep=2000)
    ok = abs(h - 2.590433) < 1e-6 and abs(st.H[-1] - 2.590433) < 5e-4
    _record("F-LEGACY", "laminar path unmoved", "ZPG fixed point 2.590433",
            f"{h:.6f}; march H[-1] {st.H[-1]:.6f}",
            "PASS" if ok else "F-FAIL -> kill 2")
    if not ok:
        raise SystemExit("F-LEGACY failed")

    # F-RERUN -- round 6's protocol verbatim, RECORDED (not re-judged)
    print("== F-RERUN: round 6's protocol, post-fix (RECORDED, not re-judged) ==")
    lead = S.march_correlation(np.array([X_TR * 1.02]),
                               C2.blasius_state(X0, ue=U_INF, rho=RHO, mu=MU,
                                                H=2.591100),
                               X0, S.flat_plate_ue(U_INF), rho=RHO, mu=MU,
                               n_substep=8000, x_tr=X_TR)
    x_s, th_s, H_s = lead.x[0], lead.theta[0], lead.H[0]
    stations = np.geomspace(x_s * 1.02, X1, 120)
    kw = dict(rho=RHO, mu=MU, n_substep=8000, x_tr=X_TR)
    A = S.march_correlation(stations, (th_s, H_s), x_s,
                            S.flat_plate_ue(U_INF), **kw)
    B = S.march_correlation(stations, (th_s, H_s * SEED_PERTURB), x_s,
                            S.flat_plate_ue(U_INF), **kw)
    d = np.abs(B.H - A.H) / A.H
    start = next((i for i in np.where(d <= ATTRACT_TOL)[0]
                  if np.all(d[i:] <= ATTRACT_TOL)), None)
    _record("E-ATTRACT", "post-fix", "<= 1 % and staying",
            f"Re_theta {A.re_theta[start]:.0f} (round 6: 578)"
            if start is not None else "never", "RECORDED")

    m = (np.arange(A.x.size) >= start) & (A.re_theta <= RE_HI)
    rows, n_out = [], 0
    for i in np.where(m)[0]:
        ret, cf = A.re_theta[i], A.cf[i]
        cfa, cfb = cf_coles_fernholz(ret), cf_power_law(ret)
        b, s = band_at(ret)
        if s > 0.10:
            continue
        dev = max(abs(cf/cfa - 1.0), abs(cf/cfb - 1.0))
        n_out += dev > b
        rows.append({"re_theta": ret, "H": A.H[i], "cf": cf, "band": b,
                     "dev": dev, "inside": int(dev <= b)})
    _record("E-PHYS", "post-fix", "1.05 < H < 4",
            f"[{A.H.min():.4f}, {A.H.max():.4f}] (round 6: [1.2768, 1.6359])",
            "RECORDED")
    _record("E-H", "post-fix", "[1.25, 1.50]",
            f"[{A.H[m].min():.4f}, {A.H[m].max():.4f}] "
            f"(round 6: [1.3256, 1.4919])", "RECORDED")
    _record("E-CF", "post-fix", "inside the derived band at every station",
            f"{len(rows)-n_out}/{len(rows)} inside (round 6: 60/69)",
            "RECORDED")

    out = [r for r in rows if not r["inside"]]
    print(f"  stations still outside: {len(out)}")
    for r in out:
        print(f"         Re_th={r['re_theta']:8.0f} H={r['H']:.4f} "
              f"dev={r['dev']:.3%} band={r['band']:.3%} "
              f"excess={r['dev']-r['band']:+.3%}")
    with open(os.path.join(RESULTS, "window_gccfix.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # F-DELTA -- what moved, against round 6's committed CSV
    old = {float(r["re_theta"]): r for r in csv.DictReader(
        open(os.path.join(R6, "results", "window.csv")))}
    ok_r = [(r, min(old, key=lambda k: abs(k - r["re_theta"]))) for r in rows]
    dcf = [abs(r["cf"]/float(old[k]["cf"]) - 1.0) for r, k in ok_r
           if abs(k - r["re_theta"]) / r["re_theta"] < 0.02]
    dH = [abs(r["H"]/float(old[k]["H"]) - 1.0) for r, k in ok_r
          if abs(k - r["re_theta"]) / r["re_theta"] < 0.02]
    _record("F-DELTA", "how far the fix moved round 6's readings", "RECORDED",
            f"c_f {100*min(dcf):.3f}..{100*max(dcf):.3f} %; "
            f"H {100*min(dH):.3f}..{100*max(dH):.3f} %", "RECORDED")

    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n== summary ==  {time.perf_counter()-t0:.2f} s")
    for t, mm, b, meas, v in SUMMARY:
        print(f"  {v:12s} [{t}] {mm} = {meas}")
    print("\n★ Round 6's verdict is NOT re-judged here; these are the same "
          "criteria re-measured after a transcription fix, recorded alongside.")
    print("★ The fix is justified by the source, not by whether E-CF improved.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
