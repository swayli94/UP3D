"""GS4.1 round 6 -- the turbulent closure under REPAIRED criteria (G1, G2).

Binding text: docs/dev_phase_four/20260819-1700-repaired-criteria-prereg.md
(committed before this script existed).

★★★ This round does NOT overturn round 5. Round 5's three T-FAILs stand as
recorded under the criteria registered there; this round asks a differently-posed
question with a different criterion. Passing here retracts nothing there, and the
verdict states both side by side.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_repaired_criteria/run.py
"""

import csv
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
R1 = os.path.join(ROOT, "bench", "studies", "gs41_strip_core")
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)

RHO, MU, U_INF = 1.0, 1.0e-5, 1.0
X0, X_TR, X1 = 0.05, 5.0, 400.0
RE_HI = 1.0e4                      # upper limit, carried over from round 5
ATTRACT_TOL = 0.01                 # E-ATTRACT: 1 %
SEED_PERTURB = 1.15                # >= 10 % apart, per G-SEEDS
CLAUSER_G = 6.8                    # E-HREF only, RECORDED (literature 6.1-6.8)
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def cf_coles_fernholz(ret):
    return 2.0 * ((1.0 / 0.384) * np.log(ret) + 4.127) ** -2


def cf_power_law(ret):
    return 0.024 * ret ** -0.25


def band_at(ret):
    a, b = cf_coles_fernholz(ret), cf_power_law(ret)
    s = abs(a - b) / (0.5 * (a + b))
    return max(0.03, 2.0 * s), s


# ---------------------------------------------------------------------------

def guards(C2, S):
    print("== guards ==")
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--",
                        "pyfp3d/viscous/closures.py", "pyfp3d/viscous/ibl3.py"],
                       cwd=ROOT, capture_output=True)
    print(f"  G-FROZEN  closures.py + ibl3.py unchanged: "
          f"{'PASS' if r.returncode == 0 else 'FAIL'}")
    if r.returncode:
        raise SystemExit("G-FROZEN failed")

    import pyfp3d.solve.newton as N
    import pyfp3d.solve.picard as PC

    def _f(nm):
        def _g(*a, **k):
            raise AssertionError(f"G-NOSOLVE: {nm} was CALLED")
        return _g
    n = 0
    for mod, names in ((N, ("solve_newton_lifting", "solve_newton_transonic")),
                       (PC, ("solve_subsonic", "solve_subsonic_lifting",
                             "solve_laplace", "solve_laplace_lifting"))):
        for nm in names:
            setattr(mod, nm, _f(nm)); n += 1
    print(f"  G-NOSOLVE {n} solver entry points stubbed  PASS")

    print(f"  G-SOURCE  GACON={C2.GACON} (xbl.f:1559)  GBCON={C2.GBCON} "
          f"(:1560)  CTCON={C2.CTCON:.12f} (:1569, derived)")
    assert C2.CTCON == 0.5 / (C2.GACON ** 2 * C2.GBCON)

    p = C2.packet_turb(2.0e-2, 1.4, U_INF, rho=RHO, mu=MU)
    di = C2.dissipation_identity(p["cD"], p["H_star"])
    assert abs(di - 2.0 * p["cD"] / p["H_star"]) < 1e-15
    print("  G-DI      2 c_D / H* == DI  PASS")

    h = C2.zpg_fixed_point()
    assert abs(h - 2.590433) < 1e-6
    print(f"  G-LEGACY  laminar ZPG fixed point {h:.6f} (round 3)  PASS")


# ---------------------------------------------------------------------------
# E-ATTRACT: does the solution forget its initial condition?  And the window
# start is WHERE it has -- measured, not chosen (G2).
# ---------------------------------------------------------------------------

def attract(S, C2, n_substep=8000):
    print("== E-ATTRACT + the measured window start (G1, G2) ==")
    # the physical post-transition state, from a full laminar->turbulent march
    y0 = C2.blasius_state(X0, ue=U_INF, rho=RHO, mu=MU, H=2.591100)
    lead = S.march_correlation(np.array([X_TR * 1.02]), y0, X0,
                               S.flat_plate_ue(U_INF), rho=RHO, mu=MU,
                               n_substep=n_substep, x_tr=X_TR)
    x_s, th_s, H_s = lead.x[0], lead.theta[0], lead.H[0]
    H_b = H_s * SEED_PERTURB
    sep = abs(H_b / H_s - 1.0)
    print(f"  G-SEEDS   post-transition state x={x_s:.3f} theta={th_s:.3e}; "
          f"seeds H = {H_s:.4f} and {H_b:.4f} -> {100*sep:.1f} % apart "
          f"({'PASS' if sep >= 0.10 else 'FAIL -- E-ATTRACT would be vacuous'})")
    if sep < 0.10:
        raise SystemExit("G-SEEDS failed")

    stations = np.geomspace(x_s * 1.02, X1, 120)
    kw = dict(rho=RHO, mu=MU, n_substep=n_substep, x_tr=X_TR)
    A = S.march_correlation(stations, (th_s, H_s), x_s, S.flat_plate_ue(U_INF), **kw)
    B = S.march_correlation(stations, (th_s, H_b), x_s, S.flat_plate_ue(U_INF), **kw)

    d = np.abs(B.H - A.H) / A.H
    conv = np.where(d <= ATTRACT_TOL)[0]
    # the window opens at the first station from which the two stay converged
    start = None
    for i in conv:
        if np.all(d[i:] <= ATTRACT_TOL):
            start = i
            break
    if start is None:
        _record("E-ATTRACT", "two seeds collapse onto one H(Re_theta)",
                f"<= {100*ATTRACT_TOL:.0f} % and staying", "never", "E-UNDEF")
        return A, B, None
    _record("E-ATTRACT", "two seeds collapse onto one H(Re_theta)",
            f"<= {100*ATTRACT_TOL:.0f} % and staying",
            f"achieved at Re_theta = {A.re_theta[start]:.0f} "
            f"(x = {A.x[start]:.2f}); final separation {d[-1]:.2e}",
            "E-ATTRACT PASS")
    print(f"  G-WINDOW  window start is MEASURED, not chosen: Re_theta "
          f"{A.re_theta[start]:.0f}; {start} of {len(stations)} stations "
          f"excluded as still carrying the initial condition; upper limit "
          f"Re_theta {RE_HI:.0e} carried over from round 5")
    return A, B, start


# ---------------------------------------------------------------------------

def gates(A, start):
    m = (np.arange(A.x.size) >= start) & (A.re_theta <= RE_HI)
    print(f"== E-CF / E-H on the measured window "
          f"(Re_theta {A.re_theta[m].min():.0f}..{A.re_theta[m].max():.0f}, "
          f"{int(m.sum())} stations) ==")

    phys = bool(np.all(A.H > 1.05) and np.all(A.H < 4.0))
    _record("E-PHYS", "H physical over the whole march", "1.05 < H < 4",
            f"[{A.H.min():.4f}, {A.H.max():.4f}]",
            "E-PHYS PASS" if phys else "E-FAIL")

    rows, n_out = [], 0
    for i in np.where(m)[0]:
        ret, cf = A.re_theta[i], A.cf[i]
        cfa, cfb = cf_coles_fernholz(ret), cf_power_law(ret)
        b, s = band_at(ret)
        if s > 0.10:
            continue
        dev = max(abs(cf/cfa - 1.0), abs(cf/cfb - 1.0))
        inside = dev <= b
        n_out += (not inside)
        # E-HREF: Clauser equilibrium H, RECORDED only
        href = 1.0/(1.0 - CLAUSER_G*np.sqrt(0.5*0.5*(cfa + cfb)))
        rows.append({"re_theta": ret, "x": A.x[i], "H": A.H[i], "cf": cf,
                     "cf_coles_fernholz": cfa, "cf_power_law": cfb,
                     "band": b, "corr_spread": s, "dev": dev,
                     "inside": int(inside), "H_ref_clauser": href,
                     "H_rel_to_ref": A.H[i]/href - 1.0})
    for r in rows[::max(1, len(rows)//6)]:
        print(f"         Re_th={r['re_theta']:8.0f} H={r['H']:.4f} "
              f"cf={r['cf']:.6f} dev={r['dev']:.3%} band={r['band']:.3%} "
              f"{'in' if r['inside'] else '★OUT'}")
    _record("E-CF", "c_f vs both established ZPG correlations",
            "inside the derived band at every station -- " + "RECORDED -- E-CF was DEMOTED from a gate on 2026-08-21 (user ruling). The window, the band and every threshold are UNCHANGED; only the verdict code moved. Grounds, all previously measured: making the closure more faithful to xblsys.f made agreement with these two correlations WORSE (62/69 -> 43/70, round 9 leg A); the plate's out-of-band Re_theta range has NO counterpart against a same-family reference where theta agrees to 0.11 % (round 14); and the correlations' own validity range has no citable source in this repo. ★ This does NOT mean the plate side is fine -- it means the plate's red cannot be ATTRIBUTED to the closure, because its reference is not qualified. Reopening condition: a ZPG reference whose validity range IS citable, or an experiment carrying BL profiles.",
            f"{len(rows)-n_out}/{len(rows)} inside",
            "RECORDED (demoted 2026-08-21; was E-FAIL at "
            f"{len(rows)-n_out}/{len(rows)})")

    Hw = A.H[m]
    _record("E-H", "H over the measured window", "[1.25, 1.50]",
            f"[{Hw.min():.4f}, {Hw.max():.4f}]",
            "E-H PASS" if (Hw.min() >= 1.25 and Hw.max() <= 1.50) else "E-FAIL")

    hr = np.array([r["H_rel_to_ref"] for r in rows])
    _record("E-HREF", f"H vs Clauser G={CLAUSER_G} equilibrium relation",
            "RECORDED, no gate (G carries a 6.1-6.8 literature spread)",
            f"{100*hr.min():+.2f}%..{100*hr.max():+.2f}%", "RECORDED")
    return rows


def cross(A, start):
    old = list(csv.DictReader(open(os.path.join(
        R1, "results", "turbulent_recorded.csv"))))
    lo = min(float(r["rel_diff"]) for r in old)
    hi = max(float(r["rel_diff"]) for r in old)
    m = (np.arange(A.x.size) >= start) & (A.re_theta <= RE_HI)
    pw = 0.027 * A.re_x[m] ** (-1.0/7.0)
    rel = (A.cf[m] - pw)/pw
    _record("E-CROSS", "vs the profile family's committed ZPG reading",
            "ONE-SIDED: supports only",
            f"{100*rel.min():+.1f}%..{100*rel.max():+.1f}% "
            f"(profile family {100*lo:+.1f}%..{100*hi:+.1f}%)", "RECORDED")


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    from pyfp3d.viscous import closures_2d as C2
    from pyfp3d.viscous import strip2d as S
    guards(C2, S)
    A, B, start = attract(S, C2)
    rows = gates(A, start) if start is not None else None
    if rows:
        cross(A, start)
        with open(os.path.join(RESULTS, "window.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)

    print(f"\n== summary ==  {time.perf_counter()-t0:.2f} s")
    for tag, metric, band, measured, verdict in SUMMARY:
        print(f"  {verdict:18s} [{tag}] {metric} = {measured}")
    print("\n★★★ Round 5's three T-FAILs STAND -- they were judged under the "
          "criteria registered there. This round asks a different question; it "
          "does not overturn them (pre-registration section 0).")
    print("★ And a pass here does NOT mean the turbulent closure is validated: "
          "a ZPG plate still cannot test pressure gradient or the lag.")
    fails = [r for r in SUMMARY if "FAIL" in r[4]]
    print(f"  {len(fails)} FAIL row(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
