"""GS4.1 round 11 -- XFOIL's two-point scheme, so C3 becomes a reading.

Binding text: phases/p4/docs/dev_phase_four/20260821-0100-two-point-scheme-prereg.md +
addendum #1 (20260821-0130). Both committed before any code.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_two_point/run.py
"""

import csv
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "bench", "studies", "gs41_lag_xfoil"))

import run as R9                                                   # noqa: E402
import scheme as TP                                                # noqa: E402

RHO, MU, BAND = R9.RHO, R9.MU, 0.02
XSRC = os.environ.get("GS41_XFOIL_SRC", os.path.join(
    "/tmp/claude-1000/-home-lrz-codes-UP3D/"
    "3cd995ca-fd78-494f-80eb-3e7fd575323a/scratchpad/xfsrc/Xfoil699src/src"))
SUMMARY = []

#: G-WHOLE2 -- every plain assignment in BLDIF (xblsys.f:1554-1979), classified.
#: Unclassified => FAIL. Same instrument as round 9's A-WHOLE, moved to BLDIF.
BLDIF_TABLE = {
    ("XLOG", 1571): ("N/A-SIMI", "ITYP=0 similarity station; we seed from XFOIL"),
    ("ULOG", 1572): ("N/A-SIMI", "ITYP=0"),
    ("TLOG", 1573): ("N/A-SIMI", "ITYP=0"),
    ("HLOG", 1574): ("N/A-SIMI", "ITYP=0"),
    ("DDLOG", 1575): ("N/A-SIMI", "ITYP=0"),
    ("XLOG", 1578): ("HERE", "residuals(): XLOG"),
    ("ULOG", 1579): ("HERE", "residuals(): ULOG"),
    ("TLOG", 1580): ("HERE", "residuals(): TLOG"),
    ("HLOG", 1581): ("HERE", "residuals(): HLOG"),
    ("DDLOG", 1586): ("N/A-DERIV", "multiplies only the Z_ derivative rows"),
    ("HUPWT", 1601): ("HERE", "scheme.HUPWT"),
    ("HDCON", 1607): ("HERE", "scheme.HDCON_WALL"),
    ("HDCON", 1613): ("N/A-WAKE", "ITYP=3"),
    ("ARG", 1620): ("HERE", "residuals(): ARG"),
    ("HL", 1621): ("HERE", "residuals(): HL"),
    ("HLSQ", 1630): ("HERE", "residuals(): HLSQ, cap 15"),
    ("EHH", 1631): ("HERE", "residuals(): the exp"),
    ("UPW", 1632): ("HERE", "residuals(): UPW"),
    ("REZC", 1666): ("N/A-LAM", "amplification equation; transition is imposed"),
    ("SA", 1690): ("HERE", "lag block"),
    ("CQA", 1691): ("HERE", "lag block"),
    ("CFA", 1692): ("HERE", "lag block"),
    ("HKA", 1693): ("HERE", "★ UPWIND-weighted HKA, distinct from BLMID's mean"),
    ("USA", 1695): ("HERE", "lag block"),
    ("RTA", 1696): ("HERE", "lag block"),
    ("DEA", 1697): ("HERE", "lag block"),
    ("DA", 1698): ("HERE", "lag block"),
    ("ALD", 1703): ("N/A-WAKE", "ALD = DLCON in the wake"),
    ("ALD", 1705): ("HERE", "ALD = 1 on a wall layer"),
    ("GCC", 1710): ("HERE", "GCCON on the turbulent wall"),
    ("HKC", 1711): ("HERE", "HKC"),
    ("HKC", 1715): ("HERE", "the 0.01 floor"),
    ("GCC", 1720): ("N/A-LAM", "laminar branch"),
    ("HKC", 1721): ("N/A-LAM", "laminar branch"),
    ("HR", 1726): ("HERE", "HR"),
    ("UQ", 1730): ("HERE", "UQ"),
    ("SCC", 1759): ("HERE", "SCC"),
    ("SLOG", 1766): ("HERE", "SLOG"),
    ("DXI", 1767): ("HERE", "DXI"),
    ("REZC", 1769): ("HERE", "the lag residual"),
    ("HA", 1846): ("HERE", "momentum HA"),
    ("MA", 1847): ("N/A-M0", "edge Mach = 0"),
    ("XA", 1848): ("HERE", "momentum XA"),
    ("TA", 1849): ("HERE", "momentum TA"),
    ("HWA", 1850): ("N/A-M0", "no wall blowing => DW = 0"),
    ("CFX", 1853): ("HERE", "momentum CFX, with BLMID's midpoint CFM"),
    ("BTMP", 1865): ("HERE", "momentum BTMP"),
    ("REZT", 1867): ("HERE", "the momentum residual"),
    ("XOT1", 1905): ("HERE", "X1/T1"),
    ("XOT2", 1906): ("HERE", "X2/T2"),
    ("HA", 1908): ("HERE", "energy HA"),
    ("HSA", 1909): ("HERE", "energy HSA"),
    ("HCA", 1910): ("N/A-M0", "density thickness HC = MSQ*(...) = 0"),
    ("HWA", 1911): ("N/A-M0", "no wall blowing"),
    ("DIX", 1913): ("HERE", "energy DIX"),
    ("CFX", 1914): ("HERE", "energy CFX"),
    ("BTMP", 1918): ("HERE", "energy BTMP"),
    ("REZH", 1920): ("HERE", "the kinetic-energy residual"),
}


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def g_whole2():
    path = os.path.join(XSRC, "xblsys.f")
    if not os.path.exists(path):
        _record("G-WHOLE2", "every plain assignment in BLDIF is classified",
                "no unclassified pair", f"source absent at {path}", "UNDEFINED")
        return None
    src = open(path).read().splitlines()
    pat = re.compile(r"^\s{6,9}([A-Z][A-Z0-9]*)\s*=")
    rows, unc = [], []
    for i in range(1554, 1979):
        l = src[i - 1]
        if l[:1] in ("C", "c", "#", "*"):
            continue
        m = pat.match(l)
        if not m:
            continue
        key = (m.group(1), i)
        cls, note = BLDIF_TABLE.get(key, ("UNCLASSIFIED", ""))
        rows.append({"var": key[0], "line": key[1], "class": cls,
                     "maps_to": note})
        if cls == "UNCLASSIFIED":
            unc.append(key)
    n = {c: sum(r["class"] == c for r in rows) for c in
         ("HERE", "N/A-SIMI", "N/A-WAKE", "N/A-LAM", "N/A-M0", "N/A-DERIV",
          "UNCLASSIFIED")}
    print("  " + "  ".join(f"{k}={v}" for k, v in n.items() if v))
    _record("G-WHOLE2", "every plain assignment in BLDIF is classified",
            "no unclassified (name, line) pair",
            f"{len(rows)} assignments, {n['UNCLASSIFIED']} unclassified"
            + (f" {unc}" if unc else ""),
            "G-WHOLE2 PASS" if not unc else "G-FAIL -> kill 1")
    return rows


def refine(s, k):
    """Subdivide each interval into k, linearly in arc length."""
    if k == 1:
        return s.copy()
    out = [s[0]]
    for a, b in zip(s[:-1], s[1:]):
        out += list(np.linspace(a, b, k + 1)[1:])
    return np.asarray(out)


#: X-CONSIST's edge velocity -- ANALYTIC, and deliberately NOT zero-pressure
#: gradient. Addendum #2: on the real u_e this criterion cannot be satisfied by
#: any correct transcription, because the two schemes consume u_e differently
#: (station endpoints against the PCHIP's derivative) and a C1 interpolant's own
#: error does not shrink when the STATIONS are refined. m != 0 keeps every ULOG
#: term alive, which was addendum #1's reason for leaving the synthetic case.
CONSIST_M, CONSIST_U0 = -0.08, 1.10


def x_consist(S, side, i_tr):
    """★★ The independent oracle, and a kill criterion. The two marchers share
    `closures_2d` and share NO discretisation code, so if the difference between
    them shrinks under refinement they are discretising the same equations. If
    it does not, my BLDIF transcription is wrong and the rest of the round could
    not be interpreted."""
    base = side.s[i_tr:]
    x0 = float(base[0])
    ue_p = lambda x: CONSIST_U0 * (x / x0) ** CONSIST_M
    ue_fn = lambda x: (CONSIST_U0 * (x / x0) ** CONSIST_M,
                       CONSIST_U0 * CONSIST_M / x0 * (x / x0) ** (CONSIST_M - 1))
    y0 = (float(side.theta[i_tr]),
          float(side.theta[i_tr] * side.H[i_tr]),
          float(np.sqrt(max(side.ctau[i_tr], 1e-12))))
    print(f"         analytic u_e = {CONSIST_U0} (x/{x0:.5f})^{CONSIST_M} -- "
          f"non-ZPG so every ULOG term is alive, and no interpolant (addendum #2)")
    diffs = []
    for k in (1, 2, 4, 8, 16):
        st = refine(base, k)
        tp = TP.march(st[1:], y0, st[0], ue_p, RHO, MU)
        rk = S.march_correlation(st[1:], (y0[0], y0[1] / y0[0], y0[2]),
                                 st[0], ue_fn, rho=RHO, mu=MU,
                                 n_substep=max(8000, 60 * st.size),
                                 x_tr=float(st[0]) * 0.999, lag=True)
        j = [int(np.argmin(np.abs(st[1:] - v))) for v in base[1:]]
        d = max(float(np.max(np.abs(tp["theta"][j] / rk.theta[j] - 1.0))),
                float(np.max(np.abs(tp["H"][j] / rk.H[j] - 1.0))))
        diffs.append((k, d))
        print(f"         refine x{k:<2d} ({st.size:4d} stations): "
              f"max |two-point - RK4| / RK4 = {d:.3e}")
    ok = all(b[1] < a[1] for a, b in zip(diffs, diffs[1:]))
    order = np.polyfit(np.log([1.0 / k for k, _ in diffs]),
                       np.log([d for _, d in diffs]), 1)[0]
    _record("X-CONSIST", f"{side.name}: two-point vs RK4 under refinement",
            "difference must DECREASE monotonically; order reported",
            f"{diffs[0][1]:.3e} -> {diffs[-1][1]:.3e} over x1..x16, "
            f"measured order {order:.2f}",
            "X-CONSIST PASS" if ok else "X-FAIL -> kill 1")
    return ok, order


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
    whole = g_whole2()

    with tempfile.TemporaryDirectory() as wd:
        surf, wake = R9.run_xfoil(wd)
    sa, sb, _ = R9.split_at_stagnation(surf)

    # G-CLOSURE: one sentinel must move BOTH marchers, or they are not running
    # the same equations.
    side, i_tr = sa, R9.transition_index(sa)
    R9.guard_xtr(side, i_tr)
    ue_fn, ue_p = R9.ue_interp(side)
    st = side.s[i_tr:]
    y3 = (float(side.theta[i_tr]), float(side.theta[i_tr] * side.H[i_tr]),
          float(np.sqrt(max(side.ctau[i_tr], 1e-12))))
    y2 = (y3[0], y3[1] / y3[0], y3[2])

    def both():
        a = TP.march(st[1:], y3, st[0], lambda x: float(ue_p(x)), RHO, MU)
        b = S.march_correlation(st[1:], y2, st[0], ue_fn, rho=RHO, mu=MU,
                                n_substep=8000, x_tr=float(st[0]) * 0.999,
                                lag=True)
        return float(a["H"][-1]), float(b.H[-1])
    on = both()
    keep = C2.CD_OUT_US
    C2.CD_OUT_US = 1.0
    try:
        off = both()
    finally:
        C2.CD_OUT_US = keep
    moved = (abs(on[0] - off[0]) > 0, abs(on[1] - off[1]) > 0)
    _record("G-CLOSURE", "one closure sentinel moves BOTH marchers",
            "CD_OUT_US 0.995 -> 1.0 must change two-point AND RK4",
            f"two-point {on[0]:.8f} -> {off[0]:.8f}; RK4 {on[1]:.8f} -> "
            f"{off[1]:.8f}",
            "G-CLOSURE PASS" if all(moved) else "G-FAIL -> kill 1")
    if not all(moved):
        raise SystemExit("G-CLOSURE failed")

    print("== X-CONSIST (the independent oracle; kill criterion) ==")
    ok, order = x_consist(S, side, i_tr)
    if not ok:
        raise SystemExit("X-CONSIST failed -- the transcription is wrong, "
                         "stopping before interpreting anything (kill 1)")

    print("== X-ANSWER / X-C3 ==")
    rows = []
    for side in (sa, sb):
        i_tr = R9.transition_index(side)
        ue_fn, ue_p = R9.ue_interp(side)
        st = side.s[i_tr:]
        y3 = (float(side.theta[i_tr]), float(side.theta[i_tr] * side.H[i_tr]),
              float(np.sqrt(max(side.ctau[i_tr], 1e-12))))
        tp = TP.march(st[1:], y3, st[0], lambda x: float(ue_p(x)), RHO, MU)
        rk = S.march_correlation(st[1:], (y3[0], y3[1] / y3[0], y3[2]), st[0],
                                 ue_fn, rho=RHO, mu=MU, n_substep=8000,
                                 x_tr=float(st[0]) * 0.999, lag=True)
        idx = np.arange(i_tr + 1, len(side.s))
        for q, a, b in (("theta", tp["theta"], rk.theta),
                        ("H", tp["H"], rk.H), ("cf", tp["cf"], rk.cf),
                        ("ctau", tp["ctau"], rk.ctau)):
            xf = getattr(side, q)[idx]
            for j, k in enumerate(idx):
                rows.append({"side": side.name, "q": q, "x": side.x[k],
                             "xfoil": xf[j], "two_point": a[j], "rk4": b[j],
                             "dev_tp": abs(a[j] / xf[j] - 1.0) if xf[j] else np.nan,
                             "dev_rk": abs(b[j] / xf[j] - 1.0) if xf[j] else np.nan,
                             "c3": abs(a[j] / b[j] - 1.0)})
    verdicts(rows, order)
    with open(os.path.join(RESULTS, "stations.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    if whole:
        with open(os.path.join(RESULTS, "bldif_index.csv"), "w",
                  newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(whole[0].keys()))
            w.writeheader(); w.writerows(whole)
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n== summary ==  {time.perf_counter()-t0:.2f} s")
    for tag, metric, band, measured, verdict in SUMMARY:
        print(f"  {verdict:22s} [{tag}] {metric} = {measured}")
    print("\n★ The two-point scheme is an INSTRUMENT with a registered expiry "
          "(this verdict). It lives in bench/, never in pyfp3d/.")
    print("★ Nothing here validates the closure: XFOIL is the same "
          "Drela-Giles family, so agreement means the TRANSCRIPTION is exact.")
    print("★ This round establishes NOTHING about laminar discretisation -- "
          "laminar intervals are not implemented (addendum #1).")
    return 0


def verdicts(rows, order):
    for side in sorted({r["side"] for r in rows}):
        print(f"== {side} ==")
        print("   quantity   two-point vs XFOIL      RK4 vs XFOIL        "
              "C3 = |two-point - RK4|")
        bad = []
        for q in ("theta", "H", "cf", "ctau"):
            z = [r for r in rows if r["side"] == side and r["q"] == q]
            dt = np.array([r["dev_tp"] for r in z])
            dr = np.array([r["dev_rk"] for r in z])
            c3 = np.array([r["c3"] for r in z])
            nt, nr = int(np.sum(dt <= BAND)), int(np.sum(dr <= BAND))
            print(f"   {q:8s}  {nt:3d}/{len(z)} in, worst "
                  f"{100*np.nanmax(dt):6.2f} %   {nr:3d}/{len(z)} in, worst "
                  f"{100*np.nanmax(dr):6.2f} %   median {100*np.median(c3):5.2f} %"
                  f", worst {100*np.nanmax(c3):6.2f} %")
            if nt != len(z):
                bad.append(q)
        _record("X-ANSWER", f"{side}: two-point on XFOIL's own stations, "
                "PER STATION", f"every station <= {100*BAND:.0f} %",
                "; ".join(
                    f"{q} {int(np.sum([r['dev_tp'] for r in rows if r['side']==side and r['q']==q] <= np.float64(BAND)))}"
                    for q in ()) or
                "; ".join(
                    f"{q} {sum(1 for r in rows if r['side']==side and r['q']==q and r['dev_tp']<=BAND)}/"
                    f"{sum(1 for r in rows if r['side']==side and r['q']==q)} inside"
                    for q in ("theta", "H", "cf", "ctau")),
                "X-ANSWER PASS -- what remains IS the discretisation, and the "
                "closure transcription is exact" if not bad
                else f"X-FAIL on {bad} -- what remains is NOT discretisation")
        c3all = np.array([r["c3"] for r in rows if r["side"] == side])
        _record("X-C3", f"{side}: the discretisation, MEASURED not bounded",
                "RECORDED (round 10 could only bound it at 0.21-0.42)",
                f"median {100*np.median(c3all):.3f} %, worst "
                f"{100*np.nanmax(c3all):.2f} %", "X-C3 RECORDED")


if __name__ == "__main__":
    sys.exit(main())
