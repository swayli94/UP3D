"""GS4.1 round 10 -- separate the three candidates behind the transition gap.

Binding text: phases/p4/docs/dev_phase_four/20260820-2000-transition-attribution-prereg.md
plus addenda #1 (the CT column is two quantities; L-TURB withdrawn) and #2
(TRDIF may be subsumed by a continuous march). All committed before this script.

★ Zero library change is EXPECTED. If a library change turns out to be needed,
that is itself a result to record (guard G-FROZEN-LIB).

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_transition_attrib/run.py
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

#: the round-9 machinery is IMPORTED, not retyped -- same XFOIL batch, same
#: column conversions, same stagnation split, same transition guard.
import run as R9                                                        # noqa: E402

RHO, MU = R9.RHO, R9.MU
BAND = 0.02                       # carried over from round 9's L-TURB, not re-picked
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def dev_table(side, st, idx, quants=("theta", "H", "cf", "ctau")):
    """Per-station |relative deviation| of a marched arm against XFOIL, on a
    FIXED index set (T-SET). Returns {q: array} aligned with idx."""
    out = {}
    for q in quants:
        xf = getattr(side, q)[idx]
        ours = getattr(st, "theta" if q == "theta" else q)[
            [int(np.argmin(np.abs(st.x - side.s[k]))) for k in idx]]
        with np.errstate(divide="ignore", invalid="ignore"):
            out[q] = np.where(xf != 0.0, np.abs(ours / xf - 1.0), np.nan)
    return out


def march_from(S, C2, side, i_seed, x_tr_s, n_substep, lag=True):
    stations = side.s[i_seed + 1:]
    ue_fn, _ = R9.ue_interp(side)
    y0 = (float(side.theta[i_seed]), float(side.H[i_seed]))
    return S.march_correlation(stations, y0, float(side.s[i_seed]), ue_fn,
                               rho=RHO, mu=MU, n_substep=n_substep,
                               x_tr=float(x_tr_s), lag=lag)


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    from pyfp3d.viscous import closures_2d as C2
    from pyfp3d.viscous import strip2d as S

    print("== guards ==")
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", "pyfp3d/"],
                       cwd=ROOT, capture_output=True)
    print(f"  G-FROZEN-LIB  pyfp3d/ unchanged this round: "
          f"{'PASS' if r.returncode == 0 else 'FAIL -- and that is a RESULT'}")
    print("  G-PROV        reference = the locally rebuilt XFOIL 6.99, used "
          "whole; round 7's three provenance facts carried in run.py's G-PROV")
    print("  G-REUSE       round 9's XFOIL batch, column conversions, stagnation "
          "split and G-XTR are IMPORTED from gs41_lag_xfoil/run.py, not retyped")

    with tempfile.TemporaryDirectory() as wd:
        surf, wake = R9.run_xfoil(wd)
    sa, sb, kstag = R9.split_at_stagnation(surf)

    rows = []
    for side in (sa, sb):
        i_tr = R9.transition_index(side)
        if not R9.guard_xtr(side, i_tr):
            raise SystemExit("G-XTR failed -- stopping (kill 1)")
        i_lam = i_tr - 1                       # the last fully laminar station
        s_lo, s_hi = float(side.s[i_lam]), float(side.s[i_tr])
        print(f"-- {side.name}: clean laminar seed at x/c {side.x[i_lam]:.5f} "
              f"(H {side.H[i_lam]:.4f}); the transition interval is "
              f"x/c {side.x[i_lam]:.5f}..{side.x[i_tr]:.5f}, arc length "
              f"{s_lo:.5f}..{s_hi:.5f}; the trip we set is x/c {R9.X_TRIP}")

        # T-SET: the comparison index set is fixed ONCE, from the widest x_tr,
        # so that moving x_tr cannot change which samples exist.
        idx = np.arange(i_tr, len(side.s))
        print(f"   T-SET  {len(idx)} stations from x/c {side.x[idx[0]]:.5f} to "
              f"{side.x[idx[-1]]:.5f}; fixed before any arm ran, so the sweep "
              f"cannot change the sample set")

        # --- T-C2: sweep x_tr across exactly one XFOIL interval -------------
        fr = np.linspace(0.0, 1.0, 11)
        for f in fr:
            s_tr = s_lo + f * (s_hi - s_lo)
            x_tr_xc = np.interp(s_tr, side.s, side.x)
            for nsub, lbl in ((4000, "fine"), (len(side.s), "one-per-interval")):
                st = march_from(S, C2, side, i_lam, s_tr, nsub)
                d = dev_table(side, st, idx)
                rows.append({"side": side.name, "arm": lbl, "f": f,
                             "x_tr": x_tr_xc, "n_substep": nsub,
                             **{f"n_out_{q}": int(np.nansum(d[q] > BAND))
                                for q in d},
                             **{f"worst_{q}": float(np.nanmax(d[q])) for q in d},
                             "dev_H_first": float(d["H"][0]),
                             "dev_theta_first": float(d["theta"][0])})
    with open(os.path.join(RESULTS, "sweep.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    verdicts(rows)
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n== summary ==  {time.perf_counter()-t0:.2f} s")
    for tag, metric, band, measured, verdict in SUMMARY:
        print(f"  {verdict:22s} [{tag}] {metric} = {measured}")
    return 0


def verdicts(rows):
    import numpy as np
    QS = ("theta", "H", "cf", "ctau")

    for side in sorted({r["side"] for r in rows}):
        fine = [r for r in rows if r["side"] == side and r["arm"] == "fine"]
        fine.sort(key=lambda r: r["f"])
        print(f"== {side} ==")
        print("   x_tr(x/c)   n_out theta/H/cf/ctau      worst H    dev at the "
              "first predicted station (theta, H)")
        for r in fine:
            print(f"   {r['x_tr']:.5f}    "
                  f"{r['n_out_theta']:3d}/{r['n_out_H']:3d}/{r['n_out_cf']:3d}/"
                  f"{r['n_out_ctau']:3d}        {100*r['worst_H']:6.2f} %    "
                  f"{100*r['dev_theta_first']:6.2f} %, "
                  f"{100*r['dev_H_first']:6.2f} %")

        # --- T-C1': can a clean seed + an interior x_tr put the transition
        # interval's own station inside the band?
        best = min(fine, key=lambda r: r["dev_H_first"])
        ok = best["dev_H_first"] <= 0.02 and best["dev_theta_first"] <= 0.02
        _record("T-C1'", f"{side}: the transition-interval station, PREDICTED "
                "from a clean laminar seed",
                "inside the round-9 band (2 %) for theta and H",
                f"best at x_tr = x/c {best['x_tr']:.5f}: theta "
                f"{100*best['dev_theta_first']:.2f} %, H "
                f"{100*best['dev_H_first']:.2f} %",
                "T-C1' PASS -- TRDIF is absorbed by the continuous march"
                if ok else "T-FAIL -- a continuous switch is NOT sufficient")

        # --- T-C2: how much of the outlier count does x_tr alone account for?
        n0 = fine[-1]                      # x_tr at the station (round 9's choice)
        tot = lambda r: sum(r[f"n_out_{q}"] for q in QS)
        lo = min(fine, key=tot)
        _record("T-C2", f"{side}: x_tr swept across ONE XFOIL interval",
                "derived range = the interval containing the trip; "
                "sample set fixed by T-SET",
                f"outliers {tot(n0)} at the station boundary -> {tot(lo)} at "
                f"x/c {lo['x_tr']:.5f} (best); the trip we set is x/c "
                f"{R9.X_TRIP}",
                "T-C2 RECORDED")

        # --- T-C3: bound the discretisation PER QUANTITY (addendum #3).
        # ★ The first execution divided the largest movement over ALL four
        # quantities by the worst deviation in H alone -- c_f and Ctau movement
        # against an H deviation, which is question 5 again.
        ratios = {}
        for q in QS:
            mv, base = 0.0, 0.0
            for a in fine:
                b = [r for r in rows if r["side"] == side
                     and r["arm"] == "one-per-interval" and r["f"] == a["f"]]
                if b:
                    mv = max(mv, abs(b[0][f"worst_{q}"] - a[f"worst_{q}"]))
                    base = max(base, a[f"worst_{q}"])
            ratios[q] = (mv, base, mv / base if base else float("nan"))
        worst_ratio = max(v[2] for v in ratios.values())
        excluded = worst_ratio < 0.1
        print("         T-C3 per quantity (movement fine -> ~one step per "
              "XFOIL interval, each against ITS OWN observed deviation):")
        for q, (mv, base, ra) in ratios.items():
            print(f"           {q:6s} movement {100*mv:6.2f} pp / observed "
                  f"{100*base:6.2f} %  = {ra:.3f}")
        print("         ★ A LOOSE bound is not a main effect: our coarse RK4 is "
              "NOT XFOIL's two-point implicit scheme, which is designed for "
              "exactly that resolution, so this over-states what the "
              "discretisation is worth. It can exclude; it cannot attribute.")
        _record("T-C3", f"{side}: discretisation, BOUNDED per quantity",
                "movement / that quantity's own observed deviation; "
                "<= 0.1 excludes (an upper bound, never 'C3 = 0')",
                "; ".join(f"{q} {v[2]:.3f}" for q, v in ratios.items()),
                "T-C3 EXCLUDED (upper bound)" if excluded
                else "T-C3 NOT EXCLUDED (the bound is loose, NOT a main effect)")


if __name__ == "__main__":
    sys.exit(main())
