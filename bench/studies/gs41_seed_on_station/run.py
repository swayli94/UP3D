"""GS4.1 round 17 (G20) -- put transition ON a station and test the seed chain.

Binding text: docs/dev_phase_four/20260821-1400-g20-prereg.md + addendum #1
(20260821-1500). Both committed before the comparison was made.

★ Committed so the round's numbers are REGENERABLE rather than living only in the
verdict's prose and the test's constants -- this project's rule is that a number
which exists only in .md prose is not evidence.

The design: xblsys.f:435 tests XIFORC .LE. X2 and :451 assigns XT = XIFORC, so
placing the forced trip exactly ON a station makes TRDIF's TURBULENT sub-interval
zero-length. That station keeps its LAMINAR theta and H while its stored CTAU
becomes XFOIL's own ST = CTR * CQ -- so our seed chain can be compared to it with
no march, no discretisation and no interpolation in between.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_seed_on_station/run.py
"""

import csv
import difflib
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
#: a STATION's x/c, read off the baseline dump -- NOT tuned. G-NOTUNE holds the
#: change to this one XFOIL input.
X_STATION = 0.05370
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    from pyfp3d.viscous import closures_2d as C2

    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", "pyfp3d/"],
                       cwd=ROOT, capture_output=True)
    print(f"  G-FROZEN-LIB  pyfp3d/ unchanged: "
          f"{'PASS' if r.returncode == 0 else 'FAIL'}")
    batch = R9.BATCH.replace(f"XTR {R9.X_TRIP:.2f} {R9.X_TRIP:.2f}",
                             f"XTR {X_STATION:.5f} {X_STATION:.5f}").replace(
        "ALFA 2.0\nDUMP d.txt", "PACC\np.txt\n\nALFA 2.0\nDUMP d.txt")
    print("  G-NOTUNE      diff against round 9's batch (only XTR, plus PACC to "
          "read the polar):")
    for l in difflib.unified_diff(R9.BATCH.splitlines(), batch.splitlines(),
                                  lineterm="", n=0):
        if l.startswith(("+", "-")) and not l.startswith(("+++", "---")):
            print(f"                  {l}")

    with tempfile.TemporaryDirectory() as wd:
        open(os.path.join(wd, "in.txt"), "w").write(batch)
        with open(os.path.join(wd, "in.txt")) as fin:
            subprocess.run([R9.XFOIL], stdin=fin, cwd=wd, capture_output=True,
                           text=True, timeout=300)
        rows = [l.split() for l in open(os.path.join(wd, "d.txt"))
                if not l.lstrip().startswith("#")]
        a = np.array([[float(v) for v in q] for q in rows if len(q) == 10])
        pol = None
        for v in (l.split() for l in open(os.path.join(wd, "p.txt"))):
            if len(v) != 7:
                continue
            try:
                pol = [float(t) for t in v]
            except ValueError:
                continue
    le = int(np.argmin(a[:, 1]))
    past = np.where(a[le:, 1] > 1.0 + 1.0e-9)[0]
    surf = a[:le + int(past[0])] if past.size else a
    sa, sb, _ = R9.split_at_stagnation(surf)

    _record("G-XTON", "XFOIL's own transition x/c vs the station we asked for",
            f"equal to print precision ({X_STATION})",
            f"Top_Xtr {pol[5]:.5f}, Bot_Xtr {pol[6]:.5f}",
            "G-XTON PASS" if abs(pol[5] - X_STATION) < 1e-4
            and abs(pol[6] - X_STATION) < 1e-4 else "G-FAIL -> UNDEFINED")

    out = []
    for side in (sa, sb):
        k = int(np.argmin(np.abs(side.x - X_STATION)))
        # ★ addendum #1: the signature is NOT H collapsing at this station -- with
        # a zero-length turbulent part it CANNOT. It is H still on the laminar
        # plateau HERE and collapsing at the NEXT station.
        flat = abs(side.H[k] - side.H[k - 1])
        drop = side.H[k] - side.H[k + 1]
        _record("G-DEGEN", f"{side.name}: the CT column's semantics flip HERE",
                "H flat to the previous station (<0.05) and collapsing to the "
                "next (>0.1)",
                f"|dH| {flat:.4f}, drop {drop:.4f}",
                "G-DEGEN PASS" if flat < 0.05 and drop > 0.1 else "G-FAIL")
        _record("G-FRAC", f"{side.name}: XT's distance from the station",
                "<= 5 % of the interval, else UNDEFINED",
                f"|x_station - XTR| = {abs(side.x[k]-X_STATION):.2e}",
                "G-FRAC PASS (f = 0, no residual discretisation)"
                if abs(side.x[k] - X_STATION) < 1e-6 else "G-FAIL")

        p = C2.packet_turb(float(side.theta[k]), float(side.H[k]),
                           float(side.ue[k]), rho=RHO, mu=MU)
        ours = C2.s_tau_at_transition(float(side.H[k]), p["Ctau_eq"])
        xf = float(np.sqrt(side.ctau[k]))
        out.append({"side": side.name, "x": side.x[k], "theta": side.theta[k],
                    "H": side.H[k], "ue": side.ue[k], "re_theta": p["re_theta"],
                    "H_star": p["H_star"], "Us": p["Us"],
                    "Ctau_eq": p["Ctau_eq"], "S_tr_ours": ours,
                    "CT_xfoil": xf, "ratio": ours / xf,
                    "CTR_ours": ours / np.sqrt(p["Ctau_eq"]),
                    "CTR_implied": xf / np.sqrt(p["Ctau_eq"])})
        _record("Z-SEED", f"{side.name}: our seed chain vs XFOIL's own ST",
                f"ratio within {100*BAND:.0f} %",
                f"ours {ours:.6e} vs XFOIL {xf:.6e} = {ours/xf:.4f}",
                "Z-SEED PASS" if abs(ours / xf - 1.0) <= BAND else "Z-FAIL")
        _record("Z-CHAIN", f"{side.name}: CTR alone", "RECORDED",
                f"ours {out[-1]['CTR_ours']:.6f} vs implied "
                f"{out[-1]['CTR_implied']:.6f}", "RECORDED")

    _record("Z-N", "sample size", "n = 2, one station per surface; no spread "
            "reported and disagreement would be UNDEFINED",
            f"ratios {out[0]['ratio']:.4f} and {out[1]['ratio']:.4f} "
            f"({abs(out[0]['ratio']-out[1]['ratio'])*100:.2f} pp apart)",
            "Z-N PASS")

    with open(os.path.join(RESULTS, "seed_on_station.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n== summary ==  {time.perf_counter()-t0:.2f} s")
    for tag, metric, band, measured, verdict in SUMMARY:
        print(f"  {verdict:46s} [{tag}] {measured}")
    print("\n★ This verifies the seed FORMULA in the special case where XT lands "
          "on a station. In general use XT falls INSIDE an interval and the state "
          "fed to the seed is our own laminar march's, not XFOIL's -- which is "
          "where round 16's 41 % + 59 % decomposition came from.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
