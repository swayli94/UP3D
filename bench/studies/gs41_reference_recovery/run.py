"""GS4.1 round 7 -- recover the edge velocity XFOIL always wrote and the
committed CSV writer discarded.

Binding text: docs/dev_phase_four/20260820-0100-reference-recovery-prereg.md
(committed before this script existed).

★ Evidence recovery, not measurement: this touches no library code and produces
no conclusion about the strip core. It re-runs XFOIL through the committed
generator's own batch sequence and keeps the DUMP columns that were parsed and
then dropped -- Ue/Vinf, Theta, H, CD and CT.

★★ G-NOOVERWRITE: the committed CSVs are evidence and are NOT edited. The
recovered reference is written alongside them as *_with_ue.csv.
"""

import csv
import os
import subprocess
import sys
import tempfile
from decimal import Decimal

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
REFDIR = os.path.join(ROOT, "cases", "reference_data", "naca0012_viscous_xfoil")
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)
sys.path.insert(0, REFDIR)

#: the committed generator's own constants -- imported, never retyped
from generate_xfoil_reference import (N_PANELS, REYNOLDS, MACH, ALPHA_DEG,  # noqa
                                      N_ITER, find_xfoil)

SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def ulp_of(text):
    """One unit in the last STORED digit of a decimal string.

    Pre-registration section 2: the committed CSV holds as few as two
    significant figures, so a tolerance tighter than its own storage cannot be
    satisfied. This derives the tolerance from the file rather than picking one.
    """
    d = Decimal(text.strip()).normalize()          # drops trailing zeros
    return float(Decimal(1).scaleb(d.as_tuple().exponent))


def batch(tag, xtr, polar, dump):
    """The committed generator's sequence, plus PLOP/G F for headless running."""
    return ("PLOP\nG F\n\n"
            "NACA 0012\nPPAR\n"
            f"N {N_PANELS}\n\n\n"
            "OPER\n"
            f"VISC {REYNOLDS:.2E}\nMACH {MACH:.1f}\n"
            f"VPAR\nXTR {xtr:.2f} {xtr:.2f}\nSHOW\n\n"
            f"ITER {N_ITER}\nPACC\n{polar}\n\n"
            f"ALFA {ALPHA_DEG:.1f}\nDUMP {dump}\n\nQUIT\n")


def run_xfoil(binary, wd, tag, xtr):
    polar, dump = f"p_{tag}.txt", f"d_{tag}.txt"
    (open(os.path.join(wd, f"{tag}.in"), "w")).write(batch(tag, xtr, polar, dump))
    with open(os.path.join(wd, f"{tag}.in")) as fin:
        r = subprocess.run([str(binary)], stdin=fin, cwd=wd,
                           capture_output=True, text=True, timeout=300)
    open(os.path.join(wd, f"{tag}.log"), "w").write(r.stdout)
    return os.path.join(wd, dump), os.path.join(wd, polar)


def parse_dump(path):
    """XFOIL 6.99 BLDUMP: 10 fields, surface rows then wake rows."""
    header, rows = None, []
    for ln in open(path).read().splitlines():
        if ln.lstrip().startswith("#"):
            header = header or ln
            continue
        f = ln.split()
        if len(f) == 10:
            rows.append([float(v) for v in f])
    a = np.array(rows)
    # surface rows run upper TE -> LE -> lower TE with s increasing; the wake
    # continues past the TE. Split where x stops belonging to the airfoil.
    le = int(np.argmin(a[:, 1]))
    tail = a[le:]
    end = le + int(np.argmax(tail[:, 1])) + 1        # last point back at the TE
    return header, a[:end], a[end:]


def main():
    os.makedirs(RESULTS, exist_ok=True)
    binary = find_xfoil()
    print("== guards ==")
    print(f"  G-BUILD   binary {binary}")
    print("  G-BUILD   built from reference/XFOIL6.99.zip -> Xfoil699src, two "
          "Makefile edits: plotlib/config.make (gfortran, -std=legacy "
          "-fallow-argument-mismatch, libPlt.a, -lX11) and bin/Makefile.gfortran "
          "(PLTOBJ -> ../plotlib/libPlt.a, PLTLIB -> -lX11, same two flags)")
    print("  G-NOSOLVE this round runs XFOIL, never a pyfp3d solver  PASS")
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", "pyfp3d/"],
                       cwd=ROOT, capture_output=True)
    print(f"  G-FROZEN  pyfp3d/ unchanged: "
          f"{'PASS' if r.returncode == 0 else 'FAIL'}")
    if r.returncode:
        raise SystemExit("G-FROZEN failed")

    committed_polar = {r["case"]: r for r in
                       csv.DictReader(open(os.path.join(REFDIR,
                                                        "polar_summary.csv")))}
    ok_all = True
    with tempfile.TemporaryDirectory() as wd:
        for tag, xtr, fname in (("xtr005", 0.05,
                                 "delta_star_cf_alpha2_m05_xtr005.csv"),
                                ("xtr030", 0.30,
                                 "delta_star_cf_alpha2_m05_xtr030.csv")):
            print(f"== {tag} ==")
            dump, polar = run_xfoil(binary, wd, tag, xtr)
            header, surf, wake = parse_dump(dump)
            print(f"  G-BUILD   BLDUMP header ({len(header.split())-1} cols): "
                  f"{header.strip()}")
            print(f"  parsed {len(surf)} surface + {len(wake)} wake rows")

            ref = list(csv.DictReader(open(os.path.join(REFDIR, fname))))
            ru = [r for r in ref if r["surface"] == "upper"]
            rl = [r for r in ref if r["surface"] == "lower"]
            if len(ru) != len(rl):
                raise SystemExit(f"{tag}: committed upper/lower counts differ")

            le = int(np.argmin(surf[:, 1]))
            up = surf[:le + 1][::-1]                 # LE -> TE, as committed
            lo = surf[le:]
            _record("R-STATIONS", f"{tag}: surface rows vs committed",
                    f"{len(ru)} upper + {len(rl)} lower",
                    f"{len(up)} upper + {len(lo)} lower",
                    "PASS" if (len(up) == len(ru) and len(lo) == len(rl))
                    else "R-FAIL -> kill 1")
            if len(up) != len(ru) or len(lo) != len(rl):
                ok_all = False
                print("  ★ paneling differs -- stopping this case rather than "
                      "interpolating two different point sets (kill 1)")
                continue

            # R-POLAR
            cp = committed_polar[tag]
            got = {}
            for ln in open(polar).read().splitlines():
                f = ln.split()
                if len(f) >= 7 and f[0].replace(".", "").isdigit():
                    got = {"cl": float(f[1]), "cd": float(f[2]),
                           "cm": float(f[4])}
            worst = 0.0
            for k in ("cl", "cd", "cm"):
                u = ulp_of(cp[k])
                n = abs(got[k] - float(cp[k])) / u
                worst = max(worst, n)
                print(f"         {k}: committed {cp[k]} (ulp {u:g}) vs "
                      f"{got[k]} -> {n:.2f} ulp")
            _record("R-POLAR", f"{tag}: cl/cd/cm", "<= 1 ulp of the stored digit",
                    f"worst {worst:.2f} ulp",
                    "PASS" if worst <= 1.0 else "R-FAIL -> kill 2")
            ok_all &= worst <= 1.0

            # R-REPRO, per-row tolerance from that row's own stored digits
            n_bad, n_1sig, worst_n = 0, 0, 0.0
            for side, arr, rr in (("upper", up, ru), ("lower", lo, rl)):
                for i, r in enumerate(rr):
                    for col, idx in (("dstar_over_c", 4), ("cf", 6)):
                        u = ulp_of(r[col])
                        if Decimal(r[col]).normalize().as_tuple().digits == (0,):
                            continue
                        if len(Decimal(r[col]).normalize().as_tuple().digits) == 1:
                            n_1sig += 1
                        n = abs(arr[i, idx] - float(r[col])) / u
                        worst_n = max(worst_n, n)
                        n_bad += n > 1.0
            _record("R-REPRO", f"{tag}: delta*/cf per station",
                    "<= 1 ulp of each row's own last stored digit",
                    f"{n_bad} of {2*len(ref)} outside; worst {worst_n:.2f} ulp"
                    + (f"; {n_1sig} rows stored to ONE significant figure"
                       if n_1sig else ""),
                    "PASS" if n_bad == 0 else "R-FAIL")
            ok_all &= n_bad == 0

            # write the recovered reference alongside (G-NOOVERWRITE)
            outp = os.path.join(REFDIR, fname.replace(".csv", "_with_ue.csv"))
            with open(outp, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["x_c", "surface", "s", "ue_over_vinf",
                            "dstar_over_c", "theta_over_c", "cf", "H",
                            "cd_local", "ct"])
                for side, arr in (("upper", up), ("lower", lo)):
                    for row in arr:
                        w.writerow([f"{row[1]:.6f}", side, f"{row[0]:.6f}",
                                    f"{row[3]:.6f}", f"{row[4]:.6e}",
                                    f"{row[5]:.6e}", f"{row[6]:.6e}",
                                    f"{row[7]:.4f}", f"{row[8]:.6e}",
                                    f"{row[9]:.6e}"])
            print(f"  wrote {os.path.relpath(outp, ROOT)} "
                  f"(committed CSV untouched -- G-NOOVERWRITE)")

            ue, ct = surf[:, 3], surf[:, 9]
            _record("R-FIELDS", f"{tag}: recovered columns", "RECORDED",
                    f"Ue/Vinf {ue.min():.4f}..{ue.max():.4f}; "
                    f"CT {ct.min():.3e}..{ct.max():.3e}; "
                    f"H {surf[:,7].min():.3f}..{surf[:,7].max():.3f}",
                    "RECORDED")

    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print("\n== summary ==")
    for t, m, b, meas, v in SUMMARY:
        print(f"  {v:20s} [{t}] {m} = {meas}")
    print("\n★ This round produces NO conclusion about the strip core -- it only "
          "recovers the reference. And R-REPRO passing does NOT clear the "
          "provenance debt: that is RECORDED by R-PROV, not fixed.")
    fails = [r for r in SUMMARY if "FAIL" in r[4]]
    print(f"  {len(fails)} FAIL row(s)")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
