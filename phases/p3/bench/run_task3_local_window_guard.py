"""Guard for the ONE new generator capability task 3 needs: a local chordwise refinement window.

Pre-registration §5.3 requires this knob to pass a guard before any reading built on it counts.

★ The guard is NOT bit-identity, and that is a deliberate design choice with a cost, recorded here
rather than glossed. `airfoil_surface_distribution(..., local_window, local_factor)` warps the
arclength parameter by a density function and keeps the STATION COUNT FIXED, because a leg that
refined the shock band by ADDING stations would change the global DOF count too and could then not
be separated from RG -- the exact confound phase two could never escape. Holding the count fixed
buys that separation and pays for it by making the region OUTSIDE the window coarser: points are
moved in, not created. So what this guard certifies is (i) the count really is unchanged, (ii) the
window really is denser, and (iii) HOW MUCH coarser the outside got -- the residual confound, as a
measured bound in the phase-two sense ("a bound is not cleanliness"), not as a claim of purity.

Outputs (TRACKED): bench/gate_results/task3_local_window_guard.csv
"""

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))

from pyfp3d.meshgen.planar import naca0012_coordinates            # noqa: E402
from pyfp3d.meshgen.structured import airfoil_surface_distribution  # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_local_window_guard.csv")
WINDOW, FACTOR, N_STATIONS = (0.45, 0.75), 3.0, 160


def spacings(coords, window):
    """median spacing inside / outside the chordwise window, and the counts."""
    x = coords[:, 0]
    d = np.linalg.norm(np.diff(coords, axis=0), axis=1)
    xm = 0.5 * (x[:-1] + x[1:])
    m = (xm >= window[0]) & (xm <= window[1])
    return (float(np.median(d[m])), int(m.sum()),
            float(np.median(d[~m])), int((~m).sum()))


def main():
    ref = naca0012_coordinates(n_half=401)
    base = airfoil_surface_distribution(ref, N_STATIONS, le_cluster=0.8)
    loc = airfoil_surface_distribution(ref, N_STATIONS, le_cluster=0.8,
                                       local_window=WINDOW, local_factor=FACTOR)
    din_b, nin_b, dout_b, nout_b = spacings(base, WINDOW)
    din_l, nin_l, dout_l, nout_l = spacings(loc, WINDOW)
    r_in, r_out = din_l / din_b, dout_l / dout_b

    count_ok = len(loc) == len(base) == N_STATIONS
    dense_ok = r_in < 1.0
    #: ★ vacuity guard, the phase-two lesson: a ratio that cannot move is not a PASS
    moved_ok = abs(r_in - 1.0) > 0.05

    print(f"local chordwise window {WINDOW} factor {FACTOR}   n_stations {N_STATIONS}\n")
    print(f"  stations      base {len(base)}  local {len(loc)}   "
          f"-> {'✓ unchanged' if count_ok else '★ CHANGED'}")
    print(f"  in-window     spacing {din_b:.6f} -> {din_l:.6f}   ratio {r_in:.3f}  "
          f"(points {nin_b} -> {nin_l})")
    print(f"  out-of-window spacing {dout_b:.6f} -> {dout_l:.6f}   ratio {r_out:.3f}  "
          f"(points {nout_b} -> {nout_l})")
    print(f"\n  ✓ count fixed: {count_ok}   ✓ window denser: {dense_ok}   "
          f"✓ effect non-vacuous: {moved_ok}")
    print(f"  ★ residual confound, MEASURED not assumed: the outside coarsens by "
          f"{100*(r_out-1):.1f} %, against the window's {100*(1-r_in):.1f} % refinement "
          f"= a {(1-r_in)/(r_out-1):.1f}x ratio.")
    print("    That is a BOUND on the confound, not cleanliness. If a later reading lands")
    print("    marginal, it must be resolved with a dedicated control leg, not argued away.")

    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    row = dict(window=str(WINDOW), factor=FACTOR, n_stations=N_STATIONS,
               n_base=len(base), n_local=len(loc), count_fixed=count_ok,
               d_in_base=round(din_b, 8), d_in_local=round(din_l, 8), ratio_in=round(r_in, 6),
               d_out_base=round(dout_b, 8), d_out_local=round(dout_l, 8),
               ratio_out=round(r_out, 6), n_in_base=nin_b, n_in_local=nin_l,
               confound_ratio=round((1 - r_in) / (r_out - 1), 3),
               window_denser=dense_ok, non_vacuous=moved_ok,
               verdict="PASS" if (count_ok and dense_ok and moved_ok) else "FAIL")
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row)); w.writeheader(); w.writerow(row)
    print(f"\nwrote {CSV}   verdict {row['verdict']}")
    return 0 if row["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
