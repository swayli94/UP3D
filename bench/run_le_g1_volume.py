"""G1's missing half: did the `h_far` arm change the NEAR-BODY VOLUME, not just the far cells?

The third factorial concluded (criterion P) that the LE deficit is dominated by the
bulk/far mesh, on the strength of a far-field-only arm moving the LE band by 0.098.
Its guard G1 checked the SURFACE mesh -- wall triangle count and LE-band tangential
spacing -- and found only a 0.343 % change. But G1 never checked the near-body VOLUME.

That matters for what the conclusion may SAY. If h_far 2.4 -> 1.8 also re-graded the
cells near the wing, then the "far-only" arm was not far-only inside the volume, and P's
wording has to become "the BULK MESH (including near-body grading) controls the LE band"
rather than "the FAR FIELD controls it" -- two statements pointing at different next
steps. The conclusion survives either way; the attribution does not.

Prior that makes the grading explanation the more likely one: B5 (2026-07-12) swept the
far-field RADIUS over R in {15, 30, 60, 120} c and measured option a (Dirichlet+vortex,
which is what these legs use) to be DOMAIN-ROBUST -- Gamma within 0.45 % / 1.09 % of the
truth across the whole range. If the domain and its BC are that insensitive, a far-field
SPACING change moving the LE band by 0.098 is more easily explained by mesh grading than
by the boundary condition.

No solves: mesh statistics only.

Outputs (TRACKED): bench/gate_results/le_g1_volume.csv
"""

import csv
import os
import sys

os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from run_le_factorial import le_face_count                          # noqa: E402
from run_le_factorial2 import _build                                # noqa: E402
from run_le_response import le_geometry                             # noqa: E402

CSV = os.path.join(HERE, "gate_results", "le_g1_volume.csv")
#: the third factorial's four legs verbatim
LEGS = (("G00_base", 0.010, 2.4), ("G10_faronly", 0.010, 1.8),
        ("G01_leonly", 0.0075, 2.4), ("G11_both", 0.0075, 1.8))
#: radial shells (in MAC) over which the volume grading is compared. The wing sits near
#: the origin, so these walk from "hugging the body" out toward the far field, and a
#: grading change shows up as a shifted cell-size profile rather than a single number.
SHELLS = ((0.0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0), (8.0, 1e9))


def cell_sizes(mc):
    """Per-tet size (cube root of volume) and centroid radius."""
    p = mc.nodes[mc.elements]
    v = np.abs(np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                         np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))) / 6.0
    h = np.cbrt(np.maximum(v, 1e-300))
    c = p.mean(axis=1)
    return h, np.linalg.norm(c, axis=1)


def main():
    print("G1's missing half: what did the h_far arm change INSIDE the volume?\n")
    rows, ref = [], None
    for tag, h_le, h_far in LEGS:
        mc, _wc = cut_wake(_build(h_le, h_far))
        ht, hn, aniso = le_geometry(mc)
        h, r = cell_sizes(mc)
        row = dict(tag=tag, h_le=h_le, h_far=h_far, n_tet=len(mc.elements),
                   n_wall=len(mc.boundary_faces["wall"]),
                   n_le_faces=le_face_count(mc),
                   le_ht=round(ht, 9), le_hn=round(hn, 9),
                   le_aniso=round(aniso, 5))
        for lo, hi in SHELLS:
            m = (r >= lo) & (r < hi)
            k = f"h_med_r{lo:g}_{'inf' if hi > 1e8 else f'{hi:g}'}"
            row[k] = round(float(np.median(h[m])), 8) if m.any() else None
            row[k.replace("h_med", "n_cells")] = int(m.sum())
        rows.append(row)
        if tag == "G00_base":
            ref = row
        print(f"  {tag:14} tets {row['n_tet']:>7}  LE h_t {ht:.6f}  LE h_n {hn:.6f}  "
              f"aniso {aniso:.3f}")
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}),
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")

    print("\n=== the far-only arm (G10) against the base, per radial shell ===")
    g10 = next(r for r in rows if r["tag"] == "G10_faronly")
    print(f"  {'shell (MAC)':>14}{'base h_med':>13}{'far-only':>12}{'change':>10}"
          f"{'base cells':>12}{'far cells':>11}")
    worst_near = 0.0
    for lo, hi in SHELLS:
        k = f"h_med_r{lo:g}_{'inf' if hi > 1e8 else f'{hi:g}'}"
        a, b = ref[k], g10[k]
        if a is None or b is None:
            continue
        d = 100.0 * (b - a) / a
        if hi <= 2.0:
            worst_near = max(worst_near, abs(d))
        print(f"  {f'{lo:g}-{hi:g}' if hi < 1e8 else f'{lo:g}+':>14}{a:>13.6f}"
              f"{b:>12.6f}{d:>9.2f}%{ref[k.replace('h_med','n_cells')]:>12}"
              f"{g10[k.replace('h_med','n_cells')]:>11}")
    d_hn = 100.0 * (g10["le_hn"] - ref["le_hn"]) / ref["le_hn"]
    d_ht = 100.0 * (g10["le_ht"] - ref["le_ht"]) / ref["le_ht"]
    print(f"\n  LE band: h_t {d_ht:+.3f} %   h_n {d_hn:+.3f} %   "
          f"(G1 checked h_t only)")
    print(f"  worst near-body (r < 2 MAC) median-size change: {worst_near:.2f} %")
    print("\n=== reading (threshold fixed before running: 1 % on the near-body shells) ===")
    if worst_near < 1.0 and abs(d_hn) < 1.0:
        print("  NEAR-BODY UNCHANGED -> P's wording stands: the FAR FIELD controls the")
        print("  LE band. The r_far branch is then the one to test, against B5's prior")
        print("  that option a is domain-robust (so the prediction is 'no change').")
    else:
        print("  ★ NEAR-BODY MOVED -> P's conclusion stands but its WORDING must change:")
        print("  the arm was not far-only inside the volume, so what is attributed is the")
        print("  BULK MESH INCLUDING NEAR-BODY GRADING, not the far field as such. An")
        print("  erratum goes on the third factorial's verdict.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
