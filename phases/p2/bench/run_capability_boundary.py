"""Locate each cell's capability boundary to +-0.01 in Mach instead of +-0.03.

The matrix walked a nine-rung ladder whose spacing near the interesting end is 0.02-0.03,
so every envelope it reported is located only to that precision -- a limitation the
matrix pre-registration recorded rather than hid. This round pays it down. No new meshes
are needed: refining the MACH ladder only needs more Mach points on the meshes that
already exist. (I had wrongly said this was blocked on the broken gmsh; it never was.
gmsh gates new mesh LEVELS, not new Mach points.)

Protocol, per cell, from the committed matrix rows:
  - a cell with a CLEAN rung and a stop rung: two interior points at 1/3 and 2/3 of the
    bracket. Two points rather than a bisection ON PURPOSE -- a bisection would assume the
    CLEAN/not-CLEAN transition is monotone in Mach, and nothing here establishes that.
    Two independent interior readings can show non-monotonicity if it exists.
  - ls_wb_medium, which has NO valid rung: M0.40 and M0.45, below the ladder's floor.
    This asks whether the cell has any valid Mach at all at its conventional alpha.
  - the two cells that exhausted the ladder (conf/ls wing coarse, both reaching M0.84):
    M0.86 and M0.88, because "ladder exhausted" is not a located boundary.

Recipes, classify(), MMAX_LIMIT and measure() are all imported from the matrix runner, so
these rows are directly comparable to the matrix rows and to each other.

Cheapest cell first: this is hours long, and an interruption should leave whole cells
finished. conf_naca_fine is last on its own (~26 min per point).

Outputs (TRACKED): bench/gate_results/capability_boundary.csv
"""

import csv
import os
import sys

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

HERE = os.path.dirname(os.path.abspath(__file__))
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import run_capability_matrix as cap                                 # noqa: E402

MATRIX = os.path.join(_GATE, "capability_matrix.csv")
cap.CSV = os.path.join(_GATE, "capability_boundary.csv")

BELOW_FLOOR = (0.45, 0.40)          # for a cell with no valid rung
ABOVE_LADDER = (0.86, 0.88)         # for a cell that exhausted the ladder


def plan():
    rows = list(csv.DictReader(open(MATRIX)))
    out = {}
    for cell in dict.fromkeys(r["cell"] for r in rows):
        rr = [r for r in rows if r["cell"] == cell]
        clean = [float(r["m_inf"]) for r in rr if r["status"] == "CLEAN"]
        stop = [r for r in rr if r["status"] != "CLEAN"]
        if not clean:
            out[cell] = list(BELOW_FLOOR)
        elif not stop:
            out[cell] = list(ABOVE_LADDER)
        else:
            lo, hi = max(clean), float(stop[0]["m_inf"])
            out[cell] = [round(lo + (hi - lo) * f, 3) for f in (0.34, 0.67)]
    return out


def done():
    if not os.path.exists(cap.CSV):
        return set()
    return {(r["cell"], round(float(r["m_inf"]), 4))
            for r in csv.DictReader(open(cap.CSV))}


def main():
    todo = plan()
    meta = {c[0]: c for c in cap.CELLS}
    already = done()
    only = [c for c in os.environ.get("PYFP3D_CAP_CELLS", "").split(",") if c]
    #: cheapest first, by the matrix's own measured wall time for the cell's last point
    rows = list(csv.DictReader(open(MATRIX)))
    cost = {}
    for cell in todo:
        w = [float(r["wall_s"]) for r in rows if r["cell"] == cell]
        cost[cell] = max(w) if w else 0.0
    order = sorted(todo, key=lambda c: cost[c])

    summary = {}
    for cell in order:
        if only and cell not in only:
            continue
        _, path, geom, mdir, level, alpha, fn = meta[cell]
        print(f"\n=== {cell} ({level}) — refine at {todo[cell]} "
              f"(cell cost ~{cost[cell]:.0f}s/pt) ===", flush=True)
        got = []
        for m in todo[cell]:
            if (cell, round(m, 4)) in already and not only:
                continue
            row, payload = cap.measure(cell, path, geom, mdir, level, alpha, fn, m)
            cap.append_row(row)
            print(f"  M{m:<6} {row['status']:18s} conv={row.get('converged')} "
                  f"m_att={row.get('m_attained')} M_max={row.get('m_max')} "
                  f"cl_p={row.get('cl_p')} |R|={row.get('res_final')} "
                  f"d10={row.get('descent10')} ({row.get('wall_s')}s)", flush=True)
            if payload is not None and row["status"].startswith("CLEAN"):
                cap.save_cp(cell, m, geom, payload)
            got.append((m, row["status"]))
        summary[cell] = got

    print("\n=== BOUNDARY REFINEMENT ===")
    for cell, got in summary.items():
        cl = [m for m, s in got if s == "CLEAN"]
        print(f"  {cell:18s} " + "  ".join(f"M{m}:{s}" for m, s in got)
              + (f"   -> envelope now >= M{max(cl)}" if cl else "   -> no new clean point"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
