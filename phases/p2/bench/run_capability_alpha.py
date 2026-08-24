"""The capability matrix's alpha axis: can lowering alpha buy back an M_max-limited level?

Pre-registered in phases/p2/docs/dev_phase_two/20260803-1100-alpha-axis-prereg.md, committed before
this file was written.

The 13-cell matrix walked M_inf only and held alpha at each case's convention (M6 3.06 deg,
NACA 1.25 deg) -- half of the relaxation the user actually granted ("lower the freestream
Mach AND the angle of attack"). Five cells were stopped by the M_max < 1.4 criterion rather
than by convergence, and one of them -- ls_wb_medium -- has NO valid point at all: the
ladder's floor M_inf = 0.50 already reads M_max = 1.5025 while converging cleanly. There is
no lower Mach rung to retreat to, so alpha is the only remaining lever for that cell.

Everything is reused from run_capability_matrix rather than copied -- the same recipes, the
same classify(), the same MMAX_LIMIT, the same measure()/save_cp() -- so these rows are
directly comparable to the matrix rows. Two module globals are rebound below (the output
CSV path, and nothing else), which is the whole reason append_row/save_cp read them from
module scope.

Pre-registered readings: P1 the alpha sensitivity S = -dM_max/dalpha should be much larger
on the wing-body cells than on the 2.5-D NACA cell (if it is not, the "singularity vs
genuine transonic pocket" distinction between the two kinds of CLEAN_OVER_MMAX collapses
and must be dropped); P2 ls_wb_medium at alpha = 0 should be clean below M_max 1.0, which
is B23's own claim about the junction pocket being lift-coupled -- if it is not, B23's
attribution does NOT extrapolate here and that is a new finding, not a number to explain
away.

Outputs (TRACKED): bench/gate_results/capability_alpha.csv
                   bench/gate_results/capability/<cell>_a<alpha>_M<M>_cp.csv
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

MATRIX_CSV = os.path.join(_GATE, "capability_matrix.csv")
#: ★ rebound so append_row (which reads cap.CSV from module scope) writes the alpha
#: rows to their own file instead of appending into the matrix.
cap.CSV = os.path.join(_GATE, "capability_alpha.csv")

#: descending; the sweep stops at the FIRST clean rung (pre-registration Sec 4).
#: The final 0.0 is the mechanism check for P2, not a deliverable condition -- lift is
#: ~zero there, so what is read is M_max, not cl.
ALPHA_LADDER = {cap.ALPHA_M6: (2.00, 1.00, 0.50, 0.00),
                cap.ALPHA_NACA: (0.80, 0.50, 0.25, 0.00)}

#: cells whose alpha = 0 reading is required by pre-registration P2 regardless of where
#: the descending ladder stops. ls_wb_medium is the cell with zero valid points in the
#: matrix, and B23's lift-coupling claim is exactly a claim about its alpha = 0 state.
P2_CELLS = {"ls_wb_medium"}


def first_over_mmax():
    """Per cell: the Mach of the first CLEAN_OVER_MMAX row in the matrix, i.e. the level
    alpha is being asked to buy back. Data-driven off the committed CSV so this cannot
    drift out of step with the matrix."""
    if not os.path.exists(MATRIX_CSV):
        raise SystemExit(f"matrix CSV missing: {MATRIX_CSV}")
    rows = list(csv.DictReader(open(MATRIX_CSV)))
    out = {}
    for r in rows:
        if r["status"] == "CLEAN_OVER_MMAX" and r["cell"] not in out:
            out[r["cell"]] = (float(r["m_inf"]), float(r["m_max"]))
    return out


def done():
    """(cell, alpha) pairs already measured, so a relaunch costs nothing."""
    if not os.path.exists(cap.CSV):
        return set()
    return {(r["cell"], round(float(r["alpha"]), 4))
            for r in csv.DictReader(open(cap.CSV))}


def main():
    targets = first_over_mmax()
    if not targets:
        print("no CLEAN_OVER_MMAX cells in the matrix -- nothing for alpha to buy back")
        return 0
    meta = {c[0]: c for c in cap.CELLS}
    already = done()
    if already:
        print(f"resuming -- {len(already)} (cell, alpha) points already measured",
              flush=True)
    only = [c for c in os.environ.get("PYFP3D_CAP_CELLS", "").split(",") if c]

    summary = {}
    #: cheapest first: this sweep is minutes, but an interruption should still leave
    #: whole cells finished rather than all of them half-done.
    order = sorted(targets, key=lambda c: cap._LEVEL_ORDER[meta[c][4]])
    for cell in order:
        if only and cell not in only:
            continue
        _, path, geom, mdir, level, alpha0, fn = meta[cell]
        m, m_max0 = targets[cell]
        print(f"\n=== {cell}  ({path}, {geom}, {level}) — hold M{m}, "
              f"was M_max {m_max0:.4f} at alpha {alpha0} ===", flush=True)
        hit = None
        for a in ALPHA_LADDER[alpha0]:
            if (cell, round(a, 4)) in already and not only:
                continue
            row, payload = cap.measure(cell, path, geom, mdir, level, a, fn, m)
            cap.append_row(row)
            st = row["status"]
            print(f"  alpha {a:<5} {st:16s} conv={row.get('converged')} "
                  f"lim/flr={row.get('n_limited')}/{row.get('n_floored')} "
                  f"M_max={row.get('m_max')} cl_p={row.get('cl_p')} "
                  f"|R|={row.get('res_final')} ({row.get('wall_s')}s)"
                  f"{'  ' + row['note'] if row.get('note') else ''}", flush=True)
            if st in ("MESH_MISSING", "ERROR"):
                break
            if payload is not None and st.startswith("CLEAN"):
                cap.save_cp(f"{cell}_a{str(a).replace('.', '')}", m, geom, payload)
            if st == "CLEAN":
                hit = a
                print(f"  -> alpha-relaxed envelope: M{m} at alpha {a}", flush=True)
                break
        #: ★ P2 needs the alpha = 0 reading SPECIFICALLY, and the loop above breaks at
        #: the first clean rung -- so on a cell that goes clean early, alpha = 0 would
        #: never be measured and the pre-registered mechanism check would silently not
        #: happen. Run it explicitly. (My first draft "resolved" this in a comment
        #: claiming the P2 cell has no clean rung to break on; that is an assumption
        #: about the result, made before measuring it.)
        if cell in P2_CELLS and (cell, 0.0) not in done():
            row, payload = cap.measure(cell, path, geom, mdir, level, 0.0, fn, m)
            cap.append_row(row)
            print(f"  alpha 0.0   {row['status']:16s} [P2 mechanism check] "
                  f"conv={row.get('converged')} M_max={row.get('m_max')} "
                  f"({row.get('wall_s')}s)", flush=True)
            if payload is not None and row["status"].startswith("CLEAN"):
                cap.save_cp(f"{cell}_a00", m, geom, payload)
        summary[cell] = (m, hit)

    print("\n=== ALPHA-RELAXED ENVELOPE (M_max < 1.4 kept, alpha lowered) ===")
    for cell, (m, hit) in summary.items():
        if hit is None:
            print(f"  {cell:18s} M{m}: NOT recovered at any tested alpha")
        else:
            print(f"  {cell:18s} M{m}: clean at alpha {hit} "
                  f"(convention {meta[cell][5]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
