"""Product metric M1, evaluated honestly (phase two GS1.5).

M1 (docs/dev_phase_two/roadmap.md 2): NACA0012 M0.80 / alpha 1.25 --
  (a) upper-surface shock at 0.61 +- 0.02 (the project's own Euler-anchored
      reference band, cases/reference_data/naca0012_m080/),
  (b) cl agreeing within 3 % between two mesh levels,
  (c) cl varying by less than 3 % over upwind_c in [1, 3].

Solved at the TARGET Mach directly (no Mach ramp -- GS1.2b measured the ramp to
reproduce the same answers where it works, to be 2.8x slower, and to fail
outright at the fine level).

This script exists so the FAIL is an artifact rather than a sentence: exit code
0 only if all three criteria hold. Decision D3 keeps M1 in the metric table with
its measured failure while M1a locks the in-envelope capability.

    python bench/run_m1_gate.py        # ~6 min
"""
import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                      # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve              # noqa: E402
from pyfp3d.post.shock import shock_report                     # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients        # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting           # noqa: E402

M_INF, ALPHA = 0.80, 1.25
SHOCK_REF, SHOCK_TOL = 0.61, 0.02
CS = (1.0, 1.5, 3.0)
LEVELS = ("coarse", "medium")


def main():
    # GS1b.3: --entropy re-runs the SAME gate with the entropy-corrected
    # density. The three criteria, the condition, the C sweep and the recipe are
    # untouched -- only the density law changes, so the FAIL recorded on
    # 2026-07-29 stays directly comparable.
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--entropy", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    print(f"entropy_correction = {a.entropy}")
    rows = []
    for level in LEVELS:
        path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
        if not path.exists():
            print(f"skip {level}: mesh missing")
            continue
        mc, wc = cut_wake(read_mesh(path))
        dz = float(np.ptp(mc.nodes[:, 2]))
        for C in CS:
            t0 = time.perf_counter()
            try:
                r = solve_newton_lifting(
                    mc, wc, m_inf=M_INF, alpha_deg=ALPHA, upwind_c=C,
                    m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
                    precond="direct", direct_refactor_every=4,
                    n_newton_max=80, entropy_correction=a.entropy)
                err = ""
            except Exception as exc:                           # noqa: BLE001
                rows.append(dict(level=level, C=C, converged=False,
                                 error=type(exc).__name__))
                print(f"  {level} C={C}: FAILED {type(exc).__name__}",
                      flush=True)
                continue
            rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                            m_inf=M_INF), M_INF)
            f = wall_force_coefficients(mc.nodes, mc.elements,
                                        mc.boundary_faces["wall"], r["phi"],
                                        alpha_deg=ALPHA, s_ref=dz,
                                        m_inf=M_INF)
            rows.append(dict(level=level, C=C, converged=r["converged"],
                             clamped=r.get("clamped"),
                             cl_p=round(f["cl"], 6),
                             x_shock=rep["upper"].get("x_shock"),
                             m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                             res_final=r["residual_history"][-1],
                             wall_s=round(time.perf_counter() - t0, 1),
                             error=err))
            print(f"  {level:7s} C={C:<4} conv={str(r['converged']):5s} "
                  f"cl={rows[-1]['cl_p']} x_shock={rows[-1]['x_shock']} "
                  f"M_max={rows[-1]['m_max']}", flush=True)

    # ★ `bench/results/` is gitignored (it holds the big bitcheck npz dumps), so
    # gate evidence goes to `bench/gate_results/`, which is TRACKED. Found
    # 2026-07-29: the GS1.5 close-out round file claimed the M1 FAIL artifact was
    # committed and it never was -- discipline #3 says a number living only in a
    # .md is not evidence, and that was exactly the situation.
    out_dir = HERE / "gate_results"
    out_dir.mkdir(exist_ok=True)
    out_csv = out_dir / ("m1_gate_entropy.csv" if a.entropy else "m1_gate.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)

    print("\n=== M1 criteria ===")
    ok = True
    # (a) shock band, default C, both levels
    for level in LEVELS:
        m = [r for r in rows if r.get("level") == level and r.get("C") == 1.5]
        if not m or m[0].get("x_shock") is None:
            print(f"  (a) {level}: no result -> FAIL")
            ok = False
            continue
        x, conv = m[0]["x_shock"], m[0]["converged"]
        good = conv and abs(x - SHOCK_REF) <= SHOCK_TOL
        print(f"  (a) {level:7s} x_shock {x:.4f} vs {SHOCK_REF}+-{SHOCK_TOL}"
              f"  conv={conv}  -> {'PASS' if good else 'FAIL'}")
        ok &= bool(good)
    # (b) two-level cl agreement at the default C
    cl = {}
    for level in LEVELS:
        m = [r for r in rows if r.get("level") == level and r.get("C") == 1.5
             and r.get("cl_p") is not None]
        if m:
            cl[level] = (m[0]["cl_p"], m[0]["converged"])
    if len(cl) == 2:
        d = (cl["medium"][0] - cl["coarse"][0]) / abs(cl["coarse"][0])
        good = abs(d) < 0.03 and all(c[1] for c in cl.values())
        print(f"  (b) cl coarse {cl['coarse'][0]:.4f} -> medium "
              f"{cl['medium'][0]:.4f} = {100 * d:+.1f} % (< 3 %)"
              f"  converged={[c[1] for c in cl.values()]}"
              f"  -> {'PASS' if good else 'FAIL'}")
        ok &= bool(good)
    else:
        print("  (b) FAIL: fewer than two levels produced a result")
        ok = False
    # (c) dissipation sensitivity per level
    for level in LEVELS:
        sub = [r for r in rows if r.get("level") == level
               and r.get("cl_p") is not None and r.get("converged")]
        if len(sub) < 2:
            print(f"  (c) {level}: {len(sub)} converged leg(s) -> FAIL")
            ok = False
            continue
        v = [r["cl_p"] for r in sub]
        spread = (max(v) - min(v)) / min(v)
        good = spread < 0.03
        print(f"  (c) {level:7s} cl over C in [1,3]: {min(v):.4f}..{max(v):.4f}"
              f" = {100 * spread:.1f} % (< 3 %) -> {'PASS' if good else 'FAIL'}")
        ok &= bool(good)

    print(f"\nM1: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
