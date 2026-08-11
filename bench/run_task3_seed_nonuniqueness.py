"""The decisive control for task 3 item 1: does the SEED, at a FIXED mesh, move the answer?

★★ This exists because my own registered criterion C2 could not see the thing that mattered.
C2 asked for "the same mode on two independent runs" -- but two runs at a fixed seed and a fixed
thread count are bit-identical BY CONSTRUCTION, so C2 certifies reproducibility and says nothing
about solution UNIQUENESS. Phase two had already registered the soft-shift defect (a converged,
0-clamp solve landing on a DIFFERENT solution, with nothing firing on it); this script asks whether
that defect dominates the refinement comparison the four legs were built to make.

It also scopes the effect, because a defect without a scope is not a capability statement: the same
probe runs in-envelope (M0.72) and at alpha = 0, i.e. with the shock but without the lift coupling.

Outputs (TRACKED): bench/gate_results/task3_seed_nonuniqueness.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
sys.path.insert(0, HERE)

from pyfp3d.post.section_cut import wall_cp_curve                 # noqa: E402
from pyfp3d.post.shock import shock_report                        # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients           # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting              # noqa: E402
from run_le14_common_root import classify_failure                # noqa: E402
from run_task3_refinement_paradox import SOLVE_KW, build          # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_seed_nonuniqueness.csv")
SEEDS = (0, 5, 12)
#: (tag, m_inf, alpha, mesh) -- the mesh is IDENTICAL across seeds within a block. A tuple is
#: hybrid-generator arguments; a string is a committed .msh path.
#: ★★ The last two blocks are the FAMILY control. Without them the reading would be scoped to the
#: mesh family this round happens to have built, while the claim it bears on (the phase-two M1
#: readings, capability boundary section 3.6) was measured on the UNSTRUCTURED family. Annotating
#: that document from a hybrid-only measurement would be a mis-scoped guard -- the exact phase-two
#: failure of checking the skin and concluding about the volume.
BLOCKS = (("M080_a125_R0", 0.80, 1.25, (160, 0.004, None, 1.0)),
          ("M080_a125_RG", 0.80, 1.25, (320, 0.002, None, 1.0)),
          ("M072_a125_R0", 0.72, 1.25, (160, 0.004, None, 1.0)),
          ("M080_a000_R0", 0.80, 0.00, (160, 0.004, None, 1.0)),
          ("M080_a125_unstr_coarse", 0.80, 1.25, "cases/meshes/naca0012_2.5d/coarse.msh"),
          ("M080_a125_unstr_medium", 0.80, 1.25, "cases/meshes/naca0012_2.5d/medium.msh"))


def main():
    rows = []
    print(f"{'block':16}{'seed':>5}{'conv':>6}{'|R|':>10}{'lim/flr':>10}{'M_max':>8}"
          f"{'x_shock':>9}{'cl':>11}{'s':>5}")
    for tag, m_inf, alpha, args in BLOCKS:
        if isinstance(args, str):
            from pyfp3d.mesh.reader import read_mesh
            from pyfp3d.mesh.wake_cut import cut_wake
            mc, wc = cut_wake(read_mesh(os.path.join(os.path.dirname(HERE), args)))
        else:
            (mc, wc), _ = build(*args)
        z = float(np.unique(mc.nodes[:, 2]).mean())
        s_ref = float(np.ptp(mc.nodes[:, 2]))
        for seed in SEEDS:
            kw = dict(SOLVE_KW); kw["n_picard_seed"] = seed
            t = time.perf_counter()
            r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=alpha, **kw)
            w = time.perf_counter() - t
            h = np.asarray(r.get("residual_history", []), dtype=float)
            rep = shock_report(wall_cp_curve(mc, r["phi"], z=z, m_inf=m_inf), m_inf)
            f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"],
                                        r["phi"], alpha_deg=alpha, s_ref=s_ref, m_inf=m_inf)
            xs = rep["upper"].get("x_shock")
            nlim, nflr = int(r.get("n_limited") or 0), int(r.get("n_floored") or 0)
            #: ★ C1 finally has non-converged legs to answer for. The four registered legs all
            #: converged, so C1 was satisfied VACUOUSLY there; here it is answered for real.
            if r["converged"]:
                mode = ""
            else:
                mode = classify_failure(
                    h, np.asarray(r.get("clamp_history", []), dtype=float),
                    np.asarray(r.get("F_history", []), dtype=float),
                    int(r.get("n_gmres_stalled") or 0), str(r.get("accept_reason")),
                    nlim, nflr)[0]
            rows.append(dict(block=tag, m_inf=m_inf, alpha=alpha, seed=seed, n_tets=len(mc.elements),
                             converged=bool(r["converged"]),
                             res_final=float(h[-1]) if len(h) else None,
                             n_limited=nlim, n_floored=nflr,
                             m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                             x_shock=None if xs is None else round(float(xs), 6),
                             cl=round(float(f["cl"]), 6), n_newton=len(h),
                             failure_mode=mode, accept_reason=str(r.get("accept_reason")),
                             threads=os.environ["NUMBA_NUM_THREADS"], wall_s=round(w, 1)))
            print(f"{tag:16}{seed:>5}{str(r['converged']):>6}{rows[-1]['res_final']:>10.2e}"
                  f"{f'{nlim}/{nflr}':>10}{rows[-1]['m_max']:>8.4f}"
                  f"{('-' if xs is None else f'{xs:.4f}'):>9}{rows[-1]['cl']:>11.6f}{w:>5.0f}  {mode or '-'}",
                  flush=True)

    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")

    print("\n=== C1: every non-converged leg carries a MODE (answered here, not vacuously) ===")
    bad = [r for r in rows if not r["converged"]]
    if not bad:
        print("  VACUOUS -- no non-converged leg")
    for r in bad:
        print(f"  {r['block']:24} seed {r['seed']:>2}  mode={r['failure_mode']:20} "
              f"clamps={r['n_limited']}/{r['n_floored']}  |R|={r['res_final']:.2e}  "
              f"accept_reason={r['accept_reason']}")
    print(f"  -> {len(bad)} non-converged legs, "
          f"{sum(1 for r in bad if r['failure_mode'])} classified")

    print("\n=== seed-induced spread at a FIXED mesh (this is the reading) ===")
    for tag, _, _, _ in BLOCKS:
        rs = [r for r in rows if r["block"] == tag]
        #: ★★ DEFECT FIXED 2026-08-12: these two lines had NO `converged` filter, so the first
        #: published run's spreads were computed over states that included CLAMPED GARBAGE (e.g. an
        #: unstructured-coarse leg at x_shock 0.8818 with 2600/580 clamps). That inflated the
        #: unstructured coarse x_shock spread 0.0122 -> 0.2745 (22x) and manufactured a 0.3060
        #: spread for unstructured medium, where NOTHING converges and the spread is UNDEFINED.
        #: The hybrid numbers were unaffected (3/3 converged there), so the round's core reading
        #: survives -- but its "family-independent, and the unstructured family is worse" extension
        #: was REFUTED by this fix. A spread over non-solutions is not a spread over solutions.
        good = [r for r in rs if r["converged"]]
        xs = [r["x_shock"] for r in good if r["x_shock"] is not None]
        cls = [r["cl"] for r in good]
        if len(cls) < 2:
            print(f"  {tag:16} converged {len(good)}/{len(rs)}   spread UNDEFINED "
                  f"(fewer than two converged seeds -- NOT 'a small spread')")
            continue
        conv = sum(r["converged"] for r in rs)
        print(f"  {tag:16} converged {conv}/{len(rs)}   x_shock spread "
              f"{(max(xs)-min(xs)) if xs else float('nan'):.4f} c   "
              f"cl spread {max(cls)-min(cls):.6f} "
              f"({100*(max(cls)-min(cls))/max(1e-12, abs(np.mean(cls))):.1f} % of mean)")
    print("\n  ★ M1's criterion (b)/(c) needs x_shock to +-0.0055 c. Compare that with the")
    print("    seed-induced spread above at M0.80/alpha1.25: the CHOICE OF SEED, at a fixed")
    print("    mesh, is the dominant term -- so a refinement comparison at that condition is")
    print("    unreadable until solution selection is pinned.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
