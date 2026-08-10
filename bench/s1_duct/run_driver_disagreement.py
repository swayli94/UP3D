"""GS1b.7 F1 + F2: why do the two conforming drivers put the shock 0.054c apart?

Same mesh, same upwind_c, same m_crit, same condition, no entropy correction:
solve_transonic_lifting (Picard continuation, what P4/G4.1 locks) reads x_shock
0.6041 and solve_newton_lifting (coupled Newton, what G8.1/M1 lock) reads 0.6581.
That is 1.35x the entropy effect and wider than the Euler anchor tolerance, so it
has to be attributed before any shock position is a deliverable.

F1 reads both solutions side by side INCLUDING the final residual -- if one is orders
of magnitude looser, it is a non-converged state and the gap is a convergence
artifact rather than a model difference.

F2 is the decisive probe: CROSS-SEED. Hand the Newton the Picard's phi and see where
it goes, and hand the Picard the Newton's phi likewise.

    both land on one state  =>  one solution; the other leg had not converged
    both stay put           =>  genuinely TWO solutions (GS1.1 proved this
                                discretisation admits them in the duct)

Outputs: results/gs1b_7_driver_disagreement.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                          # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                         # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve                 # noqa: E402
from pyfp3d.post.shock import shock_report                        # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients           # noqa: E402
from pyfp3d.solve.continuation import solve_transonic_lifting     # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting              # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

LEVEL, M_INF, ALPHA, C = "coarse", 0.80, 1.25, 1.5


def read_out(mc, dz, phi, gam, tag, res, nit, extra=""):
    rep = shock_report(wall_cp_curve(mc, phi, z=0.5 * dz, m_inf=M_INF), M_INF)
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], phi,
                                alpha_deg=ALPHA, s_ref=dz, m_inf=M_INF)
    row = dict(case=tag, x_shock=round(float(rep["upper"]["x_shock"]), 5),
               cl_p=round(f["cl"], 6),
               gamma=round(float(np.atleast_1d(gam)[0]), 8),
               res_final=res, n_iter=nit, note=extra)
    print(f"  {tag:34s} x_shock={row['x_shock']:.4f} cl={row['cl_p']:.6f} "
          f"gamma={row['gamma']:.6f} |R|={res:.2e} n={nit} {extra}", flush=True)
    return row


def main():
    mc, wc = cut_wake(read_mesh(REPO / f"cases/meshes/naca0012_2.5d/{LEVEL}.msh"))
    dz = float(np.ptp(mc.nodes[:, 2]))
    rows = []
    print(f"NACA0012 {LEVEL} M{M_INF} alpha{ALPHA}, C={C}, entropy OFF\n")

    # ---------------- F1: both drivers, side by side ----------------
    t0 = time.perf_counter()
    pic = solve_transonic_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                  max_gamma_evals=12, n_picard_eval=800)
    t_pic = time.perf_counter() - t0
    rows.append(read_out(mc, dz, pic["phi"], pic["gamma"],
                         "F1 Picard continuation",
                         float(pic.get("residual_history", [np.nan])[-1]),
                         pic.get("n_picard_total", pic.get("n_picard", -1)),
                         f"({t_pic:.0f}s) kutta_conv={pic.get('kutta_converged')}"))

    t0 = time.perf_counter()
    new = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, upwind_c=C,
                               m_crit=0.95, freeze_tol=1e-6,
                               freeze_refresh_max=8, precond="direct",
                               direct_refactor_every=4, n_newton_max=80)
    t_new = time.perf_counter() - t0
    rows.append(read_out(mc, dz, new["phi"], new["gamma"], "F1 coupled Newton",
                         new["residual_history"][-1], new["n_newton"],
                         f"({t_new:.0f}s) conv={new['converged']} "
                         f"clamped={new['clamped']}"))

    # ---------------- F2: cross-seed ----------------
    print("\n  F2 cross-seed (the decisive probe):", flush=True)
    t0 = time.perf_counter()
    n_from_p = solve_newton_lifting(
        mc, wc, m_inf=M_INF, alpha_deg=ALPHA, upwind_c=C, m_crit=0.95,
        freeze_tol=1e-6, freeze_refresh_max=8, precond="direct",
        direct_refactor_every=4, n_newton_max=80,
        phi_init=np.asarray(pic["phi"], dtype=np.float64),
        gamma_init=np.atleast_1d(np.asarray(pic["gamma"], dtype=np.float64)),
        n_picard_seed=0)
    rows.append(read_out(mc, dz, n_from_p["phi"], n_from_p["gamma"],
                         "F2 Newton seeded FROM Picard",
                         n_from_p["residual_history"][-1], n_from_p["n_newton"],
                         f"({time.perf_counter()-t0:.0f}s) "
                         f"conv={n_from_p['converged']}"))

    # ★ F1b, cheaper and more decisive than a re-solve: evaluate the NEWTON's own
    # residual AT the Picard state. If |R(phi_Picard)| is large, the Picard state is
    # not a solution of the same discrete equations at all, and the 0.054c is not a
    # disagreement between two solutions -- it is one leg not being converged.
    # `solve_transonic_lifting` takes no phi_init, so this replaces the second
    # cross-seed leg (recorded, not skipped silently).
    from pyfp3d.solve.newton import NewtonWorkspace
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA)
    ws.set_mach(M_INF)
    for tag, sol in (("Picard", pic), ("Newton", new)):
        pf = np.asarray(sol["phi"], dtype=np.float64)[:ws.n_red][ws.free].copy()
        g = np.atleast_1d(np.asarray(sol["gamma"], dtype=np.float64))
        R, F, _ = ws.eval_residual(pf, g, C, 0.95, 3.0, 0.05)
        rn, fn = float(np.max(np.abs(R))), float(np.max(np.abs(F)))
        rows.append(dict(case=f"F1b Newton residual AT the {tag} state",
                         x_shock=None, cl_p=None, gamma=round(float(g[0]), 8),
                         res_final=rn, n_iter=None,
                         note=f"|F_kutta|={fn:.3e}"))
        print(f"  Newton |R| at the {tag:6s} state: {rn:.3e}   "
              f"|F_kutta| = {fn:.3e}", flush=True)

    with open(OUT / "gs1b_7_driver_disagreement.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1b_7_driver_disagreement.csv")

    print("\n=== reading ===")
    p, n, nfp = rows[0], rows[1], rows[2]
    rp = next(r for r in rows if "AT the Picard" in r["case"])
    rnw = next(r for r in rows if "AT the Newton" in r["case"])
    print(f"  F1 gap: {abs(p['x_shock'] - n['x_shock']):.4f} c   "
          f"residuals {p['res_final']:.2e} (Picard) vs {n['res_final']:.2e} "
          f"(Newton)")
    print(f"  F2 Newton from Picard's phi -> {nfp['x_shock']:.4f} "
          f"(Newton's own {n['x_shock']:.4f}, Picard's {p['x_shock']:.4f})")
    print(f"  F1b Newton |R|: at the Picard state {rp['res_final']:.2e}, "
          f"at its own {rnw['res_final']:.2e}")
    stays = abs(nfp["x_shock"] - p["x_shock"]) < 0.01
    if rp["res_final"] > 1e3 * max(rnw["res_final"], 1e-14):
        print("  => the PICARD state is NOT a solution of the same discrete "
              "equations (its residual in the Newton's own norm is orders above "
              "the Newton's) -- the 0.054c is a CONVERGENCE gap, not two "
              "solutions.")
    elif stays:
        print("  => TWO SOLUTIONS: seeded from the Picard state the Newton stays "
              "there, and both states carry small residuals.")
    else:
        print("  => read the table: the Newton moved off the Picard state; check "
              "whether the Picard residual was already small.")


if __name__ == "__main__":
    main()
