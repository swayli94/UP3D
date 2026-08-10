"""GS1b.2 Q3: convert the dissipation lever (Q2) into a DENSITY lever, so it can
be compared with the entropy lever measured in Q1.

Q2 measured M_inf_max(C) at fixed h: raising the dissipation constant pushes the
fold up. But `upwind_c` is not a physical quantity -- to predict what an entropy
correction (which changes the density by a KNOWN factor, Q1) would buy, the C
axis has to be re-expressed in density units:

    for each C, at a COMMON freestream Mach on all branches, measure how much
    density the upwinding actually removes:  d(C) = max |rho_tilde - rho| / rho
    over the elements (and the mean over the supersonic zone),

then the transmission ratio is

    T = dM_inf_max / dd   [Mach per unit relative density deficit]

and Q1's entropy deficit (5.1 % at M1 = 1.30, 10.1 % at M1 = 1.40) times T is
the predicted ΔM_inf_max of the entropy correction.

Honest caveats, recorded before the numbers (they bound what this can claim):
  1. The artificial density acts at the UPSTREAM donor element (a retarded
     density over the whole supersonic zone); the entropy factor acts on the
     DOWNSTREAM side of the shock only. Same units, different support -- so T is
     a SCALE estimate, not a prediction.
  2. Artificial density is dissipative (it also smears the shock); the entropy
     correction is not. If anything that makes T an OPTIMISTIC bound for how
     much fold travel a given density change buys, since part of C's effect is
     smearing rather than the density level itself.
Both are quoted in the round file's verdict, not buried here.

Outputs: results/gs1b_2_q3_transmission.csv
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pyfp3d.physics.isentropic import density_field              # noqa: E402
from pyfp3d.physics.isentropic import limit_q2_field            # noqa: E402
from pyfp3d.kernels.upwind import UpwindOperator                   # noqa: E402
from pyfp3d.mesh.reader import read_mesh                           # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                          # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve                  # noqa: E402
from pyfp3d.post.shock import shock_report                         # noqa: E402
from pyfp3d.kernels.jacobian import PicardOperator               # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting               # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

LEVEL = "medium"
ALPHA = 1.25
#: the common Mach: inside every branch measured in Q2 (the smallest
#: M_inf_max was 0.7925 at C = 1.5), so all four C legs are read on a
#: converged, unclamped state at the SAME freestream condition.
M_COMMON = 0.7875
C_LIST = (1.5, 2.0, 3.0, 4.0)
#: Q2's measured branch ends (>= means the walk hard-stopped at 0.86 still alive)
M_MAX_Q2 = {1.5: 0.7925, 2.0: 0.8000, 3.0: 0.8600, 4.0: 0.8600}
M_MAX_IS_BOUND = {1.5: False, 2.0: False, 3.0: True, 4.0: True}


def density_deficit(mc, phi, m_inf, C):
    """(rho_tilde - rho)/rho statistics on a solved state."""
    op = PicardOperator(mc.nodes, mc.elements)
    upw = UpwindOperator(mc.nodes, mc.elements, weighted=False)
    u_inf = 1.0
    grad, q2 = op.velocities(phi)
    grad = grad.copy()
    q2n = q2 / u_inf ** 2
    q2l = limit_q2_field(q2n, m_inf, 1.02, 1.4)
    rho = density_field(q2l, m_inf, 1.4)
    rt = upw.rho_tilde(grad, q2l, rho, m_inf, C, 0.95, 1.4, 0.05).copy()
    rel = (rt - rho) / rho
    suplike = upw.nu > 0.0
    return dict(
        d_max=float(-rel.min()),                     # deepest deficit
        d_mean_sup=float(-rel[suplike].mean()) if suplike.any() else 0.0,
        n_sup=int(suplike.sum()), nu_max=float(upw.nu.max()))


def main():
    path = REPO / f"cases/meshes/naca0012_2.5d/{LEVEL}.msh"
    mc, wc = cut_wake(read_mesh(path))
    dz = float(np.ptp(mc.nodes[:, 2]))
    rows = []
    print(f"{LEVEL}: {len(mc.nodes)} nodes, common M_inf = {M_COMMON}\n")
    for C in C_LIST:
        r = solve_newton_lifting(mc, wc, m_inf=M_COMMON, alpha_deg=ALPHA,
                                 upwind_c=C, m_crit=0.95, freeze_tol=1e-6,
                                 freeze_refresh_max=8, precond="direct",
                                 direct_refactor_every=4, n_newton_max=80)
        ok = bool(r["converged"]) and not r.get("clamped", False)
        st = density_deficit(mc, r["phi"], M_COMMON, C)
        rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                         m_inf=M_COMMON), M_COMMON)
        rows.append(dict(
            level=LEVEL, m_inf=M_COMMON, upwind_c=C, usable=ok,
            gamma=round(float(r["gamma"][0]), 8),
            m_max=round(float(np.sqrt(r["mach2_max"])), 5),
            x_shock=rep["upper"].get("x_shock"),
            d_max_pct=round(100.0 * st["d_max"], 4),
            d_mean_sup_pct=round(100.0 * st["d_mean_sup"], 4),
            n_sup=st["n_sup"], nu_max=round(st["nu_max"], 5),
            m_inf_max_q2=M_MAX_Q2[C], q2_is_lower_bound=M_MAX_IS_BOUND[C]))
        print(f"  C={C:<4} usable={ok} M_max={rows[-1]['m_max']:.4f} "
              f"x_sh={rows[-1]['x_shock']} "
              f"deficit max={rows[-1]['d_max_pct']:.3f}% "
              f"mean_sup={rows[-1]['d_mean_sup_pct']:.3f}% "
              f"(n_sup={st['n_sup']})  M_inf_max={M_MAX_Q2[C]}"
              f"{'+' if M_MAX_IS_BOUND[C] else ''}", flush=True)

    with open(OUT / "gs1b_2_q3_transmission.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1b_2_q3_transmission.csv")

    print("\n=== transmission ratio T = dM_inf_max / d(deficit) ===")
    base = rows[0]
    for r in rows[1:]:
        dd = (r["d_mean_sup_pct"] - base["d_mean_sup_pct"]) / 100.0
        dm = r["m_inf_max_q2"] - base["m_inf_max_q2"]
        tag = " (LOWER BOUND: branch did not terminate)" \
            if r["q2_is_lower_bound"] else ""
        if abs(dd) > 1e-12:
            print(f"  C {base['upwind_c']} -> {r['upwind_c']}: "
                  f"d(deficit_mean_sup) = {100*dd:+.3f} pp, "
                  f"dM_inf_max = {dm:+.4f}  =>  T = {dm/dd:8.3f} "
                  f"Mach per unit deficit{tag}")
        else:
            print(f"  C {base['upwind_c']} -> {r['upwind_c']}: "
                  f"deficit unchanged, cannot form T")

    print("\n=== Q3: predicted entropy-correction travel ===")
    print("  Q1 measured entropy density deficit: 5.08 % (M1=1.30) .. "
          "10.13 % (M1=1.40)")
    for r in rows[1:]:
        dd = (r["d_mean_sup_pct"] - base["d_mean_sup_pct"]) / 100.0
        dm = r["m_inf_max_q2"] - base["m_inf_max_q2"]
        if abs(dd) <= 1e-12:
            continue
        T = dm / dd
        print(f"  via C {base['upwind_c']}->{r['upwind_c']} (T={T:.3f}): "
              f"dM_inf_max ~ {T*0.0508:+.4f} (5.08 %) .. "
              f"{T*0.1013:+.4f} (10.13 %)"
              f"{'  [T is a lower bound]' if r['q2_is_lower_bound'] else ''}")


if __name__ == "__main__":
    main()
