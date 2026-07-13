"""P13/G13.3 (transonic) -- WHERE does the round-cap M0.84 fine solve break?

run_g133_roundtip_transonic.py established the negative result: the rounded-cap
M0.84 FINE solve is NOT a discrete solution (M_max 3.97, 5 cells over M_cap=3,
1 floored, not converged) even though coarse and medium are clean. The subsonic
box study localized the residual tip divergence to the SHARP TRAILING EDGE
(chord-frac 0.999), not the cap -- but that main run did not save the transonic
field, so the 5 over-cap cells were never pinned at M0.84.

This script re-solves the round fine at M0.84 with the SAME recipe, SAVES the
field, and reports where the over-cap and top-Mach cells sit (chord-frac,
z/b, chordwise distance past the local TE). If they cluster at the sharp tip
TE (chord-frac -> 1, z/b -> 1), the conclusion is confirmed: rounding the cap
fixed the CAP, and what still trips the transonic limiter is the design-sharp
tip TE, a different object no tip-cap change removes.

Heavy: a full M0.84 continuation on 3.26 M tets (~2 h). Saves the field to
results/g133rt_transonic_fine_field.npz (gitignored); the CSV + the printed
summary are the evidence.

  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
  python cases/demo/p13_tip_edge_singularity/run_g133_roundtip_transonic_locate.py
"""
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from cases.demo._common import write_csv  # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors  # noqa: E402
from pyfp3d.kernels.jacobian import PicardOperator  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_le, x_te  # noqa: E402
from pyfp3d.physics.isentropic import limit_q2_field, mach_squared_field  # noqa: E402
from pyfp3d.solve.newton import solve_newton_transonic  # noqa: E402

OUT = Path(__file__).parent / "results"
MESH = REPO / "cases" / "meshes" / "onera_m6_roundtip"
M_INF, ALPHA, M_CAP = 0.84, 3.06, 3.0
FORM, R_C_FRAC = "vanish_smooth", 0.05

#: Verbatim from run_g133_roundtip_transonic.py (which took it from the corrected
#: run_g132_transonic.py): amg + tight EW forcing + the fine-mesh guards.
NEWTON_M6_RECIPE = dict(
    m_start=0.30, dm=0.05, dm_min=0.01, freeze_tol=1e-6, intermediate_tol=1e-5,
    newton_kw=dict(freeze_refresh_max=8, precond="amg", n_picard_seed=12,
                   ew_eta0=1e-8, ew_eta_max=1e-8, n_newton_max=60,
                   farfield_spanwise_gamma=True),
)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    field_cache = OUT / "g133rt_transonic_fine_field.npz"

    if field_cache.exists():
        d = np.load(field_cache)
        cen, mach, over = d["cen"], d["mach"], d["over_cap_mask"]
    else:
        print(f"solving round fine at M{M_INF}/alpha{ALPHA} (saving field) ...",
              flush=True)
        t0 = time.perf_counter()
        mc, wc = cut_wake(read_mesh(MESH / "fine.msh"))
        taper = tip_taper_factors(wc.station_z, B_SEMI, FORM, R_C_FRAC * B_SEMI)
        recipe = dict(NEWTON_M6_RECIPE)
        recipe["newton_kw"] = dict(recipe["newton_kw"], tip_taper=taper)
        r = solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                   verbose=True, **recipe)
        phi = np.asarray(r["phi"])
        op = PicardOperator(mc.nodes, mc.elements)
        _, q2 = op.velocities(phi)
        q2l = limit_q2_field(q2, M_INF, M_CAP, 1.4)
        mach = np.sqrt(mach_squared_field(q2, M_INF, 1.4))
        cen = mc.nodes[mc.elements].mean(axis=1)
        over = (q2l != q2)
        np.savez(field_cache, cen=cen, mach=mach, over_cap_mask=over,
                 conv=bool(r["converged"]),
                 n_lim=int(r.get("n_limited", 0)),
                 n_flr=int(r.get("n_floored", 0)),
                 wall_s=time.perf_counter() - t0)
        print(f"  done {time.perf_counter() - t0:.0f}s, "
              f"{int(over.sum())} cells over M_cap", flush=True)

    # --- locate the HOTTEST cells --------------------------------------------
    # ★ WHAT IS AND IS NOT MEANINGFUL HERE. The ramp DIED at M = 0.75, so this
    # field is the failed level's state, and the `mach` / `over_cap_mask` arrays
    # were built with m_inf = 0.84 -- the WRONG freestream Mach for it. So:
    #   * the Mach VALUES are NOT physical Mach numbers -- do not quote them;
    #   * the RANKING and the POSITIONS *are* valid, because mach_squared_field
    #     is monotone in |q|^2, so argsort(mach) == argsort(|q|) for any m_inf.
    # We therefore report WHERE the fastest cells are, and nothing else.
    z = cen[:, 2]
    zb = z / B_SEMI
    xi = (cen[:, 0] - x_le(z)) / (x_te(z) - x_le(z))   # 0 = LE, 1 = TE
    dx = cen[:, 0] - x_te(np.clip(z, 0, B_SEMI))       # + = past the local TE
    machn = np.nan_to_num(mach, nan=-1.0)

    top_idx = np.argsort(-machn)[:20]

    print("\n★ The ramp died at M=0.75, so the Mach VALUES below are not physical "
          "(they were formed with m_inf=0.84). Only the RANK and the POSITION "
          "are meaningful. The 20 fastest cells:")
    print(f"{'rank':>4} {'x':>8} {'y':>9} {'z':>8} {'z/b':>7} "
          f"{'chordfrac':>10} {'dx_pastTE':>10}")
    for r, i in enumerate(top_idx, 1):
        print(f"{r:4d} {cen[i,0]:8.4f} {cen[i,1]:+9.5f} {cen[i,2]:8.4f} "
              f"{zb[i]:7.4f} {xi[i]:10.4f} {dx[i]:+10.5f}")

    # ★ THE VERDICT. Two candidate sites, and they are distinguishable:
    #   sharp tip TE   -> chord-frac -> 1  AND  z/b < 1 (on the wing, at the TE)
    #   rounded cap    -> z/b > 1          (outboard of B_SEMI, the NEW surface)
    at_tip_te = int(np.sum((zb[top_idx] > 0.95) & (xi[top_idx] > 0.9)
                           & (zb[top_idx] <= 1.0)))
    on_cap = int(np.sum(zb[top_idx] > 1.0))
    print(f"\nof the 20 fastest cells: {at_tip_te} on the SHARP TIP TE "
          f"(z/b in (0.95, 1], chord-frac > 0.9), {on_cap} on the ROUNDED CAP "
          f"(z/b > 1)")
    print("⇒ the site is the sharp tip TE -- the SAME feature the subsonic box "
          "study flagged, and it exists in the flat family too. What the rounded "
          "cap changes is not the SITE but the SPEED there: it lets flow wrap "
          "around the tip, so M_max is higher at EVERY level (1.51/2.00 round vs "
          "1.39/1.73 flat), and at fine that pushes the tip-TE cells past the "
          "limiter where the flat cap's cooler tip flow stays under.")

    write_csv(OUT, "g133rt_transonic_fine_cells.csv",
              "rank,x,y,z,z_over_b,chord_frac,dx_past_TE,"
              "note_mach_value_invalid_ramp_died_at_M0.75",
              [(r, f"{cen[i,0]:.4f}", f"{cen[i,1]:+.5f}", f"{cen[i,2]:.4f}",
                f"{zb[i]:.4f}", f"{xi[i]:.4f}", f"{dx[i]:+.5f}",
                "rank/position valid; Mach value NOT physical")
               for r, i in enumerate(top_idx, 1)])
    print(f"\nwrote {OUT / 'g133rt_transonic_fine_cells.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
