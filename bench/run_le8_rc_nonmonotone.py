"""LE-8: is the taper's robustness really NON-MONOTONE in r_c, or is that freeze scatter?

LE-7 measured, at M0.88 / alpha 3.06 / conforming wing-body medium:

    r_c   0 (none)  FAIL  6/6 clamps  |R| 8.00e-05   1357 s
          0.025     OK    0/0         |R| 1.660e-14   414 s
          0.0375    FAIL  0/0         |R| 4.025e-05   1485 s
          0.05      OK    0/0         |R| 1.780e-14   467 s

So "wider taper = safer" is not true as stated, and that matters for adoption: if robustness is
non-monotone in r_c then r_c = 0.025 passing at M0.88 is ONE SAMPLE, not a trend, and it cannot
be promised as generally robust.

Two candidate explanations, and they have different consequences:

  (A) FREEZE-PATH DEPENDENCE. The project has documented that an inexact iterate path freezes a
      different upwind selection, so each converged state is an exact root of a DIFFERENT
      discrete system; B32 records the conforming wing-body medium ramp needing freeze_tol
      raised 1e-6 -> 1e-5 for exactly this reason. The 0.0375 leg's 1485 s against ~415 s for
      the successes is the signature of a struggling freeze path. If varying freeze_tol flips it
      to converged, the non-monotonicity dissolves and r_c is a clean lever again.
  (B) A REAL DISCRETE EFFECT. The number of TE stations falling inside the taper band is an
      integer, so a genuine resolution/parity effect is possible. Then r_c must be mapped, not
      interpolated, and adoption needs wider coverage.

Legs:
  freeze_tol 1e-4 and 1e-6 at r_c 0.0375   -- the (A) test; the recipe default here is 1e-5
  r_c 0.030 and 0.045 at the default        -- maps the pattern between the known points

Reported alongside: how many TE stations sit inside the band at each r_c, since that is the
quantity (B) would act through, and it costs nothing to record.

Outputs (TRACKED): bench/gate_results/le8_rc_nonmonotone.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import run_capability_matrix as cap                                 # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.post.surface import planform_area                       # noqa: E402
from pyfp3d.post.unified import wall_forces                         # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)

CSV = os.path.join(HERE, "gate_results", "le8_rc_nonmonotone.csv")
MP = os.path.join(REPO, "cases", "meshes", "onera_m6_wingbody_conforming",
                  "medium.msh")
M_TARGET, ALPHA = 0.88, 3.06
#: (r_c, freeze_tol) -- default freeze_tol here is 1e-5 (the B32 wing-body value)
LEGS = [(0.0375, 1e-4), (0.0375, 1e-6), (0.030, 1e-5), (0.045, 1e-5)]
KEYS = ["r_c", "freeze_tol", "n_stations_in_band", "converged", "m_attained",
        "n_limited", "n_floored", "res_final", "cl_p", "wall_s", "note"]


def main():
    print(f"M{M_TARGET} / alpha {ALPHA} / wing-body medium")
    print("LE-7 reference: r_c 0=FAIL, 0.025=OK, 0.0375=FAIL, 0.05=OK "
          "(all at freeze_tol 1e-5)\n")
    rows = []
    for rc, ftol in LEGS:
        mc, wc = cut_wake(read_mesh(MP))
        t = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth", rc * B_SEMI)
        #: the quantity a genuine discrete effect would act through
        n_in = int(np.sum(t < 1.0))
        t0 = time.perf_counter()
        try:
            seed = solve_newton_lifting(mc, wc, m_inf=cap.WB_MSTART,
                                        alpha_deg=ALPHA, **cap.CONF_SEED_KW)
            nk = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure", tip_taper=t,
                      phi_init=seed["phi"], gamma_init=seed["gamma"],
                      n_picard_seed=0)
            r = solve_newton_transonic(mc, wc, m_inf=M_TARGET, alpha_deg=ALPHA,
                                       m_start=cap.WB_MSTART, dm=cap.DM,
                                       dm_min=0.01, freeze_tol=ftol,
                                       intermediate_tol=1e-4, newton_kw=nk)
            sref = planform_area(mc.nodes, mc.boundary_faces["wall"])
            cl = wall_forces(mc, phi=np.asarray(r["phi"]), alpha_deg=ALPHA,
                             s_ref=sref, m_inf=M_TARGET)["cl"]
            m_att = float(r.get("m_last_converged", r.get("m_final", M_TARGET)))
            conv = bool(r["converged"]) and abs(m_att - M_TARGET) < 1e-9
            wall = time.perf_counter() - t0
            print(f"  r_c {rc:<7} freeze_tol {ftol:.0e}  n_band={n_in:3d}  "
                  f"conv={conv} lim/flr={r.get('n_limited')}/"
                  f"{r.get('n_floored')} |R|={float(r['residual_history'][-1]):.3e}"
                  f"  cl_p {cl:.7f} ({wall:.0f}s)", flush=True)
            rows.append(dict(r_c=rc, freeze_tol=ftol, n_stations_in_band=n_in,
                             converged=conv, m_attained=m_att,
                             n_limited=r.get("n_limited"),
                             n_floored=r.get("n_floored"),
                             res_final=float(r["residual_history"][-1]),
                             cl_p=cl, wall_s=round(wall, 1), note=""))
        except Exception as exc:                                   # noqa: BLE001
            wall = time.perf_counter() - t0
            print(f"  r_c {rc:<7} freeze_tol {ftol:.0e}  DIED "
                  f"{type(exc).__name__}: {str(exc)[:60]} ({wall:.0f}s)",
                  flush=True)
            rows.append(dict(r_c=rc, freeze_tol=ftol, n_stations_in_band=n_in,
                             converged=False, wall_s=round(wall, 1),
                             note=f"{type(exc).__name__}: {exc}"))
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")
    print("\n=== reading ===")
    flip = [r for r in rows if r["r_c"] == 0.0375 and r["converged"]]
    if flip:
        tols = ", ".join(f"{r['freeze_tol']:.0e}" for r in flip)
        print(f"  (A) FREEZE-PATH DEPENDENCE: r_c 0.0375 converges at freeze_tol "
              f"{tols}, so its failure at 1e-5")
        print("  is a freeze-path artefact, not an r_c property. The")
        print("  non-monotonicity dissolves and r_c is a clean lever.")
    else:
        print("  (B) the r_c 0.0375 failure SURVIVES both freeze_tol changes, so it")
        print("  is not obviously a freeze artefact. Robustness is then genuinely")
        print("  non-monotone in r_c, and r_c = 0.025 passing M0.88 is one sample --")
        print("  adoption needs wider coverage before it can be promised.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
