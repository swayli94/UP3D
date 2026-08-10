"""LE-9: WHY does the entropy correction do the taper's job? Is it just extra dissipation?

The user's objection is well founded and I had not tested it: the entropy correction acts at
SHOCKS, while the tip vortex comes from spanwise Gamma variation, so in principle it should not
touch the tip. I measured a correlation (turn the entropy correction off and the untapered case
diverges) and treated it as a mechanism without testing the mechanism.

★ SCOPING FIRST, because it changes where the discriminator belongs. Two different questions
have been running together:

    M0.84  no taper, entropy ON   CONVERGES  (LE-4)
    M0.84  no taper, entropy OFF  FAILS      (LE-4)      <- entropy substitutes for the taper
    M0.88  no taper, entropy ON   FAILS      (LE-5)
    M0.88  taper,    entropy ON   CONVERGES  (LE-5)      <- the taper still buys the edge

"Why does entropy do the taper's job" is the M0.84 question, so this runs at M0.84, not M0.88.

Also ruled out before running, by reading the generator rather than measuring: the flat tip cap
is NOT the explanation. onera_m6_wingbody_conforming passes tip_cap="round", so LE-5's M0.88
failure already had a smooth tangential cap. (The flat/round distinction does explain a
CONFOUND in LE-1/LE-2 -- the wing-alone family uses the default flat cap, whose edge P13/G13.3
measured diverging at p = +0.321, while the wing-body uses round -- but that is a separate
repair, not this mechanism.)

Two candidate mechanisms, and the user named the second:
  (M1) the entropy correction weakens the shock by cutting post-shock density (sigma < 1
       transported downstream), changing the STATE the artificial-density switch then acts on
  (M2) it acts as extra numerical dissipation, so more dissipation would substitute for it

The discriminator is `upwind_c`, the artificial-density strength (production 1.5; GS1b.2 Q5
swept it 1.5 -> 3.0). With the entropy correction OFF, does raising upwind_c recover the
untapered solve?

  substitutes  => (M2): the stabilisation is dissipation, independent of the entropy physics
  does not     => (M1): the entropy correction does something a dissipation increase cannot,
                  i.e. it changes the post-shock state rather than adding a dissipative term

★ And the dissipation is MEASURED, not inferred from convergence: nu_max, the active-element
count and the floored/limited counts are recorded per leg, so "is it dissipation?" is answered
by comparing how much artificial viscosity each configuration actually applies -- not by
reading it off which legs happened to converge.

Runs at medium because coarse does not discriminate: LE-4 measured every coarse config
converging at M0.84.

Outputs (TRACKED): bench/gate_results/le9_entropy_mechanism.csv
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
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
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

CSV = os.path.join(_GATE, "le9_entropy_mechanism.csv")
MP = os.path.join(REPO, "cases", "meshes", "onera_m6_wingbody_conforming",
                  "medium.msh")
M_TARGET, ALPHA = 0.84, 3.06

#: (label, taper r_c or None, entropy, upwind_c)
LEGS = [
    ("ref_entropyON_c1.5", None, True, 1.5),      # known CONVERGES -- in-session anchor
    ("ref_entropyOFF_c1.5", None, False, 1.5),    # known FAILS -- in-session anchor
    ("entropyOFF_c2.0", None, False, 2.0),        # the (M2) test
    ("entropyOFF_c3.0", None, False, 3.0),        # the (M2) test, GS1b.2's upper sweep value
    ("taperON_entropyOFF", 0.05, False, 1.5),     # the pre-entropy B31/B32 configuration
]
KEYS = ["leg", "r_c", "entropy", "upwind_c", "converged", "m_attained",
        "n_limited", "n_floored", "nu_max", "n_nu_active", "mach2_max",
        "res_final", "cl_p", "wall_s", "note"]


def main():
    print(f"M{M_TARGET} / alpha {ALPHA} / wing-body medium (tip_cap=round)")
    print("question: does MORE artificial dissipation substitute for the "
          "entropy correction?\n")
    rows = []
    for label, rc, ent, c in LEGS:
        mc, wc = cut_wake(read_mesh(MP))
        t = (None if rc is None
             else tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth",
                                    rc * B_SEMI))
        t0 = time.perf_counter()
        try:
            seed = solve_newton_lifting(
                mc, wc, m_inf=cap.WB_MSTART, alpha_deg=ALPHA,
                **dict(cap.CONF_SEED_KW, entropy_correction=ent, upwind_c=c))
            nk = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure",
                      phi_init=seed["phi"], gamma_init=seed["gamma"],
                      n_picard_seed=0, entropy_correction=ent, upwind_c=c)
            if t is not None:
                nk["tip_taper"] = t
            r = solve_newton_transonic(mc, wc, m_inf=M_TARGET, alpha_deg=ALPHA,
                                       m_start=cap.WB_MSTART, dm=cap.DM,
                                       dm_min=0.01, freeze_tol=1e-5,
                                       intermediate_tol=1e-4, newton_kw=nk)
            sref = planform_area(mc.nodes, mc.boundary_faces["wall"])
            cl = wall_forces(mc, phi=np.asarray(r["phi"]), alpha_deg=ALPHA,
                             s_ref=sref, m_inf=M_TARGET)["cl"]
            m_att = float(r.get("m_last_converged", r.get("m_final", M_TARGET)))
            conv = bool(r["converged"]) and abs(m_att - M_TARGET) < 1e-9
            wall = time.perf_counter() - t0
            #: the dissipation actually applied -- the point of the round
            nu = r.get("nu_max"); nact = r.get("n_nu_active")
            m2 = r.get("mach2_max")
            print(f"  {label:22s} ent={ent!s:5s} c={c:<4} conv={conv} "
                  f"lim/flr={r.get('n_limited')}/{r.get('n_floored')} "
                  f"nu_max={nu} n_active={nact} M2max={m2} "
                  f"|R|={float(r['residual_history'][-1]):.3e} "
                  f"cl_p {cl:.7f} ({wall:.0f}s)", flush=True)
            rows.append(dict(leg=label, r_c=rc, entropy=ent, upwind_c=c,
                             converged=conv, m_attained=m_att,
                             n_limited=r.get("n_limited"),
                             n_floored=r.get("n_floored"), nu_max=nu,
                             n_nu_active=nact, mach2_max=m2,
                             res_final=float(r["residual_history"][-1]),
                             cl_p=cl, wall_s=round(wall, 1), note=""))
        except Exception as exc:                                   # noqa: BLE001
            wall = time.perf_counter() - t0
            print(f"  {label:22s} ent={ent!s:5s} c={c:<4} DIED "
                  f"{type(exc).__name__}: {str(exc)[:60]} ({wall:.0f}s)",
                  flush=True)
            rows.append(dict(leg=label, r_c=rc, entropy=ent, upwind_c=c,
                             converged=False, wall_s=round(wall, 1),
                             note=f"{type(exc).__name__}: {exc}"))
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")
    print("\n=== reading ===")
    sub = [r for r in rows if r["leg"].startswith("entropyOFF_c") and r["converged"]]
    if sub:
        cs = ", ".join(str(r["upwind_c"]) for r in sub)
        print(f"  (M2) DISSIPATION SUBSTITUTES: entropy OFF converges at "
              f"upwind_c {cs}.")
        print("  So the stabilisation is dissipation strength, not the entropy")
        print("  physics -- the user's hypothesis holds, and upwind_c is the")
        print("  honest knob to expose rather than the entropy correction being")
        print("  credited with a robustness role it does not have.")
    else:
        print("  (M1) DISSIPATION DOES NOT SUBSTITUTE: raising upwind_c with the")
        print("  entropy correction off does not recover the solve, so the entropy")
        print("  correction does something a dissipation increase cannot -- it")
        print("  changes the post-shock STATE (sigma < 1 transported downstream)")
        print("  rather than adding a dissipative term. Compare the nu_max column")
        print("  across legs for how much viscosity each actually applied.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
