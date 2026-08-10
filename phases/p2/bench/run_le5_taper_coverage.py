"""LE-5: does removing the tip taper cost any ENVELOPE? The adoption evidence.

LE-4 established the mechanism -- the entropy correction (S1b default) took over the taper's
robustness job, measured by a single-knob A/B: untapered M0.84 medium converges with
entropy_correction=True (0/0 clamps, |R| 2.070e-14) and does NOT without it (|R| 7.273e-06).
And every taper config, including `none`, reached M0.84 cleanly at both levels.

But one Mach and one alpha is not adoption evidence. B31/B32 justified the taper on an envelope
claim (the conforming wing-body medium ceiling M0.79 -> M0.84), so removing it has to be tested
on the envelope, at the conditions where the tip loading is HARDEST:

  ceiling legs   M 0.86, 0.88 at alpha 3.06 -- above the level B32 claimed
  alpha legs     alpha 4.0, 5.0 at M 0.84   -- B23 measured the junction pocket growing
                 SUPERLINEARLY with alpha, and the P13 tip singularity scales with loading,
                 so higher alpha is the stress direction, not lower

Both configs (`none` vs the production vanish_smooth/0.05), both mesh levels.

ADOPTION RULE, fixed before running: removing the taper is adopted only if `none` converges
wherever the production taper does. A single condition where the taper converges and `none`
does not means the taper still earns its bias somewhere, and adoption is refused (or scoped to
the conditions where it does not).

★ Concern #3 from the LE-4 adoption list is already resolved by code reading, no run needed:
NewtonWorkspace guards the B31 Gamma-pin blend behind `_taper_active` at all three sites
(newton.py:315 row scaling, :444 residual branch, :563 Jacobian diagonal), so removing the
taper returns exactly the plain pre-B31 pressure row -- the blend and kutta_weld_sign are never
touched. The comments say as much: "taper == 1 leaves K untouched" and "the untapered Kutta
residual bit-for-bit".

Outputs (TRACKED): bench/gate_results/le5_taper_coverage.csv
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

CSV = os.path.join(_GATE, "le5_taper_coverage.csv")
MDIR = "onera_m6_wingbody_conforming"
CONFIGS = [("none", 0.0), ("vanish_smooth", 0.05)]
#: (label, m_inf, alpha)
LEGS = [("ceiling_M086", 0.86, 3.06), ("ceiling_M088", 0.88, 3.06),
        ("alpha_4", 0.84, 4.00), ("alpha_5", 0.84, 5.00)]
KEYS = ["level", "form", "r_c", "leg", "m_inf", "alpha", "converged",
        "m_attained", "n_limited", "n_floored", "res_final", "cl_p", "wall_s",
        "note"]


def run_leg(mesh_path, form, r_c, m_target, alpha):
    mc, wc = cut_wake(read_mesh(mesh_path))
    t = (None if form == "none"
         else tip_taper_factors(wc.station_z, B_SEMI, form, r_c * B_SEMI))
    t0 = time.perf_counter()
    seed = solve_newton_lifting(mc, wc, m_inf=cap.WB_MSTART, alpha_deg=alpha,
                                **cap.CONF_SEED_KW)
    nk = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure",
              phi_init=seed["phi"], gamma_init=seed["gamma"], n_picard_seed=0)
    if t is not None:
        nk["tip_taper"] = t
    r = solve_newton_transonic(mc, wc, m_inf=m_target, alpha_deg=alpha,
                               m_start=cap.WB_MSTART, dm=cap.DM, dm_min=0.01,
                               freeze_tol=1e-5, intermediate_tol=1e-4,
                               newton_kw=nk)
    return mc, r, time.perf_counter() - t0


def append(row):
    head = not os.path.exists(CSV)
    with open(CSV, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
        if head:
            w.writeheader()
        w.writerow(row)


def main():
    levels = os.environ.get("PYFP3D_LE5_LEVELS", "coarse,medium").split(",")
    verdict = {}
    for level in levels:
        mp = os.path.join(REPO, "cases", "meshes", MDIR, f"{level}.msh")
        if not os.path.exists(mp):
            print(f"{level}: mesh missing"); continue
        print(f"\n########## {level} ##########", flush=True)
        for label, m, a in LEGS:
            for form, r_c in CONFIGS:
                tag = f"{form}/{r_c}"
                try:
                    mc, r, wall = run_leg(mp, form, r_c, m, a)
                    m_att = float(r.get("m_last_converged",
                                        r.get("m_final", m)))
                    conv = bool(r["converged"]) and abs(m_att - m) < 1e-9
                    sref = planform_area(mc.nodes, mc.boundary_faces["wall"])
                    cl = wall_forces(mc, phi=np.asarray(r["phi"]),
                                     alpha_deg=a, s_ref=sref, m_inf=m)["cl"]
                    print(f"  {label:14s} {tag:20s} M{m} a{a}  conv={conv} "
                          f"m_att={m_att} lim/flr={r.get('n_limited')}/"
                          f"{r.get('n_floored')} "
                          f"|R|={float(r['residual_history'][-1]):.2e} "
                          f"cl_p {cl:.7f} ({wall:.0f}s)", flush=True)
                    append(dict(level=level, form=form, r_c=r_c, leg=label,
                                m_inf=m, alpha=a, converged=conv,
                                m_attained=m_att, n_limited=r.get("n_limited"),
                                n_floored=r.get("n_floored"),
                                res_final=float(r["residual_history"][-1]),
                                cl_p=cl, wall_s=round(wall, 1), note=""))
                    verdict[(level, label, form)] = conv
                except Exception as exc:                           # noqa: BLE001
                    print(f"  {label:14s} {tag:20s} M{m} a{a}  DIED "
                          f"{type(exc).__name__}: {str(exc)[:60]}", flush=True)
                    append(dict(level=level, form=form, r_c=r_c, leg=label,
                                m_inf=m, alpha=a, converged=False,
                                note=f"{type(exc).__name__}: {exc}"))
                    verdict[(level, label, form)] = False

    print("\n=== ADOPTION RULE: does `none` converge wherever the taper does? ===")
    bad = []
    for (level, label, form), ok in verdict.items():
        if form != "none":
            continue
        t_ok = verdict.get((level, label, "vanish_smooth"))
        if t_ok and not ok:
            bad.append((level, label))
        print(f"  {level:7s} {label:14s} none={ok!s:5s} taper={t_ok!s:5s}"
              f"{'   ★ taper still earns its bias here' if (t_ok and not ok) else ''}")
    if bad:
        print(f"\n  => ADOPTION REFUSED (or must be scoped): {bad}")
    else:
        print("\n  => `none` converges wherever the taper does, on every tested "
              "condition. Adoption evidence complete for this coverage.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
