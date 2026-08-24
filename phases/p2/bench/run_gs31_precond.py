"""GS3.1 (b)+(c): direct vs amg -- answer invariance and speed, as a TRACKED artifact.

Pre-registered in phases/p2/docs/dev_phase_two/20260802-0200-gs31-prereg.md. Tolerance 1e-8
relative on cl_p / cl_KJ / LE upper RMS / pooled RMS / M_max; acceptance >= 3x wall.
An answer moving above 1e-8 is a FINDING that blocks the flip, not an acceptable error.

★ This file exists because the first version of this comparison printed to a scratch log
and nothing else. The session's scratchpad was wiped at a session boundary and T1's amg
result -- the one number that decides whether refinement studies are affordable -- was
lost before it was ever read. That is discipline #3 ("evidence needs a committed
artifact") and I had followed it for every other measurement this phase. Everything here
writes to bench/gate_results/gs31_precond.csv as it goes, one row per (leg, precond), so
a wipe costs at most the leg in flight.

Also runs the determinism control the round owes: `direct` twice on the cheapest leg. If
direct reproduces itself bitwise while direct-vs-amg differs at 1e-4, the difference is
amg's; if direct does not reproduce itself, the 1e-8 criterion was never satisfiable and
the frozen selection is path-dependent (the B15/B21 churn class).

Legs are selected with PYFP3D_GS31_LEGS (comma-separated) so an expensive leg can be run
alone; T1's DIRECT wall time is not re-paid -- it is already recorded in
bench/gate_results/tip_allscales.csv at 3950 s -- and is merged with from_prior=True.
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

from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry         # noqa: E402
from pyfp3d.mesh.reader import read_mesh, write_mesh                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared           # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.surface import (cl_kj_3d, planform_area,           # noqa: E402
                                 wall_force_coefficients)
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)
from run_le_response import SCRATCH, build                          # noqa: E402
from run_m3_budget import (ALPHA, ETAS, M6_NEWTON_KW, M_INF,        # noqa: E402
                           N_UNMASKED, band_rms, parse_experiment,
                           station_rms)
from tests.test_p8_newton import NEWTON_M6_RECIPE                   # noqa: E402

OUT = os.path.join(HERE, "gate_results")
os.makedirs(OUT, exist_ok=True)
CSV = os.path.join(OUT, "gs31_precond.csv")

TAPER_FORM, TAPER_RC_FRAC = "vanish_smooth", 0.05
#: (tag, h_wall, h_le, h_te) -- the same meshes the GS2.1 sweep used
LEGS = {
    "Tp5": (0.020, 0.010, 0.010),
    "T0": (0.015, 0.0075, 0.0075),
    "T1": (0.010, 0.00375, 0.005),
}
#: T1's direct leg is already measured and committed (tip_allscales.csv, 3950.4 s);
#: re-paying 66 minutes to reproduce a number we hold would be discipline #2's exact
#: prohibition. Merged with from_prior=True rather than silently reused.
#: ★ EMPTIED 2026-08-02. The T1 entry above was hand-copied from a 4-6 decimal
#: printout, and the script then computed "invariance" against it -- producing a
#: spurious 2e-4 floor for EVERY configuration and a printed conclusion ("the paths
#: froze different selections, 1e-8 unreachable") that was purely my rounding. Worse,
#: it made the real question unanswerable: at T1 amg_tight matches direct to the 6
#: stored digits, which is about 1e-6 relative and does NOT establish the 1e-8 the
#: criterion asks for. Saving 66 minutes cost the adoption basis. A from_prior row must
#: carry full precision or not exist.
PRIOR_DIRECT = {}
TOL, SPEEDUP = 1e-8, 3.0
FIELDS = ("cl_p", "cl_kj", "le_upper", "pooled", "m_max")


#: the amg path is an INEXACT Newton -- Eisenstat-Walker with ew_eta0 = ew_eta_max =
#: 1e-2, i.e. each linear solve is converged to only ~1 %. `direct` solves exactly. So
#: "amg_tight" clamps the forcing to 1e-10 and asks the discriminating question: does the
#: direct-vs-amg gap shrink SMOOTHLY with the linear tolerance (then it is linear-solve
#: inexactness) or does it stay put / jump (then the two paths froze DIFFERENT upwind
#: selections, and each state is an exact root of its own frozen system -- which is why
#: both report |R| ~ 1e-14).
EW_TIGHT = 1e-10
#: calibration ladder. The CRITERION is unchanged (reproduce direct to 1e-8); this
#: finds the LOOSEST forcing that still meets it, which is engineering calibration of a
#: knob, not fitting the criterion to the answer. Run with PYFP3D_GS31_EW=1e-6,1e-5.
EW_LADDER = [float(x) for x in os.environ.get("PYFP3D_GS31_EW", "").split(",") if x]


def solve_with(mc, wc, precond):
    nk = dict(M6_NEWTON_KW)
    nk["precond"] = "amg" if precond.startswith("amg") else precond
    if precond.startswith("amg"):
        nk.pop("direct_refactor_every", None)     # a direct-only knob
    if precond == "amg_tight":
        nk.update(ew_eta0=EW_TIGHT, ew_eta_max=EW_TIGHT)
    elif precond.startswith("amg_cap"):
        # ★ the SHIPPABLE form: press only the CAP (ew_eta_max), leaving ew_eta0 at
        # its 1e-2 default so the first solve -- where _ew_forcing returns eta0
        # unclamped because r_prev is None -- stays cheap. This keeps
        # Eisenstat-Walker's adaptivity (eta = 0.9 (r/r_prev)^2) and changes only how
        # LOOSE it is allowed to be, which is the knob that actually controls the
        # drift. Measured separately from amg_ew* because "ship what you measured".
        nk.update(ew_eta_max=float(precond.split("amg_cap")[1]))
    elif precond.startswith("amg_ew"):
        # fixed forcing: eta0 = eta_max pins every solve and disables the adaptivity
        e = float(precond.split("amg_ew")[1])
        nk.update(ew_eta0=e, ew_eta_max=e)
    taper = tip_taper_factors(wc.station_z, B_SEMI, TAPER_FORM,
                              TAPER_RC_FRAC * B_SEMI)
    r0 = solve_newton_lifting(mc, wc, m_inf=0.70, alpha_deg=ALPHA, **nk)
    kw = dict(NEWTON_M6_RECIPE)
    kw["newton_kw"] = dict(kw["newton_kw"], **nk, kutta_estimator="pressure",
                           tip_taper=taper, phi_init=r0["phi"],
                           gamma_init=r0["gamma"], n_picard_seed=0)
    return solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw)


def metrics(mc, wc, r, exp):
    phi = np.asarray(r["phi"])
    B, _ = precompute_element_geometry(mc.nodes, mc.elements)
    g = np.einsum("eaj,ea->ej", B, phi[mc.elements])
    m2 = mach_number_squared(np.einsum("ej,ej->e", g, g), M_INF)
    cur = {e: section_cp_curve(mc, phi, eta=e, b_semi=B_SEMI, m_inf=M_INF)
           for e in ETAS}
    acc = {}
    for e in ETAS[:N_UNMASKED]:
        for k, (ss, nn) in band_rms(cur, exp, e).items():
            a, n = acc.get(k, (0.0, 0))
            acc[k] = (a + ss, n + nn)
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], phi,
                                alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
    o = np.argsort(wc.station_z)
    return dict(
        cl_p=f["cl"],
        cl_kj=float(cl_kj_3d(np.atleast_1d(r["gamma"])[o], wc.station_z[o],
                             s_ref, B_SEMI)),
        le_upper=(acc["LE_upper"][0] / acc["LE_upper"][1]) ** 0.5,
        pooled=float(np.mean([station_rms(cur, exp, e)[0]
                              for e in ETAS[:N_UNMASKED]])),
        m_max=float(np.sqrt(m2.max())))


def append_row(row):
    """Write as we go -- a wipe then costs at most the leg in flight."""
    head = not os.path.exists(CSV)
    keys = ["leg", "precond", "repeat", "n_tets", "converged", "res_final",
            "n_newton", "wall_s", *FIELDS, "from_prior", "src", "note"]
    with open(CSV, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        if head:
            w.writeheader()
        w.writerow(row)


def main():
    exp = parse_experiment()
    only = [t for t in os.environ.get("PYFP3D_GS31_LEGS", "").split(",") if t]
    got = {}
    for tag, (hw, hle, hte) in LEGS.items():
        if only and tag not in only:
            continue
        cache = os.path.join(SCRATCH, f"gs31_{tag}.msh")
        if os.path.exists(cache):
            mesh = read_mesh(cache)
        else:
            t0 = time.perf_counter()
            mesh = build(hw, hle, hte)
            print(f"  built {tag} in {time.perf_counter()-t0:.0f}s", flush=True)
            try:
                write_mesh(mesh, cache)
            except Exception as exc:                              # noqa: BLE001
                print(f"  (not cached: {exc})")
        mc, wc = cut_wake(mesh)
        print(f"\n=== {tag}: {len(mc.elements)} tets ===", flush=True)
        got[tag] = {}
        # the determinism control: direct twice on the cheapest leg only
        # amg_tight on EVERY leg: Tp5 showed the two criteria pull against each
        # other (loose amg 2.49x but 2e-4 drift; tight amg 1e-13 but 1.23x), and
        # the speedup grows with mesh size, so which setting -- if any -- passes
        # both is decided at the LARGEST leg, not the cheapest.
        plan = [("direct", 0), ("amg", 0), ("amg_tight", 0)]
        plan += [(f"amg_ew{e:g}", 0) for e in EW_LADDER]
        plan += [(f"amg_cap{e:g}", 0)
                 for e in [float(x) for x in
                           os.environ.get("PYFP3D_GS31_CAP", "").split(",") if x]]
        if tag == "Tp5":
            plan.insert(1, ("direct", 1))
        for precond, rep in plan:
            if precond == "direct" and rep == 0 and tag in PRIOR_DIRECT:
                p = PRIOR_DIRECT[tag]
                got[tag][("direct", 0)] = p
                append_row(dict(leg=tag, precond="direct", repeat=0,
                                n_tets=len(mc.elements), converged=True,
                                wall_s=p["wall"], from_prior=True,
                                src=p["src"], **{k: p[k] for k in FIELDS}))
                print(f"  direct[0] FROM PRIOR ({p['src']}): {p['wall']:.0f}s "
                      f"cl_p={p['cl_p']:.9f}", flush=True)
                continue
            t0 = time.perf_counter()
            try:
                r = solve_with(mc, wc, precond)
            except Exception as exc:                              # noqa: BLE001
                print(f"  {precond}[{rep}] FAILED: {type(exc).__name__}: {exc}",
                      flush=True)
                append_row(dict(leg=tag, precond=precond, repeat=rep,
                                n_tets=len(mc.elements), converged=False,
                                note=f"{type(exc).__name__}: {exc}"))
                continue
            wall = time.perf_counter() - t0
            m = metrics(mc, wc, r, exp)
            m["wall"] = wall
            got[tag][(precond, rep)] = m
            append_row(dict(leg=tag, precond=precond, repeat=rep,
                            n_tets=len(mc.elements),
                            converged=bool(r["converged"]),
                            res_final=float(r["residual_history"][-1]),
                            n_newton=r["n_newton"], wall_s=round(wall, 1),
                            from_prior=False, **m))
            print(f"  {precond}[{rep}] conv={r['converged']} "
                  f"|R|={r['residual_history'][-1]:.2e} n={r['n_newton']} "
                  f"({wall:.0f}s)  cl_p={m['cl_p']:.9f} "
                  f"LEu={m['le_upper']:.6f} M_max={m['m_max']:.5f}", flush=True)

        d, a = got[tag].get(("direct", 0)), got[tag].get(("amg", 0))
        if d and a:
            worst = max(abs(a[k] - d[k]) / max(abs(d[k]), 1e-30) for k in FIELDS)
            print(f"  --- invariance: worst rel {worst:.2e} vs tol {TOL:.0e} -> "
                  f"{'PASS' if worst <= TOL else 'FAIL (blocks the flip)'}")
            for k in FIELDS:
                rel = abs(a[k] - d[k]) / max(abs(d[k]), 1e-30)
                print(f"      {k:9s} {d[k]:.9f} vs {a[k]:.9f}  rel {rel:.2e}")
            sp = d["wall"] / a["wall"]
            print(f"  --- speed {d['wall']:.0f}s -> {a['wall']:.0f}s = {sp:.2f}x "
                  f"-> {'>=3x OK' if sp >= SPEEDUP else 'below 3x, recorded'}")
        at = got[tag].get(("amg_tight", 0))
        if d and a and at:
            w_loose = max(abs(a[k]-d[k])/max(abs(d[k]),1e-30) for k in FIELDS)
            w_tight = max(abs(at[k]-d[k])/max(abs(d[k]),1e-30) for k in FIELDS)
            print(f"  --- EW forcing 1e-2 -> {EW_TIGHT:.0e}: worst rel "
                  f"{w_loose:.2e} -> {w_tight:.2e} "
                  f"({w_loose/max(w_tight,1e-30):.0f}x tighter)  "
                  f"wall {a['wall']:.0f}s -> {at['wall']:.0f}s")
            if w_tight <= TOL:
                print("      => the gap IS linear-solve inexactness; amg reaches "
                      "1e-8 with tight forcing, and the flip is possible at that "
                      "setting (speed then read from amg_tight, not amg)")
            elif w_tight > 0.1 * w_loose:
                print("      => the gap does NOT follow the linear tolerance: the "
                      "two paths froze DIFFERENT upwind selections, so 1e-8 is "
                      "unreachable by tightening and the criterion needs re-specifying "
                      "on measurement, not preference")
            else:
                print("      => partial: shrinks with the tolerance but not to 1e-8 "
                      "-- RECORDED, both mechanisms in play")

        d2 = got[tag].get(("direct", 1))
        if d and d2:
            bit = all(d[k] == d2[k] for k in FIELDS)
            print(f"  --- DETERMINISM control: direct twice "
                  f"{'BITWISE IDENTICAL' if bit else 'DIFFERS'}")
            if not bit:
                for k in FIELDS:
                    print(f"      {k:9s} {d[k]!r} vs {d2[k]!r}")
            print(f"      => the 1e-8 criterion is "
                  f"{'satisfiable in principle; the 1e-4 gap is AMG-side' if bit else 'NOT satisfiable -- the frozen selection is path-dependent (B15/B21 churn)'}")
    print(f"\nwrote {CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
