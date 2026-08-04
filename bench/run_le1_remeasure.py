"""LE-1: re-measure the LE deficit with the validated instrument, on all six combinations.

Uses bench/le_instrument.py's band_lift -- the integrated instrument that passed LE-0's
machine-precision wiring guard (bands sum to the library's own cl bit-identically) and that
carries no sampling noise of its own.

★ Registered before running, because it decides how the phase continues: THE DEFICIT MAY NOT
SURVIVE. The old instrument was measured too weak to resolve effects below ~15 %, and P11
independently measured that a wall-Cp order collapse of exactly this kind was an artefact of
refining one length scale while the bulk mesh stayed fixed -- with all scales refined, order
came back at 1.89-1.98. These _ss ladders DO refine all scales by 2x. So "the LE deficit
dissolves" is a live outcome, and it is a RESULT, not a failure: it would retire the LE phase's
premise and correct several readings from the last two rounds.

What is being read, per (geometry, wake path, alpha):
  R per band on the band-integrated lift, against the ladder's own first-order value and its
  p -> 0 falsification ceiling. Plus the band DECOMPOSITION of the total, since LE-0 already
  showed the total's non-monotonicity is bands cancelling rather than noise -- a sign the
  previous sign-blind |dCp| norm could not have seen.

Outputs (TRACKED): bench/gate_results/le1_bands.csv   (one row per geom/path/alpha/level)
                   bench/gate_results/le1_R.csv       (one row per geom/path/alpha/band)
"""

import csv
import math
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import run_capability_matrix as cap                                 # noqa: E402
from le_instrument import band_lift                                 # noqa: E402
from pyfp3d.post.surface import planform_area                       # noqa: E402
from pyfp3d.post.unified import wall_forces                         # noqa: E402

OUT = os.path.join(HERE, "gate_results")
CSV_B = os.path.join(OUT, "le1_bands.csv")
CSV_R = os.path.join(OUT, "le1_R.csv")
M_INF = 0.50
ALPHAS = (0.0, 3.06)

#: (geom, path, mesh dir, solver, chord mode, levels, h ladder)
CASES = [
    ("naca2.5d", "conforming", "naca0012_2.5d", cap.conf_wing, "flat",
     ("xcoarse", "coarse", "medium"), (0.040, 0.020, 0.010)),
    ("naca2.5d", "level-set", "naca0012_wakefree_2.5d", cap.ls_naca, "flat",
     ("xcoarse", "coarse", "medium"), (0.040, 0.020, 0.010)),
    #: ★ FLAT-cap wing-alone ladders, kept as the superseded reference. wing3d records
    #: P13/G13.3 measuring the flat cap DIVERGING under refinement (p = +0.321), so an order
    #: computed on these has a false premise -- these rows exist to be compared against the
    #: round-tip rows below, not to be quoted.
    ("m6wing_flat", "conforming", "onera_m6", cap.conf_wing, "wing3d",
     ("xcoarse_ss", "coarse_ss", "medium"), (0.060, 0.030, 0.015)),
    ("m6wing_flat", "level-set", "onera_m6_wakefree", cap.ls_wing, "wing3d",
     ("xcoarse_ss", "coarse_ss", "medium"), (0.060, 0.030, 0.015)),
    #: ★ ROUND-tip ladders (2026-08-04) -- the ones whose orders may be quoted, since the
    #: tip is tangential with no edge and does not diverge under refinement.
    ("m6wing", "conforming", "onera_m6_roundtip", cap.conf_wing, "wing3d",
     ("xcoarse", "coarse", "medium"), (0.060, 0.030, 0.015)),
    ("m6wing", "level-set", "onera_m6_wakefree", cap.ls_wing, "wing3d",
     ("xcoarse_rt", "coarse_rt", "medium_rt"), (0.060, 0.030, 0.015)),
    ("wingbody", "conforming", "onera_m6_wingbody_conforming", cap.conf_wingbody,
     "wing3d", ("xcoarse", "coarse", "medium"), (0.044, 0.030, 0.015)),
    ("wingbody", "level-set", "onera_m6_wingbody", cap.ls_wingbody, "wing3d",
     ("xcoarse", "coarse", "medium"), (0.044, 0.030, 0.015)),
]

B_KEYS = ["geom", "path", "alpha", "level", "h_wall", "n_tri", "cl_LE", "cl_MID",
          "cl_TE", "cl_OUTSIDE", "cl_total", "v1_sum_err", "v1_ref_err",
          "converged", "wall_s", "note"]
R_KEYS = ["geom", "path", "alpha", "band", "d1", "d2", "R", "R_first_order",
          "R_falsify_ceiling", "verdict", "note"]


def append(path, row, keys):
    head = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        if head:
            w.writeheader()
        w.writerow(row)


def main():
    only = os.environ.get("PYFP3D_LE1_ONLY", "")
    for geom, path, mdir, fn, cmode, levels, hs in CASES:
        if only and f"{geom}/{path}" != only:
            continue
        hx, hc, hm = hs
        r_first = (hm - hc) / (hc - hx)
        r_max = math.log(hm / hc) / math.log(hc / hx)
        print(f"\n=== {geom} / {path}   h {hs}  (p=1 -> R {r_first:.3f}, "
              f"falsified at >= {r_max:.3f}) ===", flush=True)
        for a in ALPHAS:
            vals, ok = {}, True
            for lv, h in zip(levels, hs):
                mp = os.path.join(REPO, "cases", "meshes", mdir, f"{lv}.msh")
                if not os.path.exists(mp):
                    print(f"  a{a} {lv}: mesh missing", flush=True); ok = False; break
                t0 = time.perf_counter()
                try:
                    mesh, op, r, phi, mvop = fn(mp, M_INF, a)
                except Exception as exc:                           # noqa: BLE001
                    print(f"  a{a:<5} {lv:11s} ERROR {type(exc).__name__}: "
                          f"{str(exc)[:70]}", flush=True)
                    append(CSV_B, dict(geom=geom, path=path, alpha=a, level=lv,
                                       h_wall=h, note=f"{type(exc).__name__}"), B_KEYS)
                    ok = False; break
                wall = time.perf_counter() - t0
                kw = (dict(phi=phi) if mvop is None
                      else dict(mvop=mvop, phi_ext=phi))
                m_eff = M_INF if mvop is None else r.get("m_final", M_INF)
                conv = (bool(r["converged"]) if mvop is None
                        else bool(r.get("target_reached", False)))
                sref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
                b = band_lift(mesh, alpha_deg=a, m_inf=m_eff, s_ref=sref,
                              geom=cmode, **kw)
                s = b["LE"] + b["MID"] + b["TE"] + b["OUTSIDE"]
                ref = wall_forces(mesh, alpha_deg=a, s_ref=sref, m_inf=m_eff,
                                  **kw)["cl"]
                e1, e2 = abs(s - b["TOTAL"]), abs(b["TOTAL"] - ref)
                #: the guard is not decoration -- a banding or normalisation slip would
                #: make every R below meaningless, so it is checked on EVERY row.
                if e1 > 1e-12 or e2 > 1e-12:
                    print(f"  a{a:<5} {lv:11s} ★ V1 GUARD FAILED sum {e1:.2e} "
                          f"ref {e2:.2e}", flush=True)
                print(f"  a{a:<5} {lv:11s} conv={conv} n_tri={b['n_tri']:6d}  "
                      f"LE {b['LE']:+.7f}  MID {b['MID']:+.7f}  "
                      f"TE {b['TE']:+.7f}  TOT {b['TOTAL']:+.7f}  ({wall:.0f}s)",
                      flush=True)
                vals[lv] = b
                append(CSV_B, dict(geom=geom, path=path, alpha=a, level=lv, h_wall=h,
                                   n_tri=b["n_tri"], cl_LE=b["LE"], cl_MID=b["MID"],
                                   cl_TE=b["TE"], cl_OUTSIDE=b["OUTSIDE"],
                                   cl_total=b["TOTAL"], v1_sum_err=e1,
                                   v1_ref_err=e2, converged=conv,
                                   wall_s=round(wall, 1), note=""), B_KEYS)
            if not ok or len(vals) < 3:
                continue
            lx, lc, lm = levels
            for name in ("LE", "MID", "TE", "TOTAL"):
                d1 = vals[lc][name] - vals[lx][name]
                d2 = vals[lm][name] - vals[lc][name]
                R = d2 / d1 if d1 else float("nan")
                flip = (d1 > 0) != (d2 > 0)
                verdict = ("non-monotone: coarsest level outside the asymptotic range"
                           if flip else
                           "FALSIFIED: R >= this ladder's p->0 ceiling" if R >= r_max
                           else "first-order signature"
                           if abs(R - r_first) < 0.05 * r_first
                           else "converging, order far below first" if R > r_first
                           else "converging, order at or above first")
                print(f"    -> a{a} {name:6s} d1 {d1:+.8f}  d2 {d2:+.8f}  "
                      f"R {R:+.4f}  [{verdict}]", flush=True)
                append(CSV_R, dict(geom=geom, path=path, alpha=a, band=name,
                                   d1=round(d1, 10), d2=round(d2, 10),
                                   R=round(R, 5), R_first_order=round(r_first, 4),
                                   R_falsify_ceiling=round(r_max, 4),
                                   verdict=verdict, note=""), R_KEYS)
    print(f"\nwrote {CSV_B}\nwrote {CSV_R}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
