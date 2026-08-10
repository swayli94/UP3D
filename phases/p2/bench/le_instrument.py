"""LE-0: the replacement convergence instrument -- band-integrated lift contribution.

The previous instrument (per-band L2 of section-Cp differences, sampled at 5 spanwise
stations on a fixed x/c grid) was measured too weak to answer the question it was built for:
station-to-station R scatter gave the pooled mean a 10.6-18.6 % standard error, the same size
as the 18-24 % effect, and the section-cut sawtooth put high-frequency noise into the
level-to-level difference. Both failure modes are structural, so the fix is structural.

This instrument integrates instead of sampling:

    cl_band = -(Cp * A) @ n_out / S_ref  projected on the lift direction,
              restricted to wall triangles whose centroid sits in a chordwise band

It removes every noise source the old one had at once -- no section cutting, no spanwise
station sampling, no interpolation onto a chosen grid -- and the area integration averages
the sawtooth out rather than differencing two realisations of it. The only remaining noise is
the solver's own frozen-selection scatter (~1e-5 in cl), which is 100-1000x below the
level-to-level differences.

★ It also carries a wiring guard the old instrument structurally could not: the band
contributions must sum EXACTLY to the total cl that wall_forces reports, because they are the
same integral over a partition of the same triangles. A silent banding or normalisation bug
therefore cannot survive.

Validation required before any physics is read off it (LE-0):
  V1  band contributions sum to the total cl to round-off
  V2  the 2.5-D NACA control reproduces ~first order, since that case is known healthy
  V3  determinism -- same solution, same number (the instrument itself adds no noise)

Outputs (TRACKED): bench/gate_results/le_instrument_validation.csv
"""

import csv
import math
import os
import sys

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

from pyfp3d.meshgen.wing3d import chord_at, x_le                    # noqa: E402
from pyfp3d.post.surface import _pressure_force, planform_area      # noqa: E402
from pyfp3d.post.unified import (_cp_from_q2,                       # noqa: E402
                                 _conforming_wall_state,
                                 _d11_wall_state, wall_forces)

#: S2's committed bands, so the decomposition is comparable to the existing error budget.
BANDS = (("LE", 0.0, 0.15), ("MID", 0.15, 0.85), ("TE", 0.85, 1.001))


def band_lift(mesh, *, phi=None, mvop=None, phi_ext=None, alpha_deg=0.0,
              m_inf=0.0, s_ref=1.0, wall_tag="wall", u_inf=1.0,
              geom="wing3d", smooth_passes=0):
    """Per-band lift contribution + the total, from ONE pass over the wall triangles.

    geom = "wing3d"  -> chord-normalise per spanwise station via wing3d.x_le/chord_at
    geom = "flat"    -> x/c = x (the 2.5-D airfoil already sits on x in [0, 1])
    """
    wall = np.asarray(mesh.boundary_faces[wall_tag], dtype=np.int64)
    if phi is not None:
        q2, _, area, n_out = _conforming_wall_state(mesh, phi, wall, u_inf,
                                                    smooth_passes)
    else:
        q2, _, area, n_out = _d11_wall_state(mesh, mvop, phi_ext, wall, u_inf,
                                             smooth_passes)
    cp = _cp_from_q2(q2, m_inf, 1.4)
    cen = mesh.nodes[wall].mean(axis=1)
    if geom == "wing3d":
        z = cen[:, 2]
        xc = (cen[:, 0] - np.array([x_le(zi) for zi in z])) / \
             np.array([chord_at(zi) for zi in z])
    else:
        xc = cen[:, 0]
    out = {}
    for name, lo, hi in BANDS:
        m = (xc >= lo) & (xc < hi)
        if not np.any(m):
            out[name] = 0.0
            continue
        _cf, cl, _cd = _pressure_force(cp[m], area[m], n_out[m], s_ref, alpha_deg)
        out[name] = cl
    #: everything OUTSIDE the bands, so the sum check is exact rather than approximate --
    #: a triangle whose centroid falls beyond x/c = 1 (the tip cap wraps around) would
    #: otherwise silently vanish and make the guard pass for the wrong reason.
    m_out = ~((xc >= BANDS[0][1]) & (xc < BANDS[-1][2]))
    if np.any(m_out):
        _cf, cl, _cd = _pressure_force(cp[m_out], area[m_out], n_out[m_out],
                                       s_ref, alpha_deg)
        out["OUTSIDE"] = cl
    else:
        out["OUTSIDE"] = 0.0
    _cf, cl_tot, _cd = _pressure_force(cp, area, n_out, s_ref, alpha_deg)
    out["TOTAL"] = cl_tot
    out["n_tri"] = len(wall)
    return out


def main():
    import run_capability_matrix as cap
    rows = []

    CASES = [
        ("naca2.5d", "conforming", "naca0012_2.5d", cap.conf_wing, "flat",
         ("xcoarse", "coarse", "medium"), (0.040, 0.020, 0.010), 1.25),
        ("m6wing", "conforming", "onera_m6", cap.conf_wing, "wing3d",
         ("xcoarse_ss", "coarse_ss", "medium"), (0.060, 0.030, 0.015), 3.06),
    ]
    M_INF = 0.50
    for geom, path, mdir, fn, cmode, levels, hs, alpha in CASES:
        print(f"\n=== {geom} / {path}  alpha {alpha}  M{M_INF} ===", flush=True)
        vals = {}
        for lv, h in zip(levels, hs):
            mp = os.path.join(REPO, "cases", "meshes", mdir, f"{lv}.msh")
            if not os.path.exists(mp):
                print(f"  {lv}: mesh missing"); continue
            mesh, op, r, phi, mvop = fn(mp, M_INF, alpha)
            kw = (dict(phi=phi) if mvop is None else dict(mvop=mvop, phi_ext=phi))
            m_eff = M_INF if mvop is None else r.get("m_final", M_INF)
            sref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
            b = band_lift(mesh, alpha_deg=alpha, m_inf=m_eff, s_ref=sref,
                          geom=cmode, **kw)
            #: V1 -- the guard the old instrument could not have
            s = b["LE"] + b["MID"] + b["TE"] + b["OUTSIDE"]
            ref = wall_forces(mesh, alpha_deg=alpha, s_ref=sref, m_inf=m_eff,
                              **kw)["cl"]
            err_sum = abs(s - b["TOTAL"])
            err_ref = abs(b["TOTAL"] - ref)
            print(f"  {lv:11s} n_tri={b['n_tri']:6d}  LE {b['LE']:+.8f}  "
                  f"MID {b['MID']:+.8f}  TE {b['TE']:+.8f}  OUT {b['OUTSIDE']:+.8f}",
                  flush=True)
            print(f"  {'':11s} V1 sum-vs-total {err_sum:.3e}   "
                  f"total-vs-wall_forces {err_ref:.3e}", flush=True)
            vals[lv] = b
            rows.append(dict(geom=geom, path=path, level=lv, h_wall=h, alpha=alpha,
                             n_tri=b["n_tri"], cl_LE=b["LE"], cl_MID=b["MID"],
                             cl_TE=b["TE"], cl_OUTSIDE=b["OUTSIDE"],
                             cl_total=b["TOTAL"], v1_sum_err=err_sum,
                             v1_ref_err=err_ref))
        if len(vals) < 3:
            continue
        lx, lc, lm = levels
        hx, hc, hm = hs
        r_first = (hm - hc) / (hc - hx)
        r_max = math.log(hm / hc) / math.log(hc / hx)
        print(f"  --- R per band (p=1 -> {r_first:.3f}, falsified at >= "
              f"{r_max:.3f}) ---", flush=True)
        for name in ("LE", "MID", "TE", "TOTAL"):
            d1 = vals[lc][name] - vals[lx][name]
            d2 = vals[lm][name] - vals[lc][name]
            R = d2 / d1 if d1 else float("nan")
            print(f"    {name:7s} d1 {d1:+.8f}  d2 {d2:+.8f}  R {R:+.4f}", flush=True)
            rows.append(dict(geom=geom, path=path, level="R", alpha=alpha,
                             cl_LE=name, cl_total=R, v1_sum_err=d1, v1_ref_err=d2))
    p = os.path.join(_GATE, "le_instrument_validation.csv")
    with open(p, "w", newline="") as fh:
        keys = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
