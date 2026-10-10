"""Is the 3-D grid-convergence failure lift-coupled, or is it the wall geometry?

Pre-registered in phases/p2/docs/dev_phase_two/20260803-2000-alpha-discriminator-prereg.md, committed
before this file was written. Read that first -- it records the two premises the user
corrected (distance to a VISCOUS experiment is not an accuracy measure for a full-potential
solver, and this phase is about the capability boundary and finding bugs).

2.5-D conforming sits exactly at first order (R = 0.492/0.500/0.510/0.509 against its
ladder's p=1 value of 0.500) while no 3-D combination reaches it. Both share every kernel,
the assembly, the upwinding, the Newton and the Kutta machinery, so the defect lives in what
only exists in 3-D. Alpha splits the candidates unevenly: three known 3-D-only defect
classes are LIFT-COUPLED and vanish as alpha -> 0 (the B23 junction pocket, the P13/B31 tip
sheet termination, the B19-B21 mixed-plain aux density contamination), while only the G1.6
faceted-wall error is thickness-driven and present at alpha=0.

Instrument: R = ||Cp(medium) - Cp(coarse)|| / ||Cp(coarse) - Cp(xcoarse)||, PER CHORD BAND,
as an arc-length L2 integrated by trapezoid on a dense fixed quadrature. NOT R on cl -- the
ONERA M6 uses the symmetric ONERA D section, so cl == 0 at alpha=0 and a cl-based
discriminator would be degenerate at precisely the decisive point, auto-confirming the
lift-coupled hypothesis. And NOT a pooled point-sample norm: that first attempt was invalid
and the pre-registered control caught it -- see the BANDS note below.

M 0.50 deliberately: fully subsonic, the artificial-density upwinding never arms, so any
non-convergence there is not a transonic artefact.

METRIC PRIORITY (user directive 2026-08-03), and why my first choice was wrong: I had proposed
making spurious cd the PRIMARY metric, because at subsonic with no shock d'Alembert fixes its
exact value at zero, so it measures discretisation error directly instead of through a ratio.
That reference only exists in the shock-free subsonic corner -- and the target regime is
transonic M 0.3-0.87, where cd is physically nonzero. I had picked a metric validated exactly
where the solver is easy and undefined where it matters. So:

  PRIMARY   cl (via R) and section Cp (via per-band R). Standard practice, and the two
            quantities this project's gates are written against. The banded-Cp instrument
            carries ~10 % uncertainty (see the BANDS note) -- do not read finer than that,
            and cl has no exact reference, so R is all there is.
  ANOMALY   cd. At subsonic with no shock the exact value IS zero, so cd is a bug detector
            with a known answer: if it does not fall toward zero under refinement, that is a
            located defect, not a convergence rate. Healthy baseline measured on the NACA
            control: 0.008956 -> 0.004566 -> 0.001429, order 0.97 then 1.68. Once shocks are
            present the test changes to "does cd converge to a stable value", not to zero.

Outputs (TRACKED): bench/gate_results/alpha_discriminator.csv
                   bench/gate_results/alpha_discriminator_R.csv
"""

import csv
import math
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
from pyfp3d.post.surface import planform_area                       # noqa: E402
from pyfp3d.post.unified import section_cp, wall_forces              # noqa: E402

OUT = os.path.join(_GATE)
CSV_PTS = os.path.join(OUT, "alpha_discriminator.csv")
CSV_R = os.path.join(OUT, "alpha_discriminator_R.csv")

M_INF = 0.50
#: ★★ THE FIRST INSTRUMENT WAS INVALID AND THE PRE-REGISTERED CONTROL CAUGHT IT.
#: A pooled point-sample norm over the whole chord gave R = 0.52..1.19 on the SAME three
#: solutions depending only on the sampling grid (uniform60 -> 1.19 against uniform200 ->
#: 0.57): the LE suction peak ALIASES under sparse point sampling. So the control failed its
#: "must be near first order" condition and the round was void as designed.
#:
#: I then mis-read the broken instrument as physics -- claiming the LE "dominates the norm
#: and converges slowly", and calling it independent corroboration of G1.6. It is not: with a
#: proper quadrature the LE BAND converges at R = 0.477, i.e. AT first order. What "dropping
#: x<0.05 recovers first order" actually showed was that the ALIASING is concentrated at the
#: LE, not that the LE converges slowly.
#:
#: Redesign: per-BAND (S2's committed bands), arc-length L2 via trapezoid on a DENSE fixed
#: grid. Validated the same way the first one was invalidated -- R moves 0.15-0.60 % across
#: n_quad 100..800 -- and the control now passes in all three bands (LE 0.477 / MID 0.325 /
#: TE 0.610 against this ladder's first-order 0.500 and ceiling 1.000). Per band, not pooled,
#: because the three rates genuinely differ and pooling hides which region carries the error.
BANDS = (("LE", 0.0, 0.15), ("MID", 0.15, 0.85), ("TE", 0.85, 1.01))
N_QUAD = 400

#: (geom, path, mdir, solver, (xcoarse, coarse, medium), (h_x, h_c, h_m), etas, alphas)
CASES = [
    ("m6wing", "conforming", "onera_m6", cap.conf_wing,
     ("xcoarse_ss", "coarse_ss", "medium"), (0.060, 0.030, 0.015),
     (0.20, 0.44, 0.65, 0.80, 0.90), (0.0, 1.0, 2.0, 3.06)),
    ("m6wing", "level-set", "onera_m6_wakefree", cap.ls_wing,
     ("xcoarse_ss", "coarse_ss", "medium"), (0.060, 0.030, 0.015),
     (0.20, 0.44, 0.65, 0.80, 0.90), (0.0, 1.0, 2.0, 3.06)),
    ("wingbody", "conforming", "onera_m6_wingbody_conforming", cap.conf_wingbody,
     ("xcoarse", "coarse", "medium"), (0.044, 0.030, 0.015),
     (0.20, 0.44, 0.65, 0.80, 0.90), (0.0, 1.0, 2.0, 3.06)),
    ("wingbody", "level-set", "onera_m6_wingbody", cap.ls_wingbody,
     ("xcoarse", "coarse", "medium"), (0.044, 0.030, 0.015),
     (0.20, 0.44, 0.65, 0.80, 0.90), (0.0, 1.0, 2.0, 3.06)),
    #: ★ CONTROL. Must come back healthy at both alphas; if it does not, the instrument or
    #: the sampling grid is at fault and the 3-D numbers must NOT be interpreted.
    ("naca2.5d", "conforming", "naca0012_2.5d", cap.conf_wing,
     ("xcoarse", "coarse", "medium"), (0.040, 0.020, 0.010),
     (None,), (0.0, 1.25)),
]


def sections(mesh, kw, m_eff, etas):
    """Raw section-Cp dicts per station -- no resampling here, so the quadrature choice stays
    inside band_diff_norm where it is validated."""
    out = {}
    for eta in etas:
        if eta is None:                      # 2.5-D: the mid-span slice
            z = 0.5 * float(np.ptp(mesh.nodes[:, 2]))
            out[eta] = section_cp(mesh, z=z, m_inf=m_eff, **kw)
        else:
            out[eta] = section_cp(mesh, eta=eta, b_semi=cap.B_SEMI, m_inf=m_eff, **kw)
    return out


def band_diff_norm(sa, sb, lo, hi):
    """L2 norm of (Cp_a - Cp_b) over one chord band, trapezoid-integrated on a dense FIXED
    grid, both sides and all stations summed in quadrature."""
    xs = np.linspace(lo, min(hi, 1.0), N_QUAD)
    tot = 0.0
    for eta in sa:
        for side in ("upper", "lower"):
            v = []
            for s in (sa, sb):
                x, cp = (np.asarray(s[eta][f"x_{side}"]),
                         np.asarray(s[eta][f"cp_{side}"]))
                o = np.argsort(x)
                v.append(np.interp(xs, x[o], cp[o]))
            tot += float(np.trapz((v[0] - v[1]) ** 2, xs))
    return math.sqrt(tot)


def append(path, row, keys):
    head = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        if head:
            w.writeheader()
        w.writerow(row)


PT_KEYS = ["geom", "path", "alpha", "level", "h_wall", "n_nodes", "status",
           "converged", "res_final", "m_max", "cl_p", "cd_p", "sawtooth",
           "wall_s", "note"]
R_KEYS = ["geom", "path", "alpha", "band", "d1_cp", "d2_cp", "R_field",
           "R_first_order", "R_falsify_ceiling", "verdict", "d1_cl", "d2_cl",
           "R_cl", "note"]


def main():
    only = os.environ.get("PYFP3D_ALPHA_ONLY", "")
    for geom, path, mdir, fn, levels, hs, etas, alphas in CASES:
        if only and f"{geom}/{path}" != only:
            continue
        hx, hc, hm = hs
        r_first = (hm - hc) / (hc - hx)
        r_max = math.log(hm / hc) / math.log(hc / hx)
        print(f"\n=== {geom} / {path}  levels {levels}  h {hs} "
              f"(p=1 -> R {r_first:.3f}, falsified at >= {r_max:.3f}) ===", flush=True)
        for a in alphas:
            vecs, cls, cds, ok = {}, {}, {}, True
            for lv, h in zip(levels, hs):
                mesh_path = os.path.join(REPO, "cases", "meshes", mdir, f"{lv}.msh")
                if not os.path.exists(mesh_path):
                    print(f"  a{a} {lv}: mesh missing", flush=True); ok = False; break
                t0 = time.perf_counter()
                try:
                    mesh, op, r, phi, mvop = fn(mesh_path, M_INF, a)
                except Exception as exc:                           # noqa: BLE001
                    print(f"  a{a:<5} {lv:11s} ERROR {type(exc).__name__}: "
                          f"{str(exc)[:70]}", flush=True)
                    append(CSV_PTS, dict(geom=geom, path=path, alpha=a, level=lv,
                                         h_wall=h, status="ERROR",
                                         note=f"{type(exc).__name__}: {exc}"), PT_KEYS)
                    ok = False; break
                wall = time.perf_counter() - t0
                if mvop is None:
                    kw, m_eff = dict(phi=phi), M_INF
                    conv = bool(r["converged"])
                    res = float(r["residual_history"][-1])
                else:
                    mf = r.get("m_final", M_INF)
                    kw, m_eff = dict(mvop=mvop, phi_ext=phi), mf
                    conv = bool(r.get("target_reached", False))
                    res = float(r["levels"][-1]["residual_norm"])
                try:
                    v = sections(mesh, kw, m_eff, etas)
                except Exception as exc:                           # noqa: BLE001
                    print(f"  a{a:<5} {lv:11s} CP ERROR {type(exc).__name__}: "
                          f"{str(exc)[:60]}", flush=True); ok = False; break
                #: cl recorded as the SECONDARY reading only -- it is identically zero at
                #: alpha=0 on these symmetric sections, which is exactly why the binding
                #: instrument is the Cp field norm.
                try:
                    row, _pl = cap._postprocess(f"{geom}_{path}", path, geom, lv, a,
                                                M_INF, wall, mesh, op, r, phi, mvop)
                    clp, mmax = row["cl_p"], row["m_max"]
                except Exception:                                  # noqa: BLE001
                    clp, mmax = float("nan"), float("nan")
                #: cd = the ANOMALY criterion (exact 0 here: subsonic, no shock). Not a
                #: convergence measure -- see the metric-priority note in the docstring.
                try:
                    sref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
                    f = wall_forces(mesh, alpha_deg=a, s_ref=sref, m_inf=m_eff, **kw)
                    cdp = float(f["cd_pressure"])
                except Exception:                                  # noqa: BLE001
                    cdp = float("nan")
                #: high-frequency content of the raw section Cp, so the sawtooth the
                #: instrument has to survive is recorded per level rather than assumed.
                try:
                    x0 = np.asarray(v[etas[0]]["x_upper"])
                    c0 = np.asarray(v[etas[0]]["cp_upper"])[np.argsort(x0)]
                    saw = float(np.sqrt(np.mean(np.diff(c0, 2) ** 2)))
                except Exception:                                  # noqa: BLE001
                    saw = float("nan")
                vecs[lv], cls[lv], cds[lv] = v, clp, cdp
                print(f"  a{a:<5} {lv:11s} conv={conv} |R|={res:.2e} "
                      f"M_max={mmax} cl_p={clp} cd_p={cdp:.6f} saw={saw:.4f} "
                      f"({wall:.0f}s)", flush=True)
                append(CSV_PTS, dict(geom=geom, path=path, alpha=a, level=lv, h_wall=h,
                                     n_nodes=len(mesh.nodes),
                                     status="OK" if conv else "NOT_CONVERGED",
                                     converged=conv, res_final=res, m_max=mmax,
                                     cl_p=clp, cd_p=cdp, sawtooth=saw,
                                     wall_s=round(wall, 1), note=""), PT_KEYS)
            if not ok or len(vecs) < 3:
                append(CSV_R, dict(geom=geom, path=path, alpha=a, band="-",
                                   note="incomplete ladder"), R_KEYS)
                continue
            lx, lc, lm = levels
            d1c, d2c = cls[lc] - cls[lx], cls[lm] - cls[lc]
            Rc = d2c / d1c if d1c else float("nan")
            for bname, lo, hi in BANDS:
                d1 = band_diff_norm(vecs[lc], vecs[lx], lo, hi)
                d2 = band_diff_norm(vecs[lm], vecs[lc], lo, hi)
                R = d2 / d1 if d1 else float("nan")
                verdict = ("FALSIFIED: R >= this ladder's p->0 ceiling" if R >= r_max
                           else "first-order signature"
                           if abs(R - r_first) < 0.05 * r_first
                           else "converging, order far below first" if R > r_first
                           else "converging, order at or above first")
                print(f"  -> a{a} {bname:3s}: d1 {d1:.6f}  d2 {d2:.6f}  "
                      f"R {R:.4f}  [{verdict}]", flush=True)
                append(CSV_R, dict(geom=geom, path=path, alpha=a, band=bname,
                                   d1_cp=round(d1, 8), d2_cp=round(d2, 8),
                                   R_field=round(R, 5),
                                   R_first_order=round(r_first, 4),
                                   R_falsify_ceiling=round(r_max, 4),
                                   verdict=verdict, d1_cl=round(d1c, 8),
                                   d2_cl=round(d2c, 8), R_cl=round(Rc, 5),
                                   note=""), R_KEYS)
            #: cl R is PRIMARY at alpha > 0 and meaningless at alpha = 0 (cl == 0 by
            #: symmetry). cd is the anomaly flag, never a rate.
            cd_fall = all(cds[levels[i + 1]] < cds[levels[i]] for i in range(2))
            print(f"     cl R {Rc:.3f} "
                  f"{'(PRIMARY)' if a > 0 else '(void at alpha=0: cl == 0 by symmetry)'}"
                  f"   cd {cds[lx]:.6f} -> {cds[lc]:.6f} -> {cds[lm]:.6f} "
                  f"{'OK falling toward 0' if cd_fall else 'ANOMALY: not falling'}",
                  flush=True)
    print(f"\nwrote {CSV_PTS}\nwrote {CSV_R}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
