r"""Two probes written for the 2026-09-05 independent gate audit.

    python bench/studies/gate_audit_20260905/run_probes.py [--probe both|A|B]

Both are CHEAP (probe A ~30 s, probe B ~4 min at 8 threads) and both write a CSV
under `results/`.  Neither touches `pyfp3d/`.

**Probe A -- mesh-realization noise on an EXACT case.**
C06's lifting cylinder has a closed-form answer, so its `cl_rel` is a true error.
`r_far` is swept over [20, 55] chords at FIXED wall spacing: every value is far
enough that the vortex-corrected far field is converged, so `r_far` is an INERT
knob and the only thing it changes is which unstructured triangulation gmsh
returns.  The spread of `cl_rel` over that sweep is therefore the error bar that
a single committed mesh does not carry.  ★ The probe was written to test the
opposite hypothesis (that r_far = 15 truncation dominates C06's residual ~1 %);
it refutes it -- see the note.

**Probe B -- Steinhoff-Jameson non-uniqueness.**
A conservative full-potential model is documented to admit multiple solutions
for a symmetric section at zero incidence in a narrow transonic band near
M ~ 0.85 (Steinhoff & Jameson, AIAA-81-1019; Salas, Jameson & Melnik, NASA
TP-2385, which attributes it to the approximate shock treatment of the
conservative potential model rather than to physics).  pyFP3D's declared
envelope is M_inf 0.3-0.87, and `tests/D/test_D05` gates |cl| at alpha = 0 --
but only at M 0.72 and M 0.75, i.e. BELOW the documented band.

The probe holds mesh, recipe, Mach, incidence, entropy setting and thread count
fixed and varies ONE thing, `n_picard_seed`, which changes nothing but the
initial guess.  It also runs cold vs an upward warm-start chain, and an
entropy-correction A/B.

★ Read only the legs with `converged == 1`: the non-converged ones carry 0
clamps, so they are not the clamping mode, but this probe does not read
`accept_reason` and therefore does not name their mode.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
OUT = Path(__file__).resolve().parent / "results"

from pyfp3d.mesh.reader import read_mesh                        # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                       # noqa: E402
from pyfp3d.meshgen.extrude import extrude_single_layer         # noqa: E402
from pyfp3d.meshgen.planar import airfoil_wake_2d               # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve            # noqa: E402
from pyfp3d.post.shock import shock_report                      # noqa: E402
from pyfp3d.post.surface import (wall_force_coefficients,       # noqa: E402
                                 wall_tangential_gradient_quadratic)
from pyfp3d.solve.newton import solve_newton_lifting            # noqa: E402
from pyfp3d.solve.picard import solve_laplace_lifting           # noqa: E402


def _write(name, rows, fields):
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / name).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  -> {OUT / name}  ({len(rows)} rows)")


# ----------------------------------------------------------------- probe A ---
#: C06's own constants, imported rather than retyped so the probe cannot drift
#: from the gate it is measuring.
from tests.C.test_C06_lifting_cylinder import (A, GAMMA, U_INF,   # noqa: E402
                                               XC, _circle_polyline)


def probe_a(r_fars=(20., 25., 30., 35., 40., 45., 50., 55.),
            hs=(0.04, 0.02)):
    rows = []
    for h in hs:
        for r_far in r_fars:
            p2, tri, eg, ig = airfoil_wake_2d(
                _circle_polyline(max(80, int(round(2.0 / h)))),
                model_name=f"audit_c06_{r_far:g}", r_far=r_far, h_wall=h,
                h_far=150.0 * h, h_wake=3.0 * h,
                dist_min=0.1, dist_max=6.0, wake_dist_max=1.5)
            dz = 2.0 * h
            mc, wc = cut_wake(extrude_single_layer(p2, tri, eg, ig, dz=dz,
                                                   name="audit_c06"))
            n_st = int(np.unique(wc.te_station).size)
            r = solve_laplace_lifting(mc, wc, alpha_deg=0.0, u_inf=U_INF,
                                      gamma_fixed=np.full(n_st, GAMMA))
            phi = np.asarray(r["phi"])
            wf = mc.boundary_faces["wall"]
            wall = np.unique(wf)
            grad = wall_tangential_gradient_quadratic(mc.nodes, wf, phi)
            cp = 1.0 - np.sum(grad[wall] ** 2, axis=1) / U_INF ** 2
            th = np.arctan2(mc.nodes[wall, 1], mc.nodes[wall, 0] - XC)
            err = cp - (1.0 - (2.0 * np.sin(th)
                               + GAMMA / (2 * np.pi * A * U_INF)) ** 2)
            f = wall_force_coefficients(mc.nodes, mc.elements, wf, phi,
                                        alpha_deg=0.0, u_inf=U_INF,
                                        s_ref=dz, m_inf=0.0)
            cl_exact = 2.0 * GAMMA / (U_INF * 2 * A)
            rows.append(dict(
                h_wall=h, r_far=r_far, n_nodes=int(mc.nodes.shape[0]),
                n_wall=int(wall.size), cl=f"{float(f['cl']):.9e}",
                cl_exact=f"{cl_exact:.9e}",
                cl_rel=f"{abs(float(f['cl']) - cl_exact) / cl_exact:.9e}",
                cp_rms=f"{float(np.sqrt(np.mean(err ** 2))):.9e}",
                cp_max=f"{float(np.max(np.abs(err))):.9e}"))
            print(f"  h={h} r_far={r_far:5.1f}  cl_rel = "
                  f"{100 * float(rows[-1]['cl_rel']):7.4f} %", flush=True)
    _write("farfield_realization_c06.csv", rows,
           list(rows[0].keys()))
    for h in hs:
        v = np.array([100 * float(r["cl_rel"]) for r in rows
                      if float(r["h_wall"]) == h])
        print(f"  ==> h={h}: cl_rel% min {v.min():.4f} max {v.max():.4f} "
              f"spread {np.ptp(v):.4f} pp  std {v.std(ddof=1):.4f} pp")


# ----------------------------------------------------------------- probe B ---
from tests.D.test_D05_euler_naca0012 import RECIPE               # noqa: E402

#: the alpha = 0 legs D05 already gates, for scale: |cl| there is ~1.2e-3.
ALPHA0_BASELINE_CL = 1.2e-3


def _case(level):
    mc, wc = cut_wake(read_mesh(REPO / "cases" / "meshes" / "naca0012_2.5d"
                                / f"{level}.msh"))
    return mc, wc, float(np.ptp(mc.nodes[:, 2])), float(np.mean(mc.nodes[:, 2]))


def _solve(mc, wc, dz, zmid, m, *, seed=None, phi0=None, g0=None, ent=True):
    kw = dict(RECIPE)
    kw["entropy_correction"] = ent
    if seed is not None:
        kw["n_picard_seed"] = seed
    if phi0 is not None:
        kw["phi_init"], kw["gamma_init"], kw["n_picard_seed"] = phi0, g0, 0
    t0 = time.time()
    r = solve_newton_lifting(mc, wc, m_inf=m, alpha_deg=0.0, **kw)
    phi = np.asarray(r["phi"])
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], phi, alpha_deg=0.0,
                                u_inf=1.0, s_ref=dz, m_inf=m)
    rep = shock_report(section_cp_curve(mc, phi, z=zmid, smooth_passes=1,
                                        m_inf=m), m)
    xs = [float(rep[s]["x_shock"]) if rep[s].get("has_shock") else float("nan")
          for s in ("upper", "lower")]
    return dict(
        mach=m, cl=float(f["cl"]), cd_pressure=float(f["cd_pressure"]),
        x_shock_upper=xs[0], x_shock_lower=xs[1], dx_shock=xs[0] - xs[1],
        converged=int(bool(r.get("converged"))),
        residual=float(np.asarray(r["residual_history"], float)[-1]),
        n_limited=int(r.get("n_limited", 0)),
        n_floored=int(r.get("n_floored", 0)),
        wall_s=round(time.time() - t0, 2),
    ), phi, np.asarray(r["gamma"])


def probe_b(levels=("coarse", "medium"),
            machs=(0.78, 0.80, 0.82, 0.83, 0.84, 0.85, 0.86),
            seeds=(0, 2, 5, 8, 12), seed_mach=0.86):
    rows = []
    fields = ["level", "n_nodes", "arm", "n_picard_seed", "entropy_correction",
              "mach", "cl", "cl_over_alpha0_baseline", "cd_pressure",
              "x_shock_upper", "x_shock_lower", "dx_shock", "converged",
              "residual", "n_limited", "n_floored", "n_threads", "wall_s"]
    nthreads = os.environ.get("NUMBA_NUM_THREADS", "?")

    def add(level, n, arm, seed, ent, d):
        rows.append(dict(level=level, n_nodes=n, arm=arm, n_picard_seed=seed,
                         entropy_correction=int(ent),
                         cl_over_alpha0_baseline=round(
                             abs(d["cl"]) / ALPHA0_BASELINE_CL, 2),
                         n_threads=nthreads, **d))
        print(f"  {level:6s} {arm:12s} seed={str(seed):>4s} ent={int(ent)} "
              f"M={d['mach']:.2f} cl={d['cl']:+11.6f} dx={d['dx_shock']:+.4f} "
              f"conv={d['converged']} |R|={d['residual']:.1e} "
              f"lim/flr={d['n_limited']}/{d['n_floored']}", flush=True)

    for level in levels:
        mc, wc, dz, zmid = _case(level)
        n = int(mc.nodes.shape[0])
        for m in machs:                                     # cold, entropy ON
            d, _, _ = _solve(mc, wc, dz, zmid, m)
            add(level, n, "cold", RECIPE["n_picard_seed"], True, d)
        phi = g = None                                      # upward warm chain
        for m in machs:
            d, phi, g = _solve(mc, wc, dz, zmid, m, phi0=phi, g0=g)
            add(level, n, "ramp", 0, True, d)
        for m in machs:                                    # cold, entropy OFF
            d, _, _ = _solve(mc, wc, dz, zmid, m, ent=False)
            add(level, n, "cold_no_entropy", RECIPE["n_picard_seed"], False, d)
        for s in seeds:              # ★ the single-variable leg: guess only
            d, _, _ = _solve(mc, wc, dz, zmid, seed_mach, seed=s)
            add(level, n, "seed_sweep", s, True, d)

    _write("nonuniqueness_naca0012_alpha0.csv", rows, fields)
    conv = [r for r in rows if r["arm"] == "seed_sweep" and r["converged"]]
    for level in levels:
        c = [r for r in conv if r["level"] == level]
        if len(c) >= 2:
            cl = np.array([r["cl"] for r in c])
            print(f"  ==> {level} M{seed_mach} alpha=0, {len(c)} CONVERGED legs "
                  f"differing only in the initial guess: cl in "
                  f"[{cl.min():+.6f}, {cl.max():+.6f}], span {np.ptp(cl):.6f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", choices=("both", "A", "B"), default="both")
    a = ap.parse_args()
    if a.probe in ("both", "A"):
        print("probe A -- mesh-realization noise on C06's exact case")
        probe_a()
    if a.probe in ("both", "B"):
        print("probe B -- alpha=0 transonic non-uniqueness")
        probe_b()
