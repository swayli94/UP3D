"""GS1.3: is the Kutta-closure amplification (phase one's 1/(1-b) ~ 14x)
physical or an implementation artefact?

Three measurements (pre-registered in
phases/p2/docs/dev_phase_two/20260728-1930-s1-gamma-amplification.md):

  A  estimator cross-check: sweep the artificial-dissipation constant C with
     BOTH Kutta estimators -- "probe" (potential-jump at a wall probe node) and
     "pressure" (pressure equality over wall-adjacent control volumes, P14).
     They are completely different algebra for the same physical condition, so
     agreement means the amplification belongs to the flow problem.
  B  known-perturbation amplification: on the probe path the closure residual is
     F = tip_taper * T(phi) - Gamma, so tip_taper = 1 + eps is a clean known
     perturbation of the Kutta target. At convergence Gamma = (1+eps) T(Gamma),
     hence A = (dGamma/Gamma)/eps = 1/(1-b). Two eps values check linearity.
  C  subsonic control: the same A at M 0.5 (no shock). If A(transonic) >>
     A(subsonic) the amplification is the transonic feedback, not the algebra.

Outputs: results/gs1_3_gamma_amp.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                      # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve              # noqa: E402
from pyfp3d.post.shock import shock_report                     # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients        # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,         # noqa: E402
                                 solve_newton_transonic)

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA = 1.25
M_TRANS = 0.7875
M_SUB = 0.50
LEVELS = ("coarse", "medium")
CS = (1.0, 1.5, 2.0, 3.0)
EPS = (1.0e-3, 1.0e-2)


def _post(mc, r, m_inf, dz):
    rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=m_inf),
                       m_inf)
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], r["phi"],
                                alpha_deg=ALPHA, s_ref=dz, m_inf=m_inf)
    return rep, f


def solve(mc, wc, m_inf, C, estimator, taper=None):
    """One solve; transonic goes through the Mach ramp, subsonic direct."""
    nk = dict(freeze_refresh_max=8, precond="direct", n_newton_max=60,
              kutta_estimator=estimator)
    if taper is not None:
        nk["tip_taper"] = np.full(wc.n_stations, taper, dtype=np.float64)
    if m_inf > 0.70:
        return solve_newton_transonic(
            mc, wc, m_inf=m_inf, alpha_deg=ALPHA, m_start=0.70, dm=0.025,
            dm_min=0.003, upwind_c=C, freeze_tol=1e-6, newton_kw=nk)
    return solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=ALPHA,
                                upwind_c=C, m_crit=0.95, freeze_tol=1e-6, **nk)


def part_a(rows):
    print("\n=== A: estimator cross-check (dcl/dx_shock) ===", flush=True)
    for level in LEVELS:
        mc, wc = cut_wake(read_mesh(
            REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"))
        dz = float(np.ptp(mc.nodes[:, 2]))
        for est in ("probe", "pressure"):
            for C in CS:
                t0 = time.perf_counter()
                try:
                    r = solve(mc, wc, M_TRANS, C, est)
                except Exception as exc:                       # noqa: BLE001
                    rows.append(dict(part="A", level=level, estimator=est,
                                     C=C, converged=False,
                                     error=type(exc).__name__))
                    print("  ", rows[-1], flush=True)
                    continue
                rep, f = _post(mc, r, M_TRANS, dz)
                done = [lr["m"] for lr in r["level_results"]
                        if lr["converged"]]
                row = dict(part="A", level=level, estimator=est, C=C,
                           converged=r["converged"],
                           m_reached=max(done) if done else None,
                           cl_p=round(f["cl"], 6),
                           gamma=round(float(r["gamma"][0]), 8),
                           x_shock=rep["upper"].get("x_shock"),
                           m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                           n_floored=r["n_floored"], n_limited=r["n_limited"],
                           wall_s=round(time.perf_counter() - t0, 1), error="")
                print("  ", {k: row[k] for k in
                             ("level", "estimator", "C", "converged", "cl_p",
                              "x_shock", "m_max")}, flush=True)
                rows.append(row)


def part_bc(rows):
    print("\n=== B/C: known-perturbation amplification A = (dGamma/Gamma)/eps "
          "===", flush=True)
    for level in LEVELS:
        mc, wc = cut_wake(read_mesh(
            REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"))
        dz = float(np.ptp(mc.nodes[:, 2]))
        for m_inf, tag in ((M_SUB, "subsonic"), (M_TRANS, "transonic")):
            try:
                r0 = solve(mc, wc, m_inf, 1.5, "probe")
            except Exception as exc:                           # noqa: BLE001
                print(f"  {level} {tag}: baseline FAILED "
                      f"{type(exc).__name__}", flush=True)
                continue
            g0 = float(r0["gamma"][0])
            _, f0 = _post(mc, r0, m_inf, dz)
            print(f"  {level} {tag}: baseline gamma {g0:.8f} "
                  f"cl_p {f0['cl']:.6f} conv={r0['converged']}", flush=True)
            for eps in EPS:
                try:
                    r1 = solve(mc, wc, m_inf, 1.5, "probe", taper=1.0 + eps)
                except Exception as exc:                       # noqa: BLE001
                    rows.append(dict(part="BC", level=level, regime=tag,
                                     eps=eps, converged=False,
                                     error=type(exc).__name__))
                    print("  ", rows[-1], flush=True)
                    continue
                g1 = float(r1["gamma"][0])
                _, f1 = _post(mc, r1, m_inf, dz)
                amp = ((g1 - g0) / g0) / eps
                row = dict(part="BC", level=level, regime=tag, eps=eps,
                           converged=r1["converged"],
                           gamma0=round(g0, 8), gamma1=round(g1, 8),
                           dgamma_rel=round((g1 - g0) / g0, 8),
                           amplification=round(amp, 4),
                           cl0=round(f0["cl"], 6), cl1=round(f1["cl"], 6),
                           dcl_rel=round((f1["cl"] - f0["cl"]) / f0["cl"], 6),
                           error="")
                print(f"    eps {eps:<8g} dGamma/Gamma {row['dgamma_rel']:+.6e} "
                      f"-> A = {amp:8.3f}   (dcl/cl {row['dcl_rel']:+.4e})",
                      flush=True)
                rows.append(row)


def report(rows):
    print("\n=== A headline: dcl/dx_shock per estimator ===")
    for level in LEVELS:
        for est in ("probe", "pressure"):
            sub = [r for r in rows if r.get("part") == "A"
                   and r.get("level") == level and r.get("estimator") == est
                   and r.get("x_shock") and r.get("converged")]
            if len(sub) < 2:
                print(f"  {level:7s} {est:9s} too few converged legs")
                continue
            xs = np.array([r["x_shock"] for r in sub])
            cl = np.array([r["cl_p"] for r in sub])
            slope = np.polyfit(xs, cl, 1)[0]
            print(f"  {level:7s} {est:9s} n={len(sub)}  "
                  f"x_shock {xs.min():.4f}..{xs.max():.4f}  "
                  f"cl {cl.min():.4f}..{cl.max():.4f}  "
                  f"dcl/dx_shock = {slope:6.2f} per chord")
    print("\n=== B/C headline: amplification A = 1/(1-b) ===")
    for level in LEVELS:
        for tag in ("subsonic", "transonic"):
            sub = [r for r in rows if r.get("part") == "BC"
                   and r.get("level") == level and r.get("regime") == tag
                   and r.get("amplification") is not None]
            if not sub:
                continue
            vals = ", ".join(f"eps={r['eps']:g}: A={r['amplification']:.3f}"
                             for r in sub)
            print(f"  {level:7s} {tag:10s} {vals}")


def main():
    rows = []
    part_a(rows)
    part_bc(rows)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(OUT / "gs1_3_gamma_amp.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1_3_gamma_amp.csv")
    report(rows)


if __name__ == "__main__":
    main()
