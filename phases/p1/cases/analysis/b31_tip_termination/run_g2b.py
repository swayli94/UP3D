"""GB31.2b -- COARSE validation of the pressure-estimator tip-taper port
(the B31 Gamma-pin row blend in solve/newton.py; pre-reg GB31.2b).

Three legs on the wing-body CONFORMING coarse mesh at M0.5 alpha=3.06
(subsonic, one strict solve_newton_lifting level each, the production
level kwargs of wb30.CONF_LEVEL_KW minus the ramp-only Mach chain):

  leg P   pressure, no taper    baseline (production recipe, untapered)
  leg R   probe + taper         the proven G13.2/GB31.2a treatment with
                                the probe estimator (reference profile)
  leg Q   pressure + taper      the new B31 row blend (this gate)

Same taper array on the tapered legs: tip_taper_factors(wc.station_z,
B_SEMI, "vanish_smooth", 0.05*B_SEMI) -- the F3/G13.2 parameters.

Expectation (pre-registered in the GB31.2b task brief): pressure+taper
~ probe+taper within estimator delta, tip-most station Gamma ~ unloaded,
all legs converged. The blend unloads HARDER mid-taper than the probe's
exact Gamma_b = t*Gamma_Kutta (Schur-slope note in
tests/test_b31_pressure_taper.py) -- the gamma_tail rows make the
difference visible station by station.

Run:  python cases/analysis/b31_tip_termination/run_g2b.py
Artifact: results/g2b_coarse_pressure_taper.csv (sections "leg" and
"gamma_tail", the wb31 section-tagged row pattern).
"""

import csv
import time

import numpy as np

from wb31 import OUT  # noqa: F401  (sets sys.path for wb30/pyfp3d)
import wb30  # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI  # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting  # noqa: E402

FORM, R_C_FRAC = "vanish_smooth", 0.05
M_INF = 0.5
N_TAIL = 6                       # stations reported in the gamma_tail rows

LEGS = {
    "P": dict(estimator="pressure", taper=False),   # production baseline
    "R": dict(estimator="probe", taper=True),       # G13.2 reference
    "Q": dict(estimator="pressure", taper=True),    # the B31 port
}


def run_leg(mc, wc, tag, estimator, taper):
    """One cold-started strict M0.5 level (default 5-step Picard seed, so
    every leg freezes sigma/weld at a comparable seed state)."""
    kw = dict(wb30.CONF_LEVEL_KW)
    kw["kutta_estimator"] = estimator
    if taper:
        kw["tip_taper"] = tip_taper_factors(wc.station_z, B_SEMI, FORM,
                                            R_C_FRAC * B_SEMI)
    t0 = time.perf_counter()
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=wb30.ALPHA,
                             verbose=True, **kw)
    wall_s = time.perf_counter() - t0
    clamps = r["clamp_history"]
    phi, gamma = r["phi"], np.asarray(r["gamma"], dtype=float).ravel()
    m2 = wb30.mach2_conf(r["workspace"], phi, gamma, M_INF)
    cents = mc.nodes[mc.elements].mean(axis=1)
    tip = cents[:, 2] > 0.95 * B_SEMI
    rec = dict(tag=tag, estimator=estimator, taper=bool(taper),
               converged=bool(r["converged"]),
               accept_reason=str(r["accept_reason"]),
               res=float(r["residual_history"][-1]),
               n_newton=int(r["n_newton"]),
               n_limited=int(clamps[-1][0]) if clamps else 0,
               n_floored=int(clamps[-1][1]) if clamps else 0,
               cl_p=wb30.cl_p_conf(mc, phi, M_INF),
               gamma_mean=float(gamma.mean()),
               gamma_tip=float(gamma[-1]),
               tip_box_mmax=float(np.sqrt(m2[tip].max())),
               wall_s=wall_s,
               _gamma=gamma)
    print(f"  [{tag}] converged={rec['converged']} "
          f"res={rec['res']:.2e} n_newton={rec['n_newton']} "
          f"cl_p={rec['cl_p']:.4f} g_mean={rec['gamma_mean']:.4f} "
          f"g_tip={rec['gamma_tip']:.2e} "
          f"tipM={rec['tip_box_mmax']:.3f} ({wall_s:.0f}s)", flush=True)
    return rec


def main():
    mc, wc = wb30.build_conf("coarse")
    z = np.asarray(wc.station_z, dtype=float)
    f = tip_taper_factors(z, B_SEMI, FORM, R_C_FRAC * B_SEMI)
    print(f"coarse: n_st={wc.n_stations}, taper tail={f[-N_TAIL:]}",
          flush=True)
    rows = []
    for tag, spec in LEGS.items():
        print(f"=== GB31.2b leg {tag}: {spec['estimator']} "
              f"taper={spec['taper']} ===", flush=True)
        rec = run_leg(mc, wc, tag, **spec)
        row = {k: v for k, v in rec.items() if not k.startswith("_")}
        row.update(section="leg", m=M_INF)
        rows.append(row)
        for j in range(wc.n_stations - N_TAIL, wc.n_stations):
            rows.append(dict(section="gamma_tail", leg=tag, station=j,
                             z=float(z[j]), taper_f=float(f[j]),
                             gamma=float(rec["_gamma"][j])))
    keys = ["section", "leg", "tag", "m", "estimator", "taper",
            "converged", "accept_reason", "res", "n_newton", "n_limited",
            "n_floored", "cl_p", "gamma_mean", "gamma_tip",
            "tip_box_mmax", "wall_s", "station", "z", "taper_f", "gamma"]
    with open(OUT / "g2b_coarse_pressure_taper.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT / 'g2b_coarse_pressure_taper.csv'}")


if __name__ == "__main__":
    main()
