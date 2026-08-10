"""GB31.2b medium gate -- PRODUCTION recipe (pressure estimator) + tip taper
at the conforming dying levels (pre-reg GB31.2b; the 2b library port is in
pyfp3d/solve/newton.py, default-off).

Legs (all strict, kutta_estimator="pressure", taper = vanish_smooth
r_c=0.05*b_semi -- the GB31.2a-proven cure parameters):

  G  0.82 pressure+taper   COST CONTROL: same M, same estimator, +/-taper
                           vs the B30 pressure 0.82 converged state
                           (cl_p 0.2722) -- the F3-style cost readout
  E  0.83 pressure+taper   CURE TEST at the B30 dying level (pressure
                           no-taper died: res 8.1e-6, 0+1..2 tip clamps)
  F  0.84 pressure+taper   CAP RUNG (historical reached target)

Seeds: G from g2_conf_m0_81, E from g2_conf_m0_82, F chains E.
Acceptance = production semantics.  Guardrails: corridor corrM <= 1.3,
clamps (expect 0+0 per 2a), weld-sign flips recorded (2b hazard note:
frozen at the first residual).

Leg F2 (diagnostic, added after F failed with sigma_flips=3): same
0.84 pressure+taper solve but seeded from the HEALTHY probe+taper 0.84
converged state (g2a_m0_84_taper.npz) instead of chaining E.  The weld
sign s_j is frozen at the first residual evaluation; F's seed (E's 0.83
state re-evaluated at 0.84) carried 3 tip stations whose Gamma-
sensitivity sign differed from the mesh median, so those pins may pull
with the wrong orientation for the whole solve.  F2 converging =>
seed/sign-freeze interaction (fixable); F2 failing too => structural
estimator interaction at 0.84.

Decision (pre-reg):  E converges with clean guardrails -> **GB31.2b ✓**,
production-adoption candidate (user's call); E dies -> the probe/pressure
estimator interaction is itself ceiling-relevant -> ◐ report.

Run:  python cases/analysis/b31_tip_termination/run_g2c.py
Artifacts: results/g2c_medium_pressure_taper.csv (+ .png)
"""

import csv
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb31 import OUT, classify_conf
import wb30  # noqa: E402  (importable after wb31's sys.path hook)
from pyfp3d.constraints.wake import tip_taper_factors  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI  # noqa: E402
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting  # noqa: E402

FORM, R_C_FRAC = "vanish_smooth", 0.05
B30_OUT = wb30.OUT
UP = wb30.UPWIND_DEFAULT
LEGS = [("G", 0.82, "g2_conf_m0_81.npz", "g2c_m0_82_ptaper"),
        ("E", 0.83, "g2_conf_m0_82.npz", "g2c_m0_83_ptaper"),
        ("F", 0.84, "g2c_m0_83_ptaper.npz", "g2c_m0_84_ptaper"),
        ("F2", 0.84, "g2a_m0_84_taper.npz", "g2c_m0_84_ptaper_f2")]
# committed B30 pressure reference for the cost readout
REF_082 = dict(cl_p=0.2722, source="b30 g2_conf_m0_82 (pressure, no taper)")


def run_level(mc, wc, m_inf, seed_file, tag):
    """One strict PRODUCTION CONF level (pressure) + tip_taper."""
    cache = OUT / f"{tag}.npz"
    if cache.exists():
        return wb30._cache_load(cache)
    seed = OUT / seed_file
    if not seed.exists():
        seed = B30_OUT / seed_file
    d = np.load(seed, allow_pickle=True)
    phi0, gamma0 = d["phi"], d["gamma"]
    kw = dict(wb30.CONF_LEVEL_KW)     # keeps kutta_estimator="pressure"
    kw["tip_taper"] = tip_taper_factors(wc.station_z, B_SEMI, FORM,
                                        R_C_FRAC * B_SEMI)
    t0 = time.perf_counter()
    r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=wb30.ALPHA,
                             phi_init=phi0, gamma_init=gamma0,
                             n_picard_seed=0, verbose=True, **kw)
    wall_s = time.perf_counter() - t0
    clamps = r["clamp_history"]
    rec = dict(m_inf=float(m_inf), tag=tag,
               converged=bool(r["converged"]),
               accept_reason=str(r["accept_reason"]),
               residual_norm=float(r["residual_history"][-1]),
               n_newton=int(r["n_newton"]),
               n_limited=int(clamps[-1][0]) if clamps else 0,
               n_floored=int(clamps[-1][1]) if clamps else 0,
               gamma=float(np.mean(r["gamma"])),
               kutta_sigma_sign_flips=int(r.get("kutta_sigma_sign_flips",
                                                0)),
               weld_sign=r.get("kutta_weld_sign", None),
               step_nlim=[int(c[0]) for c in clamps],
               step_nflr=[int(c[1]) for c in clamps],
               wall_s=wall_s)
    wb30._cache_save(cache, rec, r["phi"], r["gamma"])
    rec["phi"], rec["gamma"] = r["phi"], r["gamma"]
    print(f"  [{tag}] converged={rec['converged']} "
          f"res={rec['residual_norm']:.2e} nlim={rec['n_limited']} "
          f"nflr={rec['n_floored']} n_newton={rec['n_newton']} "
          f"({wall_s:.0f}s)", flush=True)
    return rec


def measure(mc, wc, rec):
    ws = NewtonWorkspace(mc, wc, wb30.ALPHA, 1.0, 1.4, (0.25, 0.0), True,
                         kutta_estimator="pressure")
    ws.set_mach(rec["m_inf"])
    phi, gamma, m_inf = rec["phi"], rec["gamma"], rec["m_inf"]
    lim, flr = wb30.clamp_masks_conf(ws, phi, gamma, m_inf,
                                     upwind_c=UP["upwind_c"],
                                     m_crit=UP["m_crit"],
                                     m_cap=UP["m_cap"],
                                     rho_floor=UP["rho_floor"])
    m2 = wb30.mach2_conf(ws, phi, gamma, m_inf)
    cents = mc.nodes[mc.elements].mean(axis=1)
    cls = classify_conf(mc, wc, cents)
    clamp = lim | flr
    corr_m, _ = wb30.corridor_mmax(m2, cents)
    g = np.asarray(gamma, dtype=float).ravel()
    return dict(cl_p=wb30.cl_p_conf(mc, phi, m_inf),
                mmax=float(np.sqrt(m2.max())), corr_mmax=corr_m,
                mask_lim=int(lim.sum()), mask_flr=int(flr.sum()),
                cap_wall_clamps=int(np.count_nonzero(clamp
                                                     & (cls == "cap_wall"))),
                gamma_tip=float(g[-1]),
                gamma_root=float(g[0]))


def main():
    mc, wc = wb30.build_conf("medium")
    rows = {}
    for leg, m_inf, seed, tag in LEGS:
        if leg == "F" and not rows["E"]["converged"]:
            print("  [leg F] skipped -- E (0.83) did not converge",
                  flush=True)
            break
        print(f"=== GB31.2b leg {leg}: {m_inf} pressure+taper ===",
              flush=True)
        rows[leg] = run_level(mc, wc, m_inf, seed, tag)
    out_rows = []
    meas = {}
    for leg, rec in rows.items():
        row = dict(leg=leg, tag=rec["tag"], m=rec["m_inf"],
                   converged=rec["converged"],
                   accept_reason=rec["accept_reason"],
                   res=rec["residual_norm"], n_newton=rec["n_newton"],
                   n_limited=rec["n_limited"], n_floored=rec["n_floored"],
                   wall_s=round(rec["wall_s"], 1),
                   sigma_flips=rec["kutta_sigma_sign_flips"])
        meas[leg] = measure(mc, wc, rec)
        row.update(meas[leg])
        out_rows.append(row)
    keys = ["leg", "tag", "m", "converged", "accept_reason", "res",
            "n_newton", "n_limited", "n_floored", "wall_s", "sigma_flips",
            "cl_p", "mmax", "corr_mmax", "mask_lim", "mask_flr",
            "cap_wall_clamps", "gamma_tip", "gamma_root"]
    with open(OUT / "g2c_medium_pressure_taper.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {OUT / 'g2c_medium_pressure_taper.csv'}")
    _fig(wc, rows)
    print("\n--- verdict hint (pre-reg GB31.2b) ---")
    g = rows.get("G")
    e = rows.get("E")
    if g is not None and g["converged"]:
        cost = 100.0 * (meas["G"]["cl_p"] - REF_082["cl_p"]) \
            / REF_082["cl_p"]
        print(f"  leg G cost control: cl_p {meas['G']['cl_p']:.4f} vs "
              f"{REF_082['cl_p']} ({REF_082['source']}) = {cost:+.2f}% "
              f"(F3 band -1.1..-1.6%)")
    if e is not None and e["converged"] and e["n_limited"] == 0 \
            and e["n_floored"] == 0:
        f_ok = rows.get("F", {}).get("converged")
        print(f"  leg E (0.83 dying level) CONVERGED strict with the "
              f"PRODUCTION estimator + taper, 0 clamps; cap rung 0.84 "
              f"converged={f_ok} -> **GB31.2b ✓** production-adoption "
              f"candidate (user's call)")
    elif e is not None:
        print(f"  leg E died under the production estimator "
              f"(res={e['residual_norm']:.1e}, "
              f"clamps={e['n_limited']}+{e['n_floored']}) -> ◐ "
              f"estimator interaction is ceiling-relevant, report to user")


def _fig(wc, rows):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    labels = [f"{k}: M{r['m_inf']}" for k, r in rows.items()]
    res = [max(r["residual_norm"], 1e-16) for r in rows.values()]
    colors = ["tab:green" if r["converged"] else "tab:red"
              for r in rows.values()]
    ax.bar(range(len(rows)), res, color=colors)
    ax.set_yscale("log")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("final |R|_inf")
    ax.set_title("GB31.2b medium gate: pressure+taper (green=converged)")
    ax.grid(alpha=0.3, axis="y")
    ax = axes[1]
    z = np.asarray(wc.station_z, dtype=float)
    n_tail = max(4, int(0.15 * len(z)))
    for k, r in rows.items():
        if not r["converged"]:
            continue
        g = np.asarray(r["gamma"], dtype=float).ravel()
        ax.plot(z[-n_tail:], g[-n_tail:], "o-", ms=3,
                label=f"{k}: M{r['m_inf']}")
    ax.set_xlabel("station z (tip end)")
    ax.set_ylabel("Γ")
    ax.set_title("tip Γ(z), converged pressure+taper legs")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "g2c_medium_pressure_taper.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT / 'g2c_medium_pressure_taper.png'}")


if __name__ == "__main__":
    main()
