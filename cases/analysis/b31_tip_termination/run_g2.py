"""GB31.2a -- conforming tip-taper FACTORIAL at the 0.83 dying level
(zero library change; pre-reg GB31.2a).

The production CONF recipe runs `kutta_estimator="pressure"` WITHOUT the
G13.2 tip taper (F4: taper is incompatible with the pressure estimator).
This gate isolates the TAPER factor with the probe estimator:

  leg A  0.82 probe no-taper   estimator-delta control (0.82 converges
                               with pressure; RECORDED cl_p/Gamma delta)
  leg B  0.83 probe no-taper   factorial control (does it die like the
                               pressure 0.83?)
  leg C  0.83 probe + taper    treatment (vanish_smooth, r_c=0.05*b_semi,
                               the proven F3/G13.2 parameters)
  leg D  0.84 probe + taper    ONLY IF leg C converged AND leg B died --
                               the pre-reg climb rung (dm=0.0125 from 0.83
                               rounds onto the <=0.84 cap => single rung)

Seeds follow the B30 G2 warm-start chain (0.82 from g2_conf_m0_81, 0.83
from g2_conf_m0_82; leg D chains leg C).  Acceptance = production
semantics.  Clamp masks are estimator-independent (upwind/density only)
so the mask ws keeps the pressure-style workspace for all legs.

Decision tree (pre-locked):
  B died AND C converged (clamps vs 0+2 baseline shrunk) -> **✓ taper is
      the cause** -> climb (leg D) -> GB31.2b (pressure-estimator port)
  C converged but B also converged / ceiling unmoved -> **◐** cost table
  C died or tip clamps unshrunk -> **✗** CONF-side C-class closes

Run:  python cases/analysis/b31_tip_termination/run_g2.py
Artifacts: results/g2_conf_taper.csv, results/g2_conf_taper.png
"""

import csv
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb31 import OUT, classify_conf, gamma_profile_conf
import wb30  # noqa: E402  (importable after wb31's sys.path hook)
from pyfp3d.constraints.wake import tip_taper_factors  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI  # noqa: E402
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting  # noqa: E402

FORM, R_C_FRAC = "vanish_smooth", 0.05
B30_OUT = wb30.OUT
UP = wb30.UPWIND_DEFAULT


def run_level(mc, wc, m_inf, seed_file, taper, tag):
    """One strict CONF level with the PROBE estimator (optionally tapered).
    Cached to results/{tag}.npz (wb30 cache pattern)."""
    cache = OUT / f"{tag}.npz"
    if cache.exists():
        return wb30._cache_load(cache)
    seed = OUT / seed_file
    if not seed.exists():
        seed = B30_OUT / seed_file
    d = np.load(seed, allow_pickle=True)
    phi0, gamma0 = d["phi"], d["gamma"]
    kw = dict(wb30.CONF_LEVEL_KW)
    kw.pop("kutta_estimator")               # -> the default probe estimator
    if taper:
        kw["tip_taper"] = tip_taper_factors(wc.station_z, B_SEMI, FORM,
                                            R_C_FRAC * B_SEMI)
    t0 = time.perf_counter()
    r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=wb30.ALPHA,
                             phi_init=phi0, gamma_init=gamma0,
                             n_picard_seed=0, verbose=True, **kw)
    wall_s = time.perf_counter() - t0
    clamps = r["clamp_history"]
    rec = dict(m_inf=float(m_inf), taper=bool(taper), tag=tag,
               converged=bool(r["converged"]),
               accept_reason=str(r["accept_reason"]),
               residual_norm=float(r["residual_history"][-1]),
               n_newton=int(r["n_newton"]),
               n_limited=int(clamps[-1][0]) if clamps else 0,
               n_floored=int(clamps[-1][1]) if clamps else 0,
               gamma=float(np.mean(r["gamma"])),
               step_nlim=[int(c[0]) for c in clamps],
               step_nflr=[int(c[1]) for c in clamps],
               step_res=[float(x) for x in r["residual_history"]],
               wall_s=wall_s)
    wb30._cache_save(cache, rec, r["phi"], r["gamma"])
    rec["phi"], rec["gamma"] = r["phi"], r["gamma"]
    print(f"  [{tag}] converged={rec['converged']} "
          f"res={rec['residual_norm']:.2e} nlim={rec['n_limited']} "
          f"nflr={rec['n_floored']} ({wall_s:.0f}s)", flush=True)
    return rec


def measure(mc, wc, rec):
    """Cost + guardrail + tip-clamp ownership on a returned state."""
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
    return dict(cl_p=wb30.cl_p_conf(mc, phi, m_inf),
                mmax=float(np.sqrt(m2.max())), corr_mmax=corr_m,
                mask_lim=int(lim.sum()), mask_flr=int(flr.sum()),
                cap_wall_clamps=int(np.count_nonzero(clamp
                                                     & (cls == "cap_wall"))),
                gamma_tip=float(np.asarray(gamma).ravel()[-1]))


def main():
    mc, wc = wb30.build_conf("medium")
    rows = {}
    print("=== GB31.2a leg A: 0.82 probe no-taper (estimator control) ===",
          flush=True)
    rows["A"] = run_level(mc, wc, 0.82, "g2_conf_m0_81.npz", False,
                          "g2a_m0_82_notaper")
    print("=== GB31.2a leg B: 0.83 probe no-taper (factorial control) ===",
          flush=True)
    rows["B"] = run_level(mc, wc, 0.83, "g2_conf_m0_82.npz", False,
                          "g2a_m0_83_notaper")
    print("=== GB31.2a leg C: 0.83 probe + taper (treatment) ===",
          flush=True)
    rows["C"] = run_level(mc, wc, 0.83, "g2_conf_m0_82.npz", True,
                          "g2a_m0_83_taper")
    cure = (not rows["B"]["converged"]) and rows["C"]["converged"]
    if cure:
        print("=== GB31.2a leg D: 0.84 probe + taper (climb cap rung) ===",
              flush=True)
        rows["D"] = run_level(mc, wc, 0.84, "g2a_m0_83_taper.npz", True,
                              "g2a_m0_84_taper")
    out_rows = []
    for leg, rec in rows.items():
        row = dict(leg=leg, tag=rec.get("tag", ""), m=rec["m_inf"],
                   taper=rec["taper"], converged=rec["converged"],
                   accept_reason=rec["accept_reason"],
                   res=rec["residual_norm"], n_newton=rec["n_newton"],
                   n_limited=rec["n_limited"], n_floored=rec["n_floored"],
                   wall_s=round(rec["wall_s"], 1))
        row.update(measure(mc, wc, rec))
        out_rows.append(row)
    keys = ["leg", "tag", "m", "taper", "converged", "accept_reason",
            "res", "n_newton", "n_limited", "n_floored", "wall_s", "cl_p",
            "mmax", "corr_mmax", "mask_lim", "mask_flr", "cap_wall_clamps",
            "gamma_tip"]
    with open(OUT / "g2_conf_taper.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {OUT / 'g2_conf_taper.csv'}")
    _fig(wc, rows)
    print("\n--- verdict hint (pre-reg GB31.2a) ---")
    b, c = rows["B"], rows["C"]
    if cure:
        d_ok = rows.get("D", {}).get("converged")
        print(f"  control died, treatment CONVERGED at 0.83 -> taper is the "
              f"cause ✓; climb cap rung 0.84 converged={d_ok} -> "
              f"GB31.2b (pressure-estimator port) TRIGGERED")
    elif c["converged"]:
        print(f"  BOTH control and treatment converged at 0.83 -> estimator "
              f"effect, taper not isolated -> ◐ cost table to the user")
    else:
        print(f"  treatment died at 0.83 (res={c['residual_norm']:.1e}, "
              f"clamps={c['n_limited']}+{c['n_floored']}) -> ✗ CONF-side "
              f"C-class closes, B10 roll-up remains")


def _fig(wc, rows):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    labels = [f"{k}: M{r['m_inf']}\ntaper={r['taper']}"
              for k, r in rows.items()]
    conv = [r["converged"] for r in rows.values()]
    res = [max(r["residual_norm"], 1e-16) for r in rows.values()]
    colors = ["tab:green" if c else "tab:red" for c in conv]
    ax.bar(range(len(rows)), res, color=colors)
    ax.set_yscale("log")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("final |R|_inf")
    ax.set_title("GB31.2a factorial legs (green=converged)")
    ax.grid(alpha=0.3, axis="y")
    ax = axes[1]
    z = np.asarray(wc.station_z, dtype=float)
    f = tip_taper_factors(z, B_SEMI, FORM, R_C_FRAC * B_SEMI)
    n_tail = max(4, int(0.15 * len(z)))
    for k, r in rows.items():
        if not r["converged"]:
            continue
        g = np.asarray(r["gamma"], dtype=float).ravel()
        ax.plot(z[-n_tail:], g[-n_tail:], "o-", ms=3,
                label=f"{k}: M{r['m_inf']} taper={r['taper']}")
    ax2 = ax.twinx()
    ax2.plot(z[-n_tail:], f[-n_tail:], "k--", lw=1.0,
             label="taper F(z)")
    ax2.set_ylabel("F(z)", color="k")
    ax2.set_ylim(-0.02, 1.05)
    ax.set_xlabel("station z (tip end)")
    ax.set_ylabel("Γ")
    ax.set_title("tip Γ(z), converged legs + taper factor")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "g2_conf_taper.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT / 'g2_conf_taper.png'}")


if __name__ == "__main__":
    main()
