"""GB32.1 medium re-validation -- the B32 per-step weld-sign refresh at the
conforming dying/cap levels (pre-reg GB32.1; the fix is in
pyfp3d/solve/newton.py::_refresh_weld_sign, taper-active pressure path only).

Legs (all strict, kutta_estimator="pressure", taper = vanish_smooth
r_c=0.05*b_semi, identical recipes/seeds to b31 run_g2c.py):

  G-re  0.82   re-solved post-fix; sigma_flips == 0 seed => the refresh
               must NEVER write: phi/gamma BIT-IDENTICAL to the b31
               g2c_m0_82_ptaper cache, kutta_weld_sign_updates == 0
  E-re  0.83   same bit-identity gate vs g2c_m0_83_ptaper
  F-fix 0.84   chained from the b31 E state (the seed whose first-eval
               freeze caught 3 transient flip stations -> 97+63 limit
               cycle pre-fix): target = strict converged, 0 clamps,
               cl_p ~= 0.276 (the F2 healthy-seed value), corrM <= 1.3

Decision (pre-reg):  G/E bit-identical AND F-fix converges clean ->
**GB32.1 ✓**, proceed to GB32.2 adoption; F converges but G/E drift ->
side-effect ◐, cost table to the user; F still dies -> rollback ✗.

Run:  python cases/analysis/b32_tip_taper_adoption/run_g1.py
Artifacts: results/g1_weld_fix.csv (+ .png), results/*.npz caches
"""

import csv
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_B31 = Path(__file__).resolve().parent.parent / "b31_tip_termination"
sys.path.insert(0, str(_B31))

import wb31  # noqa: E402  (sys.path hook for the b30 helpers)
import wb30  # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI  # noqa: E402
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting  # noqa: E402

OUT = Path(__file__).parent / "results"
OUT.mkdir(exist_ok=True)
B31_OUT = wb31.OUT
B30_OUT = wb30.OUT
FORM, R_C_FRAC = "vanish_smooth", 0.05
UP = wb30.UPWIND_DEFAULT
# (leg, m_inf, seed, tag, pre-fix reference cache for bit-identity)
LEGS = [("G-re", 0.82, B30_OUT / "g2_conf_m0_81.npz", "g1_m0_82_ptaper_b32",
         B31_OUT / "g2c_m0_82_ptaper.npz"),
        ("E-re", 0.83, B30_OUT / "g2_conf_m0_82.npz", "g1_m0_83_ptaper_b32",
         B31_OUT / "g2c_m0_83_ptaper.npz"),
        ("F-fix", 0.84, B31_OUT / "g2c_m0_83_ptaper.npz",
         "g1_m0_84_ptaper_b32", None)]
F2_REF = dict(cl_p=0.27609, source="b31 g2c F2 (healthy-seed 0.84)")


def run_level(mc, wc, m_inf, seed_path, tag):
    """One strict PRODUCTION CONF level (pressure) + tip_taper, post-fix."""
    cache = OUT / f"{tag}.npz"
    if cache.exists():
        return wb30._cache_load(cache)
    d = np.load(seed_path, allow_pickle=True)
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
               sigma_flips=int(r.get("kutta_sigma_sign_flips", 0)),
               weld_updates=int(r.get("kutta_weld_sign_updates", 0)),
               wall_s=wall_s)
    wb30._cache_save(cache, rec, r["phi"], r["gamma"])
    rec["phi"], rec["gamma"] = r["phi"], r["gamma"]
    print(f"  [{tag}] converged={rec['converged']} "
          f"res={rec['residual_norm']:.2e} nlim={rec['n_limited']} "
          f"nflr={rec['n_floored']} n_newton={rec['n_newton']} "
          f"weld_updates={rec['weld_updates']} ({wall_s:.0f}s)",
          flush=True)
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
    corr_m, _ = wb30.corridor_mmax(m2, cents)
    return dict(cl_p=wb30.cl_p_conf(mc, phi, m_inf),
                mmax=float(np.sqrt(m2.max())), corr_mmax=corr_m,
                mask_lim=int(lim.sum()), mask_flr=int(flr.sum()))


def bit_check(rec, ref_path):
    """phi/gamma bit-identity vs the pre-fix b31 cache."""
    if ref_path is None or not rec["converged"]:
        return ""
    d = np.load(ref_path, allow_pickle=True)
    same = (np.array_equal(d["phi"], rec["phi"])
            and np.array_equal(d["gamma"], rec["gamma"]))
    return bool(same)


def main():
    mc, wc = wb30.build_conf("medium")
    rows = []
    for leg, m_inf, seed, tag, ref in LEGS:
        print(f"=== GB32.1 leg {leg}: {m_inf} pressure+taper post-fix ===",
              flush=True)
        rec = run_level(mc, wc, m_inf, seed, tag)
        row = dict(leg=leg, tag=tag, m=rec["m_inf"],
                   converged=rec["converged"], res=rec["residual_norm"],
                   n_newton=rec["n_newton"], n_limited=rec["n_limited"],
                   n_floored=rec["n_floored"],
                   sigma_flips=rec["sigma_flips"],
                   weld_updates=rec["weld_updates"],
                   wall_s=round(rec["wall_s"], 1),
                   bit_identical_prefix=bit_check(rec, ref))
        row.update(measure(mc, wc, rec))
        rows.append(row)
    keys = ["leg", "tag", "m", "converged", "res", "n_newton", "n_limited",
            "n_floored", "sigma_flips", "weld_updates", "wall_s",
            "bit_identical_prefix", "cl_p", "mmax", "corr_mmax",
            "mask_lim", "mask_flr"]
    with open(OUT / "g1_weld_fix.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT / 'g1_weld_fix.csv'}")
    _fig(rows)
    print("\n--- verdict hint (pre-reg GB32.1) ---")
    ge = [r for r in rows if r["leg"] in ("G-re", "E-re")]
    f = next((r for r in rows if r["leg"] == "F-fix"), None)
    ge_ok = all(r["bit_identical_prefix"] is True
                and r["weld_updates"] == 0 for r in ge)
    f_ok = (f is not None and f["converged"] and f["n_limited"] == 0
            and f["n_floored"] == 0 and f["corr_mmax"] <= 1.3)
    print(f"  G/E bit-identity + zero refresh writes: {ge_ok}")
    if f is not None:
        cost = 100.0 * (f["cl_p"] - F2_REF["cl_p"]) / F2_REF["cl_p"]
        print(f"  F-fix: converged={f['converged']} clamps="
              f"{f['n_limited']}+{f['n_floored']} corrM={f['corr_mmax']:.3f} "
              f"cl_p={f['cl_p']:.5f} vs {F2_REF['source']} = {cost:+.2f}% "
              f"weld_updates={f['weld_updates']}")
    if ge_ok and f_ok:
        print("  -> **GB32.1 ✓** chained-seed 0.84 cured by the refresh; "
              "proceed to GB32.2 adoption")
    elif f_ok:
        print("  -> ◐ F fixed but G/E paths drifted (refresh side effect) "
              "-- cost table to the user")
    else:
        print("  -> ✗ F still fails -- rollback candidate, report to user")


def _fig(rows):
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    labels = [f"{r['leg']}\nM{r['m']}" for r in rows]
    res = [max(r["res"], 1e-16) for r in rows]
    colors = ["tab:green" if r["converged"] else "tab:red" for r in rows]
    ax.bar(range(len(rows)), res, color=colors)
    for i, r in enumerate(rows):
        ax.annotate(f"upd={r['weld_updates']}", (i, res[i]),
                    ha="center", va="bottom", fontsize=8)
    ax.set_yscale("log")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("final |R|_inf")
    ax.set_title("GB32.1: per-step weld-sign refresh "
                 "(green=converged; upd=applied sign flips)")
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(OUT / "g1_weld_fix.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT / 'g1_weld_fix.png'}")


if __name__ == "__main__":
    main()
