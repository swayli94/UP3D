"""GB31.1 -- tip-termination atlas (cache-only; zero solves, zero lib change).

Inputs = committed B30 caches (cases/analysis/b30_transonic_ceiling/results):
  LS  : g2_ls_m0_7875.npz   (the c=1.5 dying state, 5 lim + 1 flr)
  CONF: g2_conf_m0_83.npz   (the c=1.5 dying state, 0 lim + 2 flr)

Deliverables (pre-reg GB31.1):
  results/g1_tip_atlas.csv  -- sections clamp_cells / clamp_owner /
                               ring_profile / straddler / gamma
  results/g1_tip_atlas.png  -- (a) LS ring delta(q), (b) CONF Gamma(z) +
                               taper overlay, (c) clamp ownership pivot
PASS = every dying-level clamp cell gets an unambiguous ownership class
AND the atlas lands on disk.

Run:  python cases/analysis/b31_tip_termination/run_g1.py
"""

import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb31 import (B_SEMI, OUT, TIP_TE, clamp_cell_rows, classify_conf,
                  classify_ls, gamma_profile_conf, ownership_summary,
                  ring_profile_ls, straddler_census_ls, tip_distances)
import wb30  # noqa: E402  (importable after wb31's sys.path hook)
from pyfp3d.solve.newton import NewtonWorkspace  # noqa: E402

UP = wb30.UPWIND_DEFAULT
B30_OUT = wb30.OUT
M_DYING = {"LS": 0.7875, "CONF": 0.83}


def leg_ls():
    mesh = wb30.load_mesh(wb30.LS_MESH_DIR / "medium.msh")
    wls, cm, mvop = wb30.build_ls_flat(mesh)
    m_inf = M_DYING["LS"]
    d = np.load(B30_OUT / "g2_ls_m0_7875.npz", allow_pickle=True)
    phi = d["phi"]
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    cls = classify_ls(mesh, cm, cents)
    lim, flr = wb30.clamp_masks_ls(mesh, mvop, phi, m_inf,
                                   upwind_c=UP["upwind_c"],
                                   m_crit=UP["m_crit"], m_cap=UP["m_cap"],
                                   rho_floor=UP["rho_floor"])
    m2 = wb30.mach2_ls(mvop, phi, m_inf)
    cells = clamp_cell_rows("LS", dict(limited=lim, floored=flr),
                            cls, m2, cents)
    ring = ring_profile_ls(mesh, wls, cm, mvop, phi, cents, m_inf)
    strad = straddler_census_ls(mesh, cm, mvop, phi, cents, m_inf, cls)
    print(f"  [LS] clamps: {int(lim.sum())} lim + {int(flr.sum())} flr; "
          f"owners: {[(r['cls'], r['count']) for r in ownership_summary(cells)]}",
          flush=True)
    if ring:
        last = ring[-1]
        print(f"  [LS] last-ring |delta|: mean={last['delta_mean']:.4f} "
              f"max={last['delta_max']:.4f} (F1 anchor ~0.026); "
              f"tip-box straddler census: "
              f"{[(r['cls'], round(r['ratio'], 2)) for r in strad]}",
              flush=True)
    return cells, ring, strad


def leg_conf():
    mc, wc = wb30.build_conf("medium")
    m_inf = M_DYING["CONF"]
    ws = NewtonWorkspace(mc, wc, wb30.ALPHA, 1.0, 1.4, (0.25, 0.0), True,
                         kutta_estimator="pressure")
    ws.set_mach(m_inf)
    d = np.load(B30_OUT / "g2_conf_m0_83.npz", allow_pickle=True)
    phi, gamma = d["phi"], d["gamma"]
    cents = mc.nodes[mc.elements].mean(axis=1)
    cls = classify_conf(mc, wc, cents)
    lim, flr = wb30.clamp_masks_conf(ws, phi, gamma, m_inf,
                                     upwind_c=UP["upwind_c"],
                                     m_crit=UP["m_crit"],
                                     m_cap=UP["m_cap"],
                                     rho_floor=UP["rho_floor"])
    m2 = wb30.mach2_conf(ws, phi, gamma, m_inf)
    cells = clamp_cell_rows("CONF", dict(limited=lim, floored=flr),
                            cls, m2, cents)
    gprof = gamma_profile_conf(wc, gamma)
    tip_most = gprof[-1]
    print(f"  [CONF] clamps: {int(lim.sum())} lim + {int(flr.sum())} flr; "
          f"owners: {[(r['cls'], r['count']) for r in ownership_summary(cells)]}",
          flush=True)
    print(f"  [CONF] tip-most station z={tip_most['z']:.4f} "
          f"Gamma_last={tip_most['gamma']:.5f} "
          f"(taper F={tip_most['taper_f']:.3f} -> "
          f"Gamma_eff={tip_most['gamma_eff']:.5f})", flush=True)
    return cells, gprof


def _write(rows):
    keys = ["section", "leg", "clamp", "cls", "count", "elem", "x", "y",
            "z", "dist_tip_edge", "dist_tip_te", "mach", "q_lo", "q_hi",
            "z_mid", "n_elem", "delta_mean", "delta_max", "m2_main_max",
            "m2_side_max", "ratio", "station", "gamma", "taper_f",
            "gamma_eff"]
    with open(OUT / "g1_tip_atlas.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT / 'g1_tip_atlas.csv'}")


def _fig(ring, gprof, owners, cells):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    ax = axes[0]
    if ring:
        z = [r["z_mid"] for r in ring]
        ax.semilogy(z, [r["delta_mean"] for r in ring], "o-",
                    label="ring |δ| mean")
        ax.semilogy(z, [r["delta_max"] for r in ring], "s--",
                    label="ring |δ| max")
        for r in cells:
            if r["leg"] == "LS":
                ax.axvline(r["z"], color="r", alpha=0.35, lw=0.8)
        ax.axhline(0.026, color="k", ls=":", lw=1.0,
                   label="F1 anchor |δ|≈0.026")
        ax.set_xlabel("z of ring mid")
        ax.set_title("LS termination-ring δ(q), last 15% span\n"
                     "(red = LS clamp-cell z)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    ax = axes[1]
    if gprof:
        z = [r["z"] for r in gprof]
        ax.plot(z, [r["gamma"] for r in gprof], "o-", ms=3,
                label="Γ(z) dying state")
        ax.plot(z, [r["gamma_eff"] for r in gprof], "s--", ms=3,
                label="F(z)·Γ(z) — vanish_smooth r_c=0.05·b")
        ax.axvline(B_SEMI, color="k", ls=":", lw=1.0, label="B_SEMI")
        ax.set_xlim(z[0] + 0.8 * (z[-1] - z[0]), B_SEMI * 1.005)
        ax.set_xlabel("station z (tip end)")
        ax.set_title("CONF tip Γ(z) + GB31.2 taper overlay")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    ax = axes[2]
    labels = sorted({(r["leg"], r["cls"]) for r in owners})
    counts = [sum(r["count"] for r in owners
                  if (r["leg"], r["cls"]) == k) for k in labels]
    ax.barh([f"{k[0]}\n{k[1]}" for k in labels], counts, color="tab:red")
    ax.set_title("dying-level clamp-cell ownership\n"
                 f"(tip edge z={B_SEMI}, tip TE {tuple(np.round(TIP_TE, 3))})")
    ax.set_xlabel("cells")
    ax.grid(alpha=0.3, axis="x")
    fig.suptitle("GB31.1 tip-termination atlas (B30 dying states)")
    fig.tight_layout()
    fig.savefig(OUT / "g1_tip_atlas.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT / 'g1_tip_atlas.png'}")


def main():
    print("=== GB31.1 LS leg (medium, m=0.7875 dying state) ===", flush=True)
    cells_ls, ring, strad = leg_ls()
    print("=== GB31.1 CONF leg (medium, m=0.83 dying state) ===", flush=True)
    cells_conf, gprof = leg_conf()
    cells = cells_ls + cells_conf
    owners = ownership_summary(cells)
    _write(cells + owners + ring + strad + gprof)
    _fig(ring, gprof, owners, cells)
    unnamed = [r for r in cells if r["cls"] in ("", None)]
    print(f"\nGB31.1 PASS check: {len(cells)} clamp cells classified, "
          f"{len(unnamed)} unclassified -> "
          f"{'PASS' if not unnamed else 'FAIL'}")


if __name__ == "__main__":
    main()
