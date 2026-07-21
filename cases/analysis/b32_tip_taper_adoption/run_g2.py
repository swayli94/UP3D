"""GB32.2 adoption cost table + guardrails (pre-reg GB32.2).

Reads the demo's NEW conf_taper_* caches (post-adoption) against the
pre-adoption references (old conf_* caches where they exist; the
historical B9/CL_M05 anchor 0.2173 for the medium M0.50 leg, which the
pre-B32 demo never re-solved), and computes per-leg:

  cl_p cost %        (expected ~ -1% at M0.5, ~ -3% at 0.79; pre-reg
                      trigger: worse than -4% on any leg -> ◐ to user)
  corridor corrM     (guardrail <= 1.3, wb30.corridor_mmax)
  clamp counts       (honesty: reported, nonzero is not "clean")

Run AFTER the gated demo re-solve:
  python cases/analysis/b32_tip_taper_adoption/run_g2.py
Artifacts: results/g2_adoption_cost.csv (+ .png)
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_B31 = Path(__file__).resolve().parent.parent / "b31_tip_termination"
sys.path.insert(0, str(_B31))

import wb31  # noqa: E402,F401  (sys.path hook for the b30 helpers)
import wb30  # noqa: E402
from pyfp3d.solve.newton import NewtonWorkspace  # noqa: E402

OUT = Path(__file__).parent / "results"
OUT.mkdir(exist_ok=True)
DEMO = wb30.REPO_ROOT / "cases/demo/b18_wingbody_transonic/results"
UP = wb30.UPWIND_DEFAULT
# (label, level, target_tag, pre-adoption reference cl_p or None,
#  reference source)
LEGS = [("medium M0.50", "medium", "05", 0.2173, "B9/CL_M05 historical"),
        ("medium M0.65", "medium", "065", None, "old conf_medium_065"),
        ("medium M0.75", "medium", "075", None, "old conf_medium_075"),
        ("medium M0.79", "medium", "079", None, "old conf_medium_079"),
        ("medium climb", "medium", "084", None, "(new record, no ref)"),
        ("coarse M0.60", "coarse", "06", None, "old conf_coarse_06"),
        ("coarse M0.84", "coarse", "084", None, "old conf_coarse_084")]


def load(path):
    if not path.exists():
        return None
    return np.load(path, allow_pickle=True)


def guardrails(mc, wc, phi, gamma, m_inf):
    ws = NewtonWorkspace(mc, wc, wb30.ALPHA, 1.0, 1.4, (0.25, 0.0), True,
                         kutta_estimator="pressure")
    ws.set_mach(m_inf)
    lim, flr = wb30.clamp_masks_conf(ws, phi, gamma, m_inf,
                                     upwind_c=UP["upwind_c"],
                                     m_crit=UP["m_crit"],
                                     m_cap=UP["m_cap"],
                                     rho_floor=UP["rho_floor"])
    m2 = wb30.mach2_conf(ws, phi, gamma, m_inf)
    cents = mc.nodes[mc.elements].mean(axis=1)
    corr_m, _ = wb30.corridor_mmax(m2, cents)
    return dict(corr_mmax=corr_m, mmax_field=float(np.sqrt(m2.max())),
                mask_lim=int(lim.sum()), mask_flr=int(flr.sum()))


def main():
    meshes = {}
    rows = []
    for label, level, tag, ref_clp, ref_src in LEGS:
        new = load(DEMO / f"conf_taper_{level}_{tag}.npz")
        if new is None:
            print(f"  [skip] {label}: conf_taper_{level}_{tag}.npz missing")
            continue
        clp_new, m_new = float(new["clp"]), float(new["m"])
        if ref_clp is None and "old" in ref_src:
            old = load(DEMO / f"conf_{level}_{tag}.npz")
            ref_clp = float(old["clp"]) if old is not None else float("nan")
        cost = (100.0 * (clp_new - ref_clp) / ref_clp
                if ref_clp == ref_clp and ref_clp else float("nan"))
        if level not in meshes:
            meshes[level] = wb30.build_conf(level)
        mc, wc = meshes[level]
        g = guardrails(mc, wc, new["phi"], new["gamma"], m_new)
        rows.append(dict(leg=label, m=m_new, cl_p_new=round(clp_new, 6),
                         cl_p_ref=(round(ref_clp, 6)
                                   if ref_clp is not None
                                   and ref_clp == ref_clp else ""),
                         ref=ref_src, cost_pct=(round(cost, 2)
                                                if cost == cost else ""),
                         reached=bool(new["reached"]),
                         nlim=int(new["nlim"]), nflr=int(new["nflr"]),
                         mmax_cache=round(float(new["mmax"]), 4),
                         **{k: (round(v, 4) if isinstance(v, float) else v)
                            for k, v in g.items()}))
        print(f"  [{label}] m={m_new:.4g} cl_p {clp_new:.4f} vs "
              f"{ref_clp if ref_clp is not None and ref_clp==ref_clp else 'n/a'} "
              f"({ref_src}) = {cost:+.2f}% | corrM {g['corr_mmax']:.3f} "
              f"clamps {g['mask_lim']}+{g['mask_flr']}", flush=True)
    keys = ["leg", "m", "cl_p_new", "cl_p_ref", "ref", "cost_pct",
            "reached", "nlim", "nflr", "mmax_cache", "corr_mmax",
            "mmax_field", "mask_lim", "mask_flr"]
    with open(OUT / "g2_adoption_cost.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT / 'g2_adoption_cost.csv'}")
    _fig(rows)
    print("\n--- verdict hint (pre-reg GB32.2) ---")
    bad_cost = [r for r in rows if isinstance(r["cost_pct"], float)
                and r["cost_pct"] < -4.0]
    bad_corr = [r for r in rows if r["corr_mmax"] > 1.3]
    bad_clamp = [r for r in rows if r["mask_lim"] + r["mask_flr"] > 0]
    print(f"  legs beyond the -4% cost trigger: "
          f"{[r['leg'] for r in bad_cost] or 'none'}")
    print(f"  legs beyond corrM 1.3: {[r['leg'] for r in bad_corr] or 'none'}")
    print(f"  legs with nonzero clamps: "
          f"{[(r['leg'], r['mask_lim'], r['mask_flr']) for r in bad_clamp] or 'none'}")
    if not (bad_cost or bad_corr):
        print("  -> guardrails green; adoption cost table landed "
              "(demo checks.csv is the PASS/FAIL of record)")


def _fig(rows):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    labeled = [r for r in rows if isinstance(r["cost_pct"], float)]
    ax.bar(range(len(labeled)), [r["cost_pct"] for r in labeled],
           color="tab:orange")
    ax.axhline(-4.0, ls="--", color="tab:red", lw=1,
           label="-4% ◐ trigger")
    ax.axhline(0.0, color="k", lw=0.5)
    ax.set_xticks(range(len(labeled)))
    ax.set_xticklabels([r["leg"] for r in labeled], fontsize=8,
                       rotation=20)
    ax.set_ylabel("cl_p cost % (taper vs pre-adoption)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    ax = axes[1]
    ax.bar(range(len(rows)), [r["corr_mmax"] for r in rows],
           color="tab:blue")
    ax.axhline(1.3, ls="--", color="tab:red", lw=1, label="corrM guard 1.3")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([r["leg"] for r in rows], fontsize=8, rotation=20)
    ax.set_ylabel("corridor corrM")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.suptitle("GB32.2 adoption: cl_p cost + corridor guardrails")
    fig.tight_layout()
    fig.savefig(OUT / "g2_adoption_cost.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT / 'g2_adoption_cost.png'}")


if __name__ == "__main__":
    main()
