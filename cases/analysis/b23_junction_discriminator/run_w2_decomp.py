"""W2 decomposition: is the fuselage spurious lift pocket-band contamination?

GB9.4 anchor: |cl_fus|/cl_p_wing = 0.16 (conf) / 0.20 (LS) at medium; the LS
value GROWS with refinement. D1/D5 localized the supersonic pocket to the
wake-sheet inboard free edge along the fuselage waterline (z ~= z_junc,
x >~ 1.1, behind the wing). Question: how much of cl_fus comes from wall
triangles inside that band vs the rest of the fuselage skin?

Method: reuse the cached D1 solves (no re-solve). For each (level, alpha),
take the fuselage wall triangles, compute per-triangle Cp via the same
`_d11_wall_state` + `_cp_from_q2` core `wall_forces` uses, then integrate the
pressure lift (same `_pressure_force` core) over triangle subsets:

  * all          -- the whole fuselage (must reproduce d1 cl_fus exactly)
  * band         -- |z - z_junc| < bw  AND  x > x_band0 (pocket corridor)
  * outside      -- the complement of band
  * poles        -- nose/tail polar caps (GB9.6 suspects), x < x_nose+h or
                    x > x_tail-h band ends; reported for attribution

Bandwidth sweep bw in {0.5, 1, 2, 4} x h_body gives the sensitivity curve;
if cl_fus(outside) stays small across the sweep while cl_fus(band) carries
most of the total, W2 is pocket-band contamination (same root as W1).

Artifacts: results/w2_decomposition.csv, results/w2_decomposition.png
"""

import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb_common import (ALPHA_REF, LS_MESH_DIR, OUT, Z_JUNC, build_ls,
                       load_mesh)
from pyfp3d.meshgen.fuselage import FuselageParams
from pyfp3d.post.surface import _cp_from_q2, _pressure_force, planform_area
from pyfp3d.post.surface_ls import _d11_wall_state

M_INF = 0.5
FUS = FuselageParams()

# pocket corridor anchors (D1/D5): supersonic elements live at x in
# [1.13, 2.37], z within +-0.02 of z_junc; band starts behind the wing.
X_BAND0 = 1.0
H_BODY = 0.03          # medium h_body scale for the bandwidth ladder
BW_LADDER = (0.015, 0.03, 0.06, 0.12)

# polar-cap margins (GB9.6: Cp scatter max lives at the nose/tail poles)
X_NOSE = FUS.x_center - FUS.length / 2.0
X_TAIL = FUS.x_center + FUS.length / 2.0
POLE_MARGIN = 0.10

CASES = (("medium", (0.0, 1.0, 2.0, ALPHA_REF)),
         ("coarse", (0.0, ALPHA_REF)))


def fuselage_cl_parts(mesh, mvop, phi_ext, alpha_deg, s_ref):
    wall = np.asarray(mesh.boundary_faces["fuselage"], dtype=np.int64)
    q2, _, area, n_out = _d11_wall_state(mesh, mvop, phi_ext, wall, 1.0)
    cp = _cp_from_q2(q2, M_INF, 1.4)
    cents = mesh.nodes[wall].mean(axis=1)

    def cl_of(mask):
        _, cl, _ = _pressure_force(cp[mask], area[mask], n_out[mask],
                                   s_ref, alpha_deg)
        return cl

    allm = np.ones(len(wall), dtype=bool)
    rows = {}
    rows["all"] = cl_of(allm)
    # polar caps (nose/tail)
    pole = (cents[:, 0] < X_NOSE + POLE_MARGIN) | \
           (cents[:, 0] > X_TAIL - POLE_MARGIN)
    rows["poles"] = cl_of(pole)
    rows["no_poles"] = cl_of(~pole)
    # pocket corridor at each bandwidth
    for bw in BW_LADDER:
        band = (np.abs(cents[:, 2] - Z_JUNC) < bw) & (cents[:, 0] > X_BAND0)
        rows[f"band_bw{bw:g}"] = cl_of(band)
        rows[f"out_bw{bw:g}"] = cl_of(~band)
        rows[f"areafrac_bw{bw:g}"] = float(area[band].sum() / area.sum())
    return rows


def main():
    levels = sys.argv[1:] or [lv for lv, _ in CASES]
    alphas = dict(CASES)
    out_rows = []
    for level in levels:
        mp = LS_MESH_DIR / f"{level}.msh"
        if not mp.exists():
            print(f"[skip] {level}: mesh absent")
            continue
        mesh = load_mesh(mp)
        s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
        print(f"=== W2 {level}: {len(mesh.elements)} tets ===", flush=True)
        for a in alphas[level]:
            cache = OUT / f"d1_{level}_a{a:.2f}_m{M_INF:.2f}.npz"
            if not cache.exists():
                print(f"  [skip] a={a}: no cached solve ({cache.name})")
                continue
            d = np.load(cache, allow_pickle=True)
            _, _, mvop = build_ls(mesh, a)
            parts = fuselage_cl_parts(mesh, mvop, d["phi_ext"], a, s_ref)
            rec = dict(level=level, alpha=a, **parts)
            out_rows.append(rec)
            print(f"  [a={a:>4}] cl_fus(all)={parts['all']:+.5f} "
                  f"band(bw=0.06)={parts['band_bw0.06']:+.5f} "
                  f"out={parts['out_bw0.06']:+.5f} "
                  f"poles={parts['poles']:+.5f}", flush=True)

    if not out_rows:
        print("nothing measured")
        sys.exit(0)
    keys = ["level", "alpha", "all", "poles", "no_poles"]
    for bw in BW_LADDER:
        keys += [f"band_bw{bw:g}", f"out_bw{bw:g}", f"areafrac_bw{bw:g}"]
    with open(OUT / "w2_decomposition.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {OUT/'w2_decomposition.csv'}")
    _write_fig(out_rows)


def _write_fig(rows):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    for level, col in (("coarse", "tab:blue"), ("medium", "tab:red")):
        rr = sorted((r for r in rows if r["level"] == level),
                    key=lambda r: r["alpha"])
        if not rr:
            continue
        a = [r["alpha"] for r in rr]
        ax[0].plot(a, [r["all"] for r in rr], "o-", color=col,
                   label=f"{level} all")
        ax[0].plot(a, [r["band_bw0.06"] for r in rr], "s--", color=col,
                   label=f"{level} band (|z-z_j|<0.06, x>1)")
        ax[0].plot(a, [r["out_bw0.06"] for r in rr], "^:", color=col,
                   label=f"{level} outside")
        ax[0].plot(a, [r["no_poles"] for r in rr], "d-.", color=col,
                   alpha=0.5, label=f"{level} minus poles")
    ax[0].axhline(0.0, color="0.5", lw=0.8)
    ax[0].set_xlabel("alpha (deg)")
    ax[0].set_ylabel("cl_fus contribution")
    ax[0].set_title("fuselage lift = pocket band + rest?")
    ax[0].grid(alpha=0.3)
    ax[0].legend(fontsize=7)

    # bandwidth sensitivity at the reference alpha (medium)
    rr = [r for r in rows if r["level"] == "medium"
          and abs(r["alpha"] - ALPHA_REF) < 1e-9]
    if rr:
        r = rr[0]
        bws = list(BW_LADDER)
        ax[1].plot(bws, [r[f"band_bw{b:g}"] for b in bws], "s-",
                   color="tab:red", label="band cl_fus")
        ax[1].plot(bws, [r[f"out_bw{b:g}"] for b in bws], "^-",
                   color="tab:green", label="outside cl_fus")
        ax[1].plot(bws, [r[f"areafrac_bw{b:g}"] for b in bws], "d:",
                   color="0.4", label="band area fraction")
        ax[1].axhline(r["all"], color="tab:red", ls="--", lw=0.8,
                      label=f"all = {r['all']:.4f}")
        ax[1].set_xscale("log")
        ax[1].set_xlabel("pocket-band half-width bw (around z_junc, x>1)")
        ax[1].set_title(f"bandwidth sensitivity (medium, a={ALPHA_REF})")
        ax[1].grid(alpha=0.3)
        ax[1].legend(fontsize=8)
    fig.suptitle("W2 decomposition: where does the spurious fuselage lift "
                 "come from? (cached D1 solves, M0.5)")
    fig.tight_layout()
    fig.savefig(OUT / "w2_decomposition.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT/'w2_decomposition.png'}")


if __name__ == "__main__":
    main()
