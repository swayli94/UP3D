"""D2 -- junction-only refinement A/B at FIXED h_wall = 0.015 (medium family).

Question (PRE_REGISTRATION.md): is the pocket a function of the JUNCTION-REGION
mesh, and which way does it move?

  * Mmax grows as h_junction shrinks  -> corner-singularity signature
  * Mmax shrinks                      -> under-resolution
  * Mmax flat                         -> not junction-mesh-local

Only h_junction varies; every other size law stays at the medium family's
(self-similar policy). The baseline medium (h_junction = 0.0075) mesh already
exists in cases/meshes/onera_m6_wingbody/ and anchors the ladder.

Run:  python cases/analysis/wb_junction_discriminator/run_d2.py [--mesh-only]
Artifacts: meshes/junc_*.msh (gitignored), results/d2_summary.csv,
           results/d2_meshes.csv, results/d2_junction.png
"""

import argparse
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb_common import (ALPHA_REF, FUS, LS_MESH_DIR, OUT, load_mesh, measure,
                       solve_ls)

HERE = OUT.parent
MESH_DIR = HERE / "meshes"
MESH_DIR.mkdir(exist_ok=True)

#: h_junction variants (in units of h_wall = 0.015). 0.5 = the baseline
#: medium mesh (reused, not regenerated).
H_WALL = 0.015
VARIANTS = {"j100": 1.00 * H_WALL,
            "j050": 0.50 * H_WALL,          # baseline medium (LS_MESH_DIR)
            "j025": 0.25 * H_WALL}


def gen_meshes():
    """Generate the non-baseline variants; return {tag: path}."""
    from pyfp3d.mesh.reader import write_mesh
    from pyfp3d.meshgen.wingbody import onera_m6_wingbody_mesh

    paths = {"j050": LS_MESH_DIR / "medium.msh"}
    stats = []
    for tag, h_j in VARIANTS.items():
        if tag == "j050":
            continue
        mp = MESH_DIR / f"{tag}.msh"
        paths[tag] = mp
        if mp.exists():
            print(f"[mesh {tag}] exists, skipped", flush=True)
            continue
        print(f"[mesh {tag}] h_junction={h_j:.5f} ...", flush=True)
        mesh = onera_m6_wingbody_mesh(
            h_wall=H_WALL, fuselage=FUS, tip_cap="round",
            h_far=200.0 * H_WALL, h_wake=3.0 * H_WALL, h_edge=0.5 * H_WALL,
            h_tip=0.25 * H_WALL, h_junction=h_j,
            h_body=2.0 * H_WALL, h_body_tip=0.25 * H_WALL)
        write_mesh(mesh, str(mp))
        print(f"[mesh {tag}] {len(mesh.elements)} tets", flush=True)
        stats.append((tag, h_j, len(mesh.nodes), len(mesh.elements)))
    for tag in ("j050",):
        m = load_mesh(paths[tag])
        stats.append((tag, VARIANTS[tag], len(m.nodes), len(m.elements)))
    with open(OUT / "d2_meshes.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tag", "h_junction", "n_nodes", "n_tets"])
        w.writerows(stats)
    return paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh-only", action="store_true")
    args = ap.parse_args()

    paths = gen_meshes()
    if args.mesh_only:
        return

    rows = []
    for tag in ("j100", "j050", "j025"):
        mp = paths[tag]
        if not mp.exists():
            print(f"[skip] {tag}: mesh absent")
            continue
        mesh = load_mesh(mp)
        print(f"=== D2 {tag}: {len(mesh.elements)} tets, "
              f"h_junction={VARIANTS[tag]:.5f} ===", flush=True)
        # the j050 baseline IS the D1 medium case (same mesh, alpha, solver
        # settings) -- seed its cache instead of re-solving (~15 min saved)
        if tag == "j050":
            seed_src = OUT / f"d1_medium_a{ALPHA_REF:.2f}_m0.50.npz"
            seed_dst = OUT / f"d2_j050_a{ALPHA_REF:.2f}_m0.50.npz"
            if seed_src.exists() and not seed_dst.exists():
                import shutil
                shutil.copyfile(seed_src, seed_dst)
                print("  [cache] seeded j050 from the D1 medium solve",
                      flush=True)
        rec, mvop = solve_ls(mesh, f"d2_{tag}", ALPHA_REF)
        m, top = measure(mesh, mvop, rec["phi_ext"], ALPHA_REF)
        rows.append(dict(tag=tag, h_junction=VARIANTS[tag],
                         n_tets=len(mesh.elements), converged=rec["conv"],
                         n_outer=rec["n"], wall_s=rec["wall_s"], **m))
        print(f"  [{tag}] Mmax={m['mmax']:.2f} at "
              f"(x={m['argmax_x']:.3f}, z={m['argmax_z']:.3f}, "
              f"|r-R|={m['argmax_dist_fus']:.4f}) n_sup={m['n_sup']} "
              f"cl_fus/wing={m['cl_fus_over_wing']:.3f}", flush=True)
        _write_topk(tag, top)

    keys = ["tag", "h_junction", "n_tets", "converged", "n_outer", "wall_s",
            "mmax", "m2_max", "argmax_elem", "argmax_x", "argmax_y", "argmax_z",
            "argmax_dist_fus", "argmax_z_minus_zjunc", "argmax_dist_te_junc",
            "n_sup", "sup_x_min", "sup_x_max", "sup_y_min", "sup_y_max",
            "sup_z_min", "sup_z_max", "cl_p", "cl_fus", "cl_fus_over_wing"]
    with open(OUT / "d2_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT/'d2_summary.csv'}")
    _write_fig(rows)


def _write_topk(tag, top):
    p = OUT / f"d2_topk_{tag}.csv"
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(top[0]))
        w.writeheader()
        w.writerows(top)


def _write_fig(rows):
    rr = sorted(rows, key=lambda r: r["h_junction"], reverse=True)
    h = [r["h_junction"] for r in rr]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
    ax[0].semilogx(h, [r["mmax"] for r in rr], "o-", color="tab:red")
    ax[0].axhline(1.0, color="0.5", ls="--", lw=0.8)
    ax[0].set_ylabel("Mmax"); ax[0].set_title("pocket strength vs h_junction")
    ax[1].semilogx(h, [r["n_sup"] for r in rr], "o-", color="tab:red")
    ax[1].set_ylabel("# supersonic elements"); ax[1].set_title("pocket extent")
    ax[2].semilogx(h, [r["cl_fus_over_wing"] for r in rr], "o-", color="tab:red")
    ax[2].set_ylabel("|cl_fus| / cl_p_wing"); ax[2].set_title("fuselage lift")
    for a in ax:
        a.set_xlabel("h_junction (m, fixed h_wall=0.015)")
        a.grid(alpha=0.3, which="both")
        a.invert_xaxis()
    fig.suptitle("D2: junction-only refinement at fixed medium family "
                 "(LS Picard M0.5 alpha=3.06)")
    fig.tight_layout()
    fig.savefig(OUT / "d2_junction.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT/'d2_junction.png'}")


if __name__ == "__main__":
    main()
