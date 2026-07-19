"""D3 -- junction fairing (fillet) A/B on the canonical geometry.

Pre-registered question (PRE_REGISTRATION.md): with the mesh policy fixed,
does filling the reentrant corner kill the crease maximum?

  * sharp crease (D4 baseline, same h)  -> max |q| / Mmax
  * sphere-pipe fairing radius 0.03/0.05 -> max |q| / Mmax

The fairing is the _fillet_junction sphere chain (OCC's own fillet fails on
the LE/TE cusps -- recorded in the wingbody.py docstring). A fairing that
kills the maximum => the crease itself was the cause AND a physically
legitimate mitigation exists (real aircraft have fairings).

Run:  python cases/analysis/wb_junction_discriminator/run_d3.py [--mesh-only]
Artifacts: meshes/d3_*.msh (gitignored), results/d3_summary.csv,
           results/d3_fillet.png
"""

import argparse
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from run_d4 import ALPHAS, FUS_CANON, M_INF, R_FAR, measure_crease, solve_legs
from wb_common import OUT, load_mesh

MESH_DIR = OUT.parent / "meshes"
MESH_DIR.mkdir(exist_ok=True)

#: (tag, h_wall, fillet radius). Baseline sharp meshes come from D4
#: (d4_h02 / d4_h04); the fairing legs reuse the same size laws.
H_REF = 0.02
LEGS = {"f030": 0.030, "f050": 0.050}
H_COARSE_CHECK = 0.04


def gen_meshes():
    from pyfp3d.mesh.reader import write_mesh
    from pyfp3d.meshgen.wingbody import onera_m6_wingbody_mesh

    paths, stats = {}, []
    for tag, rad in LEGS.items():
        for h in (H_REF,):
            mp = MESH_DIR / f"d3_{tag}_{h:.3f}.msh"
            paths[(tag, h)] = mp
            if mp.exists():
                print(f"[mesh {tag} h={h}] exists, skipped", flush=True)
                continue
            print(f"[mesh {tag} h={h}] fairing r={rad} ...", flush=True)
            mesh = onera_m6_wingbody_mesh(
                h_wall=h, fuselage=FUS_CANON, r_far=R_FAR,
                name=f"d3_{tag}", h_far=25.0 * h, h_wake=3.0 * h,
                h_edge=0.5 * h, h_tip=0.25 * h, h_junction=0.5 * h,
                h_body=2.0 * h, h_body_tip=0.25 * h, tip_cap="round",
                junction_fillet=rad)
            write_mesh(mesh, str(mp))
            print(f"[mesh {tag} h={h}] {len(mesh.elements)} tets", flush=True)
    for (tag, h), mp in paths.items():
        m = load_mesh(mp)
        stats.append((tag, h, LEGS[tag], len(m.nodes), len(m.elements)))
    with open(OUT / "d3_meshes.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tag", "h_wall", "fillet_radius", "n_nodes", "n_tets"])
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
    for tag, rad in LEGS.items():
        h = H_REF
        mesh = load_mesh(paths[(tag, h)])
        for alpha in ALPHAS:
            r = solve_legs(mesh, f"d3_{tag}", alpha)
            m_lap = measure_crease(mesh, r["phi_lap"], 0.0)
            m_sub = measure_crease(mesh, r["phi_sub"], M_INF)
            rows.append(dict(tag=tag, fillet_radius=rad, h_wall=h, alpha=alpha,
                             n_tets=len(mesh.elements),
                             qmax_lap=m_lap["qmax"],
                             lap_x=m_lap["qmax_x"], lap_z=m_lap["qmax_z"],
                             lap_dist_fus=m_lap["qmax_dist_fus"],
                             qmax_sub=m_sub["qmax"], mmax_sub=m_sub["mmax"],
                             sub_x=m_sub["qmax_x"], sub_z=m_sub["qmax_z"],
                             sub_dist_fus=m_sub["qmax_dist_fus"],
                             conv_sub=r["conv_sub"]))
            print(f"  [{tag} a={alpha}] |q|max lap={m_lap['qmax']:.3f} "
                  f"sub={m_sub['qmax']:.3f} Mmax={m_sub['mmax']:.2f} at "
                  f"(x={m_sub['qmax_x']:.3f}, z={m_sub['qmax_z']:.3f}, "
                  f"|r-R|={m_sub['qmax_dist_fus']:.4f})", flush=True)

    with open(OUT / "d3_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT/'d3_summary.csv'}")
    _write_fig(rows)


def _write_fig(rows):
    """Fairing legs vs the D4 sharp-crease baseline at the same h."""
    import csv as _csv
    base = {}
    p = OUT / "d4_summary.csv"
    if p.exists():
        with open(p) as fh:
            for r in _csv.DictReader(fh):
                if abs(float(r["h_wall"]) - H_REF) < 1e-9:
                    base[float(r["alpha"])] = r

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    labels = ["sharp (D4 h02)"] + [f"fairing {LEGS[t]}" for t in LEGS]
    for j, alpha in enumerate(ALPHAS):
        qlap = [float(base[alpha]["qmax_lap"]) if base.get(alpha) else np.nan]
        mmax = [float(base[alpha]["mmax_sub"]) if base.get(alpha) else np.nan]
        for tag in LEGS:
            r = next((r for r in rows
                      if r["tag"] == tag and r["alpha"] == alpha), None)
            qlap.append(r["qmax_lap"] if r else np.nan)
            mmax.append(r["mmax_sub"] if r else np.nan)
        x = np.arange(len(labels))
        a1 = ax[j]
        a1.bar(x - 0.2, qlap, width=0.4, label="max |grad phi| (Laplace)",
               color="tab:blue")
        a1.bar(x + 0.2, mmax, width=0.4, label="Mmax (M0.5)",
               color="tab:red")
        a1.set_xticks(x); a1.set_xticklabels(labels, fontsize=8, rotation=12)
        a1.set_title(f"alpha = {alpha}")
        a1.grid(alpha=0.3, axis="y")
        a1.legend(fontsize=8)
    fig.suptitle("D3: does the fairing kill the crease maximum? "
                 "(canonical geometry, h = 0.02)")
    fig.tight_layout()
    fig.savefig(OUT / "d3_fillet.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT/'d3_fillet.png'}")


if __name__ == "__main__":
    main()
