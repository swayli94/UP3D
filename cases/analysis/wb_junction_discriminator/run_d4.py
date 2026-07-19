"""D4 -- canonical crease case: smooth cylinder + symmetric plate, NO wake,
h-ladder. The cleanest possible reproduction attempt of the junction class.

Geometry: FuselageParams with a short blunt body (constant-radius cylinder
with rounded ends) + the M6 symmetric wing section sticking out, meshed by
the SAME onera_m6_wingbody_mesh generator (so the junction construction --
and only the junction construction -- is shared with the real wing-body).
No wake level set, no Kutta, no lifting machinery: solves are single-valued
Laplace (incompressible) and solve_subsonic (M0.5), both at alpha = 0 with
an alpha = 3.06 cross-flow leg.

Pre-registered question (PRE_REGISTRATION.md): does max |q| / Mmax at the
crease GROW with refinement? For a 270 deg reentrant corner the continuous
potential-flow velocity is singular, |q| ~ r^(pi/theta - 1) = r^(-1/3), so
the discrete max should grow ~ h^(-1/3) (x1.26 per halving) -- the
"refinement-worsening" signature is then the CONTINUOUS solution's own
character, not a discretization bug. A filleted corner (bounded solution)
must instead saturate. The fitted log-log slope is recorded as evidence.

Run:  python cases/analysis/wb_junction_discriminator/run_d4.py [--mesh-only]
Artifacts: meshes/d4_*.msh (gitignored), results/d4_meshes.csv,
           results/d4_summary.csv, results/d4_crease.png
"""

import argparse
import csv
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb_common import FUS, OUT, REPO_ROOT, load_mesh

MESH_DIR = OUT.parent / "meshes"
MESH_DIR.mkdir(exist_ok=True)

#: Canonical blunt body: 2.5 long, short rounded nose/tail, r = 0.15 (the
#: same r_f, so the junction station z_junc = r_f matches the real case).
from pyfp3d.meshgen.fuselage import FuselageParams
FUS_CANON = FuselageParams(r_f=0.15, length=2.5, x_center=0.403,
                           l_nose=0.3, l_tail=0.3, r_tail=0.05)

R_FAR = 6.0
LEVELS = {"h04": 0.04, "h02": 0.02, "h01": 0.01}
ALPHAS = (0.0, 3.06)
M_INF = 0.5


def gen_meshes():
    from pyfp3d.mesh.reader import write_mesh
    from pyfp3d.meshgen.wingbody import onera_m6_wingbody_mesh

    paths = {}
    stats = []
    for tag, h in LEVELS.items():
        mp = MESH_DIR / f"d4_{tag}.msh"
        paths[tag] = mp
        if mp.exists():
            print(f"[mesh {tag}] exists, skipped", flush=True)
            continue
        print(f"[mesh {tag}] h_wall={h} ...", flush=True)
        mesh = onera_m6_wingbody_mesh(
            h_wall=h, fuselage=FUS_CANON, r_far=R_FAR, name=f"d4_{tag}",
            h_far=25.0 * h, h_wake=3.0 * h, h_edge=0.5 * h,
            h_tip=0.25 * h, h_junction=0.5 * h,
            h_body=2.0 * h, h_body_tip=0.25 * h, tip_cap="round")
        write_mesh(mesh, str(mp))
        print(f"[mesh {tag}] {len(mesh.elements)} tets", flush=True)
    for tag in LEVELS:
        m = load_mesh(paths[tag])
        stats.append((tag, LEVELS[tag], len(m.nodes), len(m.elements)))
    with open(OUT / "d4_meshes.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tag", "h_wall", "n_nodes", "n_tets"])
        w.writerows(stats)
    return paths


def solve_legs(mesh, tag, alpha):
    """Laplace + subsonic M0.5 at the given alpha (cached). Returns dict."""
    from pyfp3d.kernels.jacobian import PicardOperator
    from pyfp3d.physics.isentropic import mach_squared_field
    from pyfp3d.solve.picard import solve_laplace, solve_subsonic

    cache = OUT / f"d4_{tag}_a{alpha:.2f}.npz"
    a = np.radians(alpha)
    far = np.unique(mesh.boundary_faces["farfield"])
    xyz = mesh.nodes[far]
    phi_inf = xyz[:, 0] * np.cos(a) + xyz[:, 1] * np.sin(a)   # u_inf = 1
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        print(f"  [d4 {tag} a={alpha}] CACHED", flush=True)
        return dict(phi_lap=d["phi_lap"], phi_sub=d["phi_sub"],
                    conv_sub=bool(d["conv_sub"]))
    op = PicardOperator(mesh.nodes, mesh.elements)
    t0 = time.perf_counter()
    r_lap = solve_laplace(mesh.nodes, mesh.elements, far, phi_inf)
    r_sub = solve_subsonic(mesh.nodes, mesh.elements, far, phi_inf, M_INF)
    print(f"  [d4 {tag} a={alpha}] lap cg={r_lap['n_cg_iterations']} "
          f"res={r_lap['residual_norm']:.1e} | sub conv={r_sub['converged']} "
          f"picard={r_sub['n_picard']} ({time.perf_counter()-t0:.0f}s)",
          flush=True)
    np.savez(cache, phi_lap=r_lap["phi"], phi_sub=r_sub["phi"],
             conv_sub=r_sub["converged"])
    return dict(phi_lap=r_lap["phi"], phi_sub=r_sub["phi"],
                conv_sub=bool(r_sub["converged"]))


def measure_crease(mesh, phi, m_inf):
    """Elementwise |grad phi| and M^2 maxima + localization of the max."""
    from pyfp3d.kernels.jacobian import PicardOperator
    from pyfp3d.meshgen.fuselage import radius_at
    from pyfp3d.meshgen.wingbody import junction_z
    from pyfp3d.physics.isentropic import mach_squared_field

    op = PicardOperator(mesh.nodes, mesh.elements)
    grad, q2 = op.velocities(phi)
    q = np.sqrt(np.maximum(q2, 0.0)).copy()          # u_inf = 1
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    i0 = int(np.argmax(q))
    z_junc = junction_z(FUS_CANON)
    c = cents[i0]
    rec = dict(qmax=float(q[i0]),
               qmax_x=float(c[0]), qmax_y=float(c[1]), qmax_z=float(c[2]),
               qmax_dist_fus=float(abs(np.hypot(c[1], c[2])
                                       - radius_at(FUS_CANON, float(c[0])))),
               qmax_z_minus_zjunc=float(c[2] - z_junc))
    if m_inf > 0.0:
        m2 = mach_squared_field(q2, m_inf, 1.4)
        rec["mmax"] = float(np.sqrt(m2.max()))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh-only", action="store_true")
    args = ap.parse_args()
    paths = gen_meshes()
    if args.mesh_only:
        return

    rows = []
    for tag, h in LEVELS.items():
        mesh = load_mesh(paths[tag])
        for alpha in ALPHAS:
            r = solve_legs(mesh, tag, alpha)
            m_lap = measure_crease(mesh, r["phi_lap"], 0.0)
            m_sub = measure_crease(mesh, r["phi_sub"], M_INF)
            rows.append(dict(tag=tag, h_wall=h, alpha=alpha,
                             n_tets=len(mesh.elements),
                             qmax_lap=m_lap["qmax"],
                             lap_x=m_lap["qmax_x"], lap_z=m_lap["qmax_z"],
                             lap_dist_fus=m_lap["qmax_dist_fus"],
                             lap_z_minus_zjunc=m_lap["qmax_z_minus_zjunc"],
                             qmax_sub=m_sub["qmax"], mmax_sub=m_sub["mmax"],
                             sub_x=m_sub["qmax_x"], sub_z=m_sub["qmax_z"],
                             sub_dist_fus=m_sub["qmax_dist_fus"],
                             sub_z_minus_zjunc=m_sub["qmax_z_minus_zjunc"],
                             conv_sub=r["conv_sub"]))
            print(f"  [{tag} a={alpha}] |q|max lap={m_lap['qmax']:.3f} "
                  f"sub={m_sub['qmax']:.3f} Mmax={m_sub['mmax']:.2f} at "
                  f"(x={m_sub['qmax_x']:.3f}, z={m_sub['qmax_z']:.3f}, "
                  f"|r-R|={m_sub['qmax_dist_fus']:.4f})", flush=True)

    with open(OUT / "d4_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT/'d4_summary.csv'}")
    _write_fig(rows)


def _write_fig(rows):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for alpha, col in ((0.0, "tab:blue"), (3.06, "tab:red")):
        rr = sorted((r for r in rows if r["alpha"] == alpha),
                    key=lambda r: r["h_wall"], reverse=True)
        h = np.array([r["h_wall"] for r in rr])
        ql = np.array([r["qmax_lap"] for r in rr])
        mm = np.array([r["mmax_sub"] for r in rr])
        ax[0].loglog(h, ql, "o-", color=col, label=f"a={alpha}")
        ax[1].loglog(h, mm, "o-", color=col, label=f"a={alpha}")
        # slope annotation between successive pairs
        for i in range(len(h) - 1):
            s = np.log(ql[i + 1] / ql[i]) / np.log(h[i + 1] / h[i])
            xm = np.sqrt(h[i] * h[i + 1])
            ax[0].annotate(f"slope {s:+.2f}", (xm, np.sqrt(ql[i]*ql[i+1])),
                           fontsize=8, color=col)
    # reference h^(-1/3) ray (the 270-deg-corner singularity exponent)
    href = np.array([min(LEVELS.values()), max(LEVELS.values())])
    ax[0].loglog(href, 2.2 * (href / href[0]) ** (-1/3), "k--", lw=0.9,
                 label="h^(-1/3) reference")
    ax[0].set_xlabel("h_wall"); ax[0].set_ylabel("max |grad phi| (u_inf=1)")
    ax[0].set_title("Laplace: crease max vs h (singularity exponent)")
    ax[1].axhline(1.0, color="0.5", ls="--", lw=0.8)
    ax[1].set_xlabel("h_wall"); ax[1].set_ylabel("Mmax (M0.5 subsonic)")
    ax[1].set_title("compressible: pocket strength vs h")
    for r in rows:
        ax[2].plot(r["sub_x"], r["sub_z"], "o", ms=9, fillstyle="none",
                   mew=1.5, color="tab:red" if r["alpha"] else "tab:blue",
                   label=f"{r['tag']} a={r['alpha']}")
    ax[2].axvspan(0.0, 0.806, color="tab:green", alpha=0.10,
                  label="plate chord (crease)")
    ax[2].axhline(0.15, color="0.5", ls="--", lw=0.8, label="z_junc")
    ax[2].set_xlabel("x of max element"); ax[2].set_ylabel("z of max element")
    ax[2].set_title("where is the max?")
    for a in ax:
        a.grid(alpha=0.3, which="both")
        a.legend(fontsize=7)
    fig.suptitle("D4 canonical crease: cylinder + symmetric plate, no wake")
    fig.tight_layout()
    fig.savefig(OUT / "d4_crease.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT/'d4_crease.png'}")


if __name__ == "__main__":
    main()
