"""
P11 demo -- curved wall-adjacent elements: the measured NEGATIVE, and the
re-attribution of G1.6 (roadmap track_p.md Sec P11; design.md Sec 5.1.3).

What this shows, honestly:
  1. G11.3 (PASS): the machinery is sound -- a planar wall + planar
     projection produces a stiffness delta that is EXACTLY zero (curved and
     flat share one code path), and the flat-geometry quadrature stiffness
     matches the independent P1 assembly to machine precision (locked in
     tests/test_p11_curved_walls.py).
  2. G11.4 (RECORDED): the curved layer really does move the domain -- edge
     midpoints shift O(h^2) onto the sphere, the chordal sliver volume is
     removed, and the result is quadrature-rule-insensitive (deg2 vs deg3).
  3. G11.1 (NOT MET, XFAIL): with a geometrically correct curved wall layer
     the medium-sphere Cp error moves 11.56% -> 11.33% -- the SAME value as
     the G1.4 boundary-data oracle ceiling, and nowhere near the < 2% gate.
     Mechanism (pre-registered risk, fired): the mapped-P1 basis on
     quadratic geometry loses exact linear reproduction at O(h) -- the same
     order as the facet-normal error the curving removes.
  4. G11.2 (NEGATIVE, XFAIL -- and the premise refuted): the clean h_min
     sweep REPLICATES the P1-era order collapse (0.88/0.56/0.42, now with a
     committed script), and curving does NOT restore it (0.80/0.50/0.39).
     The collapse was never geometric:
     - E6 control: an icosphere-extruded structured shell with the SAME
       flat-facet wall converges at ~2nd order (1.86/1.98) and reaches
       2.1% max Cp at h~0.027;
     - E8 discriminator: at fixed h_min = 0.03, refining ONLY the far mesh
       (h_max 3.0 -> 1.0) drops the wall phi error 3.2x and moves its argmax
       from r = 1.53 (the coarsening transition zone) to the wall; the
       medium -> far-refined-h03 order is ~1.9. The "decreasing order" was
       the fixed-bulk-mesh pollution floor of a single-variable h_min sweep.
  => G1.6's 11.6% at medium is essentially the intrinsic P1-field capability
     at h = 0.08 (structured control interpolates to ~11%), with the wall
     geometric crime contributing only ~0.2 pp. See the P11 close-out for
     the route fork (isoparametric P2 layer vs Option C re-spec = user).

Standalone + self-checking:  python cases/demo/p11_curved_walls/run_demo.py
Outputs: cases/demo/p11_curved_walls/results/{*.png, *.csv}
Sweep meshes (h05/h03/h02/h03_far10.msh) are LOCAL gitignored caches under
cases/meshes/sphere_shell/, regenerated on demand (~8 min); with them in
place the demo runs in ~5 min (16 threads). Exit code 0 iff every non-XFAIL
check passes.
"""

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, GRID, INK, MESH_DIR, MUTED, REPO_ROOT,
    S1_BLUE, S2_AQUA, S3_YELLOW, S5_VIOLET,
    CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402

from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.post.surface import wall_tangential_gradient_quadratic  # noqa: E402
from pyfp3d.solve.curved_wall import (  # noqa: E402
    assemble_curved_stiffness_delta, curved_volumes, curved_wall_geometry,
    plane_closest_point_normal,
)
from pyfp3d.solve.picard import solve_laplace  # noqa: E402
from pyfp3d.solve.wall_correction import sphere_closest_point_normal  # noqa: E402
from tests.mesh_utils import (  # noqa: E402
    generate_sphere_shell_mesh, generate_structured_cube_mesh, icosphere,
)
from tests.test_p11_curved_walls import _boundary_tris_on_plane  # noqa: E402

sys.path.insert(0, str(REPO_ROOT / "cases" / "meshes" / "sphere_shell"))
from generate_sphere_shell import generate_sphere_shell  # noqa: E402

OUT = Path(__file__).parent / "results"
SPHERE = MESH_DIR / "sphere_shell"
A = 1.0

# Local gitignored sweep-mesh caches (regenerated on demand, like the M6 .msh)
SWEEP_MESHES = {
    "h05": dict(h_min=0.05, h_max=3.0, r_out=20.0, dist_min=0.3, dist_max=10.0),
    "h03": dict(h_min=0.03, h_max=3.0, r_out=20.0, dist_min=0.3, dist_max=10.0),
    "h02": dict(h_min=0.02, h_max=3.0, r_out=20.0, dist_min=0.3, dist_max=10.0),
    "h03_far10": dict(h_min=0.03, h_max=1.0, r_out=20.0, dist_min=0.3, dist_max=10.0),
}


def ensure_mesh(name):
    path = SPHERE / f"{name}.msh"
    if not path.exists():
        t0 = time.time()
        generate_sphere_shell(path, **SWEEP_MESHES[name])
        print(f"  generated {path.name} in {time.time() - t0:.0f} s")
    return path


def run_sphere_case(mesh_path, curved, rule="deg2"):
    mesh = read_mesh(str(mesh_path))
    nodes, elements = mesh.nodes, mesh.elements
    wall_faces = mesh.boundary_faces["wall"]
    wall_nodes = np.unique(wall_faces)
    farfield_nodes = np.unique(mesh.boundary_faces["farfield"])
    r = np.linalg.norm(nodes, axis=1)
    phi_exact = nodes[:, 0] * (1.0 + 0.5 * A**3 / r**3)

    delta, geo = None, None
    if curved:
        geo = curved_wall_geometry(
            nodes, elements, wall_faces,
            lambda p: sphere_closest_point_normal(p, radius=A),
        )
        delta = assemble_curved_stiffness_delta(len(nodes), elements, geo, rule=rule)
    res = solve_laplace(
        nodes, elements, farfield_nodes, phi_exact[farfield_nodes],
        stiffness_delta=delta, rtol=1e-11, maxiter=3000,
    )
    phi = res["phi"]
    phi_err = np.abs(phi - phi_exact)

    grad_wall = wall_tangential_gradient_quadratic(nodes, wall_faces, phi)
    q2 = np.sum(grad_wall[wall_nodes] ** 2, axis=1)
    cp = 1.0 - q2
    cos_t = np.clip(nodes[wall_nodes, 0] / r[wall_nodes], -1.0, 1.0)
    cp_exact = 1.0 - 2.25 * (1.0 - cos_t**2)
    err = np.abs(cp - cp_exact)

    out = dict(
        n=len(nodes), cp_max=float(err.max()), cp_mean=float(err.mean()),
        phi_w=float(phi_err[wall_nodes].max()), phi_all=float(phi_err.max()),
        r_argmax=float(r[np.argmax(phi_err)]), res=res["residual_norm"],
        cos_t=cos_t, cp=cp, cp_ex=cp_exact,
    )
    if geo is not None:
        out.update(
            n_curved=len(geo["curved_tets"]), max_offset=geo["max_offset"],
            sliver=float(curved_volumes(geo["geom10_flat"]).sum()
                         - curved_volumes(geo["geom10"]).sum()),
        )
    return out


def demo_null_and_geometry(checks):
    print("\n[G11.3] planar null test + [G11.4] geometric fidelity")
    nodes, elements = generate_structured_cube_mesh(6)
    wall = _boundary_tris_on_plane(nodes, elements, axis=2, value=0.0)
    geo = curved_wall_geometry(
        nodes, elements, wall, plane_closest_point_normal((0, 0, 0), (0, 0, 1))
    )
    dA = assemble_curved_stiffness_delta(len(nodes), elements, geo)
    dmax = float(np.abs(dA.data).max()) if dA.nnz else 0.0
    checks.add("G11.3", "planar wall + planar projection -> dA", f"{dmax:.1e}",
               "== 0 exactly (bitwise)", dmax == 0.0)

    rows = []
    for level in ("coarse", "medium"):
        mesh = read_mesh(str(SPHERE / f"{level}.msh"))
        geo = curved_wall_geometry(
            mesh.nodes, mesh.elements, mesh.boundary_faces["wall"],
            lambda p: sphere_closest_point_normal(p, radius=A),
        )
        sliver = float(curved_volumes(geo["geom10_flat"]).sum()
                       - curved_volumes(geo["geom10"]).sum())
        rows.append((level, len(geo["curved_tets"]), f"{geo['max_offset']:.4e}",
                     f"{geo['mean_offset']:.4e}", f"{sliver:.5f}"))
        checks.add("G11.4", f"{level}: max midpoint offset (O(h^2))",
                   f"{geo['max_offset']:.2e}", "recorded", True)
        checks.add("G11.4", f"{level}: chordal sliver volume removed",
                   f"{sliver:.4f}", "> 0 (domain moved onto the sphere)", sliver > 0)
    write_csv(OUT, "g11_4_geometry.csv",
              "level,n_curved_tets,max_offset,mean_offset,sliver_volume_removed", rows)


def demo_g11_1_ab(checks):
    print("\n[G11.1] sphere coarse/medium: flat vs curved (the A/B)")
    rows = []
    results = {}
    for level in ("coarse", "medium"):
        f = run_sphere_case(SPHERE / f"{level}.msh", curved=False)
        c = run_sphere_case(SPHERE / f"{level}.msh", curved=True)
        c3 = run_sphere_case(SPHERE / f"{level}.msh", curved=True, rule="deg3")
        results[level] = (f, c)
        rows.append((level, f["n"], f"{f['cp_max']:.4f}", f"{c['cp_max']:.4f}",
                     f"{c3['cp_max']:.4f}", f"{f['cp_mean']:.4f}", f"{c['cp_mean']:.4f}",
                     c["n_curved"], f"{c['max_offset']:.3e}"))
        print(f"  {level}: flat {f['cp_max']:.4f} -> curved {c['cp_max']:.4f} "
              f"(deg3 {c3['cp_max']:.4f})")
        checks.add("G11.4", f"{level}: quadrature A/B |deg2 - deg3| on max Cp err",
                   f"{abs(c['cp_max'] - c3['cp_max']):.1e}", "< 1e-4 (recorded)",
                   abs(c["cp_max"] - c3["cp_max"]) < 1e-4)
    write_csv(OUT, "g11_1_ab.csv",
              "level,n_nodes,cp_max_flat,cp_max_curved,cp_max_curved_deg3,"
              "cp_mean_flat,cp_mean_curved,n_curved_tets,max_offset", rows)

    f, c = results["medium"]
    checks.add("G11.1", "medium curved max |Cp err|", f"{c['cp_max']:.4f}",
               "< 0.02 (the G1.6 criterion)", c["cp_max"] < 0.02, xfail=True,
               note="NOT MET; equals the G1.4 oracle ceiling 0.1133 -- "
                    "superparametric O(h) consistency replaces the O(h) normal error")
    checks.add("G11.1", "curved-vs-flat medium gain", f"{f['cp_max'] - c['cp_max']:.4f}",
               "recorded (~0.002 = the geometric-crime share)", True)

    # meridian overlay: the "nothing moved" visual
    apply_style()
    fig, ax = plt.subplots(figsize=(8, 5.5))
    th = np.linspace(0, 180, 300)
    ax.plot(th, 1.0 - 2.25 * np.sin(np.radians(th)) ** 2, "-", color=INK,
            lw=2, label="exact $1 - \\frac{9}{4}\\sin^2\\theta$")
    for (res, color, label) in ((f, S1_BLUE, "flat facets (11.56%)"),
                                (c, CRITICAL, "curved wall layer (11.33%)")):
        t = np.degrees(np.arccos(res["cos_t"]))
        o = np.argsort(t)
        ax.plot(t[o], res["cp"][o], ".", ms=2.5, alpha=0.45, color=color, label=label)
    ax.set_xlabel("theta (deg)"); ax.set_ylabel("$C_p$")
    ax.set_title("G11.1: medium sphere -- curved wall elements move almost nothing")
    ax.legend(loc="lower center"); ax.grid(True, color=GRID, alpha=0.5)
    finish(fig, OUT, "g11_1_meridian_ab.png")
    return results


def demo_g11_2_sweep(checks):
    print("\n[G11.2] clean h_min sweep 0.08/0.05/0.03/0.02: flat vs curved")
    hs = np.array([0.08, 0.05, 0.03, 0.02])
    paths = [SPHERE / "medium.msh", ensure_mesh("h05"), ensure_mesh("h03"),
             ensure_mesh("h02")]
    data = {}
    rows = []
    for label in ("flat", "curved"):
        runs = [run_sphere_case(p, curved=(label == "curved")) for p in paths]
        pw = np.array([r["phi_w"] for r in runs])
        cpm = np.array([r["cp_max"] for r in runs])
        orders = np.log(pw[:-1] / pw[1:]) / np.log(hs[:-1] / hs[1:])
        data[label] = (runs, pw, cpm, orders)
        print(f"  {label}: phi_w orders {np.round(orders, 2)}  cp_max {np.round(cpm, 4)}")
        for h, r in zip(hs, runs):
            rows.append((label, h, r["n"], f"{r['phi_w']:.6e}", f"{r['phi_all']:.6e}",
                         f"{r['cp_max']:.4f}", f"{r['cp_mean']:.4f}", f"{r['r_argmax']:.2f}"))
    write_csv(OUT, "g11_2_sweep.csv",
              "path,h_min,n_nodes,phi_err_wall_max,phi_err_all_max,cp_max,cp_mean,"
              "phi_err_argmax_radius", rows)

    orders_flat = data["flat"][3]
    orders_curved = data["curved"][3]
    baseline = np.array([0.88, 0.56, 0.42])
    checks.add("G11.2", "flat phi_w orders replicate the P1-era evidence",
               np.array2string(np.round(orders_flat, 2)),
               "== 0.88/0.56/0.42 (+-0.05)",
               bool(np.all(np.abs(orders_flat - baseline) < 0.05)),
               note="the root-cause sweep now exists as a committed script")
    checks.add("G11.2", "curved phi_w per-pair orders",
               np.array2string(np.round(orders_curved, 2)),
               ">= 1.5 each (pre-registered)",
               bool(np.all(orders_curved >= 1.5)), xfail=True,
               note="NEGATIVE: curving does not restore the order -- "
                    "the collapse is not geometric (see E6/E8)")
    return hs, data


def demo_e6_control(checks):
    print("\n[E6] structured-shell control (icosphere family, SAME flat-facet wall)")
    rows = []
    pw, cpm, hs = [], [], []
    for subdiv, n_layers in ((3, 24), (4, 48), (5, 96)):
        nodes, elements, wall_nodes, farfield_nodes = generate_sphere_shell_mesh(
            subdivisions=subdiv, n_layers=n_layers, r_inner=1.0, r_outer=25.0,
            grading=1.5)
        _, faces = icosphere(subdiv)
        r = np.linalg.norm(nodes, axis=1)
        phi_exact = nodes[:, 0] * (1.0 + 0.5 * A**3 / r**3)
        res = solve_laplace(nodes, elements, farfield_nodes,
                            phi_exact[farfield_nodes], rtol=1e-11, maxiter=3000)
        phi = res["phi"]
        grad_wall = wall_tangential_gradient_quadratic(
            nodes, faces.astype(np.int64), phi)
        q2 = np.sum(grad_wall[wall_nodes] ** 2, axis=1)
        cp = 1.0 - q2
        cos_t = np.clip(nodes[wall_nodes, 0] / r[wall_nodes], -1, 1)
        err = np.abs(cp - (1.0 - 2.25 * (1.0 - cos_t**2)))
        h = float(np.median(np.linalg.norm(
            nodes[faces[:, 0]] - nodes[faces[:, 1]], axis=1)))
        pw.append(np.abs(phi - phi_exact)[wall_nodes].max())
        cpm.append(err.max()); hs.append(h)
        rows.append((f"s{subdiv}", h, len(nodes), f"{pw[-1]:.6e}",
                     f"{cpm[-1]:.4f}", f"{err.mean():.4f}"))
        print(f"  s{subdiv}: h~{h:.3f} n={len(nodes)} phi_w={pw[-1]:.2e} cp_max={cpm[-1]:.4f}")
    pw = np.array(pw); hs_a = np.array(hs)
    orders = np.log(pw[:-1] / pw[1:]) / np.log(hs_a[:-1] / hs_a[1:])
    write_csv(OUT, "e6_ico_control.csv",
              "level,h_wall_median,n_nodes,phi_err_wall_max,cp_max,cp_mean", rows)
    checks.add("E6", "structured shell phi_w orders (flat facets!)",
               np.array2string(np.round(orders, 2)), ">= 1.6 each",
               bool(np.all(orders >= 1.6)),
               note="the same 'variational crime', ~2nd order anyway "
                    "=> the crime is not what limits the gmsh family")
    checks.add("E6", "structured shell max Cp at h~0.027", f"{cpm[-1]:.4f}",
               "recorded (~2.1% -- P1 intrinsic capability)", True)
    return hs, pw, cpm, orders


def demo_e8_floor(checks):
    print("\n[E8] bulk-floor discriminator: refine ONLY the far mesh at h_min=0.03")
    base = run_sphere_case(ensure_mesh("h03"), curved=False)
    far = run_sphere_case(ensure_mesh("h03_far10"), curved=False)
    medium = run_sphere_case(SPHERE / "medium.msh", curved=False)
    drop = base["phi_w"] / far["phi_w"]
    order_bulk = float(np.log(medium["phi_w"] / far["phi_w"]) / np.log(0.08 / 0.03))
    rows = [("h03 (h_max=3.0)", base["n"], f"{base['phi_w']:.6e}",
             f"{base['r_argmax']:.2f}", f"{base['cp_max']:.4f}", f"{base['cp_mean']:.4f}"),
            ("h03_far10 (h_max=1.0)", far["n"], f"{far['phi_w']:.6e}",
             f"{far['r_argmax']:.2f}", f"{far['cp_max']:.4f}", f"{far['cp_mean']:.4f}")]
    write_csv(OUT, "e8_bulk_floor.csv",
              "mesh,n_nodes,phi_err_wall_max,phi_err_argmax_radius,cp_max,cp_mean", rows)
    print(f"  phi_w {base['phi_w']:.2e} -> {far['phi_w']:.2e} ({drop:.2f}x), "
          f"argmax r {base['r_argmax']:.2f} -> {far['r_argmax']:.2f}, "
          f"order(medium->far-refined h03) = {order_bulk:.2f}")
    checks.add("E8", "phi_w drop from far-mesh-only refinement", f"{drop:.2f}x",
               ">= 2x (the sweep floor is bulk pollution)", drop >= 2.0)
    checks.add("E8", "h03 phi-err argmax radius (h_max=3.0)",
               f"{base['r_argmax']:.2f}", "> 1.2 (transition zone, not the wall)",
               base["r_argmax"] > 1.2)
    checks.add("E8", "order medium -> far-refined h03", f"{order_bulk:.2f}",
               ">= 1.6 (2nd order restored WITHOUT curved elements)",
               order_bulk >= 1.6)
    return base, far


def make_summary_figure(hs5, sweep, ico):
    hs_ico, pw_ico, cpm_ico, _ = ico
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))

    ax = axes[0]
    _, pw_f, cpm_f, _ = sweep["flat"]
    _, pw_c, cpm_c, _ = sweep["curved"]
    ax.loglog(hs5, cpm_f, "o-", color=S1_BLUE, label="gmsh family, flat")
    ax.loglog(hs5, cpm_c, "s--", color=CRITICAL, label="gmsh family, curved layer")
    ax.loglog(hs_ico, cpm_ico, "^-", color=S2_AQUA, label="structured shell, flat")
    ax.axhline(0.02, color=S3_YELLOW, lw=1.5)
    ax.text(0.021, 0.021, "G1.6 gate: 2%", color=INK, fontsize=9)
    ref = 0.115 * (np.array([0.11, 0.02]) / 0.08) ** 1.5
    ax.loglog([0.11, 0.02], ref, ":", color=MUTED)
    ax.text(0.045, 0.028, "$O(h^{1.5})$", color=MUTED, fontsize=9)
    ax.set_xlabel("wall h"); ax.set_ylabel("max |Cp err|")
    ax.set_title("Cp error: curving moves nothing;\nstructure + resolution move everything")
    ax.legend(fontsize=9); ax.grid(True, which="both", color=GRID, alpha=0.5)

    ax = axes[1]
    ax.loglog(hs5, pw_f, "o-", color=S1_BLUE, label="gmsh h_min sweep (h_max fixed): flat")
    ax.loglog(hs5, pw_c, "s--", color=CRITICAL, label="same, curved layer")
    ax.loglog(hs_ico, pw_ico, "^-", color=S2_AQUA, label="structured shell (all scales refine)")
    for i in range(3):
        o = np.log(pw_f[i] / pw_f[i + 1]) / np.log(hs5[i] / hs5[i + 1])
        ax.annotate(f"{o:.2f}", xy=(np.sqrt(hs5[i] * hs5[i + 1]),
                    np.sqrt(pw_f[i] * pw_f[i + 1])), color=S1_BLUE, fontsize=8)
    ax.set_xlabel("wall h"); ax.set_ylabel("max wall |phi err|")
    ax.set_title("the 'order collapse' is the fixed-bulk floor,\nnot a wall variational crime")
    ax.legend(fontsize=8, loc="upper left"); ax.grid(True, which="both", color=GRID, alpha=0.5)
    finish(fig, OUT, "p11_negative_and_reattribution.png")


def main():
    t0 = time.time()
    checks = CheckList("P11 curved wall-adjacent elements (sphere leg)")
    demo_null_and_geometry(checks)
    demo_g11_1_ab(checks)
    hs5, sweep = demo_g11_2_sweep(checks)
    ico = demo_e6_control(checks)
    demo_e8_floor(checks)
    make_summary_figure(hs5, sweep, ico)

    write_csv(OUT, "summary.csv", "metric,value", [
        ("medium_cp_max_flat", f"{sweep['flat'][0][0]['cp_max']:.4f}"),
        ("medium_cp_max_curved", "see g11_1_ab.csv"),
        ("g14_oracle_ceiling_2026_07_06", "0.1133"),
        ("flat_sweep_phi_orders", np.array2string(np.round(sweep["flat"][3], 2)).replace(",", ";")),
        ("curved_sweep_phi_orders", np.array2string(np.round(sweep["curved"][3], 2)).replace(",", ";")),
        ("ico_phi_orders", np.array2string(np.round(ico[3], 2)).replace(",", ";")),
        ("runtime_s", f"{time.time() - t0:.0f}"),
    ])
    code = checks.report(OUT)
    print(f"\ntotal {time.time() - t0:.0f} s")
    return code


if __name__ == "__main__":
    sys.exit(main())
