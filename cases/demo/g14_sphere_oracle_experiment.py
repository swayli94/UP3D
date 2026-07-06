"""
Gate G1.4 (formerly G1.2-a): sphere oracle ceiling experiment for the
Option A true-normal weak-flux correction (design.md §5.1, roadmap G1.4).

Assembles the correction RHS on the sphere wall facets from the EXACT
analytic gradient and EXACT true normals (no lagged iteration), solves
once, and measures the medium-mesh Cp error -- that number is the accuracy
ceiling of the first-order boundary-data correction. Two variants:

  t-form:    b_i = <N_i, grad(phi_exact) . t>,  t = n_facet - (n_facet.n)n
             (the Option A formula as specified)
  full-flux: b_i = <N_i, grad(phi_exact) . n_facet>  (7-point quadrature;
             the entire natural-BC consistency defect, an upper bound on
             what ANY boundary-data correction can restore)

MEASURED OUTCOME (2026-07-06, see artifacts/G1.4/summary.csv):
medium-mesh max |Cp err| = 0.1156 uncorrected, 0.1164 t-form, 0.1133
full-flux -- the ceiling is ~11.3%, far above the 5% DP1 branch point.
Same mechanism as the G1.3 cylinder pre-study: with body-fitted wall
vertices the exact solution's net flux through each facet is ~zero
(divergence theorem over the facet/true-surface sliver), so there is
almost no boundary-data defect to correct; the sphere's Cp error lives in
the domain perturbation + P1 approximation, not in the BC data.

Headless (roadmap §0.1): writes artifacts/G1.4/{cp_overlay.png,summary.csv}.
Runtime ~1 min. Usage:  python cases/demo/g14_sphere_oracle_experiment.py
"""

from pathlib import Path

import numpy as np

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.post.surface import wall_tangential_gradient_quadratic
from pyfp3d.solve.picard import solve_laplace
from pyfp3d.solve.wall_correction import (
    assemble_wall_flux_correction_rhs,
    sphere_closest_point_normal,
    wall_correction_geometry,
)

REPO = Path(__file__).resolve().parents[2]
A = 1.0  # sphere radius

# 7-point Gauss rule on the triangle (degree 5) for the full-flux variant
_A1, _B1 = 0.0597158717, 0.4701420641
_A2, _B2 = 0.7974269853, 0.1012865073
BARY7 = np.array([
    [1 / 3, 1 / 3, 1 / 3],
    [_A1, _B1, _B1], [_B1, _A1, _B1], [_B1, _B1, _A1],
    [_A2, _B2, _B2], [_B2, _A2, _B2], [_B2, _B2, _A2],
])
W7 = np.array([0.225] + [0.1323941527] * 3 + [0.1259391805] * 3)


def sphere_phi_exact(nodes):
    r = np.linalg.norm(nodes, axis=1)
    return nodes[:, 0] * (1.0 + 0.5 * A**3 / r**3)


def sphere_grad_exact(p):
    x, y, z = p[:, 0], p[:, 1], p[:, 2]
    r = np.linalg.norm(p, axis=1)
    c = 0.5 * A**3
    g = np.zeros_like(p)
    g[:, 0] = 1.0 + c / r**3 - 3 * c * x * x / r**5
    g[:, 1] = -3 * c * x * y / r**5
    g[:, 2] = -3 * c * x * z / r**5
    return g


def full_flux_rhs(mesh, wall_faces, geometry):
    tri = mesh.nodes[wall_faces]
    qp = np.einsum("qi,fik->fqk", BARY7, tri)
    g = sphere_grad_exact(qp.reshape(-1, 3)).reshape(qp.shape)
    flux = np.einsum("fqk,fk->fq", g, geometry["n_facet"])
    rhs = np.zeros(len(mesh.nodes))
    for v in range(3):
        contrib = geometry["area"][:, None] * W7[None, :] * BARY7[None, :, v] * flux
        np.add.at(rhs, wall_faces[:, v], contrib.sum(axis=1))
    return rhs


def cp_error(mesh, phi, wall_faces, wall_nodes, cp_exact):
    grad = wall_tangential_gradient_quadratic(mesh.nodes, wall_faces, phi)
    cp = 1.0 - np.sum(grad[wall_nodes] ** 2, axis=1)
    return cp, np.abs(cp - cp_exact)


def run_level(level):
    mesh = read_mesh(REPO / "cases" / "meshes" / "sphere_shell" / f"{level}.msh")
    wall_faces = mesh.boundary_faces["wall"]
    wall_nodes = np.unique(wall_faces)
    farfield_nodes = np.unique(mesh.boundary_faces["farfield"])
    phi_ex = sphere_phi_exact(mesh.nodes)
    cos_t = mesh.nodes[wall_nodes, 0] / np.linalg.norm(mesh.nodes[wall_nodes], axis=1)
    cp_exact = 1.0 - 2.25 * (1.0 - cos_t**2)

    geometry = wall_correction_geometry(
        mesh.nodes, mesh.elements, wall_faces, sphere_closest_point_normal
    )
    grad_qp = sphere_grad_exact(geometry["qp"].reshape(-1, 3)).reshape(
        geometry["qp"].shape
    )
    rhs_t = assemble_wall_flux_correction_rhs(
        len(mesh.nodes), wall_faces, geometry, grad_qp
    )
    rhs_full = full_flux_rhs(mesh, wall_faces, geometry)

    out = {"level": level, "wall_nodes": wall_nodes, "cos_t": cos_t,
           "cp_exact": cp_exact, "rhs_t_max": np.abs(rhs_t).max(),
           "rhs_full_max": np.abs(rhs_full).max()}
    for name, rhs in [("uncorrected", None), ("t_form", rhs_t),
                      ("full_flux", rhs_full)]:
        res = solve_laplace(mesh.nodes, mesh.elements, farfield_nodes,
                            phi_ex[farfield_nodes], body_source_rhs=rhs,
                            rtol=1e-11, maxiter=3000)
        cp, err = cp_error(mesh, res["phi"], wall_faces, wall_nodes, cp_exact)
        out[name] = {"cp": cp, "err": err,
                     "phi_wall_err": np.abs(res["phi"] - phi_ex)[wall_nodes].max()}
        print(f"  {level:7s} {name:11s} max|Cp err| = {err.max():.4f}  "
              f"mean = {err.mean():.4f}  wall phi err = "
              f"{out[name]['phi_wall_err']:.3e}")
    return out


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = REPO / "artifacts" / "G1.4"
    out_dir.mkdir(parents=True, exist_ok=True)

    results = [run_level(level) for level in ("coarse", "medium")]

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
    for ax, r in zip(axes, results):
        theta = np.degrees(np.arccos(np.clip(r["cos_t"], -1, 1)))
        order = np.argsort(theta)
        t_line = np.linspace(0, 180, 200)
        ax.plot(t_line, 1.0 - 2.25 * np.sin(np.radians(t_line)) ** 2,
                "k-", linewidth=2, label="exact")
        for name, style in [("uncorrected", "."), ("t_form", "x"),
                            ("full_flux", "+")]:
            ax.plot(theta[order], r[name]["cp"][order], style, markersize=3,
                    alpha=0.5, label=name)
        ax.set_xlabel("theta (deg)")
        ax.set_title(f"{r['level']}: oracle correction ceiling")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("Cp")
    fig.suptitle("G1.4 sphere oracle: boundary-data corrections cannot close "
                 "the gap (ceiling ~11.3% vs <2% target)")
    fig.savefig(out_dir / "cp_overlay.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    with open(out_dir / "summary.csv", "w") as f:
        f.write("level,variant,max_cp_err,mean_cp_err,max_wall_phi_err,max_rhs\n")
        for r in results:
            for name, rhs_key in [("uncorrected", None), ("t_form", "rhs_t_max"),
                                  ("full_flux", "rhs_full_max")]:
                rhs_val = "" if rhs_key is None else f"{r[rhs_key]:.3e}"
                f.write(f"{r['level']},{name},{r[name]['err'].max():.6e},"
                        f"{r[name]['err'].mean():.6e},"
                        f"{r[name]['phi_wall_err']:.6e},{rhs_val}\n")
    print(f"\nArtifacts written to {out_dir}")


if __name__ == "__main__":
    main()
