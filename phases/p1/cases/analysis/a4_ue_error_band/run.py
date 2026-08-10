"""A4 — wall edge-velocity (u_e) error-band study (Track V input-quality prerequisite).

Purpose (from the 2026-07-20 wing-body/Track-V review §3.3 risk 1 + the audit
20260722-0335): a boundary-layer method (Track V IBL3) consumes the inviscid
WALL edge velocity u_e = |grad_tangential phi| and its arc derivative du_e/ds as
its input. The P1 wall-gradient recovery has an O(h) error that is the SAME
class as the V6 floor / the G1.6 re-attribution (intrinsic P1 max-norm
capability). If that input error band is not quantified FIRST, the V1/V3 gates
would mix "viscous-model error" with "inviscid-input error" and the M6/wing-body
VII comparison would be un-attributable.

This study measures the u_e error band on the two cases with an ANALYTIC
surface-speed ground truth (so the number is exact, not cross-model):
  * cylinder 2.5-D (M0): u_e = 2 U sin(theta),  Cp = 1 - 4 sin^2(theta)
  * sphere shell (M0):   u_e = 1.5 U sin(theta), Cp = 1 - 2.25 sin^2(theta)
both under the exact far-field Dirichlet potential + natural (zero-flux) wall BC
(the same setup tests/test_m0_cylinder + test_laplace_sphere use), for the two
committed mesh levels and BOTH wall-recovery schemes (linear
wall_tangential_gradient vs quadratic wall_tangential_gradient_quadratic).

It also records the SHARP-TE structural constraint on the airfoil geometry IBL
actually targets (NACA0012): the quadratic recovery RAISES on the sub-6-deg TE
wedge, so at a sharp TE -- the single most BL-sensitive station -- u_e is
available only from the linear/smoothed path, i.e. the smooth-wall band below
is the OPTIMISTIC bound there.

Outputs (committed):
  results/ue_bands.csv      -- u_e abs/rel error by region x level x scheme
  results/te_constraint.csv -- NACA0012 quadratic-recovery availability
  results/ue_error_band.png -- u_e(theta) numeric-vs-exact + error panels

Run: python cases/analysis/a4_ue_error_band/run.py   (~30 s; cheap M0 solves)
"""
import os
import numpy as np

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.post.surface import (
    wall_tangential_gradient,
    wall_tangential_gradient_quadratic,
    wall_crease_angles,
)
from pyfp3d.solve.picard import solve_laplace

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
MESHES = os.path.join(HERE, "..", "..", "meshes")

# ---- region bands (theta from the +x stagnation point, degrees) -------------
# LE/stagnation: u_e -> 0, du_e/ds MAX (the IBL initial-condition band);
# shoulder/peak: u_e is largest; overall = every wall node.
BANDS = {
    "LE_stag_0-20_160-180": lambda th: (th <= 20.0) | (th >= 160.0),
    "shoulder_peak_70-110": lambda th: (th >= 70.0) & (th <= 110.0),
    "overall": lambda th: np.ones_like(th, dtype=bool),
}


def _solve_exact_dirichlet(mesh, phi_exact):
    ff = np.unique(mesh.boundary_faces["farfield"])
    res = solve_laplace(mesh.nodes, mesh.elements, ff, phi_exact[ff],
                        rtol=1e-11, maxiter=3000)
    return res["phi"], res["residual_norm"]


def _ue_case(name, mesh_path, phi_exact_fn, ue_exact_fn, coef):
    mesh = read_mesh(mesh_path)
    nodes = mesh.nodes
    wall_faces = mesh.boundary_faces["wall"]
    wall = np.unique(wall_faces)

    phi_exact = phi_exact_fn(nodes)
    phi, resnorm = _solve_exact_dirichlet(mesh, phi_exact)

    r = np.linalg.norm(nodes[:, :ndim_of(name)], axis=1)
    cos_t = np.clip(nodes[wall, 0] / r[wall], -1.0, 1.0)
    theta = np.degrees(np.arccos(cos_t))
    ue_exact = ue_exact_fn(nodes, r)[wall]

    rows = []
    per_scheme_ue = {}
    for scheme, fn in (("linear", wall_tangential_gradient),
                       ("quadratic", wall_tangential_gradient_quadratic)):
        try:
            grad = fn(nodes, wall_faces, phi)
        except ValueError as e:
            per_scheme_ue[scheme] = None
            rows.append(dict(case=name, level=level_of(mesh_path), scheme=scheme,
                             region="ALL", n=0, abs_max=np.nan, abs_rms=np.nan,
                             rel_max_pct=np.nan, note=f"RAISED: {e}"))
            continue
        ue = np.linalg.norm(grad[wall], axis=1)
        per_scheme_ue[scheme] = (theta, ue, ue_exact)
        err = np.abs(ue - ue_exact)
        for region, mask_fn in BANDS.items():
            m = mask_fn(theta)
            if not np.any(m):
                continue
            # relative error only where u_e is non-tiny (peak band): guards the
            # stagnation 0/0; overall rel uses the max |u_e| as the scale.
            scale = max(np.abs(ue_exact[m]).max(), 1e-12)
            rows.append(dict(
                case=name, level=level_of(mesh_path), scheme=scheme,
                region=region, n=int(m.sum()),
                abs_max=float(err[m].max()), abs_rms=float(np.sqrt((err[m]**2).mean())),
                rel_max_pct=float(100.0 * err[m].max() / scale), note=""))
    return rows, per_scheme_ue, resnorm


def ndim_of(name):
    return 2 if "cylinder" in name else 3


def level_of(path):
    return os.path.splitext(os.path.basename(path))[0]


# ---- exact fields -----------------------------------------------------------
def cyl_phi(nodes, a=1.0):
    r2 = nodes[:, 0] ** 2 + nodes[:, 1] ** 2
    return nodes[:, 0] * (1.0 + a ** 2 / r2)


def cyl_ue(nodes, r):
    # surface speed 2 U sin(theta); on the wall sin = |y|/r
    sin_t = np.sqrt(np.clip(1.0 - (nodes[:, 0] / r) ** 2, 0.0, 1.0))
    return 2.0 * sin_t


def sph_phi(nodes, a=1.0):
    r = np.linalg.norm(nodes, axis=1)
    return nodes[:, 0] * (1.0 + 0.5 * a ** 3 / r ** 3)


def sph_ue(nodes, r):
    sin_t = np.sqrt(np.clip(1.0 - (nodes[:, 0] / r) ** 2, 0.0, 1.0))
    return 1.5 * sin_t


def main():
    os.makedirs(RESULTS, exist_ok=True)
    all_rows = []
    plot_data = {}
    cases = [
        ("cylinder", "cylinder_2.5d", cyl_phi, cyl_ue, 2.0),
        ("sphere", "sphere_shell", sph_phi, sph_ue, 1.5),
    ]
    for name, folder, phi_fn, ue_fn, coef in cases:
        for level in ("coarse", "medium"):
            mesh_path = os.path.join(MESHES, folder, f"{level}.msh")
            if not os.path.exists(mesh_path):
                continue
            rows, per_scheme, resnorm = _ue_case(name, mesh_path, phi_fn, ue_fn, coef)
            all_rows.extend(rows)
            plot_data[(name, level)] = per_scheme
            r_over = [r for r in rows if r["region"] == "overall" and r["scheme"] == "quadratic"]
            tag = f"{name} {level}: res={resnorm:.1e}"
            if r_over:
                tag += (f"  quad overall abs_max={r_over[0]['abs_max']:.4f}"
                        f" rel_max={r_over[0]['rel_max_pct']:.1f}%")
            print(tag)

    # ---- NACA0012 sharp-TE structural constraint ---------------------------
    # IBL's most sensitive station is the TE. The quadratic recovery's
    # _wall_vertex_normals guard RAISES when a vertex's incident face normals
    # nearly cancel (|sum area*n|/sum area < 0.05 ~ a sub-6-deg wedge). Measure
    # whether NACA0012's TE trips it, and the actual TE wedge angle, so the
    # constraint is stated with a threshold, not as a blanket claim.
    te_rows = []
    for lvl in ("coarse", "medium"):
        naca = os.path.join(MESHES, "naca0012_2.5d", f"{lvl}.msh")
        if not os.path.exists(naca):
            continue
        mesh = read_mesh(naca)
        wf = mesh.boundary_faces["wall"]
        phi, resnorm = _solve_exact_dirichlet(mesh, mesh.nodes[:, 0].copy())
        ang, _ = wall_crease_angles(mesh.nodes, mesh.elements, wf)
        # turning angle across the TE edge ~ max; wedge = 180 - turning.
        te_turning = float(np.nanmax(ang))
        te_wedge = 180.0 - te_turning
        lin_ok = quad_ok = True
        lin_note = quad_note = "ok"
        try:
            wall_tangential_gradient(mesh.nodes, wf, phi)
        except ValueError as e:
            lin_ok, lin_note = False, str(e)[:50]
        try:
            wall_tangential_gradient_quadratic(mesh.nodes, wf, phi)
        except ValueError as e:
            quad_ok, quad_note = False, str(e)[:50]
        te_rows.append(dict(case="naca0012", level=lvl,
                            te_wedge_deg=round(te_wedge, 2),
                            linear_available=lin_ok, quadratic_available=quad_ok,
                            linear_note=lin_note, quadratic_note=quad_note))
        print(f"naca0012 {lvl}: TE wedge={te_wedge:.1f} deg  "
              f"linear_avail={lin_ok} quadratic_avail={quad_ok} ({quad_note})")

    # ---- write CSVs --------------------------------------------------------
    _write_csv(os.path.join(RESULTS, "ue_bands.csv"),
               ["case", "level", "scheme", "region", "n",
                "abs_max", "abs_rms", "rel_max_pct", "note"], all_rows)
    _write_csv(os.path.join(RESULTS, "te_constraint.csv"),
               ["case", "level", "te_wedge_deg", "linear_available",
                "quadratic_available", "linear_note", "quadratic_note"], te_rows)

    _plot(plot_data)
    print(f"\nwrote {RESULTS}/ue_bands.csv, te_constraint.csv, ue_error_band.png")


def _write_csv(path, cols, rows):
    with open(path, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r.get(c, "")) for c in cols) + "\n")


def _plot(plot_data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keys = [k for k in plot_data if plot_data[k]]
    fig, axes = plt.subplots(2, len(keys), figsize=(4.2 * len(keys), 7))
    if len(keys) == 1:
        axes = axes.reshape(2, 1)
    for j, key in enumerate(keys):
        name, level = key
        per = plot_data[key]
        ax_u, ax_e = axes[0, j], axes[1, j]
        for scheme, style in (("linear", "C0."), ("quadratic", "C1.")):
            d = per.get(scheme)
            if d is None:
                continue
            th, ue, ue_ex = d
            o = np.argsort(th)
            ax_u.plot(th[o], ue[o], style, ms=2.5, alpha=0.5, label=scheme)
            ax_e.semilogy(th[o], np.abs(ue - ue_ex)[o] + 1e-12, style, ms=2.5,
                          alpha=0.5, label=scheme)
        # exact
        d = per.get("quadratic") or per.get("linear")
        th, _, ue_ex = d
        o = np.argsort(th)
        ax_u.plot(th[o], ue_ex[o], "k-", lw=1.5, label="exact")
        ax_u.set_title(f"{name} {level}: u_e(theta)")
        ax_u.set_xlabel("theta (deg)"); ax_u.set_ylabel("u_e / U")
        ax_u.legend(fontsize=8); ax_u.grid(alpha=0.3)
        ax_e.set_title(f"{name} {level}: |u_e error|")
        ax_e.set_xlabel("theta (deg)"); ax_e.set_ylabel("|u_e - u_e,exact|")
        ax_e.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "ue_error_band.png"), dpi=140,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
