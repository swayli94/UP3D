"""
Track B / B9 gate GB9.6 (RECORDED, no pass/fail): the fuselage surface-Cp
discretization-error guardrail.

The fuselage is exactly the smooth-body geometry class where the flat-facet
natural-BC "variational crime" bites (G1.6: sphere Cp 11.6% at medium,
saturating ~3.6% under h-refinement). The wing lift gates were exonerated
only because lift is circulation-dominated -- an argument that does NOT
cover body SURFACE pressures. So before any B9 wing-body surface-Cp is
read, this quantifies the error on an ISOLATED body of revolution (the M2
fuselage, wing off), subsonic, non-lifting.

Error metric (honest statement of choice -- there is no analytic Cp for the
splined body):
  * PRIMARY = azimuthal Cp scatter. At alpha = 0 the exact solution is
    axisymmetric, so the exact Cp depends only on x. The per-x-bin standard
    deviation of the per-triangle Cp is therefore a REFERENCE-FREE, absolute
    readout of the discretization error -- no truth needed.
  * SECONDARY = h-sweep self-convergence: the binned meridian Cp(x) per
    level against the finest, max/median over the mid-body. Stated plainly
    as SELF-convergence (this error class partially HIDES from Richardson --
    G1.6 saturates -- so the primary metric is the headline).

h_body in {0.060, 0.030, 0.015} = the wing-body coarse/medium/fine fuselage
skin resolutions (h_body = 2 h_wall), so the recorded error speaks for the
resolutions B9 actually reads fuselage Cp at.

Run:  python cases/analysis/b9_fuselage_guardrail/run_guardrail.py
Artifacts: guardrail_cp_sweep.png, guardrail_summary.csv, checks.csv
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.meshgen.fuselage import FuselageParams, add_fuselage_solid, radius_at
from pyfp3d.meshgen.wing3d import MAC, _collect_3d
from pyfp3d.post.unified import wall_cp
from pyfp3d.solve.picard import solve_subsonic

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

M_INF = 0.5
R_FAR = 25.0 * MAC
H_BODY = {"coarse": 0.060, "medium": 0.030, "fine": 0.015}
FUS = FuselageParams()


def fuselage_alone_mesh(h_body: float, n_profile: int = 120):
    """Isolated body of revolution in a half-ball, wing OFF. Groups
    wall/farfield/symmetry (so the default post-processing extractors work).
    Skin at h_body, radius-driven at the two tips (the M2 body-sizing lesson),
    graded to h_far in the field."""
    import gmsh

    p = FUS
    h_far = 200.0 * h_body / 2.0     # keep h_far/h_body ~ the M2 far gradient
    xc = 0.5 * (p.x_nose_tip + p.x_tail_tip)
    r_far = R_FAR
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("fuselage_alone")
        occ = gmsh.model.occ
        body = add_fuselage_solid(occ, p, n_profile=n_profile)
        ball = occ.addSphere(xc, 0.0, 0.0, r_far)
        below = occ.addBox(xc - 2 * r_far, -2 * r_far, -2 * r_far,
                           4 * r_far, 4 * r_far, 2 * r_far)
        half, _ = occ.cut([(3, ball)], [(3, below)])
        fluid, _ = occ.cut(half, body)
        occ.synchronize()

        tol = 1e-3
        groups = {"wall": [], "farfield": [], "symmetry": []}
        for dim, tag in gmsh.model.getEntities(2):
            bb = gmsh.model.getBoundingBox(dim, tag)
            extent = max(bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2])
            if (bb[5] - bb[2]) < tol and abs(bb[5]) < tol:
                groups["symmetry"].append(tag)
            elif extent > 1.2 * r_far:
                groups["farfield"].append(tag)
            else:
                groups["wall"].append(tag)
        for g, tags in groups.items():
            assert tags, f"no '{g}' faces"

        field = gmsh.model.mesh.field
        f = field.add("Distance")
        field.setNumbers(f, "SurfacesList", groups["wall"])
        field.setNumber(f, "Sampling", 200)
        t = field.add("Threshold")
        field.setNumber(t, "InField", f)
        field.setNumber(t, "SizeMin", h_body)
        field.setNumber(t, "SizeMax", h_far)
        field.setNumber(t, "DistMin", 0.05)
        field.setNumber(t, "DistMax", 0.55 * r_far)
        # radius-driven tip refinement (facet angle ~ h / R)
        tips = []
        for x_tip, ramp in ((p.x_nose_tip, 0.5 * p.l_nose),
                            (p.x_tail_tip, p.l_tail + p.r_tail)):
            b = field.add("Ball")
            field.setNumber(b, "XCenter", x_tip)
            field.setNumber(b, "YCenter", 0.0)
            field.setNumber(b, "ZCenter", 0.0)
            field.setNumber(b, "Radius", 0.2 * p.r_f)
            field.setNumber(b, "Thickness", ramp)
            field.setNumber(b, "VIn", 0.25 * h_body)
            field.setNumber(b, "VOut", h_body)
            tips.append(b)
        tmin = field.add("Min")
        field.setNumbers(tmin, "FieldsList", tips)
        fmax = field.add("Max")
        field.setNumbers(fmax, "FieldsList", [t, tmin])
        fmin = field.add("Min")
        field.setNumbers(fmin, "FieldsList", [fmax])
        field.setAsBackgroundMesh(fmin)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.model.mesh.generate(3)
        mesh = _collect_3d(groups, name="fuselage_alone")
        return mesh, xc, r_far
    finally:
        gmsh.finalize()


def solve_level(level: str):
    mesh, xc, r_far = fuselage_alone_mesh(H_BODY[level])
    ff = np.unique(mesh.boundary_faces["farfield"])
    # freestream potential u_inf * x on the far field (alpha = 0, non-lifting)
    dvals = mesh.nodes[ff, 0]
    r = solve_subsonic(mesh.nodes, mesh.elements, ff, dvals, m_inf=M_INF)
    cp = wall_cp(mesh, phi=r["phi"], m_inf=M_INF, wall_tag="wall")
    return mesh, r, cp, xc


def azimuthal_scatter(cp, p, nbins=40):
    """Per-x-bin std/IQR of per-triangle Cp on the body (reference-free at
    alpha=0). Restricted to the mid-body (nose_end..body_end) so the tip
    caps' own resolution does not dominate the readout."""
    x = cp["x"]
    c = cp["cp"]
    lo, hi = p.x_nose_end, p.x_body_end
    sel = (x >= lo) & (x <= hi)
    bins = np.linspace(lo, hi, nbins + 1)
    idx = np.clip(np.digitize(x[sel], bins) - 1, 0, nbins - 1)
    xs, std, iqr, med = [], [], [], []
    for i in range(nbins):
        m = idx == i
        if m.sum() < 4:
            continue
        ci = c[sel][m]
        xs.append(0.5 * (bins[i] + bins[i + 1]))
        std.append(float(np.std(ci)))
        iqr.append(float(np.subtract(*np.percentile(ci, [75, 25]))))
        med.append(float(np.median(ci)))
    return (np.array(xs), np.array(std), np.array(iqr), np.array(med))


def meridian_cp(cp, p, nbins=60):
    """Binned meridian Cp(x) over the whole body (for the h-sweep delta)."""
    x, c = cp["x"], cp["cp"]
    bins = np.linspace(p.x_nose_tip, p.x_tail_tip, nbins + 1)
    idx = np.clip(np.digitize(x, bins) - 1, 0, nbins - 1)
    xm = 0.5 * (bins[:-1] + bins[1:])
    med = np.array([np.median(c[idx == i]) if (idx == i).any() else np.nan
                    for i in range(nbins)])
    return xm, med


def main():
    levels = list(H_BODY)
    data = {}
    for lv in levels:
        mesh, r, cp, xc = solve_level(lv)
        xs, std, iqr, med = azimuthal_scatter(cp, FUS)
        xm, mer = meridian_cp(cp, FUS)
        data[lv] = dict(n_tets=len(mesh.elements), converged=r["converged"],
                        xs=xs, std=std, iqr=iqr, med=med, xm=xm, mer=mer,
                        scatter_med=float(np.median(std)),
                        scatter_max=float(np.max(std)))
        print(f"[{lv}] n_tets={len(mesh.elements)} conv={r['converged']} "
              f"azimuthal-std med={np.median(std):.4f} max={np.max(std):.4f}",
              flush=True)

    # h-sweep self-convergence vs the finest level
    fine = data[levels[-1]]
    for lv in levels[:-1]:
        d = data[lv]
        merp = np.interp(fine["xm"], d["xm"], d["mer"])
        mid = (fine["xm"] >= FUS.x_nose_end) & (fine["xm"] <= FUS.x_body_end)
        delta = np.abs(merp[mid] - fine["mer"][mid])
        d["selfconv_max"] = float(np.nanmax(delta))
        d["selfconv_med"] = float(np.nanmedian(delta))

    # --- figures ---------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for lv, col in zip(levels, ("tab:red", "tab:blue", "tab:green")):
        d = data[lv]
        ax[0].plot(d["xm"], d["mer"], "-", color=col, lw=1.3,
                   label=f"{lv} (h_body={H_BODY[lv]}, {d['n_tets']} tets)")
        ax[1].plot(d["xs"], d["std"], "o-", color=col, ms=3, lw=1.0,
                   label=f"{lv} scatter med={d['scatter_med']:.3f}")
    ax[0].axvspan(FUS.x_nose_end, FUS.x_body_end, color="0.85", alpha=0.5)
    ax[0].set_xlabel("x"); ax[0].set_ylabel("meridian median Cp")
    ax[0].set_title("Fuselage-alone Cp(x), M0.5, alpha=0\n(grey = mid-body)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].set_xlabel("x"); ax[1].set_ylabel("azimuthal Cp std (per x-bin)")
    ax[1].set_title("PRIMARY metric: azimuthal Cp scatter\n"
                    "(reference-free discretization error at alpha=0)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.suptitle("GB9.6 fuselage surface-Cp guardrail (RECORDED, no pass/fail)")
    fig.tight_layout()
    fig.savefig(OUT / "guardrail_cp_sweep.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # --- summary CSV -----------------------------------------------------
    with open(OUT / "guardrail_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["level", "h_body", "n_tets", "converged",
                    "azimuthal_scatter_med", "azimuthal_scatter_max",
                    "selfconv_med_vs_fine", "selfconv_max_vs_fine"])
        for lv in levels:
            d = data[lv]
            w.writerow([lv, H_BODY[lv], d["n_tets"], d["converged"],
                        f"{d['scatter_med']:.5f}", f"{d['scatter_max']:.5f}",
                        f"{d.get('selfconv_med', float('nan')):.5f}",
                        f"{d.get('selfconv_max', float('nan')):.5f}"])

    # --- checks.csv (RECORDED only) --------------------------------------
    with open(OUT / "checks.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gate", "name", "value", "criterion", "status"])
        for lv in levels:
            d = data[lv]
            w.writerow(["GB9.6", f"{lv}_converged", d["converged"],
                        "solve converges", "PASS" if d["converged"] else "FAIL"])
        w.writerow(["GB9.6", "azimuthal_scatter_recorded",
                    "; ".join(f"{lv}={data[lv]['scatter_med']:.4f}"
                              for lv in levels),
                    "RECORDED (no pass/fail; caveat carried until P11/Option C)",
                    "RECORDED"])
    conv = all(data[lv]["converged"] for lv in levels)
    print(f"\nGB9.6 RECORDED. body surface-Cp azimuthal scatter (medium ="
          f" {data['medium']['scatter_med']:.4f}) is the caveat every "
          f"wing-body surface-pressure claim carries until P11/Option C.")
    sys.exit(0 if conv else 1)


if __name__ == "__main__":
    main()
