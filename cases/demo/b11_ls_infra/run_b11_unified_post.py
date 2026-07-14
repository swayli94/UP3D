"""
Track B / B11 demo -- the unified post-processing dispatch is bit-identical.

`post/surface.py` (conforming) and `post/surface_ls.py` (level-set) now share
private cores (`_cp_from_q2`, `_pressure_force`, `_wall_plane_crossings`,
`_section_curve_dict`, `_d11_wall_state`) under one keyword-dispatched upper
layer, `post/unified.py` (`wall_cp` / `wall_forces` / `section_cp`, selecting
the path by `phi=` vs `mvop=,phi_ext=`). This demo solves ONE level-set case on
NACA0012, extracts its wall Cp through BOTH the unified entry point and the
legacy per-path functions on BOTH dispatch forms (the conforming form is fed the
main-potential slice of the same solution), and self-checks that every unified
output equals the legacy one to the last bit (max|Δcp| = 0.0).

Outputs: results/cp_unified_overlay.png, results/summary.csv.
Exit 0 iff the identity holds. Standalone; runtime ~30 s.
"""

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, INK, MESH_DIR, S1_BLUE, S2_AQUA,
    CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402

from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.post import unified  # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve  # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients  # noqa: E402
from pyfp3d.post.surface_ls import (  # noqa: E402
    section_cp_curve_levelset, wall_cp_levelset,
)
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import (  # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

OUT = Path(__file__).resolve().parent / "results"
NACA_DIR = MESH_DIR / "naca0012_2.5d"
LEVEL = "coarse"
ALPHA = 2.0
M_INF = 0.5


def build(mesh):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
                       direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return wls, cm, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                        levelset=wls)


def main():
    apply_style()
    checks = CheckList("Track B / B11 -- unified post-processing dispatch")
    path = NACA_DIR / f"{LEVEL}.msh"
    if not path.exists():
        print(f"SKIP: {path} not generated (gitignored)")
        return 0
    mesh = read_mesh(path)
    wls, cm, mvop = build(mesh)

    t0 = time.time()
    r = solve_multivalued_lifting(mvop, mesh, M_INF, alpha_deg=ALPHA)
    phi_ext = r["phi_ext"]
    phi_main = mvop.main_potential(phi_ext)
    print(f"  solved: gamma={r['gamma']:+.5f} ({r['n_outer']} outers, "
          f"{time.time()-t0:.0f}s)")

    z = 0.5 * (float(mesh.nodes[:, 2].min()) + float(mesh.nodes[:, 2].max()))
    rows = []

    # --- level-set path: unified vs legacy ---------------------------------
    u_cp = unified.wall_cp(mesh, mvop=mvop, phi_ext=phi_ext, m_inf=M_INF)
    l_cp = wall_cp_levelset(mesh, mvop, phi_ext, m_inf=M_INF)
    d_ls_cp = float(np.max(np.abs(u_cp["cp"] - l_cp["cp"])))
    u_sec = unified.section_cp(mesh, mvop=mvop, phi_ext=phi_ext, z=z, m_inf=M_INF)
    l_sec = section_cp_curve_levelset(mesh, mvop, phi_ext, z=z, m_inf=M_INF)
    d_ls_sec = float(np.max(np.abs(u_sec["cp_lower"] - l_sec["cp_lower"])))
    rows.append(("levelset", "wall_cp", d_ls_cp))
    rows.append(("levelset", "section_cp_lower", d_ls_sec))
    checks.add("G11.2", "LS wall_cp unified==legacy", d_ls_cp, "== 0.0",
               d_ls_cp == 0.0)
    checks.add("G11.2", "LS section_cp unified==legacy", d_ls_sec, "== 0.0",
               d_ls_sec == 0.0)

    # --- conforming path: unified vs legacy (fed the main-potential slice) --
    u_fc = unified.wall_forces(mesh, phi=phi_main, alpha_deg=ALPHA, m_inf=M_INF)
    l_fc = wall_force_coefficients(mesh.nodes, mesh.elements,
                                   mesh.boundary_faces["wall"], phi_main,
                                   alpha_deg=ALPHA, m_inf=M_INF)
    d_cf_cp = float(np.max(np.abs(u_fc["cp_tri"] - l_fc["cp_tri"])))
    d_cf_cl = abs(u_fc["cl"] - l_fc["cl"])
    u_csec = unified.section_cp(mesh, phi=phi_main, z=z, m_inf=M_INF)
    l_csec = section_cp_curve(mesh, phi_main, z=z, m_inf=M_INF)
    d_cf_sec = float(np.max(np.abs(u_csec["cp_upper"] - l_csec["cp_upper"])))
    rows.append(("conforming", "wall_forces_cp", d_cf_cp))
    rows.append(("conforming", "wall_forces_cl", d_cf_cl))
    rows.append(("conforming", "section_cp_upper", d_cf_sec))
    checks.add("G11.2", "conforming wall_forces unified==legacy",
               max(d_cf_cp, d_cf_cl), "== 0.0", d_cf_cp == 0.0 and d_cf_cl == 0.0)
    checks.add("G11.2", "conforming section_cp unified==legacy", d_cf_sec,
               "== 0.0", d_cf_sec == 0.0)

    # --- overlay plot ------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(u_sec["x_upper"], u_sec["cp_upper"], color=S1_BLUE, lw=1.8,
            label="unified upper")
    ax.plot(u_sec["x_lower"], u_sec["cp_lower"], color=S2_AQUA, lw=1.8,
            label="unified lower")
    ax.plot(l_sec["x_upper"], l_sec["cp_upper"], color=INK, lw=0.0,
            marker="o", ms=3, mfc="none", label="legacy (overlaid)")
    ax.plot(l_sec["x_lower"], l_sec["cp_lower"], color=INK, lw=0.0,
            marker="o", ms=3, mfc="none")
    ax.invert_yaxis()
    ax.set_xlabel("x/c")
    ax.set_ylabel("Cp")
    ax.set_title(f"B11 unified section Cp -- NACA0012 M{M_INF} a{ALPHA} "
                 f"(unified == legacy, max|Δcp| = {max(d_ls_sec, d_cf_sec):.1e})")
    ax.legend(frameon=False)
    finish(fig, OUT, "cp_unified_overlay.png")

    write_csv(OUT, "summary.csv", "path,quantity,max_abs_diff", rows)
    return checks.report(OUT, "checks_unified_post.csv")


if __name__ == "__main__":
    sys.exit(main())
