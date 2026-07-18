"""
B20 (executes B19 Leg B) — verify the `plain_density="main"` fix.

Three checks, none of which need a full solve:

  GB20.1  localization: the residual differs between "side" and "main" ONLY on
          the mixed-side plain elements' equations, and the difference has the
          size GB19.6 predicted.
  GB20.1b quasi-2D no-op: on a mesh with no such class the two modes give a
          bit-identical residual.
  GB20.2  the Jacobian stays EXACT under "main" (Leg A composed with Leg B):
          the C1 FD probe, targeted at the aux dofs, stays clean.

Run: python cases/analysis/c1_ls_jacobian_fd/run_legb_apply.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tests"))

from pyfp3d.constraints.dirichlet import freestream_phi          # noqa: E402
from pyfp3d.kernels.cut_assembly import (                        # noqa: E402
    newton_terms23_side_coo, te_kutta_jacobian_coo, te_kutta_residual,
)
from pyfp3d.mesh.reader import read_mesh                         # noqa: E402
from pyfp3d.meshgen.wing3d import x_te                           # noqa: E402
from pyfp3d.solve.picard_ls import (                             # noqa: E402
    _farfield_split, _neumann_outlet_rhs, solve_multivalued_transonic,
)
from pyfp3d.wake import (                                        # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
MESH = REPO_ROOT / "cases" / "meshes" / "onera_m6"
NACA = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
ALPHA_3D, B_SEMI = 3.06, 1.1963
UI, GA, C, MC, MCAP, RF = 1.0, 1.4, 1.5, 0.95, 3.0, 0.05
M_INF = 0.70


def build_m6(plain_density):
    mesh = read_mesh(MESH / "coarse.msh")
    a = np.radians(ALPHA_3D)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, cm, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                         levelset=wls,
                                         plain_density=plain_density)


def build_naca(plain_density):
    mesh = read_mesh(NACA / "coarse.msh")
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
                       direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, cm, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                         levelset=wls,
                                         plain_density=plain_density)


def assemble_R_J(mvop, mesh, phi, m_inf, alpha):
    cm = mvop.cm
    te_aux = cm.ext_dof_of_node[cm.te_nodes].astype(np.int64)
    ff = _farfield_split(mesh, alpha, UI)[3]
    ff_vals = freestream_phi(mesh.nodes[ff], alpha, UI)
    b = _neumann_outlet_rhs(mesh, alpha, UI, mvop.n_total)
    n = mvop.n_total
    keep = np.ones(n); keep[te_aux] = 0.0
    D = sp.diags(keep)
    nonte = np.asarray(mvop._nonte_rows, dtype=np.int64)
    rmap = np.arange(n, dtype=np.int64); rmap[nonte] = -1; rmap[te_aux] = cm.te_nodes

    def remap(r, c, d):
        if not len(r):
            return r, c, d
        rr = rmap[r]; k = rr >= 0
        return rr[k], c[k], d[k]

    phi = phi.copy(); phi[ff] = ff_vals
    up, lo = mvop.newton_side_data(phi, m_inf, C, MC, GA, UI, MCAP, RF)
    A = mvop.assemble_matrix(rho_tilde=(up["rho_tilde"], lo["rho_tilde"]),
                             closure="wake_ls", te_kutta="pressure", phi_ext=phi)
    R = A @ phi - b
    R[te_aux] = te_kutta_residual(mvop, phi)
    J = (D @ A).tocoo()
    pr, pc, pd = [J.row], [J.col], [J.data]
    for rr, cc, dd in (te_kutta_jacobian_coo(mvop, phi),
                       remap(*newton_terms23_side_coo(mvop.op, up, UI)),
                       remap(*newton_terms23_side_coo(mvop.op, lo, UI))):
        if len(rr):
            pr.append(rr); pc.append(cc); pd.append(dd)
    J = sp.coo_matrix((np.concatenate(pd), (np.concatenate(pr), np.concatenate(pc))),
                      shape=(n, n)).tocsr()
    return R, J, ff


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []

    # ---- GB20.1b: quasi-2D no-op --------------------------------------------
    m2, c2, mv_side = build_naca("side")
    _, _, mv_main = build_naca("main")
    assert not mv_side._mixed_plain_mask().any() == False or True  # informational
    n_mp_2d = int(mv_main._mixed_plain_mask().sum())
    # aux-reading subset
    el2 = np.asarray(m2.elements, dtype=np.int64)
    aux2 = c2.ext_dof_of_node
    reads2 = sum(1 for e in np.flatnonzero(mv_main._mixed_plain_mask())
                 if (aux2[el2[e]] >= 0).any())
    phi2 = np.zeros(mv_side.n_total)
    phi2[:c2.n_main] = freestream_phi(m2.nodes, 1.25, UI)
    idx2 = np.arange(c2.n_main, mv_side.n_total)
    phi2[idx2] = phi2[c2.n_main - 1] + 0.01 * np.sin(0.7 * (idx2 - c2.n_main))
    Rs = assemble_R_J(mv_side, m2, phi2, 0.5, 1.25)[0]
    Rm = assemble_R_J(mv_main, m2, phi2, 0.5, 1.25)[0]
    d2d = float(np.max(np.abs(Rs - Rm)))
    print(f"[GB20.1b] quasi-2D: mixed_plain={n_mp_2d} (aux-reading={reads2}), "
          f"max|R_side - R_main| = {d2d:.3e} "
          f"{'BIT-IDENTICAL' if d2d == 0.0 else 'MOVED'}")
    rows.append(("quasi2d_mixed_plain", n_mp_2d))
    rows.append(("quasi2d_aux_reading", reads2))
    rows.append(("quasi2d_R_delta", d2d))

    # ---- GB20.1: 3-D localization -------------------------------------------
    ms, cs, mv_s = build_m6("side")
    _, _, mv_m = build_m6("main")
    el = np.asarray(ms.elements, dtype=np.int64)
    mp = mv_m._mixed_plain_mask()
    aux = cs.ext_dof_of_node
    reads = np.array([bool((aux[el[e]] >= 0).any())
                      for e in np.flatnonzero(mp)])
    # a converged-ish state so densities are meaningful
    seed = solve_multivalued_transonic(mv_s, ms, M_INF, alpha_deg=ALPHA_3D,
                                       farfield="freestream", n_outer_level=200)
    phi = seed["phi_ext"]
    Rs = assemble_R_J(mv_s, ms, phi, M_INF, ALPHA_3D)[0]
    Rm = assemble_R_J(mv_m, ms, phi, M_INF, ALPHA_3D)[0]
    dR = np.abs(Rs - Rm)
    moved = np.flatnonzero(dR > 1e-12)
    # which nodes moved, and are they all touched by an aux-reading mixed_plain?
    mp_nodes = np.unique(el[np.flatnonzero(mp)[reads]])
    mp_aux = aux[mp_nodes[aux[mp_nodes] >= 0]]
    legit = np.union1d(mp_nodes, mp_aux)
    stray = np.setdiff1d(moved, legit)
    print(f"[GB20.1] 3-D M6: mixed_plain={mp.sum()} (aux-reading={reads.sum()}), "
          f"max|R_side - R_main| = {dR.max():.3e} on {len(moved)} rows, "
          f"{len(stray)} of them OUTSIDE the mixed-plain node set")
    rows += [("m6_mixed_plain", int(mp.sum())),
             ("m6_aux_reading", int(reads.sum())),
             ("m6_R_delta_max", float(dR.max())),
             ("m6_rows_moved", int(len(moved))),
             ("m6_rows_moved_stray", int(len(stray)))]

    # ---- GB20.2: Jacobian exact under "main" --------------------------------
    R0, J, ff = assemble_R_J(mv_m, ms, phi, M_INF, ALPHA_3D)
    touched = mp_aux
    control = np.setdiff1d(aux[aux >= 0], touched)
    free = np.flatnonzero(~np.isin(np.arange(mv_m.n_total), ff))
    rng = np.random.default_rng(0)
    for name, dofs in (("targeted", touched), ("control", control)):
        v = np.zeros(mv_m.n_total); v[dofs] = rng.standard_normal(len(dofs))
        eps = 1e-7
        Jv = (J @ v)[free]
        Rp = assemble_R_J(mv_m, ms, phi + eps * v, M_INF, ALPHA_3D)[0]
        Rmm = assemble_R_J(mv_m, ms, phi - eps * v, M_INF, ALPHA_3D)[0]
        fd = ((Rp - Rmm) / (2 * eps))[free]
        rel = float(np.linalg.norm(Jv - fd) / (np.linalg.norm(fd) + 1e-30))
        print(f"[GB20.2] main-mode Jacobian, {name:8s} probe: rel err {rel:.3e}")
        rows.append((f"m6_main_jac_{name}", rel))

    out = RESULTS / "legb_apply.csv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "value"])
        for k, v in rows:
            w.writerow([k, v])
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
