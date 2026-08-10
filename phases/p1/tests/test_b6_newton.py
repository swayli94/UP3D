"""
Track B / B6-Newton: the level-set Newton solver (design_track_b.md §5.5).

The post-B6 re-derivation that closes the FP-fold cases the B6 Picard reaches
only as a bounded stall (medium M0.7875). Structure is simpler than the P8
conforming Newton: no Gamma DOF (implicit Kutta -> the TE jump is a solution
mode, no delta-Gamma elimination / Woodbury), wake-LS rows linear (constant
Jacobian), Terms 1-3 reused per side through the DOF indirection.

Fast tests: the exact-Jacobian FD check (the correctness gate) + a coarse
M0.70 Newton to a genuine discrete solution with a terminal quadratic drop.
The heavy M0.80 dual-mesh and medium M0.7875 fold runs are gated
(PYFP3D_TRANSONIC_GATES=1) / live in the demo.
"""

import os
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from pyfp3d.constraints.dirichlet import freestream_phi
from pyfp3d.kernels.cut_assembly import (
    newton_terms23_side_coo,
    te_kutta_jacobian_coo,
    te_kutta_residual,
)
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.solve.newton_ls import solve_multivalued_newton
from pyfp3d.solve.picard_ls import (
    _farfield_split,
    _neumann_outlet_rhs,
    solve_multivalued_transonic,
)
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
M3_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_wakefree_2.5d"
GATES = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"
AL, UI, GA, C, MC, MCAP, RF = 1.25, 1.0, 1.4, 1.5, 0.95, 3.0, 0.05


def _load(directory, level):
    path = directory / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    return read_mesh(path)


def _build(mesh):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def _assemble_R_J(mvop, mesh, phi, M):
    """The Newton residual and exact Jacobian at phi (mirrors newton_ls, the
    piece the FD check exercises)."""
    cm = mvop.cm
    te_aux = cm.ext_dof_of_node[cm.te_nodes].astype(np.int64)
    ff_nodes = _farfield_split(mesh, AL, UI)[3]
    ff_vals = freestream_phi(mesh.nodes[ff_nodes], AL, UI)
    b = _neumann_outlet_rhs(mesh, AL, UI, mvop.n_total)
    n_total = mvop.n_total
    keep_row = np.ones(n_total); keep_row[te_aux] = 0.0
    D = sp.diags(keep_row)
    nonte = np.asarray(mvop._nonte_rows, dtype=np.int64)
    row_map = np.arange(n_total, dtype=np.int64)
    row_map[nonte] = -1
    row_map[te_aux] = cm.te_nodes

    def remap(r, c, d):
        if not len(r):
            return r, c, d
        rr = row_map[r]; k = rr >= 0
        return rr[k], c[k], d[k]

    phi = phi.copy(); phi[ff_nodes] = ff_vals
    up, lo = mvop.newton_side_data(phi, M, C, MC, GA, UI, MCAP, RF)
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
    J = sp.coo_matrix((np.concatenate(pd),
                       (np.concatenate(pr), np.concatenate(pc))),
                      shape=(n_total, n_total)).tocsr()
    return R, J, ff_nodes


def test_newton_jacobian_fd():
    """The correctness gate: the exact Jacobian J @ v matches the central
    finite difference of the residual along v to ~1e-8 at a transonic state
    (supersonic pocket present -> Terms 2/3 active). The row-map (drop non-TE
    aux, reroute TE aux -> TE main) is what makes the linear wake-LS rows come
    out exact -- omitting it was FD-caught at 1e-4."""
    mesh, mvop = _build(_load(M0_DIR, "coarse"))
    M = 0.70
    seed = solve_multivalued_transonic(mvop, mesh, M, alpha_deg=AL,
                                       farfield="neumann", n_outer_level=300)
    phi0 = seed["phi_ext"]
    R0, J, ff = _assemble_R_J(mvop, mesh, phi0, M)
    assert mvop.n_nu_active > 50, "want an active pocket to exercise Terms 2/3"

    free = np.flatnonzero(~np.isin(np.arange(mvop.n_total), ff))
    rng = np.random.default_rng(0)
    v = np.zeros(mvop.n_total); v[free] = rng.standard_normal(len(free))
    Jv = (J @ v)[free]
    eps = 1e-7
    Rp = _assemble_R_J(mvop, mesh, phi0 + eps * v, M)[0]
    Rm = _assemble_R_J(mvop, mesh, phi0 - eps * v, M)[0]
    fd = ((Rp - Rm) / (2 * eps))[free]
    rel = np.linalg.norm(Jv - fd) / (np.linalg.norm(fd) + 1e-30)
    assert rel < 1e-6, f"Jacobian FD rel err {rel:.2e}"


def test_newton_converges_m070_coarse():
    """A coarse M0.70 Newton (seeded) reaches a GENUINE discrete solution:
    machine residual, 0 limited / 0 floored, with a terminal quadratic drop --
    the property the transonic Picard's bounded stall lacks."""
    mesh, mvop = _build(_load(M0_DIR, "coarse"))
    seed = solve_multivalued_transonic(mvop, mesh, 0.65, alpha_deg=AL,
                                       farfield="neumann", n_outer_level=400)
    r = solve_multivalued_newton(mvop, mesh, 0.70, alpha_deg=AL,
                                 farfield="neumann", phi_init=seed["phi_ext"],
                                 gamma_init=seed["gamma"], n_seed=30,
                                 n_newton_max=40)
    assert r["converged"], f"did not converge: {r['residual_history'][-3:]}"
    assert r["n_limited"] == 0 and r["n_floored"] == 0
    assert r["residual_history"][-1] < 1e-9
    # terminal quadratic: some late consecutive pair drops by > 2 digits
    h = np.array(r["residual_history"])
    ratios = h[1:] / (h[:-1] + 1e-300)
    assert ratios.min() < 1e-2, "expected a quadratic terminal drop"


@pytest.mark.skipif(not GATES, reason="PYFP3D_TRANSONIC_GATES != 1")
@pytest.mark.parametrize("mesh_dir", [M0_DIR, M3_DIR],
                         ids=["m0-embedded", "m3-wakefree"])
def test_g_b6_newton_m080_coarse(mesh_dir):
    """Gated: coarse M0.80 LS Newton reaches a machine-converged discrete
    solution, 0 lim/flr, shock in the Picard-quality band. Measured Gamma:
    M0 0.2124 (-7.4% of conforming Newton), M3 0.2322 (+1.2%)."""
    mesh, mvop = _build(_load(mesh_dir, "coarse"))
    seed = solve_multivalued_transonic(mvop, mesh, 0.75, alpha_deg=AL,
                                       farfield="neumann", n_outer_level=600)
    r = solve_multivalued_newton(mvop, mesh, 0.80, alpha_deg=AL,
                                 farfield="neumann", phi_init=seed["phi_ext"],
                                 gamma_init=seed["gamma"], n_seed=40,
                                 n_newton_max=40)
    assert r["converged"]
    assert r["residual_history"][-1] < 1e-8
    assert r["n_limited"] == 0 and r["n_floored"] == 0
    assert abs(r["gamma"] / 0.2295 - 1.0) < 0.10
