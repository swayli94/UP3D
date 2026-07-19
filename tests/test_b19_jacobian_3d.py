"""
Track B / B19: the LS-Newton Jacobian is exact on MIXED-SIDE PLAIN elements.

Why this file exists at all is the point. `tests/test_b6_newton.py` already
carries an exact-Jacobian FD gate -- but it runs on the quasi-2D NACA mesh,
where the broken element class is INERT: it has 129 mixed-side plain elements
(uncut tets whose four nodes straddle the wake level set, beyond the spanwise
tip clip or where the sheet extension passes upstream of the TE) but 0 of them
touch a cut node -- the B19 erratum; the real invariant is "no mixed-plain
element READS an aux", not "the class is empty". AUX-TOUCHING mixed-plain
elements exist only in 3-D. So a real Jacobian defect lived through B6, B7,
B12-B18 and every M6 convergence gate, and was found by an external code
review (2026-07-17 C1), not by the suite.

The defect: `mass_conservation_coo` scatters a plain element onto MAIN dofs,
while `newton_side_data` builds its gradient -- hence its density, hence the
residual -- from the SIDE field, which reads cut nodes' AUX dofs
(`side_potentials`). Terms 2/3 used the scatter map for their columns too, so
the aux sensitivity landed on a main column: J != dR/dphi there. Measured
pre-fix on M6 coarse: rel err 1.146e-01 on probe directions along those aux
dofs vs 6.33e-10 elsewhere, and eps-INDEPENDENT across three decades (a
missing term, not FD noise). Full study: `cases/analysis/c1_ls_jacobian_fd/`.

The fix (B19 Leg A) splits the two roles: rows from `dofvec` (where R lands),
columns from `readvec` (what the side field reads -- `side_potentials`' own
per-node rule applied to every element, not just cut ones).

These tests are the standing guard so the class cannot regress unseen again.
The heavy 3-D FD gate is gated; the two structural checks are not, because
they are cheap and they are the ones that encode the invariant.
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
from pyfp3d.meshgen.wing3d import x_te
from pyfp3d.solve.picard_ls import (
    _farfield_split,
    _neumann_outlet_rhs,
    solve_multivalued_transonic,
)
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M1_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6"
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
GATES = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"

ALPHA_3D, B_SEMI = 3.06, 1.1963
UI, GA, C, MC, MCAP, RF = 1.0, 1.4, 1.5, 0.95, 3.0, 0.05


def _load(directory, level):
    path = directory / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    return read_mesh(path)


def _build_3d(mesh):
    a = np.radians(ALPHA_3D)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return cm, MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def _mixed_side_plain(mesh, cm):
    """Uncut elements whose nodes straddle the level set -- the broken class."""
    el = np.asarray(mesh.elements, dtype=np.int64)
    special = np.zeros(len(el), dtype=bool)
    special[cm.cut_elems] = True
    special[cm.te_lower_elems] = True
    side = cm.node_side[el]
    return np.flatnonzero((side.min(axis=1) != side.max(axis=1)) & ~special)


def test_read_map_generalizes_the_scatter_map():
    """The read map is not a NEW rule: on cut elements it must reproduce
    dofs_upper/dofs_lower exactly. That equivalence is what makes B19 a
    generalization of the existing per-element special case into the per-node
    rule it should always have been -- and it is why the fix cannot change
    anything on a mesh whose only special elements are cut ones."""
    mesh = _load(M0_DIR, "coarse")
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)

    ru, rl = mvop._side_readvecs()
    assert np.array_equal(ru[cm.cut_elems], cm.dofs_upper)
    assert np.array_equal(rl[cm.cut_elems], cm.dofs_lower)


def test_quasi_2d_blind_spot_is_aux_reading_elements_not_mixed_ones():
    """The structural blind spot, stated CORRECTLY and asserted.

    ★ Measured 2026-07-18, and it corrects both the C1 write-up and this
    author's first version of this test: the quasi-2D mesh is NOT free of
    mixed-side plain elements -- coarse has **129** of them. What it is free
    of is the subset that can actually reach an aux DOF: **0** of those 129
    touch a cut node. A mixed-side plain element only mis-scatters when its
    side field reads an aux value, i.e. when one of its nodes is a cut node,
    so THAT is the invariant which makes B6's 2.5-D FD gate blind to C1 --
    not the coarser "no mixed elements" claim.
    """
    mesh = _load(M0_DIR, "coarse")
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    el = np.asarray(mesh.elements, dtype=np.int64)
    msp = _mixed_side_plain(mesh, cm)
    assert len(msp) > 0, "premise check: the 2.5-D mesh does have mixed elements"
    aux_reading = [e for e in msp if (cm.ext_dof_of_node[el[e]] >= 0).any()]
    assert len(aux_reading) == 0, (
        f"{len(aux_reading)} quasi-2D mixed-side plain elements now touch a cut "
        "node: test_b6_newton's FD gate would cover the B19 class and this "
        "file's premise needs revisiting")


@pytest.mark.skipif(not GATES, reason="3-D FD gate; PYFP3D_TRANSONIC_GATES=1")
def test_jacobian_fd_on_mixed_side_plain_elements():
    """GB19.1: the exact Jacobian matches the central FD of the residual along
    probe directions supported on the aux dofs of cut nodes that touch a
    mixed-side plain element -- the directions that measured 1.146e-01 before
    the fix. A control direction on the remaining aux dofs must stay clean, so
    a pass cannot come from the probe simply missing the region."""
    mesh = _load(M1_DIR, "coarse")
    cm, mvop = _build_3d(mesh)
    el = np.asarray(mesh.elements, dtype=np.int64)

    msp = _mixed_side_plain(mesh, cm)
    assert len(msp) > 0, "M6 coarse must contain the class this gate exists for"

    m_inf = 0.70
    seed = solve_multivalued_transonic(mvop, mesh, m_inf, alpha_deg=ALPHA_3D,
                                       farfield="freestream", n_outer_level=200)
    phi0 = seed["phi_ext"]
    _, J, ff = _assemble_R_J(mvop, mesh, phi0, m_inf)
    assert mvop.n_nu_active > 50, "want an active pocket so Terms 2/3 are live"

    aux_of = cm.ext_dof_of_node
    msp_nodes = np.unique(el[msp])
    touched = aux_of[msp_nodes[aux_of[msp_nodes] >= 0]]
    control = np.setdiff1d(aux_of[aux_of >= 0], touched)
    assert len(touched) > 0 and len(control) > 0

    free = np.flatnonzero(~np.isin(np.arange(mvop.n_total), ff))
    rng = np.random.default_rng(0)
    for name, dofs in (("targeted", touched), ("control", control)):
        v = np.zeros(mvop.n_total)
        v[dofs] = rng.standard_normal(len(dofs))
        eps = 1e-7
        Jv = (J @ v)[free]
        Rp = _assemble_R_J(mvop, mesh, phi0 + eps * v, m_inf)[0]
        Rm = _assemble_R_J(mvop, mesh, phi0 - eps * v, m_inf)[0]
        fd = ((Rp - Rm) / (2 * eps))[free]
        rel = np.linalg.norm(Jv - fd) / (np.linalg.norm(fd) + 1e-30)
        assert rel < 1e-6, f"{name} probe: Jacobian FD rel err {rel:.2e}"


def _assemble_R_J(mvop, mesh, phi, m_inf):
    """Same objects and row map as tests/test_b6_newton.py::_assemble_R_J, so
    this gate measures the SHIPPED Newton system (3-D alpha/mesh)."""
    cm = mvop.cm
    te_aux = cm.ext_dof_of_node[cm.te_nodes].astype(np.int64)
    ff_nodes = _farfield_split(mesh, ALPHA_3D, UI)[3]
    ff_vals = freestream_phi(mesh.nodes[ff_nodes], ALPHA_3D, UI)
    b = _neumann_outlet_rhs(mesh, ALPHA_3D, UI, mvop.n_total)
    n_total = mvop.n_total
    keep_row = np.ones(n_total)
    keep_row[te_aux] = 0.0
    D = sp.diags(keep_row)
    nonte = np.asarray(mvop._nonte_rows, dtype=np.int64)
    row_map = np.arange(n_total, dtype=np.int64)
    row_map[nonte] = -1
    row_map[te_aux] = cm.te_nodes

    def remap(r, c, d):
        if not len(r):
            return r, c, d
        rr = row_map[r]
        k = rr >= 0
        return rr[k], c[k], d[k]

    phi = phi.copy()
    phi[ff_nodes] = ff_vals
    up, lo = mvop.newton_side_data(phi, m_inf, C, MC, GA, UI, MCAP, RF)
    A = mvop.assemble_matrix(rho_tilde=(up["rho_tilde"], lo["rho_tilde"]),
                             closure="wake_ls", te_kutta="pressure",
                             phi_ext=phi)
    R = A @ phi - b
    R[te_aux] = te_kutta_residual(mvop, phi)
    J = (D @ A).tocoo()
    pr, pc, pd = [J.row], [J.col], [J.data]
    for rr, cc, dd in (te_kutta_jacobian_coo(mvop, phi),
                       remap(*newton_terms23_side_coo(mvop.op, up, UI)),
                       remap(*newton_terms23_side_coo(mvop.op, lo, UI))):
        if len(rr):
            pr.append(rr)
            pc.append(cc)
            pd.append(dd)
    J = sp.coo_matrix((np.concatenate(pd),
                       (np.concatenate(pr), np.concatenate(pc))),
                      shape=(n_total, n_total)).tocsr()
    return R, J, ff_nodes


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
