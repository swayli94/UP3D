"""
P2 wake-cut unit and gate tests: node duplication, topology asserts,
G2.1 (freestream preservation on the cut mesh) and G2.2 (prescribed-jump
consistency). Roadmap P2; design.md Sec 4.

The synthetic case is a quasi-2D "flat-plate wake" built with the M0
extruder (no Gmsh needed): rectangle [0,4]x[-1,1], wall = left edge,
farfield = the other three edges, wake = the interior line y = 0 from the
wall to the farfield. Same cut topology as the airfoil (TE nodes on the
wall, far-end wake nodes on the Dirichlet boundary) at a few hundred tets.
"""

import numpy as np
import pytest

from pyfp3d.kernels.residual import assemble_stiffness_matrix
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import assert_wake_topology, cut_wake
from pyfp3d.meshgen.extrude import extrude_single_layer
from pyfp3d.constraints.wake import WakeConstraint, kutta_targets
from pyfp3d.solve.picard import solve_laplace_lifting


def synthetic_wake_mesh(nx: int = 16, ny: int = 8, dz: float = 0.25):
    """Structured quasi-2D strip with an interior wake sheet along y = 0.

    ny must be even so the wake line lies on a grid line.
    """
    assert ny % 2 == 0
    xs = np.linspace(0.0, 4.0, nx + 1)
    ys = np.linspace(-1.0, 1.0, ny + 1)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    points2d = np.stack([X.ravel(), Y.ravel()], axis=1)

    def nid(i, j):
        return i * (ny + 1) + j

    tris = []
    for i in range(nx):
        for j in range(ny):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            tris.append((a, b, c))
            tris.append((a, c, d))
    triangles = np.asarray(tris, dtype=np.int64)

    jmid = ny // 2
    wall = np.array([[nid(0, j), nid(0, j + 1)] for j in range(ny)])
    far = np.concatenate([
        [[nid(nx, j), nid(nx, j + 1)] for j in range(ny)],
        [[nid(i, 0), nid(i + 1, 0)] for i in range(nx)],
        [[nid(i, ny), nid(i + 1, ny)] for i in range(nx)],
    ])
    wake = np.array([[nid(i, jmid), nid(i + 1, jmid)] for i in range(nx)])

    return extrude_single_layer(
        points2d, triangles,
        edge_groups={"wall": wall, "farfield": far},
        interior_edge_groups={"wake": wake},
        dz=dz, name="synthetic_wake",
    )


@pytest.fixture(scope="module")
def synthetic_cut():
    mesh = synthetic_wake_mesh()
    mesh_cut, wc = cut_wake(mesh)
    return mesh, mesh_cut, wc


def test_cut_counts_and_stations(synthetic_cut):
    mesh, mesh_cut, wc = synthetic_cut
    n_wake_nodes = len(np.unique(mesh.boundary_faces["wake"]))
    # Every sheet node (TE included) is duplicated -- see wake_cut docstring.
    assert len(wc.slave_nodes) == n_wake_nodes
    assert len(mesh_cut.nodes) == len(mesh.nodes) + n_wake_nodes
    assert np.isin(wc.te_nodes, wc.master_nodes).all()
    # Quasi-2D extrusion: one station, Gamma is a single scalar (M0 spec).
    assert wc.n_stations == 1
    assert len(wc.te_nodes) == 2  # one TE node per z-plane
    # Probes: one node off the TE, off the wake, one per side and TE node.
    assert len(wc.kutta_upper) == len(wc.te_nodes)
    nodes = mesh_cut.nodes
    assert np.all(nodes[wc.kutta_upper, 1] > 0)
    assert np.all(nodes[wc.kutta_lower, 1] < 0)


def test_plus_side_attachment(synthetic_cut):
    """Elements above the sheet reference slaves, below reference masters."""
    _, mesh_cut, wc = synthetic_cut
    el = np.asarray(mesh_cut.elements, dtype=np.int64)
    centroids = mesh_cut.nodes[el].mean(axis=1)
    has_slave = np.isin(el, wc.slave_nodes).any(axis=1)
    has_master = np.isin(el, wc.master_nodes).any(axis=1)
    assert np.all(centroids[has_slave, 1] > 0)
    assert np.all(centroids[has_master, 1] < 0)


def test_topology_assert_fires_on_broken_cut(synthetic_cut):
    """Re-pointing one slave reference back to its master must be caught."""
    _, mesh_cut, wc = synthetic_cut
    broken = read_back = mesh_cut  # alias; we operate on copies below
    import copy
    broken = copy.deepcopy(read_back)
    el = np.asarray(broken.elements, dtype=np.int64)
    e, k = np.argwhere(np.isin(el, wc.slave_nodes))[0]
    slave = el[e, k]
    master = wc.master_nodes[np.where(wc.slave_nodes == slave)[0][0]]
    el[e, k] = master
    broken.elements = el.astype(np.int32)
    with pytest.raises(AssertionError):
        assert_wake_topology(broken, wc)


def test_g21_freestream_on_cut_mesh(synthetic_cut):
    """G2.1 = V0-with-cut: phi = x, Gamma = 0 -> reduced residual machine
    zero at every free node except wall nodes (whose rows carry the
    physical through-wall flux of phi = x, exactly as in V0)."""
    _, mesh_cut, wc = synthetic_cut
    A = assemble_stiffness_matrix(mesh_cut.nodes, mesh_cut.elements)
    con = WakeConstraint(A, wc)
    R = con.A_reduced @ mesh_cut.nodes[: con.n_reduced, 0]

    to_master = np.arange(len(mesh_cut.nodes))
    to_master[wc.slave_nodes] = wc.master_nodes
    check = np.ones(con.n_reduced, dtype=bool)
    for tag in ("wall", "farfield"):
        check[np.unique(to_master[np.unique(mesh_cut.boundary_faces[tag])])] = False
    assert np.max(np.abs(R[check])) < 1e-12
    # The cut machinery itself: folded wake-master rows are interior-like.
    masters_int = wc.master_nodes[check[wc.master_nodes]]
    assert np.max(np.abs(R[masters_int])) < 1e-13


def test_g22_prescribed_jump(synthetic_cut):
    """G2.2: Gamma = const prescribed -> [phi] reproduced exactly at every
    wake node pair and the residual clean at the cut."""
    _, mesh_cut, wc = synthetic_cut
    gamma = 0.37
    r = solve_laplace_lifting(mesh_cut, wc, alpha_deg=0.0, gamma_fixed=gamma)
    jump = r["phi"][wc.slave_nodes] - r["phi"][wc.master_nodes]
    np.testing.assert_allclose(jump, gamma, rtol=0, atol=1e-13)
    assert r["residual_norm"] < 1e-8


def test_g21_g22_on_naca_coarse(mesh_dir):
    """The same two gates on the real M0 NACA0012 coarse mesh."""
    mesh = read_mesh(mesh_dir / "naca0012_2.5d" / "coarse.msh")
    mesh_cut, wc = cut_wake(mesh)
    assert wc.n_stations == 1

    A = assemble_stiffness_matrix(mesh_cut.nodes, mesh_cut.elements)
    con = WakeConstraint(A, wc)
    R = con.A_reduced @ mesh_cut.nodes[: con.n_reduced, 0]
    to_master = np.arange(len(mesh_cut.nodes))
    to_master[wc.slave_nodes] = wc.master_nodes
    check = np.ones(con.n_reduced, dtype=bool)
    for tag in ("wall", "farfield"):
        check[np.unique(to_master[np.unique(mesh_cut.boundary_faces[tag])])] = False
    assert np.max(np.abs(R[check])) < 1e-12  # G2.1

    gamma = 0.3
    r = solve_laplace_lifting(mesh_cut, wc, alpha_deg=0.0, gamma_fixed=gamma)
    jump = r["phi"][wc.slave_nodes] - r["phi"][wc.master_nodes]
    np.testing.assert_allclose(jump, gamma, rtol=0, atol=1e-12)  # G2.2
    assert r["residual_norm"] < 1e-8

    # Kutta probe plumbing: the prescribed circulation flow's probe jump
    # is within a few percent of Gamma (smooth-flow jump one node off TE).
    t = kutta_targets(r["phi"], wc)
    assert abs(t[0] - gamma) < 0.15 * gamma


def test_topology_asserts_all_wake_meshes(mesh_dir):
    """Agent-rules hard rule 7: wake-cut topology asserts on every mesh in
    cases/meshes/ that carries a wake sheet (meshes without a 'wake' tag
    have nothing to cut and are skipped by construction)."""
    ran = 0
    for msh in sorted(mesh_dir.glob("*/*.msh")):
        try:
            mesh = read_mesh(msh)
        except ValueError:
            # Surface-only assets (e.g. cessna/cessna_surface.msh) carry no
            # tets and are not solver meshes; hard rule 7 covers volume
            # meshes with a wake sheet.
            continue
        if "wake" not in mesh.boundary_faces:
            continue
        mesh_cut, wc = cut_wake(mesh)  # runs assert_wake_topology internally
        assert_wake_topology(mesh_cut, wc)
        ran += 1
    assert ran >= 2, "expected at least the naca0012_2.5d coarse+medium family"
