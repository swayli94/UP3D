"""GV5.5 TE-outflow row replacement tests (a band-(a) suite leg).

The V1 variant of bench/studies/v5_5_te_floor/PRE_REGISTRATION.md:
first-order extrapolation rows at the TE outflow nodes -- row 0
(x-momentum, the delta carrier) replaced by U[i,0] - U[up,0], row 2
(kinetic-energy, the shape carrier) by U[i,1] - U[up,1]; the matching
J_e rows zeroed; default OFF = legacy bit-identical. Synthetic small
plates + FD, no heavy compute; the te_outflow_pairs helper is checked
on the committed 2.5-D NACA strip case (mesh read only).
"""

import numpy as np
import pytest

from pyfp3d.viscous import closures as C
from pyfp3d.viscous.ibl3 import IBL3Solver
from pyfp3d.viscous.surface_mesh import (
    SurfaceMesh,
    structured_rectangle_surface,
)

Q0 = 1.0
RHO0 = 1.0
MU0 = 1.0e-5


def _plate_mesh(nx=8, nz=4, x0=0.2, x1=1.2, z0=-0.3, z1=0.3):
    xyz, tris = structured_rectangle_surface(x0, x1, z0, z1, nx, nz)
    return SurfaceMesh.from_wall_faces(xyz, tris)


def _outflow_pairs(sm, x1=1.2):
    """(te_node, upstream_node) pairs on the plate's outflow edge: for
    each node at x = x1, the element-adjacent node with the largest
    x < x1 (deterministic)."""
    nbrs = [set() for _ in range(sm.n_node)]
    for t in range(sm.n_tri):
        a, b, c = (int(v) for v in sm.triangles[t])
        nbrs[a].update((b, c))
        nbrs[b].update((a, c))
        nbrs[c].update((a, b))
    pairs = []
    for i in np.where(np.abs(sm.xyz[:, 0] - x1) < 1.0e-12)[0]:
        cand = [j for j in nbrs[int(i)] if sm.xyz[j, 0] < x1 - 1.0e-12]
        cand.sort(key=lambda j: (-sm.xyz[j, 0], j))
        pairs.append((int(i), cand[0]))
    return np.asarray(pairs, dtype=np.int64)


def _lam_seed_field(xyz):
    n = xyz.shape[0]
    U = np.zeros((n, 6))
    for i in range(n):
        U[i] = C.blasius_seed(max(xyz[i, 0], 1.0e-3), q=Q0, rho=RHO0,
                              mu=MU0)
    return U


def _make_solver(sm, flags, te=False, pairs=None):
    n = sm.xyz.shape[0]
    u_e = np.zeros((n, 3))
    u_e[:, 0] = Q0
    inflow = np.abs(sm.xyz[:, 0] - 0.2) < 1.0e-12
    st = C.blasius_seed(0.2, q=Q0, rho=RHO0, mu=MU0)
    return IBL3Solver(
        sm, u_e, RHO0, MU0, 0.0, flags, inflow, st,
        te_pairs=pairs, te_extrapolate=te,
    )


def _physical_state(sm, rng):
    U = _lam_seed_field(sm.xyz)
    U[:, 1] += 0.5 * np.sin(3.0 * sm.xyz[:, 0])
    U[:, 2] = 0.3 * np.cos(5.0 * sm.xyz[:, 2])
    U[:, 3] = 0.1 * np.sin(4.0 * sm.xyz[:, 0])
    U[:, 4] = 0.02
    U[:, 5] = 0.005
    return U


# ---------------------------------------------------------------------------
# default-off bit identity + the row-replacement structure
# ---------------------------------------------------------------------------


def test_default_off_bit_identical():
    sm = _plate_mesh()
    flags = np.zeros(sm.n_node, dtype=np.int64)
    pairs = _outflow_pairs(sm)
    s_off = _make_solver(sm, flags)
    s_kw = _make_solver(sm, flags, te=False, pairs=pairs)
    U = _lam_seed_field(sm.xyz)
    R1, J1 = s_off.residual_jacobian(U)
    R2, J2 = s_kw.residual_jacobian(U)
    assert np.array_equal(R1, R2)
    assert np.array_equal(J1.data, J2.data)


def test_row_replacement_residual_rows():
    sm = _plate_mesh()
    flags = np.zeros(sm.n_node, dtype=np.int64)
    pairs = _outflow_pairs(sm)
    s_off = _make_solver(sm, flags)
    s_on = _make_solver(sm, flags, te=True, pairs=pairs)
    U = _physical_state(sm, None)
    R_off = s_off.residual(U)
    R_on = s_on.residual(U)
    te_i = pairs[:, 0]
    for i, up in pairs:
        assert R_on[i, 0] == U[i, 0] - U[up, 0]
        assert R_on[i, 2] == U[i, 1] - U[up, 1]
    # every other row is untouched: rows 1,3,4,5 at the TE nodes and
    # every row elsewhere
    mask = np.ones_like(R_on, dtype=bool)
    mask[te_i, 0] = False
    mask[te_i, 2] = False
    assert np.array_equal(R_on[mask], R_off[mask])


def test_jacobian_rows_structure():
    sm = _plate_mesh()
    flags = np.zeros(sm.n_node, dtype=np.int64)
    pairs = _outflow_pairs(sm)
    s_on = _make_solver(sm, flags, te=True, pairs=pairs)
    U = _physical_state(sm, None)
    _, J = s_on.residual_jacobian(U)
    J = J.tocsr()
    for i, up in pairs:
        for r, c_diag, c_off in ((6 * i + 0, 6 * i + 0, 6 * up + 0),
                                 (6 * i + 2, 6 * i + 1, 6 * up + 1)):
            row = J.indices[J.indptr[r]:J.indptr[r + 1]]
            dat = J.data[J.indptr[r]:J.indptr[r + 1]]
            got = dict(zip(row.tolist(), dat.tolist()))
            # exactly two nonzero entries: +1 on the constrained column,
            # -1 on the upstream partner; the rest zeroed by the
            # row replacement (the CSR pattern itself is denser)
            nz = {c: v for c, v in got.items() if v != 0.0}
            assert nz == {c_diag: 1.0, c_off: -1.0}


# ---------------------------------------------------------------------------
# FD exactness with the flag ON (the committed FD discipline)
# ---------------------------------------------------------------------------


def test_jacobian_fd_with_te_rows():
    sm = _plate_mesh(nx=5, nz=3)
    n = sm.n_node
    rng = np.random.default_rng(1)
    flags = (rng.random(n) > 0.5).astype(np.int64)
    pairs = _outflow_pairs(sm)
    solver = _make_solver(sm, flags, te=True, pairs=pairs)
    U = _physical_state(sm, rng)
    U[:, 4] = np.where(flags == 1, 0.02, C.CTAU_LAM)
    R0, J = solver.residual_jacobian(U)
    assert np.all(np.isfinite(R0))
    assert np.all(np.isfinite(J.data))
    eps = 1.0e-6
    for _ in range(4):
        v = rng.standard_normal((n, 6))
        v /= np.max(np.abs(v))
        Rp = solver.residual(U + eps * v)
        Rm = solver.residual(U - eps * v)
        fd = (Rp - Rm).ravel() / (2.0 * eps)
        an = J @ v.ravel()
        err = np.max(np.abs(fd - an)) / max(np.max(np.abs(an)), 1.0e-12)
        assert err < 1.0e-5, f"Jacobian FD mismatch (te rows on): {err:.3e}"


def test_edge_jacobian_rows_zeroed():
    sm = _plate_mesh()
    flags = np.zeros(sm.n_node, dtype=np.int64)
    pairs = _outflow_pairs(sm)
    s_on = _make_solver(sm, flags, te=True, pairs=pairs)
    U = _physical_state(sm, None)
    _, Je = s_on.residual_edge_jacobian(U)
    Je = Je.tocsr()
    for i, _ in pairs:
        for r in (6 * i + 0, 6 * i + 2):
            seg = Je.data[Je.indptr[r]:Je.indptr[r + 1]]
            assert np.all(seg == 0.0)


# ---------------------------------------------------------------------------
# guards
# ---------------------------------------------------------------------------


def test_out_of_pattern_guard_raises():
    sm = _plate_mesh()
    flags = np.zeros(sm.n_node, dtype=np.int64)
    far = int(np.argmax(sm.xyz[:, 0]))  # node at x1 ...
    near = int(np.argmin(sm.xyz[:, 0]))  # ... paired with the x0 node
    with pytest.raises(ValueError, match="out of the CSR pattern"):
        _make_solver(sm, flags, te=True,
                     pairs=np.asarray([[far, near]], dtype=np.int64))


def test_empty_pairs_guard_raises():
    sm = _plate_mesh()
    flags = np.zeros(sm.n_node, dtype=np.int64)
    with pytest.raises(ValueError, match="empty te_pairs"):
        _make_solver(sm, flags, te=True, pairs=None)


# ---------------------------------------------------------------------------
# the quasi-2-D lock + convergence with the flag ON (a smoke, not the gate)
# ---------------------------------------------------------------------------


def test_lock_and_convergence_with_te_rows():
    sm = _plate_mesh(nx=10, nz=5)
    flags = np.zeros(sm.n_node, dtype=np.int64)
    pairs = _outflow_pairs(sm)
    solver = _make_solver(sm, flags, te=True, pairs=pairs)
    U0 = _lam_seed_field(sm.xyz)
    U, info = solver.solve(U0, tol=1.0e-9, max_iter=60)
    assert info["converged"]
    assert np.all(np.isfinite(U))
    # the quasi-2-D lock is retained (the replaced rows involve only the
    # delta / A columns -- no crossflow injection); the laminar Ctau pins
    # sit at CTAU_LAM exactly as in the legacy path
    inflow = np.abs(sm.xyz[:, 0] - 0.2) < 1.0e-12
    assert np.max(np.abs(U[~inflow, 2])) < 1.0e-12
    assert np.max(np.abs(U[~inflow, 3])) < 1.0e-12
    assert np.allclose(U[~inflow, 5], C.CTAU_LAM, atol=1.0e-12)


# ---------------------------------------------------------------------------
# te_outflow_pairs on the committed 2.5-D NACA strip case
# ---------------------------------------------------------------------------


def test_te_outflow_pairs_on_naca_strip():
    from pyfp3d.viscous.coupling import te_outflow_pairs
    from tests.v5_state import build_naca_case

    mc, wc, cfg, case = build_naca_case()
    st = case.stations
    pairs = te_outflow_pairs(case)
    te_row = int(np.argmax(st.xc))
    te_nodes = set(np.where(st.station_of == te_row)[0].tolist())
    # every TE-station node is covered exactly once, upstream partners
    # sit outside the TE station on the same side, and the pair shares
    # an element (the CSR in-pattern discipline)
    assert sorted(pairs[:, 0].tolist()) == sorted(te_nodes)
    assert len(pairs) >= 2  # the wake-cut upper/lower TE copies
    in_elem = set()
    for t in range(case.sm.n_tri):
        a, b, c = (int(v) for v in case.sm.triangles[t])
        for u, v in ((a, b), (b, c), (c, a)):
            in_elem.add((u, v))
            in_elem.add((v, u))
    for i, up in pairs:
        assert st.station_of[up] != te_row
        assert st.side_node[up] == st.side_node[i]
        assert (int(i), int(up)) in in_elem
    # the upstream partners are the immediate upstream stations
    for i, up in pairs:
        assert st.xc[st.station_of[up]] > 0.9
        assert st.xc[st.station_of[up]] < st.xc[te_row]
