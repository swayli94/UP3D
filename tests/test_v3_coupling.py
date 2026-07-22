"""Track V V3 loose coupling (binding: docs/roadmap/track_v.md "V3 -- Loose
coupling (2.5-D ladder + fuselage smoke)"; design: docs/design_track_v.md;
module under test: pyfp3d/viscous/coupling.py).

Covers (ungated, lightweight):
  - build_airfoil_case wiring on the coarse wake-cut NACA0012 strip: the
    station table is a closed loop (TE duplication collapses to one row),
    stagnation = min-x station, per-side transition flags at x_tr, LE band
    == inflow candidates == x/c <= le_band_x, per-node TE side split;
  - volume<->surface scatter/gather roundtrip via volume_node_of;
  - delta* == 0 -> m_dot == 0 -> EXACT zero wall RHS (the loose loop's
    first outer iteration must be a no-op on the FP side);
  - build_closed_body_case wiring on a synthetic icosphere shell (nose-pole
    inflow candidates, both-end stagnation bands, x_tr_frac flags);
  - coarse-mesh smoke: run_loose_coupling for two outer iterations stays
    finite, drives a nonzero transpiration RHS, and the IBL state is
    physical (delta* >= 0 everywhere after the floor, quasi-2D crossflow
    small -- the D13 SIII.D.1 local-basis lock on a curved wall).

Runs in both lanes: default JIT and PYFP3D_NOJIT=1 (the smoke leg is
skipped under NOJIT -- pure-Python numba fallback is too slow for the two
outer FP+IBL solves; the gate evidence runs JIT anyway).
"""

import os
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    build_airfoil_case,
    build_closed_body_case,
    make_picard_lifting_driver,
    run_loose_coupling,
    station_average,
)
from pyfp3d.viscous.transpiration import (
    assemble_transpiration_rhs,
    transpiration_from_delta_star,
)

from .mesh_utils import generate_sphere_shell_mesh, icosphere

REPO_ROOT = Path(__file__).parent.parent
NACA_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"

M_INF, ALPHA, RE = 0.5, 2.0, 3.0e6
NOJIT = os.environ.get("PYFP3D_NOJIT", "0") == "1"


@pytest.fixture(scope="module")
def airfoil_case():
    mc, wc = cut_wake(read_mesh(str(NACA_DIR / "coarse.msh")))
    cfg = CouplingConfig(re_chord=RE, m_inf=M_INF, alpha_deg=ALPHA)
    case = build_airfoil_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg
    )
    return mc, wc, cfg, case


# ---------------------------------------------------------------------------
# airfoil case wiring
# ---------------------------------------------------------------------------


def test_station_table_closed_loop(airfoil_case):
    """The wake-cut strip groups TE duplicates into ONE station row, so the
    station graph is a closed loop: order covers every station, starts at
    the min-x (stagnation) station, and both loop neighbours of the
    stagnation exist (le_nbrs)."""
    _, _, _, case = airfoil_case
    st = case.stations
    assert st is not None
    assert len(st.order) == len(st.xc)
    assert st.order[0] == st.stag_row
    assert st.xc[st.stag_row] == pytest.approx(np.min(st.xc))
    assert len(st.le_nbrs) == 2
    # every surface node maps to a valid station
    assert np.all((st.station_of >= 0) & (st.station_of < len(st.xc)))


def test_te_copies_one_row_two_sides(airfoil_case):
    """TE = max-x station row carries side 0 (split per node instead); its
    nodes split into upper/lower copies by side_node."""
    _, _, _, case = airfoil_case
    st = case.stations
    te_row = int(st.order[np.argmax(st.xy[st.order, 0])])
    assert st.side[te_row] == 0
    te_nodes = np.where(st.station_of == te_row)[0]
    assert len(te_nodes) >= 4  # 2 sides x >=2 span lines
    sides = set(np.sign(st.side_node[te_nodes]).tolist())
    assert sides == {-1, 1}
    # the two span lines of one side share (x, y) but are distinct nodes
    assert len(np.unique(st.side_node[te_nodes])) == 2


def test_transition_flags_and_le_band(airfoil_case):
    """D-TR: turbulent exactly on stations with x/c >= x_tr of their side;
    the stagnation station stays laminar; LE band == x/c <= le_band_x ==
    inflow candidates."""
    _, _, cfg, case = airfoil_case
    st = case.stations
    xc_n = st.xc[st.station_of]
    expect = np.zeros(case.sm.n_node, dtype=np.int64)
    expect[(st.side_node == 1) & (xc_n >= cfg.x_tr_upper)] = 1
    expect[(st.side_node == -1) & (xc_n >= cfg.x_tr_lower)] = 1
    np.testing.assert_array_equal(case.turbulent_flags, expect)
    assert np.all(case.turbulent_flags[st.station_of == st.stag_row] == 0)
    le = st.xc[st.station_of] <= cfg.le_band_x
    np.testing.assert_array_equal(case.le_band_surf, le)
    np.testing.assert_array_equal(case.inflow_candidates, le)
    assert np.all(case.seed_fetch >= 0.0)
    assert case.inflow_fetch > 0.0


def test_scatter_gather_roundtrip(airfoil_case):
    """volume -> surface gather and surface -> volume scatter through
    volume_node_of are mutually consistent on wall nodes."""
    mc, _, _, case = airfoil_case
    rng = np.random.default_rng(0)
    f_vol = rng.standard_normal(len(mc.nodes))
    f_surf = f_vol[case.sm.volume_node_of]
    back = np.zeros(len(mc.nodes))
    back[case.sm.volume_node_of] = f_surf
    wall_nodes = np.unique(case.wall_faces.reshape(-1))
    np.testing.assert_array_equal(back[wall_nodes], f_vol[wall_nodes])


def test_zero_delta_star_exact_zero_rhs(airfoil_case):
    """delta* == 0 -> transpiration m_dot == 0 -> the assembled wall RHS is
    the exact zero vector (first-iteration no-op identity, GV2.1(b)
    discipline carried through the V3 loop)."""
    mc, _, cfg, case = airfoil_case
    sm = case.sm
    ue = np.zeros((sm.n_node, 3))
    ue[:, 0] = 1.0
    rho = np.ones(sm.n_node)
    m_surf = transpiration_from_delta_star(sm, rho, ue, np.zeros(sm.n_node))
    assert np.array_equal(m_surf, np.zeros(sm.n_node))
    m_vol = np.zeros(len(mc.nodes))
    m_vol[sm.volume_node_of] = m_surf
    rhs = assemble_transpiration_rhs(mc.nodes, case.wall_faces, m_vol)
    assert np.array_equal(rhs, np.zeros(len(mc.nodes)))


# ---------------------------------------------------------------------------
# closed-body case wiring (synthetic icosphere shell)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def sphere_case():
    nodes, elements, wall_nodes, _ = generate_sphere_shell_mesh(
        subdivisions=1, n_layers=2, r_inner=1.0, r_outer=2.0, grading=1.5
    )
    _, faces = icosphere(1)
    cfg = CouplingConfig(re_chord=1.0e6, m_inf=0.3, alpha_deg=0.0)
    case = build_closed_body_case(nodes, elements, faces, cfg, x_tr_frac=0.3)
    return case


def test_closed_body_wiring(sphere_case):
    """Nose-pole inflow candidates (pole + one ring), transition at
    x_tr_frac of the body length, stagnation bands at BOTH ends, no
    station table."""
    case = sphere_case
    assert case.stations is None
    x = case.sm.xyz[:, 0]
    x0, x1 = float(np.min(x)), float(np.max(x))
    nose = int(np.argmin(x))
    assert case.inflow_candidates[nose]
    # all candidates sit in the nose band
    assert np.all(x[case.inflow_candidates] < x0 + 0.5 * (x1 - x0))
    # transition flags at 30% body length
    np.testing.assert_array_equal(
        case.turbulent_flags, (x >= x0 + 0.3 * (x1 - x0)).astype(np.int64)
    )
    # stagnation bands at both ends
    assert np.any(case.le_band_surf & (x < x0 + 0.05 * (x1 - x0)))
    assert np.any(case.le_band_surf & (x > x1 - 0.05 * (x1 - x0)))
    assert case.inflow_fetch > 0.0
    # the TAIL stagnation band is the Dirichlet outflow closure (a
    # closed surface has no natural outflow; the Newton Jacobian is
    # exactly singular without it -- GV3.3 2026-07-22). The pin covers
    # the whole band: a narrow pole+ring pin leaves the tail-cone BL
    # free to separate (Goldstein singularity; same diag), and the
    # pin's frozen seed delta* is masked out of the transpiration
    # source by run_loose_coupling.
    assert case.outflow_pin_surf is not None
    expect = x > x1 - 0.05 * (x1 - x0)
    np.testing.assert_array_equal(case.outflow_pin_surf, expect)


# ---------------------------------------------------------------------------
# coarse smoke: two outer iterations end to end
# ---------------------------------------------------------------------------


@pytest.mark.skipif(NOJIT, reason="two FP+IBL solves are JIT-lane only")
def test_smoke_two_outer_iterations(airfoil_case):
    """End-to-end smoke on the coarse strip (M0.5, alpha 2, Re 3e6, x_tr
    0.05/0.05): the loop converges (GV3.2 criterion, |ds change| < 1e-3)
    in well under 10 outer iterations at omega = 1.0, the transpiration
    RHS is nonzero, delta* is nonnegative and O(1e-3), and the crossflow
    lock on the curved wall holds (max|DS2|, max|CF2| << ds*/cf
    magnitudes -- the D13 SIII.D.1 local-basis fix; pre-fix values were
    25.9 / 0.15)."""
    mc, wc, cfg, case = airfoil_case
    cfg = CouplingConfig(
        re_chord=RE, m_inf=M_INF, alpha_deg=ALPHA, n_outer_max=10,
    )
    driver = make_picard_lifting_driver(mc, wc, M_INF, ALPHA)
    res = run_loose_coupling(driver, case, cfg)

    assert res.converged
    assert res.n_outer <= 10
    assert np.all(np.isfinite(res.phi))
    assert np.all(np.isfinite(res.U))
    assert np.max(np.abs(res.m_dot)) > 0.0
    assert np.all(res.delta_star >= 0.0)
    ds_max = float(np.max(res.outs[:, C.OUT_DS1]))
    assert 1.0e-4 < ds_max < 5.0e-2
    # crossflow lock on the curved wall (local-basis fix)
    assert np.max(np.abs(res.outs[:, C.OUT_DS2])) < 0.5 * ds_max
    cf_max = float(np.max(np.abs(res.outs[:, C.OUT_CF1])))
    assert np.max(np.abs(res.outs[:, C.OUT_CF2])) < 0.5 * cf_max
    # the inflow BC pins the whole stagnation band (not a single node):
    # a single pinned station leaves the nose split under-anchored (V3
    # debug 2026-07-22 -- the near-singular near-separation trap)
    assert res.history[1]["inflow_n_pinned"] >= 4
    # station averages are finite and positive downstream of transition
    ds_st = station_average(case, res.outs[:, C.OUT_DS1])
    assert np.all(np.isfinite(ds_st))
    assert np.max(ds_st) > 0.0
