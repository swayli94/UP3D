"""Track V V5 wing case wiring (binding: docs/roadmap/track_v.md GV5.0,
2026-07-22; module under test: pyfp3d/viscous/coupling.py::build_wing_case).

Covers (ungated, lightweight; on the gitignored M6 coarse mesh -- skips
when absent, M1 precedent):
  - local-x/c wiring: LE band (x/c <= inflow_band_x) == inflow candidates
    MINUS the tip mask; the tip band (z > z_tip*(1-frac), frac = 0.05 =
    the production tip_taper r_c = 0.05*b_semi, B32) is the pin+mask band
    (outflow_pin_surf) and is excluded from the inflow candidates;
  - per-side forced-transition flags at x_tr of the LOCAL chord (the
    taper makes a global-x transition cut wrong outboard);
  - the TE-duplicated wall gives distinct upper/lower TE copies (natural
    outflow on both sides, mid-span);
  - volume<->surface scatter/gather roundtrip + the delta* == 0 -> exact
    zero wall RHS identity (the loop's first-iteration no-op, GV2.1(b)).

No coupled solve here (the GV5.0 gate evidence lives in
cases/analysis/v5_m6_bridge/).
"""

import os
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at, x_le, x_te
from pyfp3d.viscous.coupling import CouplingConfig, build_wing_case
from pyfp3d.viscous.transpiration import (
    assemble_transpiration_rhs,
    transpiration_from_delta_star,
)

REPO_ROOT = Path(__file__).parent.parent
M6_COARSE = REPO_ROOT / "cases" / "meshes" / "onera_m6" / "coarse.msh"

M_INF, ALPHA, RE = 0.5, 3.06, 11.72e6 / 0.64607  # per meter (Re_MAC / MAC)
TIP_FRAC = 0.05


@pytest.fixture(scope="module")
def wing_case():
    if not M6_COARSE.exists():
        pytest.skip("onera_m6/coarse.msh not generated (gitignored)")
    mc, wc = cut_wake(read_mesh(str(M6_COARSE)))
    cfg = CouplingConfig(re_chord=RE, m_inf=M_INF, alpha_deg=ALPHA)
    case = build_wing_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg,
        x_le=x_le, chord_at=chord_at, tip_mask_frac=TIP_FRAC,
    )
    return mc, wc, cfg, case


def _local_xc(case):
    x, z = case.sm.xyz[:, 0], case.sm.xyz[:, 2]
    zc = np.clip(z, 0.0, B_SEMI)
    return (x - x_le(zc)) / chord_at(zc)


def test_tip_mask_band(wing_case):
    """outflow_pin_surf == the tip mask band z > 0.95*z_tip, and it is a
    thin resolved band (not empty, not the whole wing)."""
    _, _, _, case = wing_case
    z = case.sm.xyz[:, 2]
    z_tip = float(np.max(z))
    assert z_tip == pytest.approx(B_SEMI, rel=1.0e-9)  # flat-cap family
    expect = z > z_tip * (1.0 - TIP_FRAC)
    np.testing.assert_array_equal(case.outflow_pin_surf, expect)
    n_tip = int(case.outflow_pin_surf.sum())
    assert 10 <= n_tip <= 0.2 * case.sm.n_node


def test_inflow_band_excludes_tip(wing_case):
    """LE-band inflow candidates sit at x/c <= inflow_band_x AND outside
    the tip mask (the tip pin wins the overlap); laminar Blasius seed
    discipline is the airfoil's stagnation-band one."""
    _, _, cfg, case = wing_case
    xc_n = _local_xc(case)
    assert np.all(xc_n[case.inflow_candidates] <= cfg.inflow_band_x)
    assert not np.any(case.inflow_candidates & case.outflow_pin_surf)
    assert np.any(case.inflow_candidates)
    # the A4 recovery zone covers BOTH the LE band and the tip mask
    le = (xc_n < cfg.le_band_x) | case.outflow_pin_surf
    np.testing.assert_array_equal(case.le_band_surf, le)
    assert np.all(case.seed_fetch > 0.0)
    assert case.inflow_fetch > 0.0


def test_transition_flags_local_chord(wing_case):
    """D-TR per side at x_tr of the LOCAL chord; the LE stagnation band
    stays laminar."""
    _, _, cfg, case = wing_case
    xc_n = _local_xc(case)
    cent_y = case.sm.xyz[case.sm.triangles].mean(axis=1)[:, 1]
    ysum = np.zeros(case.sm.n_node)
    ycnt = np.zeros(case.sm.n_node)
    np.add.at(ysum, case.sm.triangles.reshape(-1), np.repeat(cent_y, 3))
    np.add.at(ycnt, case.sm.triangles.reshape(-1), 1.0)
    side = np.where(ysum / np.maximum(ycnt, 1.0) >= 0.0, 1, -1)
    expect = np.zeros(case.sm.n_node, dtype=np.int64)
    expect[(side == 1) & (xc_n >= cfg.x_tr_upper)] = 1
    expect[(side == -1) & (xc_n >= cfg.x_tr_lower)] = 1
    np.testing.assert_array_equal(case.turbulent_flags, expect)
    assert np.all(case.turbulent_flags[case.inflow_candidates] == 0)


def test_te_copies_distinct_midspan(wing_case):
    """Natural outflow on both TE lines: at a mid-span station the two TE
    copies (same (x, z), TE-duplicated by the wake cut) are distinct
    surface nodes with opposite side signs."""
    _, _, _, case = wing_case
    xyz = case.sm.xyz
    z0 = 0.5 * B_SEMI
    xe = x_te(z0)
    near = (np.abs(xyz[:, 2] - z0) < 0.02 * B_SEMI) & (
        np.abs(xyz[:, 0] - xe) < 0.02 * chord_at(z0))
    sel = np.where(near)[0]
    assert len(sel) >= 2
    ys = xyz[sel, 1]
    assert np.any(ys > 0.0) and np.any(ys < 0.0), (
        "upper and lower TE copies must both exist at mid-span"
    )


def test_scatter_gather_and_zero_rhs(wing_case):
    """volume->surface gather / surface->volume scatter consistency on
    wall nodes, and delta* == 0 -> EXACT zero transpiration RHS (the
    loop's first-iteration no-op, GV2.1(b) discipline)."""
    mc, _, _, case = wing_case
    rng = np.random.default_rng(0)
    f_vol = rng.standard_normal(len(mc.nodes))
    f_surf = f_vol[case.sm.volume_node_of]
    back = np.zeros(len(mc.nodes))
    back[case.sm.volume_node_of] = f_surf
    wall_nodes = np.unique(case.wall_faces.reshape(-1))
    np.testing.assert_array_equal(back[wall_nodes], f_vol[wall_nodes])

    sm = case.sm
    ue = np.zeros((sm.n_node, 3))
    ue[:, 0] = 1.0
    m_surf = transpiration_from_delta_star(
        sm, np.ones(sm.n_node), ue, np.zeros(sm.n_node))
    assert np.array_equal(m_surf, np.zeros(sm.n_node))
    m_vol = np.zeros(len(mc.nodes))
    m_vol[sm.volume_node_of] = m_surf
    rhs = assemble_transpiration_rhs(mc.nodes, case.wall_faces, m_vol)
    assert np.array_equal(rhs, np.zeros(len(mc.nodes)))
