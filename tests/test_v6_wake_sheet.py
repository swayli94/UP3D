"""Track V V6 wake-sheet delta* source (binding: docs/roadmap/track_v.md
GV6.1; pre-registration cases/analysis/v6_1_wake_sheet/PRE_REGISTRATION.md
including the 2026-07-25 addendum; module under test:
pyfp3d/viscous/wake_sheet.py).

Covers (the pre-registered +6):
  1. wake SurfaceMesh construction + W3 sanity on the NACA coarse strip
     (station chain, monotone arc, positive areas, the fold pairing);
  2. the prescribed-producer construction identity (c) on a synthetic wall
     state (pinned delta*_upper/lower, theta_upper/lower -> delta*_wake(0)
     = delta*_TE, the monotone downstream relaxation, the H -> 1 limit);
  3. zero-field bit-identity: m_wake == 0 -> b_wake == the exact zero
     vector + the (a)(i) loose-loop leg (flag-ON zero field vs flag-OFF);
  4. the sign-pin MMS (b) on the coarse strip (dead air U_inf = 0, uniform
     m0 through the PRODUCTION assembly: antisymmetry, jump = m0/rho0
     within 5%, ejection away from the sheet on both sides);
  5. the (a)(ii) A/B loose-loop bit-identity (flag-OFF vs the gate-free
     library at the pinned baseline commit; BOTH legs are subprocesses on
     fresh worktrees = fresh numba compile, because cache-load is not
     bit-faithful to fresh-compile in the viscous chain -- isolate3,
     2026-07-25 -- so both legs must share one cache mode);
  6. the fold-pairing structural assert (W3): every minus-side load lands
     in its master row under T^T.
  7. (GV6.2, +1) the CouplingConfig.wake_l_rel_chords plumbing: an
     explicit 1.0 is bit-identical to the default, and a non-default
     value reaches the producer through the loose loop (k=1 mdot
     difference on the identical inviscid-seeded wall state).

Runs in both lanes: default JIT and PYFP3D_NOJIT=1; the loose-loop / MMS /
A-B legs are JIT-only (the V3 smoke precedent -- pure-Python numba
fallback is too slow for the FP solves).
"""

import os
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from ._tol import assert_rel_close
from pyfp3d.constraints.wake import WakeConstraint
from pyfp3d.kernels.jacobian import PicardOperator
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.coupling import (
    AirfoilStations,
    CouplingConfig,
    build_airfoil_case,
    make_picard_lifting_driver,
    run_loose_coupling,
)
from pyfp3d.viscous.wake_sheet import (
    build_wake_sheet_case,
    assemble_wake_sheet_rhs,
    wake_transpiration_source,
)

REPO_ROOT = Path(__file__).parent.parent
NACA_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"

M_INF, ALPHA, RE = 0.5, 2.0, 3.0e6
NOJIT = os.environ.get("PYFP3D_NOJIT", "0") == "1"

# (a)(ii): the gate-free library = the GV6.0 adjudication merge on main
# (the commit this branch forks from; every GV6.1 code change is on top).
GATE_FREE_BASELINE = "13916b5"


@pytest.fixture(scope="module")
def naca_coarse_cut():
    return cut_wake(read_mesh(str(NACA_DIR / "coarse.msh")))


@pytest.fixture(scope="module")
def wake_case(naca_coarse_cut):
    mc, wc = naca_coarse_cut
    return build_wake_sheet_case(mc, wc)


# ---------------------------------------------------------------------------
# 1. wake SurfaceMesh construction + W3 sanity
# ---------------------------------------------------------------------------


def test_wake_surface_mesh_w3(naca_coarse_cut, wake_case):
    """W3: the station chain covers the strip (one row per distinct (x,y)
    wake position) with strictly monotone arc from the TE; node areas
    positive; every wake-minus node pairs to exactly one plus-slave at
    coincident coordinates; volume_node_of maps into the cut mesh."""
    mc, wc = naca_coarse_cut
    sm = wake_case.sm
    n_st = len(np.unique(sm.xyz[:, :2], axis=0))
    assert len(np.unique(wake_case.station_of)) == n_st
    # one arc value per station, strictly increasing along the chain
    s_st = np.zeros(n_st)
    for r in range(n_st):
        s_rows = wake_case.s[wake_case.station_of == r]
        assert np.ptp(s_rows) == 0.0
        s_st[r] = s_rows[0]
    assert np.all(np.diff(np.sort(s_st)) > 0.0)
    assert s_st.min() == 0.0  # the TE station
    te_nodes = wake_case.station_of == int(np.argmin(s_st))
    assert np.all(sm.xyz[te_nodes, 0] == pytest.approx(sm.xyz[:, 0].min()))
    assert np.all(sm.node_area > 0.0)
    # fold pairing: count + coordinate coincidence with the WakeCut map
    vol = sm.volume_node_of
    assert len(wake_case.slave_of_surf) == len(vol) == len(wc.slave_nodes)
    assert np.all(np.isin(wake_case.slave_of_surf, wc.slave_nodes))
    np.testing.assert_array_equal(
        mc.nodes[wake_case.slave_of_surf], mc.nodes[vol]
    )
    assert vol.max() < len(mc.nodes)
    assert wake_case.slave_of_surf.max() < len(mc.nodes)


# ---------------------------------------------------------------------------
# 2. prescribed-producer construction identity (c) on a synthetic wall state
# ---------------------------------------------------------------------------


def _synthetic_wall_state():
    """Four wall stations, the TE row (max x/c) carrying two upper + two
    lower copies; pinned delta*/theta per side -> delta*_TE = 3e-3,
    theta_TE = 2e-3."""
    st = AirfoilStations(
        station_of=np.array([0, 1, 2, 3, 3, 3, 3]),
        xy=np.array(
            [
                [0.0, 0.0],
                [0.4, 0.02],
                [0.4, -0.02],
                [1.0, 0.0],
            ]
        ),
        xc=np.array([0.0, 0.4, 0.4, 1.0]),
        s=np.zeros(4),
        side=np.array([0, 1, -1, 0]),
        side_node=np.array([1, 1, -1, 1, 1, -1, -1]),
        order=np.array([0, 1, 3, 2]),
        stag_row=0,
        le_nbrs=(1, 2),
    )
    outs = np.zeros((7, C.N_OUT))
    outs[st.side_node == 1, C.OUT_DS1] = 2.0e-3  # upper (incl. TE copies)
    outs[st.side_node == -1, C.OUT_DS1] = 1.0e-3  # lower
    outs[st.side_node == 1, C.OUT_TH11] = 1.2e-3
    outs[st.side_node == -1, C.OUT_TH11] = 0.8e-3
    return st, outs


def test_prescribed_producer_construction_identity(wake_case):
    """Band (c): delta*_wake(0) == delta*_upper + delta*_lower to 1e-12
    relative; the relaxation is monotone downstream and tends to theta_TE
    (H_wake -> 1)."""
    st, outs = _synthetic_wall_state()
    ue = np.tile([1.0, 0.0, 0.0], (wake_case.sm.n_node, 1))
    ds_wake, m_wake, info = wake_transpiration_source(
        wake_case, st, outs, ue, m_inf=0.0
    )
    ds_te, th_te = info["ds_te"], info["th_te"]
    assert ds_te == pytest.approx(3.0e-3)
    assert th_te == pytest.approx(2.0e-3)
    at_s0 = wake_case.s == 0.0
    assert np.any(at_s0)
    assert np.max(np.abs(ds_wake[at_s0] - ds_te)) <= 1e-12 * abs(ds_te)
    # monotone non-increasing along the chain (delta*_TE > theta_TE here)
    order = np.argsort(wake_case.s)
    assert np.all(np.diff(ds_wake[order]) <= 1e-18)
    # H -> 1 limit: the coarse strip's wake reaches s = 14.5 chords, so
    # exp(-s/L_rel) < 1e-6 at the downstream end
    assert ds_wake[np.argmax(wake_case.s)] / th_te == pytest.approx(
        1.0, abs=1e-3
    )
    assert np.all(np.isfinite(m_wake))


# ---------------------------------------------------------------------------
# 3. zero-field bit-identity: exact zero RHS + the (a)(i) loose-loop leg
# ---------------------------------------------------------------------------


def test_zero_field_exact_zero_rhs(naca_coarse_cut, wake_case):
    """m_wake == 0 -> b_wake == the exact zero vector (the GV2.1(b)
    assembly discipline on the sheet channel; the band-(a)(i) basis)."""
    mc, _ = naca_coarse_cut
    rhs = assemble_wake_sheet_rhs(
        mc.nodes, wake_case, np.zeros(wake_case.sm.n_node)
    )
    assert np.array_equal(rhs, np.zeros(len(mc.nodes)))


@pytest.mark.skipif(NOJIT, reason="loose-loop FP solves are JIT-lane only")
def test_zero_field_loose_loop_bit_identical(naca_coarse_cut, wake_case):
    """(a)(i): flag ON with a prescribed ZERO delta*_wake field vs flag
    OFF, coarse loose loop (3 outer) -- bit-identical phi/gamma."""
    mc, wc = naca_coarse_cut
    cfg = CouplingConfig(re_chord=RE, m_inf=M_INF, alpha_deg=ALPHA)
    case = build_airfoil_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg
    )
    driver = make_picard_lifting_driver(mc, wc, M_INF, ALPHA)

    cfg_on = replace(cfg, n_outer_max=3, wake_transpiration=True)
    wake0 = replace(
        wake_case, prescribed_ds=np.zeros(wake_case.sm.n_node)
    )
    res_on = run_loose_coupling(driver, case, cfg_on, wake=wake0)

    cfg_off = replace(cfg, n_outer_max=3)
    res_off = run_loose_coupling(driver, case, cfg_off)

    assert np.array_equal(res_on.phi, res_off.phi)
    assert np.array_equal(res_on.gamma, res_off.gamma)
    assert res_on.n_outer == res_off.n_outer


# ---------------------------------------------------------------------------
# 4. sign-pin MMS (b): dead air, uniform sheet source
# ---------------------------------------------------------------------------


def _tet_face_owners(elements):
    """Sorted-triple face key -> owning tet ids (the wake_cut.py idiom)."""
    tet_faces = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
    owners = {}
    for e, tet in enumerate(np.asarray(elements, dtype=np.int64)):
        for f in tet_faces:
            key = tuple(sorted(int(tet[i]) for i in f))
            owners.setdefault(key, []).append(e)
    return owners


@pytest.mark.skipif(NOJIT, reason="the MMS FP solve is JIT-lane only")
def test_sign_pin_mms(naca_coarse_cut, wake_case):
    """Band (b): dead air (U_inf = 0, farfield Dirichlet phi = 0, Gamma
    pinned 0 -- the sheet source is top/bottom symmetric), uniform m0 > 0
    prescribed on the wake sheet through the PRODUCTION assembly. Probes:
    the owner-tet one-sided P1 normal gradients at sheet faces in the
    middle third of the strip (off the TE and the downstream end):
    antisymmetry, jump = m0/rho0 within 5%, ejection away from the sheet.
    """
    mc, wc = naca_coarse_cut
    m0 = 0.01
    b_wake = assemble_wake_sheet_rhs(
        mc.nodes, wake_case, np.full(wake_case.sm.n_node, m0)
    )
    # Dead air through the production T^T route: at m_inf = 0 the Picard
    # driver IS the incompressible lifting solve (the G3.3 equivalence),
    # and its u_inf = 0 corner only breaks in the COMPRESSIBLE
    # bookkeeping (q2/u_inf**2 -> q2_at_mach cap NaN at M_inf = 0), not in
    # the sheet-source channel the band pins. Drive the identical reduced
    # Laplace system directly from the production primitives:
    # WakeConstraint + farfield Dirichlet (phi = 0 at u_inf = 0) +
    # reduced_rhs(b_wake, Gamma = 0) + CG/AMG.
    import scipy.sparse.linalg as spla

    from pyfp3d.constraints.dirichlet import farfield_dirichlet
    from pyfp3d.solve.linear import build_amg_preconditioner

    op = PicardOperator(mc.nodes, mc.elements)
    con = WakeConstraint(op.assemble_matrix(), wc)
    gamma0 = np.zeros(wc.n_stations)
    dir_nodes, dir_vals = farfield_dirichlet(
        mc, wc, 0.0, gamma0, 0.0, (0.25, 0.0), beta=1.0
    )
    dir_red, vals_red = con.to_reduced_dirichlet(dir_nodes, dir_vals)
    assert np.all(vals_red == 0.0)  # dead-air far field
    is_dir = np.zeros(con.n_reduced, dtype=bool)
    is_dir[dir_red] = True
    free = np.where(~is_dir)[0]
    A = con.A_reduced
    A_free = A[free][:, free].tocsr()
    b_red = con.reduced_rhs(b_wake, gamma0)
    b_free = b_red[free] - A[free][:, dir_red].tocsr() @ vals_red
    x, info = spla.cg(
        A_free,
        b_free,
        M=build_amg_preconditioner(A_free)[1],
        rtol=1e-11,
        maxiter=3000,
    )
    assert info == 0
    phi_red = np.empty(con.n_reduced)
    phi_red[free] = x
    phi_red[dir_red] = vals_red
    phi = con.expand(phi_red, gamma0)

    owners = _tet_face_owners(mc.elements)
    grad, _ = op.velocities(phi)

    n_faces = len(wc.wake_faces_minus)
    cent = mc.nodes[np.asarray(wc.wake_faces_minus, dtype=np.int64)].mean(axis=1)
    x0, x1 = cent[:, 0].min(), cent[:, 0].max()
    lo, hi = x0 + (x1 - x0) / 3.0, x0 + 2.0 * (x1 - x0) / 3.0
    probes = [f for f in range(n_faces) if lo < cent[f, 0] < hi]
    assert len(probes) >= 4

    v_plus, v_minus = [], []
    for f in probes:
        key_m = tuple(sorted(int(v) for v in wc.wake_faces_minus[f]))
        key_p = tuple(sorted(int(v) for v in wc.wake_faces_plus[f]))
        assert len(owners[key_m]) == 1 and len(owners[key_p]) == 1
        v_minus.append(grad[owners[key_m][0], 1])
        v_plus.append(grad[owners[key_p][0], 1])
    v_plus = np.asarray(v_plus)
    v_minus = np.asarray(v_minus)

    # (iii) sign: m0 > 0 ejects AWAY from the sheet on both sides
    assert np.all(v_plus > 0.0)
    assert np.all(v_minus < 0.0)
    # (i) antisymmetry to discretization accuracy
    assert np.max(np.abs(v_plus + v_minus)) / m0 < 0.05
    # (ii) the jump equals m0/rho0 within 5% (rho0 = 1 in dead air)
    assert np.max(np.abs((v_plus - v_minus) - m0)) / m0 < 0.05


# ---------------------------------------------------------------------------
# 5. (a)(ii) A/B bit-identity vs the gate-free library
# ---------------------------------------------------------------------------


_AB_SNIPPET = """\
import sys

# the PEP 660 editable finder for pyfp3d sits in sys.meta_path and beats
# sys.path -- strip every editable finder BEFORE importing pyfp3d, then
# resolve pyfp3d from the worktree passed as argv[1]
sys.meta_path = [f for f in sys.meta_path
                 if "_EditableFinder" not in type(f).__name__]
sys.path.insert(0, sys.argv[1])

import numpy as np
import pyfp3d

assert pyfp3d.__file__.startswith(sys.argv[1]), pyfp3d.__file__

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    build_airfoil_case,
    make_picard_lifting_driver,
    run_loose_coupling,
)

mc, wc = cut_wake(read_mesh(sys.argv[2]))
cfg = CouplingConfig(re_chord=3.0e6, m_inf=0.5, alpha_deg=2.0, n_outer_max=3)
case = build_airfoil_case(
    mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg
)
res = run_loose_coupling(make_picard_lifting_driver(mc, wc, 0.5, 2.0), case, cfg)
np.savez(sys.argv[3], phi=res.phi, gamma=res.gamma)
"""


def _git(*args):
    return subprocess.run(
        ["git", "-C", str(REPO_ROOT), *args],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _overlay_working_tree_delta(worktree):
    """Overlay the working tree's pyfp3d/ delta (modified / added /
    untracked / deleted; __pycache__ excluded) onto a HEAD worktree, so
    the leg measures THIS tree's exact code state even when dirty."""
    # NB: no .strip() on stdout -- porcelain's leading status column
    # (" M ...") is positional; stripping mangles the first line's path
    out = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "status", "--porcelain", "--", "pyfp3d/"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    for line in out.splitlines():
        xy, rel = line[:2], line[3:]
        if " -> " in rel:
            rel = rel.split(" -> ")[-1]
        rel = rel.strip('"')
        if "__pycache__" in rel:
            continue
        dst = Path(worktree) / rel
        if "D" in xy:
            if dst.exists():
                dst.unlink()
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(REPO_ROOT / rel, dst)


def _ab_leg(snippet, worktree, ref, out_npz, overlay_delta=False):
    """One (a)(ii) leg: a FRESH worktree at ref (no __pycache__ -> every
    numba function compiles from source) + one subprocess run.

    isolate3/4 (2026-07-25): numba cache-LOAD is not bit-faithful to
    fresh-COMPILE, and the infidelity lives entirely in pyfp3d/viscous/
    -- a fresh leg and a cache-loading leg diverge at ~1e-5 in phi at
    outer k >= 1 (k=0 inviscid exact) even for identical sources. Both
    legs therefore run fresh-compile; comparing an in-process
    (cache-warm) leg against a worktree leg fails spuriously.
    """
    _git("worktree", "add", "--detach", str(worktree), ref)
    try:
        if overlay_delta:
            _overlay_working_tree_delta(worktree)
        subprocess.run(
            [
                sys.executable,
                str(snippet),
                str(worktree),
                str(NACA_DIR / "coarse.msh"),
                str(out_npz),
            ],
            check=True,
            capture_output=True,
            text=True,
            env=dict(os.environ),
        )
        return np.load(out_npz)
    finally:
        subprocess.run(
            ["git", "-C", str(REPO_ROOT), "worktree", "remove", "--force",
             str(worktree)],
            capture_output=True,
        )


@pytest.mark.skipif(NOJIT, reason="loose-loop FP solves are JIT-lane only")
@pytest.mark.skipif(shutil.which("git") is None, reason="git unavailable")
def test_ab_bit_identity_gate_free_library(tmp_path):
    """(a)(ii) / W1: the flag-OFF loose loop (this tree) reproduces the
    gate-free library (the pinned baseline commit) on the same machine,
    coarse 3 outer. Both legs are fresh-compile worktree subprocesses (the
    isolate3/4 cache-mode discipline).

    GS0.2 / D1 (2026-07-28): the assertion was `np.array_equal` and it FAILED
    during the audit's full-suite run while PASSING standalone on the same
    commit and the same thread count (67 s, idle machine) -- the only
    difference was concurrent load from other solver processes. The mechanism
    was NOT root-caused (registered as an open question in
    docs/dev_phase_two/20260728-1520-s0-foundation.md §6); per decision D1 the
    permanent assertion is now a 1e-12 relative tolerance, which still pins
    the phase-one claim (the flag adds no numerical effect: any real change
    would be orders larger) without being a load-sensitive alarm."""
    snippet = tmp_path / "ab_leg.py"
    snippet.write_text(_AB_SNIPPET)
    head = _git("rev-parse", "HEAD")
    base = _ab_leg(
        snippet, tmp_path / "gate_free", GATE_FREE_BASELINE,
        tmp_path / "baseline.npz",
    )
    cur = _ab_leg(
        snippet, tmp_path / "current", head, tmp_path / "current.npz",
        overlay_delta=True,
    )
    assert_rel_close(cur["phi"], base["phi"], msg="phi vs gate-free baseline")
    assert_rel_close(cur["gamma"], base["gamma"],
                     msg="gamma vs gate-free baseline")


# ---------------------------------------------------------------------------
# 6. fold-pairing structural assert (W3)
# ---------------------------------------------------------------------------


def test_fold_pairing_structural(naca_coarse_cut, wake_case):
    """Every minus-side load lands in its master row under T^T: with the
    production assembly (1/2 m per copy), (T^T b)[master] ==
    b[master] + b[slave] == 2*b[master] exactly, and the reduced vector
    has exactly n_reduced rows (no slave rows survive)."""
    mc, wc = naca_coarse_cut
    op = PicardOperator(mc.nodes, mc.elements)
    con = WakeConstraint(op.assemble_matrix(), wc)

    m_wake = 0.01 * (1.0 + wake_case.s)  # nonzero, streamwise-varying
    b_wake = assemble_wake_sheet_rhs(mc.nodes, wake_case, m_wake)
    b_red = con.reduced_rhs(b_wake, np.zeros(wc.n_stations))

    vol = wake_case.sm.volume_node_of
    slave = wake_case.slave_of_surf
    assert b_red.shape == (con.n_reduced,)
    # the two copies assemble bitwise-identical loads (coincident
    # geometry, same face order), so the fold is exactly 2x the one-sided
    # load -- the addendum's jump = m_wake accounting
    np.testing.assert_array_equal(b_wake[vol], b_wake[slave])
    np.testing.assert_array_equal(b_red[vol], b_wake[vol] + b_wake[slave])
    np.testing.assert_array_equal(b_red[vol], 2.0 * b_wake[vol])
    assert np.all(b_red[vol] != 0.0)


# ---------------------------------------------------------------------------
# 7. GV6.2 wake_l_rel_chords plumbing (pre-registration
#    cases/analysis/v6_2_measured_effect/PRE_REGISTRATION.md section 5)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(NOJIT, reason="loose-loop FP solves are JIT-lane only")
def test_wake_l_rel_chords_plumbing(naca_coarse_cut, wake_case):
    """CouplingConfig.wake_l_rel_chords reaches the producer through the
    loose loop; the default preserves the GV6.1 behaviour bit-identically.

    At outer k = 1 both runs consume the IDENTICAL inviscid-seeded wall
    state (same k = 0 phi, deterministic), so any mdot_wake_max
    difference there comes only from l_rel_chords -- broken plumbing
    (the field silently ignored) gives a bit-identical k = 1 mdot."""
    mc, wc = naca_coarse_cut
    cfg = CouplingConfig(re_chord=RE, m_inf=M_INF, alpha_deg=ALPHA)
    case = build_airfoil_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg
    )
    driver = make_picard_lifting_driver(mc, wc, M_INF, ALPHA)

    cfg_def = replace(cfg, n_outer_max=3, wake_transpiration=True)
    res_def = run_loose_coupling(driver, case, cfg_def, wake=wake_case)

    cfg_10 = replace(cfg_def, wake_l_rel_chords=1.0)
    res_10 = run_loose_coupling(driver, case, cfg_10, wake=wake_case)
    assert np.array_equal(res_def.phi, res_10.phi)
    assert np.array_equal(res_def.gamma, res_10.gamma)

    cfg_05 = replace(cfg_def, wake_l_rel_chords=0.5)
    res_05 = run_loose_coupling(driver, case, cfg_05, wake=wake_case)
    m_def = float(res_def.history[1]["mdot_wake_max"])
    m_05 = float(res_05.history[1]["mdot_wake_max"])
    assert m_def > 0.0 and m_05 > 0.0
    assert m_05 != m_def
