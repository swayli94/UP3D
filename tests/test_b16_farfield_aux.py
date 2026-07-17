"""
Track B / B16 gates: LS Newton far-field BC generalisation -- pin the
far-field-BOUNDARY aux DOFs on a Dirichlet far field (solve/picard_ls.py
`farfield_aux_dofs`, solve/newton_ls.py `farfield_aux="pin"`;
design_track_b.md §16, roadmap §B16).

Root cause (B9 recorded follow-up, measured GB16.1): a wake level set has no
outflow clip, so the sheet reaches the far-field boundary and the outer nodes
it crosses each carry an aux DOF governed only by a near-singular wake-LS row
on a giant outer tet. At a converged freestream Picard state those aux hold
garbage (|jump| 22-53 vs Gamma ~ 0.06); the Picard fixed point tolerates it but
the Newton residual reads it as an O(1) inconsistency (8 far-field MAIN rows
|R| ~ 84), which is why every committed LS Newton recipe churns on the
wing-body. Pinning the far-field aux to the branch value their host carries
(freestream: the same phi_inf; vortex: phi_inf +/- gamma by side) removes it.

  * GB16.2 inertness: farfield="neumann" is byte-identical pin vs legacy (the
    neumann outer aux are constrained by wake-LS there; no committed neumann
    run is perturbed). B12/B15 neumann anchors are the gated tier / demo.
  * GB16.5 Schur compat: the pinned aux leave a contiguous free-aux tail; the
    Schur split constructs with n_aux_expected = n_ext - n_pinned and still
    fails loudly on a genuine mis-count.
  * (gated) GB16.3 escape: the wing-body freestream Newton, which churns under
    "legacy", converges to machine residual with 0 limited under "pin".

The heavy wing-body M0.5 solves + the B12/B15 anchors live in the B16 demo
(cases/demo/b16_farfield_aux) under PYFP3D_TRANSONIC_GATES; the ungated tests
here run on the committed coarse 2.5D mesh (cheap).
"""

import inspect
import os
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.constraints.dirichlet import freestream_phi
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.solve.newton_ls import solve_multivalued_newton
from pyfp3d.solve.picard_ls import farfield_aux_dofs
from pyfp3d.solve.schur_ls import SchurReducedSystem
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

GATES = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"
REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
WB_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6_wingbody"
ALPHA = 2.0
WB_ALPHA, WB_M = 3.06, 0.5


def _naca(level="coarse"):
    path = M0_DIR / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    return read_mesh(path)


def _naca_mvop(mesh):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def _wingbody(level="coarse"):
    """The B9 wing-body LS mesh + operator (skips if the gitignored mesh is
    absent: regenerate with
    `python cases/meshes/onera_m6_wingbody/generate_onera_m6_wingbody.py
    --levels coarse medium`)."""
    path = WB_DIR / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored); regenerate with "
                    "cases/meshes/onera_m6_wingbody/"
                    "generate_onera_m6_wingbody.py --levels coarse medium")
    from pyfp3d.meshgen.fuselage import FuselageParams
    from pyfp3d.meshgen.wingbody import te_polyline
    mesh = read_mesh(str(path))
    a = np.radians(WB_ALPHA)
    wls = WakeLevelSet(te_polyline(FuselageParams()),
                       direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


# ---------------------------------------------------------------------------
# 1. the helper vs a brute-force set
# ---------------------------------------------------------------------------
def test_farfield_aux_dofs_bruteforce():
    """`farfield_aux_dofs` == {ext_dof_of_node[n] for n in far-field if it has
    an aux}, hosts aligned, and non-empty on the naca 2.5D wake (the sheet
    reaches the downstream far field)."""
    mesh = _naca()
    mvop = _naca_mvop(mesh)
    hosts, aux = farfield_aux_dofs(mesh, mvop.cm)
    ff = np.unique(mesh.boundary_faces["farfield"])
    want = {int(mvop.cm.ext_dof_of_node[n]) for n in ff
            if mvop.cm.ext_dof_of_node[n] >= 0}
    assert set(int(a) for a in aux) == want
    assert np.all(aux >= mvop.n_main) and np.all(hosts < mvop.n_main)
    assert np.all(mvop.cm.ext_dof_of_node[hosts] == aux)
    assert aux.size > 0                                   # sheet reaches outflow


# ---------------------------------------------------------------------------
# 2. the value rules, locked to machine (freestream jump=0, vortex jump=gamma)
# ---------------------------------------------------------------------------
def test_pin_value_rules_jump():
    """freestream pin => [phi]=0 on the ring; vortex pin => [phi]=gamma
    (side-signed), both exact by construction."""
    mesh = _naca()
    mvop = _naca_mvop(mesh)
    hosts, aux = farfield_aux_dofs(mesh, mvop.cm)
    side = mvop.cm.node_side[hosts]

    # freestream: main = aux = phi_inf(host) -> jump = side*(main-aux) = 0
    phi = np.zeros(mvop.n_total)
    phi[hosts] = freestream_phi(mesh.nodes[hosts], ALPHA, 1.0)
    phi[aux] = freestream_phi(mesh.nodes[hosts], ALPHA, 1.0)
    assert np.allclose(mvop.node_jump(phi, hosts), 0.0, atol=1e-14)

    # vortex: aux = main - side*gamma -> jump = side*(main-aux) = gamma
    gamma = 0.137
    phi = np.zeros(mvop.n_total)
    phi[hosts] = np.linspace(-1.0, 1.0, hosts.size)       # arbitrary main
    phi[aux] = phi[hosts] - side * gamma
    assert np.allclose(mvop.node_jump(phi, hosts), gamma, atol=1e-13)


# ---------------------------------------------------------------------------
# 3. GB16.2 cheap tier: neumann is byte-identical pin vs legacy
# ---------------------------------------------------------------------------
def test_neumann_bit_identical_pin_vs_legacy():
    """pin is inert on farfield='neumann' (the outer aux are constrained by the
    wake-LS rows there); the two states are byte-identical without needing to
    converge."""
    mesh = _naca()
    common = dict(m_inf=0.7, alpha_deg=ALPHA, farfield="neumann",
                  n_seed=5, n_newton_max=3)
    a = solve_multivalued_newton(mesh=mesh, mvop=_naca_mvop(mesh),
                                 farfield_aux="pin", **common)
    b = solve_multivalued_newton(mesh=mesh, mvop=_naca_mvop(mesh),
                                 farfield_aux="legacy", **common)
    assert np.array_equal(a["phi_ext"], b["phi_ext"])


# ---------------------------------------------------------------------------
# 4. the knob: default "pin", validated
# ---------------------------------------------------------------------------
def test_farfield_aux_knob():
    p = inspect.signature(solve_multivalued_newton).parameters
    assert p["farfield_aux"].default == "pin"
    mesh = _naca()
    with pytest.raises(ValueError, match="farfield_aux"):
        solve_multivalued_newton(mesh=mesh, mvop=_naca_mvop(mesh), m_inf=0.5,
                                 farfield_aux="bogus")


# ---------------------------------------------------------------------------
# 5. GB16.5: the Schur split accepts a pinned-aux free set, rejects mis-counts
# ---------------------------------------------------------------------------
def test_schur_split_with_pinned_aux():
    """Reducing the free set by the pinned far-field aux leaves a contiguous
    aux tail; SchurReducedSystem constructs with n_aux_expected = n_ext -
    n_pinned and raises on the un-adjusted (legacy) count."""
    mesh = _naca()
    mvop = _naca_mvop(mesh)
    hosts, aux = farfield_aux_dofs(mesh, mvop.cm)
    assert aux.size > 0
    # a freestream-state fused matrix (the raw material for the split)
    phi = np.zeros(mvop.n_total)
    phi[:mvop.n_main] = freestream_phi(mesh.nodes, ALPHA, 1.0)
    cut = np.flatnonzero(mvop.cm.ext_dof_of_node >= 0)
    phi[mvop.cm.ext_dof_of_node[cut]] = phi[cut]
    A = mvop.assemble_matrix(closure="wake_ls", te_kutta="pressure",
                             phi_ext=phi).tocsr()
    ff = np.unique(mesh.boundary_faces["farfield"])
    is_dir = np.zeros(mvop.n_total, dtype=bool)
    is_dir[ff] = True
    is_dir[aux] = True                                    # B16 pin
    free = np.flatnonzero(~is_dir)
    A_free = A[free][:, free]
    n_pin = aux.size
    schur = SchurReducedSystem(A_free, free, mvop.n_main,
                               n_aux_expected=mvop.n_ext - n_pin)
    assert schur.n_aux == mvop.n_ext - n_pin
    assert np.all(schur.main_free < mvop.n_main)
    with pytest.raises(ValueError, match="aux"):           # legacy count now wrong
        SchurReducedSystem(A_free, free, mvop.n_main,
                           n_aux_expected=mvop.n_ext)


# ---------------------------------------------------------------------------
# 6. vortex pin runs finite and applies the constant-jump rule
# ---------------------------------------------------------------------------
def test_vortex_pin_smoke_constant_jump():
    """A vortex-BC pin Newton step produces a finite state whose far-field ring
    jumps are all equal (the single-gamma value rule) -- the wiring proof."""
    mesh = _naca()
    mvop = _naca_mvop(mesh)
    hosts, aux = farfield_aux_dofs(mesh, mvop.cm)
    r = solve_multivalued_newton(mesh=mesh, mvop=mvop, m_inf=0.5,
                                 alpha_deg=ALPHA, farfield="vortex",
                                 farfield_aux="pin", n_seed=3, n_newton_max=1)
    assert np.all(np.isfinite(r["phi_ext"]))
    jumps = mvop.node_jump(r["phi_ext"], hosts)
    assert np.allclose(jumps, jumps[0], atol=1e-10)        # constant on the ring


# ---------------------------------------------------------------------------
# 7. wing-body: the pin wiring on the real mesh (cheap, ungated)
# ---------------------------------------------------------------------------
def test_wingbody_pin_wiring():
    """On the B9 wing-body coarse mesh the far-field aux are exactly the 4
    outflow-ring aux, and a single freestream pin Newton step drives their jump
    to 0 (freestream value rule) -- no full solve needed."""
    mesh, mvop = _wingbody("coarse")
    hosts, aux = farfield_aux_dofs(mesh, mvop.cm)
    assert aux.size >= 1                                   # sheet reaches outflow
    assert np.all(mesh.nodes[hosts, 0] > 5.0)             # genuinely far downstream
    r = solve_multivalued_newton(mesh=mesh, mvop=mvop, m_inf=WB_M,
                                 alpha_deg=WB_ALPHA, farfield="freestream",
                                 farfield_aux="pin", n_seed=3, n_newton_max=1)
    assert np.all(np.isfinite(r["phi_ext"]))
    assert np.allclose(mvop.node_jump(r["phi_ext"], hosts), 0.0, atol=1e-10)


# ---------------------------------------------------------------------------
# 8. GB16.1 premise lock: the legacy path leaves garbage far-field aux
# ---------------------------------------------------------------------------
def test_wingbody_legacy_junk_aux_present():
    """The B16 premise: under a legacy freestream Picard/Newton warm start the
    outer aux carry O(10) garbage jumps. If this ever stops reproducing, B16's
    root cause has changed and the phase should be re-examined."""
    mesh, mvop = _wingbody("coarse")
    npz = REPO_ROOT / "cases/demo/b9_wingbody/results/ls_coarse.npz"
    if not npz.exists():
        pytest.skip("b9 ls_coarse.npz cache absent (run cases/demo/b9_wingbody)")
    phi = np.load(npz)["phi_ext"]
    cut = np.flatnonzero(mvop.cm.ext_dof_of_node >= 0)
    x = mesh.nodes[cut, 0]
    jump = mvop.node_jump(phi, cut)
    assert np.abs(jump[x >= 10.0]).max() > 10.0            # garbage on the outer ring


# ---------------------------------------------------------------------------
# 9. GB16.3 (gated): the wing-body freestream Newton converges under pin
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not GATES, reason="wing-body M0.5 Newton is ~1 min "
                                      "(PYFP3D_TRANSONIC_GATES=1)")
def test_wingbody_freestream_pin_converges():
    """GB16.3 coarse: legacy churns (residual O(1), thousands of limited
    cells); pin reaches machine residual with 0 limited. The far-field BC layer
    is the fix. (A handful of floored cells may remain at the wing-fuselage
    junction -- the pre-existing B8 mixed-plain / G1.6 fuselage-Cp class, GB9.4;
    it is orthogonal to the BC fix and NOT asserted here.)"""
    mesh, mvop = _wingbody("coarse")
    kw = dict(m_inf=WB_M, alpha_deg=WB_ALPHA, farfield="freestream",
              n_seed=35, n_newton_max=40, tol_residual=1e-8, freeze_tol=1e-6)
    legacy = solve_multivalued_newton(mesh=mesh, mvop=_wingbody("coarse")[1],
                                      farfield_aux="legacy", **kw)
    pin = solve_multivalued_newton(mesh=mesh, mvop=mvop,
                                   farfield_aux="pin", **kw)
    assert legacy["residual_history"][-1] > 1.0           # legacy churns
    assert pin["residual_history"][-1] < 1e-6             # pin -> machine
    assert pin["n_limited"] == 0                          # limiter quiet
    hosts, _ = farfield_aux_dofs(mesh, mvop.cm)
    assert np.abs(mvop.node_jump(pin["phi_ext"], hosts)).max() < 1e-6
