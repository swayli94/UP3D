"""
Track B / B14 gates: `precond="schur"` -- Schur-eliminated aux block + AMG on
the SPD Picard main block (solve/schur_ls.py; design_track_b.md §5.3, roadmap
§B14).

Motivation (B11 measured): the fused wake_ls matrix is nonsymmetric only
because of the aux rows, and the SPD spring surrogate's jump==0 prior kills
the circulation mode (gamma 0.0033 vs 0.139) while ILU diverges at 2.5D medium
(gamma = -136.99). B14 eliminates the aux thin-strip block EXACTLY per step
(`lu_aa = splu(J_aa)`) and runs GMRES on the reduced main-free operator with
AMG on the SPD single-valued Picard block -- no springs, no full-size
factorization.

  * GB14.1 diagnostic-first: J_aa factors and its measured 1-norm condition
    estimate is finite ("measure, don't assume" -- the pre-registered gate).
  * GB14.2 correctness: schur reaches the spsolve gamma to |dgamma| < 1e-8 on
    the coarse lifting + Newton solves, 0 stalls / 0 fallbacks; one linear
    solve agrees with spsolve and satisfies the aux rows to ~machine (the
    honesty identity: reduced rtol IS the full-system main-row residual).
  * GB14.5 inertness: the default precond=None paths are byte-identical and
    the new counter stays 0.
  * (gated) GB14.3 discriminating tier: 2.5D MEDIUM lifting -- where ILU
    DIVERGED -- converges to the spsolve gamma. "Passing there is what a real
    escape means" (roadmap B14).

All ungated tests run on the committed coarse 2.5D mesh (cheap). M6-scale
capability + timing (GB14.4) is the demo's job (cases/demo/b14_schur_precond),
not pytest.
"""

import inspect
import os
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse.linalg as spla

from pyfp3d.constraints.dirichlet import freestream_phi
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.solve.newton_ls import solve_multivalued_newton
from pyfp3d.solve.picard_ls import solve_multivalued_lifting
from pyfp3d.solve.schur_ls import (
    SchurReducedSystem,
    jaa_diagnostic,
    main_block_preconditioner,
)
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

GATES = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"
REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
ALPHA = 2.0

NEWTON_COMMON = dict(m_inf=0.7, alpha_deg=ALPHA, farfield="neumann",
                     n_seed=20, n_newton_max=25)


def _mesh(level="coarse"):
    path = M0_DIR / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    return read_mesh(path)


def _mvop(mesh):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def _fused_free_system(mesh, mvop):
    """The freestream-state wake_ls matrix, free-reduced (Dirichlet = the
    far field), plus the free index set -- the raw material for the split."""
    phi_ext = np.zeros(mvop.n_total)
    phi_ext[:mvop.n_main] = freestream_phi(mesh.nodes, ALPHA, 1.0)
    cut_nodes = np.flatnonzero(mvop.cm.ext_dof_of_node >= 0)
    phi_ext[mvop.cm.ext_dof_of_node[cut_nodes]] = phi_ext[cut_nodes]
    A = mvop.assemble_matrix(closure="wake_ls", te_kutta="pressure",
                             phi_ext=phi_ext).tocsr()
    ff = np.unique(mesh.boundary_faces["farfield"])
    is_dir = np.zeros(mvop.n_total, dtype=bool)
    is_dir[ff] = True
    free = np.flatnonzero(~is_dir)
    return A[free][:, free], free


# ---------------------------------------------------------------------------
# Wiring / inertness (GB14.5).
# ---------------------------------------------------------------------------

def test_schur_accepted_and_bogus_rejected():
    p = inspect.signature(solve_multivalued_newton).parameters
    assert p["precond"].default is None
    p = inspect.signature(solve_multivalued_lifting).parameters
    assert p["precond"].default is None
    mesh = _mesh()
    with pytest.raises(ValueError, match="schur"):
        solve_multivalued_lifting(_mvop(mesh), mesh, 0.5, alpha_deg=ALPHA,
                                  precond="bogus")
    with pytest.raises(ValueError, match="schur"):
        solve_multivalued_newton(mesh=mesh, mvop=_mvop(mesh),
                                 precond="bogus", **NEWTON_COMMON)


def test_default_paths_bit_identical():
    """GB14.5: the B14 edits leave precond=None byte-identical, with the new
    counter inert at 0."""
    mesh = _mesh()
    a = solve_multivalued_lifting(_mvop(mesh), mesh, 0.5, alpha_deg=ALPHA)
    b = solve_multivalued_lifting(_mvop(mesh), mesh, 0.5, alpha_deg=ALPHA)
    assert np.array_equal(a["phi_ext"], b["phi_ext"])
    assert a["n_schur_fallback"] == 0
    a = solve_multivalued_newton(mesh=mesh, mvop=_mvop(mesh), **NEWTON_COMMON)
    b = solve_multivalued_newton(mesh=mesh, mvop=_mvop(mesh), **NEWTON_COMMON)
    assert np.array_equal(a["phi_ext"], b["phi_ext"])
    assert a["n_schur_fallback"] == 0 and a["n_gmres_total"] == 0


# ---------------------------------------------------------------------------
# The split itself + the honesty identity (GB14.2 linear part).
# ---------------------------------------------------------------------------

def test_schur_block_split_composes():
    """The four blocks recompose the free matrix exactly, and the split sits
    where the DOF layout says it must (aux = the contiguous tail)."""
    mesh = _mesh()
    mvop = _mvop(mesh)
    A_free, free = _fused_free_system(mesh, mvop)
    schur = SchurReducedSystem(A_free, free, mvop.n_main,
                               n_aux_expected=mvop.n_ext)
    assert schur.n_aux == mvop.n_ext
    assert np.all(schur.main_free < mvop.n_main)
    rng = np.random.default_rng(0)
    x = rng.standard_normal(free.size)
    blocks = np.concatenate([
        schur.J_mm @ x[:schur.n_mf] + schur.J_ma @ x[schur.n_mf:],
        schur.J_am @ x[:schur.n_mf] + schur.J_aa @ x[schur.n_mf:]])
    full = A_free @ x
    assert np.linalg.norm(full - blocks) <= 1e-14 * np.linalg.norm(full)
    # a wrong n_aux expectation fails loudly, never silently mis-splits
    with pytest.raises(ValueError, match="aux"):
        SchurReducedSystem(A_free, free, mvop.n_main,
                           n_aux_expected=mvop.n_ext + 1)


def test_schur_linear_exactness_and_aux_rows():
    """One linear solve at tight rtol agrees with spsolve, and the aux rows
    are satisfied to ~machine by back-substitution -- the honesty identity
    that makes the reduced rtol THE full-system main-row residual."""
    mesh = _mesh()
    mvop = _mvop(mesh)
    A_free, free = _fused_free_system(mesh, mvop)
    schur = SchurReducedSystem(A_free, free, mvop.n_main,
                               n_aux_expected=mvop.n_ext)
    rng = np.random.default_rng(1)
    b = rng.standard_normal(free.size)
    M = main_block_preconditioner(mvop, None, schur.main_free)
    x, _, info = schur.solve(b, M=M, rtol=1e-12, maxiter=50)
    assert info == 0
    x_ref = spla.spsolve(A_free.tocsc(), b)
    nb = np.linalg.norm(b)
    assert np.linalg.norm(x - x_ref) <= 1e-8 * np.linalg.norm(x_ref)
    assert np.linalg.norm((A_free @ x - b)[schur.n_mf:]) <= 1e-10 * nb


def test_jaa_diagnostic_finite():
    """GB14.1 pre-registered diagnostic: the aux strip factors and its 1-norm
    condition estimate is finite -- measured, not assumed (the constant-jump
    null vector mixes main+aux columns, so J_aa is generically nonsingular)."""
    mesh = _mesh()
    mvop = _mvop(mesh)
    A_free, free = _fused_free_system(mesh, mvop)
    schur = SchurReducedSystem(A_free, free, mvop.n_main,
                               n_aux_expected=mvop.n_ext)   # splu succeeded
    d = jaa_diagnostic(schur)
    assert d["n_aux"] == mvop.n_ext
    assert np.isfinite(d["cond1"]) and d["cond1"] < 1e14


# ---------------------------------------------------------------------------
# GB14.2: nonlinear correctness vs the spsolve reference (coarse, cheap).
# ---------------------------------------------------------------------------

def test_schur_lifting_coarse_matches():
    """Coarse M0.5 lifting -- the exact case where the B11 spring surrogate
    STALLED to gamma 0.0033: schur must land on the spsolve gamma."""
    mesh = _mesh()
    ref = solve_multivalued_lifting(_mvop(mesh), mesh, 0.5, alpha_deg=ALPHA)
    s = solve_multivalued_lifting(_mvop(mesh), mesh, 0.5, alpha_deg=ALPHA,
                                  precond="schur")
    assert ref["converged"] and s["converged"]
    assert abs(s["gamma"] - ref["gamma"]) < 1e-8
    assert s["n_gmres_stalled"] == 0 and s["n_schur_fallback"] == 0


def test_schur_newton_coarse_matches():
    """Coarse M0.7 Newton (supersonic pocket: Terms 2/3 live in J_mm and are
    invisible to the Term-1 AMG -- schur must still converge unstalled)."""
    mesh = _mesh()
    ref = solve_multivalued_newton(mesh=mesh, mvop=_mvop(mesh),
                                   **NEWTON_COMMON)
    s = solve_multivalued_newton(mesh=mesh, mvop=_mvop(mesh),
                                 precond="schur", **NEWTON_COMMON)
    assert ref["converged"] and s["converged"]
    assert abs(s["gamma"] - ref["gamma"]) < 1e-8
    assert s["n_limited"] == 0 and s["n_floored"] == 0
    assert s["n_gmres_stalled"] == 0 and s["n_schur_fallback"] == 0


def test_schur_seed_precond_smoke():
    """seed_precond="schur" reaches the Picard seed through _seed_from_picard
    (wiring check, deliberately tiny budgets)."""
    mesh = _mesh()
    r = solve_multivalued_newton(mesh=mesh, mvop=_mvop(mesh), m_inf=0.5,
                                 alpha_deg=ALPHA, seed_precond="schur",
                                 n_seed=5, n_newton_max=3)
    assert np.all(np.isfinite(r["phi_ext"]))


# ---------------------------------------------------------------------------
# GB14.3 (gated): the discriminating tier -- 2.5D MEDIUM lifting, ILU's grave.
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not GATES, reason="set PYFP3D_TRANSONIC_GATES=1")
def test_schur_medium_lifting_escape():
    """2.5D medium lifting: ILU diverges here (gamma = -136.99, B11/B12) --
    the roadmap's pre-registered discriminating tier. schur must converge to
    the spsolve gamma (0.14137632). Passing here is what 'a real escape'
    means."""
    mesh = _mesh("medium")
    ref = solve_multivalued_lifting(_mvop(mesh), mesh, 0.5, alpha_deg=ALPHA)
    s = solve_multivalued_lifting(_mvop(mesh), mesh, 0.5, alpha_deg=ALPHA,
                                  precond="schur")
    assert ref["converged"] and s["converged"]
    assert abs(s["gamma"] - ref["gamma"]) < 1e-8
    assert s["n_gmres_stalled"] == 0 and s["n_schur_fallback"] == 0
