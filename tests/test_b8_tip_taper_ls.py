"""
Track B / B8 gate: level-set tip-edge desingularization (row-blend tip taper).

docs/roadmap.md Track B B8; docs/agent-rules.md P13/G13.2 finding (8).

P13/G13.2 fixed the tip/wake-edge singularity on the CONFORMING path with a
spanwise loading taper Gamma_eff(z) = F(z)*Gamma_Kutta(z). That cannot be ported
literally to the level-set path: the LS path has NO Gamma DOF and its TE Kutta
row s.(q_u - q_l) = 0 is HOMOGENEOUS (RHS = 0), so scaling it by F is an
algebraic no-op. The fix is a convex BLEND of the pressure-equality row with
B2's continuity weld, per TE node:

    F_i * [ s.(q_u - q_l) ]  +  (1 - F_i) * [ phi_aux - phi_main ]  = 0

F=1 inboard -> pure pressure Kutta (bit-identical). F=0 at the tip -> pure weld
-> jump = 0 at that node -> the tip is unloaded (the LS analogue of
Gamma_eff -> 0). The blend is NOT a no-op because the weld row is not
proportional to the (homogeneous) pressure row.

These fast tests pin the ALGEBRA + the two limits + the Newton consistency on
the cheap 2.5D NACA mesh; the physical mechanism probe (off-body tip-peak Mach
under refinement) and the two-path physics A/B live in the (gated, 3D M6) demo
`cases/demo/b8_tip_taper_ls/`.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.constraints.wake import tip_taper_factors
from pyfp3d.kernels.cut_assembly import te_kutta_residual, te_weld_coo
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.solve.picard_ls import solve_multivalued_lifting
from pyfp3d.solve.newton_ls import solve_multivalued_newton
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
M6_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6"
ALPHA = 2.0


def _load(directory, level):
    path = directory / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} missing")
    return read_mesh(path)


def _setup(mesh):
    """Chord-plane wake (design.md §4)."""
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(
        mesh.nodes, mesh.elements, wls,
        wall_nodes=np.unique(mesh.boundary_faces["wall"]),
    )
    return wls, cm, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                        levelset=wls)


def _phi_seed(mvop):
    """A deterministic extended state to linearize the TE Kutta row (the
    matrix-identity tests only need a fixed phi, not a converged one)."""
    rng = np.random.default_rng(0)
    return rng.standard_normal(mvop.n_total)


class TestWeldHelper:
    """The (1 - F) weld term (te_weld_coo)."""

    def test_all_ones_is_empty(self):
        _, cm, _ = _setup(_load(M0_DIR, "coarse"))
        r, c, d = te_weld_coo(cm, np.ones(len(cm.te_nodes)))
        assert len(r) == len(c) == len(d) == 0

    def test_zeros_is_full_weld(self):
        _, cm, _ = _setup(_load(M0_DIR, "coarse"))
        n = len(cm.te_nodes)
        r, c, d = te_weld_coo(cm, np.zeros(n))
        # two entries per TE node: +1 at aux, -1 at main
        assert len(r) == 2 * n
        aux = cm.ext_dof_of_node[cm.te_nodes]
        assert set(r.tolist()) == set(aux.tolist())
        # aux-diagonal +1, aux->main -1
        assert np.isclose(d[:n], 1.0).all()
        assert np.isclose(d[n:], -1.0).all()
        assert set(c[n:].tolist()) == set(np.asarray(cm.te_nodes).tolist())

    def test_partial_only_blended_nodes(self):
        _, cm, _ = _setup(_load(M0_DIR, "coarse"))
        n = len(cm.te_nodes)
        vals = [0.0, 0.5, 0.9][:n]      # blend the first min(3, n) nodes
        k = len(vals)
        F = np.ones(n)
        F[:k] = vals
        r, _, d = te_weld_coo(cm, F)
        assert len(r) == 2 * k
        # weld coefficient is (1 - F)
        assert np.isclose(sorted(d[:k]), sorted(1.0 - np.array(vals))).all()

    def test_wrong_shape_raises(self):
        _, cm, _ = _setup(_load(M0_DIR, "coarse"))
        with pytest.raises(ValueError):
            te_weld_coo(cm, np.ones(len(cm.te_nodes) + 1))


class TestBitIdentity:
    """tip_taper=None and F≡1 must both reproduce the untapered matrix
    BITWISE — the default path is untouched."""

    @pytest.mark.parametrize("level", ["coarse", "medium"])
    def test_none_equals_ones(self, level):
        _, cm, mvop = _setup(_load(M0_DIR, level))
        phi = _phi_seed(mvop)
        A0 = mvop.assemble_matrix(closure="wake_ls", te_kutta="pressure",
                                  phi_ext=phi, tip_taper=None)
        A1 = mvop.assemble_matrix(closure="wake_ls", te_kutta="pressure",
                                  phi_ext=phi,
                                  tip_taper=np.ones(len(cm.te_nodes)))
        d = (A0 - A1)
        d.eliminate_zeros()
        assert d.nnz == 0

    def test_solve_none_equals_ones(self):
        mesh = _load(M0_DIR, "coarse")
        _, cm, mvop = _setup(mesh)
        g_none = solve_multivalued_lifting(mvop, mesh, 0.0, alpha_deg=ALPHA)
        _, _, mvop2 = _setup(mesh)
        g_ones = solve_multivalued_lifting(
            mvop2, mesh, 0.0, alpha_deg=ALPHA,
            tip_taper=np.ones(len(cm.te_nodes)))
        assert g_none["gamma"] == pytest.approx(g_ones["gamma"], abs=1e-12)


class TestBlendIsNotANoOp:
    """★ finding (8): scaling the HOMOGENEOUS pressure row by F is a no-op;
    the blend is not, because the weld row is not proportional to it."""

    def test_scaling_noop_but_weld_bites(self):
        mesh = _load(M0_DIR, "coarse")
        _, cm, mvop = _setup(mesh)
        r = solve_multivalued_lifting(mvop, mesh, 0.0, alpha_deg=ALPHA)
        phi = r["phi_ext"]
        # at the converged untapered solution the pressure row is ~satisfied
        # (row.x ≈ 0), so scaling it by any F leaves x a solution -> no-op.
        kutta = np.max(np.abs(te_kutta_residual(mvop, phi)))
        # the weld row is FAR from satisfied there: there is a real jump.
        te_aux = cm.ext_dof_of_node[cm.te_nodes]
        weld = np.max(np.abs(phi[te_aux] - phi[np.asarray(cm.te_nodes)]))
        assert kutta < 1e-3
        assert weld > 0.05                      # ~ the circulation-scale jump
        assert weld > 100 * kutta               # the blend genuinely changes x

    def test_blend_changes_gamma(self):
        # The blend genuinely changes the solution (direction is NOT asserted:
        # a UNIFORM F on this tip-less 2.5D mesh is outside the intended
        # regime -- the taper is meant to be F≈1 in the bulk and ->0 only near
        # a real tip; uniform F is non-monotone, see the demo. Here we only
        # pin that it is not a no-op, unlike scaling the homogeneous row).
        mesh = _load(M0_DIR, "coarse")
        _, cm, mvop = _setup(mesh)
        g0 = solve_multivalued_lifting(mvop, mesh, 0.0, alpha_deg=ALPHA)["gamma"]
        _, _, mvop2 = _setup(mesh)
        gb = solve_multivalued_lifting(
            mvop2, mesh, 0.0, alpha_deg=ALPHA,
            tip_taper=np.full(len(cm.te_nodes), 0.5))["gamma"]
        assert abs(gb - g0) > 0.02 * abs(g0)


class TestFZeroReduction:
    """F≡0 welds every TE node -> jump pinned to 0 -> non-lifting. The LS
    analogue of a Laplace α=0 solve; also validates the unconditional mass
    reroute (the main row's total-fan balance IS the single-valued equation
    for a welded node)."""

    def test_all_welded_is_non_lifting(self):
        mesh = _load(M0_DIR, "coarse")
        _, cm, mvop = _setup(mesh)
        r = solve_multivalued_lifting(
            mvop, mesh, 0.0, alpha_deg=ALPHA,
            tip_taper=np.zeros(len(cm.te_nodes)))
        assert abs(r["gamma"]) < 1e-6
        assert abs(r["cl_kj"]) < 1e-3


class TestTaperFactorsOn3DTip:
    """F(z) is built from the per-TE-node spanwise arclength on the LS path,
    reusing the conforming tip_taper_factors. Structural only (no solve).
    Uses the real swept ONERA M6 TE (B7's level set), where there IS a tip."""

    def test_f_vanishes_at_the_m6_tip(self):
        from pyfp3d.meshgen.wing3d import B_SEMI, x_te
        mesh = _load(M6_DIR, "coarse")
        te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
        a = np.radians(ALPHA)
        wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
        cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                           wall_nodes=np.unique(mesh.boundary_faces["wall"]))
        assert len(cm.te_nodes) > 20, "M6 must have many spanwise TE nodes"
        q = cm.q[cm.te_nodes]
        F = tip_taper_factors(q, cm.span_length, "vanish_smooth",
                              r_c=0.05 * cm.span_length)
        assert (F >= 0.0).all() and (F <= 1.0).all()
        # inboard (small q) -> 1; at the tip (q -> span_length) -> ~0
        assert F[np.argmax(q)] < 0.2
        assert F[np.argmin(q)] == pytest.approx(1.0, abs=1e-9)


class TestNewtonBlend:
    """The blended residual + Jacobian on the LS Newton path (subsonic, coarse
    -> cheap). The point is CONSISTENCY: a taper still converges to a machine
    solution -- an inconsistent residual/Jacobian pair would break the
    terminal quadratic drop. F≡1 stays bit-identical."""

    def test_newton_none_bit_identical(self):
        mesh = _load(M0_DIR, "coarse")
        _, cm, mvop = _setup(mesh)
        a = solve_multivalued_newton(mvop, mesh, 0.0, alpha_deg=ALPHA,
                                     n_newton_max=12, verbose=False)
        _, _, mvop2 = _setup(mesh)
        b = solve_multivalued_newton(mvop2, mesh, 0.0, alpha_deg=ALPHA,
                                     n_newton_max=12,
                                     tip_taper=np.ones(len(cm.te_nodes)))
        assert a["gamma"] == pytest.approx(b["gamma"], abs=1e-12)

    def test_newton_taper_converges(self):
        # a tapered blend converges to a machine solution (consistency of the
        # blended residual + Jacobian) and changes Gamma. Direction/locality
        # of the tip unloading is a 3D-tip property -> the demo, not here.
        mesh = _load(M0_DIR, "coarse")
        _, cm, mvop = _setup(mesh)
        base = solve_multivalued_newton(mvop, mesh, 0.0, alpha_deg=ALPHA,
                                        n_newton_max=20)
        _, _, mvop2 = _setup(mesh)
        tap = solve_multivalued_newton(
            mvop2, mesh, 0.0, alpha_deg=ALPHA, n_newton_max=20,
            tip_taper=np.full(len(cm.te_nodes), 0.7))
        assert tap["converged"]
        assert tap["residual_history"][-1] < 1e-9
        assert abs(tap["gamma"] - base["gamma"]) > 1e-4


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
