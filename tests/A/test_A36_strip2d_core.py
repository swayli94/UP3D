"""GS4.1 round 1: locks for the 2-D strip core (`pyfp3d/viscous/strip2d.py`).

These run on the ungated cadence and need no solve. They exist because phase 3
already paid for the opposite arrangement once: `meshgen/structured.py` shipped
with its only assertions living in a bench script, i.e. on no cadence at all.

The headline readings of the round are anchored here on purpose. The strip
core's verdict is that the closure family's flat-plate fixed point sits
+4.5 % (H) and +6.9 % (c_f) from Blasius; if a later closure change moves
those, it must be loud rather than silently re-baselining a recorded FAIL.
"""

import numpy as np
import pytest

from pyfp3d.viscous import closures as C
from pyfp3d.viscous import strip2d as S

RHO, MU, U = 1.0, 1.0e-5, 1.0

# Blasius, from the classical similarity solution (bench study integrates the
# ODE; these are the values it reproduces to six digits).
H_BLASIUS = 2.591100
CF_BLASIUS = 0.664115


class TestClosureContract:
    """Properties of `closures.py` that the marching FORM depends on."""

    def test_laminar_thicknesses_are_re_d_independent(self):
        """The implicit ODE `M y' = F` is exact under a pressure gradient only
        because laminar theta/theta* do not depend on the edge speed."""
        st = (2.0e-2, 7.883, 0.0, 0.0, C.CTAU_LAM, 0.0)
        o0, _, de0 = C.closure_scalar(st, q=1.0, rho=RHO, mu=MU)
        o1, _, _ = C.closure_scalar(st, q=1.7, rho=RHO, mu=MU)
        assert o1[C.OUT_TH11] == pytest.approx(o0[C.OUT_TH11], rel=1e-14)
        assert o1[C.OUT_THS1] == pytest.approx(o0[C.OUT_THS1], rel=1e-14)
        assert de0[C.OUT_TH11, 0] == 0.0
        assert de0[C.OUT_THS1, 0] == 0.0

    def test_turbulent_thicknesses_do_depend_on_re_d(self):
        """The counterpart: the explicit d/d(re_d) term in `_rhs` is not dead
        code, so a turbulent strip under a pressure gradient needs it."""
        st = np.array([2.0e-2, 60.0, 0.0, 0.0, 1.0e-3, 0.0])
        _, _, de = C.closure_scalar(st, q=1.0, rho=RHO, mu=MU, turbulent=True)
        assert abs(de[C.OUT_TH11, 0]) > 0.0


class TestSimilarityFixedPoint:
    """The algebraic (march-free, discretization-free) similar state."""

    def test_flat_plate_fixed_point_anchored(self):
        A, H, cf = S.similarity_fixed_point(m=0.0, rho=RHO, mu=MU)
        assert A == pytest.approx(8.028406562, rel=1e-6)
        assert H == pytest.approx(2.707931, rel=1e-6)
        assert cf == pytest.approx(0.710301160, rel=1e-6)

    def test_matches_the_independently_established_gv11_value(self):
        """GV1.1 established H* = 2.7083 by three independent constructions;
        this algebra is a fourth and must agree with it."""
        _, H, _ = S.similarity_fixed_point(m=0.0, rho=RHO, mu=MU)
        assert H == pytest.approx(2.707931, abs=5e-5)

    def test_recorded_gap_to_blasius(self):
        """The round's FAIL, as a number: H inside +-5 %, c_f outside."""
        _, H, cf = S.similarity_fixed_point(m=0.0, rho=RHO, mu=MU)
        assert abs(H / H_BLASIUS - 1.0) == pytest.approx(0.0450893, abs=1e-4)
        assert abs(cf / CF_BLASIUS - 1.0) == pytest.approx(0.06945, abs=1e-4)

    def test_invariant_under_the_dimensional_scales(self):
        """A self-similar state cannot depend on rho, mu or u_e."""
        base = S.similarity_fixed_point(m=0.0, rho=1.0, mu=1.0e-5, ue=1.0)
        for kw in ({"rho": 1.3}, {"mu": 4.0e-6}, {"ue": 55.0}):
            got = S.similarity_fixed_point(m=0.0, **{**{"rho": 1.0,
                                                        "mu": 1.0e-5,
                                                        "ue": 1.0}, **kw})
            assert got == pytest.approx(base, rel=1e-9)

    @pytest.mark.parametrize("m,expect_H", [(0.0, 2.707931),
                                            (1.0 / 23.0, 2.469240),
                                            (1.0, 2.087926)])
    def test_wedge_fixed_points(self, m, expect_H):
        _, H, _ = S.similarity_fixed_point(m=m, rho=RHO, mu=MU)
        assert H == pytest.approx(expect_H, rel=1e-6)


class TestSeeding:
    def test_seed_round_trips_theta_and_H(self):
        y = S.similar_seed(1.0e-3, 2.4, ue=U, rho=RHO, mu=MU)
        out, _, _ = C.closure_scalar((y[0], y[1], 0.0, 0.0, C.CTAU_LAM, 0.0),
                                     q=U, rho=RHO, mu=MU)
        assert out[C.OUT_TH11] == pytest.approx(1.0e-3, rel=1e-10)
        assert out[C.OUT_H1] == pytest.approx(2.4, rel=1e-9)

    def test_seed_stays_on_the_branch_connected_to_blasius(self):
        """H(A) folds, so a bare Newton can land on the mirror root."""
        for target in (2.2, 2.45, H_BLASIUS, 2.9):
            y = S.similar_seed(1.0e-3, target, ue=U, rho=RHO, mu=MU)
            assert S._H_BRANCH_LO <= y[1] <= S._H_BRANCH_HI

    def test_unrepresentable_H_raises_rather_than_extrapolating(self):
        with pytest.raises(ValueError, match="outside the physical A-branch"):
            S.similar_seed(1.0e-3, 9.9, ue=U, rho=RHO, mu=MU)


class TestMarch:
    def test_flat_plate_march_reaches_the_algebraic_fixed_point(self):
        """The march and the algebra are independent constructions; a march
        that lands elsewhere is a marching defect, not a closure property."""
        A_alg, H_alg, cf_alg = S.similarity_fixed_point(m=0.0, rho=RHO, mu=MU)
        y0 = S.similar_seed(0.664 * 0.01 / np.sqrt(1.0e3), H_BLASIUS,
                            ue=U, rho=RHO, mu=MU)
        st = S.march(np.array([1.0, 100.0]), y0, 0.01, S.flat_plate_ue(U),
                     rho=RHO, mu=MU, n_substep=2000)
        assert st.H[-1] == pytest.approx(H_alg, rel=2e-4)
        assert st.A[-1] == pytest.approx(A_alg, rel=2e-4)
        assert st.cf[-1] * np.sqrt(st.re_x[-1]) == pytest.approx(cf_alg,
                                                                 rel=2e-4)

    @pytest.mark.parametrize("m", [0.0, 1.0 / 23.0, 1.0])
    def test_pressure_gradient_terms_hold_the_similar_state(self, m):
        """Seeded exactly at the similar state, a correct march must not move.

        This is the guard on the pressure-gradient terms specifically: they are
        the part the GV1.1 reference arm never had.
        """
        A, H, _ = S.similarity_fixed_point(m=m, rho=RHO, mu=MU)
        ue_fn = S.falkner_skan_ue(m)
        x0 = 0.2
        ue0 = ue_fn(x0)[0]
        out, _, _ = C.closure_scalar((1.0, A, 0.0, 0.0, C.CTAU_LAM, 0.0),
                                     q=ue0, rho=RHO, mu=MU)
        t = out[C.OUT_TH11]
        d1 = 0.5 * (1.0 - m) + (H + 2.0) * m
        delta0 = np.sqrt(A * MU / (RHO * ue0) * x0 / (t * d1))
        y0 = S.similar_seed(t * delta0, H, ue=ue0, rho=RHO, mu=MU)
        st = S.march(np.array([1.0, 5.0, 25.0]), y0, x0, ue_fn,
                     rho=RHO, mu=MU, n_substep=400)
        assert np.max(np.abs(st.H - H)) / H < 1.0e-8

    def test_self_convergence_is_fourth_order(self):
        """RK4: the discretization error is what refinement moves. The error
        against Blasius is NOT -- that is the model floor (round-1 V-FAIL)."""
        y0 = S.similar_seed(0.664 * 0.01 / np.sqrt(1.0e3), H_BLASIUS,
                            ue=U, rho=RHO, mu=MU)
        stations = np.geomspace(0.1, 100.0, 7)
        runs = {n: S.march(stations, y0, 0.01, S.flat_plate_ue(U), rho=RHO,
                           mu=MU, n_substep=n) for n in (250, 500, 1000)}
        ref = S.march(stations, y0, 0.01, S.flat_plate_ue(U), rho=RHO, mu=MU,
                      n_substep=16000)
        errs = [np.abs(runs[n].H - ref.H).max() for n in (250, 500, 1000)]
        assert errs[0] > errs[1] > errs[2]
        order = -np.polyfit(np.log([250, 500, 1000]), np.log(errs), 1)[0]
        assert order > 3.5, f"RK4 order collapsed to {order:.2f}"

    def test_error_against_blasius_does_not_move_under_refinement(self):
        """The V-ORDER FAIL, locked: refinement cannot buy a model floor."""
        y0 = S.similar_seed(0.664 * 0.01 / np.sqrt(1.0e3), H_BLASIUS,
                            ue=U, rho=RHO, mu=MU)
        stations = np.geomspace(0.1, 100.0, 7)
        e = [np.abs(S.march(stations, y0, 0.01, S.flat_plate_ue(U), rho=RHO,
                            mu=MU, n_substep=n).H - H_BLASIUS).max()
             for n in (250, 1000)]
        assert e[0] == pytest.approx(e[1], rel=1e-8)
        assert e[0] / H_BLASIUS > 0.04

    def test_march_rejects_ill_posed_inputs(self):
        y0 = np.array([1.0e-3, 8.0])
        ue = S.flat_plate_ue(U)
        with pytest.raises(ValueError, match="strictly increasing"):
            S.march(np.array([2.0, 1.0]), y0, 0.1, ue)
        with pytest.raises(ValueError, match="downstream of x_start"):
            S.march(np.array([0.05, 1.0]), y0, 0.1, ue)
        with pytest.raises(ValueError, match="x_start must be positive"):
            S.march(np.array([1.0, 2.0]), y0, 0.0, ue)


class TestFrozenAndReuse:
    """G-REUSE / G-FROZEN as suite-cadence locks, not only bench-script ones."""

    def test_no_closure_constant_literals_in_the_strip_core(self):
        src = open(S.__file__).read()
        for name in ("KAPPA", "B_SPALDING", "A1_BRADSHAW", "C_L_DEFAULT",
                     "RECOVERY_R", "GAMMA_AIR"):
            assert repr(getattr(C, name)) not in src, (
                f"{name}'s value is inlined in strip2d.py -- the strip core "
                "must CALL the closure, not restate it")

    def test_strip_core_does_not_touch_the_3d_solver(self):
        """Checks the IMPORTS, not the prose: the module docstring names
        ibl3.py in order to state independence from it, and a bare substring
        search cannot tell that apart from a dependency.
        """
        import ast
        tree = ast.parse(open(S.__file__).read())
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(a.name for a in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.add(node.module or "")
                imported.update(f"{node.module}.{a.name}" for a in node.names)
        assert not any("ibl3" in m for m in imported), imported
        assert not any(m.startswith("pyfp3d.solve") for m in imported), imported
