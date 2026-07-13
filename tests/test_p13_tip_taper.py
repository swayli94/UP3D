"""
P13/G13.2 -- tip-edge desingularization by a spanwise loading taper.

The tip-edge singularity (G13.1: peak Mach ~ h^-p with p = 0.59, a 1/sqrt(r)
flat-plate-edge singularity) is driven by the TRAILING vorticity, and its
DISCRETE strength is set by the circulation retained at the OUTERMOST TE
station: that station sheds Gamma_last as a concentrated vortex over the last
cell (the free-edge nodes are single-valued, so the jump falls to 0 in one
element), inducing a velocity ~ Gamma_last / h. With Gamma_last ~ h^q the edge
peak grows as h^(q-1), i.e.

    p  ~  1 - q        ==>       the criterion is  q >= 1

(measured: the untapered baseline has q = 0.44 -> p_pred 0.56 vs p_meas 0.52).
`tip_taper_factors` supplies F(z) with Gamma_eff = F * Gamma_Kutta; writing
F ~ u^s near the tip (u = z_tip - z) gives Gamma_eff ~ u^(1/2 + s) for the
near-elliptic Gamma ~ sqrt(u), which is what drives q.

These tests lock:
  * the taper forms and their near-tip exponents s (the theory ladder),
  * that the DEFAULT (tip_taper=None) path is BIT-IDENTICAL,
  * that the Newton Kutta row K is scaled by the SAME F_j as the residual
    (taper is phi-independent, so dF/dphi picks up exactly that factor --
    getting this wrong silently breaks quadratic convergence),
  * that the taper reaches the transonic driver through newton_kw.

The quantitative regularization evidence (p -> 0, M6 fine becoming a genuine
discrete solution) lives in cases/demo/p13_tip_edge_singularity/.
"""
import numpy as np
import pytest
import scipy.sparse as sp

from pyfp3d.constraints.wake import kutta_targets, tip_taper_factors
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting

FORMS = ("none", "tanh_half", "vanish_sqrt", "vanish_linear", "vanish_smooth")
# near-tip exponent s of F ~ u^s  (F(tip)=1/2 for tanh_half => s = 0)
S_EXPONENT = {"none": 0.0, "tanh_half": 0.0, "vanish_sqrt": 0.5,
              "vanish_linear": 1.0, "vanish_smooth": 2.0}


# --------------------------------------------------------------------------
# 1. the taper function itself (no mesh needed -- always runs)
# --------------------------------------------------------------------------
def test_taper_range_and_far_field_is_untouched():
    z_tip = 2.0
    z = np.linspace(0.0, z_tip, 200)
    for form in FORMS:
        F = tip_taper_factors(z, z_tip, form, r_c=0.2)
        assert F.shape == z.shape
        assert np.all(F >= 0.0) and np.all(F <= 1.0 + 1e-12)
        # far inboard of the taper the wing must be (essentially) untouched
        assert F[0] == pytest.approx(1.0, abs=1e-6), form
        # monotone non-decreasing inboard (F rises from the tip to 1)
        assert np.all(np.diff(F) <= 1e-12), form


def test_tanh_tails_reach_inboard_but_the_vanish_forms_have_COMPACT_support():
    """★ The structural defect of the tanh form. tanh never reaches 1, so its
    tails unload the wing arbitrarily far inboard: on the M6 it depresses
    F below 0.99 over 57 of 83 TE stations (down to eta ~ 0.77), costing ~7x
    the lift of a compact taper for the same regularization and even breaking
    TE pressure closure at eta = 0.90 -- where there is no singularity to fix.
    The vanish_* forms have COMPACT support: F == 1 exactly outside r_c."""
    z_tip, r_c = 1.0, 0.1
    z = np.linspace(0.0, z_tip, 400)
    F_tanh = tip_taper_factors(z, z_tip, "tanh_half", r_c=r_c)
    # the tanh is strictly < 1 EVERYWHERE -- unbounded support
    assert np.all(F_tanh < 1.0)
    # and it is still meaningfully depressed several r_c inboard of the tip
    far = z < z_tip - 3.0 * r_c
    assert F_tanh[far].min() < 1.0

    for form in ("vanish_sqrt", "vanish_linear", "vanish_smooth"):
        F = tip_taper_factors(z, z_tip, form, r_c=r_c)
        outside = z <= z_tip - r_c
        assert np.array_equal(F[outside], np.ones(outside.sum())), (
            f"{form} must be EXACTLY 1 outside the taper (compact support)")


def test_none_is_exactly_unity():
    z = np.linspace(0.0, 1.0, 17)
    assert np.array_equal(tip_taper_factors(z, 1.0, "none"), np.ones(17))


def test_tip_values_separate_the_forms():
    """F(tip): 1/2 for the tanh (s=0 -- the marginal case), 0 for the
    vanishing forms (s > 0)."""
    z_tip = 1.0
    F = {f: float(tip_taper_factors(np.array([z_tip]), z_tip, f, r_c=0.1)[0])
         for f in FORMS}
    assert F["none"] == pytest.approx(1.0)
    assert F["tanh_half"] == pytest.approx(0.5)          # <-- s = 0
    for f in ("vanish_sqrt", "vanish_linear", "vanish_smooth"):
        assert F[f] == pytest.approx(0.0, abs=1e-12), f


@pytest.mark.parametrize("form", ["tanh_half", "vanish_sqrt",
                                  "vanish_linear", "vanish_smooth"])
def test_near_tip_exponent_of_the_tapered_circulation(form):
    """Gamma ~ sqrt(u) (near-elliptic) tapered by F ~ u^s must give
    Gamma_eff ~ u^(1/2 + s). This is the ladder the G13.2 probe walked:
    only s > 1/2 sends the sheet-edge strength gamma_eff = dGamma_eff/du to
    zero -- and the tanh (s = 0) does NOT, which is why it only regularizes
    marginally (q ~ 1) and at ~7x the lift cost."""
    z_tip, r_c = 1.0, 0.2
    u = np.array([1e-6, 1e-5])                 # deep inside the taper
    F = tip_taper_factors(z_tip - u, z_tip, form, r_c=r_c)
    g_eff = np.sqrt(u) * F
    alpha = np.log(g_eff[1] / g_eff[0]) / np.log(u[1] / u[0])
    assert alpha == pytest.approx(0.5 + S_EXPONENT[form], abs=0.05), (
        f"{form}: Gamma_eff ~ u^{alpha:.3f}, expected u^{0.5 + S_EXPONENT[form]}")


def test_bad_inputs_raise():
    z = np.linspace(0.0, 1.0, 5)
    with pytest.raises(ValueError):
        tip_taper_factors(z, 1.0, "not_a_form", r_c=0.1)
    with pytest.raises(ValueError):
        tip_taper_factors(z, 1.0, "vanish_linear", r_c=0.0)


# --------------------------------------------------------------------------
# 2. solver wiring (quasi-2D NACA: 1 station, but exercises every code path)
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def naca_coarse():
    from .conftest import REPO_ROOT
    mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
                     / "coarse.msh")
    return cut_wake(mesh)


def _args():
    return dict(m_inf=0.5, alpha_deg=2.0, upwind_c=1.5, precond="amg",
                tol_residual=1e-10)


def test_default_path_is_bit_identical(naca_coarse):
    """tip_taper=None must reproduce the untapered solver BIT-for-BIT (the
    guard that P13 cannot regress P2-P10)."""
    mc, wc = naca_coarse
    a = solve_newton_lifting(mc, wc, **_args())
    b = solve_newton_lifting(mc, wc, tip_taper=np.ones(wc.n_stations), **_args())
    assert np.array_equal(np.asarray(a["phi"]), np.asarray(b["phi"]))
    assert np.array_equal(np.asarray(a["gamma"]), np.asarray(b["gamma"]))


def test_kutta_row_is_scaled_by_the_taper(naca_coarse):
    """K = dF/dphi must carry the SAME F_j as the residual, since
    F_j = taper_j * mean(probe jumps)_j - Gamma_j and the taper does not
    depend on phi. If K is left unscaled the Jacobian is inconsistent and
    quadratic convergence dies silently."""
    mc, wc = naca_coarse
    taper = np.full(wc.n_stations, 0.37)
    ws0 = NewtonWorkspace(mc, wc, alpha_deg=2.0)
    ws1 = NewtonWorkspace(mc, wc, alpha_deg=2.0, tip_taper=taper)
    K0 = (sp.diags(taper) @ ws0.K).toarray()
    assert np.allclose(ws1.K.toarray(), K0, rtol=0, atol=1e-15)


def test_residual_uses_the_tapered_target(naca_coarse):
    """F = taper * kutta_targets - Gamma (the definition of Gamma_eff)."""
    mc, wc = naca_coarse
    taper = np.full(wc.n_stations, 0.6)
    ws = NewtonWorkspace(mc, wc, alpha_deg=2.0, tip_taper=taper)
    ws.set_mach(0.5)
    phi_free = np.zeros(ws.n_free)
    gamma = np.full(wc.n_stations, 0.05)
    _, F, state = ws.eval_residual(phi_free, gamma, 1.5, 0.95, 3.0, 0.05)
    expect = taper * kutta_targets(state["phi_cut"], wc) - gamma
    assert np.allclose(F, expect, rtol=0, atol=1e-15)


def test_taper_shape_is_validated(naca_coarse):
    mc, wc = naca_coarse
    with pytest.raises(ValueError):
        NewtonWorkspace(mc, wc, tip_taper=np.ones(wc.n_stations + 3))


def test_taper_reaches_the_transonic_driver(naca_coarse):
    """solve_newton_transonic has no tip_taper argument of its own: it must
    arrive through newton_kw and survive into the workspace that is REUSED
    across Mach levels. Tested behaviourally (the driver pops `workspace`
    from its result): a taper of 0.8 must visibly suppress the accepted
    circulation."""
    from pyfp3d.solve.newton import solve_newton_transonic
    mc, wc = naca_coarse
    kw = dict(upwind_c=1.5, precond="amg", tol_residual=1e-9)
    args = dict(m_inf=0.5, alpha_deg=2.0, m_start=0.5, dm=0.05)

    t = 0.8
    base = solve_newton_transonic(mc, wc, newton_kw=dict(kw), **args)
    tap = solve_newton_transonic(
        mc, wc, newton_kw=dict(tip_taper=np.full(wc.n_stations, t), **kw),
        **args)
    assert base["converged"] and tap["converged"]
    gb = float(np.asarray(base["gamma"])[0])
    gt = float(np.asarray(tap["gamma"])[0])

    # ★ The taper is AMPLIFIED, not applied. Gamma is a FIXED POINT of
    #   Gamma = t * Gamma_Kutta(Gamma), and the Kutta map has slope b ~ 0.93
    #   (P2 measurement). Linearising about the untapered fixed point Gamma*:
    #
    #       Gamma / Gamma*  =  t (1 - b) / (1 - t b)
    #
    #   which for t = 0.8, b = 0.93 gives 0.22 -- NOT 0.8. Because b is close
    #   to 1 the gain 1/(1 - t b) is large, so even a mild taper suppresses
    #   the converged circulation hard. This is why the measured Gamma_last
    #   decay exponents q (~3.3) far exceed the naive s + 1/2, and it is the
    #   reason the taper must be kept COMPACT and local to the tip.
    b = 0.93
    expect = t * (1.0 - b) / (1.0 - t * b)
    assert gt / gb == pytest.approx(expect, rel=0.25), (
        f"taper did not reach the driver, or the Kutta fixed-point gain "
        f"changed: measured {gt / gb:.3f}, expected ~{expect:.3f}")


# --------------------------------------------------------------------------
# 3. 3D: the taper must actually unload the TIP and leave the ROOT alone
#    (ONERA M6 .msh are gitignored -- regenerate with generate_onera_m6.py)
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def m6_coarse():
    from .conftest import REPO_ROOT
    p = REPO_ROOT / "cases" / "meshes" / "onera_m6" / "coarse.msh"
    if not p.exists():
        pytest.skip("onera_m6/coarse.msh not generated (gitignored)")
    return cut_wake(read_mesh(p))


def test_taper_is_local_to_the_tip_on_m6(m6_coarse):
    """The shipped model (vanish_smooth, r_c = 0.05 b) must be a TIP model:
    inboard circulation essentially untouched, tip stations unloaded. A form
    whose tails reach inboard (the tanh does, to eta ~ 0.77) silently re-rigs
    the whole wing -- measured in the G13.2 physics demo."""
    from pyfp3d.meshgen.wing3d import B_SEMI
    mc, wc = m6_coarse
    args = dict(m_inf=0.5, alpha_deg=3.06, upwind_c=1.5, precond="amg",
                tol_residual=1e-9, farfield_spanwise_gamma=True)
    taper = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth",
                              0.05 * B_SEMI)
    # the taper itself is local
    eta = np.asarray(wc.station_z) / B_SEMI
    assert np.all(taper[eta < 0.90] > 0.999), "taper must not reach eta < 0.90"
    assert 3 <= np.count_nonzero(taper < 0.99) <= 20, "taper must be resolved"

    base = solve_newton_lifting(mc, wc, **args)
    tap = solve_newton_lifting(mc, wc, tip_taper=taper, **args)
    assert base["converged"] and tap["converged"]
    gb, gt = np.asarray(base["gamma"]), np.asarray(tap["gamma"])

    # root/mid-span loading is preserved ...
    inb = eta < 0.90
    rel = np.abs(gt[inb] - gb[inb]) / np.max(np.abs(gb))
    assert rel.max() < 0.05, f"inboard loading disturbed by {rel.max():.3f}"
    # ... while the OUTERMOST station is strongly unloaded (that IS the model)
    j = int(np.argmax(wc.station_z))
    assert gt[j] < 0.5 * gb[j], "the tip station must be unloaded"
