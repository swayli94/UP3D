"""
Track B / B17 gates: the far-field aux pin carries jump = GAMMA, not 0.

B16 pinned the far-field-boundary aux (the outer nodes a wake level set crosses
on its way out, having no outflow clip) to cure the near-singular outer wake-LS
rows that make the LS Newton churn on the wing-body. On farfield="freestream"
B16 pinned them to jump=0. ★ B17 MEASURED (cases/demo/b17_farfield_pin_gamma)
that jump=0 REMOVES the wake circulation the outflow physically carries: the
medium ONERA-M6 wing-body cl_p drops 0.2165 -> 0.1690 (a 22% resolution-dependent
error, the GB16.4 "non-convergence" that was actually a BC-modelling error --
both Picard-pin and Newton-pin converge to the SAME wrong 0.169). The coarse
"match to conforming" was a coincidence (jump=0 there cancelled the coarse
legacy's outer-tet garbage). The fix is farfield_aux="pin_gamma": aux = host
phi_inf - side*gamma (jump=gamma), the same Dirichlet conditioning cure with the
physical ring value, and cl_p stays monotone-convergent to conforming.

This file locks the WIRING (cheap, ungated, on the committed coarse 2.5D NACA
mesh whose wake reaches the far field); the heavy wing-body M0.5 lift triangle
lives in the demo under PYFP3D_TRANSONIC_GATES.

  * The Picard knob (solve_multivalued_lifting farfield_aux) mirrors the Newton
    one and defaults to "legacy" so every committed Picard run is bit-identical.
  * pin_gamma imposes jump=gamma on the ring (Newton, refreshed per step);
    pin imposes jump=0 (the B16 reproduction).
  * Both are freestream-only; vortex/neumann pin raises / is inert.
"""

import inspect
import os
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.constraints.dirichlet import freestream_phi
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.solve.newton_ls import solve_multivalued_newton
from pyfp3d.solve.picard_ls import farfield_aux_dofs, solve_multivalued_lifting
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

GATES = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"
REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
ALPHA = 2.0


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


# ---------------------------------------------------------------------------
# 1. the Picard knob: default "pin_gamma" (the B17 fix), validated
# ---------------------------------------------------------------------------
def test_picard_farfield_aux_knob_default_pin_gamma():
    p = inspect.signature(solve_multivalued_lifting).parameters
    assert p["farfield_aux"].default == "pin_gamma"
    mesh = _naca()
    with pytest.raises(ValueError, match="farfield_aux"):
        solve_multivalued_lifting(mesh=mesh, mvop=_naca_mvop(mesh), m_inf=0.5,
                                  farfield_aux="bogus")


# ---------------------------------------------------------------------------
# 2. the default (pin_gamma) is INERT on the non-freestream far fields -- so the
#    changed default leaves every committed 2.5D vortex/neumann run untouched
# ---------------------------------------------------------------------------
def test_picard_pin_inert_on_non_freestream():
    """pin/pin_gamma act only on farfield='freestream'; on vortex (the Picard
    default farfield, used by every committed 2.5D NACA case) and neumann they
    are inert, byte-identical to legacy. This is what makes the default flip to
    'pin_gamma' safe."""
    mesh = _naca()
    for ff in ("vortex", "neumann"):
        common = dict(m_inf=0.5, alpha_deg=ALPHA, farfield=ff,
                      n_outer_max=4, tol_residual=1e-7)
        leg = solve_multivalued_lifting(mesh=mesh, mvop=_naca_mvop(mesh),
                                        farfield_aux="legacy", **common)
        for mode in ("pin", "pin_gamma"):
            got = solve_multivalued_lifting(mesh=mesh, mvop=_naca_mvop(mesh),
                                            farfield_aux=mode, **common)
            assert np.array_equal(leg["phi_ext"], got["phi_ext"]), (ff, mode)


def test_newton_pin_gamma_is_default():
    """B17 flipped the Newton default 'pin' (jump=0) -> 'pin_gamma' (jump=gamma).
    'pin' is kept as an explicit diagnostic value."""
    p = inspect.signature(solve_multivalued_newton).parameters
    assert p["farfield_aux"].default == "pin_gamma"
    mesh = _naca()
    with pytest.raises(ValueError, match="farfield_aux"):
        solve_multivalued_newton(mesh=mesh, mvop=_naca_mvop(mesh), m_inf=0.5,
                                 farfield_aux="bogus")


# ---------------------------------------------------------------------------
# 3. the RING VALUE: pin_gamma imposes jump=gamma, pin imposes jump=0
# ---------------------------------------------------------------------------
def test_newton_freestream_pin_gamma_ring_jump():
    """The freestream pin_gamma ring carries jump=gamma (side-signed, constant on
    the ring, nonzero); the plain 'pin' ring carries jump=0. The ring is written
    with the gamma at the top of each step (frozen within the step), so on a
    short un-converged run it equals that step's gamma_mean, not the final one."""
    mesh = _naca()
    mvop = _naca_mvop(mesh)
    hosts, aux = farfield_aux_dofs(mesh, mvop.cm)
    assert aux.size > 0

    r_g = solve_multivalued_newton(mesh=mesh, mvop=_naca_mvop(mesh), m_inf=0.5,
                                   alpha_deg=ALPHA, farfield="freestream",
                                   farfield_aux="pin_gamma", n_seed=5,
                                   n_newton_max=4)
    jump = mvop.node_jump(r_g["phi_ext"], hosts)
    assert np.allclose(jump, jump[0], atol=1e-12)          # constant on the ring
    assert abs(jump[0]) > 1e-3                              # carries circulation
    # equals the gamma the LAST step wrote (step_records[-1] gamma at entry)
    g_written = r_g["step_records"][-1]["gamma_mean"]
    assert np.allclose(jump, g_written, atol=1e-9)

    r_0 = solve_multivalued_newton(mesh=mesh, mvop=_naca_mvop(mesh), m_inf=0.5,
                                   alpha_deg=ALPHA, farfield="freestream",
                                   farfield_aux="pin", n_seed=5, n_newton_max=4)
    assert np.allclose(mvop.node_jump(r_0["phi_ext"], hosts), 0.0, atol=1e-10)


def test_picard_pin_gamma_ring_jump():
    """A settled Picard pin_gamma imposes jump=gamma on the far-field ring
    (constant, nonzero, ~ the extracted circulation); plain pin imposes jump=0.
    The ring is written at each outer's top with the previous outer's gamma, so a
    non-machine-converged run matches the extracted gamma only to the one-outer
    lag -- a loose tolerance keeps the WIRING assertion robust."""
    mesh = _naca()
    mvop = _naca_mvop(mesh)
    hosts, aux = farfield_aux_dofs(mesh, mvop.cm)

    r_g = solve_multivalued_lifting(mesh=mesh, mvop=_naca_mvop(mesh), m_inf=0.5,
                                    alpha_deg=ALPHA, farfield="freestream",
                                    farfield_aux="pin_gamma", n_outer_max=60,
                                    tol_residual=1e-7)
    gamma = float(np.mean(mvop.te_jump(r_g["phi_ext"])))
    jump = mvop.node_jump(r_g["phi_ext"], hosts)
    assert np.allclose(jump, jump[0], atol=1e-12)          # constant on the ring
    assert abs(jump[0]) > 1e-3                              # carries circulation
    assert np.allclose(jump, gamma, atol=2e-3)             # ring ~ extracted gamma

    r_0 = solve_multivalued_lifting(mesh=mesh, mvop=_naca_mvop(mesh), m_inf=0.5,
                                    alpha_deg=ALPHA, farfield="freestream",
                                    farfield_aux="pin", n_outer_max=60,
                                    tol_residual=1e-7)
    assert np.allclose(mvop.node_jump(r_0["phi_ext"], hosts), 0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# 4. the TE-aux disjointness guard (a pin must never overwrite a Kutta row)
# ---------------------------------------------------------------------------
def test_picard_pin_te_aux_guard_present():
    """The Picard pin path carries the same TE-aux disjointness RuntimeError as
    the Newton path. On the NACA 2.5D wake the TE aux and far-field aux are
    disjoint, so a normal pin run does NOT raise -- assert the guard code exists
    and the normal path is clean."""
    import pyfp3d.solve.picard_ls as m
    src = inspect.getsource(m.solve_multivalued_lifting)
    assert "a TE aux DOF lies on the far-field boundary" in src
    mesh = _naca()
    # a real pin_gamma outer completes without the guard firing
    r = solve_multivalued_lifting(mesh=mesh, mvop=_naca_mvop(mesh), m_inf=0.5,
                                  alpha_deg=ALPHA, farfield="freestream",
                                  farfield_aux="pin_gamma", n_outer_max=2)
    assert np.all(np.isfinite(r["phi_ext"]))
