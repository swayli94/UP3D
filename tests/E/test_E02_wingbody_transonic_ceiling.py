"""The conforming WING-BODY transonic capability, locked (gated).

WHY THIS FILE EXISTS. Phase 3's first task deletes the level-set route (ruling D5),
and `cases/demo/b18_wingbody_transonic/` -- whose headline check IS the two-path
cross-model comparison -- is archived with it. Checking before archiving showed the
capability that demo also recorded, "conforming production reaches M0.84 at medium",
had NO live test lock anywhere:

  - tests/A/test_A28_wingbody_topology.py stops at M0.5 (junction structure + a
    Laplace lifting check),
  - tests/A/test_A26_gamma_pin_row_blend.py locks the Gamma-pin blend ALGEBRA and its FD
    exactness on a NACA 2.5-D case plus an M6 taper profile -- not the wing-body,
  - tests/test_b18_wingbody_transonic.py was 4/4 level-set and is archived.

So the capability boundary asserted M0.84 while nothing would go red if it broke --
and phase 3's next move (hexahedral body-fitted meshing) is exactly the kind of change
that could break it. This file closes that hole BEFORE the archive step, which is why
this round is not purely subtractive.

WHAT IT LOCKS, AND WHAT IT DOES NOT. There is no external truth for wing-body lift in
this project (`cases/reference_data/` has M6 wing Cp and a shock band, nothing for the
wing-body), so this is a DRIFT LOCK, not a correctness claim. It says: the production
conforming chain still climbs to M0.84 on the medium wing-body, still lands with ZERO
clamps, and still reports the lift it reported when the capability was recorded.

★ Every anchor is the value b18's committed evidence recorded, not a fresh reading:
    phases/p1/cases/demo/b18_wingbody_transonic/results/checks.csv  (GB18.1 row: "B32 taper
        climb: m_reached=0.84 (cl_p 0.2738; clamps 0+0)")
    phases/p1/cases/demo/b18_wingbody_transonic/results/cl_vs_mach.csv  (0.84 conforming medium
        -> cl_p 0.2738)
and the RECIPE is that demo's `conf_ramp` verbatim -- seed at M0.70 with the probe
Kutta, then the Mach continuation with pressure Kutta + the production tip taper. The
taper's -1.3 % cl_p bias rides on this number, as it does on every wing-body lift here.

COST: one full continuation ramp on the medium wing-body, ~35 min. Gated, and the mesh
is gitignored -- regenerate with
cases/meshes/onera_m6_wingbody_conforming/generate_onera_m6_wingbody_conforming.py.
"""

import os
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.constraints.wake import tip_taper_factors
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.wing3d import B_SEMI
from pyfp3d.post.surface import planform_area
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting, solve_newton_transonic
from tests.conftest import REPO_ROOT

REPO_ROOT = REPO_ROOT
MESH = REPO_ROOT / "cases" / "meshes" / "onera_m6_wingbody_conforming" / "medium.msh"
GATES = os.environ.get("PYFP3D_TRANSONIC_GATES") == "1"

ALPHA = 3.06
M_START, M_TARGET = 0.70, 0.84
#: b18's conf_ramp recipe, verbatim
SEED_KW = dict(farfield_spanwise_gamma=True, precond="direct",
               direct_refactor_every=1000, n_newton_max=60)
RAMP_NK = dict(freeze_refresh_max=8, precond="direct",
               direct_refactor_every=1000, n_newton_max=80,
               farfield_spanwise_gamma=True)
TAPER = ("vanish_smooth", 0.05)

#: committed anchors (b18 results/checks.csv + cl_vs_mach.csv, B32 refresh)
CL_P = 0.2738
CL_P_RTOL = 0.02          # 2 %: well outside run-to-run scatter, far inside the
#                           ~13 % transonic rise this chain measures (0.2143 -> 0.2738)


@pytest.mark.skipif(not GATES, reason="heavy gated wing-body ramp (~35 min)")
def test_conforming_wingbody_medium_reaches_m084():
    """The production conforming chain climbs the medium wing-body to M0.84 with
    zero clamps, at the lift b18 recorded."""
    if not MESH.exists():
        pytest.skip(f"{MESH.name} not generated (gitignored); run "
                    "cases/meshes/onera_m6_wingbody_conforming/"
                    "generate_onera_m6_wingbody_conforming.py")
    mc, wc = cut_wake(read_mesh(MESH))
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    seed = solve_newton_lifting(mc, wc, m_inf=M_START, alpha_deg=ALPHA, **SEED_KW)
    taper = tip_taper_factors(wc.station_z, B_SEMI, TAPER[0], TAPER[1] * B_SEMI)
    r = solve_newton_transonic(
        mc, wc, m_inf=M_TARGET, alpha_deg=ALPHA, m_start=M_START, dm=0.05,
        dm_min=0.01, freeze_tol=1e-5, intermediate_tol=1e-4,
        newton_kw=dict(RAMP_NK, kutta_estimator="pressure", tip_taper=taper,
                       phi_init=seed["phi"], gamma_init=seed["gamma"],
                       n_picard_seed=0))
    m_reached = r["level_history"][-1][0]
    assert r["converged"], f"ramp not converged at m={m_reached:.4g}"
    assert abs(m_reached - M_TARGET) < 1e-9, (
        f"reached only M{m_reached:.4g}, not M{M_TARGET} -- the wing-body "
        f"transonic ceiling MOVED DOWN (b18 recorded M0.84)")
    #: ★ the clamp clause is the physical one: GS1.4 says a clamped state is not a
    #: solution, so "reached M0.84" only means something with 0/0.
    assert r["n_limited"] == 0 and r["n_floored"] == 0, (
        f"{r['n_limited']} limited / {r['n_floored']} floored at M{m_reached:.4g}; "
        f"b18 recorded 0+0, and a clamped state is not a solution")
    #: ★ GS4.0 R2: was post.unified.wall_forces, which phase 3 had collapsed onto
    #: its conforming half and which its own docstring declared np.array_equal to
    #: this function. VERIFIED before the swap on the M6 wing M0.8395 state --
    #: cl 0.2775363765023681 both ways, cf/cp_tri array_equal -- then unified.py
    #: was deleted, because a dispatch layer over one path is a second name for
    #: the same computation.
    cl_p = float(wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
        alpha_deg=ALPHA, s_ref=s_ref, m_inf=m_reached)["cl"])
    assert abs(cl_p / CL_P - 1.0) < CL_P_RTOL, (
        f"cl_p {cl_p:.4f} vs the committed {CL_P:.4f} "
        f"({100 * (cl_p / CL_P - 1.0):+.2f} %, tol {100 * CL_P_RTOL:.0f} %)")


if __name__ == "__main__":
    pytest.main([__file__, "-x", "-q"])
