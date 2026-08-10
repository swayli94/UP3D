"""N1 probe (Kimi inspection 2026-07-19): does the SHIPPED freeze_side_state
capture a different (upstream, branch) selection than a B20-patched capture
would, on M6 coarse at a seeded M0.70 transonic state?

Rationale: freeze_side_state (multivalued.py:739-745) computes q2l/rho on the
UNPATCHED side field, while newton_side_data (:667-668) applies
_apply_main_density first. If the captured selection differs, the
"frozen reproduces live bitwise at the freeze point" invariant is broken on
3-D meshes, and the B15 frozen Newton finish iterates on a wrong selection.
"""
from pathlib import Path

import numpy as np

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.meshgen.wing3d import x_te
from pyfp3d.physics.isentropic import density_field, limit_q2_field
from pyfp3d.solve.picard_ls import solve_multivalued_transonic
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO = Path("/home/lrz/code/UP3D_kimi")
ALPHA_3D, B_SEMI = 3.06, 1.1963
M_INF, C, MC, GA, UI, MCAP, RF = 0.70, 1.5, 0.95, 1.4, 1.0, 3.0, 0.05

mesh = read_mesh(REPO / "cases" / "meshes" / "onera_m6" / "coarse.msh")
a = np.radians(ALPHA_3D)
te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                   wall_nodes=np.unique(mesh.boundary_faces["wall"]))
mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)

seed = solve_multivalued_transonic(mvop, mesh, M_INF, alpha_deg=ALPHA_3D,
                                   farfield="freestream", n_outer_level=200)
phi = seed["phi_ext"]

# (1) shipped capture (computes q2l/rho on the UNPATCHED side field)
frozen_shipped = mvop.freeze_side_state(phi, M_INF, C, MC, GA, UI, MCAP, RF)

# (2) patched capture: identical, but with _apply_main_density applied first
upw_u, upw_l = mvop._side_upwind()
phi_up, phi_lo = mvop.side_potentials(phi)
patched = []
for phi_s, upw in ((phi_up, upw_u), (phi_lo, upw_l)):
    grad, q2 = mvop.op.velocities(phi_s)
    grad, q2 = mvop._apply_main_density(phi, grad, q2)
    q2l = limit_q2_field(q2 / UI**2, M_INF, MCAP, GA)
    rho = density_field(q2l, M_INF, GA)
    patched.append(upw.freeze_upwind_state(grad, q2l, rho, M_INF, C, MC, GA, RF))

mp = mvop._mixed_plain_mask()
el = np.asarray(mesh.elements, dtype=np.int64)
aux_touch = mp & (cm.ext_dof_of_node[el] >= 0).any(axis=1)
print(f"mixed-side plain elements: {int(mp.sum())}; aux-touching: {int(aux_touch.sum())}")
for name, (s, p) in zip(("upper", "lower"), zip(frozen_shipped, patched)):
    us, bs = s
    up_, bp = p
    d_b = bs != bp
    d_u = us != up_
    print(f"{name}: branch differs {int(d_b.sum())} "
          f"(mp {int((d_b & mp).sum())}, aux-touching mp {int((d_b & aux_touch).sum())}); "
          f"upstream differs {int(d_u.sum())} "
          f"(mp {int((d_u & mp).sum())}, aux-touching mp {int((d_u & aux_touch).sum())})")
