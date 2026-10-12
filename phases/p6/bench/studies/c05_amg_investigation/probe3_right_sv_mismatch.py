"""Q1 probe part 3: correct soft-mode direction + preconditioner-model mismatch.

Probe 2's Cop = J^{-T} J^{-1} = (J J^T)^{-1} returns the LEFT singular vector
(residual-space direction). The direction GMRES must invert is the RIGHT
singular vector, the top eigenvector of (J^T J)^{-1}, matvec x ->
J^{-1}(J^{-T} x) via the same splu. Both are measured here, plus the
preconditioner-model mismatch ||J_ff - A_ff||_F/||J_ff||_F and where the
mismatch lives (row-localized).

Appends metrics to freeze_point_metrics.csv. Does not touch the GMRES CSV/PNG
(probe 2's long-run evidence stands).
"""
import csv
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, "/home/lrz/codes/UP3D")
os.chdir("/home/lrz/codes/UP3D")

from tests import _nozzle_case as NZ                      # noqa: E402
from pyfp3d.solve.newton import solve_newton_nonlifting   # noqa: E402
from pyfp3d.solve.linear import (build_amg_preconditioner,  # noqa: E402
                                 build_ilu_preconditioner)

OUT = "/home/lrz/codes/UP3D/bench/studies/c05_amg_investigation"
M_INF, UPWIND_C, NX = 0.80, 1.5, 200

ex = NZ.exact_solution(M_INF)
mesh = NZ.nozzle_mesh(NX, max(6, NX // 16))
dir_nodes = np.unique(np.concatenate([
    mesh.boundary_faces["inlet"].ravel(),
    mesh.boundary_faces["outlet"].ravel()]))
dir_vals = ex["phi_of_x"](mesh.nodes[dir_nodes, 0])
phi0 = ex["phi_of_x"](mesh.nodes[:, 0])

r = solve_newton_nonlifting(
    mesh, m_inf=M_INF, upwind_c=UPWIND_C, dirichlet_nodes=dir_nodes,
    dirichlet_values=dir_vals, phi_init=phi0, entropy_correction=True,
    precond="amg", n_newton_max=40)
ws = r["workspace"]
phi = np.asarray(r["phi"], dtype=np.float64)
R_free, _, state = ws.eval_residual(phi[ws.free], np.zeros(0), UPWIND_C,
                                    0.95, 3.0, 0.05)
J_ff, _ = ws.assemble_coupled(state, UPWIND_C, 0.95, 0.05)
J_ff = J_ff.tocsr()
rhs = -R_free
n = J_ff.shape[0]
A_ff = ws.op.assemble_matrix(state["rho_t"])[ws.free][:, ws.free].tocsr()
lu = spla.splu(J_ff.tocsc())
_, M_amg = build_amg_preconditioner(A_ff)
M_ilu = build_ilu_preconditioner(J_ff)

metrics = {}
rng = np.random.default_rng(0)


def soft_directions():
    """(u_left, sigma_min_l), (v_right, sigma_min_r) via LU-preconditioned
    eigsh on the two inverse Gram operators (both SPD)."""
    JTJ_inv = spla.LinearOperator(
        (n, n), matvec=lambda x: lu.solve(lu.solve(x, trans="T")))
    lam_r, V_r = spla.eigsh(JTJ_inv, k=1, which="LA", maxiter=2000,
                            v0=rng.standard_normal(n))
    JJT_inv = spla.LinearOperator(
        (n, n), matvec=lambda x: lu.solve(lu.solve(x), trans="T"))
    lam_l, V_l = spla.eigsh(JJT_inv, k=1, which="LA", maxiter=2000,
                            v0=rng.standard_normal(n))
    return (V_l[:, 0], float(1.0 / np.sqrt(lam_l[0])),
            V_r[:, 0], float(1.0 / np.sqrt(lam_r[0])))


u_l, smin_l, v_r, smin_r = soft_directions()
metrics["sigma_min_left_eigsh"] = smin_l
metrics["sigma_min_right_eigsh"] = smin_r
print(f"sigma_min: left-op {smin_l:.4e}  right-op {smin_r:.4e}", flush=True)

x_nodes = mesh.nodes[:, 0]
for tag, vec in (("left_u", u_l), ("right_v", v_r)):
    va = np.abs(vec)
    i_peak = int(np.argmax(va))
    supp = va >= 0.01 * va.max()
    metrics[f"soft_{tag}_x_at_peak"] = float(x_nodes[ws.free[i_peak]])
    metrics[f"soft_{tag}_support_frac"] = float(np.count_nonzero(supp) / n)
    metrics[f"soft_{tag}_x_min"] = float(x_nodes[ws.free[supp]].min())
    metrics[f"soft_{tag}_x_max"] = float(x_nodes[ws.free[supp]].max())
    print(f"  {tag}: peak@x={metrics[f'soft_{tag}_x_at_peak']:.3f} "
          f"support [{metrics[f'soft_{tag}_x_min']:.3f},"
          f"{metrics[f'soft_{tag}_x_max']:.3f}] "
          f"({100*metrics[f'soft_{tag}_support_frac']:.2f}% dofs)", flush=True)

# the inversion test GMRES actually faces: M ~= J^{-1} applied to J v_right
jv = J_ff @ v_r
print(f"  |J v_right| = {np.linalg.norm(jv):.4e} (= sigma_min)", flush=True)
for tag, M in (("amg", M_amg), ("ilu", M_ilu)):
    err = float(np.linalg.norm(M @ jv - v_r) / np.linalg.norm(v_r))
    metrics[f"soft_right_inversion_err_{tag}"] = err
    print(f"  inversion err on right SV [{tag}]: {err:.4e}", flush=True)
metrics["rhs_loading_on_right_sv"] = float(abs(v_r @ rhs) / np.linalg.norm(rhs))
print(f"  rhs loading on right SV: {metrics['rhs_loading_on_right_sv']:.4e}",
      flush=True)

# preconditioner-model mismatch: how far is the matrix AMG was given from the
# matrix GMRES runs on, and where
D = (J_ff - A_ff).tocsr()
fro_J = spla.norm(J_ff)
metrics["model_mismatch_fro_ratio"] = float(spla.norm(D) / fro_J)
dA = np.abs(A_ff.diagonal())
dD = np.abs(D.diagonal())
metrics["diag_shift_max_rel"] = float((dD / np.maximum(dA, 1e-30)).max())
rownorm = np.asarray(sp.linalg.norm(D, axis=1)).ravel() \
    if hasattr(sp.linalg, "norm") else np.sqrt(np.asarray(D.multiply(D)
                                               .sum(axis=1)).ravel())
top = np.argsort(rownorm)[-200:]
metrics["mismatch_top200_x_min"] = float(x_nodes[ws.free[top]].min())
metrics["mismatch_top200_x_max"] = float(x_nodes[ws.free[top]].max())
# row-mismatch profile by x-bin for the README
nb = 40
edges = np.linspace(0.0, NZ.LENGTH, nb + 1)
xc = x_nodes[ws.free]
prof = np.zeros(nb)
for b in range(nb):
    m = (xc >= edges[b]) & (xc < edges[b + 1])
    if m.any():
        prof[b] = rownorm[m].mean()
with open(os.path.join(OUT, "model_mismatch_profile.csv"), "w",
          newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["x_bin_left", "x_bin_right", "mean_row_norm_J_minus_A"])
    for b in range(nb):
        w.writerow([f"{edges[b]:.4f}", f"{edges[b+1]:.4f}", f"{prof[b]:.16e}"])
print(f"  ||J-A||/||J|| = {metrics['model_mismatch_fro_ratio']:.4f}, "
      f"top-200 mismatch rows x in [{metrics['mismatch_top200_x_min']:.2f},"
      f"{metrics['mismatch_top200_x_max']:.2f}]", flush=True)

old = {}
with open(os.path.join(OUT, "freeze_point_metrics.csv")) as fh:
    for row in csv.reader(fh):
        if row and row[0] != "metric":
            old[row[0]] = row[1]
old.update({k: v for k, v in metrics.items()})
with open(os.path.join(OUT, "freeze_point_metrics.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["metric", "value"])
    for k, v in old.items():
        w.writerow([k, v])
print("DONE", flush=True)
