"""Q1 probe part 2: long GMRES (production budget 600 inner iters) + cond2 +
soft mode via LU-preconditioned eigsh (ARPACK svds-SM failed to converge in
probe 1; C = J^{-T} J^{-1} is SPD and its top eigenpair gives sigma_min and
the right singular / soft-mode direction robustly).

Rewrites gmres_freeze_point.csv, extends freeze_point_metrics.csv, regenerates
c05_amg_gmres.png. newton_histories.csv from probe 1 stands (deterministic).
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
NY = max(6, NX // 16)

ex = NZ.exact_solution(M_INF)
mesh = NZ.nozzle_mesh(NX, NY)
dir_nodes = np.unique(np.concatenate([
    mesh.boundary_faces["inlet"].ravel(),
    mesh.boundary_faces["outlet"].ravel()]))
dir_vals = ex["phi_of_x"](mesh.nodes[dir_nodes, 0])
phi0 = ex["phi_of_x"](mesh.nodes[:, 0])

print("== reproduce amg leg to the freeze point ==", flush=True)
r = solve_newton_nonlifting(
    mesh, m_inf=M_INF, upwind_c=UPWIND_C, dirichlet_nodes=dir_nodes,
    dirichlet_values=dir_vals, phi_init=phi0, entropy_correction=True,
    precond="amg", n_newton_max=40)
ws = r["workspace"]
phi = np.asarray(r["phi"], dtype=np.float64)
print(f"  amg leg: n={r['n_newton']} |R|_end={r['residual_history'][-1]:.6e}",
      flush=True)

R_free, _, state = ws.eval_residual(phi[ws.free], np.zeros(0), UPWIND_C,
                                    0.95, 3.0, 0.05)
J_ff, _ = ws.assemble_coupled(state, UPWIND_C, 0.95, 0.05)
J_ff = J_ff.tocsr()
rhs = -R_free
n = J_ff.shape[0]

A_pic = ws.op.assemble_matrix(state["rho_t"])
A_ff = A_pic[ws.free][:, ws.free].tocsr()
lu = spla.splu(J_ff.tocsc())           # probe 1 proved this succeeds
x_ref = lu.solve(rhs)
ml, M_amg = build_amg_preconditioner(A_ff)
M_ilu = build_ilu_preconditioner(J_ff)

metrics = {}

# ---- conditioning: sigma_max via power iteration on J^T J (deterministic),
# ---- sigma_min via eigsh on C = J^{-T} J^{-1} (SPD, LU-solve matvecs)
print("== conditioning ==", flush=True)
rng = np.random.default_rng(0)
v = rng.standard_normal(n)
for _ in range(200):
    w = J_ff.T @ (J_ff @ v)
    v = w / np.linalg.norm(w)
smax = float(np.sqrt(v @ (J_ff.T @ (J_ff @ v))))
metrics["sigma_max_power"] = smax

Jinv = spla.LinearOperator((n, n), matvec=lu.solve,
                           rmatvec=lambda x: lu.solve(x, trans="T"))
Cop = spla.LinearOperator(
    (n, n), matvec=lambda x: Jinv.rmatvec(Jinv.matvec(x)))
lam, V = spla.eigsh(Cop, k=1, which="LA", maxiter=2000,
                    v0=rng.standard_normal(n))
smin = float(1.0 / np.sqrt(lam[0]))
v_soft = V[:, 0]
metrics["sigma_min_eigsh"] = smin
metrics["cond2_Jff"] = smax / smin
print(f"  sigma_max={smax:.4e} sigma_min={smin:.4e} cond2={smax/smin:.4e}",
      flush=True)

# soft mode support + how much the stalled RHS loads it
v_abs = np.abs(v_soft)
x_nodes = mesh.nodes[:, 0]
i_peak = int(np.argmax(v_abs))
thr = 0.01 * v_abs.max()
supp = v_abs >= thr
metrics["soft_mode_x_at_peak"] = float(x_nodes[ws.free[i_peak]])
metrics["soft_mode_support_frac"] = float(np.count_nonzero(supp) / n)
metrics["soft_mode_x_min"] = float(x_nodes[ws.free[supp]].min())
metrics["soft_mode_x_max"] = float(x_nodes[ws.free[supp]].max())
metrics["rhs_soft_mode_loading"] = float(abs(v_soft @ rhs) / np.linalg.norm(rhs))
jv = J_ff @ v_soft
for tag, M in (("amg", M_amg), ("ilu", M_ilu)):
    err = float(np.linalg.norm(M @ jv - v_soft) / np.linalg.norm(v_soft))
    metrics[f"soft_mode_inversion_err_{tag}"] = err
print(f"  soft mode: peak@x={metrics['soft_mode_x_at_peak']:.3f} support "
      f"x in [{metrics['soft_mode_x_min']:.3f},{metrics['soft_mode_x_max']:.3f}]"
      f" ({100*metrics['soft_mode_support_frac']:.2f}% of dofs) rhs_loading="
      f"{metrics['rhs_soft_mode_loading']:.4f}", flush=True)
print(f"  inversion err: amg={metrics['soft_mode_inversion_err_amg']:.4e} "
      f"ilu={metrics['soft_mode_inversion_err_ilu']:.4e}", flush=True)

# ---- long GMRES: production per-step budget (restart 60 x 10 cycles)
print("== long GMRES at the freeze point ==", flush=True)
RTOL, RESTART, NINNER = 1e-12, 60, 600


def gmres_traced(A, b, M, tag):
    seq = []

    def cb(xk):
        rk = b - A @ xk
        tr = float(np.linalg.norm(rk))
        pr = float(np.linalg.norm(M @ rk)) if M is not None else tr
        seq.append((tr, pr))

    x, info = spla.gmres(A, b, M=M, rtol=RTOL, atol=0.0, restart=RESTART,
                         maxiter=NINNER, callback=cb, callback_type="x")
    err = float(np.linalg.norm(x - x_ref) / np.linalg.norm(x_ref))
    rel = seq[-1][0] / seq[0][0]
    print(f"  {tag:14s}: info={info} iters={len(seq)} rel_true={rel:.4e} "
          f"head={[f'{s[0]:.2e}' for s in seq[:4]]} "
          f"tail={[f'{s[0]:.2e}' for s in seq[-3:]]} x_err={err:.3e}",
          flush=True)
    metrics[f"gmres600_{tag}_iters"] = len(seq)
    metrics[f"gmres600_{tag}_info"] = int(info)
    metrics[f"gmres600_{tag}_rel_true_res"] = rel
    metrics[f"gmres600_{tag}_x_err_vs_direct"] = err
    return seq


legs = [("none", J_ff, None), ("amg", J_ff, M_amg), ("ilu", J_ff, M_ilu),
        ("amg_on_picard", A_ff, M_amg)]
seqs = {tag: gmres_traced(A, rhs, M, tag) for tag, A, M in legs}

with open(os.path.join(OUT, "gmres_freeze_point.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["precond", "iter", "true_res", "precond_res"])
    for tag, _, _ in legs:
        for i, (tr, pr) in enumerate(seqs[tag]):
            w.writerow([tag, i, f"{tr:.16e}", f"{pr:.16e}"])

# merge with probe-1 metrics (they were measured on the same freeze point)
old = {}
with open(os.path.join(OUT, "freeze_point_metrics.csv")) as fh:
    for row in csv.reader(fh):
        if row and row[0] != "metric":
            old[row[0]] = row[1]
old.update({k: v for k, v in metrics.items()})
old.pop("svd_sigma_max", None)
with open(os.path.join(OUT, "freeze_point_metrics.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["metric", "value"])
    for k, v in old.items():
        w.writerow([k, v])

# ---- figure: Newton histories (probe-1 CSV) + long GMRES
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

nh = {"amg": [], "ilu": [], "direct": []}
with open(os.path.join(OUT, "newton_histories.csv")) as fh:
    for row in csv.reader(fh):
        if row and row[0] in nh:
            nh[row[0]].append(float(row[2]))

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for pc, mk in (("amg", "o-"), ("ilu", "s-"), ("direct", "^-")):
    ax[0].semilogy(range(len(nh[pc])), nh[pc], mk, ms=3, lw=1,
                   label=f"precond={pc}")
ax[0].set_xlabel("Newton step")
ax[0].set_ylabel(r"$|R|_\infty$")
ax[0].set_title(f"C05 nozzle nx={NX}: production legs (W2.2b recipe)")
ax[0].grid(alpha=0.3, which="both")
ax[0].legend(fontsize=8)

style = {"none": ("0.5", ":", "none (J_ff)"),
         "amg": ("r", "-", "amg on A_ff (J_ff)"),
         "ilu": ("b", "-", "ilu on J_ff (J_ff)"),
         "amg_on_picard": ("g", "--", "amg on A_ff (A_ff control)")}
for tag, _, _ in legs:
    tr = [s[0] for s in seqs[tag]]
    c, ls, lab = style[tag]
    ax[1].semilogy(range(len(tr)), np.asarray(tr) / tr[0], ls, color=c, lw=1.0,
                   label=lab)
pr_amg = [s[1] for s in seqs["amg"]]
ax[1].semilogy(range(len(pr_amg)),
               np.asarray(pr_amg) / [s[0] for s in seqs["amg"]][0],
               "m--", lw=0.9, label="amg: PRECONDITIONED res")
for cyc in range(1, NINNER // RESTART):
    ax[1].axvline(cyc * RESTART, color="0.85", lw=0.5, zorder=0)
ax[1].set_xlabel(f"GMRES inner iteration (restart={RESTART}, vertical lines "
                 "= restarts)")
ax[1].set_ylabel(r"$\|r_k\|/\|r_0\|$")
ax[1].set_title(f"freeze point $|R|_\\infty$"
                f"={float(old['freeze_r_inf']):.2e}, "
                f"cond2={metrics['cond2_Jff']:.2e}")
ax[1].grid(alpha=0.3, which="both")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "c05_amg_gmres.png"), dpi=140)
print("DONE", flush=True)
