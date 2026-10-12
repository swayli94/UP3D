"""Q1 probe (pre-registration 20261012-1000, Q1): why does precond="amg" stall
at |R|~3.9e-6 on the C05 nozzle production leg while direct/ilu converge to
machine zero.

Three parts, one freeze point:
  1. Reproduce the production legs (amg 40-step stall; ilu/direct machine zero).
  2. At the amg leg's terminal (stalled) state, extract J_ff via the production
     workspace (eval_residual + assemble_coupled; wc=None so the GMRES operator
     IS J_ff, no low-rank term) and the Picard block A_ff the AMG hierarchy is
     built on. Kill: assemble_coupled raising or splu(J_ff) failing.
  3. GMRES with per-iteration residual recording on J_ff x = -R_free for three
     preconditioners (none / amg-on-A_ff / ilu-on-J_ff) plus the control
     GMRES+amg on A_ff itself; direct splu as reference. Discriminator:
     amg flat in the first iterations => (1) preconditioner never approximates
     the shock block; slow linear decay => (2) conditioning; preconditioned
     residual claims convergence while the true residual stalls => (3).

Writes CSVs + PNG + metrics into bench/studies/c05_amg_investigation/.
Env: conda up3d, PYTHONNOUSERSITE=1, NUMBA/OMP/OPENBLAS threads 16.
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
os.makedirs(OUT, exist_ok=True)

M_INF, UPWIND_C, NX = 0.80, 1.5, 200
NY = max(6, NX // 16)
M_CRIT, M_CAP, RHO_FLOOR = 0.95, 3.0, 0.05

print("== setup: exact solution + mesh ==", flush=True)
ex = NZ.exact_solution(M_INF)
mesh = NZ.nozzle_mesh(NX, NY)
dir_nodes = np.unique(np.concatenate([
    mesh.boundary_faces["inlet"].ravel(),
    mesh.boundary_faces["outlet"].ravel()]))
dir_vals = ex["phi_of_x"](mesh.nodes[dir_nodes, 0])
phi0 = ex["phi_of_x"](mesh.nodes[:, 0])
print(f"mesh: {len(mesh.nodes)} nodes, {len(mesh.elements)} tets, "
      f"{len(dir_nodes)} dirichlet", flush=True)

# ---------------------------------------------------------------- part 1
print("== part 1: production legs (amg / ilu / direct) ==", flush=True)
hist, results = {}, {}
for pc in ("amg", "ilu", "direct"):
    r = solve_newton_nonlifting(
        mesh, m_inf=M_INF, upwind_c=UPWIND_C, dirichlet_nodes=dir_nodes,
        dirichlet_values=dir_vals, phi_init=phi0, entropy_correction=True,
        precond=pc, n_newton_max=40)
    hist[pc] = [float(v) for v in r["residual_history"]]
    results[pc] = r
    print(f"  {pc:6s}: converged={r['converged']} reason={r['accept_reason']} "
          f"n_newton={r['n_newton']} |R|_end={hist[pc][-1]:.3e} "
          f"n_gmres_total={r.get('n_gmres_total')}", flush=True)

with open(os.path.join(OUT, "newton_histories.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["leg", "step", "residual_inf"])
    for pc in ("amg", "ilu", "direct"):
        for i, v in enumerate(hist[pc]):
            w.writerow([pc, i, f"{v:.16e}"])

# ---------------------------------------------------------------- part 2
print("== part 2: freeze-point J_ff extraction (amg terminal state) ==",
      flush=True)
ws = results["amg"]["workspace"]
phi = np.asarray(results["amg"]["phi"], dtype=np.float64)
R_free, _, state = ws.eval_residual(phi[ws.free], np.zeros(0), UPWIND_C,
                                    M_CRIT, M_CAP, RHO_FLOOR)
r_inf_freeze = float(np.max(np.abs(R_free)))
try:
    J_ff, B = ws.assemble_coupled(state, UPWIND_C, M_CRIT, RHO_FLOOR)
except Exception as exc:  # kill criterion
    print(f"KILL: assemble_coupled raised {exc!r}", flush=True)
    sys.exit(2)
J_ff = J_ff.tocsr()
rhs = -R_free
print(f"  freeze |R|_inf={r_inf_freeze:.6e}  n_free={ws.n_free}  "
      f"nnz(J_ff)={J_ff.nnz}  finite={np.all(np.isfinite(J_ff.data))}",
      flush=True)

# the AMG block exactly as the driver builds it
A_pic = ws.op.assemble_matrix(state["rho_t"])
A_red = (ws.con.T.T @ (A_pic @ ws.con.T)).tocsr()
A_ff = A_red[ws.free][:, ws.free].tocsr()

try:
    lu_direct = spla.splu(J_ff.tocsc())
    x_ref = lu_direct.solve(rhs)
    assert np.all(np.isfinite(x_ref))
except Exception as exc:  # kill criterion
    print(f"KILL: splu(J_ff) failed: {exc!r}", flush=True)
    sys.exit(3)
print("  splu(J_ff) OK (no exact singularity)", flush=True)


def mmatrix_metrics(A, name):
    d = A.diagonal()
    off = A.tocsr().copy()
    off.setdiag(0)
    off.eliminate_zeros()
    abs_off = np.abs(off)
    row_off = np.asarray(abs_off.sum(axis=1)).ravel()
    weak = row_off > np.abs(d)
    At = (A - A.T).tocsr()
    num = spla.norm(At)
    den = spla.norm(A)
    return {
        f"{name}_n": A.shape[0],
        f"{name}_nnz": A.nnz,
        f"{name}_diag_min": float(d.min()),
        f"{name}_diag_max": float(d.max()),
        f"{name}_n_neg_diag": int((d < 0).sum()),
        f"{name}_frac_pos_offdiag": float((off.data > 0).sum()
                                          / max(1, off.data.size)),
        f"{name}_frac_weak_diagdom_rows": float(np.count_nonzero(weak)
                                                / A.shape[0]),
        f"{name}_skew_fro_ratio": float(num / den),
    }


metrics = {"freeze_r_inf": r_inf_freeze, "n_free": int(ws.n_free)}
metrics.update(mmatrix_metrics(A_ff, "A_ff_picard"))
metrics.update(mmatrix_metrics(J_ff, "J_ff"))
metrics["rho_t_min"] = float(np.min(state["rho_t"]))
metrics["rho_t_max"] = float(np.max(state["rho_t"]))
print(f"  A_ff: diag[{metrics['A_ff_picard_diag_min']:.3e},"
      f"{metrics['A_ff_picard_diag_max']:.3e}] pos_offdiag="
      f"{metrics['A_ff_picard_frac_pos_offdiag']:.4f} skew="
      f"{metrics['A_ff_picard_skew_fro_ratio']:.3e}", flush=True)
print(f"  J_ff: diag[{metrics['J_ff_diag_min']:.3e},"
      f"{metrics['J_ff_diag_max']:.3e}] neg_diag={metrics['J_ff_n_neg_diag']} "
      f"pos_offdiag={metrics['J_ff_frac_pos_offdiag']:.4f} skew="
      f"{metrics['J_ff_skew_fro_ratio']:.3e}", flush=True)

# pyamg hierarchy metrics
ml, M_amg = build_amg_preconditioner(A_ff)
lv_n = [lv.A.shape[0] for lv in ml.levels]
lv_nnz = [lv.A.nnz for lv in ml.levels]
metrics["amg_n_levels"] = len(ml.levels)
metrics["amg_level_n"] = ";".join(map(str, lv_n))
metrics["amg_level_nnz"] = ";".join(map(str, lv_nnz))
metrics["amg_operator_complexity"] = float(sum(lv_nnz) / lv_nnz[0])
# coarse-level size IS the level-0 aggregate count; B kept per level by SA
metrics["amg_n_aggregates_lvl0"] = (int(ml.levels[1].A.shape[0])
                                    if len(ml.levels) > 1
                                    else int(ml.levels[0].A.shape[0]))
B0 = getattr(ml.levels[0], "B", None)
metrics["amg_nearnull_dim"] = (int(B0.shape[1]) if B0 is not None else -1)
print(f"  amg: levels={lv_n} opcx={metrics['amg_operator_complexity']:.3f} "
      f"agg0={metrics['amg_n_aggregates_lvl0']} B_dim="
      f"{metrics['amg_nearnull_dim']}", flush=True)

M_ilu = build_ilu_preconditioner(J_ff)

# conditioning + soft mode (real-arithmetic SVD path)
print("== conditioning / soft mode ==", flush=True)
try:
    smax = float(spla.svds(J_ff, k=1, which="LM",
                           return_singular_vectors=False)[0])
    u_s, s_s, vt_s = spla.svds(J_ff, k=6, which="SM", maxiter=2000)
    smin = float(s_s[0])
    v_soft = vt_s[0].copy()
    metrics["svd_sigma_max"] = smax
    metrics["svd_sigma_min"] = smin
    metrics["svd_cond2"] = smax / smin
    metrics["svd_sigma_head"] = ";".join(f"{v:.4e}" for v in s_s[:4])
    # soft-mode support: map free index -> node -> x coordinate
    v_abs = np.abs(v_soft)
    node_of_free = ws.free
    x_nodes = mesh.nodes[:, 0]
    i_peak = int(np.argmax(v_abs))
    metrics["soft_mode_x_at_peak"] = float(x_nodes[node_of_free[i_peak]])
    # top-1% support spread
    thr = 0.01 * v_abs.max()
    supp = v_abs >= thr
    metrics["soft_mode_support_frac"] = float(np.count_nonzero(supp)
                                              / len(v_abs))
    metrics["soft_mode_x_min"] = float(x_nodes[node_of_free[supp]].min())
    metrics["soft_mode_x_max"] = float(x_nodes[node_of_free[supp]].max())
    # inversion quality: M ~= J_ff^{-1} should map J_ff v back to v
    jv = J_ff @ v_soft
    for tag, M in (("amg", M_amg), ("ilu", M_ilu)):
        err = float(np.linalg.norm(M @ jv - v_soft) / np.linalg.norm(v_soft))
        metrics[f"soft_mode_inversion_err_{tag}"] = err
        print(f"  soft mode: sigma_min={smin:.4e} cond2={smax/smin:.3e} "
              f"peak@x={metrics['soft_mode_x_at_peak']:.3f} "
              f"inv_err[{tag}]={err:.4e}", flush=True)
except Exception as exc:
    print(f"  svds failed (non-fatal): {exc!r}", flush=True)
    metrics["svd_sigma_max"] = float("nan")

# ---------------------------------------------------------------- part 3
print("== part 3: GMRES per-iteration residuals at the freeze point ==",
      flush=True)
RTOL, RESTART, MAXCYC = 1e-12, 60, 40


def gmres_traced(A, b, M, tag):
    seq = []

    def cb(xk):
        rk = b - A @ xk
        tr = float(np.linalg.norm(rk))
        pr = float(np.linalg.norm(M @ rk)) if M is not None else tr
        seq.append((tr, pr))

    x, info = spla.gmres(A, b, M=M, rtol=RTOL, atol=0.0, restart=RESTART,
                         maxiter=MAXCYC, callback=cb, callback_type="x")
    err = float(np.linalg.norm(x - x_ref) / np.linalg.norm(x_ref))
    print(f"  {tag:14s}: info={info} iters={len(seq)} "
          f"true_res {seq[0][0]:.3e} -> {seq[-1][0]:.3e} "
          f"(rel {seq[-1][0]/seq[0][0]:.3e}) x_err={err:.3e}", flush=True)
    metrics[f"gmres_{tag}_iters"] = len(seq)
    metrics[f"gmres_{tag}_info"] = int(info)
    metrics[f"gmres_{tag}_rel_true_res"] = seq[-1][0] / seq[0][0]
    metrics[f"gmres_{tag}_x_err_vs_direct"] = err
    return x, info, seq


legs = [
    ("none", J_ff, rhs, None),
    ("amg", J_ff, rhs, M_amg),
    ("ilu", J_ff, rhs, M_ilu),
    ("amg_on_picard", A_ff, rhs, M_amg),   # control: hierarchy health
]
gmres_seqs = {}
for tag, A, b, M in legs:
    _, _, gmres_seqs[tag] = gmres_traced(A, b, M, tag)

with open(os.path.join(OUT, "gmres_freeze_point.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["precond", "iter", "true_res", "precond_res"])
    for tag, _, _, _ in legs:
        for i, (tr, pr) in enumerate(gmres_seqs[tag]):
            w.writerow([tag, i, f"{tr:.16e}", f"{pr:.16e}"])

with open(os.path.join(OUT, "freeze_point_metrics.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["metric", "value"])
    for k, v in metrics.items():
        w.writerow([k, v])

# ---------------------------------------------------------------- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for pc, mk in (("amg", "o-"), ("ilu", "s-"), ("direct", "^-")):
    ax[0].semilogy(range(len(hist[pc])), hist[pc], mk, ms=3, lw=1,
                   label=f"precond={pc}")
ax[0].set_xlabel("Newton step")
ax[0].set_ylabel(r"$|R|_\infty$")
ax[0].set_title(f"C05 nozzle nx={NX}: production legs")
ax[0].grid(alpha=0.3, which="both")
ax[0].legend(fontsize=8)

style = {"none": ("0.5", ":", "none (J_ff)"),
         "amg": ("r", "-", "amg on A_ff (J_ff)"),
         "ilu": ("b", "-", "ilu on J_ff (J_ff)"),
         "amg_on_picard": ("g", "--", "amg on A_ff (A_ff control)")}
for tag, _, _, _ in legs:
    tr = [s[0] for s in gmres_seqs[tag]]
    c, ls, lab = style[tag]
    ax[1].semilogy(range(len(tr)), np.asarray(tr) / tr[0], ls, color=c, lw=1.2,
                   label=lab)
pr_amg = [s[1] for s in gmres_seqs["amg"]]
ax[1].semilogy(range(len(pr_amg)),
               np.asarray(pr_amg) / [s[0] for s in gmres_seqs["amg"]][0],
               "m--", lw=0.9, label="amg: PRECONDITIONED res")
ax[1].set_xlabel(f"GMRES inner iteration (restart={RESTART})")
ax[1].set_ylabel(r"$\|r_k\|/\|r_0\|$")
ax[1].set_title(f"freeze point $|R|_\\infty$={r_inf_freeze:.2e}: "
                "preconditioner comparison")
ax[1].grid(alpha=0.3, which="both")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "c05_amg_gmres.png"), dpi=140)
print(f"artifacts written to {OUT}", flush=True)
print("DONE", flush=True)
