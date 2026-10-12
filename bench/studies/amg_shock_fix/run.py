r"""「amg 能否修」研究（phase 7 第 4 轮，S3）。预注册：
`docs/dev_phase_seven/20261012-2230-amg-shock-fix-prereg.md`。

机制（p6 第 47 轮已闭合）：amg 的层级建在不含 Term-3 激波耦合的 Picard 块 A_ff 上，
对 J_ff 的激波中心近零模态（σ_min 2.75e-6，右奇异矢峰在激波 x=11.4）无反演。
本轮测两个可修方向 + 代价：

- **方向 1**：AMG 直接建在完整 J_ff 上（非对称、略密）——pyamg 建不建得动、
  建了收不收敛、层级贵多少。
- **方向 2**：AMG(A_ff) 打底 + 激波方向低秩校正（加性 deflation：
  M⁻¹r = M_amg r + W (ZᵀJW)⁻¹ Zᵀ r，W/Z = J_ff 的右/左软奇异矢，k ≤ 6）。

测量台 = p6 探针同一冻结点（C05 生产腿 amg 末态，n_free=5174）+ 真实生产 Newton 腿
（C05 喷管 nx=200，本脚本自写的**简化 Newton 环**：无 EW/线搜索/冻结 —— 与生产
驱动器的差异如实写进 README；对照臂 = 同一环上的 direct）。
**不改库、不改默认值。**

用法：`python run.py`（求解 + CSV + PNG）；环境：conda up3d、PYTHONNOUSERSITE=1、
NUMBA_NUM_THREADS=16；BLAS 由 pyfp3d 自身钉 1（断言）。
"""
import csv
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, REPO)
os.chdir(REPO)

from tests import _nozzle_case as NZ                       # noqa: E402
from pyfp3d import blas_thread_count                       # noqa: E402
from pyfp3d.solve.linear import (build_amg_preconditioner,  # noqa: E402
                                 build_ilu_preconditioner)
from pyfp3d.solve.newton import solve_newton_nonlifting    # noqa: E402

import pyamg                                               # noqa: E402

M_INF, UPWIND_C, NX = 0.80, 1.5, 200
NY = max(6, NX // 16)
M_CRIT, M_CAP, RHO_FLOOR = 0.95, 3.0, 0.05
RTOL, RESTART, MAXCYC = 1e-12, 60, 10     # 10×60 = 600 内迭代 = 生产单步预算
NEWTON_TOL, NEWTON_MAX = 1e-11, 40


def gmres_traced(A, b, M, x_ref, budget_inner=600):
    """GMRES，逐内迭代记真残差；返回 (x, info, seq, n_inner)。"""
    seq = []

    def cb(xk):
        seq.append(float(np.linalg.norm(b - A @ xk)))

    x, info = spla.gmres(A, b, M=M, rtol=RTOL, atol=0.0, restart=RESTART,
                         maxiter=MAXCYC, callback=cb, callback_type="x")
    err = float(np.linalg.norm(x - x_ref) / np.linalg.norm(x_ref))
    return x, info, seq, len(seq), err


def amg_metrics(ml):
    lv_n = [lv.A.shape[0] for lv in ml.levels]
    lv_nnz = [lv.A.nnz for lv in ml.levels]
    return dict(n_levels=len(lv_n), levels=";".join(map(str, lv_n)),
                nnz=";".join(map(str, lv_nnz)),
                opcx=float(sum(lv_nnz) / lv_nnz[0]))


def soft_modes(J, k):
    """J 的底部 k 个奇异三要素（σ, u, v）。svds-SM 在 J_ff 上不收敛
    （p6 probe1/probe2 实测）⇒ 走 LU 预条件 eigsh：C = J⁻ᵀJ⁻¹ = (JJᵀ)⁻¹ 是 SPD，
    顶部特征对给 σ_min 一侧；u = 特征矢，v = Jᵀu/σ。v0 钉死（可复算）。"""
    n = J.shape[0]
    lu = spla.splu(J.tocsc())
    Jinv = spla.LinearOperator((n, n), matvec=lu.solve,
                               rmatvec=lambda x: lu.solve(x, trans="T"))
    Cop = spla.LinearOperator((n, n),
                              matvec=lambda x: Jinv.rmatvec(Jinv.matvec(x)))
    lam, U = spla.eigsh(Cop, k=k, which="LA", maxiter=2000,
                        v0=np.random.default_rng(0).standard_normal(n))
    sig = 1.0 / np.sqrt(lam)
    V = np.column_stack([J.T @ U[:, i] / sig[i] for i in range(k)])
    return sig, U, V


def make_deflated(M_amg, J, k, W6, Z6):
    """M⁻¹ r = M_amg r + W (ZᵀJW)⁻¹ Zᵀ r（W/Z = 前 k 个右/左软奇异矢）。"""
    W, Z = W6[:, :k], Z6[:, :k]
    E = np.asarray(Z.T @ J @ W)
    E_inv = np.linalg.inv(E)

    def mv(r):
        return M_amg @ r + W @ (E_inv @ (Z.T @ r))

    return spla.LinearOperator(J.shape, matvec=mv)


def main():
    assert blas_thread_count() == 1, "证据纪律：BLAS 必须是库钉死的 1"
    print("== setup：精确解 + 网格 + 生产 amg 腿（复现停滞，取冻结点）==",
          flush=True)
    ex = NZ.exact_solution(M_INF)
    mesh = NZ.nozzle_mesh(NX, NY)
    dir_nodes = np.unique(np.concatenate([
        mesh.boundary_faces["inlet"].ravel(),
        mesh.boundary_faces["outlet"].ravel()]))
    dir_vals = ex["phi_of_x"](mesh.nodes[dir_nodes, 0])
    phi0 = ex["phi_of_x"](mesh.nodes[:, 0])

    r_amg = solve_newton_nonlifting(
        mesh, m_inf=M_INF, upwind_c=UPWIND_C, dirichlet_nodes=dir_nodes,
        dirichlet_values=dir_vals, phi_init=phi0, entropy_correction=True,
        precond="amg", n_newton_max=NEWTON_MAX)
    print(f"  生产 amg 腿：converged={r_amg['converged']} "
          f"|R|_end={r_amg['residual_history'][-1]:.3e} "
          f"n_newton={r_amg['n_newton']}", flush=True)

    ws = r_amg["workspace"]
    phi = np.asarray(r_amg["phi"], dtype=np.float64)
    R_free, _, state = ws.eval_residual(phi[ws.free], np.zeros(0), UPWIND_C,
                                        M_CRIT, M_CAP, RHO_FLOOR)
    J_ff = ws.assemble_coupled(state, UPWIND_C, M_CRIT, RHO_FLOOR)[0].tocsr()
    A_pic = ws.op.assemble_matrix(state["rho_t"])
    A_ff = (ws.con.T.T @ (A_pic @ ws.con.T)).tocsr()[ws.free][:, ws.free].tocsr()
    rhs = -R_free
    x_ref = spla.splu(J_ff.tocsc()).solve(rhs)
    print(f"  冻结点 |R|_inf={np.max(np.abs(R_free)):.3e} n_free={J_ff.shape[0]} "
          f"nnz(J_ff)={J_ff.nnz} nnz(A_ff)={A_ff.nnz}", flush=True)

    cost = {}

    # ---------------------------------------------------------- 方向 1
    print("== C1 方向 1：AMG 建在完整 J_ff 上 ==", flush=True)
    M_amgA = build_amg_preconditioner(A_ff)[1]
    t = []
    for _ in range(3):
        t0 = time.perf_counter()
        ml_A = pyamg.smoothed_aggregation_solver(A_ff)
        t.append(time.perf_counter() - t0)
    cost["amg_A_setup_s"] = float(np.median(t))
    m_A = amg_metrics(ml_A)
    try:
        t = []
        for _ in range(3):
            t0 = time.perf_counter()
            ml_J = pyamg.smoothed_aggregation_solver(J_ff)
            t.append(time.perf_counter() - t0)
        cost["amg_J_setup_s"] = float(np.median(t))
        m_J = amg_metrics(ml_J)
        finite = bool(np.all([np.all(np.isfinite(lv.A.data))
                              for lv in ml_J.levels]))
        print(f"  pyamg(J_ff) 建成：levels={m_J['levels']} opcx={m_J['opcx']:.3f} "
              f"finite={finite} setup={cost['amg_J_setup_s']*1e3:.1f} ms "
              f"(A_ff {cost['amg_A_setup_s']*1e3:.1f} ms)", flush=True)
        M_amgJ = ml_J.aspreconditioner()
        _, info1, seq1, n1, err1 = gmres_traced(J_ff, rhs, M_amgJ, x_ref)
        print(f"  GMRES+amg(J_ff)：iters={n1} info={info1} "
              f"rel_true={seq1[0] and seq1[-1]/seq1[0]:.3e} x_err={err1:.3e}",
              flush=True)
    except Exception as exc:  # kill：方向 1 建不动
        print(f"  KILL-C1：pyamg(J_ff) 抛异常：{exc!r}", flush=True)
        m_J, seq1, n1, err1 = None, None, -1, float("nan")
        cost["amg_J_setup_s"] = float("nan")

    # ---------------------------------------------------------- 方向 2
    print("== C2 方向 2：AMG(A_ff) 打底 + 软模态低秩校正 ==", flush=True)
    t0 = time.perf_counter()
    s6, U6, V6 = soft_modes(J_ff, 20)
    cost["softmodes20_s"] = time.perf_counter() - t0
    #: ★ eigsh(which="LA") 的 σ 降序、σ_min 在最后一列 ⇒ 反转成升序，
    #:   之后 [:, :k] 才是「最软的 k 个」（第一版没反转，扫的是软集里**最大**的
    #:   σ —— 换 k=20 时集合的头会动，两次扫描不可比；已在判定里登记）。
    s6, U6, V6 = s6[::-1], U6[:, ::-1], V6[:, ::-1]
    W6, Z6 = V6, U6
    print(f"  软模态 k=20（LU 预条件 eigsh）：σ 升序头="
          f"{';'.join(f'{v:.3e}' for v in s6[:4])} … σ_20={s6[-1]:.3e} "
          f"cost={cost['softmodes20_s']*1e3:.0f} ms", flush=True)
    defl_rows = []
    for k in (1, 2, 4, 6, 12, 20):
        M_d = make_deflated(M_amgA, J_ff, k, W6, Z6)
        _, info, seq, n, err = gmres_traced(J_ff, rhs, M_d, x_ref)
        rel = seq[-1] / seq[0]
        print(f"  defl k={k}：iters={n} info={info} rel_true={rel:.3e} "
              f"x_err={err:.3e}", flush=True)
        defl_rows.append((k, n, info, rel, err))

    # 对照（与 p6 探针同口径）：amg(A_ff) 与 ilu(J_ff)
    _, _, seq_amg, n_amg, err_amg = gmres_traced(J_ff, rhs, M_amgA, x_ref)
    M_ilu = build_ilu_preconditioner(J_ff)
    _, _, seq_ilu, n_ilu, err_ilu = gmres_traced(J_ff, rhs, M_ilu, x_ref)
    print(f"  对照：amg(A_ff) iters={n_amg} rel={seq_amg[-1]/seq_amg[0]:.3e}；"
          f"ilu(J_ff) iters={n_ilu} rel={seq_ilu[-1]/seq_ilu[0]:.3e}", flush=True)

    # ---------------------------------------------------------- C3 端到端
    print("== C3 端到端：简化 Newton 环（无 EW/线搜索/冻结，对照臂 direct）==",
          flush=True)

    def a_ff_of(st):
        A_pic = ws.op.assemble_matrix(st["rho_t"])
        return (ws.con.T.T @ (A_pic @ ws.con.T)).tocsr()[ws.free][:, ws.free].tocsr()

    def defl_preconditioner(J, st, k=4):
        """每步：amg(A_ff(st)) 打底 + 当前 J 的 k 个软模态校正（eigsh 钉 v0）。"""
        M_amg = build_amg_preconditioner(a_ff_of(st))[1]
        sig, U, V = soft_modes(J, k)
        return make_deflated(M_amg, J, k, V, U), "gmres"

    variants = [("ctrl_direct", lambda J, st: (None, "direct"))]
    if m_J is not None:
        variants.append(("amg_on_Jff",
                         lambda J, st: (pyamg.smoothed_aggregation_solver(J)
                                        .aspreconditioner(), "gmres")))
    variants.append(("defl_k4", defl_preconditioner))

    legs = []
    for tag, build_M in variants:
        phi_f = phi0[ws.free].copy()
        hist, t0 = [], time.perf_counter()
        for step in range(NEWTON_MAX):
            R_f, _, st = ws.eval_residual(phi_f, np.zeros(0), UPWIND_C,
                                          M_CRIT, M_CAP, RHO_FLOOR)
            rn = float(np.max(np.abs(R_f)))
            hist.append(rn)
            if rn < NEWTON_TOL:
                break
            J = ws.assemble_coupled(st, UPWIND_C, M_CRIT, RHO_FLOOR)[0].tocsr()
            M, kind = build_M(J, st)
            if kind == "direct":
                d = spla.splu(J.tocsc()).solve(-R_f)
            else:
                d, info = spla.gmres(J, -R_f, M=M, rtol=1e-10, atol=0.0,
                                     restart=RESTART, maxiter=MAXCYC)
            phi_f = phi_f + d
        wall = time.perf_counter() - t0
        ok = bool(hist[-1] < NEWTON_TOL)
        print(f"  {tag:18s}: steps={len(hist)} |R|_end={hist[-1]:.3e} "
              f"wall={wall:.2f} s converged={ok}", flush=True)
        legs.append(dict(tag=tag, hist=hist, wall=wall, converged=ok,
                         steps=len(hist)))

    # ---------------------------------------------------------- 落盘
    with open(os.path.join(HERE, "freeze_point_gmres.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["variant", "iter", "true_res"])
        for tag, seq in (("amg_on_Jff", seq1), ("amg_on_Aff", seq_amg),
                         ("ilu_on_Jff", seq_ilu)):
            if seq is None:
                continue
            for i, v in enumerate(seq):
                w.writerow([tag, i, f"{v:.16e}"])
    with open(os.path.join(HERE, "deflation_sweep.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["k", "iters", "info", "rel_true_res", "x_err_vs_direct"])
        for k, n, info, rel, err in defl_rows:
            w.writerow([k, n, info, f"{rel:.6e}", f"{err:.6e}"])
    with open(os.path.join(HERE, "newton_legs.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["variant", "step", "residual_inf"])
        for leg in legs:
            for i, v in enumerate(leg["hist"]):
                w.writerow([leg["tag"], i, f"{v:.16e}"])
    with open(os.path.join(HERE, "cost_table.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "value"])
        for k, v in cost.items():
            w.writerow([k, v])
        for leg in legs:
            w.writerow([f"newton_{leg['tag']}_wall_s", leg["wall"]])
            w.writerow([f"newton_{leg['tag']}_steps", leg["steps"]])
            w.writerow([f"newton_{leg['tag']}_converged", int(leg["converged"])])

    # ---------------------------------------------------------- 图
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ax[0].semilogy(range(len(seq_amg)), np.asarray(seq_amg) / seq_amg[0],
                   "r-", lw=1.2, label="amg(A_ff) - production shape (flat)")
    if seq1 is not None:
        ax[0].semilogy(range(len(seq1)), np.asarray(seq1) / seq1[0],
                       color="0.4", ls=":", lw=1.2,
                       label="amg(J_ff) - direction 1 (measured negative)")
    ax[0].semilogy(range(len(seq_ilu)), np.asarray(seq_ilu) / seq_ilu[0],
                   "b-", lw=1.2, label="ilu(J_ff) (control)")
    ks = [r[0] for r in defl_rows]
    rels = [r[3] for r in defl_rows]
    ax[0].set_xlabel("GMRES restart cycles (x60 inner)")
    ax[0].set_ylabel(r"$\|r_k\|/\|r_0\|$")
    ax[0].set_title(f"freeze point |R|_inf={np.max(np.abs(R_free)):.2e}: controls")
    ax[0].grid(alpha=.3, which="both"), ax[0].legend(fontsize=8)
    ax[1].semilogy(ks, rels, "go-", label="amg(A_ff) + soft-mode correction")
    ax[1].axhline(1e-6, color="r", ls="--", lw=1, label="1e-6")
    ax[1].set_xlabel("k (softest modes, ascending sigma)")
    ax[1].set_ylabel("rel true residual @600 inner")
    ax[1].set_title("direction 2: deflation sweep (freeze point)")
    ax[1].grid(alpha=.3, which="both"), ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "amg_shock_fix.png"), dpi=140)
    plt.close(fig)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
