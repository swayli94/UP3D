r"""hex 族 NACA 气泡缺陷的**残余拆分探针**（phase 7 第 3 轮顺带项；时间盒 = 两臂）。

残余（p6 第 47 轮判定 §2.2 如实登记）：缺陷在 hex 生成器（拓扑/分布/拼缝）还是在
该拓扑上被触发的求解器熵冻结，(a)(b)(c) 三件仪器分不开 —— hex 各腿 sigma 冻结搅动
|dσ| 0.19–0.21 vs tet 0.067。本探针用**动配方的两条臂**试分（hex medium，
M0.80/α1.25，D05 生产配方）：

- **臂 E0**：`entropy_correction=False`（其余同 RECIPE）—— 关掉整个熵修正；
- **臂 F0**：`freeze_tol=None`（其余同 RECIPE）—— 只关冻结，熵修正照常跑。

G-DOMAIN（看数之前写死，预注册 20261012-2135 §1.2）：

| 结果 | 读法 |
|---|---|
| E0 气泡消失且亏损闭合 | 熵修正/冻结链路被牵连 ⇒ 残余偏向求解器侧 |
| E0 气泡仍在且亏损仍在 | 缺陷不依赖熵修正 ⇒ 残余偏向生成器/拓扑侧 |
| F0 气泡消失 | 冻结交互（锁住错误结构）被单独牵连 |
| 臂不收敛 / 气泡移动但仍在 / 两臂矛盾 | RECORDED，残余保持未分开，不升级仪器 |

生产臂（RECIPE 原样）**不重算**：锚 = 已提交的 D15 门证据
`cases/gates/D15_hex_naca_boundary/summary.csv` 的 `hex|medium` 行（证据纪律 2）。

用法：`python run.py`（两臂求解 + arms.csv）；`python run.py --analysis-only` 只读 CSV。
环境：conda up3d；NUMBA_NUM_THREADS=16；BLAS 由 pyfp3d 自身钉 1（断言）。
"""
import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "tests", "D"))
sys.path.insert(0, REPO)

#: 生产配方与气泡计/分段仪器**从 D15 门文件导入**，不抄第二份 —— 探针与门量的是同一个量
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "d15", os.path.join(REPO, "tests", "D", "test_D15_hex_naca_boundary.py"))
D15 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(D15)


def run_arm(name, **override):
    """一条臂：D15 的 `_one` 同构，但配方带 override。hex medium。"""
    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.mesh.wake_cut import cut_wake
    from pyfp3d.post.section_cut import section_cp_curve
    from pyfp3d.post.shock import shock_report
    from pyfp3d.post.surface import wall_force_coefficients
    from pyfp3d.solve.newton import solve_newton_lifting

    kw = dict(D15.RECIPE)
    kw.update(override)
    mc, wc = cut_wake(read_mesh(os.path.join(
        REPO, "cases", "meshes", "naca0012_hex_2.5d", "medium.msh")))
    dz = float(np.ptp(mc.nodes[:, 2]))
    r = solve_newton_lifting(mc, wc, m_inf=D15.M_INF, alpha_deg=D15.ALPHA, **kw)
    phi = np.asarray(r["phi"])
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], phi,
                                alpha_deg=D15.ALPHA, u_inf=1.0, s_ref=dz,
                                m_inf=D15.M_INF)
    cur = section_cp_curve(mc, phi, z=float(np.mean(mc.nodes[:, 2])),
                           smooth_passes=1, m_inf=D15.M_INF)
    up = shock_report(cur, D15.M_INF)["upper"]
    frz = r.get("sigma_freeze_report") or {}
    hist = np.asarray(r["residual_history"], float)
    bubble_cp, bubble_x = D15._bubble(cur)
    return dict(arm=name, cl=float(f["cl"]),
                x_shock=(float(up["x_shock"]) if up.get("has_shock")
                         else float("nan")),
                bubble_max_cp=bubble_cp, bubble_x_at_max=bubble_x,
                converged=bool(r.get("converged")),
                residual=float(hist[-1]), n_newton=int(len(hist)),
                n_limited=int(r.get("n_limited", 0)),
                n_floored=int(r.get("n_floored", 0)),
                sigma_frozen=bool(frz.get("frozen_in_transient", False)),
                blas_threads=int(__import__("pyfp3d").blas_thread_count()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis-only", action="store_true")
    args = ap.parse_args()
    path = os.path.join(HERE, "arms.csv")

    from pyfp3d import blas_thread_count
    assert blas_thread_count() == 1, "证据纪律：BLAS 必须是库钉死的 1"

    if not args.analysis_only:
        rows = [run_arm("E0_entropy_off", entropy_correction=False),
                run_arm("F0_freeze_off", freeze_tol=None)]
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["arm", "blas_threads", "cl", "x_shock", "bubble_max_cp",
                        "bubble_x_at_max", "converged", "residual", "n_newton",
                        "n_limited", "n_floored", "sigma_frozen"])
            for r in rows:
                w.writerow([r["arm"], r["blas_threads"],
                            f"{r['cl']:.9e}", f"{r['x_shock']:.9e}",
                            f"{r['bubble_max_cp']:.9e}", f"{r['bubble_x_at_max']:.9e}",
                            int(r["converged"]), f"{r['residual']:.9e}",
                            r["n_newton"], r["n_limited"], r["n_floored"],
                            int(r["sigma_frozen"])])

    rows = {r["arm"]: r for r in csv.DictReader(open(path))}
    #: 锚 = 已提交 D15 证据的 hex medium 行（不重算，纪律 2）
    anchor = {r["leg"]: r for r in csv.DictReader(open(os.path.join(
        REPO, "cases", "gates", "D15_hex_naca_boundary", "summary.csv")))}
    a = anchor["hex|medium"]
    print(f"anchor (committed D15 hex|medium): cl={a['cl']} "
          f"bubble={a['bubble_max_cp']} x_shock={a['x_shock']}")
    for name, r in rows.items():
        print(f"{name}: cl={float(r['cl']):.6f} bubble={float(r['bubble_max_cp']):+.4f} "
              f"@x={float(r['bubble_x_at_max']):.3f} x_shock={float(r['x_shock']):.4f} "
              f"converged={r['converged']} |R|={float(r['residual']):.2e} "
              f"n_newton={r['n_newton']} lim/flr={r['n_limited']}/{r['n_floored']} "
              f"sigma_frozen={r['sigma_frozen']}")


if __name__ == "__main__":
    main()
