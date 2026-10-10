r"""D 类门的**我方**三级网格收敛表 —— `grid_convergence_ours.csv` 的写入器（W2.1a）。

★ 为什么需要它：参照侧（CFL3D）每个量有四档 + 观察阶 + 三分判定，我方此前
全两档 —— 一个网格一个数，没有数值不确定度可读（合并分析
`docs/dev_phase_six/20260906-0100-gate-audit-synthesis.md` §5 W2.1 / H8）。
xcoarse 档进来后（每族三级），我方腿也能读一个**三级**的阶。

★ 仪器是参照侧的，**不重写**：`implied_order` import 自
`cases/reference_data/cfl3d/generate_m6_reference.py`（F08 用合成数据验过牙），
import 方式照 `tests/F/test_F08_cfl3d_reference.py`（sys.path 拼接 —— 那个
目录不是包，只有拼进去才能 import）。判定词汇与参照侧 `asymptotic` 列
逐字一致（见各 VERDICT_* 常量）。

★★ **三档只能给一个数**：看不到 `unstable` —— 参照侧 131 个量里 41 个的判定
翻转是加了第四档才看见的（`generate_cfl3d_reference.py` 的
`derive_grid_convergence` docstring）。所以本表是 **RECORDED，不进值锁**；
判据化（|S − D| vs √(u_num² + u_D²)）是 W2.1b 的事，不在本轮。
"""
import csv
import os
import sys

import numpy as np

#: 三分判定词汇，与参照侧逐字一致（`generate_cfl3d_reference.py` 的
#: `asymptotic` / `error_bar` 列措辞）。单三元组没有 `unstable` 可分。
VERDICT_GROWING = "no (deltas are GROWING)"
VERDICT_YES = "yes (3-rung)"
VERDICT_UNCHECKED = "yes (3-rung, unchecked)"
VERDICT_UNDEFINED = "undefined"

#: 进表的量 = 门 `MEASURED` 里的物理列（`blas_threads` / `n_threads` /
#: `residual` 这类出处与诊断列不算）。
QUANTITIES = ("cl", "cd", "cp_rms_upper", "cp_rms_lower", "x_shock")


def _implied_order(h, vals):
    """参照侧估计器本身（`generate_m6_reference.implied_order`）。"""
    root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    cfl3d = os.path.join(root, "cases", "reference_data", "cfl3d")
    if cfl3d not in sys.path:
        sys.path.insert(0, cfl3d)
    import generate_m6_reference as G
    return G.implied_order(np.asarray(h, float), np.asarray(vals, float), 0)


def write_grid_convergence_ours(path, runs, cases, levels, mesh_dir,
                                quantities=QUANTITIES):
    """每 (case, quantity) 一行：三级读数、相邻差、比值、implied_order、三分判定。

    `levels` 必须就是 (xcoarse, coarse, medium) —— 列名写死这三级。
    h 由三级的 `n_nodes`（`<mesh_dir>/<level>_stats.csv`）按 2-D 的
    h ∝ N^(−1/2) 构造（归一化无关，估计器只用比值；3-D 用 N^(−1/3)，别混）。
    任何一腿不收敛 ⇒ `basis=not_converged` 且 ratio/order/verdict 留空
    （不收敛腿的阶无意义 —— 参照侧对读不出的量同样不给判定）。
    """
    assert tuple(levels) == ("xcoarse", "coarse", "medium"), levels

    def f9(v):
        return f"{v:.9e}" if np.isfinite(v) else "nan"

    n = {}
    for lv in levels:
        with open(os.path.join(mesh_dir, f"{lv}_stats.csv")) as fh:
            n[lv] = float(next(r[1] for r in csv.reader(fh)
                               if r and r[0] == "n_nodes"))
    h = np.array([n[lv] for lv in levels]) ** -0.5

    with open(path, "w", newline="") as fh:
        fh.write("# 我方腿的三级网格收敛读数（W2.1a，RECORDED，不进值锁）。\n")
        fh.write("# ★ 三档只能给一个数，看不到 unstable：参照侧 131 个量里 41 个的\n")
        fh.write("#   判定翻转是加第四档才看见的 —— 本表的 yes 都带 (3-rung) 限定。\n")
        fh.write("# h ∝ N^(−1/2)，N = 三级 n_nodes（" + " / ".join(
            f"{lv} {int(n[lv])}" for lv in levels) + "）。\n")
        w = csv.writer(fh)
        w.writerow(["case", "quantity", "xcoarse", "coarse", "medium",
                    "delta_xc_c", "delta_c_m", "ratio", "implied_order",
                    "verdict", "basis"])
        for case in cases:
            for q in quantities:
                vals = [float(runs[(case, lv)][q]) for lv in levels]
                conv = all(bool(runs[(case, lv)]["converged"]) for lv in levels)
                basis = "xc-c-m" if conv else "not_converged"
                xc, c, m = vals
                d_xc_c, d_c_m = c - xc, m - c
                if not conv:
                    ratio_s = order_s = verdict = ""
                else:
                    ratio, p = _implied_order(h, vals)
                    if not np.isfinite(ratio):
                        #: 估计器的两个哨兵，参照侧两者都判 undefined：
                        #: 输入带 nan（如 x_shock 未检出）→ nan；
                        #: delta_xc_c == 0（估计器返回 inf）→ 差值精确为零，无可比。
                        ratio_s = ("nan" if ratio != ratio else
                                   "UNDEFINED (coarse delta = 0)")
                        order_s, verdict = "", VERDICT_UNDEFINED
                    else:
                        ratio_s = f"{ratio:.3f}"
                        order_s = "" if p != p else f"{p:.2f}"
                        verdict = (VERDICT_GROWING if ratio >= 1.0 else
                                   VERDICT_UNCHECKED if ratio < 0.2 else
                                   VERDICT_YES)
                w.writerow([case, q, f9(xc), f9(c), f9(m),
                            f9(d_xc_c) if not np.isfinite(d_xc_c)
                            else f"{d_xc_c:+.9e}",
                            f9(d_c_m) if not np.isfinite(d_c_m)
                            else f"{d_c_m:+.9e}",
                            ratio_s, order_s, verdict, basis])
