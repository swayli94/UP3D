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
`derive_grid_convergence` docstring）。所以本表是 **RECORDED，不进值锁**。

★★★ **W2.1b（2026-10-11）**：不确定度**接进判据** —— D05/D06 的 `|S − D|` 腿
改写为 `|S − D| ≤ min(标定带, u_val)`，`u_val = √(u_num² + u_D²)`
（ASME V&V 20），`u_num = |S_coarse − S_medium|`（与参照侧 `error_bar` 的
「最细区间差」口径逐字一致；**不用 Richardson 外推** —— 阶在多半腿上读不出来，
外推会把不可判写成可判）。判定函数（`ladder_readout` / `u_num_from_ladder` /
`u_val` / `band_for` / `regime` / `leg_uncertainty`）就在下面，**输入全部可注入**
（G-TEETH 用扰动副本喂，零求解）。本表随之加 `u_num / u_D / u_val / regime`
四列。**min 语义：任何腿的实际带永不变宽。**
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


def mesh_node_counts(mesh_dir, levels=("xcoarse", "coarse", "medium")):
    """三级的 `n_nodes`（dict）—— 读 `<mesh_dir>/<level>_stats.csv`。"""
    n = {}
    for lv in levels:
        with open(os.path.join(mesh_dir, f"{lv}_stats.csv")) as fh:
            n[lv] = float(next(r[1] for r in csv.reader(fh)
                               if r and r[0] == "n_nodes"))
    return n


def mesh_h(mesh_dir, levels=("xcoarse", "coarse", "medium")):
    """h ∝ N^(−1/2)，N = 三级 `n_nodes`。归一化无关（估计器只用比值）；
    **2-D 用 −1/2**，3-D 用 −1/3，别混。"""
    n = mesh_node_counts(mesh_dir, levels)
    return np.array([n[lv] for lv in levels]) ** -0.5


def ladder_readout(vals, h):
    """(verdict, ratio, p) —— 估计器调用 + 三分词汇映射；writer 与门共用同一台仪器。

    三档看不到 unstable（那要两个三元组 = 第四档）。估计器的两个哨兵都归
    `undefined`：输入带 nan（如 x_shock 未检出）→ nan；`delta_xc_c == 0`
    （估计器返回 inf）→ 差值精确为零，无可比。
    """
    ratio, p = _implied_order(h, vals)
    if not np.isfinite(ratio):
        return VERDICT_UNDEFINED, ratio, p
    verdict = (VERDICT_GROWING if ratio >= 1.0 else
               VERDICT_UNCHECKED if ratio < 0.2 else VERDICT_YES)
    return verdict, ratio, p


def u_num_from_ladder(vals, verdict):
    """u_num = |S_coarse − S_medium| —— 我侧最细一档差，与参照侧 `error_bar`
    的「最细区间差」口径逐字一致（**不用 Richardson 外推**：阶在多半腿上读不出，
    外推会把不可判写成可判）。梯子判定不在渐近族（verdict 不以 "yes" 开头）
    ⇒ None（UNDEFINED：这条腿没有误差棒）。
    """
    if not verdict.startswith("yes"):
        return None
    _xc, c, m = (float(v) for v in vals)
    return abs(m - c)


def u_val(u_num, u_D):
    """ASME V&V 20：u_val² = u_num² + u_D²。任一缺 → None（UNDEFINED）。"""
    if u_num is None or u_D is None or not np.isfinite(u_D):
        return None
    return float(np.hypot(u_num, u_D))


def band_for(calibrated_band, u_val_):
    """**min 语义：带永不变宽**；u_val 无定义 → 标定带。"""
    if u_val_ is None:
        return calibrated_band
    return min(calibrated_band, u_val_)


def regime(sd, calibrated_band, u_val_):
    """`grid_convergence_ours.csv` 的 `regime` 列（预注册 §2.3 的四词词汇）。

    sd = |S − D|（与带同单位）。u_val 或 sd 缺 → `no-bar`（没有误差棒的锁）；
    u_val > 标定带 → `floor`（组合不确定度比标定还宽 ⇒ 标定带仍管事，
    但要说出来）；带收紧到 u_val：sd ≤ u_val → `in-uncertainty`，
    sd > u_val → `resolved-model-error {n}×`（**判定，不是红** —— 模型误差
    被分辨出来是 V&V 的正经产出）。
    """
    if u_val_ is None or sd is None:
        return "no-bar"
    if calibrated_band is not None and u_val_ > calibrated_band:
        return "floor"
    if sd <= u_val_:
        return "in-uncertainty"
    return f"resolved-model-error {sd / u_val_:.1f}×"


def leg_uncertainty(vals3, h, u_D, calibrated_band, sd=None):
    """一条腿的全部不确定度量 → dict(verdict/u_num/u_D/u_val/band/regime)。

    门文件只断言，仪器全在这里；**全部输入可注入**（G-TEETH 拿扰动副本喂，
    零求解）。`band` = min(标定带, u_val)（永不变宽）；`regime` 的 floor
    判定用**标定带本身**。
    """
    verdict, _ratio, _p = ladder_readout(vals3, h)
    u_n = u_num_from_ladder(vals3, verdict)
    uv = u_val(u_n, u_D)
    return dict(verdict=verdict, u_num=u_n, u_D=u_D, u_val=uv,
                band=(None if calibrated_band is None
                      else band_for(calibrated_band, uv)),
                regime=regime(sd, calibrated_band, uv))


def write_grid_convergence_ours(path, runs, cases, levels, mesh_dir,
                                quantities=QUANTITIES, u_d_of=None,
                                sd_of=None, band_of=None):
    """每 (case, quantity) 一行：三级读数、相邻差、比值、implied_order、三分判定，
    以及（W2.1b）`u_num / u_D / u_val / regime` 四列。

    `levels` 必须就是 (xcoarse, coarse, medium) —— 列名写死这三级。
    任何一腿不收敛 ⇒ `basis=not_converged` 且 ratio/order/verdict 留空
    （不收敛腿的阶无意义 —— 参照侧对读不出的量同样不给判定）。
    `u_d_of` / `sd_of` / `band_of` 是 (case, quantity) → float|None 的注入钩子：
    u_D = 参照侧散布，sd = 该腿判据的 |S − D|，band = 标定带（**原生单位**）；
    缺省或返回 None ⇒ 对映列留空、regime 退 `no-bar`。
    """
    assert tuple(levels) == ("xcoarse", "coarse", "medium"), levels

    def f9(v):
        return f"{v:.9e}" if np.isfinite(v) else "nan"

    def f9o(v):  # None → 留空（没测到/无定义），与 nan（测了但没有）分开
        return "" if v is None else f9(v)

    n = mesh_node_counts(mesh_dir, levels)
    h = np.array([n[lv] for lv in levels]) ** -0.5

    with open(path, "w", newline="") as fh:
        fh.write("# 我方腿的三级网格收敛读数（W2.1a RECORDED 表；W2.1b 加 "
                 "u_num/u_D/u_val/regime 四列，仍不进值锁）。\n")
        fh.write("# ★ 三档只能给一个数，看不到 unstable：参照侧 131 个量里 41 个的\n")
        fh.write("#   判定翻转是加第四档才看见的 —— 本表的 yes 都带 (3-rung) 限定。\n")
        fh.write("# u_num = |S_coarse − S_medium|（最细区间差，与参照侧 error_bar "
                 "同口径）；u_D = 参照侧散布；\n")
        fh.write("# u_val = √(u_num² + u_D²)；regime 词汇见 tests/D/_grid_convergence.py。\n")
        fh.write("# h ∝ N^(−1/2)，N = 三级 n_nodes（" + " / ".join(
            f"{lv} {int(n[lv])}" for lv in levels) + "）。\n")
        w = csv.writer(fh)
        w.writerow(["case", "quantity", "xcoarse", "coarse", "medium",
                    "delta_xc_c", "delta_c_m", "ratio", "implied_order",
                    "verdict", "basis", "u_num", "u_D", "u_val", "regime"])
        for case in cases:
            for q in quantities:
                vals = [float(runs[(case, lv)][q]) for lv in levels]
                conv = all(bool(runs[(case, lv)]["converged"]) for lv in levels)
                basis = "xc-c-m" if conv else "not_converged"
                xc, c, m = vals
                d_xc_c, d_c_m = c - xc, m - c
                if not conv:
                    ratio_s = order_s = verdict = ""
                    u_n = None
                else:
                    verdict, ratio, p = ladder_readout(vals, h)
                    u_n = u_num_from_ladder(vals, verdict)
                    if not np.isfinite(ratio):
                        ratio_s = ("nan" if ratio != ratio else
                                   "UNDEFINED (coarse delta = 0)")
                        order_s = ""
                    else:
                        ratio_s = f"{ratio:.3f}"
                        order_s = "" if p != p else f"{p:.2f}"
                ud = u_d_of(case, q) if u_d_of else None
                uv = u_val(u_n, ud)
                sd = sd_of(case, q) if sd_of else None
                bd = band_of(case, q) if band_of else None
                reg = regime(sd, bd, uv)
                w.writerow([case, q, f9(xc), f9(c), f9(m),
                            f9(d_xc_c) if not np.isfinite(d_xc_c)
                            else f"{d_xc_c:+.9e}",
                            f9(d_c_m) if not np.isfinite(d_c_m)
                            else f"{d_c_m:+.9e}",
                            ratio_s, order_s, verdict, basis,
                            f9o(u_n), f9o(ud), f9o(uv), reg])
