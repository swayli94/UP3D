r"""读 CFL3D 参考数据集的 **Cp 分布** —— D05..D10 共用。

★★★ **为什么图必须是 Cp 分布，而不是力系数柱状图**（使用者裁决 2026-09-05）：
一个 cl 是**一个数**，两条曲线差 2 % 可能来自激波位置、前缘峰、尾缘卸载中的
任何一处，甚至来自互相抵消的两处。**Cp 分布显示差在哪里**，才是验证计算是否
正确的东西。柱状图保留，但只作为 Cp 图的附注。

★ 一份读取器而不是四份拷贝：四个数据集的 Cp schema 有两种
（2-D 每例一文件 `x_c,y_c,cp,mach_local,surface`；3-D 一个长表带
`level,turb,eta_requested,...`），都在这里处理。
"""
import csv
import os

import numpy as np


def read_2d_cp(ref_dir, case, level="L3", turb="none"):
    """2-D 数据集：``cp_<case>_<turb>_<level>.csv`` -> {surface: (x, cp)}。

    ★ 文件名里的 `turb` 段对 Euler 是字面量 ``none``，对 RANS 是 ``sst``/``sa``
    —— 那是数据集自己的命名，不是我们的约定，所以由调用方传进来。
    """
    p = os.path.join(ref_dir, f"cp_{case}_{turb}_{level}.csv")
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"{p}: the reference Cp for {case}/{turb}/{level} is missing -- "
            f"check the dataset's own file naming before assuming a layout")
    by = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            by.setdefault(r["surface"], []).append(
                (float(r["x_c"]), float(r["cp"])))
    out = {}
    for s, v in by.items():
        v.sort()
        out[s] = (np.array([a for a, _ in v]), np.array([b for _, b in v]))
    return out


def read_3d_cp(ref_dir, level, eta, turb="none"):
    """3-D 数据集：`cp_stations.csv` 的长表 -> {surface: (x, cp)}。"""
    p = os.path.join(ref_dir, "cp_stations.csv")
    by = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            if (r["level"] != level or r.get("turb", "none") != turb
                    or abs(float(r["eta_requested"]) - eta) > 1e-9):
                continue
            by.setdefault(r["surface"], []).append(
                (float(r["x_c"]), float(r["cp"])))
    if not by:
        raise LookupError(
            f"{p}: no rows for level={level!r} turb={turb!r} eta={eta} -- "
            f"the ladder may not have reached that rung")
    out = {}
    for s, v in by.items():
        v.sort()
        out[s] = (np.array([a for a, _ in v]), np.array([b for _, b in v]))
    return out


def cp_rms(x_ref, cp_ref, x_us, cp_us):
    """把**我们的**曲线插到**参考的** x 上再取 RMS。

    ★ 方向是刻意的：参考点是被比较的基准，插值它会把它平滑掉。
    只在两者 x 范围的交集上比较 —— 范围外的外插值不是测量。
    """
    m = (x_ref >= np.min(x_us)) & (x_ref <= np.max(x_us))
    if m.sum() < 3:
        return float("nan"), 0
    d = np.interp(x_ref[m], x_us, cp_us) - cp_ref[m]
    return float(np.sqrt(np.mean(d ** 2))), int(m.sum())


def band_from_two(a, b):
    """两个湍流模型的 Cp -> (x, lo, hi)，即一条**带**。

    ★★ D08/D10 的参考是两个模型，所以 Cp 图上正确的画法是**阴影带**，
    不是两条曲线里挑一条。带宽就是这道门在 Cp 上的分辨底噪。
    """
    xa, ca = a
    xb, cb = b
    x = np.union1d(xa, xb)
    lo_a, lo_b = np.interp(x, xa, ca), np.interp(x, xb, cb)
    return x, np.minimum(lo_a, lo_b), np.maximum(lo_a, lo_b)
