r"""hex 族 NACA0012 M0.80/α1.25 的 cl −24 % 归因（预注册 Q2：
phases/p6/docs/dev_phase_six/20261012-1000-research-amg-hex-prereg.md）。

实测背景（W2.6，20261012-1000-w25-w26-verdict.md §4.2）：hex medium 全收敛下
cl 0.2582 vs tet medium committed 0.341340（D05），x_shock 一致（差 0.010 c），
几何无辜（hex 壁面对精确剖面 max 差 2.15e-5）。候选：① 分辨率分布；
② 系统性缺陷；③ tet 侧参考值本身有偏。

判别方法（预注册写死）：
  (a) 两族 medium 的截面 Cp 逐面比较（升力差集中在哪一段）；
  (b) hex 加细到 fine，看差是否收缩 —— 收缩 ⇒ ①，不收缩 ⇒ ②；
  (c) hex coarse→medium→fine 三级形态 vs tet（W2.1a 的读数习惯：
      参照侧 implied_order 仪器 + yes/no/unstable 三分词汇，
      h ∝ N^(−1/2)，2-D 用 −1/2）。

生产配方与 D05 逐字一致（tests/D/test_D05_euler_naca0012.py 的 RECIPE）。
kill criterion：hex fine 不收敛 ⇒ 停，记录不收敛本身，不换网格。

用法：
  python run.py --legs hex:coarse            # 单腿（冒烟）
  python run.py                              # 全部 4 腿 + 分析 + 图
环境：conda up3d；NUMBA_NUM_THREADS=16；BLAS 由 pyfp3d 自身钉 1
（脚本断言 blas_thread_count()==1，把出处写进 CSV）。
"""
import argparse
import csv
import importlib.util
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))

#: 生产配方，与 `tests/D/test_D05_euler_naca0012.py` 的 RECIPE 逐字一致
RECIPE = dict(upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=400,
              n_picard_seed=5)
M_INF, ALPHA = 0.80, 1.25

TET_DIR = os.path.join(REPO, "cases", "meshes", "naca0012_2.5d")
HEX_DIR = os.path.join(REPO, "cases", "meshes", "naca0012_hex_2.5d")
D05_DIR = os.path.join(REPO, "cases", "gates", "D05_euler_naca0012")
CFL3D_REF = os.path.join(REPO, "cases", "reference_data", "cfl3d",
                         "euler_naca0012")

#: 新算的腿（hex 三级 + tet medium 一条，后者为 Cp 曲线 —— committed 证据里
#: 只有 cl/cd/x_shock 数值与 PNG，没有 Cp 数据列；其 cl 必须与 committed
#: 0.341340 五位相合，作锚）。
FRESH_LEGS = (("hex", "coarse"), ("hex", "medium"), ("hex", "fine"),
              ("tet", "medium"))

#: ★ Cp 分段（事前写死，看数之前）：LE 峰 [0,0.05]、前中段 [0.05,0.45]、
#: 激波区 [0.45,0.75]（两族 x_shock 都在 ~0.60）、尾缘卸载 [0.75,1.0]。
SEGMENTS = (("LE_peak", 0.00, 0.05), ("front", 0.05, 0.45),
            ("shock", 0.45, 0.75), ("TE_unload", 0.75, 1.00))


def _f9(v):
    return f"{float(v):.9e}"


def mesh_stats_wide(mesh_dir, level):
    """hex 族 stats 是宽表（一行）；tet 族是长表（metric,value）。各读各的。"""
    path = os.path.join(mesh_dir, f"{level}_stats.csv")
    with open(path) as fh:
        rows = list(csv.reader(fh))
    if rows and rows[0] and rows[0][0] == "metric":
        return {r[0]: r[1] for r in rows[1:] if r}
    return dict(zip(rows[0], rows[1]))


def run_leg(family, level):
    """与 D05 `_one()` 同构的一腿：求解 + 力系数 + 截面 Cp + 激波报告。"""
    from pyfp3d import blas_thread_count
    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.mesh.wake_cut import cut_wake
    from pyfp3d.post.section_cut import section_cp_curve
    from pyfp3d.post.shock import shock_report
    from pyfp3d.post.surface import wall_force_coefficients
    from pyfp3d.solve.newton import solve_newton_lifting

    mesh_dir = TET_DIR if family == "tet" else HEX_DIR
    st = mesh_stats_wide(mesh_dir, level)
    mc, wc = cut_wake(read_mesh(os.path.join(mesh_dir, f"{level}.msh")))
    dz = float(np.ptp(mc.nodes[:, 2]))
    t0 = time.time()
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **RECIPE)
    wall = time.time() - t0
    phi = np.asarray(r["phi"])
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], phi,
                                alpha_deg=ALPHA, u_inf=1.0, s_ref=dz,
                                m_inf=M_INF)
    cur = section_cp_curve(mc, phi, z=float(np.mean(mc.nodes[:, 2])),
                           smooth_passes=1, m_inf=M_INF)
    up = shock_report(cur, M_INF)["upper"]
    frz = r.get("sigma_freeze_report") or {}
    hist = np.asarray(r["residual_history"], float)
    row = dict(family=family, level=level, source="fresh",
               n_nodes=int(float(st["n_nodes"])),
               n_tets=int(float(st.get("n_tets") or st.get("n_elements"))),
               dz=dz,
               cl=float(f["cl"]), cd=float(f["cd_pressure"]),
               x_shock=(float(up["x_shock"]) if up.get("has_shock")
                        else float("nan")),
               converged=bool(r.get("converged")),
               residual=float(hist[-1]), n_newton=int(len(hist)),
               n_limited=int(r.get("n_limited", 0)),
               n_floored=int(r.get("n_floored", 0)),
               sigma_frozen=bool(frz.get("frozen_in_transient", False)),
               sigma_churn=bool(frz.get("selection_churn", False)),
               wall_time_s=wall,
               n_threads=int(os.environ.get("NUMBA_NUM_THREADS", "0") or 0),
               blas_threads=blas_thread_count())
    # 残差历史全量落盘（kill criterion 记录：fine 不收敛的形状在此）
    with open(os.path.join(HERE, f"residual_history_{family}_{level}.csv"),
              "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "residual"])
        for i, res in enumerate(hist):
            w.writerow([i, _f9(res)])
    return row, cur


def tet_committed_ladder():
    """tet 三级 committed 值（D05 summary.csv，证据纪律：不重算已提交 artifact）。"""
    out = {}
    with open(os.path.join(D05_DIR, "summary.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["case"] == "M0.80_a1.25":
                out[r["level"]] = dict(
                    family="tet", level=r["level"], source="committed",
                    cl=float(r["cl"]), cd=float(r["cd"]),
                    x_shock=float(r["x_shock"]),
                    converged=r["converged"] == "1",
                    residual=float(r["residual"]),
                    n_threads=int(r["n_threads"]),
                    blas_threads=int(r["blas_threads"]))
    return out


def _sorted_unique(x, y):
    o = np.argsort(x)
    xs, ys = np.asarray(x)[o], np.asarray(y)[o]
    keep = np.concatenate(([True], np.diff(xs) > 0))
    return xs[keep], ys[keep]


def segment_cl(cur, x0, x1, n=4001):
    """一段对 cl 的贡献：∫(cp_lower − cp_upper) dx（2-D 截面，弦长 1）。"""
    xg = np.linspace(x0, x1, n)
    tot = 0.0
    for side, sgn in (("lower", +1.0), ("upper", -1.0)):
        xs, ys = _sorted_unique(cur[f"x_{side}"], cur[f"cp_{side}"])
        tot += sgn * float(np.trapezoid(np.interp(xg, xs, ys), xg))
    return tot


def total_cl_from_cp(cur):
    return segment_cl(cur, 0.0, 1.0)


def load_grid_convergence_instrument():
    """W2.1a 的仪器（tests/D/_grid_convergence.py），按 bench 惯例
    spec_from_file_location 加载，不走包 import。"""
    spec = importlib.util.spec_from_file_location(
        "_grid_convergence",
        os.path.join(REPO, "tests", "D", "_grid_convergence.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--legs", default=",".join(f"{f}:{l}" for f, l in FRESH_LEGS),
                    help="逗号分隔的 family:level 列表")
    ap.add_argument("--analysis-only", action="store_true",
                    help="不求解，从 ladder.csv/cp_curves.csv 重出分析")
    args = ap.parse_args()
    wanted = [tuple(s.split(":")) for s in args.legs.split(",") if s]

    os.makedirs(HERE, exist_ok=True)
    ladder_path = os.path.join(HERE, "ladder.csv")
    cp_path = os.path.join(HERE, "cp_curves.csv")

    from pyfp3d import blas_thread_count
    assert blas_thread_count() == 1, (
        f"BLAS 未钉 1（实测 {blas_thread_count()}）——证据纪律要求 BLAS=1")

    rows, curves = [], {}
    if not args.analysis_only:
        for family, level in wanted:
            print(f"[run] {family}:{level} ...", flush=True)
            row, cur = run_leg(family, level)
            rows.append(row)
            curves[(family, level)] = cur
            print(f"[done] {family}:{level} cl={row['cl']:.6f} "
                  f"cd={row['cd']:.6f} x_shock={row['x_shock']:.4f} "
                  f"converged={row['converged']} |R|={row['residual']:.2e} "
                  f"n_newton={row['n_newton']} wall={row['wall_time_s']:.0f}s",
                  flush=True)
            # ★ 每腿落地一次：中途中断也保留已算证据
            _append_ladder(ladder_path, row)
            _append_cp(cp_path, family, level, cur)
    else:
        rows = list(csv.DictReader(open(ladder_path)))
        curves = _read_cp(cp_path)

    if len(wanted) == len(FRESH_LEGS) or args.analysis_only:
        analyze(rows, curves)


def _append_ladder(path, row):
    cols = ["family", "level", "source", "n_nodes", "n_tets", "dz",
            "cl", "cd", "x_shock", "converged", "residual", "n_newton",
            "n_limited", "n_floored", "sigma_frozen", "sigma_churn",
            "wall_time_s", "n_threads", "blas_threads"]
    new = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(cols)
        w.writerow([row[c] if c not in ("cl", "cd", "x_shock", "residual",
                                        "dz", "wall_time_s")
                    else _f9(row[c]) for c in cols])


def _append_cp(path, family, level, cur):
    new = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["family", "level", "side", "x", "cp"])
        for side in ("upper", "lower"):
            xs, ys = _sorted_unique(cur[f"x_{side}"], cur[f"cp_{side}"])
            for x, y in zip(xs, ys):
                w.writerow([family, level, side, _f9(x), _f9(y)])


def _read_cp(path):
    curves = {}
    for r in csv.DictReader(open(path)):
        key = (r["family"], r["level"])
        curves.setdefault(key, {})
        curves[key].setdefault(f"x_{r['side']}", []).append(float(r["x"]))
        curves[key].setdefault(f"cp_{r['side']}", []).append(float(r["cp"]))
    return {k: {kk: np.asarray(vv) for kk, vv in v.items()}
            for k, v in curves.items()}


def analyze(rows, curves):
    gc = load_grid_convergence_instrument()
    tet_com = tet_committed_ladder()

    def row_of(family, level):
        for r in rows:
            if r["family"] == family and r["level"] == level:
                return r
        return None

    def fget(r, k):
        return r[k] if isinstance(r[k], float) else float(r[k])

    # —— 三级形态（W2.1a 读数习惯：h ∝ N^(−1/2)，三分词汇）——
    hex_levels = ("coarse", "medium", "fine")
    hex_n = {lv: float(mesh_stats_wide(HEX_DIR, lv)["n_nodes"])
             for lv in hex_levels}
    hex_h = np.array([hex_n[lv] for lv in hex_levels]) ** -0.5
    tet_levels = ("xcoarse", "coarse", "medium")
    tet_n = {lv: float(mesh_stats_wide(TET_DIR, lv)["n_nodes"])
             for lv in tet_levels}
    tet_h = np.array([tet_n[lv] for lv in tet_levels]) ** -0.5

    hex_vals = {q: [fget(row_of("hex", lv), q) for lv in hex_levels]
                for q in ("cl", "cd", "x_shock")}
    tet_vals = {q: [tet_com[lv][q] for lv in tet_levels]
                for q in ("cl", "cd", "x_shock")}

    with open(os.path.join(HERE, "grid_convergence_hex.csv"), "w",
              newline="") as fh:
        fh.write("# hex 族 coarse/medium/fine 三级读数（对照 tet 族 committed "
                 "xcoarse/coarse/medium，后者抄自 D05 grid_convergence_ours.csv）。\n")
        fh.write("# 仪器 = tests/D/_grid_convergence.py 的 ladder_readout "
                 "（参照侧 implied_order + 三分词汇）；h ∝ N^(−1/2)。\n")
        fh.write("# hex N: " + " / ".join(f"{lv} {int(hex_n[lv])}"
                                          for lv in hex_levels) + "；")
        fh.write("tet N: " + " / ".join(f"{lv} {int(tet_n[lv])}"
                                        for lv in tet_levels) + "\n")
        fh.write("# ★ 收敛口径：strict = 门的 converged 旗标（fine 的 cl/cd "
                 "被判 not_converged）；physical = |R| ≤ 1e-6 的物理解"
                 "（fine |R|=2.25e-7，进物理族——收敛旗标本身是线程依赖噪声，"
                 "见 D05 docstring）。\n")
        w = csv.writer(fh)
        w.writerow(["family", "quantity", "level1", "level2", "level3",
                    "delta_12", "delta_23", "ratio", "implied_order",
                    "verdict", "basis", "u_num_l23"])

        def is_conv(r):
            c = r["converged"]
            return bool(c) if isinstance(c, bool) else c in ("True", "1")

        # hex 三级的 strict（门旗标）与 physical（|R| ≤ 1e-6）收敛状态
        hex_conv_strict = all(is_conv(row_of("hex", lv)) for lv in hex_levels)
        hex_conv_phys = all(
            is_conv(row_of("hex", lv))
            or float(row_of("hex", lv)["residual"]) <= 1e-6
            for lv in hex_levels)
        for fam, vals3, h3 in (("hex", hex_vals, hex_h),
                               ("tet", tet_vals, tet_h)):
            for q, vv in vals3.items():
                # 收敛口径（参照 D05 write_grid_convergence_ours：任一腿不收敛
                # ⇒ basis=not_converged，阶/判定留空）。tet committed 三级全收敛。
                # hex fine 的 converged 旗标 False（|R| 2.25e-7 撞 400 步帽）：
                # cl/cd 按 strict 判 not_converged；x_shock 单列 physical 口径
                # （旗标本身是线程依赖噪声，D05 docstring）。
                if fam == "tet" or hex_conv_strict:
                    verdict, ratio, p = gc.ladder_readout(vv, h3)
                    u_n = gc.u_num_from_ladder(vv, verdict)
                    basis = "xc-c-m" if fam == "tet" else "c-m-f"
                elif q == "x_shock" and hex_conv_phys:
                    verdict, ratio, p = gc.ladder_readout(vv, h3)
                    u_n = gc.u_num_from_ladder(vv, verdict)
                    basis = "c-m-f (physical)"
                else:
                    verdict, ratio, p = "", float("nan"), float("nan")
                    u_n = None
                    basis = "not_converged"
                w.writerow([fam, q, *[_f9(v) for v in vv],
                            f"{vv[1]-vv[0]:+.9e}", f"{vv[2]-vv[1]:+.9e}",
                            (f"{ratio:.3f}" if np.isfinite(ratio) else ""),
                            ("" if p != p else f"{p:.2f}"), verdict,
                            basis, ("" if u_n is None else _f9(u_n))])

    # —— Cp 分段归因（(a)，段边界事前写死）——
    seg_rows = []
    for (fam, lv), cur in curves.items():
        for name, x0, x1 in SEGMENTS:
            seg_rows.append(dict(family=fam, level=lv, segment=name,
                                 x0=x0, x1=x1, cl_seg=segment_cl(cur, x0, x1)))
        seg_rows.append(dict(family=fam, level=lv, segment="TOTAL",
                             x0=0.0, x1=1.0, cl_seg=total_cl_from_cp(cur)))
    with open(os.path.join(HERE, "cp_segment_cl.csv"), "w", newline="") as fh:
        fh.write("# 分段 cl 贡献 ∫(cp_lower − cp_upper) dx；段边界事前写死："
                 + ", ".join(f"{n} [{a},{b}]" for n, a, b in SEGMENTS) + "\n")
        w = csv.writer(fh)
        w.writerow(["family", "level", "segment", "x0", "x1", "cl_seg"])
        for r in seg_rows:
            w.writerow([r["family"], r["level"], r["segment"], r["x0"],
                        r["x1"], _f9(r["cl_seg"])])
    # 差表：hex medium/fine − tet medium，逐段
    with open(os.path.join(HERE, "cp_segment_delta.csv"), "w", newline="") as fh:
        fh.write("# 逐段 Δcl_seg = hex − tet_medium（两族 medium 对表 + "
                 "fine 对表；负 = hex 在该段少升力）\n")
        w = csv.writer(fh)
        w.writerow(["segment", "tet_medium", "hex_medium", "hex_fine",
                    "d_hex_medium", "d_hex_fine"])
        for name, _a, _b in (*SEGMENTS, ("TOTAL", 0, 1)):
            g = {(r["family"], r["level"]): r["cl_seg"] for r in seg_rows
                 if r["segment"] == name}
            tm = g.get(("tet", "medium"), float("nan"))
            hm = g.get(("hex", "medium"), float("nan"))
            hf = g.get(("hex", "fine"), float("nan"))
            w.writerow([name, _f9(tm), _f9(hm), _f9(hf),
                        _f9(hm - tm), _f9(hf - tm)])

    # —— 气泡计：上翼面 x∈[0.2,0.55] 的 max Cp（正 = 中部减速气泡存在）——
    with open(os.path.join(HERE, "bubble_metrics.csv"), "w", newline="") as fh:
        fh.write("# 上翼面 [0.2,0.55] 内 max Cp 及其位置：正 cp = 局部已减速到"
                 "亚声速以下（气泡/双激波结构）；hex 族该结构随加细增强并前移。\n")
        w = csv.writer(fh)
        w.writerow(["family", "level", "max_cp_upper_20_55", "x_at_max"])
        for (fam, lv), cur in sorted(curves.items()):
            xs, ys = _sorted_unique(cur["x_upper"], cur["cp_upper"])
            m = (xs >= 0.2) & (xs <= 0.55)
            i = int(np.argmax(ys[m]))
            w.writerow([fam, lv, _f9(ys[m][i]), _f9(xs[m][i])])

    _figures(curves, tet_com, rows)


def _figures(curves, tet_com, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    spec = importlib.util.spec_from_file_location(
        "_cfl3d_cp", os.path.join(REPO, "tests", "D", "_cfl3d_cp.py"))
    c3 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(c3)

    def fget(r, k):
        return r[k] if isinstance(r[k], float) else float(r[k])

    def one(family, level):
        return next(r for r in rows
                    if r["family"] == family and r["level"] == level)

    rc = c3.read_2d_cp(CFL3D_REF, "n0012_m0800_a1.25", "L4", "none")

    # 图 1：(a) 两族 medium 的 Cp 对照 + 分段着色
    fig, ax = plt.subplots(figsize=(10, 6.2))
    for side, ls in (("upper", "-"), ("lower", "--")):
        if side in rc:
            ax.plot(rc[side][0], rc[side][1], ls, color="#c0392b", lw=1.5,
                    label="CFL3D Euler L4" if side == "upper" else None)
    for (fam, lv), sty, lab in (
            (("tet", "medium"), dict(lw=1.5, color="#1a56db"), "tet medium"),
            (("hex", "medium"), dict(lw=1.5, color="#e08a00"), "hex medium"),
            (("hex", "fine"), dict(lw=1.2, ls="-.", color="#2a9d3f"),
             "hex fine")):
        cur = curves.get((fam, lv))
        if cur is None:
            continue
        for side in ("upper", "lower"):
            xs, ys = _sorted_unique(cur[f"x_{side}"], cur[f"cp_{side}"])
            ax.plot(xs, ys, label=lab if side == "upper" else None, **sty)
    for name, x0, x1 in SEGMENTS:
        ax.axvspan(x0, x1, alpha=0.05)
        ax.text((x0 + x1) / 2, 1.02, name, ha="center",
                fontsize=8, color="#666")
    ax.invert_yaxis()
    ax.grid(alpha=.3)
    ax.set_xlabel("x/c")
    ax.set_ylabel("$C_p$")
    cl_tm = fget(one("tet", "medium"), "cl")
    cl_hm = fget(one("hex", "medium"), "cl")
    cl_hf = fget(one("hex", "fine"), "cl")
    ax.set_title(f"Q2(a): Cp, M0.80/α1.25 — cl: tet med {cl_tm:.4f} / "
                 f"hex med {cl_hm:.4f} / hex fine {cl_hf:.4f} "
                 f"(CFL3D L4 0.3489)", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "cp_tet_vs_hex_medium.png"), dpi=140)
    plt.close(fig)

    # 图 2：(c) hex 三级形态 vs tet medium
    fig, ax = plt.subplots(figsize=(10, 6.2))
    for side, ls in (("upper", "-"), ("lower", "--")):
        if side in rc:
            ax.plot(rc[side][0], rc[side][1], ls, color="#c0392b", lw=1.4,
                    label="CFL3D Euler L4" if side == "upper" else None)
    styles = {("hex", "coarse"): dict(lw=0.9, ls=":", color="#9aa0a6"),
              ("hex", "medium"): dict(lw=1.3, ls="--", color="#e08a00"),
              ("hex", "fine"): dict(lw=1.5, color="#2a9d3f"),
              ("tet", "medium"): dict(lw=1.3, ls="--", color="#1a56db")}
    for key, sty in styles.items():
        cur = curves.get(key)
        if cur is None:
            continue
        for side in ("upper", "lower"):
            xs, ys = _sorted_unique(cur[f"x_{side}"], cur[f"cp_{side}"])
            ax.plot(xs, ys, label=" ".join(key) if side == "upper" else None,
                    **sty)
    ax.invert_yaxis()
    ax.grid(alpha=.3)
    ax.set_xlabel("x/c")
    ax.set_ylabel("$C_p$")
    ax.set_title("Q2(c): hex coarse→medium→fine morphology vs tet medium, "
                 "M0.80/α1.25", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "cp_hex_ladder_vs_tet.png"), dpi=140)
    plt.close(fig)

    # 图 3：cl 梯子（hex 三级 fresh + tet 三级 committed + CFL3D L4 线）
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    hex_lv = ("coarse", "medium", "fine")
    tet_lv = ("xcoarse", "coarse", "medium")
    hn = [float(mesh_stats_wide(HEX_DIR, lv)["n_nodes"]) for lv in hex_lv]
    tn = [float(mesh_stats_wide(TET_DIR, lv)["n_nodes"]) for lv in tet_lv]
    hcl = [fget(one("hex", lv), "cl") for lv in hex_lv]
    tcl = [tet_com[lv]["cl"] for lv in tet_lv]
    ax.plot(np.array(hn) ** -0.5, hcl, "o-", color="#2a9d3f",
            label="hex (fresh)")
    ax.plot(np.array(tn) ** -0.5, tcl, "s--", color="#1a56db",
            label="tet (committed D05)")
    ax.axhline(0.348881, color="#c0392b", lw=1.2, ls=":",
               label="CFL3D Euler L4 0.3489")
    for h, v, lv in zip(np.array(hn) ** -0.5, hcl, hex_lv):
        ax.annotate(f"hex {lv}\n{v:.4f}", (h, v), fontsize=7,
                    textcoords="offset points", xytext=(6, -12))
    for h, v, lv in zip(np.array(tn) ** -0.5, tcl, tet_lv):
        ax.annotate(f"tet {lv}\n{v:.4f}", (h, v), fontsize=7,
                    textcoords="offset points", xytext=(6, 6))
    ax.set_xlabel("$h \\propto N^{-1/2}$")
    ax.set_ylabel("$c_l$")
    ax.grid(alpha=.3)
    ax.set_title("Q2(b)+(c): cl ladder, M0.80/α1.25", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "cl_ladder.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
