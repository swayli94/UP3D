"""C10 — hex 结构族 vs tet 族的**圆柱交叉门**（C 类；预注册
docs/dev_phase_six/20261012-1300-c10-hex-cylinder-prereg.md）。

★★★ **攻击点**：若 hex 生成器（`pyfp3d/meshgen/structured.py`）有一个**系统
缺陷**，所有用它的门共享该缺陷。同几何（半径 1 的圆柱）两拓扑对照：

- **非升力腿**（Γ=0）：hex medium 与已提交 tet（`cylinder_2.5d/medium.msh`，
  同几何同 dz）各解 Laplace 极限（`solve_newton_nonlifting`，m_inf=0），
  壁面 Cp 对**精确解**（仪器是 C02 的 `tests.mesh_utils.cylinder_grad_exact`）；
- **升力腿**（Γ=0.6）：hex（2026-10-12 补的 wake 组 —— O-grid 的接缝就是
  wake，圆柱没有尖尾缘但环量通道必须有分支切口）vs C06 的内存 tet 族，
  各对**自己的**精确 cl = 2Γ/(U·c)。★ 口径：两族半径不同（hex a=1、c=2；
  C06 族 a=0.5、c=1），所以各比各的相对误差、共用一条带 —— 写明免得有人
  拿绝对 cl 直接对。

★ **带 = max(实测) + 两族实现噪声合成**（预注册的失败标准：带不许写在
实测差之外而不过噪声地板）。实测表（2026-10-12，r_far 惰性旋钮变体，
生产求解路径）：

非升力（Γ=0，m_inf=0 的 Laplace 极限）：

| 实现 | cp_rms | cp_max | |cl| |
|---|---|---|---|
| tet 已提交 medium | 0.0321 | 0.0590 | 2.4e-5 |
| tet r_far=15 / 25 | 0.0390 / 0.0293 | 0.0704 / 0.0518 | ≤ 6e-5 |
| hex 已提交 medium | 0.0365 | 0.0630 | 5.0e-5 |
| hex r_far=15 / 25 | 0.0411 / 0.0343 | 0.0717 / 0.0581 | ≈ 5e-5 |

升力（Γ=0.6，cl 对各自精确解的相对误差 %）：

| 实现 | cl_rel |
|---|---|
| hex medium / r15 / r25 | 0.42 / 0.58 / 0.33 % |
| tet（C06 族 medium，W2.5 实测跨实现） | [0.95, 1.14] % |

★ 注意读法：两族的**误差水平**本来就不同（hex 壁面法向更细），所以
「交叉」断言的不是两族误差相等，而是各自的误差**今天落在哪里** ——
一族出系统缺陷时它自己会炸出带。
"""
import os

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import (wall_force_coefficients,
                                 wall_tangential_gradient_quadratic)
from pyfp3d.solve.newton import solve_newton_nonlifting
from pyfp3d.solve.picard import solve_laplace_lifting
from tests._gate_evidence import (assert_matches_committed, blas_thread_count,
                                  fmt)
from tests.conftest import REPO_ROOT, gate_figures_enabled
from tests.mesh_utils import cylinder_grad_exact

#: Γ=0.6 与 C06 同；hex 族半径 1 ⇒ 弦长 2 ⇒ cl 精确值 0.6；C06 族 cl 精确值 1.2。
GAMMA = 0.6
CL_EXACT_HEX = 2.0 * GAMMA / 2.0
CL_EXACT_TET = 2.0 * GAMMA / 1.0

# —— 判据（标定；推导见模块 docstring 的实测表与下面每行的注释）——
NL_CP_RMS_MAX = 0.06     # max 实测 0.0411 + tet 极差 0.0097 + hex 极差 0.0068 ≈ 0.058
NL_CP_MAX_MAX = 0.11     # max 实测 0.0717 + 极差合成 0.0186 + 0.0136 ≈ 0.104
NL_CL_MAX = 2.0e-3       # Γ=0 的 |cl|：实测 ≤ 6e-5；带按「wake/几何类缺陷」量级写
NL_CD_MAX = 2.0e-3       # 同上（d'Alembert 对照；实测 ≤ 1.7e-4）
NL_CROSS_CP_RMS = 0.017  # 两族 cp_rms 差：实测 0.0044 + 极差合成 0.0097 + 0.0068
LIFT_CL_REL_MAX = 0.016  # max 实测 1.14 % + tet 极差 0.19pp + hex 极差 0.24pp ≈ 1.57 %
LIFT_CROSS_CL_REL = 0.012  # 两族 cl_rel 差：实测偏置差 0.72pp + 噪声合成 0.43pp ≈ 0.0115

MESH_TET = os.path.join(str(REPO_ROOT), "cases", "meshes", "cylinder_2.5d",
                        "medium.msh")
MESH_HEX = os.path.join(str(REPO_ROOT), "cases", "meshes", "cylinder_hex_2.5d",
                        "medium.msh")


def _cp_error(mesh, phi):
    """壁面 Cp 对精确解（非升力圆柱）：rms / max / 壁面 theta。"""
    wf = mesh.boundary_faces["wall"]
    wall = np.unique(wf)
    grad = wall_tangential_gradient_quadratic(mesh.nodes, wf, phi)
    cp = 1.0 - np.sum(grad[wall] ** 2, axis=1)
    cp_ex = 1.0 - np.sum(cylinder_grad_exact(mesh.nodes[wall]) ** 2, axis=1)
    err = cp - cp_ex
    th = np.arctan2(mesh.nodes[wall, 1], mesh.nodes[wall, 0])
    return dict(cp_rms=float(np.sqrt(np.mean(err ** 2))),
                cp_max=float(np.abs(err).max()),
                theta=th, cp=cp, cp_exact=cp_ex)


def _solve_nonlifting(mesh):
    """Γ=0 的 Laplace 极限（W2.2a 的非升力入口；m_inf=0 逐位等于 Laplace）。"""
    r = solve_newton_nonlifting(mesh, m_inf=0.0, precond="direct")
    phi = np.asarray(r["phi"])
    d = _cp_error(mesh, phi)
    f = wall_force_coefficients(mesh.nodes, mesh.elements,
                                mesh.boundary_faces["wall"], phi, alpha_deg=0.0,
                                u_inf=1.0,
                                #: s_ref 是**面积**（弦×展宽）；半径 1 的圆柱弦 2
                                s_ref=2.0 * float(np.ptp(mesh.nodes[:, 2])),
                                m_inf=0.0)
    d.update(cl=float(f["cl"]), cd=abs(float(f["cd_pressure"])),
             converged=bool(r.get("converged", True)),
             n_nodes=len(mesh.nodes), n_tets=len(mesh.elements))
    return d


def _solve_lifting(mesh, cl_exact):
    """Γ 规定的 Laplace 升力解（C06 的求解路径原样：wake 切开、Kutta 环短路）。"""
    mc, wc = cut_wake(mesh)
    n_station = int(np.unique(wc.te_station).size)
    r = solve_laplace_lifting(mc, wc, alpha_deg=0.0, u_inf=1.0,
                              gamma_fixed=np.full(n_station, GAMMA))
    #: s_ref 是**面积**（弦×展宽）；cl_exact = 2Γ/(U·c) ⇒ c = 2Γ/cl_exact
    #:（hex 族半径 1 ⇒ 弦 2；C06 族弦 1）—— 两族口径不同，各用各的。
    chord = 2.0 * GAMMA / cl_exact
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"],
                                np.asarray(r["phi"]), alpha_deg=0.0, u_inf=1.0,
                                s_ref=chord * float(np.ptp(mc.nodes[:, 2])),
                                m_inf=0.0)
    cl = float(f["cl"])
    return dict(cl=cl, cl_rel=abs(cl - cl_exact) / cl_exact,
                cd=abs(float(f["cd_pressure"])),
                n_kutta=int(r["n_kutta_updates"]),
                residual=float(r["residual_norm"]),
                converged=bool(r.get("converged", True)), n_nodes=len(mesh.nodes),
                n_tets=len(mesh.elements))


def _tet_lifting_reference():
    """C06 的内存 tet 族（半径 0.5、弦 1）—— 升力腿的 tet 参照由它提供：
    已提交的 cylinder_2.5d 族没有 wake 组，承载不了环量。复用 C06 自己的
    `_one_level`（同一网格、同一求解调用）。"""
    from tests.C.test_C06_lifting_cylinder import _one_level
    d = _one_level(0.02)
    return dict(cl=d["cl"], cl_rel=d["cl_rel"], converged=True,
                n_nodes=d["n_nodes"])


@pytest.fixture(scope="module")
def runs():
    out = {}
    out["nl_hex"] = _solve_nonlifting(read_mesh(MESH_HEX))
    out["nl_tet"] = _solve_nonlifting(read_mesh(MESH_TET))
    out["lift_hex"] = _solve_lifting(read_mesh(MESH_HEX), CL_EXACT_HEX)
    out["lift_tet"] = _tet_lifting_reference()
    return out


class TestNonliftingCrosscheck:
    """Γ=0：两族在同一精确解（Cp = 1 − 4 sin²θ）上对照。红了 = 一族的
    生成器有系统缺陷（或壁面恢复/求解路径坏了）—— 不是容差漂移。"""

    def test_both_families_cp_error_inside_the_noise_band(self, runs):
        for fam in ("nl_hex", "nl_tet"):
            d = runs[fam]
            assert d["converged"], f"{fam} 的 Γ=0 解未收敛"
            assert d["cp_rms"] <= NL_CP_RMS_MAX, (
                f"{fam}: Cp rms {d['cp_rms']:.4f} > {NL_CP_RMS_MAX}"
                f"（带 = max 实测 0.0411 + 两族极差合成）—— 若只有 hex 超带，"
                "查 hex 生成器；若 tet 也超，查壁面恢复或求解器")
            assert d["cp_max"] <= NL_CP_MAX_MAX, (
                f"{fam}: Cp max {d['cp_max']:.4f} > {NL_CP_MAX_MAX}"
                "（带推导同上）—— 局部炸点 = 系统缺陷")
            assert abs(d["cl"]) <= NL_CL_MAX, (
                f"{fam}: Γ=0 的 cl {d['cl']:+.2e} 超出 |cl| 带 {NL_CL_MAX}"
                " —— 对称外形出升力 = 组/wake 标记错了")
            assert d["cd"] <= NL_CD_MAX, (
                f"{fam}: Γ=0 的 |cd| {d['cd']:.2e} > {NL_CD_MAX}"
                "（d'Alembert 对照，实测 ≤ 1.7e-4）")

    def test_the_two_families_agree_with_each_other(self, runs):
        """两族 cp_rms 的差必须在两族实现噪声合成之内 —— 这才是「交叉」本身。"""
        d = abs(runs["nl_hex"]["cp_rms"] - runs["nl_tet"]["cp_rms"])
        assert d <= NL_CROSS_CP_RMS, (
            f"hex vs tet 的 cp_rms 差 {d:.4f} > 噪声合成带 {NL_CROSS_CP_RMS}"
            f"（实测 0.0044 = 0.0365 − 0.0321）—— 一族偏出去了")


class TestLiftingCrosscheck:
    """Γ=0.6：两族各对自己的精确 cl。注意口径：hex 弦长 2、C06 族弦长 1，
    各比各的相对误差、共用一条带。"""

    def test_both_families_cl_rel_inside_the_shared_band(self, runs):
        for fam in ("lift_hex", "lift_tet"):
            d = runs[fam]
            assert d["converged"], f"{fam} 的升力解未收敛"
            assert d["cl_rel"] <= LIFT_CL_REL_MAX, (
                f"{fam}: cl_rel {100 * d['cl_rel']:.3f} % > 共用带 "
                f"{100 * LIFT_CL_REL_MAX:.1f} %（带 = max 实测 + 两族极差合成）"
                " —— 一族有系统缺陷（wake 路由 / 组标记 / 几何）")
        if runs["lift_hex"]["n_kutta"] != 0:
            raise AssertionError("gamma_fixed 没有短路 Kutta 环（hex 腿）")

    def test_the_cross_family_bias_difference_is_registered_not_gated(self, runs):
        """两族的 cl_rel 水平本来就不同（hex 壁面法向更细）—— 差值按实测
        登记，带 = 偏置差实测 + 两族噪声合成；超出 = 一族的偏置漂移了。"""
        d = abs(runs["lift_hex"]["cl_rel"] - runs["lift_tet"]["cl_rel"])
        assert d <= LIFT_CROSS_CL_REL, (
            f"两族 cl_rel 差 {100 * d:.3f}pp > 带 {100 * LIFT_CROSS_CL_REL:.2f}pp"
            f"（实测偏置差 0.72pp + 噪声合成 0.43pp）—— 一族的误差水平漂移了")


class TestCommittedEvidenceIsLoadBearing:
    r"""★★★ 新鲜计算 vs 已提交 `summary.csv`（`tests/_gate_evidence.py` 的机器）。
    ★ `blas_threads` 出处列进锁（第 34 轮先例：答案只随库实际 BLAS 线程数动）。"""

    MEASURED = ("cp_rms", "cp_max", "cl", "cd", "cl_rel", "blas_threads")

    def test_matches_committed_summary(self, runs, gate_evidence_dir):
        fresh = {}
        for leg, d in runs.items():
            row = {k: (fmt(d[k]) if k in d and isinstance(d[k], float)
                       else None) for k in ("cp_rms", "cp_max", "cl", "cd",
                                            "cl_rel")}
            row = {k: v for k, v in row.items() if v is not None}
            row["blas_threads"] = fmt(blas_thread_count())
            fresh[leg] = row
        n = assert_matches_committed(
            gate_evidence_dir, fresh, self.MEASURED,
            key_of=lambda r: r["leg"],
            refresh_hint="PYFP3D_GATE_FIGURES=1 pytest "
                         "tests/C/test_C10_hex_cylinder_crosscheck.py")
        #: ★ 行数守卫(2026-10-12 实测坑):导出腿少写一条腿时,逐列守卫
        #:   照绿(缺行的列在别的行里已被比过)。n = 2 条非升力腿 × 5 列
        #:   (cp_rms, cp_max, cl, cd, blas_threads) + lift_hex 4 列
        #:   (cl, cd, cl_rel, blas_threads) + lift_tet 3 列(cl, cl_rel,
        #:   blas_threads —— C06._one_level 不回 cd)= 17,一条腿丢了 n 就掉。
        assert n >= 17, (f"只比了 {n} 个数（2 nl 腿 × 5 列 + lift_hex × 4 列"
                         f" + lift_tet × 3 列 = 17）—— 有腿的行丢了")


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图/CSV 证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_evidence(runs, gate_evidence_dir):
    """写 `summary.csv` + 一张 Cp(θ) 对照图（两族 + 精确解，两条腿各一幅）。"""
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["leg", "blas_threads", "n_nodes", "n_tets", "converged",
                    "cp_rms", "cp_max", "cl", "cd", "cl_rel", "n_kutta",
                    "residual"])
        for leg in ("nl_hex", "nl_tet", "lift_hex", "lift_tet"):
            d = runs[leg]
            w.writerow([leg, blas_thread_count(), d["n_nodes"],
                        d.get("n_tets", "-"),
                        int(d["converged"]),
                        fmt(d["cp_rms"]) if "cp_rms" in d else "-",
                        fmt(d["cp_max"]) if "cp_max" in d else "-",
                        fmt(d["cl"]), fmt(d["cd"]) if "cd" in d else "-",
                        fmt(d["cl_rel"]) if "cl_rel" in d else "-",
                        d.get("n_kutta", "-"),
                        f"{d['residual']:.3e}" if "residual" in d else "-"])

    # 非升力腿的 Cp(θ) 重算一次（要曲线不要数）：hex 与 tet 各解一次画上去。
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    for k, (mesh_path, name) in enumerate(((MESH_HEX, "hex medium"),
                                           (MESH_TET, "tet medium"))):
        mesh = read_mesh(mesh_path)
        r = solve_newton_nonlifting(mesh, m_inf=0.0, precond="direct")
        d = _cp_error(mesh, np.asarray(r["phi"]))
        o = np.argsort(d["theta"])
        ax[0].plot(d["theta"][o], d["cp"][o], ".", ms=2.5,
                   label=f"{name} (rms {d['cp_rms']:.4f})")
    th = np.linspace(-np.pi, np.pi, 721)
    ax[0].plot(th, 1.0 - 4.0 * np.sin(th) ** 2, "k-", lw=1.6,
               label="exact: 1 − 4 sin²θ")
    ax[0].invert_yaxis()
    ax[0].set_xlabel("theta [rad]"), ax[0].set_ylabel("Cp")
    ax[0].set_title("Γ=0：两族在同一精确解上（攻击点 = 一族的系统缺陷）")
    ax[0].legend(fontsize=8), ax[0].grid(alpha=.3)

    names = ["lift_hex", "lift_tet"]
    vals = [100 * runs[n]["cl_rel"] for n in names]
    ax[1].bar(names, vals, color=["#1a56db", "#9aa0a6"])
    ax[1].axhline(100 * LIFT_CL_REL_MAX, color="#c0392b", ls="--",
                  label=f"shared band {100 * LIFT_CL_REL_MAX:.1f} %")
    ax[1].set_ylabel("cl_rel vs exact [%]")
    ax[1].set_title("Γ=0.6：各对各自精确解（hex 弦 2 / C06 族弦 1，口径在 docstring）")
    ax[1].legend(fontsize=8), ax[1].grid(alpha=.3, axis="y")
    fig.suptitle("C10  hex vs tet cylinder crosscheck", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir),
                             "c10_hex_cylinder_crosscheck.png"), dpi=130)
    plt.close(fig)
