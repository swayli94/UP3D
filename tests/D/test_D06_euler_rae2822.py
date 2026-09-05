r"""D06 — pyFP3D 无粘 vs CFL3D Euler — RAE2822（2D-5..2D-6）。

Case 7 M0.725/α2.55 · Case 9/10 M0.730/α3.19（`euler_rae2822/`，2 工况 × 3 档）。
α 一律**实验值不修正**（裁决二）。本门读参考的 **L3**。

---

## ★★★ 实测（先测后登记），以及第一次测量为什么是错的

第一次在 **6 线程**上测，得到"M0.730 medium 不收敛"，我据此写了一条
`assert not converged`。**那条断言自己红了**（测试跑在 8 线程），才把问题揪出来。
项目纪律 #1 是**上限 16，含 BLAS/OMP**。@16 线程重测：

| 工况（medium） | pyFP3D cl | CFL3D cl | Δcl | dx | @6t | @8t | **@16t** |
|---|---|---|---|---|---|---|---|
| M0.725/α2.55 | 1.016760 | 0.982077 | **+3.53 %** | +0.0309 | ✓ | ✓ | ✓ |
| M0.730/α3.19 | 1.113866 | 1.110904 | **+0.27 %** | +0.0051 | ✗ | ✓ | **✓** |

M0.730 medium 逐线程：4 线程 ✗（|R| 5.5e-06，80 步封顶）· 6 线程 ✗（2.9e-06）·
8 线程 ✓（**2.96e-13**）· 12 线程 ✓（3.1e-13），而 **cl 在四者间只差 0.15 %**。

★★★ **收敛旗标是算例与线程数的联合性质，不是算例的性质**，而**答案是稳定的**。
⇒ 本门**不对 `converged` 设任何断言**，只把它连同线程数写进证据；
**设门依据是分歧的大小**。

★ **误差跨零翻号**：两例都是 coarse 负、medium 正（−6.84 → +3.53、
−32.8 → +0.27）。跨零的"一致"更像**穿越**而非收缩 —— 与 D07 上 cl 三档
ratio 0.153、第四档翻成 3.161 是同一个签名。所以只在**幅度**上设带，
并明写它不单调收缩。

★ **方向不设门，因为它在两个翼型上相反**：NACA0012 上 pyFP3D 的激波在 CFL3D
**上游**（D05，4/5），RAE2822 上在**下游**（+0.031 / +0.005）。

## ★ 与 D05 相反的一点：这份参考紧到能分辨

RAE2822 参考自身的 L2→L3 激波差是 **0.00088 / 0.00100**，pyFP3D 与它差
0.0309 / 0.0051 ⇒ **35× / 5×，可分辨**。D05 的 M0.75 是 0.26×，不可分辨。
**同一条判据，两个相反的结论 —— 参照物的分辨率是逐例的性质。**
"""
import os

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.shock import shock_report
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting
from tests._gate_evidence import assert_matches_committed, fmt
from tests.conftest import REPO_ROOT, gate_figures_enabled

REF_DIR = os.path.join(str(REPO_ROOT), "cases", "reference_data", "cfl3d",
                       "euler_rae2822")
RECIPE = dict(upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80,
              n_picard_seed=5)
CASES = (
    ("M0.725_a2.55", "rae2822_m0725_a2.55", 0.725, 2.55, "gate"),
    ("M0.730_a3.19", "rae2822_m0730_a3.19", 0.730, 3.19, "gate"),
)
LEVELS = ("coarse", "medium")

# —— 判据（标定，不是实测值）——
CL_BAND = 0.06              # M0.725 medium |Δcl|：实测 0.0353
CL_BAND_M0730 = 0.03        # M0.730 medium |Δcl|：实测 0.0027
RESOLVABLE_FRAC = 1.0       # 与参考自身 L2→L3 差比较的可分辨下限


def _read_reference(path=None):
    """读已提交的 CFL3D 参考（L3）。★ 行为验证见 `TestReferenceIsLoadBearing`。"""
    import csv
    p = path or os.path.join(REF_DIR, "forces.csv")
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            if r["level"] == "L3":
                out[r["case"]] = dict(cl=float(r["cl"]), cd=float(r["cd"]))
    if not out:
        raise RuntimeError(f"{p}: no L3 rows -- the reference layout changed")
    return out


def _reference_uncertainty():
    import csv
    out = {}
    with open(os.path.join(REF_DIR, "grid_convergence.csv")) as fh:
        for r in csv.DictReader(fh):
            out[(r["case"], r["quantity"])] = abs(float(r["delta_L2_L3"]))
    return out


def _reference_shock(case):
    import csv
    with open(os.path.join(REF_DIR, "shock.csv")) as fh:
        for r in csv.DictReader(fh):
            if (r["case"] == case and r["level"] == "L3"
                    and r["surface"] == "upper" and r.get("x_shock")):
                return float(r["x_shock"])
    return float("nan")


def _one(level, m_inf, alpha):
    mc, wc = cut_wake(read_mesh(os.path.join(
        str(REPO_ROOT), "cases", "meshes", "rae2822_2.5d", f"{level}.msh")))
    dz = float(np.ptp(mc.nodes[:, 2]))
    r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=alpha, **RECIPE)
    phi = np.asarray(r["phi"])
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], phi,
                                alpha_deg=alpha, u_inf=1.0, s_ref=dz,
                                m_inf=m_inf)
    #: ★★ RAE2822 是**有弯度**的：上下面必须用局部外法向判，不能用 y > 0。
    #:   那正是 D12 修掉的库缺陷（37 个点被错分）。`section_cp_curve` 已修，
    #:   这里记下来是因为任何新写的截面处理都会再踩一次。
    cur = section_cp_curve(mc, phi, z=float(np.mean(mc.nodes[:, 2])),
                           smooth_passes=1, m_inf=m_inf)
    up = shock_report(cur, m_inf)["upper"]
    frz = r.get("sigma_freeze_report") or {}
    return dict(
        cl=float(f["cl"]), cd=float(f["cd_pressure"]),
        x_shock=(float(up["x_shock"]) if up.get("has_shock") else float("nan")),
        converged=bool(r.get("converged")),
        residual=float(np.asarray(r["residual_history"], float)[-1]),
        n_limited=int(r.get("n_limited", 0)),
        n_floored=int(r.get("n_floored", 0)),
        sigma_frozen=bool(frz.get("frozen_in_transient", False)),
    )


@pytest.fixture(scope="module")
def runs():
    return {(nm, lv): _one(lv, m, a) for nm, _c, m, a, _k in CASES
            for lv in LEVELS}


@pytest.fixture(scope="module")
def ref():
    return _read_reference()


class TestWhatIsGateable:
    def test_m0725_lift_magnitude(self, runs, ref):
        """幅度带，**不是**收缩判据：误差跨零（coarse −6.84 % → medium +3.53 %），
        所以不能要求单调收缩，只能要求 medium 的幅度在带内。"""
        got = runs[("M0.725_a2.55", "medium")]["cl"]
        want = ref["rae2822_m0725_a2.55"]["cl"]
        rel = abs(got - want) / abs(want)
        assert rel <= CL_BAND, (
            f"M0.725 cl {got:.6f} vs CFL3D Euler {want:.6f} = {rel*100:.2f} % "
            f"> {CL_BAND*100:.0f} %")

    def test_m0730_lift_magnitude(self, runs, ref):
        """@16 线程 Δcl = +0.27 %。★ 同样**不断言 `converged`**。"""
        got = runs[("M0.730_a3.19", "medium")]["cl"]
        want = ref["rae2822_m0730_a3.19"]["cl"]
        rel = abs(got - want) / abs(want)
        assert rel <= CL_BAND_M0730, (
            f"M0.730 cl {got:.6f} vs CFL3D Euler {want:.6f} = {rel*100:.2f} % "
            f"> {CL_BAND_M0730*100:.0f} %")

    def test_shock_difference_is_resolvable_against_this_reference(self, runs):
        """★★ 与 D05 的 M0.75 相反：这里差是参考自身不确定度的 35×，**可分辨**。

        锁住的是"可分辨"这个前提本身。若参考日后放松到不可分辨，本条会红，
        那时上面那条激波比较就该改成 RECORDED —— 与 D05 的镜像断言配对。"""
        unc = _reference_uncertainty()[("rae2822_m0725_a2.55", "x_shock_upper")]
        diff = abs(runs[("M0.725_a2.55", "medium")]["x_shock"]
                   - _reference_shock("rae2822_m0725_a2.55"))
        assert diff > RESOLVABLE_FRAC * unc, (
            f"the shock difference {diff:.6f} has fallen to or below the "
            f"reference's own uncertainty {unc:.6f} -- it is no longer "
            f"resolvable and must be re-specified as RECORDED")


class TestWhatIsOnlyRecorded:
    def test_convergence_flag_is_recorded_not_gated(self, runs):
        """★★★ 这里曾经有一条 `assert not converged`。

        它把一个**环境依赖的结果、而且是坏的那一侧**焊进了门里：M0.730 medium
        在 4/6 线程不收敛、8/12/16 线程收敛，而 cl 在四者间只差 0.15 %。
        **是那条断言自己红了**才暴露问题。⇒ 现在只断言旗标**可读**，
        它的值连同线程数写进证据。"""
        for nm, _c, _m, _a, _k in CASES:
            r = runs[(nm, "medium")]
            assert isinstance(r["converged"], bool)
            assert np.isfinite(r["residual"])

    def test_direction_of_the_shock_offset_is_not_gated(self, runs):
        """★ 记录方向，**不设门**：NACA0012 上 pyFP3D 的激波在上游，
        RAE2822 上在下游。没有一致方向 ⇒ 没有可登记的符号。
        本条只断言两个符号都被算出来了、可读。"""
        for nm, _c, _m, _a, _k in CASES:
            d = runs[(nm, "medium")]["x_shock"]
            assert np.isfinite(d), f"{nm}: no shock detected at medium"


class TestReferenceIsLoadBearing:
    def test_gate_actually_reads_the_reference(self, tmp_path):
        """扰动副本必须传播 —— F06 的教训，行为验证。"""
        import csv
        src = os.path.join(REF_DIR, "forces.csv")
        with open(src) as fh:
            rows = list(csv.DictReader(fh))
            cols = list(rows[0].keys())
        for r in rows:
            if r["case"] == "rae2822_m0725_a2.55" and r["level"] == "L3":
                r["cl"] = f"{float(r['cl']) + 0.4321:.6f}"
        dst = tmp_path / "forces.csv"
        with open(dst, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        base = _read_reference()["rae2822_m0725_a2.55"]["cl"]
        pert = _read_reference(str(dst))["rae2822_m0725_a2.55"]["cl"]
        assert abs((pert - base) - 0.4321) < 1e-9, (
            "the reader did not follow a perturbed reference")


class TestCommittedEvidenceIsLoadBearing:
    MEASURED = ("cl", "cd", "x_shock", "residual")

    def test_matches_committed_summary(self, runs, gate_evidence_dir):
        fresh = {f"{nm}|{lv}": {k: fmt(runs[(nm, lv)][k]) for k in self.MEASURED}
                 for nm, _c, _m, _a, _k in CASES for lv in LEVELS}
        assert_matches_committed(
            gate_evidence_dir, fresh, self.MEASURED,
            key_of=lambda r: f"{r['case']}|{r['level']}",
            refresh_hint="PYFP3D_GATE_FIGURES=1 pytest "
                         "tests/D/test_D06_euler_rae2822.py")


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图/CSV 证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_evidence(runs, ref, gate_evidence_dir):
    """★ 刷新是两遍：先带标志刷新，再不带标志验证。"""
    import csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    unc = _reference_uncertainty()
    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        # ★★ 线程数进证据：收敛旗标是"算例 x 线程数"的联合性质，
        #    不带线程数的旗标是无意义的（项目里 gated count 那条同理）。
        nthr = os.environ.get("NUMBA_NUM_THREADS", "unset")
        w.writerow(["case", "level", "kind", "n_threads", "mach", "alpha_deg", "cl", "cd",
                    "x_shock", "residual", "converged", "n_limited",
                    "n_floored", "sigma_frozen", "cl_ref_cfl3d", "d_cl_rel",
                    "x_shock_ref", "d_x_shock", "ref_shock_uncertainty",
                    "dx_over_ref_uncertainty"])
        for nm, cid, m, a, kind in CASES:
            xr = _reference_shock(cid)
            u = unc.get((cid, "x_shock_upper"), float("nan"))
            for lv in LEVELS:
                d = runs[(nm, lv)]
                dx = d["x_shock"] - xr
                w.writerow([nm, lv, kind, nthr, m, a, fmt(d["cl"]), fmt(d["cd"]),
                            fmt(d["x_shock"]), fmt(d["residual"]),
                            int(d["converged"]), d["n_limited"],
                            d["n_floored"], int(d["sigma_frozen"]),
                            f"{ref[cid]['cl']:.6f}",
                            f"{(d['cl'] - ref[cid]['cl']) / abs(ref[cid]['cl']):.6e}",
                            f"{xr:.6f}", f"{dx:+.6f}", f"{u:.6f}",
                            f"{abs(dx) / u:.2f}"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for k, (nm, cid, m, a, kind) in enumerate(CASES):
        ax = axes[k]
        d = runs[(nm, "medium")]
        ax.bar(["pyFP3D", "CFL3D Euler"], [d["cl"], ref[cid]["cl"]],
               color=["#1a56db", "#c0392b"])
        lab = "GATED" if kind == "gate" else "RECORDED"
        conv = "" if d["converged"] else "  NOT CONVERGED"
        ax.set_title(f'{nm}  [{lab}]{conv}\ncl {d["cl"]:.6f} vs '
                     f'{ref[cid]["cl"]:.6f}', fontsize=9,
                     color=("black" if d["converged"] else "#c0392b"))
        ax.grid(alpha=.3, axis="y")
    fig.suptitle("D06  pyFP3D inviscid vs CFL3D Euler on RAE2822 -- the "
                 "0.38 % agreement is on a leg that DID NOT CONVERGE",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "d06_vs_cfl3d_euler.png"),
                dpi=130)
    plt.close(fig)
