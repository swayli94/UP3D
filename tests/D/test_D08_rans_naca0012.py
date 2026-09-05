r"""D08 — pyFP3D FP+IBL vs CFL3D RANS — NACA0012（R-1..R-4）。

参考 `rans_naca0012/`：**每工况 SST + SA**，所以参考是一条**带**而不是一个点，
`|SST − SA|` 就是这道门的分辨底噪（数据请求的原话）。

---

## ★★ 判据的形状：**落在带内**，不是"与某个模型差 X %"

后者要武断地挑一个湍流模型当真值。前者用的正是两模型设计的用意。

| 工况（trip 匹配）| pyFP3D coarse | pyFP3D medium | RANS 带 [SST, SA] | 带宽 |
|---|---|---|---|---|
| M0.50/α2.0 **xtr 0.05** | 0.275298（+7.2 %）| 0.283039（**+10.2 %**）| [0.2500, 0.2569] | 2.7 % |
| M0.50/α2.0 **xtr 0.30** | 0.256161（**−1.9 %**）| 0.281915（**+6.1 %**）| [0.2612, 0.2658] | 1.8 % |

★★★ **两条腿都在加密时向上跑，而且 xtr 0.30 从带下方穿到带上方。**
⇒ 没有一条腿"落在带内"，而且**加密使它变差** —— 与 D13 记录的
"cl 离 XFOIL 从 +2.3 % 走到 +5.6 %" 同一签名（那一轮把根因指到湍流闭合的
**H 族**：H 偏小 ⇒ 摩擦偏高 ⇒ 边界层过于饱满）。本门**记录**这个方向，
把带内包含设成**将来**的判据。

## ★ 三条工况被排除，理由是**口径**不是难度

`rans_naca0012/` 另有 M0.352/α12.86、M0.778/α2.03、M0.803/α−0.1，
CFL3D 侧都是**全湍流**（`x_tr` 空、`x_tr_actual` 空），而 pyFP3D 的 IBL 用
**固定转捩点**。⇒ **口径不同，不是同一个量**。
★ 而且参考数据自己的 `note` 就写着 near-stall 那条
"NEAR STALL -- widest model spread, deliberately NOT gate material"。
**参考的作者已经做过一次设门/记录的判断，读它比自己重下判断更可靠。**
"""
import os

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.viscous.coupling import (CouplingConfig, build_airfoil_case,
                                     make_picard_lifting_driver,
                                     run_loose_coupling)
from tests._gate_evidence import assert_matches_committed, fmt
from tests.conftest import REPO_ROOT, gate_figures_enabled

REF_DIR = os.path.join(str(REPO_ROOT), "cases", "reference_data", "cfl3d",
                       "rans_naca0012")
#: 只有 trip 匹配的两条 —— 见 docstring 的口径段
CASES = (("xtr005", "n0012_m0500_a2.00_xtr005", 0.50, 2.00, 3.0e6, 0.05),
         ("xtr030", "n0012_m0500_a2.00_xtr030", 0.50, 2.00, 3.0e6, 0.30))
LEVELS = ("coarse", "medium")

# —— 判据（标定）——
ABOVE_BAND_MAX = 0.20     # 记录带：medium 高出 RANS 带的上限，实测 0.102 / 0.061
BAND_WIDTH_MAX = 0.10     # 参考带宽的合理上限（实测 0.027 / 0.018）


def _read_band(path=None):
    import csv
    p = path or os.path.join(REF_DIR, "forces.csv")
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            if r["level"] == "L3":
                out.setdefault(r["case"], {})[r["turb_model"]] = dict(
                    cl=float(r["cl"]), cd=float(r["cd"]),
                    cdv=float(r["cd_friction"]))
    if not out:
        raise RuntimeError(f"{p}: no L3 rows -- the reference layout changed")
    return out


def _transition_caliber(path=None):
    """每工况的 trip 口径 —— 用来证明被排除的工况是**口径**问题。"""
    import csv
    p = path or os.path.join(REF_DIR, "forces.csv")
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            if r["level"] == "L3":
                out[r["case"]] = dict(x_tr=r["x_tr"], note=r["note"])
    return out


def _one(level, m, a, re_c, xtr):
    mc, wc = cut_wake(read_mesh(os.path.join(
        str(REPO_ROOT), "cases", "meshes", "naca0012_2.5d", f"{level}.msh")))
    cfg = CouplingConfig(re_chord=re_c, m_inf=m, alpha_deg=a,
                         x_tr_upper=xtr, x_tr_lower=xtr, n_outer_max=10)
    case = build_airfoil_case(mc.nodes, mc.elements,
                              mc.boundary_faces["wall"], cfg)
    res = run_loose_coupling(make_picard_lifting_driver(mc, wc, m, a),
                             case, cfg)
    dz = float(np.ptp(mc.nodes[:, 2]))
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"],
                                np.asarray(res.phi), alpha_deg=a, u_inf=1.0,
                                s_ref=dz, m_inf=m)
    return dict(cl=float(f["cl"]), n_outer=int(res.n_outer),
                converged=bool(res.converged))


@pytest.fixture(scope="module")
def runs():
    return {(nm, lv): _one(lv, m, a, re_c, xtr)
            for nm, _c, m, a, re_c, xtr in CASES for lv in LEVELS}


@pytest.fixture(scope="module")
def band():
    return _read_band()


def _rel_to_band(cl, lo, hi):
    """0 表示落在带内；正/负表示高出/低于带，按带边相对。"""
    if lo <= cl <= hi:
        return 0.0
    return (cl - hi) / abs(hi) if cl > hi else (cl - lo) / abs(lo)


class TestWhatIsGateable:
    r"""★ 目前**没有一条腿落在带内**，所以能设门的是"分歧仍在记录带内"
    加上参考侧那条带本身是良态的。"""

    def test_reference_band_is_well_posed(self, band):
        """★ 带宽必须有限且非零 —— 一个塌成一点的带会让下面的判据变成
        "与某个模型比"，那正是本门刻意避开的形状。"""
        for _nm, cid, *_ in CASES:
            d = band[cid]
            w = abs(d["sst"]["cl"] - d["sa"]["cl"]) / abs(d["sa"]["cl"])
            assert 0.001 < w < BAND_WIDTH_MAX, (
                f"{cid}: band width {w*100:.2f} % is not usable as a "
                f"resolution floor")

    @pytest.mark.parametrize("name", ["xtr005", "xtr030"])
    def test_medium_stays_within_the_recording_band(self, runs, band, name):
        """记录带（20 %），不是一致性判据。实测 +10.2 % / +6.1 %。"""
        cid = next(c for n, c, *_ in CASES if n == name)
        d = band[cid]
        lo, hi = sorted((d["sst"]["cl"], d["sa"]["cl"]))
        r = _rel_to_band(runs[(name, "medium")]["cl"], lo, hi)
        assert abs(r) <= ABOVE_BAND_MAX, (
            f"{name}: medium cl {runs[(name,'medium')]['cl']:.6f} is "
            f"{r*100:+.1f} % from the RANS band [{lo:.4f}, {hi:.4f}], beyond "
            f"the recording band {ABOVE_BAND_MAX*100:.0f} %")

    def test_excluded_cases_are_excluded_for_CALIBER_not_difficulty(self):
        """★★ 三条被排除的工况在参考侧是**全湍流**，与我们的固定转捩不是
        同一个量。本条锁住那个事实：如果参考哪天给了它们一个 trip，
        它们就该被纳入。"""
        cal = _transition_caliber()
        for cid in ("n0012_m0352_a12.86", "n0012_m0778_a2.03",
                    "n0012_m0803_am0.10"):
            assert cal[cid]["x_tr"] == "", (
                f"{cid} now carries a trip ({cal[cid]['x_tr']}) -- its "
                f"transition caliber matches ours and it should be added to "
                f"CASES")
        assert "NOT gate material" in cal["n0012_m0352_a12.86"]["note"], (
            "the reference no longer flags the near-stall case as non-gate "
            "material -- re-read its note before using it")


class TestWhatIsOnlyRecorded:
    def test_refinement_moves_away_from_the_band(self, runs, band):
        """★★★ 两条腿都在加密时向上跑，xtr 0.30 还从带下穿到带上。

        红了说明方向反转（加密开始靠近），那是**好消息**，应把带内包含
        升级为真判据。"""
        moved_up = 0
        for nm, cid, *_ in CASES:
            d = band[cid]
            lo, hi = sorted((d["sst"]["cl"], d["sa"]["cl"]))
            c = _rel_to_band(runs[(nm, "coarse")]["cl"], lo, hi)
            m = _rel_to_band(runs[(nm, "medium")]["cl"], lo, hi)
            if m > c:
                moved_up += 1
        assert moved_up == len(CASES), (
            f"only {moved_up} of {len(CASES)} legs still move AWAY from the "
            f"band under refinement -- an improvement; re-specify the "
            f"band-containment criterion as a gate")


class TestReferenceIsLoadBearing:
    def test_gate_actually_reads_the_reference(self, tmp_path):
        import csv
        src = os.path.join(REF_DIR, "forces.csv")
        with open(src) as fh:
            rows = list(csv.DictReader(fh))
            cols = list(rows[0].keys())
        for r in rows:
            if (r["case"] == "n0012_m0500_a2.00_xtr005"
                    and r["level"] == "L3" and r["turb_model"] == "sst"):
                r["cl"] = f"{float(r['cl']) + 0.3141:.6f}"
        dst = tmp_path / "forces.csv"
        with open(dst, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        base = _read_band()["n0012_m0500_a2.00_xtr005"]["sst"]["cl"]
        pert = _read_band(str(dst))["n0012_m0500_a2.00_xtr005"]["sst"]["cl"]
        assert abs((pert - base) - 0.3141) < 1e-9


class TestCommittedEvidenceIsLoadBearing:
    MEASURED = ("cl",)

    def test_matches_committed_summary(self, runs, gate_evidence_dir):
        fresh = {f"{nm}|{lv}": {"cl": fmt(runs[(nm, lv)]["cl"])}
                 for nm, *_ in CASES for lv in LEVELS}
        assert_matches_committed(
            gate_evidence_dir, fresh, self.MEASURED,
            key_of=lambda r: f"{r['case']}|{r['level']}",
            refresh_hint="PYFP3D_GATE_FIGURES=1 pytest "
                         "tests/D/test_D08_rans_naca0012.py")


@pytest.mark.skipif(not gate_figures_enabled(),
                    reason="图/CSV 证据是 opt-in：PYFP3D_GATE_FIGURES=1")
def test_export_evidence(runs, band, gate_evidence_dir):
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nthr = os.environ.get("NUMBA_NUM_THREADS", "unset")
    with open(os.path.join(str(gate_evidence_dir), "summary.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "level", "n_threads", "x_tr", "cl", "cl_sst",
                    "cl_sa", "band_lo", "band_hi", "rel_to_band",
                    "n_outer", "converged"])
        for nm, cid, _m, _a, _re, xtr in CASES:
            d = band[cid]
            lo, hi = sorted((d["sst"]["cl"], d["sa"]["cl"]))
            for lv in LEVELS:
                r = runs[(nm, lv)]
                w.writerow([nm, lv, nthr, xtr, fmt(r["cl"]),
                            f'{d["sst"]["cl"]:.6f}', f'{d["sa"]["cl"]:.6f}',
                            f"{lo:.6f}", f"{hi:.6f}",
                            f'{_rel_to_band(r["cl"], lo, hi):+.6e}',
                            r["n_outer"], int(r["converged"])])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for k, (nm, cid, *_rest) in enumerate(CASES):
        ax = axes[k]
        d = band[cid]
        lo, hi = sorted((d["sst"]["cl"], d["sa"]["cl"]))
        ax.axhspan(lo, hi, color="#c0392b", alpha=.2,
                   label="CFL3D RANS band (SST-SA)")
        ax.plot(LEVELS, [runs[(nm, lv)]["cl"] for lv in LEVELS], "o-",
                color="#1a56db", label="pyFP3D FP+IBL")
        ax.set_title(f"{nm}  (trip {_rest[-1]})\nrefinement moves AWAY",
                     fontsize=9, color="#c0392b")
        ax.set_ylabel("cl"), ax.grid(alpha=.3), ax.legend(fontsize=8)
    fig.suptitle("D08  pyFP3D FP+IBL vs the CFL3D RANS BAND -- the reference "
                 "is two turbulence models, so the criterion is containment",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(str(gate_evidence_dir), "d08_vs_rans_band.png"),
                dpi=130)
    plt.close(fig)
