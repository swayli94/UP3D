r"""D09 — pyFP3D FP+IBL vs CFL3D RANS — RAE2822（R-5..R-6）。

工况 `rans_rae2822/`：M0.725/α2.55 与 M0.730/α3.19，Re 6.5e6，trip 0.03，
**每工况 SST + SA 两个湍流模型**。α 一律实验值不修正（裁决二）。

---

## ★★★ 这道门记录的是一个**能力边界**：耦合环在 RAE2822 上跑不起来

四条腿（2 工况 × coarse/medium）**全部**在第 **0** 次外迭代就死：

    RuntimeWarning: invalid value encountered in sqrt
      pyfp3d/viscous/coupling.py:801  mach_e = np.sqrt(mach_squared_field(...))
    → NaN 传进 numba 闭包 → SystemError: CPUDispatcher(closure_all)

**`SystemError` 只是症状**（numba 对 NaN 输入的表现），根因是 `mach_squared_field`
返回**负值**。实测（coarse, M0.725/α2.55）：

| | RAE2822（有弯度） | NACA0012（对称，同一套机制） |
|---|---|---|
| M² < 0 的站位 | **9 / 206** | **0 / 206** |
| min M² | **−63.90** | — |
| max q² | **23.60**（\|u_e\| = 4.86 倍来流）| **1.98** |
| 位置 x/c | 0.227–0.287 与 0.744–0.764，**全在上表面** | — |

激波在 x/c ≈ 0.60–0.63，所以坏点分布在激波**上游**（超声速平台）和**下游**两处。

★★ **定位到，但不归因。** D06 在**同一张网格**上测得无粘解是健康的
（收敛，cl 0.914939，激波 0.6049）⇒ FP 解没问题，问题在
`run_loose_coupling` 的**站位边缘速度提取**。

★★★ **一个被测量否掉的假设，写在这里免得有人重走。** 我先怀疑是分侧错误 ——
D12 测出 `section_cut.py` 用全局 y 基准在 RAE2822 上错分 37 个点，而
`coupling.py:334` 用的是同一族规则（`ysum/ycnt >= y_stag`），且该修复**没有
回移**（纪律 #9）。**证伪**：那 9 个点的 y 分别是 +0.054…+0.059（x/c 0.25）
和 +0.040…+0.042（x/c 0.75），**正是 RAE2822 上表面该有的位置**；而"一个站位
的节点全在一侧"在**对称的 NACA0012 上同样是 101/102**，即它是正常结构而不是
缺陷。分侧是对的。（`coupling.py:334` 的全局 y 基准仍然是一个**未被这次测量
触发**的隐患，单独记着。）

★ 与 GV5.2 的关系：那一轮在 RAE2822 上用**松耦合 + Newton**协议测过同样两个
工况，记录 "P2 loop runaway at k = 4, mdot_max = 1.59"。**那是第 4 次外迭代**，
本轮是**第 0 次**就死 —— 协议不同（这里是 Picard 驱动 + 固定 trip 0.03），
所以两者是**不同的读数，不是矛盾**。

## ★ 参考侧仍然是完整可用的，只是无人可比

`rans_rae2822/` 的两个工况各有 SST/SA，|SST−SA| 在 cl 上 **9.3 % / 11.0 %** ——
这道门将来能用的**分辨底噪**已经在那里了。
"""
import os

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.viscous.coupling import (CouplingConfig, build_airfoil_case,
                                     make_picard_lifting_driver,
                                     run_loose_coupling)
from tests.conftest import REPO_ROOT

REF_DIR = os.path.join(str(REPO_ROOT), "cases", "reference_data", "cfl3d",
                       "rans_rae2822")
CASES = (("M0.725_a2.55", "rae2822_m0725_a2.55", 0.725, 2.55),
         ("M0.730_a3.19", "rae2822_m0730_a3.19", 0.730, 3.19))
RE, X_TR = 6.5e6, 0.03
#: 实测 9/206；带一点余量，因为它记录的是"仍然坏着"，不是一个要收紧的目标
MAX_BAD_STATIONS = 20


def _read_band(path=None):
    """每工况的 (SST, SA) cl —— 参考是一条**带**，不是一个点。"""
    import csv
    p = path or os.path.join(REF_DIR, "forces.csv")
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            if r["level"] == "L3":
                out.setdefault(r["case"], {})[r["turb_model"]] = float(r["cl"])
    if not out:
        raise RuntimeError(f"{p}: no L3 rows -- the reference layout changed")
    return out


def _probe(fam, m, a, xtr, level="coarse"):
    """跑一次耦合环，返回它**死在哪里**（而不是让它把测试炸掉）。"""
    import pyfp3d.viscous.coupling as CP
    grab = {}
    orig = CP.mach_squared_field

    def spy(q2, m_inf, gamma):
        v = orig(q2, m_inf, gamma)
        grab.setdefault("q2", np.array(q2))
        grab.setdefault("m2", np.array(v))
        return v

    CP.mach_squared_field = spy
    try:
        mc, wc = cut_wake(read_mesh(os.path.join(
            str(REPO_ROOT), "cases", "meshes", fam, f"{level}.msh")))
        cfg = CouplingConfig(re_chord=RE, m_inf=m, alpha_deg=a,
                             x_tr_upper=xtr, x_tr_lower=xtr, n_outer_max=10)
        case = build_airfoil_case(mc.nodes, mc.elements,
                                  mc.boundary_faces["wall"], cfg)
        err = None
        try:
            with np.errstate(invalid="ignore"):
                res = run_loose_coupling(
                    make_picard_lifting_driver(mc, wc, m, a), case, cfg)
            ok = bool(res.converged)
        except Exception as e:              # noqa: BLE001 -- the point is to catch it
            ok, err = False, type(e).__name__
    finally:
        CP.mach_squared_field = orig
    m2 = grab.get("m2")
    return dict(ok=ok, err=err,
                n_stations=(0 if m2 is None else int(m2.size)),
                n_bad=(0 if m2 is None else int((m2 < 0).sum())),
                min_m2=(float("nan") if m2 is None else float(m2.min())),
                max_q2=(float("nan") if "q2" not in grab
                        else float(grab["q2"].max())))


@pytest.fixture(scope="module")
def probes():
    return {"rae": _probe("rae2822_2.5d", 0.725, 2.55, X_TR),
            "naca_control": _probe("naca0012_2.5d", 0.50, 2.00, 0.05)}


class TestTheCapabilityLimitIsRecorded:
    r"""★ 记录的是**边界**，不是一个数值判据。每条断言都写明"红了意味着什么"。"""

    def test_rae2822_loop_still_fails_at_the_first_outer_iteration(self, probes):
        """★★ 红了是**好消息**：耦合环在 RAE2822 上跑通了，
        这道门应当改写成"落在 SST–SA 带内"的真判据（`_read_band` 已备好）。"""
        p = probes["rae"]
        assert not p["ok"], (
            "the RAE2822 loose loop now RUNS -- this is an improvement, not a "
            "regression.  Re-specify D09 as a band-containment gate against "
            "the SST/SA reference (|SST-SA| is 9.3 % / 11.0 % in cl).")

    def test_the_failure_is_a_negative_mach_squared_not_something_else(self,
                                                                      probes):
        """★★★ 锁的是**根因**，不是症状。`SystemError: CPUDispatcher` 是 numba
        遇到 NaN 的表现；真正的事实是 `mach_squared_field` 返回负值。
        如果哪天它以别的方式失败，本条会红，那说明故障换了，需要重新诊断。"""
        p = probes["rae"]
        assert p["n_bad"] > 0, (
            "the RAE2822 loop no longer produces M^2 < 0 -- the recorded "
            "failure mode has changed and needs re-diagnosing")
        assert p["n_bad"] <= MAX_BAD_STATIONS, (
            f"{p['n_bad']} stations with M^2 < 0, was 9 of 206 -- the failure "
            f"has grown substantially")
        assert p["max_q2"] > 5.0, (
            f"max q2 {p['max_q2']:.3f}: the recorded blow-up (23.6, i.e. "
            f"|u_e| 4.9x freestream) is gone but the loop still fails -- "
            f"a different mechanism")

    def test_the_symmetric_section_is_clean_on_the_same_machinery(self, probes):
        """★★ 这条是**判别项**：同一套机制在对称翼型上完全干净
        （0/206，max q² 1.98）。没有它，上面的失败可能只是"IBL 一概不能用"。"""
        c = probes["naca_control"]
        assert c["n_bad"] == 0, (
            f"the NACA0012 control now also produces {c['n_bad']} stations "
            f"with M^2 < 0 -- the failure is no longer camber-specific and "
            f"D08's gates are in doubt too")
        assert c["max_q2"] < 5.0


class TestReferenceIsLoadBearing:
    def test_the_band_is_read_from_the_committed_file(self, tmp_path):
        """★ 参考虽然暂时无人可比，读取器仍要行为验证 —— 等 D09 能设门时
        它就是那条判据的地基。"""
        import csv
        src = os.path.join(REF_DIR, "forces.csv")
        with open(src) as fh:
            rows = list(csv.DictReader(fh))
            cols = list(rows[0].keys())
        for r in rows:
            if (r["case"] == "rae2822_m0725_a2.55" and r["level"] == "L3"
                    and r["turb_model"] == "sst"):
                r["cl"] = f"{float(r['cl']) + 0.2468:.6f}"
        dst = tmp_path / "forces.csv"
        with open(dst, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        base = _read_band()["rae2822_m0725_a2.55"]["sst"]
        pert = _read_band(str(dst))["rae2822_m0725_a2.55"]["sst"]
        assert abs((pert - base) - 0.2468) < 1e-9

    def test_the_reference_band_width_is_the_future_resolution_floor(self):
        """★ 记录（不设门）：|SST − SA| 是这道门将来的分辨底噪。"""
        b = _read_band()
        for cid, d in b.items():
            w = abs(d["sst"] - d["sa"]) / abs(d["sa"])
            assert 0.0 < w < 0.5, f"{cid}: implausible band width {w*100:.1f} %"
