r"""D10 — pyFP3D FP+IBL vs CFL3D RANS — ONERA M6（3D-2）。

参考 `rans_onera_m6/`：M 0.8395 / α 3.06（TEST 2308 逐字），三档 × SST/SA，
Re 14.62e6 **根弦**。本门读 **L3**。

---

## ★★★ 这道门记录一个能力边界：三维粘性耦合在举力机翼上跑不起来

新鲜实测（coarse，M 0.8395 / α 3.06，trip 0.05 双面，tip_mask 0.05）：

    RuntimeError: FP driver did not converge at outer iter 1
                  (mdot_max = 6.218e-02)

**第 1 次外迭代**就退出。这与 GV5.0 记录的机制一致 ——
"coarse: root-upper-TE separation patch drives δ*↔ṁ↔u_e runaway (ṁ_max ×12.4
over k; the GV3.3-stern/Veldman class, first measurement on a lifting wing)"。

★★ **为什么重跑 coarse 而不是读归档。** GV5.3 在**同一工况**（TEST 2308 逐字）
下跑过并提交了粘性 Cp（`phases/p2/bench/studies/v5_3_m6_cp/results/`），
读它会是**跨时间**比较 —— phase-2 的代码状态对今天的 CFL3D，
即"两个数是不是同一个东西"家族的第四种。所以 coarse 在**今天的代码**上重跑；
**medium 不重跑**（GV5.3 记录 12479 s ≈ 3.5 h），它的归档读数
（cl_KJ 0.2685 → 0.2658，粘性使升力下降 1.03 %，"input-limited"）
在此**引用并标注跨时间**，不作为本门的断言依据。

## ★ 参考侧完整可用，等着将来

| 量 | SST L3 | SA L3 | \|SST−SA\| |
|---|---|---|---|
| cl | 0.272204 | 0.280566 | **3.07 %** |
| cd | 0.017492 | 0.018373 | 5.04 % |
| cd_friction | 0.005162 | 0.005747 | 11.3 % |

那条 **3.07 % 的带**就是这道门将来的分辨底噪。★ 顺带一提，它比参考自身的
网格差（cl 上 0.35 %）大 **9.5 倍** —— 见 `rans_onera_m6/` 的说明：
**主导不确定度是湍流模型，不是网格**。

## ★ 雷诺数口径（查过，不是假定）

GV5.0/GV5.3 用 Re_MAC 11.72e6 表达为**每米**（其网格是 NASA 米制），
CFL3D 参考用**根弦上的** 14.62e6（其网格归一到根弦 = 1）。
11.72e6 × (0.8059 / 0.64607) = **14.62e6** ⇒ **同一个物理雷诺数，两种参考长度。**
"""
import os

import numpy as np
import pytest

from tests.conftest import REPO_ROOT

REF_DIR = os.path.join(str(REPO_ROOT), "cases", "reference_data", "cfl3d",
                       "rans_onera_m6")
M_INF, ALPHA = 0.8395, 3.06
MAC, RE_MAC = 0.64607, 11.72e6
X_TR, TIP_FRAC = 0.05, 0.05


def _read_band(path=None):
    import csv
    p = path or os.path.join(REF_DIR, "forces.csv")
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            if r["level"] == "L3":
                out[r["turb"]] = dict(cl=float(r["cl"]), cd=float(r["cd"]),
                                      cdv=float(r["cd_friction"]))
    if "sst" not in out or "sa" not in out:
        raise RuntimeError(f"{p}: L3 needs both turbulence models")
    return out


@pytest.fixture(scope="module")
def probe():
    """跑一次三维耦合环，返回它**怎么失败的**，而不是让它炸掉测试。"""
    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.mesh.wake_cut import cut_wake
    from pyfp3d.meshgen.wing3d import chord_at, x_le
    from pyfp3d.viscous.coupling import (CouplingConfig, build_wing_case,
                                         make_picard_lifting_driver,
                                         run_loose_coupling)
    mc, wc = cut_wake(read_mesh(os.path.join(
        str(REPO_ROOT), "cases", "meshes", "onera_m6", "coarse.msh")))
    cfg = CouplingConfig(re_chord=RE_MAC / MAC, m_inf=M_INF, alpha_deg=ALPHA,
                         x_tr_upper=X_TR, x_tr_lower=X_TR)
    case = build_wing_case(mc.nodes, mc.elements, mc.boundary_faces["wall"],
                           cfg, x_le=x_le, chord_at=chord_at,
                           tip_mask_frac=TIP_FRAC)
    try:
        with np.errstate(invalid="ignore"):
            res = run_loose_coupling(
                make_picard_lifting_driver(mc, wc, M_INF, ALPHA), case, cfg)
        return dict(ok=True, err=None, msg="", n_outer=int(res.n_outer))
    except Exception as e:                  # noqa: BLE001 -- the point
        return dict(ok=False, err=type(e).__name__, msg=str(e), n_outer=None)


class TestTheCapabilityLimitIsRecorded:
    def test_the_3d_viscous_loop_still_fails(self, probe):
        """★★ 红了是**好消息**：三维粘性耦合跑通了，这道门应改写成
        "落在 SST–SA 带内"的真判据（`_read_band` 已备好，带宽 3.07 %）。"""
        assert not probe["ok"], (
            "the 3-D viscous loop now RUNS on the lifting wing -- an "
            "improvement, not a regression.  Re-specify D10 as a "
            "band-containment gate against rans_onera_m6/ (|SST-SA| 3.07 % "
            "in cl).")

    def test_it_fails_in_the_recorded_way(self, probe):
        """★★★ 锁的是**失败的形状**，不只是"它失败了"。

        实测：`RuntimeError: FP driver did not converge at outer iter 1
        (mdot_max = 6.218e-02)` —— 即 δ*↔ṁ↔u_e 反馈在**第一次**外迭代就把
        FP 驱动推离收敛（GV5.0 的 GV3.3-stern/Veldman 类）。
        换一种失败方式 = 故障变了 = 需要重新诊断，本条会红。"""
        assert probe["err"] == "RuntimeError", (
            f"the failure changed type: {probe['err']}: {probe['msg'][:90]}")
        assert "did not converge at outer iter" in probe["msg"], (
            f"the failure changed shape: {probe['msg'][:120]}")
        n = probe["msg"].split("outer iter")[1].split("(")[0].strip()
        assert int(n) <= 3, (
            f"the loop now survives to outer iter {n} (was 1) -- it is "
            f"getting further; re-read GV5.0's runaway description")


class TestReferenceIsLoadBearing:
    def test_gate_actually_reads_the_reference(self, tmp_path):
        import csv
        src = os.path.join(REF_DIR, "forces.csv")
        with open(src) as fh:
            rows = list(csv.DictReader(fh))
            cols = list(rows[0].keys())
        for r in rows:
            if r["level"] == "L3" and r["turb"] == "sst":
                r["cl"] = f"{float(r['cl']) + 0.1111:.6f}"
        dst = tmp_path / "forces.csv"
        with open(dst, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        base = _read_band()["sst"]["cl"]
        pert = _read_band(str(dst))["sst"]["cl"]
        assert abs((pert - base) - 0.1111) < 1e-9

    def test_the_band_is_the_future_resolution_floor(self):
        """★ 记录：|SST − SA| = 3.07 % in cl，是参考自身网格差（0.35 %）的
        9.5 倍 ⇒ 将来设门时**主导不确定度是湍流模型，不是网格**。"""
        b = _read_band()
        w = abs(b["sst"]["cl"] - b["sa"]["cl"]) / abs(b["sa"]["cl"])
        assert 0.01 < w < 0.10, f"implausible band width {w*100:.2f} %"

    def test_reynolds_caliber_matches_between_the_two_sides(self):
        """★★ pyFP3D 侧用 Re_MAC 每米、CFL3D 侧用根弦上的 Re。
        本条把那次换算锁住 —— 它是"两个数是不是同一个东西"的一次显式检查。"""
        import csv
        with open(os.path.join(REF_DIR, "forces.csv")) as fh:
            re_root = float(next(csv.DictReader(fh))["re_root_chord"])
        assert abs(RE_MAC * (0.8059 / MAC) / re_root - 1.0) < 0.01, (
            f"Re caliber mismatch: our Re_MAC {RE_MAC:.4g} on MAC {MAC} "
            f"rescales to {RE_MAC*(0.8059/MAC):.4g} on the root chord, but "
            f"the reference says {re_root:.4g}")
