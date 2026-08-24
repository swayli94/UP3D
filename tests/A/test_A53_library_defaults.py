"""A53 — 库默认值门（占位）。

承接台账 §1 的实测表：ew_eta0/ew_eta_max=1e-06、upwind_c=1.5、m_crit=0.95、
m_cap=3.0、rho_floor=0.05、precond='amg'、entropy_correction=True、
kutta_estimator='probe'、tip_taper=None、gamma_target=None、tip_cap='round'。
★ 现在**没有任何东西**检查它们有没有被静默改掉

★★ **占位（裁决 ④「暂时没有数据的先占位」）** —— 本文件目前只有一条 skip，
它的存在是为了让**门号被占住、并且缺什么写在这里**，而不是散落在计划文档里。
★ 判据形状（裁决一/三）：C 类判**收敛阶 + 绝对误差不能太大**；
D 类对 **Euler 设门**（无粘）· 对 **RANS 设门**（有粘耦合）· 对**实验记录偏置**。
"""
import pytest


@pytest.mark.skip(reason=(
    "PLACEHOLDER A53 — 门已编号、内容未实现。"
    " 缺什么见本文件 docstring。裁决 ④（2026-08-24）：暂时没有数据的先占位。"))
def test_library_defaults_placeholder():
    """占位：实现时把这条 skip 换成真判据，并按 G-TEETH 实测它会红。"""
