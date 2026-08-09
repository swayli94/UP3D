# level-set 删除清单(phase three 的第一件工作)

依据裁决 **D5**(2026-08-09,使用者):**放弃 level-set 路线,未来只开发 conforming**;
phase two **冻结不删**,**删除是 phase three 的第一件工作**。
本文件把面点清,让删除变成**对着清单勾**,而不是考古。

★ **口径**:下面每个"是 LS"的判定都来自 **AST 扫真实 `import`**,不是名字子串 `grep` ——
后者会把注释算成引用(`newton.py` 提到 `newton_ls` 只是 backport 纪律的注释)。
`onera_m6_roundtip` 曾被我按名字误判为 LS 专用,**核实消费者后确认它是 conforming 的**(见 §4)。

## 1. 库:9 个文件 / **4624 行**,只有**一个**集成点

| 行数 | 文件 | 是什么 |
|---|---|---|
| 1196 | `pyfp3d/solve/newton_ls.py` | LS 的耦合 Newton 驱动器(含 freeze 爬坡) |
| 1044 | `pyfp3d/solve/picard_ls.py` | LS 的 Picard 孪生 |
| 907 | `pyfp3d/wake/multivalued.py` | 多值算子(aux DOF、TE 跳跃、element_mach2) |
| 556 | `pyfp3d/kernels/cut_assembly.py` | 切割单元装配核 |
| 333 | `pyfp3d/wake/cut_elements.py` | 切割映射(含 `inboard_clip`、`outboard_fringe`) |
| 225 | `pyfp3d/wake/levelset.py` | 尾迹水平集(含 `sheet_direction`) |
| 186 | `pyfp3d/post/surface_ls.py` | LS 侧后处理 |
| 165 | `pyfp3d/solve/schur_ls.py` | B14 的 Schur+AMG 预条件子(**只为 LS 造的**) |
| 12 | `pyfp3d/wake/__init__.py` | 导出面 |
| **4624** | | |

★★ **唯一的集成点:`pyfp3d/post/unified.py`** —— 全库只有它 `import` 进 LS 侧
(`surface_ls` 的 `_d11_wall_state`、`section_cp_curve_levelset`)。它的 docstring 明说自己是
**两条路径的统一入口**,所以删 LS 之后它**塌缩成 conforming 那一半**,不是删掉。

⚠ **更正 roadmap**:GS5.1 写的是"约 **2500** 行",**实测 4624 行(仅库)**,差 **1.85 倍**。

## 2. 测试:18 个整删 / 8 个摘腿

**纯 LS,可整文件删(18 个 / 4297 行)**

```
test_b2_multivalued  test_b6_newton  test_b6_transonic  test_b8_span_blend
test_b8_tip_taper_ls  test_b9_wingbody_ls  test_b11_linear_ls  test_b11_post_unified
test_b12_lagged_lu_ls  test_b13_lagged_picard  test_b14_schur_ls  test_b15_ls_newton_freeze
test_b16_farfield_aux  test_b17_farfield_pin_gamma  test_b19_jacobian_3d
test_b22_ls_3d_anchors  test_b31_tip_fringe  test_m2_wingbody
```

**混合,只删 LS 腿(8 个 / 2248 行)** —— 删这些**不是纯减法**,见 §5:

| 文件 | LS 腿承载什么 |
|---|---|
| `test_a1_instrumentation.py` | 四驱动器的计时/成本对比(A1 的结论) |
| `test_b1_cut_elements.py` | 切割单元的结构锁(B25 的 `inboard_clip` 等) |
| `test_b3_lifting.py` / `test_b45_farfield.py` | B3/B5 的远场 A/B(**B5 的域稳健性结论在这里**) |
| `test_b4_te_control_volume.py` | TE 控制体两路径对照 |
| `test_b7_onera_m6.py` | M6 LS 3-D 机构 + 跨声速门(3 条已 strict-xfail) |
| `test_b18_wingbody_transonic.py` | 翼身跨声速的两路径对照 |
| ★ `test_v2_newton_rhs_channel.py` | **Track V(粘性)** 的 transpiration 通道走 LS 的 `b_base` |

## 3. bench / demo

| | 纯 LS(整删) | 混合(摘腿) |
|---|---|---|
| **bench** | 5 个 / 592 行:`run_b22_reanchor`、`run_b7_reanchor`、`run_capability_budget`、`run_le13_roundtip_envelope`、`run_le15_roundtip_mode` | 1 个 / 567 行:**`run_capability_matrix`**(13 个 cell 里含 LS 的) |
| **cases/demo** | **14 个 / 4599 行**:b3、b6、b7、b8(×4)、b11(×3)、b12、b13、b14、b15、b16、b17、m6_medium_ls_workflow | 8 个 / 4013 行:b9、b18、b4p5、**p13**、**p14** 等(跨模型对比) |

## 4. 网格:**只有两个族**是 LS 专用

| 族 | 判定 | 消费者 |
|---|---|---|
| `cases/meshes/naca0012_wakefree_2.5d` | **LS 专用,可删** | 全是 b1/b2/b3/b45/b6 |
| `cases/meshes/onera_m6_wakefree` | **LS 专用,可删** | 全是 b1/b2/b7/b22 |
| `cases/meshes/onera_m6_roundtip` | ★ **不是 LS —— 留下** | `test_m5_round_tip`、M5 demo、P13 的圆角脚本(**都是 conforming**) |

★ `.msh` 本身 gitignored,要删的是**生成器脚本 + 已提交的 stats CSV / PNG**。

## 5. ★★ 删除时会踩到的三件事(这不是纯减法)

1. **混合文件里的 LS 腿承载着交叉验证** —— P14 的跨模型一致(**0.17% / 0.36%**)、
   A1 的四驱动器成本对比、B18 的两路径对照、**B5 的远场域稳健性**(Γ 在 15→120c 内
   0.45%/1.09%)。删掉之后这些结论**只剩已提交的 CSV/PNG,不再可重跑**。
   ⇒ 这正是 **D5 记下的代价**;phase three 的能力边界必须在**没有交叉验证**的前提下陈述。
2. ★ **`test_v2_newton_rhs_channel.py` 是粘性侧的测试**,却走 LS 的 `b_base` 通道 ⇒
   删 LS 会**碰到 Track V 的接线**,phase three 需要为 V2 的 transpiration 通道
   **在 conforming 上补一条等价测试**,否则粘性通道会失去一条锁。
3. **`post/unified.py` 是塌缩不是删除** —— 它是设计好的两路径统一入口(§1)。

## 6. 建议的删除顺序(降低中途破坏的风险)

1. 先删 **demo 与 bench 的纯 LS 文件**(不影响任何测试)
2. 再删 **18 个纯 LS 测试文件**,跑常开全套 + 快层确认基线只减不乱
3. 再**摘 8 个混合测试的 LS 腿**,每摘一个跑一次(它们含结构锁与跨路径对照)
4. ★ **先给 V2 的 transpiration 通道补 conforming 等价测试**,再动 `b_base`
5. 然后**塌缩 `post/unified.py`**
6. 最后删 **9 个库文件 + 2 个网格族生成器**,跑门控全套(**3:04:44 @8 线程**)
7. 收口:六面 + 快层(**10.7 min**),并更正 roadmap 的"约 2500 行"

## 7. 删完之后应当只剩的一条尾迹路线

**conforming**:`mesh/wake_cut.py`(尾迹切割)+ `constraints/te_pressure.py`(P14 压强相等 Kutta)
+ `solve/newton.py` / `solve/picard.py`,配 `tip_taper`(B31/B32,带 **−1.3% cl** 的模型偏差)。
