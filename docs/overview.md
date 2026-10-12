# pyFP3D 总览

**pyFP3D**：三维非结构网格**全势流（full-potential）跨声速**外流求解器，Python + Numba。
每个节点一个标量未知量 φ，Galerkin P1 四面体，超声速区人工密度迎风，尾迹切口 + Kutta
条件产生升力。目标场景：机翼/翼身组合体，M∞ 0.3–0.87（★ α=0 的**实测**可信上界只有
**M 0.78–0.80**，见能力边界），工作站尺度（1–3M 节点分钟级）。

> **本文件是人读的高层总览，不是权威进度源。**
> 当前阶段指针与工作纪律 → [agent-rules.md](agent-rules.md)（经 CLAUDE.md `@include` 注入）；
> 理论与数值 → [design.md](design.md)；粘性耦合线 → [design_track_v.md](design_track_v.md)；
> **现行能力声明（实测边界）** →
> [phases/p6/docs/dev_phase_six/PHASE_SIX_CAPABILITY_BOUNDARY.md](../phases/p6/docs/dev_phase_six/PHASE_SIX_CAPABILITY_BOUNDARY.md)（2026-10-11）；
> 已结束阶段的归档 → [phases/README.md](../phases/README.md)。

## 当前状态（2026-10-12，phase 6 收口实测）

- 未门控套件 **825 passed / 0 failed / 33 skipped / 2 xfailed**；门控全集 812/21/3/1
  （含 1 XPASS）；快层 5/5。**BLAS 由库运行时钉成单线程**（证据 CSV 带 `blas_threads`
  出处列）；证据生成/验证必须在**暖缓存态**（重编译后先白跑一遍再测）。
- 测试按目录分层：`tests/A`+`B` 未门控（机制/回归）、`tests/C`+`E` 快层（每次收口跑）、
  `tests/D` 门控全集（外部参照/相位边界）、`tests/F` 仪器（仓库不变量）。
- 三条实测能力边界（每条都有门）：① α=0 跨声速非唯一（D14）；② **有粘耦合在举力机翼上
  跑不起来**（D09/D10）；③ 生产配置在精确解上的已登记代价（C05 等）。
- 外部参照：CFL3D Euler/RANS 数据集在 `cases/reference_data/cfl3d/`（D05–D08、F08）；
  XFOIL 参照（D13）；实验数据只记偏置不设门（D11/D12）。

## 仓库地图（活的部分）

| 路径 | 内容 |
|---|---|
| `pyfp3d/` | 求解器库：`kernels/`（numba 内核）、`mesh/`、`meshgen/`、`physics/`、`constraints/`（Kutta/Dirichlet）、`solve/`（Picard/Newton/线性解）、`viscous/`（IBL 耦合）、`post/` |
| `tests/` | A–F 六层测试 + 共享 helper（`_gate_evidence.py` 证据锁等） |
| `cases/` | `gates/`（C/D 门的已提交证据 PNG/CSV + INDEX.md）、`reference_data/`（外部真值，永不编辑）、`meshes/`、`demo/` |
| `bench/` | 按角色命名的指标/门禁脚本（见 `bench/README_ROLES.md`）+ `studies/`（只剩被活测试加载的 7 个）+ `gate_results/` |
| `docs/` | 四份纲领文档（本文件、design、design_track_v、agent-rules）；阶段工作文档随各自阶段在 `phases/pN/docs/` |
| `phases/` | phase 1–6 归档，索引见 [phases/README.md](../phases/README.md) |
| `tools/` | 外部工具链（CFL3D、XFOIL、cgrid），gitignored，不入仓库 |

## 各阶段做了什么（详见 phases/README.md）

- **phase 1**（07-06→07-25）：从零建求解器。五条 track：P 求解（Picard→Newton、Kutta、
  P14 壁压相等估计器）、M 网格、B level-set 尾迹、V 粘性起步（V1–V6）、A 校验。
- **phase 2**（07-28→08-10）：纲领重建（五条产品指标、S0–S6）；激波判据台（Laval 喷管）；
  熵修正实现；第一版能力边界。M1 如实 FAIL——可信包线实测 M_max ≲ 1.2。
- **phase 3**（08-10→08-16）：删 level-set 路线（裁决 D5，−4624 行）；hex 剖分路线判死
  （壁面误差由 O(h) 刻面几何决定）；实测 M0.80/α1.25 离散解不唯一。
- **phase 4**（08-18→08-20）：GS4.1 二维 strip IBL / 剪应力滞后（粘性重启第一块）。
- **phase 5**（08-21→08-23）：M1 二维无粘标定（网格收敛阶梯、产品指标第一次设门）。
- **phase 6**（08-24→10-12）：全面梳理 + 门重建。门号治理（556 条断言归类、tests/A..F
  重编号、证据锁）、CFL3D/XFOIL 三层对标、能力边界换代、BLAS/缓存确定论、两个归因研究
  （C05 amg、hex 族 cl −24%）。套件 595/12/2 → 831/33/2。

## 下一阶段（phase 7）候选工作项

以下全部是 phase 6 收口时**登记在册的移交项**（出处：能力边界 §9 + 收口判定 §5），
**尚未排期，开工顺序 = 使用者裁决**：

1. **hex 族 NACA 能力边界门**（「缺陷仍在」式，红=好消息）：归因已完成（hex 收敛到 tet
   从未出现的中部气泡/双激波结构，cl −24%），门建在归因之后更准；残余 = 生成器缺陷与
   熵冻结交互的归属未分开。
2. **C05 预条件子**：写文书（激波腿用 direct/ilu，amg 只适用无激波/亚声速）+ 研究
   「amg 能否修」（AMG 建在完整 J_ff 上，或 AMG 打底 + 激波方向低秩校正）。不改默认值。
3. **M3b 速度指标的锁**（需要相对基准的形状设计）。
4. **W3.1 G-TEETH 定量化**、**W3.3 ZPG 湍流平板锚**（需先有可引用关联式）、
   **W3.4 非线性 MMS**（需体积源口）、**W3.6 E01/E03 升格评估**。
5. 开放研究项：xtr030 多吸引子（可能是其他 bug 或不足引发）、D05 M0.50 cl 阶偏低、
   激波定位（S2）、翼尖自由边奇点（P13）、tet medium 在 r_far 20/30 重生成变体上不收敛、
   D06 cp_rms_lower ratio 345.9 待查。
6. 数据受限项：W3.5 翼身 CFL3D Euler 数据（M0.50+M0.84）暂时得不到。

## 文档地图

| 文档 | 职能 | 何时更新 |
|---|---|---|
| [design.md](design.md) | 理论与数值参考（方程、离散、内核规则、求解策略、验证阶梯、风险登记） | 方法/勘误变化时 |
| [design_track_v.md](design_track_v.md) | 粘性耦合线的设计、现状与边界 | Track V 推进时 |
| [agent-rules.md](agent-rules.md) | 当前阶段指针 + 工作纪律（经 CLAUDE.md `@include`） | 阶段变化时 |
| overview.md（本文件） | 人读总览 + 仓库/文档地图 | 顺手刷新 |
| `phases/p6/docs/dev_phase_six/PHASE_SIX_CAPABILITY_BOUNDARY.md` | **现行能力声明** | 下一阶段产出新边界时换代 |
| `phases/pN/docs/dev_phase_*/progress.md` | 各阶段轮次台账（历史） | 只进不改 |
| [PROJECT_STRUCTURE.md](../PROJECT_STRUCTURE.md) | 目录树与逐模块状态 | 结构变化时 |
