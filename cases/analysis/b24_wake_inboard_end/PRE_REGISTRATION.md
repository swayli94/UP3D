# B24 尾流内侧自由端修相 — 预注册（PRE-REGISTRATION）

> **日期**：2026-07-19（在任何测量运行之前写就）
> **分支**：`kimi/b24-wake-inboard-end`（叠在 `kimi/wingbody-junction-discriminator`
> 之上——B23 未合入，本相位的控制侧数据全部来自 B23 已提交证据）
> **执行令**：B23 裁决（`cases/analysis/b23_junction_discriminator/VERDICT.md`）——
> "W1 袋 = level-set 尾流片内侧自由边奇点（P13 自由边类）；根治 = Track B 尾流
> 模型层的'内侧自由端'处理，首选 TE polyline 内端沿机身水线延伸，判据 =
> 袋消失/移位 + 翼面载荷不动。"
> **用户批准**：2026-07-19（"要的"）

## 假说（B23 证据链，一句话）

level-set 尾流片的 TE polyline 在机身水线 z_junc 处戛然而止 → 内侧自由边 =
终止于机身旁的涡片边缘 → 交叉流速度奇点（与翼尖自由边同族）沿机身水线印在
流场里 → 虚假超声速袋 + 机身袋带升力污染。**让近场不存在自由边**（polyline
内端沿机身水线延伸、自由边推到远场）应使袋消失或移位到远场，且翼面载荷不动。

## 机制界面核查（写于跑前，代码只读确认）

1. **TE 节点/γ 不被延伸污染**：`CutElementMap` 的 TE 节点 = (|s|<tol ∧ |d|<tol)
   ∩ `wall_nodes`，而 `wall_nodes` 只含"wall"组（翼面蒙皮）——延伸段水线节点在
   "fuselage"组，被过滤排除（`pyfp3d/wake/cut_elements.py:127-134`）。
   γ = `mean(mvop.te_jump)`（标量，`picard_ls.py:766`）key 在 `cm.te_nodes`
   （`multivalued.py:368`）→ 延伸后 TE 集合与 γ 定义不变。实验里**验证
   n_te_nodes 不变**（coarse 76 / medium 150，M2 锁定值）。
2. **polyline 多段支持**：`WakeLevelSet` 原生多段（`levelset.py:138-140` 明确
   支持 curved/kinked TE per-panel）。延伸只是更长的 polyline，不动 levelset/
   CutElementMap/MultivaluedOperator/solver 任何代码。
3. **延伸片不穿机身**：水线段 (x, 0, R(x)) 沿 d=(cosα,sinα,0) 扫出的点
   半径 √(R(x)²+t²sin²α) ≥ R(x) ≥ R(x+t·cosα)（等半径段与收缩段均成立）——
   片恒在机身外。片对机身表面为二次方抬升（贴面擦行）→ R1 风险在案。
4. **远场 aux 有处置**：延伸细条（z≈R_tail）把更多 sheet aux DOF 送到远场边界，
   `farfield_aux="pin_gamma"`（B17，现默认）正是为 outflow 跳量=γ 而设；
   "TE aux 落远场→RuntimeError"护栏不触发（TE 集合不变，见 1）。
5. **2×2 条件数**：水线段与 d 夹角≈α，`det = L²sin²α`（cond ~1/sin²α ≈ 350
   @ α=3.06°）——双精度下可接受；`update_direction` 的 1e-12 guard 不触发。
6. **零新网格**：level set 是求解侧几何——A/B 复用现有 coarse/medium 网格；
   **A 侧（对照）数据已全部提交**（B23 D1/D2b/W2），本相位只需跑 B 侧。

## 实验设计

### 几何变体（B 侧 = 处理）

`te_polyline(p)` 加可选延伸（默认 None = bit-identical，现状 2 点线）：

- **B1「贴水线」**：polyline 从远场端 (x_far, 0, R_tail) → 机尾 (x_tail, 0,
  R_tail) → 沿水线 (x, 0, R(x)) 上行至交界 TE 点 (x_te(z_junc), 0, z_junc)
  → 现有翼 TE 至翼尖。水线采样间距 ≈ h_body。q=0 落在远场端 →
  **内侧自由边移出近场**。片覆盖翼后 z∈[z_junc, B]（现状）+ 机身后
  z∈[R_tail, z_junc]（新）+ z≈R_tail 细条至远场（新）。
- **B3「离锥」**（仅当 R1 触发时启用）：同 B1 但 z = R(x)+δ（δ ≈ 1×h_body），
  片离开机身蒙皮，消除贴面擦行；拓扑等价，物理近似稍粗。

x_far 取远场边界内侧（与现有翼 sheet 的 outflow 行为一致）；精确值在实现时
按网格远场半径记录于结果 CSV。

**排除的变体（跑前记录）**：z≡z_junc 直线延伸**不**覆盖 z<z_junc 内岸区，
自由边仍在原处——拓扑上无效，不作为腿。

### 求解腿（同一配方，LS Picard M0.5，freestream + pin_gamma）

| 腿 | 网格 | α | 说明 |
|---|---|---|---|
| A 侧 | medium/coarse | 0,1,2,3.06 / 0,3.06 | **已完成**（B23 D1 committed） |
| B1-medium | medium | 0, 2, 3.06 | α=1 跳过（D1 已示过渡） |
| B1-coarse | coarse | 0, 3.06 | 加密趋势锚 |
| B3-* | 同上 | 同上 | 仅 R1 触发时 |

α=0 腿是**自检**：延伸片不得在零升力下引入新污染源（对称性）。
求解缓存 `results/*.npz`（gitignored）；committed 证据 = CSV/PNG。

### 度量（每条腿；复用 b23 的 measure / 区域拆分 / W2 分解核心）

1. 交界走廊 Mmax（z<0.5 ∧ x>0.8）+ sheet 走廊（|s|<0.03）Mmax —— **主度量**
2. n_sup（全场 + 走廊）；袋位置 top-200（q, |s|, x, z）
3. cl_p、γ（收敛值）
4. cl_fus 分解：袋带（|z−z_junc|<0.06, x>1）/ 带外 / 极点（W2 核心原样复用）
5. 收敛性：res、n_outer、nlim/nflr
6. 翼尖区 Mmax（z>0.8）—— 翼尖护栏（延伸不得影响翼尖奇点类）
7. 远场 sheet-aux DOF 计数；n_te_nodes（必须 = 76/150）
8. 袋峰 x 位置 vs x_tail（=2.42）—— "移位到远场" 判据用

## 判定（预注册）

- **E1（主）**：B 侧 medium α=3.06 交界走廊 Mmax ≤ 1.3 且走廊 n_sup = 0
  → **袋消失**；或袋峰移出 x > x_tail（远场）→ **袋移位**。任一 = 机制确认。
- **E2（护栏）**：cl_p 相对 A 侧 |Δ| ≤ 2%；|Δγ| ≤ 5%；翼尖 Mmax 比值 ∈ [0.5, 2]
  （实现混沌容差，B23-D2b 已标定该类散布）；α=0 无袋（走廊 Mmax ≤ 1.3）；
  n_outer ≤ 1.5× A 侧；nlim+nflr ≤ A 侧 + 10；n_te_nodes 不变。
- **E3（W2）**：袋带 cl_fus 分量塌缩（|band| ≤ 0.3× A 侧同态）；带外分量
  |Δ| ≤ 20%（物理 carryover 不动）。

### 裁决树

1. **E1 ✓ ∧ E2 ✓** → 机制确认 + 修法在手 → 进入收尾流程：demo 化
   （cases/demo/b24_*）、GB9.4 重设落地、GB20.5/B18 重锚、Track B 文档。
2. **E1 ✓ ∧ E2 ✗** → 机制确认但修法有载荷副作用 → 记录偏差来源（γ/几何），
   判"有效但需二期修载荷"，**不得**直接改 gate。
3. **袋移位到近场新位置（x_te < x ≤ x_tail 仍有走廊袋）** → "自由边迁移"——
   记录新位置，启用 B3；B3 仍迁移 → 记"延伸类不足"，退回 P13 rescope 路线
   （B23 VERDICT §b-2）。
4. **E1 全 ✗（袋仍在原处同量级）** → 自由边假说**被证伪**——记"不可判别"，
   B23 裁决需翻案，重新归因（重大发现，单独报告）。
5. R1 触发（贴面 cut 元导致不收敛/残差爆发 > A 侧 10×）→ 记 R1，转 B3 腿。

## 风险登记

| # | 风险 | 概率 | 影响 | 缓解 |
|---|------|------|------|------|
| R1 | 片贴机身二次方抬升 → 壁面擦边 cut 元/sliver，残差行为恶化 | 高 | 中 | B3 离锥变体；判定树 5 |
| R2 | 延伸 cut 元集合的 aux 块条件数（B16 类近奇异） | 中 | 高 | pin_gamma 已是 B16 的修；度量 7 监控；必要时 extent 裁剪细条 |
| R3 | TE/γ 定义被延伸污染 | **已排除**（机制核查 1） | — | n_te_nodes 不变作为运行期断言 |
| R4 | z≈R_tail 远场细条网格欠分辨（h_wake 场未覆盖） | 中 | 低 | 首轮接受（对流跳出）；若 R2 发病再 meshgen 加 wake field（独立腿） |
| R5 | 水线段与 d 近平行 → evaluate 2×2 病态 | 低 | 中 | cond ~350 可接受；若报错，延伸段实现侧显式处理 |
| R6 | α=0 延伸片引入新污染源 | 低 | 高 | α=0 自检腿（判定树 E2 条款） |
| R7 | B1 物理性：贴水线片 = 机身水线涡脱落的粗模型，可能过拟合袋而动载荷 | 中 | 中 | E2 载荷护栏；Track V 前的临时模型地位在收尾文档里写明 |

## 工程纪律

- 线程：NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16。
- `te_polyline` 默认路径 bit-identical（m2 测试 23/23 必须保持绿 + 现状
  2 点线返回值不变断言）。
- 不改任何 kernel/assembly；改动面 = `meshgen/wingbody.py: te_polyline` 可选参数
  + 本目录 harness（复用 b23 wb_common，extend 参数透传）。
- 首跑前确认无陈旧 npz；日志不得出现 silent reuse。
- 所有数字以 committed CSV 为准。

## 时间与产出

- 实现（te_polyline 延伸 + harness 透传 + 自检）：~0.5 天
- 求解：B1-medium ×3（~15 min each）+ B1-coarse ×2（~3 min each）≈ 1 h（+ B3 备用）
- 分析 + VERDICT + 收尾建议：~0.5 天
- 产出：`PRE_REGISTRATION.md`（本文件）、`run_e1.py`、`results/e1_summary.csv`、
  `results/e1_pocket_map.png`、`results/e1_w2.csv`、`VERDICT.md`

---

## 附录 A（2026-07-19，首轮执行后补记）：R5 发病 + 新风险 R8 与 levelset 修复

预注册后首轮执行暴露两个几何问题，均在"改动面"承诺内处理，如实登记：

**R5 已发病（低概率 → 确定）**：α=0 时延伸段全在 y=0 平面沿 x 向，与
d=(1,0,0) 精确平行，`update_direction` 的 1e-12·L guard 触发。
处理（预注册 R5 缓解条款"延伸段实现侧显式处理"）：延伸段加 0.2° 下倾
（`tilt_deg=0.2`），y = −tan(0.2°)·(x−x_j)，机身处位移 ≤2.8e-3 ≪ h_body，
2×2 cond ~8e4 → ε 误差 ~1e-11 ≪ 1e-6·h 容差；交点 y=0 与翼 TE 自动连续。

**R8（新风险，预注册未预见）：corner 偷窃**。α>0 时水线片（z≈0.15 幕帘
以 sinα 上升）的**后向延长面**穿过翼根后缘下游下方区域，`evaluate()` 的
面板选择度量 s²+excess² 不含 d，把该区域节点按"更近平面"判给水线面板
（d≈−0.17，翼根后缘 TE 节点随之失去全部 d>0 crossing → 无 cut 元 →
CutElementMap A3 断言爆炸，B1 medium α=2 首轮即挂）。coarse 复现：
missing=2（正交界 + 第一个外侧 TE 节点），α=0 无恙。
处理：**`WakeLevelSet.evaluate`/`surface_normals` 面板选择度量加
min(0,d)² 后退罚项**——垂直落点在片后（d<0）说明该点真正最近的是
面板 TE 边而非片本身，罚项让"片真的在"的面板胜出。单段 TE 上 argmin
退化为唯一面板 ⇒ committed 行为 bit-identical（b1/b2/v0/m2 82 绿 +
LS 套件 65 绿验证）；新增回归测试
`test_panel_selection_prefers_on_sheet_foot`（修复前 d=−0.159 失败，
修复后过）。修复后 coarse/medium α∈{0,2,3.06} missing=0，n_te 不变。
B3 离锥（δ=0.06 coarse 试）**加重** R8（抬升幕帘 z 向侵入翼根展向范围，
missing 增至 4–7），B3 腿保留但首选 B1。

改动面追加：`pyfp3d/wake/levelset.py` 面板选择度量（非 kernel/assembly，
单段 bit-identical）+ `tests/test_b1_cut_elements.py` 回归测试。
