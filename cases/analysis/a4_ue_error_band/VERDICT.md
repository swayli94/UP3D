# A4 — Wall edge-velocity (u_e) error-band study — VERDICT

> **日期**：2026-07-22 · **性质**：Track-V 输入质量前置测量（A-track 小项，非 gate）·
> **触发**：审计 `20260722-0335-b28-b32-audit-pre-trackv.md` §5.2.1 + 2026-07-20
> wingbody-trackv-review §3.3 风险 1（IBL 消费无粘壁面 u_e，必须先量化其 O(h) 误差带，
> 否则 V-gate 会把"粘性模型误差"与"无粘输入误差"混在一起）。
> **方法**：在两个有**解析壁面速度真值**的 M0 算例上直接测 u_e = |∇_t φ| 误差：
> 圆柱 2.5-D（u_e = 2U sinθ）与球壳（u_e = 1.5U sinθ），精确远场 Dirichlet + 自然壁面 BC
> （同 `test_m0_cylinder` / `test_laplace_sphere` 的设置），coarse+medium × 线性/二次两种
> 壁面恢复。真值解析 ⇒ 无 chasing 风险，无需预注册 A/B。

## 1. 结论（一句话）

**中网格光滑壁面 IBL 输入 u_e 误差带 ≈ 峰值处 2.5% 相对 / 0.04·U∞ max-norm / 0.012·U∞
rms，O(h) 可改善；在 LE/驻点带（u_e→0、du_e/ds 最大 = IBL 初值最敏感处）相对误差最大
（medium 4–7%，coarse 6–12%）。恢复方案（线性 vs 二次）无普适赢家，按区/算例差 ~1%。**
锐 TE 的二次恢复 guard 仅在楔形角 <~6° 触发——**NACA0012 的 16° TE 不触发**（审计原文
"锐 TE 只能用线性"过强，据此更正）。

## 2. 测量（`results/ue_bands.csv`，quadratic overall）

| 算例 | coarse abs_max | medium abs_max | O(h)? | medium rel（峰值尺度） | medium abs_rms |
|---|---|---|---|---|---|
| 圆柱 2.5-D | 0.0842 | 0.0451 | ✓（≈×0.5） | 2.3 % | 0.0128 |
| 球壳（3-D） | 0.0805 | 0.0392 | ✓（≈×0.5） | 2.6 % | 0.0114 |

- **O(h) 确认**：coarse→medium（h 减半）误差减半——这是 P1 壁面梯度的固有地板，与 V6 地板 /
  G1.6 重归因同类（P1 场在 h≈0.04–0.08 的 max-norm 能力）。
- **分区结构（IBL 关键）**：
  - **LE/驻点带**（θ∈[0,20]∪[160,180]，u_e→0）：相对误差最大——圆柱 quad 12.3%(coarse)/
    6.6%(medium)、球 5.4%/4.2%，因 O(h) 绝对误差 ÷ 小 u_e。**这正是 IBL 布初值、需要
    du_e/ds 的带，也是 u_e 最不可信的带。**
  - **肩/峰带**（θ∈[70,110]，u_e 最大）：最准——圆柱 quad rel 1.1%/0.5%、球 5.4%/2.6%。
- **恢复方案权衡**（无普适赢家）：圆柱 LE 带**线性远好于二次**（medium 0.020 vs 0.045）；
  圆柱峰带**二次好于线性**；球壳二次全域 ~30% 更好。⇒ IBL 实现应按区选，或 LE 带用线性。

## 3. TE 结构约束（`results/te_constraint.csv`，更正审计）

| 算例 | TE 楔形角 | 线性恢复可用 | 二次恢复可用 |
|---|---|---|---|
| NACA0012 coarse | 16.3° | ✓ | ✓ |
| NACA0012 medium | 16.4° | ✓ | ✓ |

`wall_tangential_gradient_quadratic` 的 `_wall_vertex_normals` guard（|Σ area·n̂|/Σarea<0.05
≈ 楔形 <~6°）在 NACA0012 的 16° TE **不触发**。故二次恢复在 NACA0012 上可用；guard 只对
**薄/尖点 TE（<~6°，如超临界/RAE 类翼型或 M6 foilmod 极锐 TE）**才成为天花板——后续 GV
gate 若用薄 TE 翼型需重测。

## 4. Track V 预算建议

1. **GV1.x/GV3.x 必须把这条 ~2.5%（medium，峰值）无粘输入误差带单列**，与粘性模型误差分开
   归因。GV3.3 的"CL 下移 ≈0.02 向实验"方向检验（0.02/0.27 ≈ 7%）**安全高于**此地板 ⇒ 该
   gate 成立；但 LE/驻点带的**紧** δ*/u_e 逐点对标会**输入受限**（那里相对误差 4–7% @ medium）。
2. **du_e/ds（IBL 初值）在 LE/驻点带最不可信**——V1 布 IBL 种子点应避开或对该带的 u_e 用
   线性恢复 + 光滑（`smooth_wall_tangential_gradients`）。
3. **恢复方案按区选**（LE 带线性、峰带二次）可省 ~1% 输入误差——V1 数据布局第一设计点。
4. **TE guard 不阻塞 NACA0012**（16° > 6°）；薄 TE 翼型的 GV gate 需另行核实二次恢复可用性。

## 5. 复算

`python cases/analysis/a4_ue_error_band/run.py`（~30 s，M0 Laplace 解便宜，无缓存）。
证据：`results/ue_bands.csv` / `te_constraint.csv` / `ue_error_band.png`。
