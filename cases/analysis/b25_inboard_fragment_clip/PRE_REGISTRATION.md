# PRE-REGISTRATION — B25 尾流内侧 fragment 裁切（LS 片 = conforming 拓扑）

> **日期**：2026-07-19（v2 重写）
> **分支**：`kimi/b25-inboard-fragment-clip`（本地叠在 B24 分支上，未 push 链）
> **状态**：预注册，待用户裁决后实施
> **前置**：B23 VERDICT（袋 = 尾流内侧自由边奇点）、B24 VERDICT（(b)-1 直纹
> 延伸关闭）、`pyfp3d/meshgen/wingbody.py` docstring（conforming 片拓扑原文）、
> `pyfp3d/mesh/wake_cut.py` docstring（自由边/边界复制规则原文）

---

## 0. 作废记录（v1 焊合版，不采用）

v1（`b25_inboard_edge_taper`，commit a72224f 后重写）拟把 B8 翼尖行混合焊合
rescope 到翼根（T1 终止环焊合 / T2 TE 焊合）。用户质询后裁决**不采用**：

1. **不物理**：根部 Γ 全翼最大，焊合 = 断言根区不脱落涡量 = 拿真实根升力
   换数值正则化（v1 预注册 Q1 自报的第一条）。
2. **B8 oracle**：翼尖同款焊合完美卸载 Γ(tip) 后峰照样发散——奇点住在片
   的终止方式（函数空间），不住在环量里；约束侧处理先验无效。
3. **答非所问**：袋的病是片拓扑不完整（涡在流体内部终止，违反
   Helmholtz），焊合不补拓扑，只把"涡去哪儿"偷换成"涡不存在"。

## 1. 物理推导链（为什么是这一条）

连续模型要求：涡面的每条涡线有去处（贴物面 / 出域边界 / 闭合）。
**conforming 片满足它的方式**（wingbody.py docstring 原文）：

> the sheet is extended inboard to below the symmetry plane (z_lo) and
> downstream through the far-field sphere, then fragment+embed'ed. The
> fragment trims it automatically to: exposed wing TE -> fuselage waterline
> (the y=0 top meridian, z=R(x), junction TE -> tail tip) -> symmetry edge
> (z=0 aft of the body) -> sphere arc -> free tip edge.

即：conforming 片 = **尾流面被机身裁切（fragment）**——片的内侧边界是
尾流面与机身表面的**交线**（水线→尾锥→对称面→远场球），侧边界永远落在
物面或域边界上，从不在流体内部。机身侧面的片就是尾流面本身（不是从水线
扫出的另一块面）。

**B24 的认错对象**：直纹延伸把水线当成片的**前缘**（d=0），从水线扫出
幕帘——物理上水线是片的**侧边界**（翼根涡线对流的轨迹），不是脱落缘。
幕帘在机尾脱体释放，根子在此。

**本 phase 的 LS 等价物**：LS 尾流面（直纹面，随 α 重定向）不动；
把 CutElementMap 切元判定的**内侧 clip 从"交界站 q≥0"换成机身裁切**：

```
on_sheet = (d_cross > 0) & (q_cross <= span_length) & inboard_ok(x)
inboard_ok(x,y,z):
    x <= x_end_body :  y^2 + z^2 >= R(x)^2   # 机身外侧（迹线=水线，在物面上）
    x >  x_end_body :  z >= 0                # 机尾后：到对称面（域边界）
```

- 片在机身旁向内延伸直到撞上机身表面（迹线在物面 ✓），机尾后延伸到
  对称面 z=0（域边界 ✓），下游出远场（现状已有，B16）✓——片的边界
  全部落在物面/域边界上，**拓扑与 conforming 逐点同构**，自由边不再
  存在于流体内部（Helmholtz 满足）。
- s, d, q 语义、直纹面、法向、TE/Kutta/γ 定义**全部不变**（前缘仍只有
  翼面 TE；水线段不进 polyline）；B24 的 R8 修复、R5 tilt、延伸 polyline
  设施本路线**不需要**（单段 polyline 保持）。
- 翼尖 clip（q ≤ span_length）不变——翼尖自由边两条路径共存，不在本
  phase 范围。

## 2. 实验设计

### 2.1 改动面（预注册承诺）

| 文件 | 改动 |
|---|---|
| `pyfp3d/wake/cut_elements.py` | `on_sheet` 的可选 `inboard_clip` 参数：`None`（默认）= 现状 q≥0，**bit-identical**；传入 callable 时替换内侧 clip。q≥0 逻辑、tip clip、d 逻辑不动 |
| `pyfp3d/meshgen/fuselage.py` 或 harness 内 | 机身裁切函数 `inboard_ok`（`radius_at` 已有；x_end_body = 机身末端，实现时确认 `radius_at` 域外行为） |
| `tests/test_b1_cut_elements.py` | 单测：默认 None bit-identical；机身内点拒绝（y²+z²<R²）；机身旁 z>R 接受；机尾后 z≥0 接受、z<0 拒绝 |
| `cases/analysis/b25_inboard_fragment_clip/` | harness（wb25.py 复用 wb24/wb_common + run_f1.py） |

不动：kernel/assembly 求解层、levelset.py、网格、B10 shelved 机器、
committed A 侧缓存。

### 2.2 腿（LS Picard M0.5、freestream+pin_gamma；b23 配方）

**本路线无旋钮**（拓扑修正，无宽度/形式参数）——干净的 A/C 对照：

| 腿 | level | α | 说明 |
|---|---|---|---|
| A | coarse/medium | 0, 2, 3.06 | 对照 = committed b23 D1 缓存（q≥0 clip），同码测量 |
| C | coarse | 0, 2, 3.06 | fragment clip |
| C | medium | 0, 2, 3.06 | fragment clip |

α=0 腿 = 自检（新 cut 集在机身旁/对称面出现，但 α=0 时 Γ≈0，物理上
应近似惰性；判据 cl_p |Δ| ≤ 1%、无袋、收敛正常）。

### 2.3 度量（复用 b24 measure_e1 + 新增）

1. 走廊 Mmax / n_sup / 峰位置 (x,z)（E1 主度量，b23 一致定义）
2. 峰元 q 站位（vs q=0 交界边 / vs 机身迹线）
3. cl_p、γ=mean(te_jump)、根区 te_jump 剖面失真（A 侧对照）
4. cl_fus 分解（band/out，W2）
5. 收敛：n_outer/res/nlim/nflr
6. 新增 census：n_cut（含机身旁/对称面条带分项）、n_aux on symmetry 边界、
   n_aux on farfield、n_te_nodes 断言（76/150 不变）
7. 条带 sliver 指标：机身旁 cut 元的最小二面角/体积分位数（贴壁病态监测）

## 3. 判定（预注册）

- **C-A（治愈）**：medium α=3.06 走廊 Mmax ≤ 1.3 且 n_sup=0；cl_p ≤ 2%、
  |Δγ| ≤ 5%、根区剖面失真 ≤ 5%；α=0 惰性；收敛正常（n_outer ≤ 1.5× A，
  res ≤ 1e-7，无精确奇异）；带外 carryover |Δ| ≤ 20% → **拓扑修复成立**
  → demo 化 + GB9.4/GB20.5 重设 + Track B 文档收尾。
- **C-B（迁移/贴壁发病）**：袋消失但条带与机身迹线擦行处出现新峰或
  收敛病态（sliver 指标同步恶化）→ 记"拓扑对、贴壁实现层发病"，
  转 sliver 处理腿（cut 容差/局部策略，独立预注册）。
- **C-C（依旧）**：袋在原处同量级 → 自由边假说第三次被挑战（概率极低，
  B23/B24 已双重证实）→ 重新归因，单独报告。
- **C-D（α=0 污染）**：自检腿非惰性（|Δcl_p| > 1% 或新袋）→ 记 R 系
  发病，判"裁切引入新污染源"，回滚分析。

## 4. 风险登记

| # | 风险 | 概率 | 影响 | 缓解 |
|---|---|---|---|---|
| S1 | 迹线擦行 sliver cut 元（片贴机身水面相交，α>0 迹线移到 y>0 侧） | 高 | 中 | 度量 7 监测；C-B 出口；conforming fragment 是同病网格版（已共存） |
| S2 | 对称面到达段 aux DOFs（新边界 aux 类别） | 中 | 中 | B16/B17 远场 aux 同类先例；census 度量 6；必要时 BC 侧重定向适配（wake_cut.py 规则：边界面重指向 + 侧） |
| S3 | 分支切割单连通性（conforming 注释：片不到对称面则 branch cut 失败） | 中 | 高 | 本设计机尾后片到 z=0 对称面，拓扑满足；若 aux 块全局奇异 → B16 pin_gamma 机制 |
| S4 | 裁切函数仅支持解析机身（radius_at），一般几何（NURBS）不适用 | 确定 | 低 | 本 phase 只声明 wingbody 解析机身；一般化留文档注记 |
| S5 | 机身旁条的跳量强度由 wake LS 条件对流决定，与机身 wall BC 交互（根涡贴着壁面对流） | 中 | 中 | E1/E3 度量监测 band/out；物理上这正是 conforming 的对应物，对照两路径 |
| S6 | 默认路径回归 | 低 | 高 | None ⇒ bit-identical + b1/b2/m2/v0 全绿才跑 |

## 5. 工程纪律

- 线程：NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16。
- 默认 None ⇒ 现状 bit-identical（b1/b2/m2/v0 + B 套件绿是开跑前置）。
- 结论以 committed CSV/PNG 为准；npz/msh gitignored；日志本地。
- 不 push；push/PR 先获用户确认。

## 6. 时间与产出

- 实现（clip 参数 + 裁切函数 + 单测 + harness）：~0.5–1 天
- 求解：coarse ~1 min × 3 + medium ~15 min × 3 ≈ 1 h（A 侧全缓存）
- 分析 + VERDICT：~0.5 天
- 产出：本文件、`wb25.py`、`run_f1.py`、`results/f1_summary.csv`、
  `results/f1_pocket.png`、`VERDICT.md`

## 7. 与其他路线的关系（备忘）

- 若 C-A 成立：LS 路径翼身全包线可信 → conforming 退居交叉验证；
  Track V（IBL 耦合）的片拓扑前置就位。
- 若 C-B/C-C：残余路线 = B10 roll-up 类（函数空间，SHELVED，独立重大
  决策）；交付维持 B23 §(a) 遮蔽 + conforming。
- B24 设施去向：`te_polyline(extend=)` 保留入库（本路线不消费；Track V
  或需）；R8 levelset 修复是多段 TE 的永久正确性条件，保留。
