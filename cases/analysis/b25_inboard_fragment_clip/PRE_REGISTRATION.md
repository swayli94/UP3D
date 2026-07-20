# PRE-REGISTRATION — B25 尾流内侧 fragment 裁切（LS 片 = conforming 拓扑）

> **日期**：2026-07-19（v2 重写；v2.1 = Claude 审计修订，见附录 A）
> **分支**：`kimi/b25-inboard-fragment-clip`（本地叠在 B24 分支上，未 push 链）
> **状态**：预注册 v2.1（含附录 A 审计修订），实施中
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

---

## 附录 A — v2.1 审计修订（2026-07-19，Claude 审计，用户转发批准）

实施前对 §2.3/§3 的三处判据修订。v2 原文保留未改，本条与原文冲突时
以本条为准。

**A1（修订 §3 C-A 的 cl_p 护栏，方向条款）**：cl_p 相对 A 侧 |Δ| ≤ 2%，
**或**位于 A 侧与 conforming oracle（medium 0.2173，B18 记录）之间——
向 oracle 移动 = 记录性通过（袋污染若压低根区升力，根治后 cl_p 上移
是最佳结果，+2.6% 的跨模型残差不应误报失败）。|Δγ| ≤ 5% 同加方向
条款（向 conforming γ 移动从宽，背离从严）。

**A2（修订 §2.2 α=0 自检与 §3 C-D，绝对阈）**：α=0 的 cl_p 在 ~5e-5
量级，相对判据噪声主导（B24 实测 A 5.31e-5 vs B1 5.39e-5 = 1.5%
相对差，完全惰性）。改为**绝对判据**：|Δcl_p| ≤ 1e-4 且走廊无袋
（corr Mmax 不升、峰位点逐侧一致）才算惰性。——冒烟实测已确认此
修订必要：coarse α=0 d_cl_p = +5e-5（相对 +10.1%，按 v2 原文会误触
C-D；按 A2 绝对阈惰性成立）。

**A3（拓宽 §3 C-B + 新增度量 8，条带跳量锚定病）**：C-B 出口从仅
"贴壁 sliver 病"拓宽为"贴壁 sliver 病 **或 aux 锚定病**"。物理：
fragment clip 新增的条带（尾锥区 + 机尾后对称面条）的 wake-LS 对流
特征线上溯终止于机身表面——既不是 TE（Kutta 锚）也不是远场
（pin_gamma 锚），对流闭合在那里没有流入数据；conforming 无对流
闭合故存在性证明不延伸到此点。B16 同类病（无锚 aux → cond1
O(1e19)、Picard 吸收近奇异行后"收敛但跳量垃圾"coarse |jump| 53.4）
是前科。新增：

> **度量 8（条带跳量锚定监视）**：条带 cut 元 aux 节点的
> |jump| = |phi_u - phi_l| 的 max / p95，及 max|jump| / |γ| 比值
> （B16 诊断手法原样复用）。判读：O(γ) 量级为健康；出现
> |jump| >> |γ|（量级偏离）即 aux 锚定病 → C-B 出口。

**A4（措辞修订 §1，非阻塞）**："拓扑与 conforming 逐点同构"降为
"**拓扑等价**"——α>0 时片与机身的真实交线不在 y=0 水线上
（conforming 由 OCC fragment 算真交线，LS 用解析旋成面裁切），离散
上 h 尺度等价。同节几何预判：直纹面在等半径段基本潜入机身内部、
尾锥区才钻出，新增 cut 元集中在**尾锥 + 机尾后**——census（度量 6）
读数时"等半径段无新 cut 元"是预期，不是 clip 未生效。冒烟实测符合：
coarse 新增条带 ~500（body，尾锥带 z∈[R(x), r_f)）+ ~100–350
（sym，机尾后 z∈[0, r_f)），等半径段为零。

**A5（修订 §7 后续路线，非阻塞）**："若 C-A 成立 → LS 翼身全包线
可信"超前——A/C 只跑 Picard M0.5。C-A 成立后的**具名后续实验 =
重跑 B18 LS Newton freeze-ramp**（现天花板 ~M0.55 是带着袋测的；若
袋是限制器，根治后天花板应抬升），全包线结论等它落地。
