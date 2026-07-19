# PRE-REGISTRATION — B25 内侧自由边约束侧处理（P13/B8 rescope 到翼根）

> **日期**：2026-07-19
> **分支**：`kimi/b25-inboard-edge-taper`（本地叠在 `kimi/b24-wake-inboard-end` 上，
> 与 B23/B24 同一条未 push 链）
> **状态**：预注册，待用户裁决后实施
> **前置**：B23 VERDICT §(b)-2（本 phase 的出处）、B24 VERDICT（(b)-1 水线延伸
> 已关闭）、roadmap/track_b.md B8 close-out（翼尖约束侧处理实测死亡的 oracle）

---

## 1. 背景与出处

B23 裁决：翼身交界袋 = level-set 尾流**内侧自由边**奇点（W1），根治候选
§(b)-1（几何移除自由边）/ §(b)-2（P13 自由边处理 rescope 到内侧边）。
B24 执行 (b)-1：自由边推到远场后袋随之消失（机制二次证实），但贴面/离锥
两变体都换来同等或更糟的机尾释放奇点——**(b)-1 关闭**。

本 phase = §(b)-2 的约束侧部分：不移动片的几何，而是用 B8 为翼尖自由边
建造的**行混合焊合**机制，把内侧自由边"焊成零强度边"。

## 2. Oracle：B8 在翼尖实测了同类处理死亡（必须前置声明）

B8 close-out（track_b.md，committed 证据）：

1. TE 行混合焊合（`F·[压力 Kutta 行] + (1−F)·[φ_aux−φ_main]=0`）能完美卸载
   翼尖环量（Γ_last ~ h^4.73）、完全局部、零收敛代价——**但翼尖峰照样发散**
   （p +1.37…+1.58）。
2. 机理钉死：峰元在**几何尖外侧的 beyond_tip 未切元**里（spanwise clip 拒绝
   切），是**片终止方式（函数空间）的奇点，不是脱落环量的奇点**；
   span_blend（终止环焊合）够不到它。
3. 原话："Both constraint-side routes now measured dead; any further cure
   must change the FUNCTION SPACE at the termination."

**B24 补充证据（同向）**：袋峰跟着自由边走（延伸→机尾）——终止类奇点。

**先验结论**：约束侧处理在内侧边的期望结局是死亡（B8 机理迁移）。本实验
仍然要做，因为内侧边与翼尖有三个结构差异，迁移不是逻辑必然：

| 差异 | 翼尖（B8） | 内侧边（本 phase） |
|---|---|---|
| 环量 | Γ(tip)→0  emergent，卸载无成本 | **Γ(root) = 全翼最大**，焊合 = 直接卸根升力 |
| 峰的位置 | beyond-tip 未切元（开旷流场） | 交界 z≈z_junc，**贴着机身的流体狭条**（W2 压上蒙皮） |
| 焊合可达性 | 峰在 clip 之外，焊合够不到 | 待判别（本实验 Q2 度量直接回答） |

## 3. 实验设计

### 3.1 处理机制（全部现有机器，anchor 翻转）

现状 2 点 TE polyline：q=0 在**交界**（内侧端），q=L 在翼尖。翼尖 anchor 的
F(u) 族（`tip_taper_factors(station, z_tip, form, r_c)`，u = max(z_tip−z, 0)）
经镜像即得根部 anchor：`F_root(q) = tip_taper_factors(L−q, L, form, r)`。

- **T1 终止环焊合**（B8 span_blend 的 root rescope）：`MultivaluedOperator`
  的 `span_blend=(form, r_blend)` 元组加可选第三元 `anchor ∈ {"tip","root"}`；
  2 元组 ⇒ "tip"（**bit-identical**），"root" ⇒ 镜像 station。效果：q<r 的
  非 TE cut 节点的 wake-LS 行被 `w·[LS 行] + (1−w)·s·[φ_aux−φ_main]` 替代，
  w→0 时跳量沿 r 光滑焊到零（而不是跨一个元 Heaviside 跳零）。
- **T2 TE 焊合**（B8 tip_taper 的 root rescope）：**零库改动**——harness 直接
  算 `tip_taper = tip_taper_factors(L−q_te, L, form, r)` 传给
  `solve_multivalued_lifting(tip_taper=...)`（既有参数）。效果：q<r 的 TE
  节点压力 Kutta 行被焊合 ⇒ 根站跳量=0 ⇒ 内侧边不脱落环量。

form 固定 **vanish_linear**（B8  bracket 中唯一 s>1/2 的最小偏置正则形）。

### 3.2 腿（LS Picard M0.5、freestream+pin_gamma，b23 配方；A 侧 = committed
b23 D1 缓存，同码测量）

| 腿 | level | α | 处理 | 目的 |
|---|---|---|---|---|
| T1-c | coarse | 3.06 | span_blend root, r ∈ {0.05, 0.10, 0.20} | 终止环焊合宽度扫描（cheap） |
| T1-c0 | coarse | 0.0 | span_blend root, r=0.10 | α=0 自检（应近似惰性） |
| T1-m | medium | 3.06, 2.0 | span_blend root, r = coarse 最优者 | medium 确认 |
| T2-c | coarse | 3.06 | tip_taper root, r ∈ {0.05, 0.10, 0.20} | 环量路线宽度扫描 |

r 单位 = m（展向绝对宽度；semispan 0.85，r=0.20 ≈ 24% 展长，故意过宽以
暴露单调性）。

### 3.3 度量（复用 b24 measure_e1 + 新增）

1. 走廊 Mmax / n_sup / 峰位置 (x, z)（E1 主度量，b23 一致定义）
2. **峰元 q 站位 vs 边（q=0）**：峰在 q<0 的 beyond-edge 未切区（B8 机理：
   焊合够不到）还是 q≥0 的切元条内（焊合可达）——B8 迁移判别器
3. cl_p、γ = mean(te_jump)、**根区 γ 剖面失真**（max |Δte_jump| over
   TE 节点，A 侧对照）——Q1 成本度量
4. cl_fus 分解（band / out，W2）
5. 收敛性：n_outer / res / nlim / nflr / n_span_blended
6. n_te_nodes 不变断言（76/150，沿用）

## 4. 判定（预注册）

- **T-A（治愈，oracle 被击败）**：medium α=3.06 走廊 Mmax ≤ 1.3 且 n_sup=0，
  且 cl_p ≤ 2% A 侧、根区 γ 失真 ≤ 5%、α=0 自检惰性 → 修法在手 → demo 化 +
  GB9.4 重设流程。
- **T-B（权衡）**：袋强随 r 单调下降但 r=0.20 仍未消失，且 cl_p 成本随 r
  增长 → 出权衡曲线（袋强 vs r vs cl_p），**用户仲裁**可接受折中。
- **T-C（oracle 迁移，先验结局）**：任意 r 袋强在 A 侧 [0.7, 1.3] 带内，
  或峰元坐实 q<0 beyond-edge → 约束侧在内侧边同样死亡 → **(b)-2 约束侧
  部分关闭**；残余根治路线 = 函数空间/终止重构（B10 free-wake roll-up，
  目前 SHELVED——独立的重大路线决策，不在本 phase）；交付维持 B23 §(a)
  遮蔽规则。

## 5. 风险登记

| # | 风险 | 概率 | 影响 | 缓解 |
|---|---|---|---|---|
| Q1 | 根部 Γ 最大，焊合一阶卸升力（不同于翼尖零成本） | 高 | 高 | 度量 3 显式报 cl_p/根区 γ；判定 T-B 权衡 |
| Q2 | 峰在 beyond-edge 未切区，焊合够不到（B8 机理迁移） | 高 | 中 | 度量 2 直接判别；判定 T-C |
| Q3 | 焊合行与 pin_gamma / 远场 aux 交互（B16 类近奇异） | 中 | 中 | 监控 res/n_outer + n_span_blended；发病则缩 r 或如实记 |
| Q4 | anchor 翻转破默认路径 | 低 | 高 | 2 元组 ⇒ tip bit-identical；单测锁定；B 套件全绿才跑 |
| Q5 | α=0 焊合引入污染 | 低 | 中 | T1-c0 自检腿（跳量本 ≈0，应近似惰性；判据 cl_p \|Δ\| ≤ 1%） |
| Q6 | 焊合改变 implicit Kutta 全局平衡，γ 重调（B8：sheet 侧 δ-pin 曾把全局 Γ 重调 ~10×） | 中 | 中 | mean-γ 与根区剖面分别报，不混用 |

## 6. 工程纪律

- 线程：NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16。
- 改动面 = `pyfp3d/wake/multivalued.py`（span_blend 第三元 anchor，2 元组
  bit-identical）+ 本目录 harness（T2 零库改动）。不动 kernel/assembly、
  不动网格、不动 B10 shelved 机器。
- 单测：anchor 翻转值正确（镜像）、默认 bit-identical、n_span_blended 计数；
  跑前 test_b8_span_blend / b1 / b2 / m2 / v0 全绿。
- 不改任何 committed A 侧缓存；B 侧新缓存独立 tag；结论以 committed CSV/PNG 为准。

## 7. 时间与产出

- 实现（anchor + harness + 单测）：~0.5 天
- 求解：coarse ~1 min × 8 + medium ~15 min × 3 ≈ 1 h
- 分析 + VERDICT：~0.5 天
- 产出：本文件、`wb25.py`、`run_t1.py`、`results/t1_summary.csv`、
  `results/t1_trade.png`、`VERDICT.md`
