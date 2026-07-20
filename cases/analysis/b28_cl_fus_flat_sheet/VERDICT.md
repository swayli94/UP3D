# VERDICT — B28 cl_fus 带外 flat-vs-tilted 解耦 + GB9.4 重设

> **日期**：2026-07-20
> **分支**：`kimi/b28-cl-fus-gate-respec`（基于 origin/main 0f27900）
> **预注册**：`PRE_REGISTRATION.md`（跑前写就；本文件按其判定树执行）
> **判定**：**F1（位置唯一因子）成立** —— C-vs-oracle 的带外差由尾流片
> **位置**（平 y=0 vs 斜随 α）唯一决定；片位置匹配后 LS 与 conforming
> 的机身升力在 15% 容差带内一致（实测 7.3%）。GB9.4 按预注册 §3.4-F1
> 重设落地（B9 demo LS 腿改 flat-fragment + 带外跨模型一致硬 gate）。
> **证据**：`results/f2_summary.csv`、`results/f2_decomposition.png`、
> `results/w2_conf_coarse.csv`（均由 `run_f2.py` / `run_w2_conf_coarse.py`
> 重跑自缓存重新生成）；求解缓存 `results/*.npz`（gitignored 同 B25）。

---

## 1. 三腿结果（F = fragment clip + 平片 sheet_direction=(1,0,0)，唯一变量 = 片位置）

| 腿 | conv | n_outer | res | cl_p | cl_fus | band | **out** | poles |
|---|---|---|---|---|---|---|---|---|
| F medium α=3.06（决定性） | ✓ strict | 60 | 9.0e-8 | 0.21512 | 0.03216 | −0.00043 | **0.03259** | 0.00065 |
| F medium α=2.0 | ✓ strict | 66 | 6.5e-8 | 0.14060 | 0.02104 | −0.00025 | **0.02129** | 0.00043 |
| F coarse α=3.06 | ✗（注 1） | 220 | 6.1e-7 | 0.20593 | 0.02893 | −0.00045 | **0.02938** | 0.00064 |

对照锚（committed，只读）：C medium α=3.06 out **0.05035**；conforming
oracle out **0.03514**（medium，`b25/w2_conf.csv`）/ **0.03456**（coarse，
本 phase `w2_conf_coarse.csv` 补锚，加密不敏感 ✓）；A out 0.02138。

## 2. 预注册判定逐项核对（§3.3 判定树，决定性腿）

- **|F_out − oracle| ≤ 0.15·|oracle|**：0.03259 vs 0.03514 = **7.25%** ✓
  —— 且落在带外量自身加密噪声（8–10%）之内，即**测量分辨极限内的相等**。
- **|F_out − C_out| > 0.15·|C_out|**：0.03259 vs 0.05035 = **35.3%** ✓。
- **位置因子 r = |C_out − F_out| / |C_out − oracle| = 1.17**：位置解释了
  全部间隙并轻微过冲（过冲量在噪声带内）。
- **趋势腿同向**：α=2.0 时 F out 0.02129 vs C 0.03466（−38.6%），与决
  定性腿同型；coarse 端点 F out 0.02938 vs C 0.05575（−47.3%）、距
  coarse oracle 0.03456 为 14.9%（贴边但同向）。

⇒ **BRANCH F1：带外差是片位置单独造成的。** "斜片 LS 机身升力虚高"
不是误差，是**片位置模型敏感性**；B25 注 2 的 flat-vs-tilted 口头归因
被本实验证实，且夹逼关系（A 低、C 高夹 oracle）分解为：A 低 = 袋压场
压制 + 斜片位置，C 高 = 斜片位置单独。

**副判定（带内语义）**：F band −0.0004 ≈ 0 = oracle（0.0005）类；C 的
带内高载（0.0202）**同样是位置效应**（斜片骑上机身子午面 y>0 侧对带
内三角面加载），不是 fragment 拓扑的近片载荷。带内/带外两个分量在片
位置匹配后都与 oracle 一致。

## 3. 护栏（对照 C 同 α 腿，预注册 §3.2）

| 项 | 阈 | medium 3.06 | medium 2.0 | coarse 3.06 | 判 |
|---|---|---|---|---|---|
| Δcl_p vs C | ≤2% | +1.22% | +0.87% | +0.51% | ✓ |
| Δγ vs C | ≤5% | +1.21% | +0.84% | +0.68% | ✓ |
| root te_jump 失真 vs C | ≤5% | 2.10% | 1.35% | 0.39% | ✓ |
| nlim / nflr | 0 / 0 | 0/0 | 0/0 | 0/0 | ✓ |
| sliver min 二面角 | ≥ C−5° | 11.04°（=C 同一 tet） | 同左 | 20.69° | ✓ |
| strip aux jump | O(γ) | 1.18×γ | 1.18×γ | 1.12×γ | ✓ |
| n_te_nodes | 150 / 76 断言 | 150 ✓ | 150 ✓ | 76 ✓ | ✓ |
| 走廊 corrM / n_sup_corr | ≤1.3 / 0 | 0.64 / 0 | 0.63 / 0 | 0.60 / 0 | ✓ 袋保持治愈 |
| tip Mmax | P13 监视 | 3.24 | 1.04 | 2.49 | 记录（C 同态 4.77/1.08/1.77，实现散布类） |

翼面解（cl_p/γ/根区剖面）在片位置更换下不动 ⇒ 片位置只重分配**机身**
载荷，不动翼面载荷——与"机身是近尾流压力探针"的物理图像一致。

## 4. 注 1（coarse 收敛性，R1 实发与处置）

平片几何 Picard 收敛比 ~0.94/层（斜片 ~0.77；coarse 远场巨型外单元
病态放大此效应，B16/B17 的 coarse 老问题类）。处置 = 预算 80→220 +
phi_init 暖启动续算（库配方不动）：medium 两腿 60/66 层 strict 收敛
（res ≤ 9e-8），**medium 的慢速率未显影**；coarse 腿 220 层 res
6.1e-7 未过 1e-7 门（比率 0.97/层，沿程几何下降），其分解量在 res
2.7e-5 → 6.1e-7 间 5 位有效数字不动（out 0.02938 不变）——按 B25
注 1 同例记为"意图通过、字面未过"的 RECORDED 态；判定树只用严格
收敛的 medium 决定性腿，coarse 仅作同向佐证。

成本实录：medium 4435 s + 4619 s（平片 cut 集 34242 ≈ 斜片 3×，单层
更贵），coarse 184+234 s，CutElementMap 构建 ~6 min；合计 ~2.7 h，
超预注册 §3.6 预算（40 min）——超支原因 = R1 处置（收敛速率）+
cut 集放大，记录在案。

## 5. GB9.4 重设落地（预注册 §3.4-F1）

- **B9 demo LS 腿**改采 flat-fragment（`sheet_direction=(1,0,0)` +
  `inboard_clip`；legacy farfield 配方不动保持 B9 committed 对照链；
  n_outer_max 80→220 同 R1 处置；缓存改版 `ls_flat_*.npz`）。
- **GB9.4 重设**为两款：
  (i) 硬 gate（medium）：`|cl_fus_out(conf) − cl_fus_out(LS)| ≤
  0.15·|cl_fus_out(conf)|`（coarse RECORDED）；
  (ii) RECORDED：两路径 band/out/poles/总量全分解（`cross_model_m05.csv`
  增列）。
  旧 ≤5% 阈值以 erratum 退役（物理 carryover 不再是误差，B23 §(c)）。
- **demo 重跑实测**（committed `cases/demo/b9_wingbody/results/checks.csv` /
  `cross_model_m05.csv`，8/8 PASS）：
  - **GB9.4 medium PASS**：out-band conf 0.0351 / LS 0.0376，gap **7.0%**
    ≤ 15%；band 0.0005 / 0.0008；total 0.0356 / 0.0384。
  - GB9.4 coarse RECORDED：out 0.0346 / 0.0380（gap 10.0%）；
    band 0.0002 / 0.0017；total 0.0347 / 0.0397。
  - **GB9.5 未翻界、无需重锚**：medium cl_p 0.5% / cl_kj 0.3%
    （wobble 条款未触发）；cl_fus/wing：conf 0.164、LS 0.176
    （旧斜片 LS 为 0.205）。
  - LS 收敛：medium Picard 59 层（1044 s）、coarse 73 层（35 s），
    均 conv=True，远未及 220 层预算。
- 工程纪律：knob 默认 None 逐位不变（`TestSheetDirection` 4 测试 +
  b1/m2/v0 73 绿）；demo conf 腿走 committed 缓存未动；GB9.5 若因
  配置更换翻出 <1% 带按预注册 wobble 条款记 erratum 重锚。

## 6. 结论与后续（待用户裁决）

1. **GB9.4 的"机身虚假升力"标签正式退役**：机身升力 = 物理 carryover
   基线（~10% of cl_p，两路径共有、加密不敏感）+ 袋压强印记（B25 已
   愈）+ 片位置敏感性（本 phase 定量，非误差）。B26 §4 的"C 侧
   out-band ×2"监视项就此关闭——位置敏感性，不是病灶。
2. **M2 接线嫌疑排除**：F1 意味着带外差不在拓扑/根涡强度；M2 台账
   "交界最内 TE 节点 CV fan 含机身面"不再是 cl_fus 嫌疑（该项作为
   M2 自身的遗留保留在 track_m 台账）。
3. **哪个片位置更物理**（平 y=0 顺 conforming 近似 vs 斜随来流）——
   B25 注 2 的老问题，本 phase 给出新证据：平片 = 跨模型一致配置。
   是否把 flat-fragment 推为翼身 LS 的**生产配置**（含 B18 跨声速
   demo / B26 天花板链路的 GB18.5/GB27.4 重锚）留用户裁决；本 phase
   只动 B9 demo。
4. 若采纳平片为生产配置：B18 demo 的 LS 腿、B26 的 (b) 类天花板
   归因（翼尖 P13 + 高 M Newton）不受影响（片位置不碰翼面解），但
   GB18.5 的 cl_fus 记录值会改。
