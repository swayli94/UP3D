# VERDICT — B31 翼尖片终止重定规格（C 类根治）+ LS step 语义配套评估

> **日期**：2026-07-22
> **分支**：`kimi/b30-transonic-ceiling-attribution`（用户裁决：B31 在
> B30 分支内开展）
> **预注册**：`PRE_REGISTRATION.md`（跑前写就，含 GB31.1–31.5 判定树；
> 本文件按其执行）
> **判定**：**CONF 侧 C 类根治成立**——GB31.2a 因子实验判定 taper 是
> 因（同种子同配方 ±taper：0.83 死→治），GB31.2b 生产估计器
> （pressure）+ taper 移植在 medium 门禁**治愈濒死级 0.83 并在健康
> 种子下 strict 收敛 0.84**（全收敛腿 0 钳制）；0.84 链式种子失败经
> F2 探针确诊为**种子相关的 weld 符号冻结隐患**（可修复，非估计器
> 结构缺陷）。**LS 侧 C 类关闭（阴性定论）**——C1 尖外渐隐局部治愈
> 但 inboard 回灌 −19.5%/−33.2%（F5 式全局 re-level，护栏破 20–30×），
> C3 片过尖延伸 coarse 探针直接发散（mmax 7.61，片伸到 q≈14 触及远
> 场）；预注册候选耗尽，余路 = B10 roll-up rescope。**GB31.4 关闭**
> （step 接受非两翼约束，复检钩未触发）。**采用裁决全部交用户**
> （§8）；生产配方与 demo 锚不变。
> **证据**：`results/g1_tip_atlas.csv`+`.png`、`results/g2_conf_taper.csv`
> +`.png`、`results/g2b_coarse_pressure_taper.csv`、
> `results/g2c_medium_pressure_taper.csv`+`.png`、
> `results/g3_probe_coarse.csv`、`results/g4_step_semantics.md`（由
> `run_g1/g2/g2b/g2c/g3.py` 重跑自缓存可再生成）；求解缓存
> `results/*.npz`（本分支新算，按 b23/b25/b28/b30 分析缓存惯例一并
> 提交；demo 侧种子缓存仍按 2026-07-19 用户指令不提交）。

---

## 1. GB31.1 — 钳制归属图谱：**PASS（机制定位）**

两翼濒死级全部 **8 个钳制单元 100% = cap_wall**（`g1_tip_atlas.csv`）：

| 腿 | 钳制 | cls | z | 距片尖边 | 峰 M |
|---|---|---|---|---|---|
| LS 0.7875 | 5 lim + 1 flr | cap_wall ×6 | 1.197–1.199 | 0.003–0.012 | 3.98 @ z=1.198 |
| CONF 0.83 | 2 flr | cap_wall ×2 | 1.197–1.198 | 0.003–0.010 | 1.11 |

- 钳制贴在**圆帽壁的片尖边**上（距尖边 ≤0.012），不是自由片侧边——
  B8/B30 的"翼尖"标签精确化为**片终止缘**。
- LS 终止环 δ(q)（`ring_profile`）：末环 mean 0.0451 / 各环 max
  0.0774；m2_main 从内侧 0.64 单调升到末环 **3.54**——尖侧马赫
  过冲长在终止环的 **main 侧**。
- straddler 单元 side/main M² 比：1.54 与 **7.68**（beyond_tip 者
  m2_main 4.22 vs m2_side 32.41——跨终止单元的片侧缘携带极端速度）。
- CONF Γ 尾站（station 148, z=1.1927）：Γ_last=0.0078，现状 taper
  F=0.0104 → Γ_eff=8.1e-5（**末站 Γ 已基本被 taper 卸载**，但生产
  pressure 估计器不经此路径——这正是 2b 要移植的接口）。

## 2. GB31.2a — CONF taper 因果探针：**✓（taper 是因）**

因子实验（`g2_conf_taper.csv`，同种子同配方 ±taper，probe 估计器）：

| 腿 | m | taper | conv | res | 步 | 钳制 | wall | cl_p | mmax | corrM |
|---|---|---|---|---|---|---|---|---|---|---|
| A 对照 | 0.82 | 无 | ✗ | 4.2e-6 | 80 | 0+2 | 2025s | 0.2550 | 2.67 | 1.195 |
| B 对照 | 0.83 | 无 | ✗ | 3.7e-6 | 80 | 0+3 | 6025s | 0.2592 | 2.85 | 1.207 |
| C 治疗 | 0.83 | ✓ | **✓ tol** | 2.3e-14 | 10 | **0+0** | 60s | 0.2554 | 1.82 | 1.207 |
| D 治疗 | 0.84 | ✓ | **✓ tol** | 1.9e-14 | 9 | **0+0** | 57s | 0.2603 | 1.88 | 1.218 |

- 0.83 配对（B vs C）同 M 同估计器同种子，唯差 taper：**死→治**，
  治疗腿零钳制、mmax 2.85→1.82、速度快 100× → **taper 是因**，
  天花板 0.82→≥0.84（probe 下），触发 2b。
- **U1 勘误落盘**：对照腿显示 probe 估计器本身**弱于**生产
  pressure（0.82 即死，而 B30 pressure 0.82 收敛）——估计器切换
  效应已按预注册 §4 U1 记录，不混入归因；这正是 2b 必须在
  pressure 下重验的原因。

## 3. GB31.2b — pressure+taper 生产移植：**✓（采用候选，交用户）**

**库移植**（`pyfp3d/solve/newton.py`，default-off 逐位一致）：
`F_j = t_j·σ_j·F_raw_j + (1−t_j)·s_j·Γ_j`，t_j=taper 因子；t=1 处与
生产 pressure 行逐位一致。weld 符号 s_j = 首次残差评估时
sign(diag D0)（记录为 `kutta_weld_sign`，翻号计数
`kutta_sigma_sign_flips`）。测试 54+15+35 全绿
（`tests/test_b31_pressure_taper.py` 新增）。

**coarse 验证**（`g2b_coarse_pressure_taper.csv`）：cl 代价 −1.05%
（G13.2 带内），尖 Γ 卸载至 2%，尖峰 M 0.994→0.803，0 翻号。

**medium 门禁**（`g2c_medium_pressure_taper.csv`，strict 生产配方
pressure + vanish_smooth r_c=0.05·b_semi）：

| 腿 | m | 种子 | conv | res | 步 | 钳制 | σ翻号 | cl_p | mmax | corrM |
|---|---|---|---|---|---|---|---|---|---|---|
| G 代价对照 | 0.82 | b30 0.81 | ✓ tol | 6.3e-11 | 8 | 0+0 | 0 | 0.2640 | 2.02 | 1.134 |
| E 治愈级 | 0.83 | b30 0.82 | **✓ tol** | 2.4e-14 | 11 | **0+0** | 0 | 0.2684 | 2.10 | 1.157 |
| F 爬升级 | 0.84 | 链 E | ✗ | 4.9e-5 | 80 | 97+63 | **3** | (0.2750) | 9.95 | — |
| F2 诊断 | 0.84 | 健康 2a-D | **✓ tol** | 2.1e-14 | 13 | **0+0** | 1 | 0.2761 | 2.82 | 1.180 |

- **E（B30 濒死级）被生产配方治愈**：strict、0 钳制、corrM 护栏
  ≤1.3 ✓ → 判定树 **GB31.2b ✓**，采用候选交用户。
- **代价（honesty 标注）**：G 腿 cl_p 0.2640 vs B30 pressure 0.82 的
  0.2722 = **−3.00%**，超出 F3 单翼带（−1.1..−1.6%）——翼身构型
  的 taper 升力代价高于单翼，采用即接受此代价（或先做代价分解）。
- **0.84 诊断（F 死→F2 生）**：F 从链式种子（E 的 0.83 态）出发时
  首次残差冻结到 **3 个翻号站**的 weld 符号，尖 Γ 被反向 pin 到
  −1.5e-4，求解落入 97 lim + 63 flr 极限环；F2 换健康种子（2a-D
  的 0.84 收敛态，σ翻号=1）后**同配方同库代码 strict 收敛、
  0 钳制**。→ 失败是**种子相关的符号冻结隐患**（移植设计时已
  登记的 hazard），**不是估计器结构缺陷**；修复提名见 §8-②。

## 4. GB31.3 — LS 片终止函数空间 re-spec：**C1 ✗ / C3 ✗ → LS 侧 C 类关闭（阴性）**

coarse 探针（`g3_probe_coarse.csv`；baseline 5 步 0 钳制收敛，
mmax 1.298；库改动 default-off 逐位一致，19/19 新单测 + 129 回归绿）：

| 候选 | conv | 步 | 钳制 | 尖峰 mmax | inboard 漂移 | 判定 |
|---|---|---|---|---|---|---|
| C1 fringe_fade（尖外渐隐） | ✓ | 5 | 0+0 | **0.87**（治愈） | cl **−19.51%** / Γ **−33.24%** | **✗** 护栏破 20–30× |
| C1 消融 fringe_only（仅移位） | ✓ | 6 | 0+0 | 2.18（新边更差） | cl +1.96% | （消融对照） |
| C3 片过尖延伸（no-clip） | **✗** | 80 | 21+25 | **7.61**（尖部） | cl +23.2%（非解态） | **✗** 直接发散 |

- **C1 ✗（F5 式全局 re-level 重演）**：尖峰确实治愈（mmax→0.87，
  δ 0.037→0.0013 光滑衰减），但 inboard q<0.95 span 的 cl/Γ 被
  回灌 −19.5%/−33.2%，三条通道已定位：q=span 共享 aux DOF 排流、
  椭圆主场耦合、pin_gamma 标量反馈。渐隐虽全部在尖**外**，"尖外
  不回灌"的机制论证**不成立**（预注册 U2 风险实发）。
- **C3 ✗（发散 + 远场侵入）**：no-clip 延伸使 beyond-span 切割单元
  44→**220**、触及远场 6→**26**、片伸到 q≈14.2；pin-vs-local 语义
  核查一致（mismatch ≤2e-5，排除实现 bug）——是真实的动力学
  失败：探针 80 步 21 lim + 25 flr、尖部 mmax 7.61。
- **medium 门禁按 triage 规则豁免**：候选须在 coarse（baseline
  5 步收敛处）至少收敛才值得花 medium 解；两候选均未过。用户可
  推翻此豁免要求补跑 LS 0.7875 medium 门禁。
- 预注册候选阶梯（C1→C3）耗尽且无效应/更差 → **LS 侧 C 类关闭
  （阴性定论）**，余路 = B10 roll-up rescope（§8-③，交用户）。

## 5. GB31.4 — LS step 语义：**关闭（证据性；复检钩未触发）**

`results/g4_step_semantics.md`：LS 已有 best-of-tried 回退线搜索
（newton_ls.py:899-922），"移植"框架正式退役；两翼濒死级死法均非
step 接受病理（LS = active set 卡死 + freeze 窗语义；CONF = 光滑
Newton 停滞）。**复检钩评估（GB31.3 之后）**：C1 收敛步态正常
（5–6 步 0 钳制）；C3 死于钳制计数增长 + 残差平台（与 baseline
死级同类），无 trial 全非有限 / λ 打满地板 / merit-inf 两表脱同步
证据 → **维持关闭**。

## 6. 成本与护栏

- 求解成本（缓存前）：GB31.2a 4 解（对照 2025+6025s，治疗 60+57s）；
  2b coarse + medium 4 解（108+90+475+86s）；GB31.3 coarse 探针
  4 变体（≤1025s/个）；GB31.1/4 ≈ 0（缓存/文献）。均在预注册 §6
  预算内。
- 护栏：走廊 corrM ≤ 1.3 ✓（全部收敛腿 1.13–1.22）；honesty ✓
  （F 腿连同 97+63 钳制与 σ翻号=3 全报；F2 标注为换种子诊断探针，
  不计为爬升链收敛；2a 对照腿记录 probe 估计器弱点；未收敛态不
  计天花板）。
- 生产配方与 b18 demo 锚**不变**（采用裁决待用户，锚刷新条款在
  用户采纳 ✓ 后才执行）。
- 库改动全部 default-off 逐位一致：`pyfp3d/solve/newton.py`
  （pressure+taper blend）、`pyfp3d/wake/cut_elements.py` +
  `pyfp3d/wake/multivalued.py`（outboard_fringe 哨兵）；新增测试
  `tests/test_b31_pressure_taper.py`、`tests/test_b31_tip_fringe.py`；
  `tests/test_p14_te_pressure.py` 一处配套更新。

## 7. 预注册判定树终态

- GB31.1 **PASS**（钳制 = 片终止缘 cap_wall，机制定位）→
  GB31.2a **✓**（taper 是因，0.83 死→治配对干净）→
  GB31.2b **✓**（生产估计器 + taper：0.83 治愈、0.84 健康种子
  收敛；代价 −3.00% 超 F3 带如实标注；0.84 链式种子失败确诊
  为可修复的符号冻结隐患）→
  GB31.3 **C1 ✗ / C3 ✗**（LS 侧 C 类关闭，阴性定论）→
  GB31.4 **关闭**（复检钩未触发）→ GB31.5 = 本文件。

## 8. 出口路由（全部交用户裁决）

1. **CONF taper 采用**（本 phase 提名的主出口）：生产配方
   `CONF_RAMP_NK` 引入 tip_taper（vanish_smooth，r_c=0.05·b_semi），
   配套 demo 锚刷新。收益 = 濒死级 0.83 治愈、天花板 0.82→≥0.84；
   代价 = cl_p −3.00%（翼身，超 F3 单翼带）——采纳即接受代价，
   或先立项做代价分解（taper 半径/形式扫描）。
2. **weld 符号冻结修复**（F2 确诊的小改，建议作为 ① 的配套前提）：
   `kutta_weld_sign` 由首次残差冻结改为逐步刷新（或翻号检测时
   重冻结），消除链式种子在 0.84 触发的极限环。库改动、自带
   测试，工作量小。
3. **LS 侧余路 = B10 roll-up / 显式尖涡 rescope**（预注册 §5 路由）：
   LS 片终止函数空间 re-spec 候选已耗尽（C1 回灌、C3 发散），
   若仍要抬 LS 天花板需模型级根治，另立 phase。
4. **`freeze_max_clamped` 重定规格**：仍为用户保留项，本 phase
   未动（LS 0.80@c=2.0 一单元之差的窗语义问题原样保留）。
5. **merit 语义升级**：GB31.4 证据性关闭，**不提名**。

**demo 锚刷新：不执行**（等待用户对 ① 的采纳决定）。
