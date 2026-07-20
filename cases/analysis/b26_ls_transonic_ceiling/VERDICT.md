# VERDICT — B26 袋愈后 LS 跨声速天花板重测

> **日期**：2026-07-20
> **分支**：`kimi/b25-inboard-fragment-clip`（叠在 B25 上，未 push）
> **判定**：**B26-A（天花板抬升）成立** —— 袋是 LS 翼身跨声速天花板的
> 限制器；袋愈后 medium 天花板 0.50 → **0.7625**，coarse 0.82 → **0.84
> （reached）**，且两档 C 侧死因均从 (a) 类袋拒级变为 (b) 类高 M Newton
> 停滞（峰在翼尖，P13 类，非袋）。
> **证据**：`results/g1_summary.csv`、`results/g1_levels.csv`、
> `results/g1_peaks.csv`、`results/g1_ceiling.png`（均由
> `run_g1.py` 无参全量重跑自缓存重新生成）。

---

## 1. 四腿结果（同码 A/C 对照，B18 配方逐字冻结）

| 腿 | m_last | 死于 | 失效分类 | 死因细节 | 峰位 | wall |
|---|---|---|---|---|---|---|
| A coarse | 0.82 | 0.84 | (b)+dm | Newton 停滞 res 6.3e-4、lim/flr 4/9 振荡、γ 冻结 | 机身旁条带 M3.28（q≈0） | 136 s |
| C coarse | **0.84 reached** | — | — | 0.84 strict 收敛 res 6.9e-11（freeze 接住 6/1） | 翼尖 M14.3（P13 类） | 85 s |
| A medium | 0.50 | 0.5125 | **(a)+dm** | 0.55 loose res 1.1e-13 但 8/3（>freeze_max_clamped=8）拒级；0.525/0.5125 strict 重试同拒 | 机身旁条带 M6.17 + 交界贴面 M3.53（dist_fus=0.005） | 1560 s |
| C medium | **0.7625** | 0.775 | (b)+dm | 0.80 loose res 8.2e-6（22/9）→ 0.775 strict 两次 res ~2e-6（8/7）停滞；0.7625 strict 收敛 res 2.6e-11 | 翼尖 M4.18（z=1.20，P13 类）；走廊 corrM=1.07 干净 | 2333 s |

名义阶梯 0.50→0.84 dm=0.05（8 rung）；A/C 共用同一阶梯；ramp 诚实
停止。合计求解 ~69 min（预注册 T2 预算内）。

## 2. 预注册判定逐项核对（§3 B26-A）

- **C medium m_last_converged ≥ 0.60**：0.7625 ✓ —— 爬过 B18 死亡点
  0.50 共 5 个 loose rung（0.55/0.60/0.65/0.70/0.75，全部 0–1 clamp）
  + 0.7625 strict（res 2.6e-11，1/2，freeze 接受）；warm-start 链完整。
- **收敛级无 (a) 类失效** ✓ —— C medium 全部收敛级 clamp ≤ 1/2 且经
  freeze 接受；濒死级是 (b) 类 Newton 停滞（res ~2e-6 不紧致），不是
  strict 门拒级。
- **cl_p 与 conforming 锚同趋势** ✓ —— 收敛级 cl_kj 单调
  0.1289@0.50 → 0.1491@0.7625（+16% 跨声速抬升）；cl_p 0.2542@0.84
  （C coarse 收敛态）与 0.2475@0.775（C medium 濒死态，res 2e-6 近收
  敛）对照 conforming 锚 0.2617@0.84 / 0.2579@0.79，同趋势、低
  2–4%（LS cl_p 对 cl_kj 的既有口径差，B17 同类）。
- coarse 侧同向佐证：A 死于 0.84（条带袋 M3.28 在 q≈0 卡住 Newton），
  C reached 0.84（走廊 corrM 3.28→1.10，袋消失）。

## 3. 独立发现（预注册 T1 条款）：A 侧与 committed 锚的背离 = B21/B22 效应

A 侧重测（主对照，P14 同码纪律）与 B18/GB20.5 committed 锚显著背离，
按预注册 2.2 条款记为 **B21/B22 freeze-capture 修复的独立发现**，不影
响 B26 判定：

- **A coarse**：锚 = 死于 0.55 第一级（Mmax 1.31）；现码 = 爬到 0.82。
  coarse 天花板主要是 freeze-capture bug，袋在 coarse 只在 0.84 才显影。
- **A medium**：锚 = 死于 0.50（GB20.5 main：res 1.1e-13、3/3、
  Mmax 5.22、807 s）；现码 = **0.50 收敛**（res 9.5e-11、同为 3/3、
  **Mmax 同为 5.22**）——袋在 0.50 依旧存在，但修复后的 freeze 现在能
  捕获这 3/3 个单元让 strict 门通过；袋在 0.55 爆发（Mmax 13.1，
  8/3 > freeze_max_clamped=8）重新封死 ramp，死因仍是 (a) 类。
- 即：B18 的"0.50 死亡"是 capture bug 与袋的叠加；B21/B22 剥掉了前一
  层，袋的真实杀伤线在 A medium 上是 0.55（0.525/0.5125 重试同拒）。

## 4. 物理结论

1. **袋是天花板限制器（medium 决定性）**：同一配方、同一网格、同一
   代码，唯一变量 inboard_clip → medium 天花板 0.50 → 0.7625
   （+5 rung），coarse 0.82 → 0.84（A 死在袋位点，C 袋消并收敛）。
2. **新天花板是物理/Newton 类而非离散误差**：C 两侧濒死级峰都在
   **翼尖**（z≈1.20，P13 已知独立奇点），走廊全程干净
   （corrM ≤ 1.10，峰在尾端条带 x≈2.40，q≈−0.13）；死因 (b) 类与
   conforming medium 0.80+ 停滞同类。LS 天花板现在与 conforming 天花
   板同址：coarse 0.84 = conforming coarse；medium 0.7625 ≈
   conforming medium 0.79。
3. **guardrail 全绿**：条带 aux |jump|/γ = 1.17–1.18（B25 的 1.16 一
   致，高 M 下无锚定病）；sliver 最小二面角 11.0°（C medium，健康）；
   n_te 76/150 不变（M2 锁断言通过）；n_aux_symmetry 68/150（coarse/
   medium，预注册度量 7 记录）；cl_p 收敛态对 conforming 锚 −2.9%。
4. **监视项（不阻塞判定）**：C 侧 cl_fus ≈ 0.076–0.078，约为 A 侧
   （0.038–0.048）两倍，增量主要在 out-band 分量（0.057–0.068 vs
   0.022–0.033）——条带贴体区的机身带外升力，GB9.4 类，留给 P11/曲
   面壁元路线。

## 5. 工程纪律核对

- 库代码零改动（A = 默认 clip 同码复测；C = B25 已入库
  `inboard_clip`）；LS_RAMP_KW 逐字冻结（`run_demo.py:90-92` 原
  样）；α=3.06、dm=0.05、dm_min=0.01、upwind 默认。
- 前置门禁 b1/b2/m2/v0 全绿（86/86，24 s）；线程上限 16。
- 证据 = committed CSV/PNG；npz 缓存本地（gitignore 政策同 B25）。
- 未 push；push/PR 待用户确认。

## 6. 后续（预注册 §3 B26-A 出口，待用户裁决）

1. **B18 demo 刷新**（LS 腿复活：coarse 0.84 reached、medium
   0.7625）+ GB9.4/GB20.5 重设 + Track B 文档收尾。
2. 新天花板归因：(b) 类高 M Newton 停滞（翼尖 P13 奇点 + 激波系
   lim/flr 振荡），与 conforming medium 0.80+ 停滞是否同一机制——
   若是，LS 与 conforming 在翼身上已同极限，跨声速 LS 包线评估开题
   （B25 §7 落地）。
3. cl_fus out-band 增量（C 侧 ×2）→ P11/曲面壁元路线输入。
4. A medium 0.55 的袋爆发（Mmax 13.1，超过 freeze_max_clamped=8 可捕
   获上限）→ 若需继续抬 A 侧，freeze_max_clamped 与袋的交互值得记录
   （C 侧已不需要）。
