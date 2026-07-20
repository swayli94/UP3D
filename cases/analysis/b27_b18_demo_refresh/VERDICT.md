# VERDICT — B27 B18 demo 刷新（LS 腿复活）+ Track B 文档收尾

> **日期**：2026-07-20
> **分支**：`kimi/b25-inboard-fragment-clip`（叠在 B25/B26 上，未 push 链）
> **判定**：**B27 落地** —— GB27.1/GB27.2/GB27.3(0.65) **PASS**，
> GB27.3(0.75)/GB27.4/GB27.5 **RECORDED**（均按预注册 §2.4 预期）。
> B18 demo 门面从"LS junction-limited（closed-negative）"重述为
> "袋愈后 LS 天花板与 conforming 同址"，B26-A 由此落账。
> **证据**：`results/g27_consistency.csv`（336/336 bit-identical）+
> 刷新的 `cases/demo/b18_wingbody_transonic/results/`（checks.csv 8/8
> PASS、cl_vs_mach.csv、cross_model.csv、3 PNG）。

---

## 1. 六腿结果（全量重解，无缓存起步，总求解 ~1 h 39 min，T2 预算内）

| # | 腿 | 结果 | 对锚 |
|---|---|---|---|
| L1 | conf coarse 0.84（0.70 起）+ coarse 0.60 cross | reached 0.84，cl_p **0.2617**，Mmax 2.15；0.60 cross conf **0.2178** | bit = B18 committed ✓ |
| L2 | conf medium 0.65（0.60）+ 0.79（0.70） | **0.2321** / **0.2579** | bit = B18 committed ✓ |
| L3 | conf medium 0.75（0.70）**新点** | reached strict，cl_p **0.2483**，Mmax 2.27，res 8.3e-11 | cl(M) 单调 0.2173/0.2321/0.2483/0.2579 ✓ |
| L4 | LS A/C ceiling probe coarse 0.50→0.84 | A：m_last **0.82** 死 0.84（b+dm，条带袋 M3.28 q≈0 卡 Newton）；C：**reached 0.84**，cl_p 0.2542 | bit = B26 committed ✓ |
| L5 | LS A/C ceiling probe medium 0.50→0.84 | A：m_last **0.50** 死 0.5125（a+dm，Mmax 6.17）；C：m_last **0.7625** 死 0.775（b+dm，翼尖 M4.18），濒死态 cl_p 0.2475 | bit = B26 committed ✓ |
| L6 | LS+clip medium 定向 0.50→0.65、0.50→0.75 | 两腿均 reached（strict 末级，0/0 与 0/1 clamp）：cl_p **0.2266** / **0.2421** | 新测量（跨模型第二/三横坐标） |

各腿耗时：conf coarse ~17 s/级量级；LS A/C coarse 2×~2 min；conf medium
0.65/0.75/0.79 ≈ 6+6+6 min；LS A medium 26 min、C medium 39 min
（B26 实测 1560 s/2333 s 同量级）；定向腿 ~5+6 min。

## 2. 预注册 gates 逐项核对

- **GB27.1（PASS）**：conforming 腿复现 committed 锚，`cl_vs_mach` 4 位
  小数全一致（0.2173/0.2321/0.2579/0.2617 + cross 0.2178 + M_max 2.15 +
  reached 标志）。B21/B22 对 conforming 路径惰性成立——`g27_consistency.csv`
  conforming 8/8 bit，无 T1 漂移。
- **GB27.2（PASS）**：LS ceiling A/C 复现 B26 committed——C coarse
  reached 0.84；C medium m_last=0.7625、死 0.775 (b) 类；A medium
  m_last=0.50、死 0.5125 (a) 类；A coarse m_last=0.82。summary 36 项 +
  levels 292 项全 bit-identical（同码同 16 线程，零 T5 漂移）。
- **GB27.3（0.65 PASS / 0.75 RECORDED）**：新跨模型——0.65：
  conf 0.2321 vs LS+clip 0.2266，**|gap| 2.4% ≤ 5%** PASS；0.75：
  conf 0.2483 vs LS+clip 0.2421，|gap| **2.5%**（RECORDED，不设阈值）。
  ★ 跨模型 gap 现在全 Mach 平整：M0.5 2.6% / M0.65 2.4% / M0.75 2.5%
  ——同一条 ~2.5% 带 = B17 已知的 LS cl_p↔cl_kj 口径差，非物理分歧。
- **GB27.4（RECORDED）**：GB18.4 重答入 demo（袋 = B23 inboard 自由边
  奇点，C 侧已愈；残余限制器 = 翼尖 P13 类 + 高 M Newton，峰位引 B26
  committed g1_peaks：A 濒死峰 M6.17 @ 交界条带 x=2.12 + 贴面 M3.53
  dist_fus=0.005，C 濒死峰 M4.18 @ 翼尖 z=1.197，走廊 corrM 1.07
  干净）；GB18.5 刷新（conf cl_fus 0.0423 = 16% @0.79 live；C 侧新天花板
  cl_fus 0.0781 / band 0.0216 / out-band 0.0565，引 committed
  g1_summary——out-band ×2 为 P11/曲面壁元输入）。
- **GB27.5（RECORDED）**：T1 勘误已入 demo docstring——A 侧复测爬过
  B18 committed 锚（死 0.50/0.55）是 **B21/B22 freeze-capture 修复效应**
  （B26 §3 独立发现），非物理漂移；袋真实杀伤线 A medium = 0.55
  （Mmax 13.1 > freeze_max_clamped=8）。

## 3. 独立观察（不阻塞判定）

1. **coarse 0.60 cross 的 LS 侧换了口径**：旧 committed 行是 A 侧
   0.2174（gap 0.2%）；本次 C 侧 0.2133（gap **2.1%**，increment 49%，
   under-resolved 注记保留）。读法：A 侧粗网格上的"近 0 差"是袋污染态
   恰好落在 conforming 附近的巧合；愈后 C 侧落回与全表一致的 ~2–2.6%
   口径差带（§2 GB27.3 ★ 行）。旧 0.2% 行就此退役。
2. **`b18_sections_conf_medium.png` 在旧 committed demo 里就是空图**：
   `section_cp_curve` 早已改返回 dict（x_upper/cp_upper/…），旧代码按
   tuple 解包被 `except: pass` 静默吞掉。B27 顺手修复（上下翼面两条
   曲线），不属于 B21/B22 回归——是更早就存在的 API 漂移。
3. **conf medium 0.75 新点** Mmax 2.27、res 8.3e-11 strict 收敛，cl(M)
   单调性成立，0.80+ 停滞叙事不变。

## 4. 工程纪律核对

- 库代码零改动；A 侧 = 默认 clip 同码复测；C 侧 = B25 已入库
  `inboard_clip`；LS_RAMP_KW / CONF 配方逐字冻结。
- 前置门 b1/b2/m2/v0 全绿（86/86，24 s）；线程上限 16；
  `PYFP3D_TRANSONIC_GATES=1` 门控惯例。
- 证据 = committed CSV/PNG；npz/msh gitignored。
- 求解 ~1 h 39 min（预注册 T2 预算 2–2.5 h 内）。
- 未 push；push/PR 待用户确认。

## 5. 文档落点（本 phase 收尾项）

- `cases/demo/b18_wingbody_transonic/run_demo.py` docstring + checks 注
  = 预注册 §2.4 叙事逐字落位。
- `docs/demo_report.md` B18 行刷新；`docs/design_track_b.md` 新 §22
  （B26/B27）+ §18 erratum 指针；`docs/agent-rules.md` B26/B27 条目 +
  next-phase 优先级；`docs/overview.md` B-track 表 GB20.5 故事勘误
  （袋归因：刻面 G1.6 → 自由边奇点已愈；G1.6 退居 cl_fus 嫌疑）；
  `docs/roadmap.md` B-track 行 B18 描述刷新。
- 遗留 inspection 债（预注册 T4 条款，本 phase 不扩 scope）：
  `docs/demo_report/track_b.md` 的 B18 详章（~L1756 起）与 L8 头部仍讲
  旧"交界受限"故事；`docs/agent-rules.md` / `docs/roadmap/track_b.md`
  缺 B23/B24/B25 条目（本分支未 push 链的既存债）。均留待合入时统一收尾。

## 6. 后续（预注册 §6 出口）

1. **LS 翼身全包线评估开题**（B25 §7 落地）：LS+clip 与 conforming
   天花板同址（coarse 0.84 = 0.84；medium 0.7625 ≈ 0.79），conforming
   退居交叉验证，Track V 片拓扑前置全部就位。
2. **(b) 类天花板归因**（下一 phase 候选）：LS+clip medium 0.775 与
   conforming medium 0.80+ 停滞若同机制（翼尖 P13 奇点 + 激波系
   lim/flr 振荡），后续投入转向两路径共享的 Newton 鲁棒性。
3. **cl_fus out-band ×2**（C 侧 0.057 vs A 侧 0.022）→ P11/曲面壁元
   路线输入监视项。
