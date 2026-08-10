# 翼身交界判别实验 — 预注册（PRE-REGISTRATION）

> **日期**：2026-07-19（在任何测量运行之前写就）
> **分支**：`kimi/wingbody-junction-discriminator`（Kimi 工作副本；不动 main）
> **执行令**：P11 close-out（commit 3ab214f，roadmap track_p.md §P11）——
> "GB9.4/GB20.5 need their own discriminator (junction = crease geometry ≠
> smooth-wall faceting) before any cross-attribution is quoted again."
> **用户目标**：交界/机身表面数据可信（本判别 → 证据裁决），随后实现 LS 工作流。

## 判别对象（两个伤口，基线锚 = committed 证据，不重算）

| 伤口 | 基线（B20 re-baseline 2026-07-19） | 随加密行为 |
|---|---|---|
| **W1** LS 翼身交界虚假超声速袋 | medium M0.5: Mmax **5.22** genuine（nlim 3/nflr 3，res 1.1e-13 收敛态）；coarse Mmax 1.31 @ ~M0.55 | **恶化**（1.31 → 5.22） |
| **W2** 机身虚假升力 | GB9.4：\|cl_fus\|/cl_p_wing = conf **0.16** / LS **0.20** @ medium | LS **增长**（0.164 → 0.205），conf 持平 |

辅助锚：GB9.6 孤立机身护栏（`cases/analysis/b9_fuselage_guardrail/`）——
方位角 Cp 散布 median 随加密衰减 0.0036/0.0022/0.0010，但 max 增长
0.042/0.096/0.117（鼻/尾极点）。光滑机身蒙皮中位数正常 ⇒ 全局场阶不缺；
最大值长在几何退化点 ⇒ 局部奇点嫌疑。

P11 已测死、本 campaign **不复查**的路线：mapped-P1 + 曲面几何（GLM 路线 B；
medium Cp 增益 0.23pp = oracle ceiling）。"G1.6 类"标签对翼身**禁用**（球锚已溶解）。

## 四个实验（D1–D4）

### D1 — α 扫描：袋是否需要升力/尾流？（现有网格，只算）

- 设置：LS Picard（B9 实测翼身配方 `farfield="freestream"`，B17 仲裁
  `farfield_aux="pin_gamma"`），M0.5；medium α ∈ {0, 1, 2, 3.06}，coarse α ∈ {0, 3.06}。
- 度量（每个状态）：Mmax = sqrt(max element_mach2)；袋位置（argmax 元素质心的
  x/y/z、到机身面距离 |√(y²+z²)−R(x)|、到交界 TE 点距离、所在 z 相对 z_junc）；
  |cl_fus|/cl_p_wing；收敛性（res、n_outer、nlim/nflr）。
- **判定（预注册）**：
  - α = 0 已有 Mmax ≫ 1 的袋 → **几何驱动**（升力/尾流/Kutta 均不在场）；
  - α = 0 无袋（Mmax = O(1) 合理值），袋随 α 出现并增长 → **升力/尾流耦合驱动**
    （交界 TE 内端 = level-set polyline 内端 + 隐式 Kutta 最内侧站位）；
  - 自检：α = 0 时 cl_fus ≈ 0（对称性），否则测量通道本身可疑。

### D2 — 交界局部加密：袋是不是交界区网格的函数？（固定 h_wall，只变 h_junction）

- 网格：medium 族固定 h_wall = 0.015，h_junction ∈ {0.015 (1.0 h_wall)、
  0.0075 (0.5，= 基线 medium)、0.00375 (0.25)、0.001875 (0.125)}；其余尺寸律不变
  （自相似族策略）。基线 medium 网格已有，新生成 3 个变体。
- 求解：M0.5，α = 3.06，同 D1 配方。
- **判定（预注册）**：
  - Mmax 随 h_junction↓ **增长** → 角点奇点"越分辨越尖"签名（加密捕捉到更多
    奇异峰）→ **折痕奇点驱动**；
  - Mmax **下降** → 欠分辨驱动（须与 coarse→medium 全场恶化 1.31→5.22 对照解读：
    全场加密变了一切，本实验只变交界区）；
  - Mmax **不动** → 非交界区网格因素（远场/全场/求解器层）。

### D3 — 整流片（fillet）A/B：杀死折痕能否杀死袋？

- D3a（规范案例腿）：在 D4 的规范几何上做折痕倒圆 A/B（几何简单，OCC fillet 可行）。
- D3b（真实翼身腿）：对真实翼身流体域的交界棱边试 `occ.fillet`
  （半径 ~0.03–0.05 = 0.2–0.33 r_f）。OCC 倒圆若在合理尝试内拒绝，
  **记录为"低成本不可行"并退回 D3a + D2 证据**，不陷入几何攻坚战。
- 求解：M0.5，α = 3.06（D3a 加 α = 0），同 D1 配方；网格策略与未倒圆变体一致。
- **判定（预注册）**：固定网格策略下袋**消失/显著减弱**（Mmax 降到 O(1) 量级）
  → 锐折痕确认 + 缓解手段在手（真实飞机本有整流罩，物理合法）；
  袋**依旧** → 袋不是锐折痕本身产生的。

### D4 — 规范折痕案例 h 阶梯：折痕类能否在干净几何上复现？

- 几何：圆柱（r = 0.15，沿 x）+ 对称翼型平板（M6 剖面，无弯度）从圆柱面伸出；
  **无尾流、非升力**（α = 0；可加压差腿 α > 0 作对照）。小远场域控成本。
- 求解：`solve_subsonic`（M0.5，可压）+ Laplace（不可压）两条腿；h 阶梯 3 级。
- 度量：Mmax(h)、|q|max/|q∞|(h)（Laplace）、袋位置。
- **判定（预注册）**：规范几何复现"随加密恶化"→ **折痕类确认**
  （与 M6 翼身/尾流/机身长度无关的几何类）；不复现 → 真实翼身的袋需要
  翼身特有因素（尾流内端/近场/升力），矛头转回 D1 的判定。

## 裁决树（第三步 = 拿证据后执行，预注册）

1. **D1 α=0 有袋 ∧ D4 复现 ∧ D3 杀死** → **折痕几何驱动**。
   路线 = 几何/网格层处理（fillet 或交界分级尺寸策略），天–周级；P2 无指征。
2. **D1 α=0 无袋 ∧ 袋随 α 增长** → **升力/尾流耦合驱动**。
   路线 = 修交界 TE 内端处理（Track B wake 层：level-set polyline 内端 +
   最内 Kutta 站位在折痕角点处的行为），非几何相位。
3. **D2 袋随加密增 ∧ D3 杀不死** → **角点场阶不足**。
   P2（P11 close-out 三岔口选项 b）重回桌面，预期管理 = "expectation, not promise"。
4. 组合/矛盾/均不显著 → 记录"不可判别"，列下一轮判别实验，不强行归因。

W2（机身升力）的读取：同一批测量的 |cl_fus|/cl_p 随 α / h_junction / fillet
的行为；结合 GB9.6 极点证据区分"交界贡献"与"极点贡献"。若 W2 与 W1 的归因
不同步（例如袋是折痕的、升力是极点的），分别记录，不捆绑。

## 工程纪律

- 线程：NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16。
- 求解缓存：本目录 `results/*.npz`（gitignored）；committed 证据 = CSV/PNG。
- B20 demo-cache 教训：首跑前确认无陈旧 npz；日志中不得出现 silent reuse。
- 不改任何 kernel/assembly；若 meshgen 加 fillet 参数，默认行为 bit-identical。
- 所有数字以 committed CSV 为准；图 = 各 leg PNG。
