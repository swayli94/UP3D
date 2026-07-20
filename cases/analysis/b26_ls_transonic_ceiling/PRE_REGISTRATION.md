# PRE-REGISTRATION — B26 袋愈后 LS 跨声速天花板重测（B18 LS Newton freeze-ramp 重跑）

> **日期**：2026-07-19
> **分支**：`kimi/b25-inboard-fragment-clip`（本地叠在 B25 上，未 push 链）
> **状态**：预注册，待用户裁决后实施
> **前置**：B18 committed 结果（LS 翼身 Newton 天花板 M0.50/0.55）、
> GB20.5（袋归因 G1.6）、B20 勘误（袋 = 真实收敛场）、B25 VERDICT
> （袋愈：corrM 14.66→0.63，Picard 腿 nlim=nflr=0）

---

## 1. 物理推导链（为什么是这一条）

**B18 的失效模式（committed `b18.../checks.csv` + GB20.5 同码重跑
`c1.../legb_b18_hypothesis.csv`）**：LS Newton freeze-ramp 在翼身上
**第一级就死**——coarse 死于 m=0.55（Mmax=1.31），medium 死于 m=0.50
（Mmax=5.22，nlim 3 / nflr 3）。失效**不是残差发散**：post-B20 下该级
Newton 残差收敛到 1.1e-13，但交界袋处的 floored/limited 单元使
strict 接受门（要求 0 limited / 0 floored，track_b.md GB16.3 注）否决
该级 → ramp 死。即：**天花板 = 袋触发限幅/截断计数，接受门拒级**。

**B25 F1（本分支已 verdict）**：袋 = 尾流内侧自由边奇点；`inboard_clip`
把片内侧边界移到机身表面/对称面（conforming fragment 拓扑）后，
medium α=3.06 Picard：corrM 14.66→0.63、走廊 n_sup 88→0、
**nlim=nflr=0**、条带跳量 1.16×γ（无锚定病）、无 sliver。

**假设**：若袋是天花板的限制器（接受门所拒的 nlim/nflr 全部来自袋），
则同一 ramp 配方 + `inboard_clip` 应爬过 B18 死亡点。**零结果**（天花
板不动）同样有价值：说明限制器住在别处（G1.6 刻面几何、wake-LS 高
M 条件数、upwind 机制），重新归因。

**与 P11 的区隔**（机制相邻的并行路线）：GB20.5 把袋归因于 G1.6 刻面
几何误差（P11 曲面壁元攻击的根因）；B25 已证明拓扑自由边是袋的驱动
奇点并治愈之。B26 测的是**操作后果**（天花板动没动），不裁决两种归
因的份额；若 B26-B（天花板不动），G1.6/P11 升为首席嫌疑。

**与 B18 docstring 的区隔**：`run_demo.py:12-19` 的数字是初版遗留，
与 committed checks.csv（post-B20 刷新）不一致——一律以 checks.csv
为准。

## 2. 实验设计

### 2.1 改动面（预注册承诺）

**库代码零改动**。`inboard_clip` 已在 B25 入库（`cut_elements.py`，
默认 None bit-identical；`fuselage.make_inboard_clip`）。

| 文件 | 改动 |
|---|---|
| `cases/analysis/b26_ls_transonic_ceiling/` | 新建：本文件、`wb26.py`（复用 `wb25.build_ls_clip` + `wb_common`）、`run_g1.py` |
| `results/g1_summary.csv`、`results/g1_ceiling.png`、`VERDICT.md` | 产出 |

B18 配方**逐字冻结**（`run_demo.py:90-92`）：

```python
LS_RAMP_KW = dict(farfield="freestream", farfield_aux="pin_gamma",
                  freeze_tol=1e-4, freeze_max_clamped=8,
                  intermediate_tol=1e-3, n_seed=30,
                  direct_refactor_every=1000, n_newton_max=80)
```

求解器 `solve_multivalued_newton_transonic`（newton_ls.py:1000），
α=3.06 固定，dm=0.05，dm_min=0.01（默认）。upwind 默认
（c=1.5, m_crit=0.95, m_cap=3.0, rho_floor=0.05）。

### 2.2 腿（A/C 同码对照，A 兼作 B18 复现门）

两侧共用同一标称阶梯 **m_start=0.50 → m_target=0.84，dm=0.05**
（0.84 = conforming coarse 天花板，物理上限参照）。A 侧
`inboard_clip=None`（bit-identical，即 B18 配方的当前码复测）；C 侧
fragment clip。ramp 诚实停止（失败级减半至 dm_min 耗尽）。

| 腿 | level | 说明 |
|---|---|---|
| A | coarse | 同码对照 + B18 复现门（对 checks.csv 行 3） |
| C | coarse | fragment clip（冒烟 + 趋势） |
| A | medium | 同码对照 + GB20.5 锚（legb_b18_hypothesis.csv） |
| C | medium | **主判定腿** |

执行顺序：A coarse → C coarse（冒烟，先验设施）→ A medium → C medium。

**A 侧复现门**：A medium 应复现 committed 锚——死于 m=0.50 级、
Mmax≈5.22、nlim/nflr≈3/3（checks.csv 行 5；GB20.5 main 腿 res
1.1e-13 同址）。注意 B21/B22 动过 Newton freeze 捕获（n1_freeze_fix
sweep），当前码与 B18 提交时点不严格同码：A 重测是**主对照**（P14
纪律：同码测量），committed CSV 是**历史锚**；若 A 重测与锚显著背
离，记为 B21/B22 效应的独立发现，不判 B26 失败。

### 2.3 度量（每级 × 每侧）

1. **主度量**：`m_last_converged` / `m_final` / `target_reached`
   （ramp 诚实字段，B18 同款天花板定义）
2. 级细节：`converged`、n_newton、res_final、n_limited、n_floored、
   Mmax（全场）、wall_s
3. 袋监视（B25 f1 定义）：走廊 corrM / n_sup / 峰位 (x,z,q)
4. 峰位归因（b23 topk 机器）：dist_fus_surface、z−z_junc、
   dist_te_junc——区分"交界条带峰"与"翼面激波峰"（跨声速下走廊
   合法含翼面激波，判据看**位点**不看有无）
5. cl_p（收敛级）vs conforming 锚（cl_vs_mach.csv：0.2173@0.5、
   0.2321@0.65、0.2579@0.79、0.2617@0.84）
6. 度量 8（B25 附录 A3）：条带 aux |jump| max / 对 γ 比值——高 M 下
   锚定病监视
7. census：n_cut、n_aux_symmetry、n_te_nodes 断言（76/150）
8. sliver：条带最小二面角 / 体积分位数
9. 失效分类（预注册）：(a) strict 门 nlim/nflr>0 且峰在交界位点；
   (b) Newton 不收敛（res 停滞 / revert 耗尽）；(c) dm 减半耗尽

## 3. 判定（预注册）

- **B26-A（天花板抬升）**：C medium `m_last_converged` ≥ 0.60（爬过
  B18 死亡点 0.50 至少两级且 warm-start 链完整），收敛级无 (a) 类失
  效，cl_p 与 conforming 锚同趋势 → **袋是天花板限制器成立** →
  B18 demo 刷新 + GB9.4/GB20.5 重设 + Track B 文档收尾。
- **B26-B（天花板不动）**：C medium 死于与 A 同级且同 (a) 类特征
  （nlim/nflr 在交界位点）→ 袋非（唯一）限制器；按 G1.6 刻面
  （GB20.5，P11 路线）/ wake-LS 高 M 条件数 / upwind-条带交互排序
  重新归因，单独报告。
- **B26-C（新失效模式）**：C 爬过 A 死亡点但以 (b)/(c) 类或新位点
  死亡 → 记录失效分类与位点，单独归因（不在本 phase 判 C-A 延伸）。

## 4. 风险登记

| # | 风险 | 概率 | 影响 | 缓解 |
|---|---|---|---|---|
| T1 | A 重测与 committed 锚背离（B21/B22 freeze 修复的时滞效应） | 中 | 中 | 2.2 复现门条款：A 重测为主对照，背离记独立发现 |
| T2 | C 爬升远 → medium 成本失控（7 级 × ~15–25 min ≈ 2–3 h） | 高 | 低 | coarse 先行冒烟；medium C 后台跑；n_newton_max=80 不变 |
| T3 | 高 M 下 upwind/ν 与新条带交互（条带贴体区 ν 活性） | 中 | 中 | 度量 6/8 + n_nu_active 记录；失效分类 (a)/(b) 区分 |
| T4 | G1.6 刻面误差独立于袋封住天花板（GB20.5 归因份额） | 中 | 高 | B26-B 出口已预注册；P11 区隔见 §1 |
| T5 | seed Picard 在 m_start=0.50 行为漂移（n_seed=30） | 低 | 低 | B25 F1 同点已收敛干净（nlim=nflr=0），作先验 |

## 5. 工程纪律

- 线程：NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16。
- 默认路径 bit-identical（A 侧即 B18 配方当前码复测）；b1/b2/m2/v0
  绿为开跑前置（B25 已绿，本 phase 零库改动，跑前复核一次即可）。
- 结论以 committed CSV/PNG 为准；npz/msh gitignored；日志本地。
- 不 push；push/PR 先获用户确认。
- 重型运行同 B18 惯例显式启动（参照 PYFP3D_TRANSONIC_GATES 精神）。

## 6. 时间与产出

- 实现（wb26.py + run_g1.py，复用 wb25/wb_common/measure 机器）：~0.5 天
- 求解：coarse A/C ~0.5 h；medium A ~15 min（早死型，GB20.5 实测
  783–807 s）；medium C 0.5–3 h（取决于爬升距离，T2）
- 分析 + VERDICT：~0.5 天
- 产出：本文件、`wb26.py`、`run_g1.py`、`results/g1_summary.csv`、
  `results/g1_ceiling.png`、`VERDICT.md`

## 7. 与其他路线的关系（备忘）

- B26-A 成立：B18 demo 刷新（LS 腿复活）→ GB9.4/GB20.5 重设 → LS
  翼身全包线评估才真正开题（B25 §7 的原话落地）；conforming 退居
  交叉验证；Track V 片拓扑前置全部就位。
- B26-B/C：残余路线 = G1.6/P11（几何）、wake-LS 条件数（C1 遗产）、
  B10 roll-up 类（SHELVED，独立重大决策）；交付维持 conforming。
- 依赖：B25 的 `inboard_clip`（已入库）；不阻塞 P11，两者机制相邻
  但判据正交（B26 测操作后果，P11 修几何根因）。
