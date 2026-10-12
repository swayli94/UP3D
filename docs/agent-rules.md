# pyFP3D Agent 工作规则

> 本文件经 `CLAUDE.md` 的 `@docs/agent-rules.md` 注入每个 session。只放**现行**规则与指针；
> 阶段前史已随归档移走（phase 1 叙事在 `phases/p1/docs/`，phase 6 收口时的旧版全文在
> `phases/p6/docs/agent-rules.md`——历史快照，勿据其行事）。

## 当前阶段

**phase 6（全面梳理 + 门重建）已收口归档（2026-10-12）；下一阶段尚未开工，开工内容 = 使用者裁决。**

- **现行能力声明**（实测边界，引用能力时以它为准）：
  `phases/p6/docs/dev_phase_six/PHASE_SIX_CAPABILITY_BOUNDARY.md`（2026-10-11）。
- **套件基线**（2026-10-12 实测）：未门控 **825 passed / 0 failed / 33 skipped / 2 xfailed**；
  门控全集 812/21/3/1（1 XPASS）；快层 5/5。**BLAS 由库运行时钉成单线程**（读库不读 env）。
- 下一阶段候选工作项登记在 `docs/overview.md` 的「下一阶段」节（全部未排期）。
- 阶段工作面（轮次台账、prereg/verdict）的最近实例：`phases/p6/docs/dev_phase_six/`；
  轮次文件模板：`phases/p2/docs/dev_phase_two/_TEMPLATE.md`。

★ 本块只写指针与实测基线——会烂的散文一律进当阶段的 progress.md，不进这里。

## 工作纪律（每条都有实测出处；完整叙事见归档）

1. **线程上限 16，含 BLAS/OMP**：`NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
   OPENBLAS_NUM_THREADS=16`；漏掉 BLAS 上限在这台 16C/32T 合租机上代价 ~33%（P8 记录）。
   有重活在跑就先查负载、降到 8 线程。★ phase 6 起 **BLAS 由 `pyfp3d` 运行时钉单线程**
   （H30；pytest 里设 env 不生效——OpenBLAS 只在 numpy 首载时读一次）。
2. **不重算贵重已提交产物**（P4 heavy demo ~40 min、P5 medium 45–75 min、conforming fine
   M0.5 ~34 min 等）。已提交的 CSV/PNG 是权威；只有当真改动会移动那些数字、且你会提交
   刷新时才重跑。例行修改在便宜的 coarse 路径上验证。
3. **证据需要已提交的产物**——只存在于 .md 散文里的数字不是证据（2026-07-13 审计）。
   跑不起第二次的计算，正是要提交它 CSV 的理由。门证据写到被跟踪的 `cases/gates/<gate_id>/`，
   并由 `tests/_gate_evidence.py` 的新鲜复算锁守（gitignored 的目录里的东西不是证据）。
4. **折叠区纪律**（NACA medium M≈0.79，dcl/dM ≈ 6–10）：只用单网格回归锁，不做网格收敛
   声明；保留 Mach 爬坡（G10.3 判定）；折叠附近**忌**松的中间容差（G10.2b）。
5. **M6 fine（~45 万自由度）conforming Newton 配方**：`precond="amg"` + 紧 EW forcing
   （η=1e-8）+ `m_start=0.30, n_picard_seed=12`——medium 配方的 `direct` 是 4h39m/26GB 的
   splu 陷阱（P9；24 GB RSS 是标志）。★ phase 6 边界：amg 在**带激波**的腿上会失效
   （C05 喷管腿停在 |R|~3.9e-6，机制已归因：AMG 层级建在不含激波耦合的 Picard 块上）——
   含激波腿的预条件子选择见能力边界 §9 登记。
6. **门控测试**：重型跨声速/Newton/M6 门只在 `PYFP3D_TRANSONIC_GATES=1` 下跑；M6 `.msh`
   被 gitignore（用 `cases/meshes/onera_m6/generate_onera_m6.py` 重新生成，~30 s）。
7. **不重开已关闭的负结果**（没有新证据时）：G1.6 的修复路线（h-加密/恢复微调/Nitsche/
   边界数据修正）已排除；B8 的约束侧翼尖修法实测死路；展向 Γ 光滑化是死路；
   「0.019 缺口是分辨率」保持 *strongly indicated, NOT earned*。完整清单见
   `phases/p6/docs/dev_phase_six/20260824-0300-carried-forward.md`。
8. **两条路径 ⇒ 回移检查**（A3，2026-07-18 实测）：落在一条求解路径上的修复必须显式对照
   另一条并记录答案（level-set 线已在 phase 3 删除，但纪律对任何未来的第二路线成立）。
   同源教训：测试的界不得比它声称锁的门文更松——B3 的测试曾用 2% 断言一条 0.3% 的门（7×），
   直到 A3 实测 0.1441% 并收紧。
9. **收口刷新是五个面，不是两个**：阶段台账、agent-rules 当前阶段块 + 基线、
   overview.md、**PROJECT_STRUCTURE.md（页脚和目录树都算）**、`phases/README.md` 的阶段行。
   2026-07-17 审计发现的 17 处一致性缺陷，大多数就在最后两个面。
10. **再基线需要勘误清单，不只是新证据**（B22）：重新生成证据移动了已提交数字时，grep
    被移动的值横跨文档面逐处更正/标注，**同一个提交里**完成；五面流程只管新章节。
11. **逐位 A/B ⇒ 双腿同一 numba 缓存模式**（GV6.1 实测）：numba cache-load 与
    fresh-compile 在 `pyfp3d/viscous/` 上不逐位一致。★ phase 6 加码（第 43/44 轮）：
    **证据生成与证据验证必须在暖缓存态**——IBL 松耦合门（D08/D13）的答案随缓存出处态漂移
    1e-3…3.6e-2；重编译/清 `__pycache__` 后先白跑一遍再刷新或验证证据；新 clone 的第一跑
    会飘、第二跑绿——这不是缺陷，是编译扰动测量（依据：
    `phases/p6/docs/dev_phase_six/20261011-1800-process-history-determinism.md`）。
12. **内核/装配改动后先跑 `pytest tests/A/test_A01_freestream_preservation.py`**，
    再跑更大范围；全套 `pytest tests/` 的基线见顶部「当前阶段」块。

## 轮次制度（phase 2 起沿用，phase 6 定型）

1. **预注册在代码之前提交**：判据逐条 + kill criterion **事前写死**，每条判据写
   **G-DOMAIN**（相反结果落在哪、窗口是否覆盖定义域）。
2. **G-TEETH**：每条门必须实测它在被保护量坏掉时会红——「一条没有被验证过会红的门，
   不算交付」。**G-CADENCE**：门必须说明跑在哪个周期（未门控套件 / 快层 / 门控全集）。
3. **判定只有 PASS / FAIL / RECORDED / UNDEFINED**（收敛腿 < 2 报 UNDEFINED，不报方向）。
   **kill criterion 照章**：触发时按预注册处置，**不许放宽带去换绿**（UNDEFINED 是合格答案）。
4. **台账一行**进当阶段的 `progress.md`（每轮一行，新行在表上部）；测试账目必须闭合
   （如「829+2=831」）。
5. 预注册被实施推翻要**留勘误块**；再基线走纪律 10 的勘误清单。
6. 审计结论必须独立复核——不把另一个 agent 的报告当面值收下，逐条复测。

## 归档纪律（细则见 [../phases/README.md](../phases/README.md)）

- 阶段结束后，该阶段的轮次文档、一次性研究脚本、阶段性数据 `git mv` 进 `phases/pN/`
  （**跟踪，不 gitignore**——否则新 clone 里每个旧结论都失去提交产物）。
- **测试不归档**：还在检验活代码的测试留在 `tests/`；把活的门指进归档比没有门更糟。
- 被活代码/活测试引用的东西留在原位（判据 = 传递闭包，注意路径可以是逐段拼出来的）。
- 搬迁必修两类断裂（根相对路径字面量 + 相对 markdown 链接），顺序：先取基线 →
  `git mv` → 证明仪器有牙 → 改写 → 复验回基线。F07（`tests/F/test_F07_repo_invariants.py`）
  会机械断言活代码里内嵌的仓库路径必须存在。
