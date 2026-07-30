# GS1b.8：补场残差闸门 + G4.1 改走 Newton 路径 —— 预注册

轮次文件名：`20260730-0400-s1b-residual-gate-and-g41.md`
阶段 / 门：`S1b` / `GS1b.8`
分支：`claude/phase-two-s0-s1`
裁决依据：使用者 2026-07-30 —— "改成走 Newton 路径求解; 补残差闸门"

> **§1–§3 在动任何库代码之前写下并单独提交**（roadmap 原则 8）。§4 起执行后填。

---

## 1. 目的与前情

GS1b.7 测出：`solve_transonic_lifting` 的终态在 coarse M0.80 上 `|R| = 2.20e-04`，
比 Newton 的 5.46e-12 高 **8 个数量级**；Newton 从它出发 6 步就走掉；
Newton 的残差在它身上就是 2.198e-04 ⇒ **它不是同一套离散方程的解**。
而它报 `converged = True`。

★ **要说清一件事**：这不是"有人忘了"。`continuation.py` 里那个定义是**有名字的** ——
"**P4 engineering-converged regime**"（物理 M_max + Kutta 外循环 < `tol_gamma`），
`test_p4_transonic.py` 的 docstring 也写着。**定义是公开的；定义的代价没有被量过** ——
没人测过"这个态离真解有多远"（8 个数量级）以及"激波位置因此偏多少"（0.054 c，
是 M1 判据 (a) 容差 ±0.02 的 2.7 倍）。本轮补的是**代价的可见性**，不是补一个疏漏。

## 2. 要做的两件（按裁决）

### 2.1 补场残差闸门

`solve_transonic_lifting` 新增 `tol_residual`，并把上报拆成**两个不会互相冒充**的标志：

| 字段 | 含义 |
|---|---|
| `converged` | **真收敛**：physical **且** Kutta < `tol_gamma` **且** `|R| < tol_residual` |
| `engineering_converged` | 旧语义（physical 且 Kutta），**显式命名**，供历史对照 |
| `residual_final` / `not_converged_reason` | 达到的场残差，以及没过的是哪一条 |

默认 `tol_residual` 取一个**真的收敛要求**（1e-8）。这会让跨声速调用如实报
`converged=False` —— 那正是事实；把默认放宽到平台之上（例如 2e-4）等于把这一轮
要暴露的东西再藏一次。

### 2.2 G4.1 改走 Newton 路径

`tests/test_p4_transonic.py` 的 G4.1 改用 **Newton 路径**求解（收敛到 1e-12、
异构种子可复现、不钳制 —— GS1b.7 的三条证据）。Picard 的旧断言**不删**，改成
**显式标注"Picard 路径的历史回归，不是物理门"**并断言 `engineering_converged`。

★ **参考带不动**（判据 E5）。已知的后果，先写下来：Newton 的**等熵**收敛值是
**0.6581**，落在 `0.62 ± 0.03` **之外**（而熵修正 ON 是 **0.6186**，带内）。
所以按测量，G4.1 在**当前默认（等熵）**下会**不过**。处理方式沿用项目既有先例
（G1.6 的 strict xfail / M1 保留为记录在案的 FAIL）：**标为 strict xfail 并在
reason 里写明"熵修正 ON 时进带（0.6186）"** —— 不移动靶子，也不让套件带一个红。

## 3. 失败判据

| # | 判据 | 判定 |
|---|---|---|
| **E1** | 闸门真的会触发：跨声速 Picard 调用报 `converged=False`，且 `residual_final` 与 `not_converged_reason` 在返回值里 | 不触发 ⇒ 闸门是装饰 |
| **E2** | **健康路径不被打断**：亚声速 Picard 调用仍报 `converged=True`（A1 的 M0.5 四处断言） | 亚声速也红 ⇒ 闸门定得太严，重定 |
| **E3** | G4.1 走 Newton 后的读数**如实**对参考带定性（xfail 还是 pass 由测量决定） | 为了让它过而改带或改容差 ⇒ 整轮作废 |
| **E4** | **影响面逐个归因**：全套 + 门控套件跑一遍，每一个变红的都写明原因 | 只说"应该没影响" ⇒ 不算完成 |
| **E5** | 不改任何参考带、容差、工况、网格 | 改了 ⇒ 作废 |
| **E6** | 重型 demo 的 checklist 会翻（P4 的 G4.1/G4.3 行）—— **本轮只记录不重跑**（成本：~40 min），并在 demo 文档里留 erratum | 静默 ⇒ 不算完成 |

## 4. 改了什么

（执行后填）

## 5. 结果

（执行后填）

## 6. 判定与下一步

（执行后填）
