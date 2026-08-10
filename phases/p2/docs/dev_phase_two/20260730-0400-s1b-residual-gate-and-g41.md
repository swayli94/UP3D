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

| 文件 | 改动 |
|---|---|
| `pyfp3d/solve/continuation.py` | 新增 `tol_residual`（默认 **1e-8**）；把上报**拆成两个不会互相冒充的标志**：`converged`（physical **且** Kutta **且** `\|R\| < tol_residual`）与 `engineering_converged`（旧语义，显式命名）；新增 `residual_final` 与 `not_converged_reason`（后者直接写明"这是 Picard 的激波平台，此态不是离散方程的解，见 GS1b.7"） |
| `tests/test_p4_transonic.py` | ★ G4.1 **改走 Newton 路径**（新 `_transonic_case_newton`）；旧 Picard 断言**保留**但改名为 `test_p4_picard_path_historical_regression`，显式标注**不是物理门**，只断言 `engineering_converged` + 激波质量 + 一条对 0.6041 的**漂移锁**，并断言它**不是**真收敛（若哪天它真收敛了，就该把它升回物理门）；`_assert_g41` 不再要求 Picard 的 `converged` 语义（由调用方各自断言） |
| `cases/demo/p4_transonic/run_demo.py` | G4.1 的两行改读 `engineering_converged`，**新增一行记录真实场残差**；docstring 顶部加 **erratum**（committed 产物早于闸门，"converged"那一条已被取代，其余仍成立；重跑与 GS1b.6 的 B3 一起做） |

## 5. 结果

### 5.1 E1 = PASS：闸门真的触发，且说明原因

```
跨声速 M0.80：converged=False  engineering_converged=True  residual_final=2.198e-04
  reason: field residual 2.198e-04 >= tol_residual 1.0e-08 (the Picard shock
          plateau -- this state is NOT a solution of the discrete equations)
```

### 5.2 E2 = PASS：健康路径不被打断

```
亚声速 M0.50：converged=True  engineering_converged=True  residual_final=3.858e-09
```
★ 但要记一笔**余量很窄**：3.858e-09 对 1e-8 只有 **2.6 倍**。所以这个默认值不是
"随便挑一个宽的"——它刚好把亚声速放过、把跨声速平台挡住。若以后有亚声速算例落在
1e-8 之上，**应当先查那个算例**，而不是放宽闸门。

### 5.3 E3 = PASS：G4.1 走 Newton 后的定性由测量决定

| | x_shock | 参考带 0.62 ± 0.03 |
|---|---|---|
| Newton **等熵**（当前默认） | **0.6581** | **带外**（+0.038） |
| Newton **熵修正 ON** | **0.6186** | **带内**（−0.0014） |

⇒ 按测量，G4.1 在当前默认下**不过**，标为 **strict xfail**，reason 里写明"熵修正 ON
时进带（0.6186）"。**参考带一个字没动**（E5）。这与项目既有先例一致
（G1.6 的 strict xfail、M1 保留为记录在案的 FAIL）：不移动靶子，也不让套件带一个红。

★ 顺带一条值得记的读数：参考带本身含一个**有文献依据的 0–2% 弦长等熵后移余量**
（Holst PAS 2000），而 Newton 收敛态的实测后移是 **+0.043 c（4.3%）** ——
**比文献记的余量还大一倍**。也就是说我们这套离散（C = 1.5）的后移偏差超出那条带的设计范围。
这条独立于熵修正，登记为观察。

### 5.4 E4 = 影响面

**全套（非门控）：`676 passed / 28 skipped / 3 xfailed`，0 失败**
（本轮之前是 676 / 28 / **2**）。多出来的那个 xfailed 就是新的 Newton 路径 G4.1；
passed 数不变是因为一个通过的测试（旧 `test_g41_transonic_coarse_smoke`）被替换成
另一个通过的测试（Picard 历史回归）。⇒ **残差闸门没有打断任何非门控路径**
（含 A1 的四处 M0.5 `converged` 断言）。

**门控套件：`2 failed / 700 passed / 2 skipped / 3 xfailed`，6589 s。**逐个归因：

| 变红的 | 归因 |
|---|---|
| `test_b9_wingbody_conforming::test_laplace_lifting_loads_the_junction` | **既有失败**，`main` 上同样红（GS1b.6 §5.4 已归因），与本轮无关 |
| `test_p4_transonic::test_g43_robustness_sweep` | ★ **本轮闸门造成，且是预期后果**：G4.3 断言"10/10 全部收敛"，而 `converged` 现在要求场残差，跨声速点都在 Picard 的激波平台上（\|R\| ~ 2e-04） |

**G4.3 的处置**：它问的是"**驱动器**能不能在包线上跑下来"，不是"答案对不对"。
所以断言改成它**实际测到的** `engineering_converged`，并把**每个点的真实残差**写进
summary CSV（原来那个"收敛"是靠一个名字比含义强的标志暗示的）。答案质量的门在 Newton 路径上。
Newton 路径版的这个 sweep 登记为待办（还要 10 次耦合求解），本轮不做。
再跑：`tests/test_p4_transonic.py` 门控下 **3 passed / 1 xfailed**。

### 5.4.1 ★★ 我自己漏了一处，验证时才发现，而 medium 的读数是这条线最强的证据

改完 coarse 之后我以为完了。跑门控时看到 `test_g41_transonic_medium_gate` **通过**——
它**还在 Picard 路径上**断言那条参考带，**正是我刚刚批评的"在未收敛态上立物理门"**，
而且它"通过"只是因为 Picard 的 medium 激波恰好落在带内。已同样改到 Newton 路径。

改完之后 medium 的实测（Newton，M0.80，参考带 0.62 ± 0.03）：

| | converged | \|R\| | x_shock | 对带 |
|---|---|---|---|---|
| **等熵 OFF**（当前默认） | **False** | 6.9e-06 | 0.8819 | 带外 **+0.262** |
| **熵修正 ON** | **True** | **2.6e-13** | **0.6146** | **带内 −0.0054** |

⇒ **在 M1 的工况、medium 网格上：等熵的耦合 Newton 根本不收敛**
（与 M1 门记录的"medium：0 条腿收敛"一致），**而熵修正把它变成一个收敛解、并且落进参考带。**
这是熵修正到目前为止最强的单条证据 —— 它不只是把数字挪近，它让一个解不出来的工况解出来了。

★ 因此 medium 那条 xfail 与 coarse **原因不同**（是**不收敛**，不是带外），
xfail 的 reason 已按实测改写。两条都会在翻默认那天变成 pass。

### 5.4.2 summary writer 改成路径无关

`_write_g41_summary` 原来写死 `kutta_mismatch` / `n_picard_total`（Picard 才有的键），
Newton 的返回值没有 ⇒ 会 KeyError。改成"写返回值里实际有的键"，并加写
`residual_final` / `engineering_converged` / `n_newton`。

### 5.5 E6 = demo 的处置

P4 重型 demo（~40 min）**本轮不重跑**（预注册如此）。已做：脚本的 G4.1 两行改读
`engineering_converged`、新增一行记录真实残差、docstring 顶部写入 erratum
（committed 的 CSV/PNG 早于闸门："converged"那一条已被取代，其余仍成立；
重跑与 GS1b.6 的 B3 一起做）。

## 6. 判定与下一步

**两件裁决事项都已实现并测量；E1/E2/E3/E5/E6 达成，E4 的门控腿在跑。**

| 判据 | 结果 |
|---|---|
| **E1** 闸门触发 | **PASS**（跨声速 `converged=False` + `residual_final` + 具体原因） |
| **E2** 健康路径不断 | **PASS**（亚声速 `converged=True`，余量 2.6× —— 已记为需要注意的窄余量） |
| **E3** G4.1 定性由测量决定 | **PASS**（等熵 0.6581 带外 ⇒ strict xfail；熵修正 ON 0.6186 带内，写进 reason） |
| **E4** 影响面逐个归因 | **PASS** —— 非门控 **676 / 28 / 3，0 失败**；门控 **2 failed / 700 passed**，两个都已归因（1 个既有、1 个本轮闸门的预期后果，已按"驱动器健壮性"重新定性） |
| **E5** 不改带/容差/工况/网格 | 遵守 |
| **E6** demo 记录不重跑 | **PASS**（脚本改读 `engineering_converged` + 新增真实残差行 + docstring erratum） |

**这一轮的三条可引用结论**：

1. `solve_transonic_lifting` 现在**不会再把非解报成收敛**；旧语义仍可读，但必须叫它自己的名字
   （`engineering_converged`）—— 两个概念不再能互相冒充；
2. **G4.1 现在是一条建立在"真的是解"的状态上的物理门**，而它在当前默认（等熵）下
   **如实不过**（0.6581 带外），熵修正 ON 时进带（0.6186）；
3. ★ 参考带自带的等熵后移余量是 **0–2% 弦长**（Holst），而 Newton 收敛态的实测后移是
   **4.3%** —— **超出那条带的设计范围一倍**。这条独立于熵修正，是对"我们这套离散在 C = 1.5
   下的后移偏差"的一个新观察，登记备查。

**下一步**：

1. 门控腿跑完补齐 E4（每一个变红的都归因）；
2. 然后回到 **GS1b.6 的 B3**（重型 demo 重算），现在有了明确的基准：**只用 Newton 路径**，
   且 P4 的 committed 产物与 B3 一起重跑（E6 已埋好 erratum）；
3. 默认值的裁决输入已经齐了：包线内几乎不动（M1a 三级 +0.91% → +0.67% 仍过）、
   包线外方向正确且路径一致（≈ −0.04 c）、影响面 4 个锁、而**其中 G4.1 已经先一步
   改成了"ON 时才过"的形状** —— 也就是说翻默认会让 G4.1 从 xfail 变 pass，
   剩下三个 P8 锁需要与翻默认同一个 commit 重新锚定。


（执行后填）
