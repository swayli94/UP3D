# 20261012-2135 — **预注册：第 3 轮 = hex 族 NACA 能力边界门 D15（S2，星标工程项之二）+ 残余拆分探针**

**性质**：预注册，先于代码提交。**本轮建新门 + 一个时间盒探针；不改库、不改既有判据。**
**由来**：能力边界 §9 ★ 移交项第一条（hex 族 NACA 能力边界门，「缺陷仍在」式，
红 = 好消息）+ phase 7 立项 §2 S2 / §3 轮 3。归因已完成（p6 第 47 轮）：
hex 族收敛到 tet 从未出现的上翼面中部气泡/双激波结构，cl −24 %
（`phases/p6/bench/studies/hex_cl_attribution/`）；残余 = 生成器缺陷 vs 熵冻结交互
未分开（hex 各腿 sigma 冻结搅动 |dσ| 0.19–0.21 vs tet 0.067）。

---

## 1. 目的与判据（事前写死）

### 1.1 门本体（D15，`tests/D/test_D15_hex_naca_boundary.py`）

形状照抄 D14 族：**断言「缺陷仍在」，红 = 好消息**（修好了才红，红了走纪律 10/11
勘误，不许改回去换绿）。工况 M0.80/α1.25（M1 靶子），配方从
`tests/D/test_D05_euler_naca0012.py` **导入 RECIPE**（不重打常数，D14 同款做法）。

腿：默认套件跑 **hex coarse + hex medium + tet medium 锚**；**hex fine 进
`PYFP3D_TRANSONIC_GATES=1` 门控腿**（它撞 400 步帽、是最贵的一条）。
网格全是已提交的（`cases/meshes/naca0012_hex_2.5d/` 三级 .msh committed；
tet = `naca0012_2.5d`）⇒ 无网格实现噪声。

| # | 判据（全部对**新鲜求解**断言，数值期望来自已提交研究） | G-DOMAIN |
|---|---|---|
| A1 | **气泡仍在**：hex medium 上翼面 [0.2,0.55] max Cp > 0（研究实测 +0.199 @x=0.457） | 相反 = 气泡消失 ⇒ **红 = 好消息**，处置写进断言消息（重新定型为正常门） |
| A2 | **tet 对照无气泡**：tet medium 同窗 max Cp < 0（研究实测 −0.891） | 若 tet 也长出气泡 ⇒ 现象不再是 hex 专属 ⇒ 红，重新归因 |
| A3 | **升力亏损仍在**：cl_tet_medium − cl_hex_medium > 0.05（研究实测 0.083） | 亏损闭合 ⇒ 红 = 好消息，走勘误 |
| A4 | **指纹随加细增强**：hex 气泡计 coarse < medium（−0.297 → +0.199）；门控腿再加 medium < fine（+0.199 → +0.383，fine 按 physical 口径 \|R\| ≤ 1e-6 读，不断言 converged 旗标——旗标是线程依赖噪声，D05 教训） | 指纹变平/反转 ⇒ 结构性质变了 ⇒ 红，重读 |
| A5 | **激波位置两族一致**：\|x_shock_hex_medium − x_shock_tet_medium\| ≤ 0.02（研究实测 0.010 c）—— 把缺陷约束在 cl/结构而不是激波位置上 | 超带 ⇒ 缺陷形状变了 ⇒ 红，重读 |
| A6 | **亏损定位**：分段 ∫(cp_l−cp_u)dx 差里，front+shock 两段占 \|总差\| ≥ 90 %，LE/TE 各占 ≤ 5 %（研究实测：front −0.0341、shock −0.0509、LE −0.0003、TE +0.0026，总 −0.0827；段边界 [0,.05]/[.05,.45]/[.45,.75]/[.75,1] 事前写死，与研究同） | 亏损搬到 LE/TE ⇒ 归因形状变了 ⇒ 红，重读 |
| A7 | **证据锁**：`cases/gates/D15_hex_naca_boundary/summary.csv` 进 git，`assert_matches_committed`（rel_tol 1e-6）锁住。**锁哪些列/腿由 8t vs 16t 的 A/B 实测决定**（D14 做法：只锁跨线程稳定的，不稳的写进 CSV 但不进值锁）；`residual` 不进值锁（D14 已登记的教训）；`blas_threads` 出处列进锁 | 若 A/B 实测全部不稳 ⇒ 只锁 provenance 列，物理断言照跑，如实登记 |
| A8 | **G-TEETH 三条实测**：T1 改 committed CSV 一个被锁值 ⇒ 锁红；T2 把 tet medium 的新鲜数据喂给 A1/A3 的判别函数 ⇒ 判别返回「缺陷不存在」（证明判据不是恒真）；T3 把 DEFICIT_MIN 临时抬到 0.20 ⇒ A3 红。扰动全部还原，实测结果进判定文件 | 任何一条不红 ⇒ 该判据不算交付（G-TEETH 常设要求） |
| A9 | **G-CADENCE**：docstring 写明——默认套件跑 A1–A3、A5、A6、A7（hex coarse/medium + tet 锚）与证据锁；门控全集加 A4 的 fine 腿与导图腿 | — |
| A10 | `cases/gates/INDEX.md` 加 D15 条目；`tests/D` 门号占用 = D15（实测下一个空号） | — |

### 1.2 顺带探针（时间盒：**两条腿**，各一次 hex medium 求解）

目的：尝试拆开残余（生成器缺陷 vs 熵冻结交互）。两臂，只读气泡计/cl/x_shock：

- **臂 E0**：`entropy_correction=False`（其余同 RECIPE）。
- **臂 F0**：`freeze_tol=None`（关闭冻结，其余同 RECIPE）。

**G-DOMAIN（看数之前写死）**：

| 结果 | 读法 |
|---|---|
| E0 气泡消失（≤ 0）且亏损闭合 | 熵修正/冻结链路被牵连 ⇒ 残余偏向求解器侧 |
| E0 气泡仍在且亏损仍在 | 缺陷不依赖熵修正 ⇒ 残余偏向生成器/拓扑侧 |
| F0 气泡消失 | 冻结交互（锁住错误结构）被单独牵连 |
| 臂不收敛 / 气泡移动但仍在 / 两臂答案矛盾 | **RECORDED，残余保持未分开**，不升级仪器 |

**时间盒 kill criterion**：两臂总计不超过 2 次求解；跑完若分不开 ⇒ 照章 RECORDED，
**不许**为分开它而加第三臂（那需要动配方的探针轮，另立轮次）。
探针证据落 `bench/studies/hex_residual_split/`（新一次性研究：run.py + README + CSV），
并登记进 `bench/studies/README.md` 的一次性名单（保持归档二分的语义）。

### 1.3 不做

- 不修 hex 生成器、不动熵修正/冻结的任何默认值（研究项 S3/S4 之前不改行为）。
- 不碰 tet 族任何既有门（D05/D03/M1）的判据与证据；tet medium 锚腿的新鲜值必须与
  committed D05 值五位相合（0.341340），否则停下归因（锚断了）。
- 不给 hex 族「修好之后」的门定型（那是 A1 变红之后的事，断言消息里写明）。

## 2. 收口

- 默认全套：831 + Δ 闭合（D15 新增条目逐条记账；fine 腿与导图腿在门控全集）。
- 证据生成/验证在**暖缓存态**（纪律 11）；证据 CSV 带 `blas_threads` 出处列。
- 线程 A/B：8t vs 16t 各跑一遍三条默认腿，结果决定 A7 的锁名单并进判定文件。
