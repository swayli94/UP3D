# 20261104-0640 — **判定（第 11 轮）：迁移机 12 红的四臂出处测量 —— 不是单一根因；8 门与 numba JIT 代码生成同量级敏感，C07 不是，C08/D08/D13 另随缓存态移动**

轮次文件名：`20261104-0640-baseline-provenance-verdict.md`
阶段 / 门：phase 7 / 无门（调查）—— devflow 采用审计阻塞 1
登记：[20261104-0524-devflow-adoption-audit.md](20261104-0524-devflow-adoption-audit.md) §3.4（ISA 臂的冷/冷、暖/暖配对）、§12（四臂顺序）
测量代码：`080f291`（detached worktree）；证据：`bench/studies/baseline_provenance/`（仪器 + 四臂原始产物 + `classify.py` 与其三份 CSV）

---

## 1. 目的

迁移机上默认全套有 12 条红（审计 §3.3），全部是证据锁 / 散布锁 / 一条缺陷记录断言，库与测试零改动。要回答的是事先登记的两个问题：

1. 答案对 **numba JIT 代码生成**（ISA）敏感吗？敏感的量级够不够解释与 committed 证据的缺口？
2. **冷/暖缓存态**各自贡献多少？

★ **登记状态如实写**：四臂与四组配对在启动（05:16）前已写进审计报告，但报告当时在 gitignored 的 `bench/results/`，提交（`12f6e19`）晚于启动；
**每门的汇总统计与分类标签是事后描述**（`classify.py` 的 `label_post_hoc` 列），没有预注册阈值，只用来描述，**不决定处置**。

## 2. 改了什么

库 / 测试 / 默认值 / 参考数据 / committed gate 证据**零改动**。新增 `bench/studies/baseline_provenance/`（只记录不改判定的 pytest 插件、runner、比较与分类脚本、四臂产物）。

## 3. 怎么测的

- 11 个门文件（C05 C07 C08 D05 D06 D08 D11 D12 D13 D14 D15），`080f291`，独立 detached worktree，16 线程，BLAS 运行时 1 线程，`-p no:cacheprovider`。
- 四臂：`native_cold`（空 `NUMBA_CACHE_DIR`）→ `native_warm`（同目录第二跑）→ `generic_cold`（`NUMBA_CPU_NAME=generic`，另一空目录；feature 串随之为空）→ `generic_warm`。墙钟 1009 / 1045 / 1066 / 927 s。
- 比较对象：committed `summary.csv`（旧开发机生成）；参照全部是**防漂移的历史值**，不是外部真值。
- 范围：`NUMBA_CPU_NAME` 只改 numba/LLVM JIT，不改 OpenBLAS / SciPy / PyAMG 等预编译库。旧机的 codegen 从未被记录，所以 generic 臂只能检验**敏感性**，不能复现旧机。

## 4. 结果

**失败节点**（`compare_summary.md`）：`native_warm` 恰为原 12 项。`native_cold` 11（D08 绿），`generic_cold` 13（D08 红，另加 C08 `test_m_crit_is_the_cheaper_lever` 红），`generic_warm` 11（D08 绿）。

**逐门**（`classification.csv`；gap = native 暖对 committed 的相对差；ISA = generic 暖对 native 暖；cache = native 冷对 native 暖）：

| 门 | 值数 | gap > 1e-6 | gap 最大 | ISA 最大 | cache 最大 | 事后描述 |
|---|---|---|---|---|---|---|
| C05 | 24 | 13 | 2.60e-02 | 1.98e-01 | 3.6e-12 | 与 codegen 同量级敏感 |
| C07 | 24 | 8 | 4.53e-04 | **6.66e-09** | 0 | **numba JIT 解释不了** |
| C08 | 30 | 4 | 5.44e-04 | 1.46e-03 | **4.27e-04** | 随缓存态移动 + codegen |
| D05 | 108 | 7 | 1.02e-04 | 4.28e-05 | 5.5e-10 | 与 codegen 同量级敏感 |
| D06 | 36 | 15 | 3.55e-02 | 2.01e+01 | 1.0e-09 | 与 codegen 敏感（见 F5） |
| D08 | 28 | 20 | 2.22e-01 | 1.42e-01 | **2.65e-01** | 随缓存态移动 + codegen |
| D11 | 36 | 5 | 1.47e-04 | 6.05e-05 | 3.0e-11 | 与 codegen 同量级敏感 |
| D12 | 28 | 14 | 9.03e-02 | 2.43e+01 | 1.2e-09 | 与 codegen 敏感（见 F5） |
| D14 | 22 | 10 | 2.62e-04 | 2.87e-04 | 3.4e-10 | 与 codegen 同量级敏感 |
| D15 | 15 | 1 | 1.14e-06 | 4.89e-06 | 0 | 与 codegen 同量级敏感 |
| D13 asym（npz） | — | — | — | 4.73e-02（abs） | **9.13e-02**（abs） | 四臂 coarse 最大不对称度 1.294 / 1.202 / 1.152 / 1.155，限 1.15 |

★ D05/D06/D14/D15 的门在调用比对函数前已把值舍入到 10 位有效数字，这四门的 “0” 意为 10 位内相同。

**发现：**

- **F1 当前机暖态可复现，且与缓存的编译历史无关（显示精度内）。** 早上的定向重跑用默认 `__pycache__`，本轮 `native_warm` 用全新缓存目录：10 个可对照值 `.9e` 全等，D13 不对称度 `repr` 全等（`1.2022292349421009`）；1 个（D13 `dstar`）因键格式未对上（`premise_b.csv`）。⇒ 裁决 A 前提 (b) 得到**显示精度级**支持；全列位级仍未证。
- **F2 冷/暖缓存对 8 门不是原因。** C05 C07 D05 D06 D11 D12 D14 D15 的 native 冷对暖 ≤ 3.6e-12（舍入门 ≤ 1.0e-09），它们的红不是冷编译扰动。
- **F3 numba JIT 代码生成足以在同量级上移动其中 7 门**（C05 D05 D11 D14 D15 及 D06/D12）：ISA 扰动最大值是 gap 最大值的 0.4×–7.6×（D06/D12 见 F5）。这与“跨机 codegen 差异”假设**一致**，但**不是证明**：旧机的 codegen 不可知。
- **F4 C07 不是 JIT 代码生成。** 差值只在 medium 一层（8 列，gap 2e-6…4.5e-4），而 generic 与 native 只差 ≤ 6.7e-9。候选转向 JIT 之外：预编译库（OpenBLAS 按 CPU 选 kernel、SciPy/PyAMG）、包或 Python 补丁版本。
- **F5 D06/D12 共用的 RAE2822 coarse M0.73 α3.19 腿站在刀刃上。** cl：committed 0.7438、native 0.7652（冷暖相同）、generic 冷 0.7279、**generic 暖 2.0185**（不物理，另一个根）。同一份代码，换 codegen + 换缓存态就跳根；而且 native 与 committed 之间就差 **2.9 %**（原日志只报第一个越界列，看不出来）。这是 phase 7 G-ATTRACTOR 意义上的发现，登记，不在本轮处置。
- **F6 D08 的绿/红随缓存态翻转。** `xtr030|medium` 的 `cp_rms_lower_sst` 对 committed 的绝对差：native 冷 3.98e-4（绿）、native 暖 5.43e-3（红）、generic 冷 5.23e-3（红）、generic 暖 2.74e-3（绿），带 3.89e-3。当前机实测的冷暖散布（相对 2.65e-1）大于旧机标定的散布带。
- **F7 D13 不对称度四臂全部 > 1.15**（1.152–1.294），自身跨状态变化 0.14，大于该量登记的散布带 0.056；历史 0.90–1.12 一次也没有回到。
- **F8 numba 统计（裁决 A 前提 (a) 的试点）。** 暖臂 61 个 dispatcher 上 **0 miss / 34 hit**，冷臂 **50 miss / 0 hit**（`numba_cache.csv`）；61 = 库内 `cache=True` 处数。冷臂多出的 16 次 miss 全是被调用方的辅助函数，随调用方一起编译，暖态下从不单独调度。⇒ 冷态信号可靠；但 **G-TEETH（强制测量进程内编译必须红）未做**，前提 (a) 未证。

**全套测试**：本轮未跑全套（只跑 11 个门文件 ×4）；当前机全套仍是审计 §3.3 的 12 failed / 827 passed / 35 skipped / 2 xfailed。

## 5. 判定

**RECORDED** —— 12 红**不是单一根因**：
- 8 门（C05 C07 D05 D06 D11 D12 D14 D15）与冷/暖缓存无关，其中 7 门的缺口与 numba JIT codegen 扰动同量级（**一致，未证明**），**C07 被 JIT 解释不了**；
- C08、D08、D13 还随缓存态移动，D08 在同一台机器上绿红翻转；
- D06/D12 的一条跨声速腿在 codegen 扰动下跳到不物理的根。

一句话：迁移机 12 红 = 跨平台出处差（7 门与 JIT codegen 同量级、C07 指向 JIT 之外）叠加 IBL/C08 的缓存态敏感；暖态可复现、与缓存编译历史无关（显示精度内）。

## 6. 下一步

- **处置需要使用者裁决**（审计 §3.4 第 3 条、纪律 10）：(i) 生成当前机的平台专用证据（证据 CSV 先补有效 codegen 出处列，带完整勘误清单）；(ii) 保持红色、登记为跨平台证据差；(iii) 其他。本轮不做任何再基线。adoption Change 在裁决前仍不 `start`。
- 便宜的补充臂（只登记，未跑）：C07 单文件 × `OPENBLAS_CORETYPE` 两档，检验 F4 的预编译库候选。
- F5（D06/D12 刀刃腿）是独立的 G-ATTRACTOR 研究问题，不在 adoption 线上处置。
- 裁决 A 前提 (a) 的 G-TEETH 与全列位级的前提 (b) 留给 adoption Change tasks 4.4。
- detached worktree `/home/lrz/code/up3d-baseline-prov` 可在本轮提交后 `git worktree remove`。
