# 20261010-1100 — **预注册:H30 层 2 返工 —— 把 BLAS 固定做成运行时生效,F09 升级为读库判据**

**性质**:预注册(判据与 kill criterion 事前写死),随代码提交。
**对象**:`pyfp3d/__init__.py`(加运行时固定)、`tests/F/test_F09_blas_determinism.py`
(重写判据载体)、`bench/run_capability_locks.py`(快层显式走逃生舱)。
**由来**:[Kimi 独立审计](20261010-1038-progress-audit-kimi.md) §4.2 复测确认
[Claude 审计](20261009-2200-progress-audit.md) §4.2:env 固定只在 numpy 未导入时有效,
而 pytest 里 numpy 先导入 ⇒ **固定从未生效**;证据实际锁在**默认 24 档**(10-10 实测,
4/4 复现)。

---

## 1. 目的

**让「BLAS 单线程固定」在两种 import 顺序下都真实生效,并让 F09 能看见它失效。**

**预先写下的失败标准**:
- 运行时固定后,numpy 先导入的子进程里 `openblas_get_num_threads()` 仍 ≠ 1 ⇒ 本轮 FAIL
  (不许退回"只在 pyfp3d-first 顺序下有效"就宣布完成 —— 那正是现状)。
- F09 的新前提断言若读的是 env 而不是库 ⇒ 本轮 FAIL(同一族错误不得换个位置再犯)。

## 2. 改动内容(事前登记)

1. **`pyfp3d/__init__.py`**:保留现有 env 强制覆盖(它覆盖 numpy 未导入的情形 + 逃生舱语义),
   **追加运行时固定**:若 `numpy` 已在 `sys.modules`(即 env 已来不及),从 `/proc/self/maps`
   找到已加载的 `libopenblas*`(实测 numpy 2.4.6 conda 版经 `libopenblas.so.0` 链接),
   `ctypes` 调 `openblas_set_num_threads(1)`。**找不到已加载的库就 `warnings.warn`,不静默**
   (env 路径仍在,属于"尽力 + 响亮",不是"假装成功")。逃生舱语义不变:
   `PYFP3D_ALLOW_BLAS_THREADS=1` 时整个固定(含运行时)跳过。
   ★ 仍然**不 import numpy**(延迟初始化语义不变);numpy 先导入的情形由 maps 路径覆盖。
2. **`tests/F/test_F09` 重写**:
   - 探针子进程改为**两种 import 顺序**(`pyfp3d`-first 与 `numpy`-first)× 外部 BLAS {8, 16}
     = 4 条腿,载体不变(medium M0.80/α1.25,`n_newton_max=80`,实测 BLAS 敏感);
   - **前提断言改读库的实际线程数**(ctypes + `/proc/self/maps`,与库修复同一份小工具,
     探针内联,不进库 API)—— 不再读 `os.environ`;
   - 判据不变:同顺序内、跨外部 BLAS 档,`cl` 与 `resid` **逐位相同**;
   - 顺手修掉审计发现的自相矛盾注释(:38「强制覆盖」vs :108「setdefault」)。
3. **`bench/run_capability_locks.py`**:`os.environ.setdefault("PYFP3D_ALLOW_BLAS_THREADS", "1")`,
   注释写明理由 —— 快层的 G8.2 墙钟预算与位级锚是在 **BLAS=8** 标定的;运行时固定优先级
   高于 env,不设逃生舱会让快层静默滑到 BLAS=1(标定作废且归因困难)。
   **这是本轮的处置决策:快层保 8 线程、显式逃生舱,而不是接受 BLAS=1 重刷
   `capability_locks.csv`**(后者把一次仪器修复和一次证据再基线混在一个提交里,违反
   「一个提交一件事」;若将来要把快层也归一到 BLAS=1,单开一轮)。

## 3. 判据(事前写死)

| # | 判据 | 期望 |
|---|---|---|
| C1 | ctypes 探针:`numpy` 先导入 + 外部 `OPENBLAS_NUM_THREADS=16`,import pyfp3d 后**库实际线程数 == 1** | 修复前 16,修复后 1 |
| C2 | F09 四腿:同顺序内跨外部档 `cl`/`resid` 逐位相同;两顺序间也逐位相同(都是 BLAS=1) | 全等 |
| C3 | F09 前提断言:四腿的库实际线程数 == 1(读库,非 env) | 全 1 |
| C4 | **G-TEETH**:`PYFP3D_ALLOW_BLAS_THREADS=1` 跑 F09 ⇒ **必须红**(前提断言与判据都红) | 红 |
| C5 | 全套未门控:除已知 6 条证据锁外**无新增红**;6 条锁的失配量级预期**变小**(从"距 24 档"变成"距 1 档",10-09 实测 D05 ≈3.1e-05)—— **本轮不重刷证据**(下一轮,带纪律 11 勘误清单) | 无新增红 |
| C6 | 全套墙钟变化记录(固定生效后 BLAS=1,预期略涨;H30 记过四层 +4.5 %) | 只记录,不设带 |

**G-DOMAIN**:相反结果 = numpy 先导入时库线程数仍为外部值 ⇒ 运行时固定无效 ⇒ FAIL,
退回调查(CD 路径 / maps 解析),**不得**改判据去适应。

**G-CADENCE**:F09 在默认周期(未门控)跑;四腿成本 ≈ 2× 现状(F09 现为 2 腿)。

## 4. kill criterion

若 `/proc/self/maps` 路线在本机以外的环境不可用(非 Linux、或 maps 里找不到
libopenblas):**不伪造绿** —— 警告 + 退回 env-only(与现状等价,不更差),在本文件登记,
并给 F09 的前提断言标"该环境下不可验"。本轮在本机(24 核 Linux)必须 C1–C5 全过。

## 5. 明确不做

- 不重刷 D05/D06/D11/D12 证据(下一轮,走纪律 11)。
- 不动 D08/D13 的锁形状(等 D08/D13 的 BLAS=24 探针结果 + 使用者裁决)。
- 不动 M1 门预算(再下一轮)。
- 不改任何判据的数值容差。
