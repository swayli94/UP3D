# 20261010-1430 — **预注册:M1 门预算 80 → 400 重跑,产品指标按实测更新(F06 锁随动)**

**性质**:预注册(判据与 kill criterion 事前写死),随代码提交。
**对象**:`bench/run_m1_gate.py`(预算、出处列、诚实的线程报告)、
`pyfp3d/__init__.py`(`blas_thread_count()` 提为公共 API)、
`tests/_gate_evidence.py`(改为委托)、`bench/gate_results/m1_gate_default{,_3seed}.csv`
(**重刷**)、`tests/F/test_F06_m1_spec.py`(锁的实测值随动)、
`docs/dev_phase_six/progress.md`(产品指标表)。
**约束(事前重申)**:M1 的**规格**不动 —— `M1B_TOL = M1C_TOL = 0.20`、种子 {0,5}、
C ∈ {1.0, 1.5, 3.0}、NACA0012 M0.80/α1.25、无 Mach 坡道;**不把** `bench/usability.py`
或 x_shock 带接进 M1(立项 §5.6,使用者裁决)。

---

## 1. 目的

H30 已归因 `n_newton_max = 80` 对 M1 靶子腿(M0.80/α1.25 medium)不够,答案门
(D03/D05/D06/D11/D12)早已改 400,而 **M1 门自己是产品指标门、断言的是答案**,
按 H30 预注册自己的划分标准属于「答案门」—— 但它仍跑 80。本轮把它对齐到 400
并重跑,回答:**M1c 在 seed 0 上的覆盖失败(medium 只有 C=1.5 一条收敛腿)有多少
是同一个预算缺陷。**

**预先写下的失败标准**:
- 若重跑后我不报告新旧 CSV 的逐腿对照就更新 F06 的锁 ⇒ 失败
  (产品指标的移动必须有逐腿证据,不许"绿了就行")。
- 若趁机动 `M1B_TOL`/`M1C_TOL`/种子/工况 ⇒ 失败(规格是裁决,不是这轮的对象)。

## 2. 改动内容(事前登记)

1. `bench/run_m1_gate.py:187`:`n_newton_max` 80 → **400**(与答案门一致;
   `tests/D/test_D05:73` 与 `test_D11:62` 的「与 bench/run_m1_gate.py 逐字一致」
   注释**由此恢复为真** —— 顺带核销一处文档失配)。
2. 同文件:CSV 行加 **`blas_threads`**(读库);开头的线程报告行除 env 值外
   **加印库实际 BLAS 线程数** —— 「环境变量不是库的状态」,这份脚本恰恰是把
   per-leg 墙钟写进证据的脚本,它的出处行必须是真话。
3. `blas_thread_count()` 从 `tests/_gate_evidence.py` 提为
   `pyfp3d/__init__.py` 的**公共 API**(`_gate_evidence` 改为委托,六门 import 不变);
   bench 脚本从库取,不走 bench←tests 依赖。
4. 重跑 **default(2 seeds)与 --third-seed(3 seeds)**,覆盖两个已提交 CSV
   (再基线;旧 CSV 已存 `/tmp` 备份,勘误清单见 §4-C4)。
5. `tests/F/test_F06_m1_spec.py` 的锁**按新实测更新**:
   `test_measured_verdicts_on_the_committed_data` 的逐种子判定、
   `test_the_margins_are_recorded` 的两个余量(0.178563 / 0.16333)—— 新值必须是
   **从新 CSV 经 `_criteria` 实算**的,断言消息里写明「2026-10-10,预算 400 + BLAS=1」。

## 3. 判据(事前写死)

| # | 判据 | 期望 |
|---|---|---|
| C1 | 新 CSV 的出处列为真:`blas_threads` 全为 1,`n_threads` 语义不变 | 逐行检查 |
| C2 | 预算假设的回答:已提交 CSV 里 `n_newton=80` 的那些 medium 腿,在 400 下的 `converged`/`accept_reason` 逐腿列出 | 记录,不论方向 |
| C3 | `_criteria` 在新 CSV 上重算:M1b/M1c 逐种子判定;**M1c seed-0 是否仍是覆盖失败** | 记录;F06 锁按此更新 |
| C4 | 勘误:grep 旧 CSV 被引用的值(F06 两个余量、docstring 里的 M1 数字、progress.md 产品指标表),逐处更正或标注 | 清单闭合 |
| C5 | F06 单跑全绿(锁更新后);受影响面(A27/A30/A06 的机制测试**不经过** run_m1_gate,必须不受影响) | F06 绿;三者复跑绿 |
| C6 | 全套收口(标准纪律):无新增红 | 0 failed |

**G-DOMAIN**:C2 的相反结果 = 400 下仍不收敛 ⇒ 预算假设被否,M1c seed-0 的覆盖
失败需要**别的**解释 ⇒ 如实登记(这本身就是有价值的结果,不许为了避免它而把预算
再加大 —— 400 是答案门的既定值,不是可以随手拧的旋钮)。

**G-CADENCE**:run_m1_gate.py 不在测试收集里;它的锁是 F06(默认周期,零求解)。

## 4. kill criterion

- 若重跑中某条腿把**收敛腿变成了不收敛**(例如 seed-5 C=1.5 medium 在 400 下反而失败),
  那不是「预算不够」能解释的 ⇒ **停**,先归因(候选:400 步走进了另一个吸引子 ——
  与 D06/D08 一族),不许直接刷 CSV 了事。
- 若新 CSV 的 medium 腿数/种子结构与旧 CSV 对不上(腿丢了)⇒ 停,查脚本改动。

## 5. 明确不做

- 不动 M1 规格(容差、种子、C 集合、工况、无坡道)。
- 不接 usability / x_shock 带进 M1(裁决)。
- 不动 `m1_gate.csv` / `m1_gate_default_8t.csv` / `m1_gate_entropy.csv`(历史记录)。
- 不动 A27/A30/A06 的 80 步预算(那是机制测试的**测试对象**,不是答案)。
- 不重跑 M1a/M3 等其他产品指标。
