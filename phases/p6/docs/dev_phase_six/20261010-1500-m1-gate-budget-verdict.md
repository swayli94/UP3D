# 20261010-1500 — **判定:M1 门预算 80→400 重跑 —— 覆盖失败不是预算缺陷(假设被否),产品指标判定不变**

轮次文件名:`20261010-1500-m1-gate-budget-verdict.md`
阶段 / 门:无门(产品指标实测)—— 预注册:
[20261010-1430-m1-gate-budget-prereg.md](20261010-1430-m1-gate-budget-prereg.md)
提交:分支 `claude/phase-six-audit`(本文件随判定提交)

---

## 1. 目的

回答:**M1c 在 seed 0 上的覆盖失败(medium 只有 C=1.5 一条收敛腿)有多少是
`n_newton_max=80` 的预算缺陷。** 判据与 kill criterion 见预注册,事前写死,未改。

## 2. 改了什么

- `bench/run_m1_gate.py`:`n_newton_max` 80 → **400**(与答案门对齐;D05:73/D11:62
  的「逐字一致」注释由此恢复为真);CSV 行加 **`blas_threads`**(读库);
  线程报告行加印**库实际 BLAS 线程数**(「环境变量不是库的状态」)。
- `pyfp3d/__init__.py`:`blas_thread_count()` 提为**公共 API**
  (`_loaded_openblas_paths()` 公共化,运行时固定的探针改走它);
  `tests/_gate_evidence.py` 的同义函数改为委托(六门 import 不变)。
- `bench/gate_results/m1_gate_default{,_3seed}.csv`:重刷(预算 400 + BLAS=1)。
- `tests/F/test_F06_m1_spec.py`:**未改** —— 见 §4.1,锁的判定与余量全部在带内。
- `docs/dev_phase_six/progress.md`:产品指标表两行补「本阶段变化」。
- `docs/dev_phase_six/20261009-2200-progress-audit.md`:§4.4 加勘误标注
  (其「线索,未验证」已闭合)。

## 3. 怎么测的

```bash
python bench/run_m1_gate.py                # default,2 seeds × 3 C × 2 levels = 12 腿
python bench/run_m1_gate.py --third-seed   # + seed 12(诊断,不投票)= 18 腿
python -m pytest tests/F/test_F06_m1_spec.py tests/A/test_A27_gamma_target.py \
                 tests/A/test_A30_seed_fallback_contract.py tests/A/test_A06_entropy_correction.py
python -m pytest tests/ -q -rfEsx          # C6 收口全套(16 线程纪律)
```

## 4. 结果

### 4.1 判据逐条

| 判据 | 结果 |
|---|---|
| C1 出处列为真 | ✅ `blas_threads` 全为 1(新旧两 CSV 逐行);线程报告行同时印 env 与库实际 |
| C2 预算假设逐腿回答 | ✅ 见 §4.2 —— **假设被否** |
| C3 `_criteria` 新 CSV 重算 | ✅ seed 0:M1b PASS(−16.33 %)/ M1c FAIL-覆盖;seed 5:M1b PASS(−16.15 %)/ M1c PASS(17.86 %)—— **与已提交判定逐字相同** |
| C4 勘误 | ✅ F06 锁的两个余量(0.178563 / 0.16333)在新 CSV 上复算仍在带内(|Δ| ≤ 4.3e-05 < 1e-4),**锁无需更新**;Claude 审计 §4.4 已标注;progress.md 产品指标表已补 |
| C5 F06 + 机制测试 | ✅ **64 passed + 1 skipped**(F06/A27/A30/A06;80 步预算的机制测试不受影响,正是预注册的边界) |
| C6 收口全套 | 见基线行 |

### 4.2 ★★★ 预算假设被否(G-DOMAIN 的相反结果,如实记录)

| 腿 | 旧(80 预算) | 新(400 预算) | 判读 |
|---|---|---|---|
| medium seed-0 C=3.0 | `limit_cycle` @80,cl 0.519284 | **`limit_cycle` @400**,cl 0.519279,\|R\| 4.7e-07 | **不是预算** —— 400 步全烧完仍是极限环,答案不动 |
| medium seed-0 C=1.0 | `clamping` @47,1135/0 | `clamping` @47 | 钳制失败,与预算无关 |
| medium seed-5 C=1.0 | `clamping` @70,cl 0.63779 | `clamping` @78,cl **−0.02045** | 钳制失败;**失败态本身跨档换了吸引子**(记录,不投票) |
| medium C=1.5(两种子) | 收敛 @62,cl 0.34134 | 收敛 @**203**,cl **0.34134** | 答案不变;步数 62→203 = BLAS=1 走的路径不同(H30 附录 C.2 的 b1@400 = 203 步,**逐位咬合**) |
| coarse 全部 6 腿 | — | **逐位相同**(cl 6 位、步数同) | 粗网格对 BLAS 稳健 |

⇒ **M1c seed-0 的覆盖失败没有一个百分点是预算缺陷**:C=1.0 是钳制、C=3.0 是真极限环。
Claude 审计 §4.4 的线索「4 条 medium 腿停在 80 步、0 钳制」被实测闭合 —— 那条腿的
「80 步封顶」是极限环的**症状**而不是预算的**症状**(与 H30 附录 C.3 的镜像教训同族:
`budget_limited` 与 limit cycle 互相误判过两次,这是第三次,而这次事先写了判别式)。
**M1c seed-0 的覆盖失败现在需要一个非预算的解释**(极限环本身 = 未冻结选择集抖动;
钳制 = m_cap 物理),登记进裁决项。
★ 3seed 记录的变化更大(它产生于 2026-08-16 的另一代码/线程状态):其 seed-0 C=1.5
medium(clamping 失败)与 seed-5 C=1.5 medium(limit_cycle @80)在新配置上**都收敛**
(@203,cl 0.34134)—— 即旧 3seed 记录的失败两腿是**当时环境的吸引子**,不是缺陷。

### 4.3 产品指标(登记,非目标)

判定**逐字不变**:M1b PASS(两种子,−16.33 / −16.15 %);M1c PASS(seed 5,17.86 %)/
FAIL-覆盖(seed 0)。(a) 的 RECORDED 读数:x_shock coarse 0.6195/0.6073、medium
0.6006,均在 committed 0.62 ± 0.03 带内(与 GS1b.11 的 published 值逐位一致)。
**移动的只是归因,不是判定** —— 已在 progress.md 产品指标表登记。

## 5. 判定

**PASS** —— C1–C5 全过(C3 的结果是「判定不变 + 覆盖失败归因闭合」,这正是
预注册 G-DOMAIN 写下的相反结果分支);kill criterion 未触发(无收敛腿转为不收敛;
新旧 CSV 腿数/种子结构一致:12/12、18/18)。

一句话:**M1 门对齐 400 步重跑,判定逐字不变;M1c seed-0 的覆盖失败实测不是预算
缺陷(极限环在 400 步仍是极限环),审计线索被否,F06 锁因余量全部在带内而未动。**

## 6. 下一步

1. **§34 轮 §4.5 的 6 项历史文档失配立案一轮**(L3→L4 参考列更正会连锁
   CL_BAND_M0730 校准叙事,需自己的预注册)。
2. **W2.1(我方误差棒)** —— 仪器与产品指标基线都已干净。
3. 裁决项(追加):M1c seed-0 的覆盖失败的**非预算**处置 —— C=3.0 的真极限环
   (未冻结选择集抖动,frozen-selection 路线已有登记)与 C=1.0 的钳制(m_cap 物理),
   是否立项;`layer_attribution_D.csv` 删/留(仍开着)。

## 基线行

**803 passed / 32 skipped / 2 xfailed / 0 failed @1591.29 s @16 线程纪律(实际 BLAS=1)**
—— 与第 34 轮逐项相同;本轮的 `pyfp3d/__init__.py` 重构(公共 API 提取)未移动任何答案。
