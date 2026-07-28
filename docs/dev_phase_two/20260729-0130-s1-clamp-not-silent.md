# GS1.4：钳制不得静默

轮次文件名：`20260729-0130-s1-clamp-not-silent.md`
阶段 / 门：`S1` / `GS1.4`
分支：`claude/phase-two-s0-s1`

---

## 1. 目的与失败判据（执行前写下）

钳制（速度限制器 `m_cap`、人工密度下限 `rho_floor`）是代数覆盖：**它绑定的地方，
那个单元离散的已经不是全速势方程**，所以"残差为零"不再说明任何流动性质。
GS1.1 已经量到后果：同一套边界数据下，钳制能驻留一个**机器零残差的伪解**
（|R| 9.2e-13，激波离唯一精确位置 40 个单元，72 个单元被钉在 `rho_floor` 上），
而物理支路对下限四个数量级完全不敏感。

本轮把这件事变成**契约**：任何驱动器都不得把钳制态报成 converged。

| # | 判据 | 判定 |
|---|---|---|
| 1 | 审计出**所有**会把钳制态报成收敛的驱动器，逐个修 | 漏掉一个 = FAIL |
| 2 | 契约要有测试锁住（不是锁数字，是锁行为），且**探针必须真的触发钳制**（否则测试是空的） | 空测试 = FAIL |
| 3 | 伪解回归进 `bench/`（roadmap 对 GS1.4 的原文要求） | 缺 = FAIL |
| 4 | **影响面必须量出来**：全套测试跑一遍，逐个说明因此失败的用例 | 只说"应该没影响" = FAIL |
| 5 | 默认关闭/无关路径逐位不变（`bitcheck` 10/10） | 不满足 = FAIL |

**范围**：不加新旋钮；不改任何数值配方；只改"什么算收敛"的判定与上报。

## 2. 改了什么

### 2.1 审计结果：只有一个缺口

| 驱动器 | 钳制时是否拒绝 converged | 处置 |
|---|---|---|
| `solve_newton_lifting`（`newton.py`） | ✓ 早已拒绝（收敛判定里带 `n_limited == 0 and n_floored == 0`） | 只加 `clamped` 汇总标志 |
| `solve_multivalued_newton`（`newton_ls.py`） | ✓ 同上 | 不动 |
| `solve_transonic_lifting`（`continuation.py`） | ✓ 通过 `physical` 标志 | 不动 |
| `solve_subsonic`（`picard.py`，非升力） | 不适用（该路径不用限制器/下限，构造上不可能钳制） | 不动 |
| **`solve_subsonic_lifting`（`picard.py`）** | **✗ 会报 converged=True** | **修** |

被修的这一个正是**其他所有路径的 Picard 种子**——一个钳制的种子会把伪支路
静默地交给调用者。

### 2.2 具体改动（库）

* `pyfp3d/solve/picard.py::solve_subsonic_lifting`：收敛判定加
  `and not clamped`（`clamped = n_limited > 0 or n_floored > 0`），并在返回字典里
  加 `clamped` 标志；
* `pyfp3d/solve/newton.py::solve_newton_lifting`：只加 `clamped` 汇总标志
  （行为本来就已正确，但两个计数器太容易被忽略——本轮的 GS1.7 扫描就又一次
  从 `M_max = 3.0` 的钳制态里读了数）。

### 2.3 测试与回归（新增）

* `tests/test_s1_clamp_not_silent.py`（4 个常开测试，~20 s）：锁**契约**不锁数字 ——
  两个驱动器都上报 `clamped` 且与计数器一致；用 `m_cap = 0.55`（低于 M0.5 工况自身
  的吸力峰 M ≈ 0.68）**强制触发**钳制，驱动器必须永不报收敛；干净解与钳制解必须
  可区分（Γ 差 > 1e-4、φ 不相等）。
  ★ 第一版探针用 `m_cap = 1.02` **没有钳制到任何单元**（q²(M=1.02) = 3.62 远高于
  该工况的 q²_max ≈ 1.78），测试是空的 —— 判据 2 的"探针必须真触发"就是为这个写的，
  它当场生效了。
* `bench/s1_duct/regress_floor_spurious.py`（判据 3）：把 GS1.1 的发现变成带断言的
  常驻回归。实测：**物理支路对 `rho_floor` 四个数量级散布 0.00e+00**（激波
  11.93155 一位不变，`n_floored` 恒为 0）；**伪支路在出厂下限 0.05 下仍然复现**
  （|R| 9.2e-13、x_shock 7.98、72 个 floored）。退出码 0。

## 3. 怎么测的

```bash
PYTHONPATH=. python -m pytest tests/test_s1_clamp_not_silent.py -q     # ~20 s
PYTHONPATH=. python bench/s1_duct/regress_floor_spurious.py            # ~1 min
python bench/bitcheck.py --diff bench/results/bit_a.npz bench/results/bit_gs14.npz
PYTHONPATH=. python -m pytest tests/ -q                                # 影响面
```

## 4. 结果

| 判据 | 结果 |
|---|---|
| 1 审计全覆盖 | **PASS** —— 5 个驱动器逐个查明，缺口只有 1 个并已修（§2.1） |
| 2 契约测试 + 探针有效 | **PASS** —— 4/4 通过；空探针问题被判据 2 当场逮住并修正 |
| 3 伪解回归 | **PASS** —— `regress_floor_spurious.py` 退出码 0；物理支路散布 **0.00e+00**，伪支路复现 |
| 4 影响面 | 见下 |
| 5 逐位不变 | **PASS** —— `bitcheck` 10/10（改动只影响"是否标记为收敛"，不改任何算术） |

**全套测试**：`1 failed / 657 passed / 27 skipped / 2 xfailed`，939 s @8 线程。
唯一的失败是 `test_v6_wake_sheet::test_ab_bit_identity_gate_free_library`，
**归因链如下（每一步都是测量，不是推断）**：

1. 我的第一个猜测是"GS1.4 改了 `solve_subsonic_lifting` 的 break 条件，导致 Picard
   多跑几步" —— **撤回 picard 改动后同样失败**，猜测被否证；
2. 换到 `main` 的 worktree 上跑同一个测试 —— **`main` 上也失败**，与本分支无关，
   是既有失败；
3. 决定性检查：把**同一个 commit** 跑两遍（两个全新 worktree、同机同线程数）——
   `max relative 1.024，6104/6106 个节点不同`，**与该测试报告的差异完全同量级**。

⇒ **这个测试在比较一个不可复现的计算**：松耦合粘性回路在 coarse 上不是 run-to-run
确定性的。（对照：`bitcheck` 的无粘探针 10/10 逐位一致 —— 非确定性局限在
`pyfp3d/viscous/`，与第一阶段纪律 #12 的记载一致。）

处置：**该测试退役**（`@pytest.mark.skip`，reason 里写明完整归因链与证据脚本），
两遍测量固化为 `bench/s1_duct/check_loose_loop_determinism.py`。
GV6.1 (a)(ii) 的**证据本身不受影响**（committed 在 `cases/analysis/v6_1_wake_sheet/`），
退役的是一个**声称的东西它做不到**的活守卫 —— 这正是原则 1 的要求。

退役后：`tests/test_v6_wake_sheet.py` 7 passed / 1 skipped。

## 5. 判定

**GS1.4 = PASS（5/5 判据）。**

三条可引用结论：

1. 钳制契约的缺口**只有一个**（`solve_subsonic_lifting`，即所有路径的 Picard 种子），
   已修；两个主驱动器新增 `clamped` 汇总标志。
2. 伪解回归常驻化：物理支路对 `rho_floor` 四个数量级**散布 0.00e+00**，
   伪支路在出厂下限下仍然复现（|R| 9.2e-13、72 个 floored）。
3. ★ 意外收获（判据 4 逼出来的）：`test_ab_bit_identity_gate_free_library` 的既有失败
   与本轮无关，**根因是松耦合粘性回路不可复现**（同 commit 两遍差 max rel 1.024）。
   该测试退役并留下证据脚本。

★ 我在归因上先猜错了一次（以为是自己的改动），是**撤回改动重跑**这个动作纠正了它 ——
判据 4 要求"逐个说明失败原因"而不是"应该没影响"，这条纪律当场生效了两次
（另一次是判据 2 逮住了空探针）。

## 6. 下一步

**GS1.5 = S1 收口**：按 GS1.3b/GS1.7 的结论，M1 的工况处在 fold 邻域，
不是靠加密、换闭合或熵修正能达到的。收口轮要做的是：
① 如实跑一次 M1 门并记录 FAIL；② 把 S1 的整体结论（能力边界 + 已排除的路线 +
仍开放的路线）写成一份可交给使用者裁决的总结；③ 更新 roadmap 的 S1 状态。
