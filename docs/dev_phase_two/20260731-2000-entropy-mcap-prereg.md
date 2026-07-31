# 预登记:熵修正不得从**被限幅的单元**里读激波前马赫(3-D m_cap 逃逸)

轮次文件名:`20260731-2000-entropy-mcap-prereg.md`
分支:`claude/phase-two-s2`
**本文件在改任何库代码之前提交**(roadmap 原则 8)。
上游证据:`20260731-1800-m3-budget.md` §2 + `bench/gate_results/m3_budget.csv`

---

## 1. 诊断(已测,不是猜)

M6 medium、M0.8395、P14 配方、`entropy_correction=True`:
`σ_min = 0.0`、`m1_max = 2.9999999999999996`、57 个 floored 单元、|R| 2.49e-06 未收敛
—— 与 G8.2 那条 strict xfail 的病征逐位吻合。同样的代码在 M6 **coarse** 上完全健康
(σ_min 0.9713、m1_max 1.3376、0 floored、严格收敛)。

**机制(读代码 + 那个数字定位):**

* `NewtonWorkspace.eval_residual` 里 `q2l = limit_q2_field(q2n, m_inf, m_cap, γ)`,
  `m_cap` 默认 **3.0**;`lim = (q2l == q2n)`(True = **未**被限幅)。
* `refresh_sigma` 把 **`state["q2l"]`**(**限幅后**的场)交给 `EntropyOperator.sigma`。
* 于是在被限幅的单元上,膝点行走取到的"激波前马赫"**就是限幅值 3.0**
  ⇒ `s = σ_RH(3.0) ≈ 0.328`,而不是一个物理值;
* 再经 `transport_sigma` 沿 donor 链**连乘** ⇒ 0.328^k → **σ → 0** ⇒ 密度塌陷
  ⇒ 更多单元被限幅 ⇒ **正反馈**。这就是 coarse/medium 的分水岭:coarse 没有被限幅的单元。

⚠ 这是**默认路径上的缺陷**(熵修正 2026-07-31 已翻为默认 ON,是我翻的),
优先级高于 S2 的新工作。

## 2. 修法(事前写定,并说明为什么这是**正确性**而不是调参)

`limit_q2_field` 一旦限幅,就等于系统已经宣布"**这个单元的速度不是物理值**"。
从一个非物理速度算熵产,按定义无意义。所以:

**(F1)** `shock_factor_sweep` 接受 `lim` 掩码;**行走在被限幅的单元处停止**;
且若后激波单元本身或其 donor 被限幅,则 **`s = 1.0`(不施加修正)**。
**(F2)** `EntropyOperator.sigma(..., lim=None)`;`lim=None` ⇒ 全 True
(**向后兼容**,且是下面 A/B 的对照腿)。
**(F3)** `refresh_sigma` 传入 `state["lim"]`。

**不做**的事(记录,免得日后当成遗漏):不给 σ 加下限(floor)。
若 (F1) 之后 σ→0 仍出现,那是**另一个**机制,应单独测量,而不是用一个 floor 盖住。

## 3. 事前登记的判据

| 编号 | 内容 | 判据 |
|---|---|---|
| **E1** | 位相同性:**无被限幅单元**的算例,修后必须与修前**逐位相同** | 2.5-D M1a coarse/medium/fine + G4.1 coarse 的 cl_p 全部 **bitwise 相等**;任一不等 ⇒ 先查原因,不得接受 |
| **E2** | M6 **medium** ON 腿:R2 有效性闸门**转为 PASS** | σ_min **> 0**、m1_max **< m_cap**、严格收敛(\|R\| ≤ 1e-8) |
| **E3** | M6 medium ON 腿的 cl_p 回到 P14 锚点邻域 | \|cl_p − 0.277628\| / 0.277628 **< 5%**(当前 ON 腿是 0.2110,低 24%) |
| **E4** | M6 **coarse** ON 腿:**逐位不变**(它本无被限幅单元,是 (F1) 的空对照) | cl_p 与 σ_min bitwise 相等 |
| **E5** | 修后 M6 medium 的 R3 读数(激波带能否被熵修正改善)= **RECORDED** | 无 pass/fail;这是 M3 第二项的第一次可读测量 |

**失败即负结果**:若 E2 不过(σ 仍塌陷),则 (F1) 被**证伪**,机制另有其物 ——
按 §2 的"不做"条款,**不得**用 floor 盖过去,而应把 E2 的失败作为本轮结论记录。

## 4. 影响范围与成本

改动仅在 `pyfp3d/kernels/entropy.py`(+`lim` 参数)与 `pyfp3d/solve/newton.py`
(传 `state["lim"]`);**Picard 路径**(`picard.py`)也读 σ ⇒ **backport check 必做**(纪律 #9),
答案写进本轮结论。成本:E1 约 2 min,E2/E3 两腿 medium 约 9 min,E4 约 20 s。
