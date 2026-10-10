# 20261010-1230 — **判定:H30 层 2 返工 —— BLAS 运行时固定生效,F09 升级为读库判据;同批修三个仪器缺陷**

轮次文件名:`20261010-1230-h30-layer2-runtime-verdict.md`
阶段 / 门:无门(工程 / 仪器修复)—— 预注册:
[20261010-1100-h30-layer2-runtime-prereg.md](20261010-1100-h30-layer2-runtime-prereg.md)
提交:分支 `claude/phase-six-audit`(本文件随判定提交)

---

## 1. 目的

让「BLAS 单线程固定」在两种 import 顺序下都真实生效,并让 F09 能看见它失效。
判据 C1–C6 与 kill criterion 见预注册,**事前写死,未事后修改**。

## 2. 改了什么

- `pyfp3d/__init__.py` — env 强制覆盖保留;**追加运行时固定**:numpy 已在
  `sys.modules` 时,从 `/proc/self/maps` 找已加载的 `libopenblas*`,
  `openblas_set_num_threads(1)` 并回读验证;找不到/压不住就 `warnings.warn`,
  不静默。逃生舱 `PYFP3D_ALLOW_BLAS_THREADS=1` 语义不变(整个固定跳过)。
  仍不 import numpy。
- `tests/F/test_F09_blas_determinism.py` — **重写**:探针跑**两种 import 顺序**
  (pyfp3d-first 与 numpy-first = pytest 的真实顺序)× 外部 BLAS {8,16} = 4 腿;
  **前提断言改读库的实际线程数**(ctypes + /proc/self/maps,不读 os.environ);
  判据不变(同顺序跨档逐位相同)+ 新增「两顺序互相逐位相同」;
  顺手修掉自相矛盾的机制注释(:38 强制覆盖 vs :108 setdefault,实际为强制覆盖)。
- `bench/run_capability_locks.py` — 快层子进程显式走逃生舱(BLAS=8 标定保留),
  **逃生舱只设在子进程 env 上**(教训见 §4.4)。
- 同批仪器修复(独立提交):`bench/run_m3_budget.py:261` 死键
  (`sigma_converged`/`entropy_converged` 都不存在 → 真键 `sigma_transport_converged`,
  R2 守门员 transport 支路原来按构造恒真);`tests/F/test_F05:108` 裸
  `import usability` → `from bench import usability`(原来靠 D11/D12 收集期的
  sys.path 污染才能解析,**单跑即红**,实测复现);D01:101 / D02:105 /
  D02:174+213(CSV+PNG)/ D03:344 共 **5 处无守卫写盘**补 `gate_figures_enabled()`。

## 3. 怎么测的

```bash
conda activate up3d; export PYTHONNOUSERSITE=1
# C1 ctypes 探针(4 种组合):numpy-先/后 × 外部设 16/不设
# C2/C3 F09 正常路径:        python -m pytest tests/F/test_F09_blas_determinism.py
# C4 G-TEETH:               PYFP3D_ALLOW_BLAS_THREADS=1 python -m pytest tests/F/test_F09...
# 泄漏复现:                  python -m pytest tests/F/test_F01_provenance.py tests/F/test_F09...
# C5/C6 全套(16 线程纪律):  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
#                            python -m pytest tests/ -q -rfEsx -p no:cacheprovider
# 写盘守卫 G-TEETH:默认跑 D01+D02 后 git status 必须干净
```

## 4. 结果

### 4.1 判据逐条

| 判据 | 结果 |
|---|---|
| C1 numpy 先导入 + 外部 16/不设 ⇒ 库实际线程 | **==1,==1**(修复前 16/24);pyfp3d-first 不退化(1);逃生舱放行(16) |
| C2 F09 四腿跨档逐位相同 + 两顺序互相同 | **4 passed @161.64 s** |
| C3 前提断言读库:四腿 `openblas_get_num_threads()==1` | 绿(含在 4 passed 内) |
| C4 G-TEETH(逃生舱) | **红得正是地方**:`test_the_library_is_actually_pinned` 红(库跑 8/16)+ `test_the_answer_is_bit_identical` 红(cl 0.34011432026323 vs 0.34133984793989103,与文档记载的敏感值同族) |
| C5 全套无新增红 | 见 §4.3/§4.4(第一次全套抓到一个我自己的泄漏缺陷,修后复跑) |
| C6 墙钟 | 全套 2172.95 s vs 修复前 2062.73 s(+5.3 %,与 H30 记的 +4.5 % 同量级;两次负载不同,只记录) |
| 写盘守卫 G-TEETH | 默认跑 D01+D02 后 `git status` 干净(6 passed + 5 skipped);F05 单跑 **9 passed**(修复前 1 failed) |

### 4.2 ★ 全套下证据锁的失配切换到 BLAS=1 签名(固定生效的间接实测)

| 门 · 列 | 修复前(=Claude 的"外部 BLAS 16"列) | 修复后全套 | Claude 的"外部 BLAS 1"列 |
|---|---|---|---|
| D05 cp_rms_upper | 3.34e-05 | **3.11e-05** | 3.11e-05 |
| D06 cp_rms_upper | 7.02e-05 | **3.33e-03** | 3.33e-03 |
| D11 cl_p | 2.05e-06 | **1.75e-06** | 1.75e-06 |
| D12 mach_max | 3.37e-06 | **1.55e-04** | 1.55e-04 |

⇒ 套件进程现在**真的跑在 BLAS=1**(与另一份审计独立测的 BLAS=1 读数逐位一致 —
— 两份审计的数字在这一步互相咬合)。这 6 条锁仍红是**预期**:证据锁在默认 24 档,
重刷属下一轮(带纪律 11 勘误清单)。
★ D08/D13 的"IBL 不可复现"说法被 24 档实验**否掉**(两条锁在 24 档全绿,
1625 s)—— 6 条红是同一个 BLAS 根因;IBL 松耦合的已知不可复现性仍是
`76b70c0` 登记的开放议题,但它不是本轮红色的原因。

### 4.3 全套(第一次,修复泄漏前):**8 failed / 795 passed / 32 skipped / 2 xfailed @2172.95 s**

预期的 6 条证据锁 + **F09 两条新红**。

### 4.4 ★★ 第一次全套抓到了我自己 —— 泄漏链与修法

我把快层的逃生舱写成了 `run_capability_locks.py` 的**模块级**
`os.environ.setdefault("PYFP3D_ALLOW_BLAS_THREADS", "1")`。全套里 F01 在收集期
import 这个 bench 模块 ⇒ 逃生舱泄漏进整个 pytest 进程的 `os.environ` ⇒
F09 的探针子进程继承它 ⇒ 固定被放行 ⇒ F09 红(§4.3)。
证据锁**没**被污染:pyfp3d 在 A 类收集期已完成固定,**库状态不随后续 env 变动
回退**——所以 §4.2 的 BLAS=1 签名有效。
**修法**:逃生舱只设在快层**子进程的 env**(`main()` 里显式构造),模块 import
零副作用。复验:import 后 `os.environ` 无该变量;**F01+F09 同进程 27 passed**
(正是原泄漏场景)。
★ 这与 F05 的 sys.path 污染同族:**bench 模块的模块级副作用会被收集期 import
带进套件进程**。登记:值得一条 F 类不变量(bench 模块 import 不得改
`os.environ` / `sys.path`),本轮只修实例,不立新门。

### 4.5 全套(收口,修复泄漏后)

**797 passed / 6 failed / 32 skipped / 2 xfailed @1590.15 s @16 线程纪律(load ~4)**。
账目闭合:837 collected = 797 + 6 + 32 + 2;与修复前基线 835 的差 = F09 由 2 条
扩到 4 条(+2 collected,+2 passed)。**6 条 failed 正好是、也只是那 6 条证据锁**
(失配已在 §4.2 切到 BLAS=1 签名);F09 四条在全套内转绿 ⇒ 泄漏修复生效。
★ 本次全套后 `git status` 无任何 `cases/` 改动 —— 写盘守卫在全套尺度上兑现
(此前每跑必改写 D01/D02 的 3 个已跟踪 CSV)。

## 5. 判定

**PASS** —— C1–C4 全过;C5 收口全套无预期外红(6 条证据锁 = 已归因、待重刷);
C6 记录。kill criterion 未触发。

一句话:**BLAS 固定从「写给空气的 env」变成「钉进库的运行时调用」,F09 换成
读库判据并覆盖 pytest 的真实 import 顺序;全套把 6 条证据锁的失配切到 BLAS=1
签名(固定生效的间接实测),下一轮在真 BLAS=1 下重刷证据。**

## 6. 下一步

1. **在真 BLAS=1 下重刷 D05/D06/D11/D12 证据**(`PYFP3D_GATE_FIGURES=1`),
   走纪律 11 勘误清单;证据 CSV 加**库实际 BLAS 线程数**列(从库读)。
   先查 D06 为什么 BLAS 敏感度大两个数量级(3.33e-03 vs D05 3.11e-05)。
   D08/D13 同批重刷(根因相同,已实测)。
2. M1 门 `n_newton_max` 80 → 400 重跑(M1c seed-0 覆盖失败的归因);
   顺带更正 D05:73 / D11:62 已过期的「逐字一致」注释。
3. ~~快层跑一次~~ **本轮已跑:5/5 绿 @550 s(9.2 min),load 3.3**;逃生舱路径
   端到端生效(子进程 BLAS=8 标定保留),`capability_locks.csv` 已随本轮重刷
   (仅墙钟列动;这也是 H30 三层之后快层的第一次复跑,顺带收掉 W0.5 的尾巴)。
4. 登记未修:bench 模块 import 副作用的 F 类不变量(§4.4);`n_threads` 列的
   语义(记的是 numba,不是 BLAS)随重刷一并改。
5. 之后回 W2.1(我方误差棒)。

**基线行**(收口全套):**797 passed / 6 failed / 32 skipped / 2 xfailed,
0 error @1590.15 s @NUMBA=16/OMP=16/OPENBLAS=16(实际 BLAS=1,固定生效)**;
6 条 failed 全部是证据锁(D05 D06 D08 D11 D12 D13),失配 = BLAS=1 签名,
属已归因而非回归;重刷证据见下一步 1。
