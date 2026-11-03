# devflow_adoption_audit —— openspec-devflow 采用评估的支撑产物（phase 7 第 10 轮）

轮次文件：[`20261104-0524-devflow-adoption-audit.md`](../../../docs/dev_phase_seven/20261104-0524-devflow-adoption-audit.md)

**性质：一次性审计 / 调查；不改生产库、测试、默认值、committed gate 证据或参考数据。**
无测试加载本目录；phase 7 期间留原位，收口归档。

## 文件

| 文件 | 是什么 | 怎么复现 |
|---|---|---|
| `exp_extra_inputs.py` / `_output.txt` | devflow v0.2.0（`cc0a9f1`）对**被 Git 忽略**的 Case `extraInputs` 与 checkout 内 raw 的指纹/脏检查行为，三臂 + 对照臂 | 见脚本 docstring（需 `openspec-devflow-dev` 环境与 devflow checkout） |
| `exp_junit_xfail.py` / `_output.txt` | pytest 的 JUnit 报告怎样编码 skip 与 strict xfail（devflow 判定读的计数） | `conda run -n up3d python -I exp_junit_xfail.py` |
| `targeted_12red_20261104_console.txt` | 迁移机上 12 红所在 11 个门文件的定向重跑原始 console 日志（`080f291`，16 线程，暖缓存，1051 s） | 同轮次文件 §3.3 的命令；**全套那一跑的日志未保留** |

## 读数（2026-11-04）

```text
arm1 extra-input hash changed: True | dirty after ignored edit: []
arm1 ignored in-checkout raw: hash changed: True | dirty: []
arm2 (tracked control) dirty after edit: ['cases/meshes/m6.msh']
arm3 absent ignored extra-input: {'cases/meshes/m6.msh': None} | dirty: []

console: 1 passed, 1 skipped, 1 xfailed
testsuite counts: errors="0" failures="0" skipped="2" tests="3"
skipped types: ['pytest.skip', 'pytest.xfail']
```

⇒ 被忽略文件经 `extraInputs`/raw 进入指纹但**不报 dirty**，缺文件时得 `None` 而**不报错**；
pytest 把 xfail 计入 `skipped`，所以 devflow 的“全部 skip ⇒ `skip`”规则会把只含 xfail 的 Check 判 `skip`。

★ `targeted_12red_20261104_console.txt` 的断言消息只有 `.9e` 精度且每个测试只报第一个越界列 ——
它**不是**位级证据。全精度比较由 baseline-provenance 四臂运行承担（见轮次文件 §12）。
