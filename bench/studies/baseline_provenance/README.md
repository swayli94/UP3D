# baseline_provenance —— 迁移机 12 红的四臂出处测量（phase 7 第 11 轮）

判定：[`20261104-0640-baseline-provenance-verdict.md`](../../../docs/dev_phase_seven/20261104-0640-baseline-provenance-verdict.md)
登记（四臂与四组配对）：[`20261104-0524-devflow-adoption-audit.md`](../../../docs/dev_phase_seven/20261104-0524-devflow-adoption-audit.md) §3.4 / §12

**性质：一次性测量；不改生产库、测试、默认值、committed gate 证据或参考数据。** 无测试加载本目录。

## 运行

- 代码：`080f291`，在 detached worktree `/home/lrz/code/up3d-baseline-prov` 中运行（避开主 checkout 的并发修改）；
  每臂 `head.txt` / `status_*.txt` / `pyfp3d_origin.txt` 证明 HEAD、干净工作树与 `pyfp3d` 导入来源。
- 11 个门文件（12 红所在），16 线程，`-p no:cacheprovider`；四臂顺序：
  `native_cold`（空 `NUMBA_CACHE_DIR`）→ `native_warm`（同目录第二跑）→ `generic_cold`（`NUMBA_CPU_NAME=generic`，另一个空目录）→ `generic_warm`。
- 2026-11-04 05:16–06:24（runner `run_arms.sh`；`runner_progress.txt`）。`run_arms.sh` 与 `prov_plugin.py` 里的绝对路径指向当时的
  `bench/results/baseline_provenance_20261104/`（gitignored 运行目录），保留原样作为出处。

## 仪器

- `prov_plugin.py`：pytest 插件，**只记录不改判定**（dry-check 时故意的红仍红）。包住 `tests._gate_evidence` 的两个比对函数，
  把 `fresh` 值以 `repr` 与同键的 committed 值一并写入 `fresh.jsonl`；D13 的 `runs` fixture 写入 `d13_runs.npz`；
  会话结束写 61 个 pyfp3d dispatcher 的 numba `cache_hits/misses`（`numba_stats.json`）。
  ★ D05/D06/D14/D15 的门在调用比对函数前已把值舍入到 10 位有效数字，所以这四个门的“逐位相同”= 10 位有效数字相同。
- `codegen_env.py`：每臂的有效 codegen（CPU 名、feature 串摘要、numba/llvmlite/LLVM、配置覆盖、BLAS）。
- `compare.py`：四组配对比较 → `compare_summary.md`、`compare_values.csv`（运行目录里自动生成后复制过来）。
- `classify.py`：**判定里引用的全部数字**由它从本目录算出 → `classification.csv`、`premise_b.csv`、`numba_cache.csv`。
  复现：`conda run -n up3d python -I classify.py`。★ 其分类标签是**事后描述**，未预注册。

## 产物

`arms/<arm>/`：`junit.xml`、`fresh.jsonl`、`d13_runs.npz`、`numba_stats.json`、`codegen.json`、`pytest_console.txt`
（原 `pytest.log`，`*.log` 被全局 gitignore 故改名）及出处小文件。numba 缓存目录（各 4.6 MB）不提交。
