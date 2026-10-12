# amg_shock_fix —— 「amg 能否修」（phase 7 第 4 轮，S3）

预注册：[`docs/dev_phase_seven/20261012-2230-amg-shock-fix-prereg.md`](../../../docs/dev_phase_seven/20261012-2230-amg-shock-fix-prereg.md)。
机制前提（p6 第 47 轮已闭合，`phases/p6/bench/studies/c05_amg_investigation/`）：
amg 的层级建在不含 Term-3 激波耦合的 Picard 块 A_ff 上，对 J_ff 的激波中心近零模态
（σ_min = 2.746e-6，右奇异矢峰在激波 x=11.4）无反演 ⇒ GMRES 600 内迭代真残差
纹丝不动 ⇒ Newton 停滞（C05 生产腿 40 步停 |R| 3.9e-6）。
**性质：研究/测量，不改库、不改判据、不改默认值；本目录随 phase 7 归档。**

## 答案（一句话）

**能修，走方向 2（AMG(A_ff) 打底 + 激波方向低秩校正）；方向 1（AMG 直接建在
完整 J_ff 上）实测死。**

## 实测（2026-10-12，BLAS=1 库钉，pyamg 钉种，eigsh v0 钉 rng(0)）

冻结点 = 生产 amg 腿末态（|R|∞ = 3.914e-6，n_free = 5174，GMRES 预算 =
生产单步预算 600 内迭代）：

| 预条件子 | rel 真残差 @600 内迭代 | 读法 |
|---|---|---|
| amg(A_ff)（生产形状） | 0.968 | 平 —— 机制现场 |
| **amg(J_ff)（方向 1）** | **0.732**（层级建成：4 层、opcx 1.066、有限） | **实测负**：层级建得动、建得便宜（9.2 ms vs 6.9 ms），但非对称 J_ff 上它不做预条件子 |
| **amg(A_ff) + 软模态校正（方向 2）** | k=1 1.11 → k=2 0.52 → k=4 0.135 → k=6 0.0139 → **k=12 6.9e-9** → k=20 3.3e-10（360 内迭代干净收敛） | **可修，单调剂量-响应** |
| ilu(J_ff)（对照） | 2 循环收敛 | 参照 |

端到端（本脚本内的**简化 Newton 环**：全步长、无 EW/线搜索/冻结 —— 与生产驱动器的
差异就在这三个省略，对照臂 direct 在同一环上）：

| 臂 | Newton 步 | 末 \|R\| | 墙钟 |
|---|---|---|---|
| ctrl_direct | 7 | 2.7e-14 | 0.13 s |
| amg_on_Jff | 40（帽） | 2.6e-07 | 18.5 s |
| **defl_k4**（每步重估 4 个软模态） | **8** | **2.2e-12** | 3.4 s |

## 代价（C4，口径如实）

- 层级建造：amg(J_ff) 9.2 ms vs amg(A_ff) 6.9 ms（×1.25；J_ff nnz 63254 vs 54178）。
- 方向 2 的成本驱动 = **每步软模态估计**：本研究用 LU 预条件 eigsh（k=20 一次
  45 ms；每步 k=4 重估）。★ 在 C05 这个尺度（5174 自由度）splu 本就瞬时，
  **defl_k4 比 direct 慢 25× 是预期形状，不是反对意见**；
  **生产尺度（M6 fine 45 万）的代价本轮测不了、不外推**（预注册 C4 写明）。
- ★ 采用前必须答的下一问（登记给将来）：软模态估计器要换成不带 LU 的
  （如 AMG 预条件的子空间迭代）或跨步复用 —— 那是一轮新研究，不在本轮时间盒。

## ★ 仪器勘误（我自己的，执行中抓出）

1. `svds(which="SM")` 在 J_ff 上不收敛（ARPACK 51741 迭代，0/6）—— p6 probe1
   踩过同一个坑（probe2 改 LU 预条件 eigsh），我第一版照样踩。⇒ 本研究软模态
   一律走 `soft_modes()`（C=(JJᵀ)⁻¹ 的 eigsh，v0 钉死）。
2. 第一版 k 扫描把 eigsh(which="LA") 返回的**降序** σ 当升序用，「前 k 个」
   实际是软集里**最大**的 k 个；且换 k 总集时集合的头会动 ⇒ 两次扫描不可比
   （k=6 总集时 k=4 读数 0.0345 vs k=20 总集时 0.588）。已修为升序取最软 k；
   修正后 k=6 读数（0.0139）与修正前那一列**逐位一致**（两次都含 σ_min 时一致），
   交叉验证通过。

## 复算

```
conda activate up3d
PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
    python bench/studies/amg_shock_fix/run.py     # 全管线 ~40 s + CSV + PNG
```

## 文件

| 文件 | 内容 |
|---|---|
| `run.py` | 全管线（生产 amg 腿复现 → 冻结点提取 → 两方向 → 端到端三环） |
| `freeze_point_gmres.csv` | 对照三腿的逐循环真残差 |
| `deflation_sweep.csv` | 方向 2 的 k 扫描（最软 k，升序） |
| `newton_legs.csv` | 端到端三环的逐 Newton 步 \|R\|∞ |
| `cost_table.csv` | 层级建造 / 软模态估计 / 端到端墙钟 |
| `amg_shock_fix.png` | 左：冻结点对照；右：deflation 剂量-响应 |
