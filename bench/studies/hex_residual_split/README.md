# hex_residual_split —— hex 族 NACA 气泡缺陷的**残余拆分探针**（phase 7 第 3 轮顺带项）

预注册：[`docs/dev_phase_seven/20261012-2135-d15-hex-naca-gate-prereg.md`](../../../docs/dev_phase_seven/20261012-2135-d15-hex-naca-gate-prereg.md) §1.2（时间盒 = 两臂）。
残余出处：p6 第 47 轮判定 §2.2 —— 缺陷在 hex 生成器还是在该拓扑上被触发的熵冻结，
当时的三件仪器分不开（hex 各腿 sigma 冻结搅动 |dσ| 0.19–0.21 vs tet 0.067）。
**性质：研究/测量，不改库、不改判据；本目录随 phase 7 归档。**

## 复算

```
conda activate up3d
NUMBA_NUM_THREADS=16 python bench/studies/hex_residual_split/run.py     # 两臂求解 + arms.csv（~1 min）
python bench/studies/hex_residual_split/run.py --analysis-only          # 只读 CSV 重打印
```

- 两臂都在 hex **medium** 上、M0.80/α1.25、D05 生产配方，只动一个旋钮：
  **E0** = `entropy_correction=False`；**F0** = `freeze_tol=None`（只关 N5 迎风选择冻结；
  熵修正照常跑，且实测它自己的 sigma 冻结仍介入，`sigma_frozen=1`）。
- 锚（生产臂）**不重算**：读已提交的 `cases/gates/D15_hex_naca_boundary/summary.csv`
  的 `hex|medium` 行（证据纪律 2）。仪器（气泡计/分段）从 D15 门文件导入，不抄第二份。
- BLAS 由库钉 1（脚本断言）。

## 实测（2026-10-12，BLAS=1，16 线程）

| 臂 | cl | 气泡计 | x_shock | 收敛 |
|---|---|---|---|---|
| 生产（committed D15 锚） | 0.258223 | **+0.199** @0.457 | 0.6094 | ✓ 31 步 |
| **E0**（熵修正 OFF） | 0.359584 | **−0.646** @0.543（消失） | 0.5590 | ✓ 32 步，0 钳制 |
| **F0**（N5 冻结 OFF） | 0.327155 | **−0.255** @0.543（消失） | 0.5255 | ✓ 14 步，0 钳制 |

## 读法（按预注册的 G-DOMAIN 表）

**两臂都把气泡拆没了** —— E0 命中「熵修正/冻结链路被牵连 ⇒ 偏求解器侧」，
F0 命中「冻结交互单独牵连」。合起来的结论：

1. **气泡不是网格几何的后果** —— 几何无辜此前已实测（壁面对精确剖面 max 差
   2.15e-5，W2.6）；现在两臂各只动一个求解器旋钮就把气泡拆掉。
   **残余从「分不开」偏向求解器侧的冻结/熵修正路径。**
2. ★ **两臂不是干净析因**：E0 摘掉的是整个熵修正（连 sigma 场一起），F0 只摘
   N5 选择冻结而 sigma 冻结仍介入 —— 「sigma 机制是否放大」没被分开，如实登记。
3. ★ 两臂的落点**互不相同也都不等于 tet**（cl：E0 0.3596 **越过** tet 0.3413 与
   CFL3D L4 0.3489；F0 0.3272 在中间；x_shock：0.5590 / 0.5255 / 0.5991 三处）——
   与立项的中心风险同族：**这条腿站在多吸引子邻域，生产配方只是其中一个分支**。
   深一层的问题（hex 拓扑为什么把瞬态带进这个分支，tet 为什么不）登记进
   phase 7 的 S5（xtr030 多吸引子）同一族，不在本轮升级仪器。

## 文件

| 文件 | 内容 |
|---|---|
| `run.py` | 两臂求解 + 读数打印（锚读 committed 门证据） |
| `arms.csv` | 两臂全状态：cl/x_shock/气泡计/收敛/钳制/sigma 冻结/出处列 |
