# pyFP3D 开发现状（快照）

> **快照日期：2026-07-11。** 本文件是给人读的高层概览，**不是**权威进度源：
> 阶段/gate 状态以 [roadmap.md](roadmap.md)（含进度台账）为准，当前阶段细节以
> [agent-rules.md](agent-rules.md) 的 "Current phase" 行为准，证据在
> [demo_report.md](demo_report.md)。若本文件与它们冲突，以它们为准。
> 每次 gate 关闭后按 CLAUDE.md 工作流更新 roadmap/agent-rules 时，顺手刷新本文件。

## 一句话状态

求解器主线（Track P）P0–P8 全部关闭，P10 的两个效率/可行性 gate（G10.2、G10.3）
已关；**当前正在跑 P9 判别阶段**（网格收敛 vs 精度地板归因，纯证据阶段、不改数值），
其结论（G9.3）决定是否开启 P11 曲面壁元这一大工程。

## 阶段一览

| 阶段 | 状态 | 要点 |
|------|------|------|
| P0–P3 | ✓ | 基础设施、Laplace、尾迹割缝/Kutta、亚声速压缩（P1 仅 G1.6 未关，见下） |
| P4 | ✓（带勘误） | 人工密度跨声速 Picard；★2026-07-11 勘误：Picard 态**不是**离散方程的解（Newton 残差 2.2e-4），P4 gate 降级为 Picard 质量/鲁棒性 gate（warm-start 引擎） |
| P5 | ✓（带注记） | ONERA M6 M0.84/α3.06；同类欠收敛注记（真解 cl 高 +5.8%/+7.9%，激波位置不变） |
| P6 | ✓ | 表面 Cp 锯齿 = 壁面梯度**恢复**伪迹（非通量），法向门控平滑修复 |
| P7 | ✓ | 冻结选择下的 ∂ρ̃/∂φ 精确灵敏度，FD 验证 3e-10~6e-9 |
| P8 | ✓ | 全耦合 (φ_red, Γ) Newton：精确 Jacobian + 直接步 + 停滞自适应冻结；G8.2 M6 medium 249 s（G10.2 后 ~145 s）、G8.3 全套件 302 s |
| P9 | **运行中** | 三点网格研究 + Richardson 外推，预注册判别带（cl_KJ∞ ≥0.283 分辨率主导 / ≤0.278 地板确认）；M6 coarse/medium 已缓存（干净收敛 0/0），fine（~450k dof）解算进行中 |
| P10 | ◐ | G10.2 ✓（分裂裁决：M6 +41% 提速已晋升；折叠区禁用松容差）；G10.3 ✓（裁决：**保留 Mach ramp**）；G10.1（非升力 Newton 入口）未动、无顺序约束 |
| P11 | 未开 | 曲面/等参壁元（承接 G1.6 + V6<1%），**是否开启取决于 G9.3 裁决（用户仲裁）** |
| P12 | 未开 | Backlog：离散伴随、VII transpiration、混合单元/BO 标定 |
| M0/M1 | ✓ | 准 2D 挤出网格族 + M6 后掠尾迹网格族（.msh gitignored，脚本再生） |
| M2 / Track B / Track V | 未开 | 翼身组合体宜与 Track B（水平集嵌入尾迹）同排；Track V = VII 粘性耦合，已设计未动工 |

## 长期挂起项（勿反复重提）

- **G1.6 球面 Cp <2%**（strict xfail，11.6%）：已定根因 = 平坦面片壁上的自然边界条件
  （变分罪）；h 加密、恢复调参、Nitsche、边界数据修正**均已用证据排除**。
  唯一在案路线 = Option C 重定义 + P11 曲面壁元。
- **V6 <1%**（CL_p 对 CL_KJ 的 O(h) 离散地板）：挂到 P11，P9 正在检验其归因。
- 后掠 TE Kutta 探针跨站共享（P5 记录的鲁棒性隐患，未修）。

## 当前最大困难（按层级）

1. **模型极限（不是 bug）：FP 非唯一性折叠。** NACA medium 在 M≈0.79 附近进入
   全势方程的折叠区（dcl/dM ≈ 6–10）：同工况不同网格不可比、warm-start 有陷阱
   （G10.2b 负结果）、M0.80 无可达孤立解、无 ramp 直接求解 class C（G10.3）。
   对策已固化为纪律：折叠区只做单网格回归锁、Mach ramp 保留、松容差禁用。
2. **精度地板：平坦 P1 壁面。** 0.019 的 M6 升力缺口（cl_KJ 0.2692 vs
   Tranair/KRATOS 0.288）归因于尖 TE/LE 的 P1 壁面梯度——但这一归因是推断级，
   P9 正是为把它变成证据级（或推翻）而设。修复（P11 曲面壁元）是公认的大工程。
3. **技术未知数（进行中）：M6 fine ~450k dof 的 Newton 解算** 的鲁棒性与性能
   （真 3D LU 填充是已量化的坑，lagged-LU 每层一次分解 + AMG 兜底，预算 ≤2 h）。

## 运行/成本纪律（易踩的坑）

- 16 线程封顶且**必须同时盖住 BLAS/OMP**（`NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
  OPENBLAS_NUM_THREADS=16`；漏了 BLAS 会慢 ~33% 并挂掉 G8.2 断言）。
- 贵重工件不随手重算（P4 heavy demo ~40 min、G4.1 medium ~17 min、P5 medium
  45–75 min、M6 fine 网格数分钟）；已提交的 CSV/PNG 是权威。
- fine 网格与解算 npz 一律 gitignored，本地缓存、demo 缺则重算。
- 回归基线：**184 passed + 8 skipped + 2 xfailed**（重 gate 走
  `PYFP3D_TRANSONIC_GATES=1`）。内核/装配改动后先跑 `tests/test_v0_freestream.py`。
