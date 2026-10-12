# 20261012-1400 — **判定:C10 hex 圆柱交叉门 —— cylinder_hex 补 wake + 两拓扑对照落地;裁决 c 立项登记**

轮次文件名:`20261012-1400-c10-hex-cylinder-verdict.md`
预注册:[20261012-1300-c10-hex-cylinder-prereg.md](20261012-1300-c10-hex-cylinder-prereg.md)
(使用者裁决 b「cylinder_hex 补 wake 组 + 非升力对照」与 c「tet medium r_far 20/30
不收敛立项研究」的执行轮)

---

## 1. 改了什么

- **`pyfp3d/meshgen/structured.py::cylinder_o_grid_2d`**:返回 4 元组 → 5 元组
  (`airfoil_hybrid_2d` 的既有约定),新增 `interior_edge_groups` 带 **wake** ——
  即 O-grid 的接缝(θ = 0,壁面点 (radius, 0) 沿 +x 排到远场)。圆柱没有尖尾缘,
  但环量通道必须有分支切口。**只加组:每条接缝边本来就是四边形竖边,零节点/
  零几何变动。**
- **`cases/meshes/cylinder_hex_2.5d/generate_cylinder_hex.py`**:解包并传
  `interior_edge_groups`;stats 加 `n_tris_wake` 列;inspection PNG 画 wake(绿)。
- **三级网格重生成**(coarse/medium/fine)+ 三个 `.msh.manifest.json`(既有
  manifest 机制,首随本轮)。
- **新门 `tests/C/test_C10_hex_cylinder_crosscheck.py`**:非升力腿(Γ=0,
  hex vs tet 各对精确 Cp)+ 升力腿(Γ=0.6,hex vs C06 内存 tet 族,各对各自
  精确 cl —— 口径:hex 弦 2 ⇒ 0.6,C06 族弦 1 ⇒ 1.2)。证据锁带
  `blas_threads` 出处列与**行数守卫**(见 §5 意外 a)。
- **证据**:`cases/gates/C10_hex_cylinder_crosscheck/`(summary.csv 4 行 + PNG);
  `cases/gates/INDEX.md` 加 C10 节。
- **裁决 c 登记**:`PHASE_SIX_CAPABILITY_BOUNDARY.md` §9 开放研究项加一行
  (tet medium 在 r_far 20/30 不收敛,M0.80 腿脆弱性,与 §5 敏感腿同族)。
- **波及修复**:`cylinder_o_grid_2d` 签名变更波及两个 phase-3 bench 脚本的
  解包(`run_hex_g0_single_var.py` ×2 行、`run_hex_q2_wall_accuracy.py` ×1 行,
  `_interior` 丢弃,语义不变 —— 见 §5 意外 b)。

## 2. 判据逐条(预注册 §3)

| 判据 | 结果 |
|---|---|
| C1 重生成:既有组几何差 + wake 存在 + 两次同参数 sha256 一致 | ✅ 既有组几何差 **0.0**(kill 阈值 1e-12);wake 面 62/72/82(三级);两次 sha256 逐位一致(coarse `a25c5af9…` / medium `bdd93a29…` / fine `c6de34cf…`);`cut_wake` 通过(1 TE 站位,74 主从对) |
| C2 C10 两腿在带内 + G-TEETH 两条 | ✅ 门 **5 passed / 1 skipped**(导出 opt-in);G-TEETH 见 §3 |
| C3 全套收口 820+Δ 闭合 | **留给主线收口**(本代理只跑单文件:F01 **23 passed**、A14+A18 **19 passed**、C10 5/1) |
| C4 旧 sha 引用(manifest/F01)不破 | ✅ F01 绿 |

## 3. G-TEETH(两条都实测)

1. **改坏 committed 值**:nl_hex 的 cp_rms +27 % ⇒ 锁红(指名列与相对差);
   还原 ⇒ 绿。
2. **改坏几何**:`LEVELS["medium"]` n_theta 96→90 重生成 ⇒ 锁红
   (cp_rms 相对差 5.18e-2 > 1e-6)。★ 注意分工:此时 4 条判据带测试仍绿
   (90 在噪声带内)—— **带内扰动正是锁的职责**,不是判据失灵。
   还原参数重生成:sha256 逐位回 `bdd93a291f…` ⇒ 门绿。

## 4. 两腿读数(committed `summary.csv`;带推导在门 docstring)

非升力(Γ=0,m_inf=0 Laplace 极限,`solve_newton_nonlifting` direct):

| 腿 | cp_rms | cp_max | cl | cd |
|---|---|---|---|---|
| nl_hex | 0.03654 | 0.06300 | −2.5e-5 | 7.1e-6 |
| nl_tet | 0.03210 | 0.05899 | +1.2e-5 | 1.4e-5 |

带:cp_rms ≤ 0.06、cp_max ≤ 0.11、|cl| ≤ 2e-3、cd ≤ 2e-3;两族交叉差
**0.0044** ≤ 0.017。

升力(Γ=0.6,`solve_laplace_lifting` gamma_fixed;hex `n_kutta=0` 短路成立):

| 腿 | cl | cl_exact | cl_rel |
|---|---|---|---|
| lift_hex | 0.597450 | 0.6(弦 2) | **0.425 %** |
| lift_tet(C06 族) | 1.186373 | 1.2(弦 1) | **1.136 %** |

共用带 ≤ 1.6 %;两族偏置差 **0.71 pp** ≤ 1.2 pp(登记带)。
**G-DOMAIN 的相反结果没有出现**:hex 族在非升力圆柱上落在带内,
NACA 的 −24 % 嫌疑**没有**在圆柱上复现 —— 嫌疑更偏向「翼型 hex 拓扑上
被触发的求解/分辨率问题」,生成器的圆柱半干净。

## 5. 常量零改动 + 意外清单

`git diff` 内**无任何既有判据常量变动**;新增常量是 C10 自己的标定带
(新文件,推导注释在门 docstring)。stats CSV 除新增 `n_tris_wake` 列外
逐位不动(几何不动的又一证据)。

意外:
a. **证据锁行数守卫初版常量写错**(18 → 17):`lift_tet` 行没有 `cd` 列
   (C06._one_level 不回 cd),锁的机器对缺失列跳过不计。修守卫常量,
   不是修数据 —— 守卫的防丢腿功能不变。
b. **签名波及**:phase-3 两个 bench 脚本按 4 元组解包。g0 已修并**跑通**
   (G0 PASS,单变量逐位性质在加 wake 后仍成立 —— 这正是 wake 只加组
   不动几何的独立验证);q2 已修解包,但该脚本另有**既有的 REPO 路径漂移**
   (`phases/p3/cases/...` 不存在,本轮之前就跑不到解包行),不修 —— 修
   路径会让它重写 tracked 的 `bench/gate_results/hex_q2_wall_accuracy.csv`,
   属重刷其他证据,超出本轮授权。**登记**:q2 脚本需要一轮自己的布局修复。
c. 主线交接原话让 §9 指向 `bench/studies/hex_cl_attribution/README.md`;
   核实该 README **没有** r_far 记录,真实出处是
   `20261012-1000-w25-w26-verdict.md` §4.2(顺带发现)与 §6(裁决 c)——
   按真实出处写。

## 6. 判定

**PASS**。裁决 b 落地(wake 补进、两拓扑对照门绿、G-TEETH 实测),裁决 c
登记进能力边界 §9。kill criterion 一条都没触发(几何差 0.0、hex 圆柱腿
全部收敛)。

## 7. 下一步(留给主线 / 裁决层)

1. 全套收口(预注册 C3):820+Δ 闭合,Δ 预期 +6(C10 的 5 门 + 1 导出腿)。
2. commit 全部本轮文件(含三个新 manifest)。
3. §9 工程清单的「W2.6(hex 族求解交叉门)」措辞:圆柱腿已由本门交付,
   NACA 腿停在归因(hex_cl_attribution,cl −24 %,缺陷仍在侧)—— 是否
   改写由主线定夺。
4. q2 bench 脚本的 REPO 路径漂移(§5b)待一轮自己的修复。
