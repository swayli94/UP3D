# 20261012-1300 — **预注册:cylinder_hex 补 wake 组 + 非升力对照(裁决 b);tet medium 不收敛立项登记(裁决 c)**

**性质**:预注册(判据与 kill criterion 事前写死),随代码提交。
**对象**:`cases/meshes/cylinder_hex_2.5d/generate_cylinder_hex.py`(补 wake 组)、
该族三级网格(**重生成,纪律 11 勘误**)、新建
`tests/C/test_C10_hex_cylinder_crosscheck.py`。
**由来**:使用者 2026-10-12 裁决 b「cylinder_hex 补 wake 组,同时包括非升力对照」
与 c「tet medium 在 r_far 20/30 不收敛立项研究」。
**已核实**:没有任何测试读 cylinder_hex 的已提交 .msh(A18 读的是生成器函数),
重生成不破任何测试;重生成仍走纪律 11 勘误清单。

---

## 1. 目的

把 W2.6 停掉的圆柱腿补上:cylinder_hex 族补 wake 组(让它能承载升力环量),
并把**非升力对照**(Γ=0,有精确 Cp)一并纳入 —— 两拓扑( hex 结构族 vs
cylinder_2.5d tet 族)在同一精确解上对照,攻击的是「生成器有系统缺陷时
所有门共享该缺陷」这个单点。

**预先写下的失败标准**:
- 重生成的网格若**不是**逐位/结构可控的(边界组变了但几何变了多少没量)⇒ FAIL;
- 对照门的带若写在**当前实测差**之外而没过噪声地板 ⇒ FAIL(带必须 ≥
  两族实现噪声合成,算式写明);
- 补 wake 若动了既有 wall/farfield 组的几何 ⇒ FAIL(只允许**加**组)。

## 2. 改动内容(事前登记)

1. **生成器补 wake 组**:cylinder_hex 的生成器加 wake 片(参照
   naca0012_hex 族的做法 —— 读它的 generate_naca0012_hex.py 怎么写 wake),
   既有 wall/farfield/symmetry 组几何不动;重生成三级(coarse/medium/fine),
   每组生成报告:节点/单元数、边界组计数、与旧网格的几何差(应有 wake 面、
   其余逐位或近逐位)。
2. **新门 C10**(`tests/C/test_C10_hex_cylinder_crosscheck.py`):
   - **非升力对照腿**:hex 与 tet(cylinder_2.5d medium)各解 Γ=0 圆柱,
     cl/Cp vs 精确解(C02 已有的精确解),带 = 两族实现噪声合成
     (C06 的实测 + 本轮实测,算式写明);
   - **升力腿**(wake 组就位后):hex vs tet 的 cl vs 精确解,同上带;
   - 证据进 `cases/gates/C10_hex_cylinder_crosscheck/`(新目录 + INDEX.md),
     **blas_threads 出处列进锁**。
   - G-TEETH:改坏 committed 值 ⇒ 红;改坏几何(生成器参数)⇒ 红。
3. **裁决 c 登记**:tet medium 在 r_far 20/30 不收敛(M0.80 腿脆弱性)写进
   能力边界文书的开放研究项 + 本轮判定文件。

## 3. 判据(事前写死)

| # | 判据 | 期望 |
|---|---|---|
| C1 | 重生成的网格:既有组几何差报告 + wake 组存在;重生成可复现(同参数两次 sha256 一致) | 实测 |
| C2 | C10 两腿在带内;G-TEETH 两条实测 | 绿 + 红 |
| C3 | 全套收口:820+Δ / 0 / 32 / 2 闭合 | 0 failed |
| C4 | 纪律 11:重生成网格的旧 sha 引用(manifest/F01 机制)不破 | F01 绿 |

**G-DOMAIN**:C2 的相反结果 = hex 族在非升力圆柱上也偏差超带 ⇒
与 NACA 的 −24 % 同族(系统缺陷)⇒ 红是**发现**,登记,不许放宽带。

## 4. kill criterion

- 补 wake 若让既有组的几何动超过 1e-12 ⇒ 停,生成器改动错了;
- 若 hex 圆柱腿不收敛 ⇒ 停,记录(与 hex fine 的 NACA 不收敛同族登记)。

## 5. 明确不做

- 不动 naca0012_hex 族(它已有 wake);不动 C02 的 tet 圆柱判据;
- 不建 NACA 的 hex 门(归因已完成,建门形状等裁决:能力边界门 vs 生成器审查);
- 不改任何既有判据带。
