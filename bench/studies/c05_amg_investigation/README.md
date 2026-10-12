# c05_amg_investigation —— Q1:C05 喷管生产腿上 `precond="amg"` 为什么 40 步停在 |R|~3.9e-6

预注册:`docs/dev_phase_six/20261012-1000-research-amg-hex-prereg.md` §Q1。
性质:**研究/测量**,不建门、不改判据、不改库(库与判据零改动,本目录为唯一新增)。

## 结论:机制 ① —— 预条件子没近似到激波块(GMRES+amg 从第一步起就不动)

更准确地说(比预注册候选 ① 的字面机制更锐):**AMG 层级本身是健康的,
失败在于它被建在了不含激波耦合的 Picard Term-1 块上。** 证据链:

1. **Newton 层复现**(`newton_histories.csv`,W2.2b 配方复算):
   amg 40 步停在 |R|∞=3.906e-6(尾部 37→39 步 3.885e-6→3.906e-6,微升,
   与 C05 模块 docstring「降到 ~4e-6 后单调回升」一致),每步烧满
   600 次 GMRES 内迭代(n_gmres_total=24000=40×600,即**每步都没达到 EW η**);
   ilu 10 步到 9.2e-15(全程仅 58 次 GMRES),direct 10 步 9.2e-15。
2. **同一冻结点**(amg 腿末态,|R|∞=3.9136e-6,n_free=5174,wc=None ⇒
   GMRES 算子 ≡ J_ff,无低秩项)上取 J_ff 成功、splu 正常(无精确奇异,
   **未触发 kill**):
   - **判别式直接读数**(`gmres_freeze_point.csv`,rtol=1e-12,restart=60,
     600 内迭代 = 生产单步预算):**GMRES+amg 真残差相对降 0.994 —— 全程平的,
     从头部(2.25e-5→2.24e-5)到尾部不动**;预条件残差同样平
     (1.0166e-3→1.0153e-3 / 40 步)。对照:GMRES+ilu **2 次内迭代**收敛
     (x 与 direct 差 3.1e-13);**GMRES+amg 作用于 A_ff 本身也是 2 次收敛**
     (层级健康对照)。
   - ⇒ ① 成立(前几步就不动);② 排除(不是缓慢线性降;且同一 J_ff 上 ilu
     2 步收敛,cond2=2.50e5 远非病态);③ 排除(预条件残差与真残差一样平,
     判据没被骗)。
3. **机制闭合**(`freeze_point_metrics.csv`、`model_mismatch_profile.csv`):
   - J_ff 含一个近零软模态:σ_min=2.7462e-6(左/右两个逆 Gram 算子各自
     独立量到同值),σ_max=0.686,cond2=2.50e5。**右奇异矢(待解方向)峰值在
     x=11.4 = 激波处**;左奇异矢(残差方向)峰值在 x=6.4 = 喉部声速点。
   - 该软模态的预条件子反演误差 ||M(Jv)−v||/||v||:**amg = 0.950(完全反不动)**,
     ilu = 0.602(也不精确,但 ILU 把谱压紧后 Krylov 一两步就收拾掉单个软方向)。
   - AMG 被喂的矩阵 A_ff:对称(skew 1.2e-17)、对角全正、rho_t∈[0.678,1.273],
     pyamg 建出 4 层(opcx 1.079,335 聚合,B 维 1)——**「Picard 块失去
     M-矩阵性 ⇒ 聚合误读」这一子机制被证伪**:A_ff 确有 24.9% 正非对角元
     (非 M-矩阵),但层级在 A_ff 上工作得很好(对照腿 2 步收敛)。
   - 模型错配 ||J_ff−A_ff||_F/||J_ff||_F=0.0954,错配最大的前 200 行集中在
     **x∈[4.6,11.9]** —— 声速点到激波的超音速段(逐 x-bin 分布见
     `model_mismatch_profile.csv`)。全局 9.5% 的 Frobenius 差集中压在
     软模态上,这就是「看起来差不多的块」做不出预条件子的原因。

**一句话**:喷管生产腿的 amg 停滞 = ①;AMG 建在不含 Term-3 迎风/激波耦合的
Picard 块 A_ff 上,对 J_ff 的激波中心近零模态(σ_min=2.7e-6,峰在 x=11.4)几乎
无反演能力,GMRES 在生产预算(600 内迭代/步)内残差纹丝不动,Newton 步因此
恒近零步,|R| 停在 3.9e-6;ilu 直接分解 J_ff 本身,2 次内迭代解决。

**附带记录**(非本轮问题,但测量中在场):三条腿(含收敛的 ilu/direct)都打出
sigma frozen-in-transient 警告(amg 腿 last |dσ|=4.8e-3,tail n_shock
69;68;72;71;…),即熵修正在末态被冻结于瞬态 —— 它不影响本结论(ilu/direct 带着
同一个冻结 σ 收敛到机器零),但它是生产腿公有的已知瞬态冻结,登记于此。

## 文件

| 文件 | 内容 |
|---|---|
| `newton_histories.csv` | 三条生产腿的逐 Newton 步 |R|∞(leg,step,residual_inf) |
| `gmres_freeze_point.csv` | 冻结点上四条 GMRES 腿的逐内迭代真残差+预条件残差(precond,iter,true_res,precond_res) |
| `freeze_point_metrics.csv` | 全部标量读数(矩阵性、pyamg 层级、σ、软模态、反演误差、GMRES 汇总) |
| `model_mismatch_profile.csv` | ||(J_ff−A_ff) 行|| 的逐 x-bin 分布(错配定位) |
| `c05_amg_gmres.png` | 左:三条 Newton 腿;右:冻结点 GMRES 逐次残差(含 amg 预条件残差与 A_ff 对照) |
| `probe1_reproduce_and_extract.py` | 复现三腿 + 冻结点 J_ff/A_ff 提取 + 矩阵性 + 40 迭代 GMRES |
| `probe2_long_gmres_softmode.py` | 600 迭代长 GMRES + cond2 + 软模态(v1 的 svds-SM 不收敛,改 LU 预条件 eigsh) |
| `probe3_right_sv_mismatch.py` | 右奇异矢方向的反演误差 + ||J−A|| 错配定位 |

## 复算(判据 C3)

```
cd /home/lrz/codes/UP3D
PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
  conda run -n up3d --no-capture-output python bench/studies/c05_amg_investigation/probe1_reproduce_and_extract.py
# 依次再跑 probe2、probe3(各自 ~1–2 min;探针内 sys.path 钉死仓库根)
```

环境:conda `up3d`(python 3.11.15 / numpy 2.4.6 / scipy 1.17.1 / pyamg 5.3.0),
分支 `claude/phase-six-audit`。pyamg 层级经 `build_amg_preconditioner` 的钉种构造,
测量可复算(本轮 probe1/probe2 两次独立重跑 amg 腿,末态 |R| 分别为
3.906e-06 与 3.905879e-06,一致;提取点残差 3.913594e-06 与驱动器末步记录值的
小差来自末次 σ 刷新后的重估,两条链路均已在探针注释里写明)。
