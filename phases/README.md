# `phases/` —— phase 1 与 phase 2 的归档

**这里放的是已结束阶段的材料。** phase 3 的工作面在仓库的原位置
(`pyfp3d/`、`tests/`、`bench/`、`cases/`、`docs/`)。

整理于 **2026-08-10**,基准提交 **`d224223`**(phase two 的 PR #25 合并点)。

## 分类口径(量出来的,不是凭印象)

相位归属 = **文件首次提交的日期**,边界 **2026-07-28**(phase two 纲领建立日)。
按此,整理前的 1606 个被跟踪文件里 phase 1 约 1029、phase 2 约 524。

## ★ 什么留在原位:被保留代码引用的东西

规则是**可测的**:**如果一个仍然存在的测试或脚本读它,它就是 phase 3 仍然需要的,于是留在原位。**
按这条规则,`cases/demo` 的 33 个子目录里 **9 个留下**、`cases/analysis` 的 33 个里 **10 个留下** ——
例如 `cases/demo/p11_curved_walls/`,`tests/test_laplace_sphere.py` 要读它来跑 **G1.6 Option C 的活判据**。

留在原位的还有:`pyfp3d/` 全部、`cases/reference_data/`(外部真值,永不编辑)、`cases/meshes/`、
73 个测试文件(**全部 conforming 锚点**)、`bench/` 的 7 个指标脚本闭包 + `bench/gate_results/`
(能力边界按路径引用那些 CSV)、以及 `docs/` 的 overview / design / design_track_v / agent-rules /
inspection 与 6 个 phase-3 要用的 phase-two 文件。

## 目录

| | 内容 |
|---|---|
| `p1/` | phase 1(Track P/M/B/V/A):24 个 demo + 23 个 analysis 证据链、phase-1 的 roadmap 与 demo_report、`design_track_b.md`、`analysis/`、`archive/`、**18 个纯 level-set 测试** |
| `p2/` | phase 2:45 个 bench 脚本 + `s1_duct/`、**85 个轮次文件** |

## ★★ 归档脚本的可运行性:如实说

移动会打断脚本找路的方式,已**机械修好并验证**的有三类:
1. **58 处 `parents[N]`** 深度(两种写法:`Path(__file__).resolve().parents[k]` 与
   `HERE.parents[k]`,二者相差 1);
2. **22 处指向已移动目录的路径字面量**(只改真的移动了的,`cases/meshes/` 与留下的 demo 不动);
3. **41 个归档 bench 脚本的 `gate_results` 指向** —— 那个目录**留在** `bench/`,所以归档脚本
   是**横向**去取,不是往自己下面找;每一个都验证过解析结果落在含 CSV 的目录上。

**没有逐一验证到"能真的跑起来"。** 残留的少数交叉引用依赖被 import 的兄弟脚本来设置 `sys.path`,
静态检查看不透。所以口径是:

- **证据是已提交的 CSV / PNG**(项目纪律:"没有提交产物的说法不算证据"),它们**都还在,且仍被跟踪**;
- **脚本是出处(provenance)**;若要一个**保证能跑**的 phase 1/2 树,用一条命令:

```bash
git worktree add ../up3d-prereorg d224223     # 整理之前的完整状态
```

## 为什么是跟踪而不是 gitignore

最初的要求是把归档 gitignore。**改为跟踪**(使用者 2026-08-10 裁决),理由是:
gitignore 意味着**新 clone 里没有这些文件** —— 而那会让 phase 1/2 的每一个结论在 HEAD 上
**失去提交产物**,与本项目的头号纪律直接冲突。目录整洁的效果不变,证据仍在。
