"""生产求解配方 —— 算例层的常量，不是库的一部分。

★ 使用者裁决 2026-08-24（选项 a）：`NEWTON_M6_RECIPE` 原先住在
`tests/test_p8_newton.py` 里，而 `bench/` 与 `cases/demo/` 从那里 import 它 ——
**方向反了**：一个生产配方住在测试文件里，十几个脚本依赖一个测试文件的名字。
配方是**基准/算例层**的东西（放进 `pyfp3d/` 会让库知道 ONERA M6），所以落在这里。

谁定义，谁都从这里读：`tests/`、`bench/`、`cases/demo/` 一律 `from bench.recipes import ...`。
"""

NEWTON_M6_RECIPE = dict(
    dm=0.05, dm_min=0.01, freeze_tol=1e-6, intermediate_tol=1e-5,
    newton_kw=dict(freeze_refresh_max=8, precond="amg", n_newton_max=60,
                   n_picard_seed=0, farfield_spanwise_gamma=True),
)
