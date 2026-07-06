# Per-phase demo cases

One self-contained, self-checking demo per completed roadmap phase. Each
`run_demo.py` regenerates every figure/CSV in its `results/` folder and
verifies the phase's acceptance numbers (exit code 0 = all checks pass;
the G1.6 open gate is an expected XFAIL). Conclusions and analysis live in
[docs/demo_report.md](../../docs/demo_report.md).

| Folder | Phase | Run | Runtime |
|---|---|---|---|
| `p0_infrastructure/` | P0 mesh infrastructure (G0.1–G0.4) | `python cases/demo/p0_infrastructure/run_demo.py` | ~5 s |
| `p1_laplace/` | P1 Laplace solver (V0, G1.1–G1.2; G1.4/DP1 negative results; G1.6 open) | `python cases/demo/p1_laplace/run_demo.py` | ~15 s |
| `p2_kutta_lifting/` | P2 wake cut + Kutta (G2.1–G2.5) | `python cases/demo/p2_kutta_lifting/run_demo.py` | ~20 s |
| `m0_meshgen/` | M0 quasi-2D meshing pipeline | `python cases/demo/m0_meshgen/run_demo.py` | ~10 s |

`_common.py` holds the shared chart style and the `CheckList` acceptance
recorder. `results/` contents (PNG + `summary.csv` + `checks.csv`) are
committed so the report's figures stay valid without rerunning.
