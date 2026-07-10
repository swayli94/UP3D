# Per-phase demo cases

One self-contained, self-checking demo per completed roadmap phase. Each
`run_demo.py` regenerates every figure/CSV in its `results/` folder and
verifies the phase's acceptance numbers (exit code 0 = all checks pass;
the G1.6 open gate is an expected XFAIL). Conclusions and analysis live in
[docs/demo_report.md](../../docs/demo_report.md).

Heavy transonic parts run only under `PYFP3D_TRANSONIC_GATES=1` (skipped by
default); runtimes below are the default (light) path.

| Folder | Phase | Run | Runtime |
|---|---|---|---|
| `p0_infrastructure/` | P0 mesh infrastructure (G0.1–G0.4) | `python cases/demo/p0_infrastructure/run_demo.py` | ~5 s |
| `p1_laplace/` | P1 Laplace solver (V0, G1.1–G1.2; G1.4/DP1 negative results; G1.6 open) | `python cases/demo/p1_laplace/run_demo.py` | ~15 s |
| `p2_kutta_lifting/` | P2 wake cut + Kutta (G2.1–G2.5) | `python cases/demo/p2_kutta_lifting/run_demo.py` | ~20 s |
| `m0_meshgen/` | M0 quasi-2D meshing pipeline | `python cases/demo/m0_meshgen/run_demo.py` | ~10 s |
| `p3_subsonic/` | P3 subsonic compressible (G3.1–G3.3) | `python cases/demo/p3_subsonic/run_demo.py` | ~1 min |
| `p4_transonic/` | P4 transonic artificial density (G4.1–G4.3) | `python cases/demo/p4_transonic/run_demo.py` | ~3 min (+heavy) |
| `m1_wing_mesh/` | M1 swept-wing mesh (ONERA M6) | `python cases/demo/m1_wing_mesh/run_demo.py` | ~30 s |
| `p5_onera_m6/` | P5 3D validation (ONERA M6, G5.1–G5.2) | `python cases/demo/p5_onera_m6/run_demo.py` | heavy |
| `p6_surface_recovery/` | P6 surface-Cp recovery / sawtooth removal (G6.1–G6.4) | `python cases/demo/p6_surface_recovery/run_demo.py` | ~6 min (+heavy M6) |
| `p7_diff_flux/` | P7 frozen-selection ∂ρ̃/∂φ of the walk flux (G7.3) | `python cases/demo/p7_diff_flux/run_demo.py` | ~10 s (+heavy converged-field check) |

Future phases (P8 Newton, P9 curved wall elements) will add `p8_newton/`,
`p9_curved_walls/` when they run.

`_common.py` holds the shared chart style and the `CheckList` acceptance
recorder. `results/` contents (PNG + `summary.csv` + `checks.csv`) are
committed so the report's figures stay valid without rerunning.
