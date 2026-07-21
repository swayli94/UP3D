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
| `p8_newton/` | P8 fully-coupled (φ_red, Γ) Newton (G8.1–G8.3) | `python cases/demo/p8_newton/run_demo.py` | ~2 min (+heavy) |
| `p8_capability/` | P8 capability assessment (not a gate; 36/36) | `python cases/demo/p8_capability/run_demo.py` | heavy |
| `p9_grid_discrimination/` | P9 grid convergence / accuracy-gap discrimination (G9.1–G9.3) | `python cases/demo/p9_grid_discrimination/run_demo.py` | heavy |
| `p10_newton_usability/` | P10 continuation usability (G10.2/G10.3) | `python cases/demo/p10_newton_usability/run_demo.py` | heavy |
| `p11_curved_walls/` | **P11** curved wall elements: measured NEGATIVE + G1.6 re-attribution (14 PASS + 2 XFAIL) | `python cases/demo/p11_curved_walls/run_demo.py` | ~4 min warm; +~8 min first run (regenerates the gitignored sphere sweep meshes) |
| `p13_tip_edge_singularity/` | **P13** tip / wake-edge singularity — see the four scripts below | (below) | (below) |
| `b3_levelset_lifting/` | Track B B3 lifting solve, implicit Kutta | `python cases/demo/b3_levelset_lifting/run_demo.py` | ~1 min |
| `b4p5_farfield/` | Track B **B5** far-field A/B (dir name kept from the old numbering) | `python cases/demo/b4p5_farfield/run_demo.py` | ~2 min |
| `b6_transonic/` | Track B B6 transonic on the level-set path | `python cases/demo/b6_transonic/run_demo.py` | heavy |
| `b7_onera_m6/` | Track B B7 ONERA M6 3D dual-mesh (35/35) | `python cases/demo/b7_onera_m6/run_demo.py` | heavy |
| `b8_tip_taper_ls/` | Track B B8 tip taper on the level-set path (characterized-not-cured) | `python cases/demo/b8_tip_taper_ls/run_demo.py` | heavy |
| `b9_wingbody/` | Track B **B9** wing-body cross-model (LS Picard vs conforming Newton; **8 PASS** post-B28 GB9.4 re-spec — out-band cross-model ≤15%, medium 7.0%) | `python cases/demo/b9_wingbody/run_demo.py` | heavy |
| `b11_ls_infra/` | Track B B11 unified post infrastructure (bitwise legacy locks) | `python cases/demo/b11_ls_infra/run_b11_unified_post.py` | ~1 min |
| `b12_lagged_lu/` | Track B B12 lagged-LU on the LS Newton | `python cases/demo/b12_lagged_lu/run_demo.py` | heavy |
| `b13_lagged_picard/` | Track B B13 lagged-LU on the LS Picard | `python cases/demo/b13_lagged_picard/run_demo.py` | heavy |
| `b14_schur_precond/` | Track B **B14** Schur-eliminated aux + AMG (7/7 incl. gated M6) | `python cases/demo/b14_schur_precond/run_demo.py` | heavy |
| `b15_ls_newton_ramp/` | Track B **B15** LS Newton freeze-ramp (the strict transonic LS route) | `python cases/demo/b15_ls_newton_ramp/run_demo.py` | heavy |
| `b16_farfield_aux/` | Track B **B16** far-field aux pin — the LS-Newton churn fix (9 PASS + 1 XFAIL, resolved by B17) | `python cases/demo/b16_farfield_aux/run_demo.py` | heavy |
| `b17_farfield_pin_gamma/` | Track B **B17** `farfield_aux="pin_gamma"` (jump=γ; 3 coarse PASS + gated medium) | `python cases/demo/b17_farfield_pin_gamma/run_demo.py` | heavy |
| `b18_wingbody_transonic/` | Track B **B18** wing-body transonic M0.84 — **8 gates PASS** (B27 refresh / B29 flat / **B32 taper**: conforming coarse M0.84 + **medium climb reaches M0.84** cl_p 0.2738 0 clamps (taper); LS+`inboard_clip`+flat reaches the SAME ceiling site: coarse 0.84 / medium 0.775; cross-model M0.65 0.3% PASS + M0.75 0.2%) | `python cases/demo/b18_wingbody_transonic/run_demo.py` | heavy |
| `m5_round_tip/` | Track M **M5** round tip cap (9/9) | `python cases/demo/m5_round_tip/run_demo.py` | ~2 min |
| `m6_medium_ls_workflow/` | Track M/B M6 medium level-set workflow | `python cases/demo/m6_medium_ls_workflow/run_demo.py` | heavy |
| `p14_pressure_kutta/` | **P14** pressure-equality Kutta estimator (28 checks) | `python cases/demo/p14_pressure_kutta/run_demo.py` | heavy |

*(Table completed in A3 2026-07-18: it had stopped at `b7_onera_m6` and was
missing 13 directories that exist on disk. When a phase closes, add its row
here — see CLAUDE.md workflow step 5.)*

**P13 (`p13_tip_edge_singularity/`) is the one folder with several scripts** — the
phase has three gates and each carries its own acceptance set:

| Script | Gate | Checks |
|---|---|---|
| `run_demo.py` | G13.1 — characterize the tip/wake-edge singularity | 10 passed |
| `run_taper_probe.py` | G13.2 — the spanwise loading taper regularizes it | 16 passed + 2 xfailed |
| `run_taper_physics.py` | G13.2 — what the taper costs, and that it stays LOCAL | 10 passed + 1 xfailed |
| `run_g132_transonic.py` | G13.2 — the transonic clause: is G9.1's fine mesh a discrete solution now? (**heavy, ~45 min**) | see below |
| `run_g133_ladder.py` | G13.3 — what still blocks 3D grid convergence (the flat tip cap) | 5 passed |

The P13 xfails are *documented negative results*, not breakage: `tanh_half` sits
exactly on the regularization criterion's knife edge **and** is disqualified for
unbounded support, and `vanish_smooth` at r_c = 0.03 is under-regularized — which
is why the shipped r_c is 0.05.

`_common.py` holds the shared chart style and the `CheckList` acceptance
recorder. **When several scripts share one `results/` folder they must pass
distinct `CheckList.report(out, fname)` names** — they all defaulted to
`checks.csv` and silently overwrote each other's committed evidence (found and
fixed 2026-07-13; hence `checks_g131.csv`, `checks_g132_probe.csv`, …).
`results/` contents (PNG + CSVs) are committed so the report's figures stay valid
without rerunning; the heavy `*.npz` solution caches are gitignored.
