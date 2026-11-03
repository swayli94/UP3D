# baseline-provenance arms: comparison

source: `/home/lrz/code/UP3D/bench/results/baseline_provenance_20261104`; values compared at full precision (`repr`).
Scope: `NUMBA_CPU_NAME=generic` changes numba/LLVM JIT codegen only.

## Arms

| arm | tests | failures | errors | skipped | numba hits | numba misses | dispatchers | wall |
|---|---|---|---|---|---|---|---|---|
| native_cold | 120 | 11 | 0 | 12 | 0 | 50 | 61 | 1009 s |
| native_warm | 120 | 12 | 0 | 12 | 34 | 0 | 61 | 1045 s |
| generic_cold | 120 | 13 | 0 | 12 | 0 | 50 | 61 | 1066 s |
| generic_warm | 120 | 11 | 0 | 12 | 34 | 0 | 61 | 927 s |

## Failed nodes per arm

| node | native_cold | native_warm | generic_cold | generic_warm |
|---|---|---|---|---|
| `tests.C.test_C05_nozzle_quasi1d.TestCommittedEvidenceIsLoadBearing::test_fresh_run_reproduces_the_committed_summary` | FAIL | FAIL | FAIL | FAIL |
| `tests.C.test_C07_karman_trefftz.TestCommittedEvidenceIsLoadBearing::test_fresh_run_reproduces_the_committed_summary` | FAIL | FAIL | FAIL | FAIL |
| `tests.C.test_C08_ringleb.TestArtificialDissipationPricedAgainstTruth::test_m_crit_is_the_cheaper_lever` | pass | pass | FAIL | pass |
| `tests.C.test_C08_ringleb.TestCommittedEvidenceIsLoadBearing::test_fresh_run_reproduces_the_committed_summary` | FAIL | FAIL | FAIL | FAIL |
| `tests.D.test_D05_euler_naca0012.TestCommittedEvidenceIsLoadBearing::test_matches_committed_summary` | FAIL | FAIL | FAIL | FAIL |
| `tests.D.test_D06_euler_rae2822.TestCommittedEvidenceIsLoadBearing::test_matches_committed_summary` | FAIL | FAIL | FAIL | FAIL |
| `tests.D.test_D08_rans_naca0012.TestCommittedEvidenceIsLoadBearing::test_matches_committed_summary` | pass | FAIL | FAIL | pass |
| `tests.D.test_D11_naca0012_experiment_bias.TestCommittedEvidenceIsLoadBearing::test_fresh_run_reproduces_the_committed_summary` | FAIL | FAIL | FAIL | FAIL |
| `tests.D.test_D12_rae2822_experiment_bias.TestCommittedEvidenceIsLoadBearing::test_fresh_run_reproduces_the_committed_summary` | FAIL | FAIL | FAIL | FAIL |
| `tests.D.test_D13_ibl_vs_xfoil.TestCommittedEvidenceIsLoadBearing::test_fresh_run_reproduces_the_committed_summary` | FAIL | FAIL | FAIL | FAIL |
| `tests.D.test_D13_ibl_vs_xfoil.TestRecorded::test_surface_asymmetry_is_missing` | FAIL | FAIL | FAIL | FAIL |
| `tests.D.test_D14_alpha0_transonic_nonuniqueness.TestCommittedEvidenceIsLoadBearing::test_matches_committed_summary` | FAIL | FAIL | FAIL | FAIL |
| `tests.D.test_D15_hex_naca_boundary.TestCommittedEvidenceIsLoadBearing::test_matches_committed_summary` | FAIL | FAIL | FAIL | FAIL |

## Pairwise (all locked values)

| pairing | values | bit-identical | max rel diff |
|---|---|---|---|
| ISA cold (native_cold vs generic_cold) | 351 | 147 | 2.010e-01 |
| ISA warm (native_warm vs generic_warm) | 351 | 146 | 2.428e+01 |
| cache native (native_cold vs native_warm) | 351 | 235 | 2.093e-01 |
| cache generic (generic_cold vs generic_warm) | 351 | 236 | 2.787e+01 |

## Per gate: max rel diff vs committed, per arm

| gate | values | native_cold | native_warm | generic_cold | generic_warm |
|---|---|---|---|---|---|
| cases/gates/C05_nozzle_quasi1d | 24 | 2.602e-02 | 2.602e-02 | 1.768e-01 | 1.768e-01 |
| cases/gates/C07_karman_trefftz | 24 | 4.526e-04 | 4.526e-04 | 4.526e-04 | 4.526e-04 |
| cases/gates/C08_ringleb | 30 | 9.708e-04 | 5.439e-04 | 9.054e-02 | 2.001e-03 |
| cases/gates/D05_euler_naca0012 | 108 | 1.023e-04 | 1.023e-04 | 5.954e-05 | 5.954e-05 |
| cases/gates/D06_euler_rae2822 | 36 | 3.549e-02 | 3.549e-02 | 2.132e-02 | 1.935e+01 |
| cases/gates/D08_rans_naca0012 | 28 | 2.488e-02 | 2.222e-01 | 2.140e-01 | 1.141e-01 |
| cases/gates/D11_naca0012_experiment_bias | 36 | 1.469e-04 | 1.469e-04 | 8.639e-05 | 8.639e-05 |
| cases/gates/D12_rae2822_experiment_bias | 28 | 9.031e-02 | 9.031e-02 | 4.519e-02 | 2.657e+01 |
| cases/gates/D14_alpha0_transonic_nonuniqueness | 22 | 2.619e-04 | 2.619e-04 | 5.384e-05 | 5.384e-05 |
| cases/gates/D15_hex_naca_boundary | 15 | 1.135e-06 | 1.135e-06 | 3.751e-06 | 3.751e-06 |

## Per gate: pairwise bit-identical / max rel diff

| gate | ISA cold | ISA warm | cache native | cache generic |
|---|---|---|---|---|
| cases/gates/C05_nozzle_quasi1d | 4/24, 1.98e-01 | 4/24, 1.98e-01 | 6/24, 3.61e-12 | 6/24, 1.60e-12 |
| cases/gates/C07_karman_trefftz | 1/24, 6.66e-09 | 1/24, 6.66e-09 | 24/24, 0.00e+00 | 24/24, 0.00e+00 |
| cases/gates/C08_ringleb | 0/30, 9.02e-02 | 0/30, 1.46e-03 | 6/30, 4.27e-04 | 6/30, 8.30e-02 |
| cases/gates/D05_euler_naca0012 | 82/108, 4.28e-05 | 82/108, 4.28e-05 | 107/108, 5.50e-10 | 108/108, 0.00e+00 |
| cases/gates/D06_euler_rae2822 | 21/36, 5.61e-02 | 21/36, 2.01e+01 | 32/36, 1.01e-09 | 31/36, 1.90e+01 |
| cases/gates/D08_rans_naca0012 | 8/28, 2.01e-01 | 8/28, 1.42e-01 | 8/28, 2.09e-01 | 8/28, 1.30e-01 |
| cases/gates/D11_naca0012_experiment_bias | 8/36, 6.05e-05 | 8/36, 6.05e-05 | 8/36, 3.01e-11 | 8/36, 1.19e-11 |
| cases/gates/D12_rae2822_experiment_bias | 9/28, 1.24e-01 | 9/28, 2.43e+01 | 9/28, 1.15e-09 | 9/28, 2.79e+01 |
| cases/gates/D14_alpha0_transonic_nonuniqueness | 5/22, 2.87e-04 | 4/22, 2.87e-04 | 20/22, 3.38e-10 | 21/22, 1.63e-10 |
| cases/gates/D15_hex_naca_boundary | 9/15, 4.89e-06 | 9/15, 4.89e-06 | 15/15, 0.00e+00 | 15/15, 0.00e+00 |

## D13 `runs` fixture (npz), pairwise

| pairing | arrays | bit-identical | max abs diff |
|---|---|---|---|
| ISA cold | 38 | 20 | 1.431e-01 |
| ISA warm | 38 | 19 | 1.000e+00 |
| cache native | 38 | 19 | 1.000e+00 |
| cache generic | 38 | 20 | 1.410e-01 |
- native_cold: coarse max asym 1.2935390829359936; medium max asym 1.0367791077921389
- native_warm: coarse max asym 1.2022292349421009; medium max asym 1.036114571261331
- generic_cold: coarse max asym 1.1523819079972915; medium max asym 1.037377629891366
- generic_warm: coarse max asym 1.1549252558970544; medium max asym 1.0342227075361436
