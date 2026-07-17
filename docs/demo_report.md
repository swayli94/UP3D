# Phase demo report — evidence for completed phases (index)

**Since 2026-07-15 the per-phase evidence sections live in one file per track
under [demo_report/](demo_report/)** — split verbatim from this file:
[demo_report/track_p.md](demo_report/track_p.md) (solver),
[demo_report/track_m.md](demo_report/track_m.md) (meshing),
[demo_report/track_b.md](demo_report/track_b.md) (level-set wake; includes the
2026-07-15 M6-medium LS workflow demo). References of the form "demo_report §P4"
resolve via the table below. When a phase closes, add its demo section to the
matching track file and a row here.

**Scope.** One self-contained demo case per completed roadmap phase, designed as
*evidence that the phase's functionality works, is numerically stable, and
physically sensible* — not merely that tests pass. Each demo is a standalone
script with built-in acceptance checks against the roadmap gate criteria; its
figures and CSVs are committed under `cases/demo/<phase>/results/` so the
report's numbers and images always correspond to a reproducible run.

**Reproduce.** `python cases/demo/<phase>/run_demo.py` (headless, matplotlib
Agg). Exit code 0 = all checks pass; `results/checks.csv` holds the measured
value, criterion, and PASS/FAIL/XFAIL status per check.

**Honesty rule.** Negative results are documented as evidence, not hidden as
gaps: G1.6 is shown as a strict XFAIL with its root cause; P4's same-day
re-open/re-close lives in §P4's addenda; B8 closed as characterized-NOT-cured;
G13.3-transonic is recorded as NEGATIVE. **A claim without a committed artifact
is not evidence** (the 2026-07-13 audit rule).

| Phase | Demo | Checks | Verdict | Report section in |
|---|---|---|---|---|
| P0 mesh infrastructure | `cases/demo/p0_infrastructure/` | 4 PASS | closed, reproduced | [track_p](demo_report/track_p.md) |
| P1 Laplace solver | `cases/demo/p1_laplace/` | 9 PASS + 1 XFAIL (G1.6) | closed gates reproduced; G1.6 open by design | [track_p](demo_report/track_p.md) |
| P2 wake cut + Kutta | `cases/demo/p2_kutta_lifting/` | 11 PASS | closed, reproduced | [track_p](demo_report/track_p.md) |
| M0 quasi-2D meshing | `cases/demo/m0_meshgen/` | 6 PASS | closed, reproduced | [track_m](demo_report/track_m.md) |
| P3 subsonic compressible | `cases/demo/p3_subsonic/` | 14 PASS | closed, reproduced | [track_p](demo_report/track_p.md) |
| P4 transonic artificial density | `cases/demo/p4_transonic/` | 10 PASS | closed, reproduced (re-closed 2026-07-07; Picard-quality per the 2026-07-11 erratum) | [track_p](demo_report/track_p.md) |
| M1 swept-wing meshing (ONERA M6) | `cases/demo/m1_wing_mesh/` | 13 PASS | closed, reproduced (M1b self-similar ladder 2026-07-13) | [track_m](demo_report/track_m.md) |
| P5 3D validation (ONERA M6) | `cases/demo/p5_onera_m6/` | 16 PASS | closed 2026-07-08 (V6 < 1% deferred; see P13 for its O(h) closure trend) | [track_p](demo_report/track_p.md) |
| P6 surface-pressure recovery | `cases/demo/p6_surface_recovery/` | 6 PASS (incl. gated M6) | closed 2026-07-08 (sawtooth = recovery artifact) | [track_p](demo_report/track_p.md) |
| P7 differentiable walk flux | `cases/demo/p7_diff_flux/` | 7 PASS (incl. gated converged-field) | closed 2026-07-10 (FD 3–5e-10) | [track_p](demo_report/track_p.md) |
| P8 fully-coupled Newton | `cases/demo/p8_newton/` | 15 PASS (parts 2–3 gated) | closed 2026-07-11 (G8.1 + G8.2 + G8.3) | [track_p](demo_report/track_p.md) |
| P8 capability assessment | `cases/demo/p8_capability/` | 36 PASS (full matrix gated) | **evaluation demo, not a gate** (2026-07-11) | [track_p](demo_report/track_p.md) |
| P10 (partial) G10.2 continuation tolerance | `cases/demo/p10_newton_usability/` | split A/B verdict | G10.2 + G10.3 closed 2026-07-11; phase stays open (G10.1) | [track_p](demo_report/track_p.md) |
| P9 grid-convergence & accuracy-gap discrimination | `cases/demo/p9_grid_discrimination/` | 11 PASS + 3 XFAIL | closed 2026-07-11 | [track_p](demo_report/track_p.md) |
| **Track B** B1 cut-element identification | `tests/test_b1_cut_elements.py` (test-only) | 34 PASS | closed 2026-07-11 | [track_b](demo_report/track_b.md) |
| **Track B** B2 multivalued assembly | `tests/test_b2_multivalued.py` (test-only) | 17 PASS | closed 2026-07-11 | [track_b](demo_report/track_b.md) |
| **Track B** B3 + B4 lifting + TE Kutta | `cases/demo/b3_levelset_lifting/` | 13 demo PASS (+6, +8 tests) | closed 2026-07-12 | [track_b](demo_report/track_b.md) |
| **Track B** B5 far-field A/B | `cases/demo/b4p5_farfield/` | 9 demo PASS (+10 tests) | closed 2026-07-12 | [track_b](demo_report/track_b.md) |
| **Track B** B6 transonic (level-set) + LS Newton | `cases/demo/b6_transonic/` | 14 demo PASS (+9, +2 tests; +2, +2 Newton) | ◐ coarse gate ✓ 2026-07-12; the medium quantitative item closed by B15/GB15.4 | [track_b](demo_report/track_b.md) |
| **Track B** B7 ONERA M6 3D gate | `cases/demo/b7_onera_m6/` | 35 PASS | closed 2026-07-12 (M_max re-read 2026-07-14: honest main-field 1.392) | [track_b](demo_report/track_b.md) |
| P13 G13.1 tip/wake-edge characterization | `cases/demo/p13_tip_edge_singularity/` | 10 PASS | closed 2026-07-13 (1/√r, dΓ/dz-driven) | [track_p](demo_report/track_p.md) |
| P13 G13.2 spanwise loading taper | `run_taper_probe.py` + `run_taper_physics.py` | PASS | conforming fix closed 2026-07-13 | [track_p](demo_report/track_p.md) |
| P13 G13.3 ladder + third singularity | `run_g133_ladder.py` | 5/5 | flat-cap wall edge located (p=+0.32) → M5 | [track_p](demo_report/track_p.md) |
| **Track M** M5 rounded tip cap | `cases/demo/m5_round_tip/` | 9/9 | closed 2026-07-13 (seam crease O(h)) | [track_m](demo_report/track_m.md) |
| P13 G13.3 subsonic Richardson | `run_g133_roundtip.py` | 9/9 | earned 2026-07-13: p=2.31, cl_KJ(h→0)=0.2050 | [track_p](demo_report/track_p.md) |
| P13 G13.3 transonic | `run_g133_roundtip_transonic.py` + `_locate.py` | 5/5 | **NEGATIVE** 2026-07-13/14: round fine never reaches M0.84 (sharp tip-TE, amplified) | [track_p](demo_report/track_p.md) |
| **Track B** B8 tip-edge desingularization | `cases/demo/b8_tip_taper_ls/` + re-spec demos | 12/12 + 8/8 | closed 2026-07-14 **characterized-not-cured** (metric artifact + both cures negative) | [track_b](demo_report/track_b.md) |
| **Track B** B11 LS infrastructure | `cases/demo/b11_ls_infra/` | PASS | closed 2026-07-14 (unified post + ILU escape; AMG stalls on lifting) | [track_b](demo_report/track_b.md) |
| **Track B** B12 + B13 lagged-LU | `cases/demo/b12_lagged_lu/` + `b13_lagged_picard/` | 6/6 + 6/6 | closed 2026-07-14 (Newton 2.18×; lifting 6.55×) | [track_b](demo_report/track_b.md) |
| M6 medium LS workflow | `cases/demo/m6_medium_ls_workflow/` | 10/10 | demo, not a gate (2026-07-15): sub+transonic at conforming-comparable cost | [track_b](demo_report/track_b.md) |
| **Track B** B15 LS Newton ramp + freeze | `cases/demo/b15_ls_newton_ramp/` | 19/19 incl. gated M6 | closed 2026-07-15 (plateau gone, 3.5×; four errata) | [track_b](demo_report/track_b.md) |
| **Track A** A1 solver bottleneck study | `cases/analysis/a1_solver_bottleneck/` | 11/11 (2.5-D) + 4/4 (gated 3-D) | closed 2026-07-16: GA1.1–GA1.5; 3-D Newton is precond-bound (the 2.5-D seed headline does not transfer); GA1.5 reproduces G8.2/B15 digit-for-digit and found 4 harness defects | [track_a](demo_report/track_a.md) |
| **Track A** A2 TE/Kutta fidelity | `cases/analysis/a2_te_kutta_fidelity/` | 22/22 (zero-solve) + 4/4 (gated intervention) | closed 2026-07-17: GA2.1–GA2.5; S1 Γ(z) jitter = a probe-difference Kutta-target measurement artifact (fixed-Γ discriminator D=7.33/25.70 coarse/medium), not flow content; S2 TE Cp jump = potential-jump Kutta form error (34×/133× vs level-set) + P1 recovery artifact; fix routed to P14 (no `pyfp3d/` edits) | [track_a](demo_report/track_a.md) |
| **Track B** B9 wing-body cross-model | `cases/demo/b9_wingbody/` (+ guardrail `cases/analysis/b9_fuselage_guardrail/`) | 7 PASS + 1 XFAIL | closed 2026-07-17 (RE-SPEC'D, user-approved): LS (Picard) + conforming (NEW capability, Newton) on the M2 wing-body, M0.5 coarse+medium. ★ the two wake models AGREE to **cl_p 0.4% / cl_kj 0.6%** at medium (conf 0.2173/0.2188 vs LS 0.2165/0.2175, GB9.5 PASS; coarse 12.8% = resolution) — the wing-body analogue of P14's cross-model. GB9.1/9.2/9.3 ✓; **GB9.4 XFAIL** (fuselage lift 16-20%, resolution/model-sensitive ⇒ G1.6 fuselage-Cp error, band NOT moved); GB9.6 RECORDED (azimuthal Cp scatter median 0.0036/0.0022/0.0010, max grows at the poles). ★ LS uses PICARD: the committed LS Newton recipes all diverge on the wing-body — `neumann` is unbounded under the fuselage blockage (te_aux perfect 1.8e-8, 8 far-field fluid rows |R|≈84 in the never-exercised freestream-Newton path) | [track_b](demo_report/track_b.md) |
| **Track B** B16 LS Newton far-field aux pin | `cases/demo/b16_farfield_aux/` | 9 PASS + 1 XFAIL | closed 2026-07-18 (churn fix; lift-convergence OPEN) (NEW, user-directed; executes the B9 recorded follow-up): the wing-body LS-Newton churn is a **near-singular far-field aux block** — a wake sheet with no outflow clip leaves the outer nodes it crosses on near-singular wake-LS rows; at the freestream Picard state they hold garbage (|jump| **53.4** at x≥10 vs Γ̄ 0.0586), which Picard absorbs but Newton reads as the **8 far-field MAIN rows max\|R\|=84.457** (aux-block cond1 **O(1e19)→8.70e6** pinned). `farfield_aux="pin"` (default, mode-adaptive) ⇒ **coarse** freestream Newton **res 5.88e-14, 0 limited** (lift matches conforming 0.1%) where legacy churns at **7.95**/3690 limited. ★★ **GB16.4 XFAIL — UNRESOLVED non-convergence (user-flagged):** the {Newton-pin, LS-Picard, conforming} lift triangle does NOT close and flips with resolution — medium Newton-pin **0.1690** STALLS at res 7e-6, 22% below the Picard≈conforming pair (0.2165/0.2173, B9's 0.4%) ⇒ at least one path is not converged (medium Newton-pin non-converged, or B9's LS-Picard≈conforming a non-converged coincidence); analysis deferred. ⚠ the proposal's weld-vs-wake_ls mechanism was self-corrected. neumann byte-identical | [track_b](demo_report/track_b.md) |
| **P14** pressure-equality Kutta estimator | `cases/demo/p14_pressure_kutta/` (+ diagnostic `cases/analysis/p14_te_pressure_diag/` 20/20) | 28 PASS incl. gated M0.84 | 2026-07-17 CLOSED (opened + closed same day): A2's fix built and confirmed — S1+S2 both gone in one estimator swap — M0.84 Γ(z) roughness 0.0970→**0.0043** / 0.0365→**0.0024** (at/below the LS band), all-station raw TE Cp gap 0.2206→**0.0040** / 0.1585→**0.0024** (55×/67×; a same-day erratum corrected the first write-up's wrong-metric 0.318/0.228 baseline — A2 measured the gap two ways). ★ **Cross-model (V14.6)**: conforming-pressure 0.2776/0.2823 vs level-set 0.2772/0.2813 — the two independent wake models agree to **0.17%/0.36%** where the probe path was 4.5%/4.3% below LS ⇒ the conforming-vs-LS lift disagreement *was* the Kutta form. ★ **G14.7 ✓ re-specced at close** (user-arbitrated): the estimator swap moves cl_KJ +4.85% off the *probe-path* G8.2 locks and ONTO the level-set oracle (0.15%/0.34%), closing **69% of P9's 0.019 gap** — the move is the finding, re-locked against level-set. Discriminator D 7.33→**1.80** (A2's inconclusive zone, recorded not rounded). ★ **V14.7 self-correction**: the TE Cp *spike* was asserted untouched (A2's "shared recovery artifact") — **measured false**: 0.1143→**0.0533**, below level-set's 0.0743; the conforming excess over LS was Kutta-form error too, with a ~0.05 shared recovery floor remaining | [track_p](demo_report/track_p.md) |

> Track-P renumber (2026-07-08, then 2026-07-11 ×2): P6 = surface recovery;
> P7 = differentiable flux (Newton prereq); P8 = fully-coupled Newton;
> P9 = grid-convergence discrimination; P10 = Newton generality/continuation;
> P11 = curved wall elements; P12 = backlog; P13 (2026-07-13, appended) =
> tip/wake-edge singularity.
> **Track-B renumber (2026-07-12 ×2, 2026-07-13):** a new **B4** (TE control
> volume) was inserted, then the half-integer IDs were regularized away — the
> far-field gate is now **B5** (was B3.5, then B4.5; its demo dir keeps the old
> `b4p5_` name), transonic is **B6**, ONERA M6 3D is **B7**; 2026-07-13 inserted
> B8 = LS tip taper (old B8 multi-wake → B9, old B9 curved wake → B10).
> See roadmap/track_b.md for the full mapping.

---

## Cross-phase summary

- **Functionality**: every closed gate's headline number is reproduced from
  scratch by the demos (MMS slope 1.96, CG 8→11→14, cl −0.82%, Γ-lift
  cross-check 0.015%, spanwise ratio 2.05, cylinder slope 1.02).
- **Numerical stability**: machine-zero consistency checks (V0, G2.1, G2.2,
  G2.5a) hold on the largest committed meshes, the linear solver is
  mesh-independent, and the Kutta outer loop converges in 2 updates.
- **Physical soundness**: stagnation/suction structure on the sphere,
  smooth TE flow, and three mutually independent lift routes agreeing —
  physics cross-checks, not code self-consistency.
- **Open item**: G1.6 remains open (11.6% vs 2%); its Option C acceptance
  re-spec is the open P1 task, with Option A conclusively ruled out by the
  oracle demos.
- **P3 additions**: compressibility validated against the classical
  correction band (sphere 0.32% vs PG; airfoil cl inside [PG, KT]);
  M∞ → 0 is bit-identical to the P1/P2 Laplace drivers; assembly is
  colored-`prange` with precomputed geometry (~160× hot reassembly) and
  every solve is bit-reproducible run-to-run (seeded AMG setup).

---

