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

- P0 mesh infrastructure — `cases/demo/p0_infrastructure/` — 4 PASS — closed, reproduced — [track_p](demo_report/track_p.md)
- P1 Laplace solver — `cases/demo/p1_laplace/` — 9 PASS + 1 XFAIL (G1.6) — closed gates reproduced; G1.6 open by design —
  [track_p](demo_report/track_p.md)
- P2 wake cut + Kutta — `cases/demo/p2_kutta_lifting/` — 11 PASS — closed, reproduced — [track_p](demo_report/track_p.md)
- M0 quasi-2D meshing — `cases/demo/m0_meshgen/` — 6 PASS — closed, reproduced — [track_m](demo_report/track_m.md)
- P3 subsonic compressible — `cases/demo/p3_subsonic/` — 14 PASS — closed, reproduced — [track_p](demo_report/track_p.md)
- P4 transonic artificial density — `cases/demo/p4_transonic/` — 10 PASS — closed, reproduced (re-closed 2026-07-07;
  Picard-quality per the 2026-07-11 erratum) — [track_p](demo_report/track_p.md)
- M1 swept-wing meshing (ONERA M6) — `cases/demo/m1_wing_mesh/` — 13 PASS — closed, reproduced (M1b self-similar ladder 2026-07-13) —
  [track_m](demo_report/track_m.md)
- P5 3D validation (ONERA M6) — `cases/demo/p5_onera_m6/` — 16 PASS — closed 2026-07-08 (V6 < 1% deferred;
  see P13 for its O(h) closure trend) — [track_p](demo_report/track_p.md)
- P6 surface-pressure recovery — `cases/demo/p6_surface_recovery/` — 6 PASS (incl. gated M6) — closed 2026-07-08 (sawtooth = recovery
  artifact) — [track_p](demo_report/track_p.md)
- P7 differentiable walk flux — `cases/demo/p7_diff_flux/` — 7 PASS (incl. gated converged-field) — closed 2026-07-10 (FD 3–5e-10) —
  [track_p](demo_report/track_p.md)
- P8 fully-coupled Newton — `cases/demo/p8_newton/` — 15 PASS (parts 2–3 gated) — closed 2026-07-11 (G8.1 + G8.2 + G8.3) —
  [track_p](demo_report/track_p.md)
- P8 capability assessment — `cases/demo/p8_capability/` — 36 PASS (full matrix gated) — **evaluation demo, not a gate** (2026-07-11) —
  [track_p](demo_report/track_p.md)
- P10 (partial) G10.2 continuation tolerance — `cases/demo/p10_newton_usability/` — split A/B verdict —
  G10.2 + G10.3 closed 2026-07-11; phase stays open (G10.1) — [track_p](demo_report/track_p.md)
- P9 grid-convergence & accuracy-gap discrimination — `cases/demo/p9_grid_discrimination/` — 11 PASS + 3 XFAIL —
  closed 2026-07-11 — [track_p](demo_report/track_p.md)
- **Track B** B1 cut-element identification — `tests/test_b1_cut_elements.py` (test-only) — 34 PASS —
  closed 2026-07-11 — [track_b](demo_report/track_b.md)
- **Track B** B2 multivalued assembly — `tests/test_b2_multivalued.py` (test-only) — 17 PASS — closed 2026-07-11 —
  [track_b](demo_report/track_b.md)
- **Track B** B3 + B4 lifting + TE Kutta — `cases/demo/b3_levelset_lifting/` — 13 demo PASS (+6, +8 tests) —
  closed 2026-07-12 — [track_b](demo_report/track_b.md)
- **Track B** B5 far-field A/B — `cases/demo/b4p5_farfield/` — 9 demo PASS (+10 tests) — closed 2026-07-12 —
  [track_b](demo_report/track_b.md)
- **Track B** B6 transonic (level-set) + LS Newton — `cases/demo/b6_transonic/` — 14 demo PASS (+9, +2 tests; +2, +2 Newton) —
  ◐ coarse gate ✓ 2026-07-12; the medium quantitative item closed by B15/GB15.4 — [track_b](demo_report/track_b.md)
- **Track B** B7 ONERA M6 3D gate — `cases/demo/b7_onera_m6/` — 35 PASS — closed 2026-07-12 (M_max re-read 2026-07-14:
  honest main-field 1.392) — [track_b](demo_report/track_b.md)
- P13 G13.1 tip/wake-edge characterization — `cases/demo/p13_tip_edge_singularity/` — 10 PASS — closed 2026-07-13 (1/√r, dΓ/dz-driven) —
  [track_p](demo_report/track_p.md)
- P13 G13.2 spanwise loading taper — `run_taper_probe.py` + `run_taper_physics.py` — PASS — conforming fix closed 2026-07-13 —
  [track_p](demo_report/track_p.md)
- P13 G13.3 ladder + third singularity — `run_g133_ladder.py` — 5/5 — flat-cap wall edge located (p=+0.32) → M5 —
  [track_p](demo_report/track_p.md)
- **Track M** M5 rounded tip cap — `cases/demo/m5_round_tip/` — 9/9 — closed 2026-07-13 (seam crease O(h)) —
  [track_m](demo_report/track_m.md)
- P13 G13.3 subsonic Richardson — `run_g133_roundtip.py` — 9/9 — earned 2026-07-13: p=2.31, cl_KJ(h→0)=0.2050 —
  [track_p](demo_report/track_p.md)
- P13 G13.3 transonic — `run_g133_roundtip_transonic.py` + `_locate.py` — 5/5 — **NEGATIVE** 2026-07-13/14:
  round fine never reaches M0.84 (sharp tip-TE, amplified) — [track_p](demo_report/track_p.md)
- **Track B** B8 tip-edge desingularization — `cases/demo/b8_tip_taper_ls/` + re-spec demos — 12/12 + 8/8 —
  closed 2026-07-14 **characterized-not-cured** (metric artifact + both cures negative) — [track_b](demo_report/track_b.md)
- **Track B** B11 LS infrastructure — `cases/demo/b11_ls_infra/` — PASS — closed 2026-07-14 (unified post + ILU escape;
  AMG stalls on lifting) — [track_b](demo_report/track_b.md)
- **Track B** B12 + B13 lagged-LU — `cases/demo/b12_lagged_lu/` + `b13_lagged_picard/` — 6/6 + 6/6 —
  closed 2026-07-14 (Newton 2.18×; lifting 6.55×) — [track_b](demo_report/track_b.md)
- M6 medium LS workflow — `cases/demo/m6_medium_ls_workflow/` — 10/10 — demo, not a gate (2026-07-15):
  sub+transonic at conforming-comparable cost — [track_b](demo_report/track_b.md)
- **Track B** B15 LS Newton ramp + freeze — `cases/demo/b15_ls_newton_ramp/` — 19/19 at close-out; 17/20 under the first B20 re-baseline;
  **20/20 refreshed by B22 after the B21 fix (2026-07-19)** — closed 2026-07-15 (plateau gone; four errata). ★ Erratum trail:
  the post-B20 stall (GB20.7 "real capability loss") was actually the B21/N1 freeze-capture patch gap —
  the refreshed committed record is γ **0.088343**, M_max **2.4818**, res 9.0e-14, 0 lim/1 flr, freeze armed 6/6 levels with 0 reverts,
  **511 s = 4.51×** vs the committed Picard; see the B21/B22 sections — [track_b](demo_report/track_b.md)
- **Track B** B14 Schur-eliminated aux + AMG structural preconditioner — `cases/demo/b14_schur_precond/` — 7/7 incl.
  gated M6 coarse+medium; 5/7 under the first B20 re-baseline; **7/7 refreshed by B22 (2026-07-19)** —
  closed 2026-07-17 ( `precond="schur"` , the A1 precond bottleneck gone). B22-refreshed post-B21 numbers:
  medium ramp lagged 505 s vs schur **345 s = 1.47×**, precond share 45.8% → **1.8%**, γ **0.088343** both arms (\|Δγ\| 8.6e-13), 0
  fallbacks; anchor constants re-pinned to the B21 state (γ 0.088343 / M_max 2.4818). ⚠ this index row was MISSING from B14's close-out
  until B22 added it (the D9 close-out-debt class) — [track_b](demo_report/track_b.md)
- **Track A** A1 solver bottleneck study — `cases/analysis/a1_solver_bottleneck/` — 11/11 (2.5-D) + 4/4 (gated 3-D) —
  closed 2026-07-16: GA1.1–GA1.5; 3-D Newton is precond-bound (the 2.5-D seed headline does not transfer);
  GA1.5 reproduces G8.2/B15 digit-for-digit and found 4 harness defects — [track_a](demo_report/track_a.md)
- **Track A** A2 TE/Kutta fidelity — `cases/analysis/a2_te_kutta_fidelity/` — 22/22 (zero-solve) + 4/4 (gated intervention) —
  closed 2026-07-17: GA2.1–GA2.5; S1 Γ(z) jitter = a probe-difference Kutta-target measurement artifact (fixed-Γ discriminator D=7.33/25.70
  coarse/medium), not flow content; S2 TE Cp jump = potential-jump Kutta form error (34×/133× vs level-set) + P1 recovery artifact;
  fix routed to P14 (no `pyfp3d/` edits) — [track_a](demo_report/track_a.md)
- **Track A** A3 inspection response / C1 Jacobian check — `cases/analysis/c1_ls_jacobian_fd/` — 3 probes + a 3-decade eps sweep
  (measurement, not pass/fail) — closed 2026-07-18: GA3.1–GA3.6; ★★ **C1 CONFIRMED** — the LS Newton Jacobian is not dR/dφ on mixed-side
  plain elements (targeted 1.146e-01 vs control 6.33e-10, and **eps-independent**: 1.532e-01 at 1e-6/1e-7/1e-8, max/min 1.00 ⇒ a missing
  term, not FD noise) ⇒ **quasi**-Newton in 3-D; R untouched so every converged state and gate number stands. RECORDED not fixed.
  Plus C2/C3 backported, reader C4/C5, C6/C7/P1/T1/T2/F0, 17/17 docs findings — [track_a](demo_report/track_a.md)
- **Track B** B9 wing-body cross-model — `cases/demo/b9_wingbody/` (+ guardrail `cases/analysis/b9_fuselage_guardrail/` ) —
  7 PASS + 1 XFAIL — closed 2026-07-17 (RE-SPEC'D, user-approved): LS (Picard) + conforming (NEW capability, Newton) on the M2 wing-body,
  M0.5 coarse+medium. ★ the two wake models AGREE to **cl_p 0.4% / cl_kj 0.6%** at medium (conf 0.2173/0.2188 vs LS 0.2165/0.2175, GB9.5
  PASS; coarse 12.8% = resolution) — the wing-body analogue of P14's cross-model. GB9.1/9.2/9.3 ✓; **GB9.4 XFAIL** (fuselage lift 16-20%,
  resolution/model-sensitive ⇒ G1.6 fuselage-Cp error, band NOT moved); GB9.6 RECORDED (azimuthal Cp scatter median 0.0036/0.0022/0.0010,
  max grows at the poles). ★ LS uses PICARD: the committed LS Newton recipes all diverge on the wing-body —
  `neumann` is unbounded under the fuselage blockage (te_aux perfect 1.8e-8, 8 far-field fluid rows — R —
  ≈84 in the never-exercised freestream-Newton path) — [track_b](demo_report/track_b.md)
- **Track B** B16 LS Newton far-field aux pin — `cases/demo/b16_farfield_aux/` — 9 PASS + 1 XFAIL — closed 2026-07-18 (churn fix;
  lift-convergence OPEN) (NEW, user-directed; executes the B9 recorded follow-up): the wing-body LS-Newton churn is a **near-singular
  far-field aux block** — a wake sheet with no outflow clip leaves the outer nodes it crosses on near-singular wake-LS rows;
  at the freestream Picard state they hold garbage ( — jump — **53.4** at x≥10 vs Γ̄ 0.0586), which Picard absorbs but Newton reads as the
  **8 far-field MAIN rows max\|R\|=84.457** (aux-block cond1 **O(1e19)→8.70e6** pinned). `farfield_aux="pin"` (default, mode-adaptive) ⇒
  **coarse** freestream Newton **res 5.88e-14, 0 limited** (lift matches conforming 0.1%) where legacy churns at **7.95**/3690 limited.
  ★★ **GB16.4 XFAIL — UNRESOLVED non-convergence (user-flagged):** the {Newton-pin, LS-Picard, conforming} lift triangle does NOT close and
  flips with resolution — medium Newton-pin **0.1690** STALLS at res 7e-6, 22% below the Picard≈conforming pair (0.2165/0.2173, B9's 0.4%) ⇒
  at least one path is not converged (medium Newton-pin non-converged, or B9's LS-Picard≈conforming a non-converged coincidence);
  analysis deferred. ⚠ the proposal's weld-vs-wake_ls mechanism was self-corrected. neumann byte-identical.
  **★ GB16.4 RESOLVED by B17 (below).** — [track_b](demo_report/track_b.md)
- **Track B** B17 far-field aux pin_gamma (resolves GB16.4) — `cases/demo/b17_farfield_pin_gamma/` — 3 coarse PASS + gated medium —
  closed 2026-07-18 (NEW, user-directed; executes the GB16.4 open follow-up): **GB16.4 was a BC-modelling error, not a non-convergence.**
  B16's freestream pin forced the outflow wake jump to **0**, removing the circulation the wake physically carries out (−22% at medium;
  the coarse "match" was a jump=0/legacy-garbage cancellation). **Decisive:** an independent Picard-pin converges cleanly (res 7.5e-8) to
  the SAME medium **0.1691** the Newton-pin "stalls" at (0.1690) ⇒ both solvers agree per-BC. Fix = `farfield_aux="pin_gamma"` (aux=host
  φ∞−side·γ, jump→γ, the new default on both solvers): triangle closes MONOTONE to conforming — coarse **0.2087**, medium **0.2117**
  (Picard)/**0.2115** (Newton; **0.2114** post-B20, and the medium Newton trajectory now converges to \|R\|~1e-13 —
  the junction churn was the B20 mixed-plain contamination), both agree 0.1%. GB17.1 ring jump 53→0→γ + **B9 coarse-12.8%-was-contamination
  erratum**; GB17.2 cl_p≡cl_KJ move together (not a post artifact; plotted sectional cl is Γ-based);
  GB17.6 vortex brackets from +2.5% (does not close the gap). ★ the far-field conditioning, the outflow circulation, and the junction churn
  are three orthogonal issues B16 conflated — [track_b](demo_report/track_b.md)
- **Track B** B18 wing-body transonic (M0.84) — **B27 refresh 2026-07-20** — `cases/demo/b18_wingbody_transonic/` — 8 gates PASS —
  closed 2026-07-18 (NEW, user-directed; executes the GB16.6 debt); **refreshed 2026-07-20 (B27)**: the old "LS junction-limited
  (closed-negative)" story is **RETIRED** — B25/B26 measured that the junction pocket is the B23 inboard free-edge singularity and that
  clipping the wake sheet to the conforming fragment topology ( `inboard_clip` , B25) heals it, so **the LS ceiling is now co-located with
  the conforming ceiling**. **Conforming** unchanged the reference: coarse **M0.84 cl_p 0.2617**, medium **M0.79 cl_p 0.2579** strict, cl(M)
  0.2173/0.2321/**0.2483**/0.2579 @ M0.50/0.65/**0.75 (new)**/0.79. **Level-set**: the A side (no clip) is still pocket-limited (medium dies
  0.5125 class (a): the pocket erupts at 0.55, Mmax 13.1 > freeze_max_clamped=8; coarse dies 0.84 class (b)) —
  and its climbing past the committed B18 anchors (died 0.50/0.55) is the **B21/B22 freeze-capture fix** (B26 T1 finding), not physics
  drift; the **C side (+clip) REACHES coarse M0.84 (cl_p 0.2542) and medium M0.7625** (dies 0.775 class (b), dying peak at the wing TIP =
  P13 class, junction corridor clean) — the residual limiter is the same high-M Newton/shock class as conforming's.
  **Cross-model upgraded**: M0.5 (2.6%, B9/B17) + **M0.65 medium 2.4% (PASS ≤5%)** + M0.75 medium 2.5% (recorded) + coarse M0.6 (2.1%
  C-side, under-resolved — the old 0.2% A-side row retired as a fortuitous pocket-state agreement); all gaps sit in the ~2.5% B17 cl_p/cl_kj
  convention band. Consistency (B27 gates): conforming legs bit-reproduce the committed B18 anchors (GB27.1), LS A/C legs bit-reproduce the
  committed B26 anchors (GB27.2) — 336/336 in `cases/analysis/b27_b18_demo_refresh/results/g27_consistency.csv` .
  ★ the pre-B27 sections PNG was silently EMPTY (section_cp_curve tuple→dict API drift, swallowed by `except: pass` ) —
  fixed in the refresh. GB18.1/18.2/18.3(M0.65) PASS + GB18.3(M0.75)/18.4/18.5 RECORDED; cl_fus: conf 0.0423 @M0.79, C-side 0.0781 with
  out-band 0.0565 (×2, P11 input). **B29 refresh (same day, user-adjudicated B28 §6)**: flat-fragment adopted as the LS production
  config — C side = clip + flat sheet ( `sheet_direction=(1,0,0)` ); M0.5 LS anchors 0.2115/0.2184; medium ceiling **0.7625→0.775**
  (live peak M3.98 @ wing tip); cross-model gaps **0.5/1.1/1.1 %** (M0.5/0.65/0.75, was 2.6/2.4/2.5); GB18.5 live flat cl_fus
  **0.0382** (band −0.0006 / out 0.0388) @0.7875 vs conf 0.0423 — the ×2 out-band reading retired; 8/8 PASS —
  [track_b](demo_report/track_b.md)
- **Track B** B19 LS-Newton Jacobian exactness (3-D) — `cases/analysis/c1_ls_jacobian_fd/` — GB19.1–19.6 ✓ (GB19.4 recorded NEGATIVE) —
  closed 2026-07-18: TWO Leg-A defects (DOF column maps + gradient factors), targeted probe **1.146e-01 → 1.33e-08** with the
  ε-discriminator flipped (ε-independent → ~1/ε); R bit-identical ( `git stash` A/B per fix) so no converged result moves; GB19.4:
  NO convergence gain (+3.6% wall) — the plateau is B15 selection churn; Leg B/GB19.6 measured the mixed-plain side-density contamination
  (spurious supersonic q² 3.2229 vs 1.3379, 45.3% ρ error on 252 elements) and routed it to B20 — [track_b](demo_report/track_b.md)
- **Track B** B20 mixed-plain main-field density (adopted permanently) — `cases/analysis/c1_ls_jacobian_fd/` + re-baselined
  B7/B9/B15/B16/B17/B18 demo CSVs — GB20.1–20.6 ✓, GB20.7 answered 2026-07-19 (superseded by B21 same day) —
  closed 2026-07-18, knob removed (user-arbitrated), re-baselined 2026-07-19: every moved 3-D number went B20's way (B7 M_max 1.453→1.392;
  B16 legacy limited 3690→11; B18 wing-body res 6.8e-5→1.1e-13; M6 coarse ramp 0.7875-stall→M0.84 converged) EXCEPT the M6-medium ramp,
  which stalled at M0.6625 — GB20.7 called it "a real capability loss (freeze_tol axis)"; **B21 then found the true mechanism:
  the stall was B20's own patch gap in `freeze_side_state` (N1), and fixing it restores M0.84** — [track_b](demo_report/track_b.md)
- **Track B** B21 freeze-capture alignment (N1) — `cases/analysis/c1_ls_jacobian_fd/` ( `run_n1_freeze_capture.py` ,
  `n1_freeze_fix_sweep.csv` ) — GB21.1–21.3 ✓ — closed 2026-07-19 (executes the 2026-07-19 Kimi-inspection N1 finding):
  `freeze_side_state` captured the frozen selection on the UNPATCHED side field — the one B20 consumer the patch missed (probe:
  83+9 selection differences vs live, all aux-touching mixed-plain). One-line fix ⇒ the committed M6-medium recipe reaches **M0.84 again**
  (γ **0.088343**, M_max 2.4818, res 9.0e-14, 0 lim/1 flr, **515 s** — faster and cleaner than pre-B20's 657 s/3 clamped),
  freeze_tol-insensitive across 1e-3/1e-5. **GB20.7 overturned**; "contamination was a stabiliser" retired;
  GB15.4 capability clause stands (small numeric re-baseline; B15/B14 demo refresh = recorded follow-up).
  3-D capture lock in `test_b15_ls_newton_freeze.py` , verified failing pre-fix — [track_b](demo_report/track_b.md)
- **Track B** B22 evidence refresh + 3-D anchor locks + re-baseline process rule — refreshed `cases/demo/b15_ls_newton_ramp/` +
  `cases/demo/b14_schur_precond/` ; locks `tests/test_b22_ls_3d_anchors.py` — B15 demo **20/20**, B14 demo **7/7**;
  2 gated anchor tests — closed 2026-07-19 (executes B21's recorded follow-up + the Kimi N3 finding + the §2 process recommendation):
  LS npz caches deleted, zero `cached` lines, both demos re-solved green on the B21 state (medium γ 0.088343 / M_max 2.4818 / 511 s;
  coarse ramp γ 0.084931 / M_max 1.3684); **the N3 gap is closed** — gated anchor locks now re-solve the committed M6 coarse+medium ramps
  and assert m_final/γ/M_max/clamps absolutely, so the next silent re-baseline fails the suite; CLAUDE.md workflow step 5 + agent-rules
  discipline #11 gained the re-baseline erratum checklist; next-phase priority analysis in
  `docs/analysis/next_phase_priorities_2026-07-19.md` (recommends P11) — [track_b](demo_report/track_b.md)
- **Track B** B23 wing-body junction discriminator — `cases/analysis/b23_junction_discriminator/` — pre-registered D-campaign, committed
  CSV/PNG — closed 2026-07-19: the junction pocket is **lift/wake-coupled, NOT geometric** — α=0 clean at both levels (Mmax 0.66/0.64,
  cl_fus ≈ 0 self-check), pocket appears with α and grows superlinearly (medium α=3.06 corrM **14.66**, 104 supersonic elements, peak BEHIND
  the fuselage x=2.13 @ z≈z_junc). Attribution: the wake sheet's inboard boundary ends at the junction station (q≥0) —
  the **P13-class free-edge singularity transplanted inboard**. Routes: (b)-1 → B24 (closed negative), (b)-2 → B25 (the cure).
  P11 close-out input delivered — [track_b](demo_report/track_b.md)
- **Track B** B24 wake inboard-end waterline extension (negative) — `cases/analysis/b24_wake_inboard_end/` —
  pre-registered E1, committed CSV/PNG — closed 2026-07-19: the pocket **follows the free edge** (B1 flush moves the peak past x_tail every
  leg — mechanism re-confirmed), but both extension variants trade the singularity for equal-or-worse forms (B1 medium α=3.06 corrM **78.56
  NON-converged**; B3 offset-cone migrates back near-field every leg). Decision-tree exit 3: "extension-class insufficient", **the (b)-1
  route is CLOSED** → back to (b)-2 = B25 — [track_b](demo_report/track_b.md)
- **Track B** B25 inboard fragment clip (the CURE) — `cases/analysis/b25_inboard_fragment_clip/` + `meshgen/fuselage.py:make_inboard_clip` +
  `wake/cut_elements.py` — pre-registered F1 v2.1, all primary criteria decisive — closed 2026-07-19:
  the sheet's inboard boundary moves to the fuselage surface / symmetry plane (= conforming fragment topology);
  **default None ⇒ bit-identical**. Medium α=3.06: corridor corrM **14.66 → 0.63**, n_sup **88 → 0**, cl_p **+0.38 %** (within [A, oracle
  0.2173]), \|Δγ\| 0.37 %, root te_jump 0.28 %, α=0 inert, 56 outer ≤ 1.5×A, sliver min dihedral 11.0°.
  One secondary guardrail (out-of-band cl_fus carryover +135 %) oracle-attributed to the flat-vs-tilted sheet model difference —
  recorded non-blocking, P11 watch item. Chain closed: B23 → B24 → B25. Tests `TestInboardFragmentClip` (4) + M2/B1 locks —
  [track_b](demo_report/track_b.md)
- **Track B** B26 post-cure LS transonic ceiling re-measured — `cases/analysis/b26_ls_transonic_ceiling/` —
  pre-registered G1, B26-A (ceiling lift) — closed 2026-07-20: same-code A/C (default vs `inboard_clip` ), B18 recipe frozen verbatim, ~69
  min. **C medium m_last 0.50 → 0.7625** (five loose rungs + 0.7625 strict res 2.6e-11; dies 0.775);
  **C coarse 0.82 → 0.84 REACHED** (strict res 6.9e-11). Death cause flips (a)-pocket-rejection → **(b)-class high-M Newton stall with the
  peak at the WING TIP** (medium M4.18 @ z=1.20, corridor corrM 1.07 clean) = the conforming 0.80+ stall class, no longer the junction.
  cl_p same-trend vs conforming (0.2542@0.84 vs 0.2617; 2–4 % low = B17 convention gap). Supersedes the B18 LS-ceiling story —
  [track_b](demo_report/track_b.md)
- **Track B** B27 B18 demo refresh (LS legs resurrected) — `cases/analysis/b27_b18_demo_refresh/` + refreshed
  `cases/demo/b18_wingbody_transonic/` — GB27.1/27.2/27.3(0.65) PASS; 27.3(0.75)/27.4/27.5 RECORDED —
  closed 2026-07-20 (full re-solve ~1 h 39 min, no cache): conforming legs bit-reproduce the committed B18 anchors;
  LS A/C legs bit-reproduce B26 committed (**336/336 bit-identical**); checks.csv **8/8 PASS**. Transonic cross-model now exists —
  **M0.65 2.4 % PASS (≤5 %)**, M0.75 2.5 % RECORDED; gap flat across Mach (2.6/2.4/2.5 %) = one ~2.5 % B17 convention band. Independent:
  old coarse-0.60-cross 0.2 % row retired (cured C-side 0.2133, gap 2.1 %); the silently-empty sections PNG (older `section_cp_curve` API
  drift) fixed; conf medium **0.75 NEW point 0.2483** strict. Façade: "junction-limited" → "**post-cure LS ceiling co-located with
  conforming** (coarse 0.84 = 0.84; medium 0.7625 ≈ 0.79)"; Track V sheet-topology prerequisites in place —
  [track_b](demo_report/track_b.md)
- **Track B** B29 flat-fragment = the wing-body LS production config — `cases/demo/b18_wingbody_transonic/` (+
  `cases/analysis/b28_cl_fus_flat_sheet/VERDICT.md` §6, user-adjudicated) — closed 2026-07-20 (same branch as B28; no
  `pyfp3d/` change): B18 LS production side C = `inboard_clip` + `sheet_direction=(1,0,0)` (NEW `ls_flat_*` caches; the B26
  tilted `ls_C_*` stale; A side tilted kept as the historical comparison; conforming bit-reproduce). Medium ceiling 0.7625 →
  **0.775** (dies 0.7875, a+dm; live peak M3.98 @ wing tip z=1.20 — GB18.4's C side now live); cross-model **M0.5 0.5 % /
  M0.65 1.1 % PASS (≤5 %) / M0.75 1.1 %** (was 2.6/2.4/2.5 %); GB18.5 live flat decomposition cl_fus **0.0382** (band −0.0006,
  out-band 0.0388, poles 0.0007) @0.7875 vs conf 0.0423 @0.79 — the B26 tilted ×2 out-band reading retired (B28 position
  sensitivity); M0.5 LS anchors re-pinned 0.2087/0.2117 → 0.2115/0.2184; checks.csv **8/8 PASS**; `test_b9_wingbody_ls` on
  the production wiring (5/5) — [track_b](demo_report/track_b.md)
- **P11** curved wall-adjacent elements (sphere leg) — `cases/demo/p11_curved_walls/` — 14 PASS + 2 XFAIL (G11.1/G11.2 recorded negatives) —
  2026-07-19 CLOSED (opened + closed same day, user-directed): the DP1 curved-element route measured **NEGATIVE** —
  a verified curved wall layer (planar null test ΔA ≡ 0 bitwise; quadrature = P1 reference to 1.3e-15;
  deg2/deg3 A/B 5.5e-9) moves medium sphere Cp only **11.56%→11.33%** (= the G1.4 oracle ceiling); the pre-registered superparametric risk
  fired (mapped-P1 linear-reproduction deviation **O(h)**, max 0.138 coarse). ★★ **G1.6 re-attributed**:
  the h_min sweep's order collapse (0.88/0.56/0.42, replicated exactly by the committed script) is the **fixed-bulk-mesh pollution floor**
  (E8: far-mesh-only refinement at h_min=0.03 → φ_wall 3.17× lower, argmax r=1.53→wall, order 1.89);
  a structured icosphere shell with the SAME flat facets converges at order 1.67/1.98 and hits 2.14% max Cp at h≈0.036 ⇒ medium's 11.6% ≈
  **intrinsic P1-field capability at h=0.08** (geometric-crime share ≈0.2 pp). G1.6 xfail stays; route fork (Option C re-spec with measured
  passing form / isoparametric P2 layer / accept) = user's call — [track_p](demo_report/track_p.md)
- **P14** pressure-equality Kutta estimator — `cases/demo/p14_pressure_kutta/` (+ diagnostic `cases/analysis/p14_te_pressure_diag/` 20/20) —
  28 PASS incl. gated M0.84 — 2026-07-17 CLOSED (opened + closed same day): A2's fix built and confirmed —
  S1+S2 both gone in one estimator swap — M0.84 Γ(z) roughness 0.0970→**0.0043** / 0.0365→**0.0024** (at/below the LS band), all-station raw
  TE Cp gap 0.2206→**0.0040** / 0.1585→**0.0024** (55×/67×; a same-day erratum corrected the first write-up's wrong-metric 0.318/0.228
  baseline — A2 measured the gap two ways). ★ **Cross-model (V14.6)**: conforming-pressure 0.2776/0.2823 vs level-set 0.2772/0.2813 —
  the two independent wake models agree to **0.17%/0.36%** where the probe path was 4.5%/4.3% below LS ⇒ the conforming-vs-LS lift
  disagreement *was* the Kutta form. ★ **G14.7 ✓ re-specced at close** (user-arbitrated): the estimator swap moves cl_KJ +4.85% off the
  *probe-path* G8.2 locks and ONTO the level-set oracle (0.15%/0.34%), closing **69% of P9's 0.019 gap** —
  the move is the finding, re-locked against level-set. Discriminator D 7.33→**1.80** (A2's inconclusive zone, recorded not rounded).
  ★ **V14.7 self-correction**: the TE Cp *spike* was asserted untouched (A2's "shared recovery artifact") —
  **measured false**: 0.1143→**0.0533**, below level-set's 0.0743; the conforming excess over LS was Kutta-form error too, with a ~0.05
  shared recovery floor remaining — [track_p](demo_report/track_p.md)

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

