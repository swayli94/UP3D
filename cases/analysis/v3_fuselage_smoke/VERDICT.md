# GV3.3 VERDICT — honest FAIL (0 PASS / 2 FAIL / 7 RECORDED)

Date: 2026-07-22. Branch: kimi/track-v3-loose-coupling.
Binding text: docs/roadmap/track_v.md GV3.3; pre-registration:
PRE_REGISTRATION.md (committed before the first execution). Evidence:
results/ (gv3_3_summary.csv is the machine-readable verdict; this file
is the narrative).

**Headline: the closed-body coupling machinery is now STABLE — after
three debug rounds (exactly-singular Newton at the aft pole →
transpiration-sink runaway → Goldstein separation crash) the final
scheme (tail-band Dirichlet pin + transpiration masking on the pinned
band + FP non-convergence guard) runs all 10 outer iterations without a
numerical event, and the MID-BODY axisymmetry is excellent
(σ/μ(δ*) 0.018–0.068, crossflow ratio ~1e-6). Both binding bands still
FAIL, concentrated in two localized zones: the immediate post-trip ring
(x/L 0.20–0.33) and the tail cone (x/L ≥ 0.82). The loose loop does not
converge — the tail-cone transpiration keeps growing (ṁ_max ×5.7 over
k = 5→10), the measured loose-coupling instability that motivates V4.**

## Case conditions (as pre-registered)

FuselageParams BoR (L = 4.0295), full-2π revolve, full R_FAR = 25 MAC
sphere, groups {wall, farfield}; coarse mesh (26.5k nodes / 139k tets /
12k wall triangles; stats in cases/meshes/fuselage_bor/results/). M 0.3,
α = 0, Re = 3.0e6 per body length (re_chord 7.443e5). x_tr = 0.05,
stag_band_frac = 0.05. Driver: compressible Picard (non-lifting) + V2
body_source_rhs via run_loose_coupling; ω = 1.0, n_outer_max = 10,
tol_ds 1e-3.

## Debug history (three rounds, all measured — docs/temp/v3_bor_diag.py)

1. **First smoke — crash.** On a closed surface the BL characteristics
   converge to the aft pole: there is no natural outflow boundary (an
   airfoil has the TE), and the IBL Newton Jacobian is exactly singular
   there ("Matrix is exactly singular"; garbage state then crashed the
   closure with a ZeroDivisionError).
2. **Round 1 — 5 % tail band pinned to turbulent seeds.** Singularity
   removed, but the frozen fat-seed δ* over the tail-cone convergence
   zone acts as a transpiration SINK: ṁ_sink −1.05 (k=1) → −3.3 (k=3) →
   −5.7 (k=4), tail u_e 1.29 → 1.51 → 1.80, positive feedback gain > 1,
   FP non-converged at k=4 (φ ~ 2e23), blow-up at k=5 (q ~ 1e25). A
   Picard-type loose loop cannot converge for ANY ω when the feedback
   gain exceeds 1 (the Veldman lesson — exactly the V4 motivation).
3. **Round 2 — narrow pin (pole + first ring only).** The tail-cone BL,
   left free, genuinely separates and the Newton march hits the
   Goldstein separation singularity head-on — crash even earlier (k≈2).
4. **Round 3 (FINAL) — tail-band pin + transpiration masking + FP
   guard.** The pin covers the last 5 % of body length (blocks the
   Goldstein march); the pinned band's frozen-seed δ* is masked out of
   the transpiration source (boundary data generates no ṁ — the tiny
   net-source imbalance is absorbed by the Dirichlet far field); a
   non-converged FP solve now raises instead of feeding garbage forward.
   No numerical event in 10 outer iterations. Airfoil path untouched
   (outflow_pin_surf = None there); changes confined to coupling.py.

## Assertions

- **FAIL — GV3.3(a) azimuthal σ/μ(δ*): worst 0.5533 at x/L = 0.940**
  (band ≤ 0.15 at every window station). 12 of 63 window stations
  exceed the band, in two localized zones:
  - post-trip ring x/L ∈ [0.20, 0.33]: σ/μ up to 0.389 at 0.285. δ* is
    small there (~3e-3) and the coarse mesh's azimuthal imprint through
    the instantaneous trip (x_tr = 0.05) dominates the ratio;
  - tail cone x/L ≥ 0.82: σ/μ rises 0.08 → 0.55 toward the pin edge —
    the growing free/pinned δ* discontinuity under meridian convergence
    (see loop record below).
  The mid-body x/L ∈ [0.34, 0.82] — the actual axisymmetric-collapse
  test — passes everywhere: σ/μ 0.018–0.068 (43 stations).
- **FAIL — GV3.3(b) crossflow: max|B|/max|A| = 0.2631,
  max|Cτ2|/max|Cτ1| = 0.2295** (band ≤ 0.05). Both maxima sit in the
  tail cone; over the body the ratio is ~1e-6–1e-4 (panel
  gv3_3_panels_coarse.png, bottom-right), consistent with the V1/V2
  local-basis regression levels (crossflow leakage 1.8e-4/1.6e-3, suite
  green). Not a basis defect — a tail-zone state defect.

## RECORDED

- (b) max|Ψ| = 2.367e3; max|DS2|/|DS1| = 1.657e-2/5.609e-2;
  max|CF2|/|CF1| = 4.623/4.515.
- (c) 653 inflow-pinned nodes (nose ring + tail band) / 1330 stagnation
  (A4 u_e-recovery) band nodes.
- (d) wall |ΔCp| on/off: max 1.6560 (at the pin-band edge, where ṁ is
  masked) / mean 0.0147. 20 nodes with δ* < 0 (clipped to the physical
  floor; counted per iteration in the history). Station-mean H(x):
  max 3.19 at x/L = 0.030 (nose laminar zone) vs 1.41 on the cylinder —
  **no tail H rise** in the meridian mean: the pin+mask does what the
  pre-registration's (d) escape hatch describes (indicated tail
  separation masked, not chased; δ* sign change absent inside the
  window, so no window exclusion was needed).
- **Loose loop: 10/10 outer iterations, ω = 1.0, NOT converged** — and
  not benignly: ṁ_max 0.26 → 0.27 (k=1–5, nearly stationary) then
  0.35 → 0.51 → 0.71 → 1.02 → 1.47 (k=6–10), ds_change_rel 0.06 → 0.74
  over the same range, IBL residual floor degrading 3e-5 → 4.7e-4
  (gv3_3_history_coarse.csv). With the artificial seed-sink masked, the
  residual instability is the REAL tail-cone feedback: meridian
  convergence amplifies any azimuthal δ* asymmetry into ṁ, which the FP
  solve returns as u_e distortion at the cone. Gain < 1 early, → 1
  late. This is the measured, case-level justification for V4
  (quasi-simultaneous) — the counterweight to GV3.2's 4–5-iteration
  convergence on the airfoil.
- Re per body length 3.0e6 (re_chord 7.4432e5).

## Medium level

Not executed. The pre-registration marks medium "RECORDED if
generated"; the coarse verdict is decisive (both FAILs are
state/physics-localized, not resolution artifacts the medium would
clear) and a 1M-tet loose loop on the same drifting tail would cost
> 1 h for a RECORDED line. The mesh exists
(cases/meshes/fuselage_bor/medium.msh, gitignored — regenerate with
generate_fuselage_bor.py) and `--levels medium` runs it on request.

## Conclusion

**GV3.3: honest FAIL (0 PASS / 2 FAIL / 7 RECORDED), exit 1 as
designed.** The smoke's engineering purpose is met: the closed-body
coupling runs end-to-end without a numerical event, mid-body
axisymmetry is demonstrated on a genuinely 3-D unstructured surface,
and the two FAILs + the loop record pin the remaining defect precisely
at the tail cone — loose coupling cannot close the stern of a closed
body, which is the strongest measured argument FOR V4 in the V4
skip/proceed decision (against GV3.2's airfoil-side 4–5 iteration
convergence). Per the pre-registration no fourth fix round was
attempted; the decision input goes to the user.
