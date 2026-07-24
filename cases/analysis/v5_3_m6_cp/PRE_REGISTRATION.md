# PRE-REGISTRATION — GV5.3: M6 wing direction+magnitude check vs committed Cp

Committed BEFORE the first code change (Track V discipline; user-sequenced
2026-07-24: GV5.1d → GV5.5 → GV5.2–5.4). Binding text:
`docs/roadmap/track_v.md` GV5.3 (re-anchored 2026-07-22, user-directed: the
committed `cases/reference_data/onera_m6_experiment/` holds **Cp only** — no
experimental CL value is committed, so this gate does NOT use the external
"experiment ≈ 0.26–0.27" figure).

## 1. Question

On the ONERA M6 wing at the committed experimental condition (TEST 2308:
**M 0.8395, α 3.06**, Re_MAC 11.72e6), does the loose VII loop move the
solution in the physically expected direction by more than the committed
input-noise floor?

- **(a)** Does the viscous CL move **DOWN** from the converged same-mesh
  inviscid baseline by more than the A4 input floor (2.5 %, medium)?
- **(b)** Does the viscous wall Cp at the 7 committed span stations match
  the committed experiment **better** than the same-mesh inviscid baseline
  (Cp RMS-to-experiment decreases)?

This is a **direction+magnitude** gate on a 3-D transonic lifting wing —
the first live read of whether the committed IBL model adds physical
content beyond the inviscid solve at a condition with a shock. It is NOT
a shock-position or absolute-Cp validation gate (no such band is
registered).

## 2. Conditions, meshes, data

- ONERA M6, conforming path: `cases/meshes/onera_m6/{coarse,medium}.msh`
  (gitignored; regenerate with `cases/meshes/onera_m6/generate_onera_m6.py`).
  **Medium binding**, coarse recorded (the committed inviscid anchors and
  the A4 floor are medium-quoted; GV5.0 recommended the medium anchor).
- Condition = the dataset-labeled TEST 2308 verbatim: **M 0.8395,
  α 3.06**, Re_MAC = 11.72e6 (the bridge convention; the experiment title
  reads 0.117E+08), forced transition x_tr/c = 0.05 both sides (the
  committed GV5.0 convention; no experimental transition data committed).
- Experiment: `cases/reference_data/onera_m6_experiment/experiment-Cp.dat`
  — Tecplot, 7 zones, stations Y/b = 0.20 / 0.44 / 0.65 / 0.80 / 0.90 /
  0.96 / 0.99; per-point (X/L, Y/b, Z/L, Cp) with X/L = local x/c.
  **Ground truth — never modified.**
- Inviscid anchors (committed, M 0.84, P14 pressure-Kutta,
  `cases/demo/p14_pressure_kutta/results/m084_pressure.csv`): coarse
  cl_p 0.2628 / cl_KJ 0.2688; medium cl_p 0.2776 / cl_KJ 0.2823. The
  ΔM = 0.0005 between the dataset label (0.8395) and the anchor run
  (0.84) is absorbed by the wiring-guard tolerance (§6 W1).

## 3. Protocol (the committed recipes verbatim)

- Loose VII loop = the GV3.1 recipe verbatim (ω = 1.0, ≤ 10 outer,
  tol_ds = 1e-3) driving `viscous/coupling.py::build_wing_case` (the
  GV5.0 wing case: LE-band laminar pin per local x/c, both TE natural
  outflow, root symmetry natural, **tip band z > 0.95·b_semi pinned +
  ṁ-masked** = the production tip_taper radius, scope guard).
- FP driver = the **P14 transonic recipe verbatim** (no re-invention):
  the k = 0 cold start seeds M 0.70 with a probe-Kutta Newton solve
  (`solve_newton_lifting`, M6_NEWTON_KW) then ramps to M 0.8395 with
  `solve_newton_transonic` (**NEWTON_M6_RECIPE imported from
  `tests/test_p8_newton.py`**, pressure Kutta, n_picard_seed = 0); outer
  iterations k ≥ 1 warm-start `solve_newton_lifting`
  (kutta_estimator="pressure", n_picard_seed=0) from the previous state.
  No tip_taper on the FP side (P14 verbatim — the tip mask lives on the
  IBL side per the scope guard); this makes the k = 0 solve directly
  comparable to the committed P14 anchors (§6 W1).
- Metrics are read on the k = 0 (inviscid) and terminal (viscous) states
  of the SAME run — the "same-mesh inviscid baseline" of the binding
  text is the loop's own k = 0 state.

## 4. Metrics and bands

### (a) ΔCL direction gate (PASS/FAIL, medium binding, coarse recorded)

Δcl = cl(terminal) − cl(k = 0), quoted for BOTH committed estimators
(cl_p pressure integral, cl_KJ Γ integration; `post/surface.py`
idioms, S_ref = planform area of the wall tag):

- Δcl < 0 AND |Δcl| / cl(k=0) > 2.5 % (the A4 medium input floor) ⇒
  **PASS** (direction confirmed beyond input noise).
- |Δcl| / cl(k=0) ≤ 2.5 % ⇒ **RECORDED, flagged input-limited** (the
  pre-registered "smaller move" clause — not a FAIL).
- Δcl > 0 AND |Δcl| / cl(k=0) > 2.5 % ⇒ **FAIL** (wrong direction
  beyond the input floor).

Both estimators are quoted; the gate reads on cl_KJ (the Γ-integrated
estimator the binding text anchors, 0.2823), cl_p is the consistency
read (a disagreement between estimators larger than the A4 floor is
itself RECORDED).

### (b) Per-station Cp RMS direction gate (PASS/FAIL, medium binding)

Computed Cp: `post/section_cut.py::section_cp_curve` (the committed P5
idiom) at each η = Y/b station on the k = 0 and terminal phi, x/c
normalized by the cut loop's own extent, compressible Cp at the run
Mach. Experiment side assignment: Z/L sign (positive = upper/suction
side — verified at execution by the LE stagnation point and the
suction-side match, printed in the CSV; a flip would be caught by the
k = 0 wiring read). Computed Cp interpolated linearly onto each
experiment point's (x/c, side); per-station RMS over all zone points
(both sides pooled).

- Binding stations = the 5 **unmasked** stations η = 0.20 / 0.44 / 0.65 /
  0.80 / 0.90. **PASS** = RMS_viscous < RMS_inviscid at ≥ 4 of 5
  stations AND the pooled (all unmasked-zone points) RMS decreases.
  Otherwise **FAIL**.
- The 2 **tip-masked** stations η = 0.96 / 0.99 are RECORDED-only: the
  IBL is pinned + ṁ-masked there by construction, so viscous ≈ inviscid
  and |ΔRMS| ≈ 0 is expected; a large shift is a wiring anomaly to
  investigate, not a gate event.
- The A4 input band is quoted as the scale annotation: medium 2.5 %
  peak-relative u_e (LE band 4–7 %) ⇒ a Cp-scale uncertainty of order
  2·|Cp|·2.5 % ≈ 0.05·|Cp| at peak suction; a per-station |ΔRMS|
  smaller than 0.05 in Cp units is annotated **input-limited** (still
  counted in the direction tally, flagged).

### (c) Convergence and guards (RECORDED)

n_outer / converged per level, ds_change_rel, mdot_max history,
ds_neg_floored, IBL n_iter / final residual per outer, crossflow
max|B|/max|A| (the GV5.0 read), FP-solve stats per call (converged /
n_newton / rescue path taken / accept_reason), per-level wall time
(8-thread session constraint — flagged non-comparable).

## 5. Execution mechanics (pre-registered, no band/protocol change)

- **FP rescue chain** (the GV5.2 addenda #2/#3 semantics verbatim,
  pre-registered here up-front against the M ≥ 0.725 shock-cell
  plateaus): each FP call tries (i) the strict solve (target |R| ≤
  1e-10); on non-convergence (ii) a `solve_newton_transonic`
  continuation warm-start (m_start = 0.80 → target, NEWTON_M6_RECIPE);
  on a stall at a |R| ≲ 1e-9 plateau (iii) honesty-guarded acceptance
  (Kutta constraint converged, zero limiter/floor activity at
  acceptance, accept_reason logged per call; the per-attempt path is
  recorded in run.log and the counts in summary.csv). Strict stays the
  first choice at every call.
- **Hard raise** on a non-finite state or an FP failure the chain
  cannot absorb (GV3.3 discipline: loud, not silent-RECORDED).
- **mdot runaway** (the GV3.3 transpiration-runaway class: the
  loud-fail guard in `run_loose_coupling`) ⇒ the level's metrics are
  RECORDED at the last finite state as a **recipe limit** (the GV5.2 §6
  precedent), the run continues with the remaining levels, and the gate
  verdict reflects the remaining readable level(s) with the runaway
  level RECORDED. If the runaway hits the BINDING level, (a) and (b)
  read RECORDED-not-gated and the outcome is "gates unreadable —
  recipe limit RECORDED", NOT PASS.
- **Loose non-convergence** (not converged in ≤ 10 outer, no runaway):
  the level's (a)/(b) gates still read at the terminal state (the GV5.2
  precedent: the terminal state carries the loop's physical content),
  and the non-convergence is itself RECORDED under (c).
- **Wiring guards** (a fire = RuntimeError = recipe error, NOT a gate
  verdict): W1 = the k = 0 inviscid cl_p / cl_KJ within 1 % (relative)
  of the committed P14 anchors (coarse 0.2628 / 0.2688, medium 0.2776 /
  0.2823; the 1 % absorbs the ΔM = 0.0005 dataset-vs-anchor label and
  fixed-point scatter — a miss means the driver recipe drifted, stop);
  W2 = the experiment-zone side mapping sanity (each zone's max-Cp
  point sits at x/c < 0.05 = the LE stagnation point, and the pooled
  k = 0 RMS computed with the chosen mapping is below the RMS computed
  with the flipped mapping); W3 = section extraction guards from
  `section_cp_curve` (min_points_per_side) propagate as errors.
- Runner exits 1 iff any PASS/FAIL gate reads FAIL (honest-FAIL
  discipline; RECORDED-only outcomes exit 0).

## 6. Artifacts

`cases/analysis/v5_3_m6_cp/`: this pre-registration, run.py, VERDICT.md
(filled at execution), `results/summary.csv`,
`results/history_{coarse,medium}.csv` (the loose-loop per-outer table),
`results/cp_stations_{k0,final}_{level}.csv` (per-station computed Cp
curves), `results/cp_rms_{level}.csv` (per-station RMS vs experiment,
inviscid / viscous / Δ + input-limited flags),
`results/cp_overlay_{level}.png` (7-station Cp overlays: experiment vs
inviscid vs viscous), `results/cl_history_{level}.png`. run.log stays
in the working tree (`*.log` is gitignored — the GV5.2 precedent).

## 7. Priors and expected reads (context, not bands)

- GV5.0 (M 0.5, same wing case + loop): the loop did NOT converge ≤ 10
  outer at either level (coarse: tip-adjacent separation-patch runaway;
  medium: bounded δ* limit cycle 2–12 %/k); ΔCL DOWN both estimators
  (medium −2.4 % / −2.1 %, input-limited). At M 0.8395 the suction and
  δ* grow, so the viscous signal should be LARGER than at M 0.5 — but
  the loop's contraction at a transonic condition is untested (GV5.2:
  the 2.5-D loose loop failed to contract at M ≥ 0.725).
- The A4 band (medium 2.5 % peak-rel u_e, LE band 4–7 %) is the quoted
  input-noise floor for both gates.
- Cost estimate @8 threads: k = 0 ramp ≈ 10 s (coarse) / ≈ 5 min
  (medium, P14 measured 288 s @16 threads); each outer = one warm FP
  solve + one IBL solve; the GV5.0 medium outer cost was contention-
  polluted — quoted flagged. Expected total ≈ 30–90 min for both
  levels.
