# GV6.2 PRE-REGISTRATION — measured wake-IBL on/off effect vs the A4 band

**Committed before execution** (Track-V discipline). Date: 2026-07-25.
Branch `kimi/track-v6-gv6-2`.

**Gate text** (`docs/roadmap/track_v.md`, V6): *"**GV6.2 measured
effect**: wake-IBL on/off cl (and TE-region Cp) delta on the GV3.1 case,
direction-checked against XFOIL's wake modelling; RECORDED with the A4
input band quoted."*

**Inputs** (committed, never re-computed):

- GV6.0 ruling (user 2026-07-25,
  `cases/analysis/v6_0_design_adjudication/DESIGN_ADJUDICATION.md` §6):
  Option A + producer (i); producer (ii) — the solved wake IBL — opens
  **only if the GV6.2 measured effect is significant vs the A4 input
  band**. The significance adjudication is the **user's**; this gate
  RECORDS, it does not rule.
- GV6.1 ✓ CLOSED (6 PASS / 0 FAIL / 7 RECORDED,
  `cases/analysis/v6_1_wake_sheet/VERDICT.md`): the non-binding on/off
  smoke measured **Δ-cl +0.00015** (cl ≈ 0.2825, ≈ +0.05 %), TE-region
  (x/c > 0.9) max |ΔCp| **0.00250**, outer counts 3 vs 3, ṁ_wake max
  1.8e-2; δ*_TE 6.62e-3 c / θ_TE 4.48e-3 c; L_rel = 1.0 c pinned MODEL
  CHOICE. Caveat carried: the committed GV3.1 medium fixed point is not
  cross-run reproducible (the v5_tight_coupling scatter caveat + the
  GV6.1 numba cache-mode finding); only in-session A/B pairings are
  meaningful.
- A4 (`cases/analysis/a4_ue_error_band/VERDICT.md`): the medium
  smooth-wall u_e input band ≈ **2.5 % peak-relative** / 0.04·U∞
  max-norm (LE/stagnation band 4–7 %, not exercised here).

## 1. Scope

- **No FP-solver behavior change. One additive library plumbing change**
  (needed by band (c)): `CouplingConfig.wake_l_rel_chords: float =
  WAKE_L_REL_CHORDS` (default 1.0) threaded to the
  `wake_transpiration_source` call (`viscous/coupling.py:916-924`, which
  today always uses the library default). Default = bit-identical; the
  GV6.1 (a) flag-OFF bit-identity discipline is untouched. +1 plumbing
  test (§5). Nothing else in `pyfp3d/` changes (guard G4).
- The gate case = **the committed GV3.1 medium recipe verbatim**: M 0.5,
  α 2°, Re 3e6, x_tr 0.05 both surfaces, loose Picard leg,
  `CouplingConfig` defaults (n_outer_max 10, tol_ds 1e-3, ω = 1.0),
  medium mesh `cases/meshes/naca0012_2.5d/medium.msh` (986 wake nodes /
  493 stations / s_max 14.5 c).
- **All bands RECORDED** — the gate text's verdict type. The runner
  exits 0 unless a §4 guard fires (guards are recipe-error raisers, not
  bands; honest-anomaly readings are recorded, never tuned).

## 2. The XFOIL wake-reference sourcing question — **RULED 2026-07-25
## (user): Option A** (posed for user ruling before band (b) executes)

The gate text's direction check needs XFOIL's *wake* solution downstream
of the TE. The committed reference CSVs
(`cases/reference_data/naca0012_viscous_xfoil/`) carry **surface rows
only** — the XFOIL DUMP wake rows (8 fields: s, x, y, u_e/V∞, δ*, θ,
C_f, H) were discarded at generation time
(`generate_xfoil_reference.py:168-176`). Options:

- **Option A (RECOMMENDED) — analysis-local sourcing.** This gate's
  runner drives the pinned XFOIL binary itself (the gitignored build
  under `tools/xfoil/`, build recipe in the reference README) at the
  committed reference conditions verbatim (NACA 0012, 280 panels, M 0.5,
  Re 3e6, α 2°, xtr 0.05/0.05, Ncrit 9, ITER 200 — the
  `generate_xfoil_reference.py` batch script's xtr005 leg), parses the
  DUMP wake rows, and commits them as a gate artifact
  `results/xfoil_wake.csv` next to the VERDICT.
  `cases/reference_data/` is **untouched**; the polar reproduction guard
  G3 pins the run to the committed reference.
- **Option B — extend the committed reference.** Amend
  `generate_xfoil_reference.py` to also emit `wake_*.csv` into
  `cases/reference_data/naca0012_viscous_xfoil/` and regenerate (README
  records why). The proper long-term home, but it modifies the
  ground-truth directory (agent-rules hard rule 6 protects the existing
  CSVs; *adding* files is a ground-truth change needing the user's
  explicit ruling).
- **Option C — no wake rows.** Direction-check against the committed
  surface TE values + XFOIL's documented wake-model behaviour only. No
  empirical wake profile; arguably does not satisfy "direction-checked
  against XFOIL's wake modelling".

**Viability smoke (disclosed, NOT committed, no number enters any
band):** 2026-07-25, this session, scratch dir — the pinned binary ran
rc = 0 at the committed xtr005 conditions; the saved polar reproduces
the committed `polar_summary.csv` xtr005 row to every printed digit
(cl 0.2691 / cd 0.00926 / cm 0.0011); the DUMP carries 33 wake rows
spanning x/c ∈ [1.0001, 2.0]. Band (b)'s numbers are recorded fresh at
execution from the runner's own XFOIL invocation.

## 3. Bands (all RECORDED)

- **(a) On/off measured effect — RECORDED.** In-process A/B on the §1
  case: the OFF leg then the ON leg (L_rel = 1.0 c) in one process from
  the same initial state (identical JIT/cache state, threads, seed —
  the GV6.1 cache-mode finding cannot enter an in-process pairing).
  Recorded: cl_off / cl_on (final), Δ-cl (abs + rel), per-outer cl
  histories, outer counts, TE-region (x/c > 0.9, both surfaces) max
  |ΔCp| and its location, ṁ_wake max, δ*_TE / θ_TE (final), wall times
  (flagged non-comparable). **A4 quoting (pinned formulas)**: rel Δ-cl
  vs the A4 medium peak band 0.025; max |ΔCp| vs the first-order
  propagation δCp_A4 = 2·(u_e/U∞)_TE·0.025 with (u_e/U∞)_TE the OFF
  leg's TE-region max recovered u_e/U∞ (recorded; Cp = 1 − (u_e/U∞)²).
  Context anchors: GV3.1's inviscid-cl 2.6 % below XFOIL ≈ the
  A4-consistent cl floor; the GV6.1 smoke's +0.00015 / 0.00250 as the
  prior reading. The significance read (above/below the A4 input band)
  is recorded; the producer-(ii) opening decision stays with the user
  (GV6.0 ruling).
- **(b) XFOIL wake direction check — RECORDED** (executes only after
  the §2 ruling; the text below assumes Option A, B is equivalent with
  the reference path as the source). Comparisons, all recorded:
  (i) near-wake relaxation **direction**: the sign of dδ*/dx over
  XFOIL's wake rows vs the producer's monotone δ*_TE → θ_TE relaxation
  (mass-transpiration relaxation, H → 1 far-wake limit) — agreement or
  honest anomaly (recorded as found, never tuned);
  (ii) the residual fraction (δ*−θ)/(δ*_TE−θ_TE) at x/c = 2.0 from
  XFOIL's rows vs the producer's e^(−1) = 0.368 at L_rel = 1.0 c;
  (iii) the TE anchor: our δ*_TE/c (ON leg, final outer) vs XFOIL's
  first-wake-station δ*/c and vs the committed surface-TE sum
  δ*_upper(TE) + δ*_lower(TE) = 0.01190 c (xtr005) — quoted with the
  GV3.1 caveat (our wall δ* runs low vs XFOIL, the GV3.1 FAIL finding);
  (iv) downstream θ: XFOIL's wake θ evolution vs the producer's
  conserved-θ construction — model-form difference, recorded.
  Our δ*_wake(s) profile is the ON leg's final-outer producer state;
  the comparison window is XFOIL's wake extent (x/c ∈ [1, 2]); our
  strip extends to 14.5 c (recorded, no XFOIL counterpart).
- **(c) L_rel sensitivity sweep — RECORDED.** The ON leg re-run at
  L_rel ∈ {0.5, 2.0}·c (the (a) legs give 1.0 c; the OFF leg serves all
  three), same in-process pairing discipline: Δ-cl and TE-region max
  |ΔCp| per L_rel. **L_rel = 1.0 c stays the pinned MODEL CHOICE
  regardless of the readings** — the sweep records sensitivity, it is
  not a tuning instrument (the GV6.1 registration).

## 4. Guards (recipe-error raisers)

- **G1**: every A/B pair runs in-process, same initial state, threads,
  seed (the GV6.1 numba cache-mode finding: no cross-compile comparison
  is used anywhere in this gate).
- **G2**: flag-OFF = legacy bit-identical stands (GV6.1 (a) gated it;
  this gate adds no code path that can reach the legacy solve — the
  only library delta is the §1 additive config field with the default
  preserving today's call).
- **G3**: the runner's XFOIL invocation reproduces the committed
  `polar_summary.csv` xtr005 row (cl/cd/cm to the printed digits);
  mismatch = harness error (wrong binary/conditions), raise before any
  wake row is read.
- **G4**: the gate's `pyfp3d/` diff touches only `CouplingConfig` (one
  field), the one call site (`viscous/coupling.py`), and
  `tests/` — asserted by inspecting the committed diff before the
  VERDICT is written.

## 5. Tests (+1)

`tests/test_v6_wake_sheet.py`: the `wake_l_rel_chords` plumbing — a
non-default value reaches the producer (the produced δ*_wake(s) matches
the pinned formula `θ_TE + (δ*_TE − θ_TE)·exp(−s/(L·c))` at L = the
configured value), and the default preserves the GV6.1 behaviour. The
existing v6 suite (7/7 JIT) stands unchanged.

## 6. Execution discipline

- Thread cap **8** (NUMBA/OMP/OPENBLAS — the standing temporary
  user-directed session constraint, the GV6.1 §7 discipline), identical
  across every leg; wall times flagged non-comparable. If the user
  lifts the cap before execution, an addendum re-pins before the run.
- Committed artifacts are never re-computed; GV6.1 (d)'s numbers are
  quoted as the prior smoke, not re-used as this gate's measurement —
  band (a) is its own in-session pairing.
- `cases/reference_data/` untouched under Option A (Option B would
  require the explicit §2 ruling first).
- Addenda: any change to pinned bands/guards/recipe/constants requires
  an addendum section appended to this file, committed BEFORE the
  (re-)execution it affects.
- Runner exit 0 = guards clean (all bands RECORDED); VERDICT.md +
  results/ (summary.csv, histories, te_cp.csv, wake_profiles.csv,
  xfoil_wake.csv per the §2 ruling, gv6_2.png) committed.

## 7. Explicitly OUT of scope (registered follow-ups)

- Producer (ii) — the solved wake IBL (opens only on the user's GV6.2
  significance adjudication, GV6.0 ruling).
- The LS sheet-source leg; the Newton/tight driver wiring of b_wake
  (GV6.1 §8).
- Coarse cross-check; the xtr030 XFOIL variant; any operating point
  other than the committed GV3.1 recipe.
- Geometric wake relaxation (forbidden by the V6 design constraint).
- Any tuning of L_rel (band (c) records sensitivity only).
