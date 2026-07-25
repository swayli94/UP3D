# GV6.2 VERDICT — measured wake-IBL on/off effect vs the A4 band

- Date: 2026-07-25 · Branch: `kimi/track-v6-gv6-2`
- Binding text: `docs/roadmap/track_v.md` GV6.2
- Pre-registration: [PRE_REGISTRATION.md](PRE_REGISTRATION.md)
  (committed before execution, `ab89483`; the XFOIL wake-reference
  sourcing **RULED Option A 2026-07-25 (user)** — analysis-local
  sourcing, `cases/reference_data/` untouched)
- Evidence: `results/summary.csv`, `results/history_{off,on_l05,
  on_l10,on_l20}.csv`, `results/te_cp.csv`, `results/wake_profiles.csv`,
  `results/xfoil_wake.csv`, `results/gv6_2.png` — regenerate with
  `python cases/analysis/v6_2_measured_effect/run.py` (exit 0 = guards
  clean; all bands RECORDED per the gate text)
- **VERDICT: GV6.2 EXECUTED — 0 PASS / 0 FAIL / 24 RECORDED.**
  Guards G1–G4 clean (G3 XFOIL polar reproduction verified in-line;
  G4 the `pyfp3d/` diff is exactly the registered plumbing: 1 file,
  +6 lines). Test lane: `tests/test_v6_wake_sheet.py` 8/8 (7 GV6.1 +
  the new GV6.2 plumbing test, JIT). **The measured on/off effect sits
  at 0.022× (cl) / 0.051× (Cp) of the A4 input band — NOT significant
  vs the A4 input band; the producer-(ii) opening decision is the
  user's (GV6.0 ruling), this gate RECORDS.**
- **Post-verdict adjudication 2026-07-25 (user)**: the recorded
  reading accepted — the significance condition is NOT met →
  **producer (ii) NOT opened; V6 ✓ CLOSED** (the GV6.0 ruling clause,
  `docs/roadmap/track_v.md` V6).

## Result table

| gate | metric (band) | measured | verdict |
|---|---|---|---|
| (a) | cl flag OFF / ON (final, L_rel = 1.0 c) | 0.28246 / 0.28261 | RECORDED |
| (a) | on/off Δ-cl (abs + rel; in-process A/B) | +0.00015 (+0.0547 %) | RECORDED |
| (a) | A4 quoting: \|rel Δ-cl\| vs 0.025 | 0.0547 % vs 2.50 % (ratio 0.022) | RECORDED |
| (a) | outer count ON vs OFF | 3 vs 3 | RECORDED |
| (a) | ṁ_wake max over outers | 1.800e-02 | RECORDED |
| (a) | δ*_TE / θ_TE (ON final) | 6.61795e-03 / 4.47577e-03 | RECORDED |
| (a) | TE-region (x/c > 0.9) max \|ΔCp\| + location | 0.00250 at x/c 1.0000 upper | RECORDED |
| (a) | A4 quoting: max\|ΔCp\| vs δCp_A4 = 2·(u_e/U∞)_TE·0.025 | 0.00250 vs 0.0493 ((u_e/U∞)_TE = 0.9861; ratio 0.051) | RECORDED |
| (a) | wall per leg (8-thread, NON-COMPARABLE) | 69 s / 69 s | RECORDED |
| (c) | L_rel = 0.5 c: Δ-cl / TE max \|ΔCp\| | +0.00017 (+0.0615 %) / 0.00173 | RECORDED |
| (c) | L_rel = 2.0 c: Δ-cl / TE max \|ΔCp\| | +0.00013 (+0.0464 %) / 0.00363 | RECORDED |
| (c) | L_rel = 1.0 c stays the pinned MODEL CHOICE | unchanged | RECORDED |
| (b) | XFOIL wake rows sourced (Option A) | 33 rows, x/c ∈ [1.0001, 1.9999] | RECORDED |
| (b) | (i) XFOIL wake δ* direction | 32/32 negative steps (monotone decreasing) | RECORDED |
| (b) | (ii) residual fraction at x/c = 2: XFOIL vs e^(−1) | 0.1105 vs 0.3650 | RECORDED |
| (b) | (ii′) XFOIL effective relaxation length (derived) | 0.454 c vs pinned 1.0 c | RECORDED |
| (b) | (iii) TE anchor: ours vs XFOIL first-wake vs committed surface-TE sum | 0.00662 vs 0.01440 vs 0.01191 | RECORDED |
| (b) | (iv) downstream θ: XFOIL θ(2)/θ(TE) vs conserved-θ | 0.6948 vs 1.0 (by construction) | RECORDED |
| (b) | XFOIL wake H at TE / at x/c = 2 | 2.0143 / 1.0760 | RECORDED |

## Per-band analysis

### (a) the measured effect is two decades of margin below the A4 input band

The in-process A/B on the committed GV3.1 medium recipe (OFF then ON,
identical initial state/threads/JIT state — the GV6.1 cache-mode
finding cannot enter an in-process pairing): **Δ-cl = +0.00015
(+0.0547 %)** on cl ≈ 0.2825, both legs converged in 3 outers; the
TE-region max |ΔCp| = **0.00250** at x/c = 1.0 upper. Quoted against
the pinned A4 propagation: 0.0547 % vs the 2.5 % medium peak u_e band
(ratio **0.022**); 0.00250 vs δCp_A4 = 2·0.9861·0.025 = 0.0493 (ratio
**0.051**). The wake-sheet correction moves the converged state by
~1/46 (cl) resp. ~1/20 (Cp) of the inviscid-input error floor — the
effect is real (deterministic in-session A/B, reproduced digit-for-digit
against the GV6.1 (d) smoke's +0.00015 / 0.00250 on the same 8-thread
env) but **not significant vs the A4 input band**. This is the
registered input to the user's GV6.0 producer-(ii) adjudication ("the
solved wake IBL opens only if the GV6.2 measured effect is significant
vs the A4 input band"); the decision itself is NOT this gate's.

### (b) XFOIL direction check: direction agrees, the rate and the TE
### anchor do not (model-form, recorded)

The runner drove the pinned XFOIL binary at the committed xtr005
conditions (the committed generator's batch script verbatim); G3
reproduced the committed polar to the printed digits (cl 0.2691 /
cd 0.00926 / cm 0.0011). The 33 wake rows (x/c ∈ [1.0001, 1.9999],
`xfoil_wake.csv`):

- **(i) direction AGREEMENT**: XFOIL's wake δ* decreases monotonically
  (32/32 negative steps) — the same relaxation direction as the
  producer's δ*_TE → θ_TE construction (H → 1: XFOIL 2.0143 → 1.0760).
- **(ii) rate DISAGREEMENT (recorded, not tuned)**: XFOIL's residual
  fraction at x/c = 2 is 0.1105 vs the producer's e^(−1) = 0.3650 —
  XFOIL relaxes ≈ 2.2× faster, an effective relaxation length
  **0.454 c** vs the pinned 1.0 c. Of the (c) sweep the 0.5 c leg is
  the closest to XFOIL's rate (`gv6_2.png` middle panel) — recorded;
  L_rel = 1.0 c stays the pinned MODEL CHOICE (the sweep records
  sensitivity, it is not a tuning instrument).
- **(iii) TE anchor low (the GV3.1 caveat)**: our δ*_TE = 0.00662 c vs
  XFOIL's first-wake 0.01440 c and the committed surface-TE sum
  0.01191 c — our wall δ* runs at ≈ 0.46–0.56 of XFOIL at the TE,
  consistent with GV3.1's registered δ* H-family offset in APG (the
  GV3.1 FAIL finding). The wake sheet source therefore acts on a
  smaller defect than XFOIL's — one more reason the (a) effect is
  small, quoted with the caveat, not chased.
- **(iv) downstream θ**: XFOIL's wake θ falls to 0.695 of its TE value
  at x/c = 2 (wake dissipation continues) vs the producer's
  conserved-θ construction — a model-form difference, recorded.

### (c) the insignificance is L-robust

Over the 4× sweep L_rel ∈ {0.5, 1.0, 2.0} c the on/off Δ-cl reads
+0.00017 / +0.00015 / +0.00013 (+0.046…+0.062 %) and the TE max |ΔCp|
0.00173 / 0.00250 / 0.00363 — the measured effect varies by ≈ ±15 %
around the 1.0 c value and stays 20–60× below the A4 band at every
L_rel. The (a) reading does not depend on the relaxation-length
choice within this family.

## Guards

- **G1** (in-process A/B pairing): every on/off pair ran in one
  process, same initial state, threads, JIT/cache state — no
  cross-compile comparison anywhere in the gate (the GV6.1 numba
  cache-mode finding).
- **G2** (flag-OFF legacy bit-identity): untouched — the only library
  delta is additive with a default that preserves the GV6.1 call;
  the plumbing test asserts explicit 1.0 ≡ default bit-identical
  (phi/gamma) on the coarse loop.
- **G3** (XFOIL polar reproduction): verified in-line before any wake
  row was read (cl 0.2691 / cd 0.00926 / cm 0.0011 == the committed
  `polar_summary.csv` xtr005 row to the printed digits).
- **G4** (diff scope): `git diff --stat -- pyfp3d/ tests/` =
  `pyfp3d/viscous/coupling.py` +6 (the import, the
  `wake_l_rel_chords` field, the one call-site argument) and
  `tests/test_v6_wake_sheet.py` +42 (the plumbing test + docstring) —
  exactly the registered scope, nothing else.

## Numerical settings (as run)

- Thread cap **8** (NUMBA/OMP/OPENBLAS; the standing temporary
  user-directed session constraint, PRE_REGISTRATION §6 — wall times
  non-comparable).
- Mesh: committed `cases/meshes/naca0012_2.5d/medium.msh` (986 wake
  nodes / 493 stations / s_max 14.5 c); the committed GV3.1 recipe
  verbatim (M 0.5 / α 2° / Re 3e6 / xtr 0.05 both surfaces,
  CouplingConfig defaults: n_outer_max 10, tol_ds 1e-3, ω = 1.0, loose
  Picard leg); all four legs converged at k = 3.
- Producer constants: L_rel ∈ {0.5, 1.0, 2.0} c swept, **1.0 c the
  pinned MODEL CHOICE, unchanged**; W2 TE-continuity asserted inside
  every producer call (never fired).
- XFOIL 6.99 (the gitignored pinned build under `tools/xfoil/`),
  batch per the committed generator (NACA 0012, 280 panels, M 0.5,
  Re 3e6, α 2°, xtr 0.05/0.05, Ncrit 9, ITER 200; converged in ~5
  Newton iterations).
- Wall per loose leg 69 s (8-thread session; NON-COMPARABLE).
