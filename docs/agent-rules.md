# pyFP3D Agent Rules

Current phase: **B9 ✓ CLOSED 2026-07-17 (RE-SPEC'D, user-approved): wing-body
cross-model validation — BOTH wake models on the M2 wing-body (M∞0.5, α3.06°,
coarse+medium).** ★ HEADLINE: the level-set (Picard) and conforming (NEW
capability, P14 Newton) wing-body lifts AGREE to **cl_p 0.4% / cl_kj 0.6%** at
medium (conf 0.2173/0.2188 vs LS 0.2165/0.2175; coarse 12.8% = resolution) —
the wing-body analogue of P14's wing-alone cross-model. GB9.1/9.2/9.3/9.5 ✓,
GB9.6 RECORDED, **GB9.4 XFAIL** (fuselage lift 16-20%, resolution/model-sensitive
⇒ G1.6 fuselage-Cp error, band NOT moved). Demo `cases/demo/b9_wingbody/`
7 PASS + 1 XFAIL; guardrail `cases/analysis/b9_fuselage_guardrail/`.
- ★ **Conforming wing-body is the NEW capability** —
  `onera_m6_wingbody_mesh(embed_wake=True)`: fuselage built as TWO π-revolves
  (`add_fuselage_solid_split`, else the waterline-imprinted single revolve
  surface is unmeshable), through-body sheet, Netgen OFF; `cut_wake` /
  constraints / P14 pressure-Kutta ALL unchanged; embed_wake=False bit-identical.
- ★ **LS uses PICARD, not Newton — measured (user-questioned).** The committed
  LS Newton recipes (lagged-LU, `precond="schur"` B14, N5 freeze, B15 ramp) ALL
  diverge/churn on the subsonic wing-body: the FAILURE IS THE FAR-FIELD BC, not
  the solver. `neumann` is unbounded under the fuselage blockage (res 1e43);
  the `freestream` Newton path (never exercised — all committed LS Newton use
  neumann) is locally inconsistent (te_aux perfect 1.8e-8, 8 far-field fluid
  rows |R|≈84). LS Picard + `freestream` converges cleanly (res 3e-7). A
  `newton_ls.py` freestream Dirichlet bug was fixed in passing; the outer-row
  residual inconsistency is a recorded follow-up. **Next phase = user's call.**

**P14 results (evidence: [demo_report/track_p.md](demo_report/track_p.md) §P14).**
S1 and S2 both die in one estimator swap: M0.84 Γ(z) roughness 0.0970 →
**0.0043** (coarse) / 0.0365 → **0.0024** (medium) — at/below the level-set
band; all-station raw TE Cp gap 0.2206 → **0.0040** / 0.1585 → **0.0024**
(**55×/67×**), on G14.6's PRIMARY clause (raw recovery), fallback unused.
★ **Metric-baseline trap (I fell in it; erratum same-day):** A2 measured the TE
gap TWO ways — *section-last-point* (0.318/0.228 conf vs 0.009/0.002 LS = the
34×/133× headline) and the *all-station sweep* (0.2206/0.1585 conf). Quote the
one your pipeline actually ran. And `a2_te_gap.csv`'s LS rows are ≈0 because
they read the LS's OWN control volumes (its own constraint residual — A2's
"cannot be used as A/B"): to compare TE gaps across paths you must re-measure
the LS wall through the same section sweep (V14.6 does; LS = 0.0047 medium).
Wiring scope (user-arbitrated): coupled Newton + `solve_laplace_lifting` only —
`solve_subsonic_lifting`'s inner secant and `continuation.py` stay probe-based.
- ★ **The cross-MODEL result (V14.6, `cross_model_medium_m084.csv`) — the
  phase's strongest evidence.** The level-set path has ALWAYS used
  pressure-equality Kutta (B4). Medium M0.84: conforming-pressure
  **0.2776/0.2823** vs level-set **0.2772/0.2813** = agreement to
  **0.17%/0.36%**, where the probe path sat **4.5%/4.3% BELOW** LS — from a
  different wake model, DOF space, and mesh family (`onera_m6_wakefree`). The
  long-standing conforming-vs-LS lift disagreement WAS the Kutta form. Caveats:
  cross-model, NOT a same-mesh A/B; the LS state carries 1 lim/2 flr (B15
  caveat) vs 0/0; and "both agree" ≠ "both right" (a shared model error like
  the rigid planar wake is common to both by construction).
- ★ **G14.7 ✓ CLOSED — re-specced to the level-set oracle; the lift move is
  the finding (user-arbitrated 2026-07-17).** The gate opened against the G8.2
  **probe** locks; it XFAILed as written (band not moved after the fact,
  cl_KJ +4.85%), and the pre-registered mechanism fired exactly: the closures
  agree pointwise to the probe's own O(h) reading bias (cross-read 0.79% at
  medium M0.84 — a shifted closure, not a wandered solution), which the Kutta
  map's b ≈ 0.93 amplifies 1/(1−b) ≈ 14× into Γ, so the lift MUST move — and it
  moves ONTO the level-set answer (0.15%/0.34%, V14.6). **User verdict: accept
  the move, re-lock against the level-set oracle** (`< 1%`, PASS). Direction
  recorded: |cl_KJ − 0.288| 0.0188 → 0.0057, **69% of P9's "0.019 gap" was
  Kutta-estimator bias** — P9 could not see it (both its meshes shared the
  estimator, common mode to its Richardson). Closing G14.7 asserts the two
  paths AGREE, NOT that the M6 fine converges (it is not a discrete solution)
  NOR that "the 0.019 gap is resolution" (still *strongly indicated, NOT
  earned*).
- ★ **CORRECTION (measured 2026-07-17, V14.7) — the TE Cp SPIKE drops too, and
  A2's S2 decomposition needs a nuance.** P14's own earlier write-ups asserted
  the spike was "untouched, a wake-model-independent P1 recovery artifact" —
  that was a PREDICTION from A2's reasoning, never measured on the pressure
  path. Measured medium M0.84 raw: probe **0.1143** → pressure **0.0533**
  (2.1×), and now BELOW the level-set path's 0.0743. A2 was right that a
  spike is shared (~0.05 residual = the genuine recovery floor); it was wrong
  that the conforming EXCESS over LS was recovery — that part was Kutta-form
  error too (a wrong Kutta gives a genuinely wrong TE flow, and the
  common-mode spike metric cannot separate that from a recovery artifact).
  Tell: P6 smoothing no longer helps on the pressure path (0.0533→0.0660→
  0.0626 over 0/1/2 passes; A2 measured 0.147→0.081 on the probe path).
  **General lesson: do not carry a prior phase's attribution into a new
  measurement as if it were a result — measure it.**
- **Recorded honestly, not rounded away:** the fixed-Γ discriminator on the new
  estimator is **D = 1.80** (probe 7.33) — inside A2's INCONCLUSIVE zone
  (confirm > 3 / refute < 1.5), so the estimator is 4× cleaner but not a
  perfect measurement operator.
- **Recipe (measured):** seed the pressure Newton from the probe solution — the
  quadratic row's basin is smaller (M6 medium M0.5: cold Picard-5 seed wanders
  to cl +16% and fail-fasts at 29 steps/417 s; probe-seeded = 3 quadratic
  steps/26 s, faster than the probe path itself). The M0.84 ramp seeds level 0
  the same way. tip_taper + pressure raises NotImplementedError (the B8 blend
  is not re-derived).

B9 scope guards (user-arbitrated 2026-07-14; RE-SPEC'D 2026-07-17,
user-approved): **subsonic M∞ 0.5 ONLY** (M0.84 excluded — the round-cap
transonic fine is still a non-converged limit cycle, G13.3 transonic
NEGATIVE); α 3.06° (the committed M6 subsonic convention). Meshes = the M2
wing-body family (`cases/meshes/onera_m6_wingbody/`, wake-free, for the LS
leg) **plus the NEW wake-embedded conforming variant**
(`cases/meshes/onera_m6_wingbody_conforming/`) — coarse + medium only, both
gitignored, regenerate before running (~5 min wake-free; conforming TBD).
The 2026-07-16 body/far-field re-spec stands (5 root chords, wing centered,
2-diameter ellipsoid nose, graded skin, R_FAR = 25 MAC, ★ needs
`Mesh.OptimizeNetgen` — sliver lottery 0.31/4.80/2.63 without it; wing
untouched, TE nodes bit-identical 76/150 so B9's Kutta stations do not move).
**`wall_tag` stays `"wall"` on BOTH paths** — widening it to include
`fuselage` would mint spurious Kutta stations along the sheet–body waterline.
The old M2 open verification item (innermost TE node's wall-adjacent CV fan
touches fuselage wall faces; the upper/lower CVs must take only wing-side
elements — `multivalued.py::_build_te_control_volumes`, and the conforming
analogue `te_pressure.py::TEControlVolumes`) is now **gate GB9.3**, both
paths. Cross-model gate GB9.5 pre-registered at < 1% (medium, cl_p(wing) +
exposed-span cl_KJ, same-extractor discipline); GB9.4 fuselage-no-lift ≤ 5%
of wing cl_p at medium; GB9.6 = the kept 2026-07-14 fuselage-Cp guardrail
(RECORDED, no pass/fail).

## Track status (one line each; authority = docs/roadmap/*.md ledgers)

- **Track P** ([track_p.md](roadmap/track_p.md)): P0–P9 ✓ (P1: G1.6 strict
  xfail) · P10 ◐ (G10.1 open) · P11 conditional-not-opened · P13 ◐ (G13.3
  transonic NEGATIVE-open) · **P14 ✓ CLOSED 2026-07-17** (pressure-equality
  Kutta estimator, from A2 — S1 jitter + S2 TE Cp gap both gone, and the
  conforming path now matches level-set on lift; G14.1–G14.7 ✓).
- **Track M** ([track_m.md](roadmap/track_m.md)): M0–M5 ✓; M2 solver leg CLOSED
  by B9 2026-07-17 (both wake models run on the wing-body; the conforming leg
  added a wake-embedded variant `onera_m6_wingbody_conforming/`).
- **Track B** ([track_b.md](roadmap/track_b.md)): B1–B8, B11–B15 ✓ ·
  B6 ◐ (medium quantitative closed by GB15.4) · **B14 ✓ CLOSED 2026-07-17**
  (`precond="schur"` Schur-eliminated-aux + AMG(SPD Picard main block),
  `pyfp3d/solve/schur_ls.py`; the A1 precond bottleneck is GONE — M6 medium
  M0.84 42.6% → 2.6%, ramp 1.43× / subsonic 2.08×, γ = the committed GB15.4;
  ★ SLOWER at small scale, the fine memory-bounded route stays the unbuilt
  designed use-case) · B10 shelved · **B9 ✓ CLOSED 2026-07-17 (RE-SPEC'D)** —
  wing-body cross-model: LS (Picard) + conforming (NEW capability, Newton)
  agree 0.4%/0.6% at medium M0.5; GB9.4 fuselage-lift XFAIL ⇒ G1.6; LS Newton
  diverges on the wing-body = the neumann far-field blockage, not the solver.
- **Track V** ([track_v.md](roadmap/track_v.md)): designed, zero implementation.
- **Track A** ([track_a.md](roadmap/track_a.md)): created 2026-07-15 · **A1 ✓**
  (2026-07-16, GA1.1–GA1.5; 4-driver timing instrumentation + cost benchmark) ·
  **A2 ✓** (CLOSED 2026-07-17, GA2.1–GA2.5; TE/Kutta fidelity attribution —
  **S1 conforming Γ(z) jitter = a per-station probe-difference Kutta-target
  measurement artifact, NOT flow content** (fixed-Γ discriminator D=7.33/25.70
  coarse/medium; closure |F|/|Γ| ≤ 0.6%); **S2 TE Cp jump = potential-jump
  Kutta form error** (34×/133× vs level-set) **+ P1 recovery artifact** (both
  paths); fix routed to **P14** designed-not-started; no `pyfp3d/` edits;
  **B9 stays NEXT**).
  ★ **In 3-D both Newton paths are PRECONDITIONER-bound (~40% of
  wall, lagged LU already on); the "Picard seed is the cost" result is 2.5-D
  ONLY** — quote a dominant phase with its mesh. ★ **"G8.2 = 250 s" is a
  superseded number**: P10's promoted recipe made it ~145 s.

Human-readable snapshot + document map: [overview.md](overview.md). The pre-split
narrative history of this file:
[archive/agent-rules-2026-07-15.md](archive/agent-rules-2026-07-15.md) (archive,
not a spec; its GB15.3 timings are pre-CSV — trust the committed CSVs).

## Session disciplines (each one measured; citation in the named record)

1. **Thread cap 16 incl. BLAS/OMP** — `NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
   OPENBLAS_NUM_THREADS=16`; missing the BLAS caps costs ~33% on this 16C/32T
   shared box and fails the G8.2 300 s assert (P8 record). Check load before
   heavy runs; drop to 8 threads when another heavy job is in flight.
2. **Do not recompute expensive committed artifacts** (P4 heavy demo ~40 min,
   G4.1 medium ~17 min, P5 medium 45–75 min, conforming fine M0.5 ~34 min,
   G13.2 fine transonic ~45 min). The committed CSV/PNG is authoritative;
   re-run only when a real change moves those numbers AND you will commit the
   refresh. Verify routine edits on the cheap coarse path.
3. **Evidence needs a committed artifact** — a number that exists only in .md
   prose is not evidence (2026-07-13 audit; CLAUDE.md doc-map rule). If a run
   is too expensive to repeat, that is the reason to commit its CSV.
4. **Fold-zone discipline** (NACA medium M≈0.79, dcl/dM ≈ 6–10): single-mesh
   regression locks only, never grid-convergence claims; keep the Mach ramp
   (G10.3 verdict); loose intermediate tolerances CONTRAINDICATED near folds
   on the conforming path (G10.2b) — the LS ramp's freeze-armed loose levels
   are the designed exception (B15; the LS path has no Γ DOF).
5. **M6 fine (~450k dofs) conforming Newton: `precond="amg"` + tight EW
   forcing (η=1e-8) + `m_start=0.30, n_picard_seed=12`** — the medium recipe's
   `precond="direct"` is the 4h39m/26GB splu trap (P9; 24 GB RSS is the tell).
   The LS path has no fine escape yet (that is B14, not built); LS medium uses
   B12/B13 lagged-LU (`direct_refactor_every`).
6. **LS transonic recipes**: Picard — `tol_residual=1e-5` (above the shock
   plateau) or every level burns its full budget; strict convergence — the B15
   Newton ramp (`solve_multivalued_newton_transonic`), with `freeze_tol` ABOVE
   the Mach-rising churn floor and the `freeze_max_clamped` caveat quoted
   honestly (the converged M0.84 state carries 3 clamped cells; see
   design_track_b.md §14 + the roadmap B15 corrections block).
7. **Gated tests**: heavy transonic/Newton/M6 gates run only under
   `PYFP3D_TRANSONIC_GATES=1`; M6 `.msh` are gitignored (regenerate ~30 s via
   `cases/meshes/onera_m6/generate_onera_m6.py`).
8. **Do not re-open closed negatives** without new evidence: G1.6 fix routes
   (h-refinement / recovery tweaks / Nitsche / boundary-data corrections) are
   ruled out; B8's constraint-side tip cures are measured dead; spanwise-Γ
   smoothing is a dead route; "the 0.019 gap is resolution" stays *strongly
   indicated, NOT earned* (2026-07-14 wording arbitration); B15 did NOT revive
   G9.1.

Baseline: **442 passed + 20 skipped + 2 xfailed** (2026-07-17, B9 wing-body,
+13 passed / +1 skipped = `tests/test_b9_wingbody_{conforming,ls}.py`
(conforming 8 + 1 gated skip, LS 5), measured 1084.20 s @16 threads).
Previous: 429 + 19 + 2 (2026-07-17, B14 Schur+AMG, +8/+1 =
`tests/test_b14_schur_ls.py`; 1043.37 s @16 threads); 421 + 18 + 2 (P14 tier 1+2, +15 =
`tests/test_p14_te_pressure.py`; 1015.17 s @8 threads); 406 (= the 399 M2
number once A1's 7 tests landed), 973.59 s @8 threads, 2026-07-16; 396,
988.73 s @16 threads, 2026-07-15; lineage in [overview.md](overview.md). After
any kernel/assembly change run `pytest tests/test_v0_freestream.py` first
(CLAUDE.md hard rule 1).
