# pyFP3D Agent Rules

Current phase: **B20 ✓ CLOSED 2026-07-18 (NEW, user-directed; executes B19 Leg
B): mixed-side plain elements can now take their density from the MAIN field,
and the Leg B hypothesis is answered — as a SPLIT.** New `plain_density` knob on
`MultivaluedOperator`: `"side"` (**default, bit-identical to every committed
result**) / `"main"`. ★ **The reporting layer had already made this call** —
`element_mach2` has defaulted to `mixed_plain="main"` since 2026-07-14 — so B20
makes the ASSEMBLY agree with the DIAGNOSTIC. ★★ **A workspace-aliasing bug was
caught by measuring an unexpected result, not explaining it:**
`PicardOperator.velocities` returns VIEWS into a shared buffer, so recomputing
the main gradient inside the density path overwrote the caller's side values in
place — quasi-2D "moved" 0.77, subsonic Γ TRIPLED, the Jacobian degraded: all
ONE bug (2940 elements clobbered vs the 129 in the mask), fixed with `.copy()`.
*Accepting the plausible "Γ tripled ⇒ Leg B has a big effect" story would have
recorded an aliasing bug as a physics finding.* **GB20.1 ✓** quasi-2D
bit-identical, M6 R moves on 164 of ~12k rows. **GB20.2 ✓** the Jacobian stays
EXACT under main (8.07e-09 / 6.29e-10) — Leg A ∘ Leg B compose. **GB20.3 ✓**
2.5-D subsonic Γ **+0.0000 %** (the class is 3-D only ⇒ every committed
quasi-2D lock untouched). **GB20.4 ✓** M6 coarse ramp→M0.84: side **m 0.7875
NOT converged** vs main **M0.84 CONVERGED**. ★★ **GB20.5 RECORDED — the
hypothesis SPLITS.** B18 medium wing-body @M0.5: side res 6.8e-5 / **82
clamped** / Mmax 3.920 (a CLAMPED non-converged number) vs main **res 1.1e-13 /
6 clamped** / Mmax 5.220 (a genuine converged solution). ⇒ the **CONVERGENCE**
pathology was largely the contamination (churn/clamps cured), but the junction
**POCKET is REAL** and B19's literal hypothesis is **REFUTED** — removing the
contamination UNCLAMPED it, revealing a genuine M≈5.2 spike at subsonic
freestream = the **G1.6/GB9.4 faceted-geometry** error. Main still cannot pass
M0.5. **Do not repeat "the pocket is mixed-plain contamination" — measured
false.** ★★ **ADOPTED PERMANENTLY — the knob is REMOVED (user-arbitrated
2026-07-18).** It existed only to make the A/B measurable; leaving an internal
inconsistency behind a switch would mean shipping a defect as an option. The
decisive argument needs no physics: the element's stiffness is contracted with
the MAIN field while its density came from a DIFFERENT one — **one equation,
two velocity fields** — in an element no wake jump even passes through.
**Accepted cost, now PAID and measured (2026-07-19 re-baseline):** suite
**465+22+2 unchanged** and gated 3-D LS **67/67 green** — *no test lock breaks*;
★ the 3-D numbers are not locked by tests at all, they live in demo evidence
(open process gap). Every moved number went B20's way — B7 M_max 1.453→**1.392**
and tip Γ −0.0003→**−0.0000**; B16 legacy limited **3690→11**, pin floored
**3→0** (those 3 were recorded as "B8/G1.6 class, not chased" — they were this
contamination), medium clamps 42/40→**0/0**; B18 wing-body residual
6.8e-5→**1.1e-13**; M6 coarse ramp 0.7875-not-converged→**M0.84 converged**.
B9's cross-model headline (LS 0.2165/0.2175) did not move one digit.
★ **ONE regression: M6 medium M0.84** — the ramp reaches **M0.6625 (2/5)**,
γ 0.088338→0.071909 (B15 17/20, B14 5/7; every failure traces to this one case).
★★ **The old number's validation was COMMON-MODE:** B15 compared M_max against
the LS *Picard* (2.4549) — both LS solvers read the same contaminated density,
so it never checked correctness; the conforming record is **1.995**. And the new
1.5822 is at M0.6625, not M0.84 — not like-for-like. ★★ **GB20.7 ANSWERED (2026-07-19): a REAL capability loss, not a recipe
mismatch.** `freeze_tol` swept 1e-3→1e-6: the ceiling moves only 0.6625→**0.6750**
(1e-5/1e-6, and a level starts converging on `tol` instead of escaping via
`assignment_cycle`) — the hypothesis was directionally right, quantitatively
insufficient. **Post-B20 LS transonic envelope (M6 wake-free): coarse M0.84
converged (was 0.7875 not-converged), medium ≈M0.675 (was M0.84).**
★ **Synthesis: the contamination was an unintended STABILISER** — same pattern
as GB20.5's wing-body (converges beautifully, climbs less far). The old M0.84
medium state carried M_max 2.45–2.49, ~25 % above the conforming reference
**1.995**, and was validated only against the equally-contaminated LS Picard.
**State it as a trade:** the old code went further into states we now distrust;
the new code stops earlier and what it produces is clean (0/0 clamps, |R|~1e-13).
⇒ **GB15.4's "reaches M0.84" is now a NEGATIVE and GB14.4's "== committed
GB15.4 physics" is superseded** — both need a re-spec against the new envelope
(G14.7 precedent), left OPEN for the user since it redefines a committed
capability claim. Bound: only `freeze_tol` was varied.
★ **Demo-cache trap:** heavy demos reuse gitignored `results/*.npz` — delete the
LS caches and verify zero `cached` lines, or a re-run is a no-op (it cost me one
false "B7 unchanged" result). The junction pocket is NOT fixed by any of this
(G1.6).

**B19 ✓ CLOSED 2026-07-18 (NEW, user-directed; executes the
A3/GA3.6 C1 finding as TWO deliberately separated legs): the level-set Newton
Jacobian is now EXACT in 3-D, and the residual's own asymmetry is measured and
routed.** ★★ **Leg A was TWO defects, not one.** (1) **DOF maps** — Terms 2/3
used the mass-conservation SCATTER map for rows AND columns, but columns must
follow `side_potentials`' per-node READ map (they coincide on cut elements —
`readvec` reproduces `dofs_upper`/`dofs_lower`, asserted — and diverge on
mixed-side plain elements, a 3-D-only class). (2) **Gradient factors, the same
duality one level down** — the residual is `ρ̃(grad of the READ field)·V·(grad
of the SCATTER field · B_a)`, so the ROW factor needs `grad_row` while the
COLUMN factor keeps the side gradient; the code used the side gradient for both.
★ **Fixing (1) alone left 1.4697e-02 — 8× better but STILL ε-independent ⇒
recorded PARTIAL rather than rounded into a pass, which is what forced (2)
out.** ★ A block isolation (`|FD23−J23| ≡ |FD23−J2|` exactly ⇒ Term 3 emits
nothing there) found (2) — **the first reading of the row classification pointed
at Term 3 and was WRONG**; acting on that inference would have put a new bug
into correct code. **GB19.1 ✓ targeted probe 1.145684e-01 → 1.333699e-08**
(control 6.33e-10 unchanged); ★★ **the ε discriminator FLIPPED** — pre-fix
1.532e-01 at every ε (spread 1.00 = missing term), post-fix 1.6e-09/2.1e-08/
2.2e-07 (spread 131.5, ~1/ε = FD roundoff). **GB19.2 ✓ max|ΔR| = 0.000e+00**
bit-identical, `git stash` A/B after EACH fix ⇒ no converged LS result moves.
★ **GB19.4 ✓ RECORDED NEGATIVE — NO convergence gain**: γ 0.07212068 identical
to 8 dp, M_max 1.134235 to 6 dp, same 40 steps, same plateau, **+3.6 % wall**.
The plateau is the **B15 selection-churn limit cycle** and an exact derivative
cannot fix a discontinuous selection ⇒ **B19 must NOT be credited with a
convergence improvement**; it buys correctness (a Newton, not a quasi-Newton).
**GB19.5 ✓** `tests/test_b19_jacobian_3d.py` (+3) closes the blind spot;
★★ **ERRATUM — the blind spot was mis-stated** by C1 and by my own first test:
quasi-2-D has **129** mixed-side plain elements, not zero; it has **0** that
READ an aux (0 of 129 touch a cut node) — that is the real invariant.
★★ **GB19.6 (Leg B) — the residual asymmetry is NOT benign.** Only **252**
elements (**0.19 %** of volume) read an aux, but there max|ρ_side − ρ_main| =
**0.4474 (45.3 %)** and **the SIDE field reads q² 3.2229 (M≈1.80, at the M_cap
limiter) where the MAIN field reads 1.3379** — a **spurious supersonic state**
the artificial-density switch then acts on, so the contamination reaches the
solver's density, not just a diagnostic. Third bite from one element class
(B8 metric ×5 → Jacobian → residual). **HYPOTHESIS not result:** candidate
contributor to B18's refinement-worsening wing-body pocket (M_max 3.96);
wing-alone measurement, no causal link; named test recorded. **NOT adopted** —
changing the density source changes R ⇒ its own phase. Evidence
`cases/analysis/c1_ls_jacobian_fd/`.

**A3 ✓ CLOSED 2026-07-18 (NEW, user-directed): response to the
2026-07-17 independent Kimi inspection — docs consistency (17 findings) +
cross-path hardening + the C1 Jacobian verification.** ★★ **HEADLINE: C1 is
REAL and now verified — the LS Newton Jacobian is NOT the derivative of its
residual on mixed-side plain elements** (ONERA M6 coarse M0.70: targeted probe
‖Jv−FD‖/‖FD‖ = **1.146e-01** vs control **6.33e-10**, eight orders apart; and
**eps-INDEPENDENT** — 1.532e-01 at eps 1e-6/1e-7/1e-8, max/min **1.00** across
three decades ⇒ a missing term, NOT FD noise). The class is 3378 elements, not
just the 428 beyond-tip ones. **Bounded consequence: R is untouched, so every
converged LS state, γ, cl and gate number STANDS — what degrades is the
convergence rate (the LS Newton is a quasi-Newton in 3-D).** B6's FD gate could
not see it (quasi-2D has no such elements); B7/B15's M6 gates are convergence
gates, not FD gates. **RECORDED, NOT FIXED** — the side-aware column-mapping fix
is a shipped-kernel change that would move committed step-count trajectories ⇒
its own phase, user's call. Evidence `cases/analysis/c1_ls_jacobian_fd/`.
★ **Backported two B15-era LS fixes the review found un-backported for three
phases** (C2 selection-epoch fail-fast scoping, C3 `freeze_max_reverts` disarm)
to `solve/newton.py` — both dormant on every committed run. ★ **Reader hardened
(C4/C5):** unnamed physical surface groups were SILENTLY dropped, whose chain
ends at Γ(root) pinned to 0 with no error on an imported mesh; +3 tests,
verified failing before the fix. ★ **T1: the B3 gate is now actually locked** —
the test asserted 2% against a 0.3% gate; measured **0.1441%**, tightened.
★ **NO `pyfp3d/` numerics change** — every edit is dormant-path, default-inert,
or a measured no-op (C7b: 0 affected nodes on all six committed families).
★ Close-out ritual extended to **five surfaces + a backport check**
(CLAUDE.md step 5, disciplines #9/#10) — the audit's 17 findings were mostly
close-out debt. Suite **463 + 21 + 2** (+3 reader). Response report:
[inspection/20260718-response-to-kimi-inspection.md](inspection/20260718-response-to-kimi-inspection.md).

**B18 ✓ CLOSED 2026-07-18 (NEW, user-directed; appended after
B17; executes the GB16.6 debt): wing-body transonic (M0.84) — conforming reaches
it, level-set is junction-limited.** ★★ The wing-body transonic capability is
**asymmetric, and that is the finding.** **Conforming** (Newton + pressure Kutta,
Mach continuation) IS the wing-body transonic path: coarse reaches **M0.84 (cl_p
0.2617)**, medium reaches **M0.79 strict (cl_p 0.2579)** with a clean transonic
rise cl_p(M) = **0.2173/0.2321/0.2579** at M0.50/0.65/0.79 (medium M0.80+ stalls —
NOT slivers, the medium mesh is clean; a sharper shock/junction interaction,
recorded not chased; ★ the conforming wing-body medium ramp needs `freeze_tol`
raised to the wing-body churn floor 1e-6→1e-5, the B17 lesson). **Level-set** (B15
freeze-ramp + B17 pin_gamma) does NOT reach transonic on the wing-body: the
wing-fuselage junction spurious supersonic pocket (**G1.6/GB9.4/B8 mixed-plain**,
M²≈1.27 already at M0.5) **WORSENS with refinement** — coarse ceiling ~M0.575
(Mmax 1.44), medium dies at the FIRST transonic level ~M0.5 (Mmax artifact 3.96,
nlim 43/nflr 40); the direct analogue of GB9.4's fuselage-lift-grows-with-
refinement, a closed-negative discretization error (discipline #8), characterized
not chased. ⇒ **no common transonic Mach at medium** (LS can't leave 0.5), so the
trustworthy cross-model stays **M0.5 (2.6%)**; a coarse M0.60 transonic
cross-model was skipped (LS coarse ceiling 0.575 < 0.60). GB18.1 PASS + GB18.2–5
RECORDED. ★ **repays the GB16.6 evidence debt** (spec'd RECORDED but never
implemented; B18 executes it as a negative). ★ **NO `pyfp3d/` numerics change** —
pure demo/tests/docs on existing `solve_newton_transonic` +
`solve_multivalued_newton_transonic`. fine excluded (G13.3). Tests
`tests/test_b18_wingbody_transonic.py` (4, ungated); demo
`cases/demo/b18_wingbody_transonic/` (7 gates: 1 PASS + 6 RECORDED).

**B17 ✓ CLOSED 2026-07-18 (resolves GB16.4): the far-field aux pin must carry
jump=γ, not 0.** ★★
HEADLINE: **GB16.4 was NOT a non-convergence — it was a BC-modelling error in
B16's freestream pin.** B16 pinned the outflow wake jump to **0**, which REMOVES
the circulation the wake physically carries out (medium cl_p 0.2165→0.1690, a
−22% resolution-dependent error; the coarse "match to conforming" was a
coincidence — jump=0 there cancelled the coarse legacy's outer-tet garbage).
**Decisive discriminator:** giving the **Picard** driver the same freestream pin
(new `farfield_aux` knob on `solve_multivalued_lifting`) makes medium Picard-pin
converge cleanly (res 7.5e-8) to cl_p **0.1691** — matching the "stalled"
Newton-pin **0.1690** to 0.1% ⇒ two independent solvers on the same value ⇒ a
genuine BC-determined state, NOT a Newton stall. Fix = `farfield_aux="pin_gamma"`
(aux = host φ∞ − side·γ, jump→γ, refreshed with the live γ — the new default on
**both** solvers): the triangle closes MONOTONE to conforming (cl_p wing) —
coarse conf 0.2089 / legacy 0.1853 / pin0 0.2086 / **pin_gamma 0.2087**; medium
conf 0.2173 / legacy 0.2165 / pin0 0.1690 / **pin_gamma 0.2117 (Picard) = 0.2115
(Newton)**, both solvers agreeing 0.1%. GB17.1–17.4 ✓, GB17.5/17.6 RECORDED;
demo `cases/demo/b17_farfield_pin_gamma/` (3 coarse PASS + gated medium), tests
`tests/test_b17_farfield_pin_gamma.py` (6).
- ★ **B16 conflated two orthogonal issues:** the far-field near-singular
  **conditioning** (the pin cures it, jump value irrelevant — cond1 O(1e19)→8.7e6
  either way) and the outflow **circulation** (needs jump=γ). A third, pre-existing
  issue — the wing-fuselage-junction churn (medium Newton-pin_gamma still carries
  nlim 42/nflr 40, res 5.5e-5, the **G1.6/GB9.4** class) — survives but limits only
  the residual floor, not the lift (γ stable 0.06420, cl_p 0.2115 correct).
- ★ **Post-processing is NOT the cause (user's suspicion checked, GB17.2):** cl_p
  (surface-pressure integral) and cl_KJ (circulation integral) move together
  (~22% both) ⇒ a real flow-state change. The "section Cp looks aligned yet cl_p
  differs 22%" is a Cp-axis scale illusion; the plotted spanwise sectional cl(z)
  is the **Γ-based** `2Γ/(u·c)`, and per-station ∫Cp differs 24–44% while Cp
  curves differ only ~0.03–0.05 on a ±1 axis.
- ★ **Defaults (user-arbitrated 2026-07-18):** `pin_gamma` is the new default on
  BOTH `solve_multivalued_newton` (was `"pin"`) and `solve_multivalued_lifting`
  (was free/legacy). It acts ONLY on `farfield="freestream"`, inert (bit-identical
  to legacy) on vortex/neumann ⇒ every committed 2.5D vortex/neumann Picard run +
  every neumann Newton anchor byte-untouched. B9/B16 **freestream** Picard demos
  pinned to explicit `farfield_aux="legacy"`; B16 jump=0 reproduces with explicit
  `"pin"`. **B9 erratum:** coarse 12.8% was far-field contamination, not resolution
  (its medium legacy≈conforming headline stands). `"pin"` (jump=0) kept as the
  diagnostic value.
- ★ **vortex evaluation (GB17.6, user-requested):** `farfield="vortex"` does NOT
  close the 2.6% residual — it BRACKETS conforming from the other side (medium
  +2.5%) and its free far-field aux churn at coarse (res 3.2, needs its own pin).
  freestream pin_gamma stays recommended.

**B16 ✓ CLOSED 2026-07-18 (churn fix; GB16.4 RESOLVED BY B17 above):** the
wing-body LS-Newton churn is a **near-singular far-field aux block** — a wake
level set has no outflow clip, so the sheet reaches the far field and the outer
nodes it crosses carry aux DOFs governed only by near-singular wake-LS rows on
giant outer tets; at the freestream Picard state they hold garbage (coarse
|jump| **53.4** at x≥10 vs Γ̄ 0.0586), which Picard absorbs but the Newton
residual reads as the **8 far-field MAIN rows max|R| = 84.457** (aux-block cond1
**O(1e19) → 8.70e6** pinned). The pin fixes the CONDITIONING (coarse freestream
Newton res 5.88e-14, 0 limited vs legacy 7.95/3690) — but B16's freestream jump=0
value was wrong for the lift; see B17 above. GB16.1/16.2/16.3-coarse/16.5 ✓,
GB16.6 RECORDED; demo `cases/demo/b16_farfield_aux/` (9 PASS + 1 XFAIL, the XFAIL
now resolved by B17).
- ★ **The proposal's mechanism was WRONG** (self-corrected): it blamed the
  Picard weld (`closure="continuity"`) vs Newton `wake_ls`, but the lifting
  Picard uses `wake_ls` too — the difference is fixed-point ABSORPTION of the
  near-singular rows vs Newton's residual reading them, NOT the closure.
- ★ **GB16.3 coarse honest limit:** pin carries `n_flr=3` at the wing-fuselage
  junction = the **B8 mixed-plain / G1.6 fuselage-Cp** class (M²_side 7.32 vs
  M²_main 0.29, same root as GB9.4's xfail), orthogonal to the BC fix, not chased.
- ★ **Pinning the 4 far-field-boundary aux also cured the 4 INTERIOR junk aux**
  (their wake-LS rows now anchor to clean Dirichlet data) — the R2 risk void.

**B9 ✓ CLOSED 2026-07-17 (RE-SPEC'D, user-approved): wing-body cross-model
validation** — the level-set (Picard) and conforming (NEW capability, P14
Newton) wing-body lifts AGREE to **cl_p 0.4% / cl_kj 0.6%** at medium (conf
0.2173/0.2188 vs LS 0.2165/0.2175; coarse 12.8% = resolution). GB9.1/9.2/9.3/9.5
✓, GB9.6 RECORDED, **GB9.4 XFAIL** (fuselage lift 16-20% ⇒ G1.6 fuselage-Cp
error, band NOT moved). Demo `cases/demo/b9_wingbody/` 7 PASS + 1 XFAIL. The B9
LS-Newton follow-up ("neumann res 1e43; freestream Newton 8 rows |R|≈84") is now
**closed by B16** — the freestream Newton path works with the aux pin.
- ★ **Conforming wing-body is the NEW capability** —
  `onera_m6_wingbody_mesh(embed_wake=True)`: fuselage built as TWO π-revolves
  (`add_fuselage_solid_split`, else the waterline-imprinted single revolve
  surface is unmeshable), through-body sheet, Netgen OFF; `cut_wake` /
  constraints / P14 pressure-Kutta ALL unchanged; embed_wake=False bit-identical.

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
gitignored, regenerate before running (~4–5 min wake-free; conforming TBD).
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
  M0.84 43.6% → 2.6% (same-session A/B; A1's earlier independent read was
  42.6%), ramp 1.43× / subsonic 2.08×, γ = the committed GB15.4;
  ★ SLOWER at small scale, the fine memory-bounded route stays the unbuilt
  designed use-case) · B10 shelved · **B9 ✓ CLOSED 2026-07-17 (RE-SPEC'D)** —
  wing-body cross-model: LS (Picard) + conforming (NEW capability, Newton)
  agree 0.4%/0.6% at medium M0.5; GB9.4 fuselage-lift XFAIL ⇒ G1.6; LS Newton
  diverges on the wing-body = the neumann far-field blockage, not the solver ·
  **B16 ✓ CLOSED 2026-07-18 (churn fix; GB16.4 resolved by B17)** (NEW, executes
  the B9 follow-up) — the wing-body LS-Newton churn is a near-singular far-field
  aux block (cond1 O(1e19), 8 rows |R|≈84 reproduced); the pin drops it to 8.7e6
  ⇒ COARSE freestream Newton res 5.88e-14 / 0 limited where legacy churns at
  7.95; neumann byte-identical · **B17 ✓ CLOSED 2026-07-18 (NEW, resolves
  GB16.4)** — GB16.4 was a BC-modelling error, not a non-convergence: B16's
  freestream pin forced the outflow wake jump to 0, removing the physical
  circulation (medium −22%); an independent Picard-pin converges to the same
  0.169 the Newton-pin "stalls" at (both solvers agree per-BC).
  `farfield_aux="pin_gamma"` (jump→γ, new default on both solvers) closes the
  triangle MONOTONE to conforming (coarse 0.2087, medium 0.2117 Picard / 0.2115
  Newton); acts only on freestream, inert on vortex/neumann; vortex brackets from
  +2.5%; B9 coarse-12.8% erratum'd as far-field contamination · **B18 ✓ CLOSED
  2026-07-18 (NEW, executes the GB16.6 debt)** — wing-body transonic M0.84:
  CONFORMING reaches it (coarse M0.84 0.2617, medium M0.79 0.2579, clean cl(M)
  rise 0.2173/0.2321/0.2579), LEVEL-SET junction-limited (coarse ~M0.575, medium
  dies ~M0.5 — the G1.6/GB9.4 junction pocket WORSENS with refinement, Mmax
  1.4→4.0, closed-negative); no common transonic Mach at medium ⇒ cross-model
  stays M0.5 (2.6%); GB18.1 PASS + GB18.2–5 RECORDED; no `pyfp3d/` change.
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
9. **Two paths ⇒ backport check** (A3, 2026-07-18, measured): a fix landing on
   the conforming path or the level-set path must be explicitly checked against
   the other, with the answer recorded in the phase entry. Two B15 LS
   robustness fixes (selection-epoch fail-fast scoping, `freeze_max_reverts`
   disarm) sat un-backported to `newton.py` for three phases until an external
   review found them. Also: a test's bound must not be looser than the gate
   text it claims to lock — B3's test asserted 2% against a 0.3% gate
   (7×) until A3 measured 0.1441% and tightened it.
10. **Close-out refresh is FIVE surfaces, not two** — track ledger,
   agent-rules current-phase + baseline, overview.md, **PROJECT_STRUCTURE.md
   (footer AND the directory trees)**, and the `cases/*/README.md` row. The
   2026-07-17 audit found 17 consistency defects, most of them exactly the
   last two surfaces. Full checklist in CLAUDE.md workflow step 5.

Baseline: **465 passed + 22 skipped + 2 xfailed** (2026-07-18, B19 LS-Newton
Jacobian exactness, +2 passed / +1 skipped = `tests/test_b19_jacobian_3d.py`
(2 structural locks + 1 gated 3-D FD gate); **measured 1101.50 s @16 threads**).
Previous: 463 + 21 + 2 (2026-07-18, A3 inspection response, +3 passed =
`tests/test_mesh_reader_roundtrip.py`'s unnamed-physical-group locks);
460 + 21 + 2 (2026-07-18, B18 wing-body
transonic, +4 passed = `tests/test_b18_wingbody_transonic.py` (4 ungated on the
committed 2.5D NACA mesh; the wing-body transonic ramps live in the gated demo));
456 + 21 + 2 (2026-07-18, B17 far-field pin_gamma, +6 passed =
`tests/test_b17_farfield_pin_gamma.py`, 1097.11 s @16 threads);
450 + 21 + 2 (2026-07-17, B16 far-field aux pin, +8 passed / +1 skipped
= `tests/test_b16_farfield_aux.py`); 442 + 20 + 2 (2026-07-17, B9 wing-body, +13/+1 =
`tests/test_b9_wingbody_{conforming,ls}.py`, 1084.20 s @16 threads);
429 + 19 + 2 (2026-07-17, B14 Schur+AMG, +8/+1 =
`tests/test_b14_schur_ls.py`; 1043.37 s @16 threads); 421 + 18 + 2 (P14 tier 1+2, +15 =
`tests/test_p14_te_pressure.py`; 1015.17 s @8 threads); 406 (= the 399 M2
number once A1's 7 tests landed), 973.59 s @8 threads, 2026-07-16; 396,
988.73 s @16 threads, 2026-07-15; lineage in [overview.md](overview.md). After
any kernel/assembly change run `pytest tests/test_v0_freestream.py` first
(CLAUDE.md hard rule 1).
