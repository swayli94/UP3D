# pyFP3D Agent Rules

Current phase: **B29 ✓ CLOSED 2026-07-20 (NEWEST; same branch as B28; no
`pyfp3d/` change): flat-fragment adopted as the wing-body LS PRODUCTION
config (user-adjudicated, B28 VERDICT §6).** B18 demo LS side C = B25 clip +
B28 flat sheet (`sheet_direction=(1,0,0)`; NEW `ls_flat_*` caches, the B26
tilted `ls_C_*` stale; side A tilted kept as the historical pocket
comparison; conforming legs bit-reproduce). M0.5 LS anchors re-pinned
0.2087/0.2117 → **0.2115/0.2184**. Flat C side: coarse 0.84 reached (cl_p
0.2551); **medium ceiling 0.7625 → 0.775** (dies 0.7875, cls a+dm; live
dying peak M3.98 @ wing TIP z=1.20 — GB18.4's C side now measured live);
cross-model gaps **2.6→0.5 % (M0.5), 2.4→1.1 % (M0.65 PASS ≤5 %), 2.5→1.1 %
(M0.75)**. **GB18.5 live flat decomposition: cl_fus 0.0382 (band −0.0006 /
out 0.0388 / poles 0.0007) @0.7875 vs conf 0.0423 @0.79** — the B26 tilted
"×2 out-band" reading (P11 watch item) retired per B28 (position
sensitivity, not a lesion). checks.csv **8/8 PASS**;
`tests/test_b9_wingbody_ls.py` switched to the production wiring (flat+clip,
5/5). Evidence: `cases/demo/b18_wingbody_transonic/results/` (B29-refresh
checks.csv / cross_model.csv / cl_vs_mach.csv / PNGs).

**B27 ✓ CLOSED 2026-07-20 (stacked on the unpushed
B25/B26 chain; B18 demo refresh + Track B doc close-out; no `pyfp3d/` change):
the pocket-healed level-set reaches the SAME ceiling site as conforming — the
"LS junction-limited (closed-negative)" B18 story is RETIRED.** ★ **GB27.1
PASS** — the conforming legs BIT-reproduce the committed B18 anchors (cl_p
0.2617/0.2173/0.2321/0.2579, cross 0.2178, M_max 2.15): B21/B22 were inert on
the conforming path. ★ **GB27.2 PASS** — the LS A/C ceiling legs BIT-reproduce
the committed B26 anchors (C coarse reached 0.84; C medium m_last 0.7625 dies
0.775 (b); A medium m_last 0.50 dies 0.5125 (a); A coarse m_last 0.82):
336/336 rows bit-identical in
`cases/analysis/b27_b18_demo_refresh/results/g27_consistency.csv`. ★ **GB27.3 —
cross-model UPGRADED**: M0.5 (2.6%, B9/B17) + **M0.65 medium 2.4% PASS (≤5%)**
+ M0.75 medium 2.5% RECORDED — all gaps in the ~2.5% B17 cl_p/cl_kj convention
band; new conf medium 0.75 point (cl_p 0.2483, monotone cl(M)). ★
GB27.4/27.5 RECORDED — GB18.4 re-answered (pocket = B23 inboard free-edge
singularity, healed by B25's `inboard_clip`; the residual limiter = wing-tip
P13 + high-M Newton, same class as conforming's) and GB18.5 refreshed (C-side
cl_fus 0.0781, out-band 0.0565 ≈ ×2 the A side → P11/curved-wall input); the
A-side-vs-B18-anchor divergence is the **B21/B22 freeze-capture effect**
(erratum carried in the demo docstring). Demo 8/8 PASS on a ~1 h 39 min full
re-solve (no caches). ★ the pre-B27 `b18_sections_conf_medium.png` was
silently EMPTY (`section_cp_curve` tuple→dict API drift, swallowed by
`except: pass`) — fixed in the refresh. Evidence:
`cases/demo/b18_wingbody_transonic/results/` +
`cases/analysis/b27_b18_demo_refresh/` (PRE_REGISTRATION.md / VERDICT.md /
g27_consistency.csv); Track-B narrative re-based in
[design_track_b.md](design_track_b.md) §22 (+ §18 erratum pointer). **Next
phase = user's call** — recommended candidate: **(b)-class ceiling
attribution** (LS+clip medium 0.775 vs conforming medium 0.80+: same
mechanism? if yes the spend moves to the Newton robustness both paths share);
the LS wing-body full-envelope evaluation is now OPEN (conforming demoted to
cross-validation; Track-V sheet-topology prerequisites all in place); cl_fus
out-band ×2 → P11 watch item.

**B26 ✓ CLOSED 2026-07-20 (B26-A; the pocket-healed LS transonic ceiling
re-measurement; no `pyfp3d/` change): the junction pocket WAS the LS wing-body
transonic ceiling limiter.** Same code, same mesh, same frozen B18 recipe —
the only variable is `inboard_clip` (B25's conforming-topology inboard sheet):
medium ceiling **0.50 → 0.7625**, coarse **0.82 → 0.84 (reached)**; the C
side's death class flips from (a) pocket-rejection to **(b) high-M Newton
stall** (dying peak at the wing TIP z≈1.20 = the P13 class; junction corridor
corrM ≤1.10 clean). ★ **T1 independent finding**: the A-side re-run diverges
from the committed B18 anchors (died 0.50/0.55) — that is the **B21/B22
freeze-capture fix**, not physics drift (A medium 0.50 now CONVERGES with the
same 3/3 clamps and the same Mmax 5.22; the pocket's true kill line on the A
medium side is 0.55, Mmax 13.1 > freeze_max_clamped=8). ★ watch item: C-side
cl_fus ≈ 0.076–0.078 ≈ ×2 the A side, the excess in the out-band component
(0.057–0.068) → P11/curved-wall input. Evidence:
`cases/analysis/b26_ls_transonic_ceiling/` (g1_summary.csv / g1_levels.csv /
g1_peaks.csv / g1_ceiling.png + VERDICT.md).

**P11 ✓ CLOSED 2026-07-19 (user-directed, opened +
closed same day; the sphere leg of the curved-wall-element route): the DP1
curved-element route is a measured NEGATIVE, and G1.6 is RE-ATTRIBUTED.**
★★ **G11.1 NOT MET**: a *verified* curved wall-adjacent layer
(`pyfp3d/solve/curved_wall.py` — tet10 geometry via `closest_point_normal`
midpoint projection, mapped-P1 field, ΔA delta assembly; opt-in
`stiffness_delta` on `solve_laplace`, default bit-identical; planar null test
ΔA ≡ 0 bitwise, quadrature = P1 reference to 1.3e-15) moves the medium sphere
only **11.56% → 11.33%** = the G1.4 boundary-data oracle ceiling. The
pre-registered risk FIRED and is the mechanism: mapped-P1 on quadratic
geometry loses linear reproduction at **O(h)** (measured 0.138 max coarse) —
the same order as the O(h) facet-normal error it removes. New dead route: do
not re-propose superparametric (mapped-P1) curved wall elements. ★★ **G11.2
negative + premise REFUTED — the re-attribution headline**: the committed
sweep script replicates the P1-era φ order collapse exactly (flat
**0.88/0.56/0.42**; curved 0.80/0.50/0.39 — curving restores nothing), and
the controls show the collapse was never geometric: a structured icosphere
shell with the SAME flat facets converges at order **1.67/1.98** (2.14% max
Cp at h≈0.036), and at fixed h_min=0.03 refining ONLY the far mesh drops
wall-φ error **3.17×** (argmax r=1.53 → wall; order restored to **1.89**) —
the collapse is the **fixed-bulk-mesh pollution floor** of a single-variable
h_min sweep. ⇒ **the medium mesh's 11.6% ≈ the intrinsic P1-field max-norm
capability at h=0.08** (structured control ≈11% there); geometric-crime share
≈0.2 pp; recovery ≈0.2–0.5 pp (unchanged). The 2%-max-at-medium bar demands
O(h²) wall velocity at h=0.08 — beyond ANY P1-field method on ANY mesh.
**G1.6's strict xfail STAYS**; PROJECT_STRUCTURE "Known gaps" + the G1.6 row
carry the erratum. ★ **Route fork recorded (= user's call)**: (a) Option C
re-spec with a measured PASSING form (all-scales-refined order ≥1.8 + mean-Cp
< 1% at h_min 0.03 — measured 1.89–1.98 / 0.60%); (b) isoparametric P2 wall
layer (the only route to the literal criterion); (c) accept as a permanent
recorded limitation. ★ **Wing-body caveat**: the "three wounds, one G1.6 root
cause" table (next_phase_priorities) lost its sphere anchor — GB9.4/GB20.5
worsen-with-refinement is qualitatively different behaviour and needs its own
discriminator before "G1.6 class" is quoted again. ★ Backport check: N/A —
opt-in Laplace path only, `newton.py`/`newton_ls.py` untouched. Evidence:
demo `cases/demo/p11_curved_walls/` (14 PASS + 2 XFAIL, ~4 min warm; sweep
meshes gitignored, ~8 min first run), tests `tests/test_p11_curved_walls.py`
(8, ungated), sweep/control CSVs committed. **Next phase = user's call**
(the G1.6 route fork above, or the standing LS-fine / Track-V ordering from
the priorities analysis — whose sphere-anchor caveat is now annotated in
place).

**B22 ✓ CLOSED 2026-07-19 (executes B21's recorded
follow-up + the Kimi-inspection N3/§2/§5 items; no `pyfp3d/` change): the
B21-state evidence is refreshed and the 3-D LS numbers are finally
test-locked.** ★ **GB22.1** B15 demo **20/20** (was 17/20 under B20; caches
deleted, zero `cached` lines): medium γ 0.088343 / M_max 2.4818 / |R|
9.048e-14 / 0 lim/1 flr / 6/6 levels / **511 s = 4.51×** vs the committed
Picard — bit-consistent with `n1_freeze_fix_sweep.csv`. ★ **GB22.2** B14 demo
**7/7** (was 5/7): medium lagged 505 s vs schur **345 s = 1.47×**, precond →
**1.8 %**, γ 0.088343 both arms; anchors re-pinned; ★ the M6 COARSE ramp also
moved under B21 (γ 0.0848 → **0.084931**, M_max 1.3684) — the freeze-capture
patch touches every 3-D freeze-armed ramp; ⚠ B14's demo_report index row had
been missing since close-out (found + added). ★★ **GB22.3 — N3 CLOSED**:
`tests/test_b22_ls_3d_anchors.py` (+2 gated) RE-SOLVES the committed M6
coarse (~35 s) + medium (~9 min) ramps and asserts m_final/γ/M_max/clamps
absolutely (γ rtol 1e-4 = 20× the measured run-to-run spread, four orders
below a B20-sized move) — the alarm that did not exist when the suite stayed
green through two re-baselines in two days. ★ **GB22.4** the re-baseline
erratum checklist is now process (CLAUDE.md workflow step 5 ★ clause +
discipline #11). ★ **GB22.5 RECORDED** — next-phase priority analysis
(`docs/analysis/next_phase_priorities_2026-07-19.md`): **P11 first** (G1.6
owns three refinement-worsening wounds; the wing-body line has no numerical
suspects left after B20/B21), LS fine second (B14-ready, but no blocked
question needs it), Track V after P11 (V1 2.5-D ladder parallelizable), plus
a cheap M_max same-family cross-check (2.4818-vs-1.995 is currently a
cross-FAMILY comparison — do not quote it as a defect). **Next phase =
user's call.**

**B21 ✓ CLOSED 2026-07-19 (NEW; executes the 2026-07-19 Kimi
second-round inspection's N1 finding + the D1–D10 doc-errata wave):
`freeze_side_state` was the ONE consumer of the side q²/ρ path the B20 patch
missed** — it captured the frozen (upstream, branch) selection on the
UNPATCHED side field while `newton_side_data` runs `_apply_main_density`
first, so on 3-D meshes every armed freeze locked a selection the live system
would not make (probe: **83 upstream + 9 branch** differences at M6 coarse
seeded M0.70, all aux-touching mixed-plain;
`docs/inspection/20260719-n1-freeze-probe.py`). ★★ **The one-line fix RESTORES
the M6-medium M0.84 ramp**: the committed recipe reaches **M0.84** with γ
**0.088343** (pre-B20 0.088338), M_max 2.4818, res 9.0e-14, **0 lim/1 flr**,
**515 s** (pre-B20 657 s / 3 clamped) — faster AND cleaner, and freeze_tol
1e-3/1e-5 agree to 5e-7 (`n1_freeze_fix_sweep.csv`). ⇒ **GB20.7's "real
capability loss" verdict is OVERTURNED** — the loss was B20's own patch gap,
not the fix's intrinsic cost — and the "contamination was an unintended
stabiliser" synthesis is RETIRED. Key mechanism fact: pre-B20 the capture and
the live sweep were CONSISTENTLY unpatched; **B20's partial patch CREATED the
inconsistency**; GB20.7's freeze_tol sweep could not see it in principle (it
varies WHEN the freeze arms, N1 is WHAT it captures — the Kimi §8.3 point).
**GB15.4's capability clause STANDS** with a small numeric re-baseline
(B15/B14 demo refresh = recorded follow-up; post-B21 LS envelope: M6 coarse
AND medium reach M0.84). ★ Test lock
`test_b15_ls_newton_freeze.py::test_freeze_capture_matches_live_density_3d`
(gated, premise-asserted, **verified FAILING pre-fix**) — the 2.5-D bitwise
lock is structurally blind here (0 aux-touching mixed-plain in quasi-2-D).
★ **Fourth bite from the mixed-plain class** (B8 diagnostic → B19 Jacobian →
GB19.6 residual → B21 freeze capture) ⇒ design_track_b.md §21 standing rule:
**every consumer of the side-field q²/ρ must explicitly decide and record its
density source.** ★ Same commit wave: the Kimi D1–D10 doc errata (B18
junction-pocket attribution corrected to G1.6 across all faces; B15 success
narrative bracketed by the two-reversal erratum trail; ledger/status
contradictions; cond1 **9.1e18** not 6.36e18; the 42/40 re-baseline row
relabeled B17-not-B16; B17 medium Newton **0.2114**; B16 post-B20 number
shifts disclosed; demo_report B19/B20/B21 rows + design_track_b **§19–§21**
so the numerics spec describes the SHIPPED discretization; M2 ledger ✓;
track_a A3 header + P2-backlog pointer; PROJECT_STRUCTURE tree +
`cases/analysis/`; N4 stale comments). The N3 process gap (3-D LS numbers
unlocked) was closed the same day by B22's gated anchor locks (above).

**B20 ✓ CLOSED 2026-07-18 (NEW, user-directed; executes B19 Leg
B): mixed-side plain elements can now take their density from the MAIN field,
and the Leg B hypothesis is answered — as a SPLIT.** A TEMPORARY `plain_density`
knob (`"side"` = the old behaviour / `"main"`) was built solely to make the A/B
below measurable, and was **REMOVED on adoption** (★★ block below) — `main` is
now hard-coded. ★ **The reporting layer had already made this call** —
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
contamination); B17 pin_gamma medium clamps 42/40→**0/0** (⚠ erratum
2026-07-19: earlier drafts mislabeled this row "B16 medium pin" — B16's own
medium pin was already 0/0 pre-B20) and its Newton cl_p 0.2115→**0.2114**;
B18 wing-body residual
6.8e-5→**1.1e-13**; M6 coarse ramp 0.7875-not-converged→**M0.84 converged**.
B9's cross-model headline (LS 0.2165/0.2175) did not move one digit.
★ **ONE regression: M6 medium M0.84** — the ramp reaches **M0.6625 (2/5)**,
γ 0.088338→0.071909 (B15 17/20, B14 5/7; every failure traces to this one case).
★★ **The old number's validation was COMMON-MODE:** B15 compared M_max against
the LS *Picard* (2.4549) — both LS solvers read the same contaminated density,
so it never checked correctness; the conforming record is **1.995**. And the new
1.5822 is at M0.6625, not M0.84 — not like-for-like. ★★ **GB20.7 ANSWERED (2026-07-19) — and OVERTURNED BY B21 the same day (see
the current-phase block above).** The sweep itself stands: `freeze_tol`
1e-3→1e-6 moved the ceiling only 0.6625→**0.6750** ⇒ freeze_tol was a
contributor, not the cause. But its "REAL capability loss" verdict, the
"contamination was an unintended STABILISER" synthesis, and the "GB15.4 is now
a NEGATIVE / GB14.4 superseded ⇒ re-spec = user's call" consequence are all
**RETIRED**: the actual mechanism was N1 (the un-patched freeze capture),
which the freeze_tol axis cannot see in principle, and fixing it restores
M0.84 outright. What survives: the old M_max 2.45–2.49 was validated only
COMMON-MODE against the equally-contaminated LS Picard (conforming record
**1.995**, different mesh family — still an open cross-model question, B21's
2.4818 does not close it).
★ **Demo-cache trap:** heavy demos reuse cached `results/*.npz` (committed or
not — since 2026-07-19 analysis-chain npz are committable) — delete the
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
it, level-set is junction-limited.** ★★★ **ERRATUM 2026-07-20 (B23–B27): the
"junction-limited (closed-negative)" story is RETIRED — the pocket = the B23
inboard free-edge singularity, healed by B25's `inboard_clip`; the post-cure LS
ceiling is co-located with conforming (B26: medium 0.7625 / coarse 0.84
reached) and the demo is refreshed to 8/8 PASS (B27, 336/336 bit-identical;
cross-model M0.65 2.4% PASS / M0.75 2.5%). The conforming legs below stand —
bit-reproduced by B27.** ★★ The wing-body transonic capability is
**asymmetric, and that is the finding.** **Conforming** (Newton + pressure Kutta,
Mach continuation) IS the wing-body transonic path: coarse reaches **M0.84 (cl_p
0.2590)**, medium reaches **M0.79 strict (cl_p 0.2545)** with a clean transonic
rise cl_p(M) = **0.2143/0.2290/0.2545** at M0.50/0.65/0.79 — **B32 (2026-07-22):
the production recipe carries tip_taper (vanish_smooth 0.05·b_semi) and the climb
REACHES M0.84 at medium (cl_p 0.2738, 0 clamps)** — the old "M0.80+ stall" was the
wing-tip sheet-termination singularity (B31 attribution), cured by the taper at
≈ −1.3 % cl_p cost (★ the conforming wing-body medium ramp needs `freeze_tol`
raised to the wing-body churn floor 1e-6→1e-5, the B17 lesson). **Level-set** (B15
freeze-ramp + B17 pin_gamma) does NOT reach transonic on the wing-body: the
wing-fuselage junction spurious supersonic pocket (M²≈1.27 already at M0.5)
**WORSENS with refinement** — at close-out coarse ceiling ~M0.575 (Mmax 1.44),
medium dies at the FIRST transonic level ~M0.5 (Mmax artifact 3.96,
nlim 43/nflr 40); the direct analogue of GB9.4's fuselage-lift-grows-with-
refinement, a closed-negative discretization error (discipline #8), characterized
not chased. ★★ **Erratum by B20/GB20.5 (2026-07-19): the "B8 mixed-plain"
attribution is measured FALSE** — the pocket is the **G1.6 faceted-geometry**
error (removing the contamination converged the medium case res 6.8e-5→1.1e-13
and unclamped a GENUINE Mmax 5.22; the 3.96 was a clamped artifact); post-B20
coarse ceiling **~M0.55 (Mmax 1.31)**. ⇒ **no common transonic Mach at medium**
(LS can't leave 0.5), so the trustworthy cross-model stays **M0.5 (2.6%)**; the
coarse M0.60 cross-model, originally skipped, EXISTS in the re-baselined
artifact (conf 0.2178 vs LS 0.2174 = 0.2%). GB18.1 PASS + GB18.2–5
RECORDED. ★ **repays the GB16.6 evidence debt** (spec'd RECORDED but never
implemented; B18 executes it as a negative). ★ **NO `pyfp3d/` numerics change** —
pure demo/tests/docs on existing `solve_newton_transonic` +
`solve_multivalued_newton_transonic`. fine excluded (G13.3). Tests
`tests/test_b18_wingbody_transonic.py` (4, ungated); demo
`cases/demo/b18_wingbody_transonic/` (7 gates at the 2026-07-18 close-out:
1 PASS + 6 RECORDED — **superseded 2026-07-20 by the B27 refresh: checks.csv
8/8 PASS**, see the erratum above).

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
conf 0.2173 / legacy 0.2165 / pin0 0.1690 / **pin_gamma 0.2117 (Picard) = 0.2114
(Newton; 0.2115 pre-B20)**, both solvers agreeing 0.1%. GB17.1–17.4 ✓,
GB17.5/17.6 RECORDED;
demo `cases/demo/b17_farfield_pin_gamma/` (3 coarse PASS + gated medium), tests
`tests/test_b17_farfield_pin_gamma.py` (6).
- ★ **B16 conflated two orthogonal issues:** the far-field near-singular
  **conditioning** (the pin cures it, jump value irrelevant — cond1 O(1e19)→8.7e6
  either way) and the outflow **circulation** (needs jump=γ). A third issue at
  close-out — the wing-fuselage-junction churn (medium Newton-pin_gamma carried
  nlim 42/nflr 40, res 5.5e-5) — limited only the residual floor, not the lift
  (γ stable 0.06420). ★ Post-B20 erratum (2026-07-19): that churn WAS the
  mixed-plain contamination — the re-baselined trajectory converges to
  |R| ~1e-13, γ 0.064201, cl_p 0.2114, clamps gone (GB20.5).
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
error, band NOT moved) — **corrected by B28 (2026-07-20)**: cl_fus = physical
carryover + wake-sheet POSITION sensitivity, NOT the G1.6 error; GB9.4
re-spec'd per B23 §(c) to out-band cross-model ≤ 15% and PASSES at 7.0%
(medium), demo now 8/8. Demo `cases/demo/b9_wingbody/` 7 PASS + 1 XFAIL → **8 PASS post-B28**. The B9
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
  the same way. tip_taper + pressure runs the B31 Gamma-pin row blend
  (solve/newton.py NewtonWorkspace docstring; the pin carries the row's own
  frozen slope SIGN — measured diag D > 0 on the conforming meshes, so an
  unsigned weld would amplify instead of unload).

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
  xfail, root cause RE-ATTRIBUTED by P11) · P10 ◐ (G10.1 open) · **P11 ✓
  CLOSED 2026-07-19** (curved wall elements measured NEGATIVE; G1.6 =
  intrinsic P1 capability at h=0.08, not the wall crime; route fork = user's
  call — see the current-phase block) · P13 ◐ (G13.3
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
  agree 0.4%/0.6% at medium M0.5; GB9.4 fuselage-lift XFAIL ⇒ G1.6 (**B28
  2026-07-20 corrects**: sheet-position sensitivity not an error; gate
  re-spec'd out-band cross-model ≤15%, medium 7.0% PASS, demo 8/8); LS Newton
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
  rise 0.2173/0.2321/0.2579), LEVEL-SET junction-limited (post-B20: coarse
  ~M0.55/Mmax 1.31, medium dies ~M0.5 with a GENUINE unclamped Mmax 5.22 —
  the junction pocket WORSENS with refinement; GB20.5 corrected the
  attribution to the G1.6/GB9.4 faceted geometry, NOT mixed-plain,
  closed-negative); no common transonic Mach at medium ⇒ cross-model stays
  M0.5 (2.6%) + a post-B20 coarse M0.6 point (0.2%); GB18.1 PASS + GB18.2–5
  RECORDED; no `pyfp3d/` change — ★★ **ERRATUM 2026-07-20: the
  "junction-limited (closed-negative)" story is RETIRED (pocket = B23 free-edge
  singularity, healed by B25's `inboard_clip`; LS ceiling now co-located with
  conforming) — see B26/B27** · **B19 ✓ CLOSED 2026-07-18** — LS-Newton
  Jacobian EXACT in 3-D (probe 1.146e-01→1.33e-08, ε-discriminator flipped; R
  bit-identical, NO convergence gain; Leg B routed the mixed-plain density
  contamination to B20) · **B20 ✓ CLOSED 2026-07-18, ADOPTED PERMANENTLY
  (knob removed), re-baselined 2026-07-19** — mixed-plain main-field density;
  the apparent M6-medium regression was B20's own patch gap (resolved by B21)
  · **B21 ✓ CLOSED 2026-07-19** — N1 freeze-capture alignment RESTORES the
  M6-medium M0.84 ramp (γ 0.088343, res 9e-14, 515 s); GB20.7 overturned;
  3-D capture lock added (verified failing pre-fix) · **B22 ✓ CLOSED
  2026-07-19** — B15 demo 20/20 + B14 7/7 refreshed on the B21 state (coarse
  ramp γ 0.084931 disclosed); **N3 closed** via gated absolute anchor locks
  (`test_b22_ls_3d_anchors.py`); re-baseline erratum checklist = process;
  next-phase analysis recommends P11 (user's call) · **B23 ✓ CLOSED
  2026-07-19** — junction discriminator: the pocket is lift/wake-coupled
  (α=0 clean at both levels, grows superlinearly with α), attribution = the
  wake inboard FREE-EDGE singularity (NOT G1.6 faceting); P11 close-out
  input · **B24 ✓ CLOSED 2026-07-19 (negative)** — the pocket follows the
  free edge, but the waterline-extension variants trade the singularity for
  equal-or-worse forms (B1 corrM 78.56 non-converged); the (b)-1 route is
  CLOSED · **B25 ✓ CLOSED 2026-07-19 (the cure)** — `inboard_clip` moves
  the sheet's inboard boundary to the fuselage surface / symmetry plane (=
  conforming fragment topology): medium α=3.06 junction pocket corrM
  **14.66→0.63**, guardrails clean; default None bit-identical · **B26 ✓ CLOSED
  2026-07-20 (B26-A)** — the junction pocket WAS the LS wing-body transonic
  ceiling limiter: healed by B25's `inboard_clip`; medium 0.50→0.7625, coarse
  0.82→0.84 reached; death class flips (a)→(b) (wing-tip P13; corridor
  clean); T1: the A-side anchor divergence = the B21/B22 freeze-capture
  effect · **B27 ✓ CLOSED 2026-07-20** — B18 demo refresh: LS+clip reaches
  the conforming ceiling site (GB27.1/27.2 bit-consistency 336/336;
  cross-model upgraded M0.5 + M0.65 2.4% PASS + M0.75 2.5%); the
  "junction-limited" story RETIRED; next = (b)-class ceiling attribution
  (user's call).
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
   last two surfaces. Full checklist in CLAUDE.md workflow step 5. Format
   note: the track ledgers and the roadmap.md / overview.md status blocks are
   wrapped bullet lists, not pipe tables (converted 2026-07-20 for
   readability) — append bullets and keep physical lines ≤ ~140 chars.
11. **A re-baseline needs an erratum checklist, not just new evidence** (B22,
   from the 2026-07-19 inspection — its D1/D2/D7/D8 stale-doc findings were
   all products of this gap). When regenerated evidence moves committed
   numbers, grep the moved values across docs/ and correct or annotate EVERY
   old-section quote in the same commit; the five-surface ritual only covers
   new sections. Full wording in CLAUDE.md workflow step 5.

Baseline: **479 passed + 25 skipped + 2 xfailed** (2026-07-20, B25 inboard
fragment clip, +6 passed = `tests/test_b1_cut_elements.py::TestInboardFragmentClip`
(4) + the same file's foot-preference lock (1) + `tests/test_m2_wingbody.py`'s
waterline-extension lock (1); **measured 1100.63 s @16 threads**).
Previous: 473 + 25 + 2 (2026-07-19, P11 curved
walls, +8 passed = the ungated `tests/test_p11_curved_walls.py`;
measured 1124.94 s @16 threads);
465 + 25 + 2 (2026-07-19, B22 3-D LS
anchor locks, +2 skipped = the gated `tests/test_b22_ls_3d_anchors.py`;
measured 1127.38 s @16 threads);
465 + 23 + 2 (2026-07-19, B21
freeze-capture alignment, +1 skipped = the gated 3-D freeze-capture lock
`test_freeze_capture_matches_live_density_3d` in
`tests/test_b15_ls_newton_freeze.py`; measured 1105.87 s @16 threads);
465 + 22 + 2 (2026-07-18, B19 LS-Newton
Jacobian exactness, +2 passed / +1 skipped = `tests/test_b19_jacobian_3d.py`
(2 structural locks + 1 gated 3-D FD gate); measured 1101.50 s @16 threads);
463 + 21 + 2 (2026-07-18, A3 inspection response, +3 passed =
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
