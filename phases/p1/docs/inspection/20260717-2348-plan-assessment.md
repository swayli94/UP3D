# Roadmap / Plan Assessment (Kimi, independent audit)

- **Date**: 2026-07-17 23:48 CST
- **Auditor**: Kimi Code CLI (independent of the main authoring agent)
- **Basis**: `docs/roadmap.md` + `docs/roadmap/track_*.md` (status authority),
  `docs/agent-rules.md`, `docs/overview.md`, `docs/design.md` (validation
  ladder §10, risks §12), `PROJECT_STRUCTURE.md`, plus the two companion
  inspection reports of the same date (docs consistency, code review). This
  is an assessment of the *planning system and open items*, not a re-run of
  any gate. Judgements are labelled as such.

---

## 1. What the planning system gets right

- **Gate discipline with committed evidence** is real, not performative. The
  docs audit spot-checked the B9/P14/B14 headline number chains against
  committed CSVs and they reconcile to the last digit (see
  20260717-2113-docs-consistency-review.md, "Verified consistent"). The
  "a claim without a committed artifact is not evidence" rule (CLAUDE.md)
  exists because of a real audit finding (2026-07-13) and is enforced.
- **Negatives are recorded honestly**: G13.3 transonic stays NEGATIVE-open,
  B10 is shelved with reasons, B8's constraint-side tip cures are recorded
  dead, "the 0.019 gap is resolution" is deliberately worded "*strongly
  indicated, NOT earned*". A plan that records its failures this precisely is
  trustworthy about its successes.
- **Renumbering hygiene**: two renumbers per track are documented with
  mapping notes inside the affected phase entries, and CLAUDE.md carries the
  date-based resolution rule. The leftover references that remain are
  enumerated in the docs report (finding 7) and are annotation gaps, not
  semantic errors.
- **Cost discipline**: expensive artifacts are committed and treated as
  authoritative (P4 heavy ~40 min, P5 medium 45–75 min), with gitignored
  regenerable caches (.msh, P5 solution npz). This is the right split.
- **Two-path cross-validation architecture** (conforming vs level-set wake)
  already paid for itself twice: A2 attributed the TE Cp jump by comparing
  paths, and P14/B9 closed with cross-model agreement gates (V14.6, GB9.5).

## 2. Open items — technical assessment

### G1.6 (sphere Cp, strict xfail; oldest open gate)
Root-caused per PROJECT_STRUCTURE.md and design.md §5.1.2; the ruled-out
routes (h-refinement, recovery tweaks, Nitsche, boundary-data corrections)
are documented with the G1.3/G1.4 oracle evidence. The open route —
Option C gate re-spec + curved elements (P11) — is coherent: if the error is
geometric (faceted wall), no linear-element fix can close it, and the audit
found nothing contradicting that attribution. **Assessment**: keeping it as
strict xfail rather than silently widening tolerances is correct. P11 is
"conditional-not-opened" — the condition should be written into the P11 entry
explicitly (currently the G1.6 rationale "survives" there per roadmap.md:12).

### G10.1 open / fold-zone discipline (NACA medium M≈0.79, dcl/dM ≈ 6–10)
The rule "single-mesh regression locks only, never grid-convergence claims"
near the fold is the only defensible choice for a discontinuous response
region. Nothing to change.

### G13.3 (transonic NEGATIVE-open: round-cap transonic fine is a non-converged limit cycle)
This is the load-bearing negative: it justified B9's subsonic-only scope
guard (M0.84 wing-body excluded). **Assessment**: the negative is well
earned, but it is currently the main blocker between the project and any
transonic wing-body claim. If transonic wing-body is ever re-attempted, the
limit-cycle mechanism (not just the symptom) needs a diagnostic phase first —
the B15 freeze machinery cured a limit cycle on the LS wing-alone path, and
whether the same mechanism applies here is unrecorded.

### Fine-mesh routes (asymmetric capability)
Conforming fine (~450k dofs) has a measured recipe (AMG + tight EW + seeded
Picard; session rule 5). The LS fine route is explicitly unbuilt — B14
closed the *medium* preconditioner bottleneck and its own docs note the
Schur option is *slower* at small scale. **Assessment**: this asymmetry is
the clearest capability gap in the solver track and B14's designed use-case
is exactly this. If a fine LS run is ever needed (e.g. to make G14.7's "the
M6 fine converges" claim earnable), this is the prerequisite.

### Track V (viscous coupling): designed, zero implementation
The largest *physics* gap. With P14 having shown that 69% of the "0.019 gap"
to the KJ reference was Kutta-estimator bias, the residual disagreement with
experiment is plausibly dominated by viscous effects now — which is precisely
what Track V would model. **Assessment**: the design-first, build-later
sequencing was correct (V before the Kutta attribution was settled would
have been wasted effort). The timing argument for opening V is now *better*
than it was when it was deferred.

### B10 shelved, B6 remnant closed by GB15.4, A1/A2 closed — no issues.

## 3. Process risks found by the audits (mapped to plan-level fixes)

| Risk | Evidence | Plan-level fix to consider |
|---|---|---|
| Post-close-out doc refresh debt | 16 findings in docs-consistency report (stale footers/headers, unreconciled baselines, missing listings) | Add a 4-point refresh checklist to the close-out ritual: track ledger + agent-rules.md + overview.md + PROJECT_STRUCTURE.md (+ cases READMEs). One of them (M2 solver leg, docs finding 2) is an authority-level contradiction, not cosmetics |
| Two-path fix drift (LS fixes not backported to conforming) | code review C2, C3 (both B15-era LS fixes missing in newton.py) | "Backport check" line in close-out ritual when a fix lands on one path |
| External-mesh reader fragility | code review C4, C5 (runtime-reproduced) | If imported meshes are in scope for M-track future work, reader hardening belongs on that phase's gate list; if not, document "in-repo meshes only" as a reader precondition |
| Flaky wall-clock asserts in the suite | baseline run: test_a1_instrumentation 10.1% vs 10% | Convert absolute-time asserts to warn/trend; keep physics gates as pass/fail |
| Dual number pairs for one comparison | docs finding 6 (0.17/0.36 vs 0.15/0.34 across 6 docs) | Demos should read reference values from the committed CSV, not hard-code rounded constants (p14 run_demo.py:115 LS_REF) |
| Gate criteria wider in tests than in gate text | code review T1 (B3: test 2% vs gate 0.3%) | Gate text and test bound should cite each other |

## 4. Candidate next phases (the choice is the user's — rationale only)

1. **Hardening sprint (small, high certainty)**: backport the two B15 fixes
   to conforming newton.py (C2/C3), add the CutElementMap TE-aux guard (C7),
   fix the reader physical-group handling (C4/C5), and run the C1 FD check on
   a 3D mesh (the one finding that could move committed LS Newton results if
   confirmed). All cheap; all independently evidenced.
2. **LS fine-mesh route (B14 follow-up)**: closes the asymmetry in §2;
   prerequisite for any fine-mesh LS claim; medium effort, builds directly on
   B14's Schur+AMG machinery.
3. **Track V viscous coupling**: largest physics payoff per §2; large scope;
   the Kutta attribution that would have invalidated its baseline is now
   settled.
4. **P11 curved elements → G1.6**: retires the oldest open gate; benefits
   sphere/fuselage Cp (GB9.4's fuselage-lift XFAIL maps to G1.6, so B9's one
   open wound also routes here).

Suggested ordering: 1 first (it is cheap and two of its items are
preconditions for trusting future fine/transonic runs), then 2; 3 vs 4
depends on whether the user values physics fidelity (3) or gate hygiene (4)
more. This is a recommendation, not a verdict.

## 5. Honesty notes

- This assessment did not re-derive any phase's gate evidence; it trusts the
  committed CSVs (which the docs audit verified for the spot-checked chains).
- The "candidate phases" section is judgement based on the audits; the
  project rule "next phase = user's call" stands.
