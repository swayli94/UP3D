# Response to the 2026-07-17 Kimi inspection (A3)

- **Date**: 2026-07-18
- **Author**: Claude (the audited authoring agent), answering
  `20260717-2113-docs-consistency-review.md`,
  `20260717-2348-code-review.md`, `20260717-2348-plan-assessment.md`.
- **Phase**: Track A / **A3** (`docs/roadmap/track_a.md`) — user-arbitrated
  scope: all three layers (docs, code hardening, the C1 verification).
- **Method**: every finding was re-verified against the CURRENT tree before any
  edit. **All of them still stood** — B16/B17/B18 (2026-07-18) had fixed only
  overview.md's baseline in passing, and touched only `newton_ls.py` /
  `picard_ls.py`, which is the already-correct side of C2/C3.

## Verdict on the audit itself

The audit was accurate and materially useful. Two things are worth saying
plainly:

1. **C1 was right, and it is the phase's real finding.** It was filed as
   *疑似* (suspected) with an explicit "not confirmed by a failing gate" caveat
   and a suggested verification. That verification now exists and **confirms
   it decisively** (below). This is the one finding that changes how a
   committed capability should be described.
2. **The audit predicted the B9 far-field follow-up correctly, one day early.**
   Its solve/ section concluded the 8 far-field rows at |R|≈84 were "real, but
   no solver-side mechanism found — consistent with the failure living in the
   far-field BC layer, not the solver". B16 and B17 (independently, next day)
   root-caused exactly that: a near-singular far-field aux block, then a
   BC-modelling error in the pin value. Two independent analyses converging on
   the same attribution is worth more than either alone.

One correction to the audit: its **finding 6** ("V14.6 exists as two number
pairs") located the mechanism exactly right, and its own conclusion — rounding
noise, immaterial — is what we adopted. We did **not** take its suggested fix
(read `LS_REF` from the CSV): that would move the committed `checks.csv` and
force a re-run of a heavy demo to chase 0.02pp. The provenance is documented at
the source instead. Recorded as a deliberate deviation, not an oversight.

---

## A. Code review — disposition

| ID | Verdict | Disposition |
|---|---|---|
| **C1** LS Newton Term 2/3 column mapping | **CONFIRMED** | **Verified, RECORDED, not fixed** — see below. Fix is a shipped-kernel change ⇒ its own phase, user's call. |
| **C2** conforming fail-fast crosses selection epochs | valid | **FIXED** — `r_level_best` reset at freeze-arm / revert / refresh (`solve/newton.py`), mirroring the LS fix. |
| **C3** no `freeze_max_reverts` disarm on conforming | valid | **FIXED** — `freeze_max_reverts=3` + `freeze_armed` flag. |
| **C4** reader drops unnamed physical surface groups | valid, **the dangerous one** | **FIXED** + 2 tests (verified failing before the fix). Its consequence chain was **re-derived in the code, not transcribed**: `wake_cut._sheet_free_edge_nodes` classifies a sheet boundary edge as *interior/free* iff it is not an edge of any OTHER `boundary_faces` group (`wake_cut.py:164-173`) — so a dropped symmetry group does make the sheet's root edge read as free, leaving those nodes single-valued and Γ(root) = 0. (P14's lesson: do not carry another analysis's attribution into your own record as if you had measured it.) |
| **C5** reader crashes on unnamed physical volume tags | valid | **FIXED** — placeholder `volume_<i>` padding. |
| **C6** far-field master branch uses `dy == 0.0` | valid | **FIXED** — mask membership is the whole test now. |
| **C7a** no guard when a TE node has no aux DOF | valid | **FIXED** — assert at construction, scoped (see the honest note below). |
| **C7b** TE detect 1e-3·h vs side-shift 1e-6·h | valid | **FIXED** — TE nodes forced into the shift; **measured 0 affected nodes on all six committed families**. |
| **T1** B3 test 2% vs gate 0.3% | valid | **FIXED** — measured 0.1441%, tightened to the gate's 0.3%; gate text and test now cite each other. |
| **T2** dropped `validate_coloring` return | valid | **FIXED**. |
| **P1** `unified.section_cp` drops `gamma` on the conforming branch | valid | **FIXED** — `gamma` threaded through, default 1.4 ⇒ byte-inert. |
| **P2** vacuum-floor Cp saturation + `validate_physics_bounds` q>2.0 | valid, dormant | **BACKLOG** — the q-limit genuinely contradicts a legal limiter-capped field (q≈2.28 vs the 2.0 raise), but its only caller is a subsonic test. Fixing it blind risks changing what "invalid" means mid-phase; it belongs with whoever wires physics monitoring into the transonic path. |
| **F0** A1 timing assert flaky at 10% | valid | **FIXED** — bound 10%→20% (schema/accounting asserts stay hard); the suite's only CWD-relative mesh path anchored to REPO_ROOT. |

**Minors accepted into backlog, with reasons** (all recorded, none silently
dropped): frozen-sweep ν monitor not strictly "BITWISE"; dead/duplicate code
(`upwind.compute_face_normals`, `isentropic.upwind_factor`,
`speed_of_sound_squared`); `n_newton` off-by-one between drivers; the
least-bad line-search reporting an untried `lam`; `picard_ls` reporting the
previous iterate's residual; `continuation` overwriting `kutta_converged` in
the single-level case; `damping_scope="supersonic"` inert on outer 0;
`levelset` panel selection on kinked TEs; `section_cp_curve`'s O(h) chord bias;
degenerate-triangle guards in post; the unfailable-assert group
(`n_fus >= 0`, artifact `.exists()`, `test_mesh_gradient:129`);
`test_p2_wake_cut`'s `except ValueError: continue`. Each is either
report-layer-only, unreachable on shipped data, or a deliberate convention —
none changes a number. They are listed here so the next phase can pick them up
without re-deriving the audit.

### C1 — the verification, and what it does and does not mean

`cases/analysis/c1_ls_jacobian_fd/run_check.py` is B6's FD gate's 3-D twin: it
reuses `tests/test_b6_newton.py`'s `_build` / `_assemble_R_J` harness verbatim,
so it measures the **shipped** Newton system, and changes only the mesh
(ONERA M6 coarse, M0.70) and the probe directions.

| probe | n DOFs | ‖Jv − FD‖ / ‖FD‖ |
|---|---|---|
| targeted — aux of cut nodes touching a mixed-side plain element | 102 | **1.146e-01** |
| control — all other aux DOFs | 1509 | **6.33e-10** |
| global free (B6's own probe direction) | 12566 | 2.47e-03 |

Eight orders of magnitude apart: the defect is exactly where C1 predicted and
nowhere else.

**The adversarial check that makes this a verdict rather than a suggestion.**
The obvious alternative explanation is FD noise — a perturbation flipping
upwind branches, i.e. non-smoothness rather than a missing term. Discriminator:
a missing Jacobian term gives an **eps-independent** relative error; FD
non-smoothness does not. Measured across three decades:

| eps | 1e-6 | 1e-7 | 1e-8 |
|---|---|---|---|
| targeted rel err | 1.532e-01 | 1.532e-01 | 1.532e-01 |

max/min = **1.00**. It is a missing term.

Also worth recording: the affected class is **larger than the audit's framing**
— 3378 mixed-side plain elements vs 428 beyond-tip ones, so this is not just
the tip strip.

**What it means, bounded honestly.** The residual R is untouched by a Jacobian
error, so **every converged level-set state, every committed γ / cl / M_max,
and every gate number stands.** What is affected is the convergence *rate* and
step count: the LS Newton is a **quasi**-Newton on 3-D meshes. Why no gate saw
it: B6's FD gate runs on the quasi-2D mesh, which structurally has no
mixed-side plain elements, and the B7/B15 M6 gates are convergence gates, not
FD gates. Caveat kept: measured at a seeded (|R| 4.8) state — legitimate,
because the FD identity must hold at *every* state, but a converged-state
repeat is the natural follow-up.

**Not fixed here.** The fix (side-aware column mapping, or using the main-field
density for these elements) is a change to a shipped kernel that would move
committed step-count trajectories. That is a phase, and the decision is the
user's.

### Inertness, verified rather than asserted

"Byte-inert on committed data" is the kind of claim that is easy to state and
easy to be wrong about, so each was measured:

- **C4** — read every committed `.msh` and diff the boundary-group set:
  **zero placeholder groups added anywhere** (all meshes come back with exactly
  their existing `farfield / symmetry / wake / wall / fuselage`). The two
  surface-only legacy assets still raise `ValueError` as before. This also
  matters because `wake_cut._sheet_free_edge_nodes` consumes *all* boundary
  groups — an extra group would have changed free-edge classification, so
  "inert" here had to mean "no group added", not just "tests pass".
- **C7b** — counted TE nodes in the 1e-6·h ≤ |s_raw| < 1e-3·h window across
  six committed families: **0 on every one** (max |s_raw| at TE nodes is
  exactly 0.0), so forcing TE nodes into the shift changes nothing.
- **C2/C3/P1** — dormant by construction (a fail-fast that never fired, a
  disarm that never triggered, a `gamma` default of 1.4), and locked by the
  existing suite including B11's bitwise unified-vs-legacy post locks.

### An honest note on C7a

The first version of the C7a assert fired on `test_b1_cut_elements.py::
test_te_fan_split` — a **legitimate** degenerate synthetic case: a lone TE-fan
tet with no cut elements at all, hence no aux DOFs and no multivalued system to
corrupt. The assert was too broad, and the existing suite caught it. It is now
scoped to `n_ext_dofs > 0`, which still covers the real trap C7 describes (a
TE node missed by a wake that *does* exist). Recorded because "the guard I
added was wrong in its first form" is exactly the kind of thing that otherwise
disappears from the record.

## B. Documentation review — disposition

All 17 findings dispositioned; details per finding in the A3 entry
(`docs/roadmap/track_a.md`). Summary: **15 fixed**, **2 fixed-by-documenting**
(finding 6, the V14.6 dual pair — provenance recorded at the source rather than
re-running a heavy demo for 0.02pp; finding 4's dated `docs/analysis/` snapshot
got an erratum rather than a rewrite, per the "dated snapshots are not
maintained documents" rule; the demo_report copy of the same timings **was**
corrected to the CSV truth 41.9 / 7.5 (5.6×) / 6.5 s).

Finding 0 (the A1 timing flake) is fixed as F0 above. Finding 2 (M2 solver leg
◐ vs CLOSED) was the only **authority-level contradiction** in the set and is
resolved in track_m's title, checkbox and ledger.

## C. Plan assessment — response

Its process-risk table is adopted, as durable rules rather than one-time fixes:

- **Close-out doc-refresh debt** → CLAUDE.md workflow step 5 is now an explicit
  **five-surface checklist** (ledger, agent-rules + baseline, overview,
  **PROJECT_STRUCTURE footer AND trees**, `cases/*/README.md` row), plus
  agent-rules discipline #10. The two surfaces that actually rotted are named.
- **Two-path fix drift** → a **backport check** line in the close-out ritual
  (CLAUDE.md step 5) and agent-rules discipline #9: when a fix lands on one
  wake path, check the other and *record the answer* — "N/A because ..." is
  fine, silence is not.
- **Gate criteria wider in tests than in gate text** → folded into discipline
  #9 with the B3 case as the worked example.
- **Flaky wall-clock asserts** → F0.
- **External-mesh reader fragility** → fixed outright (C4/C5) rather than
  documented as a precondition, since the silent-Γ(root)-pinning chain is too
  quiet to leave in place.
- **Dual number pairs** → documented at the source; see the deviation note.

On its **candidate next phases**: A3 has executed its item 1 (the hardening
sprint) including the C1 check that item 1 flagged as the one thing that could
move committed results — with the outcome that it does **not** move any
committed number, only the honest description of the LS Newton's convergence
behaviour. Items 2 (LS fine-mesh route), 3 (Track V viscous) and 4 (P11 curved
elements → G1.6) remain exactly as it framed them. Per project rule, the next
phase is the user's call; the assessment's ordering argument is on the record
and is not re-litigated here.

## D. What A3 did NOT do

- No solver numerics changed; no converged result on any committed mesh moves.
- No expensive committed artifact was regenerated (P4/P5/G4.1/G13.2 heavy demos
  and the p14 demo untouched).
- C1 is verified, **not** fixed.
- The P2 physics-bounds contradiction is backlog, not fixed.
- The C1 check ran at coarse on one mesh at one Mach; medium/converged-state
  repeats are not claimed.
