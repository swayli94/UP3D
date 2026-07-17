# Documentation Consistency Review (Kimi, independent audit)

- **Date**: 2026-07-17 21:13 CST
- **Auditor**: Kimi Code CLI (independent of the main authoring agent)
- **Scope**: cross-document consistency of `README.md`, `CLAUDE.md`,
  `PROJECT_STRUCTURE.md`, `docs/agent-rules.md`, `docs/roadmap.md`,
  `docs/roadmap/track_*.md`, `docs/overview.md`, `docs/demo_report.md`,
  `docs/demo_report/track_*.md`, `docs/analysis/`, plus existence checks of
  referenced paths and spot-checks of claimed numbers against committed CSVs
  and `pytest --collect-only` counts.
- **Method note**: this audit did NOT re-run any solver; wall-clock numbers are
  taken from documents as-is. Every "verified consistent" claim below was
  checked against a committed CSV or an actual pytest collection, not merely
  cross-referenced between documents.
- **Regression baseline check (this machine, same day)**: full
  `pytest tests/ -q` at 16-thread cap (NUMBA/OMP/OPENBLAS) gave
  **1 failed, 441 passed, 20 skipped, 2 xfailed in 1129.12 s**.
  The single failure is `tests/test_a1_instrumentation.py::test_conforming_picard`
  — a timing-schema assert (`unaccounted time 1.590s is 10.1% of 15.760s wall`
  vs a 10% threshold). This is a marginal, load-sensitive instrumentation
  assertion, not a physics regression; with it passing, the documented
  baseline of 442+20+2 is confirmed. **Finding 0 below.**

---

## Verified consistent (checked against committed artifacts)

- **B9 number chain**: medium cl_p 0.2173/0.2188 vs LS 0.2165/0.2175 = 0.4%/0.6%
  (`cross_model_m05.csv` computes 0.37%/0.60%); GB9.4 XFAIL 16%/20%
  (0.0356/0.2173=16.4%, 0.0444/0.2165=20.5%); GB9.6 medians
  0.0036/0.0022/0.0010 (guardrail `checks.csv`); demo 7 PASS + 1 XFAIL; test
  counts via `pytest --collect-only`: conforming 9 (8+1 gated), LS 5.
- **P14 number chain**: roughness 0.0970→0.0043 / 0.0365→0.0024, TE gap
  0.2206→0.0040 / 0.1585→0.0024 (55×/67×), spike 0.1143→0.0533 vs LS 0.0743,
  D=1.80, +4.85%/69% — all match
  `cases/demo/p14_pressure_kutta/results/checks.csv` (28 rows); diagnostic
  `diag_checks.csv` 20/20; test_p14 has 15 tests.
- **B14**: 2.08×/1.43×, γ 0.088338, M_max 2.4938, "fine route unbuilt"
  consistent across track_b entry, ledger, and demo_report; test_b14 has 9.
- **Baseline arithmetic**: 406+15=421, 421+8/+1=429+19, 429+13/+1=442+20
  reconciles step by step; M1 16 and M2 23 tests verified via collection.
- **Mechanism spot-checks**: `PYFP3D_NOJIT` (`pyfp3d/__init__.py`),
  `PYFP3D_TRANSONIC_GATES` (40+ uses), `PYFP3D_P5_RESOLVE` (p5 demo),
  `PYFP3D_ARTIFACTS_DIR` (conftest) all exist in code.
- **Other numbers**: G8.3 CI 301.66 s (track_p:854 + demo_report/track_p:900);
  "G8.2 = 250 s superseded → 145 s" consistent across agent-rules / track_a /
  demo_report_a; V6 trend 6.30%→3.29%→1.41% (track_p:1375); G13.1 ×5-metric
  artifact +0.62; LS-vs-conf −7.4%; cl_KJ(h→0)=0.2050; discussion_notes
  deletion 0e4895a; STATUS.md removal; design_track_b §14/§15 exist.
- **Referenced paths exist**: `cases/demo/b9_wingbody/`,
  `cases/meshes/onera_m6_wingbody_conforming/`,
  `cases/analysis/b9_fuselage_guardrail/`, both track_a.md files, etc.
- **docs/temp/**: empty directory, no git history, covered by the `.gitignore`
  `temp/` pattern, referenced nowhere — no content, hence no conflict.

## Findings

### 0. [test robustness] A1 instrumentation test fails marginally at 10% threshold
- **Location**: `tests/test_a1_instrumentation.py:68` (`_assert_schema`),
  failing at :118 (`test_conforming_picard`).
- **Evidence**: full-suite run 2026-07-17 (16-thread cap, idle 32-core box):
  `unaccounted time 1.590s is 10.1% of 15.760s wall` vs assert `< 10%`.
- **Why it matters**: the documented baseline "442 passed" is reproducible only
  when this timing assert happens to pass; a hard 10% accounting threshold on
  wall-clock timing is inherently flaky on a shared machine. Consider a
  looser bound, a median-of-N, or xfail-strict on this specific assert.

### 1. [status contradiction] PROJECT_STRUCTURE.md footer stale, pre-dates B14/B9 close-out
- **Location**: `PROJECT_STRUCTURE.md:907-924` (footer, "Last updated: 2026-07-17").
- **Detail**: footer claims 2026-07-17 and includes P14 ✓ and Track A "A1, A2 ✓",
  yet still says "**B9 (wing-body LS solve, M∞0.5) = NEXT**", Track B row omits
  B14 ("B1–B8, B11–B13, B15 ✓"), and baseline reads "421 passed + 18 skipped +
  2 xfailed".
- **Evidence**: agent-rules.md:3 "B9 ✓ CLOSED 2026-07-17", :124 "B14 ✓ CLOSED
  2026-07-17", :194 baseline "442+20+2"; roadmap/track_b.md:986/991 same.
  The footer was refreshed between P14 (421) and B14/B9 (429/442) on the same
  day, contradicting its own date.

### 2. [status contradiction] M2 solver-leg closure inconsistent across authoritative track files
- **Location**: `docs/roadmap/track_m.md:46` (M2 title "◐"), `:164`
  ("[ ] Solver leg = Track B / B9" unchecked), `:247` (ledger M2 ◐) vs
  `docs/roadmap/track_b.md:986` (B9 ledger tail "**Closes the solver leg of
  Track M's M2**") and `docs/agent-rules.md:121-123` ("M2 solver leg CLOSED
  by B9 2026-07-17").
- **Detail**: the B9 close-out (e215ff7) ticked GB9.1–9.6 but did not tick the
  M2 solver leg in track_m.md. Per the document map ("track-file ledgers are
  the sole authority"), track_m (◐) conflicts with track_b/agent-rules
  (closed); roadmap.md:13, overview.md:23, PROJECT_STRUCTURE.md:919 all still
  show M2 ◐.

### 3. [number inconsistency] overview.md regression baseline + lineage missing B14/B9
- **Location**: `docs/overview.md:48-63`.
- **Detail**: "current baseline **421 passed + 18 skipped + 2 xfailed** (P14)";
  lineage stops at 421, missing 429 (B14) and 442 (B9), while `CLAUDE.md:72`
  says "the full lineage lives in docs/overview.md". Also the header says
  "snapshot date: 2026-07-15" (:3) yet the body contains 2026-07-17 B9/P14
  content — self-contradictory.

### 4. [number inconsistency] GB15.3 timings: demo_report still carries superseded pre-CSV values
- **Location**: `docs/demo_report/track_b.md:1022-1024` ("Picard 44.0 s /
  Newton 8.1 s (5.4×) / 6.8 s (6.5×)").
- **Evidence**: committed CSVs say **41.9 s / 7.5 s / 5.6× / 6.5 s**
  (`cases/demo/b15_ls_newton_ramp/results/checks.csv`, `summary.csv`);
  `docs/roadmap/track_b.md:992` ledger already has the corrected 41.9/7.5/5.6×.
  agent-rules.md:152 explicitly warns the archive's GB15.3 timings are pre-CSV,
  yet the same pre-CSV numbers remain in demo_report/track_b.md.
  `docs/analysis/capability_review_2026-07-15.md:29-31` also has 44.0/8.1 and
  writes "5.2×" — 44.0/8.1 = 5.43, internally inconsistent (minor; analysis
  docs are dated snapshots).

### 5. [broken reference] PROJECT_STRUCTURE.md "Known gaps" section does not exist
- **Location**: `CLAUDE.md:50-51` ('PROJECT_STRUCTURE.md — …and **"Known
  gaps"**: read it before touching the G1.6 sphere-Cp problem');
  `PROJECT_STRUCTURE.md:293, 408, 824, 913` four internal references to
  'see "Known gaps"'.
- **Evidence**: none of the 18 headings in PROJECT_STRUCTURE.md is "Known
  gaps" (the G1.6 content actually lives inside "### ✓ P1 gates G1.1 and G1.2
  closed…", :569-648). CLAUDE.md points at a non-existent named section.

### 6. [number inconsistency] V14.6 cross-model agreement exists as two number pairs (0.17%/0.36% vs 0.15%/0.34%) — mechanism located
- **Location**: 0.17%/0.36% in `agent-rules.md:44`, `overview.md:22`,
  `demo_report.md:65`, `demo_report/track_p.md:1746`, `roadmap/track_p.md:1644`,
  `track_b.md:335,438`; 0.15%/0.34% in `agent-rules.md:57-58`, `roadmap.md:12`,
  `track_p.md:1613,1721`, `PROJECT_STRUCTURE.md:917`,
  `demo_report/track_p.md:1774`, `demo_report.md:65` (both pairs in the same line).
- **Evidence (which is true)**: committed `cross_model_medium_m084.csv` has only
  4 decimals (0.2776/0.2823 vs 0.2772/0.2813), face value 0.14%/0.36%. In
  `run_demo.py`: the **G14.7 row** (:507-512) uses hard-coded rounded constants
  `LS_REF={0.2772, 0.2813}` (:115) → 0.15%/0.34%; the **V14.6 row** (:709-718)
  uses the A1 cache's **full-precision** LS values `float(d["cl_p"])` →
  0.17%/0.36%. Both pairs even appear in the committed `checks.csv` (:21 vs :27).
  **Conclusion**: 0.17%/0.36% is the full-precision-reference value;
  0.15%/0.34% is computed against pre-rounded LS_REF constants; the difference
  is rounding noise (~0.02pp, immaterial to the <1% conclusion), but the two
  pairs are used interchangeably across 6 documents with no note of provenance.

### 7. [renumbering leftover] track_p.md live text still uses "P9" in its old meaning (curved elements)
- **Location**: `docs/roadmap/track_p.md:415` (P5 title "V6 < 1% deferred to
  **P9 curved elements**"), `:513-516`, `:585-586` (G6.3 "→ **P9** (curved
  elements)"), `:748`.
- **Detail**: after the 2026-07-11 renumbering, P9 = grid-convergence
  discriminator (closed 2026-07-11) and curved wall elements = P11. The
  authoritative tracker's open item "V6 → P9" now points at a closed phase.
  Mapping notes exist at :867 and in the P11 entry, and the P5 entry is
  2026-07-08 historical text (CLAUDE.md's convention covers this), but these
  four spots carry no local annotation; moreover P13 re-attributed V6
  (overview.md:70-71 "not a P11 problem"), making the old deferral doubly stale.

### 8. [number inconsistency] roadmap.md renumbering counts wrong
- **Location**: `docs/roadmap.md:29-30`: "renumbered **twice** on Track P
  (2026-07-08 and 2026-07-11) and **twice** on Track B (2026-07-12 and
  2026-07-13)".
- **Evidence**: actually 3+3 — `track_p.md:867` "two same-day insertions"
  (07-11 ×2, plus 07-08 = 3); `track_b.md:10-26` "TWO renumbers landed the
  same day" (07-12 ×2, plus 07-13 = 3). demo_report.md:67,72 and CLAUDE.md:101
  state it correctly.

### 9. [document-map gaps] three file listings omit Track A files
- **Location**: `CLAUDE.md:32-34` (demo_report section lists only
  "`track_p.md`, `track_m.md`, `track_b.md`"); `docs/demo_report.md:5-8`
  (same enumeration); `docs/overview.md:33` (document map writes
  "[roadmap/track_{p,m,b,v}.md]").
- **Evidence**: `docs/demo_report/track_a.md` and `docs/roadmap/track_a.md`
  both exist and these same documents link them in their own Track A rows —
  the listings lag the 2026-07-15 track creation.

### 10. [stale headers] four file headers / ledger blurbs not refreshed after close-outs
- `docs/demo_report/track_a.md:7`: header still "(A1 ◐, opened 2026-07-15)" —
  A1 ✓ 2026-07-16, A2 ✓ 2026-07-17.
- `docs/demo_report/track_b.md:8`: header enumeration "B1–B5 ✓ B7 ✓ …
  B11–B15 ✓" omits **B9 ✓** (the same file has the B9 close-out at :1214).
- `docs/roadmap/track_b.md:966`: ledger status line "…B1 CLOSED 2026-07-11…;
  **next = B2**" (B2 long closed).
- `docs/roadmap/track_a.md:257`: ledger status line "**B9 stays NEXT**" (B9
  closed; the :239 scope guard is a historical arbitration record and fine,
  the status line is not).

### 11. [stale listings] two READMEs under cases/ seriously out of date
- `cases/demo/README.md:12-32`: demo table stops at b7_onera_m6, missing
  **10 current demo dirs**: m5_round_tip, b8_tip_taper_ls, b11_ls_infra,
  b12_lagged_lu, b13_lagged_picard, b14_schur_precond, b15_ls_newton_ramp,
  m6_medium_ls_workflow, p14_pressure_kutta, b9_wingbody (all on disk).
- `cases/analysis/README.md:14-16`: lists only a1_solver_bottleneck, missing
  a2_te_kutta_fidelity, b9_fuselage_guardrail, p14_te_pressure_diag.

### 12. [stale listings] PROJECT_STRUCTURE.md directory trees outdated in several places
- `pyfp3d/` tree missing `meshgen/fuselage.py`, `meshgen/wingbody.py` (M2,
  07-13/16), `solve/timing.py` (A1, 07-16) — yet the same tree has 07-17's
  `schur_ls.py` and `te_pressure.py`: half-refreshed state.
- `cases/meshes/` list missing `onera_m6_wingbody/`,
  `onera_m6_wingbody_conforming/`, `cessna/`, `nl7301_2element_2.5d/`,
  `zeroebwb/` (the last three are git-tracked legacy assets;
  tests/test_p2_wake_cut.py:182 still references cessna).
- `cases/demo/` list stops at m1 (same as finding 11); `tests/` list stops at
  test_p13, missing all test_b*/test_a1/test_m2/test_m5/test_p14/test_b9
  (57 actual test files).
- `:393` lists "cases/test_*.py [Deprecated]" but no .py files exist under
  that path anymore.

### 13. [number inconsistency] track_a.md GA1.2 reconciliation (398/405) superseded same-day by final account (399/406)
- **Location**: `docs/roadmap/track_a.md:61-65`: "re-spec'd the M2 body later
  the same day (**+2 tests, 396 → 398**) … post-A1 number as **405**".
- **Evidence**: `overview.md:58-61` and `agent-rules.md:197-200` both say
  **+3 → 399 → 406** (and 406+15=421 builds on 406). Context: track_a was
  written at the 15-MAC version (22 tests, 396+2=398), and the same-day
  far-field enlargement added one more lock (23 tests, 399) — correct when
  written, stale the same day, with no erratum in the gate text.

### 14. [number mixing] "42.6% → 2.6%" splices endpoints from two different measurements
- **Location**: `agent-rules.md:127`, `roadmap.md:14`, `track_b.md:8` header.
- **Evidence**: the same-session A/B precond share is **43.6% → 2.6%**
  (track_b.md:747 table, demo_report/track_b.md:1195); **42.6%** is A1's
  earlier independent measurement (demo_report/track_a.md:176; track_b.md:749
  prose also cites them separately). One point off, but strictly "42.6% → 2.6%"
  is no single A/B's numbers.

### 15. [process] docs/references/ is gitignored and untracked, yet listed in the overview document map
- `overview.md:41` lists `docs/references/` (López dissertation PDF) as part of
  the document system; `.gitignore` contains `docs/references`, and
  `git ls-files docs/references/` is empty — a fresh clone lacks the PDF while
  design_track_b / track_b cite it heavily by section ("López eq. 3.33–3.34",
  "López p.57").

### 16. [minor] scattered numbers/paths
- M2 mesh regeneration cost: `track_m.md:169` "~4 min" vs `agent-rules.md:98`
  "~5 min wake-free".
- `PROJECT_STRUCTURE.md:835` Quick Start hard-codes `cd /home/lrz/code/UP3D`
  (the Claude-side working-copy path; this copy is UP3D_kimi, and other
  machines differ again).

## Scope boundaries of this audit

- No solver/test-suite re-runs beyond the one full pytest pass documented in
  the header (plus `pytest --collect-only` counts); wall-clock numbers taken
  from documents as-is.
- `docs/analysis/` capability reviews are dated snapshots with their own
  disclaimers; not audited line-by-line beyond finding 4. `docs/archive/` not
  held to currency standards per the rules.
