# pyFP3D — Project Instructions for Claude Code

3D unstructured-mesh **full-potential** transonic flow solver (Python + Numba):
one scalar φ per node, Galerkin P1 tets, artificial-density upwinding in supersonic
zones, wake cut + Kutta condition for lift. Target: wings at M∞ 0.3–0.87,
workstation-scale (minutes for 1–3 M nodes).

## Document map (read the relevant one before coding)

Docs were split by track on 2026-07-15; the old monolith paths remain as thin
indexes, so historical references like "roadmap.md Track B" or "demo_report §P4"
resolve through one hop.

- [docs/overview.md](docs/overview.md) — human-readable snapshot: per-track
  status table, **document map** (which file is authoritative for what, and when
  to update it), regression-baseline lineage, long-standing open items.
- [docs/roadmap.md](docs/roadmap.md) — **active tracker index**: working rules,
  gate-ID/renumbering conventions, one-line status per track. The phase entries,
  gate checklists and progress ledgers live in **[docs/roadmap/](docs/roadmap/)**
  (`track_p.md` P0–P14 solver, `track_m.md` M0–M5 meshing, `track_b.md` B1–B32
  level-set wake — **B16/B17 far-field aux pin + `pin_gamma`, B18 wing-body
  transonic, B19 LS-Jacobian exactness, all ✓ CLOSED 2026-07-18; B20
  mixed-plain main-field density ADOPTED PERMANENTLY + re-baselined and B21
  N1 freeze-capture fix (restores the M6-medium M0.84 ramp; GB20.7's
  "capability loss" verdict overturned), both ✓ CLOSED 2026-07-19; B23
  junction discriminator (the pocket = the wake inboard free-edge
  singularity), B24 waterline-extension route closed (negative) and B25
  `inboard_clip` CURES the pocket (corrM 14.66→0.63, default None
  bit-identical), all ✓ CLOSED 2026-07-19; B26 post-cure LS ceiling
  re-measured = the conforming site (medium 0.7625 / coarse 0.84 reached)
  and B27 B18 demo refresh (checks 8/8 PASS, 336/336 bit-identical;
  transonic cross-model M0.65 2.4% PASS / M0.75 2.5%), both ✓ CLOSED
  2026-07-20 — the B18 "junction-limited" story is RETIRED; B28 cl_fus
  decoupling + GB9.4 re-spec (the "fuselage spurious lift" label retired;
  out-band cross-model ≤15%, medium 7.0% PASS) and B29 flat-fragment adopted
  as the wing-body LS production config, both ✓ CLOSED 2026-07-20; B30
  (b)-class ceiling attribution (conforming stall and LS+clip death = the
  SAME wing-tip P13 free-edge singularity + high-M Newton, not a wake-model
  pocket) ✓ CLOSED 2026-07-21; B31 C-class wing-tip cure (production
  pressure+taper cures the conforming 0.83 dying level via the FD-verified
  Gamma-pin row blend; LS-side C-class closed negative) and B32 conforming
  tip_taper adopted (wing-body medium ceiling M0.79 → **M0.84 reached**,
  cl_p 0.2738, 0 clamps; weld-sign per-step refresh rolled back as
  ill-posed), both ✓ CLOSED 2026-07-22; **P11
  curved wall elements ✓ CLOSED 2026-07-19 in track_p — measured NEGATIVE,
  G1.6 re-attributed to intrinsic P1 capability at h=0.08 (not the wall
  variational crime), route fork = user's call** —
  `track_v.md` V1–V4 viscous, designed-not-started, `track_a.md` A1–A3
  verification & analysis; **A3 ✓ CLOSED 2026-07-18** = the response to the
  2026-07-17 independent inspection (docs/inspection/; the 2026-07-19
  second-round inspection's N1/D1–D10 findings were executed by B21 + the
  errata wave). **Next phase = user's call.**)
  "What phase are we in" and "what gate is open" live there, nowhere
  else. Track B numerics live in a separate spec,
  [docs/design_track_b.md](docs/design_track_b.md) (it supersedes DN1).
- [docs/design.md](docs/design.md) — theory & numerics reference: equations (§2–§3),
  wake/Kutta (§4), BCs (§5), discretization (§6), Numba kernel rules (§7), solver
  strategy (§8), V0–V6 validation ladder (§10), risks/mitigations (§12); §11 is a
  pointer to roadmap.md + docs/roadmap/ since 2026-07-15.
- [docs/demo_report.md](docs/demo_report.md) — **evidence dossier index** (per-
  phase directory table); the evidence sections live in
  **[docs/demo_report/](docs/demo_report/)** (`track_p.md`, `track_m.md`,
  `track_b.md`, `track_a.md`): one self-checking demo per phase under `cases/demo/<phase>/`
  with committed figures + measured gate numbers. When a phase closes, add its
  demo section to the matching track file and a row to the index.
  **A claim without a committed artifact is not evidence.** The 2026-07-13 audit
  found the P13 M0.84 transonic result (cl_KJ 0.2866 ⇒ "the 0.019 gap is
  resolution" ⇒ "P11's lift case is refuted") existing as *prose only* — no
  script, no CSV, no cached `.npz` — after a P11 ledger status had already been
  changed on its strength. If a run is too expensive to repeat, that is the
  reason to commit its CSV, not a reason to skip it.
- [docs/analysis/](docs/analysis/) — analysis/review reports (capability
  reviews etc.), dated snapshots, non-normative. [docs/archive/](docs/archive/)
  — historical archives (e.g. the pre-2026-07-15 agent-rules narrative);
  never a coding spec (rule 11). `docs/discussion_notes/` was **DELETED
  2026-07-14** (commit 0e4895a; history via
  `git show 8aa4aee:docs/discussion_notes/<file>`). The rule stands: plan
  against **roadmap.md/roadmap-track gates + design.md numerics only**.
- [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) — layout, per-module status, and
  **"Known gaps"**: read it before touching the G1.6 sphere-Cp problem (formerly
  G1.2; P1 gates renumbered 2026-07-06, mapping in roadmap.md) — it is already
  root-caused, and the G1.3/G1.4 oracles ruled out boundary-data corrections
  (design.md §5.1.2); do not re-propose h-refinement, recovery tweaks, Nitsche,
  or flux-data corrections. Open route: Option C gate re-spec + curved elements.
- `cases/reference_data/` — ground truth, never edit.

## Hard rules and current phase

@docs/agent-rules.md

## Workflow

1. Before coding: find the open gate in the current phase's docs/roadmap/ entry
   and plan against its acceptance criterion. Every visual gate needs a headless artifact
   (`artifacts/<gate_id>/*.png` + `summary.csv`; matplotlib `Agg`, PyVista
   off-screen — never GUI-only checks).
2. After any kernel or assembly change, run the primary regression first:
   `pytest tests/test_v0_freestream.py`
3. Full suite: `pytest tests/` — current baseline **479 passed + 25 skipped +
   2 xfailed** (2026-07-20, B25 inboard fragment clip, +6 passed =
   `tests/test_b1_cut_elements.py::TestInboardFragmentClip` (4) + the same
   file's foot-preference lock (1) + `tests/test_m2_wingbody.py`'s
   waterline-extension lock (1);
   measured 1100.63 s @16 threads;
   the full lineage lives in [docs/overview.md](docs/overview.md), do not
   re-grow it here). Skip
   semantics: the M6 `.msh` are gitignored — 16 M1 tests skip until
   `cases/meshes/onera_m6/generate_onera_m6.py` runs (~30 s); the wake-free
   families likewise (M3 medium ~40 s, M4 ~12 s); the heavy transonic/Newton
   gates (P4 medium + G4.3 sweep, P5, P8 G8.1/G8.2 + FD pocket, B6 M0.80
   dual-mesh + LS-Newton, B7 M6 3D) only run under `PYFP3D_TRANSONIC_GATES=1`
   and make up most of the skips. G8.3's CI reference is 301.66 s.
4. Numba debugging: `PYFP3D_NOJIT=1` swaps `@njit` for identity — print/pdb work.
5. **When a phase closes — the refresh checklist** (extended in A3 after the
   2026-07-17 audit found 17 consistency defects, most of them close-out debt):
   tick the gate in the phase's `docs/roadmap/track_*.md` entry, then update
   **all five** surfaces, because each has gone stale at least once by being
   "obvious enough to skip":
   1. that track file's **progress ledger** (bullet entry + track-status line;
      the track ledgers are wrapped bullet lists, not pipe tables, since
      2026-07-20 — append new phases as bullets),
   2. the **"Current phase"** block in docs/agent-rules.md **and its baseline line**,
   3. **docs/overview.md** (status bullet list + the regression-baseline lineage),
   4. **PROJECT_STRUCTURE.md** — the footer one-liner AND any directory tree
      the phase added files to (this is the one that silently rots),
   5. the **`cases/demo/README.md` table row / `cases/analysis/README.md` bullet**
      for the new demo or study.
   Keep the commit phase-scoped.
   ★ **Backport check.** This codebase has TWO wake paths (conforming
   `newton.py` / level-set `newton_ls.py`, and their Picard twins). When a fix
   lands on one, explicitly check whether the other needs it and record the
   answer in the phase entry — "N/A because ..." is a fine answer, silence is
   not. Two B15-era LS robustness fixes sat un-backported for three phases
   until an external review found them (A3 / kimi C2, C3).
   ★ **Re-baseline erratum checklist** (added by B22 after the 2026-07-19
   inspection: its D1/D2/D7/D8 findings were ALL products of this rule not
   existing). The five-surface list governs NEW sections; it does not catch
   OLD sections quoting numbers a re-baseline just superseded. So any commit
   that regenerates committed evidence must carry, in the phase entry, a
   checklist of every doc location that quotes the old numbers (grep the
   moved values — e.g. `grep -rn "0.2115\|2.4938" docs/`), each one either
   corrected in place or annotated "(pre-X value; superseded, see Y)". A
   number left standing silently is a future audit finding.
6. **Cost caution — do not recompute expensive artifacts casually.** Some
   evidence is committed precisely because regenerating it is slow: the P4 heavy
   demo figures (`cases/demo/p4_transonic/run_demo.py` under
   `PYFP3D_TRANSONIC_GATES=1`) are ~40 min of Picard; the medium G4.1 gate is
   ~17 min; the P5 M6 medium from-scratch continuation+polish is ~45–75 min
   (its solution npz `cases/demo/p5_onera_m6/results/medium_solution.npz` is a
   LOCAL gitignored cache like the .msh — the demo re-solves it when absent, or
   under `PYFP3D_P5_RESOLVE=1`; the committed PNG/CSV are the evidence); the
   ONERA M6 medium/fine `.msh` are minutes to regenerate. Treat the
   committed baseline as authoritative and only rerun the heavy part when a real
   solver/mesh/reference change would move those numbers AND you will commit the
   refresh. For routine edits, verify on the cheap coarse path. Prefer reading a
   committed CSV/PNG over recomputing it.

Gate IDs are `G<phase>.<n>` per Track P numbering (Track V: `GV<phase>.<n>`;
phases V1–V4 are distinct from the validation-case IDs V0–V6). Track P was
renumbered 2026-07-08 and 2026-07-11, Track B 2026-07-12 (×2) and 2026-07-13 —
docs before those dates use the then-current IDs (e.g. pre-2026-07-11 "P9
curved walls"/"P10 backlog" read as P11/P12). The one-line convention summary
and per-phase mapping notes live in docs/roadmap.md and the affected
docs/roadmap/ entries.
