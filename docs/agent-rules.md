# pyFP3D Agent Rules

Current phase: **B15 ✓ CLOSED 2026-07-15 (LS Newton transonic ramp + N5
freeze-selection; the M6-medium Picard plateau is gone, 38.4 → 11.0 min = 3.51×)
⇒ B9 (wing-body LS solve, M∞ 0.5) = NEXT** (user-arbitrated 2026-07-14).

B9 scope guards (user-arbitrated): **subsonic M∞ 0.5 ONLY** (M0.84 excluded —
the round-cap transonic fine is still a non-converged limit cycle, G13.3
transonic NEGATIVE); mesh = the delivered M2 wing-body family
(`cases/meshes/onera_m6_wingbody/`, wake-free, LS TE polyline endpoints exact);
**open verification item**: the B4 TE control volumes are wall-adjacent, so the
innermost TE node's fan touches fuselage wall faces — verify the upper/lower CVs
take only wing-side elements (`multivalued.py::_build_te_control_volumes`;
recorded in Track M M2).

## Track status (one line each; authority = docs/roadmap/*.md ledgers)

- **Track P** ([track_p.md](roadmap/track_p.md)): P0–P9 ✓ (P1: G1.6 strict
  xfail) · P10 ◐ (G10.1 open) · P11 conditional-not-opened · P13 ◐ (G13.3
  transonic NEGATIVE-open).
- **Track M** ([track_m.md](roadmap/track_m.md)): M0–M5 ✓ except M2 ◐ (mesh ✓,
  solver leg = B9).
- **Track B** ([track_b.md](roadmap/track_b.md)): B1–B8, B11–B13, B15 ✓ ·
  B6 ◐ (medium quantitative closed by GB15.4) · B14 designed-not-scheduled ·
  B10 shelved · **B9 next**.
- **Track V** ([track_v.md](roadmap/track_v.md)): designed, zero implementation.

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

Baseline: **396 passed + 18 skipped + 2 xfailed** (measured 988.73 s
@16 threads, 2026-07-15; lineage in [overview.md](overview.md)). After any kernel/assembly change run
`pytest tests/test_v0_freestream.py` first (CLAUDE.md hard rule 1).
