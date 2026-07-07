# pyFP3D Agent Rules

Current phase: **P4 closed (2026-07-07)** (P3 same day; P2 + M0 2026-07-06; ledger has details) — transonic artificial density delivered: `kernels/upwind.py` (design.md §3 + four evidence-forced hardenings, each from a measured failure: **multi-hop upstream walk** — single-hop reach on the prism-split meshes is ~0.37 of the streamwise extent, starving ν below the (M²−1)/M² bound, blow-up at M0.75 for every ω/C; **q² limiter** M_cap = 3 — vacuum-bound transients otherwise become spurious dead-cell attractors; **shock-point operator** ν = max(ν_e, ν_up); all exact bitwise no-ops below M_crit), pseudo-transient diag(m/Δτ) in `solve/picard.py` (Δτ ≈ 1e-3–3e-3 FIXED — SER ramping re-ignites the limit cycle), and `solve/continuation.py` (Mach continuation + **Γ as an outer per-station secant around frozen-Γ density-converged solves** — nested exact Kutta runs away (Γ 0.115 → 4.99) and damped interleaving limit-cycles at transonic since the target map slope crosses 1), plus `post/shock.py` monitors and `cases/reference_data/naca0012_m080/` (Euler anchor + documented FP shift band; README says why no open FP table). Gates: G4.1 upper shock x/c 0.599 coarse / MEDIUM_G41_PLACEHOLDER medium vs 0.62 ± 0.03 band, monotone, ≤ 3 stations, no expansion shock; G4.2 **bitwise** subcritical no-op at M0.5 (ν = 0.0 exactly); G4.3 all 10 sweep cases with one `TRANSONIC_DEFAULTS` set, smooth trends, α = 0 cl ≈ −7e-4. **Transonic convergence semantics** (documented in picard/continuation docstrings): engineering-converged — physical M_max, zero limited/floored cells, Kutta |F| below eval noise, cl drift < 1e-3/hundreds of iterations; the sharp-shock residual is a bounded slowly-decaying tail, NOT 1e-10 — **P6 Newton is the designed cure; do not chase the tail with more relaxation, it was all tried and measured** (see demo_report §P4). Heavy gates (`G4.1 medium`, `G4.3` sweep) run under `PYFP3D_TRANSONIC_GATES=1`, excluded from the default suite (123 passed + 2 skipped + 2 xfailed, ~4.5 min — the always-on coarse transonic smoke is ~170 s of it). **Next priority: P5 (ONERA M6, needs Track M1 swept-wake meshing) and/or P6 (Newton — also fixes P4's residual tail and gate runtimes).** Also: P1 still open on one item — G1.6 sphere-Cp (strict xfail) awaits its **Option C acceptance re-spec**; curved/isoparametric wall elements separately scoped (design.md §5.1.2); do not propose h-refinement, recovery tweaks, Nitsche, or boundary-data corrections for G1.6 — ruled out with evidence — per [docs/roadmap.md](roadmap.md).

Authoritative docs:
- [docs/design.md](design.md): theory, discretization, architecture, and validation formulas.
- [docs/roadmap.md](roadmap.md): phases, gates, and progress ledger.

Hard rules:
1. After any kernel or assembly change, run `pytest tests/test_v0_freestream.py` first.
2. A phase closes only when its medium-mesh gate and the full coarse regression suite are green.
3. Numba kernels use SoA arrays only, no Python objects inside `@njit`, colored assembly for `prange`, zero allocation in hot loops, `cache=True`, and the `PYFP3D_NOJIT` switch.
4. All sparse linear algebra stays in SciPy/PyAMG.
5. Physics scalars live only in `physics/isentropic.py`.
6. Never edit files under `cases/reference_data/`.
7. Wake-cut changes must rerun topology asserts on all meshes in `cases/meshes/`.
8. Mesh generation stays vanilla Gmsh; node duplication happens only in the solver preprocessor.
9. Validate against full-potential references before Euler.
10. Keep commits phase-scoped and update the v2 progress ledger when a gate closes.
