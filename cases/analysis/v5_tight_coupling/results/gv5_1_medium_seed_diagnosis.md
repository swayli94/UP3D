# GV5.1 medium seed cross-check FAIL — diagnosis (2026-07-23)

Context: the amended-protocol medium leg (Addendum 2 seed = the loose
loop's converged state, regenerated with the committed GV3.1 recipe)
stopped at its recipe cross-check (summary.csv rows 30–32, verdict
FAIL, exit 2 before any pack/FD/polish work):

| quantity | HEAD regen | committed GV3.1 medium | delta |
|---|---|---|---|
| n_outer | 3 | 5 | −2 |
| converged | True | True | — |
| cl final | 0.28137437 | 0.27190697 | 9.467e-3 |
| ds_max final | 3.454219e-3 | 6.840905e-3 | rel 0.495 |

Band: converged; |dcl_k0| ≤ 1e-8; |dcl_final| ≤ 1e-3;
|dds_max| ≤ 1e-2 rel → FAIL on the last two. The coarse cross-check
PASSED the same band (n_outer 4 vs 4, dcl_final 1.52e-4, dds_max rel
2.7e-3). This file records the four-hypothesis discrimination.

## H1 — recipe fidelity: CLEAN

Line-by-line, `run_amended_leg("medium")`
(cases/analysis/v5_tight_coupling/run.py) vs
`run_case("medium", M_INF, "picard")`
(cases/analysis/v3_loose_coupling/run.py:135-162, 376-378):

- mesh: identical path `cases/meshes/naca0012_2.5d/medium.msh` through
  `read_mesh` + `cut_wake` (v3 run.py:60,138; v5 run.py phase 1).
- config: identical constructor `CouplingConfig(re_chord=3.0e6,
  m_inf=0.5, alpha_deg=2.0)` — every default (x_tr 0.05/0.05,
  le_band_x 0.05, inflow_band_x 0.02, omega 1.0, n_outer_max 10,
  tol_ds 1e-3, ibl_tol 1e-9, ibl_max_iter 100, n_smooth_passes 2,
  eps_diff 0.005, eps_diff_s 0.02) is the shared dataclass's, no
  override in either caller.
- case: identical `build_airfoil_case(mc.nodes, mc.elements,
  mc.boundary_faces["wall"], cfg)`.
- driver: identical `make_picard_lifting_driver(mc, wc, 0.5, 2.0)`,
  no kwargs; `run_loose_coupling(driver, case, cfg, probe=probe)`.
- runtime evidence: inflow pin counts measured on HEAD — medium 14 =
  committed 14, coarse 6 = committed 6; k=0 inviscid cl matches to
  1.3e-9 (below).
- v3 run.py unchanged since the artifacts' commit (no commits in
  6303b55..HEAD on cases/analysis/v3_loose_coupling/).

No recipe difference found.

## H2 — early trajectory: divergence starts AT k=1, not k=0

Per-outer-k, committed vs HEAD regen
(gv5_1_seed_reproducibility_medium.csv):

| k | cl committed | cl regen | ds_max committed | ds_max regen | ibl_iter c/r |
|---|---|---|---|---|---|
| 0 | 0.284371833975345 | 0.2843718326662933 | — | — | — |
| 1 | 0.26223950 | 0.28096765 | 6.8228e-3 | 3.4342e-3 | 100/100 |
| 2 | 0.27364211 | 0.28135485 | 6.8217e-3 | 3.4525e-3 | 100/24 |
| 3 | 0.27318587 | 0.28137437 | 6.8283e-3 | 3.4542e-3 | 21/16 |
| 4 | 0.27181456 | — | 6.8410e-3 | — | 46/— |
| 5 | 0.27190697 | — | 6.8409e-3 | — | 23/— |

k=0 (pure inviscid FP) matches to 1.3e-9 — the mesh, the FP code, and
the recipe are the same. The divergence appears at the FIRST IBL solve
(k=1), immediately at factor ~2 in delta* (3.43e-3 vs 6.82e-3). A
wiring error would show at k=0 or in the pin counts; it does not.

## H3 — environment/code sensitivity at the IBL floor: PROVEN

Controls:

- **HEAD is deterministic**: two independent HEAD medium regens are
  bit-identical (cl = 0.28137437374745844, ds_max =
  0.0034542185646497105, n_outer = 3, same per-k IBL iteration counts
  100/24/16). Same for coarse (|dphi| = |dU| = 0).
- **The IBL code path of the committed artifacts == c2dc325**:
  6303b55..c2dc325 touches only additive files (f263424: coupling.py
  +87 = build_wing_case; c2dc325: new tight.py + tests). ibl3.py /
  closures.py unchanged.
- **closures state columns are bit-identical c2dc325 ↔ HEAD**
  (probe: 256 random states, laminar+turbulent: outs max diff 0.0;
  the Stage-2 claim verified at the residual level). douts differ at
  max 3.6e-12 (derivative-stack reordering).

Three medium trajectories from the SAME k=0 state (c2dc325 and HEAD
k=0 cl are bit-identical = 0.2843718326662933; committed k=0 differs
at 1.3e-9 → the committed artifacts predate an environment shift):

| source | k=1 cl / ds_max | converged state | n_outer |
|---|---|---|---|
| committed (6303b55, old env) | 0.26224 / 6.8228e-3 | cl 0.27191, ds 6.8409e-3 | 5 |
| c2dc325, this env | 0.27793 / 4.9240e-3 | cl 0.22170, ds 9.7279e-3 | 6 |
| HEAD, this env | 0.28097 / 3.4342e-3 | cl 0.28137, ds 3.4542e-3 | 3 |

- Same IBL code, env shift only (committed vs c2dc325): k=1 delta*
  6.82e-3 → 4.92e-3 = **27.8% apart**; the converged fixed point moves
  from cl 0.272/ds 6.84e-3 to cl 0.222/ds 9.73e-3.
- Same env, IBL refactor only (c2dc325 vs HEAD, closures outs
  bit-identical, douts at 3.6e-12, k=0 bit-identical): k=1 delta*
  4.92e-3 → 3.43e-3 = **30.3% apart**.
- c2dc325's own trajectory confirms the floor's chaos: k=3 excursion
  to cl 0.2200/ds 9.56e-3 (ds_change 28.6%), then re-converges at
  k=6 to a THIRD fixed point.

Mechanism (established at Stage 3, da27e95, on coarse): the IBL
Newton does not converge at these states — it stops at the 100-iter
cap with the residual pinned at its ~1e-6 floor and cond(J_BL,BL) ~
4e10. The stopped state lies on a near-null manifold; any
perturbation (environment rounding ~1e-16, derivative-stack
reordering ~3.6e-12) is amplified along the manifold to O(0.3) in
delta*. The medium k=1 sits deep in this basin (both the committed
run and c2dc325 hit the 100-iter cap at k=1/k=2 with 28% outer
changes); the coarse k=1 happens to sit on a well-conditioned section
(committed 4.5727e-3 vs HEAD 4.5791e-3 = 0.14% apart), which is why
the coarse cross-check passed.

## H4 — mesh integrity: CLEAN

`medium.msh` md5 = 9bfd428f675ffff83648031c3e326cc8 — identical in
the working tree, at HEAD, and at 6303b55 (the artifacts' commit).
Not gitignored; last mesh-dir commits ac4ac9d / c40cb61 predate V3.

## Verdict

**Hypothesis 3 — plateau/floor drift, NOT a recipe error.** The
committed GV3.1 medium fixed point is not reproducible from any
current checkout in this environment: an environment shift alone
(same code) moves the medium converged state by ~28% in delta*, and
the Stage-2 IBL refactors (closures outs bit-identical) move it a
further ~30%. The coarse fixed point is well-conditioned and
reproduces to 1.5e-4 in cl. The recipe, mesh, pins, and inviscid path
are all verified identical.

## Recommendation (per instruction (b): reported, NOT acted on)

The medium band question goes back to the user/parent. The medium leg
stopped before the polish; no band was widened. Options on the table:

1. Accept the HEAD-regen medium seed as the amended-protocol seed
   (Addendum 2's "regenerated" language; the committed numbers stay
   the comparison baseline, with the 9.5e-3 cl / 0.50 rel delta*
   seed-vs-committed offset recorded honestly in the compare rows).
2. Re-adjudicate the medium cross-check band around
   fixed-point-agnostic criteria (converged + k=0 exact + pin counts
   + trajectory-class membership), given the measured non-
   reproducibility of the committed medium fixed point.
3. Record GV5.1 medium as blocked-by-seed-nonreproducibility with
   this mechanism, and take the coarse amended leg (already complete:
   FD PASS at seed 2.246e-8 and at the endpoint 2.244e-8; polish
   crawled on the same IBL floor, N_polish = 10 not converged) as the
   GV5.1 evidence.
