# Project Structure

This document describes the directory layout and initialization status of pyFP3D.

## Directory Tree

```
pyfp3d/                    # Main package
├── __init__.py           # Package entry point, NOJIT mode flag
├── mesh/                 # Mesh I/O, topology, metrics, coloring
│   ├── __init__.py
│   ├── reader.py         # [P0] meshio → SoA arrays + boundary tags
│   ├── metrics.py        # [P0] volumes, gradients, face adjacency; ✓ [P3] adds
│   │                       #   precompute_element_geometry (B_e/V_e once per mesh);
│   │                       #   ✓ [P4] precompute_face_normals (adjacency-ordered outward)
│   ├── coloring.py       # [P0] element graph coloring; ✓ [P3] greedy loop numba-jitted
│   │                       #   (same visit order => identical assignment) + color_partition_csr
│   └── wake_cut.py       # ✓ [P2] wake-sheet node duplication (per-node flood-fill side
│                           #   classification, no planarity assumption; TE nodes ARE
│                           #   duplicated -- see docstring + roadmap P2 assert re-spec;
│                           #   (x,y)-grouped spanwise stations -> quasi-2D collapses to
│                           #   one scalar Γ; Kutta probe nodes; preprocess topology asserts)
│                           #   ✓ [M1] 3D swept-wake support: sheet interior FREE edges
│                           #   (tip edge) stay single-valued => Γ(tip)=0 discretely, tip
│                           #   TE corner excluded from Kutta stations; off-plane Kutta
│                           #   probe fallback; both paths exactly inert on quasi-2D meshes
├── meshgen/              # ✓ [M0/M1] Mesh generation (roadmap Track M)
│   ├── __init__.py
│   ├── extrude.py        # ✓ single-layer extrusion; globally consistent
│   │                       #   prism→3-tet split (min-global-index diagonal rule);
│   │                       #   assert_quad_split_consistency (M0 preprocessor assert)
│   ├── planar.py         # ✓ vanilla-Gmsh 2D builders: cylinder annulus, NACA0012
│   │                       #   with wake line embedded via gmsh.model.mesh.embed
│   │                       #   (gmsh imported lazily; solver tests don't need it);
│   │                       #   ✓ [M3/Track B] embed_wake=False + a size-field-ONLY ±6° corridor
│   │                       #   fan -> the wake-free family (nothing in the topology knows the
│   │                       #   wake exists; default True keeps the M0 path untouched)
│   ├── structured.py     # ✓ [phase 3 / task 2] structured / hybrid layer-block generator.
│   │                       #   ★ Its load-bearing property is BIT-IDENTICAL single-variable
│   │                       #   knobs (G0 4/4): grading is ANCHORED AT THE WALL with a DERIVED
│   │                       #   layer count (growing r_far only APPENDS layers), tangential
│   │                       #   spacing depends on n_theta alone, and it uses NO gmsh size field
│   │                       #   -- phase two's four gmsh knobs were all measured out of scope.
│   │                       #   airfoil_surface_distribution also carries a LOCAL chordwise
│   │                       #   refinement window at FIXED station count (so a shock-band leg is
│   │                       #   not confounded with a global DOF change; the cost is a measured
│   │                       #   12.4 % coarsening outside -- a BOUND, not cleanliness).
│   │                       #   Locked by tests/test_meshgen_structured.py
│   ├── fuselage.py       # ✓ [M2] simplified axisymmetric fuselage as ONE splined body of
│   │                       #   revolution (fusing primitives leaves C0 seams = spurious edges);
│   │                       #   rule-driven 5*C_ROOT length, 2-diameter ellipsoid nose, graded
│   │                       #   skin with a local-RADIUS-driven tip law (a revolve facets at h/R);
│   │                       #   add_fuselage_solid_split = the TWO-pi-revolve variant B9's
│   │                       #   conforming wing-body needs (the single revolve is unmeshable once
│   │                       #   the waterline is imprinted); ✓ [B25] make_inboard_clip builds the
│   │                       #   inboard-fragment clip polygon (fuselage surface / symmetry plane)
│   ├── wingbody.py       # ✓ [M2] wing + fuselage fused into a half model; wake-free for the
│   │                       #   level-set path, and (B9) a wake-embedded conforming variant;
│   │                       #   wing3d.py byte-untouched -> M1/M4/M5 stay bit-identical
│   └── wing3d.py         # ✓ [M1] ONERA M6 half wing: OCC two-section ruled loft (exact
│                           #   straight taper, sharp foilmod TE), spherical far field
│                           #   15 MAC, chord-plane wake sheet swept from the TE --
│                           #   occ.fragment stitches the shared TE edge, then
│                           #   gmsh.model.mesh.embed makes the volume conform (fragment
│                           #   alone does NOT); axis convention chord x / lift y / span z;
│                           #   ✓ [M4/Track B] embed_wake=False -> the sheet is built but neither
│                           #   fragmented nor embedded (it feeds the Distance size field only),
│                           #   so the tets never conform to it and no `wake` group exists;
│                           #   ★ ✓ [M5] tip_cap="round" (default "flat" = bit-identical): closes
│                           #   the wing with the HALF BODY OF REVOLUTION swept by the tip section
│                           #   about its own chord line (OCC revolve of the tip half-face about an
│                           #   edge OF that face, fused onto the loft) -- removes the flat cap's
│                           #   SHARP CONVEX EDGE, the P13/G13.3 wall singularity. The cap radius
│                           #   vanishes at the LE and TE, so it degenerates to a point at each and
│                           #   the TE line / wake sheet / tip TE corner / Kutta stations / B_SEMI
│                           #   are all UNCHANGED (that is what makes the M1 A/B controlled, and
│                           #   why no solver change was needed); h_tip sizes the cap;
│                           #   geometry helpers B_SEMI / C_ROOT / x_te / x_le / chord_at /
│                           #   TIP_CAP_RADIUS
├── constraints/          # ✓ [P2] Constraint machinery
│   ├── __init__.py
│   ├── wake.py           # ✓ master–slave elimination (A_red = TᵀAT once; Γ RHS-only via
│   │                       #   precomputed per-station vectors); kutta_targets() (per-station
│   │                       #   mean of probe jumps); ✓ [P8] reduce_operator(A) → (TᵀAT, TᵀAG)
│   │                       #   pure/no-mutation (update_matrix delegates, bit-identical) —
│   │                       #   on the Newton J its H column IS the exact wake-jump ∂R_red/∂Γ
│   │                       #   ✓ [P13/G13.2] tip_taper_factors(station_z, z_tip, form, r_c)
│   │                       #   — the spanwise loading taper F(z) that desingularizes the wake's
│   │                       #   free tip edge; geometry-only (independent of φ), so the Newton
│   │                       #   Kutta row is just scaled by F. Shipped: "vanish_smooth"
│   │                       #   (smoothstep, COMPACT support), r_c = 0.05·b_semi. Consumed via
│   │                       #   solve_newton_lifting(tip_taper=…), default None = bit-identical
│   ├── te_pressure.py    # ✓ [P14] probe-free Kutta estimator: TEControlVolumes — the B4
│   │                       #   wall-adjacent upper(slave)/lower(master) TE fans on the
│   │                       #   conforming cut mesh (EXACT wall-face ownership), per-side
│   │                       #   volume-weighted P1 recovery, pressure-equality residual
│   │                       #   |q_u|²−|q_l|² per station, frozen-mean implied targets (the
│   │                       #   secant drivers' drop-in), exact Newton rows (K_p, D=∂F/∂Γ).
│   │                       #   Self-contained on (mesh_cut, wc); imported ONLY when
│   │                       #   kutta_estimator="pressure" (probe default untouched)
│   └── dirichlet.py      # ✓ far-field freestream + 2D vortex correction (branch cut ON
│                           #   the wake sheet; eliminated ⁺-side far-field wake nodes
│                           #   automatically consistent); ✓ [P3] Prandtl-Glauert scaling
│                           #   (beta stretches only the atan2 argument; beta=1 bit-exact);
│                           #   ✓ [P5] spanwise_gamma=True: Γ(z)-tapered vortex (per-station
│                           #   interpolant, 0 at/beyond the sheet tip — removes the spurious
│                           #   branch-ray jump beyond the tip; default False bit-identical)
├── physics/              # Physics constants and constitutive relations
│   ├── __init__.py
│   └── isentropic.py     # ✓ [P0] ρ(q²), M(q²), a(q²), Cp, etc.; [P2] adds
│                           #   pressure_coefficient_incompressible (Bernoulli limit);
│                           #   ✓ [P3] density_field/mach_squared_field array sweeps
│                           #   (ρ ≡ 1.0 BITWISE at M∞=0 -- the G3.3 anchor);
│                           #   ✓ [P7] mach_squared_derivative_wrt_q_sq (dM²/dq², the
│                           #   ∂ν/∂q² chain factor of the frozen-walk sensitivity)
├── kernels/              # Element-wise assembly kernels (Numba-jitted)
│   ├── __init__.py
│   ├── gradient.py       # ✓ [P3] prange velocity sweep: grad(phi)_e + q²_e from B_e,
│   │                       #   zero-alloc (outputs preallocated by PicardOperator)
│   ├── jacobian.py       # ✓ [P3] Picard matrix (6.2): symbolic CSR pattern, elem_to_csr
│   │                       #   scatter map, colored-prange data kernel, PicardOperator
│   │                       #   per-mesh workspace (B_e/V_e/coloring/pattern/buffers once);
│   │                       #   accumulation order fixed by color sequence => bit-deterministic
│   │                       #   across calls/threads; ✓ [P8/N2] assemble_newton_jacobian:
│   │                       #   exact (6.3) at frozen selection — Terms 1+2 fused on the
│   │                       #   SHARED Picard pattern (Term-2 footprint = Term-1; Term 2
│   │                       #   added only when s_e≠0 so masked elements reduce to A bitwise
│   │                       #   — a fused expression would FMA-contract differently), Term 3
│   │                       #   (rows e → cols u(e), graph-dist ≤4) as ACTIVE-SET COO
│   │                       #   (16 entries per s_u≠0 element, rebuilt per Newton step ⇒
│   │                       #   selection churn can't corrupt a reused pattern); records
│   │                       #   newton_nnz/n_term3_active (the N2 measurement); JVP
│   │                       #   FD-verified ~1e-10 (test_p8_jacobian.py)
│   ├── residual.py       # [P1] serial reference kernels (KEPT, regression-tested against)
│   │                       #   + ✓ [P3] assemble_residual_colored; assemble_stiffness_matrix
│   │                       #   now delegates to the fast path (P1/P2 drivers share it) → [P8] Newton
│   ├── upwind.py         # ✓ [P4] artificial density (3.1)-(3.2): MULTI-HOP upstream walk
│                           #   (single-hop reaches only ~1/3 extent on prism-split meshes --
│                           #   measured dissipation starvation), shock-point operator
│                           #   nu = max(nu_e, nu_up), rho_tilde floor, UpwindOperator workspace;
│                           #   exact bitwise no-op below M_crit (G4.2);
│                           #   ✓ [P6] opt-in streamline-Gaussian kernel flux (weighted=True,
│                           #   mode="kernel"; a Picard-speed path — NOT the sawtooth fix,
│                           #   which is the P6 recovery smoothing in post/surface.py);
│                           #   ✓ [P7] rho_tilde_sensitivities_sweep + UpwindOperator.
│                           #   rho_tilde_sensitivities: exact branch-wise (s_e, s_u) =
│                           #   ∂ρ̃/∂q² of the WALK flux at FROZEN u(e) (López B.3–B.6;
│                           #   floor branch → 0), FD-verified to ~4e-10 — the P8 Newton
│                           #   prerequisite; ✓ [P8/N5] classify_upwind_branches +
│                           #   rho_tilde_frozen_sweep/_sensitivities_sweep +
│                           #   freeze_upwind_state/rho_tilde_frozen: the flux at a FROZEN
│                           #   (u(e), branch) assignment — bitwise = live sweep at the
│                           #   freeze state, smooth within the assignment (no max-tie kink),
│                           #   floor-free on branches 0–2 (driver reverts on divergence) —
│                           #   the Newton finish phase on tie-degenerate prism meshes
│                           #   Term-2/Term-3 physics factor (forward path byte-identical)
│   ├── __init__.py       #   exports WakeLevelSet / CutElementMap / MultivaluedOperator
├── viscous/              # ✓ [Track V / V1] IBL3 (Drela 2013 integral boundary layer,
│   │                       #   design_track_v.md) — standalone prescribed-u_e stage shipped;
│   │                       #   GV1.1 9 PASS / 2 FAIL, V1 ✓ CLOSED 2026-07-22 (VERDICT
│   │                       #   bench/studies/v1_ibl3_standalone/VERDICT.md); V1 does NOT touch
│   │                       #   solve/ (pure additive package) — V2 adds the solve/ RHS
│   │                       #   channels, V3 the loose-coupling driver (coupling.py)
│   ├── __init__.py
│   ├── surface_mesh.py   # ✓ compact wall-surface DOF numbering + per-node local basis;
│   │                       #   wake-slot reservation + master-map hook in the data layout
│   ├── closures.py       # ✓ laminar + turbulent closure packet, analytic state derivatives
│   │                       #   (N_OUT=30 incl. stress-flux integrals), safety floors,
│   │                       #   blasius_seed; family fixed point H*≈2.7083 (≠ Blasius 2.59)
│   ├── ibl3.py           # ✓ 6-equation surface Galerkin P1 FE: strong-form divergence +
│   │                       #   only-diffusion-by-parts (D13 (74)), colored prange assembly,
│   │                       #   analytic CSR Jacobian, physical-density PTC (F_pt merit);
│   │                       #   [V5/GV5.5] default-OFF te_extrapolate TE-outflow row
│   │                       #   replacement (rows 6i+0/6i+2, exact J rows, in-pattern guard)
│   ├── transpiration.py  # ✓ [V2] δ*→ṁ = ∇_Γ·(ρ_e u_e δ*) (gradN strong-form + node_area
│   │                       #   lumping) + wall-RHS Galerkin assembly (wall_correction template,
│   │                       #   b = −load(ṁ) blowing-positive, sign pinned by GV2.1(a)) +
│   │                       #   per-zone u_e extraction per A4; GV2.1 23 PASS / 0 FAIL,
│   │                       #   V2 ✓ CLOSED 2026-07-22 (bench/studies/v2_transpiration_channel/)
│   ├── coupling.py       # ✓ [V3] loose viscous–inviscid coupling: CouplingCase builders
│   │                       #   (airfoil: LE-band Dirichlet pinning + station chain; closed
│   │                       #   body: nose+tail stagnation-band pinning, tail-ṁ masking;
│   │                       #   wing [V5]: local-x/c LE pin + TE outflow + tip-band pin/ṁ mask) +
│   │                       #   run_loose_coupling outer loop (IBL→δ*→ṁ→FP→u_e);
│   │                       #   [V5/GV5.5] CouplingConfig.te_extrapolate (default OFF) +
│   │                       #   te_outflow_pairs(case) TE↔upstream pairing;
│   │                       #   GV3.1/3.2 bench/studies/v3_loose_coupling/,
│   │                       #   GV3.3 phases/p1/cases/analysis/v3_fuselage_smoke/,
│   │                       #   GV5.0 bench/studies/v5_m6_bridge/ (✓ EXECUTED 2026-07-23);
│   │                       #   [V6/GV6.1] CouplingConfig.wake_transpiration (default
│   │                       #   OFF = legacy bit-identical) + run_loose_coupling wake=
│   │                       #   hook: b_wake added into the wall body_source_rhs;
│   │                       #   [V6/GV6.2] CouplingConfig.wake_l_rel_chords (default
│   │                       #   1.0 c pinned = bit-identical) threaded to the producer
│   ├── wake_sheet.py     # ✓ [V6/GV6.1] conforming wake-sheet δ* source: wake
│   │                       #   SurfaceMesh (V1 layout, 2nd instance) + open-chain
│   │                       #   station table + master→slave fold pairing (W3 build
│   │                       #   asserts); producer (i) prescribed δ*_wake (TE
│   │                       #   confluence + straight-wake relaxation, L_rel = 1.0c
│   │                       #   pinned, W2 1e-12 assert every call); two-sided
│   │                       #   u_e,wake; b_wake = assemble over wake± with ½ṁ per
│   │                       #   coincident face copy (2026-07-25 addendum — the Tᵀ
│   │                       #   fold sums both copies, so ½ realizes [ρv_n] = ṁ);
│   │                       #   GV6.1 6 PASS / 0 FAIL / 7 RECORDED, ✓ CLOSED
│   │                       #   2026-07-25 (bench/studies/v6_1_wake_sheet/)
│   ├── tight.py          # ✓ [V5] fixed/linear operators of the augmented tight system:
│   │                       #   W wall-load / S scatter / P pin-mask / L surface-divergence /
│   │                       #   D closure-row / G per-zone u_e-recovery operators + the
│   │                       #   Jacobian assemblies J_φ,BL / J_BL,φ (edge-data chain) /
│   │                       #   J_φφ augmentation (dṁ/dφ through ρ_e·u_e)
│   ├── tight_driver.py   # ✓ [V5] TightPack + block/augmented residual + Jacobian +
│   │                       #   newton_tight (splu + P8/P14 backtracking; probe guard =
│   │                       #   the IBL halving-on-nonfinite idiom); GV5.1 ✓ EXECUTED
│   │                       #   2026-07-23 (bench/studies/v5_tight_coupling/)
│   └── strip2d.py        # ★ [GS4.1 phase 4, round 1 2026-08-18] 2-D chordwise STRIP,
│                           #   marched along the streamwise coordinate: momentum +
│                           #   kinetic-energy integrals in conservative form, unknowns
│                           #   (delta, A[, Ctau1]) = the closure's OWN state, implicit
│                           #   M y' = F via its analytic state Jacobian, RK4 with
│                           #   substeps budgeted by share of log(x). NEVER assembles a
│                           #   global system -- that is the mechanism it is cheaper than
│                           #   ibl3.py by. CALLS closures.py (no closure formula or
│                           #   constant of its own, source-asserted); ibl3.py untouched
│                           #   (roadmap freeze). similarity_fixed_point() solves the
│                           #   family's self-similar state ALGEBRAICALLY (no march, no
│                           #   discretization) = the guard separating a marching defect
│                           #   from a closure property.
│                           #   ★★ ROUND-1 VERDICT = V-FAIL, cause in the CLOSURE, kill 1
│                           #   fired: the flat-plate fixed point is H 2.708292 (+4.52 %,
│                           #   inside ±5 %) but c_f√Re_x 0.710235 (+6.94 %, OUTSIDE) --
│                           #   the MIRROR of the pre-registered prediction. ★★ Laminar c_f
│                           #   had never been compared to an analytic value AT ALL: the
│                           #   number sat in phase 1's own committed CSV
│                           #   (v1_ibl3_standalone/results/gv1_1a_march.csv, cf_march),
│                           #   +15.8 %..+7.75 % off Blasius across its window and
│                           #   converging on this round's 0.710235 -- but GV1.1's gate (a)
│                           #   gated H and its gate (b) gated TURBULENT c_f against the
│                           #   closure's own march (common mode: two arms, one closure, so
│                           #   blind to a shared model error). Data on disk, no criterion
│                           #   covering it. ★ Do NOT put gate (b)'s 0.07 % next to a
│                           #   laminar number -- that is a cross-regime comparison.
│                           #   Discretization is NOT the cause (self-convergence order
│                           #   4.04, error 1.3e-13 = 11 decades below the gap).
│                           #   F-SIMILAR PASS 5.789 % (2/3 wedges); turbulent RECORDED.
│                           #   Evidence bench/studies/gs41_strip_core/ (4.58 s);
│                           #   locks tests/test_gs41_strip2d.py (21). ★ Whether to change
│                           #   the closure is the NEXT round + user adjudication.
├── solve/                # Linear and nonlinear solvers
│   ├── __init__.py
│   ├── linear.py         # [P1] Dirichlet elimination + CG/PyAMG preconditioner (done);
│   │                       #   ✓ [P3] build_amg_preconditioner pins pyamg's spectral-radius
│   │                       #   RNG seed -- ALL solves bit-reproducible run-to-run (G3.3);
│   │                       #   ✓ [P8/N3] solve_gmres (restarted, preconditioned, auto-retry
│   │                       #   at 2× restart) + build_ilu_preconditioner — the Newton
│   │                       #   Jacobian is NONSYMMETRIC in supersonic zones (Term 3), CG
│   │                       #   does not apply; CG/AMG paths untouched
│   ├── picard.py         # [P1] Laplace driver (done); ✓ [P2] solve_laplace_lifting():
│   │                       #   Kutta outer loop, matrix+AMG built once, Γ updates RHS-only,
│   │                       #   secant (Aitken) acceleration -> 2 updates on the linear driver;
│   │                       #   ✓ [P3] solve_subsonic (non-lifting density Picard) +
│   │                       #   solve_subsonic_lifting (NESTED: outer ρ update, inner P2
│   │                       #   secant Kutta at frozen ρ -- interleaved Γ steps rejected with
│   │                       #   measurements; AMG reuse every 4 outers; opt-in forcing-term
│   │                       #   inexact inner solves, default off; M∞=0 bitwise ≡ P2);
│   │                       #   ✓ [P4] upwind_c/m_crit (ρ̃ in matrix+residual), q² limiter
│   │                       #   (m_cap), pseudo-transient diag(m/Δτ), omega_rho, kutta_per_outer,
│   │                       #   phi_init/gamma_init continuation seeds
│   │                       #   ✓ [V2] body_source_rhs threaded through solve_subsonic /
│   │                       #   solve_subsonic_lifting (the transpiration wall RHS; lifting
│   │                       #   rides the reduced_rhs Tᵀ reduction; None = bit-identical)
│   ├── wall_correction.py # ✓ [P1/G1.3] true-normal weak-flux correction RHS (Option A);
│   │                       #   assembly-verified; correction itself RULED OUT by the
│   │                       #   G1.3/G1.4 oracles (design.md §5.1.2) -- kept as reusable
│   │                       #   facet-integral infrastructure (e.g. Gap-SBM terms)
│   ├── curved_wall.py     # ✓ [P11] curved wall-adjacent elements (tet10 geometry via
│   │                       #   closest_point_normal midpoint projection, mapped-P1 field,
│   │                       #   ΔA delta assembly; opt-in via solve_laplace stiffness_delta).
│   │                       #   Route measured NEGATIVE for G1.6 (superparametric O(h)
│   │                       #   consistency; see "Known gaps" P11 block) -- kept as
│   │                       #   curved-geometry infra + evidence machinery
│   ├── continuation.py   # ✓ [P4] Mach continuation + transonic Γ closure: outer per-station
│   │                       #   SECANT on F(Γ) = kutta_target(density-converged φ at fixed Γ) − Γ
│   │                       #   around frozen-Γ pseudo-time solves (nested/interleaved Kutta both
│   │                       #   fail at transonic -- measured; see module docstring);
│   │                       #   ✓ [P5] n_kutta_polish: fixed-Γ Kutta-closure polish after the
│   │                       #   continuation (secant-free damped fixed point, omega_rho_polish)
│   │                       #   -- the 3D secant leaves the stiffest station under-circulated
│   │                       #   and DIVERGES if pushed (INVESTIGATION_kutta_closure.md);
│   │                       #   default 0 = P4 path bit-identical
│   ├── timing.py         # ✓ [A1] the canonical timings schema shared by all four nonlinear
│   │                       #   drivers (seed/assembly/precond/linsolve/residual/kutta/other/wall)
│   │                       #   + step_records; new_timings/snapshot/step_delta helpers
│   └── newton.py         # ✓ [P8/N4] fully-coupled (φ_red, Γ) Newton driver (design.md §8.1):
│                           #   NewtonWorkspace (free/dir split, Kutta row K, affine far-field
│                           #   basis vals0_red + V_red·Γ via unit-Γ probing), ONE shared
│                           #   eval_residual path, exact δΓ elimination (J_ff + B·K) δφ =
│                           #   −R − B·F with B = J_red[free,dir]@V_red + H_J[free,:] (the
│                           #   easy-to-miss far-field vortex column INCLUDED, FD-guarded),
│                           #   GMRES + Term-1-AMG + Eisenstat–Walker, safety-only line search,
│                           #   optional consistent ptc_dtau;
│                           #   ✓ [P8/N5, G8.1] transonic robustness chain: precond="direct"
│                           #   exact steps (splu+Woodbury; the refining shock-position soft
│                           #   mode stalls η-accurate Krylov steps), stall-adaptive FREEZE of
│                           #   the upwind assignment + active-set refresh (two-cycle
│                           #   acceptance, honest residual_unfrozen floor), freeze-revert /
│                           #   level-fail-fast / best-of-tried line-search safety nets;
│                           #   solve_newton_transonic = upward Mach continuation with dm
│                           #   halving (recipe: tests/test_p8_newton NEWTON_TRANSONIC_RECIPE);
│                           #   ✓ [P8/N6, G8.2] lagged-LU direct mode: direct_refactor_every
│                           #   (default 1 = bit-identical) + direct_reuse_rtol — stale-LU-
│                           #   preconditioned GMRES on the fresh coupled operator, refactor
│                           #   fallback on GMRES failure (true-3D splu fill is ~100× the 2.5D
│                           #   cost; M6 medium 1606 s → 249 s; M6 recipe:
│                           #   tests/test_p8_newton NEWTON_M6_RECIPE, dm 0.05 + spanwise Γ);
│                           #   ✓ [P10/G10.2] level-adaptive intermediate continuation
│                           #   tolerance: solve_newton_transonic(intermediate_tol=…) opt-in
│                           #   (default None bit-identical) — loose acceptance of ORIGINAL-
│                           #   SCHEDULE intermediate levels (tol_residual_loose after ≥1
│                           #   step / 1e3 rel-drop / stall-accept, freeze off; dm-halving
│                           #   retries + final level stay strict). A/B: M6 medium +40.3%
│                           #   locks-intact (promoted into NEWTON_M6_RECIPE); fold-zone NACA
│                           #   medium NEGATIVE (untracked Γ seed) — contraindicated near
│                           #   folds, NEWTON_TRANSONIC_RECIPE unchanged;
│                           #   ✓ [V2] NewtonWorkspace(external_rhs=…) lagged external-RHS
│                           #   channel (R_free −= (Tᵀb_ext)[free]; Jacobian untouched —
│                           #   GV2.1(c) bit-invariant + FD-exact; None = bit-identical)
└── post/                 # Post-processing
    ├── __init__.py
    ├── vtk_out.py        # [P0] Write .vtu for ParaView; also the PNG/CSV gate-artifact helpers
    │                       #      (export_error_heatmap, export_matplotlib_plot) live here, not
    │                       #      in a separate artifacts.py
    ├── shock.py          # ✓ [P4] shock monitors: Cp* sonic-crossing detection, shock x/c,
    │                       #      jump-width cell count, monotonicity + expansion-shock detectors
    ├── surface.py        # [P1] nodal_gradient_recovery() (volume-weighted, for interior fields)
    │                       #      ✓ [P3] m_inf param: exact isentropic Cp (2.5) in the force integral
    │                       #      and wall_tangential_gradient() (surface-only, for wall Cp --
    │                       #      the accurate one; see "Known gaps" for why it still isn't
    │                       #      accurate *enough* to close G1.6); ✓ [P2] adds triangle-wise
    │                       #      wall force integration (owner-tet-oriented outward normals,
    │                       #      no nodal averaging across the sharp TE) and KJ sectional cl
    ├── section_cut.py    # ✓ [G1.3→P2] z = const section extraction: degenerate single-layer
    │                       #      path + [P2] general marching-tets interpolation path and
    │                       #      wall_cp_curve() sectional Cp(x/c) upper/lower split;
    │                       #      ✓ [P5] section_cp_curve() derives the LOCAL chord/x_le from
    │                       #      the cut itself (swept, tapered planform)
                            #      node carries TWO values, so wall triangles must be told WHICH
                            #      copy to read (★ D11, by the outward normal's lift-axis sign:
                            #      n_y > 0 = upper). Reading phi_main on both surfaces makes the
                            #      pressure integral junk (measured cl_pressure = −3.35 vs 0.28).
                            #      wall_cp_levelset / surface_curve_levelset / cl_pressure_levelset
                            #      (2.5D, normalised by the span extent);
                            #      ✓ [B7] section_cp_curve_levelset (the D11 per-side plane cut —
                            #      ★ its UPPER surface is BIT-IDENTICAL to section_cp_curve fed
                            #      main_potential, so every gate shock metric is unaffected; the
                            #      LOWER surface is where D11 bites) + cl_pressure_3d_levelset
                            #      (planform-area normalisation, pairs with cl_kj_3d for V6);
                            #      ✓ [B11] shares surface.py/section_cut.py private cores
                            #      (_cp_from_q2, _pressure_force, _wall_plane_crossings,
                            #      _d11_wall_state) — the three near-duplicate blocks collapsed
    └── unified.py        # ✓ [Track B/B11] the unified upper post-processing layer over BOTH
                            #      paths: wall_cp / wall_forces / section_cp, keyword-dispatched by
                            #      phi= (conforming) vs mvop=,phi_ext= (level-set); outputs are
                            #      np.array_equal to the legacy per-path functions (which are kept)

cases/                     # Test cases and reference data
├── meshes/               # Mesh families (coarse/medium/fine)
│   ├── sphere_shell/     # [P1] Gmsh sphere-shell case for gate G1.6 (coarse/medium generated);
│   │                       #   [P3] reused for G3.1 (compressible-vs-PG same-mesh comparison);
│   │                       #   [P11] sweep meshes h05/h03/h02/h03_far10.msh gitignored,
│   │                       #   regenerated by cases/demo/p11_curved_walls/run_demo.py (~8 min)
│   ├── cylinder_2.5d/    # ✓ [M0] Single-layer extruded cylinder-flow test case
│   │                       #   (generate_cylinder.py; coarse 6.9k / medium 17.3k tets
│   │                       #   committed; analytic Cp = 1 - 4 sin^2(theta) validation);
│   │                       #   G1.3 pre-study ran here (fine.msh added, 50.2k tets); found
│   │                       #   recovery-dominated (~76%), NOT sphere-crime-dominated --
│   │                       #   de-designated as G1.6 testbed (design.md §5.1.2)
│   ├── naca0012_2.5d/    # ✓ [M0] Single-layer extruded NACA0012 + embedded wake sheet
│   │                       #   (generate_naca0012.py, one parameter h_wall per level;
│   │                       #   coarse 16.4k / medium 61.8k tets committed, fine on demand)
│   ├── cylinder_hex_2.5d/  # ✓ [phase 3 / task 2] structured O-grid cylinder, hex->tet via
│   │                       #   extrude_single_layer; coarse/medium/fine .msh COMMITTED (matching
│   │                       #   the sibling unstructured family) + stats CSV + layer PNG with a
│   │                       #   NEAR-WALL ZOOM (the near-wall layer is the whole claim; a
│   │                       #   full-domain view cannot show it)
│   ├── naca0012_hex_2.5d/  # ✓ [phase 3 / task 2] HYBRID architecture (user-specified, the
│   │                       #   production one): finite-height structured layer block hugging the
│   │                       #   wall, truncated, gmsh unstructured beyond. TE marching direction
│   │                       #   pinned to +-y (a sharp TE's averaged normal points +x, which would
│   │                       #   march down the wake line) and the TE de-clustered. AR max 5.9
│   │                       #   (the v1 O-grid was 1305). coarse/medium/fine committed + 3-panel
│   │                       #   PNG (full / TE / LE zoom)
│   ├── rae2822_2.5d/     # ✓ [V5/GV5.2] RAE2822 family via the point-set airfoil path
│   │                       #   (meshgen/planar.py::airfoil_wake_2d + load_airfoil_ordinates
│   │                       #   PCHIP resample; rae2822.dat = the Cook/AGARD-AR-138 Table 6.1
│   │                       #   ordinates, positive-down lower column, sha256 in the generator
│   │                       #   header; coarse 5560 nodes / medium 20790, .msh + stats + layer
│   │                       #   PNGs committed)
│   ├── onera_m6_wingbody/          # ✓ [M2] wing-body half model, WAKE-FREE (level-set path);
│   │                       #   coarse/medium, .msh gitignored (~4-5 min regen), stats CSVs +
│   │                       #   inspection PNG committed
│   ├── onera_m6_wingbody_conforming/  # ✓ [B9] the same body with the wake sheet EMBEDDED
│   │                       #   (conforming path; Netgen OFF -- it segfaults on this geometry)
│   ├── cessna/           # legacy git-tracked surface asset (referenced by
│   │                       #   tests/test_p2_wake_cut.py); not part of any gate ladder
│   ├── nl7301_2element_2.5d/  # legacy git-tracked two-element asset; no active gate
│   ├── zeroebwb/         # legacy git-tracked BWB asset; no active gate
│   ├── onera_m6/         # ✓ [M1] ONERA M6 swept/tapered half wing + embedded wake sheet
│   │                       #   (generate_onera_m6.py, one parameter h_wall, 2x ladder:
│   │                       #   coarse 55.5k / medium 350.7k / fine 2513k tets; .msh files
│   │                       #   gitignored (large) -- regenerate coarse+medium ~30 s; the
│   │                       #   stats CSVs + inspection PNGs are the committed evidence;
│   │                       #   M1 tests skip when the meshes are absent)
│   ├── onera_m6_roundtip/ # ★ ✓ [M5] the ROUNDED-TIP M6 family (generate_onera_m6_roundtip.py;
│   │                       #   wing3d.py tip_cap="round"). The flat cap's sharp convex edge was
│   │                       #   the LAST thing blocking a 3D Richardson (P13/G13.3: its box peak
│   │                       #   Mach DIVERGES, p = +0.321). Gate metric = the SEAM CREASE ANGLE on
│   │                       #   the tip section (post/surface.py::wall_crease_angles): flat
│   │                       #   91.9° -> 92.1° under refinement (a real edge -- refinement resolves
│   │                       #   it and removes nothing) vs round 46.8° -> 25.0° (O(h) faceting of a
│   │                       #   smooth surface). Self-similar ladder from the start (no h_far
│   │                       #   clamp) with h_tip = 0.25 h_wall -- WITHOUT that the 22 mm cap would
│   │                       #   be one element wide at coarse and discretize back into a flat one.
│   │                       #   59.4k / 448k / fine tets (×1.28 of M1 at equal h_wall, level-
│   │                       #   independent); .msh gitignored, stats CSVs + PNGs committed
│   ├── naca0012_wakefree_2.5d/  # ✓ [M3] the WAKE-FREE ("O-grid analogue") NACA family — the
│   │                       #   other half of Track B's DUAL-MESH rule: no wake surface is
│   │                       #   embedded, nothing in the topology knows the wake exists, so the
│   │                       #   level set makes GENERIC cuts through generic elements (the actual
│   │                       #   workflow target). planar.py embed_wake=False + a size-field-ONLY
│   │                       #   ±6° corridor fan covering the α-sweep envelope; coarse committed,
│   │                       #   medium/fine gitignored (~40 s regen)
│   └── onera_m6_wakefree/ # ✓ [M4] the wake-free ONERA M6 family (the 3D half of the dual-mesh
│                           #   rule; wing3d.py embed_wake=False -- the chord-plane sheet feeds
│                           #   only the Distance size field, is never fragmented/embedded, and
│                           #   no `wake` tag exists). ★ Sized to land within 6–9% of M1's tet
│                           #   count at equal h_wall — that equal-sizing property is what makes
│                           #   the B7 A/B against P5/P8 a CONTROLLED comparison. No α-wedge
│                           #   corridor in 3D (the wedge volume scales with span: a ±3° envelope
│                           #   would ~4× the tets), so 3D α re-aiming stays in the near-nominal
│                           #   band; .msh gitignored, stats CSVs committed
├── reference_data/       # Ground truth (DO NOT EDIT)
│   ├── naca0012_incompressible/  # ✓ [P2] Hess–Smith panel reference (generator script +
│   │                             #   cl_reference.csv / cp_alpha4.csv / convergence.csv +
│   │                             #   README provenance; two independent lift routes agree
│   │                             #   to 0.09%, panel-count converged)
│   ├── naca0012_m05/     # ✓ [P3] the same panel solution under Prandtl-Glauert AND
│   │                             #   Karman-Tsien corrections (G3.2 reference = PG/KT
│   │                             #   midpoint + inside-bracket assert; README provenance)
│   └── naca0012_m080/    # ✓ [P4] transonic shock reference: Euler anchor (~0.60c upper,
│                                 #   ~0.35c lower) + documented conservative-FP aft-shift band;
│                                 #   README records that no open FP table was retrievable
├── demo/                 # ✓ Per-phase evidence demos (docs/demo_report.md; one
│   ├── README.md         #   self-checking run_demo.py + committed results/ per phase)
│   ├── _common.py        #   shared chart style + CheckList acceptance recorder
│   ├── p0_infrastructure/  # G0.1-G0.4: volume/gradient exactness, coloring, VTK I/O
│   ├── p1_laplace/         # V0/G1.1/G1.2 + G1.4 oracle (absorbed) + G1.6 open-gate XFAIL
│   ├── p2_kutta_lifting/   # G2.1-G2.5: cut exactness, Kutta, cl vs panel, spanwise decay
│   ├── p3_subsonic/        # G3.1-G3.3: assembly-debt evidence, sphere-vs-PG, cl bracket,
│   │                       #   monotone nested Picard, M=0 bit-identity (14 checks)
│   ├── p4_transonic/       # G4.1-G4.3: subcritical bitwise no-op, upwind-reach evidence,
│   │                       #   coarse M0.80 shock quality vs reference band (10 checks)
│   ├── m0_meshgen/         # mesh gallery, hard-rule-7 topology matrix, cylinder convergence
│   ├── m1_wing_mesh/       # M6 wing+wake gallery, tip cut planes, ingestion/station/free-edge
│   │                       #   semantics, quality ladder, freestream-on-cut-mesh (13 checks)
│   └── ...                 # one directory per closed phase thereafter: p5-p11, p13, p14,
│                           #   m5, m6, b3-b9, b11-b18 (33 total). This tree is NOT the
│                           #   authoritative list -- cases/demo/README.md carries the full
│                           #   table with runtimes, and it is the one to update on close-out.
├── analysis/             # ✓ [Track A + phase studies] analysis/measurement working dirs
│   ├── README.md         #   (distinct from demo/: measurements and A/B studies, not
│   ├── a1_solver_bottleneck/   # [A1] 4-driver timing benchmark
│   ├── a2_te_kutta_fidelity/   # [A2] TE/Kutta fidelity attribution
│   ├── b9_fuselage_guardrail/  # [B9/GB9.6] fuselage-Cp guardrail
│   ├── b23_junction_discriminator/  # [B23] junction pocket = wake inboard free-edge singularity
│   ├── b24_wake_inboard_end/   # [B24] waterline-extension route CLOSED (negative)
│   ├── b25_inboard_fragment_clip/   # [B25] inboard_clip heals the junction pocket (14.66→0.63)
│   ├── b26_ls_transonic_ceiling/    # [B26] post-cure LS ceiling: medium 0.7625 / coarse 0.84
│   ├── b27_b18_demo_refresh/ # [B27] B18 demo refresh 8/8 PASS + 336/336 bit-identical
│   ├── b28_cl_fus_flat_sheet/ # [B28] cl_fus flat-fragment decoupling + GB9.4 re-spec
│   │                           #   (the "fuselage spurious lift" label retired; b9 demo 8/8)
│   ├── b30_transonic_ceiling/ # [B30] (b)-class ceiling attribution = SAME both paths
│   │                           #   (wing-tip P13 free-edge singularity + high-M Newton)
│   ├── b31_tip_termination/   # [B31] C-class wing-tip cure: conforming pressure+taper
│   │                           #   cures the 0.83 dying level; LS-side C1/C3 closed negative
│   ├── b32_tip_taper_adoption/ # [B32] conforming tip_taper adopted (medium ceiling 0.79→0.84);
│   │                           #   ② weld-sign per-step refresh rolled back (ill-posed)
│   ├── c1_ls_jacobian_fd/      # [A3/B19/B20/B21] LS-Jacobian FD probes, Leg-B density
│   │                           #   gap, GB20.7 recipe sweep + B21 N1 freeze-capture sweep
│   ├── p14_te_pressure_diag/   # [P14] TE pressure-Kutta diagnostics
│   ├── v1_ibl3_standalone/     # [V1/GV1.1] standalone IBL3 verification (prescribed u_e)
│   ├── v2_transpiration_channel/ # [V2/GV2.1] transpiration channel (δ*→ṁ) verification
│   ├── v3_fuselage_smoke/      # [V3/GV3.3] fuselage body-of-revolution smoke
│   ├── v3_loose_coupling/      # [V3/GV3.1/3.2] loose coupling, NACA0012 2.5-D strip
│   ├── v5_1b_scaled_newton/    # [V5/GV5.1b] scaled + damped augmented Newton
│   │                           #   (2P/0F/7R adjudicated (1P/1F/7R as executed): machinery
│   │                           #   exact, floor_reached stop works;
│   │                           #   window question reframed to an above-band seed)
│   ├── v5_1c_above_band_window/ # [V5/GV5.1c] the above-band seed: the pre-floor
│   │                           #   slope-2 window read (2P/1F/7R: NO quadratic regime
│   │                           #   above the floor — λ-capped halvings + a mid-range
│   │                           #   stall, never reaching the band)
│   ├── v5_1d_near_band_window/ # [V5/GV5.1d] the near-band seed: does a quadratic
│   │                           #   basin exist ADJACENT to the floor? (2P/1F/7R: NO —
│   │                           #   the stall extends down to 24× floor; medium's first
│   │                           #   step moves AWAY from the band)
│   ├── v5_5_te_floor/          # [V5/GV5.5] TE-band (B,δ) formulation floor-breaking
│   │                           #   (2P/1F/9R: the V1 TE-outflow row replacement does NOT
│   │                           #   break the floor — m2 5554×/245998× the floor, "worse"
│   │                           #   clause; damage peaks at the LE, flag default-OFF)
│   ├── v5_2_rae2822/           # [V5/GV5.2] RAE2822 transonic VII vs committed experiment
│   │                           #   (band (b) FAIL + the loose-recipe transonic-limit anatomy:
│   │                           #   every computed shock 0.06–0.10c aft of the bracket, 1/4
│   │                           #   legs converged; band (a) 9.46°/9.92° no fallback)
│   ├── v5_3_m6_cp/             # [V5/GV5.3] M6 wing direction+magnitude check vs the
│   │                           #   committed 7-station Cp (band (b) honest FAIL: the viscous
│   │                           #   Cp does NOT move toward experiment, 1/5 stations + pooled
│   │                           #   increase; band (a) RECORDED input-limited, Δcl_KJ −2.20 %
│   │                           #   DOWN under the A4 floor; k=0 anchored to P14 via the W1
│   │                           #   wiring guard)
│   ├── v5_4_cost/              # [V5/GV5.4] augmented-step cost on M6 medium (band (a)
│   │                           #   RECORDED 7.53× the inviscid step, above the ≤ ~2× band;
│   │                           #   band (b) honest FAIL: the block preconditioner does NOT
│   │                           #   work at medium — block-Jacobi diverges, exact-BL Schur
│   │                           #   stalls 1/4; "measure before Schur": the BL elimination
│   │                           #   is cheap (1.8 s), Krylov convergence is the bottleneck)
│   ├── v5_6_schur_prec/        # [V5/GV5.6] Schur-aware reduced-space preconditioner
│   │                           #   (the GV5.4 follow-up; band (b) honest FAIL: the
│   │                           #   corrected-AMG direction measured dead at wing scale —
│   │                           #   rung 3 rel_res 0.664 vs GV5.4's 2e-7..6e-5 stagnation;
│   │                           #   coarse rung 3 works 5/5; band (a) 7.87× RECORDED)
│   ├── v5_ibl_floor/           # [V5] IBL-floor diagnosis (GV5.1 follow-up, 14 RECORDED:
│   │                           #   raw cond mostly a scaling artifact + genuine scaled (A,Ψ)
│   │                           #   stiffness 1e5–1e7 + TE-band (B,δ) floor residual inside J's range)
│   ├── v5_m6_bridge/           # [V5/GV5.0] M6 subsonic loose-coupling bridge (RECORDED)
│   └── v5_tight_coupling/      # [V5/GV5.1] augmented-Newton exactness + convergence
│                               #   (9P/1F/36R: FD PASS both levels; quadratic tail blocked
│                               #   by the IBL floor) + the medium-seed diagnosis
│                           # (was missing from this tree until 2026-07-19 — the D9 finding)

tests/                     # Unit and gate tests
├── conftest.py           # ✓ Pytest fixtures: artifacts_dir (persistent, PYFP3D_ARTIFACTS_DIR
│                           #   overridable), mesh_dir, etc.
├── test_conftest_artifacts.py       # ✓ Regression test: gate artifacts persist in artifacts/
├── test_metrics_degenerate.py       # ✓ Regression test: degenerate-tet guard in metrics.py
├── mesh_utils.py         # ✓ [P1] Dependency-free structured-cube + sphere-shell mesh generators
├── __init__.py
├── test_v0_freestream.py # ✓ [P0/P1] Primary regression test (incl. cut-free residual check)
├── test_v1_surface_mesh.py        # ✓ [V1] surface-mesh DOF/basis/geometry + master-map hook
├── test_v1_closures.py            # ✓ [V1] closure FD-vs-analytic (both lanes), floors, seeds
├── test_v1_ibl3.py                # ✓ [V1] IBL3 Jacobian FD, bit-determinism, Newton laminar+turbulent
├── test_v2_transpiration.py       # ✓ [V2] transpiration assembly/divergence/u_e exactness,
│                                  #   GV2.1(a) coarse MMS lock, GV2.1(b) Picard legs
├── test_v2_newton_rhs_channel.py  # ✓ [V2] GV2.1(b)/(c): Newton external_rhs + LS b_base
│                                  #   bit-identity, Jacobian bit-invariance + FD under lagged ṁ
├── test_v3_coupling.py            # ✓ [V3] case-builder wiring (airfoil strip + closed body),
│                                  #   inflow/outflow pinning, 2-iteration coarse smoke
├── test_v5_above_band_seed.py         # ✓ [V5/GV5.1c] synthetic seed-helper tests (9):
│                                  #   perturbation mask + calibration bisection + triple
│                                  #   filter + regression slope + pooled verdict logic
├── test_v5_near_band_seed.py          # ✓ [V5/GV5.1d] synthetic near-band tests (7):
│                                  #   window sanity vs both floor bands + the GV5.1c
│                                  #   stall region / escalation direction / band-entry
│                                  #   read / near-band calibration / imported-helper id
├── test_v5_te_outflow.py          # ✓ [V5/GV5.5] TE-outflow row replacement tests (9):
│                                  #   default-OFF bitwise + residual/J row structure +
│                                  #   flag-ON FD + J_e zeroing + out-of-pattern guards +
│                                  #   plate smoke + te_outflow_pairs on the NACA strip
├── test_meshgen_structured.py     # ✓ [phase 3] structured/hybrid generator locks (11):
│                                  #   wall-anchored grading (append-only bit-identical prefix,
│                                  #   exact first step, radii-vs-distances kept SEPARATE on
│                                  #   purpose) + station count + LE-finer-than-TE + the local
│                                  #   window's two default bit-identities and its measured
│                                  #   outside cost + block quality drift lock (AR 5.899,
│                                  #   min_area 5.08e-06) + TE ray EXACTLY +-y.
│                                  #   ★ Added because the module had ZERO tests/ coverage: G0
│                                  #   lived only in a bench script, i.e. on no cadence
├── test_meshgen_rae2822.py        # ✓ [V5/GV5.2] RAE2822 point-set meshgen tests (5):
│                                  #   ordinate load (Cook layout) + PCHIP resample
│                                  #   convention/bounds/clustering + thickness/camber
│                                  #   signature locks + no-self-intersection + TE wedge
│                                  #   + the Cp-compare helpers on both committed layouts
├── test_v5_wing_case.py           # ✓ [V5] build_wing_case wiring on the M6 wall (LE/tip/root/TE
│                                  #   BC topology, local-x/c transition, scatter/gather + zero-RHS)
├── v5_state.py                    # ✓ [V5] shared GV5.1 builders: the 2.5-D NACA0012 strip case
│                                  #   + the loose-k1 state fixture behind all tight-gate tests
├── test_v5_tight_jacobian.py      # ✓ [V5] tight Stage 1: fixed operators + J_φ,BL FD gate
├── test_v5_tight_edge.py          # ✓ [V5] tight Stage 2: J_BL,φ = J_e·D_ue·G edge-chain FD gates
├── test_v5_tight_system.py        # ✓ [V5] tight Stage 3: full-system FD gate + smoke augmented
│                                  #   Newton (line-search probe guard exercised green)
├── test_v5_tight_scaled.py        # ✓ [V5/GV5.1b] scaled+damped path: 8 tests = scaling
│                                  #   identities + μ schedule + floor-stop + k1 smoke
├── test_v6_wake_sheet.py          # ✓ [V6/GV6.1+GV6.2] wake-sheet W3 construction +
│                                  #   producer identity (c) + zero-field (a)(i)
│                                  #   bit-identity + sign-pin MMS (b) + (a)(ii)
│                                  #   fresh-compile A/B vs the gate-free library +
│                                  #   fold-pairing structure + GV6.2 wake_l_rel_chords
│                                  #   plumbing (8)
├── test_mesh_*.py        # [P0] Gates G0.1–G0.4
├── test_mesh_adjacency.py           # ✓ [P0] Regression test for build_face_adjacency fix
├── test_mesh_reader_roundtrip.py    # ✓ [P0] Regression test for write_mesh tag-loss fix
├── test_laplace_mms.py              # ✓ [P1] Gate G1.1 -- PASSES
├── test_laplace_cg_iterations.py    # ✓ [P1] Gate G1.2 (formerly G1.3) -- PASSES
├── test_laplace_sphere.py           # ✓ [P1] Gate G1.6 (formerly G1.2) -- strict xfail, see "Known gaps"
├── test_laplace_picard.py           # ✓ [P1] Regression test for solve_laplace residual_norm fix
├── test_m0_extrude.py               # ✓ [M0] Prism-split unit tests (pure numpy, no Gmsh)
├── test_m0_cylinder.py              # ✓ [M0] Cylinder-flow validation (analytic Cp, spanwise)
├── test_wall_correction_cylinder.py # ✓ [P1] Gate G1.3 -- completed, acceptance NOT met
│                                     #   (negative result locked in; acceptance = strict xfail)
├── test_m0_naca0012.py              # ✓ [M0] NACA0012 family topology/wake-sheet/ingestion
├── test_p2_wake_cut.py              # ✓ [P2] Cut topology unit tests (synthetic strip, no Gmsh),
│                                     #   G2.1 + G2.2, assert-fires-on-broken-cut, hard-rule-7
│                                     #   sweep over every wake-tagged mesh in cases/meshes/
│                                     #   + [GV5.2] cambered-TE Kutta-probe bisector fallback
├── test_p2_kutta_naca0012.py        # ✓ [P2] Gates G2.3/G2.4/G2.5 + V2.1–V2.5 artifacts
├── test_m1_onera_m6.py              # ✓ [M1] M6 family: tags/geometry/wake-tip closure/quality,
│                                     #   swept-TE station + free-edge cut semantics, G2.1-style
│                                     #   freestream preservation on the cut coarse mesh
├── test_p3_assembly.py              # ✓ [P3] Colored-assembly rewrite: fast-vs-reference bit checks
├── test_p3_subsonic.py              # ✓ [P3] Gates G3.1 + G3.3 (incl. bit-identical Laplace limit)
├── test_p3_naca0012_m05.py          # ✓ [P3] Gate G3.2 (medium-mesh nested Picard, ~45 s)
├── test_p4_upwind.py                # ✓ [P4] Gate G4.2 (bitwise subcritical no-op) + upwind units
├── test_p4_transonic.py             # ✓ [P4] Gates G4.1/G4.3 (coarse smoke always-on; medium gate
│                                     #   + sweep behind PYFP3D_TRANSONIC_GATES=1)
├── test_p5_onera_m6.py              # ✓ [P5] 4 fast + 2 gated (G5.1/G5.2 behind
│                                     #   PYFP3D_TRANSONIC_GATES=1; polish recipe + 3% V6 bound)
├── test_p6_cp_metric.py             # ✓ [P6] shock-robust sign-alternating sawtooth metric
├── test_p6_recovery.py              # ✓ [P6] G6.1 recovery smoothing + G6.4 bit-identity
├── test_p6_weighted_flux.py         # ✓ [P6] opt-in kernel-flux invariants (no-op, determinism,
│                                     #   weighted=False restores the walk bitwise)
├── test_p7_diff_flux.py             # ✓ [P7] Gate G7.3: frozen-selection ∂ρ̃/∂φ of the walk flux
│                                     #   FD-verified (JVP vs shipped rho_tilde_sweep at frozen u;
│                                     #   all regimes + floor branch; kink-locus guard documented)
├── test_p8_jacobian.py              # ✓ [P8/N2] assembled Newton Jacobian JVP vs frozen-selection
│                                     #   residual FD (all regimes, rel ~1e-10 vs 1e-6 tol; kink
│                                     #   rows lifted from the P7 element guard); pattern-sharing,
│                                     #   limiter-mask gating, forward-path bit-guard; frozen-
│                                     #   assignment machinery (bitwise-at-freeze-state + frozen
│                                     #   JVP); + gated converged-pocket FD on the NEWTON coarse
│                                     #   M0.80 field (PYFP3D_TRANSONIC_GATES=1 — G8.1 FD clause)
├── test_p8_newton.py                # ✓ [P8/N3–N5] coupled Newton: Γ-column FD (far-field-column
                                      #   trap detector), exact Kutta row, far-field Γ-linearity,
                                      #   GMRES-vs-direct, supersonic nonsymmetry, cl/Γ match vs
                                      #   P3 Picard, terminal order p_k ~ 2, m_inf=0 single-step;
                                      #   + gated G8.1 terminal-quadratic runs (coarse M0.80,
                                      #   medium M0.7875 — re-specced case set, regression-lock
                                      #   physics bands; NEWTON_TRANSONIC_RECIPE lives here);
                                      #   ✓ [P8/N6] gated G8.2 M6 medium end-to-end < 300 s
                                      #   (NEWTON_M6_RECIPE lives here too; skips without the
                                      #   gitignored onera_m6/*.msh; carries the promoted
                                      #   G10.2 intermediate_tol=1e-5 since 2026-07-11)
├── test_p10_continuation.py         # ✓ [P10/G10.2] level-adaptive intermediate tolerance:
│                                     #   default-path accept_reason lock + subsonic-ramp
│                                     #   adaptive path (final level strict, Γ matches the
│                                     #   strict run to 1e-6, total steps not worse)
└── test_p13_tip_taper.py            # ✓ [P13/G13.2] spanwise loading taper (15): tip_taper=None
                                      #   is bit-identical; F is geometry-only; compact vs
                                      #   unbounded support; and the AMPLIFICATION law
                                      #   Γ/Γ* = F(1−b)/(1−F·b) with the P2 Kutta slope b≈0.93
                                      #   (F=0.8 ⇒ 0.21×, not 0.8× — the trap that makes r_c
                                      #   have to stay small)

(This tree covers P0–P13 only. The suite has **63 test files** as of
2026-07-19: also test_p14_te_pressure, test_p11_curved_walls,
test_a1_instrumentation, test_m2_wingbody, test_m5_round_tip,
test_b1..test_b19, test_b22_ls_3d_anchors
and the mesh/post unit files. `ls tests/test_*.py` is the authoritative list —
do not read a missing entry here as a missing test.)

artifacts/                 # Gate outputs (auto-generated, gitignored)
├── G0.1/                 # Volume conservation heatmap
├── G0.2/                 # Gradient recovery plots
├── G0.3/                 # Element coloring 3D render
└── ...

phases/                    # ★★ ARCHIVE of finished phases (reorganised 2026-08-10, base
│                           #   commit d224223). TRACKED, not gitignored -- the evidence stays
│                           #   in HEAD. Read phases/README.md first.
├── p1/                    #   phase 1: 24 demo + 23 analysis evidence chains, the phase-1
│                           #   roadmap/ and demo_report/, design_track_b.md, analysis/,
│                           #   archive/, and the 18 PURE LEVEL-SET test files
└── p2/                    #   phase 2: 45 bench scripts + s1_duct/, and 85 round files
                           #   ★ What stayed in place is what a surviving test or script
                           #   READS -- that rule is measurable, and it is why
                           #   cases/demo/p11_curved_walls/ did not move (the live G1.6 gate
                           #   reads it). Archived scripts had their repo-root depth, their
                           #   moved-target path literals and their gate_results reference
                           #   mechanically fixed and verified; they are provenance, and the
                           #   committed CSV/PNG is the evidence. For a guaranteed-runnable
                           #   phase-1/2 tree: git worktree add ../up3d-prereorg d224223

bench/                     # ✓ Measurement harness -- NOT tests. ★ The line between this and
│                           #   one-off studies was only DRAWN on 2026-08-10: before that,
│                           #   bench/ and cases/analysis/ were two homes for the same
│                           #   activity, split by which phase invented them. Merged here by
│                           #   user ruling. The criterion is NOT old-vs-new, it is WILL IT BE
│                           #   RUN AGAIN. See bench/README.md.
├── run_capability_locks.py  # ★ the FAST capability tier, ritual step 0: 5 groups, measured
│                           #   TWICE on 2026-08-11 at the same 8 threads on the same box:
│                           #   891 s and 564 s, both 5/5 green. Keep both -- a 1.6x spread with
│                           #   an identical result is the standing lesson that a wall clock is a
│                           #   CALIBRATION of the machine, not a reading of the solver.
│                           #   Pins its own thread caps (one lock asserts a wall clock) and
│                           #   prints WHAT IT DOES NOT COVER every run.
├── run_bench.py             # GS0.3 per-round drift detection vs baseline_2026-07-28.csv
├── run_m1_gate.py           # product metric M1, both seed legs; exit 1 while it FAILS
├── run_m3_budget.py         # M3's per-band error budget (LE 69.6 % / MID 20.7 %)
├── run_capability_matrix.py # 13 configurations x ladder = 78 points (M5-adjacent)
├── run_g82_anchor_check.py  # G8.2's physics anchors when the wall-clock assert masks them
├── run_le14_common_root.py  # the failure CLASSIFIER (never report a bare conv=False)
├── run_hex_g0_single_var.py # [phase 3 / task 2] G0: the 4 single-variable arms, bit-identity
├── run_hex_q2_wall_accuracy.py  # [task 2] Q2 wall Cp vs the analytic cylinder, wall DOF matched
├── run_hex_q1_donor.py      # [task 2] Q1 donor determinism D1-D4, global AND supersonic-zone;
│                           #   its premise asserts converged AND 0/0 clamps (the first run
│                           #   measured a DIVERGED field and would have published it)
├── run_task3_refinement_paradox.py  # [task 3] the 4 refinement legs + the alpha=0 discriminator
├── run_task3_seed_nonuniqueness.py  # ★ [task 3] the DECISIVE control: at a FIXED mesh, the seed
│                           #   alone moves x_shock 0.168-0.306 c = 31-56x M1's tolerance, so a
│                           #   refinement comparison at M0.80/alpha1.25 is unreadable. Also
│                           #   scopes it (M0.72 unique; alpha=0 nearly unique) and classifies
│                           #   every non-converged leg
├── run_task3_local_window_guard.py  # [task 3] the local-refinement window's guard
│   # ★★ the rest of task 3's harness, added here 2026-08-16 (29 run_task3_*/run_hex_* scripts and
│   #    32 committed CSVs in gate_results/ -- this document had listed only the first four):
├── run_task3_sigma_selection.py, run_task3_sigma_strength.py,
│   run_task3_sigma_charge_count.py, run_task3_sigma_freeze.py,
│   run_task3_soft_membership.py, run_task3_magnitude_preserving.py,
│   run_task3_capture_selection.py       # the seven sigma rounds; the TEMPORARY sigma_scale knob
│                           #   they used was DELETED at its registered expiry, with its tests
├── run_task3_m6_triage.py, run_m3_budget.py (extended)  # the staleness triage: the "LE band is
│                           #   70 % of M3's budget" number was EXPIRED and is retired
├── run_task3_le_registration.py, run_task3_tip_gate.py,
│   run_task3_le_offset_scatter.py       # the LE line: chordwise misregistration refuted, the tip
│                           #   refuted with a proxy-free gate, the error split into scatter+offset
├── run_task3_fixed_budget.py, run_task3_nearbody_axis.py  # allocation at a fixed cell budget:
│                           #   a real ~17 % lever that SATURATES; also the shared mesh builder,
│                           #   guards and search that every later task-3 script imports
├── run_task3_selfconvergence.py, run_task3_common_pointset.py,
│   run_task3_nonconv_discriminator.py, run_task3_tipstations_mechanism.py,
│   run_task3_shape_and_taper.py, run_task3_alignment_spanwise.py,
│   run_task3_nearsingular_partition.py  # ★ the self-convergence chain: d_self GROWS with
│                           #   refinement at the LE upper surface, which invalidated the
│                           #   extrapolation behind "changing the mesh cannot solve M3"
├── run_task3_s4_premise_check.py        # ★ the turn to product metrics: S4's two founding
│                           #   premises re-measured on HEAD (cost ratio, crossflow), plus the
│                           #   level trend behind ruling-grade decisions about GS4.1
├── run_gs31_precond.py      # reproduces the precond/EW decision, which is LIVE and carries
│                           #   refutation conditions -- testing them needs this script
├── run_le_response.py       # LE response measurement (imported by the above)
├── bitcheck.py              # DEVELOPMENT-time bit-identity A/B (ruling D1 -- not a test)
├── gate_results/            # ✓ committed CSVs for the recurring set; cited by path from the
│                           #   capability boundary. Archived scripts reach ACROSS to it.
└── studies/                 # ★ ONE-OFF registered studies, 20 gates, each self-contained:
                           #   PRE_REGISTRATION.md -> execution -> VERDICT.md -> results/.
                           #   The pre-registration is committed BEFORE the measurement it
                           #   governs; that ordering is the point. 23 phase-1 chains are in
                           #   phases/p1/cases/analysis/ -- the archive KEEPS the old layout
                           #   on purpose, since renaming a historical snapshot falsifies it.
                           #   ★ "20 gates" counts GATES, not directories (17 measured
                           #   2026-08-18) -- several studies carry more than one gate, so
                           #   the two numbers are not the same thing and do not reconcile.
                           #   ├── gs41_strip_core/  ★ [GS4.1 round 1, phase 4] the 2-D
                           #   │     strip core vs Blasius / Falkner-Skan. Its registration
                           #   │     lives in docs/dev_phase_four/ (prereg + 3 addenda),
                           #   │     not in the study dir -- phase 4 keeps one round file
                           #   │     per round there. V-FAIL, cause in the closure.

cases/                     # inputs and demos only, now that analysis/ merged into bench/
├── reference_data/        # ✓ external ground truth -- NEVER edit
├── meshes/                # ✓ mesh generators + committed mesh statistics. ★ Two level-set-only
│                           #   families cannot leave yet: four KEPT tests still read them, and
│                           #   they go when phase 3 amputates those LS legs.
└── demo/                  # ✓ per-phase self-checking demos with committed figures

docs/                      # only what phase 3 needs; the phase-1 plan/evidence is in phases/p1
├── dev_phase_two/          # ★ six files kept:
│   ├── PHASE_TWO_CAPABILITY_BOUNDARY.md   # ★★ READ FIRST when picking the project up
│   ├── LEVELSET_DELETION_INVENTORY.md     # ★ phase three's first task, counted by AST
│   ├── roadmap.md, progress.md            #   the plan + the 70-round ledger and its tables
│   ├── DECISION-2026-08-02-precond.md     #   a decision record with refutation conditions
│   └── _TEMPLATE.md         #   the round-file format (pre-registrations are committed BEFORE
│                           #   the measurement they govern -- that ordering is the point).
│                           #   The 85 executed round files are in phases/p2/docs/dev_phase_two/
├── dev_phase_three/        # ★★ THE ACTIVE PHASE (added here 2026-08-16 -- this document is the
│                           #   surface the close-out ritual flags as the one that silently rots,
│                           #   and it had carried only phase 3's FIRST round until now).
│                           #   68 round files, one pre-registration + one verdict per round, each
│                           #   pre-registration committed BEFORE the code it governs.
│   ├── progress.md          # ★ THE INDEX: one row per round (36 rows), plus the task ①②③ status
│   ├── 20260811-*           #   task ① delete the level-set route; task ② hexahedral body-fitted
│   │                        #   mesh (Q1 donor determinism, Q2 wall accuracy -- both FAIL, ruling
│   │                        #   D6 stopped the route)
│   ├── 20260811..0813-*     #   task ③ part 1: the refinement paradox, non-uniqueness, and seven
│   │                        #   sigma rounds (selection / strength / charge count / freeze / soft
│   │                        #   membership / magnitude-preserving / capture) -- the sigma line
│   │                        #   closed with a trade-off statement in the capability boundary
│   ├── 20260813..0815-*     #   task ③ part 2: the M6 LE line -- staleness triage, chordwise
│   │                        #   registration, tip gate, offset-vs-scatter, isoparametric-P2
│   │                        #   read-only assessment, fixed-budget allocation, self-convergence
│   ├── 20260815..0816-*     #   task ③ part 3: the non-convergence discriminator chain, then the
│   │                        #   turn to product metrics -- S4's premise check, its cost trend, and
│   │                        #   20260816-1000-gs41-initiation.md (a proposal, not an opened gate)
│   └── (the level-set deletion verdict is 20260811-0100)
├── design.md, design_track_v.md            # theory/numerics reference (current); track_b is
│                           #   archived with the level-set route (ruling D5)
├── agent-rules.md          # imported by CLAUDE.md via @; phase-one "Current phase" narrative
│                           #   is frozen, the session disciplines still apply
├── overview.md             # phase-one snapshot (frozen)
└── inspection/             # independent audits -- the 2026-07-28 full audit founded phase two

pyproject.toml            # ✓ Project metadata and dependencies
setup.py                  # ✓ Legacy setup (pyproject.toml preferred)
CLAUDE.md                 # ✓ Claude Code project instructions (doc map + workflow; imports docs/agent-rules.md)
```

## Implementation Status

### ✓ Complete (P0)
- **pyfp3d/physics/isentropic.py** — All physics scalars, numba-jitted
- **pyfp3d/mesh/reader.py** — Mesh I/O (meshio → SoA), mesh validation, tagged round-trip write
- **pyfp3d/mesh/metrics.py** — Geometry: volumes, gradients, face adjacency (Numba)
- **pyfp3d/mesh/coloring.py** — Element graph coloring for @prange
- **pyfp3d/post/vtk_out.py** — VTK writer (point fields) + PNG/CSV artifact helpers
- **tests/conftest.py** — Pytest fixtures
- **tests/test_v0_freestream.py** — Smoke tests + regression baseline ✓
- **tests/test_mesh_volume.py** — Gate G0.1 (volume conservation) ✓
- **tests/test_mesh_gradient.py** — Gate G0.2 (gradient recovery) ✓
- **tests/test_mesh_coloring.py** — Gate G0.3 (element coloring) ✓
- **tests/test_io_vtk.py** — Gate G0.4 (VTK round-trip) ✓
- **tests/test_mesh_adjacency.py**, **tests/test_mesh_reader_roundtrip.py** — regression tests for
  two bugs found by manual audit (see below)
- **pyproject.toml** — Build metadata and dependencies
- **CLAUDE.md** — Claude Code project instructions, auto-loaded each session (replaces the
  former `.copilot-instructions.md`, whose content largely duplicated design.md/roadmap.md
  and had drifted; details now live only in the authoritative docs)

Three latent bugs found by manual code audit (not caught by the existing suite, because nothing
exercised these code paths) have been fixed, each now with a regression test:
- `mesh/metrics.py::build_face_adjacency` crashed under `@njit` (reflected-list dict values are
  not valid in numba nopython mode) — rewritten around a `numba.typed.Dict` keyed by sorted face,
  storing only the packed first-owner index instead of a growing list. (`test_mesh_adjacency.py`)
- `mesh/reader.py::write_mesh` silently dropped every named boundary group (`wall`, `farfield`,
  ...), writing only a legacy `"all_triangles"` block, and `.msh` was ambiguous between meshio's
  `ansys`/`gmsh` writers (defaulted to `ansys`, discarding all tag data). Now writes each boundary
  group as its own tagged `triangle` block plus `gmsh:physical`/`field_data`, explicitly via the
  `gmsh22` writer. (`test_mesh_reader_roundtrip.py`)
- `solve/picard.py::solve_laplace` reported `residual_norm` over *all* nodes, including Dirichlet
  (far-field) rows whose natural-BC flux imbalance is O(1) and never shrinks — swamping the actual
  free-dof residual (which was already converging to ~1e-10). Now restricted to free dofs and
  correctly nets out `body_source_rhs` when present. (`test_laplace_picard.py`)

A second manual audit (2026-07-06) fixed four more latent issues, again each locked in by a
regression test:
- `tests/conftest.py::artifacts_dir` handed out a `tempfile.TemporaryDirectory`, so every gate
  PNG/CSV was deleted at teardown and the repo `artifacts/` directory stayed permanently empty —
  violating the CLAUDE.md/roadmap rule that every visual gate leaves inspectable headless
  artifacts. Now defaults to the persistent (gitignored) `artifacts/`, overridable via
  `PYFP3D_ARTIFACTS_DIR`. (`test_conftest_artifacts.py`)
- `mesh/metrics.py::element_gradients` silently returned zero gradients for |det J| < 1e-20 — an
  absolute threshold that both let coplanar tets corrupt assembly without any error and zeroed
  perfectly well-shaped tiny elements (edge ~1e-7 ⇒ det ~1e-21). Now raises `ValueError` with a
  scale-relative threshold (|det J| < 1e-12 × product of edge norms). (`test_metrics_degenerate.py`)
- `post/surface.py::_wall_vertex_normals` silently averaged cancelling normals when wall-triangle
  winding is inconsistent (or a crease is razor-sharp), producing garbage tangent planes
  downstream. Now raises with a diagnostic when |Σ area·n̂| / Σ area < 0.05 at any vertex.
  (`test_post_surface.py::test_inconsistent_wall_winding_raises`)
- `mesh/reader.py::read_mesh` rebuilt `tag_names` from `field_data` unconditionally (meshio meshes
  always *have* that attribute), clobbering the default `["bulk"]` with `[""]` for meshes carrying
  no named 3D groups. Now rebuilds only when a dim-3 entry exists.
  (`test_mesh_reader_roundtrip.py::test_untagged_volume_keeps_default_bulk_name`)

Numba pitfall hit while landing the degenerate-tet guard (captured here for the end-of-P1 skill
sedimentation checkpoint): `@njit(cache=True)` disk-cache invalidation tracks only the cached
function's *own* source, not its jitted callees — after editing `metrics.py::element_gradients`,
the stale cached `kernels/residual.py::assemble_residual` (with the old callee behavior baked in)
could still be loaded depending on compile order, making a test pass or fail depending on which
test file ran first. Symptom: behavior differs between `pytest tests/test_x.py` and the full
suite. Fix during development: delete `pyfp3d/**/__pycache__/*.nb[ic]` after editing any function
that other cached kernels call, then re-run.

Doc-only fixes in the same pass: design.md §5 far-field formula typo (`z sinα cosβ`, was garbled
"cosα-corrected"), design.md §10 V1 / roadmap G1.1 manufactured-solution wording (sin·sin·sin, not
sin·cos), `coloring.py` docstring (pure-Python preprocessing, not "@njit-compatible"; architecture
is design.md §7, not §3), `isentropic.py::validate_physics_bounds` docstring bounds aligned with
the implementation, a dead `nodal_gradient_recovery_spr` cross-reference in `post/surface.py`
replaced with the real wall-recovery functions, and stale `docs/PROJECT_STRUCTURE.md` paths (this
file lives at the repo root).

### ✓ P1 gates G1.1 and G1.2 closed; G1.3–G1.6 open (P1 renumbered 2026-07-06, see roadmap.md; below)
- **pyfp3d/kernels/residual.py** — Laplace residual (6.1) + SPD stiffness matrix (6.2) assembly
- **pyfp3d/solve/linear.py** — Dirichlet elimination (principal-submatrix reduction) + CG+PyAMG
- **pyfp3d/solve/picard.py** — `solve_laplace()` driver (P1's Picard loop degenerates to one solve)
- **pyfp3d/post/surface.py** — `nodal_gradient_recovery()` (volume-weighted, interior fields),
  `wall_tangential_gradient()` (surface-only linear recovery, for wall Cp), and
  `wall_tangential_gradient_quadratic()` (surface-only quadratic patch recovery; a real but modest
  improvement, see below)
- **tests/mesh_utils.py**, **cases/meshes/sphere_shell/** — structured-cube (MMS) and sphere-shell
  (G1.6) mesh generators; sphere-shell coarse/medium `.msh` + inspection PNGs are committed
- **tests/test_laplace_mms.py** — Gate G1.1 (MMS convergence) ✓ — L2 slope ≈ 1.94–1.96 with a
  sin·cos manufactured solution and a proper 4-point quadrature-consistent load vector. (A
  harmonic-polynomial exact solution was tried first and rejected: this codebase's structured
  Kuhn-triangulated cube reproduces harmonic quadratics to machine precision at *every* h, giving
  zero convergence-order signal — the same reason central finite differences are exact for
  quadratics.)
- **tests/test_laplace_cg_iterations.py** — Gate G1.2 (formerly G1.3; CG+AMG mesh-independence) ✓ — iterations
  8→11→14 across an 8×/level node-count increase (n=8,16,32 cube), comfortably under a 2× cap.

## Known gaps

*(This is the section CLAUDE.md, `docs/design.md` §5.1.2 and the four "see
'Known gaps'" pointers above all refer to. It carries the one long-standing
open defect — G1.6 — with its ruled-out fix routes, so nobody re-proposes
them. Added as a named heading in A3: the references had been pointing at an
unnamed block since the P1 renumbering.)*

**G1.6 (formerly G1.2; incompressible sphere Cp) is still open** — `tests/test_laplace_sphere.py::test_sphere_cp_medium_mesh`
is a `strict=True` xfail against the real <2% criterion, not a loosened threshold:
- The original `nodal_gradient_recovery` (volume-weighted average of the one-sided tets touching
  each wall node) gave ~26% max / 9% mean Cp error on the medium mesh — systematically low,
  because tangential velocity physically decays moving away from the wall and every incident tet
  sits on the inward side.
- A local least-squares ("SPR-style") extrapolation fix was tried and **rejected**: on a real
  graded mesh, some nodes' 1-ring element patches are nearly coplanar, and extrapolating a linear
  model through an ill-conditioned patch blew one node's Cp error up to 429 (found on the medium
  mesh). Don't resurrect this approach without a robust conditioning safeguard.
- The fix that shipped, `wall_tangential_gradient()`, computes the surface-tangential gradient
  directly from the wall's own P1 triangulation (physically exact for a natural-BC wall, where the
  normal component is zero by construction) instead of extrapolating the volume gradient. This
  improved medium-mesh max error to ~12%, and is a genuine, well-conditioned fix (no
  extrapolation, only interpolation/averaging) — but a mesh-refinement sweep (h_min = 0.08 → 0.04
  → 0.02 → 0.015, up to ~1.2M nodes / 7.4M tets) only reaches ~3.6% max error and visibly
  saturates rather than continuing to converge, so the remaining gap is *not* simply "refine the
  mesh more."
- **Root-caused (this session)** by isolating each candidate error source on a *clean*,
  single-variable h_min refinement sweep (h_min = 0.08 → 0.05 → 0.03 → 0.02, everything else in
  `generate_sphere_shell.py` held fixed — the earlier sweep above changed h_max/r_out/dist_max
  simultaneously with h_min, confounding the picture):
  - **Recovery scheme ruled out as the dominant cause.** An oracle test feeds the *exact* analytic
    potential straight into the recovery step, bypassing the FEM solve entirely. The recovery
    operator's own bias measured this way is 9.3%→0.7% (linear) or 0.5%→0.005% (quadratic, see
    below) across the same h_min range where the *full* FEM pipeline gives 12%→4.2% — i.e.
    recovery is a small fraction of the total error at every mesh size tested, for both schemes.
  - **Under-refined bulk/far mesh is a minor contributor, not the dominant one.** Tightening
    h_max/dist_max at a fixed h_min=0.03 helps a bit (5.6% → 4.3% max error) but plateaus by
    ~4% even with a much finer far mesh (0.5 vs 3.0 h_max, 4.5M vs 1.4M tets) — refining the bulk
    mesh alone doesn't close the gap either.
  - **Confirmed dominant cause: the volume PDE solve's own accuracy next to the wall**, not
    anything in post-processing. The raw nodal potential φ itself (not just its recovered
    gradient) has the same sub-first-order error at the wall as the derived Cp, and the
    convergence *order* measured on the clean sweep is decreasing as h shrinks (0.88 → 0.56 → 0.42
    for max nodal φ error) rather than settling to a fixed rate — ~~a signature of a genuine
    geometric/consistency error, not plain discretization error. Mechanism: the natural
    (zero-flux) BC is satisfied on the flat *polyhedral* wall-facet approximation (Γ_h), not the
    true curved sphere (Γ); solving `Δφ_h = 0` on the (slightly wrong) domain bounded by Γ_h
    pollutes the whole solution through ellipticity, not just the boundary nodes. This is a
    textbook "variational crime" for Neumann conditions on curved boundaries meshed with flat
    (non-isoparametric) elements.~~ ★★ **ERRATUM 2026-07-19 (P11, measured — see the P11 block
    below): the "PDE solve at the wall dominates" half of this bullet STANDS, but the
    variational-crime mechanism is OVERTURNED.** The decreasing order is the *fixed-bulk-mesh
    pollution floor* of a sweep that shrinks h_min while h_max stays 3.0 (the φ-error argmax at
    h_min=0.03 sits at r=1.53, in the coarsening transition zone, not at the wall; refining ONLY
    the far mesh drops wall φ error 3.17× and restores order 1.89), and a structured shell with
    the SAME flat facets converges at order ~2. The sweep numbers above are correct as measured —
    the inference drawn from them was the artifact.
  - **A direct fix was tried and rejected**: a Nitsche/penalty term added to the stiffness matrix,
    weakly forcing each wall-adjacent tet's own volumetric gradient toward zero along the *true*
    (here, analytically known) surface normal, swept over penalty strength β. Result: error and CG
    iteration count both got *worse* monotonically with increasing β (e.g. medium mesh max error
    12%→17%→40%→98%→211% for β=1→10→100→1000→1e4). Diagnosis: a P1 tet spanning from the wall
    inward necessarily has a nonzero radial gradient component representing the interior falloff
    of tangential velocity (the exact solution's normal derivative is zero only exactly *at* the
    wall, not throughout the adjacent tet's finite thickness) — that's correct FEM behavior, not a
    BC violation, so this penalty fights the physically-correct solution instead of correcting an
    inconsistency. Don't resurrect this approach; a correct fix needs to change how the boundary
    integral/geometry itself is represented (see below), not add a volumetric penalty.
  - **Implemented as a genuine (if modest) improvement**: `wall_tangential_gradient_quadratic()`
    fits a local quadratic model per wall node (in its own reconstructed tangent plane, over its
    1-ring — expanded to 2-ring, then falling back to a 2-parameter linear fit, if the patch is
    rank-deficient for the 6-parameter fit; uses `np.linalg.lstsq`'s SVD-based minimum-norm
    solution throughout, never a normal-equations solve or 3D extrapolation, so it can't hit the
    ill-conditioned blowup the earlier volume-based SPR attempt did). This is exact for a locally
    quadratic field (vs. linear recovery's exactness only for locally linear fields) and cuts the
    recovery-only oracle error by roughly 20x, but — consistent with recovery not being the
    dominant error source — only trims medium-mesh total error from ~12.0% to ~11.6%. Adopted as
    the default for the G1.6 test since it's a strict, low-risk improvement, but it does not (and
    was never going to) close the gate alone. See `tests/test_post_surface.py` for regression
    coverage locking in both facts (recovery-only accuracy, and the fact that it's still not
    enough).
- **Fix routes researched and tiered (2026-07-06, see design.md §5.1 for the full writeup)**:
  verify **Option A** first — the true-normal weak-flux correction (lagged/SBM-style: correct the
  natural-BC *data* via closest-point projection onto the true surface; RHS-only, stiffness matrix
  unchanged, AMG hierarchy fully reused). An oracle experiment (exact analytic potential + exact
  normals in the correction RHS, gate G1.4) measures the accuracy ceiling of this
  first-order normal correction, then the DP1 decision point picks the route:
  Option B (Gap-SBM gap terms) if the ceiling lands at 2–5%, or curved/isoparametric wall
  elements — now demoted to the option of last resort, taken up as its own scoped effort only if
  both A and B fall short — with G1.6 redefined per Option C's geometry-consistent-reference
  yardstick in that case. *(That last-resort route was executed 2026-07-19 as phase P11 and
  measured NEGATIVE — see the P11 block below.)* Still do **not** propose further h-refinement or recovery-scheme
  tweaks; both remain ruled out with evidence. `cases/meshes/cylinder_2.5d/` is confirmed to
  share the same root cause with quantified error (max |Cp err| 0.091 coarse → 0.045 medium,
  ~O(h)) and is the designated fix-route testbed — see design.md §5.1.1 / gate G1.3; the
  gate criterion itself stays on the sphere.

### ✓ G1.3 + G1.4 completed 2026-07-06 — negative results; DP1 "> 5%" branch taken

The Option A (true-normal weak-flux correction) verification chain ran to completion the same
day it was renumbered, and falsified the route (full evidence: roadmap G1.3/G1.4/DP1 entries,
design.md §5.1.2, `artifacts/G1.3/` + oracle results in `cases/demo/p1_laplace/results/`):

- **Delivered code**: `pyfp3d/solve/wall_correction.py` (closest_point_normal callback interface
  with analytic cylinder/sphere implementations, domain-outward facet orientation from the
  owning tet, 3-point edge-midpoint facet quadrature, RHS assembly verified against a
  hand-computed single-facet case); `pyfp3d/post/section_cut.py` (P2's final interface,
  degenerate single-layer path); `tests/test_wall_correction_cylinder.py` (10 tests);
  the sphere-oracle experiment (now absorbed into `cases/demo/p1_laplace/run_demo.py`,
  2026-07-07); cylinder `fine.msh` (50.2k tets); shared cylinder-case helpers moved into
  `tests/mesh_utils.py`.
- **Core finding**: for a harmonic potential with body-fitted wall vertices, the exact net flux
  through every flat facet is exactly zero (divergence theorem over the facet/true-surface
  sliver), so boundary-DATA corrections have (almost) nothing to correct. Measured: the full
  consistency defect is ~2e-5 (coarse cylinder, ~O(h⁴)); Option A's t-form assembles to machine
  zero on the cylinder; corrected Cp errors unchanged. Sphere ceiling: 0.1156 → 0.1133 (11.3%)
  vs the < 2% target → DP1 "> 5%" branch: Option C gate re-spec + separately-scoped curved
  elements. SBM/BDT-style corrections earn their keep on *unfitted* O(h)-gap boundaries, not
  here.
- **Cylinder re-framed**: ~76% of its Cp error at every level is the surface *recovery* on the
  quasi-2D sliver-strip wall triangulation (exact-potential-through-recovery oracle), and wall
  nodal φ converges at a healthy ~1.2 order — no sphere-style sub-first-order pathology. The
  "same variational crime" designation is retracted; the case remains the M0 end-to-end check.
- **G2.5(b) reference numbers** (by-product): solved-field spanwise floor max |w|/U∞ =
  2.88e-2 / 1.50e-2 / 7.40e-3 over coarse/medium/fine, identical with/without correction.

### ✓ P11 completed 2026-07-19 — curved wall elements NEGATIVE; G1.6 re-attributed

The DP1 "> 5%" branch's curved-element route was built and measured (phase P11, sphere leg;
full record: docs/roadmap/track_p.md §P11; demo `cases/demo/p11_curved_walls/`; tests
`tests/test_p11_curved_walls.py`; code `pyfp3d/solve/curved_wall.py` + the opt-in
`stiffness_delta` hook on `solve_laplace`, default bit-identical):

- **The route is a recorded NEGATIVE.** A verified curved wall-adjacent layer (quadratic
  tet10-style geometry via projected wall-edge midpoints, mapped-P1 field, ΔA delta assembly;
  quadrature == P1 reference to 1.3e-15 on straight geometry, planar-projection null test
  exactly zero, deg2-vs-deg3 A/B 5.5e-9) moves the medium sphere only **11.56% → 11.33%** — the
  same value as the G1.4 boundary-data oracle ceiling. The pre-registered superparametric risk
  fired: the mapped-P1 basis on quadratic geometry loses exact linear reproduction at **O(h)**
  (measured max 0.138 coarse), the same order as the O(h) facet-normal error it removes. Do not
  re-propose mapped-P1 (superparametric) curved wall elements.
- **The G1.6 re-attribution (two controls).** (i) A structured icosphere-extruded shell with
  the SAME flat-facet wall converges at wall-φ order **1.67/1.98** and reaches **2.14% max Cp
  at h≈0.036** — the "crime" geometry converges fine when the mesh is structured and all scales
  refine. (ii) The clean h_min sweep's order collapse (0.88/0.56/0.42, replicated exactly by
  the committed demo script) is the **fixed-bulk-mesh pollution floor**: at h_min=0.03 the
  φ-error argmax sits at r=1.53 (the coarsening transition zone); refining ONLY the far mesh
  (h_max 3.0→1.0) drops wall φ error **3.17×** and restores order **1.89** — no curved anything.
  ⇒ **the medium mesh's 11.6% is essentially the intrinsic P1-field max-norm capability at
  h=0.08** (structured control interpolates to ≈11% there); the wall geometric crime contributes
  ≈0.2 pp; the recovery ≈0.2–0.5 pp (unchanged oracle result).
- **What stands / what fell.** Stands: recovery not dominant; Nitsche dead; boundary-data
  corrections have nothing to correct; the gate is still unmet on the committed medium mesh
  (`test_laplace_sphere.py` strict xfail unchanged). Fell: "the dominant error is the
  flat-facet natural-BC variational crime" and "closing G1.6 needs curved/isoparametric
  elements".
- **Route fork RESOLVED 2026-07-22 (user-directed): (a) Option C re-spec ADOPTED.** The active
  G1.6 gate is now the achievable, measured criterion asserted PASSING by
  `tests/test_laplace_sphere.py::TestG16Respec` (reads P11's committed sweep, no re-solve):
  **all-scales-refined φ_w order ≥ 1.8** (E6 icosphere s4→s5 = 1.98; E8 far-refinement = 1.89)
  **+ mean Cp < 1% at h_min 0.03** (E8 h03_far10 = 0.60%). The literal 2%-max-at-medium
  `test_sphere_cp_medium_mesh` **STAYS a strict xfail = the recorded P1 limitation** (it demands
  O(h²) wall velocity at h=0.08, beyond any P1-field method on any mesh). Unchosen routes on
  record: (b) isoparametric P2 wall layer (new midside DOFs, field order raised — the only route
  to the LITERAL criterion, and the one that would also tighten Track V's u_e input band, A4);
  (c) accept as a permanent recorded limitation. The wing-body "G1.6-class" cross-attribution
  (GB9.4/GB20.5) lost its sphere anchor and needs its own discriminator before being quoted again.

### ✓ M0 mesh-side items delivered (2026-07-06)

- **pyfp3d/meshgen/extrude.py** — single-layer quasi-2D extrusion: one cell layer in z, both
  planes tagged `symmetry`; every prism split into 3 tets with the globally consistent
  min-global-index diagonal rule (Dompierre et al. "indirect" subdivision — since top-layer
  indices are bottom + n_2d, the smallest index on any lateral quad is always the smaller
  *bottom* node, so both sides of every shared quad pick the same diagonal by construction);
  tagged interior edge sheets (the wake) split by the same rule so they coincide exactly with
  tet faces. `assert_quad_split_consistency()` is the M0 preprocessor assert: single-owner tet
  faces must equal the tagged boundary set, interior-sheet faces must have exactly two owner
  tets; unit tests prove it fires on a deliberately broken split (`tests/test_m0_extrude.py`).
- **pyfp3d/meshgen/planar.py** — vanilla-Gmsh 2D builders (Distance+Threshold grading like the
  sphere case): `cylinder_annulus_2d()` and `naca0012_wake_2d()` (closed-TE NACA0012, circular
  far field r=15c centered at mid-chord, wake line TE→farfield embedded with
  `gmsh.model.mesh.embed` so triangle edges conform to it). Gmsh imported lazily.
- **cases/meshes/naca0012_2.5d/** (the M0 deliverable) — `generate_naca0012.py`, one parameter
  (h_wall) per level, everything else derived: coarse 16.4k / medium 61.8k tets committed
  (on-target vs the ~15k/60k spec; fine ~240k on demand). Per-level stats CSV + layer-inspection
  PNG written at generation time (headless, roadmap §0.1).
- **cases/meshes/cylinder_2.5d/** (extra validation case) — quasi-2D circular-cylinder flow,
  no wake (non-lifting), analytic solution phi = x(1 + a²/r²), Cp = 1 − 4 sin²θ. The simplest
  end-to-end check of pipeline + P1 solver: measured max |Cp err| 9.1% (coarse, 6.9k tets) →
  4.5% (medium, 17.3k tets) with quadratic surface recovery — same curved-wall/flat-facet
  variational crime as the G1.6 sphere, converging ~O(h), fully expected on this geometry.
- **tests/test_m0_extrude.py / test_m0_cylinder.py / test_m0_naca0012.py** — 21 tests covering
  the M0 gate items: reader ingestion, tags, quad-split consistency, wake-sheet topology (one
  connected planar interior sheet TE→farfield, both z-planes, nodes not duplicated), symmetry
  planes planar/disjoint from wall, cylinder Cp vs analytic, spanwise-gradient behavior.

**G2.5(b) finding (evidence in roadmap.md G2.5 note — the gate criterion needs re-spec before
P2 closes it):** the *interpolated* freestream phi = x is spanwise-zero to machine precision
(G2.5(a) verified), but the *solved* field is not and cannot be: a 3-tet prism split is
necessarily asymmetric under the z-mirror (on a lateral quad face, ∫N_i dS is S/3 or S/6
depending on the diagonal direction), so a z-invariant field does not satisfy the discrete
equations and the discrete minimizer carries O(h) spanwise velocity noise. Measured (non-lifting
solves): cylinder max|w|/U∞ 2.9e-2 (coarse) → 1.5e-2 (medium); NACA0012 α=0 5.4e-2 (coarse) —
clean ~O(h) decay, 10 orders above the literal 1e-12 criterion. Literal machine-zero would need
a z-mirror-symmetric subdivision (requires Steiner points, violating the 3-tet M0 spec). The
tests assert the honest behavior instead: machine-zero for the interpolant, small-and-decreasing
under refinement for solved fields.

### ✓ P2 delivered, M0 closed (2026-07-06)

- **pyfp3d/mesh/wake_cut.py** — wake-sheet node duplication: ⁺-side (upper-hint) copies
  appended after the original nodes, so the reduced dof space of the wake constraint is
  exactly the original node set; per-node flood-fill side classification (adjacency through
  non-wake faces only — no planarity assumption, ready for the M1 swept wake); spanwise
  stations group TE nodes by (x, y), so a quasi-2D extrusion collapses to the single scalar
  Γ of the M0 spec while a swept TE gets per-node stations; Kutta probe nodes (wall, one
  edge off each TE node, per side); preprocess-time topology asserts (roadmap P2 list,
  proven to fire on a deliberately broken cut).
- **THE P2 spec deviation (evidence-backed): TE nodes ARE duplicated.** The roadmap
  originally asserted the opposite; implemented first and measured: a single-valued TE node
  tapers [φ] from Γ to 0 across the first wake cell ≡ a point vortex of strength Γ parked
  at the TE → wall suction ~ (Γ/2πr)², a spurious TE force ~ Γ²/h that *diverges* under
  refinement (coarse NACA0012, Γ=0.3 prescribed: peak wall |V| 4.6 U∞, −0.27 spurious out
  of cl ≈ 0.6, from 6 triangles). With the TE doubled: cl = 0.6012 vs Kutta–Joukowski 0.6.
  Roadmap assert block re-specced accordingly; design.md §4 records the theory (the TE jump
  IS the Kutta condition).
- **pyfp3d/constraints/wake.py** — master–slave elimination: φ_full = T φ_red + g(Γ),
  A_red = TᵀAT assembled once (SPD preserved), Γ RHS-only via precomputed per-station
  vectors h_j = TᵀA g_j; folding slave rows into masters is exactly the weak mass-flux
  continuity (4.2). `kutta_targets()` = per-station mean of probe jumps (also filters the
  O(h) spanwise noise on quasi-2D meshes).
- **pyfp3d/constraints/dirichlet.py** — far-field freestream at incidence + incompressible
  2D point-vortex correction with the branch cut ON the wake sheet: an eliminated ⁺-side
  far-field wake node (master + Γ) automatically equals the upper-branch vortex value, so
  no special Dirichlet casing is needed.
- **solve/picard.py::solve_laplace_lifting()** — Kutta outer loop with the matrix, Dirichlet
  split and AMG hierarchy built once; Γ updates are RHS-only. Secant (Aitken) acceleration
  from the second update: the linear map's measured slope is b ≈ 0.93, so plain ω-relaxation
  would need O(100) updates; the secant hits the affine fixed point in 2.
- **post/surface.py** — triangle-wise wall force integration (`wall_force_coefficients`):
  in-plane tangential gradient IS the wall velocity (natural BC), Cp per triangle, outward
  normals oriented by the owning tet (no winding/star-shape assumptions), no nodal averaging
  across the sharp-TE crease; `sectional_cl_from_gamma` (KJ cross-check).
- **post/section_cut.py** — general z = const marching-tets path (linear-exact, unit-tested
  at an off-node plane) + `wall_cp_curve()`: triangle-wise sectional Cp(x/c) split
  upper/lower at the intersection-segment midpoints.
- **cases/reference_data/naca0012_incompressible/** — Hess–Smith panel reference (constant
  sources + single vortex + Kutta), same closed-TE coordinate set as the mesh so G2.3 is
  method-vs-method; cl(4°) = 0.482556 at N=800, Cp-integration vs Kutta–Joukowski agree to
  0.09%, lift slope 6.91/rad vs thickness-corrected 6.90. Provenance in its README.
- **Gates** (all green, `tests/test_p2_wake_cut.py` + `tests/test_p2_kutta_naca0012.py`,
  artifacts/G2.{1..5}/): G2.1 ‖R‖∞ = 8.4e-13 (folded wake-master rows 6.9e-16); G2.2
  [φ] − Γ < 1e-13; G2.3 medium cl_p = 0.47858 → −0.82% vs panel (coarse −3.0%), Kutta
  converged in 2 updates; G2.4 Γ-cl vs pressure-cl 0.01% (coarse 0.32%); G2.5 closed under
  the re-specced criterion (b) — p99 |w|/U∞ 4.82e-3 → 2.35e-3 (ratio 2.05 at h ratio 2),
  max recorded (LE peak-gradient region, not wake), stripe-free mid-plane heatmap.
- **M0 closed** with it: wake-cut topology asserts sweep every wake-tagged mesh in
  cases/meshes/ (hard rule 7 test) and the G2.5 acceptance link is green.

### ⛔ Track B — level-set embedded wake: **DELETED 2026-08-11** (ruling D5)

★★ **This whole section describes a route that no longer exists in the tree.** Phase 3's
first task deleted it: 9 library files / **4624 lines**, `pyfp3d/wake/` removed entirely,
`post/unified.py` collapsed onto its conforming half. `from pyfp3d.wake import ...` is a
hard error now. The text below is kept as the DESIGN RECORD of what was built and why --
read it in the past tense; the code is in `phases/p1/` and in git history, and the
measured reasons for abandoning it are in docs/dev_phase_two/roadmap.md ruling D5.
Verdict: docs/dev_phase_three/20260811-0100-ls-deletion-verdict.md.

Original heading: *✓ Track B — level-set embedded wake (B1 ✓ B2 ✓ B3 ✓ B4 ✓ B5 ✓ B7 ✓;
B6 ◐ in progress)*

A **parallel** wake representation: instead of duplicating nodes so the mesh conforms to the
wake sheet (`mesh/wake_cut.py` + `constraints/wake.py`), the sheet is a **level set** and the
elements it cuts carry two DOF copies. Purpose (user-arbitrated) is **mesh/geometry workflow
capability**, not solver speed: the mesh need not know the wake exists, so α can be re-aimed
without remeshing. The conforming path stays **byte-untouched** — nothing in `wake/`,
`kernels/cut_assembly.py`, `solve/*_ls.py` or `post/surface_ls.py` is imported by it.

The two structural payoffs, both delivered:
- **The Kutta condition is IMPLICIT** — no Γ secant, no master–slave Γ. Γ is a *solution mode*
  read off the converged TE jump. ★ Consequence measured in B6/B7: since there is no
  early-stoppable Γ outer loop, the level-set **Picard** tracks the conforming **Newton** truth
  to within a few % (M6 coarse: +0.7% on the wake-free mesh), while the conforming *Picard*
  under-circulates ~8% below it.
- **Γ(tip) = 0 falls out discretely** from the level set's spanwise clip, with no free-edge
  bookkeeping (B7: tip Γ ~3e-4 on ONERA M6).

Two findings worth knowing before touching this code (both cost real time to discover):
1. **The wake level set CANNOT pin Γ** (B4): its residual is identically zero for any spatially
   constant jump, because Σ_c ∇N_c = ∇(1) = 0 (partition of unity; measured 1.9e-16). Γ needs
   its own condition — the **nonlinear TE pressure-equality Kutta**, recovered on **wall-adjacent**
   control volumes (the full element fan gives Γ +11–15%).
2. **The conforming transonic recipe does not transplant** (B6): the P4 whole-field θ·diag
   damping is a Jacobi smoother, so it throttles the (now smooth, global) circulation mode ⇒
   damping must be localized to the supersonic rows; and near the fold the live Γ→far-field-vortex
   loop has **gain > 1** ⇒ the transonic/3D recipe is `farfield="neumann"` (the López outlet).

Authoritative: [docs/design_track_b.md](phases/p1/docs/design_track_b.md) (numerics; §11 = the B7 3D gate)
and [docs/roadmap.md](phases/p1/docs/roadmap.md) Track B (gates + ledger). Evidence: demos
`cases/demo/b3_levelset_lifting/`, `b4p5_farfield/`, `b6_transonic/`, `b7_onera_m6/`; tests
`test_b1_cut_elements` / `test_b2_multivalued` / `test_b3_lifting` / `test_b4_te_control_volume` /
`test_b45_farfield` / `test_b6_transonic` / `test_b6_newton` / `test_b7_onera_m6`.
(Note the `b4p5_farfield` / `test_b45_farfield` names predate the 2026-07-12 Track-B renumber and
are kept on purpose so the committed paths stay stable — that gate is now **B5**.)

### ⏳ Next

> **This "Implementation Status" section is a P0–P2-era historical record and is NOT the
> tracker.** "What phase are we in / what gate is open" lives ONLY in
> [docs/roadmap.md](phases/p1/docs/roadmap.md) (progress ledger) and [docs/agent-rules.md](docs/agent-rules.md)
> ("Current phase"); the per-phase evidence lives in [docs/demo_report.md](phases/p1/docs/demo_report.md).
> P3–P10 and Track B/M all closed gates *after* the text above was written — read the roadmap,
> not this list.

- **Track B → B8** (level-set tip-edge desingularization / row-blend tip taper — NEW
  2026-07-13, the LS analogue of P13/G13.2) and **B9** (multi-wake: multi-element /
  wing–body) are the open Track B gates (Track-B renumber 2026-07-13: new B8 inserted,
  old B8 multi-wake → B9, old B9 curved wake → B10);
  **B6** (transonic level-set) stays ◐ open on its medium quantitative closure.
- **G1.6 re-spec per Option C — DONE 2026-07-22 (user-directed):** the achievable, measured
  acceptance criterion (all-scales-refined φ_w order ≥1.8 + mean Cp <1% at h_min 0.03, on P11's
  E6/E8 committed sweep — the geometry-consistent all-scales reference) is asserted PASSING by
  `tests/test_laplace_sphere.py::TestG16Respec`; the literal 2%-max xfail stays = recorded P1
  limit. See "Known gaps": h-refinement, recovery tweaks, Nitsche and boundary-data corrections
  are all **ruled out with evidence** — do not re-propose them.
- ~~G1.3/G1.4 oracle experiments~~ — DONE 2026-07-06 with negative results (see the G1.3+G1.4
  section above); DP1 decided the "> 5%" branch.
- ~~P3 (subsonic compressible)~~ — long since closed (P3–P9 closed; P10 partial). Ignore the
  stale entry that used to sit here.

## Quick Start

### 1. Install dependencies
```bash
cd <repo-root>   # e.g. ~/code/UP3D
pip install -e ".[dev]"
```

### 2. Run smoke tests
```bash
pytest tests/test_v0_freestream.py -xvs
```

Expected output:
```
tests/test_v0_freestream.py::test_import_pyfp3d PASSED
tests/test_v0_freestream.py::test_import_physics PASSED
tests/test_v0_freestream.py::test_isentropic_stagnation PASSED
tests/test_v0_freestream.py::test_isentropic_freestream PASSED
tests/test_v0_freestream.py::test_pressure_coefficient_bounds PASSED
```

### 3. Run physics module directly (self-test)
```bash
python -m pyfp3d.physics.isentropic
```

Expected output:
```
=== Isentropic Physics Self-Test ===

Freestream Mach M∞ = 0.5

At stagnation (q² = 0):
  ρ = 1.000000 (expected 1.0)
  Cp = 1.000000 (expected 1.0)

At freestream (q² = 1):
  ρ = 1.000000
  M = 0.500000 (expected 0.5)
  Cp = 0.000000 (expected 0.0)

Critical speed q*² = 0.923077 where M = 1.0:
  Computed M = 1.000000 (expected 1.0)

✓ All checks passed!
```

## Next Steps (P1 Completion)

P0 (mesh I/O, metrics, coloring, VTK writer, gates G0.1–G0.4) is done. G1.1 and G1.2 (formerly
G1.3) are done; G1.3 and G1.4 completed 2026-07-06 with negative results and DP1 took the
"> 5%" branch (see above). **G1.6 Option C re-spec DONE 2026-07-22 (user-directed): the active
gate is the achievable measured criterion (`TestG16Respec` PASS); the literal 2%-max stays xfail
= recorded P1 limit.** Historical planning notes below kept for the record:

1. ~~**Run the G1.3 cylinder oracle pre-study and the G1.4 sphere oracle experiment, then follow
   the DP1 decision point**~~ — DONE 2026-07-06, negative result: the Option A ceiling is
   ≈ 11.3% (medium sphere) because on body-fitted meshes the boundary-data defect it corrects is
   (near-)zero — see the "G1.3 + G1.4 completed" section above and design.md §5.1.2. DP1 took
   the "> 5%" branch. ~~The new open item: **draft the G1.6 Option C acceptance re-spec**
   (geometry-consistent reference on the same polyhedral domain)~~ — DONE 2026-07-22
   (`TestG16Respec`, all-scales order ≥1.8 + mean-Cp <1% at h_min 0.03); the curved/isoparametric-
   element effort for physical accuracy is separately scoped (P11 measured superparametric NEGATIVE).

2. ~~**Create test meshes** (M0, parallel track)~~ — DONE 2026-07-06 (single-layer re-spec
   targets ~15k/60k/240k tets, not the older 30k/150k/700k): `pyfp3d/meshgen/` +
   `cases/meshes/naca0012_2.5d/` + `cases/meshes/cylinder_2.5d/`, topology validated by
   `tests/test_m0_*.py`. See "M0 mesh-side items delivered" above.

## References

- **Design & Theory:** [docs/design.md](docs/design.md)
- **Roadmap & Gates:** [docs/roadmap.md](phases/p1/docs/roadmap.md)
- **Agent Rules:** [docs/agent-rules.md](docs/agent-rules.md)
- **Claude Code Instructions:** [CLAUDE.md](CLAUDE.md)

---

**Last updated:** 2026-07-18  
**Status:** per-track status lives in [docs/overview.md](docs/overview.md)
(human-readable snapshot) and the per-track trackers
[docs/roadmap/](phases/p1/docs/roadmap) (authoritative; docs were split by track
2026-07-15 — docs/roadmap.md and docs/demo_report.md are now thin indexes).
One-line summary: Track P — P0–P9 ✓ (P1: G1.6 open as a `strict=True` xfail;
root cause RE-ATTRIBUTED by P11, see "Known gaps"), P10 ◐, **P11 ✓ CLOSED
(2026-07-19, opened + closed same day): curved wall elements measured
NEGATIVE (medium 11.56%→11.33% = the oracle ceiling; superparametric O(h)
risk fired); G1.6 = intrinsic P1 capability at h=0.08, the order collapse was
a confounded-sweep bulk floor; **route fork RESOLVED 2026-07-22 — Option C
re-spec ADOPTED** (achievable criterion `TestG16Respec` PASS; literal 2%-max
stays xfail = recorded P1 limit)**, P13 ◐ (G13.3 transonic NEGATIVE-open), **P14 ✓ CLOSED
(2026-07-17, opened + closed same day): pressure-equality Kutta estimator;
G14.1–G14.7 ✓; the conforming path now MATCHES level-set (cl_p/cl_KJ
0.15%/0.34%), and the +4.85% cl_KJ move off the probe locks closed 69% of
P9's 0.019 gap**; Track M — M0–M5 ✓ (M2 ✓ — its solver leg was closed by B9
on 2026-07-17, both wake models now run on the wing-body); Track B —
B1–B9, B11–B32 ✓, B6 ◐, B10 shelved (**B16/B17 far-field aux pin +
`pin_gamma`, B18 wing-body transonic, B19 LS-Jacobian exactness**, all
2026-07-18; **B20 mixed-plain main-field density ADOPTED PERMANENTLY +
re-baselined, B21 N1 freeze-capture fix restoring the M6-medium M0.84 ramp —
GB20.7's "real capability loss" verdict overturned — and B22 evidence
refresh (B15 demo 20/20, B14 7/7) + gated 3-D anchor locks closing the N3
gap**, 2026-07-19; **B23 junction discriminator (pocket = wake inboard
free-edge singularity, NOT faceted geometry) + B24 waterline extension
(route closed) + B25 `inboard_clip` — the junction pocket HEALED**,
2026-07-19; **B26 (pocket-healed LS ceiling = the conforming ceiling site:
coarse 0.84 reached / medium 0.7625) + B27 B18 demo refresh, 8/8 PASS — the
B18 "LS junction-limited" story RETIRED**, 2026-07-20; **B28 cl_fus decoupling
+ GB9.4 re-spec + B29 flat-fragment adopted as the wing-body LS production
config**, 2026-07-20; **B30 (b)-class ceiling attribution (SAME mechanism both
paths = wing-tip P13 free-edge singularity + high-M Newton) + B31 C-class
wing-tip cure (conforming pressure+taper CURES it; LS-side closed negative) +
B32 conforming tip_taper adopted — wing-body medium ceiling M0.79 → M0.84
reached**, 2026-07-21/22); Track V —
**V1 ✓ CLOSED 2026-07-22 · GV1.1 9 PASS / 2 FAIL** (IBL3 solver core
shipped: `pyfp3d/viscous/` surface_mesh/closures/ibl3, standalone prescribed-u_e;
VERDICT `bench/studies/v1_ibl3_standalone/VERDICT.md` — (a) ×2 = closure-family
fixed point, (e) first-run outflow 2h grid mode FAIL → fixed by the D-HB
streamwise-tensor stabilization ε_s=0.02 = PASS,
(b)(c)(d) PASS); **V2 ✓ CLOSED 2026-07-22 · GV2.1 23 PASS / 0 FAIL /
16 RECORDED** (transpiration coupling shipped: `pyfp3d/viscous/transpiration.py`
δ*→ṁ operator + wall-RHS channels in solve_laplace / solve_subsonic(+lifting) /
newton_lifting / ls_newton — `None` ⇒ legacy path bit-identical; VERDICT
`bench/studies/v2_transpiration_channel/VERDICT.md` — (a) MMS cylinder-blowing
convergence strict-decreasing, order 1.65/1.64 ≥ 1.0, (b) five-driver ṁ=0
bit-identity, (c) FD Jacobian 6.6e-09–7.2e-08 < 1e-5); **V3 ✓ CLOSED
2026-07-22 · GV3.1/GV3.2 2 PASS / 4 FAIL / 23 RECORDED · GV3.3 0 PASS /
2 FAIL / 7 RECORDED** (loose coupling shipped: `pyfp3d/viscous/coupling.py`
+ committed XFOIL reference `cases/reference_data/naca0012_viscous_xfoil/`
+ BoR smoke generator `cases/meshes/fuselage_bor/`; IBL3 local-basis
crossflow fix 25.9/0.15 → 1.8e-4/1.6e-3 en route; VERDICTs
`bench/studies/v3_loose_coupling/` + `phases/p1/cases/analysis/v3_fuselage_smoke/` —
Δcl PASS ratio 0.542 ∈ [0.5, 2.0], loose loop 4–5 outer iters at ω = 1.0
incl. transonic M 0.72 record (4 iters, no tuning); honest FAILs localized:
cf +44 % at the first post-trip station only (XFOIL e^N ramp vs
instantaneous switch), δ* H-family offset ≤ 27.9 % at x/c = 0.074, GV3.3
tail-cone σ/μ 0.5533 / crossflow 0.2631 FAIL + loop NOT converged =
measured stern instability — V4 skip criterion met by letter (GV3.2),
counter-evidence logged (GV3.3) — **V4 ⊘ SKIPPED 2026-07-22** (user-directed;
reopen trigger = V5 stall / pre-V5 closed-body scope)); **V5 ✓ CLOSED
2026-07-25 · GV5.0 ✓ EXECUTED (16 RECORDED / 0 FAIL)** (M6 subsonic
loose-coupling bridge, VERDICT `bench/studies/v5_m6_bridge/`: bridge answer
= the loose loop is NOT sufficient on the 3-D lifting wing — coarse
root-upper-TE separation-patch runaway ṁ_max ×12.4 (GV3.3-stern class),
medium patch refined away but bounded δ* limit cycle 2–12 %/k, tol 1e-3
never met; ΔCL DOWN both estimators (medium −2.4 % input-limited); crossflow
small max|B|/|A| ≤ 0.072; tip mask validated; `build_wing_case` +
`tests/test_v5_wing_case.py` (5) new; δ*(z) CSVs feed GV5.3's bands; medium
wall-time polluted by external load, quoted flagged) · **GV5.1 ✓ EXECUTED
2026-07-23 (9 PASS / 1 FAIL / 36 RECORDED)** (augmented tight (φ, Γ, BL)
Newton shipped: `pyfp3d/viscous/tight.py` + `tight_driver.py`,
`tests/v5_state.py` + 3 tight test files; VERDICT
`bench/studies/v5_tight_coupling/VERDICT.md`, design record
`docs/design_track_v.md` §12 — band (a) FD exactness PASS both levels,
worst sweet-spot 2.2e-8 coarse / 5.1e-9 medium; band (b) quadratic tail
HONEST FAIL = the intrinsic floor of the steady IBL residual on the
cond(J_BL,BL) ~ 4e10 near-null manifold (the standalone pseudo-time solve
stalls there too), NOT a coupling defect; band (c) N_aug ≤ 2 not met
standalone nor as polish, N_total 14/13 vs loose 4/5; finding: the
committed GV3.1 medium fixed point is NOT reproducible — IBL-floor
trajectory scatter, diagnosis committed, HEAD-regen seed user-accepted;
IBL-floor follow-up diagnosis ✓ EXECUTED 2026-07-24 (14 RECORDED,
`bench/studies/v5_ibl_floor/`: raw cond 4e10–4e13 mostly a scaling
artifact (equilibrated 2e4/7e5/1e7, sub-1e-6 → 0/0/2, no exact null
directions), genuine scaled (A, Ψ) stiffness 1e5–1e7 remains; the floor
residual lives at the TE band (B, δ) equations inside J's range; closure
floors inactive; eps_diff ×4 ≤ 6 %; the pseudo-time controller bottoms
out = a formulation floor globalization alone cannot pass); **GV5.1b ✓
EXECUTED 2026-07-24 (2 PASS / 0 FAIL / 7 RECORDED adjudicated
2026-07-24; 1P/1F/7R as executed, preserved in commit 1c55906,
`bench/studies/v5_1b_scaled_newton/`, design record
`docs/design_track_v.md` §14)**: the scaled + damped machinery is
delivered and exact (solver-internal row/column equilibration +
Levenberg damping + floor-reached stop, flags default OFF = legacy
bit-identical; `tests/test_v5_tight_scaled.py` (8), tight fleet 28
green); the medium live-seed e2 read on a non-pre-registered ≤1e-10
threshold = SuperLU pivot-order machine floor through cond ~ 1e10,
adjudicated PASS under the cond-aware read tol = max(1e-10, 10·κ₁·eps)
(~4-decade margin, VERDICT §3); the amended seeds sit INSIDE the 10× floor
band from iter 0 ⇒ no above-band window by construction — fallback:
medium floor_reached at iter 5 at the same merit, coarse still
descending below GV5.1, k=1 standalone F_BL −31 % / merit 2.3× below,
μ rejection-retries 0 (scaling the active ingredient); the window
question reframed to an above-band-seed protocol → **GV5.1c ✓ EXECUTED
2026-07-24 (2 PASS / 1 FAIL / 7 RECORDED,
`bench/studies/v5_1c_above_band_window/`, design record
`docs/design_track_v.md` §15)**: calibrated above-band seeds (δ×(1+ε),
ε = 1e4 → F_BL ≈ 1e4× the floor band) — the pre-floor slope-2 window
MEASURED: NO quadratic regime above the floor (λ = 0.5-capped halvings
p = 1.00 by construction; then a mid-range stall at F_BL ~ 1e-2, never
reaching the band; binding medium median p = 0.56 honest FAIL; μ retries
0 again; band (a) PASS with the cond-aware e2 tolerance pre-registered);
the tight-Newton obstacle is bigger than the floor — a mid-range descent
barrier 3–4 decades above it → **GV5.1d ✓ EXECUTED 2026-07-24 (2 PASS /
1 FAIL / 7 RECORDED, `bench/studies/v5_1d_near_band_window/`, design
record `docs/design_track_v.md` §16)**: near-band seeds (T1 = [1e-4,
1e-3]; coarse 5.42× / medium 35× the band) — NO quadratic basin
adjacent to the floor either (coarse crawls to 24× floor, never
entering the band; medium's first accepted step moves F_BL AWAY from
the band, then crawls to 493×; binding medium median p = 1.17 honest
FAIL; μ retries 0 a third time) — the flat/ragged merit neighborhood
extends down to within ~1.5 decades of the floor: basin hunting
exhausted (GV5.1b/1c/1d) → **GV5.5 ✓ EXECUTED 2026-07-24 (2 PASS /
1 FAIL / 9 RECORDED, `bench/studies/v5_5_te_floor/`, design record
`docs/design_track_v.md` §17)**: the standalone TE-band (B, δ)
formulation item — route (a) variant V1 = TE-outflow row replacement
(rows 6i+0/6i+2 first-order extrapolation, exact Jacobian rows,
default-OFF flag `te_extrapolate`) does NOT break the floor: the
binding m2 (original-system residual at the V1 terminal) = 5554×
(coarse) / 245998× (medium, 8-thread scatter clause) the floor — the
pre-registered "worse" clause; the damage peaks at the LE suction zone
(x_c ≈ 0.027, F_B/F_Psi rows), not the TE; band (a) FD PASS both
levels; the flag stays default-OFF (legacy bit-identical) and the
escalation ladder stays registered-not-opened (opening = user's
adjudication) → **GV5.2 ✓ EXECUTED 2026-07-25 (band (b) FAIL + the
loose-recipe transonic-limit anatomy, `bench/studies/v5_2_rae2822/`,
design record `docs/design_track_v.md` §18)**: RAE2822 P1/P2 VII vs the
committed experiment — band (a) TE wedge 9.46°/9.92° mesh-crease vs
12.91° ordinate fit, quadratic available, no fallback; **band (b)
FAIL**: every computed shock 0.06–0.10 c downstream of the bracket
(medium P1 terminal 0.6288 out of band, leg non-converged RECORDED;
medium P2 k = 4 transpiration runaway §6 RECORDED), worsening with
Mach/α, not improving with refinement; only 1/4 legs converged (coarse
P1, 7 outer); M_peak 1.365/1.306 two outside-envelope RECORDED; the
`cut_wake` Kutta probe gained a bisector-normal fallback for the
RAE2822 reflex camber (fires only where the old code raised) and the
FP driver a pre-registered rescue chain (strict → Mach continuation →
honesty-guarded stall acceptance) against the M ≥ 0.725 shock-cell
plateaus; reading: the loose displacement-thickness feedback is too
weak at M ≥ 0.725 ⇒ the next transonic-VII reads come from the
tight/augmented path, not loose-loop tuning → **GV5.3 ✓ EXECUTED
2026-07-25 (band (b) honest FAIL + band (a) RECORDED input-limited,
0P/1F/17R, `bench/studies/v5_3_m6_cp/`, design record
`docs/design_track_v.md` §19)**: M6 wing TEST 2308 M0.8395/α3.06 vs the
committed 7-station Cp — band (a): Δcl_KJ −2.20 % medium / −1.03 %
coarse, direction DOWN both estimators but under the A4 2.5 % floor
(input-limited; most of the medium move arrives with the late k = 9–10
separation-patch event); **band (b) FAIL**: the viscous Cp does NOT
move toward experiment — 1/5 unmasked stations improved, pooled RMS
0.1288 → 0.1299 increased (coarse 0/5), every |ΔRMS| < 0.05 flagged
input-limited (a direction verdict; tip-masked stations clean);
neither level converges ≤ 10 outer (the GV5.0/GV5.2 signature); the
FP rescue chain load-bearing (10/21, 10/22 stall-accepts); the k = 0
inviscid baseline anchored to the committed P14 cl_KJ 0.2823/0.2688
via the W1 wiring guard (the first-execution fire root-caused to a
driver cold-start short-circuit — the ramp never ran — fixed under
addendum #1; NOT an 8-thread branch scatter); reading: the 3-D
counterpart of GV5.2 — further reads belong to the tight/augmented
path → **GV5.4 ✓ EXECUTED 2026-07-25 (0P/1F/17R,
`bench/studies/v5_4_cost/`, design record
`docs/design_track_v.md` §20)**: augmented-step cost on the 124,216-DOF
M6 medium W2 system (the GV5.1b scaled+damped driver + an injectable
`step_solve` solve callback — library change, default None = splu
bit-identical, +2 tests `tests/test_v5_tight_scaled.py`) — band (a)
RECORDED: augmented step 22.93 s vs the in-session inviscid anchor
3.05 s/step = **7.53×, above the ≤ ~2× reference band** (recorded
either way per the registration; 4/5 steps carry capped-GMRES work;
8-thread walls not comparable to 16-thread entries); **band (b) honest
FAIL (D5)**: the block preconditioner does NOT work at medium —
rung-1 block-Jacobi diverges (rel_res 5.75e4, the φ–BL off-diagonal
coupling too strong), rung-2 exact-BL Schur converges 1/4 steps
(2.66e-8 @277 it then stalls vs rtol 1e-8 @300 cap — pure AMG-φ
cannot represent the J_hB·J_BB⁻¹·J_Bh Schur correction); "measure
before Schur" answered: splu(J_BL,BL) setup 1.8 s (elimination cheap;
Krylov convergence is the bottleneck), AMG-φ 0.2 s / ILU-BL 2.5 s;
guards W1/W2/W3 PASS (cl_p 0.26429 vs the addendum-#4 P14-locked
0.2646; FD median φ 8.7e-12 / Γ 7.3e-12 / BL 0; the IBL floor
trajectory intact); reading: wing-scale augmented Newton needs a
stronger (Schur-aware, (A,Ψ)-structured) reduced-space preconditioner
before cost reads into the ≤ ~2× band; the EW-forcing variant
registered-not-opened (user adjudication); **V5 CLOSED — all five
gates executed** → **GV5.6 ✓ CLOSED 2026-07-25 (0P/1F/17R,
`phases/p1/cases/analysis/v5_6_schur_prec/`, design record
`docs/design_track_v.md` §21)**: the GV5.4 registered Schur-aware
preconditioner follow-up executed (user-adjudicated opening; the GV5.4
system/seed/protocol verbatim, NO library change — W4 diff-scope
clean, suite baseline unchanged) — rung 3 = exact-BL Schur +
bdiag(AMG(Ŝ_φφ), M_Γ) with the explicitly assembled Ĉ =
J_hB·D_BB⁻¹·J_Bh (per-node 6×6 D_BB = the quasi-simultaneous local BL
response; n_fallback 0, t_corr 0.2 s, t_lu 1.8 s reproduced), rung 4 =
block upper-triangular full-system; **band (b) honest FAIL**: the
pre-registered hypothesis falsified in its naive form — the corrected
AMG matrix catastrophically worse than plain AMG(J_φφ) at medium
(rung 3 rel_res 0.664 @242 it vs GV5.4 rung-2's 2e-7..6e-5 stagnation
band; rung 4 0.68 → 1.06, every step capped info=5); the
coarse/medium split — rung 3 WORKS at coarse (5/5 converged, 120–209
it, RECORDED); **band (a) RECORDED** 23.88 s / 3.03 s = 7.87× (≈
GV5.4's 7.53×); guards W1/W2/W3 PASS (cl_p 0.116 % vs the P14 lock;
FD medians 9.6e-12 / 7.3e-12 / 0; the IBL-floor trajectory intact);
reading: the "AMG on a sparsified Schur" direction measured dead at
wing scale — the inter-node (A,Ψ)-structured M_BB route + the
EW-forcing variant registered-not-opened (user adjudication);
V4-reopen trigger considered, NOT invoked (stays parked); **V6 ✓
CLOSED 2026-07-25 — GV6.0 RULED (Option A: conforming-only +
producer (i); the LS leg + the solved wake IBL = recorded follow-ups)
· GV6.1 ✓ CLOSED (6 PASS / 0 FAIL / 7 RECORDED) · GV6.2 ✓ CLOSED
(0 PASS / 0 FAIL / 24 RECORDED)** (conforming wake-sheet δ* source
shipped: `pyfp3d/viscous/wake_sheet.py` + the
`CouplingConfig.wake_transpiration` default-OFF hook = legacy
bit-identical, `tests/test_v6_wake_sheet.py` (8); VERDICTs
`bench/studies/v6_1_wake_sheet/VERDICT.md` +
`bench/studies/v6_2_measured_effect/VERDICT.md` — GV6.1 (a)(i)/(a)(ii)
δ*_wake = 0 bit-identity PASS, (b) sign-pin MMS PASS (antisym
0.81 %, jump 0.44 %, empirically pins the 2026-07-25 per-face ½ṁ
addendum), (c) W2 TE-continuity every outer PASS; GV6.2 measured
on/off Δ-cl +0.00015 (+0.0547 %) / TE max |ΔCp| 0.00250 =
0.022×/0.051× the A4 input band = NOT significant, L-robust over
L_rel {0.5, 1.0, 2.0} c (the `CouplingConfig.wake_l_rel_chords`
plumbing, default = bit-identical, +1 test; 1.0 c stays pinned);
XFOIL wake Option A (user ruling, reference_data untouched, G3 polar
in-line): direction agrees, rate 0.454 c vs 1.0 c recorded-not-tuned,
TE anchor low (the GV3.1 δ* caveat); harness finding logged: numba
cache-load is NOT bit-faithful to fresh-compile in
`pyfp3d/viscous/` — bit-identity A/Bs must pin one cache mode on both
legs, `results/ab_cache_mode_isolation.csv`; producer (ii) NOT opened
(the significance condition NOT met on the GV6.2 recorded reading,
user-adjudicated 2026-07-25 per the GV6.0 clause)); Track A — A1, A2,
**A3 ✓ CLOSED 2026-07-18**, **A4
RECORDED 2026-07-22** (wall u_e error-band study = Track-V input-quality
prerequisite: medium smooth-wall band ≈2.5% peak / 0.04·U∞ max-norm / O(h),
`phases/p1/cases/analysis/a4_ue_error_band/`) (A3 = response
to the 2026-07-17 independent inspection: docs consistency + cross-path
hardening + the C1 Jacobian verification, see
[docs/inspection/](docs/inspection/); the footer's "A3 ◐" was itself one of
the close-out-debt findings, fixed 2026-07-19). Next phase = the user's call.
★★★ **PHASE THREE opened 2026-08-11. Task 1 done: the LEVEL-SET ROUTE IS DELETED**
(ruling D5) — 9 library files / **4624 lines**, `pyfp3d/wake/` gone, `post/unified.py`
collapsed onto its conforming half. ONE wake route remains: `mesh/wake_cut.py` +
`constraints/te_pressure.py` + `solve/newton.py` / `picard.py`, with `tip_taper`
(B31/B32, carrying a **−1.3 % cl** model bias).
Baselines: always-on **479 passed / 12 skipped / 2 xfailed** (2026-08-12 = 474 measured in
full @478.40 s plus 5 non-interacting soft-membership asserts; earlier it went
468 → 472 → 468 in one day -- the +4 were the temporary `sigma_scale` instrument's locks
and were deleted with it at its registered expiry). Earlier the same count (2026-08-11, 494 s @8
threads quiet; +11 = `tests/test_meshgen_structured.py`, closing the structured
generator's zero-coverage debt). After the level-set deletion it was **457 / 12 / 2**, gated
**466 / 1 / 4** (2:08:50), fast tier **5 groups**. Both accounts close item by item —
docs/dev_phase_three/20260811-0100-ls-deletion-verdict.md.
⚠ `phases/p1/` still imports the deleted modules and is NOT runnable; the archive is a
snapshot, and a working tree is `git worktree add ../up3d-prereorg d224223`.

★★ **PHASE TWO (2026-07-28 .. 2026-08-10) — read
[docs/dev_phase_two/PHASE_TWO_CAPABILITY_BOUNDARY.md](docs/dev_phase_two/PHASE_TWO_CAPABILITY_BOUNDARY.md)
first.** The whole block below this line is the PHASE-ONE footer, frozen with the
phase-one docs on 2026-07-28 (docs/dev_phase_two/roadmap.md §8); its status verdicts
are history, its technical facts stay citable as phase-one records.

Phase-two baselines, measured:
- always-on suite **700 passed + 28 skipped + 2 xfailed**, 0 failed
  (1464 s @8 threads / 1669 s @16 threads — quote the thread count, it is not a
  cosmetic detail: the same M6 solve measures 566.6 s uncapped against 113.7 s capped);
- FAST capability tier **7/7 green, 670 s** (`bench/run_capability_locks.py`);
- FULL gated set **720 passed + 2 skipped + 8 xfailed**, 0 failed (3 h 04) — the first
  all-green since the 2026-08-04 round-tip switch, against 7 failed on 08-06, one of
  which had been red for seventeen days because a retired premise was never grepped out
  of the tests.
- The 8 xfailed include **four deliberate strict xfails** (b7 ×3 + b22 medium): on the
  production round-tip meshes those level-set gates have no clamp-free state to anchor,
  and per ruling **D5** that route is ABANDONED — they are abandoned-route records, not
  open obligations, and they go away with the files in phase three.
- Phase two's own metric verdicts live in `docs/dev_phase_two/progress.md`'s tracking
  table; the one that changes what to work on next is **M3a, measured UNREACHABLE for
  an inviscid solver** (model floor 0.0516–0.0707 against a viscous experiment).

---

Default suite (PHASE ONE, frozen): **652 passed + 25 skipped + 2 xfailed** (2026-07-25, Track V
V6 GV6.2 (the measured wake-IBL on/off effect vs the A4 band: Δ-cl +0.00015
(+0.0547 %) / TE max |ΔCp| 0.00250 = 0.022×/0.051× the A4 input band = NOT
significant, L-robust over L_rel {0.5, 1.0, 2.0} c (1.0 c pinned); the XFOIL
wake direction check Option A — direction agrees, rate 0.454 c vs pinned
1.0 c recorded); full-suite measured 652 @1455.80 s **@8 threads**
(temporary 8-core session constraint, user-directed; NOT comparable to the
16-thread ledger entries); +1 vs 651 = `tests/test_v6_wake_sheet.py` (the
GV6.2 `wake_l_rel_chords` plumbing test). Previous 651: V6 GV6.1 (the
conforming wake-sheet δ* source: (a)(i)/(a)(ii) δ*_wake = 0 bit-identity
PASS — the (a)(ii) harness runs both legs fresh-compile, the numba
cache-load infidelity discipline — (b) sign-pin MMS PASS empirically
pinning the per-face ½ṁ addendum, (c) W2 TE-continuity every outer PASS;
(d) Δ-cl +0.00015 / TE-region max |ΔCp| 0.00250 RECORDED for GV6.2);
full-suite measured 651 @1606.31 s **@8 threads**; +7 vs 644 =
`tests/test_v6_wake_sheet.py` (7 new: W3 construction, the producer
identity, the zero-field RHS, (a)(i) bit-identity, the sign-pin MMS, (a)(ii)
fresh-compile A/B, the fold pairing). Previous 644:
V5 GV5.4 (the augmented-step cost on M6 medium: band (a) RECORDED — the
augmented step 22.93 s vs the inviscid step 3.05 s = 7.53×, above the
≤ ~2× band; band (b) honest FAIL — the block preconditioner does NOT work
at medium: block-Jacobi diverges, exact-BL Schur stalls 1/4); full-suite
measured 644 @1230.61 s **@8 threads** (temporary 8-core session constraint,
user-directed; NOT comparable to the 16-thread ledger entries); +2 vs 642 =
`tests/test_v5_tight_scaled.py` (2 new: the `step_solve` callback wiring +
the default-None splu bit-identity guard). Previous 642:
V5 GV5.3 (the M6 wing direction+magnitude check vs committed Cp: band (b)
honest FAIL — the viscous Cp does NOT move toward the committed 7-station
experiment (1/5 unmasked stations + a pooled increase); band (a) RECORDED
input-limited, Δcl_KJ −2.20 % DOWN under the A4 floor); full-suite measured
642 @1233.89 s **@8 threads** (temporary 8-core session constraint,
user-directed; NOT comparable to the 16-thread ledger entries); +0 vs 642 =
no test changes (the gate added no library/test code; wall re-measured).
Previous 642:
V5 GV5.2 (the RAE2822 transonic VII vs committed experiment: band (b) FAIL +
the loose-recipe transonic-limit anatomy — every computed shock 0.06–0.10 c
downstream of the bracket, 1/4 legs converged); full-suite measured 642
@1245.33 s **@8 threads** (temporary 8-core session constraint,
user-directed; NOT comparable to the 16-thread ledger entries); +6 vs 636 =
`test_meshgen_rae2822.py` (5) + `test_p2_wake_cut.py` (1,
test_kutta_probes_cambered_te). Previous 636:
V5 GV5.5 (the TE-band (B, δ) formulation item: the V1 TE-outflow row
replacement does NOT break the floor — binding m2 5554×/245998× the floor,
the "worse" clause; damage peaks at the LE suction zone, not the TE; flag
default-OFF); full-suite measured 636 @1402.11 s **@8
threads** (temporary 8-core session constraint, user-directed; NOT
comparable to the 16-thread ledger entries); +9 vs 627 =
`test_v5_te_outflow.py` (9). Previous 627:
V5 GV5.1d (the near-band window read: NO quadratic basin adjacent to the
floor either — near-band seeds stall immediately (coarse crawling to 24×
floor, medium's first step moving AWAY from the band), binding medium
median p = 1.17 honest FAIL); full-suite measured 627 @1340.77 s **@8
threads** (temporary 8-core session constraint, user-directed; NOT
comparable to the 16-thread ledger entries; wall markedly below the
GV5.1c-era 3903 s on the same thread count — machine/cache conditions
differ, quoted flagged); +7 vs 620 =
`test_v5_near_band_seed.py` (7). Previous 620:
V5 GV5.1c (the above-band window read: NO quadratic regime above the floor —
λ-capped halvings + a mid-range stall, binding medium median p = 0.56 honest
FAIL); full-suite measured 620 @3903.16 s **@8 threads** (temporary 8-core
session constraint, user-directed, machine idle; NOT comparable to the
16-thread ledger entries); +9 vs 611 =
`test_v5_above_band_seed.py` (9). Previous 611:
V5 GV5.1b (scaled+damped augmented Newton; machinery exact, band (b) window
question reframed), 611 measured @6556.77 s @16 threads (wall polluted by
co-tenant load ~70–80, quoted flagged); previous 603:
V5 GV5.1 (augmented tight (φ, Γ, U) Newton; FD exactness PASS both levels,
quadratic tail HONEST FAIL on the IBL floor), 603 measured @1537.09 s;
previous 583:
V5 GV5.0 (M6 subsonic loose-coupling bridge, RECORDED entry check),
583 measured @1218.05 s; previous 578:
V3 loose coupling + GV3.1/3.2/3.3, 578 measured @1637.39 s; previous 571:
V2 transpiration channel + GV2.1, 571 measured @1321.89 s, NOJIT lane
17/17 @163.78 s; previous 554:
V1 IBL3 core + GV1.1, 554 measured @1462.64 s, NOJIT lane 35/35; previous 519:
B28–B32 close-out +
G1.6 Option C re-spec, 516 measured @1223.39 s + 3 TestG16Respec asserts);
the 16 M1 tests
skip unless the gitignored M6 meshes are regenerated (~30 s); the gated
skips include B21's freeze-capture lock and B22's 3-D anchor locks.
