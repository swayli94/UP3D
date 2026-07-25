# GV6.1 PRE-REGISTRATION — conforming wake-sheet δ* source

**Committed before the first code change** (Track-V discipline). Date:
2026-07-25. Branch `kimi/track-v6-gv6-0`.

**Gate text** (`docs/roadmap/track_v.md`, V6): *"**GV6.1 conforming sheet
source**: δ*_wake enters via the `constraints/wake.py` reduce RHS
(Tᵀ b_wake); TE thickness continuity asserted; δ*_wake = 0
**bit-identical**."*

**Ruling basis** (GV6.0, user 2026-07-25,
`cases/analysis/v6_0_design_adjudication/DESIGN_ADJUDICATION.md` §6):
**Option A** — V6 closes conforming-only (the LS sheet-source leg is a
recorded follow-up); **producer (i)** — prescribed δ*_wake (the solved wake
IBL is follow-up producer (ii)).

**Design constraints** (registration, unchanged): TE kink absorbed by Drela
local-basis adaptation; **straight wake + mass-transpiration relaxation, no
geometric relaxation**.

---

## 1. Scope

- Conforming path only (the GV6.0 ruling). Driver scope: the **loose
  Picard leg** — `solve/picard.py::solve_subsonic_lifting`, the
  `body_source_rhs` channel (`solve/picard.py:609-617,659`). The GV6.2
  measurement case (GV3.1) runs on this leg.
- The Newton (`solve/newton.py` `external_rhs`) and tight
  (`viscous/tight.py` b(m)) driver legs are **NOT wired in this gate** —
  registered follow-up wiring (the Tᵀ b channel is identical; wired when a
  gate needs it). See §8.
- New machinery, all behind a default-OFF flag (`wake_transpiration`;
  flag-OFF = legacy bit-identical):
  1. a wake `SurfaceMesh` + an open-chain station table built on
     `wc.wake_faces_minus` (the plus side is the slave copies; they fold
     into the master rows under T);
  2. the prescribed δ*_wake producer (§2);
  3. the b_wake wiring into the loose loop (§3);
  4. the TE-continuity runtime assert (§2, W2);
  5. tests (§6).
- **No existing solver-path code is modified in behavior**: the only touch
  point is an additive, flag-gated term in the Picard `body_source_rhs`
  assembly.

## 2. The prescribed δ*_wake producer (producer (i), pinned a priori)

Updated once per loose-loop outer iteration k from the current wall-IBL
state:

- **TE confluence** (the registered TE thickness continuity — satisfied by
  CONSTRUCTION and asserted at runtime, W2):
  `δ*_TE = δ*_upper(TE) + δ*_lower(TE)`,
  `θ_TE  = θ_upper(TE)  + θ_lower(TE)`,
  read off the wall IBL state at the two TE copies (the wake-cut split
  nodes; `side_node` distinguishes them, `coupling.py:143-145`).
- **Straight-wake relaxation** (mass-transpiration relaxation; no geometric
  relaxation): along the straight-wake station arc s (s = 0 at the TE,
  downstream positive),
  `δ*_wake(s) = θ_TE + (δ*_TE − θ_TE) · exp(−s / L_rel)`,  `L_rel = 1.0·c`.
  Physics basis: with no wall shear the wake momentum thickness is
  approximately conserved (far-field drag theorem), so the wake profile
  fills out — H_wake = δ*/θ relaxes from its TE value toward 1. The
  exponential is the **MODEL CHOICE — RECORDED**: the single free constant
  L_rel is pinned a priori at 1.0 c; the L_rel sensitivity is a GV6.2
  recorded item, NOT tuned in this gate. The construction gives
  δ*_wake(0) = δ*_TE identically (the W2 assert) and δ*_wake → θ_TE
  downstream (H → 1).
- **u_e,wake / ρ_e,wake**: per wake surface node, the average of the two
  sides' FP-recovered values,
  `u_e(wake node) = ½·(ue_vol[minus node] + ue_vol[plus slave node])`,
  `ue_vol` from the loose loop's per-zone recovery (`post/surface.py`
  discipline, `coupling.py:752-763`); ρ_e from the same state via the wall
  chain's isentropic relation.
- **ṁ_wake** = `transpiration_from_delta_star(wake_sm, ρ_e, u_e, δ*_wake)`
  — reused verbatim (SurfaceMesh-generic): ṁ = div_Γ(ρ_e u_e δ*). The
  downstream end of the wake strip is left natural (no pin; the ṁ field is
  recorded as computed).

## 3. The sheet-source channel (the registered mechanism)

- `b_wake = assemble_transpiration_rhs(nodes, wake_faces, ṁ_full)` with
  `wake_faces = wake_minus ∪ wake_plus` and `ṁ_full` the full-length
  volume-node vector: minus-side nodes carry ṁ_wake (via the wake
  SurfaceMesh's `volume_node_of`); each plus-side **slave** node carries the
  same ṁ_wake value as its minus-side partner. The two sides' Galerkin
  loads fold under Tᵀ into the master row — that fold IS the weak sheet
  flux (GV6.0 survey §1.1).
- Injection: added to the existing `body_source_rhs` slot of
  `solve_subsonic_lifting`, flag-gated.
- **Sign derivation** (the MMS pin of band (b)): the transpiration load is
  the negated Galerkin load — positive ṁ = blowing out of the body
  (`transpiration.py:19-27,114`). On the wake sheet the two faces carry
  opposite outward normals; folding sums both sides' loads into the master
  row. In the continuous limit this realizes the internal jump condition
  `[ρ ∂φ/∂n] = ρ∂φ/∂n|₊ − ρ∂φ/∂n|₋ = ṁ_wake`
  across the sheet: a thickening wake (ṁ > 0) ejects fluid symmetrically
  AWAY from the sheet on both sides. Band (b) asserts exactly this.

## 4. Bands (verdicts)

- **(a) δ*_wake = 0 bit-identical — PASS/FAIL** (binding). Two legs, both
  in-session: (i) flag-ON with a prescribed ZERO δ*_wake field vs flag-OFF;
  (ii) flag-OFF vs the library without this gate's code (A/B on the same
  machine/threads/seed). Each must be bit-identical on the coarse loose
  loop (3 outer). The GV2.1(b) discipline.
- **(b) Sign-pin MMS — PASS/FAIL** (binding). On the wake-cut NACA0012
  coarse strip at U∞ = 0 (dead air, farfield Dirichlet), prescribe a
  spatially uniform ṁ₀ > 0 on the wake sheet. Assert at probe nodes
  adjacent to the sheet mid-chord of the strip (away from the TE and the
  downstream strip end):
  (i) antisymmetry v_n|₊ = −v_n|₋ to discretization accuracy;
  (ii) the jump v_n|₊ − v_n|₋ = ṁ₀/ρ₀ within **5 %** (pinned; the
      GV2.1(a) cylinder-MMS lock is 3 % for a smooth exterior field — the
      internal sheet source is a singular field, so the lock is set at
      5 %, probes kept off the singular endpoints);
  (iii) the sign convention: ṁ₀ > 0 ejects fluid away from the sheet on
      BOTH sides.
  A mismatch is a RECIPE ERROR (fix the implementation, never the
  tolerance) — the GV2.1(a) discipline.
- **(c) TE-continuity assert — PASS/FAIL** (binding, runtime): the
  constructed δ*_wake(0) equals δ*_upper(TE) + δ*_lower(TE) to 1e-12
  relative. A failure means the TE-copy wiring is wrong (recipe error).
- **(d) GV3.1 smoke — RECORDED** (non-binding): the loose loop on the
  GV3.1 case (medium, the committed GV3.1 recipe verbatim) with the flag
  ON vs OFF: on/off Δcl, TE-region Cp shift, outer-count behavior, the ṁ_wake
  field scale. RECORDED only — the measured effect with the A4 input band
  and the XFOIL direction check is **GV6.2's** band, not this gate's. The
  committed GV3.1 medium fixed point is NOT reproducible (the committed
  IBL-floor-trajectory scatter, `v5_tight_coupling/VERDICT.md` §4) — every
  cross-run comparison carries the scatter caveat; the smoke's on/off
  pairing is done in-session on the same seed (A/B discipline).

## 5. Wiring guards (recipe-error raisers)

- **W1**: in-session A/B bit-identity of the flag-OFF loose loop vs the
  library without this gate's code (coarse, 3 outer; same machine, threads,
  seed). Any deviation = wiring error. (Cross-thread / cross-commit
  reproduction of the committed GV3.1 numbers is NOT required — the
  committed scatter caveat.)
- **W2**: the TE-continuity construction identity (band (c)) holds at every
  outer iteration.
- **W3**: wake-SurfaceMesh sanity — the station count matches the wake
  strip; the arc length is strictly monotone from the TE; node areas
  positive; every wake-minus node pairs to exactly one plus-slave (the fold
  pairing); `volume_node_of` maps into the cut mesh's node table.

## 6. Tests (new; the expected suite delta +6)

`tests/test_v6_wake_sheet.py`:

1. wake SurfaceMesh construction + W3 sanity on the NACA coarse strip;
2. the prescribed-producer construction identity (c) on a synthetic wall
   state (pinned δ*_upper/lower, θ_upper/lower → δ*_wake(0), the monotone
   downstream relaxation, the H → 1 limit);
3. zero-field bit-identity: ṁ_wake ≡ 0 → b_wake ≡ exact zero vector (the
   `transpiration.py` zero discipline) + the (a)(i) loose-loop leg;
4. the sign-pin MMS (b) on the coarse strip;
5. the (a)(ii) A/B loose-loop bit-identity (flag-OFF vs gate-free library);
6. the fold-pairing assert (W3) on the strip: every minus-side load lands
   in its master row (checked structurally against T).

## 7. Execution discipline

- Temporary 8-thread session constraint
  (`NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8`,
  user-directed); all wall times quoted flagged non-comparable.
- Committed artifacts are never re-computed; `cases/reference_data/` is
  ground truth, untouched.
- Runner exit 0 = all binding bands PASS; exit 1 = any FAIL (honest-FAIL
  discipline); VERDICT.md + results/ (summary.csv, the smoke CSVs/PNGs)
  committed.
- Addenda: any change to pinned bands/guards/recipe/producer constants
  requires an addendum section appended to this file, committed BEFORE the
  (re-)execution it affects.

## 8. Explicitly OUT of scope (registered follow-ups)

- The LS sheet-source leg (GV6.0 ruling, Option A).
- The solved wake IBL (GV6.0 ruling, producer (ii)) — opened only if the
  GV6.2 measured effect is significant against the A4 input band.
- The Newton / tight driver wiring of b_wake (the Tᵀ b channel is
  identical; wired when a gate needs it).
- Geometric wake relaxation (forbidden by the design constraint).
- The L_rel sensitivity sweep (a GV6.2 recorded item).
- GV6.2 itself (the measured on/off effect with the A4 band + the XFOIL
  direction check; its own pre-registration — NOTE: the committed XFOIL
  reference CSVs carry surface rows only, the DUMP wake rows were discarded
  at generation time; how GV6.2 sources its XFOIL wake reference is a
  GV6.2-pre-registration item, likely a user ruling).

---

## Addendum 2026-07-25 (PRE-CODE design review): the per-face ½ṁ_wake recipe

Filed per §7 BEFORE the first code change (branch `kimi/track-v6-gv6-1`).
A pre-code weak-form review found §3's scatter recipe internally
inconsistent with §3's own jump claim by a factor 2. The bands and guards
are UNCHANGED; the recipe text is corrected as follows.

**The inconsistency.** §3 as registered: "minus-side nodes carry ṁ_wake …
each plus-side slave node carries the same ṁ_wake value", and "folding sums
both sides' loads into the master row … this realizes [ρ∂φ/∂n] = ṁ_wake".
Both statements cannot hold: the Tᵀ fold SUMS the Galerkin loads of the two
coincident face copies into the master row, so two copies each carrying ṁ
give the folded load −2∫Nṁ dS. The wall-transpiration calibration
(`transpiration.py:19-27,114`) is ONE face: b = −∫Nṁ dS realizes
ρ∂φ/∂n_body = ṁ (factor 1). On the sheet the folded equation constrains the
SUM of the two one-sided fluxes (the homogeneous case is the registered
[ρ dφ/dn] = 0 continuity, `constraints/wake.py:18-21`), so two full-ṁ
copies realize **[ρv_n] = 2ṁ_wake**, not ṁ_wake.

**The physics.** ṁ_wake = div_Γ(ρ_e u_e δ*_wake) is the divergence of the
TOTAL wake defect mass flux = the total per-area ejection of the sheet =
the jump strength: continuity of the defect gives [ρv_n] = ṁ_wake, split
symmetrically as v_n^± = ±ṁ_wake/(2ρ) (the source-sheet identity). Band
(b)'s pinned expectation (jump = ṁ₀/ρ₀ within 5 %) encodes exactly this.

**Corrected recipe (pinned).** Each face copy carries **½ṁ_wake**:
`ṁ_full[minus nodes] = ṁ_full[plus slaves] = 0.5·ṁ_wake`. The fold then
contributes −∫Nṁ_wake dS once and the realized jump is ṁ_wake. (Loading
the minus faces only with the full ṁ_wake yields the same folded vector;
the ½-both-copies form is chosen to keep §3's registered both-sides-loaded
mechanism.) Band (b)'s pinned 5 % jump lock is UNCHANGED — it now also
empirically pins this factor (a 2× error lands 100 % outside the lock);
bands (a)/(c), W1–W3, the producer constants (L_rel = 1.0 c), and the §8
scope are all untouched.
