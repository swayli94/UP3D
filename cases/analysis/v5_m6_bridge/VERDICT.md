# GV5.0 VERDICT — M6 subsonic loose-coupling bridge (RECORDED)

Gate: `docs/roadmap/track_v.md` **GV5.0** (2026-07-22, user-directed;
RECORDED entry check, no pass/fail). Pre-registration:
`PRE_REGISTRATION.md` (committed f263424 before the first execution).
Runner: `run.py` (regenerates every artifact in `results/`). Executed
2026-07-22/23 on `kimi/track-v5-tight-coupling`.

**Answer to the bridge question, up front:** the V3 loose loop does NOT
converge on the M6 wing within 10 outer iterations at EITHER level, in two
measured, mesh-dependent regimes:

- **coarse — separation-zone feedback runaway**: a root-upper-TE
  separation-influenced patch (H ≈ 4–5.5) drives δ*↔ṁ↔u_e feedback with
  gain > 1 (ṁ_max ×12.4 monotone over k, δ*max 3.2e-3 → 1.5e-2, cl still
  drifting at k = 10). The GV3.3 stern/Veldman class, now measured on a
  lifting wing.
- **medium — separation resolved away, bounded limit cycle**: refinement
  removes the patch entirely (0 nodes with H > 3.5 at x/c > 0.9; root-TE
  δ* healthy at 4.3e-3 with H ≈ 1.85), the runaway disappears (ṁ_max
  peaks at k = 2 then decays and settles ~0.063), but the loop still
  never meets tol_ds = 1e-3 — it sits in a bounded δ* oscillation
  (Δδ*/δ*max ≈ 2–12 % per k) while cl settles to ±0.0002.

The bulk 3-D IBL solution is sane at both levels (physical δ*(z)
distributions, small crossflow, quasi-2-D lock retained only at the LE
pin), ΔCL is DOWN on both lift estimators, and the FP side stays clean
throughout (every warm Newton solve converges; the GV3.3 guardrails never
fire). Recorded as the motivation evidence for V5's tight coupling: the
loose loop is not sufficient on the 3-D lifting wing, either in the
separation-influenced regime (coarse) or as a bounded-but-unconverged
cycle (medium).

## 1. Configuration as executed

Per the pre-registration, one logged addition, no numerics deviations:
ONERA M6 coarse/medium (flat-cap family, wake-cut at ingestion),
M∞ 0.5 / α 3.06, forced transition x_tr/c = 0.05 both sides,
Re_MAC = 11.72e6 (re_chord = 1.81405e7 per meter), production Newton
recipe (`farfield_spanwise_gamma`, `precond="direct"`, pressure Kutta
(P14), `tip_taper` vanish_smooth r_c = 0.05·b_semi (B32); cold start
probe→pressure chain), LE-band x/c ≤ 0.02 laminar pin, tip band
z > 0.95·b_semi pinned + ṁ-masked, A4 zones LE+tip, eps 0.005/0.02,
ω = 1.0, tol_ds 1e-3, ≤ 10 outer iterations. The logged addition: a
per-outer-iteration progress print in the probe (k, cl_p, cl_kj, elapsed)
— pure logging, added after the first medium run died blind at a 4 h
task timeout; no numerics touched.

IBL surfaces: coarse 2690 nodes / 5152 tris (tip-masked 142, LE-band 263);
medium 10205 nodes / 19958 tris (tip-masked 513, LE-band 1055).

**Wall-time caveat**: the medium wall time (23266 s) is contention-
polluted — two unrelated heavy jobs (another project's case studies,
~14 cores) shared the machine for the first ~3.5 h (k = 0 took 7268 s vs
the P14 demo's ~90 s for the identical recipe; post-contention outer
iterations took ~15–27 min). Quote wall times as polluted, not as solver
cost. The coarse run (1465 s) ran on an idle machine.

## 2. Recorded headline numbers

| metric | coarse | medium |
|---|---|---|
| loop converged (≤10, tol 1e-3) | False (10 iters) | False (10 iters) |
| regime | runaway (gain > 1) | bounded limit cycle |
| ds_change_rel at k=10 | 0.613 | 0.070 |
| ṁ_max trend over k | ×12.41 (0.026 → 0.328) | peak k=2, settles ~0.063 |
| Δcl_p (viscous−inviscid) | −0.0105 (−5.20 %) DOWN | −0.0050 (−2.41 %) DOWN † |
| Δcl_kj | −0.0097 (−4.75 %) DOWN | −0.0044 (−2.10 %) DOWN † |
| max\|B\|/max\|A\| (live) | 0.0630 (z/b = 0.948, tip edge) | 0.0723 (z/b = 0.145) |
| max\|Cτ2\|/max\|Cτ1\| (live) | 0.902 (TE band, small/small) | 0.126 |
| IBL residual floor | 3.3e-6 … 8.4e-6 | 1.1e-6 … 1.8e-6 |
| TE separation content (H>3.5, x/c>0.9) | 10 nodes (root upper patch + 1 lower) | 0 nodes |
| wall time (polluted, §1) | 1465 s | 23266 s |
| final FP n_newton / converged | 4 / True | 2 / True |

† medium Δcl below the A4 2.5 % input floor → flagged **input-limited**
per the pre-registered A4 discipline (a smaller move than the inviscid
u_e input band is RECORDED-only).

2.5-D reference count (GV3.2, same criterion, ω = 1.0): 4 (coarse) / 5
(medium) — the 3-D wing loop does NOT meet it at either level.

IBL inner solves hit the 100-iteration cap at every outer iteration of
both levels (~1e-6–8e-6 residual floor) — the SAME floor class as the
GV3.1 airfoil (1.1e-6, also capped), so the inner tolerance is not the
loop driver.

## 3. The instability, localized and measured

Coarse (`gv5_0_history_coarse.csv`, `gv5_0_surface_coarse.csv`): δ*max
oscillates and grows (3.2e-3 → 1.5e-2), ṁ_max grows monotone ×12.4, cl_p
drifts down (0.2017 → 0.1912, still moving at k = 10). The top-12 live δ*
nodes all sit at the TE (x/c > 0.955): 11 of 12 on the **upper surface
over the root ~8 % span** (z/b 0.012–0.084, δ* up to 1.29e-2,
H 2.6–5.5), plus one lower-TE node at z/b = 0.24 (H 17.8). The lower
surface at the root stays thin and healthy (δ* 1.4–3.2e-3, H 1.4–1.9).
Separation-indicated content (H > 3.5) is 70/2548 live nodes, of which
**58 sit at x/c ≤ 0.23** — the known near-stagnation θ→0 H-spike artifact
class just behind the LE pin (δ* there tiny and healthy; the A4 LE
input-limited zone), not separation. The genuine TE content is ~10
nodes, dominated by the root-upper patch. Mechanism: the root upper-TE BL
thickens toward mild separation under the α 3.06 suction-side APG;
thickening δ* ⇒ growing transpiration source ⇒ FP u_e distortion ⇒
further thickening (Veldman gain > 1; the GV3.3 stern class). The tip
mask works as designed — no tip-side runaway (tip-band δ* ≤ 1.1e-3,
H ≤ 2.83, frozen seeds, ṁ ≡ 0).

Medium (`gv5_0_history_medium.csv`, `gv5_0_surface_medium.csv`): the
patch is gone — zero H > 3.5 nodes at x/c > 0.9; the max-δ* nodes sit at
the root TE with healthy H ≈ 1.8–1.9 (δ* ≤ 4.3e-3). ṁ_max peaks at
k = 2 (0.121) then decays to ~0.063; cl_p settles 0.2035 ± 0.0002. What
remains unconverged is a bounded max-norm δ* oscillation
(Δδ*/δ*max 0.023–0.125, no trend) — a limit cycle, not a blowup.

## 4. The bulk 3-D solution (the bridge's positive payload)

- **δ*(z) spanwise distributions** (`gv5_0_dstar_spanwise_*.csv/.png`):
  physical at every station of both levels — smooth growth with x/c,
  upper > lower, spanwise-uniform mid-span with the recorded root
  behaviour (coarse: separation uptick; medium: healthy maximum at the
  root TE). These are the GV5.3 band pre-registration feed.
- **Crossflow, first live 3-D exercise** (`gv5_0_crossflow_*.png`):
  max|B|/max|A| = 0.063 (coarse, argmax at the tip-mask edge) / 0.072
  (medium, z/b = 0.145) — small everywhere; the coarse
  |Cτ2|/|Cτ1| = 0.90 max sits at the TE edge band (Cτ1 small there ⇒
  small/small ratio artifact; medium's 0.126 confirms). The LE pin keeps
  B ≡ 0 there (seed states). Ψ maxima (6.9e3 coarse / 4.0e3 medium) sit
  at the LE-band edge (profile twist where q → 0; same artifact class as
  the front H spikes).
- **ΔCL DOWN on both estimators** (direction as pre-registered): coarse
  cl_p −5.20 % / cl_kj −4.75 % (still drifting at k = 10 — quoted at the
  final iterate, honestly flagged); medium cl_p −2.41 % / cl_kj −2.10 %
  (settled; both under the A4 2.5 % input floor ⇒ input-limited, so the
  medium magnitude is RECORDED-only by the pre-registered discipline).
- **FP side clean**: every warm-started pressure-Kutta Newton converges
  (final: 4 steps coarse / 2 medium), no clamp/floor events, the GV3.3
  non-convergence guardrail never fires.

## 5. Consequences (recorded, not decided here)

1. GV5.0's bridge answer = the loose loop is NOT sufficient on the 3-D
   lifting wing: runaway where a separation-influenced zone exists
   (coarse), bounded-but-unconverged cycling once refinement removes it
   (medium). This is the measured motivation for V5's tight
   (augmented-Newton) coupling, and a live wing case for the V4 reopen
   trigger (track_v.md V4 note; the trigger text reads "V5's augmented
   Newton stalls, or closed-body viscous cases enter scope" — here the
   LOOSE loop is what stalls, on an OPEN wing, which strengthens the
   tight-coupling route's priority rather than reopening V4 by itself).
2. GV5.3 inherits: (a) the δ*(z) CSVs as its band pre-registration feed;
   (b) the tip mask (z > 0.95·b_semi) validated as a working mask at both
   levels; (c) the root-upper-TE zone flagged as the separation-sensitive
   region — coarse-level conclusions there are mesh-limited (the medium
   BL stays attached), so GV5.3's direction+magnitude check should anchor
   on the medium level.
3. The medium Δcl magnitude (−2.4 %) sits under the A4 input floor: a
   tighter viscous ΔCL claim on the M6 wants the O(h²) input band (the
   unchosen G1.6 route (b) per the scope guards) or the fine level — both
   out of GV5.0 scope, recorded for GV5.3's planning.

## 6. Artifact index

- `results/gv5_0_history_{coarse,medium}.csv` — per-outer-iteration loop
  record (δ*max, Δδ*, ṁ_max, IBL stats, cl_p/cl_kj/cd_p).
- `results/gv5_0_surface_{level}.csv` — per-node surface dump (geometry,
  x/c, side, δ*, H, B, Ψ, Cτ1/2, cf1, pin/mask flags).
- `results/gv5_0_dstar_spanwise_{level}.csv/.png` — the GV5.3 feed.
- `results/gv5_0_crossflow_{level}.png` — planform δ*, |B|, |Cτ2|/|Cτ1|.
- `results/gv5_0_cl_history.png` — cl_p / cl_kj vs outer iteration.
- `results/summary.csv` — all rows RECORDED, 0 FAIL.
