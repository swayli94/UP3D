# P5 medium gate — spanwise-Γ-smoothing investigation (A–E)

**Date:** 2026-07-08 **Verdict:** ❌ **DEAD ROUTE — do not re-attempt.**
Spanwise-Γ smoothing does **not** fix the medium-gate trailing-edge M>2 spike.

This is the lab-notebook record so the investigation is **not repeated**. The
one-line conclusions live in [docs/roadmap.md](../../../docs/roadmap.md) P5
OPEN block, [docs/demo_report.md](../../../docs/demo_report.md) P5 section, and
[docs/agent-rules.md](../../../docs/agent-rules.md); this file keeps the full
method + per-test numbers + caveats behind them.

## The failure being investigated

Medium (350.7k) gate FAILS: M_max 5.204, 8 floored / 4 limited cells, V6
1.73%. Location (verified, `python diagnose_medium.py`): **26/350718 cells
M>2, split 18 ON the wing at the outboard TE (x/c 1.00–1.11, z/b 0.80–0.81) +
8 at the far-field sphere** (x≈10). The single M_max=5.20 cell is far-field;
the dominant M≈2.7–3.1 spikes are the wing-TE cluster, on **well-shaped**
elements (shape-quality q≈0.65, aspect ratio ≤3.1 — not slivers). The band
sharpens under refinement (coarse 1.47 → medium 5.20).

## Hypothesis (tested here, REFUTED)

The wing-TE cluster co-locates with the steepest spanwise Γ roll-off and the
per-station Kutta Γ has a non-monotone dip there (0.076→0.043→0.048 across z/b
0.7/0.8/0.9), so we hypothesised the ~8% per-station Kutta **spanwise noise**,
carried through the wake-cut master–slave Γ-jump, drives the TE spike — and
that spanwise-smoothing the per-station Γ would remove it.

The smoother tried (Nadaraya–Watson Gaussian kernel in z + a Γ(tip)=0 anchor,
applied to the Kutta target each outer step; wired as a `gamma_smooth` hook in
`solve_transonic_lifting`, default off / no-op on single-station meshes).
**Reverted after this investigation** — not in the tree. Self-contained core,
for the record only:

```python
def smooth_spanwise_gamma(gamma, station_z, sigma_frac=0.03, b_semi=None):
    g, z = np.asarray(gamma, float), np.asarray(station_z, float)
    if g.size < 2: return g
    span = float(z.max() - z.min())
    if span <= 0: return g
    sigma = sigma_frac * span
    if b_semi is None: b_semi = z.max() + np.median(np.diff(np.sort(z)))
    zz, gg = np.r_[z, b_semi], np.r_[g, 0.0]          # tip Γ=0 anchor
    w = np.exp(-0.5 * ((zz[:, None] - zz[None, :]) / sigma) ** 2)
    return ((w @ gg) / w.sum(1))[: g.size]
```

## Tests and results (A–E)

| # | Test (mesh) | Result | Reads on the hypothesis |
|---|---|---|---|
| **A** | zero-cost correlation from cache (coarse+medium) | TE-band max Mach vs \|dΓ/dz\|: **coarse +0.25, medium +0.63**; vs Γ scatter: +0.33 / +0.53. Strengthens with refinement. | Looks supportive — but only a correlation. |
| **B** | coarse OFF vs ON σ0.03/0.05 (full solve ×3, ~15 min) | ON σ0.03: outboard-band Mmax **1.478→1.395**, Γ scatter halved, CL/shocks/V6 preserved, 0 floored/limited. σ0.05 over-smooths (CL −1.3%, Γ_tip 0.021→0.028). **Caveat: OFF field M_max = 1.97, not the cached 1.47** (see reproducibility note). | Weak positive, but band-Mmax move (~6%) is within the run-to-run noise. |
| **C** | coarse **injected-kink** at fixed Γ (×3 warm-started) | baseline Γ0 → M 1.47 (0/0); halve the z/b0.80 station → **M 8.38, 3 floored/2 limited**; smooth(kink) → **M 1.40, 0/0**. | Heals — but this is a **tautology**: smoothing removes an *injected* kink; it does not show the *natural* spike is kink-driven. |
| **D** | medium σ0.03 in the continuation (full solve, 32 min) | **WORSE**: M_max 5.20→**10.75**, floored/limited 8/4→6/20, M>2 26→58 (wing 18→**50**), V6 1.73→1.82%. Trajectory blew up at the M0.80 level (Mmax 13.08). Smoothing couples stations → the per-station secant is invalid → density solve destabilised on the stiffer medium mesh. | ❌ |
| **E** | medium **fixed-Γ** (bypasses the secant; ×3 warm-started, ~2 min each) | baseline (cached Γ) reproduces **5.204 exactly** (path check ✓). smooth σ0.03: **\|dΓ/dz\| 1.81→0.26** (flatter than the physical coarse 0.66) yet **M_max 5.341, still 14 wing M>2 cells**. smooth σ0.05: **diverges** (M_max 15.5, NaN). | ❌ **decisive refutation** |

**E is the decisive test.** Flattening the spanwise Γ gradient **7×** — to
smoother-than-physical — leaves the wing-TE spike essentially unchanged
(M_max 5.20→5.34, 18→14 cells). Spanwise Γ smoothness is **not** the lever.
The A-step +0.63 correlation was **correlation, not causation**.

## Corrected root cause & next tests (guides new work)

**Root cause:** a **trailing-edge discretization singularity at fine
resolution** — the P1 element-constant gradient at the sharp zero-thickness TE
grows as the TE tet shrinks. Same P1-singularity family as the **G1.6** LE-Cp
gate. Independent of the wake/Kutta model.

**Fix routes to try (NOT spanwise-Γ smoothing):**
1. **Wing-TE cluster (dominant, 18 cells) — TE-element-level:** curved/
   isoparametric or locally-refined TE elements; the **P6** consistent/
   differentiable artificial-density flux; or **P7** Newton for the
   sharp-feature tail.
2. **Far-field cluster (8 cells incl. the M_max=5.20 cell) — taper the
   2D-vortex far-field correction toward the tip** (`constraints/dirichlet.py`).
   Independent, cheap; removes M_max but not sufficient for the `physical`
   gate on its own.

## Reproducibility caveats (read before interpreting any medium re-run)

- **Field M_max is run-to-run noisy.** The coarse medium-recipe solve gave OFF
  field M_max 1.97 where the committed cache says 1.47, while CL/shocks/V6/Γ
  were stable. The far-field/wake tail is only *engineering-converged* (drho
  floors ~0.18, not 1e-8), so the field-wide max wanders between runs; judge
  the **wing** (band Mmax, CL, shocks, V6), not the bare field M_max.
- The medium fixed-Γ solves (E) also floor at drho≈0.18 in 400 iters. A very
  long (n_picard≈1500) fixed-Γ σ0.03 run was NOT executed; but a 7× flatter Γ
  giving the same M_max at equal iterations is already conclusive.

## Meta-lesson

Do the **cheapest falsification test first.** E (medium fixed-Γ, bypassing the
secant) is what refuted the hypothesis; it was cheaper than D (the 32-min
continuation) and should have preceded it. C's "heal" was a tautology
(smoothing an *injected* kink) and did not deserve the confidence it was given.
