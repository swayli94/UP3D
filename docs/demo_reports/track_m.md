# Phase demo report — Track M (mesh generation)

> Split verbatim from `docs/demo_report.md` on 2026-07-15 (content unchanged;
> only this header was added; sections keep their original chronological order).
> Scope, reproduce instructions and the honesty/evidence rule: see the
> [demo_report.md](../demo_report.md) index. Roadmap gates: [tracks/](../tracks/).

## M0 — quasi-2D meshing pipeline (closed; acceptance link = G2.5)

**Purpose.** Show the mesh-side evidence for M0: the pipeline (vanilla-Gmsh
planar mesh → single-layer extrusion → globally consistent min-global-index
prism→3-tet split) produces topologically sound, refinable meshes that the
solver actually converges on. M0's formal acceptance was G2.5 — the Track M
↔ Track P link demonstrated in the P2 demo above.

**Case setup.** All seven committed meshes in `cases/meshes/`
(naca0012_2.5d coarse/medium, cylinder_2.5d coarse/medium/fine,
sphere_shell coarse/medium); end-to-end solve on the cylinder family
against the analytic Cp = 1 − 4 sin²θ.

**Key figures.**

![mesh gallery](../cases/demo/m0_meshgen/results/mesh_gallery.png)
![topology asserts](../cases/demo/m0_meshgen/results/topology_asserts.png)
![mesh quality](../cases/demo/m0_meshgen/results/mesh_quality.png)
![cylinder Cp convergence](../cases/demo/m0_meshgen/results/cylinder_cp_convergence.png)

**Measured results.**

| Check | Measured | Criterion |
|---|---|---|
| topology asserts (tags, quad-split consistency, wake cut) on all 7 meshes | all pass | hard rule 7 |
| min tet volume across families | 1.3e-6 > 0 | no degenerate/inverted tets |
| isotropic sphere family max aspect ratio | 3.5 | < 5 (quasi-2D far-field anisotropy is by design) |
| cylinder solve residual, all 3 levels | ≤ 3.1e-11 | < 1e-8 |
| max wall Cp error coarse → medium → fine | 9.1% → 4.5% → 2.2% | monotone, ≥ 25%/level |
| Cp error slope vs h | 1.02 | ~O(h), documented curved-wall limit |

**Conclusion & analysis.** Every committed mesh passes the full topology
assert battery (agent-rules hard rule 7), including the wake-cut asserts on
the NACA family — the wake is a single conforming interior sheet from TE to
far field, and the prism split is globally consistent so no lateral quad is
cracked. Element quality is deliberate rather than accidental: the isotropic
sphere family stays under aspect 3.5, while the quasi-2D families carry
their large aspect ratios exactly where the single-layer design puts them
(far field: in-plane coarsening at fixed dz). The pipeline's meshes are not
just valid but *useful*: the cylinder case solves to 1e-11 residuals at
every level and its wall-Cp error falls monotonically at the expected O(h)
(slope 1.02) — first order because of the same flat-facet curved-wall
recovery limit root-caused at G1.6, i.e. a documented solver-side limit, not
a meshing defect. Combined with G2.3–G2.5 running on these meshes (P2 demo),
M0's deliverable is demonstrated end to end.

---

## M1 — swept/tapered wing meshing, ONERA M6 (closed; consumed by P5)

**Purpose.** Show the mesh-side evidence for M1: a scripted, refinable
ONERA M6 half-wing tet mesh whose chord-plane wake sheet — swept from the
sharp (foilmod zero-thickness) TE, ending exactly at the tip, reaching the
spherical far field at 15 MAC — is ingested by the P2 solver preprocessor
with the topology asserts green. The new mesh-side machinery is
`pyfp3d/meshgen/wing3d.py` (OCC ruled loft + `occ.fragment` +
`mesh.embed`); the new solver-side machinery is wake_cut.py's handling of
a swept TE (per-node stations, off-plane Kutta-probe fallback) and of the
sheet's interior FREE edge at the tip (single-valued nodes ⇒ Γ(tip) = 0
discretely).

**Case setup.** The `cases/meshes/onera_m6` family (coarse 55.5k /
medium 350.7k tets; fine 2.513M validated at generation time). The .msh
files are large and gitignored — regenerate coarse+medium with
`generate_onera_m6.py` (~30 s) before running this demo; the committed
per-level stats CSVs and inspection PNGs are the persistent evidence.
Solver axis convention: chord x, lift y, span z.

**Key figures.**

![wing + wake gallery](../cases/demo/m1_wing_mesh/results/wing_wake_gallery.png)
![tip cut planes](../cases/demo/m1_wing_mesh/results/tip_cut_planes.png)
![mesh quality](../cases/demo/m1_wing_mesh/results/mesh_quality.png)

**Measured results.**

| Check | Measured | Criterion |
|---|---|---|
| tags + P2 topology asserts through cut_wake, coarse & medium | pass | M1 gate "same asserts" |
| per-node TE stations on the swept TE | 83 (coarse) / 166 (medium) | == n_TE_nodes |
| tip free-edge nodes single-valued, at z = b | 106 / 208, none duplicated | wake-tip semantics |
| wake-tip closure: tip edge one open chain from the exact tip TE corner | pass (both levels) | no cracks / self-intersections |
| Kutta probe pairs found on the unstructured TE | 83 / 166, y>0 upper, y<0 lower | design.md (4.4) fallback |
| min dihedral coarse/medium/fine | 7.5° / 11.0° / 3.5° | ≥ 2° |
| max aspect ratio coarse/medium/fine | 9.3 / 6.9 / 6.5 | ≤ 60 |
| refinement ladder (one h_wall parameter, 2×) | 55.5k → 350.7k → 2513k tets | monotone ~2³/level |
| freestream residual on the CUT coarse mesh | 4.3e-14 | < 1e-10 (G2.1 analogue) |

**Conclusion & analysis.** The M1 gate items are all measured green: the
solver preprocessor ingests the family unchanged (same read_mesh/cut_wake
call as the 2.5D cases), the quality report is comfortably inside bounds
on all three levels, and the family is one script with one parameter. Two
findings worth recording: (1) for a sheet that ends *inside* the domain,
`occ.fragment` alone does not make the tet mesh conform — it stitches the
shared TE edge and the boundary trims, but `gmsh.model.mesh.embed` must
still be called on the trimmed sheet face; (2) the sheet's tip edge is an
interior free edge whose node stars are NOT split by the sheet, so the
duplication map must exclude them — which is also the physically correct
discrete statement Γ(tip) = 0 (the trailing jump vanishes at the tip).
Both are documented in the wake_cut.py module docstring; the free-edge
path is exactly inert on the quasi-2D meshes (their sheets have no free
edges), which the unchanged P2/M0 test battery confirms. Mesh sizes are
runtime-driven per the P4 lesson (solver wall time is the binding
constraint): coarse is the P5 development mesh, medium the gate mesh.

---

## Track M / M5 — the tip cap, rounded (`cases/demo/m5_round_tip/`, 9/9 PASS, 2026-07-13)

G13.3 named the defect; this is the fix. `meshgen/wing3d.py` grew
**`tip_cap="round"`** (default `"flat"` ⇒ every existing family bit-identical, so
the P5 / P8-G8.2 / B7 / M1 locks are untouched), which closes the wing with the
**half body of revolution swept by the tip section about its own chord line** —
`{√(y² + (z−B_SEMI)²) ≤ t(x)}`, with `t(x)` the tip section's local half
thickness. It is an OCC revolve of the tip section's upper half-face about an
edge *of* that face (the degenerate-at-the-axis kind, like a sphere from a half
disc), fused onto the loft.

★ **Why this construction and not a fillet, a loft, or curved elements.** The cap
radius **vanishes at the LE and the TE**, so the cap degenerates to a *point* at
each — which means the **TE line, the wake sheet, the tip TE corner, the Kutta
stations and B_SEMI are all unchanged**. Only the tip *wall* moves. That is what
makes the A/B against M1 controlled, and it is why the fix needed **no solver
code at all**. (A constant-radius fillet is geometrically impossible here — the
section thickness goes to zero at the TE — and P11's isoparametric *elements*
cannot regularize a sharp *edge* in the first place.) The revolved solid was
verified analytically and asserted at generation to protrude past the wing by at
most **4 × 10⁻⁶ m**, three orders below the finest edge size, and never to reach
aft of the local TE — so it can neither cut the wake sheet nor spawn slivers.

★ **The gate metric is the SEAM CREASE ANGLE, because the solver never sees the
CAD.** It only ever sees the triangulation, so "the geometry is round" is not the
claim that has to be true — "the *mesh* has no edge" is. Measure the turning
angle between the outward normals of adjacent wall triangles across the
tip-section seam at `z = B_SEMI` (the locus that *is* the sharp edge in M1), away
from the LE and TE — both sharp **by design** in either family, the TE because it
carries the Kutta condition (`post/surface.py::wall_crease_angles`):

| seam crease (p99) | coarse | medium | fine | exponent q |
|---|---|---|---|---|
| **flat cap (M1)** | **91.9°** | **92.1°** | **92.1°** | **−0.00** — does not move |
| **round cap (M5)** | 46.8° | 24.5° | **13.7°** | **−0.92** — O(h) |

(round cap, *max* rather than p99: 46.8 → 25.0 → 18.1°, q = −0.68 — see below.)

That contrast is the whole gate. A sharp convex edge creases by its own turning
angle: halving `h` **resolves** it better and removes nothing, which is exactly
why refinement made the flow singularity *worse* (p = +0.321). Facets
approximating a **smooth** surface crease by O(h · curvature) and halve when `h`
halves — the discrete statement of "there is no edge in the limit".

Two honesty notes on the metric. (i) The flat cap's **median** seam crease is
**0.00°**: the rest of the seam is already smooth, so this metric is reading the
edge and nothing else — it is not a diffuse mesh-quality number. (ii) The round
cap's **max** decays more slowly (q = −0.68) than its p99 (−0.92), because the
max is a single worst facet at the *thin* end of the seam window, where the cap
radius is smallest and the local facet size is set by `h_edge` (the TE/LE
refinement) rather than by `h_tip`. Both **decay**, which is the claim; the flat
cap's decays not at all.

**★ `h_tip = 0.25 · h_wall` is load-bearing, not a tuning knob.** The cap radius
is only `TIP_CAP_RADIUS` = 22 mm (1.9 % of the semi-span). At `h_wall` the coarse
cap would be about *one element* wide — i.e. the mesh would quietly discretize the
rounded cap back into a flat one and the geometry change would do nothing. `h_tip`
scales with the level, so the cap costs the same *fraction* at every level and the
refinement ray is preserved.

Other measured items (demo `run_demo.py`, 9/9; `tests/test_m5_round_tip.py`, 19):
cap resolved to its apex (wall `z_max` 1.218465 / 1.218466 vs the analytic apex
1.218467); quality inside the M1 bounds (min dihedral 4.05° / 5.18°, max aspect
16.2 / 6.4); `cut_wake` keeps the M1 semantics (85 TE stations; the tip TE corner
is still a free-edge node ⇒ Γ(tip) = 0 discretely); G2.1 freestream on the cut
mesh < 1e-10; **tip TE corner offset exactly 0.0** and the wake sheet still in the
chord plane, ending at B_SEMI. **Cost ×1.29 / ×1.28 tets** (59,359 / 448,197 vs
M1's `coarse_ss` 46,067 / `medium` 350,718) — level-independent, as a self-similar
ladder requires. Note the comparison is against M1's `coarse_ss`, not its shipped
`coarse`: the latter still carries the M1b `h_far` clamp and would report a
spurious ×1.07.

The **flow** consequence is P13/G13.3's, measured on the round ladder below.

---

