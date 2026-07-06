# NACA0012 M∞ = 0.80 α = 1.25° transonic reference (gate G4.1)

The canonical transonic benchmark (AGARD AR-138/AR-211 lineage; also the
GAMM 1979 workshop case): a strong upper-surface shock and a weak
lower-surface shock.

## Reference values (`shock_reference.csv`)

| quantity | value | provenance |
|---|---|---|
| upper shock x/c (Euler anchor) | 0.60–0.63 | canonical Euler solutions; e.g. [arXiv:2406.07441](https://arxiv.org/pdf/2406.07441) states the two shocks sit "at approximately 60% and 35% of the chord"; the same positions are visible in Jameson & Ou, *50 years of transonic aircraft design*, Prog. Aerosp. Sci. (2011), Fig. 2(a) ([PDF](http://aero-comlab.stanford.edu/Papers/aj_PAS_2011.pdf)) |
| lower shock x/c (Euler anchor) | 0.33–0.38 | same sources |
| FP-vs-Euler systematic shift | conservative FP shocks sit AT-to-AFT of Euler with slightly stronger jumps (isentropic shock, no entropy/vorticity production) | Holst, *Transonic flow computations using nonlinear potential methods*, Prog. Aerosp. Sci. 36 (2000), NTRS [20020078395](https://ntrs.nasa.gov/api/citations/20020078395/downloads/20020078395.pdf), §FP-vs-Euler discussion; Salas et al. comparisons cited therein |
| **G4.1 gate value: upper shock x/c = 0.62 ± 0.03** | Euler anchor midpoint 0.615 + 0–2% chord conservative-FP aft shift, rounded | derived from the two rows above |
| lower shock x/c (reported, not gated) | 0.35 ± 0.04 | Euler anchor; the lower shock is weak and its position is more mesh-sensitive |

## Provenance discipline note

A digitized *full-potential* Cp table for this exact case was not
retrievable from an open source at generation time (Holst's review
discusses the M = 0.83 nonuniqueness band and 3D cases, not this 2D
case's table; AGARD AR-211 itself is not openly downloadable). The gate
reference is therefore an **Euler anchor + documented FP shift band**,
recorded here with sources, rather than a fabricated "FP table". If a
proper digitized FP dataset is obtained later, replace
`shock_reference.csv` and record the change here (hard rule 6 applies to
these files once committed).

Measured for context (2026-07-07, pyFP3D): coarse-mesh (16.4k tets)
converged solution — upper shock 0.599, lower 0.362, both monotone, no
expansion shock, cl ≈ 0.334, M_max = 1.36.
