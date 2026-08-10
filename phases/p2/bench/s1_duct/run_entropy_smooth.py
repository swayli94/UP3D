"""RETIRED 2026-07-31 (GS1b.10). This measured the GS1b.4 SMOOTHING widths
(eps_m / f_lo / f_hi) on the donor-chain sigma, and both are gone: GS1b.4 was a
measured negative (S1 FAIL, and the widths behaved as a fitting parameter), and
GS1b.10 replaced the chain entirely with an FV upwind transport plus an additive
entropy production density -- which has no smoothing width to sweep.

Kept for the record because its CSV (results/gs1b_4_smooth.csv) is the evidence
behind that negative. The stability question it asked is now asked by
run_entropy_stability.py against the new construction.
"""

raise SystemExit(__doc__)
