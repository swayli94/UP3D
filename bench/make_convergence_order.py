"""Three-level cl_p convergence order for all six geometry x wake-path combinations.

This is what the xcoarse levels were added for. Before them, FIVE of the six combinations
had only two mesh levels, so no order was computable anywhere except conforming NACA -- the
one place a fine level already existed. Reads only committed CSVs; no solve.

Ladders differ per family BY NECESSITY, and the reason is recorded with each:
  NACA      xcoarse(0.040)/coarse(0.020)/medium(0.010)   uniform 2x; xcoarse needs the
            h_far clamp OFF or its first interval would not refine the far field at all.
  M6 wing   xcoarse_ss(0.060)/coarse_ss(0.030)/medium(0.015)  uniform 2x; the SHIPPED
            coarse is clamped (P13/M1b) and is off the refinement ray, so it cannot sit
            in a convergence ladder -- hence the _ss levels.
  wing-body xcoarse(0.044)/coarse(0.030)/medium(0.015)   ratios 1.47 then 2.00, because
            0.060 BREAKS the mesh (min dihedral 19.50 -> 0.70 deg; TIP_CAP_RADIUS = 0.0222
            is a fixed scale) and both paths must share h_wall.

Outputs (TRACKED): bench/gate_results/convergence_order.csv
"""
import csv, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "gate_results")
LAD = {("naca2.5d", "conforming"): ("xcoarse", "coarse", "medium", 0.040, 0.020, 0.010),
       ("naca2.5d", "level-set"): ("xcoarse", "coarse", "medium", 0.040, 0.020, 0.010),
       ("m6wing", "conforming"): ("xcoarse_ss", "coarse_ss", "medium", 0.060, 0.030, 0.015),
       ("m6wing", "level-set"): ("xcoarse_ss", "coarse_ss", "medium", 0.060, 0.030, 0.015),
       ("wingbody", "conforming"): ("xcoarse", "coarse", "medium", 0.044, 0.030, 0.015),
       ("wingbody", "level-set"): ("xcoarse", "coarse", "medium", 0.044, 0.030, 0.015)}


#: ★ THE PRIMARY READING IS `R = delta2/delta1`, NOT `order_p`. Read this before using p.
#:
#: p is a single-power fit to three points, and it is NOT reliable here -- not because of
#: noise (the deltas are 5e-3..1.6e-2 against a ~1e-5..3e-5 frozen-selection scatter in cl,
#: so 100-1000x above it) but because of MODEL FORM: if the error carries two components of
#: different order, a one-term fit blends them, and different functionals blend them
#: differently. Measured signature of exactly that: on the wing-body, cl_p gives p = 0.20
#: while cl_KJ gives p = 0.45 on the same states.
#:
#: p also actively MISLEADS the ranking. It reported wing-body level-set at 0.72 against
#: conforming at 0.20 -- suggesting one is 3.5x better -- when in R terms both are >= 1,
#: i.e. both FALSIFIED, with no better or worse between them.
#:
#: R needs no power law to be useful:
#:   * R is bounded above by the LADDER'S OWN p -> 0 limit,
#:       R_max = ln(h_m/h_c) / ln(h_c/h_x),
#:     so R >= R_max FALSIFIES convergence at any positive order without identifying it.
#:     ★ R_max IS LADDER-SPECIFIC: 1.000 on a uniform 2x ladder but 1.810 on the
#:     wing-body's 1.47/2.00 pair. I first wrote 1.0 as a universal threshold and thereby
#:     called the wing-body "falsified" when its R = 1.24-1.69 sits BELOW its own 1.810 --
#:     the same "criterion calibrated on one baseline, applied to a non-comparable one"
#:     mistake this file already warns about for p. On a non-uniform ladder growing
#:     differences are COMPATIBLE with a small positive order, because the second interval
#:     is the bigger refinement step.
#:   * on a uniform 2x ladder, first order means R = 0.5 exactly -- "the difference halves
#:     when h halves" is then a direct reading, not a fit.
#:   * a SIGN FLIP between delta1 and delta2 means the sequence is non-monotone: the
#:     coarsest level is outside the asymptotic range and no order exists at all.
def solve_p(d1, d2, hx, hc, hm):
    """cl(h) = cl* + C h^p on three levels. SECONDARY -- see the note above. None when no
    positive p fits, which happens for a non-monotone sequence."""
    if d1 == 0:
        return None
    R, rx, rm = d2 / d1, hx / hc, hm / hc
    f = lambda p: (rm ** p - 1) / (1 - rx ** p) - R
    lo, hi = 0.01, 8.0
    if f(lo) * f(hi) > 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if f(lo) * f(mid) <= 0 else (mid, hi)
    return 0.5 * (lo + hi)


def main():
    rows = [r for r in csv.DictReader(open(os.path.join(OUT, "capability_matrix.csv")))
            if r["status"].startswith("CLEAN")]

    def cl(path, geom, lv, m):
        for r in rows:
            if (r["path"] == path and r["geom"] == geom and r["level"] == lv
                    and abs(float(r["m_inf"]) - m) < 1e-9):
                return float(r["cl_p"])

    out = []
    for (geom, path), (lx, lc, lm, hx, hc, hm) in LAD.items():
        for m in (0.50, 0.60, 0.65, 0.70, 0.75, 0.78, 0.80, 0.84):
            a, b, c = (cl(path, geom, lx, m), cl(path, geom, lc, m),
                       cl(path, geom, lm, m))
            if None in (a, b, c):
                continue
            d1, d2 = b - a, c - b
            p = solve_p(d1, d2, hx, hc, hm)
            R = d2 / d1 if d1 else float("nan")
            flip = (d1 > 0) != (d2 > 0)
            #: the ladder's own first-order value and its p -> 0 ceiling, both computed
            #: from the h's rather than hard-coded -- see the R_max warning in the header.
            r_first = (hm ** 1 - hc ** 1) / (hc ** 1 - hx ** 1)
            r_max = math.log(hm / hc) / math.log(hc / hx)
            verdict = ("non-monotone: coarsest level outside the asymptotic range"
                       if flip else
                       f"FALSIFIED: R >= this ladder's p->0 ceiling {r_max:.3f}"
                       if R >= r_max else
                       "first-order signature (R = this ladder's p=1 value)"
                       if abs(R - r_first) < 0.05 * abs(r_first) else
                       "converging, order far below first" if R > r_first
                       else "converging, order at or above first")
            out.append(dict(geom=geom, path=path, m_inf=m, h_xcoarse=hx, h_coarse=hc,
                            h_medium=hm, cl_xcoarse=a, cl_coarse=b, cl_medium=c,
                            delta1=round(d1, 8), delta2=round(d2, 8),
                            R=round(R, 4), R_first_order=round(r_first, 4),
                            R_falsify_ceiling=round(r_max, 4), verdict=verdict,
                            order_p=(None if p is None else round(p, 4)),
                            monotone=abs(d2) < abs(d1),
                            note=("p is SECONDARY -- see the header note" )))
    p = os.path.join(OUT, "convergence_order.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)
    print(f"wrote {p} ({len(out)} rows)")
    print("\n=== reading ===")
    for (geom, path) in LAD:
        sel = [r for r in out if r["geom"] == geom and r["path"] == path]
        rs = [r["R"] for r in sel]
        print(f"  {geom:10s} {path:11s} R = "
              f"{', '.join(f'{x:.3f}' for x in rs)}   "
              f"(p=1 -> {sel[0]['R_first_order']}, falsified at >= "
              f"{sel[0]['R_falsify_ceiling']})")
        for v in dict.fromkeys(r["verdict"] for r in sel):
            print(f"                            -> {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
