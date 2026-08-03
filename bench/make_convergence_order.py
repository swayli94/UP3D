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
import csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "gate_results")
LAD = {("naca2.5d", "conforming"): ("xcoarse", "coarse", "medium", 0.040, 0.020, 0.010),
       ("naca2.5d", "level-set"): ("xcoarse", "coarse", "medium", 0.040, 0.020, 0.010),
       ("m6wing", "conforming"): ("xcoarse_ss", "coarse_ss", "medium", 0.060, 0.030, 0.015),
       ("m6wing", "level-set"): ("xcoarse_ss", "coarse_ss", "medium", 0.060, 0.030, 0.015),
       ("wingbody", "conforming"): ("xcoarse", "coarse", "medium", 0.044, 0.030, 0.015),
       ("wingbody", "level-set"): ("xcoarse", "coarse", "medium", 0.044, 0.030, 0.015)}


def solve_p(d1, d2, hx, hc, hm):
    """cl(h) = cl* + C h^p, fitted on three levels. None when no positive p fits --
    which happens for a NON-MONOTONE sequence, and that is a result (the coarsest level is
    outside the asymptotic range), not a missing number."""
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
            out.append(dict(geom=geom, path=path, m_inf=m, h_xcoarse=hx, h_coarse=hc,
                            h_medium=hm, cl_xcoarse=a, cl_coarse=b, cl_medium=c,
                            delta1=round(d1, 8), delta2=round(d2, 8),
                            order_p=(None if p is None else round(p, 4)),
                            monotone=abs(d2) < abs(d1),
                            note=("non-monotone: coarsest level outside the asymptotic "
                                  "range" if p is None else "")))
    p = os.path.join(OUT, "convergence_order.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)
    print(f"wrote {p} ({len(out)} rows)")
    print("\n=== reading ===")
    for (geom, path) in LAD:
        sel = [r for r in out if r["geom"] == geom and r["path"] == path]
        ps = [r["order_p"] for r in sel if r["order_p"] is not None]
        print(f"  {geom:10s} {path:11s} {len(sel)} Mach pts, "
              f"p = {', '.join(f'{x:.2f}' for x in ps) if ps else 'not fittable'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
