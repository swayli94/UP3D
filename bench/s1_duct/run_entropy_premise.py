"""GS1b.2 premise check, part Q1: how much stronger is the isentropic
mass-conserving jump than Rankine-Hugoniot?

That difference IS the entire travel of the entropy-correction lever: if the
isentropic jump is only a fraction of a percent off R-H at the shock strengths
where the fold sits (M1 ~ 1.3-1.4, measured in GS1b.1), then no implementation
of a non-isentropic density law can move the fold and route B is not worth
building.

Uses the shipped relations (pyfp3d/physics/isentropic.py) for the isentropic
branch, so this is a statement about THIS code, not about a textbook.

Outputs: results/gs1b_2_q1_jump.csv
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))

from pyfp3d.physics.isentropic import (GAMMA, critical_speed_squared,   # noqa
                                       density_isentropic, q2_at_mach)

M_INF = 0.80          # the normalisation of the shipped relations
M1S = (1.15, 1.20, 1.30, 1.35, 1.40, 1.50, 1.60)


def isentropic_jump(m1, m_inf=M_INF, gamma=GAMMA):
    """Downstream state of the mass-conserving ISENTROPIC jump: the other root
    of rho(q^2) q = const, exactly as the solver's own flux implies."""
    q2_1 = q2_at_mach(m1, m_inf, gamma)
    u1 = float(np.sqrt(q2_1))
    mdot = float(density_isentropic(q2_1, m_inf, gamma) * u1)
    u_star = float(np.sqrt(critical_speed_squared(m_inf, gamma)))
    lo, hi = 1e-12, u_star
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if float(density_isentropic(mid * mid, m_inf, gamma) * mid) < mdot:
            lo = mid
        else:
            hi = mid
    u2 = 0.5 * (lo + hi)
    return u1, u2, u1 / u2          # u1/u2 == rho2/rho1 by mass conservation


def rh_jump(m1, gamma=GAMMA):
    """Rankine-Hugoniot normal shock: rho2/rho1 and the total-pressure ratio."""
    m2 = m1 * m1
    rho_ratio = ((gamma + 1.0) * m2) / ((gamma - 1.0) * m2 + 2.0)
    p0_ratio = (((gamma + 1.0) * m2 / ((gamma - 1.0) * m2 + 2.0))
                ** (gamma / (gamma - 1.0))
                * ((gamma + 1.0) / (2.0 * gamma * m2 - (gamma - 1.0)))
                ** (1.0 / (gamma - 1.0)))
    return rho_ratio, p0_ratio


def corrected_jump(m1, sigma, m_inf=M_INF, gamma=GAMMA):
    """FP mass-conservation jump with the ENTROPY-CORRECTED downstream density
    rho_s = sigma * rho_isen(q), sigma = p02/p01. Returns (u2, rho2/rho1).

    ADDED 2026-07-29 to fix an exponent error of mine: the round file first
    quoted the density factor as sigma^(1/(gamma-1)), which applies the exponent
    TWICE. For a perfect gas with constant total enthalpy, rho0 = p0/(R*T0) and
    T0 is preserved across the shock, so rho02/rho01 = p02/p01 EXACTLY and the
    density factor is sigma itself. This routine checks that against R-H rather
    than trusting the algebra -- and it is exact to 5 decimals at every strength,
    which is also the 0-D unit test GS1b.3 must pass.
    """
    u1 = float(np.sqrt(q2_at_mach(m1, m_inf, gamma)))
    rho1 = float(density_isentropic(u1 * u1, m_inf, gamma))
    mdot = rho1 * u1
    lo, hi = 1e-6, u1
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if sigma * float(density_isentropic(mid * mid, m_inf, gamma)) * mid \
                < mdot:
            lo = mid
        else:
            hi = mid
    u2 = 0.5 * (lo + hi)
    rho2 = sigma * float(density_isentropic(u2 * u2, m_inf, gamma))
    return u2, rho2 / rho1


def main():
    rows = []
    print(f"{'M1':>6s} {'rho2/rho1 isen':>15s} {'rho2/rho1 R-H':>14s} "
          f"{'over-strength':>14s} {'sigma=p02/p01':>13s} "
          f"{'1-sigma':>8s} {'corrected':>10s} {'corr vs R-H':>12s}")
    for m1 in M1S:
        u1, u2, r_isen = isentropic_jump(m1)
        r_rh, p0 = rh_jump(m1)
        over = r_isen / r_rh - 1.0
        u2c, r_corr = corrected_jump(m1, p0)
        rows.append(dict(m1=m1, u1=round(u1, 6), u2_isentropic=round(u2, 6),
                         rho_ratio_isentropic=round(r_isen, 6),
                         rho_ratio_rh=round(r_rh, 6),
                         over_strength=round(over, 6),
                         sigma_p0_ratio=round(p0, 6),
                         density_reduction=round(1.0 - p0, 6),
                         u2_corrected=round(u2c, 6),
                         rho_ratio_corrected=round(r_corr, 6),
                         corrected_vs_rh=round(r_corr / r_rh - 1.0, 9)))
        print(f"{m1:6.2f} {r_isen:15.5f} {r_rh:14.5f} {100 * over:13.2f} % "
              f"{p0:13.5f} {100 * (1 - p0):7.2f}% {r_corr:10.5f} "
              f"{100 * (r_corr / r_rh - 1.0):11.4f} %")

    out = HERE / "results" / "gs1b_2_q1_jump.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}")
    band = [r for r in rows if 1.29 <= r["m1"] <= 1.41]
    print("\nQ1 verdict band (M1 = 1.30..1.40, where GS1b.1 found the fold):")
    for r in band:
        print(f"  M1 {r['m1']}: isentropic jump is "
              f"{100 * r['over_strength']:.2f} % stronger than R-H; "
              f"entropy factor sigma = p02/p01 = {r['sigma_p0_ratio']:.5f} "
              f"({100 * r['density_reduction']:.2f} % density reduction); "
              f"corrected jump vs R-H {100 * r['corrected_vs_rh']:+.5f} %")
    print("\n  => rho_s = sigma * rho_isen(q^2) with sigma = p02/p01 reproduces\n"
          "     Rankine-Hugoniot EXACTLY at every strength above (the 0-D unit\n"
          "     test for GS1b.3). NOT sigma^(1/(gamma-1)) -- that applies the\n"
          "     exponent twice; corrected round-file numbers accordingly.")


if __name__ == "__main__":
    main()
