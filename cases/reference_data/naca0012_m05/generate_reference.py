"""
Reference-data generator for gate G3.2: NACA0012 at M_inf = 0.5, subsonic
compressible, from the SAME Hess-Smith panel solution the P2 gates use
(../naca0012_incompressible/generate_panel_reference.py), corrected for
compressibility two independent ways:

  - Prandtl-Glauert (PG):  Cp_c = Cp_i / beta,  beta = sqrt(1 - M^2).
    First-order in the perturbation; known to under-predict suction-peak
    amplification on thick sections.
  - Karman-Tsien (KT):     Cp_c = Cp_i / (beta + (M^2/(1+beta)) * Cp_i/2).
    Hodograph-based; the standard subcritical engineering correction and
    the more accurate of the two for a 12%-thick section at M = 0.5
    (see e.g. Anderson, "Modern Compressible Flow", ch. 9). The gate
    compares against the KT value; the PG value is recorded for context.

cl is obtained by integrating the corrected Cp over the SAME panel
geometry and normals as the incompressible reference (dC = -Cp n dS),
so the only change vs the P2 reference is the pointwise Cp mapping --
method-vs-method on identical geometry, per the project's reference-data
discipline (roadmap "Reference data discipline").

Validity: M_inf = 0.5, alpha = 2 deg is comfortably subcritical (the
solver measures max local M^2 ~ 0.53), which is where both corrections
are meant to apply; KT-vs-PG spread is recorded in the CSV as the
correction-model uncertainty band.

Writes: cl_reference.csv, cp_alpha2_m05.csv.
Run:    python generate_reference.py
"""

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "naca0012_incompressible"))
from generate_panel_reference import _panel_geometry, solve_hess_smith  # noqa: E402

M_INF = 0.5
ALPHAS = [0.0, 2.0]
N_HALF = 401  # ~800 panels, same as the converged incompressible reference


def karman_tsien(cp_i: np.ndarray, m_inf: float) -> np.ndarray:
    beta = np.sqrt(1.0 - m_inf**2)
    return cp_i / (beta + (m_inf**2 / (1.0 + beta)) * cp_i / 2.0)


def prandtl_glauert(cp_i: np.ndarray, m_inf: float) -> np.ndarray:
    return cp_i / np.sqrt(1.0 - m_inf**2)


def corrected_cl(alpha_deg: float, n_half: int, m_inf: float):
    """cl from integrating the corrected panel Cp over the same geometry."""
    r = solve_hess_smith(alpha_deg, n_half)
    _, _, _, n_hat, _, _, length = _panel_geometry(n_half)
    a = np.deg2rad(alpha_deg)
    lift_dir = np.array([-np.sin(a), np.cos(a)])

    out = {"cl_incompressible": r["cl"], "cp_i": r["cp"], "xm": r["xm"],
           "ym": r["ym"]}
    for name, fn in (("pg", prandtl_glauert), ("kt", karman_tsien)):
        cp_c = fn(r["cp"], m_inf)
        cf = -(cp_c * length) @ n_hat
        out[f"cl_{name}"] = float(cf @ lift_dir)
        out[f"cp_{name}"] = cp_c
    return out


def main():
    out_dir = Path(__file__).parent

    rows = []
    for a in ALPHAS:
        r = corrected_cl(a, N_HALF, M_INF)
        rows.append([a, f"{r['cl_incompressible']:.6f}",
                     f"{r['cl_pg']:.6f}", f"{r['cl_kt']:.6f}"])
        print(f"alpha={a:4.1f}  cl_inc={r['cl_incompressible']:.6f}  "
              f"cl_PG={r['cl_pg']:.6f}  cl_KT={r['cl_kt']:.6f}  "
              f"(KT-PG spread {100*(r['cl_kt']/max(r['cl_pg'],1e-30)-1) if a else 0:.2f}%)")

    with open(out_dir / "cl_reference.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["alpha_deg", "cl_incompressible", "cl_pg", "cl_kt"])
        w.writerows(rows)

    r2 = corrected_cl(2.0, N_HALF, M_INF)
    upper = r2["ym"] > 0
    with open(out_dir / "cp_alpha2_m05.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_c", "cp_incompressible", "cp_pg", "cp_kt", "surface"])
        for xm, ci, cp, ck, up in zip(r2["xm"], r2["cp_i"], r2["cp_pg"],
                                      r2["cp_kt"], upper):
            w.writerow([f"{xm:.6f}", f"{ci:.6f}", f"{cp:.6f}", f"{ck:.6f}",
                        "upper" if up else "lower"])
    print(f"wrote {out_dir}/cl_reference.csv, cp_alpha2_m05.csv")


if __name__ == "__main__":
    main()
