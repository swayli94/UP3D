"""
Shock monitors for transonic surface solutions (roadmap P4 deliverable):
sonic-crossing shock detection on a sectional Cp(x/c) curve, shock
sharpness (cell count), monotonicity across the jump, and an
expansion-shock detector.

All monitors work on the OUTPUT of post/section_cut.wall_cp_curve (the
triangle-wise sectional Cp), against the exact isentropic critical
pressure coefficient Cp* = Cp(q*^2) (physics/isentropic.py) -- Cp < Cp*
means locally supersonic.

Reference: design.md Sec 3 properties ("monotone shock capture over 2-3
cells; no expansion shocks"), roadmap gate G4.1.
"""

from typing import Dict

import numpy as np

from pyfp3d.physics.isentropic import (
    critical_speed_squared,
    pressure_coefficient,
)


def cp_critical(m_inf: float, gamma: float = 1.4) -> float:
    """Cp at the sonic condition (exact isentropic): Cp* = Cp(q*^2)."""
    return float(pressure_coefficient(critical_speed_squared(m_inf, gamma),
                                      m_inf, gamma))


def _dedupe(x: np.ndarray, cp: np.ndarray):
    """Average Cp of points at (numerically) the same x/c -- the quasi-2D
    extrusion produces one point per z-layer triangle at each station."""
    order = np.argsort(x)
    x, cp = x[order], cp[order]
    xu, inv = np.unique(np.round(x, 6), return_inverse=True)
    cpu = np.zeros(len(xu))
    cnt = np.zeros(len(xu))
    np.add.at(cpu, inv, cp)
    np.add.at(cnt, inv, 1.0)
    return xu, cpu / cnt


def shock_metrics(x: np.ndarray, cp: np.ndarray, m_inf: float,
                  gamma: float = 1.4) -> Dict[str, object]:
    """
    Analyze one surface's Cp(x/c) curve for the terminating shock.

    Detection: the LAST supersonic->subsonic sonic crossing (Cp rising
    through Cp*). The shock position is the Cp* crossing, linearly
    interpolated; the shock width is measured between the last point
    below a supersonic plateau margin and the first point of subsonic
    recovery, counted in surface stations (cells).

    Returns dict:
        has_shock: bool
        x_shock: Cp* crossing x/c (nan if no shock)
        n_cells: stations spanned by the pressure jump
        monotone: Cp non-decreasing across the jump (within 1% of the
            jump magnitude -- measurement noise on triangle-wise data)
        n_supersonic: supersonic station count
        expansion_shock: True if the flow re-accelerates to supersonic
            DOWNSTREAM of the detected shock (forbidden by (3.2)'s
            compression-only dissipation)
        cp_min: suction peak
    """
    x, cp = _dedupe(np.asarray(x, float), np.asarray(cp, float))
    cps = cp_critical(m_inf, gamma)
    sup = cp < cps

    out = {"has_shock": False, "x_shock": float("nan"), "n_cells": 0,
           "monotone": True, "n_supersonic": int(sup.sum()),
           "expansion_shock": False, "cp_min": float(cp.min())}
    if not sup.any():
        return out

    # Last supersonic->subsonic transition = the terminating shock.
    trans = np.where(sup[:-1] & ~sup[1:])[0]
    if len(trans) == 0:
        return out
    i = int(trans[-1])
    t = (cps - cp[i]) / (cp[i + 1] - cp[i])
    out["has_shock"] = True
    out["x_shock"] = float(x[i] + t * (x[i + 1] - x[i]))

    # Jump extent: from the last supersonic station i, walk downstream
    # while the per-station rise stays comparable to the first jump
    # increment (>= 25%); the gentle recovery slope ends the count.
    first_step = cp[i + 1] - cp[i]
    k = i + 1
    while k + 1 < len(cp) and cp[k + 1] - cp[k] > 0.25 * first_step:
        k += 1
    jump_scale = max(abs(cp[k] - cp[i]), 1e-6)
    out["n_cells"] = int(k - i)

    # Monotone across the jump: from the last plateau point to recovery.
    seg = cp[i:k + 1]
    dips = np.diff(seg) < -0.01 * jump_scale
    out["monotone"] = not bool(dips.any())

    # Expansion shock: any supersonic station strictly downstream of the
    # recovery point k.
    out["expansion_shock"] = bool(sup[k + 1:].any())
    return out


def shock_report(curve: Dict[str, np.ndarray], m_inf: float,
                 gamma: float = 1.4) -> Dict[str, object]:
    """shock_metrics for both surfaces of a wall_cp_curve() result."""
    return {
        "upper": shock_metrics(curve["x_upper"], curve["cp_upper"], m_inf,
                               gamma),
        "lower": shock_metrics(curve["x_lower"], curve["cp_lower"], m_inf,
                               gamma),
        "cp_critical": cp_critical(m_inf, gamma),
    }
