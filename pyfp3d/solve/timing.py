"""
Cost accounting shared by the four nonlinear drivers (Track A / A1).

The drivers -- conforming Picard/Newton (`picard.py`, `newton.py`) and
level-set Picard/Newton (`picard_ls.py`, `newton_ls.py`) -- each report a
`timings` dict on the SAME schema, so their wall clock can be compared
phase by phase:

    seed      Picard warm start inside a Newton driver (0.0 for the Picard
              drivers, which have no seed).
    assembly  operator matrix assembly + constraint update/slicing + the
              per-solve right-hand side build.
    precond   AMG / ILU construction, splu factorization.
    linsolve  CG / GMRES / spsolve / LU back-solve (Woodbury included).
    residual  velocities + density + upwind + residual assembly.
    kutta     Kutta targets / TE pressure residual + the Gamma update.
    wall      perf_counter across the whole driver.
    other     wall - sum(the phases above); the unaccounted remainder.

`other` is the honesty term: it is what the phase breakdown does NOT
explain, and gate GA1.1 asserts it stays under 5% of wall. Phases are
exclusive -- a region timed as `residual` must not also contain an
assembly (the level-set residual does, and `newton_ls._System` splits it
out explicitly).

Reference: docs/roadmap/track_a.md A1.
"""

import time
from typing import Dict, Iterable, List

# Ordered so a stacked bar chart reads seed -> assembly -> ... -> other.
PHASES = ("seed", "assembly", "precond", "linsolve", "residual", "kutta")


def new_timings() -> Dict[str, float]:
    """A zeroed timings dict on the canonical schema."""
    t = {p: 0.0 for p in PHASES}
    t["wall"] = 0.0
    t["other"] = 0.0
    return t


class phase:
    """Accumulate the wall time of a `with` block into `timings[key]`.

    A fresh instance per use (so nesting the same key is impossible by
    construction). The cost is one small allocation and two perf_counter
    calls per OUTER iteration -- nanoseconds against an assembly.
    """

    __slots__ = ("_t", "_k", "_t0")

    def __init__(self, timings: Dict[str, float], key: str):
        self._t = timings
        self._k = key

    def __enter__(self):
        self._t0 = time.perf_counter()
        return self

    def __exit__(self, *exc):
        self._t[self._k] += time.perf_counter() - self._t0
        return False


def finalize(timings: Dict[str, float], wall: float) -> Dict[str, float]:
    """Stamp `wall` and derive `other`. Call once, on the way out."""
    timings["wall"] = float(wall)
    timings["other"] = float(wall - sum(timings[p] for p in PHASES))
    return timings


def sum_timings(parts: Iterable[Dict[str, float]]) -> Dict[str, float]:
    """Add up per-level timings dicts (the Mach-ramp total).

    The ramp wrappers report this as `timings_total`, alongside the
    final-level `timings` they have always returned -- reading the latter
    as if it were the whole ramp is the "timings dict = final level only"
    footgun (P8 record).
    """
    total = new_timings()
    for p in parts:
        for k in total:
            total[k] += float(p.get(k, 0.0))
    return total


def step_delta(timings: Dict[str, float],
               prev: Dict[str, float]) -> Dict[str, float]:
    """Per-phase time spent since the previous step, for `step_records`."""
    return {f"t_{p}": timings[p] - prev.get(p, 0.0) for p in PHASES}


def snapshot(timings: Dict[str, float]) -> Dict[str, float]:
    """Copy of the running phase counters (the `prev` of `step_delta`)."""
    return {p: timings[p] for p in PHASES}


def records_to_columns(records: List[Dict]) -> Dict[str, list]:
    """Transpose a `step_records` list into column arrays for CSV export."""
    keys: List[str] = []
    for r in records:
        for k in r:
            if k not in keys:
                keys.append(k)
    return {k: [r.get(k) for r in records] for k in keys}
