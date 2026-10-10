"""Shared tolerance assertions for the test suite.

Phase two GS0.2, decision D1 (phases/p2/docs/dev_phase_two/roadmap.md §7):

  * permanent tests assert a RELATIVE tolerance of 1e-12 -- four to six orders
    tighter than any physically meaningful change (0.1 % in lift, 0.001 chord
    in shock position are all >= 1e-4), yet immune to last-bit rounding;
  * `np.array_equal` (bit identity) is reserved for the one claim it can
    actually support: "this code did not execute a different arithmetic path",
    i.e. the SAME array object/buffer is returned (e.g. nu == 0 makes
    `rho_tilde` hand back `rho`).

    ★★ 2026-09-06 (gate audit W0.6 / H28): this bullet used to end "Mark those
    `@pytest.mark.identity`." **That marker never existed** -- zero uses in
    `tests/`, and never registered in `pyproject.toml`, so the sentence
    described a mechanism the repo does not have. The DISCIPLINE above is real
    and stays; the PROMISE is deleted rather than implemented, because
    implementing it means classifying **145 `array_equal` call sites across 31
    files** one by one, which is a judgement-heavy task of its own and not a
    tidy-up. Registering an unused marker would only move the falsehood from
    "the marker doesn't exist" to "nothing uses it".
    ⇒ If the marker is ever wanted, that classification IS the work.
  * "these two computations are mathematically equal" does NOT justify bit
    identity: two kernels summing the same terms in a different order differ
    by ~1 ULP. Measured on this repo (audit 2026-07-28 §6.2): 244 of 4051
    Jacobian entries differed by 2.2e-16 between `assemble_newton_jacobian`
    with zeroed sensitivities and `assemble_matrix` -- 5 of 10 suite failures
    on a clean environment were of exactly this kind, carrying no information.

Cross-change bit identity (before/after a refactor, same machine, same
environment) stays available as a DEVELOPMENT tool: `bench/bitcheck.py`.
"""

import numpy as np

REL_TOL = 1.0e-12


def rel_diff(got, ref, nan_equal: bool = True):
    """Worst relative difference between two arrays (NaN-aware).

    Returns (worst_rel, index_of_worst). Scale is max(|ref|, |got|) per entry,
    so an entry that is zero in one array and tiny in the other is judged
    against the tiny value, not against zero.
    """
    g = np.asarray(got, dtype=np.float64).ravel()
    r = np.asarray(ref, dtype=np.float64).ravel()
    if g.shape != r.shape:
        raise AssertionError(f"shape mismatch: {g.shape} vs {r.shape}")
    both_nan = np.isnan(g) & np.isnan(r)
    scale = np.maximum(np.abs(r), np.abs(g))
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.abs(g - r) / np.where(scale > 0.0, scale, 1.0)
    d[both_nan] = 0.0
    if d.size == 0:
        return 0.0, -1
    i = int(np.nanargmax(d))
    return float(d[i]), i


def assert_rel_close(got, ref, rtol: float = REL_TOL, nan_equal: bool = True,
                     msg: str = ""):
    """Assert max relative difference <= rtol (default 1e-12), NaN == NaN.

    Drop-in replacement for `assert np.array_equal(a, b)` on float arrays.
    """
    worst, i = rel_diff(got, ref, nan_equal=nan_equal)
    if not (worst <= rtol):
        g = np.asarray(got, dtype=np.float64).ravel()
        r = np.asarray(ref, dtype=np.float64).ravel()
        n_bad = int(np.count_nonzero(
            np.abs(g - r) > rtol * np.maximum(np.abs(g), np.abs(r))))
        raise AssertionError(
            f"{msg + ': ' if msg else ''}max relative difference "
            f"{worst:.3e} > rtol {rtol:.1e} "
            f"({n_bad}/{g.size} entries); worst at flat index {i}: "
            f"got {g[i]!r} vs ref {r[i]!r}"
        )
    return worst
