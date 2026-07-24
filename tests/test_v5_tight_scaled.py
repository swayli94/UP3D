"""Track V V5.1b -- the scaled + damped augmented Newton path (design:
cases/analysis/v5_1b_scaled_newton/PRE_REGISTRATION.md; module under test:
pyfp3d/viscous/tight_driver.py -- the newton_tight
scaling/lm_damping/floor_stop path plus the equilibrate_rc /
scaled_damped_step / mu-schedule / FloorStop helpers).

The GV5.1b helpers are solver-internal: the assembled F/J are
bit-identical to GV5.1 (the committed FD verdicts stand); these tests
gate only the transform + globalization layer:

  - equilibrate_rc: the one-pass row/column 2-norm equilibration is an
    exact factorization (diag(rn) @ Jsc @ diag(cn) == J to machine
    precision), unit column norms, zero-safe scales on a zero row/col;
  - scaled_damped_step at mu = 0: the unscaled path recovers the
    undamped splu step exactly, and scaling="rowcol" solves the SAME
    mathematical system (round-trip dx = C dy matches the direct solve);
  - the mu schedule helpers (mu_on_accept /3 bounded below,
    mu_on_reject x10 bounded above) and the FloorStop counter
    (FLOOR_CONSEC consecutive relative decreases below FLOOR_REL_TOL,
    reset by any larger decrease; a negative decrease counts as below);
  - the end-to-end mu schedule on a mock pack (monkeypatched
    augmented_residual/augmented_jacobian/_closure_packet): an
    always-contracting model accepts lam = 1 every iteration and mu
    decays /3 per accepted step; an always-expanding model exhausts the
    30-halving line search at every mu, climbs x10 to MU_MAX, takes the
    legacy least-bad fallback, and never decreases mu afterwards;
  - the floor-reached stop: a constant-F model (merit exactly flat)
    terminates "floor_reached" at iteration FLOOR_CONSEC with
    floor_stop=True and runs to the cap without it;
  - the k=1 smoke (JIT lane only): the new path runs on the real
    pre-registered k=1 pack (finite iterates, a recorded mu per
    iteration), and the DEFAULT legacy path reproduces the committed
    GV5.1 k1seed history (cases/analysis/v5_tight_coupling/results/
    gv5_1_newton_history_coarse_k1seed.csv) -- the bit-for-bit guard on
    the untouched legacy branch.
"""

import csv
import os
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as sla

from pyfp3d.viscous import closures as C
from pyfp3d.viscous import tight_driver as td
from tests.v5_state import (
    M_CRIT,
    M_CAP,
    RHO_FLOOR,
    UPWIND_C,
    build_k1_state,
    build_naca_case,
)

NOJIT = os.environ.get("PYFP3D_NOJIT", "0") == "1"
K1_JIT_ONLY = pytest.mark.skipif(
    NOJIT, reason="the k=1 fixture (FP + IBL solves) is JIT-lane only"
)

REPO_ROOT = Path(__file__).parent.parent
K1SEED_CSV = (
    REPO_ROOT
    / "cases"
    / "analysis"
    / "v5_tight_coupling"
    / "results"
    / "gv5_1_newton_history_coarse_k1seed.csv"
)


@pytest.fixture(scope="module")
def naca_case():
    return build_naca_case()


@pytest.fixture(scope="module")
def k1_state(naca_case):
    return build_k1_state(naca_case)


@pytest.fixture(scope="module")
def pack(k1_state):
    """The pre-registered k=1 pack (the GV5.1 smoke fixture)."""
    return td.build_tight_pack(k1_state, UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)


# ---------------------------------------------------------------------------
# equilibrate_rc
# ---------------------------------------------------------------------------


def test_equilibrate_rc_identity():
    """diag(rn) @ Jsc @ diag(cn) == J to machine precision; unit column
    norms; a zero row/column gets scale 1 (the pre-registered recipe)."""
    rng = np.random.default_rng(7)
    n = 60
    J = sp.random(n, n, density=0.15, format="lil", random_state=rng)
    J[5, :] = 0.0  # a zero row
    J[:, 9] = 0.0  # a zero column
    J = J.tocsr()

    rn, cn, Jsc = td.equilibrate_rc(J)

    Jr = (sp.diags(rn) @ Jsc @ sp.diags(cn)).tocsr()
    scale = float(np.max(np.abs(J.toarray())))
    rel = float(np.max(np.abs((Jr - J).toarray()))) / scale
    assert rel <= 1.0e-14

    assert rn[5] == 1.0
    assert cn[9] == 1.0

    cols = np.sqrt(np.asarray(Jsc.multiply(Jsc).sum(axis=0)).ravel())
    assert cols[9] == 0.0
    mask = np.ones(n, dtype=bool)
    mask[9] = False
    assert np.all(np.abs(cols[mask] - 1.0) <= 1.0e-12)


# ---------------------------------------------------------------------------
# scaled_damped_step
# ---------------------------------------------------------------------------


def _well_conditioned_j(seed=11, n=60):
    rng = np.random.default_rng(seed)
    # strictly diagonally dominant -> invertible, well conditioned
    return (
        sp.random(n, n, density=0.1, format="csr", random_state=rng)
        + 8.0 * sp.eye(n, format="csr")
    ).tocsr()


def test_mu_zero_step_equals_undamped():
    """mu = 0 on the unscaled path recovers the undamped splu step."""
    J = _well_conditioned_j()
    F = np.random.default_rng(13).standard_normal(J.shape[0])
    d = td.scaled_damped_step(J, F, 0.0)
    d_ref = -sla.splu(J.tocsc()).solve(F)
    assert np.max(np.abs(d - d_ref)) <= 1.0e-13 * np.max(np.abs(d_ref))


def test_rowcol_roundtrip():
    """scaling="rowcol" at mu = 0 solves the SAME system: dx = C dy
    matches the direct splu step (the live band-(a) identity)."""
    J = _well_conditioned_j(seed=17)
    F = np.random.default_rng(19).standard_normal(J.shape[0])
    d_rc = td.scaled_damped_step(J, F, 0.0, scaling="rowcol")
    d_ref = -sla.splu(J.tocsc()).solve(F)
    rel = np.max(np.abs(d_rc - d_ref)) / np.max(np.abs(d_ref))
    assert rel <= 1.0e-9


# ---------------------------------------------------------------------------
# the mu schedule + FloorStop helpers
# ---------------------------------------------------------------------------


def test_mu_schedule_helpers():
    assert td.mu_on_accept(1.0e-6) == pytest.approx(1.0e-6 / 3.0)
    assert td.mu_on_accept(1.0e-13) == td.MU_MIN
    assert td.mu_on_reject(1.0e-6) == pytest.approx(1.0e-5)
    assert td.mu_on_reject(50.0) == td.MU_MAX

    # a large decrease resets the consecutive-below-tol counter
    fl = td.FloorStop()
    assert fl.update(1.0e-5) is False
    assert fl.update(0.5) is False
    assert fl.update(1.0e-5) is False
    assert fl.update(1.0e-5) is False  # only 2 consecutive since reset

    fl2 = td.FloorStop()
    assert fl2.update(1.0e-5) is False
    assert fl2.update(1.0e-5) is False
    assert fl2.update(1.0e-5) is True

    # a negative relative decrease (merit rose on the least-bad
    # fallback) still counts as below FLOOR_REL_TOL
    fl3 = td.FloorStop()
    assert fl3.update(-1.0) is False
    assert fl3.update(-0.25) is False
    assert fl3.update(0.0) is True


# ---------------------------------------------------------------------------
# the end-to-end mu schedule on a mock pack
# ---------------------------------------------------------------------------


class _MockPack:
    """The minimal TightPack surface newton_tight touches: x_base(),
    split_F() (3 blocks), split_x() ((phi, gamma, U) with U n x 6)."""

    def __init__(self, x0):
        self._x0 = np.asarray(x0, dtype=np.float64)

    def x_base(self):
        return self._x0.copy()

    def split_F(self, F):
        return F[:4], F[4:8], F[8:]

    def split_x(self, x):
        return x[:4], x[4:6], x[6:].reshape(-1, 6)


@pytest.fixture
def mock_system(monkeypatch):
    """Monkeypatch the pack-facing operators; ds_at gets a zero
    closure packet (ds_change == 0, unused by these gates)."""

    def fake_closure(pack, U):
        return np.zeros((U.shape[0], C.OUT_DS1 + 1)), None

    monkeypatch.setattr(td, "_closure_packet", fake_closure)

    def install(F_fun, J_mat):
        monkeypatch.setattr(
            td, "augmented_residual", lambda pack, x: F_fun(np.asarray(x))
        )
        monkeypatch.setattr(
            td, "augmented_jacobian", lambda pack, x: J_mat
        )

    return install


def test_mu_schedule_accept_mock(mock_system):
    """An always-contracting model (true F = A(x - x*), A spd, J = A/2
    scaled to the identity by equilibrate_rc): every step accepts at
    lam = 1 and mu decays /3 per accepted step."""
    n = 12
    a = np.array(
        [1.0, 2.0, 0.5, 4.0, 1.5, 2.5, 0.75, 3.0, 1.25, 1.75, 2.25, 0.9]
    )
    x_star = np.linspace(-0.5, 0.7, n)
    x0 = x_star + np.linspace(0.3, -0.2, n)
    mock_system(lambda x: a * (x - x_star), sp.csr_matrix(np.diag(0.5 * a)))

    res = td.newton_tight(
        _MockPack(x0),
        max_iter=6,
        scaling="rowcol",
        lm_damping=True,
        floor_stop=False,
    )

    assert res["termination"] == "cap"
    assert not res["converged"]
    assert res["n_iter"] == 6
    hist = res["history"]
    assert len(hist) == 7
    for k, h in enumerate(hist[1:]):
        assert h["mu"] == pytest.approx(td.MU0 / 3.0**k, rel=1.0e-12)
        assert h["lam"] == pytest.approx(1.0)
        assert h["accepted"]
        assert h["mu_retries"] == 0
    merits = [h["merit"] for h in hist]
    assert all(m2 < m1 for m1, m2 in zip(merits, merits[1:]))


def test_mu_schedule_reject_mock(mock_system):
    """An always-expanding model (true F = B(x - x*), B = diag(-1, 1,
    ...), J = I): F^T d > 0, the merit rises along the whole ray, all 30
    halvings reject at every mu; mu climbs x10 to MU_MAX, the legacy
    least-bad fallback is taken, and mu never decreases afterwards."""
    n = 12
    b = np.ones(n)
    b[0] = -1.0
    x_star = np.linspace(-0.3, 0.4, n)
    x0 = x_star - np.eye(n)[0]
    mock_system(lambda x: b * (x - x_star), sp.eye(n, format="csr"))

    res = td.newton_tight(
        _MockPack(x0),
        max_iter=2,
        scaling=None,
        lm_damping=True,
        floor_stop=False,
    )

    assert res["termination"] == "cap"
    assert res["n_iter"] == 2
    hist = res["history"]
    # iter 1: solves at mu = 1e-6 .. 1e2 (8 reject-retries over the 8
    # decades), least-bad taken at the cap
    assert hist[1]["mu"] == pytest.approx(td.MU_MAX)
    assert hist[1]["mu_retries"] == 8
    assert not hist[1]["accepted"]
    assert hist[1]["lam"] == pytest.approx(2.0**-29)
    # iter 2: mu stays at the cap (no accept-decrease), first solve
    # already at MU_MAX -> no retries
    assert hist[2]["mu"] == pytest.approx(td.MU_MAX)
    assert hist[2]["mu_retries"] == 0
    assert not hist[2]["accepted"]
    assert hist[2]["lam"] == pytest.approx(2.0**-29)
    # the least-bad fallback RAISES the merit slightly
    assert hist[1]["merit"] > hist[0]["merit"]
    assert hist[2]["merit"] > hist[1]["merit"]


def test_floor_stop_mock(mock_system):
    """A constant-F model: the merit is exactly flat, equality accepts,
    rel_dec == 0 every step -> floor_reached at iteration FLOOR_CONSEC
    with floor_stop=True; the control without it runs to the cap."""
    n = 12
    c = np.linspace(1.0, 2.0, n)
    mock_system(lambda x: c.copy(), sp.eye(n, format="csr"))

    res = td.newton_tight(
        _MockPack(np.zeros(n)),
        max_iter=10,
        scaling="rowcol",
        lm_damping=True,
        floor_stop=True,
    )
    assert res["termination"] == "floor_reached"
    assert res["n_iter"] == td.FLOOR_CONSEC
    assert len(res["history"]) == td.FLOOR_CONSEC + 1

    res2 = td.newton_tight(
        _MockPack(np.zeros(n)),
        max_iter=5,
        scaling="rowcol",
        lm_damping=True,
        floor_stop=False,
    )
    assert res2["termination"] == "cap"
    assert res2["n_iter"] == 5


# ---------------------------------------------------------------------------
# the k=1 smoke (JIT lane only)
# ---------------------------------------------------------------------------


@K1_JIT_ONLY
def test_k1_smoke(pack):
    """The new path runs on the real k=1 pack; the DEFAULT legacy path
    reproduces the committed GV5.1 k1seed history (the bit-for-bit
    guard on the untouched legacy branch)."""
    res_new = td.newton_tight(
        pack,
        max_iter=2,
        scaling="rowcol",
        lm_damping=True,
        floor_stop=True,
    )
    assert np.all(np.isfinite(res_new["x"]))
    assert res_new["termination"] in ("converged", "cap", "floor_reached")
    assert all("mu" in h for h in res_new["history"])

    res_leg = td.newton_tight(pack, max_iter=2)
    with open(K1SEED_CSV, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for it in (1, 2):
        h = res_leg["history"][it]
        row = rows[it]
        for j, key in enumerate(("f_phi_max", "f_gamma_max", "f_bl_max")):
            ref = float(row[key])
            got = float(h["block_max"][j])
            if ref > 1.0e-14:
                assert abs(got - ref) <= 2.0e-6 * ref, (
                    f"iter {it} {key}: got {got:.9e}, ref {ref:.9e}"
                )
            else:
                # f_gamma sits at machine zero in the committed history
                assert abs(got) <= 1.0e-14
