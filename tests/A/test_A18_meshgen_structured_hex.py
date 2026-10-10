"""Locks for `pyfp3d/meshgen/structured.py` -- the hybrid layer-block generator built in task 2.

★★ Why this file exists. Route (A)'s one surviving deliverable is EXPERIMENTAL CONTROLLABILITY:
`r_far` / `h_wall_normal` / `n_theta` / `growth` are bit-identical single-variable knobs, which
phase two's four gmsh knobs were all measured NOT to be. That property is what every reading built
on this generator rests on -- and until now it was asserted only by `phases/p3/bench/run_hex_g0_single_var.py`,
a bench script that runs on no cadence. That is precisely the gap round 2b closed for the conforming
wing-body transonic capability: the boundary document stated a capability and nothing would go red
if it broke, right before the work most likely to break it.

Everything here is pure mesh arithmetic (no solves), so it belongs in the always-on tier.

Bit-identity, not "small change", is the criterion wherever the generator claims it -- the G0
verdict's own standard.
"""

import numpy as np
import pytest

from pyfp3d.meshgen.planar import naca0012_coordinates
from pyfp3d.meshgen.structured import (airfoil_layer_block_2d, airfoil_surface_distribution,
                                       layer_block_quality, n_layers_of, wall_anchored_distances,
                                       wall_anchored_radii)

#: the production configuration of `cases/meshes/naca0012_hex_2.5d/generate_naca0012_hex.py`
PROD = dict(thickness=0.08, h_wall_normal=0.004, growth=1.15, te_blend=0.05)


@pytest.fixture(scope="module")
def ref():
    return naca0012_coordinates(n_half=401)


class TestWallAnchoredGrading:
    """The grading is anchored AT THE WALL with a derived layer count -- that is what makes
    `r_far` and `h_wall_normal` single-variable rather than merely well-behaved."""

    def test_first_step_is_exactly_h_wall_normal(self):
        r = wall_anchored_radii(0.5, 20.0, 0.004, 1.15)
        assert r[0] == 0.5
        assert r[1] - r[0] == pytest.approx(0.004, rel=1e-15)

    def test_growing_r_far_only_APPENDS(self):
        """★ bit-identical prefix: extending the domain must not re-space the existing layers.
        A generator that rescaled to hit r_far exactly would fail this and would silently move
        every committed mesh when the domain size changed."""
        near = wall_anchored_radii(0.5, 20.0, 0.004, 1.15)
        far = wall_anchored_radii(0.5, 60.0, 0.004, 1.15)
        assert len(far) > len(near)
        assert np.array_equal(far[: len(near)], near)

    def test_h_wall_normal_halved_keeps_the_wall_radius_exact(self):
        a = wall_anchored_radii(0.5, 20.0, 0.004, 1.15)
        b = wall_anchored_radii(0.5, 20.0, 0.002, 1.15)
        assert a[0] == b[0] == 0.5
        assert b[1] - b[0] == pytest.approx(0.002, rel=1e-15)

    def test_monotone_and_reaches_the_far_field(self):
        r = wall_anchored_radii(0.5, 20.0, 0.004, 1.15)
        assert np.all(np.diff(r) > 0.0)
        assert r[-1] >= 20.0

    def test_distances_are_NOT_a_shifted_copy_of_radii(self):
        """★ Deliberate duplication, locked so nobody "simplifies" it away: accumulating from
        r_wall differs from accumulating from 0 and shifting, in the last bits (measured ~8.9e-16
        absolute). Collapsing the two would silently move every committed mesh, so the two
        functions stay separate and this test records WHY."""
        r = wall_anchored_radii(0.5, 20.0, 0.004, 1.15)
        d = wall_anchored_distances(19.5, 0.004, 1.15)
        assert len(r) == len(d)
        assert not np.array_equal(r, 0.5 + d)
        assert np.allclose(r, 0.5 + d, rtol=1e-14, atol=1e-14)


class TestSurfaceDistribution:
    def test_station_count_is_exact(self, ref):
        for n in (80, 160, 320):
            assert len(airfoil_surface_distribution(ref, n)) == n

    def test_le_cluster_puts_the_fine_spacing_at_the_LE_not_the_TE(self, ref):
        """The TE de-clustering measured in task 2: a plain cosine is dense at BOTH ends, which
        made the TE finer than the LE at equal station count."""
        c = airfoil_surface_distribution(ref, 160, le_cluster=0.8)
        d = np.linalg.norm(np.diff(c, axis=0), axis=1)
        xm = 0.5 * (c[:-1, 0] + c[1:, 0])
        le = np.median(d[xm < 0.05])
        te = np.median(d[xm > 0.95])
        assert le < te, f"LE spacing {le:.5f} should be finer than TE {te:.5f}"

    def test_local_window_default_is_bit_identical(self, ref):
        """Two ways of asking for "no local refinement" must both be bit-identical to the legacy
        call -- an inert-by-default knob is the project's standing requirement for a new option."""
        base = airfoil_surface_distribution(ref, 160, le_cluster=0.8)
        assert np.array_equal(base, airfoil_surface_distribution(
            ref, 160, le_cluster=0.8, local_window=None, local_factor=3.0))
        assert np.array_equal(base, airfoil_surface_distribution(
            ref, 160, le_cluster=0.8, local_window=(0.45, 0.75), local_factor=1.0))

    def test_local_window_refines_inside_and_the_COST_outside_is_locked(self, ref):
        """★ This knob's guard is NOT bit-identity, unlike the four G0 knobs, and the difference is
        recorded rather than glossed: holding the station count fixed means points are MOVED IN, so
        the outside gets coarser. The count is held fixed on purpose -- adding stations would change
        the global DOF count too, and the shock-band leg could then not be separated from the
        global-refinement leg, which is the exact confound phase two could never escape.

        Both halves are locked, so a later "improvement" that quietly changes the trade goes red.
        """
        w, n = (0.45, 0.75), 160
        base = airfoil_surface_distribution(ref, n, le_cluster=0.8)
        loc = airfoil_surface_distribution(ref, n, le_cluster=0.8, local_window=w, local_factor=3.0)
        assert len(loc) == len(base) == n

        def med(c):
            d = np.linalg.norm(np.diff(c, axis=0), axis=1)
            xm = 0.5 * (c[:-1, 0] + c[1:, 0])
            m = (xm >= w[0]) & (xm <= w[1])
            return np.median(d[m]), np.median(d[~m]), int(m.sum())

        din_b, dout_b, nin_b = med(base)
        din_l, dout_l, nin_l = med(loc)
        assert din_l / din_b < 0.6, "the window must actually refine (measured 0.512)"
        assert nin_l > nin_b, f"window point count {nin_b} -> {nin_l} must grow"
        #: the measured cost, locked as a BOUND: 1.124 at close-out
        assert 1.0 < dout_l / dout_b < 1.3, (
            f"outside-window coarsening {dout_l / dout_b:.3f} left the recorded band -- the "
            f"confound bound in the task-3 verdict was measured at 1.124 and any change to it "
            f"changes what that verdict may claim")


class TestLayerBlockQuality:
    """★ Both tests in this class were RED on their first run, and both because I indexed the point
    array by `len(surface)` instead of reading `info["n_stations"]`. The TE node is SPLIT for
    marching, so `n_stations == len(surface) + 1`, and `pts[len(surface)]` is the split TE COPY at
    the same coordinates as `pts[0]` -- a zero step, i.e. NaN. That is the same class of mistake as
    the LE probe earlier in this round: read the layout the generator reports, never re-derive it.
    """

    @staticmethod
    def _block(ref, n_stations=160):
        surf = airfoil_surface_distribution(ref, n_stations, le_cluster=0.8)
        #: the return order is (pts, tris, outer_idx, wall_idx, info) -- read, do not recall
        pts, tris, outer_idx, wall_idx, info = airfoil_layer_block_2d(surf, **PROD)
        return pts, tris, outer_idx, wall_idx, info

    def test_production_block_quality_drift_lock(self, ref):
        """A DRIFT lock, not a correctness claim: these are the numbers measured after the user's
        architecture correction (finite-height block + TE direction pinned to +-y + TE
        de-clustering), which took the aspect ratio from 1305 to 5.9. The v1 O-grid passed an
        orientation-only guard AT AR 1305, so "no inverted cells" is NOT "usable" -- hence the
        aspect bound is asserted, not merely recorded."""
        pts, tris, outer_idx, wall_idx, info = self._block(ref)
        q = layer_block_quality(pts, tris, wall_idx, info["n_stations"], n_layers_of(info))
        assert q["ray_crossings"] == 0
        assert q["aspect_max"] < 8.0, f"aspect_max {q['aspect_max']:.3f} (measured 5.899)"
        assert q["aspect_median"] < 3.0, f"aspect_median {q['aspect_median']:.3f} (measured 2.012)"
        #: v1 had 1.12e-07 here; measured 5.08e-06
        assert q["min_area"] > 1e-6, f"min_area {q['min_area']:.3e}"
        assert q["orthogonality_dev_median_deg"] < 1.0, (
            f"orthogonality median {q['orthogonality_dev_median_deg']:.3f} deg (measured 0.251)")

    def test_te_marching_direction_is_EXACTLY_plus_minus_y(self, ref):
        """At a sharp TE the averaged wall normal points ~+x, so marching along it sends the TE
        column straight down the wake line and the neighbouring rays cross just behind it. The fix
        blends to +-y over the last `te_blend` of chord and reaches PURE +-y at the TE, with the
        side taken from node ORDER because the TE has y == 0 exactly.

        Measured exact, so locked exact rather than "y-dominated": [0, +1] on the first station and
        [0, -1] on the split copy. A blend that stopped short would loosen these to non-integers.
        """
        pts, tris, outer_idx, wall_idx, info = self._block(ref)
        ns = info["n_stations"]
        assert info["te_split"] is True
        assert len(pts) == ns * n_layers_of(info), "the raw block is NOT welded; the weld is later"

        up = pts[ns] - pts[0]
        up = up / np.linalg.norm(up)
        lo = pts[2 * ns - 1] - pts[ns - 1]
        lo = lo / np.linalg.norm(lo)
        assert up[0] == 0.0 and up[1] == pytest.approx(1.0, abs=1e-15), f"upper TE ray {up}"
        assert lo[0] == 0.0 and lo[1] == pytest.approx(-1.0, abs=1e-15), f"lower TE ray {lo}"
        assert np.sign(up[1]) != np.sign(lo[1]), "the two split TE copies must march APART"
