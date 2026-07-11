"""
Level-set representation of a straight wake surface (Track B, B1).

The wake is the ruled half-surface swept from the trailing-edge polyline
along the (freestream) wake direction -- Lopez dissertation section 3.5 /
design_track_b.md section 3. It is represented implicitly: for a query
point x the level set returns

  s(x): signed offset from the wake surface plane through the spanwise
        closest TE segment ("+" = upper side),
  d(x): downstream coordinate along the wake direction, measured from the
        TE (d < 0 means upstream of the TE -- NOT on the wake half-surface),
  q(x): spanwise arclength along the TE polyline, UNCLAMPED -- q < 0 or
        q > span_length means the point is off the ends of the sheet.

Cut classification needs all three, because the sheet is a BOUNDED strip:

  - d > 0 excludes the region ahead of the leading edge, where the domain
    is connected across the wake plane's own extension (the naive sign
    test alone would "cut" it);
  - 0 <= q <= span_length excludes the region OUTBOARD OF THE WING TIP,
    where the sheet has ended: the tip edge is a free edge, the jump must
    vanish there (Gamma(tip) = 0 discretely -- the same semantics the
    conforming preprocessor gets from its free-edge rule,
    mesh/wake_cut.py). Without this clip a swept-wing level set would
    happily cut the whole wake-plane extension beyond the tip -- exactly
    the far-field branch-ray artifact P5 had to fix on the conforming
    path.

On the quasi-2D meshes q spans the extruded layer and no node falls off
the ends, so the clip is inert there.

3D-readiness (design_track_b.md D9): the TE is a POLYLINE, not a point --
a swept wing's ruled wake has a per-segment frame. The quasi-2D meshes are
the degenerate single-segment case. Sign convention: the caller orders the
TE polyline so that t_span x d_hat points to the intended "+" (upper) side;
for the quasi-2D NACA family (span +z, chord +x, lift +y) passing the TE
segment in +z order gives "+" = +y at alpha = 0.

Multi-wake: one WakeLevelSet instance per wake (independent by design).
Angle-of-attack sweeps: update_direction() re-aims the ruled surface; the
mesh never changes (the Track B workflow payoff).
"""

from typing import Optional, Tuple

import numpy as np


class WakeLevelSet:
    """Signed-offset/downstream evaluation of a straight ruled wake surface.

    Args:
        te_points: (m, 3) TE polyline (m >= 2), or (3,) / (1, 3) single TE
            point -- then a synthetic span segment along ``span_hint`` is
            used (still exercised by the same per-segment code path).
        direction: (3,) wake direction (typically the freestream); need not
            be normalized. Must not be parallel to any TE segment.
        span_hint: (3,) synthetic span axis for the single-point form.
        extent: optional downstream length of the sheet (None = to infinity;
            the far field truncates it physically). Recorded for consumers;
            evaluation itself does not clip.
    """

    def __init__(
        self,
        te_points: np.ndarray,
        direction: np.ndarray,
        span_hint: Tuple[float, float, float] = (0.0, 0.0, 1.0),
        extent: Optional[float] = None,
    ):
        pts = np.atleast_2d(np.asarray(te_points, dtype=np.float64))
        if pts.shape[1] != 3:
            raise ValueError("te_points must be (m, 3)")
        if len(pts) == 1:
            span = np.asarray(span_hint, dtype=np.float64)
            span = span / np.linalg.norm(span)
            # Long synthetic segment so clamping never activates in practice.
            pts = np.vstack([pts[0] - 1e3 * span, pts[0] + 1e3 * span])
        self.te_points = pts
        self.extent = extent
        self._seg_a = pts[:-1]                       # (nseg, 3)
        self._seg_v = pts[1:] - pts[:-1]             # (nseg, 3)
        seg_len2 = np.einsum("ij,ij->i", self._seg_v, self._seg_v)
        if np.any(seg_len2 <= 0.0):
            raise ValueError("degenerate TE segment (repeated points)")
        self._seg_len2 = seg_len2
        self._seg_len = np.sqrt(seg_len2)
        # arclength at each segment's start; total = span_length
        self._seg_q0 = np.concatenate([[0.0], np.cumsum(self._seg_len)[:-1]])
        self.span_length = float(self._seg_len.sum())
        self._d_hat = np.empty(3)
        self._seg_n = np.empty_like(self._seg_v)
        self.update_direction(direction)

    @property
    def direction(self) -> np.ndarray:
        return self._d_hat.copy()

    def update_direction(self, new_direction: np.ndarray) -> None:
        """Re-aim the ruled surface (alpha change); mesh untouched."""
        d = np.asarray(new_direction, dtype=np.float64)
        nrm = np.linalg.norm(d)
        if nrm == 0.0:
            raise ValueError("direction must be nonzero")
        self._d_hat = d / nrm
        # Per-segment upper normal: unit(t_span x d_hat).
        n = np.cross(self._seg_v, self._d_hat[None, :])
        n_norm = np.linalg.norm(n, axis=1)
        if np.any(n_norm < 1e-12 * np.sqrt(self._seg_len2)):
            raise ValueError("wake direction (near-)parallel to a TE segment")
        self._seg_n = n / n_norm[:, None]
        # Gram data of the OBLIQUE in-plane basis (v, d_hat). On a SWEPT
        # TE the span direction is not perpendicular to the wake direction,
        # so an orthogonal projection would leak the downstream distance
        # into the spanwise coordinate (measured: it pushes far-downstream
        # points past the tip and clips ~60% of the real M6 cut set).
        self._a12 = self._seg_v @ self._d_hat            # v . d_hat
        self._det = self._seg_len2 - self._a12**2        # > 0 (guarded above)

    def evaluate(
        self, points: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (s, d, q) at points (n, 3).

        s: signed offset from the wake surface, "+" = upper side.
        d: downstream coordinate from the (projected) TE along direction.
        q: spanwise arclength along the TE polyline, unclamped -- outside
           [0, span_length] the point is off the ends of the sheet (past
           the wing tip, or inboard of the root).

        Each point is decomposed in the panel's OBLIQUE ruled-surface frame
        (v_span, d_hat, n_hat): x - a = u*v + d*d_hat + s*n_hat, solved as
        a 2x2 system in the (v, d_hat) plane (n_hat is orthogonal to both).
        An orthogonal projection onto v would be wrong on a SWEPT TE, where
        v is not perpendicular to the wake direction: the downstream
        distance would leak into the spanwise coordinate. Panel selection
        uses |s| plus the out-of-panel spanwise excess, so a multi-segment
        (curved/kinked) TE works per panel while u stays unclamped -- which
        is what lets q report the beyond-the-tip excess instead of
        saturating at the tip.
        """
        x = np.atleast_2d(np.asarray(points, dtype=np.float64))
        n_pts = len(x)
        rel = x[:, None, :] - self._seg_a[None, :, :]          # (n, nseg, 3)
        b1 = np.einsum("pns,ns->pn", rel, self._seg_v)         # rel . v
        b2 = rel @ self._d_hat                                 # rel . d_hat
        u = (b1 - self._a12 * b2) / self._det                  # span fraction
        d_all = (self._seg_len2 * b2 - self._a12 * b1) / self._det
        s_all = np.einsum("pns,ns->pn", rel, self._seg_n)      # (n, nseg)

        excess = np.maximum(0.0, np.maximum(-u, u - 1.0)) * self._seg_len
        dist2 = s_all**2 + excess**2
        best = np.argmin(dist2, axis=1)                        # (n,)
        idx = np.arange(n_pts)
        s = s_all[idx, best]
        d = d_all[idx, best]
        q = self._seg_q0[best] + u[idx, best] * self._seg_len[best]
        return s, d, q
