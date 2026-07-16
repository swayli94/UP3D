"""
Track A / A2 -- shared metrics & estimators for the TE/Kutta fidelity study.

Case-local on purpose: A2 measures the committed machinery, it does not touch
`pyfp3d/` (the correctness anchor is GA2.1's reproduction of the committed B7
roughness numbers, not a unit-test suite).

Contents:
  * roughness_d2       -- the B7 jitter metric, VERBATIM (b7_onera_m6/run_demo.py
                          "roughness"): RMS second difference / curve range.
                          Callers order the curve by z first; B7's inputs were
                          already z-ordered.
  * local_fit / jitter_localfit -- windowed local-quadratic fit (tricube
                          weights) and the RMS-residual jitter metric it
                          induces. Robust to the non-uniform station spacing
                          the swept unstructured TE produces; also the
                          "diagnostic smooth Gamma" generator for the gated
                          intervention leg (T3/E fixed-Gamma methodology --
                          NOT the refuted smoothing-as-fix route).
  * wall_value_section -- z = const plane cut of the wall triangulation
                          carrying a NODAL field (phi): edge-crossing points
                          with linear interpolation. Side = triangle outward
                          normal n_y (the surface_ls convention -- geometric
                          y>0 misclassifies exactly at the TE, where both
                          surfaces collapse onto y ~ 0).
  * delta_phi_profile  -- Delta phi(z; x/c) = phi_upper - phi_lower read on
                          the WALL at x/c < 1. Not pinned by the wake
                          constraint (the wake jump equals the prescribed
                          per-station Gamma identically), so it discriminates
                          "jitter confined to the TE-adjacent enforcement
                          layer" from "circulation field genuinely rough".
  * cp_section_from_tri -- z = const cut carrying a PER-TRIANGLE field (Cp):
                          the cheap sweep variant of
                          post/section_cut._wall_section_points (one recovery,
                          many stations).
  * spike_metric / te_gap_from_curve -- the S2 (TE Cp) decomposition metrics.
  * probe_census       -- per-station Kutta-probe geometry features
                          (distances, off-plane offsets, probe sharing), the
                          H1 correlation table inputs.
  * pearson            -- correlation with NaN guard (supporting evidence
                          only; the P5 meta-lesson stands: correlation is not
                          causation, the intervention leg decides H1).
"""

import numpy as np


# --------------------------------------------------------------------------
# S1 jitter metrics
# --------------------------------------------------------------------------

def roughness_d2(g):
    """B7's committed jitter metric, verbatim: RMS 2nd difference / range.

    The caller passes the curve ordered by z (B7 fed z-ordered arrays)."""
    g = np.asarray(g, dtype=np.float64)
    d2 = g[:-2] - 2.0 * g[1:-1] + g[2:]
    return float(np.sqrt(np.mean(d2 ** 2)) / (g.max() - g.min()))


def local_fit(z, g, window_frac=0.08):
    """Windowed local-quadratic fit of g(z) (tricube weights).

    The window is a fraction of the span extent, so coarse/medium station
    counts see the same PHYSICAL smoothing length. Wide enough to leave the
    genuine tip roll-off in the fit and only station-scale content in the
    residual."""
    z = np.asarray(z, dtype=np.float64)
    g = np.asarray(g, dtype=np.float64)
    h = window_frac * (z.max() - z.min())
    fit = np.empty_like(g)
    for i, zi in enumerate(z):
        u = np.abs(z - zi) / h
        w = np.where(u < 1.0, (1.0 - u ** 3) ** 3, 0.0)
        m = w > 0
        if m.sum() < 4:                       # not enough support: fall back
            fit[i] = g[i]
            continue
        A = np.vander(z[m] - zi, 3)           # [dz^2, dz, 1]
        W = w[m]
        coef, *_ = np.linalg.lstsq(A * W[:, None], g[m] * W, rcond=None)
        fit[i] = coef[2]
    return fit


def jitter_localfit(z, g, window_frac=0.08):
    """RMS residual around the local-quadratic fit, normalised by the curve
    range (same normaliser as roughness_d2, so the two are comparable)."""
    g = np.asarray(g, dtype=np.float64)
    resid = g - local_fit(z, g, window_frac)
    return float(np.sqrt(np.mean(resid ** 2)) / (g.max() - g.min()))


# --------------------------------------------------------------------------
# Wall plane cuts
# --------------------------------------------------------------------------

def _tri_sides_ny(nodes, elements, wall_faces):
    """Per-wall-triangle upper mask from the outward normal's lift component
    (n_y > 0 = upper) -- the post/surface_ls side convention."""
    from pyfp3d.post.surface import wall_outward_normals
    n_out = wall_outward_normals(nodes, elements, wall_faces)
    return n_out[:, 1] > 0.0


def wall_value_section(nodes, wall_faces, values, z, upper_mask):
    """z = const cut of the wall triangulation carrying a NODAL field.

    Each crossed triangle edge contributes one point with linearly
    interpolated coordinates and value; the point inherits the triangle's
    side. Returns (x, val, is_upper) unordered."""
    tri = np.asarray(wall_faces, dtype=np.int64)
    s = nodes[tri, 2] - z                                  # (n_tri, 3)
    xs, vals, side = [], [], []
    crossed = np.where(~((s > 0).all(axis=1) | (s < 0).all(axis=1)))[0]
    edges = ((0, 1), (1, 2), (2, 0))
    for t in crossed:
        for a, b in edges:
            sa, sb = s[t, a], s[t, b]
            if sa == sb or (sa > 0) == (sb > 0):
                continue
            w = sa / (sa - sb)
            na, nb = tri[t, a], tri[t, b]
            xs.append((1 - w) * nodes[na, 0] + w * nodes[nb, 0])
            vals.append((1 - w) * values[na] + w * values[nb])
            side.append(upper_mask[t])
    return np.asarray(xs), np.asarray(vals), np.asarray(side, dtype=bool)


def delta_phi_profile(nodes, wall_faces, phi, z, xc_targets, upper_mask):
    """Delta phi(x/c) = phi_upper - phi_lower on the wall at station z.

    x/c is normalised by the section's own x extent (the swept, tapered
    planform normaliser section_cp_curve uses). `upper_mask` is the
    per-triangle side split (compute ONCE via _tri_sides_ny -- it walks the
    volume adjacency and is far too slow to redo per station). Returns an
    array aligned with xc_targets; NaN where either side lacks support."""
    x, v, s = wall_value_section(nodes, wall_faces, phi, z, upper_mask)
    out = np.full(len(xc_targets), np.nan)
    if len(x) < 6 or s.sum() < 3 or (~s).sum() < 3:
        return out
    x_le, chord = x.min(), x.max() - x.min()
    xc = (x - x_le) / chord
    for i, t in enumerate(np.asarray(xc_targets, dtype=np.float64)):
        vu = _interp_side(xc[s], v[s], t)
        vl = _interp_side(xc[~s], v[~s], t)
        if vu is not None and vl is not None:
            out[i] = vu - vl
    return out


def _interp_side(xc, v, t):
    o = np.argsort(xc)
    xs, vs = xc[o], v[o]
    if t < xs[0] or t > xs[-1]:
        return None
    return float(np.interp(t, xs, vs))


def cp_section_from_tri(nodes, wall_faces, tri_cp, z, upper_mask):
    """z = const cut carrying a PER-TRIANGLE field (the cheap many-station
    sweep: recovery once, cuts everywhere). Mirrors the geometry of
    post/section_cut._wall_section_points (segment midpoint, triangle's own
    constant Cp) but reuses a precomputed tri_cp. Returns (x_mid, cp, is_upper)."""
    tri = np.asarray(wall_faces, dtype=np.int64)
    s = nodes[tri, 2] - z
    xs, cps, side = [], [], []
    crossed = np.where(~((s > 0).all(axis=1) | (s < 0).all(axis=1)))[0]
    edges = ((0, 1), (1, 2), (2, 0))
    for t in crossed:
        pts = []
        for a, b in edges:
            sa, sb = s[t, a], s[t, b]
            if sa == sb or (sa > 0) == (sb > 0):
                continue
            w = sa / (sa - sb)
            pts.append((1 - w) * nodes[tri[t, a]] + w * nodes[tri[t, b]])
        if len(pts) < 2:
            continue
        mid = 0.5 * (pts[0] + pts[1])
        xs.append(mid[0])
        cps.append(tri_cp[t])
        side.append(upper_mask[t])
    return np.asarray(xs), np.asarray(cps), np.asarray(side, dtype=bool)


# --------------------------------------------------------------------------
# S2 (TE Cp) metrics
# --------------------------------------------------------------------------

def spike_metric(x, cp, fit_lo=0.85, fit_hi=0.97):
    """Last-point TE spike: |Cp(last) - quadratic extrapolation of the
    x/c in [fit_lo, fit_hi] trend|. The smooth extrapolated trend is the
    physical inviscid TE recovery; the residual is the artifact."""
    x = np.asarray(x, dtype=np.float64)
    cp = np.asarray(cp, dtype=np.float64)
    o = np.argsort(x)
    x, cp = x[o], cp[o]
    m = (x >= fit_lo) & (x <= fit_hi)
    if m.sum() < 4:
        return {"spike": np.nan, "cp_last": np.nan, "fit_at_1": np.nan}
    c = np.polyfit(x[m], cp[m], 2)
    return {"spike": float(abs(cp[-1] - np.polyval(c, x[-1]))),
            "cp_last": float(cp[-1]),
            "fit_at_1": float(np.polyval(c, 1.0))}


def te_gap_from_curve(sec):
    """|Cp_upper - Cp_lower| at each side's last (max-x/c) point of a
    section_cp_curve()-style dict."""
    iu = int(np.argmax(sec["x_upper"]))
    il = int(np.argmax(sec["x_lower"]))
    return float(abs(sec["cp_upper"][iu] - sec["cp_lower"][il]))


# --------------------------------------------------------------------------
# Kutta-probe census (H1 features)
# --------------------------------------------------------------------------

def probe_census(mc, wc):
    """Per-station probe-geometry features on the cut mesh.

    Returns dict of (n_st,) arrays: d_up/d_lo (probe distance from its TE
    node, station mean), asym (|d_up-d_lo|/(d_up+d_lo)), dz_up/dz_lo
    (off-plane offset |z_probe - z_te|), shared_up/shared_lo (probe node
    also used by another station -- the degeneracy
    INVESTIGATION_kutta_closure.md recorded as a known-robustness item)."""
    nodes = mc.nodes
    n_st = wc.n_stations
    te = wc.te_nodes
    d_up = np.zeros(n_st)
    d_lo = np.zeros(n_st)
    dz_up = np.zeros(n_st)
    dz_lo = np.zeros(n_st)
    counts = np.zeros(n_st)
    for k in range(len(te)):
        j = wc.te_station[k]
        pu, pl, t = wc.kutta_upper[k], wc.kutta_lower[k], te[k]
        d_up[j] += np.linalg.norm(nodes[pu] - nodes[t])
        d_lo[j] += np.linalg.norm(nodes[pl] - nodes[t])
        dz_up[j] += abs(nodes[pu, 2] - nodes[t, 2])
        dz_lo[j] += abs(nodes[pl, 2] - nodes[t, 2])
        counts[j] += 1
    counts[counts == 0] = 1
    d_up, d_lo = d_up / counts, d_lo / counts
    dz_up, dz_lo = dz_up / counts, dz_lo / counts

    def _shared(probes):
        st_of = {}
        for k, p in enumerate(probes):
            st_of.setdefault(int(p), set()).add(int(wc.te_station[k]))
        flag = np.zeros(n_st, dtype=bool)
        for p, sts in st_of.items():
            if len(sts) > 1:
                for j in sts:
                    flag[j] = True
        return flag

    return {"d_up": d_up, "d_lo": d_lo,
            "asym": np.abs(d_up - d_lo) / np.where(d_up + d_lo > 0,
                                                   d_up + d_lo, 1.0),
            "dz_up": dz_up, "dz_lo": dz_lo,
            "shared_up": _shared(wc.kutta_upper),
            "shared_lo": _shared(wc.kutta_lower)}


def pearson(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3 or a[m].std() == 0 or b[m].std() == 0:
        return float("nan")
    return float(np.corrcoef(a[m], b[m])[0, 1])
