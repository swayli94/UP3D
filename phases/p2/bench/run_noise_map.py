"""S2 opening probe: WHERE does the element-gradient noise live?

GS1b.10 sec 9 registered "shock localisation" as its own problem and named the shared
root with G1.6: element-wise CONSTANT P1 gradients on irregular tets. The S2 opening
round file (20260731-1200) then made a precision point -- the wall recovery is
POST-PROCESSING and fixes only the wall Cp, whereas localisation needs a VOLUME
gradient. So whether S2's anisotropic prism layers help BOTH problems or only one
depends on a fact nobody has measured: where the bad cells actually are.

This measures it, cheaply, on a SMOOTH SUBSONIC state (M 0.50, no supersonic zone at
all) so that every bit of cross-face |dM| is discretisation noise rather than physics.
For each interior face it records the jump and the two candidate explanatory
variables:

    wall distance   -- if the outliers sit in the wall layer, prism layers are one
                       fix for two problems
    element shape    -- aspect ratio and min dihedral angle; if the outliers track
                       shape ANYWHERE in the volume, the fix is reconstruction (or
                       mesh quality), and prisms will not reach the shock at x/c 0.6

★ ERRATUM, same day, written before any attribution was made. The first version of
this script used the raw cross-face |dM| as the noise measure and a pre-registered
discriminator on the top-1 % faces' median wall distance. It read 0.022 / 0.019 and
said "wall-layer object" -- and that reading is CONFOUNDED, which the numbers
themselves gave away: the all-face median wall distance is 1.84 chords (the count is
far-field dominated), while only 4.5 % / 2.6 % of the outliers are actually inside the
first wall layer. A first difference across a face is h * |grad M| plus noise, and in
a SMOOTH subsonic field |grad M| is genuinely largest near the body -- so a large |dM|
near the wall is exactly what a perfectly clean discretisation would also produce.
That measure cannot separate noise from physics anywhere.

The fix is a noise measure with the smooth part removed: for each element, fit a
linear model of M over its node-neighbour patch EXCLUDING the element itself, predict
at its centroid, and take the deviation. A clean discretisation of a smooth field
gives a deviation of order the truncation error; a bad cell deviates from what its own
neighbourhood says, which is precisely what makes a shock sensor fire spuriously. That
is also, not coincidentally, the reconstruction operator that would be the fix -- so
the measurement and the candidate cure share the machinery. The raw |dM| columns are
KEPT beside it, labelled confounded, because they are what GS1b.10 sec 9 quoted
(max |dM| 0.2186) and the two need to be comparable.

Discriminator (restated on the new measure, before the run): the top-1 % elements'
median wall distance against the all-element median, plus the fraction inside the
first wall layer. Prisms are a shared fix only if the outliers are genuinely a
wall-layer population; if they track element shape out in the volume, the fix is
reconstruction and prisms will not reach the shock at x/c 0.6.

Outputs (TRACKED, bench/results/ is gitignored):
    bench/gate_results/noise_map.csv      per-level statistics + the binned table
    bench/gate_results/noise_map.png      the noise vs wall distance, and where the
                                          top 1 % sit
"""

import csv
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
import numba                                                        # noqa: E402
import numpy as np                                                  # noqa: E402
from scipy.spatial import cKDTree                                   # noqa: E402

HERE = Path(__file__).resolve().parent
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = HERE.parent
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.metrics import (build_face_adjacency,              # noqa: E402
                                 compute_aspect_ratios,
                                 compute_min_dihedral_angles,
                                 precompute_element_geometry)
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting                # noqa: E402

OUT = _GATE
OUT.mkdir(exist_ok=True)

M_INF, ALPHA = 0.50, 1.25          # smooth: subcritical everywhere
LEVELS = ("coarse", "medium")
TOP_FRAC = 0.01                    # the outlier definition, fixed before the run
GAMMA = 1.4
#: the sensor's regulariser, as a fraction of M_inf. Its job is to kill the 0/0 in
#: the far field, where BOTH the second difference and the smooth variation vanish
#: and an unregularised ratio would rank uniform-freestream cells as the worst in
#: the mesh. Same role as the epsilon in the Jameson pressure switch.
SENSOR_EPS_FRAC = 1e-3


@numba.njit(cache=True)
def _patch_linear_dev(cent, mach, h, elements, n2e_off, n2e_idx, dev, gmag):
    """Leave-one-out node-patch linear fit, per element.

    For element e: fit M ~ a + g . (x - x_e)/h_e by least squares over the elements
    sharing a node with e, EXCLUDING e itself, then dev[e] = |M_e - a| and
    gmag[e] = |g| (the smooth variation across one cell size). Offsets are
    non-dimensionalised by h_e so all four normal-equation blocks are O(count) and
    one relative ridge conditions them uniformly -- this matters on the 2.5-D meshes,
    where the span-direction offsets are set by the single layer's thickness and not
    by the local cell size.

    Elements appearing in several of e's four node lists are accumulated several
    times. That is left deliberately: it weights a neighbour by how many nodes it
    shares with e (3 = face neighbour, down to 1 = touches one corner), which is the
    proximity weighting the fit wants anyway.
    """
    A = np.zeros((4, 4))
    b = np.zeros(4)
    v = np.zeros(4)
    for e in range(elements.shape[0]):
        A[:, :] = 0.0
        b[:] = 0.0
        cnt = 0
        inv_h = 1.0 / h[e]
        for a in range(4):
            nd = elements[e, a]
            for p in range(n2e_off[nd], n2e_off[nd + 1]):
                j = n2e_idx[p]
                if j == e:
                    continue
                v[0] = 1.0
                v[1] = (cent[j, 0] - cent[e, 0]) * inv_h
                v[2] = (cent[j, 1] - cent[e, 1]) * inv_h
                v[3] = (cent[j, 2] - cent[e, 2]) * inv_h
                for r in range(4):
                    for c in range(4):
                        A[r, c] += v[r] * v[c]
                    b[r] += v[r] * mach[j]
                cnt += 1
        if cnt < 8:
            dev[e] = np.nan          # recorded, never silently zero
            gmag[e] = np.nan
            continue
        for r in range(4):
            A[r, r] *= 1.0 + 1e-10
            if A[r, r] <= 0.0:
                A[r, r] = 1e-30
        sol = np.linalg.solve(A, b)
        ok = True
        for r in range(4):
            if not np.isfinite(sol[r]):
                ok = False
        if not ok:
            dev[e] = np.nan
            gmag[e] = np.nan
            continue
        dev[e] = abs(mach[e] - sol[0])
        gmag[e] = np.sqrt(sol[1] ** 2 + sol[2] ** 2 + sol[3] ** 2)


def patch_noise(cent, mach, V, elements, n_nodes):
    """(dev, gmag, sensor): the smooth-part-removed noise and its bounded ratio."""
    counts = np.bincount(elements.reshape(-1), minlength=n_nodes)
    off = np.zeros(n_nodes + 1, dtype=np.int64)
    np.cumsum(counts, out=off[1:])
    order = np.argsort(elements.reshape(-1), kind="stable")
    idx = (order // 4).astype(np.int64)
    dev = np.empty(len(elements))
    gmag = np.empty(len(elements))
    _patch_linear_dev(np.ascontiguousarray(cent), np.ascontiguousarray(mach),
                      np.cbrt(np.maximum(V, 1e-300)),
                      np.ascontiguousarray(elements.astype(np.int64)),
                      off, idx, dev, gmag)
    eps = SENSOR_EPS_FRAC * M_INF
    sensor = dev / (gmag + dev + eps)
    return dev, gmag, sensor


def element_mach(nodes, elements, phi, m_inf, gamma=GAMMA):
    """Element Mach from the P1 gradient (the same q2 the upwinding sees)."""
    B, V = precompute_element_geometry(nodes, elements)
    grad = np.einsum("eaj,ea->ej", B, phi[elements])
    q2 = np.einsum("ej,ej->e", grad, grad)
    # isentropic: a^2 = 1/M_inf^2 + (gamma-1)/2 * (1 - q2), M^2 = q2 / a^2
    a2 = 1.0 / m_inf ** 2 + 0.5 * (gamma - 1.0) * (1.0 - q2)
    return np.sqrt(np.maximum(q2 / np.maximum(a2, 1e-30), 0.0)), V, grad


def main():
    rows, panels = [], []
    for level in LEVELS:
        path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
        if not path.exists():
            print(f"skip {level}: mesh missing")
            continue
        mc, wc = cut_wake(read_mesh(path))
        t0 = time.perf_counter()
        r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                 upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6,
                                 freeze_refresh_max=8, precond="direct",
                                 direct_refactor_every=4, n_newton_max=80)
        wall_s = time.perf_counter() - t0
        res = float(r["residual_history"][-1])
        print(f"\n{level}: {len(mc.nodes)} nodes, {len(mc.elements)} tets -- "
              f"|R|={res:.2e} conv={r['converged']} ({wall_s:.0f}s)", flush=True)
        if not r["converged"]:
            print("  NOT converged -- recorded, the noise read would be "
                  "contaminated by the iterate")

        mach, V, _ = element_mach(mc.nodes, mc.elements, np.asarray(r["phi"]),
                                 M_INF)
        n_sup = int(np.count_nonzero(mach > 1.0))
        print(f"  M_max = {mach.max():.4f}, supersonic elements = {n_sup} "
              f"(the premise: this must be 0 for the jumps to be pure noise)")

        # ---- PRIMARY: the smooth-part-removed noise, per element -----------
        cent = mc.nodes[mc.elements].mean(axis=1)
        dev, gmag, sensor = patch_noise(cent, mach, V, mc.elements,
                                        len(mc.nodes))
        n_bad = int(np.count_nonzero(~np.isfinite(dev)))
        good = np.isfinite(dev)
        if n_bad:
            print(f"  {n_bad} elements had too small a patch to fit -- "
                  f"excluded and RECORDED, not zeroed")

        # ---- SECONDARY (confounded, kept for comparability) ----------------
        nb, _ = build_face_adjacency(mc.elements)
        e_idx = np.repeat(np.arange(len(mc.elements)), 4)
        n_idx = nb.reshape(-1)
        keep = n_idx > e_idx                       # each interior face once
        dm = np.abs(mach[e_idx[keep]] - mach[n_idx[keep]])

        # ---- the candidate explanatory variables ---------------------------
        wall_nodes = np.unique(mc.boundary_faces["wall"].reshape(-1))
        d_wall = cKDTree(mc.nodes[wall_nodes]).query(cent)[0]
        ar = compute_aspect_ratios(mc.nodes, mc.elements)
        dih = compute_min_dihedral_angles(mc.nodes, mc.elements)
        h = np.cbrt(np.maximum(V, 1e-300))

        k = max(1, int(TOP_FRAC * int(good.sum())))
        gi = np.flatnonzero(good)
        top = gi[np.argpartition(sensor[gi], -k)[-k:]]
        med = lambda a: float(np.median(a))                          # noqa: E731

        ratio_dwall = med(d_wall[top]) / med(d_wall[gi])
        ratio_ar = med(ar[top]) / med(ar[gi])
        ratio_dih = med(dih[top]) / med(dih[gi])
        # in the first wall layer = closer to the wall than its own cell size
        in_layer = float(np.mean(d_wall[top] < h[top]))
        n_fire = int(np.count_nonzero(sensor[gi] > 0.5))
        print(f"  sensor median {med(sensor[gi]):.5f}   p99 "
              f"{float(np.quantile(sensor[gi], 0.99)):.5f}   max "
              f"{float(sensor[gi].max()):.5f}   > 0.5: {n_fire} elements")
        print(f"  (confounded, for comparison) face |dM| median {med(dm):.5f} "
              f"p99 {float(np.quantile(dm, 0.99)):.5f} max {dm.max():.5f}")
        print(f"  top {100*TOP_FRAC:g}% ({k} elements) vs all-element medians:")
        print(f"    wall distance   {med(d_wall[top]):.5f} / "
              f"{med(d_wall[gi]):.5f} = {ratio_dwall:.3f}")
        print(f"    aspect ratio    {med(ar[top]):.3f} / {med(ar[gi]):.3f} "
              f"= {ratio_ar:.3f}")
        print(f"    min dihedral    {med(dih[top]):.2f} / {med(dih[gi]):.2f} "
              f"= {ratio_dih:.3f} deg")
        print(f"    inside the first wall layer: {100*in_layer:.1f} %")

        # ---- correlated or EXPLAINED? ---------------------------------------
        # A 1.4x median shift is a bias, not an attribution. Two harder reads:
        # the enrichment (of the worst-shaped 1 %, what fraction is also in the
        # noisiest 1 %? -- 1.0 means shape carries no information, 100 means it
        # carries all of it) and the rank correlation over the whole population.
        def spearman(a, b):
            ra = np.argsort(np.argsort(a)).astype(np.float64)
            rb = np.argsort(np.argsort(b)).astype(np.float64)
            ra -= ra.mean()
            rb -= rb.mean()
            return float(ra @ rb / np.sqrt((ra @ ra) * (rb @ rb)))

        worst_ar = gi[np.argpartition(ar[gi], -k)[-k:]]
        worst_dih = gi[np.argpartition(-dih[gi], -k)[-k:]]
        enrich_ar = float(np.isin(worst_ar, top).mean()) / TOP_FRAC
        enrich_dih = float(np.isin(worst_dih, top).mean()) / TOP_FRAC
        rho_ar = spearman(sensor[gi], ar[gi])
        rho_dih = spearman(sensor[gi], -dih[gi])
        rho_dw = spearman(sensor[gi], -d_wall[gi])
        rho_h = spearman(sensor[gi], h[gi])
        rho_g = spearman(sensor[gi], gmag[gi])
        print(f"  enrichment (worst-shape 1% that is also noisiest 1%): "
              f"aspect {enrich_ar:.1f}x, dihedral {enrich_dih:.1f}x "
              f"(1.0x = shape carries no information)")
        print(f"  Spearman rho vs the sensor: aspect {rho_ar:+.3f}, "
              f"1/dihedral {rho_dih:+.3f}, wall proximity {rho_dw:+.3f}, "
              f"cell size {rho_h:+.3f}, local |grad| {rho_g:+.3f}")

        # ★ a control on MY OWN metric, not on the mesh. The regulariser eps
        # forces sensor -> 0 wherever gmag << eps, i.e. across the whole far
        # field -- which manufactures a monotone sensor/|grad| and
        # sensor/wall-proximity relation out of nothing. So repeat the reads on
        # the NEAR FIELD only, where gmag >> eps and the ratio is meaningful. If
        # the correlations survive they are real; if they collapse they were mine.
        near = gi[d_wall[gi] < 1.0]
        gm_med = med(gmag[near])
        eps = SENSOR_EPS_FRAC * M_INF
        nk = max(1, int(TOP_FRAC * near.size))
        ntop = near[np.argpartition(sensor[near], -nk)[-nk:]]
        nworst_ar = near[np.argpartition(ar[near], -nk)[-nk:]]
        n_enrich = float(np.isin(nworst_ar, ntop).mean()) / TOP_FRAC
        n_rho_ar = spearman(sensor[near], ar[near])
        n_rho_g = spearman(sensor[near], gmag[near])
        n_rho_dw = spearman(sensor[near], -d_wall[near])
        print(f"  NEAR FIELD only (d_wall < 1 c, {near.size} elements, median "
              f"|grad| {gm_med:.4g} vs eps {eps:.4g} = "
              f"{gm_med/eps:.0f}x): aspect {n_rho_ar:+.3f}, local |grad| "
              f"{n_rho_g:+.3f}, wall proximity {n_rho_dw:+.3f}, "
              f"enrichment {n_enrich:.1f}x")
        # and the operationally relevant one: would spurious detections land
        # WHERE the shock goes? The M0.80 upper-surface shock sits near x/c 0.6.
        band = ((cent[ntop, 0] > 0.4) & (cent[ntop, 0] < 0.8)
                & (np.abs(cent[ntop, 1]) < 0.2))
        frac_band = float(band.mean())
        print(f"  of those near-field outliers, {100*frac_band:.1f} % sit in "
              f"the shock band 0.4 < x/c < 0.8, |y| < 0.2")

        # ★★ DUAL-MEASURE agreement, because the eps calibration decides WHICH
        # cells count as noisiest and this round has already caught it doing that
        # twice (the far field, then the low-gradient approach region upstream of
        # the leading edge, which holds 38-47 % of the near-field outliers). So
        # repeat the two decision reads on the raw absolute deviation `dev`, which
        # has NO regulariser at all -- its own confound is that dev ~ h^2 |grad^2 M|
        # and so favours high-curvature regions, a DIFFERENT bias. Where the two
        # measures agree, the conclusion does not depend on my metric choice; where
        # they disagree, nothing is concluded.
        dtop = near[np.argpartition(dev[near], -nk)[-nk:]]
        d_in_layer = float(np.mean(d_wall[dtop] < h[dtop]))
        d_rho_dw = spearman(dev[near], -d_wall[near])
        d_worst_ar = near[np.argpartition(ar[near], -nk)[-nk:]]
        d_enrich = float(np.isin(d_worst_ar, dtop).mean()) / TOP_FRAC
        print(f"  DUAL CHECK on raw dev (no eps): {100*d_in_layer:.1f} % in the "
              f"wall layer, wall proximity {d_rho_dw:+.3f}, shape enrichment "
              f"{d_enrich:.1f}x")
        agree_wall = (d_rho_dw < 0.3) == (n_rho_dw < 0.3)
        print(f"  => the two measures {'AGREE' if agree_wall else 'DISAGREE'} on "
              f"the wall-layer question")

        rows.append(dict(
            enrich_ar=round(enrich_ar, 3), enrich_dih=round(enrich_dih, 3),
            rho_ar=round(rho_ar, 4), rho_dih=round(rho_dih, 4),
            rho_dwall=round(rho_dw, 4), rho_h=round(rho_h, 4),
            rho_gmag=round(rho_g, 4),
            near_n=int(near.size), near_gmag_med=round(gm_med, 8),
            near_rho_ar=round(n_rho_ar, 4),
            near_rho_gmag=round(n_rho_g, 4),
            near_rho_dwall=round(n_rho_dw, 4),
            near_enrich_ar=round(n_enrich, 3),
            near_frac_in_shock_band=round(frac_band, 4),
            dev_frac_in_wall_layer=round(d_in_layer, 4),
            dev_rho_dwall=round(d_rho_dw, 4),
            dev_enrich_ar=round(d_enrich, 3), dual_agree_wall=bool(agree_wall),
            level=level, n_nodes=len(mc.nodes), n_tets=len(mc.elements),
            m_inf=M_INF, res_final=res, converged=bool(r["converged"]),
            m_max=round(float(mach.max()), 5), n_supersonic=n_sup,
            n_unfittable=n_bad,
            sensor_median=round(med(sensor[gi]), 6),
            sensor_p99=round(float(np.quantile(sensor[gi], 0.99)), 6),
            sensor_max=round(float(sensor[gi].max()), 6),
            n_sensor_gt_0p5=n_fire,
            dev_median=round(med(dev[gi]), 8),
            dev_max=round(float(dev[gi].max()), 8),
            n_faces=int(dm.size), dm_median_confounded=round(med(dm), 6),
            dm_p99_confounded=round(float(np.quantile(dm, 0.99)), 6),
            dm_max_confounded=round(float(dm.max()), 6),
            top_k=k,
            top_dwall_med=round(med(d_wall[top]), 6),
            all_dwall_med=round(med(d_wall[gi]), 6),
            ratio_dwall=round(ratio_dwall, 4),
            top_ar_med=round(med(ar[top]), 4), all_ar_med=round(med(ar[gi]), 4),
            ratio_ar=round(ratio_ar, 4),
            top_dih_med=round(med(dih[top]), 3),
            all_dih_med=round(med(dih[gi]), 3), ratio_dih=round(ratio_dih, 4),
            frac_top_in_wall_layer=round(in_layer, 4),
            wall_s=round(wall_s, 1)))

        # ---- binned table: is there a wall-distance decile that owns it? ----
        istop = np.zeros(len(mc.elements), dtype=bool)
        istop[top] = True
        edges = np.quantile(d_wall[gi], np.linspace(0.0, 1.0, 11))
        edges[-1] *= 1.0 + 1e-9
        bi = np.clip(np.digitize(d_wall, edges) - 1, 0, 9)
        for b in range(10):
            s = good & (bi == b)
            if not s.any():
                continue
            rows.append(dict(
                level=f"{level}:dwall_decile{b}", n_tets=int(s.sum()),
                all_dwall_med=round(med(d_wall[s]), 6),
                sensor_median=round(med(sensor[s]), 6),
                sensor_p99=round(float(np.quantile(sensor[s], 0.99)), 6),
                sensor_max=round(float(sensor[s].max()), 6),
                dev_median=round(med(dev[s]), 8),
                all_ar_med=round(med(ar[s]), 4),
                frac_of_top=round(float(istop[s].sum()) / k, 5)))
        panels.append((level, d_wall[gi], sensor[gi], cent[gi],
                       istop[gi], ar[gi]))

    if not panels:
        return 1
    with open(OUT / "noise_map.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "noise_map.csv")

    fig, axes = plt.subplots(2, len(panels), figsize=(6.4 * len(panels), 9.5),
                             squeeze=False)
    for j, (level, dw, sen, cen, istop, arr) in enumerate(panels):
        ax = axes[0][j]
        ax.semilogx(dw, sen, ".", ms=1, alpha=0.15, color="0.6",
                    label="all elements")
        ax.semilogx(dw[istop], sen[istop], ".", ms=3, color="crimson",
                    label=f"top {100*TOP_FRAC:g} %")
        ax.axhline(0.5, color="navy", lw=0.8, ls="--",
                   label="a sensor threshold of 0.5")
        ax.set_xlabel("distance to the wall")
        ax.set_ylabel("noise sensor  dev / (|grad| + dev + eps)")
        ax.set_title(f"{level}: smooth-part-removed noise vs wall distance\n"
                     f"(M{M_INF}, no supersonic element anywhere)")
        ax.legend(loc="best", fontsize=8)
        ax.grid(alpha=0.3, which="both")

        ax = axes[1][j]
        ax.plot(cen[:, 0], cen[:, 1], ".", ms=1, alpha=0.08, color="0.7")
        ax.plot(cen[istop, 0], cen[istop, 1], ".", ms=3.5, color="crimson")
        ax.set_xlim(-0.4, 1.6)
        ax.set_ylim(-0.6, 0.6)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(f"{level}: where the top {100*TOP_FRAC:g} % sit "
                     f"(x/c 0.6 is where the shock goes)")
    fig.tight_layout()
    fig.savefig(OUT / "noise_map.png", dpi=130)
    print("wrote", OUT / "noise_map.png")

    print("\n=== reading (the discriminator was restated before this run) ===")
    for r in rows:
        if ":" in str(r["level"]):
            continue
        rd, ra, il = (r["ratio_dwall"], r["ratio_ar"],
                      r["frac_top_in_wall_layer"])
        print(f"  {r['level']:7s} wall-distance ratio {rd:.3f}, "
              f"aspect-ratio ratio {ra:.3f}, {100*il:.1f} % in the wall layer, "
              f"sensor max {r['sensor_max']:.3f}, "
              f"{r['n_sensor_gt_0p5']} elements > 0.5")
        # ★ the verdict reads the NEAR-FIELD columns only. The global rho_gmag
        # and rho_dwall (about +0.8 both levels) were measured this round to be an
        # artefact of this script's own regulariser -- inside the near field they
        # collapse to about +0.15 and about 0. Attributing on them would be
        # attributing to my metric.
        print(f"          near-field: wall proximity {r['near_rho_dwall']:+.3f}, "
              f"shape enrichment {r['near_enrich_ar']:.1f}x, "
              f"{100*r['near_frac_in_shock_band']:.1f} % in the shock band")
        print(f"          (global rho_dwall {r['rho_dwall']:+.3f} / rho_gmag "
              f"{r['rho_gmag']:+.3f} are THIS SCRIPT'S eps artefact, not mesh "
              f"facts -- see the near-field control)")
        print(f"          dual check on raw dev (no eps): wall proximity "
              f"{r['dev_rho_dwall']:+.3f}, {100*r['dev_frac_in_wall_layer']:.1f} "
              f"% in the wall layer, shape enrichment "
              f"{r['dev_enrich_ar']:.1f}x -> "
              f"{'AGREE' if r['dual_agree_wall'] else 'DISAGREE'}")
        if not r["dual_agree_wall"]:
            print("    => the two noise measures DISAGREE: conclude NOTHING "
                  "about the wall-layer question from this run")
        elif r["near_rho_dwall"] > 0.3 or il > 0.5:
            print("    => WALL-LAYER population: prisms are one fix for two "
                  "problems")
        else:
            print("    => NOT a wall-layer population (no wall trend inside the "
                  "near field, and almost none of the outliers is in the first "
                  "layer) => prisms buy WALL ACCURACY ONLY; shock localisation "
                  "needs its own fix, in the volume")
            if r["near_enrich_ar"] > 3.0:
                print(f"       shape carries REAL but PARTIAL information "
                      f"({r['near_enrich_ar']:.1f}x enrichment in the near "
                      f"field, rank correlation only "
                      f"{r['near_rho_ar']:+.3f}): most badly-shaped cells are "
                      f"quiet and most noisy cells are not badly shaped, which "
                      f"is GS1b.10 sec 9's obstruction with a number on it")
    return 0


if __name__ == "__main__":
    sys.exit(main())
