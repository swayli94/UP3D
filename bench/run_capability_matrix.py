"""S5-cap: the solver CAPABILITY matrix, measured on today's HEAD.

Pre-registered in docs/dev_phase_two/20260802-2200-capability-matrix-prereg.md (plus the
same file's appendix recording the user's mesh-density axis and save-and-plot
requirement), committed before this file.

The question: for each (wake model x geometry x MESH LEVEL), how far up the Mach ladder
does the solver converge CLEANLY, and where does the envelope end?

Definitions, fixed in the pre-registration and implemented literally here:
  * CONVERGED  = driver converged AND n_limited == 0 AND n_floored == 0. A residual met
                 with clamps present is recorded CLAMPED and does NOT count -- the
                 project's clamp-not-silent contract (GS1.4 / B15).
  * ENVELOPE   = the last ladder point that is both converged and has M_max < 1.4, the
                 user's relaxation. M6 at M0.84 already reads M_max 2.1 from the tip
                 singularity, outside any full-potential model's validity, so "does it
                 converge inside M_local < 1.4" is the question actually asked.

Every recipe below is copied VERBATIM from a committed source, cited inline. Assembling
one by hand measures something else -- a mistake already made once this session.

Artifacts, per the user: section Cp curves and force coefficients are saved for every
successful point and plotted per cell. CSV rows are appended AS THEY ARE MEASURED so an
interruption never costs data; PNGs are generated per cell at the end so plotting never
blocks measurement.

Outputs (TRACKED): bench/gate_results/capability_matrix.csv
                   bench/gate_results/capability/<cell>_M<mach>_cp.csv
                   bench/gate_results/capability/<cell>.png
"""

import csv
import json
import os
import sys
import time
import traceback

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.fuselage import FuselageParams, make_inboard_clip  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te                      # noqa: E402
# ★ the UNIFIED entry points: one dispatch serving both wake paths
# (phi= for conforming, mvop=+phi_ext= for level-set), so "same extractor" is
# guaranteed by construction rather than by discipline. My first draft invented a
# `wall_cp_curve_mv` that does not exist -- checked against the source instead.
from pyfp3d.post.unified import section_cp, wall_forces as u_wall_forces  # noqa: E402
from pyfp3d.post.surface import (cl_kj_3d, planform_area,           # noqa: E402
                                 wall_force_coefficients)
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)
from tests.test_p8_newton import NEWTON_M6_RECIPE                   # noqa: E402

OUT = os.path.join(HERE, "gate_results")
ART = os.path.join(OUT, "capability")
os.makedirs(ART, exist_ok=True)
CSV = os.path.join(OUT, "capability_matrix.csv")

MACH_LADDER = (0.50, 0.60, 0.65, 0.70, 0.75, 0.78, 0.80, 0.82, 0.84)
MMAX_LIMIT = 1.4                       # the user's relaxation
LEG_BUDGET_S = 900.0                   # pre-registered
ALPHA_M6, ALPHA_NACA = 3.06, 1.25

# ---- recipes, verbatim from committed sources ----------------------------
# cases/demo/b18_wingbody_transonic/run_demo.py:139-156
CONF_SEED_KW = dict(farfield_spanwise_gamma=True, precond="direct",
                    direct_refactor_every=1000, n_newton_max=60)
CONF_RAMP_NK = dict(freeze_refresh_max=8, precond="direct",
                    direct_refactor_every=1000, n_newton_max=80,
                    farfield_spanwise_gamma=True)
CONF_TAPER = ("vanish_smooth", 0.05)
LS_RAMP_KW = dict(farfield="freestream", farfield_aux="pin_gamma",
                  freeze_tol=1e-4, freeze_max_clamped=8, intermediate_tol=1e-3,
                  n_seed=30, direct_refactor_every=1000, n_newton_max=80)
DM = 0.05
# tests/test_b22_ls_3d_anchors.py:52 (the M6-wing LS recipe)
LS_WING_KW = dict(farfield="neumann", n_seed=40, n_newton_max=80,
                  tol_residual=1e-10)
FUS = FuselageParams()


def append_row(row):
    head = not os.path.exists(CSV)
    keys = ["cell", "path", "geom", "level", "n_nodes", "n_tets", "m_inf",
            "alpha", "status", "converged", "n_limited", "n_floored",
            "res_final", "n_newton", "m_max", "cl_p", "cl_kj", "wall_s",
            # ★ added 2026-08-03 after the first matrix run: without these, DIFFERENT
            # failure modes were indistinguishable in the CSV, and one row was found
            # to be MISLABELLED. ls_naca_medium's "M0.78 NOT_CONVERGED" row actually
            # held a state at m_attained = 0.765 -- the ramp never reached the
            # requested Mach -- while res_final 1.58e-10 made it look like a
            # near-tolerance miss. A row must never be able to claim a Mach the
            # solver did not attain. (CLEAN rows were never affected: conv already
            # requires |m_final - m| < 1e-9. The damage was bounded to non-CLEAN rows,
            # which is why the envelope points stand.)
            # ★★ ERRATUM 2026-08-16 (GS4.0). Two things above went wrong when phase 3
            # deleted the level-set route, and both are now fixed:
            #  (1) this instrument was built for a defect found on the LEVEL-SET path
            #      using a field only the LEVEL-SET driver supplied, so deleting the
            #      route silently removed its input and "a row must never be able to
            #      claim a Mach the solver did not attain" became FALSE for the only
            #      surviving path (see the m_att branch in run_cell);
            #  (2) the parenthetical is wrong even as written: THIS script's `conv` is
            #      `bool(r["converged"])` and never looked at m_final -- that was the
            #      level-set driver's internal invariant, not ours. CLEAN rows happened
            #      to stay safe for a different reason (a conforming ramp that dies
            #      returns converged=False), and the CLEAN guard now checks explicitly.
            "m_attained", "accept_reason", "res_unfrozen", "f_final",
            "descent10", "note"]
    with open(CSV, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        if head:
            w.writeheader()
        w.writerow(row)


def classify(conv, n_lim, n_flr, m_max):
    """The pre-registered definitions, applied literally."""
    if not conv:
        return "NOT_CONVERGED"
    if n_lim or n_flr:
        return "CLAMPED"                      # residual met but not a solution
    return "CLEAN" if m_max < MMAX_LIMIT else "CLEAN_OVER_MMAX"


# ------------------------------------------------------------------ solvers
# ★ Each recipe has an m_start (the conforming wing's driver default is 0.70; the
# wing-body and the LS paths start at 0.50-0.60). A ladder point AT OR BELOW that is
# not a ramp case at all -- the driver rightly refuses to continue DOWNWARD
# (design.md Sec 12 risk 2) -- so it is solved single-level instead. That is correct
# usage, not a deviation from the recipe: the same newton_kw is passed either way.
# Caught by smoke-testing the runner on the cheapest cell before launching the matrix.
CONF_WING_MSTART, WB_MSTART, LS_WING_MSTART = 0.70, 0.50, 0.60


def conf_wing(mesh_path, m, alpha):
    mc, wc = cut_wake(read_mesh(mesh_path))
    if m <= CONF_WING_MSTART:
        r = solve_newton_lifting(mc, wc, m_inf=m, alpha_deg=alpha,
                                 **NEWTON_M6_RECIPE["newton_kw"])
    else:
        r = solve_newton_transonic(mc, wc, m_inf=m, alpha_deg=alpha,
                                   **NEWTON_M6_RECIPE)
    return mc, wc, r, np.asarray(r["phi"]), None


def conf_wingbody(mesh_path, m, alpha):
    mc, wc = cut_wake(read_mesh(mesh_path))
    taper = tip_taper_factors(wc.station_z, B_SEMI, CONF_TAPER[0],
                              CONF_TAPER[1] * B_SEMI)
    if m <= WB_MSTART:
        r = solve_newton_lifting(mc, wc, m_inf=m, alpha_deg=alpha,
                                 **dict(CONF_RAMP_NK,
                                        kutta_estimator="pressure",
                                        tip_taper=taper))
        return mc, wc, r, np.asarray(r["phi"]), None
    seed = solve_newton_lifting(mc, wc, m_inf=WB_MSTART, alpha_deg=alpha,
                                **CONF_SEED_KW)
    r = solve_newton_transonic(
        mc, wc, m_inf=m, alpha_deg=alpha, m_start=0.50, dm=DM, dm_min=0.01,
        freeze_tol=1e-5, intermediate_tol=1e-4,
        newton_kw=dict(CONF_RAMP_NK, kutta_estimator="pressure",
                       tip_taper=taper, phi_init=seed["phi"],
                       gamma_init=seed["gamma"], n_picard_seed=0))
    return mc, wc, r, np.asarray(r["phi"]), None


# ---------------------------------------------------------------------------
# ★ LEVEL-SET CELLS REMOVED 2026-08-10 (phase 3 task 1, ruling D5). What went:
# _ls_single / _ls_op / ls_wing / ls_wingbody plus the 10 "level-set" CONFIGS rows.
#
# ★★ ERRATUM OBLIGATION, stated here because a number depends on it: the
# capability boundary quotes this script's "13 configurations x ladder = 78 points,
# 72/78 converged, 0/78 clamps". That reading INCLUDED the level-set cells, so the
# committed CSV stays valid as a PRE-DELETION measurement while this script can no
# longer reproduce it -- going forward it measures the conforming cells only. The
# boundary is annotated accordingly; the number was not silently re-scoped.
# ---------------------------------------------------------------------------
#: (cell, path, geom, mesh dir, level, alpha, solver)
#: ★ ORDERED coarse -> medium -> fine, deliberately: this run is hours long, and if it
#: is interrupted the matrix should be COMPLETE at one refinement level rather than
#: having a few geometries finished and the rest blank. Mesh density is an axis (the
#: user's requirement), and the level-to-level comparison is only meaningful if every
#: geometry has the same levels measured.
_CELLS_UNORDERED = [
    ("conf_naca_coarse", "conforming", "naca2.5d", "naca0012_2.5d", "coarse",
     ALPHA_NACA, conf_wing),
    ("conf_naca_medium", "conforming", "naca2.5d", "naca0012_2.5d", "medium",
     ALPHA_NACA, conf_wing),
    ("conf_naca_fine", "conforming", "naca2.5d", "naca0012_2.5d", "fine",
     ALPHA_NACA, conf_wing),
    ("conf_wing_coarse", "conforming", "m6wing", "onera_m6", "coarse",
     ALPHA_M6, conf_wing),
    ("conf_wing_medium", "conforming", "m6wing", "onera_m6", "medium",
     ALPHA_M6, conf_wing),
    ("conf_wb_coarse", "conforming", "wingbody", "onera_m6_wingbody_conforming",
     "coarse", ALPHA_M6, conf_wingbody),
    ("conf_wb_medium", "conforming", "wingbody", "onera_m6_wingbody_conforming",
     "medium", ALPHA_M6, conf_wingbody),
    #: ★ xcoarse added 2026-08-03. Effort moved off `fine`: it costs 345-1561 s per flow
    #: point against coarse's 2.7-11.7 s, meaningless against the 2 CPU-min product
    #: target, and it adds a third mesh level to conforming NACA -- the ONE combination
    #: that already had three -- while five of the six geometry x wake-path combinations
    #: had only two and so admitted no convergence order at all.
    #:
    #: h_wall = 0.044, the largest BOTH wing-body families accept. One level coarser
    #: (0.060) does not lose accuracy, it BREAKS the mesh: min dihedral 19.50 -> 0.70 deg
    #: and aspect ratio 108, because TIP_CAP_RADIUS = 0.0222 is a FIXED geometric scale
    #: that h_wall 0.060 exceeds by 2.7x. Wing-body ONLY -- the NACA and M6-wing families
    #: clamp h_far (min(3.0, 150h) and min(2.5, 120h)), so their xcoarse would have a
    #: first interval that does not refine the far field AT ALL.
    ("conf_wb_xcoarse", "conforming", "wingbody",
     "onera_m6_wingbody_conforming", "xcoarse", ALPHA_M6, conf_wingbody),
    #: ★ the other two geometries' cheap levels, added after measuring each family's own
    #: h_far clamp and meshability floor rather than reusing the wing-body's numbers.
    #: NACA: h_wall 0.040 with the clamp OFF (clamped, h_far would pin at 3.0 exactly as
    #: at coarse => a first interval with NO far-field refinement). This family has no
    #: dihedral gate and correctly so -- a one-prism-layer 2.5-D extrusion runs aspect
    #: ~100 at EVERY level including the shipped coarse, so the 3-D bound does not apply.
    ("conf_naca_xcoarse", "conforming", "naca2.5d", "naca0012_2.5d", "xcoarse",
     ALPHA_NACA, conf_wing),
    #: M6 wing: the shipped `coarse` IS clamped (the P13/M1b defect), so it is OFF the
    #: refinement ray and cannot sit in a convergence ladder. The valid cheap ladder is
    #: {xcoarse_ss, coarse_ss, medium} -- hence coarse_ss cells here as well, not just
    #: xcoarse. Uniform 2x in every length.
    ("conf_wing_xcoarse_ss", "conforming", "m6wing", "onera_m6", "xcoarse_ss",
     ALPHA_M6, conf_wing),
    ("conf_wing_coarse_ss", "conforming", "m6wing", "onera_m6", "coarse_ss",
     ALPHA_M6, conf_wing),
]
_LEVEL_ORDER = {"xcoarse": 0, "xcoarse_ss": 0, "coarse": 1, "coarse_ss": 1,
                "medium": 2, "fine": 3}
CELLS = sorted(_CELLS_UNORDERED, key=lambda c: (_LEVEL_ORDER[c[4]], c[0]))


def measure(cell, path, geom, mdir, level, alpha, fn, m):
    mesh_path = os.path.join(REPO, "cases", "meshes", mdir, f"{level}.msh")
    if not os.path.exists(mesh_path):
        return dict(cell=cell, path=path, geom=geom, level=level, m_inf=m,
                    alpha=alpha, status="MESH_MISSING", note=mesh_path), None
    t0 = time.perf_counter()
    try:
        mesh, op, r, phi, mvop = fn(mesh_path, m, alpha)
    except Exception as exc:                                        # noqa: BLE001
        return dict(cell=cell, path=path, geom=geom, level=level, m_inf=m,
                    alpha=alpha, status="ERROR",
                    wall_s=round(time.perf_counter() - t0, 1),
                    note=f"{type(exc).__name__}: {exc}"), None
    wall = time.perf_counter() - t0
    try:
        return _postprocess(cell, path, geom, level, alpha, m, wall,
                            mesh, op, r, phi, mvop)
    except Exception as exc:                                        # noqa: BLE001
        # ★ post-processing gets the SAME discipline as the solve. The first version
        # left it outside the try, and a stale import in the level-set branch
        # (wall_forces lives in post.unified, not post.surface) killed the whole
        # multi-hour run at the first LS cell instead of recording one bad row.
        return dict(cell=cell, path=path, geom=geom, level=level, m_inf=m,
                    alpha=alpha, status="POSTPROC_ERROR",
                    wall_s=round(wall, 1),
                    note=f"{type(exc).__name__}: {exc}"), None


def _postprocess(cell, path, geom, level, alpha, m, wall, mesh, op, r, phi,
                 mvop):
    if mvop is None:                                   # conforming
        mc, wc = mesh, op
        s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        f = wall_force_coefficients(mc.nodes, mc.elements,
                                    mc.boundary_faces["wall"], phi,
                                    alpha_deg=alpha, s_ref=s_ref, m_inf=m)
        o = np.argsort(wc.station_z)
        clkj = (float(cl_kj_3d(np.atleast_1d(r["gamma"])[o], wc.station_z[o],
                               s_ref, B_SEMI)) if geom != "naca2.5d" else
                float("nan"))
        from pyfp3d.mesh.metrics import precompute_element_geometry
        from pyfp3d.physics.isentropic import mach_number_squared
        B, _ = precompute_element_geometry(mc.nodes, mc.elements)
        g = np.einsum("eaj,ea->ej", B, phi[mc.elements])
        m_max = float(np.sqrt(mach_number_squared(
            np.einsum("ej,ej->e", g, g), m).max()))
        n_lim, n_flr = int(r["n_limited"]), int(r["n_floored"])
        res = float(r["residual_history"][-1]); nn = int(r["n_newton"])
        conv = bool(r["converged"]); nodes, tets = len(mc.nodes), len(mc.elements)
        geom_obj = (mc, wc, None)
    else:                                              # level-set: REMOVED
        #: ★ phase 3 task 1 (ruling D5): the level-set branch of this dispatch is
        #: gone with the route, and so are the 10 "level-set" CONFIGS rows. Kept as
        #: an explicit raise rather than deleted silently, because the "path" column
        #: still exists in the committed CSV and a stale caller passing "level-set"
        #: must fail loudly instead of being scored as a conforming cell.
        raise ValueError(
            f"path={path!r}: the level-set route was deleted in phase 3 "
            f"(ruling D5). Only conforming cells are measurable now; the "
            f"pre-deletion 78-point matrix is in the committed CSV.")
    status = classify(conv, n_lim, n_flr, m_max)
    #: the Mach the solver ACTUALLY attained -- see the append_row note. Recorded
    #: unconditionally so a non-converged row can never again be read as if it held a
    #: state at the requested Mach.
    #:
    #: ★★ GS4.0 (2026-08-16), and the reason is worth stating in full. This line used to
    #: read `float(r.get("m_last_converged", r.get("m_final", m)))`. Both keys existed
    #: ONLY on the level-set driver (`newton_ls.py`), which phase 3 DELETED -- so from
    #: that day the chain fell through to `m`, i.e. m_att was IDENTICALLY the requested
    #: Mach, and the CLEAN guard below (`abs(m_att - m) > 1e-9`) could not fire. A guard
    #: written specifically so as NOT to trust the driver flags had become a guard that
    #: trusted them completely, silently. Found by the 2026-08-16 independent audit
    #: (docs/inspection/20260816-2200-independent-audit-zh.md, §7.1); the fields are now
    #: backported to the CONFORMING ramp (`newton.py::_ramp_honesty_fields`).
    #:
    #: Two branches, each stated rather than defaulted:
    #:  - RAMP driver: the keys are present and must be read. No `.get` default -- if a
    #:    future refactor drops them again this must raise, not degrade quietly. That is
    #:    the entire lesson of this defect.
    #:  - SINGLE-Mach driver (`solve_newton_lifting`, used by the subsonic cells): the
    #:    returned state is at `m` BY CONSTRUCTION -- there is no ramp to fall short of --
    #:    so m_att = m is the correct reading here, not a fallback.
    #: `m_final` (not `m_last_converged`) is the label, because this column labels the
    #: row's OWN cl_p / cl_kj / m_max, and those are computed from `r["phi"]`, which is
    #: the state at `m_final`. A ramp that died at 0.80 after converging 0.75 returns the
    #: FAILED 0.80 state, so 0.80 is what the row's numbers belong to.
    if "m_final" in r:
        m_att = r["m_final"]
    else:
        m_att = float(m)
    #: is the last-10-step residual descent still steep? This is the cheap
    #: discriminator between "ran out of iteration budget while converging" and
    #: "genuinely stalled" -- the two were indistinguishable in the first matrix, and
    #: they are completely different findings.
    #: ★ residual_history is the per-NEWTON-STEP history and exists on BOTH paths. The
    #: first version used the level-set path's `levels` list instead, which is per RAMP
    #: LEVEL -- a handful of entries -- so descent10 silently came back None on every
    #: level-set row, i.e. the diagnostic was dead exactly where it was first needed.
    hist = np.asarray(r.get("residual_history", [])).ravel()
    d10 = (float(hist[-11] / hist[-1]) if len(hist) >= 11 and hist[-1] > 0
           else float("nan"))
    fh = r.get("F_history")
    row = dict(cell=cell, path=path, geom=geom, level=level, n_nodes=nodes,
               n_tets=tets, m_inf=m, alpha=alpha, status=status,
               converged=conv, n_limited=n_lim, n_floored=n_flr,
               res_final=res, n_newton=nn, m_max=round(m_max, 5),
               cl_p=round(f["cl"], 8), cl_kj=round(clkj, 8),
               wall_s=round(wall, 1), m_attained=m_att,
               accept_reason=r.get("accept_reason"),
               res_unfrozen=r.get("residual_unfrozen"),
               f_final=(float(np.asarray(fh).ravel()[-1]) if fh is not None
                        and len(np.asarray(fh).ravel()) else None),
               descent10=(round(d10, 3) if d10 == d10 else None), note="")
    #: a CLEAN row that did not attain its own Mach would invalidate the envelope
    #: table, so make that unrepresentable rather than trusting the driver flags to
    #: stay consistent with each other.
    #: ★ GS4.0: `m_att is None` (no level ran at all) also lands here -- an unknown
    #: attained Mach on a CLEAN row is exactly the unrepresentable case, and a bare
    #: `abs(None - m)` would raise inside the reporting layer, which is how a 40-minute
    #: solve was lost once. The old parenthetical "CLEAN rows were never affected: conv
    #: already requires |m_final - m| < 1e-9" (see append_row) described the LEVEL-SET
    #: driver's internals -- this script's `conv` is just `bool(r["converged"])`, so the
    #: check has to be made here, not assumed upstream.
    if status.startswith("CLEAN") and (m_att is None or abs(m_att - m) > 1e-9):
        row["status"] = "MACH_NOT_ATTAINED"
        row["note"] = f"clean at {m_att} but {m} was requested"
    return row, (geom_obj, phi, f)


def save_cp(cell, m, geom, payload):
    """Section Cp -- the artifact the user asked for as next-phase reference.

    Uses pyfp3d.post.unified.section_cp.

    ★ GS4.0 erratum 2026-08-16: this used to read "for BOTH paths, so the
    conforming and level-set curves in this matrix are produced by the same
    dispatch". There is only ONE path since ruling D5, and `post/unified.py`
    was collapsed onto its conforming half in phase 3 -- the sentence had
    outlived the thing it described.
    """
    (obj, phi, _f) = payload
    mesh, op, mf = obj
    kw = (dict(phi=phi) if mf is None else dict(mvop=op, phi_ext=phi))
    m_eff = m if mf is None else mf
    stations = ([("z_mid", None)] if geom == "naca2.5d"
                else [(f"{e:.2f}", e) for e in (0.20, 0.44, 0.65, 0.80, 0.90)])
    rows = []
    for label, eta in stations:
        try:
            if eta is None:
                z = 0.5 * float(np.ptp(mesh.nodes[:, 2]))
                c = section_cp(mesh, z=z, m_inf=m_eff, **kw)
            else:
                c = section_cp(mesh, eta=eta, b_semi=B_SEMI, m_inf=m_eff, **kw)
        except Exception as exc:                                  # noqa: BLE001
            print(f"      (Cp at {label} failed: {type(exc).__name__}: {exc})",
                  flush=True)
            continue
        for side in ("upper", "lower"):
            for x, cp in zip(c[f"x_{side}"], c[f"cp_{side}"]):
                rows.append((label, side, float(x), float(cp)))
    if not rows:
        return
    p = os.path.join(ART, f"{cell}_M{str(m).replace('.', '')}_cp.csv")
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["station", "side", "x_over_c", "cp"])
        w.writerows(rows)


def plot_cell(cell):
    """One PNG per cell: Cp at every Mach measured. Generated AFTER the cell's
    ladder so plotting never blocks or endangers measurement."""
    import glob
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    files = sorted(glob.glob(os.path.join(ART, f"{cell}_M*_cp.csv")))
    if not files:
        return
    data = {}
    for fp in files:
        mach = os.path.basename(fp).split("_M")[-1].split("_cp")[0]
        d = {}
        for r in csv.DictReader(open(fp)):
            d.setdefault((r["station"], r["side"]), []).append(
                (float(r["x_over_c"]), float(r["cp"])))
        data[mach] = d
    stations = sorted({k[0] for d in data.values() for k in d})
    fig, axes = plt.subplots(1, len(stations), figsize=(4.2 * len(stations), 3.6),
                             squeeze=False, sharey=True)
    cmap = plt.get_cmap("viridis")
    for j, st in enumerate(stations):
        ax = axes[0][j]
        for i, (mach, d) in enumerate(sorted(data.items())):
            col = cmap(i / max(len(data) - 1, 1))
            for side, ls in (("upper", "-"), ("lower", "--")):
                pts = sorted(d.get((st, side), []))
                if pts:
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], ls,
                            color=col, lw=1.2,
                            label=f"M0.{mach[1:]}" if side == "upper" else None)
        ax.invert_yaxis(); ax.grid(alpha=0.3)
        ax.set_xlabel("x/c"); ax.set_title(f"{cell}  station {st}")
        if j == 0:
            ax.set_ylabel("Cp"); ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(ART, f"{cell}.png"), dpi=120)
    plt.close(fig)


def done_cells():
    """Cells whose ladder already TERMINATED in the CSV, so a resumed run skips them.
    A cell is done if it has a non-CLEAN row (the ladder stopped) or a row at the
    ladder's last Mach (it ran the whole way). Added after a post-processing bug
    killed a multi-hour run -- resume must be cheap."""
    if not os.path.exists(CSV):
        return set()
    seen = {}
    for r in csv.DictReader(open(CSV)):
        seen.setdefault(r["cell"], []).append(r)
    out = set()
    for c, rs in seen.items():
        if any(r["status"] != "CLEAN" for r in rs) or \
           any(abs(float(r["m_inf"]) - MACH_LADDER[-1]) < 1e-12 for r in rs):
            out.add(c)
    return out


def main():
    only = [c for c in os.environ.get("PYFP3D_CAP_CELLS", "").split(",") if c]
    already = done_cells()
    if already:
        print(f"resuming -- already terminated: {sorted(already)}", flush=True)
    summary = {}
    for cell, path, geom, mdir, level, alpha, fn in CELLS:
        if only and cell not in only:
            continue
        if cell in already and not only:
            continue
        print(f"\n=== {cell}  ({path}, {geom}, {level}, alpha {alpha}) ===",
              flush=True)
        envelope, stopped = None, None
        for m in MACH_LADDER:
            row, payload = measure(cell, path, geom, mdir, level, alpha, fn, m)
            append_row(row)
            st = row["status"]
            print(f"  M{m:<5} {st:16s} conv={row.get('converged')} "
                  f"lim/flr={row.get('n_limited')}/{row.get('n_floored')} "
                  f"M_max={row.get('m_max')} cl_p={row.get('cl_p')} "
                  f"|R|={row.get('res_final')} ({row.get('wall_s')}s)"
                  f"{'  ' + row['note'] if row.get('note') else ''}", flush=True)
            if st in ("MESH_MISSING", "ERROR"):
                stopped = st; break
            if payload is not None and st.startswith("CLEAN"):
                save_cp(cell, m, geom, payload)
            if st == "CLEAN":
                envelope = m
            else:
                stopped = st; break
            if (row.get("wall_s") or 0) > LEG_BUDGET_S:
                stopped = "OVER_BUDGET"; break
        try:
            plot_cell(cell)
        except Exception as exc:                              # noqa: BLE001
            print(f"  (plot failed: {type(exc).__name__}: {exc})")
        summary[cell] = dict(envelope=envelope, stopped=stopped)
        print(f"  -> envelope M{envelope} (stopped by {stopped})", flush=True)
    with open(os.path.join(OUT, "capability_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("\n=== ENVELOPE SUMMARY (clean convergence with M_max < 1.4) ===")
    for c, v in summary.items():
        print(f"  {c:20s} M{v['envelope']}   stopped by {v['stopped']}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc(); sys.exit(1)
