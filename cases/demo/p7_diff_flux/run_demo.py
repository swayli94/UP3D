"""
P7 demo -- frozen-selection differentiability of the walk flux (gate G7.3).

P7 is the P8 fully-coupled-Newton prerequisite: the exact Jacobian
(design.md Sec 6.3 / Lopez Appendix B) needs d(rho_tilde)/d(phi) at FROZEN
upstream selection u(e). The shipped P4 walk flux is already differentiable
there -- P7 derives the branch-wise sensitivities (s_e, s_u) =
(d rho_tilde_e / d q2_e, d rho_tilde_e / d q2_u) (kernels/upwind.py::
rho_tilde_sensitivities_sweep) and verifies them against a directional
central difference of the SHIPPED rho_tilde_sweep with u(e) held frozen.
The forward flux is untouched (byte-identical): G7.1/G7.2 hold by
construction and stay locked by the G4.2/G4.1/G4.3 tests.

  part 1 (always, seconds, no solve): constructed multi-regime field on the
      coarse NACA0012 2.5D mesh. V7.1 analytic-vs-FD scatter + per-regime
      rel-err histogram (the G7.3 evidence); V7.2 frozen-selection regime map.
  part 2 (PYFP3D_TRANSONIC_GATES=1, ~6 min): real NACA0012 M=0.80 alpha=1.25
      coarse walk solution (P4 G4.1 case). V7.3 regime map over the actual
      supersonic pocket; V7.4 FD accuracy on the converged field -- the
      building block is exact on the true P8 target state.

Kink note: rho_tilde is C0 but not C1 exactly AT the max(nu_e, nu_u) tie and
at the switch threshold M2 = M_c2 (the measure-zero locus of design.md Sec
3.1; Lopez's Newton converges because the active set freezes there). A
central difference straddling a kink inside its eps-stencil returns a branch
average, not a derivative, so such elements are excluded from the rel-err
statistic and their count is reported (a handful out of ~16k).

Headless; writes results/*.png + checks.csv; exits nonzero on unexpected FAIL.
"""

import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, INK_2, MESH_DIR, MUTED, S1_BLUE, S2_AQUA, S3_YELLOW,
    S4_ROSE, CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402

from pyfp3d.kernels.jacobian import PicardOperator  # noqa: E402
from pyfp3d.kernels.upwind import UpwindOperator, rho_tilde_sweep  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.physics.isentropic import (  # noqa: E402
    GAMMA, density_field, mach_number_squared,
)

OUT = Path(__file__).resolve().parent / "results"
M_CRIT = 0.95
UPWIND_C = 1.5
RHO_FLOOR = 0.05
EPS = 5e-7
KINK_GUARD = 3e-5
REL_TOL = 1e-6           # G7.3 acceptance

REGIME_COLORS = {"subsonic": S1_BLUE, "accelerating": S3_YELLOW,
                 "shock-point": S4_ROSE, "self-upstream": MUTED}


def _frozen_rho_tilde(q2, m_inf, upstream):
    """The SHIPPED flux (density law + rho_tilde_sweep, byte-identical
    physics) evaluated at a frozen selection."""
    rho = density_field(q2, m_inf)
    nu = np.empty_like(q2)
    rt = np.empty_like(q2)
    rho_tilde_sweep(q2, rho, upstream, m_inf, M_CRIT, UPWIND_C, GAMMA,
                    RHO_FLOOR, nu, rt)
    return rt


def _regime_masks(q2, upstream, m_inf):
    m2 = mach_number_squared(q2, m_inf, GAMMA)
    m2u = m2[upstream]
    mc2 = M_CRIT ** 2
    nu_e = UPWIND_C * np.maximum(0.0, 1.0 - mc2 / np.maximum(m2, mc2))
    nu_u = UPWIND_C * np.maximum(0.0, 1.0 - mc2 / np.maximum(m2u, mc2))
    self_up = upstream == np.arange(len(q2))
    masks = {
        "self-upstream": self_up,
        "subsonic": (~self_up) & (np.maximum(nu_e, nu_u) == 0.0),
        "accelerating": (~self_up) & (nu_e >= nu_u) & (nu_e > 0.0),
        "shock-point": (~self_up) & (nu_u > nu_e),
    }
    # FD-valid mask: exclude the eps-neighbourhood of the two C0 kinks
    # (irrelevant for self-upstream elements: their upwind jump is 0).
    near_kink = (~self_up) & (
        ((np.maximum(nu_e, nu_u) > 0.0) & (np.abs(nu_e - nu_u) < KINK_GUARD))
        | (np.abs(m2 - mc2) < KINK_GUARD) | (np.abs(m2u - mc2) < KINK_GUARD))
    return masks, ~near_kink


def _jvp_vs_fd(nodes, elements, phi, m_inf, seed=0):
    """Analytic JVP from (s_e, s_u) vs central-difference of the shipped
    flux at frozen selection, for one random direction."""
    op = PicardOperator(nodes, elements)
    upw = UpwindOperator(nodes, elements, weighted=False)
    grad, q2 = op.velocities(phi)
    grad, q2 = grad.copy(), q2.copy()
    rho = density_field(q2, m_inf)
    s_e, s_u, upstream = upw.rho_tilde_sensitivities(
        grad, q2, rho, m_inf, UPWIND_C, M_CRIT, rho_floor=RHO_FLOOR)
    s_e, s_u, upstream = s_e.copy(), s_u.copy(), upstream.copy()

    rng = np.random.default_rng(seed)
    delta = rng.standard_normal(len(phi))
    delta /= np.abs(delta).max()
    gradd, _ = op.velocities(delta)
    gradd = gradd.copy()
    dq2 = 2.0 * np.einsum("ij,ij->i", grad, gradd)
    jvp = s_e * dq2 + s_u * dq2[upstream]

    _, q2p = op.velocities(phi + EPS * delta)
    rt_p = _frozen_rho_tilde(q2p.copy(), m_inf, upstream)
    _, q2m = op.velocities(phi - EPS * delta)
    rt_m = _frozen_rho_tilde(q2m.copy(), m_inf, upstream)
    fd = (rt_p - rt_m) / (2.0 * EPS)

    masks, fd_valid = _regime_masks(q2, upstream, m_inf)
    return {"jvp": jvp, "fd": fd, "q2": q2, "upstream": upstream,
            "s_e": s_e, "s_u": s_u, "masks": masks, "fd_valid": fd_valid}


def _rel_err(res):
    scale = np.abs(res["fd"][res["fd_valid"]]).max()
    rel = np.abs(res["jvp"] - res["fd"]) / scale
    return rel, scale


def _scatter_and_hist(res, tag, title, fname, cl, gate_name):
    rel, scale = _rel_err(res)
    v = res["fd_valid"]
    max_rel = rel[v].max()
    n_excl = int((~v).sum())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4))
    for name, m in res["masks"].items():
        mm = m & v
        if mm.sum():
            ax1.plot(res["fd"][mm], res["jvp"][mm], ".", ms=3,
                     color=REGIME_COLORS[name], label=f"{name} ({mm.sum()})")
    lim = np.abs(res["fd"][v]).max()
    ax1.plot([-lim, lim], [-lim, lim], "-", lw=0.8, color=INK_2, zorder=0)
    ax1.set_xlabel("central difference of shipped flux (frozen u)")
    ax1.set_ylabel("analytic JVP from (s_e, s_u)")
    ax1.set_title("analytic vs FD, per element")
    ax1.legend()

    floor = 1e-14
    for name, m in res["masks"].items():
        mm = m & v
        if mm.sum():
            ax2.hist(np.log10(np.maximum(rel[mm], floor)), bins=40,
                     histtype="step", lw=1.6, color=REGIME_COLORS[name],
                     label=name)
    ax2.axvline(np.log10(REL_TOL), color=CRITICAL, lw=1.2, ls="--",
                label="G7.3 tol 1e-6")
    ax2.set_xlabel("log10 rel err (field-scale normalized)")
    ax2.set_ylabel("elements")
    ax2.set_title(f"max rel err {max_rel:.2e} ({n_excl} kink-adjacent excluded)")
    ax2.legend()
    fig.suptitle(title, fontsize=12, fontweight="semibold")
    finish(fig, OUT, fname)

    cl.add("G7.3", gate_name, f"{max_rel:.2e}", "< 1e-6", max_rel < REL_TOL,
           note=f"{n_excl} kink-adjacent elements excluded (C0 locus)")
    return max_rel, n_excl


def _regime_map(nodes, elements, res, dz, title, fname, xlim=None, ylim=None):
    cent = nodes[elements].mean(axis=1)
    mid = np.abs(cent[:, 2] - 0.5 * dz) < 0.3 * dz   # midspan layer
    fig, ax = plt.subplots(figsize=(9, 6))
    for name in ("subsonic", "accelerating", "shock-point", "self-upstream"):
        m = res["masks"][name] & mid
        if m.sum():
            ax.plot(cent[m, 0], cent[m, 1], ".", ms=2.2,
                    color=REGIME_COLORS[name], label=f"{name} ({m.sum()})")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    if xlim: ax.set_xlim(*xlim)
    if ylim: ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.legend(markerscale=4, loc="upper right")
    finish(fig, OUT, fname)


def part1_constructed(cl: CheckList):
    print("\n[1/2] constructed multi-regime field on NACA0012 coarse (no solve)")
    mesh = read_mesh(MESH_DIR / "naca0012_2.5d" / "coarse.msh")
    nodes, elements = mesh.nodes, mesh.elements
    dz = float(np.ptp(nodes[:, 2]))
    x = nodes[:, 0]
    rng = np.random.default_rng(42)
    # oscillating streamwise gradient -> sub/supersonic bands with both
    # accelerating and decelerating (shock-point) frozen-selection branches;
    # small nodal noise breaks the prism-split gradient degeneracies that
    # would otherwise park element pairs exactly on the max(nu_e,nu_u) tie.
    phi = 1.25 * x + 0.35 * np.sin(1.7 * x) + 1e-3 * rng.standard_normal(len(x))
    m_inf = 0.8

    res = _jvp_vs_fd(nodes, elements, phi, m_inf)
    max_rel, n_excl = _scatter_and_hist(
        res, "constructed", "V7.1 frozen-walk d(rho_tilde)/d(phi): analytic vs "
        "central difference -- constructed multi-regime field",
        "v71_fd_scatter_constructed.png", cl, "FD rel err (constructed field)")
    for name in ("subsonic", "accelerating", "shock-point"):
        cl.add("G7.3", f"regime populated: {name}", int(res["masks"][name].sum()),
               "> 0", res["masks"][name].sum() > 0)
    frac_excl = n_excl / len(res["fd_valid"])
    cl.add("G7.3", "kink-adjacent excluded fraction", f"{frac_excl:.2%}",
           "< 2% (measure-zero locus)", frac_excl < 0.02)

    _regime_map(nodes, elements, res, dz,
                "V7.2 frozen-selection regime map (midspan, constructed field)",
                "v72_regime_map_constructed.png")
    write_csv(OUT, "g73_constructed.csv", "quantity,value",
              [("max_rel_err", f"{max_rel:.3e}"), ("eps", EPS),
               ("n_kink_excluded", n_excl), ("n_elements", len(res["fd_valid"])),
               *[(f"n_{k}", int(v.sum())) for k, v in res["masks"].items()]])


def part2_converged(cl: CheckList):
    print("\n[2/2] real NACA0012 M=0.80 alpha=1.25 coarse walk solution [heavy]")
    from pyfp3d.mesh.wake_cut import cut_wake
    from pyfp3d.solve.continuation import solve_transonic_lifting

    mc, wc = cut_wake(read_mesh(MESH_DIR / "naca0012_2.5d" / "coarse.msh"))
    dz = float(np.ptp(mc.nodes[:, 2]))
    r = solve_transonic_lifting(mc, wc, m_inf=0.80, alpha_deg=1.25,
                                max_gamma_evals=12, n_picard_eval=800)
    res = _jvp_vs_fd(mc.nodes, mc.elements, r["phi"], 0.80)
    max_rel, n_excl = _scatter_and_hist(
        res, "converged", "V7.4 frozen-walk d(rho_tilde)/d(phi) on the CONVERGED "
        "G4.1 field (the P8 Newton target state)",
        "v74_fd_scatter_converged.png", cl, "FD rel err (converged G4.1 field)")
    _regime_map(mc.nodes, mc.elements, res, dz,
                "V7.3 frozen-selection regimes over the real supersonic pocket "
                "(NACA0012 M0.80 coarse, midspan)",
                "v73_regime_map_converged.png",
                xlim=(-0.3, 1.5), ylim=(-0.6, 0.8))
    n_super = int((res["masks"]["accelerating"] | res["masks"]["shock-point"]).sum())
    cl.add("G7.3", "supersonic pocket populated (converged)", n_super, "> 0",
           n_super > 0)
    write_csv(OUT, "g73_converged.csv", "quantity,value",
              [("max_rel_err", f"{max_rel:.3e}"), ("n_kink_excluded", n_excl),
               ("mach_max", f"{np.sqrt(r['mach2_max']):.4f}"),
               ("n_supersonic_frozen", n_super),
               *[(f"n_{k}", int(v.sum())) for k, v in res["masks"].items()]])


def main():
    apply_style()
    cl = CheckList("P7 differentiable walk flux at frozen selection (G7.3)")
    part1_constructed(cl)
    if os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1":
        part2_converged(cl)
    else:
        print("\n[2/2] converged-field check skipped (~6 min transonic solve); "
              "set PYFP3D_TRANSONIC_GATES=1 to run")
    sys.exit(cl.report(OUT))


if __name__ == "__main__":
    main()
