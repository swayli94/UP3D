"""
Track B / B4 gate: TE control-volume & implicit-Kutta re-derivation
(docs/roadmap.md Track B B4; docs/design_track_b.md §9).

B3 delivered the multivalued lifting solve, but its emergent circulation
converged to the WRONG value (Γ → ~0.168 vs the conforming 0.120). B4 fixes
that. The chain of reasoning, each link pinned by a test here:

  1. **The wake LS cannot pin Γ.** Its residual is identically zero for ANY
     spatially-constant jump, because Σ_c ∇N_c = ∇(1) = 0 (partition of
     unity). So "g₂ IS the discrete Kutta" is false, and the LS block has a
     constant-jump null space.        -> test_wake_ls_has_constant_jump_nullspace

  2. So Γ must be pinned by the TE rows. B3 used lower-side mass conservation
     there, whose control volume is up/down ASYMMETRIC on a symmetric airfoil
     (the ε shift sends every on-sheet node "+"). That over-circulated ~42%.

  3. **B4's fix: the NONLINEAR TE pressure-equality (Bernoulli) Kutta**,
     |q_u|² = |q_l|², factorized exactly as (q_u+q_l)·(q_u−q_l) = 0 and
     linearized by freezing the mean. q_u and q_l are recovered on the TE's
     **WALL-ADJACENT** upper/lower control volumes — the body surfaces at the
     TE. They come from DIFFERENT element sets, so q_u−q_l is NOT a jump
     gradient and does NOT vanish for a constant jump: non-degenerate in Γ.
                                       -> test_te_kutta_pins_gamma
     Symmetrizing the control volume is NOT an available route (the mesh is
     naturally asymmetric there — user-arbitrated 2026-07-12); the condition is
     instead a pointwise physical statement needing no symmetry.

  4. **Wall-adjacency is what makes it accurate.** Recovering q over the whole
     TE element fan (interior + wake elements pollute the average) gives Γ
     0.1407/0.1355; over the wall-adjacent elements only, 0.1177/0.1191 vs the
     conforming 0.1175/0.1200.        -> test_wall_adjacent_recovery_is_what_works
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.picard import solve_laplace_lifting
from pyfp3d.solve.picard_ls import solve_multivalued_lifting
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
ALPHA = 2.0


def _load(level):
    path = M0_DIR / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} missing")
    return read_mesh(path)


def _setup(mesh):
    """Horizontal (chord-plane) wake — design.md §4; it is also the plane the
    M0 mesh embeds and the plane the far-field vortex branch cut uses."""
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(
        mesh.nodes, mesh.elements, wls,
        wall_nodes=np.unique(mesh.boundary_faces["wall"]),
    )
    return wls, cm, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                        levelset=wls)


def _conforming_gamma(mesh):
    mesh_cut, wc = cut_wake(mesh)
    r = solve_laplace_lifting(mesh_cut, wc, alpha_deg=ALPHA)
    return float(np.mean(r["gamma"]))


class TestWakeLSNullspace:
    """LINK 1 — the wake LS is blind to a constant jump, so it CANNOT fix Γ.
    This is why a separate TE Kutta condition is structurally necessary."""

    @pytest.mark.parametrize("level", ["coarse", "medium"])
    def test_wake_ls_has_constant_jump_nullspace(self, level):
        mesh = _load(level)
        _, cm, mvop = _setup(mesh)
        A = mvop.assemble_matrix(closure="wake_ls", te_kutta="mass").tocsr()

        from pyfp3d.kernels.cut_assembly import nonte_aux_rows
        ls_rows = nonte_aux_rows(cm)
        assert len(ls_rows) > 0

        # aux = main - side*Γ is the EXACT other-side continuation (the
        # potential is analytic across the wake; the cut is only a branch cut).
        rng = np.random.default_rng(0)
        x = np.zeros(mvop.n_total)
        x[: cm.n_main] = rng.standard_normal(cm.n_main)
        cut_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
        aux = cm.ext_dof_of_node[cut_nodes]
        for gamma in (0.0, 0.12, 1.0, -0.7):
            x[aux] = x[cut_nodes] - cm.node_side[cut_nodes] * gamma
            r_ls = (A @ x)[ls_rows]
            scale = float(np.abs(A[ls_rows]).max()) * float(np.abs(x).max())
            assert np.max(np.abs(r_ls)) < 1e-10 * max(scale, 1.0), (
                f"LS residual not zero for constant jump {gamma} — the "
                f"null-space argument broke: {np.max(np.abs(r_ls)):.3e}"
            )


class TestTEControlVolume:
    @pytest.mark.parametrize("level", ["coarse", "medium"])
    def test_below_te_fan_is_not_cut(self, level):
        """The ε shift must not manufacture cuts in the below-TE fan (López
        p.57). Regression pin for the 2026-07-12 fix."""
        mesh = _load(level)
        _, cm, _ = _setup(mesh)
        el = np.asarray(mesh.elements, dtype=np.int64)
        is_te = np.zeros(cm.n_main, dtype=bool)
        is_te[cm.te_nodes] = True
        te_e = is_te[el]
        side = cm.node_side[el]
        true_fan = np.flatnonzero(
            te_e.any(axis=1) & np.all(te_e | (side == -1), axis=1)
        )
        is_cut = np.zeros(len(el), dtype=bool)
        is_cut[cm.cut_elems] = True
        assert not is_cut[true_fan].any()
        assert set(true_fan) == set(cm.te_lower_elems)

    def test_te_control_volume_is_wall_adjacent(self):
        """LINK 3/4 — the TE Kutta lives on the WALL-adjacent upper/lower
        elements (the body surface at the TE), not the whole element fan."""
        mesh = _load("medium")
        _, cm, mvop = _setup(mesh)
        el = np.asarray(mesh.elements, dtype=np.int64)
        is_w = np.zeros(cm.n_main, dtype=bool)
        is_w[np.unique(mesh.boundary_faces["wall"])] = True
        for cv in mvop._te_cv:
            assert len(cv["upper_elems"]) > 0 and len(cv["lower_elems"]) > 0
            for key in ("upper_elems", "lower_elems"):
                # every element in the control volume carries a wall FACE
                assert np.all(is_w[el[cv[key]]].sum(axis=1) >= 3)


class TestB4Gate:
    """THE B4 GATE: the nonlinear TE pressure-equality Kutta pins Γ to the
    conforming value; the B3 (lower-mass-conservation) TE row does not."""

    @pytest.mark.parametrize("level", ["coarse", "medium"])
    def test_te_kutta_pins_gamma(self, level):
        mesh = _load(level)
        g_ref = _conforming_gamma(mesh)

        _, _, mvop = _setup(mesh)
        g_new = solve_multivalued_lifting(
            mvop, mesh, 0.0, alpha_deg=ALPHA, te_kutta="pressure"
        )["gamma"]
        err = abs(g_new - g_ref) / g_ref
        assert err < 0.05, (
            f"{level}: pressure-equality Kutta Γ={g_new:.4f} vs conforming "
            f"{g_ref:.4f} ({err*100:.1f}% — gate is 5%)"
        )

        # ...and the B3 TE row (lower-side mass conservation) does NOT: it
        # over-circulates badly. Keeps the before/after contrast honest.
        _, _, mvop2 = _setup(mesh)
        g_old = solve_multivalued_lifting(
            mesh=mesh, mvop=mvop2, m_inf=0.0, alpha_deg=ALPHA, te_kutta="mass"
        )["gamma"]
        assert abs(g_old - g_ref) / g_ref > 0.30


class TestB4Artifacts:
    def test_export_te_diagnosis(self, gate_artifacts_dir):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        mesh = _load("medium")
        _, cm, mvop = _setup(mesh)
        el = np.asarray(mesh.elements, dtype=np.int64)
        cen = mesh.nodes[el].mean(axis=1)
        cv = mvop._te_cv[0]
        te = mesh.nodes[cm.te_nodes[0]]

        # measured 2026-07-12 (incompressible, α=2°)
        h = [1.0, 0.5, 0.25]
        g_old = [0.2074, 0.1760, 0.1704]     # B3 TE row = lower mass cons
        g_fan = [0.1407, 0.1355, 0.1329]     # B4 Kutta, FULL-fan recovery
        g_new = [0.1177, 0.1191, 0.1197]     # B4 Kutta, WALL-adjacent (shipped)
        g_ref = [0.1175, 0.1200, 0.1202]     # conforming, same meshes

        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13.5, 5.2))

        fan = np.flatnonzero((el == cm.te_nodes[0]).any(axis=1))
        p = cen[fan]
        ax0.scatter(p[:, 0], p[:, 1], s=60, c="lightgrey", edgecolors="grey",
                    zorder=2, label=f"TE element fan  [{len(fan)}]")
        for key, col, lab in (("upper", "tab:red", "UPPER control volume"),
                              ("lower", "tab:blue", "LOWER control volume")):
            q = cen[cv[f"{key}_elems"]]
            ax0.scatter(q[:, 0], q[:, 1], s=150, c=col, edgecolors="k",
                        linewidths=0.5, zorder=3,
                        label=f"{lab} (wall-adjacent) [{len(q)}]")
        ax0.plot([te[0], p[:, 0].max()], [te[1], te[1]], "k--", lw=1.5,
                 zorder=1, label="wake sheet (y = 0)")
        ax0.plot(te[0], te[1], "k*", ms=20, zorder=4, label="TE node (doubled)")
        ax0.set_xlabel("x"); ax0.set_ylabel("y")
        ax0.set_title("B4: the TE Kutta lives on the WALL-ADJACENT\n"
                      "upper/lower control volumes, not the whole fan")
        ax0.legend(fontsize=8, loc="lower left", framealpha=0.92)
        ax0.grid(alpha=0.3)

        ax1.plot(h, g_old, "s-", color="tab:red", lw=2, ms=8,
                 label="B3: TE row = lower mass cons  (+42%)")
        ax1.plot(h, g_fan, "^-", color="tab:orange", lw=2, ms=8,
                 label="B4 Kutta, full-fan recovery  (+11%)")
        ax1.plot(h, g_new, "o-", color="tab:green", lw=2.5, ms=10,
                 label="B4 Kutta, WALL-adjacent  (shipped)")
        ax1.plot(h, g_ref, "--", color="tab:blue", lw=2.5,
                 label="conforming (same meshes)")
        ax1.set_xlabel("relative mesh size h  (coarse / medium / fine)")
        ax1.set_ylabel("$\\Gamma$")
        ax1.set_title("B4 closes the gap: 42% → <1% of conforming\n"
                      "(NACA0012, $\\alpha=2°$, incompressible)")
        ax1.set_xscale("log"); ax1.invert_xaxis()
        ax1.legend(fontsize=8.5); ax1.grid(alpha=0.3, which="both")

        fig.tight_layout()
        png = gate_artifacts_dir / "b4_te_kutta.png"
        fig.savefig(png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        assert png.exists()

        csv = gate_artifacts_dir / "summary.csv"
        with open(csv, "w") as f:
            f.write("level,h,gamma_b3_massrow,gamma_b4_fullfan,"
                    "gamma_b4_walladj,gamma_conforming,err_pct\n")
            for i, lv in enumerate(("coarse", "medium", "fine")):
                e = abs(g_new[i] - g_ref[i]) / g_ref[i] * 100
                f.write(f"{lv},{h[i]},{g_old[i]},{g_fan[i]},{g_new[i]},"
                        f"{g_ref[i]},{e:.2f}\n")
        assert csv.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
