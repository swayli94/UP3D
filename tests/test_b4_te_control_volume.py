"""
Track B / B4 gate: TE control-volume & implicit-Kutta re-derivation
(docs/roadmap.md Track B B4; docs/design_track_b.md section 9 "OPEN DESIGN
PROBLEM").

B3 delivered the multivalued lifting solve, but the emergent circulation
converges to the WRONG value (NACA0012 a=2, incompressible: Gamma 0.207 ->
0.176 -> 0.170 on coarse/medium/fine vs the conforming/thin-airfoil 0.120).
This module is the DIAGNOSTIC harness for the re-derivation: it pins the two
structural facts that localize the defect, and emits the visual artifact the
re-derivation is steered by.

The two facts (design_track_b.md section 9):

  1. **The wake LS cannot pin Gamma.** Its residual is IDENTICALLY zero for
     any spatially-constant jump, because sum_c grad(N_c) = grad(1) = 0
     (partition of unity). So the LS block has a constant-jump null space,
     and "g2 IS the discrete Kutta condition" (the old design_track_b.md
     section 2.3 claim) is FALSE. Test: test_wake_ls_has_constant_jump_nullspace.

  2. **Gamma is therefore pinned by a SINGLE equation** -- the TE aux row
     (lower-side mass conservation at the trailing edge) -- and that row's
     control volume is UP/DOWN ASYMMETRIC on a symmetric airfoil (upper fan
     9 elements vs lower fan 6 on the medium NACA mesh), because the eps
     side-shift sends every on-sheet node "+". Test:
     test_te_control_volume_is_asymmetric.

Fact 2 is the standing hypothesis for the ~45% over-circulation; the B4 gate
is to re-derive the TE control volume so the discrete Kutta is consistent
(and, per the roadmap, to make the emergent Gamma converge to the conforming
value). These tests DOCUMENT the defect: they will need updating when B4
lands, which is intentional -- they are the before/after pins.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"

# Measured 2026-07-12 (incompressible, alpha = 2 deg), conforming = 0.1200.
CONFORMING_GAMMA = 0.1200
MEASURED_GAMMA = {"coarse": 0.2074, "medium": 0.1760, "fine": 0.1704}


def _load(level):
    path = M0_DIR / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} missing")
    return read_mesh(path)


def _setup(mesh):
    """Horizontal (chord-plane) wake -- design.md Sec 4; it is also the plane
    the M0 mesh embeds and the plane the far-field vortex branch cut uses."""
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(
        mesh.nodes, mesh.elements, wls,
        wall_nodes=np.unique(mesh.boundary_faces["wall"]),
    )
    return wls, cm


def _te_fan_census(mesh, cm, te_node):
    """Classify the element fan around a TE node: upper (all '+', plain),
    lower (below-TE fan, TE reference = aux), cut (sheet passes through)."""
    el = np.asarray(mesh.elements, dtype=np.int64)
    fan = np.flatnonzero((el == te_node).any(axis=1))
    is_cut = np.zeros(len(el), dtype=bool)
    is_cut[cm.cut_elems] = True
    is_low = np.zeros(len(el), dtype=bool)
    is_low[cm.te_lower_elems] = True
    upper = fan[~is_cut[fan] & ~is_low[fan]]
    return {
        "fan": fan,
        "upper": upper,
        "lower": fan[is_low[fan]],
        "cut": fan[is_cut[fan]],
    }


class TestWakeLSNullspace:
    """FACT 1: the wake LS is blind to a constant jump, so it cannot fix
    Gamma. This is the analytic partition-of-unity result, pinned
    numerically."""

    @pytest.mark.parametrize("level", ["coarse", "medium"])
    def test_wake_ls_has_constant_jump_nullspace(self, level):
        mesh = _load(level)
        wls, cm = _setup(mesh)
        mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
        A = mvop.assemble_matrix(closure="wake_ls").tocsr()

        from pyfp3d.kernels.cut_assembly import nonte_aux_rows
        ls_rows = nonte_aux_rows(cm)
        assert len(ls_rows) > 0

        # A field with a spatially CONSTANT jump: aux = main - side * Gamma.
        # (The potential is analytic across the wake -- the cut is only a
        # branch cut -- so this IS the exact "other side" continuation.)
        rng = np.random.default_rng(0)
        x = np.zeros(mvop.n_total)
        x[: cm.n_main] = rng.standard_normal(cm.n_main)  # arbitrary main field
        cut_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
        aux = cm.ext_dof_of_node[cut_nodes]
        for gamma in (0.0, 0.12, 1.0, -0.7):
            x[aux] = x[cut_nodes] - cm.node_side[cut_nodes] * gamma
            r_ls = (A @ x)[ls_rows]
            # scale-free: compare against the row norms
            scale = float(np.abs(A[ls_rows]).max()) * float(np.abs(x).max())
            assert np.max(np.abs(r_ls)) < 1e-10 * max(scale, 1.0), (
                f"LS residual not zero for a constant jump {gamma}: "
                f"{np.max(np.abs(r_ls)):.3e} -- the null-space argument broke"
            )


class TestTEControlVolume:
    """FACT 2: the TE aux row -- the ONLY equation that pins Gamma -- sits on
    an up/down asymmetric control volume."""

    @pytest.mark.parametrize("level", ["coarse", "medium"])
    def test_below_te_fan_is_not_cut(self, level):
        """The eps-shift must not manufacture cuts in the below-TE fan
        (Lopez p.57). Regression pin for the 2026-07-12 fix."""
        mesh = _load(level)
        _, cm = _setup(mesh)
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
        assert not is_cut[true_fan].any(), (
            "below-TE fan elements are being CUT -- they lie entirely at/below "
            "the sheet and only touch the TE vertex (Lopez section 3.5.4 p.57)"
        )
        assert set(true_fan) == set(cm.te_lower_elems)

    def test_te_control_volume_is_asymmetric(self):
        """THE STANDING B4 DEFECT: on a SYMMETRIC airfoil the TE upper and
        lower control volumes should match; they do not, because the eps
        shift sends every on-sheet node '+'."""
        mesh = _load("medium")
        _, cm = _setup(mesh)
        c = _te_fan_census(mesh, cm, cm.te_nodes[0])
        n_up, n_lo = len(c["upper"]), len(c["cut"]) + len(c["lower"])
        # Documented state (medium): 9 upper / 6 lower / 3 cut.
        assert len(c["upper"]) == 9 and len(c["lower"]) == 6
        assert len(c["cut"]) == 3
        assert n_up != len(c["lower"]), (
            "TE fan is symmetric now -- if B4 fixed it, update this pin"
        )


class TestB4Artifacts:
    """Headless diagnostic artifact steering the B4 re-derivation."""

    def test_export_te_diagnosis(self, gate_artifacts_dir):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        mesh = _load("medium")
        _, cm = _setup(mesh)
        el = np.asarray(mesh.elements, dtype=np.int64)
        cen = mesh.nodes[el].mean(axis=1)
        c = _te_fan_census(mesh, cm, cm.te_nodes[0])

        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13, 5.2))

        # (a) TE control volume: the element fan, coloured by class.
        fan_pts = cen[c["fan"]]
        pad = 0.35 * max(np.ptp(fan_pts[:, 0]), np.ptp(fan_pts[:, 1]))
        x0, x1 = fan_pts[:, 0].min() - pad, fan_pts[:, 0].max() + pad
        y0, y1 = fan_pts[:, 1].min() - pad, fan_pts[:, 1].max() + pad
        te = mesh.nodes[cm.te_nodes[0]]
        ax0.axhspan(y0, 0.0, color="tab:blue", alpha=0.05)
        ax0.plot([te[0], x1], [te[1], te[1]], "k--", lw=1.6, zorder=1,
                 label="wake sheet (y = 0)")
        for key, col, lab in (
            ("upper", "tab:red", "upper fan (main = $\\phi_u$)"),
            ("lower", "tab:blue", "below-TE fan (TE ref = aux = $\\phi_l$)"),
            ("cut", "tab:green", "cut (sheet through element)"),
        ):
            p = cen[c[key]]
            ax0.scatter(p[:, 0], p[:, 1], s=130, c=col, alpha=0.85, zorder=3,
                        edgecolors="k", linewidths=0.5,
                        label=f"{lab}  [{len(c[key])}]")
        ax0.plot(te[0], te[1], "k*", ms=20, zorder=4, label="TE node (doubled)")
        ax0.set_xlim(x0, x1)
        ax0.set_ylim(y0, y1)
        ax0.set_xlabel("x")
        ax0.set_ylabel("y")
        ax0.set_title("B4 defect: TE control volume is UP/DOWN ASYMMETRIC\n"
                      f"upper {len(c['upper'])} vs lower {len(c['lower'])} "
                      "on a SYMMETRIC airfoil (element centroids)")
        ax0.legend(fontsize=8, loc="lower left", framealpha=0.92)
        ax0.grid(alpha=0.3)

        # (b) the consequence: emergent Gamma converges to the WRONG value.
        lv = ["coarse", "medium", "fine"]
        g = [MEASURED_GAMMA[k] for k in lv]
        h = [1.0, 0.5, 0.25]
        ax1.plot(h, g, "o-", color="tab:red", lw=2, ms=9,
                 label="B-path implicit Kutta (emergent $\\Gamma$)")
        ax1.axhline(CONFORMING_GAMMA, color="tab:blue", ls="--", lw=2,
                    label=f"conforming / thin-airfoil $\\Gamma$ = {CONFORMING_GAMMA}")
        for hi, gi in zip(h, g):
            ax1.annotate(f"{gi:.4f}", (hi, gi), textcoords="offset points",
                         xytext=(0, 9), ha="center", fontsize=9)
        ax1.set_xlabel("relative mesh size h  (coarse / medium / fine)")
        ax1.set_ylabel("$\\Gamma$")
        ax1.set_title("Emergent circulation converges to ~0.168, NOT 0.120\n"
                      "(mesh-convergent ⇒ a METHOD defect, not discretization)")
        ax1.set_xscale("log")
        ax1.invert_xaxis()
        ax1.set_ylim(0.10, 0.22)
        ax1.legend(fontsize=9)
        ax1.grid(alpha=0.3, which="both")

        fig.tight_layout()
        png = gate_artifacts_dir / "b4_te_control_volume.png"
        fig.savefig(png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        assert png.exists()

        csv = gate_artifacts_dir / "summary.csv"
        with open(csv, "w") as f:
            f.write("quantity,value,note\n")
            f.write(f"te_fan_upper,{len(c['upper'])},plain '+' elements\n")
            f.write(f"te_fan_lower,{len(c['lower'])},below-TE fan (aux)\n")
            f.write(f"te_fan_cut,{len(c['cut'])},sheet through element\n")
            f.write(f"conforming_gamma,{CONFORMING_GAMMA},truth\n")
            for k in lv:
                err = abs(MEASURED_GAMMA[k] - CONFORMING_GAMMA) / CONFORMING_GAMMA
                f.write(f"gamma_{k},{MEASURED_GAMMA[k]},err {err*100:.1f}%\n")
        assert csv.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
