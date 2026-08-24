"""
Is the optimal recovery patch scale a LENGTH that tracks the mesh, or a
TOPOLOGICAL property pinned at the immediate neighbours?

Pre-registration: phases/p4/docs/dev_phase_four/20260817-1500-level-ladder-prereg.md
Previous round:   phases/p4/docs/dev_phase_four/20260817-1100-coherence-ladder-verdict.md (C-PEAK)

Three levels spanning 20x in element count, with the condition, geometry and
recipe held fixed and only h varying. K-FIXED (k*=1 everywhere) makes the result
transferable; K-DRIFT makes it a length scale; K-VANISH kills the candidate
because the previous round's finding would then be mesh-specific.

★ Genuinely zero-solve: every phi is read from a cache. The operators are
IMPORTED from run_gs40c_coherence_ladder rather than re-typed, so the two rounds
cannot drift apart.

★ FOUR FORBIDDEN SENTENCES: necessary-not-sufficient; zero-solve says nothing
about behaviour inside the solver; COHERENT IS NOT CORRECT; and k*=1 holding
across three LEVELS does not mean it holds across CONDITIONS -- all three sit at
M0.88/alpha 3.06.
"""

import csv
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "bench"))

import run_gs40c_coherence_ladder as L                              # noqa: E402
from pyfp3d.mesh.manifest import mesh_fingerprint, read_manifest    # noqa: E402
from pyfp3d.mesh.metrics import (build_face_adjacency,              # noqa: E402
                                 compute_edge_lengths,
                                 precompute_element_geometry)
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2             # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared           # noqa: E402

MDIR = REPO / "cases/meshes/onera_m6_wingbody_conforming"
OUT = REPO / "bench/gate_results/gs40e_level_ladder.csv"
M_INF, GAMMA = 0.88, 1.4
ARMS = ("P0", "P1", "P2", "P3", "N1", "S1")
LADDER = ("P0", "P1", "P2", "P3")

#: registered, and every one CARRIED OVER VERBATIM from the previous round
PEAK_FACTOR, FLAT_RATIO = L.PEAK_FACTOR, L.FLAT_RATIO
COMP_MIN, USABLE_MIN, MIN_CROSS = L.COMP_MIN, L.USABLE_MIN, L.MIN_CROSS
GATE_S = 20 * 60

STATES = (("xcoarse", REPO / "bench/gate_results/gs40d_levels/xcoarse.npz", "binding"),
          ("coarse", REPO / "bench/gate_results/gs40d_levels/coarse.npz", "binding"),
          ("medium", REPO / "bench/gate_results/le14_cache/rc0.05.npz", "anchor"))

#: G-REPRO: medium must reproduce the previous round's published C(P0..P3)
REPRO_C = dict(P0=0.191, P1=0.428, P2=0.284, P3=0.214)


def guard_recipe():
    """G-RECIPE: SOURCE-COMPARE the probe's recipe against run_le14_common_root.

    Not "I remember copying it". The two files must name the same constants and
    pass the same literals to the ramp; a drifted recipe would make the levels
    incomparable, which is the 5th question in its cross-recipe form.
    """
    probe = (REPO / "bench/run_gs40d_level_cost_probe.py").read_text()
    le14 = (REPO / "bench/run_le14_common_root.py").read_text()
    need_names = ("cap.CONF_SEED_KW", "cap.CONF_RAMP_NK", "cap.WB_MSTART", "cap.DM")
    need_lits = ("dm_min=0.01", "freeze_tol=1e-5", "intermediate_tol=1e-4",
                 'kutta_estimator="pressure"', "n_picard_seed=0")
    missing = [t for t in need_names + need_lits
               if t not in re.sub(r"\s+", "", probe).replace(",", ",")
               and t.replace(" ", "") not in re.sub(r"\s+", "", probe)]
    missing += [t + " [le14]" for t in need_names + need_lits
                if t.replace(" ", "") not in re.sub(r"\s+", "", le14)]
    assert not missing, f"G-RECIPE: recipe fragments not found: {missing}"
    print(f"G-RECIPE  {len(need_names + need_lits)} recipe fragments present in "
          "BOTH the probe and run_le14_common_root  PASS")


def load(tag, npz):
    d = np.load(npz)
    assert bool(d["conv"]), f"{tag}: not converged"
    assert int(d["nlim"]) == 0 and int(d["nflr"]) == 0, f"{tag}: clamped"
    sha = str(d["mesh_sha"]) if "mesh_sha" in d else None
    return np.asarray(d["phi"], dtype=np.float64), sha


def main():
    t0 = time.perf_counter()
    L.guard_zero_solve()
    print("G-Z       no solver driver imported; every phi comes from a cache  PASS")
    guard_recipe()

    rows, C_all, extra = [], {}, {}
    for tag, npz, role in STATES:
        msh = MDIR / f"{tag}.msh"
        man = read_manifest(msh)
        sha = mesh_fingerprint(msh)["sha256"]
        assert man is not None and sha == man["sha256"], f"G-MESH {tag}: mesh moved"
        phi, npz_sha = load(tag, npz)
        if npz_sha:
            assert npz_sha == sha, (
                f"G-MESH {tag}: the cached state was produced on a DIFFERENT mesh "
                f"({npz_sha[:12]} vs {sha[:12]})")
        print(f"\nG-MESH {tag:8} sha {sha[:12]}  "
              f"{'(npz agrees)' if npz_sha else '(npz predates mesh_sha; manifest only)'}  PASS")

        mc, _ = cut_wake(read_mesh(msh))
        nodes, elements = mc.nodes, mc.elements
        assert len(phi) == len(nodes), f"{tag}: phi/mesh mismatch"
        B, _ = precompute_element_geometry(nodes, elements)
        face_nb, _ = build_face_adjacency(elements)
        fa, fb = L.faces(face_nb)
        h_med = float(np.median(compute_edge_lengths(nodes, elements)))

        P, p1_elems = {}, None
        for k in (0, 1, 2, 3):
            po, pi = L.khop_patches(elements, face_nb, k)
            P[f"P{k}"] = L.patch_nodes(elements, po, pi)
            if k == 1:
                p1_elems = (po, pi)
        P["N1"] = L.node_ring_patch(elements, len(nodes))

        gA = np.empty((len(elements), 3)); q2 = np.empty(len(elements))
        element_velocity_q2(elements, B, phi, gA, q2)
        g, conds = {}, []
        for name in ("P0", "P1", "P2", "P3", "N1"):
            g[name] = L.refit(nodes, elements, phi, *P[name],
                              cond_out=conds if name == "P1" else None)
        g["S1"] = L.average_over_patch(nodes, elements, gA, *p1_elems)

        e_a = float(np.max(np.abs(g["P0"] - gA)))
        print(f"G-A    {tag:8} P0 vs kernel {e_a:.2e}  "
              f"{'PASS' if e_a < 1e-10 else '** FAIL **'}")
        assert e_a < 1e-10, "G-A: kill criterion 3"
        a = np.random.default_rng(0).normal(size=3)
        bar = max(1e-12, 10 * conds[0] * np.finfo(float).eps)
        worst = max(float(np.max(np.abs(L.refit(nodes, elements, nodes @ a, *P[n]) - a)))
                    for n in ("P1", "P2", "P3", "N1"))
        print(f"G-CONS {tag:8} worst {worst:.2e}  bar {bar:.2e} "
              f"(kappa {conds[0]:.2e}, measured)  "
              f"{'PASS' if worst < bar else '** FAIL **'}")
        assert worst < bar, "G-CONS: kill criterion 3"

        C, res = {}, {}
        for name in ARMS:
            m2 = mach_number_squared((g[name] * g[name]).sum(axis=1), M_INF, GAMMA)
            r = L.coherence(np.sqrt(np.maximum(m2, 0.0)), fa, fb)
            res[name] = r; C[name] = r["C"]
            #: physical patch radius: median |node - centroid| over the patch
            if name in P:
                cent = nodes[elements].mean(axis=1)
                off, idx = P[name]
                sample = range(0, len(elements), max(1, len(elements) // 4000))
                rad = float(np.median([np.median(np.linalg.norm(
                    nodes[idx[off[e]:off[e + 1]]] - cent[e], axis=1)) for e in sample]))
            else:
                rad = float("nan")
            rows.append(dict(level=tag, role=role, n_tets=len(elements),
                             h_median=round(h_med, 6), arm=name,
                             patch_radius=round(rad, 6),
                             radius_over_h=round(rad / h_med, 3),
                             n_cross=r["n_cross"], n_comp=r["n_comp"],
                             max_comp=r["max_comp"], n_ge50=r["n_usable"],
                             C=round(r["C"], 4)))
        C_all[tag] = C
        extra[tag] = dict(res=res, h=h_med, rows={r["arm"]: r for r in rows if r["level"] == tag})

        print(f"[{tag}]  h_med {h_med:.4f}   " + "  ".join(
            f"{n} {C[n]:.3f}" for n in ARMS))
        print("        n_cross " + "  ".join(f"{n} {res[n]['n_cross']}" for n in LADDER))

        if tag == "medium":
            bad = {k: (C[k], REPRO_C[k]) for k in LADDER
                   if abs(C[k] - REPRO_C[k]) > 1e-3}
            print(f"G-REPRO medium vs the published C(P0..P3): "
                  f"{'PASS' if not bad else '** FAIL ** ' + str(bad)}")
            assert not bad, "G-REPRO: kill criterion 2 -- instrument, not a finding"

    # ---------------------------------------------------------- verdict -----
    print("\n" + "=" * 76)
    per = {}
    for tag, _n, role in STATES:
        C = C_all[tag]; res = extra[tag]["res"]
        thin = [k for k in LADDER if res[k]["n_cross"] < MIN_CROSS]
        vals = [C[k] for k in LADDER]
        if thin or any(np.isnan(v) for v in vals):
            per[tag] = ("UNDEFINED", f"thin sonic set at {thin}")
        elif max(vals) / max(min(vals), 1e-30) < FLAT_RATIO:
            per[tag] = ("NO-PEAK", f"flat, max/min {max(vals)/min(vals):.2f}")
        elif max(C["P1"], C["P2"]) >= PEAK_FACTOR * max(C["P0"], C["P3"]):
            #: ties resolve to the SMALLEST k -- fixed in advance
            kstar = min(LADDER, key=lambda k: (-round(C[k], 6), LADDER.index(k)))
            per[tag] = (f"PEAK k*={kstar}",
                        f"{max(C['P1'], C['P2']) / max(C['P0'], C['P3']):.2f}x endpoints")
        else:
            per[tag] = ("NO-PEAK", "no interior maximum")
        print(f"{tag:8} ({role:8}) C = " + " ".join(f"{C[k]:.3f}" for k in LADDER)
              + f"   -> {per[tag][0]}   [{per[tag][1]}]")

    binding = [t for t, _n, r in STATES if r == "binding"]
    codes = [per[t][0] for t in binding]
    if all(c == "UNDEFINED" for c in codes):
        verdict = "UNDEFINED (both binding legs thin) -- kill criterion 4"
    elif any(c.startswith("NO-PEAK") for c in codes):
        verdict = "K-VANISH (a binding level has no interior maximum) -- candidate dies"
    elif any(c == "UNDEFINED" for c in codes):
        verdict = ("UNDEFINED (binding reduced to one leg; the registration "
                   "forbids promoting medium)")
    elif all(c == "PEAK k*=P1" for c in codes):
        verdict = "K-FIXED (k*=1 on every binding level -- topological, transferable)"
    else:
        verdict = "K-DRIFT (interior maximum everywhere, but k* moves with the level)"
    print("=" * 76)
    print(f"VERDICT: {verdict}")

    print("\nRECORDED (registered derived quantities):")
    for tag, _n, _r in STATES:
        C = C_all[tag]; res = extra[tag]["res"]
        ks = per[tag][0].split("k*=")[-1] if "k*=" in per[tag][0] else "P1"
        rr = extra[tag]["rows"].get(ks, {})
        gap = abs(C["P1"] - C["S1"]) / max(C["S1"], 1e-30)
        print(f"  {tag:8} h {extra[tag]['h']:.4f}  radius/h at {ks} = "
              f"{rr.get('radius_over_h')}   D-USABLE(>={USABLE_MIN}) = "
              f"{res[ks]['n_usable']}   |C_P1-C_S1|/C_S1 = {gap:.3f}")
    print("\n★ COHERENT IS NOT CORRECT, and k* fixed across LEVELS is not k* fixed\n"
          "  across CONDITIONS -- all three sit at M0.88 / alpha 3.06.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    OUT.with_suffix(".json").write_text(json.dumps(
        dict(verdict=verdict, per_level={k: v[0] for k, v in per.items()},
             wall_s=round(time.perf_counter() - t0, 1)), indent=2) + "\n")
    dt = time.perf_counter() - t0
    print(f"\nwrote {OUT.name} + .json;  {dt:.1f} s (gate {GATE_S} s)")
    assert dt < GATE_S, "kill criterion 5"


if __name__ == "__main__":
    main()
