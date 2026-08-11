"""G0 — is the structured generator genuinely SINGLE-VARIABLE? (phase 3 task 2, route A)

Pre-registered in docs/dev_phase_three/20260811-0300-hex-mesh-prereg.md §3: this check comes
BEFORE every reading, because phase two learned the hard way that a conclusion of the form
"refining X controls Y" cannot be earned from knobs that move other things. Four knobs of the
unstructured generator were caught out of scope -- `h_wall` (+41.7 % LE face count), `h_edge`
(sizes the LE *and* the TE), `h_far` (+24.04 % near-body median cell size) and even `r_far`
(+67.33 %, and it is not a size field at all) -- and two published conclusions had to be
rewritten because the guard only checked the surface.

★ So the criterion here is BIT-IDENTITY, not "small change": move one knob and the quantities
that knob has no business touching must be unchanged to the last bit. That is only askable
because `structured.py` builds the grading OUTWARD FROM THE WALL with a fixed growth rate and
DERIVES the layer count -- growing the domain appends layers and leaves earlier nodes alone.

No solves: mesh construction only.

Outputs (TRACKED): bench/gate_results/hex_g0_single_var.csv
"""

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)

from pyfp3d.meshgen.structured import cylinder_o_grid_2d       # noqa: E402

CSV = os.path.join(HERE, "gate_results", "hex_g0_single_var.csv")
BASE = dict(radius=1.0, r_far=20.0, n_theta=96, h_wall_normal=0.02, growth=1.15)

#: (name, the knob moved, what MUST stay bit-identical, what is ALLOWED to move)
ARMS = (
    ("r_far_x3", dict(r_far=60.0),
     "near-body nodes, wall edges, tangential spacing, radial sequence prefix",
     "radial layer count, node count"),
    ("h_wall_normal_half", dict(h_wall_normal=0.01),
     "wall nodes, wall edges, tangential spacing",
     "radial layer count, radial sequence"),
    ("n_theta_x2", dict(n_theta=192),
     "radial layer count, radial sequence, achieved outer radius",
     "wall node count, tangential spacing"),
    ("growth_softer", dict(growth=1.08),
     "wall nodes, wall edges, tangential spacing, FIRST radial spacing",
     "radial layer count, outer radii"),
)


def radii_of(pts, n_theta):
    """Radial coordinate of node column i = 0, i.e. one node per layer."""
    return np.hypot(pts[::n_theta, 0], pts[::n_theta, 1])


def main():
    p0, t0, e0, i0 = cylinder_o_grid_2d(**BASE)
    r0 = radii_of(p0, i0["n_theta"])
    nt0 = i0["n_theta"]
    print("G0: single-variable check on the structured O-grid (bit-identity, not 'small')\n")
    print(f"  baseline: n_theta {nt0}, n_radial {i0['n_radial']}, "
          f"h_t {i0['h_wall_tangential']:.8f}, r_achieved {i0['r_achieved']:.6f}\n")
    rows, all_ok = [], True
    for name, override, invariant, allowed in ARMS:
        p, t, e, i = cylinder_o_grid_2d(**dict(BASE, **override))
        nt = i["n_theta"]
        r = radii_of(p, nt)
        checks = {}
        #: wall = the first n_theta nodes by construction (layer j = 0)
        if nt == nt0:
            checks["wall_nodes_bit_identical"] = bool(np.array_equal(p0[:nt0], p[:nt0]))
            checks["wall_edges_bit_identical"] = bool(np.array_equal(e0["wall"], e["wall"]))
            checks["h_t_bit_identical"] = bool(
                i0["h_wall_tangential"] == i["h_wall_tangential"])
        if "n_theta" in override:
            #: the radial construction must not notice a tangential change
            checks["radial_sequence_bit_identical"] = bool(np.array_equal(r0, r))
            checks["r_achieved_bit_identical"] = bool(i0["r_achieved"] == i["r_achieved"])
        if "r_far" in override:
            #: ★ the decisive one: growing the domain must only APPEND layers
            k = min(len(r0), len(r))
            checks["radial_prefix_bit_identical"] = bool(np.array_equal(r0[:k], r[:k]))
            n_common = min(len(p0), len(p))
            checks["node_prefix_bit_identical"] = bool(
                np.array_equal(p0[:n_common], p[:n_common]))
        if "growth" in override:
            checks["first_radial_spacing_bit_identical"] = bool(
                (r0[1] - r0[0]) == (r[1] - r[0]))
        ok = all(checks.values())
        all_ok &= ok
        rows.append(dict(arm=name, knob=";".join(f"{k}={v}" for k, v in override.items()),
                         must_not_move=invariant, allowed_to_move=allowed,
                         n_radial=i["n_radial"], n_nodes=len(p),
                         h_t=round(i["h_wall_tangential"], 10),
                         r_achieved=round(i["r_achieved"], 8),
                         **{k: v for k, v in checks.items()},
                         verdict="PASS" if ok else "FAIL"))
        print(f"  {name:20} {'PASS' if ok else '★ FAIL'}  "
              f"n_radial {i['n_radial']:>3}  nodes {len(p):>6}")
        for k, v in checks.items():
            print(f"      {'✓' if v else '✗'} {k}")

    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}),
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")
    print("\n=== reading (criterion fixed in the pre-registration, §3 G0) ===")
    if all_ok:
        print("  G0 PASS — every knob is single-variable to the BIT. This is the property the")
        print("  unstructured generator could not provide (four knobs measured out of scope),")
        print("  and it is a property of the CONSTRUCTION, not of a lucky knob choice.")
        print("  ⇒ readings taken on this family may attribute an effect to the knob moved.")
    else:
        print("  ★ G0 FAIL — per the pre-registration, NO reading on this family counts until")
        print("  the generator is fixed; and if single-variable behaviour cannot be reached in")
        print("  two rounds, kill criterion 2 fires and route (B) becomes the candidate.")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
