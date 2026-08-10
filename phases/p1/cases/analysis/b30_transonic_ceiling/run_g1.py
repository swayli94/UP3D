"""GB30.1 -- baseline anchor verification (B30 PRE_REGISTRATION.md section 3).

Cache-only, no solves: rebuild the two production operators on the medium
meshes, verify the committed B29 demo caches (a) load, (b) carry the
committed anchor scalars (checks.csv / B29 ledger), and (c) are
shape-compatible warm-start seeds for GB30.2.

PASS = every quoted anchor reproduced at quoted precision AND both seeds
load with shapes matching the freshly built operators. FAIL = cache vs
committed CSV drift (pre-reg risk T1: record as an independent finding,
re-anchor, phase continues).

Artifacts: results/g1_anchor.csv
"""

import csv
import sys

import numpy as np

from wb30 import (ANCHORS, DEMO_OUT, OUT, build_conf, build_ls_flat,
                  load_mesh, LS_MESH_DIR)


def check_row(name, path, expect, n_dof):
    d = np.load(path, allow_pickle=True)
    row = {"cache": name, "exists": path.exists()}
    phi_key = "phi_ext" if "phi_ext" in d else "phi"
    phi = d[phi_key]
    row["phi_len"] = int(phi.shape[0])
    row["dof_match"] = bool(phi.shape[0] == n_dof)
    ok = row["dof_match"]
    for key, want in expect.items():
        got = d[key].item() if key in d else None
        row[f"got_{key}"] = got
        if isinstance(want, bool):
            match = bool(got) == want
        elif isinstance(want, str):
            match = str(got) == want
        else:
            match = got is not None and abs(float(got) - want) < 5e-5
        row[f"ok_{key}"] = bool(match)
        ok &= match
    row["PASS"] = bool(ok)
    print(f"  [{name}] phi_len={row['phi_len']} dof_match={row['dof_match']} "
          + " ".join(f"{k}={row[f'got_{k}']}(ok={row[f'ok_{k}']})"
                     for k in expect)
          + f" -> {'PASS' if row['PASS'] else 'FAIL'}", flush=True)
    return row


def main():
    rows = []

    print("=== GB30.1: building the LS production operator (medium) ===",
          flush=True)
    mesh_ls = load_mesh(LS_MESH_DIR / "medium.msh")
    _wls, cm, _mvop = build_ls_flat(mesh_ls)
    n_ext = mesh_ls.nodes.shape[0] + int((cm.ext_dof_of_node >= 0).sum())
    rows.append(check_row("ls_flat_medium_084",
                          DEMO_OUT / "ls_flat_medium_084.npz",
                          ANCHORS["ls_flat_medium_084"], n_ext))
    rows.append(check_row("ls_flat_coarse_084",
                          DEMO_OUT / "ls_flat_coarse_084.npz",
                          ANCHORS["ls_flat_coarse_084"],
                          # coarse seed shape not needed for G2; report only
                          n_dof=int(np.load(DEMO_OUT / "ls_flat_coarse_084.npz")
                                    ["phi_ext"].shape[0])))

    print("=== GB30.1: loading the conforming medium mesh ===", flush=True)
    mc, _wc = build_conf("medium")
    n_cut = mc.nodes.shape[0]
    rows.append(check_row("conf_medium_079",
                          DEMO_OUT / "conf_medium_079.npz",
                          ANCHORS["conf_medium_079"], n_cut))
    rows.append(check_row("conf_coarse_084",
                          DEMO_OUT / "conf_coarse_084.npz",
                          ANCHORS["conf_coarse_084"],
                          n_dof=int(np.load(DEMO_OUT / "conf_coarse_084.npz")
                                    ["phi"].shape[0])))

    keys = sorted({k for r in rows for k in r})
    with open(OUT / "g1_anchor.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT / 'g1_anchor.csv'}")
    n_pass = sum(r["PASS"] for r in rows)
    print(f"GB30.1: {n_pass}/{len(rows)} anchor rows PASS")
    sys.exit(0 if n_pass == len(rows) else 1)


if __name__ == "__main__":
    main()
