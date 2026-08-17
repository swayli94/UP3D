"""
Mesh provenance: a committed fingerprint for a file that is NOT committed.

★ WHY THIS EXISTS (GS4.0, 2026-08-16). Every 3-D production mesh in this
project is gitignored -- `onera_m6/*.msh`, `onera_m6_wingbody*/*.msh`,
`onera_m6_roundtip/*.msh`, `fuselage_bor/*.msh` -- because they are large and
regenerable. The generators ARE tracked, so in principle the mesh is derivable
from HEAD. In practice that failed, and the failure was measured:

  On 2026-08-04 `coarse.msh` and `medium.msh` were regenerated when the level
  names flipped flat -> round. **That geometry change was completely invisible
  in git.** The only surviving evidence was the .msh file's mtime. Every
  pre-08-04 anchor silently became a FLAT-CAP number while the tree looked
  unchanged, and a 2.0 % discrepancy that this produced stood for five days and
  was restated as "unexplained" in EIGHT documents before someone spent the 14
  seconds of compute needed to pin it down (phase 3 round 36: F-MESH).

The project already made "are these two numbers the same thing?" a standing
check, with cross-pipeline / cross-level / cross-provenance / cross-time /
cross-MESH-FAMILY as its five known members. The first four can be answered
from tracked artifacts. The fifth could not, because the mesh was not
identifiable. This module makes it identifiable:

  * `mesh_fingerprint(path)` -- content hash plus the counts that distinguish
    mesh FAMILIES, computed from the .msh alone. Cheap, pure, no solve.
  * `write_manifest(path, **extra)` -- drops `<name>.manifest.json` beside the
    mesh. `write_mesh` calls it, so all 10 generators that route through the
    project's writer get provenance for free.
  * `read_manifest(path)` -- for the consuming side: a script that commits a
    CSV can stamp the mesh hash into a column, and then "same mesh?" is a
    string comparison instead of an archaeology exercise.

★ The manifest is a SIDECAR and never changes the .msh bytes. It is also
excluded from its own hash by construction -- the hash is over the mesh file,
which the manifest is not part of. (The "does my guard measure itself?"
question, asked because three guards this season did exactly that.)

★ Deliberately NOT in the fingerprint: wall-clock, mtime, hostname, or the
generator's git SHA. Those make two byte-identical meshes compare unequal,
which is the opposite of the point. Provenance that is *not* reproducible goes
in `extra` (where a caller can put whatever it likes) and never in `sha256`.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict

import numpy as np

#: read in chunks so a 1.16 M-tet mesh does not need to be resident twice.
_CHUNK = 1 << 20

MANIFEST_SUFFIX = ".manifest.json"

#: bumped only if the meaning of a field changes; a consumer that finds an
#: unfamiliar version must say so rather than silently compare.
SCHEMA_VERSION = 1


def manifest_path(mesh_file: Path | str) -> Path:
    """`.../coarse.msh` -> `.../coarse.msh.manifest.json`.

    The full mesh filename is kept (not the stem) so `coarse.msh` and
    `coarse_flat.msh` cannot collide -- those two are exactly the pair whose
    confusion this module exists to prevent.
    """
    return Path(str(mesh_file) + MANIFEST_SUFFIX)


def _file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def mesh_fingerprint(mesh_file: Path | str) -> Dict:
    """Identify a mesh file: content hash plus the counts that separate families.

    The counts are not redundant with the hash. A hash answers "is this the same
    FILE?"; the counts answer "how does it differ?", which is what a reader of a
    committed CSV actually needs. Measured on the pair that motivated this
    module: `coarse.msh` (round tip) and `coarse_flat.msh` differ by 68 624 vs
    55 531 tets, and the wall triangles differ almost entirely at the tip
    (2013 vs 280) while the non-tip count agrees to 0.10 %.

    Args:
        mesh_file: path to a `.msh` file.

    Returns:
        dict with `sha256`, `bytes`, `n_nodes`, `n_tets`, and `n_tris`
        (per boundary group, sorted) -- all derived from the file only, so two
        byte-identical meshes always fingerprint equal regardless of when,
        where or by whom they were written.

    Raises:
        FileNotFoundError: the mesh does not exist. Deliberately loud: a
            fingerprint that quietly returns "unknown" would put this module in
            the same class as the defect it was written to prevent.
    """
    path = Path(mesh_file)
    if not path.exists():
        raise FileNotFoundError(f"cannot fingerprint a mesh that is not there: {path}")

    from pyfp3d.mesh.reader import read_mesh

    mesh = read_mesh(path)
    tris = {name: int(len(faces))
            for name, faces in sorted(mesh.boundary_faces.items())}
    return {
        "schema": SCHEMA_VERSION,
        "sha256": _file_sha256(path),
        "bytes": int(path.stat().st_size),
        "n_nodes": int(len(mesh.nodes)),
        "n_tets": int(len(mesh.elements)),
        "n_tris": tris,
    }


def write_manifest(mesh_file: Path | str, **extra) -> Path:
    """Write `<mesh>.manifest.json` beside the mesh; return its path.

    `extra` is merged in under the key `extra` (never at top level, so a caller
    cannot shadow a fingerprint field). Put non-reproducible provenance there --
    generator name, level, parameters, git SHA -- and note that none of it
    participates in `sha256`.
    """
    fp = mesh_fingerprint(mesh_file)
    if extra:
        fp["extra"] = {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                       for k, v in extra.items()}
    out = manifest_path(mesh_file)
    out.write_text(json.dumps(fp, indent=2, sort_keys=True) + "\n")
    return out


def read_manifest(mesh_file: Path | str) -> Dict | None:
    """Read the sidecar, or None if it is not there.

    None rather than an exception, because a missing manifest is the normal
    state for a mesh generated before GS4.0 -- the caller decides whether that
    is fatal. `mesh_fingerprint` is the loud one.
    """
    p = manifest_path(mesh_file)
    if not p.exists():
        return None
    return json.loads(p.read_text())


def describe_difference(a: Dict, b: Dict) -> str:
    """One human-readable line saying whether two fingerprints are the same mesh.

    Exists so that "are these two numbers from the same mesh?" is answerable in
    a print statement rather than by reading two JSON files side by side.
    """
    if a.get("sha256") == b.get("sha256"):
        return "SAME MESH (sha256 identical)"
    bits = [f"sha256 {a.get('sha256', '?')[:12]} vs {b.get('sha256', '?')[:12]}"]
    for k in ("n_nodes", "n_tets", "bytes"):
        if a.get(k) != b.get(k):
            bits.append(f"{k} {a.get(k)} vs {b.get(k)}")
    ta, tb = a.get("n_tris", {}) or {}, b.get("n_tris", {}) or {}
    for name in sorted(set(ta) | set(tb)):
        if ta.get(name) != tb.get(name):
            bits.append(f"n_tris[{name}] {ta.get(name)} vs {tb.get(name)}")
    return "DIFFERENT MESH -- " + "; ".join(bits)
