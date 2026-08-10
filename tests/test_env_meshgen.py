"""One test that says out loud whether mesh generation works in this environment.

Written because of a measured failure mode, not as boilerplate
(docs/dev_phase_two/20260801-0230-gmsh-broken.md). On 2026-08-01 `import gmsh` on
this machine raised

    OSError: libGLU.so.1: cannot open shared object file

and the repo's only symptoms were SIXTEEN silently skipped M1 tests plus a bare
OSError from any generator script. Both read as "fine". Meanwhile the committed
workflow (CLAUDE.md, discipline #7) says the M6 .msh files are gitignored and
regenerate in ~30 s -- so the state was "nobody can currently reproduce the meshes
on disk", which is exactly what the repo's evidence-must-be-regenerable rule exists
to prevent.

This test does not FAIL when gmsh is unavailable -- the skip-if-absent convention
stays, since a working gmsh is not required to run the solver on committed meshes.
What it does is make the skip diagnostic: it names the missing library and the
package that provides it, in one place, instead of leaving the reader to infer it
from a skip count.
"""

import pytest

#: measured dependency chain of the pypi gmsh wheel (4.15.2) on Ubuntu jammy --
#: found by adding one library at a time until gmsh.initialize() succeeded
_PROVIDERS = {
    "libGLU.so.1": "libglu1-mesa",
    "libOpenGL.so.0": "libopengl0",
    "libGLdispatch.so.0": "libglvnd0",
    "libGLX.so.0": "libglx0",
}


def meshgen_status():
    """(ok, detail). `detail` names the missing library and its package when it can."""
    try:
        import gmsh
    except Exception as exc:                                      # noqa: BLE001
        msg = str(exc)
        for lib, pkg in _PROVIDERS.items():
            if lib in msg:
                return False, (
                    f"gmsh cannot load: {lib} is missing (provided by {pkg}). "
                    f"The full chain this wheel needs is: "
                    f"{', '.join(f'{k} <- {v}' for k, v in _PROVIDERS.items())}. "
                    f"Fix: sudo apt install {' '.join(_PROVIDERS.values())} -- or, "
                    f"without sudo, apt-get download those packages, dpkg -x them "
                    f"into a directory and point LD_LIBRARY_PATH at its "
                    f"usr/lib/x86_64-linux-gnu (see "
                    f"docs/dev_phase_two/20260801-0230-gmsh-broken.md)")
        return False, f"gmsh import failed: {type(exc).__name__}: {msg}"
    try:
        gmsh.initialize()
    except Exception as exc:                                      # noqa: BLE001
        return False, f"gmsh imported but initialize() failed: {exc}"
    try:
        version = gmsh.option.getString("General.Version")
    finally:
        gmsh.finalize()
    return True, f"gmsh {version}"


def test_meshgen_available():
    """Mesh generation is available -- or the skip reason says exactly why not.

    Deliberately a skip and not a failure: the solver, its tests and the committed
    evidence all run against meshes already on disk, so a broken gmsh must not turn
    the suite red. It must not be INVISIBLE either, which is what this fixes.
    """
    ok, detail = meshgen_status()
    if not ok:
        pytest.skip(f"MESH GENERATION UNAVAILABLE -- {detail}")
    assert "." in detail          # a version string came back
