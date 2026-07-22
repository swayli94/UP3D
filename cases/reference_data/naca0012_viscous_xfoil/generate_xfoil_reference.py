"""
Reference-data generator: viscous NACA0012 by XFOIL 6.99 (external tool,
NOT a repo dependency -- the binary is gitignored under tools/xfoil/).

Provenance / method (see README.md):
  - Code: XFOIL 6.99, built from the MIT source tarball with gfortran
    (double precision, -std=legacy -fallow-argument-mismatch; src/xfoil.f
    patched IDEV=2 so plots go to PostScript instead of X11 -> batch runs
    are headless). Build recipe: README.md in this directory.
  - Geometry: NACA 0012 via XFOIL's "NACA 0012" command (245 buffer
    points), repaneled to 280 panel nodes (PPAR -> N 280; note PPAR needs
    TWO blank lines to return: the first applies the repanel, the second
    exits the menu).
  - Conditions (pinned): M_inf = 0.5, Re = 3.0e6 (chord-based),
    alpha = 2 deg, viscous mode, Ncrit = 9 default.
  - Two forced-transition variants (OPER -> VPAR -> XTR xt xb):
      xtr005: trip at x/c = 0.05 on both surfaces (mostly turbulent);
      xtr030: trip at x/c = 0.30 on both surfaces (upper surface actually
              transitions naturally at x/c = 0.2668 ahead of the trip --
              XFOIL e^N wins upstream of XTR; lower is tripped at 0.30).
  - Per run: converged Newton solution (ITER 200; converges in ~5
    iterations to rms ~1e-5), DUMP of the viscous BL profiles, and a
    one-point polar save file (PACC) for cl/cd/cm. The stdout convergence
    printout is parsed as an independent cross-check of the polar values.

XFOIL DUMP file layout (parsed here): one header line, then surface rows
with 12 fields (s, x, y, Ue/Vinf, Dstar, Theta, Cf, H, H*, P, m, K)
ordered upper TE -> LE -> lower TE, then wake rows with 8 fields
(discarded). XFOIL's chord is 1, so Dstar is already dstar/c. The upper/
lower split is taken at the minimum-x station.

Writes: delta_star_cf_alpha2_m05_xtr005.csv,
        delta_star_cf_alpha2_m05_xtr030.csv,
        polar_summary.csv.
Run:    python generate_xfoil_reference.py
"""

import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]

# Pinned run conditions (do not change; see README.md).
MACH = 0.5
REYNOLDS = 3.0e6
ALPHA_DEG = 2.0
N_PANELS = 280
N_ITER = 200
CASES = {
    "xtr005": (0.05, 0.05),
    "xtr030": (0.30, 0.30),
}

BUILD_HINT = """\
XFOIL binary not found under tools/xfoil/ (it is gitignored -- build it
from source once):

    mkdir -p tools/xfoil && cd tools/xfoil
    wget https://web.mit.edu/drela/Public/web/xfoil/xfoil6.99.tgz
    tar xzf xfoil6.99.tgz
    cd Xfoil/plotlib
    cp config.make.gfortranDP config.make
    # add '-std=legacy -fallow-argument-mismatch' to FFLAGS in config.make
    make                                   # builds libPlt_gDP.a
    cd ../bin && cp Makefile_gfortran Makefile
    # add '-std=legacy -fallow-argument-mismatch' to FFLAGS and FFLOPT
    # patch ../src/xfoil.f: 'IDEV = 1' -> 'IDEV = 2' (headless batch)
    make xfoil BINDIR=.                    # binary: Xfoil/bin/xfoil

See cases/reference_data/naca0012_viscous_xfoil/README.md for details.
"""


def find_xfoil() -> Path:
    """Locate an executable xfoil binary below tools/xfoil/."""
    tools = REPO_ROOT / "tools" / "xfoil"
    candidates = sorted(tools.rglob("xfoil")) if tools.is_dir() else []
    for cand in candidates:
        if cand.is_file() and not cand.suffix and os.access(cand, os.X_OK):
            return cand
    print(BUILD_HINT, file=sys.stderr)
    sys.exit(1)


def batch_input(tag: str, xtr_top: float, xtr_bot: float,
                polar_name: str, dump_name: str) -> str:
    """XFOIL 6.99 stdin script for one operating point.

    Menu quirks that matter: PPAR needs two blank lines to exit after a
    change; PACC ignores inline args and prompts for the polar save file
    then the polar dump file (blank = none).
    """
    return (
        "NACA 0012\n"
        "PPAR\n"
        f"N {N_PANELS}\n"
        "\n\n"
        "OPER\n"
        f"VISC {REYNOLDS:.2E}\n"
        f"MACH {MACH:.1f}\n"
        "VPAR\n"
        f"XTR {xtr_top:.2f} {xtr_bot:.2f}\n"
        "SHOW\n"
        "\n"
        f"ITER {N_ITER}\n"
        "PACC\n"
        f"{polar_name}\n"
        "\n"
        f"ALFA {ALPHA_DEG:.1f}\n"
        f"DUMP {dump_name}\n"
        "\n"
        "QUIT\n"
    )


def run_case(binary: Path, workdir: Path, tag: str,
             xtr_top: float, xtr_bot: float) -> dict:
    polar_name = f"polar_{tag}.txt"
    dump_name = f"dump_{tag}.txt"
    script = batch_input(tag, xtr_top, xtr_bot, polar_name, dump_name)
    (workdir / f"{tag}.in").write_text(script)
    proc = subprocess.run(
        [str(binary)], input=script, cwd=workdir,
        capture_output=True, text=True, timeout=600,
    )
    log = proc.stdout + proc.stderr
    (workdir / f"{tag}.log").write_text(log)
    if proc.returncode != 0:
        raise RuntimeError(f"xfoil run {tag} failed (rc={proc.returncode}); "
                           f"see {workdir}/{tag}.log")

    mver = re.search(r"XFOIL Version\s+(\S+)", log)
    version = mver.group(1) if mver else "unknown"

    # --- cl/cd/cm from the polar save file (last numeric data row) ---
    polar_path = workdir / polar_name
    data_rows = []
    for line in polar_path.read_text().splitlines():
        tok = line.split()
        if len(tok) >= 7:
            try:
                data_rows.append([float(v) for v in tok])
            except ValueError:
                pass
    if not data_rows:
        raise RuntimeError(f"no polar data row found in {polar_path}")
    alpha, cl, cd, _cdp, cm = data_rows[-1][:5]

    # --- cross-check against the stdout convergence printout ---
    mcl = re.findall(r"a =\s*([-\d.]+)\s+CL =\s*([-\d.]+)", log)
    mcd = re.findall(r"Cm =\s*([-\d.]+)\s+CD =\s*([-\d.]+)", log)
    if mcl and mcd:
        cl_out, cm_out, cd_out = (float(mcl[-1][1]), float(mcd[-1][0]),
                                  float(mcd[-1][1]))
        if (abs(cl_out - cl) > 1e-4 or abs(cd_out - cd) > 1e-5
                or abs(cm_out - cm) > 1e-4):
            print(f"  WARNING {tag}: polar ({cl:.4f}/{cd:.5f}/{cm:.4f}) "
                  f"vs stdout ({cl_out:.4f}/{cd_out:.5f}/{cm_out:.4f}) "
                  "disagree; using polar values")

    # --- BL profiles from the DUMP file (12-field surface rows only) ---
    rows = []
    for line in (workdir / dump_name).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        if len(fields) < 12:
            continue  # wake rows have 8 fields
        rows.append([float(v) for v in fields[:12]])
    d = np.array(rows)
    x, dstar, cf = d[:, 1], d[:, 4], d[:, 6]
    k = int(np.argmin(x))  # LE pivot: rows [0..k] = upper TE->LE, rest = lower
    upper = np.column_stack([x[:k + 1][::-1], dstar[:k + 1][::-1],
                             cf[:k + 1][::-1]])  # LE -> TE
    lower = np.column_stack([x[k + 1:], dstar[k + 1:], cf[k + 1:]])

    return {
        "tag": tag, "version": version, "alpha": alpha,
        "cl": cl, "cd": cd, "cm": cm,
        "upper": upper, "lower": lower,
        "dstar_te_upper": float(dstar[0]),
        "dstar_te_lower": float(dstar[-1]),
    }


def run_inviscid(binary: Path, workdir: Path) -> dict:
    """Inviscid XFOIL point at the pinned M/alpha (no VISC, no XTR): the
    reference for 'XFOIL's own viscous decrement' (GV3.1 dcl clause).

    The 6.99 inviscid analysis prints no CL to stdout in batch (the ALFA
    printout is a viscous-mode line), so the point is taken through PACC
    polar accumulation + a one-point ASEQ; the polar save file carries
    cl/cm. PACC ignores inline args and prompts for the save file first,
    then the dump file (blank = none).
    """
    polar_name = "polar_inviscid.txt"
    script = (
        "NACA 0012\n"
        "PPAR\n"
        f"N {N_PANELS}\n"
        "\n\n"
        "OPER\n"
        f"MACH {MACH:.1f}\n"
        "PACC\n"
        f"{polar_name}\n"
        "\n"
        f"ASEQ {ALPHA_DEG:.1f} {ALPHA_DEG:.1f} 1.0\n"
        "PACC\n"
        "\n"
        "QUIT\n"
    )
    (workdir / "inviscid.in").write_text(script)
    proc = subprocess.run(
        [str(binary)], input=script, cwd=workdir,
        capture_output=True, text=True, timeout=600,
    )
    log = proc.stdout + proc.stderr
    (workdir / "inviscid.log").write_text(log)
    if proc.returncode != 0:
        raise RuntimeError(f"xfoil inviscid run failed; see {workdir}")
    data_rows = []
    for line in (workdir / polar_name).read_text().splitlines():
        tok = line.split()
        if len(tok) >= 7:
            try:
                data_rows.append([float(v) for v in tok])
            except ValueError:
                pass
    if not data_rows:
        raise RuntimeError("no inviscid polar data row found")
    alpha, cl, _cd, _cdp, cm = data_rows[-1][:5]
    return {"alpha": alpha, "cl": cl, "cm": cm}


def main():
    out = Path(__file__).parent
    binary = find_xfoil()
    print(f"using xfoil binary: {binary}")

    results = []
    with tempfile.TemporaryDirectory(prefix="xfoil_ref_") as tmp:
        workdir = Path(tmp)
        for tag, (xt, xb) in CASES.items():
            r = run_case(binary, workdir, tag, xt, xb)
            results.append(r)
            print(f"{tag}: alpha={r['alpha']:.3f}  cl={r['cl']:.4f}  "
                  f"cd={r['cd']:.5f}  cm={r['cm']:.4f}  "
                  f"dstar_TE upper={r['dstar_te_upper']:.5f} "
                  f"lower={r['dstar_te_lower']:.5f}  "
                  f"(XFOIL {r['version']})")
        inv = run_inviscid(binary, workdir)
        print(f"inviscid: cl={inv['cl']:.4f} cm={inv['cm']:.4f} "
              "(the viscous-decrement reference for GV3.1 dcl)")

    with open(out / "inviscid_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "m_inf", "alpha_deg", "cl", "cm"])
        w.writerow(["inviscid", f"{MACH:.2f}", f"{ALPHA_DEG:.3f}",
                    f"{inv['cl']:.6f}", f"{inv['cm']:.6f}"])
    print(f"wrote {out / 'inviscid_summary.csv'}")

    for r in results:
        csv_path = out / f"delta_star_cf_alpha2_m05_{r['tag']}.csv"
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["x_c", "surface", "dstar_over_c", "cf"])
            for surf, arr in (("upper", r["upper"]), ("lower", r["lower"])):
                for x, dstar, cf in arr:
                    w.writerow([f"{x:.6f}", surf, f"{dstar:.6e}", f"{cf:.6e}"])
        print(f"wrote {csv_path}")

    with open(out / "polar_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "m_inf", "re", "alpha_deg",
                    "x_tr_upper", "x_tr_lower", "cl", "cd", "cm"])
        for r in results:
            xt, xb = CASES[r["tag"]]
            w.writerow([r["tag"], f"{MACH:.2f}", f"{REYNOLDS:.4e}",
                        f"{r['alpha']:.3f}", f"{xt:.2f}", f"{xb:.2f}",
                        f"{r['cl']:.6f}", f"{r['cd']:.6e}", f"{r['cm']:.6f}"])
    print(f"wrote {out / 'polar_summary.csv'}")


if __name__ == "__main__":
    main()
