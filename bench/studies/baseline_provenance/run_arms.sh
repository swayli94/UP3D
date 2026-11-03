#!/usr/bin/env bash
# baseline-provenance: four arms on the 11 gate files that hold the 12 reds.
# Order is registered before the run: native_cold -> native_warm -> generic_cold -> generic_warm.
# Runs in a detached worktree pinned at 080f291 so concurrent edits to the main checkout
# cannot contaminate it. Writes only under OUT and the two fresh cache dirs it creates.
set -uo pipefail

OUT=/home/lrz/code/UP3D/bench/results/baseline_provenance_20261104
WT=/home/lrz/code/up3d-baseline-prov
EXPECT_HEAD=080f29152f14e15a70e9e1c0a62e9b283f127f19
THREADS=${PROV_THREADS:-16}
FILES=(
  tests/C/test_C05_nozzle_quasi1d.py tests/C/test_C07_karman_trefftz.py tests/C/test_C08_ringleb.py
  tests/D/test_D05_euler_naca0012.py tests/D/test_D06_euler_rae2822.py tests/D/test_D08_rans_naca0012.py
  tests/D/test_D11_naca0012_experiment_bias.py tests/D/test_D12_rae2822_experiment_bias.py
  tests/D/test_D13_ibl_vs_xfoil.py tests/D/test_D14_alpha0_transonic_nonuniqueness.py
  tests/D/test_D15_hex_naca_boundary.py
)

echo $$ > "$OUT/runner.pid"
trap 'rc=$?; echo "EXIT rc=$rc at $(date -Is)" >> "$OUT/progress.txt"; date -Is > "$OUT/DONE"' EXIT
echo "START $(date -Is) threads=$THREADS host=$(hostname) load=$(cut -d" " -f1-3 /proc/loadavg)" >> "$OUT/progress.txt"

head=$(git -C "$WT" rev-parse HEAD)
if [ "$head" != "$EXPECT_HEAD" ] || [ -n "$(git -C "$WT" status --porcelain --untracked-files=no)" ]; then
  echo "ABORT worktree not at clean $EXPECT_HEAD (head=$head)" >> "$OUT/progress.txt"; exit 2
fi

run_arm() {
  local name=$1 cachedir=$2 cpu=$3 kind=$4
  local d="$OUT/$name"
  if [ -e "$d" ]; then echo "ABORT $name: $d exists" >> "$OUT/progress.txt"; exit 2; fi
  mkdir -p "$d"
  if [ "$kind" = cold ]; then
    if [ -e "$cachedir" ]; then echo "ABORT $name: cold arm but $cachedir exists" >> "$OUT/progress.txt"; exit 2; fi
    mkdir -p "$cachedir"
  fi
  find "$cachedir" -type f | wc -l > "$d/cache_files_before.txt"
  local envs=(PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=$THREADS OMP_NUM_THREADS=$THREADS
              OPENBLAS_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS MPLBACKEND=Agg
              NUMBA_CACHE_DIR=$cachedir PROV_OUT_DIR=$d PROV_WORKTREE=$WT PYTHONPATH=$OUT)
  [ -n "$cpu" ] && envs+=(NUMBA_CPU_NAME=$cpu)
  git -C "$WT" rev-parse HEAD > "$d/head.txt"
  git -C "$WT" status --porcelain --untracked-files=no > "$d/status_before.txt"
  cut -d" " -f1-3 /proc/loadavg > "$d/load_before.txt"
  ( cd "$WT" && conda run --no-capture-output -n up3d env "${envs[@]}" python "$OUT/codegen_env.py" ) > "$d/codegen.json" 2>&1
  echo "BEGIN $name $(date -Is) cache_files_before=$(cat "$d/cache_files_before.txt")" >> "$OUT/progress.txt"
  local t0=$(date +%s)
  ( cd "$WT" && conda run --no-capture-output -n up3d env "${envs[@]}" \
      python -m pytest "${FILES[@]}" -q -rfE --durations=15 -p no:cacheprovider -p prov_plugin \
      --junitxml="$d/junit.xml" ) > "$d/pytest.log" 2>&1
  local rc=$?
  local wall=$(( $(date +%s) - t0 ))
  find "$cachedir" -type f | wc -l > "$d/cache_files_after.txt"
  git -C "$WT" status --porcelain --untracked-files=no > "$d/status_after.txt"
  cut -d" " -f1-3 /proc/loadavg > "$d/load_after.txt"
  echo "$name rc=$rc wall=$wall end=$(date -Is) summary=\"$(tail -1 "$d/pytest.log")\"" >> "$OUT/progress.txt"
  # pytest rc 1 = some tests failed (expected: the 12 reds). Anything else is an instrument problem.
  if [ $rc -ne 0 ] && [ $rc -ne 1 ]; then echo "ABORT $name: pytest rc=$rc" >> "$OUT/progress.txt"; exit 3; fi
}

run_arm native_cold  "$OUT/cache_native"  ""        cold
run_arm native_warm  "$OUT/cache_native"  ""        warm
run_arm generic_cold "$OUT/cache_generic" "generic" cold
run_arm generic_warm "$OUT/cache_generic" "generic" warm

( cd "$WT" && conda run --no-capture-output -n up3d env PYTHONNOUSERSITE=1 python "$OUT/compare.py" "$OUT" ) \
  > "$OUT/compare.log" 2>&1 || echo "compare.py failed (data intact; re-run it by hand)" >> "$OUT/progress.txt"
echo "ALL ARMS DONE $(date -Is)" >> "$OUT/progress.txt"
