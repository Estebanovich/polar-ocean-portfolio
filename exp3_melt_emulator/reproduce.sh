#!/usr/bin/env bash
# Reproduce Experiment 3: hypotheses, figures and README from the committed JSON, optionally re-running
# the cross-validation. The final test is NEVER re-run: it was evaluated once, and figures/final_test.json
# and its start marker are never touched by this script.
#
#   bash exp3_melt_emulator/reproduce.sh --quick   # hypotheses, bootstrap, probabilistic JSON, Figures 1-4, README
#                                                  # from figures/*.json
#   bash exp3_melt_emulator/reproduce.sh           # full: re-run the noise floor, run_protocol.py --cv and the
#                                                  # post hoc profile-test diagnosis, compare each with the
#                                                  # committed record, restore it, then --quick
#
# Environment (see the repository RUNNING.md):
#   PY       Python for the analysis, numpy/matplotlib/cmocean (default: the pyocean conda env)
#   ANEMOI   Python with PyTorch and scikit-learn for the protocol (default: the anemoi conda env)
#   SOURCE_DATE_EPOCH  timestamp written into the PDFs (default 0)
# Data outside the repository (paths fixed in the code, as in the original runs):
#   /Volumes/LaCie/portfolio_runs/exp2/exp3_dataset.npz   the Exp. 2 dataset (sha256 checked before use)
#   /Volumes/LaCie/portfolio_runs/exp3/final_maps.npz     the final-test maps (sha256 in final_test.json)
# --quick needs nothing else: without the dataset or the maps, hypotheses.json and the README are still
# rebuilt from the JSON and Figures 1-4 are kept as committed.
# The full mode also needs the Exp. 2 run output (/Volumes/LaCie/portfolio_runs/exp2, for the noise floor)
# and about 11 h of wall clock on 5 cores. run_protocol.py --cv verifies the frozen configuration
# (figures/protocol_config.json), the dataset and the folds, and reuses the finished units in
# /Volumes/LaCie/portfolio_runs/exp3/cv_units: with the recorded units in place NOTHING is refitted and the
# --cv check only shows that cv_results.json reassembles identically (a notice says so); move that
# directory aside to refit every unit (the script never deletes it). It is run
# under Amendment 1 of the pre-registration (wall-clock budget 20 h instead of 8 h, protocol unchanged),
# exactly as the recorded run. No committed record is changed by the full mode: the noise floor is
# recomputed to a scratch file and compared with figures/noise_floor.json (an input of the frozen
# configuration); --cv and the profile-test diagnosis rewrite figures/cv_results.json and
# figures/profile_test_diag.json with new session metadata, so scripts/rerun_compare_restore.sh backs
# them up, compares the re-run results (analysis/compare_records.py, timing ignored) and restores them.
# The re-run copies are kept in /Volumes/LaCie/portfolio_runs/exp3/reproduce/.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="${PY:-/opt/anaconda3/envs/pyocean/bin/python}"
ANEMOI="${ANEMOI:-/opt/anaconda3/envs/anemoi/bin/python}"
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-0}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 VECLIB_NUM_THREADS=1 MKL_NUM_THREADS=1

QUICK=0
for a in "$@"; do
  case "$a" in
    --quick) QUICK=1 ;;
    -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 0 ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done

cd "$ROOT"
step() { echo; echo "=== $*"; }
A="$HERE/analysis"
F="$HERE/figures"
DATASET=/Volumes/LaCie/portfolio_runs/exp2/exp3_dataset.npz
RUNS3=/Volumes/LaCie/portfolio_runs/exp3
KEEPAWAKE=""; command -v caffeinate > /dev/null && KEEPAWAKE="caffeinate -is"   # macOS: keep the machine awake
AMENDMENT="Enmienda 1 del PD0 (Esteban 2026-11-10): presupuesto de reloj 20 h; protocolo científico sin cambios"

# The final test is evaluated once; nothing below may create, change or remove its record.
final_state() { for p in "$F/final_test.json" "$RUNS3/final_started.json"; do
  [ -e "$p" ] && shasum -a 256 "$p" || echo "absent $p"; done; }
FINAL_BEFORE="$(final_state)"

if [ "$QUICK" -eq 0 ]; then
  [ -f "$DATASET" ] || { echo "dataset $DATASET not found (mount LaCie)" >&2; exit 2; }
  step "noise floor, recomputed and compared with the frozen figures/noise_floor.json"
  mkdir -p "$RUNS3/reproduce"
  "$PY" "$A/noise_floor.py" --out "$RUNS3/reproduce/noise_floor.npz" \
    --summary "$RUNS3/reproduce/noise_floor.json" > /dev/null
  CHECKS=()
  "$PY" - "$F/noise_floor.json" "$RUNS3/reproduce/noise_floor.json" <<'EOF' || CHECKS+=("noise floor differs")
import json, sys
a, b = (json.load(open(p)) for p in sys.argv[1:])
keys = ("pooled_median_rms", "median_abs_integrated_rel_diff", "min_rms", "max_rms")
diff = {k: (a[k], b[k]) for k in keys if a[k] != b[k]}
print("noise floor identical" if not diff else f"noise floor DIFFERS: {diff}")
sys.exit(1 if diff else 0)
EOF
  step "cross-validation: run_protocol.py --cv (Amendment 1; never --final)"
  export RERUN_DIR="$RUNS3/reproduce"
  st=0; bash "$HERE/scripts/rerun_compare_restore.sh" cv "$F/cv_results.json" -- \
    $KEEPAWAKE "$ANEMOI" "$HERE/scripts/run_protocol.py" --cv --wall-limit-hours 20 \
    --amendment "$AMENDMENT" || st=$?
  [ "$st" -eq 0 ] || CHECKS+=("cv_results differs or --cv failed (exit $st)")
  if grep "START --cv" "$RUNS3/cv_progress.log" | tail -1 | grep -q " pending 0 "; then
    echo; echo "=== NOTE: every cross-validation unit already existed in $RUNS3/cv_units, so NOTHING WAS REFITTED:"
    echo "=== the check above only shows that cv_results.json reassembles identically from the recorded units."
    echo "=== For a refit from scratch, move that directory aside yourself (this script never deletes it), ~11 h."
  fi
  step "post hoc profile-test diagnosis (one M7 refit)"
  st=0; bash "$HERE/scripts/rerun_compare_restore.sh" profile "$F/profile_test_diag.json" -- \
    "$ANEMOI" "$A/profile_test_diag.py" > /dev/null || st=$?
  [ "$st" -eq 0 ] || CHECKS+=("profile_test_diag differs or failed (exit $st)")
fi

step "hypotheses (H3a-H3c) and ranking";        "$PY" "$A/hypotheses.py" > /dev/null
step "bootstrap intervals over runs";           "$PY" "$A/bootstrap.py" > /dev/null
step "coverage, CRPS and skill";                "$PY" "$A/probabilistic.py" > /dev/null
MAPS="$RUNS3/final_maps.npz"
if [ -f "$DATASET" ] && [ -f "$MAPS" ]; then
  step "Figures 1, 2 and 4";                    "$PY" "$A/figures.py" > /dev/null
  step "Figure 3: physics diagnostics";         "$PY" "$A/physics_diag.py" > /dev/null
else
  echo; echo "=== dataset or final-test maps not found: Figures 1-4 kept as committed"
fi
GEOM_NC="$HERE/../exp2_isomip_cavity/data_external/Ocean1_input_geom_v1.01.nc"
if [ -f "$GEOM_NC" ]; then
  step "Figure 0: protocol schematic";          "$PY" "$F/src/fig0_protocol.py" > /dev/null
else
  echo; echo "=== ISOMIP+ geometry not found (exp2_isomip_cavity/data_external/SOURCES.md): Figure 0 kept as committed"
fi
step "README";                                  "$PY" "$HERE/scripts/render_readme.py"

[ "$(final_state)" = "$FINAL_BEFORE" ] || { echo "the final-test record changed: this must not happen" >&2; exit 1; }
echo; echo "done: figures and metrics in exp3_melt_emulator/figures/, README rendered; final test untouched."
if [ "$QUICK" -eq 0 ] && [ "${#CHECKS[@]}" -gt 0 ]; then
  printf 'reproduction check FAILED: %s\n' "${CHECKS[@]}" >&2; exit 1
fi
