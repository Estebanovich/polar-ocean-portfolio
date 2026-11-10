#!/usr/bin/env bash
# Reproduce Experiment 2: metrics, figures, the Exp. 3 export and the README, optionally from scratch.
#
#   bash exp2_isomip_cavity/reproduce.sh --quick          # metrics, figures, Exp. 3 export and README from existing run output
#   bash exp2_isomip_cavity/reproduce.sh --quick --gifs   # ... plus the meltwater-tracer animation
#   bash exp2_isomip_cavity/reproduce.sh                  # full: every run campaign in its original order, then --quick
#
# Environment (see the repository RUNNING.md):
#   PY     Python with numpy/scipy/matplotlib/pyyaml/netCDF4 (default: the pyocean conda env)
#   RUNS   run-output root (default /Volumes/LaCie/portfolio_runs/exp2; must be mounted)
#   SOURCE_DATE_EPOCH  timestamp written into the PDFs (default 0, so a rerun gives identical bytes)
#   VECTORCRAFT_CLI  optional: path to vectorcraft-cli v0.22.0 to redraw Figure 0 (a)-(c)
#   BLENDER          optional, with VECTORCRAFT_CLI: path to blender to re-render Figure 0 (d) from the geometry file
# The full mode needs MITgcm (read-only source in /Volumes/ESTEBAN_MIT/MITgcm), MPI and gfortran, and the
# external data of data_external/SOURCES.md. It takes about 20 h of compute on a laptop (timings in
# notes/pd1_calibracion.md). Every run goes through scripts/run_campaign.py, which refuses to start while
# another model run is active and writes run_status.json; the analysis refuses any run that is not OK.
# The campaign order is the one actually used:
#   1. scripts/queue_calibration.sh  K1 calibration (Gamma_T 0.015/0.03/0.06 at 4 and 2 km) and the
#                                    verification runs at Gamma_T,cal (ruling R21: run suffix _r21)
#   2. scripts/queue_production.sh   the 28-run ensemble and the 4 km reference (np 1, 5 at a time),
#                                    the N1-N3 check, then Ocean1 at 2 km (np 5, four 5-year stages)
#   3. scripts/queue_r24.sh          the post hoc R24 exposure runs hol_T* (5 runs, 4 km)
# Each queue is wrapped in caffeinate: an idle Mac otherwise pauses the runs.
# Note: the queue scripts hardcode RUNS=/Volumes/LaCie/portfolio_runs/exp2 (RUNS above does not override
# them), and queue_production.sh takes Gamma_T,cal from cases.yaml (0.029356 at 4 km, 0.030959 at 2 km).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="${PY:-/opt/anaconda3/envs/pyocean/bin/python}"
RUNS="${RUNS:-/Volumes/LaCie/portfolio_runs/exp2}"; export EXP2_RUNS="$RUNS"
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-0}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

QUICK=0; GIFS=0
for a in "$@"; do
  case "$a" in
    --quick) QUICK=1 ;;
    --gifs) GIFS=1 ;;
    -h|--help) sed -n '2,26p' "$0"; exit 0 ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done

cd "$ROOT"
[ -d "$RUNS" ] || { echo "run root $RUNS not found (mount LaCie or set RUNS)" >&2; exit 2; }
step() { echo; echo "=== $*"; }
A="$HERE/analysis"
F="$HERE/figures"
R=_r21                                     # R21 run-directory suffix: the void runs are never read

if [ "$QUICK" -eq 0 ]; then
  for q in queue_calibration queue_production queue_r24; do
    step "campaign $q (caffeinate)"
    caffeinate -is bash "$HERE/scripts/$q.sh"
  done
fi

step "K1 calibration and verification (pd0 K1)"
for res in 4km 2km; do
  c=$([ "$res" = 4km ] && echo cal4 || echo cal2)
  "$PY" "$A/calibrate.py" --res "$res" --out "$F/calibration_${res}$R.json" \
    "$RUNS/${c}_gT0.015$R" "$RUNS/${c}_gT0.03$R" "$RUNS/${c}_gT0.06$R"
  "$PY" "$A/calibrate.py" --verify --res "$res" --out "$F/verify_${res}$R.json" "$RUNS/ocean0_${res}$R"
done
step "K1, V1, N1-N3, R15, R17 and the 2 vs 4 km resolution test"
"$PY" "$A/pd1_calibration_eval.py" --runs-root "$RUNS" > /dev/null
if [ -f "$HERE/data_external/osf_jpl/Ocean1_yosi_COM_MITgcm_JPL_new.nc" ]; then
  step "post hoc comparison with MITgcm-JPL";   "$PY" "$A/jpl_compare.py" --runs-root "$RUNS" > /dev/null
else
  echo; echo "=== MITgcm-JPL output not found (data_external/SOURCES.md): jpl_comparison.json kept as committed"
fi
step "Figure 1: Ocean1 (V2, V3, N1-N3)";       "$PY" "$A/fig1_ocean0_validation.py" --run "$RUNS/ocean1_2km$R" > /dev/null
step "ensemble table (N1-N3)";                 "$PY" "$A/ensemble.py" --runs-root "$RUNS" --out "$F/ensemble_table.json" > /dev/null
step "Figure 2 (H1) and Figure 3 (H2)"
"$PY" "$A/fig2_melt_sensitivity.py" --table "$F/ensemble_table.json" > /dev/null
"$PY" "$A/fig3_meltwater_export.py" --table "$F/ensemble_table.json" --runs-root "$RUNS" > /dev/null
step "post hoc H1: R23, R24 part 2, R24 readout"
"$PY" "$A/h1_posthoc_tf.py" --runs-root "$RUNS" > /dev/null
"$PY" "$A/h1_posthoc_tf.py" --runs-root "$RUNS" --quad > /dev/null
"$PY" "$A/h1_posthoc_tf.py" --runs-root "$RUNS" --r24 > /dev/null
step "Exp. 3 dataset";                         "$PY" "$A/export_exp3.py" --runs-root "$RUNS" --out "$RUNS/exp3_dataset.npz"

if [ "$GIFS" -eq 1 ]; then
  step "animation"; "$PY" "$A/animate_tracer.py" "$RUNS/ocean0_4km$R" --out "$F/fig3_anim_ocean0_4km$R.gif"
fi

if [ -n "${VECTORCRAFT_CLI:-}" ]; then
  step "Figure 0 (a)-(c): schematic (VectorCraft)"
  "$PY" "$F/src/fig0_schematic.py" --cli "$VECTORCRAFT_CLI"
  if [ -n "${BLENDER:-}" ] && [ -f "$HERE/data_external/Ocean1_input_geom_v1.01.nc" ]; then
    step "Figure 0 (d): geometry (Blender + VectorCraft)"
    "$PY" "$F/src/fig0_geometry.py" --all --blender "$BLENDER" --cli "$VECTORCRAFT_CLI"
  else
    echo; echo "=== Figure 0 (d) kept as committed (set BLENDER and provide the geometry file to re-render it)"
  fi
else
  echo; echo "=== Figure 0 kept as committed (set VECTORCRAFT_CLI to redraw it)"
fi

step "README"; "$PY" "$HERE/scripts/render_readme.py"
echo; echo "done: figures and metrics in exp2_isomip_cavity/figures/, README rendered."
