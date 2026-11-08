#!/usr/bin/env bash
# Reproduce Experiment 1: figures, metrics and README, optionally from scratch.
#
#   bash exp1_fjord_plumes/reproduce.sh --quick          # metrics, figures and README from existing run output
#   bash exp1_fjord_plumes/reproduce.sh --quick --gifs   # ... plus the animated GIFs (slow)
#   bash exp1_fjord_plumes/reproduce.sh                  # full: generate inputs, build, run every campaign, then --quick
#
# Environment (see the repository RUNNING.md):
#   PY               Python with numpy/scipy/matplotlib/pyyaml (default: the pyocean conda env)
#   RUNS             run-output root (default /Volumes/LaCie/portfolio_runs/exp1; must be mounted)
#   VECTORCRAFT_CLI  optional: path to vectorcraft-cli v0.22.0 to redraw the Figure 0 schematics
# The full mode needs MITgcm (read-only source in /Volumes/ESTEBAN_MIT/MITgcm), MPI and gfortran, and a
# pyOCEAN checkout for the (c) open-boundary cross-check. It takes several days on a laptop (run timings are in
# notes/); every case is run through scripts/run_campaign.py, which refuses to start while another MITgcm run
# is active and writes run_status.json, which the analysis checks before using any run.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="${PY:-/opt/anaconda3/envs/pyocean/bin/python}"
RUNS="${RUNS:-/Volumes/LaCie/portfolio_runs/exp1}"; export EXP1_RUNS="$RUNS"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

QUICK=0; GIFS=0
for a in "$@"; do
  case "$a" in
    --quick) QUICK=1 ;;
    --gifs) GIFS=1 ;;
    -h|--help) sed -n '2,16p' "$0"; exit 0 ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done

cd "$ROOT"
[ -d "$RUNS" ] || { echo "run root $RUNS not found (mount LaCie or set RUNS)" >&2; exit 2; }
step() { echo; echo "=== $*"; }
A="$HERE/analysis"

if [ "$QUICK" -eq 0 ]; then
  campaign() { step "campaign $1: $2 (np $3)"; "$PY" "$HERE/scripts/run_campaign.py" --config "$1" --cases "$2" --np "$3" --skip-ok; }
  # (a) 2-D line plume, 48 h; resolution and sensitivity studies; the 6-h animation rerun
  campaign a_line2d 'a_line2d_*' 1
  campaign conv_dx10 conv_dx10 1
  campaign conv_dx5 conv_dx5 1
  campaign conv_dx2p5 conv_dx2p5 1
  campaign sens_line2d 'sens_*' 1
  campaign a_line2d_anim a_line2d_twolayer_q0.03_anim 1
  # (b) two 3-D plumes: 10 m campaign with no-discharge backgrounds, then the pre-registered 5 m subset
  campaign b_twoplume3d 'b_*' 5
  campaign b_twoplume3d_dx5 'b_dx5_*' 5
  # (c) 3-D fjord (two-stage runs are handled by run_campaign.py)
  campaign c_fjord3d 'c_fjord_*' 5
  # (c) open-boundary cross-check with pyOCEAN (independent model)
  for c in pyo_shelf pyo_shelf_runoff; do
    step "pyOCEAN cross-check $c"; "$PY" "$HERE/crosscheck/pyocean_c/run_pyo_c.py" --case "$c"
  done
fi

# Open-boundary criteria of the 5 m subset (requirement before interpreting it; exits non-zero on failure)
for c in single bg_twolayer_nodischarge twolayer_x0_0 twolayer_x0_20 twolayer_x0_80; do
  r="$RUNS/b_dx5_${c/#single/twolayer_single}"
  step "OB criteria $(basename "$r")"; "$PY" "$A/ob_eval.py" "$r" --ob 196:199 --interior 20:180 > /dev/null
done

step "(a) Figure 1 and convergence";   "$PY" "$A/fig1_plume_validation.py"; "$PY" "$A/pd2_convergence_eval.py"
step "(b) plume excess 10 m and 5 m";  "$PY" "$A/excess_clean_bg.py"; "$PY" "$A/excess_dx5.py"; "$PY" "$A/d2_resolution.py"
step "(b) Figure 2";                   "$PY" "$A/fig2_two_plume_melt.py"
step "(c) Figure 3 and shelf amplitude"; "$PY" "$A/fig3_exchange_partition.py"; "$PY" "$A/c_shelf_amplitude.py"
step "(c) pyOCEAN cross-check";        "$PY" "$HERE/crosscheck/pyocean_c/eval_pyo_c.py"

if [ -n "${VECTORCRAFT_CLI:-}" ]; then
  step "Figure 0 schematics (VectorCraft)"
  "$PY" "$HERE/figures/src/fig0_schematic.py" --cli "$VECTORCRAFT_CLI"
  "$PY" "$HERE/figures/src/fig0_schematic3d.py" --cli "$VECTORCRAFT_CLI"
else
  echo; echo "=== Figure 0 kept as committed (set VECTORCRAFT_CLI to redraw it)"
fi

if [ "$GIFS" -eq 1 ]; then
  step "animations"
  for c in a_line2d_twolayer_q0.03_anim a_line2d_twolayer_q0.03 a_line2d_sermilik_q0.03; do "$PY" "$A/animate_a.py" "$RUNS/$c"; done
  for c in b_3d_twolayer_single b_3d_twolayer_x0_0 b_3d_twolayer_x0_80 b_3d_sermilik_x0_0; do "$PY" "$A/animate_b.py" "$RUNS/$c" --out "$HERE/figures/fig2_anim_$c.gif"; done
  for c in plume wind shelf all; do "$PY" "$A/animate_c.py" "$RUNS/c_fjord_$c" --out "$HERE/figures/fig3_anim_c_fjord_$c.gif"; done
fi

step "README"; "$PY" "$HERE/scripts/render_readme.py"
echo; echo "done: figures and metrics in exp1_fjord_plumes/figures/, README rendered."
