#!/usr/bin/env bash
# Reproduce Experiment 6: builds, runs, evaluation, figures and the README.
#
#   bash exp6_southern_carbon/reproduce.sh --check     # list the steps in order; run nothing, write nothing
#   bash exp6_southern_carbon/reproduce.sh --quick     # evaluation, figures and README from existing runs
#   bash exp6_southern_carbon/reproduce.sh             # everything, from the builds on
#
# Environment:
#   PY       Python with numpy/matplotlib/pyyaml/cmocean (default: the pyocean conda env)
#   RUNS     run-output root (default /Volumes/LaCie/portfolio_runs/exp6)
#   BUILDS   MITgcm build root (default /Volumes/ESTEBAN_MIT/portfolio_builds/exp6)
# Compute (laptop, <= 5 cores): G0 ~20 h (10 stages of 100 tracer years, 5 ranks); G0s ~6 h;
# GC, GT and the eight coarse branches ~58 h; fine F0/FT ~230 h. Every run goes through
# scripts/run_campaign.py, which refuses to start while another model process is running.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="${PY:-/opt/anaconda3/envs/pyocean/bin/python}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export EXP6_RUNS="${RUNS:-/Volumes/LaCie/portfolio_runs/exp6}"
export EXP6_BUILDS="${BUILDS:-/Volumes/ESTEBAN_MIT/portfolio_builds/exp6}"
export PYTHONPATH="$ROOT"

CHECK=0; QUICK=0
for a in "$@"; do
  case "$a" in
    --check) CHECK=1 ;;
    --quick) QUICK=1 ;;
    -h|--help) sed -n '2,15p' "$0"; exit 0 ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done

cd "$ROOT"
RC="exp6_southern_carbon/scripts/run_campaign.py"
STAGE=788400                                   # 100 tracer years at deltaTClock = 4000 s
step() { echo "=== $*"; }
run() { if [ "$CHECK" -eq 1 ]; then echo "+ $*"; else "$@"; fi; }
heavy() { if [ "$QUICK" -eq 1 ]; then echo "  (skipped in --quick)"; else run "$@"; fi; }

step "build: g50 (serial), g50p5 and f10 (MPI); run_campaign.py builds when the code changed"
heavy scripts/build_mitgcm.sh exp6_southern_carbon/config/g50 "$EXP6_BUILDS/g50"
heavy scripts/build_mitgcm.sh exp6_southern_carbon/config/g50p5 "$EXP6_BUILDS/g50p5" --mpi
heavy scripts/build_mitgcm.sh exp6_southern_carbon/config/f10 "$EXP6_BUILDS/f10" --mpi

step "spin-up G0: accelerated tracers, 10 stages of 100 tracer years (V6a)"
heavy "$PY" "$RC" --cases G0 --n-timesteps "$STAGE"
for k in 2 3 4 5 6 7 8 9 10; do heavy "$PY" "$RC" --cases G0 --n-timesteps "$STAGE" --continue; done

step "spin-up G0s: synchronous years after G0"
heavy "$PY" "$RC" --cases G0s

step "coarse campaign: GC and GT with their wind branches, one invocation (branches poll for the parent's pickups)"
heavy "$PY" "$RC" --cases 'G[CT]*' --max-concurrent 5 --wait-hours 80

step "fine campaign: F0 from GT in 1900 and 2060, then the FT control/wind members"
heavy "$PY" "$RC" --cases 'F0_*'
heavy "$PY" "$RC" --cases 'FT_*'

step "evaluation: spinup.json, validation.json, figdata.json, hypotheses.json, crosscheck.json"
run "$PY" -m exp6_southern_carbon.analysis.figdata --runs-root "$EXP6_RUNS"
run "$PY" -m exp6_southern_carbon.analysis.hypotheses --runs-root "$EXP6_RUNS"
if [ -f exp6_southern_carbon/analysis/crosscheck.py ]; then
  run "$PY" -m exp6_southern_carbon.analysis.crosscheck --runs-root "$EXP6_RUNS"
else
  echo "  (crosscheck.py not present yet: H6c pending)"
fi

step "figures: Figs. 0-5 (a figure whose JSON is missing is skipped)"
run "$PY" -m exp6_southern_carbon.analysis.figures

step "README: render from the template"
if [ -f exp6_southern_carbon/figures/crosscheck.json ]; then
  run "$PY" exp6_southern_carbon/scripts/render_readme.py
else
  run "$PY" exp6_southern_carbon/scripts/render_readme.py --allow-pending
fi
