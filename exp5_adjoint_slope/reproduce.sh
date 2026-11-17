#!/usr/bin/env bash
# Reproduce Experiment 5: builds, runs, evaluation, figures and the README.
#
#   bash exp5_adjoint_slope/reproduce.sh --check     # list the steps in order; run nothing, write nothing
#   bash exp5_adjoint_slope/reproduce.sh --quick     # evaluation, figure data, figures and README from existing runs
#   bash exp5_adjoint_slope/reproduce.sh             # everything, from the builds on
#
# Environment:
#   PY       Python with numpy/matplotlib/pyyaml/cmocean (default: the pyocean conda env)
#   RUNS     run-output root (default /Volumes/ESTEBAN_MIT/portfolio_runs/exp5)
#   BUILDS   MITgcm build root (default /Volumes/ESTEBAN_MIT/portfolio_builds/exp5)
# Compute (laptop, <= 5 cores): G0 spin-up ~10 h; F1/F1b ~2.7 h; ensemble adjoint ~45 min; forward ensembles ~30 min;
# V1-V3 ~1 h; fine ensemble FC/FP ~6 h. Every run goes through scripts/run_adjoint.py or
# scripts/run_campaign.py, which refuse to start while another model process is running.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="${PY:-/opt/anaconda3/envs/pyocean/bin/python}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_NUM_THREADS=1
export EXP5_RUNS="${RUNS:-/Volumes/ESTEBAN_MIT/portfolio_runs/exp5}"
export EXP5_BUILDS="${BUILDS:-/Volumes/ESTEBAN_MIT/portfolio_builds/exp5}"
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
RA="exp5_adjoint_slope/scripts/run_adjoint.py"
RC="exp5_adjoint_slope/scripts/run_campaign.py"
step() { echo "=== $*"; }
run() { if [ "$CHECK" -eq 1 ]; then echo "+ $*"; else "$@"; fi; }
heavy() { if [ "$QUICK" -eq 1 ]; then echo "  (skipped in --quick)"; else run "$@"; fi; }

step "build: forward (chan10, chan2p5) and Tapenade adjoint/tangent-linear (chan10_tap)"
heavy scripts/build_mitgcm.sh exp5_adjoint_slope/config/chan10 "$EXP5_BUILDS/chan10"
heavy scripts/build_mitgcm.sh exp5_adjoint_slope/config/chan2p5 "$EXP5_BUILDS/chan2p5" --mpi
for k in adj tlm fwd; do heavy scripts/build_mitgcm_tap.sh exp5_adjoint_slope/config/chan10_tap "$EXP5_BUILDS/chan10_tap_$k" "$k"; done

step "spin-up G0: coarse model, then the fine spin-up F0/F0b and the eddy spin-up F1/F1b"
heavy "$PY" "$RC" --cases G0 --np 1
heavy "$PY" "$RC" --cases 'F0*' --np 5
heavy "$PY" "$RC" --cases 'F1*' --np 5

step "adjoint gate: gradient checks on short windows (GS, GSU, GSF)"
heavy "$PY" "$RA" --case 'GS*' --max-concurrent 5

step "ensemble adjoint EA: ten 90-day windows from the G0 yearly pickups, then the PD2 patterns"
heavy "$PY" "$RA" --case 'EA*' --max-concurrent 5
heavy "$PY" -m exp5_adjoint_slope.analysis.pd2_patterns ens

step "forward ensembles: V3E, H5A, H5BC"
for g in 'V3E_*' 'H5A*' 'H5BC_*'; do heavy "$PY" "$RA" --case "$g" --max-concurrent 5; done

step "deterministic V1-V3: GM30, GV1, GV1_tlm, GV2d"
heavy "$PY" "$RA" --case GM30
heavy "$PY" -m exp5_adjoint_slope.analysis.pd2_patterns v2d
for g in 'GV1' 'GV1_tlm' 'GV2d_*'; do heavy "$PY" "$RA" --case "$g" --max-concurrent 5; done

step "fine ensemble: FC1-3 and FP1-3"
heavy "$PY" "$RC" --cases 'F[CP]*' --np 5

step "evaluation: validation.json, sensitivity.json, hypotheses.json, crosscheck.json"
run "$PY" -m exp5_adjoint_slope.analysis.pd2_eval
run "$PY" -m exp5_adjoint_slope.analysis.hypotheses
run "$PY" -m exp5_adjoint_slope.analysis.crosscheck

step "figure data: mean_state.json, diagnostics.json"
run "$PY" -m exp5_adjoint_slope.analysis.figdata

step "figures: Figs. 0-5 and the spin-up animation"
run "$PY" -m exp5_adjoint_slope.analysis.figures
run "$PY" -m exp5_adjoint_slope.analysis.animate

step "README: render from the template"
run "$PY" exp5_adjoint_slope/scripts/render_readme.py
