#!/usr/bin/env bash
# Reproduce Experiment 4: metrics, figures and the README, optionally from scratch.
#
#   bash exp4_arctic_bayiw/reproduce.sh --quick             # every JSON, figure and the README from existing run output
#   bash exp4_arctic_bayiw/reproduce.sh                     # full: inputs, builds and every run, then --quick, then
#                                                           #       the regenerated JSONs compared with the committed ones
#   bash exp4_arctic_bayiw/reproduce.sh [--quick] --dry-run # print the plan; run nothing, write nothing
#
# Environment (see the repository RUNNING.md):
#   PY       Python with numpy/scipy/matplotlib/xarray/pyyaml/gsw/pyproj (default: the pyocean conda env)
#   RUNS     run-output root (default /Volumes/LaCie/portfolio_runs/exp4; must exist). The full mode uses a NEW root
#            (default /Volumes/LaCie/portfolio_runs/exp4_rerun, created if missing), so the stored runs are never overwritten
#   ORIG_RUNS  the stored runs (default /Volumes/LaCie/portfolio_runs/exp4): full mode, build-control comparison
#   BUILDS   MITgcm build root (default /Volumes/ESTEBAN_MIT/portfolio_builds/exp4)
#   VECTORCRAFT_CLI  optional: path to vectorcraft-cli v0.22.0 to redraw Figure 0
# Both modes refuse to start while an MPI launcher or a model process is running (pfcommon.campaign, matched by
# process name). The quick mode reads the runs (about 6 to 10 h of reads from the run disk: every cached analysis is
# recomputed, see below), the thesis run in /Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan (V1 and the
# contrast with the thesis scripts) and the external data of data_external/SOURCES.md (profiles, bathymetry).
# Caches: the partition caches (partition.py) are keyed on the source of their kernel functions, not on every
# helper they call, and the library cache of thesis_contrast.py is not keyed on thesis_contrast.py/flux.py, so
# both are run with --recompute here. The full mode needs MITgcm (read-only source), MPI and gfortran, and about
# 33 h of compute on a laptop (timings in notes/pd1_config.md and figures/pd1_timing.json); every run goes through
# scripts/run_campaign.py, which writes run_status.json, and the analysis refuses any run that is not OK.
# The PD1 timing and diagnostic short runs (figures/pd1_timing.json) are a record of the PD1 decisions and are not
# rerun; pd1_timing.json is regenerated from them when they exist under RUNS.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="${PY:-/opt/anaconda3/envs/pyocean/bin/python}"
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-0}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

QUICK=0; DRY=0
for a in "$@"; do
  case "$a" in
    --quick) QUICK=1 ;;
    --dry-run) DRY=1 ;;
    -h|--help) sed -n '2,28p' "$0"; exit 0 ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done
if [ "$QUICK" -eq 1 ]; then RUNS="${RUNS:-/Volumes/LaCie/portfolio_runs/exp4}"
else RUNS="${RUNS:-/Volumes/LaCie/portfolio_runs/exp4_rerun}"; fi
BUILDS="${BUILDS:-/Volumes/ESTEBAN_MIT/portfolio_builds/exp4}"
ORIG_RUNS="${ORIG_RUNS:-/Volumes/LaCie/portfolio_runs/exp4}"      # the stored runs (full mode: build-control comparison)
export EXP4_RUNS="$RUNS" EXP4_BUILDS="$BUILDS"

cd "$ROOT"
step() { echo; echo "=== $*"; }
run() { if [ "$DRY" -eq 1 ]; then echo "+ $*"; else "$@"; fi; }

# refuse while a model run is active (process names, never a command-line pattern)
if ! "$PY" -c "import sys; sys.path.insert(0, '$ROOT'); from pfcommon import campaign as cp; p = cp.other_mitgcm_processes(); print(*p, sep='\n', file=sys.stderr); sys.exit(1 if p else 0)"; then
  echo "refusing to start: an MPI launcher or a model process is active (above)" >&2
  exit 3
fi
# the full mode writes to a new run root: create it (the quick mode needs the existing one)
[ "$QUICK" -eq 1 ] || run mkdir -p "$RUNS"
[ "$DRY" -eq 1 ] || [ -d "$RUNS" ] || { echo "run root $RUNS not found (mount the run disk or set RUNS)" >&2; exit 2; }

A="$HERE/analysis"
S="$HERE/scripts"
F="$HERE/figures"
R="$RUNS"

if [ "$QUICK" -eq 0 ]; then
  BACKUP="$RUNS/committed_figures_json_$(date +%Y%m%dT%H%M%S)"
  step "backup of the committed figures/*.json to $BACKUP"
  run mkdir -p "$BACKUP"
  run cp -p "$F"/*.json "$BACKUP/"
  step "inputs first: theory numbers and the representative profiles (gendata/make_case.py reads figures/profiles_vmf.json)"
  run "$PY" "$A/theory_numbers.py"
  campaign() { step "campaign $1: $2 (np $3)"; run "$PY" "$S/run_campaign.py" --config "$1" --cases "$2" --np "$3" --skip-ok "${@:4}"; }
  bitwise() {   # $1 $2: run directories; fails if any state2d/state3d/diag_rho record differs byte for byte
    run "$PY" -c "import sys; sys.path.insert(0, '$ROOT'); from pfcommon import campaign as cp; bad = {p: cp.compare_bitwise('$1', '$2', p) for p in ('state2d.*.data', 'state3d.*.data', 'diag_rho.*.data')}; print(bad); sys.exit(1 if any(bad.values()) else 0)"
  }
  # the order actually used: the thesis bridge (Task 3), the fjord campaign (Task 6), then R2/R3
  campaign bay_thesis R1 4
  # build control (Task 3): 2 h of R1 on its own build and on the KPP build must be bit-identical
  campaign bay_thesis R1 4 --end-time 7200 --run-suffix _short2h
  campaign bay_thesis_kpp R1b 4 --end-time 7200 --run-suffix _short2h
  step "build control: R1b_short2h against R1_short2h, and R1_short2h against the stored one (bit for bit)"
  bitwise "$R/R1b_short2h" "$R/R1_short2h"
  bitwise "$R/R1_short2h" "$ORIG_RUNS/R1_short2h"
  campaign bay_thesis_kpp R1k 4
  for c in F1 F2 F3 F1c V3 F1g F5 F4; do campaign fjord_vmf "$c" 5; done
  campaign fjord_vmf_400 V4 5
  campaign bay_thesis_kpp 'R[23]' 4
  step "restart tests (PD1): F2 and F4, 10 + pickup + 10 steps against 20"
  run "$PY" "$S/restart_test.py" F2
  run "$PY" "$S/restart_test.py" F4
  step "V2a column test (PD2, post hoc): doubly periodic slab with the F1 configuration"
  run "$ROOT/scripts/build_mitgcm.sh" "$HERE/config/fjord_vmf_column" "$BUILDS/fjord_vmf_column"
  run "$PY" "$S/column_test.py" make --f1-input "$R/F1/input" --out "$R/V2a_column/input"
  run "$PY" "$S/column_test.py" run --input "$R/V2a_column/input" --run "$R/V2a_column"
fi

step "theory numbers, representative profiles (PD0) and the bathymetry report (Task 1b)"
[ "$QUICK" -eq 0 ] || run "$PY" "$A/theory_numbers.py"        # the full mode ran it before the campaigns
run "$PY" "$A/bathy_report.py"
step "PD1: timing record and restart tests"
PD1_RUNS=$("$PY" -c "import json, os; d = json.load(open('$F/pd1_timing.json')); print(' '.join(p for p in ('$R/' + os.path.basename(v['run_dir']) for v in d.values()) if os.path.isdir(p)))")
[ -z "$PD1_RUNS" ] || run "$PY" "$A/pd1_timing.py" $PD1_RUNS
run "$PY" "$S/restart_test.py" F2 --compare-only
run "$PY" "$S/restart_test.py" F4 --compare-only
step "V1: R1 against the saved thesis run; contrast of the flux and budget methods with the thesis scripts"
run "$PY" "$A/v1_thesis.py" --run "$R/R1"
run "$PY" "$A/thesis_contrast.py" reference-nash
run "$PY" "$A/thesis_contrast.py" reference-kf
run "$PY" "$A/thesis_contrast.py" reference-diss
run "$PY" "$A/thesis_contrast.py" contrast --run "$R/R1" --recompute
step "V2-V4 (validation.json), ice-open-boundary check, PD2 diagnostics (post hoc), Figure 1"
run "$PY" "$A/validate.py" evaluate --f1 "$R/F1" --v3 "$R/V3" --v4 "$R/V4"
run "$PY" "$A/validate.py" ice-check "$R/F2" "$R/F4"
run "$PY" "$S/column_test.py" analyse --run "$R/V2a_column"
run "$PY" "$A/pd2_v2_diagnostics.py" --f1 "$R/F1" --validation "$F/validation.json"
run "$PY" "$A/pd2_v2_diagnostics.py" --f1 "$R/F1" --validation "$F/validation.json" --free-wave            # -> pd2_v2b_free_wave.json
run "$PY" "$A/pd2_v4_diagnostics.py" --f1 "$R/F1" --v4 "$R/V4"
run "$PY" "$A/pd2_v4_diagnostics.py" --f1 "$R/F1" --v4 "$R/V4" --descriptive   # -> pd2_v4_descriptive.json
run "$PY" "$A/validate.py" pd2-merge
run "$PY" "$A/fig1_validation.py"
step "H4a-H4d: band partition and Kang & Fringer budgets of every run (partition.json, budget.json)"
run "$PY" "$A/partition.py" compute F1 F2 F3 F4 F5 F1c F1g V3 V4 R1 R1k R2 R3 --recompute --runs-root "$R" --cache "$R/analysis_cache"
run "$PY" "$A/partition.py" assemble --runs-root "$R" --cache "$R/analysis_cache"
run "$PY" "$A/hypotheses.py"
step "post hoc: moving-reference strait pressure work; Figures 2-4"
run "$PY" "$A/posthoc_reference.py" --runs-root "$R"
run "$PY" "$A/figures.py" --cache "$R/analysis_cache"

if [ -n "${VECTORCRAFT_CLI:-}" ]; then
  step "Figure 0 (VectorCraft)"; run "$PY" "$F/src/fig0_schematic.py" --cli "$VECTORCRAFT_CLI"
else
  echo; echo "=== Figure 0 kept as committed (set VECTORCRAFT_CLI to redraw it)"
fi

step "README"; run "$PY" "$S/render_readme.py"

if [ "$QUICK" -eq 0 ]; then
  step "compare the regenerated figures/*.json with the committed ones ($BACKUP)"
  run "$PY" "$S/compare_json.py" "$BACKUP" "$F"
fi
echo; echo "done: metrics and figures in exp4_arctic_bayiw/figures/, README rendered."
