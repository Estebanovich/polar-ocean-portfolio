#!/bin/bash
# Exp. 2 production queue (PD1a approved by Esteban, 2026-11-09): the ensemble first, then
# Ocean1. Gamma_T comes from cases.yaml (GammaT_cal of the R21 calibration: 4 km 0.029356,
# 2 km 0.030959); no --set, no physics or configuration change. Every run directory carries
# the R21 suffix _r21; --skip-ok reuses a run only if deltaT, endTime, gammaT and melt_law
# match (so ocean0_4km_r21, the 4 km verification run, is reused as the ensemble reference).
#
#   1. ens: the 28 ens_* cases, Ocean0-like ramp profiles at 4 km, 2 years, np = 1,
#      5 at a time (run_campaign refuses to start if another model process exists)
#   2. ref: the ensemble reference ocean0_4km (expected: already OK, skipped)
#   3. ncheck: pd0 N1-N3 over the 29 ensemble runs (analysis/ensemble_n_check.py), summary
#      in the log. N1/N2 failures do not stop the queue; an N3 failure is listed as a gap
#      and is NOT rerun automatically
#   4. ocean1: ocean1_2km (AD16 Table 2: COLD initial state, WARM restoring), 20 years in
#      4 chained restart stages of 5 years, np = 5
#
# Usage (detached):  nohup exp2_isomip_cavity/scripts/queue_production.sh > LOG 2>&1 < /dev/null &
# Rules: refuses to start if an MPI launcher is running or LaCie is not mounted; thread
# variables = 1; every model run goes through run_campaign.py (which refuses to start while
# another model process exists, checks np / SIZE.h, and records the binary sha256 in each
# case.json).
set -u

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
PY=/opt/anaconda3/envs/pyocean/bin/python
RC="$REPO/exp2_isomip_cavity/scripts/run_campaign.py"
NCHK="$REPO/exp2_isomip_cavity/analysis/ensemble_n_check.py"
RUNS=/Volumes/LaCie/portfolio_runs/exp2
LOGS="$RUNS/logs/production"
RS=_r21

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

say() { echo "[queue $(date "+%Y-%m-%d %H:%M:%S")] $*" >&2; }

if ! mount | grep -q " on /Volumes/LaCie "; then say "LaCie not mounted, refusing"; exit 2; fi
BUSY=$(pgrep -fl mpirun)
if [ -n "$BUSY" ]; then say "an MPI launcher is running, refusing: $BUSY"; exit 1; fi
mkdir -p "$LOGS" || exit 1
cd "$REPO" || exit 1
say "start; repo $REPO at $(git rev-parse --short HEAD); pid $$"

step() {   # step <name> <command...>: run, log to $LOGS/q_<name>.log, report rc
  local name=$1; shift
  say "START $name: $*"
  "$@" > "$LOGS/q_$name.log" 2>&1
  local rc=$?
  say "END   $name rc=$rc"
  return $rc
}

# the selection must be exactly the 28 ensemble members (checked before anything runs)
N_ENS=$("$PY" "$RC" --config ocean0_4km --cases 'ens_*' --dry-run 2>&1 | grep -c '^ens_')
if [ "$N_ENS" != 28 ]; then say "ens_* selects $N_ENS cases, expected 28: refusing"; exit 1; fi

step ens "$PY" "$RC" --config ocean0_4km --cases 'ens_*' --np 1 --max-concurrent 5 \
  --run-suffix "$RS" --skip-ok \
  || say "ensemble step rc != 0: some case not OK (listed by ncheck below); continuing"
step ref "$PY" "$RC" --config ocean0_4km --cases ocean0_4km --np 1 --max-concurrent 1 \
  --run-suffix "$RS" --skip-ok || say "reference step rc != 0 (see q_ref.log); continuing"

step ncheck "$PY" "$NCHK" --run-suffix "$RS" --out "$LOGS/ensemble_n_check.json"
say "N1-N3: $(grep '^SUMMARY' "$LOGS/q_ncheck.log")"

step ocean1 "$PY" "$RC" --config ocean0_2km --cases ocean1_2km --np 5 --max-concurrent 1 \
  --stages 4 --run-suffix "$RS" --skip-ok || say "ocean1_2km not OK (see q_ocean1.log)"

LEFT=$(pgrep -fl mpirun)
say "QUEUE DONE; MPI launchers left: [${LEFT}]"
