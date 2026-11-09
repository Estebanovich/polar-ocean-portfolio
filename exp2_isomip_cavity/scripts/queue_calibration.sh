#!/bin/bash
# Exp. 2, Task 4 phase A: detached queue for the pd0 K1 Gamma_T calibration.
#
#   1. cal4_gT{0.015,0.03,0.06}: Ocean0 at 4 km, 2 years       (CAL4_NP / CAL4_CONC below)
#   2. cal2_gT{0.015,0.03,0.06}: Ocean0 at 2 km, 2 years, np = 5, one at a time
#   3. analysis/calibrate.py for each resolution -> figures/calibration_<res>.json, GammaT_cal
#   4. verification runs ocean0_4km and ocean0_2km with --set gammaT=<GammaT_cal>
#      (skipped for a resolution whose calibration did not give a value: PD2)
#
# The K1/V1/N1-N3 verdicts are NOT computed here: they are evaluated afterwards, mechanically
# against notes/pd0_preregistro.md. A second K1 iteration, if pd0 allows it, is not queued.
#
# Usage (detached):  nohup exp2_isomip_cavity/scripts/queue_calibration.sh > LOG 2>&1 &
# Rules: refuses to start if an MPI launcher is running or LaCie is not mounted; thread
# variables = 1; every model run goes through run_campaign.py (which itself refuses to start
# while another model process exists and checks the np / SIZE.h layout).
set -u

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
PY=/opt/anaconda3/envs/pyocean/bin/python
RC="$REPO/exp2_isomip_cavity/scripts/run_campaign.py"
CAL="$REPO/exp2_isomip_cavity/analysis/calibrate.py"
FIG="$REPO/exp2_isomip_cavity/figures"
RUNS=/Volumes/LaCie/portfolio_runs/exp2
LOGS="$RUNS/logs"
# 4 km layout of the calibration batch, chosen from the timing runs (notes/pd1_calibracion.md)
CAL4_CONFIG=${CAL4_CONFIG:-ocean0_4km}
CAL4_NP=${CAL4_NP:-1}
CAL4_CONC=${CAL4_CONC:-3}

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

say() { echo "[queue $(date '+%Y-%m-%d %H:%M:%S')] $*"; }

if ! mount | grep -q " on /Volumes/LaCie "; then say "LaCie not mounted, refusing"; exit 2; fi
BUSY=$(pgrep -fl mpirun)
if [ -n "$BUSY" ]; then say "an MPI launcher is running, refusing: $BUSY"; exit 1; fi
mkdir -p "$LOGS" "$FIG" || exit 1
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

step cal4 "$PY" "$RC" --config "$CAL4_CONFIG" --cases 'cal4_gT*' --np "$CAL4_NP" \
  --max-concurrent "$CAL4_CONC" --skip-ok
step cal2 "$PY" "$RC" --config ocean0_2km --cases 'cal2_gT*' --np 5 --skip-ok

calibrate() {   # calibrate <res> <run dirs...>: prints GammaT_cal on success
  local res=$1; shift
  if step "calibrate_$res" "$PY" "$CAL" --res "$res" --out "$FIG/calibration_$res.json" "$@"
  then tail -1 "$LOGS/q_calibrate_$res.log"; else echo ""; fi
}
G4=$(calibrate 4km "$RUNS"/cal4_gT0.015 "$RUNS"/cal4_gT0.03 "$RUNS"/cal4_gT0.06 | tail -1)
G2=$(calibrate 2km "$RUNS"/cal2_gT0.015 "$RUNS"/cal2_gT0.03 "$RUNS"/cal2_gT0.06 | tail -1)
say "GammaT_cal: 4km=[${G4}] 2km=[${G2}]"

if [ -n "$G4" ]; then
  step ver4 "$PY" "$RC" --config ocean0_4km --cases ocean0_4km --np 1 --max-concurrent 1 \
    --set "gammaT=$G4" --skip-ok
else say "no GammaT_cal at 4 km (see q_calibrate_4km.log): verification skipped, PD2"; fi
if [ -n "$G2" ]; then
  step ver2 "$PY" "$RC" --config ocean0_2km --cases ocean0_2km --np 5 \
    --set "gammaT=$G2" --skip-ok
else say "no GammaT_cal at 2 km (see q_calibrate_2km.log): verification skipped, PD2"; fi

LEFT=$(pgrep -fl mpirun)
say "QUEUE DONE; MPI launchers left: [${LEFT}]"
