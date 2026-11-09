#!/bin/bash
# Exp. 2, Task 4 phase A: detached queue for the pd0 K1 Gamma_T calibration, with the PD2
# sampling amendment (Gamma_T = 0.0075 added at both resolutions; Esteban, 2026-11-09,
# notes/pd1_calibracion.md). Re-runnable: runs already OK with the same deltaT/endTime/gammaT
# are skipped (--skip-ok), so a rerun only adds the missing cases.
#
#   1. cal4_gT{0.0075,0.015,0.03,0.06}: Ocean0 at 4 km, 2 years (CAL4_NP / CAL4_CONC below)
#   2. cal2_gT{0.0075,0.015,0.03,0.06}: Ocean0 at 2 km, 2 years, np = 5, one at a time
#   3. analysis/calibrate.py on all points of each resolution -> figures/calibration_<res>.json
#   4. verification runs ocean0_4km and ocean0_2km with --set gammaT=<GammaT_cal>
#      (skipped for a resolution whose calibration did not give a value: PD2)
#   5. K1 verdict of each verification run (calibrate.py --verify). Only if it fails, the one
#      extra iteration pd0 K1 allows: calibrate.py on all points + the verification run
#      -> figures/calibration_<res>_iter2.json, and one more run (<case>_iter2). Its K1
#      verdict is written but nothing further is queued (a second failure goes to PD2).
#
# The K1/V1/N1-N3 evaluation is done afterwards, mechanically against
# notes/pd0_preregistro.md.
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

say() { echo "[queue $(date "+%Y-%m-%d %H:%M:%S")] $*" >&2; }

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

GAMMAS="0.0075 0.015 0.03 0.06"
dirs() { local p=$1; for g in $GAMMAS; do echo "$RUNS/${p}_gT$g"; done; }

calibrate() {   # calibrate <label> <run dirs...>: prints GammaT_cal on success, else ""
  local label=$1; shift
  if step "calibrate_$label" "$PY" "$CAL" --res "${label%%_*}" \
       --out "$FIG/calibration_$label.json" "$@"
  then tail -1 "$LOGS/q_calibrate_$label.log"; else echo ""; fi
}

verify_and_iterate() {  # <res> <config> <np> <GammaT_cal> <calibration run dirs...>
  local res=$1 cfg=$2 np=$3 g=$4; shift 4
  local ver="$RUNS/$cfg"
  if [ -z "$g" ]; then
    say "no GammaT_cal at $res (see q_calibrate_$res.log): verification skipped, PD2"; return
  fi
  step "ver_$res" "$PY" "$RC" --config "$cfg" --cases "$cfg" --np "$np" --max-concurrent 1 \
    --set "gammaT=$g" --skip-ok || { say "verification run $cfg not OK: PD2"; return; }
  step "verify_$res" "$PY" "$CAL" --verify --res "$res" --out "$FIG/verify_$res.json" "$ver"
  local rc=$?
  if [ $rc -eq 0 ]; then say "K1 at $res: PASS (deep melt $(tail -1 "$LOGS/q_verify_$res.log"))"
    return; fi
  if [ $rc -ne 4 ]; then say "K1 check at $res failed to run (rc=$rc)"; return; fi
  say "K1 at $res: FAIL (deep melt $(tail -1 "$LOGS/q_verify_$res.log")): the one extra iteration"
  local g2
  g2=$(calibrate "${res}_iter2" "$@" "$ver" | tail -1)
  if [ -z "$g2" ]; then say "iteration 2 at $res: no GammaT_cal, PD2"; return; fi
  step "ver_${res}_iter2" "$PY" "$RC" --config "$cfg" --cases "$cfg" --np "$np" \
    --max-concurrent 1 --set "gammaT=$g2" --run-suffix _iter2 --skip-ok \
    || { say "iteration-2 run not OK: PD2"; return; }
  step "verify_${res}_iter2" "$PY" "$CAL" --verify --res "$res" \
    --out "$FIG/verify_${res}_iter2.json" "${ver}_iter2"
  say "K1 at $res after iteration 2: rc=$? (0 PASS, 4 FAIL -> PD2)"
}

G4=$(calibrate 4km $(dirs cal4) | tail -1)
G2=$(calibrate 2km $(dirs cal2) | tail -1)
say "GammaT_cal: 4km=[${G4}] 2km=[${G2}]"

verify_and_iterate 4km ocean0_4km 1 "$G4" $(dirs cal4)
verify_and_iterate 2km ocean0_2km 5 "$G2" $(dirs cal2)

LEFT=$(pgrep -fl mpirun)
say "QUEUE DONE; MPI launchers left: [${LEFT}]"
