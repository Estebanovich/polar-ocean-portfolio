#!/bin/bash
# Exp. 2, Task 4 phase A: detached queue for the pd0 K1 Gamma_T calibration, relaunched from
# scratch after ruling R21 (2026-11-09: AD16 Eq. 27 is u*^2 = Cd (u^2 + u_tidal^2); every
# earlier run used the wrong law and is void, notes/pd1_calibracion.md). Only the pd0 sample
# Gamma_T in {0.015, 0.03, 0.06} is run; the void PD2 amendment (0.0075) is NOT reused.
#
# Every run of this queue carries the run-directory suffix $RS (= _r21), so the void runs
# (cal4_gT*, cal2_gT*, ocean0_4km, ocean0_2km without suffix) are never picked up by
# --skip-ok nor read by calibrate.py; they stay on disk as a historical record.
# Re-runnable: runs already OK under the suffix with the same deltaT/endTime/gammaT are
# skipped (--skip-ok), so a rerun only adds the missing cases.
#
#   1. cal4_gT{0.015,0.03,0.06}$RS: Ocean0 at 4 km, 2 years (CAL4_NP / CAL4_CONC below)
#   2. cal2_gT{0.015,0.03,0.06}$RS: Ocean0 at 2 km, 2 years, np = 5, one at a time
#   3. analysis/calibrate.py on the 3 points of each resolution -> figures/calibration_<res>_r21.json
#      If 30 m/yr is not bracketed (exit 3, no extrapolation) that resolution stops there:
#      no verification run, PD2.
#   4. verification runs ocean0_4km$RS and ocean0_2km$RS with --set gammaT=<GammaT_cal>
#   5. K1 verdict of each verification run (calibrate.py --verify). Only if it fails, the one
#      extra iteration pd0 K1 allows: calibrate.py on the 4 points (3 + verification run)
#      -> figures/calibration_<res>_iter2_r21.json, and one more run (<case>$RS_iter2). Its K1
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
LOGS="$RUNS/logs/r21"                     # R21 relaunch logs (the void ones stay in logs/)
RS=_r21                                   # R21 run-directory suffix (see header)
GAMMAS="0.015 0.03 0.06"                  # pd0 K1 sample, nothing else
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

# fnmatch glob that selects exactly the pd0 sample (0.015, 0.03, 0.06; not 0.0075):
# checked against the case registry before anything runs
SEL='gT0.0[136]*'
for p in cal4 cal2; do
  got=$("$PY" "$RC" --config "$([ $p = cal4 ] && echo "$CAL4_CONFIG" || echo ocean0_2km)" \
        --cases "${p}_$SEL" --dry-run 2>&1 | grep -o "${p}_gT[0-9.]*" | sort -u | tr '\n' ' ')
  want=$(for g in $GAMMAS; do echo "${p}_gT$g"; done | sort -u | tr '\n' ' ')
  if [ "$got" != "$want" ]; then say "case selection for $p is [$got], expected [$want]"; exit 1; fi
done

step cal4 "$PY" "$RC" --config "$CAL4_CONFIG" --cases "cal4_$SEL" --np "$CAL4_NP" \
  --max-concurrent "$CAL4_CONC" --run-suffix "$RS" --skip-ok
step cal2 "$PY" "$RC" --config ocean0_2km --cases "cal2_$SEL" --np 5 --run-suffix "$RS" \
  --skip-ok

dirs() { local p=$1; for g in $GAMMAS; do echo "$RUNS/${p}_gT$g$RS"; done; }

calibrate() {   # calibrate <label> <run dirs...>: prints GammaT_cal on success, else ""
  local label=$1; shift
  if step "calibrate_$label" "$PY" "$CAL" --res "${label%%_*}" \
       --out "$FIG/calibration_$label$RS.json" "$@"
  then tail -1 "$LOGS/q_calibrate_$label.log"; else echo ""; fi
}

verify_and_iterate() {  # <res> <config> <np> <GammaT_cal> <calibration run dirs...>
  local res=$1 cfg=$2 np=$3 g=$4; shift 4
  local ver="$RUNS/$cfg$RS"
  if [ -z "$g" ]; then
    say "no GammaT_cal at $res (not bracketed or failed, see q_calibrate_$res.log):" \
        "queue stops for $res, no verification, no extrapolation: PD2"; return
  fi
  step "ver_$res" "$PY" "$RC" --config "$cfg" --cases "$cfg" --np "$np" --max-concurrent 1 \
    --set "gammaT=$g" --run-suffix "$RS" --skip-ok \
    || { say "verification run $cfg$RS not OK: PD2"; return; }
  step "verify_$res" "$PY" "$CAL" --verify --res "$res" --out "$FIG/verify_$res$RS.json" "$ver"
  local rc=$?
  if [ $rc -eq 0 ]; then say "K1 at $res: PASS (deep melt $(tail -1 "$LOGS/q_verify_$res.log"))"
    return; fi
  if [ $rc -ne 4 ]; then say "K1 check at $res failed to run (rc=$rc)"; return; fi
  say "K1 at $res: FAIL (deep melt $(tail -1 "$LOGS/q_verify_$res.log")): the one extra iteration"
  local g2
  g2=$(calibrate "${res}_iter2" "$@" "$ver" | tail -1)
  if [ -z "$g2" ]; then say "iteration 2 at $res: no GammaT_cal, PD2"; return; fi
  step "ver_${res}_iter2" "$PY" "$RC" --config "$cfg" --cases "$cfg" --np "$np" \
    --max-concurrent 1 --set "gammaT=$g2" --run-suffix "${RS}_iter2" --skip-ok \
    || { say "iteration-2 run not OK: PD2"; return; }
  step "verify_${res}_iter2" "$PY" "$CAL" --verify --res "$res" \
    --out "$FIG/verify_${res}_iter2$RS.json" "${ver}_iter2"
  say "K1 at $res after iteration 2: rc=$? (0 PASS, 4 FAIL -> PD2)"
}

G4=$(calibrate 4km $(dirs cal4) | tail -1)
G2=$(calibrate 2km $(dirs cal2) | tail -1)
say "GammaT_cal: 4km=[${G4}] 2km=[${G2}]"

verify_and_iterate 4km ocean0_4km 1 "$G4" $(dirs cal4)
verify_and_iterate 2km ocean0_2km 5 "$G2" $(dirs cal2)

LEFT=$(pgrep -fl mpirun)
say "QUEUE DONE; MPI launchers left: [${LEFT}]"
