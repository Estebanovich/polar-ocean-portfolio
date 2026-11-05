#!/bin/bash
# Usage: run_mitgcm.sh <build_dir> <input_dir> <run_dir> <np>
# Copies (never links: LaCie is exFAT) mitgcmuv + input, runs, writes run_status.json.
if [ $# -ne 4 ]; then echo "usage: $0 <build_dir> <input_dir> <run_dir> <np>" >&2; exit 1; fi
BUILD=$1; INPUT=$2; RUN=$3; NP=$4

# Refuse before touching anything if the target volume is not really mounted
# (a stale /Volumes/LaCie mountpoint directory is not enough).
case "$RUN" in
  /Volumes/LaCie*)
    VOL=$(echo "$RUN" | cut -d/ -f1-3)
    if ! mount | grep -q " on $VOL "; then
      echo "LaCie not mounted ($VOL)" >&2; exit 2
    fi;;
esac

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
ulimit -n 6000
mkdir -p "$RUN" || exit 1
cp "$BUILD/mitgcmuv" "$RUN/" || exit 1
cp -RL "$INPUT"/. "$RUN"/ || exit 1
cd "$RUN" || exit 1
rm -f STDOUT.* STDERR.* run_status.json

# Parse run length from `data` (comment lines ignored): nTimeSteps, else (endTime-startTime)/deltaT.
getpar() { grep -v '^ *#' data | grep -iE "^ *$1 *=" | head -1 | sed -E 's/^[^=]*= *//; s/[ ,].*$//; s/[dD]/e/'; }
NIT0=$(getpar nIter0); NIT0=${NIT0:-0}
NTS=$(getpar nTimeSteps); REASON=""
if ! [[ "$NTS" =~ ^[0-9]+$ ]]; then
  ET=$(getpar endTime); ST=$(getpar startTime); DT=$(getpar deltaT)
  NTS=$(awk -v e="$ET" -v s="${ST:-0}" -v d="$DT" 'BEGIN{ if (e==""||d==""||d+0<=0) {print ""; exit} printf "%d", (e-s)/d+0.5 }')
  if ! [[ "$NTS" =~ ^[0-9]+$ ]]; then NTS=0; REASON="cannot determine nTimeSteps from data"; fi
fi
[[ "$NIT0" =~ ^[0-9]+$ ]] || { NIT0=0; REASON="cannot parse nIter0"; }
printf '{"state": "RUNNING", "last_iter": %s, "n_timesteps": %s, "wall_s": 0.0}\n' "$NIT0" "$NTS" > run_status.json

T0=$(date +%s)
if [ "$NP" -gt 1 ]; then
  mpirun -np "$NP" ./mitgcmuv > output.txt 2> stderr.txt
else
  ./mitgcmuv > output.txt 2> stderr.txt
fi
T1=$(date +%s)
WALL=$((T1 - T0))

# Serial runs print the model log to stdout (no STDOUT.0000); MPI runs write STDOUT.0000.
LOG=STDOUT.0000; [ -s "$LOG" ] || LOG=output.txt
LAST=$(grep 'time_tsnumber' "$LOG" 2>/dev/null | tail -1 | sed -E 's/.*= *([0-9]+).*/\1/')
LAST=${LAST:-0}
STATE=FAILED
if [ -n "$REASON" ]; then :
elif ! grep -q "PROGRAM MAIN: Execution ended Normally" "$LOG" 2>/dev/null; then REASON="no normal-end line in model log"
elif [ "$LAST" -ne $((NIT0 + NTS)) ]; then REASON="last_iter $LAST != nIter0 + nTimeSteps ($((NIT0 + NTS)))"
else STATE=OK; fi
if [ -n "$REASON" ]; then
  printf '{"state": "%s", "last_iter": %s, "n_timesteps": %s, "wall_s": %s, "reason": "%s"}\n' "$STATE" "$LAST" "$NTS" "$WALL" "$REASON" > run_status.json
else
  printf '{"state": "%s", "last_iter": %s, "n_timesteps": %s, "wall_s": %s}\n' "$STATE" "$LAST" "$NTS" "$WALL" > run_status.json
fi

LEFT=$(pgrep -fl "mpirun|mitgcmuv")
echo "pgrep mpirun|mitgcmuv after run: [${LEFT}]" | tee -a output.txt
[ "$STATE" = OK ]
