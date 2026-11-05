#!/bin/bash
# Usage: run_mitgcm.sh <build_dir> <input_dir> <run_dir> <np>
# Copies (never links: LaCie is exFAT) mitgcmuv + input, runs, writes run_status.json.
if [ $# -ne 4 ]; then echo "usage: $0 <build_dir> <input_dir> <run_dir> <np>" >&2; exit 1; fi
BUILD=$1; INPUT=$2; RUN=$3; NP=$4

# Refuse before touching anything if the target volume is not mounted.
case "$RUN" in
  /Volumes/LaCie*)
    VOL=$(echo "$RUN" | cut -d/ -f1-3)
    if [ ! -d "$VOL" ] || [ "$VOL" != "/Volumes/LaCie" ]; then
      echo "LaCie not mounted ($VOL)" >&2; exit 2
    fi;;
esac

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
ulimit -n 6000
mkdir -p "$RUN" || exit 1
cp "$BUILD/mitgcmuv" "$RUN/" || exit 1
cp -R "$INPUT"/. "$RUN"/ || exit 1
cd "$RUN" || exit 1
rm -f STDOUT.* STDERR.* run_status.json

NTS=$(grep -i '^ *nTimeSteps *=' data | head -1 | sed -E 's/.*= *([0-9]+).*/\1/')
NTS=${NTS:-0}
echo "{\"state\": \"RUNNING\", \"last_iter\": 0, \"n_timesteps\": $NTS, \"wall_s\": 0.0}" > run_status.json

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
if grep -q "PROGRAM MAIN: Execution ended Normally" "$LOG" 2>/dev/null; then STATE=OK; else STATE=FAILED; fi
LAST=$(grep 'time_tsnumber' "$LOG" 2>/dev/null | tail -1 | sed -E 's/.*= *([0-9]+).*/\1/')
LAST=${LAST:-0}
echo "{\"state\": \"$STATE\", \"last_iter\": $LAST, \"n_timesteps\": $NTS, \"wall_s\": $WALL}" > run_status.json

LEFT=$(pgrep -fl "mpirun|mitgcmuv")
echo "pgrep mpirun|mitgcmuv after run: [${LEFT}]" | tee -a output.txt
[ "$STATE" = OK ]
