#!/bin/bash
# Exp. 2, ruling R24 (approved by Esteban 2026-11-10): detached queue for the POST HOC exposure
# test of H1's n > 2 (notes/pd2_ensamble.md, section R24, written and committed before this
# queue starts). Five 4 km cases hol_T{-0.5,0.0,0.5,1.0,1.5}: Holland 2008-like ramp
# (z_c = 150 m, d = 100 m), everything else as the ensemble core (config ocean0_4km:
# GammaT_cal,4km 0.029356, 2 years, melt_law, checkIniTemp). Re-runnable (--skip-ok).
#
# Usage (detached):  nohup caffeinate -is exp2_isomip_cavity/scripts/queue_r24.sh > LOG 2>&1 < /dev/null &
# Rules: refuses to start if an MPI launcher is running or LaCie is not mounted; thread
# variables = 1; np = 1, at most 5 concurrent; run_campaign.py refuses to start while another
# model process exists and records the binary sha256 in each case.json.
set -u

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
PY=/opt/anaconda3/envs/pyocean/bin/python
RC="$REPO/exp2_isomip_cavity/scripts/run_campaign.py"
RUNS=/Volumes/LaCie/portfolio_runs/exp2
LOGS="$RUNS/logs/r24"
RS=_r21

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

say() { echo "[r24 $(date "+%Y-%m-%d %H:%M:%S")] $*" >&2; }

if ! mount | grep -q " on /Volumes/LaCie "; then say "LaCie not mounted, refusing"; exit 2; fi
BUSY=$(pgrep -fl mpirun)
if [ -n "$BUSY" ]; then say "an MPI launcher is running, refusing: $BUSY"; exit 1; fi
mkdir -p "$LOGS" || exit 1
cd "$REPO" || exit 1
say "start; repo $REPO at $(git rev-parse --short HEAD); pid $$"

N=$("$PY" "$RC" --config ocean0_4km --cases 'hol_*' --dry-run 2>&1 | grep -c '^hol_')
if [ "$N" != 5 ]; then say "hol_* selects $N cases, expected 5: refusing"; exit 1; fi

say "START runs"
"$PY" "$RC" --config ocean0_4km --cases 'hol_*' --np 1 --max-concurrent 5 \
  --run-suffix "$RS" --skip-ok > "$LOGS/q_runs.log" 2>&1
say "END   runs rc=$?"

LEFT=$(pgrep -fl mpirun)
say "QUEUE DONE; MPI launchers left: [${LEFT}]"
