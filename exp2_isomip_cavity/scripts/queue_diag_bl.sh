#!/bin/bash
# Exp. 2, ruling R20: detached queue for the POST HOC boundary-layer diagnostic
# (notes/pd1_diagnostico_capa_limite.md). Ocean0 4 km, GammaT = 0.015, 2 years, one change
# per case against cal4_gT0.015: diag4_adv33 (DST3), diag4_hfac025 (hFacMin 0.25). Variant
# (c) is not run (already the baseline; any other layer needs new code). Does not touch
# K1, V1 or GammaT_cal.
#
# Usage (detached):  nohup exp2_isomip_cavity/scripts/queue_diag_bl.sh > LOG 2>&1 &
# Rules: refuses to start if an MPI launcher is running or LaCie is not mounted; thread
# variables = 1; np = 1, at most 3 concurrent; run_campaign.py checks for model processes.
set -u

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
PY=/opt/anaconda3/envs/pyocean/bin/python
RC="$REPO/exp2_isomip_cavity/scripts/run_campaign.py"
CAL="$REPO/exp2_isomip_cavity/analysis/calibrate.py"
FIG="$REPO/exp2_isomip_cavity/figures"
RUNS=/Volumes/LaCie/portfolio_runs/exp2
LOGS="$RUNS/logs"

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

say() { echo "[diag $(date "+%Y-%m-%d %H:%M:%S")] $*" >&2; }

if ! mount | grep -q " on /Volumes/LaCie "; then say "LaCie not mounted, refusing"; exit 2; fi
BUSY=$(pgrep -fl mpirun)
if [ -n "$BUSY" ]; then say "an MPI launcher is running, refusing: $BUSY"; exit 1; fi
mkdir -p "$LOGS" "$FIG" || exit 1
cd "$REPO" || exit 1
say "start; repo $REPO at $(git rev-parse --short HEAD); pid $$"

step() {   # step <name> <command...>: run, log to $LOGS/d_<name>.log, report rc
  local name=$1; shift
  say "START $name: $*"
  "$@" > "$LOGS/d_$name.log" 2>&1
  local rc=$?
  say "END   $name rc=$rc"
  return $rc
}

step runs "$PY" "$RC" --config ocean0_4km --cases 'diag4_*' --np 1 --max-concurrent 3 --skip-ok
# deep melt of the last 6 months of each case (the K1 helper; exit 4 only means "not 28-32")
for c in cal4_gT0.015 diag4_adv33 diag4_hfac025; do
  step "melt_$c" "$PY" "$CAL" --verify --res 4km --out "$FIG/diag_bl_$c.json" "$RUNS/$c"
  say "$c deep melt (months 19-24): $(tail -1 "$LOGS/d_melt_$c.log")"
done

LEFT=$(pgrep -fl mpirun)
say "QUEUE DONE; MPI launchers left: [${LEFT}]"
