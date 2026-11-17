#!/bin/bash
# Usage: build_mitgcm_tap.sh <config_dir> <build_dir> <adj|tlm|fwd> [--mpi]
# Builds MITgcm with <config_dir>/code_tap as -mods through genmake2 -tap (Tapenade, MITgcm-AD v2;
# Gaikwad et al. 2024 Sec. 4). Targets: adj -> mitgcmuv_tap_adj, tlm -> mitgcmuv_tap_tlm,
# fwd -> mitgcmuv (same code, no differentiation). Source tree is never modified.
# DRYRUN=1 prints the commands instead of running them.
set -euo pipefail
MITGCM="${MITGCM_ROOT:-/Volumes/ESTEBAN_MIT/MITgcm}"
TAP="${TAPENADE_HOME:-/Volumes/ESTEBAN_MIT/tools/tapenade-e89660a5921a1ccf68b4afc331b3b0910878a179}"
if [ $# -lt 3 ]; then echo "usage: $0 <config_dir> <build_dir> <adj|tlm|fwd> [--mpi]" >&2; exit 1; fi
case "$3" in
  adj) MK="make tap_adj"; EXE=mitgcmuv_tap_adj;;
  tlm) MK="make tap_tlm"; EXE=mitgcmuv_tap_tlm;;
  fwd) MK="make -j4"; EXE=mitgcmuv;;
  *) echo "unknown target '$3' (adj|tlm|fwd)" >&2; exit 1;;
esac
CONFIG=$(cd "$1" && pwd)
[ -d "$CONFIG/code_tap" ] || { echo "no $CONFIG/code_tap" >&2; exit 1; }
MPI=""
[ "${4:-}" = "--mpi" ] && MPI="-mpi"
if [ -z "${JAVA_HOME:-}" ] && [ -x /usr/libexec/java_home ]; then
  JAVA_HOME=$(/usr/libexec/java_home -v 11 2>/dev/null || true)
fi
GEN="$MITGCM/tools/genmake2 -rootdir $MITGCM -mods $CONFIG/code_tap -optfile $MITGCM/tools/build_options/darwin_arm64_gfortran -tap $MPI"
if [ "${DRYRUN:-0}" = "1" ]; then
  echo "export TAPENADE_HOME=$TAP JAVA_HOME=${JAVA_HOME:-}"
  echo "export PATH=$TAP/bin:\$PATH"
  echo "cd $2"
  echo "$GEN"
  echo "make depend"
  echo "$MK"
  echo "test -x $EXE"
  exit 0
fi
mkdir -p "$2"
BUILD=$(cd "$2" && pwd)
export TAPENADE_HOME="$TAP" JAVA_HOME
export PATH="$TAP/bin:$PATH"
cd "$BUILD"
$GEN
make depend
$MK
test -x "$EXE"
echo "build OK: $BUILD/$EXE"
