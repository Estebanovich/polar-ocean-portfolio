#!/bin/bash
# Usage: build_mitgcm.sh <config_dir> <build_dir> [--mpi]
# Builds MITgcm with <config_dir>/code as -mods. Source tree is never modified.
set -euo pipefail
MITGCM="${MITGCM_ROOT:-/Volumes/ESTEBAN_MIT/MITgcm}"
if [ $# -lt 2 ]; then echo "usage: $0 <config_dir> <build_dir> [--mpi]" >&2; exit 1; fi
CONFIG=$(cd "$1" && pwd)
mkdir -p "$2"
BUILD=$(cd "$2" && pwd)
MPI=""
[ "${3:-}" = "--mpi" ] && MPI="-mpi"
cd "$BUILD"
"$MITGCM/tools/genmake2" -rootdir "$MITGCM" -mods "$CONFIG/code" \
  -optfile "$MITGCM/tools/build_options/darwin_arm64_gfortran" $MPI
make depend
make -j4
test -x mitgcmuv
echo "build OK: $BUILD/mitgcmuv"
