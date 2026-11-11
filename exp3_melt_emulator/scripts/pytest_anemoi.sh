#!/usr/bin/env bash
# Run pytest under the anemoi interpreter, which has no pytest installed.
# A throw-away shim directory symlinks the pure-Python pytest stack from the
# pyocean env; nothing is installed. Default target: exp3_melt_emulator, added
# whenever no path is given (also when only options are given, e.g. -q).
# Interpreters: $PYOCEAN (has pytest) and $ANEMOI (runs the tests); defaults below.
set -euo pipefail

PYOCEAN="${PYOCEAN:-/opt/anaconda3/envs/pyocean/bin/python}"
ANEMOI="${ANEMOI:-/opt/anaconda3/envs/anemoi/bin/python}"

SITE=$("$PYOCEAN" -c "import pytest,os;print(os.path.dirname(os.path.dirname(pytest.__file__)))")
SHIM=$(mktemp -d)
trap 'rm -rf "$SHIM"' EXIT

for name in pytest _pytest pluggy iniconfig py.py; do
    ln -s "$SITE/$name" "$SHIM/$name"
done

cd "$(dirname "$0")/../.."
# Add the default target unless a path (file, directory or node id file::test) is given;
# option values such as "-k expr" are not mistaken for paths because they do not exist.
has_path=0
for arg in "$@"; do
    case "$arg" in
        -*) ;;
        *) [ -e "${arg%%::*}" ] && has_path=1 ;;
    esac
done
if [ "$has_path" -eq 0 ]; then
    set -- "$@" exp3_melt_emulator
fi
PYTHONPATH="$SHIM" "$ANEMOI" -m pytest -p no:cacheprovider "$@"
