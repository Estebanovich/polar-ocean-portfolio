#!/usr/bin/env bash
# Run pytest under the anemoi interpreter, which has no pytest installed.
# A throw-away shim directory symlinks the pure-Python pytest stack from the
# pyocean env; nothing is installed. Default target: exp3_melt_emulator.
set -euo pipefail

PYOCEAN=/opt/anaconda3/envs/pyocean/bin/python
ANEMOI=/opt/anaconda3/envs/anemoi/bin/python

SITE=$("$PYOCEAN" -c "import pytest,os;print(os.path.dirname(os.path.dirname(pytest.__file__)))")
SHIM=$(mktemp -d)
trap 'rm -rf "$SHIM"' EXIT

for name in pytest _pytest pluggy iniconfig py.py; do
    ln -s "$SITE/$name" "$SHIM/$name"
done

cd "$(dirname "$0")/../.."
if [ "$#" -eq 0 ]; then
    set -- exp3_melt_emulator
fi
PYTHONPATH="$SHIM" "$ANEMOI" -m pytest -p no:cacheprovider "$@"
