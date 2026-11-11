#!/usr/bin/env bash
# Re-run a command that rewrites a committed Exp. 3 record, compare the result with the
# committed record (analysis/compare_records.py, session metadata ignored) and restore the
# committed file, whatever happens. Used by reproduce.sh (full mode); it never changes a
# committed record.
#
#   rerun_compare_restore.sh cv|profile RECORD -- COMMAND [ARGS...]
# Exit status: that of the command if it failed, else 0 if identical, 1 if different.
# The re-run record is kept for inspection as $RERUN_DIR/<record name>.rerun
# (RERUN_DIR defaults to $TMPDIR), never inside the repository.
set -euo pipefail

[ "$#" -ge 4 ] && [ "$3" = "--" ] || { echo "usage: $0 cv|profile RECORD -- COMMAND [ARGS...]" >&2; exit 2; }
KIND="$1"; RECORD="$2"; shift 3
PY="${PY:-/opt/anaconda3/envs/pyocean/bin/python}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ -f "$RECORD" ] || { echo "$RECORD not found" >&2; exit 2; }

BACKUP="$(mktemp "${TMPDIR:-/tmp}/exp3_record.XXXXXX")"
cp -p "$RECORD" "$BACKUP"
restore() { cp -p "$BACKUP" "$RECORD"; rm -f "$BACKUP"; }
trap restore EXIT

status=0
"$@" || status=$?
if [ "$status" -ne 0 ]; then
  echo "re-run failed (exit $status); committed $RECORD restored" >&2
  exit "$status"
fi
RERUN="${RERUN_DIR:-${TMPDIR:-/tmp}}/$(basename "$RECORD").rerun"
cp -p "$RECORD" "$RERUN"
"$PY" "$HERE/../analysis/compare_records.py" "$KIND" "$BACKUP" "$RERUN" || status=$?
exit "$status"
