#!/opt/anaconda3/envs/pyocean/bin/python
"""Compare two sets of Exp. 4 metrics files (the committed figures/*.json against a regenerated set).

Usage: compare_json.py OLD NEW   (directories, or two files)

For every JSON file of OLD: identical (parsed values equal), differing (the paths of the values
that differ, with the largest relative difference |new - old| / |old| of the numbers), or missing
in NEW. Exit status 0 only if every file is identical. Used by reproduce.sh (full mode: the
backup of the committed JSONs against the regenerated ones).
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path


# files not compared, with the reason (reported, never silently skipped)
EXCLUDED = {"ice_ob_check.json": ("keyed by absolute run paths (validate.py ice-check stamps): a rerun under "
                                  "another run root adds keys; the check itself runs again in reproduce.sh "
                                  "and fails the script if the ice reaches the open boundaries")}


def diff(a, b, path: str = "") -> list[tuple[str, float | None]]:
    """[(path, relative difference or None for a non-numeric or structural difference)]."""
    num = (int, float)
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in list(a) + [k for k in b if k not in a]:
            p = f"{path}.{k}" if path else str(k)
            out += diff(a[k], b[k], p) if (k in a and k in b) else [(p, None)]
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [(path, None)]
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += diff(x, y, f"{path}[{i}]")
        return out
    if isinstance(a, num) and isinstance(b, num) and not isinstance(a, bool) and not isinstance(b, bool):
        if a == b or (isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b)):
            return []
        return [(path, abs(b - a) / abs(a) if a else math.inf)]
    return [] if a == b else [(path, None)]


def compare_dirs(old: Path, new: Path) -> dict:
    old, new = Path(old), Path(new)
    files = [old] if old.is_file() else sorted(old.glob("*.json"))
    rep = {"identical": [], "differ": {}, "missing_in_new": [], "excluded": {}}
    for f in files:
        if f.name in EXCLUDED:
            rep["excluded"][f.name] = EXCLUDED[f.name]
            continue
        g = new if new.is_file() else new / f.name
        if not g.exists():
            rep["missing_in_new"].append(f.name)
            continue
        d = diff(json.loads(f.read_text()), json.loads(g.read_text()))
        if not d:
            rep["identical"].append(f.name)
            continue
        rel = [r for _, r in d if r is not None]
        rep["differ"][f.name] = {"n": len(d), "max_rel": max(rel) if rel else None,
                                 "non_numeric": [p for p, r in d if r is None][:20],
                                 "first": [p for p, _ in d][:10]}
    return rep


def main(argv=None) -> int:
    a = sys.argv[1:] if argv is None else argv
    if len(a) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    rep = compare_dirs(Path(a[0]), Path(a[1]))
    print(f"identical: {len(rep['identical'])} files")
    for n, why in rep["excluded"].items():
        print(f"excluded: {n}: {why}")
    for n, d in rep["differ"].items():
        print(f"DIFFERS: {n}: {d['n']} values, max relative difference {d['max_rel']}, "
              f"non-numeric {d['non_numeric']}, first {d['first']}")
    for n in rep["missing_in_new"]:
        print(f"MISSING in new: {n}")
    return 0 if not rep["differ"] and not rep["missing_in_new"] else 1


if __name__ == "__main__":
    sys.exit(main())
