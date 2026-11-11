#!/opt/anaconda3/envs/pyocean/bin/python
r"""Compare a re-run record of Exp. 3 with the committed one, ignoring session metadata.

Used by reproduce.sh (full mode): ``run_protocol.py --cv`` and ``profile_test_diag.py``
rewrite figures/cv_results.json and figures/profile_test_diag.json with new timestamps,
wall times and session counters. The committed files are the record of the runs actually
made, so reproduce.sh backs them up, re-runs, compares with this script and restores them.

* ``cv``: the per-run ``entries`` must be identical, and the ``units`` identical once the
  timing fields (started, finished, fit_s, wall_s, amendment) are removed.
* ``profile``: the whole record must be identical once ``fit_s`` is removed.

Usage::

    compare_records.py cv|profile COMMITTED RERUN      (exit 0 identical, 1 different)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

TIMING_KEYS = frozenset({"started", "finished", "fit_s", "wall_s", "amendment"})


def strip(obj, ignore=TIMING_KEYS):
    """obj without the keys in ``ignore`` at any depth."""
    if isinstance(obj, dict):
        return {k: strip(v, ignore) for k, v in obj.items() if k not in ignore}
    if isinstance(obj, list):
        return [strip(v, ignore) for v in obj]
    return obj


def compare_cv(old: dict, new: dict) -> list[str]:
    """Differences between two cv_results.json records (empty list = same results)."""
    diffs = []
    oe = {(e["variant"], e["fold"], e["method"], e["case"]): e for e in old["entries"]}
    ne = {(e["variant"], e["fold"], e["method"], e["case"]): e for e in new["entries"]}
    if set(oe) != set(ne):
        diffs.append(f"entries: {len(set(oe) - set(ne))} missing, {len(set(ne) - set(oe))} new")
    diffs += [f"entry {'/'.join(k)} differs" for k in sorted(set(oe) & set(ne)) if oe[k] != ne[k]]
    ou = {u["unit"]: strip(u) for u in old["units"]}
    nu = {u["unit"]: strip(u) for u in new["units"]}
    if set(ou) != set(nu):
        diffs.append(f"units: {len(set(ou) - set(nu))} missing, {len(set(nu) - set(ou))} new")
    diffs += [f"unit {k} differs" for k in sorted(set(ou) & set(nu)) if ou[k] != nu[k]]
    return diffs


def compare_profile(old: dict, new: dict) -> list[str]:
    o, n = strip(old, {"fit_s"}), strip(new, {"fit_s"})
    return [f"{k} differs" for k in sorted(set(o) | set(n)) if o.get(k) != n.get(k)]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("kind", choices=("cv", "profile"))
    ap.add_argument("committed", type=Path)
    ap.add_argument("rerun", type=Path)
    a = ap.parse_args(argv)
    old, new = (json.loads(p.read_text()) for p in (a.committed, a.rerun))
    diffs = (compare_cv if a.kind == "cv" else compare_profile)(old, new)
    if diffs:
        print(f"{a.kind}: re-run DIFFERS from the committed record ({len(diffs)}): "
              + "; ".join(diffs[:10]) + (" …" if len(diffs) > 10 else ""))
        return 1
    print(f"{a.kind}: re-run identical to the committed record (session metadata ignored)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
