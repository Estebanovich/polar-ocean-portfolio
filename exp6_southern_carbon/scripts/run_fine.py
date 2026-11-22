#!/opt/anaconda3/envs/pyocean/bin/python
"""Fine-model campaign of Exp. 6 (plan Task 9, Step 5), resumable.

Usage::

    run_fine.py [--years 1900 2060] [--runs-root DIR] [--no-ft]

1. F0_<Y> for each year: interpolated start from GT (run_campaign), then stages of whole years
   (analysis.fine.f0_next_years: up to the 3-year minimum, then 1 year at a time) until V6d
   passes or the 8-year cap; the V6d result goes to <runs>/F0_<Y>/v6d.json.
2. The FT pairs, member by member (c1, w1, c2, w2, c3, w3; the order of cases.yaml), 1900 first,
   with --skip-ok, so a cut (power, disk) resumes where it stopped.
Completed F0 years are counted from the F0 run's own latest pickup (iteration 0 = 1 January of Y),
so a stage that failed mid-year is rerun from its last whole year."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
REPO = PKG.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp6_southern_carbon.analysis import fine  # noqa: E402
from exp6_southern_carbon.scripts import run_campaign as rc  # noqa: E402

SPY = rc.STEPS_G50                    # f10 has the same deltaT (1000 s): 31536 steps per year


def f0_years_done(run: Path) -> int:
    last = rc.latest_pickup(run) if run.exists() else None
    return int(last) // SPY if last else 0


def spin_up(case: str, runs_root: Path) -> int:
    run = runs_root / case
    while True:
        years = f0_years_done(run)
        series = fine.eke_series(run)[: 12 * years] if years else []
        v = fine.v6d(series)
        n = fine.f0_next_years(12 * years, v["pass"])
        if n == 0:
            (run / "v6d.json").write_text(json.dumps(dict(v, years=years, series=series), indent=1) + "\n")
            print(f"{case}: V6d {'PASS' if v['pass'] else 'FAIL'} after {years} years", flush=True)
            return 0
        print(f"{case}: {years} years done, running {n} more", flush=True)
        args = ["--cases", case, "--n-timesteps", str(n * SPY), "--runs-root", str(runs_root)]
        if years:
            args.append("--continue")
        if rc.main(args) != 0:
            return 1


def main(argv=None) -> int:
    from pfcommon.paths import EXP6_RUNS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--years", type=int, nargs="+", default=[1900, 2060])
    ap.add_argument("--runs-root", type=Path, default=EXP6_RUNS)
    ap.add_argument("--no-ft", action="store_true")
    a = ap.parse_args(argv)
    for y in a.years:
        if spin_up(f"F0_{y}", a.runs_root) != 0:
            return 1
    if a.no_ft:
        return 0
    for y in a.years:
        if rc.main(["--cases", f"FT_{y}_*", "--skip-ok", "--runs-root", str(a.runs_root)]) != 0:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
