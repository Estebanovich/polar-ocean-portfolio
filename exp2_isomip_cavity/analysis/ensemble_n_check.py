#!/opt/anaconda3/envs/pyocean/bin/python
r"""pd0 N1-N3 over the ensemble runs (Exp. 2 production queue, after the ensemble step).

Usage::

    ensemble_n_check.py [--runs-root /Volumes/LaCie/portfolio_runs/exp2] [--run-suffix _r21]
                        [--out LOG_DIR/ensemble_n_check.json]

For every ensemble member (``make_case.ensemble_cases``: the reference ``ocean0_4km`` and the
28 ``ens_*`` cases) it reports, mechanically against ``notes/pd0_preregistro.md``:

* N3: ``metrics.load_run`` succeeds (``require_ok``: OK, reached nTimeSteps, no NaN; and the
  R21 ``melt_law`` guard). A run that fails is listed as a **gap**; nothing is rerun;
* N1: ``metrics.equilibrium_drift(rd, '2yr')`` < 0.03;
* N2: ``|metrics.tracer_budget_residual(rd)|`` < 0.01.

It never stops the queue: the exit code is 0 whatever the verdicts (2 only for a usage
error). The summary is printed (one line per run, then totals) and written as JSON.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp2_isomip_cavity.analysis import metrics as mt  # noqa: E402
from exp2_isomip_cavity.gendata.make_case import ensemble_cases, load_cases  # noqa: E402
from pfcommon.paths import EXP2_RUNS as DEFAULT_RUNS  # noqa: E402

N1_MAX = 0.03          # pd0 N1
N2_MAX = 0.01          # pd0 N2


def check_run(run_dir: Path) -> dict:
    out = {"run_dir": str(run_dir)}
    try:
        rd = mt.load_run(run_dir)
    except Exception as e:                        # noqa: BLE001 - any failure is an N3 gap
        out.update(N3={"pass": False, "error": f"{type(e).__name__}: {e}"}, gap=True)
        return out
    out.update(N3={"pass": True}, gap=False)
    for key, fn, lim in (("N1", lambda: mt.equilibrium_drift(rd, "2yr"), N1_MAX),
                         ("N2", lambda: mt.tracer_budget_residual(rd), N2_MAX)):
        try:
            v = float(fn())
            out[key] = {"value": v, "pass": bool(abs(v) < lim)}
        except Exception as e:                    # noqa: BLE001 - reported, not fatal
            out[key] = {"value": None, "pass": False, "error": f"{type(e).__name__}: {e}"}
    return out


def check_ensemble(runs_root: Path, suffix: str, names: list[str] | None = None) -> dict:
    names = names if names is not None else ensemble_cases(load_cases())
    runs = {n: check_run(Path(runs_root) / f"{n}{suffix}") for n in names}
    return {
        "rule": "pd0 N1 (< 0.03), N2 (< 0.01), N3 (require_ok); gaps are not rerun",
        "runs": runs,
        "n_runs": len(runs),
        "gaps_N3": [n for n, r in runs.items() if r["gap"]],
        "fail_N1": [n for n, r in runs.items() if not r["gap"] and not r["N1"]["pass"]],
        "fail_N2": [n for n, r in runs.items() if not r["gap"] and not r["N2"]["pass"]],
    }


def _fmt(n: str, r: dict) -> str:
    if r["gap"]:
        return f"{n}: N3 FAIL (gap, not rerun): {r['N3']['error'][:160]}"
    v = lambda k: "n/a" if r[k]["value"] is None else f"{r[k]['value']:.4g}"   # noqa: E731
    ok = lambda k: "PASS" if r[k]["pass"] else "FAIL"                          # noqa: E731
    return f"{n}: N3 PASS, N1 {v('N1')} {ok('N1')}, N2 {v('N2')} {ok('N2')}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="pd0 N1-N3 over the ensemble runs")
    ap.add_argument("--runs-root", type=Path, default=DEFAULT_RUNS)
    ap.add_argument("--run-suffix", default="_r21")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    res = check_ensemble(a.runs_root, a.run_suffix)
    for n, r in res["runs"].items():
        print(_fmt(n, r))
    print(f"SUMMARY: {res['n_runs']} runs; N3 gaps {len(res['gaps_N3'])} {res['gaps_N3']}; "
          f"N1 fail {len(res['fail_N1'])} {res['fail_N1']}; "
          f"N2 fail {len(res['fail_N2'])} {res['fail_N2']}")
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(res, indent=2, default=float) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
