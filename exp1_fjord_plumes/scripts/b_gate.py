#!/opt/anaconda3/envs/pyocean/bin/python
"""Gate after the first case of the (b) campaign (ruling R30).

Usage::

    b_gate.py <run_dir> [--ref <ref_run_dir>]

Run it on the first campaign case (``b_3d_twolayer_x0_40``, 3 h, dt = 2 s)
before the remaining 13 cases are launched with ``run_campaign.py --skip-ok``.
PASS (exit 0) needs both:

(i) volume budget: ``|retained_last_hour| <= 0.5`` from
    ``analysis.b_checks.eta_retention`` (the open boundary passes at least half
    of the discharge in the last hour and does not drain more than 1.5 times
    it; notes/pd1_b_twoplume3d.md §6);
(ii) R29 output did not change the physics: at iteration 900 (t = 30 min at
    dt = 2 s) the case's ``icefrontDiag`` (15-min mean, 15-30 min) and its
    30-min snapshots THETA, SALT, UVEL, VVEL, ETAN must be byte-identical to
    ``icefrontDiag`` and the final state dumps T, S, U, V, Eta of the reference
    short run ``b_3d_twolayer_x0_40_short30m_dt2`` (same inputs and binary,
    made before R29). Snapshots are labelled with the iteration of the state
    they hold (pkg/diagnostics/diagnostics_write.F:62-72, suffix myIter-1), so
    they compare with the state dumps of the same iteration. W is not compared
    (WVEL may be diagnosed differently from the dumped W). Any difference, or
    a missing file, is a FAIL; the maximum absolute difference is reported.

FAIL (exit 1) prints every reason. Both runs must be OK (require_ok).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from exp1_fjord_plumes.analysis.b_checks import RETENTION_MAX, eta_retention  # noqa: E402
from pfcommon.mitgcm_io import read_mds  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402

REF = Path("/Volumes/LaCie/portfolio_runs/exp1/b_3d_twolayer_x0_40_short30m_dt2")
ITER = 900
PAIRS = (("icefrontDiag", "icefrontDiag"), ("THETA", "T"), ("SALT", "S"),
         ("UVEL", "U"), ("VVEL", "V"), ("ETAN", "Eta"))


def compare(run_dir: Path, ref_dir: Path, it: int = ITER) -> dict:
    """{run_name: {"identical": bool|None, "max_abs_diff": float|None, "missing": [...]}}."""
    out = {}
    for name, ref_name in PAIRS:
        a = Path(run_dir) / f"{name}.{it:010d}.data"
        b = Path(ref_dir) / f"{ref_name}.{it:010d}.data"
        missing = [str(p) for p in (a, b) if not p.exists()]
        if missing:
            out[name] = {"identical": None, "max_abs_diff": None, "missing": missing}
            continue
        same = a.read_bytes() == b.read_bytes()
        d = 0.0
        if not same:
            x = read_mds(Path(run_dir) / name, it)[0]
            y = read_mds(Path(ref_dir) / ref_name, it)[0]
            d = float(np.max(np.abs(x - y))) if x.shape == y.shape else float("inf")
        out[name] = {"identical": same, "max_abs_diff": d, "missing": []}
    return out


def gate(run_dir: Path, ref_dir: Path = REF) -> tuple[bool, list[str], dict]:
    run_dir, ref_dir = Path(run_dir), Path(ref_dir)
    reasons: list[str] = []
    require_ok(ref_dir)
    ret = eta_retention(run_dir)  # require_ok(run_dir) inside
    r = ret["retained_last_hour"]
    if not abs(r) <= RETENTION_MAX:
        reasons.append(f"(i) retained_last_hour = {r:.3f} outside [-{RETENTION_MAX}, "
                       f"{RETENTION_MAX}] (retained_cum = {ret['retained_cum']:.3f})")
    dt_run = json.loads((run_dir / "case.json").read_text())["params"]["deltaT"]
    dt_ref = json.loads((ref_dir / "case.json").read_text())["params"]["deltaT"]
    if dt_run != dt_ref:
        reasons.append(f"(ii) deltaT differs: run {dt_run}, reference {dt_ref}")
    ident = compare(run_dir, ref_dir)
    for name, v in ident.items():
        if v["missing"]:
            reasons.append(f"(ii) {name}: missing " + ", ".join(v["missing"]))
        elif not v["identical"]:
            reasons.append(f"(ii) {name} not byte-identical to the reference "
                           f"(max |diff| = {v['max_abs_diff']:.3e})")
    return not reasons, reasons, {"retention": ret, "identity": ident}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--ref", type=Path, default=REF)
    a = ap.parse_args(argv)
    ok, reasons, det = gate(a.run_dir, a.ref)
    print(json.dumps(det, indent=2, default=str))
    if ok:
        print(f"PASS: {a.run_dir.name} (retained_last_hour = "
              f"{det['retention']['retained_last_hour']:.3f}; R29 output byte-identical)")
        return 0
    print("FAIL: " + a.run_dir.name)
    for r in reasons:
        print("  - " + r)
    return 1


if __name__ == "__main__":
    sys.exit(main())
