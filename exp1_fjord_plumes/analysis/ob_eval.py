"""Open-boundary criteria 1-3 for (b) runs (pre-registered, notes/pd1_b_twoplume3d.md, 2026-11-06).

1. require_ok (run finished, reached nTimeSteps, no non-finite monitor lines);
2. at every UVEL snapshot (t > 0), max|u| on the OB columns <= 1.5 x max|u| on the interior columns;
3. |domain-mean eta| < 1e-3 m at every monitor line, and the OB transport (obc_E_uVel_Int) equals the
   expected outflow Q (Q0_total from case.json) within 1e-6 Q (1e-6 m3/s if Q = 0).
Columns are 0-based Python slices: 10 m runs use OB 96:99 and interior 10:90; 5 m runs use OB 196:199 and
interior 20:180 (the same 100-900 m).
Usage: ob_eval.py RUN_DIR [--ob A:B] [--interior A:B] [--field UVEL]   (exit code 1 if any criterion fails)
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pfcommon.mitgcm_io import read_mds  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402

RATIO_MAX, ETA_MAX = 1.5, 1e-3


def _slice(s: str) -> slice:
    return slice(*map(int, s.split(":")))


def evaluate(run: Path, ob: slice, inn: slice, field: str = "UVEL") -> dict:
    Q = float(json.loads((run / "case.json").read_text())["Q0_total"])
    res = {"run": run.name, "Q_expected": Q}
    try:
        require_ok(run)
        res["c1"] = {"pass": True}
    except RuntimeError as e:
        res["c1"] = {"pass": False, "error": str(e)[:200]}
    worst = 0.0
    for p in sorted(glob.glob(str(run / f"{field}.*.meta"))):
        it = int(p.split(".")[-2])
        if it == 0:
            continue
        a = read_mds(run / field, it)[0][0]
        if not np.isfinite(a).all():
            worst = np.inf
            break
        worst = max(worst, float(np.abs(a[:, :, ob]).max() / max(np.abs(a[:, :, inn]).max(), 1e-12)))
    res["c2"] = {"pass": bool(worst <= RATIO_MAX), "max_ratio": worst}
    log = (run / "STDOUT.0000").read_text(errors="replace").splitlines()

    def val(k):
        return [float(line.split()[-1]) for line in log if f"%MON {k} " in line and "NaN" not in line]

    eta, tr = val("dynstat_eta_mean"), val("obc_E_uVel_Int")
    tol = 1e-6 * abs(Q) if Q else 1e-6
    eta_max = max(map(abs, eta)) if eta else float("nan")
    tr_err = max(abs(x - Q) for x in tr) if tr else float("nan")
    res["c3"] = {"pass": bool(eta and tr and eta_max < ETA_MAX and tr_err <= tol),
                 "max_abs_eta_mean_m": eta_max, "max_transport_error_m3_s": tr_err}
    res["overall"] = "PASS" if all(res[c]["pass"] for c in ("c1", "c2", "c3")) else "FAIL"
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--ob", default="196:199")
    ap.add_argument("--interior", default="20:180")
    ap.add_argument("--field", default="UVEL")
    a = ap.parse_args()
    r = evaluate(a.run_dir, _slice(a.ob), _slice(a.interior), a.field)
    print(json.dumps(r, indent=1, default=float))
    sys.exit(0 if r["overall"] == "PASS" else 1)
