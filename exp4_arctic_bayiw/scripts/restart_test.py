#!/opt/anaconda3/envs/pyocean/bin/python
"""Restart test of the fjord_vmf build (Task 4, PD1): one 20-step run against 10 + pickup + 10.

Usage: restart_test.py CASE [--compare-only] [--out figures/pd1_restart.json]

The inputs of CASE (cases.yaml) are generated for 2400 s (20 steps of 120 s) with pChkptFreq =
1200 s, so the permanent pickups (model, pkg/seaice, Orlanski) are written at step 10, which is
also a monitor time (run_mitgcm.sh reads the last iteration from the monitor). The single run
and the two-stage run (pfcommon.campaign.run_stages, equal split, pickups copied) go to
<EXP4_RUNS>/<CASE>_restarttest_{single,split}. A restart is exact if the step-20 pickups and the
step-10 diagnostics of the second stage equal those of the single run byte for byte
(pfcommon.campaign.compare_bitwise). Rules: np = 5, one run at a time, no other MPI process.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
REPO = PKG.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pfcommon import campaign as cp  # noqa: E402
from pfcommon.paths import EXP4_BUILDS, EXP4_RUNS  # noqa: E402
from exp4_arctic_bayiw.gendata import make_case as mc  # noqa: E402

PATTERNS = ("pickup.0000000020.data", "pickup_seaice.0000000020.data",
            "pickup_orlanski*.0000000020.*.data", "state3d.0000000010.data",
            "surf.0000000010.data", "seaice.0000000010.data")


def compare(case: str, runs: Path = EXP4_RUNS) -> dict:
    a, b = runs / f"{case}_restarttest_single", runs / f"{case}_restarttest_split" / "stage2"
    out = {p: cp.compare_bitwise(a, b, p) for p in PATTERNS}
    out["bit_identical"] = all(v == [] for v in out.values())
    return out


def run(case: str, runs: Path = EXP4_RUNS, builds: Path = EXP4_BUILDS) -> dict:
    if cp.other_mitgcm_processes():
        raise RuntimeError("another MPI / model process is running")
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}

    def sh(cmd):
        return subprocess.run([str(c) for c in cmd], env=env).returncode

    inp = runs / f"{case}_restarttest_input"
    mc.make_case(case, inp, overrides={"t_end_h": 2400.0 / 3600.0})
    t = (inp / "data").read_text()
    (inp / "data").write_text(re.sub(r"pChkptFreq=[^,]+,", "pChkptFreq=1200.,", t))
    b = builds / "fjord_vmf"
    run_sh = REPO / "scripts" / "run_mitgcm.sh"
    s1 = cp.run_stages(runs / f"{case}_restarttest_single", inp, b, 5, 1, sh, run_sh)
    s2 = cp.run_stages(runs / f"{case}_restarttest_split", inp, b, 5, 2, sh, run_sh,
                       pickup_stems=("pickup", "pickup_seaice"))
    if cp.other_mitgcm_processes():
        raise RuntimeError("processes left after the restart test")
    return {"single": s1, "split": s2}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case")
    ap.add_argument("--compare-only", action="store_true")
    ap.add_argument("--out", type=Path, default=PKG / "figures" / "pd1_restart.json")
    a = ap.parse_args(argv)
    res = {} if a.compare_only else run(a.case)
    res.update(compare(a.case))
    allres = json.loads(a.out.read_text()) if a.out.exists() else {}
    allres[a.case] = res
    a.out.write_text(json.dumps(allres, indent=2))
    print(json.dumps(res))
    return 0 if res["bit_identical"] else 1


if __name__ == "__main__":
    sys.exit(main())
