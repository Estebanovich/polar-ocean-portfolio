#!/opt/anaconda3/envs/pyocean/bin/python
"""PD2 column test of the model's inertial frequency (POST HOC; Esteban, PD2 2026-11-13).

Why: V2a failed (F1 shelf near-inertial peak at 1.136 f, 3.1 df); the peak is blue-shifted over
the whole domain (figures/pd2_v2_diagnostics.json). This test separates (a) a model error in the
inertial frequency from (c) the bounded domain.

Rule fixed BEFORE running (PD2 decision): peak <= 1 df of f -> the model's inertial frequency is
right and the F1 blue shift is (c); peak > 1 df -> (a), investigate. It changes the
classification of V2a only, never its verdict.

Configuration: F1's own input namelists (data, data.free, data.kpp, data.diagnostics, eedata:
Delta t 120 s staggered, KPP, PD0 mixing, fitted linear EOS, f0 of 78 N, same output interval and
forced + free stages to 280 h) on a horizontally uniform, DOUBLY PERIODIC 10 x 10 x 40 domain:
no land and useOBCS = .FALSE., so MITgcm's default periodic boundaries apply in x and y (a
domain with no walls; model/src/ini_masks_etc.F masks come from the bathymetry only). Flat bottom
at the depth of the V2a shelf station; T and S the F1 initial column at that station; the wind
stress the F1 records at that station (spatially uniform in F1 too), spread uniformly. The build
is config/fjord_vmf_column: the fjord_vmf code (same MITgcm mods) with a 10 x 10, 1-rank SIZE.h.
Two stages as F1 (pfcommon.campaign.run_stages, forced data then data.free).

Usage:
    column_test.py make  --f1-input <F1>/input --out <run>/input [--depth 153.37 --ji 114 62]
    column_test.py run   --input <run>/input --run <run> [--allow-concurrent]   (np = 1; builds nothing)
    column_test.py analyse --run <run> [--out figures/pd2_column_test.json]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parents[1]
REPO = PKG.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pfcommon import campaign as cp  # noqa: E402
from pfcommon.paths import EXP4_BUILDS  # noqa: E402

N = 10                                   # cells in x and y (SIZE.h of fjord_vmf_column)
NAMELISTS = ("data", "data.free", "data.kpp", "data.diagnostics", "eedata")


def _column(path: Path, shape, ji):
    j, i = ji
    m = np.memmap(path, dtype=">f8", mode="r", shape=shape)
    return np.array(m[..., j, i])


def make_column_inputs(f1_input: Path, out: Path, depth_m: float = 153.3728, ji=(114, 62)) -> dict:
    """Write the column inputs from F1's input directory (see the module docstring)."""
    f1_input, out = Path(f1_input), Path(out)
    p = json.loads((f1_input / "case.json").read_text())["params"]
    nx, ny, nz = int(p["nx"]), int(p["ny"]), int(p["nz"])
    out.mkdir(parents=True, exist_ok=True)
    for f in NAMELISTS:
        t = (f1_input / f).read_text()
        if f in ("data", "data.free"):
            t = re.sub(r"delX=\d+\*", f"delX={N}*", t)
            t = re.sub(r"delY=\d+\*", f"delY={N}*", t)
        (out / f).write_text(t)
    pkg = (f1_input / "data.pkg").read_text()
    (out / "data.pkg").write_text(re.sub(r"useOBCS=\.TRUE\.", "useOBCS=.FALSE.", pkg))
    T = _column(f1_input / "T.init.bin", (nz, ny, nx), ji)
    S = _column(f1_input / "S.init.bin", (nz, ny, nx), ji)
    for name, col in (("T.init.bin", T), ("S.init.bin", S)):
        np.ascontiguousarray(np.broadcast_to(col[:, None, None], (nz, N, N)), ">f8").tofile(out / name)
    nrec = (f1_input / "taux.bin").stat().st_size // (8 * nx * ny)
    for name in ("taux.bin", "tauy.bin"):
        tau = _column(f1_input / name, (nrec, ny, nx), ji)
        np.ascontiguousarray(np.broadcast_to(tau[:, None, None], (nrec, N, N)), ">f8").tofile(out / name)
    np.full((N, N), -float(depth_m)).astype(">f8").tofile(out / "bathy.bin")
    info = {"source": str(f1_input), "station_ji": list(ji), "depth_m": float(depth_m), "nx": N,
            "ny": N, "nz": nz, "wind_records": int(nrec),
            "boundaries": "doubly periodic (no land, useOBCS=.FALSE.: MITgcm default periodicity)",
            "post_hoc": "PD2 column test (Esteban, 2026-11-13)"}
    (out / "case.json").write_text(json.dumps({"params": {"nx": N, "ny": N, "nz": nz, "case": "column"},
                                               "actual": info}, indent=1))
    return info


def _sh(cmd) -> int:
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
    return subprocess.run([str(c) for c in cmd], env=env).returncode


def run(inp: Path, run_dir: Path, build: Path = EXP4_BUILDS / "fjord_vmf_column", sh=_sh) -> str:
    """Forced stage (data) then free stage (data.free), as F1 (0-60 h, 60-280 h), on 1 rank."""
    inp = Path(inp)
    n1 = cp.nml_int((inp / "data").read_text(), "nTimeSteps")
    nf = cp.nml_int((inp / "data.free").read_text(), "nTimeSteps")
    texts = {1: (inp / "data").read_text(), 2: (inp / "data.free").read_text()}
    return cp.run_stages(Path(run_dir), inp, Path(build), 1, 2, sh, REPO / "scripts" / "run_mitgcm.sh",
                         plan=[(0, n1), (n1, nf)], stage_data=texts.__getitem__)


def analyse(run_dir: Path) -> dict:
    """V2a estimator (validate.v2a_inertial_peak: surface HKE, one Hann window, 0-T_fin) at the
    centre cell of the column run, and the pre-fixed classification rule."""
    from exp4_arctic_bayiw.analysis import load, validate as va
    ds = load.run(Path(run_dir))
    th = va.load_theory()
    f = th["cases"]["vmf_summer"]["f_s-1"]
    j = i = N // 2
    u = 0.5 * (ds["U"].isel(k=0, j=j, i=i).values + ds["U"].isel(k=0, j=j, i=i + 1).values)
    v = 0.5 * (ds["V"].isel(k=0, j=j, i=i).values + ds["V"].isel(k=0, j=j + 1, i=i).values)
    t = ds["time"].values.astype(float)
    r = va.v2a_inertial_peak(u, v, t, f)
    U = ds["U"].isel(k=0).values
    r["horizontal_uniformity_max_abs_dev_m_s"] = float(np.abs(U - U[:, :1, :1]).max())
    r["peak_over_f"] = r["f_peak_hz"] / r["f_hz"]
    r["record_h"] = float((t[-1] - t[0] + (t[1] - t[0])) / 3600)
    r["rule"] = ("fixed before running (PD2, 2026-11-13): peak <= 1 df of f -> the model inertial "
                 "frequency is right and the F1 blue shift is (c); > 1 df -> (a)")
    # the V2a evaluator's "pass" is renamed: this is a rule outcome, not a validation verdict
    r["within_1df"] = r.pop("pass")
    r["v2a_class"] = "(c)" if r["within_1df"] else "(a)"
    r["post_hoc"] = True
    return r


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("make")
    m.add_argument("--f1-input", type=Path, required=True)
    m.add_argument("--out", type=Path, required=True)
    m.add_argument("--depth", type=float, default=153.3728)
    m.add_argument("--ji", type=int, nargs=2, default=(114, 62))
    r = sub.add_parser("run")
    r.add_argument("--input", type=Path, required=True)
    r.add_argument("--run", type=Path, required=True)
    r.add_argument("--allow-concurrent", action="store_true",
                   help="run although other MPI/model processes exist (the declared exception: np = 1 "
                        "next to a 4-rank run, 5 ranks in total; PD2 2026-11-13)")
    a_ = sub.add_parser("analyse")
    a_.add_argument("--run", type=Path, required=True)
    a_.add_argument("--out", type=Path, default=PKG / "figures" / "pd2_column_test.json")
    a = ap.parse_args(argv)
    if a.cmd == "make":
        print(json.dumps(make_column_inputs(a.f1_input, a.out, a.depth, tuple(a.ji))))
        return 0
    if a.cmd == "run":
        busy = cp.other_mitgcm_processes()                 # same check as run_campaign / restart_test
        if busy and not a.allow_concurrent:
            raise RuntimeError("another MPI / model process is running (" + " | ".join(busy) + "); "
                               "pass --allow-concurrent only for the declared np = 1 exception")
        st = run(a.input, a.run)
        print(st)
        return 0 if st == "OK" else 1
    res = analyse(a.run)
    from exp4_arctic_bayiw.analysis import validate as va
    keep = {k: v for k, v in res.items() if k != "spectrum"}
    a.out.write_text(json.dumps(keep, indent=1, default=va.jsonable))
    print(json.dumps(keep, indent=1, default=va.jsonable))
    return 0


if __name__ == "__main__":
    sys.exit(main())
