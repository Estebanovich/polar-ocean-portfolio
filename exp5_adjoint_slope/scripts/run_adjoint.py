#!/opt/anaconda3/envs/pyocean/bin/python
"""Run an Exp. 5 Tapenade case (chan10_tap): adjoint, tangent-linear or forward, from a spin-up
pickup, with time restarting at 0.

Usage::

    run_adjoint.py --case GS [--max-concurrent K] [--dry-run]

Case keys (cases.yaml): ``kind`` (adj | tlm | fwd), ``window_days``, ``last_days`` (the J
window, data.cost lastinterval), ``ctrl_period_days`` (0 = time-constant wind controls; else the
gentim2d record spacing, and window_days must be a multiple of it, because without pkg/cal the
records are periodic over externForcingCycle = window; eesupp/src/get_periodic_interval.F),
``pickup_from`` + ``pickup_iter`` (the spin-up pickup, copied as pickup.G0endstate, read through
pickupSuff with nIter0 = 0), ``grdchk`` (list of {var, i, j, k, rec, eps: [..]}: one serial
subrun per point and eps), ``snapshots`` (true: every-step snapshots of VVEL, THETA, GM_PsiY,
GM_Kvy, GM_Kvz, GM_vbT and DFyE_TH for the J anchor). Executables: <EXP5_BUILDS>/chan10_tap_<kind>/ (scripts/build_mitgcm_tap.sh).
Each run directory gets run_status.json (OK if the model log ends normally), fc.json (the
final fc) and, with grdchk, grdchk.json."""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
REPO = PKG.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pfcommon import campaign as cp  # noqa: E402
from pfcommon.paths import EXP5_BUILDS, EXP5_RUNS  # noqa: E402
from exp5_adjoint_slope.analysis import grdchk  # noqa: E402
from exp5_adjoint_slope.gendata.make_case import PKG, load_case, write_case  # noqa: E402
from exp5_adjoint_slope.gendata.namelist import set_param  # noqa: E402

SUFF = "G0endstate"                       # pickupSuff: exactly 10 characters (read_pickup.F:89)
EXE = {"adj": "mitgcmuv_tap_adj", "tlm": "mitgcmuv_tap_tlm", "fwd": "mitgcmuv"}
_FC = re.compile(r"global fc =\s*([-+]?\d\.\d+D[-+]\d+)")


def subruns(p: dict) -> list[dict]:
    if not p.get("grdchk"):
        return [{"name": None}]
    out = []
    vars_ = [pt["var"] for pt in p["grdchk"]]
    for pt in p["grdchk"]:
        tag = pt["var"]
        if vars_.count(pt["var"]) > 1:              # same control at several points
            tag += f"_i{pt['i']}j{pt['j']}k{pt['k']}"
        for n, e in enumerate(pt["eps"]):
            out.append({**{k: v for k, v in pt.items() if k != "eps"}, "eps": e,
                        "name": f"{tag}_e{n}"})
    return out


def last_fc(log: str) -> float | None:
    hits = _FC.findall(log)
    return float(hits[-1].replace("D", "E")) if hits else None


def _edit(path: Path, group: str, pairs: dict) -> None:
    t = path.read_text()
    for k, v in pairs.items():
        t = set_param(t, group, k, v)
    path.write_text(t)


def install_xx(inp: Path, p: dict) -> None:
    """``xx_from: <dir>``: copy that directory's xx_<name>.<it>.data/.meta (perturbed controls
    for the V2/V3 forward runs) into the input and keep them: pkg/ctrl rewrites them as zeros
    unless doInitXX = .FALSE. (pkg/ctrl/ctrl_readparms.F:171, ctrl_init_ctrlvar.F:140)."""
    if not p.get("xx_from"):
        return
    src = Path(p["xx_from"])
    files = [f for f in src.glob("xx_*") if not f.name.startswith("._")]
    if not files:
        raise FileNotFoundError(f"no xx_* files in {src}")
    for f in files:
        shutil.copyfile(f, Path(inp) / f.name)
    _edit(Path(inp) / "data.ctrl", "CTRL_NML", {"doInitXX": ".FALSE."})


def edit_namelists(inp: Path, p: dict, sub: dict) -> None:
    inp = Path(inp)
    data = (inp / "data").read_text()
    dt = cp.nml_float(data, "deltaT")
    win = float(p["window_days"]) * 86400.0
    nts = win / dt
    if abs(nts - round(nts)) > 1e-9:
        raise ValueError(f"window {win} s is not a multiple of deltaT {dt} s")
    per = float(p.get("ctrl_period_days", 0)) * 86400.0
    d = {"nIter0": "0", "nTimeSteps": str(int(round(nts)))}
    if p.get("pickup_from"):
        d["pickupSuff"] = f"'{SUFF}'"
    if per > 0 and abs(win / per - round(win / per)) > 1e-9:
        raise ValueError(f"window {win} s is not a multiple of the control period {per} s")
    # One forcing record per cycle, reloaded every step: pkg/ctrl adds xx_fu/xx_fv to fu/fv each
    # step (pkg/ctrl/ctrl_map_forcing.F:119, model/src/forward_step.F:572) and only periodic
    # forcing resets fu from the file (model/src/external_fields_load.F:268). Without this the
    # wind control accumulates step after step.
    d["periodicExternalForcing"] = ".TRUE."
    d["externForcingPeriod"] = repr(win)
    d["externForcingCycle"] = repr(win)
    _edit(inp / "data", "PARM03", d)
    _edit(inp / "data.cost", "COST_NML", {"lastinterval": repr(float(p["last_days"]) * 86400.0)})
    _edit(inp / "data.ctrl", "CTRL_NML_GENARR",
          {"xx_gentim2d_period(1)": repr(per), "xx_gentim2d_period(2)": repr(per)})
    on = sub.get("name") is not None
    _edit(inp / "data.pkg", "PACKAGES", {"useGrdchk": ".TRUE." if on else ".FALSE."})
    if on:
        _edit(inp / "data.grdchk", "GRDCHK_NML", {
            "grdchkvarname": f"'{sub['var']}'", "grdchk_eps": repr(float(sub["eps"])),
            "iGloPos": str(sub["i"]), "jGloPos": str(sub["j"]), "kGloPos": str(sub["k"]),
            "recglo": str(sub.get("rec", 1)), "nbeg": "0", "nstep": "1", "nend": "0"})
    if p.get("snapshots"):
        (inp / "data.diagnostics").write_text(
            "# J anchor: every-step snapshots (negative frequency = snapshot).\n"
            " &DIAGNOSTICS_LIST\n"
            " fields(1:7,1) = 'VVEL    ','THETA   ','GM_PsiY ','GM_Kvy  ','GM_Kvz  ',"
            "'GM_vbT  ','DFyE_TH ',\n"
            " fileName(1) = 'snap3d',\n"
            f" frequency(1) = {-dt},\n"
            " &\n &DIAG_STATIS_PARMS\n &\n")
    for fname, groups in (p.get("nml") or {}).items():
        for group, pairs in groups.items():
            _edit(inp / fname, group, {k: str(v) for k, v in pairs.items()})


def find_last_pickup(case: str) -> tuple[Path, Path]:
    hits = sorted(Path(EXP5_RUNS, case).rglob("pickup.[0-9]*.data"),
                  key=lambda f: int(f.name.split(".")[1]))
    hits = [h for h in hits if not h.parent.name.startswith("input")]
    if not hits:
        raise FileNotFoundError(f"no pickup under {EXP5_RUNS / case}")
    return hits[-1], hits[-1].with_suffix(".meta")


def remove_tree(d: Path) -> None:
    """shutil.rmtree that tolerates files already gone: on exFAT, deleting <name> also deletes
    its AppleDouble twin ._<name>, which rmtree has already listed."""
    def missing_ok(func, path, exc_info):
        if not issubclass(exc_info[0], FileNotFoundError):
            raise exc_info[1]
    shutil.rmtree(d, onerror=missing_ok)


def copy_inputs(inp: Path, rd: Path) -> None:
    """Copy the case inputs into the run directory, skipping macOS AppleDouble files (._*),
    which exFAT volumes create and remove behind our back."""
    for f in inp.iterdir():
        if not f.name.startswith("._"):
            shutil.copyfile(f, rd / f.name)


def run_one(case: str, p: dict, sub: dict) -> str:
    rd = EXP5_RUNS / case if sub["name"] is None else EXP5_RUNS / case / sub["name"]
    if rd.exists():
        remove_tree(rd)
    inp = write_case(case, rd / "input")
    edit_namelists(inp, p, sub)
    install_xx(inp, p)
    if p.get("pickup_from"):
        if p.get("pickup_iter"):
            tag = f"{int(p['pickup_iter']):010d}"
            src = next(f for f in Path(EXP5_RUNS, p["pickup_from"]).rglob(f"pickup.{tag}.data")
                       if not f.parent.name.startswith("input"))
            pair = (src, src.with_suffix(".meta"))
        else:
            pair = find_last_pickup(p["pickup_from"])
        for f in pair:
            shutil.copyfile(f, inp / f"pickup.{SUFF}{f.suffix}")
        (rd / "pickup_source.txt").write_text(str(pair[0]) + "\n")
    copy_inputs(inp, rd)
    exe = EXP5_BUILDS / f"chan10_tap_{p['kind']}" / EXE[p["kind"]]
    shutil.copyfile(exe, rd / "mitgcmuv")
    os.chmod(rd / "mitgcmuv", 0o755)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    with open(rd / "output.txt", "w") as out, open(rd / "stderr.txt", "w") as err:
        rc = subprocess.run(["./mitgcmuv"], cwd=rd, stdout=out, stderr=err, env=env).returncode
    log = (rd / "output.txt").read_text(errors="replace")
    ok = rc == 0 and "PROGRAM MAIN: Execution ended Normally" in log
    bad = re.search(r"(%MON|cg2d).*(NaN|Infinity)", log)
    state = "OK" if ok and not bad else "FAILED"
    (rd / "run_status.json").write_text(json.dumps(
        {"state": state, "rc": rc, "nonfinite": bool(bad)}) + "\n")
    (rd / "fc.json").write_text(json.dumps({"fc": last_fc(log)}) + "\n")
    if sub["name"] is not None:
        (rd / "grdchk.json").write_text(json.dumps(grdchk.parse(log), indent=1) + "\n")
    return state


def jobs(pattern: str) -> list[tuple[str, dict, dict]]:
    """(case, params, subrun) for every registry case matching the glob, in registry order."""
    import fnmatch
    import yaml
    names = [c for c in yaml.safe_load((PKG / "cases.yaml").read_text())["cases"]
             if fnmatch.fnmatchcase(c, pattern)]
    if not names:
        raise KeyError(f"no case matches {pattern!r}")
    out = []
    for c in names:
        p = load_case(c)
        out += [(c, p, s) for s in subruns(p)]
    return out


def label(job: tuple[str, dict, dict]) -> str:
    return job[0] if job[2]["name"] is None else f"{job[0]}/{job[2]['name']}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True, help="case name or glob (one pool for all)")
    ap.add_argument("--max-concurrent", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    js = jobs(a.case)
    for c, p, s in js:
        print(f"{c}/{s['name'] or '.'}: kind={p['kind']} window={p['window_days']} d "
              f"last={p['last_days']} d pickup={p.get('pickup_from', '-')} "
              f"{'' if s['name'] is None else s}", file=sys.stderr)
    if a.dry_run:
        return 0
    busy = cp.other_mitgcm_processes()
    if busy:
        print("error: other MITgcm processes running:\n" + "\n".join(busy), file=sys.stderr)
        return 2
    by = {label(j): j for j in js}
    res = cp.run_pool(list(by), 1, a.max_concurrent, lambda n: run_one(*by[n]))
    print(json.dumps(res, indent=1))
    return 0 if all(v == "OK" for v in res.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
