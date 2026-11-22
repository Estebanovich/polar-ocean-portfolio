#!/opt/anaconda3/envs/pyocean/bin/python
"""Run Exp. 6 MITgcm cases (plan Task 4).

Usage::

    run_campaign.py --cases '<glob>' [--max-concurrent K] [--n-timesteps N] [--continue]
                    [--dry-run] [--skip-ok] [--runs-root DIR]

Every case matching the glob gets its inputs from gendata.make_case in <runs>/<case>/input/ and
runs with scripts/run_mitgcm.sh in <runs>/<case>/ on the np of its config (cases.yaml).
Start iteration:
  * pickup_from + branch_year (coarse branches): the parent's pickups at
    branch_iter(branch_year, iter_1850), where iter_1850 = last_iter of G0s (the first iteration
    of GC/GT/GTa, 1 January 1850). The run waits (polling) until the parent has written them.
  * pickup_from without branch_year: the parent's final pickups (its last_iter).
  * none: a cold start at iteration 0. With --continue, the case continues in its own run
    directory from its latest pickup (the G0 spin-up in 100-year stages).
  * F0_<Y> (fine spin-up): a cold start at iteration 0 from GT at branch_iter(Y), theta and the six
    tracers interpolated to the 10-km grid (gendata.interp.fine_init, theta noise member 0); with
    --continue, from its own latest pickup (V6d stages).
  * FT_<Y>_{c,w}<k>: F0_<Y>'s final pickups with theta noise of member k on the copy, the same k
    for c and w (spec Sec. 5: the pair shares its drift).
The three pickup families written by the model (pickup, pickup_ptracers, pickup_dic) are copied,
never linked (LaCie is exFAT)."""
from __future__ import annotations

import argparse
import fnmatch
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
REPO = PKG.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402
import yaml  # noqa: E402

from pfcommon import campaign as cp  # noqa: E402
from exp6_southern_carbon.gendata import make_case as mc  # noqa: E402

RUN_SH = REPO / "scripts" / "run_mitgcm.sh"
BUILD_SH = REPO / "scripts" / "build_mitgcm.sh"
PICKUP_STEMS = ("pickup", "pickup_ptracers", "pickup_dic")
STEPS_G50 = 31536                       # 365-day year at deltaT = 1000 s


def _sh(cmd: list) -> int:
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    return subprocess.run([str(c) for c in cmd], env=env).returncode


def select(pattern: str) -> list[str]:
    reg = yaml.safe_load((PKG / "cases.yaml").read_text())
    return [c for c in reg["cases"] if fnmatch.fnmatch(c, pattern)]


def branch_iter(year: int, iter_1850: int, steps_per_year: int = STEPS_G50) -> int:
    return int(iter_1850) + (int(year) - 1850) * int(steps_per_year)


def _status(runs_root: Path, case: str) -> dict:
    f = Path(runs_root) / case / "run_status.json"
    if not f.exists():
        raise FileNotFoundError(f"{f} (run {case} first)")
    return json.loads(f.read_text())


def pickup_files(runs_root: Path, case: str, itr: int) -> list[Path]:
    d = Path(runs_root) / case
    files = [d / f"{s}.{itr:010d}{e}" for s in PICKUP_STEMS for e in (".data", ".meta")]
    missing = [f.name for f in files if not f.exists()]
    if missing:
        raise FileNotFoundError(f"{case}: missing {', '.join(missing)}")
    return files


def latest_pickup(run_dir: Path) -> int | None:
    its = [int(f.name.split(".")[1]) for f in Path(run_dir).glob("pickup.*.meta") if not f.name.startswith("._")]
    return max(its) if its else None


def _wait_pickups(runs_root, parent, itr, wait_s, poll_s=60.0):
    t0 = time.time()
    while True:
        try:
            return pickup_files(runs_root, parent, itr)
        except FileNotFoundError:
            st = Path(runs_root) / parent / "run_status.json"
            failed = st.exists() and json.loads(st.read_text()).get("state") == "FAILED"
            if failed or time.time() - t0 >= wait_s:
                raise
            time.sleep(poll_s)


def member_of(case: str) -> int | None:
    m = re.fullmatch(r"FT_\d{4}_[cw](\d+)", case)
    return int(m.group(1)) if m else None


def fine_init_from_run(runs_root: Path, parent: str, itr: int, case: str) -> dict:
    """make_case `init` of an F0 case from the parent's (GT) pickups at itr."""
    from exp6_southern_carbon.analysis import io
    from exp6_southern_carbon.gendata import fields, interp
    pr = Path(runs_root) / parent
    g = io.grid(pr)
    vol_c = g["hfac"] * g["rac"][None] * g["drf"][:, None, None]
    th = io.fields(pr, "pickup", itr)["Theta"]
    ptr = io.fields(pr, "pickup_ptracers", itr)
    tr = np.stack([ptr[k] for k in sorted(ptr)])
    if tr.shape[0] != 6:
        raise ValueError(f"{parent} pickup_ptracers.{itr}: {sorted(ptr)}, expected 6 tracers")
    gc, gf = mc.grid_of(mc.load_case(parent)), mc.grid_of(mc.load_case(case))
    vol_f = interp.wet_volume(gf, fields.bathymetry(gf))
    return interp.fine_init(th, tr, gc, gf, vol_c, vol_f, member=0)


def plan(cases: list[str], runs_root: Path, n_timesteps: int | None = None, cont: bool = False,
         wait_s: float = 0.0) -> list[dict]:
    reg = yaml.safe_load((PKG / "cases.yaml").read_text())
    out = []
    for c in cases:
        p = mc.load_case(c)
        spy = mc.steps_per_year(p)
        extra: dict = {}
        par = p["pickup_from"]
        pickups: list[Path] = []
        if c.startswith("F0_"):
            n0 = 0
            if cont:
                last = latest_pickup(Path(runs_root) / c)
                if last is None:
                    raise FileNotFoundError(f"{c}: --continue but no pickup in {Path(runs_root) / c}")
                n0, pickups = last, pickup_files(runs_root, c, last)
            else:
                it = branch_iter(p["branch_year"], int(_status(runs_root, "G0s")["last_iter"]))
                pickup_files(runs_root, par, it)
                extra["fine_from"] = {"parent": par, "iter": it}
        elif par is None:
            n0 = 0
            if cont:
                last = latest_pickup(Path(runs_root) / c)
                if last is None:
                    raise FileNotFoundError(f"{c}: --continue but no pickup in {Path(runs_root) / c}")
                n0, pickups = last, pickup_files(runs_root, c, last)
        elif p["branch_year"] is not None and p["config"] != "f10":
            extra["iter_1850"] = int(_status(runs_root, "G0s")["last_iter"])
            n0 = branch_iter(p["branch_year"], extra["iter_1850"], spy)
            pickups = _wait_pickups(runs_root, par, n0, wait_s)
        else:
            n0 = int(_status(runs_root, par)["last_iter"])
            pickups = _wait_pickups(runs_root, par, n0, wait_s)
            if c in ("GC", "GT", "GTa"):
                extra["iter_1850"] = n0
        if p["pco2_mode"] == "transient" and "iter_1850" not in extra:
            extra["iter_1850"] = int(_status(runs_root, "G0s")["last_iter"])
        if n_timesteps is not None:
            nts = int(n_timesteps)
        elif p["n_years"] is not None:
            nts = int(p["n_years"]) * spy
        else:
            raise ValueError(f"{c}: n_years is null; pass --n-timesteps")
        extra.update(nIter0=n0, n_timesteps=nts)
        np_ = int(reg["configs"][p["config"]]["np"])
        out.append({"case": c, "config": p["config"], "np": np_, "nIter0": n0, "nTimeSteps": nts,
                    "pickups": pickups, "extra": extra})
    return out


def rmtree_exfat(d: Path) -> None:
    """shutil.rmtree that tolerates files vanishing mid-walk: on exFAT (LaCie) macOS deletes the
    ._ AppleDouble companion together with its file, so rmtree then fails on the ._ entry."""
    def onerror(func, path, exc):
        if not issubclass(exc[0], FileNotFoundError):
            raise exc[1]
    shutil.rmtree(d, onerror=onerror)


def check_size(case_dir: Path, size_h: Path) -> None:
    """RuntimeError unless SIZE.h matches the flat Exp. 6 case.json (nx, ny, nr). The pfcommon
    version reads the Exp. 5 schema (case.json["params"]) and cannot be used here."""
    c = json.loads((Path(case_dir) / "case.json").read_text())
    t = Path(size_h).read_text()
    sp = cp._size_param
    got = {"Nx": sp(t, "sNx") * sp(t, "nSx") * sp(t, "nPx"), "Ny": sp(t, "sNy") * sp(t, "nSy") * sp(t, "nPy"),
           "Nr": sp(t, "Nr")}
    want = {"Nx": c["nx"], "Ny": c["ny"], "Nr": c["nr"]}
    bad = [f"{k}: SIZE.h {got[k]} != case {want[k]}" for k in got if got[k] != want[k]]
    if bad:
        raise RuntimeError(f"{size_h} does not match {case_dir}: " + "; ".join(bad))


def prepare(item: dict, runs_root: Path) -> Path:
    inp = Path(runs_root) / item["case"] / "input"
    if inp.exists():
        rmtree_exfat(inp)
    extra = dict(item["extra"])
    ff = extra.pop("fine_from", None)
    if ff:
        extra["init"] = fine_init_from_run(runs_root, ff["parent"], ff["iter"], item["case"])
    if item["config"] != "f10" and mc.load_case(item["case"])["sponge"] == "anth":
        extra["cant_sup"] = {int(k): v for k, v in json.loads((Path(runs_root) / "cant_sup.json").read_text()).items()}
    mc.write_case(item["case"], inp, extra)
    for f in item["pickups"]:
        shutil.copyfile(f, inp / f.name)
    k = member_of(item["case"])
    if k is not None:
        from exp5_adjoint_slope.gendata.interp import perturb_pickup_theta
        from exp6_southern_carbon.gendata.interp import MEMBER_AMP, MEMBER_SEED
        g = mc.grid_of(mc.load_case(item["case"]))
        perturb_pickup_theta(inp / f"pickup.{item['nIter0']:010d}.data", g.nr, g.ny, g.nx,
                             k=k, amp=MEMBER_AMP, seed=MEMBER_SEED)
    check_size(inp, PKG / "config" / item["config"] / "code" / "SIZE.h")
    return inp


def main(argv=None) -> int:
    from pfcommon.paths import EXP6_BUILDS, EXP6_RUNS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--max-concurrent", type=int, default=1)
    ap.add_argument("--n-timesteps", type=int)
    ap.add_argument("--continue", dest="cont", action="store_true")
    ap.add_argument("--skip-ok", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--runs-root", type=Path, default=EXP6_RUNS)
    ap.add_argument("--wait-hours", type=float, default=48.0)
    a = ap.parse_args(argv)
    cases = select(a.cases)
    if not cases:
        print(f"error: no case matches {a.cases!r}", file=sys.stderr)
        return 2
    if a.dry_run:
        print("\n".join(cases))
        return 0
    reg = yaml.safe_load((PKG / "cases.yaml").read_text())
    nps = {int(reg["configs"][mc.load_case(c)["config"]]["np"]) for c in cases}
    try:
        if len(nps) != 1:
            raise ValueError(f"cases mix np values {sorted(nps)}; run them separately")
        cp.check_concurrency(nps.pop() if len(nps) == 1 else 1, a.max_concurrent)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    busy = cp.other_mitgcm_processes()
    if busy:
        print("error: other MITgcm processes running:\n" + "\n".join(busy), file=sys.stderr)
        return 2
    if a.skip_ok:
        cases = [c for c in cases if not cp.finished_ok(a.runs_root / c)]
    for cfg in sorted({mc.load_case(c)["config"] for c in cases}):
        np_ = int(reg["configs"][cfg]["np"])
        bdir = EXP6_BUILDS / cfg
        if not cp.build_up_to_date(PKG / "config" / cfg, bdir, np_ > 1):
            if _sh([BUILD_SH, PKG / "config" / cfg, bdir] + (["--mpi"] if np_ > 1 else [])) != 0:
                print(f"error: build of {cfg} failed", file=sys.stderr)
                return 1

    def run_one(c):
        (item,) = plan([c], a.runs_root, a.n_timesteps, a.cont, wait_s=a.wait_hours * 3600)
        inp = prepare(item, a.runs_root)
        rd = a.runs_root / c
        _sh([RUN_SH, EXP6_BUILDS / item["config"], inp, rd, item["np"]])
        st = json.loads((rd / "run_status.json").read_text())
        return st["state"]

    res = cp.run_pool(cases, int(reg["configs"][mc.load_case(cases[0])["config"]]["np"]) if cases else 1,
                      a.max_concurrent, run_one)
    print(json.dumps(res, indent=1))
    return 0 if all(v == "OK" for v in res.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
