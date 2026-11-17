#!/opt/anaconda3/envs/pyocean/bin/python
"""Run Exp. 5 forward MITgcm cases (chan10, chan2p5): build once per config, then run.

Usage::

    run_campaign.py --cases '<glob>' [--np N] [--max-concurrent K] [--set KEY=VALUE]
                    [--run-suffix S] [--dry-run]

Each case of cases.yaml matching the glob gets its inputs from ``gendata.make_case`` in
``<EXP5_RUNS>/<case><suffix>/input/``; ``n_timesteps`` (and ``n_iter0``) of the case set the
namelist. ``pickup_from: <case>`` with ``pickup_iter`` copies (never links) that case's
``pickup.<iter>.*`` into the input and sets nIter0. ``n_stages`` > 1 runs chained restart stages
(pfcommon.campaign.run_stages). Builds go to ``<EXP5_BUILDS>/<config>/`` with
scripts/build_mitgcm.sh (``--mpi`` if np > 1). np must equal nPx * nPy of the config's SIZE.h;
resource rules: pfcommon.campaign.check_concurrency (<= 5 processes). The adjoint configs
(``*_tap``) are not run here (scripts/run_adjoint.py)."""
from __future__ import annotations

import argparse
import fnmatch
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
REPO = PKG.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import yaml  # noqa: E402

from pfcommon import campaign as cp  # noqa: E402
from pfcommon.paths import EXP5_BUILDS, EXP5_RUNS  # noqa: E402
from exp5_adjoint_slope.gendata import interp  # noqa: E402
from exp5_adjoint_slope.gendata.make_case import grid_of, load_case, write_case  # noqa: E402

RUN_SH = REPO / "scripts" / "run_mitgcm.sh"
BUILD_SH = REPO / "scripts" / "build_mitgcm.sh"


def _sh(cmd: list) -> int:
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    return subprocess.run([str(c) for c in cmd], env=env).returncode


def select(pattern: str) -> list[str]:
    reg = yaml.safe_load((PKG / "cases.yaml").read_text())
    return [c for c in reg["cases"] if fnmatch.fnmatch(c, pattern)
            and not reg["cases"][c]["config"].endswith("_tap")]


MEMBER_AMP, MEMBER_SEED = 1e-3, 2026      # spec Sec. 4: 1e-3 degC in T0


def find_pickup(runs_root: Path, case: str, itr: int) -> list[Path]:
    tag = f"{itr:010d}"
    hits = sorted(Path(runs_root, case).rglob(f"pickup.{tag}.*"))
    hits = [h for h in hits if "input" not in h.parent.name]
    if {h.suffix for h in hits} != {".data", ".meta"}:
        raise FileNotFoundError(f"no pickup.{tag}.data/.meta under {runs_root / case}")
    return hits[:2]


def prepare(case: str, p: dict, run_dir: Path, runs_root: Path) -> Path:
    inp = write_case(case, run_dir / "input", overrides=p)
    text = (inp / "data").read_text()
    n0 = int(p.get("n_iter0", 0))
    if "pickup_from" in p:
        n0 = int(p["pickup_iter"])
        for f in find_pickup(runs_root, p["pickup_from"], n0):
            shutil.copyfile(f, inp / f.name)
        if "member" in p:                       # ensemble member: theta noise on the copy
            g = grid_of(p)
            interp.perturb_pickup_theta(inp / f"pickup.{n0:010d}.data", g.nr, g.ny, g.nx,
                                        k=int(p["member"]), amp=MEMBER_AMP, seed=MEMBER_SEED)
    text = cp.set_namelist_ints(text, {"nIter0": n0, "nTimeSteps": int(p["n_timesteps"])})
    (inp / "data").write_text(text)
    return inp


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--max-concurrent", type=int, default=1)
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--run-suffix", default="")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    try:
        cp.check_concurrency(a.np, a.max_concurrent)
        cases = select(a.cases)
        if not cases:
            raise ValueError(f"no forward case matches {a.cases!r}")
        sets = cp.parse_sets(a.set) if a.set else {}
        params = {c: {**load_case(c), **sets} for c in cases}
        for c, p in params.items():
            cp.check_np(PKG / "config" / p["config"] / "code" / "SIZE.h", a.np)
    except (ValueError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    for c in cases:
        p = params[c]
        print(f"{c}: config={p['config']} run={EXP5_RUNS / (c + a.run_suffix)} "
              f"nTimeSteps={p['n_timesteps']} stages={p.get('n_stages', 1)} "
              f"pickup={p.get('pickup_from', '-')}", file=sys.stderr)
    if a.dry_run:
        print("\n".join(cases))
        return 0
    busy = cp.other_mitgcm_processes()
    if busy:
        print("error: other MITgcm processes running:\n" + "\n".join(busy), file=sys.stderr)
        return 2
    mpi = a.np > 1
    for cfg in sorted({params[c]["config"] for c in cases}):
        bdir = EXP5_BUILDS / cfg
        if not cp.build_up_to_date(PKG / "config" / cfg, bdir, mpi):
            cmd = [BUILD_SH, PKG / "config" / cfg, bdir] + (["--mpi"] if mpi else [])
            if _sh(cmd) != 0:
                print(f"error: build of {cfg} failed", file=sys.stderr)
                return 1

    def run_one(c):
        p = params[c]
        rd = EXP5_RUNS / (c + a.run_suffix)
        inp = prepare(c, p, rd, EXP5_RUNS)
        cp.check_size(inp, PKG / "config" / p["config"] / "code" / "SIZE.h")
        return cp.run_stages(rd, inp, EXP5_BUILDS / p["config"], a.np, int(p.get("n_stages", 1)),
                             _sh, RUN_SH)

    res = cp.run_pool(cases, a.np, a.max_concurrent, run_one)
    print(json.dumps(res, indent=1))
    return 0 if all(v == "OK" for v in res.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
