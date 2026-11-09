"""Campaign helpers shared across experiments: case selection, process and disk checks,
a bounded concurrent pool of runs, and chained restart stages.

Taken from exp1_fjord_plumes/scripts/run_campaign.py without its Exp. 1 specifics (Stevens /
Orlanski OB pickups, the fixed two-stage wind cases). Additions:

- ``run_pool``: runs independent cases side by side. Resource rule (global constraints;
  Esteban 2026-11-05): either up to 5 serial (np = 1) runs at a time, or one np = 5 run.
- ``run_stages``: splits a run into N restart stages. Stage k restarts from the pickup that
  stage k-1 wrote at its last iteration, with nIter0 set to that iteration. Inputs and
  pickups are copied, never linked (the run disk is exFAT).

The MITgcm executable is never started here: callers pass ``sh`` (a list -> return code
callable) that runs scripts/run_mitgcm.sh, so tests can inject a fake runner.
"""
from __future__ import annotations

import fnmatch
import json
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable

EXE_NAME = "mitgcm" + "uv"
MAX_NP = 5
MAX_SERIAL_CONCURRENT = 5
# Match the process *name* only (no -f): with -f, any shell whose command line merely
# mentions the binary (an editor, a log grep) was taken for a model run.
PROC_PATTERN = f"^(mpirun|mpiexec|orterun|prterun|{EXE_NAME})$"
PGREP_ARGS = ("-l", PROC_PATTERN)
MIN_FREE_BYTES = 20 * 10**9          # always keep at least 20 GB free on the runs filesystem
FREE_FACTOR = 4                      # ... and at least 4 x the generated input dirs


# ---- selection and arguments -------------------------------------------------------------

def select_cases(cases: dict, config: str, pattern: str) -> list[str]:
    """Names in ``cases`` (registry order) of ``config`` matching the fnmatch ``pattern``.
    Raises ValueError if none; the message names matches of other configurations."""
    hits = [n for n in cases if fnmatch.fnmatchcase(n, pattern)]
    mine = [n for n in hits if cases[n]["config"] == config]
    if not mine:
        other = [n for n in hits if n not in mine]
        extra = f" (matches of other configs: {', '.join(other)})" if other else ""
        raise ValueError(f"no case of config {config!r} matches {pattern!r}{extra}")
    return mine


def parse_sets(items: list[str] | None) -> dict:
    """['KEY=VALUE', ...] -> {KEY: int or float} (``--set``)."""
    out = {}
    for it in items or []:
        if "=" not in it:
            raise ValueError(f"--set expects KEY=VALUE, got {it!r}")
        k, v = (x.strip() for x in it.split("=", 1))
        try:
            out[k] = int(v)
        except ValueError:
            out[k] = float(v)
    return out


def check_overrides(names: list[str], overrides: dict | None, spec: dict) -> dict:
    """Reject keys absent from the merged parameters of any selected case; give each value
    the type of the parameter it replaces. A null parameter (e.g. an uncalibrated gammaT)
    takes a float."""
    if not overrides:
        return overrides or {}
    out = dict(overrides)
    for name in names:
        case = spec["cases"][name]
        merged = {**spec["configs"][case["config"]], **case}
        for k, v in overrides.items():
            if k not in merged:
                raise ValueError(f"--set {k}: not a parameter of case {name}")
            cur = merged[k]
            if cur is None:
                out[k] = float(v)
            elif isinstance(cur, bool) or not isinstance(cur, (int, float)):
                raise ValueError(f"--set {k}: parameter is not numeric")
            elif isinstance(cur, int):
                if float(v) != int(float(v)):
                    raise ValueError(f"--set {k}={v}: integer parameter")
                out[k] = int(v)
            else:
                out[k] = float(v)
    return out


def check_concurrency(np_each: int, max_concurrent: int) -> None:
    """Allowed: np = 1 with 1..5 concurrent runs, or np = 5 with exactly 1. Anything else
    raises ValueError (global constraints: <= 5 MPI processes in total)."""
    if np_each == 1 and 1 <= max_concurrent <= MAX_SERIAL_CONCURRENT:
        return
    if np_each == MAX_NP and max_concurrent == 1:
        return
    raise ValueError(f"np={np_each} with max_concurrent={max_concurrent} not allowed: use np=1 "
                     f"with max_concurrent 1..{MAX_SERIAL_CONCURRENT}, or np={MAX_NP} with 1")


# ---- environment checks ------------------------------------------------------------------

def _pgrep() -> str:
    return subprocess.run(["pgrep", *PGREP_ARGS], capture_output=True, text=True).stdout


def other_mitgcm_processes(pgrep: Callable[[], str] | None = None) -> list[str]:
    """MPI launcher / model processes (by process name), minus this process."""
    out = (pgrep or _pgrep)()
    me = str(os.getpid())
    return [ln.strip() for ln in out.splitlines() if ln.strip() and ln.split()[0] != me]


def _size_param(text: str, name: str) -> int:
    m = re.search(rf"^\s*&?\s*{name}\s*=\s*(\d+)", text, re.MULTILINE)
    if not m:
        raise RuntimeError(f"SIZE.h: cannot find {name}")
    return int(m.group(1))


def check_size(case_dir: Path, size_h: Path) -> None:
    """Raise RuntimeError unless SIZE.h (sNx nSx nPx, sNy nSy nPy, Nr) matches case.json."""
    p = json.loads((Path(case_dir) / "case.json").read_text())["params"]
    t = Path(size_h).read_text()
    Nx = _size_param(t, "sNx") * _size_param(t, "nSx") * _size_param(t, "nPx")
    Ny = _size_param(t, "sNy") * _size_param(t, "nSy") * _size_param(t, "nPy")
    Nr = _size_param(t, "Nr")
    bad = [f"{k}: SIZE.h {s} != grid {g}" for k, s, g in
           (("Nx", Nx, p["nx"]), ("Ny", Ny, p["ny"]), ("Nr", Nr, p["nz"])) if s != g]
    if bad:
        raise RuntimeError(f"{size_h} does not match {case_dir}: " + "; ".join(bad))


def build_up_to_date(cfg_dir: Path, build_dir: Path, mpi: bool) -> bool:
    """Executable newer than every file in cfg_dir/code, and built with/without MPI as asked."""
    exe, mk = Path(build_dir) / EXE_NAME, Path(build_dir) / "Makefile"
    if not (exe.exists() and mk.exists()):
        return False
    if ("-DALLOW_USE_MPI" in mk.read_text()) != mpi:
        return False
    newest = max(f.stat().st_mtime for f in (Path(cfg_dir) / "code").iterdir() if f.is_file())
    return exe.stat().st_mtime > newest


def _dir_bytes(d: Path) -> int:
    return sum(f.stat().st_size for f in Path(d).rglob("*") if f.is_file())


def check_free_space(runs_root: Path, inputs, disk_usage=None) -> None:
    """Raise RuntimeError if the runs filesystem has less than max(20 GB, 4 x the total size
    of the input dirs ``inputs`` (a path or a list of paths)) free."""
    inputs = [inputs] if isinstance(inputs, (str, Path)) else list(inputs)
    free = (disk_usage or shutil.disk_usage)(str(runs_root)).free
    need = max(MIN_FREE_BYTES, FREE_FACTOR * sum(_dir_bytes(d) for d in inputs))
    if free < need:
        raise RuntimeError(f"not enough free space on {runs_root}: {free / 1e9:.1f} GB free, "
                           f"need {need / 1e9:.1f} GB (max(20 GB, 4 x input dirs)); "
                           "refusing to start the model")


# ---- run status ----------------------------------------------------------------------------

def _status(run_dir: Path) -> dict | None:
    try:
        return json.loads((Path(run_dir) / "run_status.json").read_text())
    except (OSError, json.JSONDecodeError):
        return None


def finished_ok(run_dir: Path) -> bool:
    """True if the run is OK and complete. run_status_final.json (written by run_stages)
    wins when present; otherwise run_status.json must be OK with last_iter >= n_timesteps."""
    run_dir = Path(run_dir)
    fin = run_dir / "run_status_final.json"
    try:
        if fin.exists():
            st = json.loads(fin.read_text())
            return st.get("state") == "OK" and st.get("last_iter") == st.get("expected_last_iter")
        st = json.loads((run_dir / "run_status.json").read_text())
    except (OSError, json.JSONDecodeError):
        return False
    return st.get("state") == "OK" and st.get("last_iter", -1) >= st.get("n_timesteps", 1 << 62)


def _write_final(run_dir: Path, state: str, stages: dict, expected: int,
                 reason: str | None = None) -> str:
    last = [st for st in stages.values() if st]
    fin = {"state": state, "last_iter": last[-1].get("last_iter") if last else None,
           "expected_last_iter": expected, "stages": stages}
    if reason:
        fin["reason"] = reason
    Path(run_dir).mkdir(parents=True, exist_ok=True)
    (Path(run_dir) / "run_status_final.json").write_text(json.dumps(fin, indent=2))
    return state


# ---- concurrency ---------------------------------------------------------------------------

def run_pool(cases: list[str], np_each: int, max_concurrent: int,
             run_one: Callable[[str], str]) -> dict[str, str]:
    """Run ``run_one(case)`` for every case with at most ``max_concurrent`` at a time
    (``check_concurrency``). Returns {case: state} in the order of ``cases``; a case whose
    ``run_one`` raises gets "ERROR: <message>" and the others go on."""
    check_concurrency(np_each, max_concurrent)

    def guarded(case):
        try:
            return run_one(case)
        except Exception as e:  # noqa: BLE001 - one failed case must not stop the others
            return f"ERROR: {type(e).__name__}: {e}"

    with ThreadPoolExecutor(max_workers=max_concurrent) as ex:
        futs = {c: ex.submit(guarded, c) for c in cases}
        return {c: futs[c].result() for c in cases}


# ---- namelists and stages --------------------------------------------------------------------

def _live_lines(text: str):
    for n, ln in enumerate(text.splitlines()):
        if not ln.lstrip().startswith(("#", "!")):
            yield n, ln


def _nml_raw(text: str, key: str) -> str:
    for _, ln in _live_lines(text):
        m = re.match(rf"\s*{key}\s*=\s*([^,\s/]+)", ln, re.IGNORECASE)
        if m:
            return m.group(1)
    raise RuntimeError(f"{key} not found")


def nml_int(text: str, key: str) -> int:
    return int(_nml_raw(text, key))


def nml_float(text: str, key: str) -> float:
    return float(_nml_raw(text, key).lower().replace("d", "e"))


def set_namelist_ints(text: str, values: dict[str, int]) -> str:
    """Replace ``KEY = n`` for each key (each must occur exactly once, outside comments).
    Refuses a namelist that also sets endTime/startTime, which would compete with
    nIter0 + nTimeSteps."""
    lines = text.splitlines(keepends=True)
    for k in ("endTime", "startTime"):
        if any(re.match(rf"\s*{k}\s*=", ln, re.IGNORECASE) for _, ln in _live_lines(text)):
            raise ValueError(f"namelist sets {k}; staged runs need nIter0 + nTimeSteps only")
    for key, val in values.items():
        hits = [n for n, ln in _live_lines(text) if re.match(rf"\s*{key}\s*=", ln, re.IGNORECASE)]
        if len(hits) != 1:
            raise ValueError(f"namelist must set {key} exactly once (found {len(hits)})")
        n = hits[0]
        lines[n] = re.sub(rf"({key}\s*=\s*)[^,\s/]+", rf"\g<1>{int(val)}", lines[n],
                          count=1, flags=re.IGNORECASE)
    return "".join(lines)


def stage_plan(n_iter0: int, n_total: int, n_stages: int) -> list[tuple[int, int]]:
    """[(nIter0, nTimeSteps)] of ``n_stages`` equal stages; n_total must divide evenly."""
    if n_stages < 1 or n_total % n_stages:
        raise ValueError(f"nTimeSteps {n_total} does not split into {n_stages} equal stages")
    n = n_total // n_stages
    return [(n_iter0 + k * n, n) for k in range(n_stages)]


def run_stages(run_dir: Path, inp: Path, build_dir: Path, np_: int, n_stages: int,
               sh: Callable[[list], int], run_sh, pickup_stems=("pickup",)) -> str:
    """Run case input ``inp`` in ``n_stages`` chained restart stages (ruling R2).

    Stage 1 runs in ``run_dir``, stage k > 1 in ``run_dir/stage<k>``. With more than one
    stage, each stage gets a copy ``run_dir/input_stage<k>`` of ``inp`` whose ``data`` has
    nIter0 / nTimeSteps of that stage; for k > 1 it also gets copies of every
    ``pickup*.<iter>.*`` file stage k-1 wrote at its last iteration. Each stem of
    ``pickup_stems`` must be present (``<stem>.<iter:010d>.data`` and ``.meta``). If
    ``data`` sets pChkptFreq, every stage end must fall on a multiple of it, since the
    permanent pickup is written there. Before each stage its run_status.json is removed,
    and a nonzero exit code of the runner makes the case FAILED. Writes
    ``run_dir/run_status_final.json``: OK only if every stage is OK and the last
    iteration is nIter0 + nTimeSteps of the input."""
    run_dir, inp = Path(run_dir), Path(inp)
    text = (inp / "data").read_text()
    n0, ntot = nml_int(text, "nIter0"), nml_int(text, "nTimeSteps")
    plan = stage_plan(n0, ntot, n_stages)
    expected = n0 + ntot
    if n_stages > 1:
        try:
            ck = nml_float(text, "pChkptFreq")
        except RuntimeError:
            ck = 0.0
        if ck > 0:
            dt = nml_float(text, "deltaT")
            for a, n in plan:
                t_end = (a + n) * dt
                if abs(t_end / ck - round(t_end / ck)) > 1e-9:
                    raise ValueError(f"stage end {t_end:g} s is not a multiple of pChkptFreq {ck:g} s")
        set_namelist_ints(text, {"nIter0": n0, "nTimeSteps": ntot})   # validates the namelist
    stages: dict = {}
    prev = None
    for k, (a, n) in enumerate(plan, 1):
        sdir = run_dir if k == 1 else run_dir / f"stage{k}"
        if n_stages == 1:
            sinp = inp
        else:
            sinp = run_dir / f"input_stage{k}"
            if sinp.exists():
                shutil.rmtree(sinp)
            shutil.copytree(inp, sinp, symlinks=False)                # copies, never links
            (sinp / "data").write_text(set_namelist_ints(text, {"nIter0": a, "nTimeSteps": n}))
            if prev is not None:
                for f in sorted(prev.glob(f"pickup*.{a:010d}.*")):
                    shutil.copyfile(f, sinp / f.name)                  # content copy
        stages[str(k)] = None
        # A stale status from an earlier attempt must never be read as this run's result:
        # run_mitgcm.sh can exit (LaCie not mounted, mkdir/cp failure) before removing it.
        (sdir / "run_status.json").unlink(missing_ok=True)
        rc = sh([run_sh, build_dir, sinp, sdir, np_])
        st = _status(sdir)
        stages[str(k)] = st
        if rc != 0:
            return _write_final(run_dir, "FAILED", stages, expected,
                                f"stage {k}: runner exited with code {rc}")
        if not (st and st.get("state") == "OK" and st.get("last_iter") == a + n):
            return _write_final(run_dir, "FAILED", stages, expected,
                                f"stage {k} not OK or did not reach iteration {a + n}")
        if k < n_stages:
            tag = f"{a + n:010d}"
            missing = [f"{s}.{tag}.{e}" for s in pickup_stems for e in ("data", "meta")
                       if not (sdir / f"{s}.{tag}.{e}").exists()]
            if missing:
                return _write_final(run_dir, "FAILED", stages, expected,
                                    f"stage {k}: missing pickup " + ", ".join(missing))
        prev = sdir
    return _write_final(run_dir, "OK", stages, expected)


def compare_bitwise(dir_a: Path, dir_b: Path, pattern: str) -> list[str]:
    """Names matching ``pattern`` that differ byte for byte between two run directories,
    or exist in only one (restart test: a staged run must equal the single run)."""
    a = {f.name: f for f in Path(dir_a).glob(pattern) if f.is_file()}
    b = {f.name: f for f in Path(dir_b).glob(pattern) if f.is_file()}
    if not a and not b:
        return [f"(no files match {pattern})"]
    bad = sorted(set(a) ^ set(b))
    bad += [n for n in sorted(set(a) & set(b)) if a[n].read_bytes() != b[n].read_bytes()]
    return sorted(bad)
