#!/opt/anaconda3/envs/anemoi/bin/python
r"""Exp. 3 validation protocol (notes/pd0_preregistro.md).

``--freeze-config`` (Task 7) writes ``figures/protocol_config.json`` with everything the
protocol runs must use, fixed at PD0 before any protocol run:

* the dataset sha256 (checked against the file before anything is read);
* the folds of ``folds.make_folds`` and the N1-sensitivity variants (``folds.drop_n1``) of
  every fold except ``final`` (the final test is evaluated once, without a repeat);
* every method's hyperparameters, read from the constructor defaults (``hparams()`` of
  each ``methods.REGISTRY`` factory's instance);
* the PD0 thresholds (``PD0_THRESHOLDS``) and the noise floor (figures/noise_floor.json);
* the git commit and the paths that differed from it at freeze time; any changed tracked
  path outside ``exp3_melt_emulator/notes/`` makes it refuse.

It refuses to overwrite an existing file; there is no override.

``--cv`` and ``--final`` (Task 8) first verify the frozen configuration: its sha256 must be
``CONFIG_SHA256`` (approved at PD0), every method's current ``hparams()`` must equal the
frozen ones, the dataset must have the frozen sha256 and the folds recomputed from it must
equal the frozen folds. Otherwise they refuse (``ProtocolError``).

``--cv`` fits every method on every fold except ``final`` (main folds and their N1
sensitivity repeat; M4/M5 skip ``profile_test``, per ``profile_test_methods``) and
evaluates the fold's test runs; a ``split == "test"`` run is never evaluated. A work unit
is (variant main|n1, fold, method); units run in up to ``MAX_WORKERS`` = 5 spawned worker
processes with 1 thread each, M7 units first and longest first. Every finished unit is
written atomically to ``UNITS_DIR/<unit>.json`` and logged with its wall time in
``CV_PROGRESS_LOG``; a restart skips the units already written with the same config sha.
Stop rule (PD0, frozen): when the first batch (the first ``workers`` units) has finished,
the total wall clock is estimated (``estimate_wall``); above 8 h no new unit is started,
the estimate is logged and the run ends without results (exit code 3). The limit can
be replaced only by ``--wall-limit-hours H`` together with ``--amendment "<text>"`` (a
documented PD0 amendment, e.g. Enmienda 1: 20 h); the override is recorded in the log,
in every unit written under it and in cv_results.json ("amendments"). Main-variant units
are scheduled before any n1 unit (Enmienda 1; ordering does not change results). At the end
``figures/cv_results.json`` is assembled from the unit files: one entry per variant x
fold x method x test run with every Task 3 metric, skill against M0 of the same
fold/variant, coverage (68/95 %) and CRPS for M5-M7 and the H3c spread quantities of
``predict_members`` for M6/M7, plus each unit's method ``info``.

``--final`` refuses if ``figures/final_test.json`` exists or if a previous ``--final``
started (marker ``FINAL_MARKER``); there is no override. It fits the 8 methods on the 25
training runs, predicts the 3 test runs once, and writes ``figures/final_test.json`` and
the maps (test runs, plus the in-sample maps of the 25 training runs for the Lambert
diagnostics) to ``FINAL_MAPS``, whose sha256 is recorded in the JSON.

Test-only reduced method settings (``TEST_OVERRIDES``) are used only when
``EXP3_TEST_MODE=1`` AND the dataset is synthetic (sha differs from the real one and all
cases are ``synth_*``); with the real dataset that variable makes every run refuse.

Usage (anemoi)::

    run_protocol.py --freeze-config [--out exp3_melt_emulator/figures/protocol_config.json]
    run_protocol.py --cv
    run_protocol.py --final
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import heapq
import json
import multiprocessing
import os
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp3_melt_emulator.data import DATASET_PATH, DATASET_SHA256, load_runs  # noqa: E402
from exp3_melt_emulator.folds import drop_n1, make_folds  # noqa: E402
from exp3_melt_emulator.methods import REGISTRY  # noqa: E402
from exp3_melt_emulator.metrics import (  # noqa: E402
    coverage, crps_gaussian, integrated_rel_error, rmse_gl, rmse_map, skill,
)

FIGURES = REPO / "exp3_melt_emulator" / "figures"
DEFAULT_CONFIG = FIGURES / "protocol_config.json"
DEFAULT_NOISE = FIGURES / "noise_floor.json"

# PD0 (notes/pd0_preregistro.md), approved by Esteban on 2026-11-10.
PD0_THRESHOLDS = {
    "gl_frac": 0.10,
    "hypothesis_folds": "H3a-H3c use the main folds only; drop_n1 results are sensitivity",
    "H3a": {"cases": "final test (3 runs, Tbot = 2.0)",
            "metric": "|integrated_rel_error| per run",
            "rule": "mean over the 3 runs: M7 < M6, and M7 < M6 in >= 2 of the 3 runs",
            "min_runs_better": 2,
            "inconclusive_if": "|mean difference| < noise_floor.median_abs_integrated_rel_diff"},
    "H3b": {"cases": "leave-one-run-out (25 main folds, not the drop_n1 sensitivity folds)", "metric": "mean rmse_map over folds",
            "rule": "mean_M7 <= max_ratio_M7_over_M6 * mean_M6", "max_ratio_M7_over_M6": 1.10},
    "H3c": {"cases": "zc_out + warm_out test runs + final test runs, for M6 and M7",
            "spread": "std (ddof = 0) over the 5 members of the integrated melt / "
                      "|integrated melt of the ensemble mean| (Ockenden et al. 2026)",
            "rule": "no case with spread < spread_max AND |integrated_rel_error| > error_min",
            "spread_max": 0.05, "error_min": 0.20,
            "ockenden_flag": 0.10,
            "ockenden_report": "confusion table spread > 10 % vs |integrated_rel_error| > 20 %"},
    "noise_floor_rule": "a difference between methods smaller than the pooled noise RMS "
                        "(rmse_map, rmse_gl) or than the median |N1 drift| "
                        "(integrated_rel_error) is not interpreted",
    "ranking": "reported, not hypothesised",
    "stop_rule": "--cv: if the wall-clock estimate (<= 5 concurrent processes, 1 thread "
                 "each, measured after the first batch of folds) exceeds 8 h, stop and report",
    "stop_rule_note": "limit raised from the plan's 6 h to 8 h by Esteban (2026-11-10), before "
                      "any protocol run, because full-LORO tf_kind selection puts the estimate "
                      "at about 6 h",
}

# Methods fitted per fold kind: M4 and M5 see only (Tbot, z_c, s) and cannot represent the
# d = 720 m profile of ocean0_4km; M0 is kept there because skill is relative to M0.
PROFILE_TEST_METHODS = ("M0", "M1", "M2", "M3", "M6", "M7")
# Paths that may differ from HEAD when the configuration is frozen (the PD0 note is edited
# with Esteban's decisions in the same Step-5 commit).
ALLOWED_DIRTY_PREFIX = "exp3_melt_emulator/notes/"


def git_state() -> dict:
    """HEAD commit and the tracked paths that differ from it (empty = clean)."""
    def git(*a):
        return subprocess.run(["git", *a], cwd=REPO, check=True, capture_output=True,
                              text=True).stdout
    commit = git("rev-parse", "HEAD").strip()
    dirty = [ln[3:] for ln in git("status", "--porcelain", "--untracked-files=no").splitlines()]
    return {"commit": commit, "dirty_paths": dirty}


def _sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _fold(f) -> dict:
    return {"name": f.name, "kind": f.kind, "train": list(f.train), "test": list(f.test)}


def _noise(noise_json) -> dict | None:
    if not noise_json:
        return None
    rec = json.loads(Path(noise_json).read_text())         # FileNotFoundError if missing
    return {"json": str(noise_json), "json_sha256": _sha256(noise_json),
            "npz_sha256": rec["sha256"], "pooled_median_rms": rec["pooled_median_rms"],
            "median_abs_integrated_rel_diff": rec["median_abs_integrated_rel_diff"]}


def freeze_config(out=DEFAULT_CONFIG, dataset_path=DATASET_PATH, dataset_sha=DATASET_SHA256,
                  noise_json=DEFAULT_NOISE, git_state=git_state) -> dict:
    """Write the frozen protocol configuration to ``out``; FileExistsError if it exists."""
    out = Path(out)
    if out.exists():
        raise FileExistsError(f"{out} exists: the protocol configuration is frozen")
    git = git_state()
    bad = [p for p in git["dirty_paths"] if not p.startswith(ALLOWED_DIRTY_PREFIX)]
    if bad:
        raise RuntimeError(f"uncommitted changes outside {ALLOWED_DIRTY_PREFIX}: {bad}; "
                           "commit them before freezing the protocol")
    runs = load_runs(dataset_path, dataset_sha)            # ValueError on sha mismatch
    folds = make_folds(runs)
    cfg = {
        "dataset": {"path": str(dataset_path), "sha256": dataset_sha},
        "git": git,
        "folds": [_fold(f) for f in folds],
        "folds_sensitivity_n1": [_fold(drop_n1(f, runs)) for f in folds if f.kind != "final"],
        "n1_fail_cases": [r.case for r in runs if not r.n1_pass],
        "profile_test_methods": list(PROFILE_TEST_METHODS),
        "methods": {mid: REGISTRY[mid]().hparams() for mid in sorted(REGISTRY)},
        "thresholds": PD0_THRESHOLDS,
        "noise_floor": _noise(noise_json),
    }
    text = json.dumps(cfg, indent=1) + "\n"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "x") as fh:                             # exclusive create: never overwrite
        fh.write(text)
    return cfg




# =========================================================================================
# Task 8: --cv and --final
# =========================================================================================

# sha256 of figures/protocol_config.json as approved at PD0 (notes/pd0_preregistro.md).
CONFIG_SHA256 = "56a9ffe5b57fd69da0462e63cb79cb1cbc7538e8444946e9bde5101f9a462486"
RUNS_DIR = Path("/Volumes/LaCie/portfolio_runs/exp3")
UNITS_DIR = RUNS_DIR / "cv_units"
CV_PROGRESS_LOG = RUNS_DIR / "cv_progress.log"
CV_RESULTS = FIGURES / "cv_results.json"
FINAL_TEST = FIGURES / "final_test.json"
FINAL_MAPS = RUNS_DIR / "final_maps.npz"
FINAL_MARKER = RUNS_DIR / "final_started.json"
FINAL_PROGRESS_LOG = RUNS_DIR / "final_progress.log"

MAX_WORKERS = 5            # P-cores only (repo rules; PD0 stop rule: <= 5 processes, 1 thread)
STOP_HOURS = 8.0           # PD0 decision 7
PROB_METHODS = ("M5", "M6", "M7")
MEMBER_METHODS = ("M6", "M7")
COVERAGE_LEVELS = (0.68, 0.95)
THREAD_VARS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "VECLIB_NUM_THREADS", "NUMEXPR_NUM_THREADS")
EXIT_STOPPED = 3

# Test-only reduced settings (EXP3_TEST_MODE=1 on a synthetic dataset only; test_overrides).
TEST_OVERRIDES = {
    "M6": {"n_members": 2, "seeds": [0, 1], "max_epochs": 3, "patience": 2},
    "M7": {"n_members": 2, "seeds": [0, 1], "max_epochs": 3, "patience": 2},
}

# Prior wall time (s, 1 thread) per unit for the stop-rule estimate, from the Task 6 timings
# quoted in the PD0 note ("Computo estimado"): one M7 fit with tf_kind fixed ~50 s and
# 2m + 1 fits with the full-LORO tf_kind choice; M6 ~66 s at m = 24; M4/M5 seconds.
PRIOR_M7_FIT_S = 50.0
PRIOR_S = {"M6": 66.0, "M5": 30.0, "M4": 30.0, "M0": 2.0, "M1": 2.0, "M2": 2.0, "M3": 2.0}


class ProtocolError(RuntimeError):
    """The frozen protocol cannot be run as requested (refusal)."""


def _now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


def _jsonable(x):
    """Plain-JSON copy (numpy scalars/arrays -> Python, dict keys -> str, tuples -> lists)."""
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, np.generic):
        return x.item()
    return x


def _atomic_write(path, data: bytes):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    with open(tmp, "wb") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def _log(path, line: str):
    line = f"{_now()} {line}"
    print(line, flush=True)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as fh:
        fh.write(line + "\n")


# --- verification of the frozen configuration ---------------------------------------------
def load_protocol(config_path, expected_sha):
    """Frozen config + runs, after every check; ProtocolError (ValueError for the dataset
    sha, from load_runs) otherwise. Returns (cfg, runs_by_case)."""
    raw = Path(config_path).read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    if got != expected_sha:
        raise ProtocolError(f"config sha256 mismatch for {config_path}: {got} != {expected_sha}")
    cfg = json.loads(raw)
    if set(cfg["methods"]) != set(REGISTRY):
        raise ProtocolError(f"methods: frozen {sorted(cfg['methods'])} != code {sorted(REGISTRY)}")
    for mid in sorted(REGISTRY):
        cur = json.loads(json.dumps(REGISTRY[mid]().hparams()))
        if cur != cfg["methods"][mid]:
            diff = sorted(k for k in set(cur) | set(cfg["methods"][mid])
                          if cur.get(k) != cfg["methods"][mid].get(k))
            raise ProtocolError(f"{mid}: current hparams differ from the frozen ones in {diff}")
    runs = load_runs(cfg["dataset"]["path"], cfg["dataset"]["sha256"])   # ValueError on sha
    folds = make_folds(runs)
    if cfg["folds"] != [_fold(f) for f in folds]:
        raise ProtocolError("folds recomputed from the dataset differ from the frozen folds")
    if cfg["folds_sensitivity_n1"] != [_fold(drop_n1(f, runs)) for f in folds if f.kind != "final"]:
        raise ProtocolError("N1-sensitivity folds recomputed from the dataset differ from the frozen folds")
    return cfg, {r.case: r for r in runs}


def test_overrides(cfg, runs) -> dict:
    """Reduced method settings for the test suite; {} unless EXP3_TEST_MODE=1, and refused
    (ProtocolError) when that variable is set on anything but a synthetic dataset."""
    if os.environ.get("EXP3_TEST_MODE") != "1":
        return {}
    cases = [r.case for r in (runs.values() if isinstance(runs, dict) else runs)]
    if cfg["dataset"]["sha256"] == DATASET_SHA256 or not all(c.startswith("synth_") for c in cases):
        raise ProtocolError("EXP3_TEST_MODE=1 is refused: the dataset is not the synthetic test set")
    return json.loads(json.dumps(TEST_OVERRIDES))


test_overrides.__test__ = False      # not a pytest test function


def make_method(mid: str, overrides: dict):
    kw = overrides.get(mid)
    if not kw:
        return REGISTRY[mid]()
    from exp3_melt_emulator.methods import nets
    cls = {"M6": nets.DeepEnsemble, "M7": nets.Hybrid}[mid]
    return cls(**kw)


def _thread_env():
    for v in THREAD_VARS:
        os.environ[v] = "1"


def _worker_init():
    _thread_env()
    try:
        import torch
        torch.set_num_threads(1)
    except ImportError:                     # pragma: no cover - anemoi has torch
        pass


def _pool(workers):
    if not 1 <= workers <= MAX_WORKERS:
        raise ProtocolError(f"workers must be in 1..{MAX_WORKERS}, got {workers}")
    _thread_env()                           # inherited by the spawned workers
    return ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
                               initializer=_worker_init)


# --- work units -----------------------------------------------------------------------------
def _uid(variant, fold, method):
    return f"{variant}__{fold}__{method}"


def prior_seconds(unit) -> float:
    """Prior single-thread wall time of a unit (scaled by the training-run count m)."""
    m = unit["n_train"]
    if unit["method"] == "M7":
        return (2 * m + 1) * PRIOR_M7_FIT_S
    if unit["method"] == "M6":
        return PRIOR_S["M6"] * m / 24
    return PRIOR_S[unit["method"]]


def cv_units(cfg, runs) -> list[dict]:
    """All --cv work units: every main-variant unit before any n1 unit (PD0 Enmienda 1);
    within each variant M7 first, then the rest, each group longest (prior) first."""
    units = []
    for variant, folds in (("main", cfg["folds"]), ("n1", cfg["folds_sensitivity_n1"])):
        for f in folds:
            if f["kind"] == "final":
                continue
            for c in f["test"]:
                if runs[c].split == "test":
                    raise ProtocolError(f"fold {f['name']}: test run {c} has split == 'test'")
            for c in f["train"]:
                if runs[c].split == "test" or runs[c].group == "reference":
                    raise ProtocolError(f"fold {f['name']}: {c} may not train")
            if not f["test"]:
                raise ProtocolError(f"fold {f['name']} has no test runs")
            for mid in sorted(cfg["methods"]):
                if f["kind"] == "profile_test" and mid not in cfg["profile_test_methods"]:
                    continue
                units.append({"unit": _uid(variant, f["name"], mid), "variant": variant,
                              "fold": f["name"], "kind": f["kind"], "method": mid,
                              "train": list(f["train"]), "test": list(f["test"]),
                              "n_train": len(f["train"])})
    # Scheduling order only: every unit is an independent, deterministic fit, so the order
    # changes when results become available (main folds first, PD0 Enmienda 1), never
    # the results themselves.
    units.sort(key=lambda u: (u["variant"] != "main", u["method"] != "M7", -prior_seconds(u)))
    return units


def _metrics(pred, run, gl_frac):
    msk = run.shelf_mask
    return {"rmse_map": rmse_map(pred, run), "integrated_rel_error": integrated_rel_error(pred, run),
            "rmse_gl": rmse_gl(pred, run, gl_frac),
            "integrated_pred": float(np.asarray(pred, float)[msk].sum()),
            "integrated_true": float(run.melt[msk].sum())}


def _method_info(m) -> dict:
    info = dict(getattr(m, "info", {}) or {})
    if m.name in ("M1", "M2", "M3"):
        info["coef"] = m.coef
    if m.name == "M4":
        info.update(degree=m.degree, loro_rmse=m.loro_rmse)
    if m.name == "M7":
        info.setdefault("tf_kind", m.tf_kind_)
    return _jsonable(info)


def run_unit(unit, train_runs, test_runs, overrides, gl_frac, config_sha, units_dir=None,
             in_sample_runs=(), return_maps=False, amendment=None):
    """Fit one method on ``train_runs`` and evaluate ``test_runs`` (worker entry point).
    Writes ``units_dir/<unit>.json`` atomically when ``units_dir`` is given."""
    t0 = time.time()
    started = _now()
    m = make_method(unit["method"], overrides)
    m.fit(list(train_runs))
    t_fit = time.time() - t0
    tests, maps = [], {}
    for run in test_runs:
        if unit["method"] in PROB_METHODS:
            mu, sigma = m.predict_dist(run)
        else:
            mu, sigma = m.predict(run), None
        rec = {"case": run.case, "group": run.group, "split": run.split, "Tbot": run.Tbot,
               "z_c": run.z_c, "s": run.s, "n1_pass": run.n1_pass,
               "metrics": _metrics(mu, run, gl_frac), "prob": None, "members": None}
        if sigma is not None:
            rec["prob"] = {f"coverage_{int(round(100 * lv))}": coverage(mu, sigma, run, lv)
                           for lv in COVERAGE_LEVELS}
            rec["prob"]["crps"] = crps_gaussian(mu, sigma, run)
            rec["prob"]["mean_sigma"] = float(np.asarray(sigma)[run.shelf_mask].mean())
            maps[f"{run.case}__sigma"] = sigma
        maps[f"{run.case}__mean"] = mu
        if unit["method"] in MEMBER_METHODS:
            pm = m.predict_members(run)
            integ = np.asarray(pm["integrated"])
            rec["members"] = {
                "integrated": pm["integrated"], "integrated_mixture": pm["integrated_mixture"],
                "mean_aleatoric": pm["mean_aleatoric"], "mean_epistemic": pm["mean_epistemic"],
                # H3c (PD0): std (ddof = 0) of the member integrated melt / |mixture|
                "spread_rel": float(integ.std() / abs(pm["integrated_mixture"])),
            }
            maps[f"{run.case}__members_mu"] = pm["mu"]
        tests.append(rec)
    for run in in_sample_runs:
        maps[f"{run.case}__mean_in_sample"] = m.predict(run)
    out = {"unit": unit["unit"], "variant": unit["variant"], "fold": unit["fold"],
           "kind": unit["kind"], "method": unit["method"], "train": unit["train"],
           "test": unit["test"], "config_sha256": config_sha, "test_overrides": overrides,
           "amendment": amendment,
           "started": started, "finished": _now(), "fit_s": t_fit, "wall_s": time.time() - t0,
           "hparams": _jsonable(m.hparams()), "info": _method_info(m), "results": _jsonable(tests)}
    if units_dir is not None:
        _atomic_write(Path(units_dir) / f"{unit['unit']}.json",
                      (json.dumps(out, indent=1) + "\n").encode())
    return (out, maps) if return_maps else out


# --- stop-rule estimate ---------------------------------------------------------------------
def _lpt_makespan(busy, durations, workers):
    heap = sorted(list(busy) + [0.0] * (workers - len(busy)))
    heapq.heapify(heap)
    for d in durations:                       # list order (already longest first)
        heapq.heappush(heap, heapq.heappop(heap) + d)
    return max(heap) if heap else 0.0


def estimate_wall(completed, in_flight, pending, elapsed_s, workers) -> dict:
    """Total wall-clock estimate (h) of this session: elapsed + simulated list-scheduling
    makespan of the in-flight and pending units on ``workers`` processes.

    Unit durations: priors (``prior_seconds``) scaled per method by the measured/prior ratio
    of the completed units (unmeasured methods keep the prior). For M7, two cost models are
    calibrated on the completed M7 units: 'fits' (~ 2m + 1 fits, the PD0 note's model) and
    'cells' (~ m (2m - 1): fits x training cells); the larger (conservative) estimate is the
    one compared with the limit.
    ``completed``: [(unit, wall_s)]; ``in_flight``: [(unit, running_s)]; ``pending``: [unit]."""
    def m7_cost(u, model):
        m = u["n_train"]
        return 2 * m + 1 if model == "fits" else m * (2 * m - 1)

    out = {"elapsed_h": elapsed_s / 3600, "workers": workers, "n_completed": len(completed),
           "n_in_flight": len(in_flight), "n_pending": len(pending)}
    ests = {}
    for model in ("fits", "cells"):
        ratio = {}
        for mid in {u["method"] for u, _ in completed}:
            done = [(u, w) for u, w in completed if u["method"] == mid]
            if mid == "M7":
                ratio[mid] = np.mean([w / m7_cost(u, model) for u, w in done])
            else:
                ratio[mid] = np.mean([w / prior_seconds(u) for u, w in done])

        def dur(u):
            if u["method"] == "M7":
                c = ratio.get("M7", PRIOR_M7_FIT_S if model == "fits" else None)
                if c is None:                 # no M7 measured: prior via the fits model
                    return prior_seconds(u)
                return c * m7_cost(u, model)
            return ratio.get(u["method"], 1.0) * prior_seconds(u)

        busy = [max(dur(u) - r, 0.0) for u, r in in_flight]
        span = _lpt_makespan(busy, [dur(u) for u in pending], workers)
        ests[model] = (elapsed_s + span) / 3600
        if model == "fits" and "M7" in ratio:
            out["m7_s_per_fit"] = float(ratio["M7"])
    out["estimate_h_by_model"] = ests
    out["estimate_h"] = float(max(ests.values()))
    return out


# --- --cv -----------------------------------------------------------------------------------
def _load_unit(path, uid, sha, overrides):
    rec = json.loads(Path(path).read_text())
    if rec.get("unit") != uid or rec.get("config_sha256") != sha:
        raise ProtocolError(f"{path}: written with another config (sha "
                            f"{rec.get('config_sha256')}) or unit; refusing to mix")
    if rec.get("test_overrides") != overrides:
        raise ProtocolError(f"{path}: written with other test overrides {rec.get('test_overrides')}")
    return rec


def assemble_cv(cfg, units, units_dir, sha, overrides, git, amendments=(), wall_limit_hours=STOP_HOURS) -> dict:
    """cv_results.json content from the unit files (skill against M0 of the same
    variant/fold; rmse_gl per run as in R8, the RMS over runs is taken when aggregating)."""
    recs = [_load_unit(Path(units_dir) / f"{u['unit']}.json", u["unit"], sha, overrides)
            for u in units]
    m0 = {(r["variant"], r["fold"], t["case"]): t["metrics"]
          for r in recs if r["method"] == "M0" for t in r["results"]}
    entries, unit_recs = [], []
    for r in recs:
        for t in r["results"]:
            ref = m0[(r["variant"], r["fold"], t["case"])]
            mt = t["metrics"]
            entries.append({
                "variant": r["variant"], "label": "sensitivity" if r["variant"] == "n1" else "main",
                "fold": r["fold"], "kind": r["kind"], "method": r["method"], "unit": r["unit"],
                "case": t["case"], "group": t["group"], "Tbot": t["Tbot"], "z_c": t["z_c"],
                "s": t["s"], "n1_pass": t["n1_pass"], "metrics": mt,
                "skill": {"rmse_map": skill(mt["rmse_map"], ref["rmse_map"]),
                          "rmse_gl": skill(mt["rmse_gl"], ref["rmse_gl"]),
                          "abs_integrated_rel_error": skill(abs(mt["integrated_rel_error"]),
                                                            abs(ref["integrated_rel_error"]))},
                "prob": t["prob"], "members": t["members"]})
        unit_recs.append({k: r[k] for k in ("unit", "variant", "fold", "kind", "method", "train",
                                            "test", "started", "finished", "fit_s", "wall_s",
                                            "info")})
    seen = list(amendments)
    for r in recs:                       # amendments under which any unit was written
        if r.get("amendment") and r["amendment"] not in seen:
            seen.append(r["amendment"])
    return {"config_sha256": sha, "dataset_sha256": cfg["dataset"]["sha256"],
            "created": _now(), "git": git, "test_overrides": overrides,
            "amendments": sorted(seen, key=lambda a: a["timestamp"]),
            "wall_limit_hours": wall_limit_hours,
            "gl_frac": cfg["thresholds"]["gl_frac"],
            "cpu_s_total": float(sum(u["wall_s"] for u in unit_recs)),
            "n1_fail_cases": cfg["n1_fail_cases"],
            "notes": "variant n1 = drop_n1 sensitivity repeat (never used for H3a-H3c); "
                     "skill = 1 - err / err_M0 of the same variant and fold; rmse_gl is the "
                     "per-run |GL-mean difference| (R8)",
            "entries": entries, "units": unit_recs}


def run_cv(config_path=DEFAULT_CONFIG, expected_sha=CONFIG_SHA256, units_dir=UNITS_DIR,
           progress_log=CV_PROGRESS_LOG, out=CV_RESULTS, workers=MAX_WORKERS,
           estimate_fn=estimate_wall, git_state=git_state, wall_limit_hours=None,
           amendment=None) -> dict:
    """Run (or resume) --cv; returns {"status": "complete"|"stopped", "ran": [...], ...}.

    The frozen stop rule uses STOP_HOURS (8 h). ``wall_limit_hours`` replaces it only
    together with ``amendment`` (the text of a documented PD0 amendment); the override is
    recorded in the progress log, in every unit written under it and in cv_results.json."""
    if (wall_limit_hours is None) != (amendment is None):
        raise ProtocolError("--wall-limit-hours and --amendment must be given together")
    am = None
    if amendment is not None:
        if not str(amendment).strip():
            raise ProtocolError("the amendment text is empty")
        if not float(wall_limit_hours) > 0:
            raise ProtocolError(f"wall_limit_hours must be > 0, got {wall_limit_hours}")
        am = {"text": str(amendment), "timestamp": _now(),
              "wall_limit_hours": float(wall_limit_hours),
              "frozen_wall_limit_hours": STOP_HOURS}
    stop_hours = am["wall_limit_hours"] if am else STOP_HOURS
    cfg, runs = load_protocol(config_path, expected_sha)
    overrides = test_overrides(cfg, runs)
    git = git_state()
    units = cv_units(cfg, runs)
    units_dir = Path(units_dir)
    done = {u["unit"] for u in units if (units_dir / f"{u['unit']}.json").exists()
            and _load_unit(units_dir / f"{u['unit']}.json", u["unit"], expected_sha, overrides)}
    pending = [u for u in units if u["unit"] not in done]
    gl_frac = cfg["thresholds"]["gl_frac"]
    _log(progress_log, f"START --cv config {expected_sha[:12]} commit {git['commit'][:7]} "
                       f"units {len(units)} done {len(done)} pending {len(pending)} "
                       f"workers {workers} stop limit {stop_hours:g} h"
                       + (" TEST_OVERRIDES" if overrides else ""))
    if am:
        _log(progress_log, f"AMENDMENT wall-clock limit {am['wall_limit_hours']:g} h replaces the "
                           f"frozen {STOP_HOURS:g} h for this session: {am['text']!r} "
                           f"(recorded {am['timestamp']})")
    t_start = time.time()
    ran, completed, failures = [], [], []
    estimate, stopped = None, False
    if pending:
        with _pool(min(workers, len(pending))) as pool:
            in_flight = {}
            first_batch = set()

            def submit():
                u = pending.pop(0)
                fut = pool.submit(run_unit, u, [runs[c] for c in u["train"]],
                                  [runs[c] for c in u["test"]], overrides, gl_frac,
                                  expected_sha, str(units_dir), (), False, am)
                in_flight[fut] = (u, time.time())
                return u

            while pending and len(in_flight) < workers:
                first_batch.add(submit()["unit"])
            checked = False
            while in_flight:
                fin, _ = wait(list(in_flight), return_when=FIRST_COMPLETED)
                for fut in fin:
                    u, t0 = in_flight.pop(fut)
                    first_batch.discard(u["unit"])
                    try:
                        rec = fut.result()
                    except Exception as exc:     # logged; no new units; raised at the end
                        failures.append((u["unit"], repr(exc)))
                        _log(progress_log, f"FAILED {u['unit']} after {time.time() - t0:.1f} s: {exc!r}")
                        continue
                    ran.append(u["unit"])
                    completed.append((u, rec["wall_s"]))
                    _log(progress_log, f"DONE {u['unit']} wall {rec['wall_s']:.1f} s "
                                       f"(fit {rec['fit_s']:.1f} s; n_train {u['n_train']}; "
                                       f"{len(done) + len(ran)}/{len(units)})")
                if not checked and not first_batch and not failures:
                    checked = True
                    now = time.time()
                    estimate = estimate_fn(
                        completed=list(completed),
                        in_flight=[(u, now - t0) for u, t0 in in_flight.values()],
                        pending=list(pending), elapsed_s=now - t_start, workers=workers)
                    _log(progress_log, f"ESTIMATE total wall clock {estimate['estimate_h']:.2f} h "
                                       f"(limit {stop_hours:g} h): {json.dumps(_jsonable(estimate))}")
                    if estimate["estimate_h"] > stop_hours:
                        stopped = True
                        _log(progress_log, f"STOP: estimate {estimate['estimate_h']:.2f} h > "
                                           f"{stop_hours:g} h; no new units; waiting for "
                                           f"{len(in_flight)} in flight")
                while pending and len(in_flight) < workers and not stopped and not failures:
                    submit()
    if failures:
        raise RuntimeError(f"--cv: {len(failures)} unit(s) failed: {failures}")
    if stopped or pending:
        _log(progress_log, f"END stopped by the stop rule; {len(units) - len(pending)}/{len(units)} "
                           f"units written; no cv_results.json")
        return {"status": "stopped", "ran": ran, "estimate": estimate}
    res = assemble_cv(cfg, units, units_dir, expected_sha, overrides, git,
                      [am] if am else [], stop_hours)
    res["session"] = {"wall_s": time.time() - t_start, "units_run": len(ran),
                      "units_resumed": len(done), "workers": workers, "estimate": estimate}
    _atomic_write(out, (json.dumps(res, indent=1) + "\n").encode())
    _log(progress_log, f"END complete: wrote {out} ({len(res['entries'])} entries) in "
                       f"{time.time() - t_start:.1f} s")
    return {"status": "complete", "ran": ran, "estimate": estimate, "out": str(out)}


# --- --final --------------------------------------------------------------------------------
def run_final(config_path=DEFAULT_CONFIG, expected_sha=CONFIG_SHA256, out=FINAL_TEST,
              maps_out=FINAL_MAPS, marker=FINAL_MARKER, progress_log=FINAL_PROGRESS_LOG,
              workers=MAX_WORKERS, git_state=git_state) -> dict:
    """The one-time final test (no override exists)."""
    out, marker = Path(out), Path(marker)
    if out.exists():
        raise FileExistsError(f"{out} exists: the final test has been evaluated")
    if marker.exists():
        raise ProtocolError(f"{marker} exists: a --final run started before "
                            f"({marker.read_text().strip()}); it is not repeated")
    cfg, runs = load_protocol(config_path, expected_sha)
    overrides = test_overrides(cfg, runs)
    git = git_state()
    (fold,) = [f for f in cfg["folds"] if f["kind"] == "final"]
    if not fold["test"] or any(runs[c].split != "test" for c in fold["test"]):
        raise ProtocolError("final fold: the test runs must be the split == 'test' runs")
    marker.parent.mkdir(parents=True, exist_ok=True)
    with open(marker, "x") as fh:            # from here on the final test counts as used
        fh.write(json.dumps({"started": _now(), "config_sha256": expected_sha,
                             "commit": git["commit"]}) + "\n")
    _log(progress_log, f"START --final config {expected_sha[:12]} commit {git['commit'][:7]}")
    t_start = time.time()
    gl_frac = cfg["thresholds"]["gl_frac"]
    train = [runs[c] for c in fold["train"]]
    test = [runs[c] for c in fold["test"]]
    units = [{"unit": _uid("main", "final", mid), "variant": "main", "fold": "final",
              "kind": "final", "method": mid, "train": list(fold["train"]),
              "test": list(fold["test"]), "n_train": len(train)} for mid in sorted(cfg["methods"])]
    units.sort(key=lambda u: (u["method"] != "M7", -prior_seconds(u)))
    recs, maps = {}, {}
    with _pool(min(workers, len(units))) as pool:
        futs = {pool.submit(run_unit, u, train, test, overrides, gl_frac, expected_sha, None,
                            train, True): u for u in units}
        for fut in futs:
            u = futs[fut]
            rec, mp = fut.result()
            recs[u["method"]] = rec
            maps.update({f"{u['method']}__{k}": np.asarray(v, float) for k, v in mp.items()})
            _log(progress_log, f"DONE {u['unit']} wall {rec['wall_s']:.1f} s")
    maps_out = Path(maps_out)
    maps_out.parent.mkdir(parents=True, exist_ok=True)
    tmp = maps_out.with_name(f"{maps_out.stem}.tmp.{os.getpid()}.npz")
    np.savez(tmp, **maps)
    os.replace(tmp, maps_out)
    maps_sha = _sha256(maps_out)
    ref = {t["case"]: t["metrics"] for t in recs["M0"]["results"]}
    entries = []
    for mid in sorted(recs):
        for t in recs[mid]["results"]:
            mt, r0 = t["metrics"], ref[t["case"]]
            entries.append({**{k: t[k] for k in ("case", "group", "Tbot", "z_c", "s", "n1_pass")},
                            "fold": "final", "method": mid, "metrics": mt,
                            "skill": {"rmse_map": skill(mt["rmse_map"], r0["rmse_map"]),
                                      "rmse_gl": skill(mt["rmse_gl"], r0["rmse_gl"]),
                                      "abs_integrated_rel_error": skill(
                                          abs(mt["integrated_rel_error"]),
                                          abs(r0["integrated_rel_error"]))},
                            "prob": t["prob"], "members": t["members"]})
    res = {"config_sha256": expected_sha, "dataset_sha256": cfg["dataset"]["sha256"],
           "created": _now(), "git": git, "test_overrides": overrides, "gl_frac": gl_frac,
           "train": fold["train"], "test": fold["test"],
           "maps": str(maps_out), "maps_sha256": maps_sha,
           "maps_keys": "<method>__<case>__mean | __sigma (M5-M7) | __members_mu (M6, M7) for "
                        "the test runs; <method>__<case>__mean_in_sample for the training runs",
           "wall_s": time.time() - t_start, "entries": entries,
           "units": [{k: recs[m][k] for k in ("unit", "method", "started", "finished", "fit_s",
                                              "wall_s", "info")} for m in sorted(recs)]}
    res = _jsonable(res)
    with open(out, "x") as fh:
        fh.write(json.dumps(res, indent=1) + "\n")
    _log(progress_log, f"END wrote {out} and {maps_out} (sha {maps_sha[:12]}) in "
                       f"{time.time() - t_start:.1f} s")
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    act = ap.add_mutually_exclusive_group(required=True)
    act.add_argument("--freeze-config", action="store_true",
                     help="write the frozen PD0 configuration (refuses to overwrite)")
    act.add_argument("--cv", action="store_true",
                     help="cross-validation on every non-final fold (resumable)")
    act.add_argument("--final", action="store_true",
                     help="the one-time final test (refuses if it was run before)")
    ap.add_argument("--wall-limit-hours", type=float, default=None,
                    help="--cv only, with --amendment: replaces the frozen 8 h stop limit")
    ap.add_argument("--amendment", default=None,
                    help="--cv only, with --wall-limit-hours: text of the PD0 amendment")
    ap.add_argument("--out", type=Path, default=DEFAULT_CONFIG,
                    help="--freeze-config output path")
    ap.add_argument("--dataset", type=Path, default=DATASET_PATH)
    ap.add_argument("--dataset-sha", default=DATASET_SHA256)
    ap.add_argument("--noise-json", default=str(DEFAULT_NOISE),
                    help="noise-floor summary to record ('' = none)")
    a = ap.parse_args(argv)
    if (a.wall_limit_hours is None) != (a.amendment is None):
        ap.error("--wall-limit-hours and --amendment must be given together")
    if a.amendment is not None and not a.cv:
        ap.error("--wall-limit-hours/--amendment apply to --cv only")
    if a.cv:
        st = run_cv(DEFAULT_CONFIG, CONFIG_SHA256, wall_limit_hours=a.wall_limit_hours,
                    amendment=a.amendment)
        return 0 if st["status"] == "complete" else EXIT_STOPPED
    if a.final:
        res = run_final(DEFAULT_CONFIG, CONFIG_SHA256)
        print(f"wrote {FINAL_TEST}: {len(res['entries'])} entries")
        return 0
    cfg = freeze_config(a.out, a.dataset, a.dataset_sha, a.noise_json, git_state=git_state)
    print(f"wrote {a.out}: {len(cfg['folds'])} folds, methods {sorted(cfg['methods'])}, "
          f"commit {cfg['git']['commit'][:7]}"
          + (f" (dirty: {cfg['git']['dirty_paths']})" if cfg['git']['dirty_paths'] else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
