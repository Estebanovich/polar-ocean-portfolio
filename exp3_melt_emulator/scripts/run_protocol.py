#!/opt/anaconda3/envs/anemoi/bin/python
r"""Exp. 3 validation protocol (notes/pd0_preregistro.md).

This task implements only ``--freeze-config``: it writes ``figures/protocol_config.json``
with everything the protocol runs (``--cv``, ``--final``, later tasks) must use, fixed at
PD0 before any protocol run:

* the dataset sha256 (checked against the file before anything is read);
* the folds of ``folds.make_folds`` and the N1-sensitivity variants (``folds.drop_n1``) of
  every fold except ``final`` (the final test is evaluated once, without a repeat);
* every method's hyperparameters, read from the constructor defaults (``hparams()`` of
  each ``methods.REGISTRY`` factory's instance);
* the PD0 thresholds (``PD0_THRESHOLDS``) and the noise floor (figures/noise_floor.json);
* the git commit and the paths that differed from it at freeze time.

It refuses to overwrite an existing file; there is no override.

Usage (anemoi)::

    run_protocol.py --freeze-config [--out exp3_melt_emulator/figures/protocol_config.json]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp3_melt_emulator.data import DATASET_PATH, DATASET_SHA256, load_runs  # noqa: E402
from exp3_melt_emulator.folds import drop_n1, make_folds  # noqa: E402
from exp3_melt_emulator.methods import REGISTRY  # noqa: E402

FIGURES = REPO / "exp3_melt_emulator" / "figures"
DEFAULT_CONFIG = FIGURES / "protocol_config.json"
DEFAULT_NOISE = FIGURES / "noise_floor.json"

# PD0 (notes/pd0_preregistro.md). PROPOSED values until Esteban approves PD0; the note and
# this dict must agree when --freeze-config is run.
PD0_THRESHOLDS = {
    "gl_frac": 0.10,
    "H3a": {"cases": "final test (3 runs, Tbot = 2.0)",
            "metric": "|integrated_rel_error| per run",
            "rule": "mean over the 3 runs: M7 < M6, and M7 < M6 in >= 2 of the 3 runs",
            "min_runs_better": 2,
            "inconclusive_if": "|mean difference| < noise_floor.median_abs_integrated_rel_diff"},
    "H3b": {"cases": "leave-one-run-out (25 folds)", "metric": "mean rmse_map over folds",
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
}


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
    runs = load_runs(dataset_path, dataset_sha)            # ValueError on sha mismatch
    folds = make_folds(runs)
    cfg = {
        "dataset": {"path": str(dataset_path), "sha256": dataset_sha},
        "git": git_state(),
        "folds": [_fold(f) for f in folds],
        "folds_sensitivity_n1": [_fold(drop_n1(f, runs)) for f in folds if f.kind != "final"],
        "n1_fail_cases": [r.case for r in runs if not r.n1_pass],
        "methods": {mid: REGISTRY[mid]().hparams() for mid in sorted(REGISTRY)},
        "thresholds": PD0_THRESHOLDS,
        "noise_floor": _noise(noise_json),
    }
    text = json.dumps(cfg, indent=1) + "\n"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "x") as fh:                             # exclusive create: never overwrite
        fh.write(text)
    return cfg


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--freeze-config", action="store_true",
                    help="write the frozen PD0 configuration (refuses to overwrite)")
    ap.add_argument("--out", type=Path, default=DEFAULT_CONFIG)
    ap.add_argument("--dataset", type=Path, default=DATASET_PATH)
    ap.add_argument("--dataset-sha", default=DATASET_SHA256)
    ap.add_argument("--noise-json", default=str(DEFAULT_NOISE),
                    help="noise-floor summary to record ('' = none)")
    a = ap.parse_args(argv)
    if not a.freeze_config:
        ap.error("nothing to do: only --freeze-config is implemented")
    cfg = freeze_config(a.out, a.dataset, a.dataset_sha, a.noise_json, git_state=git_state)
    print(f"wrote {a.out}: {len(cfg['folds'])} folds, methods {sorted(cfg['methods'])}, "
          f"commit {cfg['git']['commit'][:7]}"
          + (f" (dirty: {cfg['git']['dirty_paths']})" if cfg['git']['dirty_paths'] else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
