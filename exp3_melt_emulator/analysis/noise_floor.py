#!/opt/anaconda3/envs/pyocean/bin/python
r"""Internal-variability noise floor of the Exp. 3 melt maps (spec Sec. 3, ruling R4).

There are no duplicate runs, so the noise of a run's melt map is estimated by the
difference between its mean map over months 19-24 (records 18-23, the dataset target) and
over months 13-18 (records 12-17), the same windows as Exp. 2 pd0 N1
(``metrics.calendar_window``: 'last6', 'prev6'). Per run:

    RMS = sqrt( sum_shelf a (m_19-24 - m_13-18)^2 / sum_shelf a ),  a = RAC [m^2]

over the run's own shelf mask (pd0: wet columns with ice base < 0). The noise floor is the
median of the per-run RMS over the 25 training runs (group != reference, split != test):
no difference between methods smaller than it is interpreted (PD0 of Exp. 3). The
relative difference of the shelf-mean melt, (<m>_19-24 - <m>_13-18) / <m>_19-24 (the
signed pd0 N1 drift), and its median |.| are recorded alongside, as the analogous floor
for the integrated-melt error.

The maps are read directly from the Exp. 2 runs on LaCie (read-only use of
``exp2_isomip_cavity.analysis.metrics``; the Exp. 2 exporter and dataset are untouched).
Each run is checked against the Exp. 3 dataset: same shelf mask and the same months 19-24
map (to 1e-9 m/a), otherwise the script refuses.

Usage (pyocean; no torch / scikit-learn)::

    noise_floor.py [--runs-root /Volumes/LaCie/portfolio_runs/exp2] [--run-suffix _r21]
                   [--out /Volumes/LaCie/portfolio_runs/exp3/noise_floor.npz]
                   [--summary exp3_melt_emulator/figures/noise_floor.json]
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

from exp2_isomip_cavity.analysis import export_exp3 as ex  # noqa: E402
from exp2_isomip_cavity.analysis import metrics as mt  # noqa: E402
from exp3_melt_emulator.data import DATASET_PATH, DATASET_SHA256, load_runs  # noqa: E402

RUNS_ROOT = Path("/Volumes/LaCie/portfolio_runs/exp2")
RUN_SUFFIX = "_r21"                    # exp2_isomip_cavity/notes/exp3_dataset.md
DEFAULT_OUT = Path("/Volumes/LaCie/portfolio_runs/exp3/noise_floor.npz")
DEFAULT_SUMMARY = REPO / "exp3_melt_emulator" / "figures" / "noise_floor.json"
MATCH_ATOL = 1e-9                      # m/a, months 19-24 map vs the dataset melt


def _window_maps(rd: mt.RunData) -> tuple[np.ndarray, np.ndarray]:
    out = []
    for name in ("prev6", "last6"):     # pd0: months 13-18, months 19-24
        win = list(mt.calendar_window(rd, name))
        out.append(mt.melt_rate_m_per_a(rd.mean_field("shiDiag", "SHIfwFlx", win)))
    return out[0], out[1]


def noise_maps(run_dir) -> tuple[np.ndarray, np.ndarray]:
    """(mean melt map of months 13-18, of months 19-24), m a-1, full grid."""
    return _window_maps(mt.load_run(run_dir))


def area_rms(diff, mask, area) -> float:
    """Area-weighted RMS of ``diff`` over ``mask``."""
    mask = np.asarray(mask, bool)
    if not mask.any():
        raise ValueError("empty mask")
    a, d = np.asarray(area, float)[mask], np.asarray(diff, float)[mask]
    return float(np.sqrt(np.sum(a * d * d) / np.sum(a)))


def training_cases(runs) -> list[str]:
    """The training pool of folds.make_folds: group != reference and split != test."""
    return [r.case for r in runs if r.group != "reference" and r.split != "test"]


def build(runs, runs_root=RUNS_ROOT, suffix: str = RUN_SUFFIX) -> dict:
    """Noise maps and per-run RMS of the training runs of ``runs`` (dataset Runs)."""
    keep = set(training_cases(runs))
    prev_l, last_l, mask_l, per_run = [], [], [], []
    for r in (r for r in runs if r.case in keep):
        rd = mt.load_run(Path(runs_root) / f"{r.case}{suffix}")
        prev, last = _window_maps(rd)
        if not np.array_equal(rd.shelf, r.shelf_mask):
            raise ValueError(f"{r.case}: shelf mask differs from the Exp. 3 dataset")
        dev = float(np.max(np.abs(last[r.shelf_mask] - r.melt[r.shelf_mask])))
        if not dev <= MATCH_ATOL:
            raise ValueError(f"{r.case}: months 19-24 map differs from the dataset melt "
                             f"by {dev:.3g} m/a")
        rms = area_rms(last - prev, r.shelf_mask, rd.rA)
        a = rd.rA[r.shelf_mask]
        mean_last = float(np.sum(last[r.shelf_mask] * a) / np.sum(a))
        mean_prev = float(np.sum(prev[r.shelf_mask] * a) / np.sum(a))
        per_run.append({"case": r.case, "group": r.group, "Tbot": r.Tbot, "z_c": r.z_c,
                        "s": r.s, "N1_pass": r.n1_pass, "rms": rms,
                        "mean_melt_19_24": mean_last, "rms_over_mean": rms / mean_last,
                        "integrated_rel_diff": (mean_last - mean_prev) / mean_last})
        prev_l.append(np.where(r.shelf_mask, prev, np.nan))
        last_l.append(np.where(r.shelf_mask, last, np.nan))
        mask_l.append(r.shelf_mask)
    if not per_run:
        raise ValueError("no training run")
    rms = np.array([p["rms"] for p in per_run])
    summary = {
        "definition": "per run: area-weighted RMS over the run's shelf_mask of "
                      "(mean melt months 19-24 - mean melt months 13-18), m a-1; "
                      "noise floor = median over the training runs",
        "windows": {"prev6": "months 13-18 (records 12-17)",
                    "last6": "months 19-24 (records 18-23)"},
        "runs_root": str(runs_root), "run_suffix": suffix,
        "n_runs": len(per_run), "pooled_median_rms": float(np.median(rms)),
        "min_rms": float(rms.min()), "max_rms": float(rms.max()),
        "median_abs_integrated_rel_diff": float(np.median(
            [abs(p["integrated_rel_diff"]) for p in per_run])),
        "per_run": per_run,
    }
    arrays = {"case": np.array([p["case"] for p in per_run], dtype=str),
              "prev6": np.stack(prev_l), "last6": np.stack(last_l),
              "shelf_mask": np.stack(mask_l), "rms": rms}
    return {"summary": summary, "arrays": arrays}


def write_outputs(res: dict, out, summary_path, dataset_sha: str | None = None) -> str:
    """Byte-reproducible .npz (Exp. 2 writer) and the JSON summary; returns the sha256."""
    sha = ex.write_npz(res["arrays"], out)
    rec = {"file": str(out), "sha256": sha, "dataset_sha256": dataset_sha, **res["summary"]}
    Path(summary_path).parent.mkdir(parents=True, exist_ok=True)
    Path(summary_path).write_text(json.dumps(rec, indent=1) + "\n")
    return sha


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs-root", type=Path, default=RUNS_ROOT)
    ap.add_argument("--run-suffix", default=RUN_SUFFIX)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    a = ap.parse_args(argv)
    runs = load_runs(DATASET_PATH, DATASET_SHA256)          # sha checked
    res = build(runs, a.runs_root, a.run_suffix)
    sha = write_outputs(res, a.out, a.summary, DATASET_SHA256)
    s = res["summary"]
    print(f"{s['n_runs']} runs; noise RMS median {s['pooled_median_rms']:.4f} m/a "
          f"(min {s['min_rms']:.4f}, max {s['max_rms']:.4f}); {a.out} sha256 {sha}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
