#!/opt/anaconda3/envs/pyocean/bin/python
r"""Bootstrap intervals over runs for the Exp. 3 ranking (ruling R17).

pd0_preregistro.md ("Métricas") asks for every metric "por corrida, ... y con intervalos
bootstrap sobre las corridas". It does not fix the number of resamples or the seed; they
are the analyst's choice, made at analysis time (after the results were known): B = 10000,
seed 0. Nothing is refitted: the per-run metrics already in cv_results.json (both variants,
every fold kind) and final_test.json (main variant only) are resampled.

For each variant (main | n1) and fold kind with at least 2 held-out runs, the runs are
resampled with replacement B times; the same resampled run indices are used for every
method (paired), so the intervals of two methods are comparable. For each method, the
statistics are the mean over the resampled runs of rmse_map (Burgard et al. 2022 Eq. 32),
of |integrated_rel_error| (B22 Eq. 33-34, one shelf) and of the per-run rmse_gl (B22
Eq. 35, ruling R8), plus the RMS over runs of rmse_gl, the aggregate reported in the
ranking. The interval is the percentile 95 % interval (2.5 and 97.5 percentiles, numpy
linear interpolation). With 3 final-test runs those intervals are indicative only. The
hypotheses H3a-H3c are defined on the means and are not affected.

Usage::

    bootstrap.py [--cv figures/cv_results.json] [--final figures/final_test.json]
                 [--out figures/bootstrap.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
FIG_DIR = REPO / "exp3_melt_emulator" / "figures"
N_BOOT = 10000             # analyst's choice (R17), not fixed by the pre-registration
SEED = 0                   # analyst's choice (R17)
LEVEL = 0.95
STATS = ("rmse_map", "abs_integrated_rel_error", "rmse_gl", "rmse_gl_rms")
KIND_ORDER = ("loro", "zc_out", "warm_out", "profile_test", "final")


def _per_run(entries: list[dict]) -> dict[str, np.ndarray]:
    m = [e["metrics"] for e in entries]
    return {"rmse_map": np.array([x["rmse_map"] for x in m]),
            "abs_integrated_rel_error": np.abs([x["integrated_rel_error"] for x in m]),
            "rmse_gl": np.array([x["rmse_gl"] for x in m])}


def _stats(v: dict[str, np.ndarray], idx: np.ndarray) -> dict[str, np.ndarray]:
    """Statistics over the runs; idx (B, n) resampled indices, or (1, n) the original runs."""
    return {"rmse_map": v["rmse_map"][idx].mean(1),
            "abs_integrated_rel_error": v["abs_integrated_rel_error"][idx].mean(1),
            "rmse_gl": v["rmse_gl"][idx].mean(1),
            "rmse_gl_rms": np.sqrt((v["rmse_gl"][idx] ** 2).mean(1))}


def _groups(cv: dict, final: dict) -> dict[tuple[str, str], list[dict]]:
    g: dict[tuple[str, str], list[dict]] = {}
    for e in cv["entries"]:
        g.setdefault((e.get("variant", "main"), e.get("kind") or e["fold"]), []).append(e)
    for e in final["entries"]:
        g.setdefault(("main", "final"), []).append(e)
    return g


def bootstrap(cv: dict, final: dict, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    rng = np.random.default_rng(seed)
    lo_q, hi_q = 100 * (1 - LEVEL) / 2, 100 * (1 + LEVEL) / 2
    out: dict[str, dict] = {}
    groups = _groups(cv, final)
    keys = sorted(groups, key=lambda k: (k[0] != "main", k[0], KIND_ORDER.index(k[1])))
    for variant, kind in keys:
        ents = groups[(variant, kind)]
        cases = sorted({e["case"] for e in ents})
        if len(cases) < 2:
            continue
        idx = rng.integers(0, len(cases), size=(n_boot, len(cases)))   # paired over methods
        orig = np.arange(len(cases))[None, :]
        res = {}
        for method in sorted({e["method"] for e in ents}):
            by_case = {e["case"]: e for e in ents if e["method"] == method}
            if sorted(by_case) != cases:
                raise ValueError(f"{variant}/{kind}/{method}: runs differ from the other methods")
            v = _per_run([by_case[c] for c in cases])
            point, boot = _stats(v, orig), _stats(v, idx)
            res[method] = {"n_runs": len(cases), **{
                s: {"mean": float(point[s][0]),
                    "ci95": [float(np.percentile(boot[s], lo_q)), float(np.percentile(boot[s], hi_q))]}
                for s in STATS}}
        out.setdefault(variant, {})[kind] = res
    return {"variants": out}


def write(cv_path, final_path, out, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    cv = json.loads(Path(cv_path).read_text())
    final = json.loads(Path(final_path).read_text())
    rec = {"method": "percentile bootstrap over held-out runs, paired across methods; "
                     "per-run metrics resampled, no refits",
           "n_boot": n_boot, "seed": seed, "level": LEVEL,
           "choice": "B and seed are the analyst's choice made at analysis time (R17); "
                     "the pre-registration asks for bootstrap intervals over runs but does not "
                     "fix them",
           "note": "final test: 3 runs, intervals indicative only; H3a-H3c are defined on "
                   "means and do not use these intervals; kinds with < 2 runs are skipped",
           "sources": [Path(cv_path).name, Path(final_path).name],
           **bootstrap(cv, final, n_boot, seed)}
    Path(out).write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cv", type=Path, default=FIG_DIR / "cv_results.json")
    ap.add_argument("--final", type=Path, default=FIG_DIR / "final_test.json")
    ap.add_argument("--out", type=Path, default=FIG_DIR / "bootstrap.json")
    a = ap.parse_args(argv)
    rec = write(a.cv, a.final, a.out)
    print(f"wrote {a.out}: variants {sorted(rec['variants'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
