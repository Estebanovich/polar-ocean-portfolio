#!/opt/anaconda3/envs/pyocean/bin/python
r"""Summary of the pre-registered probabilistic metrics and skill of Exp. 3.

pd0_preregistro.md ("Métricas"): skill = 1 - err / err_M0 (same fold and variant), and for
M5, M6 and M7 the empirical coverage of the central 68 % and 95 % predictive intervals and
the CRPS. run_protocol.py computed them per held-out run (cv_results.json, final_test.json,
keys "skill" and "prob"); nothing is recomputed or refitted here. This script averages them
over the held-out runs of each variant and fold kind, per method (M0, the skill reference,
is left out), and writes figures/probabilistic.json for the README.

Usage::

    probabilistic.py [--cv figures/cv_results.json] [--final figures/final_test.json]
                     [--out figures/probabilistic.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
FIG_DIR = REPO / "exp3_melt_emulator" / "figures"
PROB_METHODS = ("M5", "M6", "M7")
PROB_KEYS = ("coverage_68", "coverage_95", "crps", "mean_sigma")
SKILL_KEYS = ("rmse_map", "abs_integrated_rel_error", "rmse_gl")
KIND_ORDER = ("loro", "zc_out", "warm_out", "profile_test", "final")


def summarise(cv: dict, final: dict) -> dict:
    groups: dict[tuple, list] = {}
    for e in cv["entries"]:
        groups.setdefault((e.get("variant", "main"), e.get("kind") or e["fold"], e["method"]), []).append(e)
    for e in final["entries"]:
        groups.setdefault(("main", "final", e["method"]), []).append(e)
    out: dict = {}
    for (variant, kind, method) in sorted(groups, key=lambda k: (k[0] != "main", k[0],
                                                                 KIND_ORDER.index(k[1]), k[2])):
        if method == "M0":
            continue
        ents = groups[(variant, kind, method)]
        rec = {"n_runs": len(ents)}
        for k in SKILL_KEYS:
            rec[f"skill_{k}"] = float(np.mean([e["skill"][k] for e in ents]))
        if method in PROB_METHODS:
            for k in PROB_KEYS:
                rec[k] = float(np.mean([e["prob"][k] for e in ents]))
        out.setdefault(variant, {}).setdefault(kind, {})[method] = rec
    return {"variants": out}


def write(cv_path, final_path, out) -> dict:
    cv = json.loads(Path(cv_path).read_text())
    final = json.loads(Path(final_path).read_text())
    rec = {"definition": "means over the held-out runs of the per-run values stored by "
                         "run_protocol.py: skill = 1 - err/err_M0 of the same fold and variant "
                         "(rmse_map, |integrated_rel_error|, rmse_gl); for M5-M7 the empirical "
                         "coverage of the central 68 % and 95 % predictive intervals over the "
                         "run's shelf cells, the CRPS (m a-1) and the mean predictive sigma",
           "sources": [Path(cv_path).name, Path(final_path).name], **summarise(cv, final)}
    Path(out).write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cv", type=Path, default=FIG_DIR / "cv_results.json")
    ap.add_argument("--final", type=Path, default=FIG_DIR / "final_test.json")
    ap.add_argument("--out", type=Path, default=FIG_DIR / "probabilistic.json")
    a = ap.parse_args(argv)
    write(a.cv, a.final, a.out)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
