#!/opt/anaconda3/envs/pyocean/bin/python
r"""Mechanical evaluation of the pre-registered hypotheses H3a-H3c of Exp. 3 (ruling R14).

Rules (exp3_melt_emulator/notes/pd0_preregistro.md, "Hipotesis", approved 2026-11-10);
only the main folds are used, never the drop_n1 sensitivity repeat:

* H3a (final test, 3 runs at Tbot = 2.0): PASS if mean |integrated_rel_error| of M7 < that
  of M6 AND M7 is better (smaller |error|) in >= 2 of the 3 runs; INCONCLUSIVE if the
  difference of the means is < 1.3 % (the noise floor of the integrated melt); else FAIL.
* H3b (loro, 25 main folds): PASS if mean rmse_map of M7 <= 1.10 x that of M6.
* H3c (replication of Ockenden et al. 2026, GMD 19:7979): over the test runs of zc_out (25),
  warm_out (5) and final (3) = 33 cases, separately for M6 and M7. Spread = std (ddof = 0)
  of the 5 members' integrated melt / |integrated melt of the ensemble mean|. PASS if no
  case has spread < 5 % AND |integrated_rel_error| > 20 %. Ockenden's suggested threshold
  (spread > 10 % of the integrated melt; literatura/Ockenden_2026*.md L731, Fig. 8b L571-576)
  is reported as a confusion table against error > 20 %.

The ranking of the 8 methods is reported, not a hypothesis (pd0): per fold kind and method,
mean rmse_map (Burgard et al. 2022 Eq. 32), mean |integrated_rel_error| (B22 Eq. 33-34, one
shelf) and the RMS over runs of the per-run rmse_gl (B22 Eq. 35, ruling R8). The drop_n1
repeat is tabulated the same way and labelled "sensitivity".

Usage::

    hypotheses.py [--cv figures/cv_results.json] [--final figures/final_test.json]
                  [--out figures/hypotheses.json]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
FIG_DIR = REPO / "exp3_melt_emulator" / "figures"

# pd0_preregistro.md, Hipotesis (thresholds approved by Esteban 2026-11-10)
# The two noise floors are read from figures/noise_floor.json (the single source, also
# recorded in protocol_config.json) and rounded to the precision at which the PD0 states
# them (0.157 m a-1, 1.3 %).
_NOISE = json.loads((FIG_DIR / "noise_floor.json").read_text())
H3A_FLOOR = round(_NOISE["median_abs_integrated_rel_diff"], 3)   # 1.3 %: median |N1 drift|
H3A_MIN_BETTER = 2         # M7 better in >= 2 of 3 final runs
H3B_RATIO = 1.10           # loro mean rmse_map M7 <= 1.10 x M6
H3C_SPREAD = 0.05          # violation: spread < 5 % ...
H3C_ERROR = 0.20           # ... and |integrated_rel_error| > 20 %
O26_SPREAD = 0.10          # Ockenden et al. 2026: spread > 10 % flags low confidence
RMSE_FLOOR = round(_NOISE["pooled_median_rms"], 3)   # 0.157 m a-1, rmse_map / rmse_gl
assert (H3A_FLOOR, RMSE_FLOOR) == (0.013, 0.157), "floors differ from the PD0 values"
H3C_KINDS = ("zc_out", "warm_out", "final")
KINDS = ("loro", "zc_out", "warm_out", "profile_test", "final")
METHODS = tuple(f"M{i}" for i in range(8))
EXPECTED = {"final": 3, "loro": 25, "h3c": 33}     # pd0 fold sizes
SPREAD_RTOL = 1e-9


def kind_of(e: dict) -> str:
    """Fold kind of an entry; final_test.json entries carry fold = 'final' and no kind."""
    return e.get("kind") or e["fold"]


def main_entries(cv: dict, final: dict) -> list[dict]:
    """Main-variant entries of cv_results.json plus all final_test.json entries."""
    out = [e for e in cv["entries"] if e.get("variant", "main") == "main"]
    return out + [e for e in final["entries"] if kind_of(e) == "final"]


def _select(entries, kind, method):
    return [e for e in entries if kind_of(e) == kind and e["method"] == method]


def spread_rel(e: dict) -> float:
    """H3c spread, recomputed from the members and checked against the stored value."""
    mem = e.get("members")
    if not mem:
        raise ValueError(f"{e['method']} {e['case']}: no ensemble members")
    s = float(np.std(np.asarray(mem["integrated"], float), ddof=0)
              / abs(mem["integrated_mixture"]))
    stored = mem.get("spread_rel")
    if stored is not None and not np.isclose(s, stored, rtol=SPREAD_RTOL, atol=0):
        raise ValueError(f"{e['method']} {e['case']}: stored spread_rel {stored} != "
                         f"recomputed {s}")
    return s


def h3a(entries, expected: int | None = EXPECTED["final"]) -> dict:
    by = {m: {e["case"]: abs(e["metrics"]["integrated_rel_error"])
              for e in _select(entries, "final", m)} for m in ("M6", "M7")}
    cases = sorted(by["M6"])
    if sorted(by["M7"]) != cases:
        raise ValueError("H3a: M6 and M7 final cases differ")
    if expected is not None and len(cases) != expected:
        raise ValueError(f"H3a: {len(cases)} final runs, expected {expected}")
    mean = {m: float(np.mean([by[m][c] for c in cases])) for m in by}
    better = [c for c in cases if by["M7"][c] < by["M6"][c]]
    diff = mean["M7"] - mean["M6"]
    if abs(diff) < H3A_FLOOR:
        verdict = "INCONCLUSIVE"
    elif mean["M7"] < mean["M6"] and len(better) >= H3A_MIN_BETTER:
        verdict = "PASS"
    else:
        verdict = "FAIL"
    return {"rule": "final test: mean |integrated_rel_error| M7 < M6 AND M7 better in >= 2 "
                    "of 3 runs; INCONCLUSIVE if |difference of means| < 1.3 %",
            "verdict": verdict, "cases": cases,
            "abs_ire": {m: [by[m][c] for c in cases] for m in by},
            "mean_abs_ire": mean, "difference_M7_minus_M6": diff,
            "abs_mean_difference": abs(diff), "floor": H3A_FLOOR,
            "runs_M7_better": better, "n_runs_M7_better": len(better),
            "note": "3 runs: an indication, not a conclusion (pd0)"}


def h3b(entries, expected: int | None = EXPECTED["loro"]) -> dict:
    v = {m: [e["metrics"]["rmse_map"] for e in _select(entries, "loro", m)]
         for m in ("M6", "M7")}
    for m, x in v.items():
        if expected is not None and len(x) != expected:
            raise ValueError(f"H3b: {m} has {len(x)} loro folds, expected {expected}")
    mean = {m: float(np.mean(x)) for m, x in v.items()}
    ratio = mean["M7"] / mean["M6"]
    return {"rule": "loro (main): mean rmse_map M7 <= 1.10 x M6",
            "verdict": "PASS" if mean["M7"] <= H3B_RATIO * mean["M6"] else "FAIL",
            "n_folds": len(v["M7"]), "mean_rmse_map": mean, "ratio_M7_over_M6": ratio,
            "threshold_ratio": H3B_RATIO,
            "difference_M7_minus_M6": mean["M7"] - mean["M6"],
            "difference_below_rmse_floor": abs(mean["M7"] - mean["M6"]) < RMSE_FLOOR,
            "rmse_floor": RMSE_FLOOR}


def h3c_cases(entries, method: str) -> list[dict]:
    out = []
    for kind in H3C_KINDS:
        for e in _select(entries, kind, method):
            s = spread_rel(e)
            err = abs(e["metrics"]["integrated_rel_error"])
            out.append({"case": e["case"], "kind": kind, "spread_rel": s, "abs_ire": err,
                        "n1_pass": e.get("n1_pass"),
                        "violation": bool(s < H3C_SPREAD and err > H3C_ERROR)})
    return out


def ockenden_confusion(cases: list[dict]) -> dict:
    """Spread > 10 % (Ockenden et al. 2026) against |integrated_rel_error| > 20 %."""
    hi_s = [c["spread_rel"] > O26_SPREAD for c in cases]
    hi_e = [c["abs_ire"] > H3C_ERROR for c in cases]
    n = lambda a, b: int(sum(1 for s, e in zip(hi_s, hi_e) if s == a and e == b))
    return {"spread_gt_10pct_error_gt_20pct": n(True, True),
            "spread_gt_10pct_error_le_20pct": n(True, False),
            "spread_le_10pct_error_gt_20pct": n(False, True),
            "spread_le_10pct_error_le_20pct": n(False, False),
            "cases_low_spread_high_error": [c["case"] for c, s, e in zip(cases, hi_s, hi_e)
                                            if not s and e]}


def h3c(entries, expected: int | None = EXPECTED["h3c"]) -> dict:
    res = {"rule": "zc_out (25) + warm_out (5) + final (3) test runs, per method: PASS if "
                   "no case has spread < 5 % AND |integrated_rel_error| > 20 %; spread = "
                   "std (ddof = 0) of the 5 member integrated melts / |ensemble-mean "
                   "integrated melt|",
           "spread_threshold": H3C_SPREAD, "error_threshold": H3C_ERROR,
           "ockenden_spread_threshold": O26_SPREAD}
    for m in ("M6", "M7"):
        cases = h3c_cases(entries, m)
        if expected is not None and len(cases) != expected:
            raise ValueError(f"H3c: {m} has {len(cases)} cases, expected {expected}")
        viol = [c for c in cases if c["violation"]]
        res[m] = {"verdict": "FAIL" if viol else "PASS", "n_cases": len(cases),
                  "n_violations": len(viol), "violations": viol,
                  "ockenden_10pct_confusion": ockenden_confusion(cases), "cases": cases}
    return res


def ranking_table(entries, kinds) -> dict:
    out = {}
    for kind in kinds:
        rows = {}
        for m in METHODS:
            sel = _select(entries, kind, m)
            if not sel:
                continue
            r = np.array([e["metrics"]["rmse_map"] for e in sel])
            a = np.array([abs(e["metrics"]["integrated_rel_error"]) for e in sel])
            g = np.array([e["metrics"]["rmse_gl"] for e in sel])
            rows[m] = {"n_runs": len(sel), "mean_rmse_map": float(r.mean()),
                       "median_rmse_map": float(np.median(r)),
                       "mean_abs_ire": float(a.mean()), "median_abs_ire": float(np.median(a)),
                       "rms_rmse_gl": float(np.sqrt(np.mean(g**2)))}
        if not rows:
            continue
        order = sorted(rows, key=lambda m: rows[m]["mean_rmse_map"])
        best = rows[order[0]]["mean_rmse_map"]
        out[kind] = {
            "methods": rows, "order_by_mean_rmse_map": order,
            "order_by_mean_abs_ire": sorted(rows, key=lambda m: rows[m]["mean_abs_ire"]),
            "order_by_rms_rmse_gl": sorted(rows, key=lambda m: rows[m]["rms_rmse_gl"]),
            "within_rmse_floor_of_best": [m for m in order
                                          if rows[m]["mean_rmse_map"] - best < RMSE_FLOOR]}
    return out


def evaluate(cv: dict, final: dict, expected: dict | None = EXPECTED) -> dict:
    """All verdicts and tables; ``expected`` (pd0 fold sizes) is checked unless None."""
    ex = expected or {}
    main = main_entries(cv, final)
    n1 = [e for e in cv["entries"] if e.get("variant") == "n1"]
    return {
        "H3a": h3a(main, ex.get("final")),
        "H3b": h3b(main, ex.get("loro")),
        "H3c": h3c(main, ex.get("h3c")),
        "ranking": {
            "note": "the ranking of the 8 methods is reported, not a hypothesis (pd0); "
                    "no difference < 0.157 m a-1 (rmse_map, rmse_gl) or < 1.3 % "
                    "(integrated_rel_error) is interpreted",
            "main": ranking_table(main, KINDS),
            "n1": {"label": "sensitivity",
                   "description": "drop_n1 repeat: the 5 N1-failing runs removed from "
                                  "training; never used for H3a-H3c",
                   "kinds": ranking_table(n1, KINDS[:-1])}},
        "amendments": cv.get("amendments", []),
    }


def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(cv_path: Path, final_path: Path, out: Path) -> dict:
    cv = json.loads(Path(cv_path).read_text())
    final = json.loads(Path(final_path).read_text())
    res = evaluate(cv, final)
    rec = {"source": {"cv_results": Path(cv_path).name, "cv_results_sha256": _sha(cv_path),
                      "final_test": Path(final_path).name,
                      "final_test_sha256": _sha(final_path),
                      "config_sha256": cv.get("config_sha256"),
                      "preregistration": "exp3_melt_emulator/notes/pd0_preregistro.md"},
           "verdicts": {"H3a": res["H3a"]["verdict"], "H3b": res["H3b"]["verdict"],
                        "H3c_M6": res["H3c"]["M6"]["verdict"],
                        "H3c_M7": res["H3c"]["M7"]["verdict"]},
           **res}
    Path(out).write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cv", type=Path, default=FIG_DIR / "cv_results.json")
    ap.add_argument("--final", type=Path, default=FIG_DIR / "final_test.json")
    ap.add_argument("--out", type=Path, default=FIG_DIR / "hypotheses.json")
    a = ap.parse_args(argv)
    rec = write(a.cv, a.final, a.out)
    print(json.dumps(rec["verdicts"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
