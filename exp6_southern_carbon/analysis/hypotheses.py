"""H6a and H6b (spec Sec. 7; frozen in notes/pd0_preregistro.md).

H6a (Esteban: signed, sign and detection only): dS = S_GT(2060) - S_GT(1900) > 0 (the wind
response moves toward uptake) and dS > 2 |S_GC(2060) - S_GC(1900)| (control drift, V6c).
S_nat = S_GC, S_ant = S_GT - S_GC (Lovenduski et al. 2008).
H6b: |pred - meas| <= 0.25 |meas| in every GT branch."""
from __future__ import annotations

import numpy as np


def h6a(s_gt: dict[int, float], s_gc: dict[int, float]) -> dict:
    ds = s_gt[2060] - s_gt[1900]
    drift = abs(s_gc[2060] - s_gc[1900])
    years = sorted(s_gt)
    s0 = np.sign(s_gt[1900])
    cross = next((y for y in years if y > 1900 and np.sign(s_gt[y]) != s0), None)
    return {"dS": float(ds), "drift_gc": float(drift), "pass": bool(ds > 0 and ds > 2 * drift),
            "s_nat": {y: float(s_gc[y]) for y in sorted(s_gc)},
            "s_ant": {y: float(s_gt[y] - s_gc[y]) for y in years if y in s_gc},
            "zero_crossing_year": cross}


def h6b(pred: dict[int, float], meas: dict[int, float], tol: float = 0.25) -> dict:
    within = {y: bool(abs(pred[y] - meas[y]) <= tol * abs(meas[y])) for y in sorted(meas)}
    return {"within": within, "n_within": int(sum(within.values())), "pass": bool(all(within.values())),
            "ratio": {y: float(pred[y] / meas[y]) if meas[y] else float("nan") for y in sorted(meas)}}


BRANCH_YEARS = (1900, 1980, 2020, 2060)


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return x.item()
    return x


def assemble(s_gt: dict, s_gc: dict, pred: dict, meas: dict, terms: dict) -> dict:
    """H6a and H6b verdicts with the numbers the PD0 asks to report, JSON-ready (string keys)."""
    a = h6a(s_gt, s_gc)
    a["dS_pct_of_S1900"] = float(100.0 * a["dS"] / abs(s_gt[1900])) if s_gt[1900] else float("nan")
    return _jsonable({"h6a": a, "h6b": h6b(pred, meas), "s_gt": s_gt, "s_gc": s_gc,
                      "pred": pred, "meas": meas, "terms": terms})


def evaluate(runs_root) -> dict:
    """Reads GC, GT and their branches under runs_root (plan Task 8)."""
    from pathlib import Path

    from . import collect, decomp, io
    root = Path(runs_root)
    f = {c: collect.flux_series(root / c) for c in ("GC", "GT")}
    s = {"GC": {}, "GT": {}}
    pred, meas, terms = {}, {}, {}
    for y in BRANCH_YEARS:
        for par in ("GC", "GT"):
            br = root / f"{par}_w{y}"
            c = io.case(br)
            n0, spy = int(c["nIter0"]), int(c["steps_per_year"])
            s[par][y] = collect.branch_sensitivity(collect.flux_series(br), f[par], n0, spy)
            if par == "GT":
                p = dict(collect.parent_state(root / par, n0), **collect.circulation(root / par, n0))
                b = collect.circulation(br, n0)
                terms[y] = dict(decomp.terms(p, b), parent=p, branch=b)
                terms[y].update(collect.temperature(br, root / par, n0))   # reported, not in sum_circ_bio
                pred[y] = terms[y]["sum_circ_bio"]
                terms[y]["meas_parts"] = collect.measured_parts(br, root / par, n0)   # post hoc (22 nov), reported
                meas[y] = float(sum(terms[y]["meas_parts"].values()))
    return assemble(s["GT"], s["GC"], pred, meas, terms)


def main(argv=None) -> int:
    import argparse
    import json
    from pathlib import Path

    from pfcommon.paths import EXP6_RUNS
    ap = argparse.ArgumentParser(description="Evaluate H6a and H6b; writes figures/hypotheses.json")
    ap.add_argument("--runs-root", type=Path, default=EXP6_RUNS)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "figures" / "hypotheses.json")
    a = ap.parse_args(argv)
    out = evaluate(a.runs_root)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(f"H6a {'PASS' if out['h6a']['pass'] else 'FAIL'}  H6b {'PASS' if out['h6b']['pass'] else 'FAIL'}  -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
