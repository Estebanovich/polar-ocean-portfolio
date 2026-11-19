"""Figure and validation data for the README (plan Task 10), extracted from the runs into committed JSON:

  figures/spinup.json      spin-up report of G0 (accelerated, V6a) and G0s (synchronous drift)
                           (analysis/spinup.py).
  figures/validation.json  V6b: residual of the annual channel DIC budget / surface term, GC and GT,
                           in each branch year; R and mixed-layer residuals reported without criterion.
  figures/figdata.json     geometry and forcing (Fig. 0), x-mean DIC sections and G200(y) of GT in the
                           decades centred on 1900 and 2060 (Fig. 1), annual channel flux of GC, GT
                           and the branches (Fig. 2).

Calendar: GC and GT start on 1 January 1850 at iter_1850 (their nIter0); the model year that ends at
iter_1850 + k spy is calendar year 1849 + k (year index k, analysis/budget.py).
Usage: python -m exp6_southern_carbon.analysis.figdata [--runs-root DIR] [--only spinup|validation|figdata]"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import budget, collect, gradients, io, spinup
from .hypotheses import BRANCH_YEARS

FIG = Path(__file__).resolve().parents[1] / "figures"
V6B_TOL = 0.01                                   # PD0: residual < 1 % of the surface term


# ---- pure helpers ----

def calendar_years(end_iters, iter_1850: int, steps_per_year: int) -> list[int]:
    return [1849 + (int(i) - int(iter_1850)) // int(steps_per_year) for i in end_iters]


def year_index(year: int) -> int:
    """Model-year index (1-based, analysis/budget.py) of calendar year `year` in GC/GT."""
    return int(year) - 1849


def xmean_section(a3d, hfac) -> np.ndarray:
    """Mean over wet cells along x, (nr, ny); NaN where a row is dry."""
    a, w = np.asarray(a3d, float), np.asarray(hfac, float) > 0
    n = w.sum(axis=2)
    with np.errstate(invalid="ignore"):
        s = np.where(w, a, 0.0).sum(axis=2) / n
    return np.where(n > 0, s, np.nan)


def g200_by_row(dic, z, hfac) -> np.ndarray:
    """G200 (gradients.g200, umol kg-1 m-1) of each row j over its wet surface columns."""
    wet = np.asarray(hfac, float)[0] > 0
    out = np.full(wet.shape[0], np.nan)
    for j in range(wet.shape[0]):
        if wet[j].any():
            m = np.zeros_like(wet)
            m[j] = wet[j]
            out[j] = gradients.g200(dic, z, m)
    return out


def v6b(budgets: dict) -> dict:
    """budgets: {run: {year: dic_budget(...)}} -> |residual| / |surf| per run and year, PD0 verdict."""
    frac = {r: {str(y): abs(b["residual"]) / abs(b["surf"]) for y, b in sorted(ys.items())}
            for r, ys in budgets.items()}
    mx = max(v for ys in frac.values() for v in ys.values())
    return {"frac": frac, "max_frac": float(mx), "tol": V6B_TOL, "pass": bool(mx < V6B_TOL)}


def _nan_to_none(x):
    """JSON-ready: arrays to lists, NaN/inf to null, numpy scalars to Python."""
    if isinstance(x, dict):
        return {str(k): _nan_to_none(v) for k, v in x.items()}
    if isinstance(x, np.ndarray):
        x = x.tolist()
    if isinstance(x, (list, tuple)):
        return [_nan_to_none(v) for v in x]
    if isinstance(x, (np.integer, np.bool_)):
        return x.item()
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    return x


# ---- readers (run directories) ----

def spinup_json(root: Path) -> dict:
    out = {"G0": spinup.report(root / "G0", j_edge=collect.J_EDGE)}
    if collect.finished_ok(root / "G0s"):
        out["G0s"] = spinup.report(root / "G0s", j_edge=collect.J_EDGE)
    return out


def validation_json(root: Path) -> dict:
    col, rml = {}, {}
    for par in ("GC", "GT"):
        run = root / par
        R = collect.region(run)
        col[par] = {y: budget.dic_budget(run, year_index(y)) for y in BRANCH_YEARS}
        rml[par] = {str(y): budget.dic_budget(run, year_index(y), R, "ml") for y in BRANCH_YEARS}
    return {"V6b": v6b(col), "budgets_channel": {p: {str(y): b for y, b in d.items()} for p, d in col.items()},
            "budgets_R_ml": rml}


def _geometry(run: Path) -> dict:
    from ..gendata import fields
    from ..gendata.make_case import grid_of, load_case
    g = grid_of(load_case("GT"))
    f = fields.tutorial_fields(g)
    sst = fields.sst_relax_monthly(g)
    y_km = (np.arange(g.ny) + 0.5) * g.dx_m / 1e3
    sponge = np.where(np.asarray(f["mask"]).max(axis=(0, 2)) > 0)[0]          # rbcs rows (data.rbcs)
    return {"y_km": y_km, "depth_m": -np.asarray(f["bathy"], float).mean(axis=1),
            "tau_x": np.asarray(f["taux"], float).mean(axis=1), "sst_mean": sst.mean(axis=(0, 2)),
            "sst_amp": float((sst.max(axis=0) - sst.min(axis=0)).mean() / 2),
            "sponge_y_km": [float(y_km[sponge[0]] - g.dx_m / 2e3), float(y_km[-1] + g.dx_m / 2e3)] if sponge.size else None,
            "R_edge_km": float(collect.J_EDGE * g.dx_m / 1e3)}


def _sections(run: Path) -> dict:
    g = io.grid(run)
    n0, spy = int(io.case(run)["nIter0"]), int(io.case(run)["steps_per_year"])
    R = collect.region(run)
    out = {"z_m": g["rc"], "y_km": g["yc"] / 1e3, "dic": {}, "g200_y": {}, "g200_R": {}}
    for y in (1900, 2060):
        nb = n0 + (y - 1850) * spy
        dec = [i for i in io.iters(run, "carbon_ann") if nb - 5 * spy < i <= nb + 5 * spy]
        dic = np.mean([io.fields(run, "carbon_ann", i)["TRAC01"] for i in dec], axis=0)
        out["dic"][str(y)] = gradients._umolkg(xmean_section(dic, g["hfac"]))
        out["g200_y"][str(y)] = g200_by_row(dic, g["rc"], g["hfac"])
        out["g200_R"][str(y)] = gradients.g200(dic, g["rc"], R)
    return out


def _series(root: Path) -> dict:
    i1850 = int(io.case(root / "GT")["nIter0"])
    out = {}
    for c in ["GC", "GT"] + [f"{p}_w{y}" for y in BRANCH_YEARS for p in ("GC", "GT")]:
        run = root / c
        if not collect.finished_ok(run):
            continue
        f = collect.flux_series(run)
        spy = int(io.case(run)["steps_per_year"])
        out[c] = {"years": calendar_years(list(f), i1850, spy), "flux": list(f.values())}
    return out


def figdata_json(root: Path) -> dict:
    return {"geometry": _geometry(root / "GT"), "sections": _sections(root / "GT"), "series": _series(root)}


def main(argv=None) -> int:
    import argparse

    from pfcommon.paths import EXP6_RUNS
    ap = argparse.ArgumentParser(description="Write figures/{spinup,validation,figdata}.json")
    ap.add_argument("--runs-root", type=Path, default=EXP6_RUNS)
    ap.add_argument("--only", choices=("spinup", "validation", "figdata"))
    a = ap.parse_args(argv)
    jobs = {"spinup": spinup_json, "validation": validation_json, "figdata": figdata_json}
    FIG.mkdir(parents=True, exist_ok=True)
    for name, fn in jobs.items():
        if a.only and name != a.only:
            continue
        out = FIG / f"{name}.json"
        out.write_text(json.dumps(_nan_to_none(fn(a.runs_root)), indent=1) + "\n")
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
