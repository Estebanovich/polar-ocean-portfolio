#!/opt/anaconda3/envs/pyocean/bin/python
r"""Exp. 2 ensemble table and the pre-registered H1 / H2 verdicts (notes/pd0_preregistro.md).

Usage::

    ensemble.py [--runs-root /Volumes/LaCie/portfolio_runs/exp2] [--run-suffix _r21]
                [--out figures/ensemble_table.json]

``collect(runs_root)`` returns a pandas-free table (dict of equal-length lists), one row per
ensemble member of cases.yaml (the reference ``ocean0_4km`` and the 28 ``ens_*`` cases):

    case, group, Tbot, z_c, s, split, status, N3_pass, N1, N1_pass, N2, N2_pass,
    mean_melt_shelf, mean_melt_deep300, total_flux, z_out, D_out, T_f_GL, notes

* every quantity is the pd0 one over the last 6 months (records 19-24, ``metrics.LAST6``),
  computed with ``exp2_isomip_cavity/analysis/metrics.py`` (``load_run`` enforces
  ``require_ok`` and the R21 ``melt_law`` guard before anything is read);
* a run that fails N3 (``load_run`` refuses it), that does not exist, or whose case.json
  disagrees with the registry is a **gap**: it is listed with its status and all its metric
  fields are None. Nothing is imputed or interpolated (Review Focus 5);
* split: 'test' for Tbot = 2.0 (the extrapolation group, Exp. 3 test data only), 'train'
  otherwise. Test runs never enter H1 / H2 (pd0: Tbot in {-0.5, ..., 1.5}).

H1 (pd0 "H1"), per z_c in {250, 400, 550} over the 5 core runs (group factorial, s = 1):
  1. quadratic least squares m = a Tbot^2 + b Tbot + c, pass if R^2 >= 0.98;
  2. power law m = k (Tbot - T_f,GL)^n fitted in log-log (ordinary least squares of ln m on
     ln(Tbot - T_f,GL)), pass if 1.3 <= n <= 2.0; T_f,GL = metrics.freezing_point_gl
     (AD16 Eq. 25, Table 4; Sbot = 34.7, pressure of the deepest wet ice base);
  H1 passes at a z_c iff 1 and 2 pass. A z_c with any gap is NOT_EVALUABLE.
  pd0 writes "<m>" without naming the mask, so ``h1_verdicts`` has no default for it: the
  caller must pass ``melt_key`` ('mean_melt_shelf' or 'mean_melt_deep300') explicitly.

H2 (pd0 "H2"), per z_c over the same 5 core runs: D_out = -z_out (``metrics.outflow_depth``,
tracer-transport-weighted depth of u C > 0 at x_front); pass if Spearman
rho(Tbot, D_out) <= -0.9. A z_c with any gap is NOT_EVALUABLE.
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

from exp2_isomip_cavity.analysis import metrics as mt  # noqa: E402

CORE_TBOT = (-0.5, 0.0, 0.5, 1.0, 1.5)       # pd0 H1/H2: Tbot in {-0.5, ..., 1.5} degC
CORE_ZC = (250.0, 400.0, 550.0)              # pd0 H1/H2: z_c in {250, 400, 550} m
TEST_TBOT = 2.0                              # extrapolation group: Exp. 3 test only
OCEAN0_TBOT = 1.0                            # AD16 Table 6: Ocean0 bottom temperature
H1_R2_MIN = 0.98                             # pd0 H1.1
H1_N_RANGE = (1.3, 2.0)                      # pd0 H1.2
H2_RHO_MAX = -0.9                            # pd0 H2
N1_MAX, N2_MAX = 0.03, 0.01                  # pd0 N1, N2
MELT_KEYS = ("mean_melt_shelf", "mean_melt_deep300")
RHO_DECIMALS = 12        # Spearman rho is rounded to 1e-12 before the pd0 comparison: with
                         # n = 5 and no ties rho is a multiple of 0.1 (e.g. exactly -0.9),
                         # which floating point can land a few ulp on either side of
DEFAULT_RUNS = Path("/Volumes/LaCie/portfolio_runs/exp2")
COLUMNS = ("case", "group", "Tbot", "z_c", "s", "split", "status", "N3_pass", "N1",
           "N1_pass", "N2", "N2_pass", "mean_melt_shelf", "mean_melt_deep300", "total_flux",
           "z_out", "D_out", "T_f_GL", "notes")
METRIC_COLUMNS = ("N1", "N1_pass", "N2", "N2_pass", "mean_melt_shelf", "mean_melt_deep300",
                  "total_flux", "z_out", "D_out", "T_f_GL")


# ---- registry -------------------------------------------------------------------------------

def registry_cases(cases_file: Path | None = None) -> list[dict]:
    """The ensemble members of cases.yaml in registry order (make_case.ensemble_cases):
    {case, group, Tbot, z_c, s}. The reference ocean0_4km has the Ocean0 linear profile
    (Tbot = 1.0, AD16 Table 6; no z_c). Geometry-group (Tbot, z_c) are the Latin-hypercube
    values make_case generates from the seed (cross-checked against case.json in collect)."""
    from exp2_isomip_cavity.gendata.make_case import ensemble_cases, load_cases
    spec = load_cases(cases_file)
    out = []
    for name in ensemble_cases(spec):
        c = spec["cases"][name]
        group = c.get("ensemble_group")
        if group == "reference":
            out.append({"case": name, "group": group, "Tbot": OCEAN0_TBOT, "z_c": None,
                        "s": float(c.get("s", 1.0))})
        else:
            out.append({"case": name, "group": group, "Tbot": float(c["Tbot"]),
                        "z_c": float(c["z_c"]), "s": float(c["s"])})
    return out


def split_of(Tbot) -> str:
    return "test" if Tbot is not None and abs(float(Tbot) - TEST_TBOT) < 1e-9 else "train"


# ---- collect --------------------------------------------------------------------------------

def _mismatch(spec: dict, actual: dict) -> list[str]:
    bad = []
    for k in ("Tbot", "z_c", "s"):
        if spec.get(k) is None or spec["group"] == "reference" and k != "s":
            continue
        a = actual.get(k)
        if a is None or abs(float(a) - float(spec[k])) > 1e-9:
            bad.append(f"{k}: registry {spec[k]} vs case.json {a}")
    return bad


def run_row(run_dir: Path, spec: dict) -> dict:
    """One table row. Any failure of load_run is an N3 gap (nothing else is read)."""
    row = {c: None for c in COLUMNS}
    row.update(case=spec["case"], group=spec["group"], Tbot=spec["Tbot"], z_c=spec["z_c"],
               s=spec["s"], split=split_of(spec["Tbot"]), notes=[])
    if not run_dir.exists():
        row.update(status="GAP_MISSING", N3_pass=False)
        return row
    try:
        rd = mt.load_run(run_dir)
    except Exception as e:                       # noqa: BLE001 - any refusal is an N3 gap
        row.update(status=f"GAP_N3: {type(e).__name__}: {e}"[:400], N3_pass=False)
        return row
    bad = _mismatch(spec, rd.actual)
    if bad:
        row.update(status="GAP_MISMATCH: " + "; ".join(bad), N3_pass=True)
        return row
    row.update(status="OK", N3_pass=True)
    win = list(mt.calendar_window(rd, "last6"))          # pd0: records 19-24

    def put(key, fn):
        try:
            row[key] = float(fn())
        except (ValueError, RuntimeError) as e:          # reported, not imputed
            row["notes"].append(f"{key}: {type(e).__name__}: {e}")

    put("mean_melt_shelf", lambda: mt.mean_melt(rd, win, "shelf"))
    put("mean_melt_deep300", lambda: mt.mean_melt(rd, win, "deep300"))
    put("total_flux", lambda: mt.total_melt_flux(rd, win))
    put("z_out", lambda: mt.outflow_depth(rd, win))
    if row["z_out"] is not None:
        row["D_out"] = -row["z_out"]                     # pd0 H2: D_out = -z_out
    put("N1", lambda: mt.equilibrium_drift(rd, "2yr"))
    put("N2", lambda: mt.tracer_budget_residual(rd))
    put("T_f_GL", lambda: mt.freezing_point_gl(rd)["T_f_GL_degC"])
    row["N1_pass"] = None if row["N1"] is None else bool(row["N1"] < N1_MAX)
    row["N2_pass"] = None if row["N2"] is None else bool(abs(row["N2"]) < N2_MAX)
    return row


def collect(runs_root, suffix: str = "_r21", cases: list[dict] | None = None) -> dict:
    """Ensemble table (dict of lists, one entry per member, registry order); see module doc."""
    cases = registry_cases() if cases is None else cases
    rows = [run_row(Path(runs_root) / f"{c['case']}{suffix}", c) for c in cases]
    return {k: [r[k] for r in rows] for k in COLUMNS}


def rows(table: dict) -> list[dict]:
    n = len(table.get("case", []))
    return [{k: table[k][i] for k in table} for i in range(n)]


def write_table(table: dict, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(table, indent=1) + "\n")
    return path


def read_table(path: Path) -> dict:
    return json.loads(Path(path).read_text())


# ---- fits -----------------------------------------------------------------------------------

def quadratic_fit(T, m) -> dict:
    """pd0 H1.1: least squares m = a T^2 + b T + c; R^2 = 1 - SS_res / SS_tot."""
    T, m = np.asarray(T, float), np.asarray(m, float)
    A = np.vstack([T**2, T, np.ones_like(T)]).T
    coef, *_ = np.linalg.lstsq(A, m, rcond=None)
    res = m - A @ coef
    ss_tot = float(((m - m.mean()) ** 2).sum())
    r2 = 1.0 - float((res**2).sum()) / ss_tot if ss_tot > 0 else float("nan")
    return {"a": float(coef[0]), "b": float(coef[1]), "c": float(coef[2]), "R2": r2}


def power_law_fit(T, m, T_f: float) -> dict:
    """pd0 H1.2: m = k (T - T_f)^n by ordinary least squares of ln m on ln(T - T_f)."""
    T, m = np.asarray(T, float), np.asarray(m, float)
    dT = T - T_f
    if not (m > 0).all():
        return {"n": None, "k": None, "R2_loglog": None,
                "reason": "melt not positive for every run: ln m undefined"}
    if not (dT > 0).all():
        return {"n": None, "k": None, "R2_loglog": None,
                "reason": "Tbot - T_f,GL not positive for every run"}
    x, y = np.log(dT), np.log(m)
    n, lnk = np.polyfit(x, y, 1)
    res = y - (n * x + lnk)
    r2 = 1.0 - float((res**2).sum()) / float(((y - y.mean()) ** 2).sum())
    return {"n": float(n), "k": float(np.exp(lnk)), "R2_loglog": r2}


def _ranks(v) -> np.ndarray:
    v = np.asarray(v, float)
    order = np.argsort(v, kind="mergesort")
    r = np.empty(len(v))
    r[order] = np.arange(1, len(v) + 1)
    for val in np.unique(v):                     # ties: average rank
        sel = v == val
        r[sel] = r[sel].mean()
    return r


def spearman(x, y) -> float:
    """Spearman rank correlation: Pearson correlation of the (average-tie) ranks.
    None if either variable is constant (rho undefined)."""
    rx, ry = _ranks(x), _ranks(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    den = np.sqrt((rx**2).sum() * (ry**2).sum())
    return None if den == 0 else float((rx * ry).sum() / den)


# ---- verdicts -------------------------------------------------------------------------------

def core_rows(table: dict, z_c: float) -> list[dict]:
    """The pd0 core runs of one z_c: group factorial, s = 1, Tbot in CORE_TBOT (sorted)."""
    out = [r for r in rows(table)
           if r["group"] == "factorial" and r["split"] == "train" and r["s"] is not None
           and abs(r["s"] - 1.0) < 1e-9 and r["z_c"] is not None
           and abs(r["z_c"] - z_c) < 1e-9
           and any(abs(r["Tbot"] - t) < 1e-9 for t in CORE_TBOT)]
    return sorted(out, key=lambda r: r["Tbot"])


def _core(table: dict, z_c: float, keys: tuple) -> tuple[list[dict], list[str]]:
    rs = core_rows(table, z_c)
    have = {round(r["Tbot"], 6) for r in rs}
    gaps = [r["case"] for r in rs
            if r["status"] != "OK" or any(r[k] is None for k in keys)]
    missing = [t for t in CORE_TBOT if round(t, 6) not in have]
    gaps += [f"(no registry row for Tbot = {t}, z_c = {z_c:.0f})" for t in missing]
    return rs, gaps


def _quality(rs: list[dict]) -> dict:
    return {"N1_fail": [r["case"] for r in rs if r["N1_pass"] is False],
            "N2_fail": [r["case"] for r in rs if r["N2_pass"] is False]}


def h1_verdicts(table: dict, *, melt_key: str) -> dict:
    """pd0 H1, per z_c (see module doc). ``melt_key`` is required (pd0 does not name it)."""
    if melt_key not in MELT_KEYS:
        raise ValueError(f"melt_key must be one of {MELT_KEYS}, not {melt_key!r}")
    per, tfs = {}, []
    for z in CORE_ZC:
        rs, gaps = _core(table, z, (melt_key, "T_f_GL"))
        key = f"{z:.0f}"
        if gaps:
            per[key] = {"verdict": "NOT_EVALUABLE", "gaps": gaps, "n_runs": len(rs) - len(
                [g for g in gaps if not g.startswith("(")])}
            continue
        T = [r["Tbot"] for r in rs]
        m = [r[melt_key] for r in rs]
        tf = {r["T_f_GL"] for r in rs}
        if max(tf) - min(tf) > 1e-9:
            raise RuntimeError(f"z_c {key}: T_f,GL differs between s = 1 runs: {sorted(tf)}")
        tf = rs[0]["T_f_GL"]
        tfs.append(tf)
        q = quadratic_fit(T, m)
        q["pass"] = bool(q["R2"] >= H1_R2_MIN)
        p = power_law_fit(T, m, tf)
        p["pass"] = bool(p["n"] is not None and H1_N_RANGE[0] <= p["n"] <= H1_N_RANGE[1])
        per[key] = {"verdict": "PASS" if q["pass"] and p["pass"] else "FAIL",
                    "n_runs": len(rs), "cases": [r["case"] for r in rs], "Tbot": T,
                    "melt": m, "quadratic": q, "power_law": p, "T_f_GL_degC": tf,
                    "quality": _quality(rs)}
    return {"hypothesis": "H1", "melt_key": melt_key,
            "rule": "pd0 H1: per z_c, 5 core runs (4 km, s = 1, Tbot -0.5..1.5, last 6 months); "
                    "quadratic R^2 >= 0.98 AND power-law n in [1.3, 2.0] "
                    "(m = k (Tbot - T_f,GL)^n, log-log fit)",
            "T_f_GL_degC": tfs[0] if tfs else None, "per_z_c": per}


def h2_verdicts(table: dict) -> dict:
    """pd0 H2, per z_c: Spearman rho(Tbot, D_out) <= -0.9 over the 5 core runs."""
    per = {}
    for z in CORE_ZC:
        rs, gaps = _core(table, z, ("D_out",))
        key = f"{z:.0f}"
        if gaps:
            per[key] = {"verdict": "NOT_EVALUABLE", "gaps": gaps}
            continue
        T = [r["Tbot"] for r in rs]
        D = [r["D_out"] for r in rs]
        rho = spearman(T, D)
        ok = rho is not None and round(rho, RHO_DECIMALS) <= H2_RHO_MAX
        per[key] = {"verdict": "PASS" if ok else "FAIL",
                    **({} if rho is not None else {"reason": "D_out constant: rho undefined"}),
                    "rho": rho, "n_runs": len(rs), "cases": [r["case"] for r in rs],
                    "Tbot": T, "D_out": D, "quality": _quality(rs)}
    return {"hypothesis": "H2",
            "rule": "pd0 H2: per z_c, 5 core runs (4 km, s = 1, last 6 months); "
                    "Spearman rho(Tbot, D_out) <= -0.9, D_out = -z_out at x_front (u C > 0)",
            "rho_rounding": f"rho rounded to {RHO_DECIMALS} decimals before the comparison",
            "per_z_c": per}


def n_summary(table: dict) -> dict:
    """N1-N3 over the ensemble (pd0), for the reports."""
    rs = rows(table)
    return {"n_members": len(rs),
            "gaps": {r["case"]: r["status"] for r in rs if r["status"] != "OK"},
            "N1_fail": [r["case"] for r in rs if r["N1_pass"] is False],
            "N2_fail": [r["case"] for r in rs if r["N2_pass"] is False],
            "N1_max": max((r["N1"] for r in rs if r["N1"] is not None), default=None),
            "N2_max_abs": max((abs(r["N2"]) for r in rs if r["N2"] is not None),
                              default=None)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Exp. 2 ensemble table (pd0 quantities)")
    ap.add_argument("--runs-root", type=Path, default=DEFAULT_RUNS)
    ap.add_argument("--run-suffix", default="_r21")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    tab = collect(a.runs_root, a.run_suffix)
    for r in rows(tab):
        print(f"{r['case']:28s} {r['split']:5s} {r['status'][:40]:40s} "
              f"m={r['mean_melt_shelf']} D_out={r['D_out']} N1={r['N1']} N2={r['N2']}")
    print(json.dumps(n_summary(tab)))
    if a.out:
        write_table(tab, a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
