"""PD2 evaluation (notes/pd0_preregistro.md, "Enmienda PD2"; post hoc, criteria written before the
runs were seen). Applies the frozen thresholds mechanically and writes figures/validation.json
and figures/sensitivity.json.

A. deterministic 30-day window (GM30): V1d (GV1 adjoint grdchk vs FD, best of the eps per point,
   and vs the TLM of GV1_tlm), V2d/V3d (GV2d_<u>_<p|m> against xx_v2d/pred.json, J0 = fc GM30).
B. ensemble adjoint EA21-EA30 (90 d): V2e/V3e (V3E_<y>_<p|m>, J0 = fc EA<y>), H5a (adjoint
   verdict + direct H5A runs, years 26-30), H5c (leads 60/30/0 d), coarse H5b prediction
   (H5BC_<y>_<p|m>, 180 d, <(J+ - J-)/2>).

Usage: python -m exp5_adjoint_slope.analysis.pd2_eval"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from pfcommon.paths import EXP5_RUNS

from . import adjsens, ensemble, h5a, validate

FIG = Path(__file__).resolve().parents[1] / "figures"
YEARS = range(21, 31)
H5A_YEARS = range(26, 31)
_SUB = re.compile(r"^(?P<point>.+)_(?P<eps>e\d+)$")


def v1d_points(adj: dict, tlm: dict) -> list[dict]:
    """adj, tlm: {subrun name: grdchk entry}. Per control point the adjoint subrun with the
    smallest FD error; its gradient is compared with the TLM of the same point (any eps). The
    point's rel_err is the worse of the two comparisons (V1: adjoint vs FD and vs TLM)."""
    best: dict[str, tuple[str, dict]] = {}
    for name, g in adj.items():
        m = _SUB.match(name)
        pt, eps = m["point"], m["eps"]
        if pt not in best or g["rel_err"] < best[pt][1]["rel_err"]:
            best[pt] = (eps, g)
    t_by = {_SUB.match(n)["point"]: g for n, g in tlm.items()}
    out = []
    for pt, (eps, g) in best.items():
        rt = None
        if pt in t_by:
            rt = abs(g["grad"] - t_by[pt]["grad"]) / abs(t_by[pt]["grad"])
        out.append({"point": pt, "eps_run": eps, "grad_adj": g["grad"], "fd_grad": g["fd_grad"],
                    "grad_tlm": t_by[pt]["grad"] if pt in t_by else None,
                    "rel_err_fd": float(g["rel_err"]), "rel_err_tlm": rt,
                    "rel_err": float(max(g["rel_err"], rt if rt is not None else 0.0))})
    return out


def _fc(case: str) -> float:
    return float(json.loads((EXP5_RUNS / case / "fc.json").read_text())["fc"])


def _grd(case: str) -> dict:
    d = EXP5_RUNS / case
    return {s.name: json.loads((s / "grdchk.json").read_text())[0]
            for s in sorted(d.iterdir()) if (s / "grdchk.json").exists()}


def _pm(prefix: str, years, j0=None):
    """Per-member (d+, d-) with J0 (dict year -> J0), or (J+, J-) when j0 is None."""
    jp = np.array([_fc(f"{prefix}_{y}_p") for y in years])
    jm = np.array([_fc(f"{prefix}_{y}_m") for y in years])
    if j0 is None:
        return jp, jm
    b = np.array([j0[y] for y in years])
    return jp - b, jm - b


def part_a() -> dict:
    pts = v1d_points(_grd("GV1"), _grd("GV1_tlm"))
    out = {"V1d": {**validate.v1(pts), "points": pts}}
    j0 = _fc("GM30")
    pred = json.loads((EXP5_RUNS / "xx_v2d" / "pred.json").read_text())
    for u in ("fu", "theta", "diffkr"):
        dp, dm = _fc(f"GV2d_{u}_p") - j0, _fc(f"GV2d_{u}_m") - j0
        out[f"V2d_{u}"] = {**validate.v2(dp, dm), "dj_plus": dp, "dj_minus": dm}
        out[f"V3d_{u}"] = {**validate.v3(pred[u], (dp - dm) / 2), "pred": pred[u]}
    return out


def drop_periodic_end(S: np.ndarray) -> np.ndarray:
    """Drop the last gentim2d record of (members, nrec, ny) sensitivities: without pkg/cal the
    records are periodic over the window and the last one is the end point, which must carry no
    sensitivity (final review I1)."""
    if np.any(S[:, -1] != 0):
        raise ValueError("last record of adxx_fu is not zero; it cannot be dropped")
    return S[:, :-1]


def part_b() -> tuple[dict, dict]:
    pat = json.loads((EXP5_RUNS / "xx_ens" / "patterns.json").read_text())
    yc = np.array(pat["y_m"])
    S = np.array([adjsens.tau_x_sensitivity(adjsens.load_adxx(EXP5_RUNS / f"EA{y}", "fu"), None)
                  for y in YEARS])                                   # (members, nrec, ny)
    S = drop_periodic_end(S)                                         # last record: periodic end
    nrec = S.shape[1]
    leads = np.array([adjsens.lead_time_days(m, 90.0) for m in range(1, nrec + 1)])
    tot = np.abs(S.sum(-1))                                # |sum_y dJ/dtau_x| per record
    j0 = {y: _fc(f"EA{y}") for y in YEARS}
    out = {}
    dp, dm = _pm("V3E", YEARS, j0)
    out["V2e"] = {**ensemble.v2e(dp, dm), "dj_plus": dp.tolist(), "dj_minus": dm.tolist()}
    out["V3e"] = ensemble.v3e(float(np.mean(pat["v3e"]["pred_members"])), dp, dm)
    out["V3e"]["pred_members"] = pat["v3e"]["pred_members"]
    out["V3e"]["fwd_members"] = ((dp - dm) / 2).tolist()
    out["H5a"] = h5a.evaluate(S.sum(1), yc)
    direct = {}
    for n in ("i", "ii"):
        p, m = _pm(f"H5A{n}", H5A_YEARS, j0)
        direct[n] = {"fwd": ensemble.mean_se((p - m) / 2), "v2e": ensemble.v2e(p, m),
                     "pred": float(np.mean([pat[f"h5a{n}"]["pred_members"][y - 21] for y in H5A_YEARS]))}
    out["H5a"]["direct_37.5km_years26-30"] = direct
    order = np.argsort(leads)
    out["H5c"] = {**ensemble.h5c(leads[order], tot.mean(0)[order]),
                  "sens_mean": tot.mean(0)[order].tolist(),
                  "sens_se": (tot.std(0, ddof=1) / np.sqrt(len(tot)))[order].tolist()}
    jp, jm = _pm("H5BC", YEARS)
    half = (jp - jm) / 2
    out["H5b_coarse_prediction"] = {"dj_mean": float(half.mean()),
                                    "dj_se": float(half.std(ddof=1) / np.sqrt(half.size)),
                                    "members": half.tolist()}
    sens = {"years": list(YEARS), "leads_days": leads.tolist(), "y_m": yc.tolist(),
            "dJdtaux_rows_mean": S[:, :nrec].mean(0).tolist(),
            "dJdtaux_rows_se": (S[:, :nrec].std(0, ddof=1) / np.sqrt(len(S))).tolist(),
            "J0_members_TW": [j0[y] for y in YEARS], "fp_pattern": pat["v3e"]["profile"]}
    return out, sens


def main() -> None:
    val = {"note": "post hoc (PD2 amendment); thresholds as frozen in PD0", "A": part_a()}
    b, sens = part_b()
    val["B"] = b
    FIG.mkdir(exist_ok=True)
    (FIG / "validation.json").write_text(json.dumps(val, indent=1) + "\n")
    (FIG / "sensitivity.json").write_text(json.dumps(sens, indent=1) + "\n")
    print(json.dumps(val, indent=1))


if __name__ == "__main__":
    main()
