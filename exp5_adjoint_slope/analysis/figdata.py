"""Figure data for the README (plan Task 7), extracted from the runs into committed JSON:

  figures/mean_state.json   G0 year 30: x-mean theta, S and sigma (linear EOS) sections, depth
                            profile, wind, J and its parts (G0/screen_v4.json).
  figures/diagnostics.json  why the adjoint became an ensemble (PD2): max |dJ/dtheta0| against the
                            window (GM30, D90, GB, GA; theta is unaffected by the wind-control
                            correction), the monthly J of the E0-E4 forward ensemble and the D90
                            finite-difference jump.

Runs are read from $EXP5_RUNS (G0 from the LaCie copy when not on the run disk).
Usage: python -m exp5_adjoint_slope.analysis.figdata"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pfcommon.campaign import nml_float
from pfcommon.mitgcm_io import read_bin, read_mds
from pfcommon.paths import EXP5_RUNS

from ..gendata import fields
from . import adjsens
from . import heat_transport as ht

FIG = Path(__file__).resolve().parents[1] / "figures"
LACIE = Path("/Volumes/LaCie/portfolio_runs/exp5")
ALPHA, BETA, RHO0 = 3.0e-5, 7.4e-4, 1027.5           # linear EOS of the configuration (data)
NR, NY, NX = 35, 45, 40


def _run(case: str) -> Path:
    d = EXP5_RUNS / case
    return d if (d / "run_status.json").exists() else LACIE / case


def xmean(a, hfac):
    """Mean over wet cells along x; NaN where the row is dry."""
    w = np.asarray(hfac, float) > 0
    s = np.where(w, a, 0.0).sum(-1)
    n = w.sum(-1)
    return np.where(n > 0, s / np.maximum(n, 1), np.nan)


def monthly_J(run: Path) -> list[float]:
    case = json.loads((run / "case.json").read_text()) if (run / "case.json").exists() \
        else json.loads((run / "input" / "case.json").read_text())
    j = case["j_qb"]
    hfs = read_bin(run / "hFacS.data", (NR, NY, NX))
    drf = read_bin(run / "DRF.data", (NR,))
    dxg = read_bin(run / "DXG.data", (NY, NX))[j]
    data = (run / "data").read_text()
    kw = dict(theta_f=case["theta_f"], rho0=nml_float(data, "rhoConst"),
              cp=nml_float(data, "HeatCapacity_Cp"))
    its = sorted(int(f.name.split(".")[1]) for f in run.glob("state3d.*.data")
                 if not f.name.startswith("._"))
    out = []
    for it in its:
        m = np.asarray(read_mds(run / "state3d", it)[0])
        out.append(ht.J_from_means(m[4], m[3], hfs, drf, dxg, j, dfy=m[5], vbt=m[6], **kw))
    return out


def mean_state() -> dict:
    g0 = _run("G0")
    st = g0 / "stage30"
    its = sorted(int(f.name.split(".")[1]) for f in st.glob("state3d.*.data")
                 if not f.name.startswith("._"))[-12:]
    m = np.mean([np.asarray(read_mds(st / "state3d", it)[0]) for it in its], axis=0)
    hfc = read_bin(g0 / "hFacC.data", (NR, NY, NX))
    th, sa = xmean(m[0], hfc), xmean(m[1], hfc)
    sig = RHO0 * (-ALPHA * (th - 0.0) + BETA * (sa - 35.0))      # anomaly, kg m-3 (linear EOS)
    yc = read_bin(g0 / "YC.data", (NY, NX))[:, 0]
    rc = read_bin(g0 / "RC.data", (NR,))
    case = json.loads((g0 / "case.json").read_text())
    scr = json.loads((g0 / "screen_v4.json").read_text())
    tolist = lambda a: np.where(np.isnan(a), None, a).tolist()   # noqa: E731
    return {"source": f"G0 {st.name}, 12 monthly means, iters {its[0]}-{its[-1]}",
            "y_m": yc.tolist(), "z_m": rc.tolist(), "depth_m": fields.depth(yc).tolist(),
            "tau_x": fields.wind_tau_x(yc).tolist(), "theta": tolist(th), "salt": tolist(sa),
            "sigma_anom": tolist(sig), "y_qb_m": case["y_qb_m"], "j_qb": case["j_qb"],
            "J": {k: scr[k] for k in ("J_TW", "J_resolved_TW", "J_bolus_TW", "J_redi_TW",
                                     "shelf_theta_degC", "shelf_salt")}}


def diagnostics() -> dict:
    growth = {}
    for case, win in (("GM30", 30), ("D90/xx_theta_e0", 90), ("GB", 180), ("GA", 360)):
        d = EXP5_RUNS / case
        th = np.abs(adjsens.load_adxx(d, "theta"))
        growth[str(win)] = {"case": case, "max_abs_dJdtheta0": float(th.max()),
                            "median_abs_dJdtheta0": float(np.median(th[th > 0]))}
    ens = {c: monthly_J(EXP5_RUNS / c) for c in ("E0", "E1", "E2", "E3", "E4")}
    a = np.array(list(ens.values()))
    g = json.loads((EXP5_RUNS / "D90" / "xx_theta_e0" / "grdchk.json").read_text())[0]
    return {"adjoint_growth": growth, "ensemble_J_TW": ens,
            "ensemble_std_TW": a.std(0, ddof=1).tolist(),
            "d90_fd": {k: g[k] for k in ("fc_ref", "fc_plus", "fc_minus", "grad", "fd_grad")},
            "theta_noise_degC": 1e-6}


def main() -> None:
    FIG.mkdir(exist_ok=True)
    for name, fn in (("mean_state", mean_state), ("diagnostics", diagnostics)):
        (FIG / f"{name}.json").write_text(json.dumps(fn(), indent=1) + "\n")
        print("wrote", FIG / f"{name}.json")


if __name__ == "__main__":
    main()
