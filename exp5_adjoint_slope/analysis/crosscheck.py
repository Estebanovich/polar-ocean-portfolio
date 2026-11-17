"""H5b cross-check and V5 (plan Task 8; notes/pd0_preregistro.md, PD2: the coarse prediction is
the paired forward ensemble H5BC, not GB).

H5b: in the fine (eddying) model, the ensemble-mean dJ = <J_P> - <J_C> for the FP pattern holds if
  1. detectable: |dJ| > 2 sigma_D, sigma_D = sqrt(s_C^2/N + s_P^2/N) (sample std, ddof = 1);
  2. its sign is the predicted sign;
  3. it is within a factor 2 of the prediction.
With N < 3 the result is only "indicative".
V5: the eddy fraction of the fine-model J, v'theta' about the zonal and time mean (> 50 %, or
declared otherwise; Stewart et al. 2018).

Usage: python -m exp5_adjoint_slope.analysis.crosscheck   (writes figures/crosscheck.json)"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

FIG = Path(__file__).resolve().parents[1] / "figures"
TW = 1e-12


def sigma_delta(j_control, j_pert) -> float:
    c, p = np.asarray(j_control, float), np.asarray(j_pert, float)
    return float(np.sqrt(c.var(ddof=1) / c.size + p.var(ddof=1) / p.size))


def h5b(j_control, j_pert, dj_pred: float, factor: float = 2.0) -> dict:
    c, p = np.asarray(j_control, float), np.asarray(j_pert, float)
    dj = float(p.mean() - c.mean())
    sd = sigma_delta(c, p)
    det = bool(abs(dj) > 2 * sd)
    sign = bool(np.sign(dj) == np.sign(dj_pred) and dj != 0)
    ratio = dj / dj_pred if dj_pred != 0 else float("nan")
    within = bool(sign and 1 / factor <= ratio <= factor)
    return {"dj": dj, "sigma_delta": sd, "dj_pred": float(dj_pred), "ratio": float(ratio),
            "detectable": det, "sign_ok": sign, "within_factor": within,
            "pass": bool(det and sign and within), "indicative": bool(min(c.size, p.size) < 3),
            "n": int(min(c.size, p.size))}


def _face(theta, j):
    th = np.asarray(theta, float)
    return 0.5 * (th[..., j, :] + th[..., j - 1, :])


def _w(grid, j):
    hf = np.asarray(grid["hfac_s"], float)[:, j, :]
    return hf * np.asarray(grid["drf"], float)[:, None] * np.asarray(grid["dxg"], float) * (hf > 0)


def eddy_fraction(v, theta, grid, j: int, theta_f: float = -1.9) -> float:
    """v, theta: (nt, nr, ny, nx). Fraction of the transport sum v (theta_face - theta_f) carried
    by deviations from the zonal-and-time mean at the faces of row j."""
    vj = np.asarray(v, float)[:, :, j, :]
    tf = _face(theta, j) - theta_f
    w = _w(grid, j)
    total = float((vj * tf * w).mean(0).sum())
    vm = vj.mean(axis=(0, 2), keepdims=True)
    tm = tf.mean(axis=(0, 2), keepdims=True)
    mean = float((vm * tm * w).mean(0).sum())
    return (total - mean) / total


def eddy_fraction_from_means(vvelth, vvel, theta, grid, j: int, theta_f: float = -1.9) -> float:
    """Same quantity from time means: VVELTH (v theta at the v face, per step), VVEL and THETA
    records (nt, nr, ny, nx); the mean-flow part uses the zonal-and-time mean v and theta."""
    w = _w(grid, j)
    vt = np.asarray(vvelth, float)[:, :, j, :] - theta_f * np.asarray(vvel, float)[:, :, j, :]
    total = float((vt * w).mean(0).sum())
    vm = np.asarray(vvel, float)[:, :, j, :].mean(axis=(0, 2), keepdims=True)
    tm = (_face(theta, j) - theta_f).mean(axis=(0, 2), keepdims=True)
    mean = float((vm[0] * tm[0] * w).sum())
    return (total - mean) / total


def eke(u, v, usq, vsq, wet) -> float:
    """Volume-mean (wet cells, equal weights) eddy kinetic energy (m2 s-2) of one monthly record:
    1/2 <UVELSQ - ubar^2 + VVELSQ - vbar^2>, ubar/vbar the zonal means of the monthly UVEL/VVEL
    (correction 2026-11-17, notes/pd0_preregistro.md). Arrays (nr, ny, nx)."""
    u, v = np.asarray(u, float), np.asarray(v, float)
    ub, vb = u.mean(-1, keepdims=True), v.mean(-1, keepdims=True)
    e = 0.5 * (np.asarray(usq, float) - ub ** 2 + np.asarray(vsq, float) - vb ** 2)
    return float(e[np.asarray(wet, bool)].mean())


def eke_equilibrated(series, n: int = 3, tol: float = 0.1) -> dict:
    """Equilibrium rule: mean of the last n months within tol of the mean of the n before."""
    s = np.asarray(series, float)
    last, prev = s[-n:].mean(), s[-2 * n:-n].mean()
    rel = float(abs(last - prev) / prev)
    return {"last": float(last), "prev": float(prev), "rel_change": rel, "pass": bool(rel < tol)}


def eke_series(run: Path) -> list[float]:
    """Monthly EKE of a fine run from state3d (THETA, SALT, UVEL, VVEL, UVELSQ, VVELSQ)."""
    from pfcommon.mitgcm_io import read_bin, read_mds
    case = json.loads((run / "input" / "case.json").read_text())
    nr, ny, nx = case["params"]["nz"], case["params"]["ny"], case["params"]["nx"]
    wet = read_bin(run / "hFacC.data", (nr, ny, nx)) > 0
    its = sorted(int(f.name.split(".")[1]) for f in run.glob("state3d.*.data") if not f.name.startswith("._"))
    out = []
    for it in its:
        a = np.asarray(read_mds(run / "state3d", it)[0])
        out.append(eke(a[2], a[3], a[4], a[5], wet))
    return out


def fine_J(run: Path, n_last: int = 6) -> float:
    """J (TW) of a fine run: mean of the last n_last 5-day jflux5d records (the last 30 days)."""
    from pfcommon.campaign import nml_float
    from pfcommon.mitgcm_io import read_bin, read_mds
    from . import heat_transport as ht
    case = json.loads((run / "input" / "case.json").read_text())
    nr, ny, nx = case["params"]["nz"], case["params"]["ny"], case["params"]["nx"]
    j = case["j_qb"]
    its = sorted(int(f.name.split(".")[1]) for f in run.glob("jflux5d.*.data")
                 if not f.name.startswith("._"))[-n_last:]
    m = np.mean([np.asarray(read_mds(run / "jflux5d", it)[0]) for it in its], axis=0)
    data = (run / "input" / "data").read_text()
    return ht.J_from_means(m[0], m[1], read_bin(run / "hFacS.data", (nr, ny, nx)),
                           read_bin(run / "DRF.data", (nr,)), read_bin(run / "DXG.data", (ny, nx))[j],
                           j, theta_f=case["theta_f"], rho0=nml_float(data, "rhoConst"),
                           cp=nml_float(data, "HeatCapacity_Cp"))


def v5(run: Path, n_months: int = 6) -> float:
    """Eddy fraction of J over the last n_months months of a fine control run."""
    from pfcommon.mitgcm_io import read_bin, read_mds
    case = json.loads((run / "input" / "case.json").read_text())
    nr, ny, nx = case["params"]["nz"], case["params"]["ny"], case["params"]["nx"]
    j = case["j_qb"]
    st = sorted(int(f.name.split(".")[1]) for f in run.glob("state3d.*.data")
                if not f.name.startswith("._"))[-n_months:]
    jf = sorted(int(f.name.split(".")[1]) for f in run.glob("jflux5d.*.data")
                if not f.name.startswith("._"))[-6 * n_months:]
    th = np.array([np.asarray(read_mds(run / "state3d", it)[0])[0] for it in st])
    fl = np.array([np.asarray(read_mds(run / "jflux5d", it)[0]) for it in jf])
    grid = {"hfac_s": read_bin(run / "hFacS.data", (nr, ny, nx)), "drf": read_bin(run / "DRF.data", (nr,)),
            "dxg": read_bin(run / "DXG.data", (ny, nx))[j]}
    return eddy_fraction_from_means(fl[:, 0], fl[:, 1], th, grid, j, case["theta_f"])


def main() -> None:
    from pfcommon.paths import EXP5_RUNS
    val = json.loads((FIG / "validation.json").read_text())["B"]["H5b_coarse_prediction"]
    jc = {c: fine_J(EXP5_RUNS / c) for c in ("FC1", "FC2", "FC3") if (EXP5_RUNS / c / "run_status.json").exists()}
    jp = {c: fine_J(EXP5_RUNS / c) for c in ("FP1", "FP2", "FP3") if (EXP5_RUNS / c / "run_status.json").exists()}
    spin = {c: eke_series(EXP5_RUNS / c) for c in ("F0b", "F1", "F1b") if (EXP5_RUNS / c / "run_status.json").exists()}
    out = {"EKE_spinup_m2s2": spin,
           "F1_equilibrium": eke_equilibrated(spin["F1"]) if "F1" in spin else None,
           "F1b_equilibrium": eke_equilibrated(spin["F1b"]) if "F1b" in spin else None,
           "J_control_TW": jc, "J_pert_TW": jp,
           "H5b": h5b(list(jc.values()), list(jp.values()), val["dj_mean"]),
           "dj_pred_se_TW": val["dj_se"], "V5_eddy_fraction": {c: v5(EXP5_RUNS / c) for c in jc}}
    (FIG / "crosscheck.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
