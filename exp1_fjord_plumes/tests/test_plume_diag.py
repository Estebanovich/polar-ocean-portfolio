"""Tests for the plume diagnostics of the (a) line-plume runs."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from exp1_fjord_plumes.analysis.plume_diag import (
    detach_depth,
    melt_rate_m_per_s,
    plume_profile,
)

NZ, NX, DT, DZ = 20, 12, 10.0, 10.0


def test_melt_conversion_known_value():
    # -1e-4 kg m-2 s-1 (ICFfwFlx < 0 = ice loss) -> +1.0905e-7 m/s = 3.44 m/yr
    m = melt_rate_m_per_s(np.array([-1.0e-4, 0.0, 2.0e-4]))
    assert m[0] == pytest.approx(1.0905e-7, rel=1e-4)
    assert m[0] * 365.25 * 86400 == pytest.approx(3.44, abs=0.01)
    assert m[1] == 0.0 and m[2] < 0  # freezing is negative


def _write_mds(prefix: Path, it: int, arr: np.ndarray, flds: list[str] | None = None) -> None:
    arr = np.asarray(arr, dtype=float)
    nrec = arr.shape[0] if flds and len(flds) > 1 else 1
    body = arr.reshape(nrec, NZ, NX)
    stem = prefix.parent / f"{prefix.name}.{it:010d}"
    np.ascontiguousarray(body, dtype=">f8").tofile(f"{stem}.data")
    fl = "".join(f" '{f}'" for f in (flds or []))
    Path(f"{stem}.meta").write_text(
        f" nDims = [ 3 ];\n dimList = [\n {NX}, 1, {NX},\n 1, 1, 1,\n {NZ}, 1, {NZ}\n ];\n"
        f" dataprec = [ 'float64' ];\n nrecords = [ {nrec} ];\n timeStepNumber = [ {it} ];\n"
        f" nFlds = [ {nrec} ];\n fldList = {{\n{fl}\n }};\n"
    )


def _make_run(root: Path, state: str = "OK") -> Path:
    run = root / "run"
    run.mkdir()
    (run / "run_status.json").write_text(
        json.dumps({"state": state, "last_iter": 720, "n_timesteps": 720, "wall_s": 1})
    )
    (run / "case.json").write_text(json.dumps({"params": {"deltaT": DT, "nx": NX, "nz": NZ}}))
    np.ascontiguousarray(np.full(NX, 5.0), dtype=">f8").tofile(run / "delX.bin")
    np.ascontiguousarray(np.full(NZ, DZ), dtype=">f8").tofile(run / "delZ.bin")
    # plume: w = 0.4 m/s in face column 0..1 below k=10 (z=-100 m), 0 above
    w = np.zeros((NZ, NX))
    w[10:, :2] = 0.4
    w[10:, 2] = 0.05  # below the 10 % threshold (0.04)? no: above -> inside plume
    for it in (360, 720):
        _write_mds(run / "WVEL", it, w)
        _write_mds(run / "THETA", it, np.full((NZ, NX), 1.5))
        _write_mds(run / "SALT", it, np.full((NZ, NX), 34.0))
        _write_mds(run / "UVEL", it, np.full((NZ, NX), 0.01))
    fw = np.zeros((2, NZ, NX))
    fw[0, 10:, 0] = -1.0e-4
    _write_mds(run / "icefrontDiag", 360, fw, ["ICFfwFlx", "ICFhtFlx"])
    _write_mds(run / "icefrontDiag", 720, fw, ["ICFfwFlx", "ICFhtFlx"])
    return run


def test_plume_profile_refuses_failed_run(tmp_path):
    run = _make_run(tmp_path, state="FAILED")
    with pytest.raises(RuntimeError):
        plume_profile(run, (0.0, 0.2))


def test_detach_depth_synthetic():
    dz = 10.0
    zc = -(np.arange(60) + 0.5) * dz  # surface first, z negative down
    w = np.where(zc < -150.0, 0.3 * np.sin(np.pi * (zc + 600) / 600), 0.0)
    zd = detach_depth(zc, w)
    assert abs(zd - (-150.0)) <= dz


def test_plume_profile_synthetic_run(tmp_path):
    run = _make_run(tmp_path)
    # hours here are 3600 s, steps 10 s: files at it=360 -> t = 1 h, it=720 -> 2 h
    prof = plume_profile(run, (0.0, 2.0))
    assert prof["z_detach"] == pytest.approx(-100.0, abs=DZ)
    k = 15
    assert prof["w_face"][k] == pytest.approx(0.4)
    assert prof["melt"][k] == pytest.approx(1.0905e-7, rel=1e-4)  # positive = melt
    assert prof["melt"][0] == 0.0
    # plume = face + cell 1 + cell 2 (0.05 >= 0.1*0.4 = 0.04): width 15 m
    assert prof["width"][k] == pytest.approx(15.0)
    assert np.isnan(prof["width"][0])  # above the detachment depth
    assert prof["T"][k] == pytest.approx(1.5)
    assert prof["t_snapshots_h"] == [1.0, 2.0]
    assert prof["t_means_h"] == [1.0, 2.0]


def test_fig1_verdict_uses_mitgcm_mode_R13(monkeypatch):
    from types import SimpleNamespace

    from exp1_fjord_plumes.analysis import fig1_plume_validation as f1

    z = np.linspace(-500, 0, 501)

    def fake(params, mode):  # jenkins would pass C2, mitgcm must fail it
        zd = -300.0 if mode == "jenkins" else -100.0
        return SimpleNamespace(z=z, u=np.full_like(z, 0.3), melt=np.full_like(z, 1e-5),
                               z_detach=zd, reached_top=False)

    monkeypatch.setattr(f1, "theory_solution", fake)
    zc = -(np.arange(50) + 0.5) * 10.0
    prof = {"z": zc, "z_detach": -300.0, "w_mean": np.full(50, 0.3), "melt": np.full(50, 1e-5)}
    out = f1.evaluate(prof, {"H": 500.0})
    assert out["pass"]["C2"] is False
    assert out["modes"]["jenkins"]["z_detach_diff_over_H"] == pytest.approx(0.0)


def test_fig1_C2p_verdict_uses_mitgcm_mode_R13(monkeypatch):
    """R13 for C2': same detachment in both modes; jenkins u would pass C2',
    the mitgcm-mode u (2x the model) must fail it."""
    from types import SimpleNamespace

    from exp1_fjord_plumes.analysis import fig1_plume_validation as f1

    z = np.linspace(-500, 0, 501)

    def fake(params, mode):
        u = 0.3 if mode == "jenkins" else 0.6
        return SimpleNamespace(z=z, u=np.full_like(z, u), melt=np.full_like(z, 1e-5),
                               z_detach=-300.0, reached_top=False)

    monkeypatch.setattr(f1, "theory_solution", fake)
    zc = -(np.arange(50) + 0.5) * 10.0
    prof = {"z": zc, "z_detach": -300.0, "w_mean": np.full(50, 0.3), "melt": np.full(50, 1e-5)}
    out = f1.evaluate(prof, {"H": 500.0})
    assert out["modes"]["jenkins"]["u_ratio"] == pytest.approx(1.0)
    assert out["modes"]["mitgcm"]["u_ratio"] == pytest.approx(0.5)
    assert out["pass"]["C2p"] is False
