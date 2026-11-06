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
    w[10:, 2] = 0.05  # just above the 10 % threshold (0.04): inside the plume
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
    # R22: edge interpolated between centres of cell 2 (x = 12.5, w = 0.05) and
    # cell 3 (x = 17.5, w = 0): 12.5 + 5 (0.05 - 0.04)/0.05 = 13.5 m
    assert prof["width"][k] == pytest.approx(13.5)
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
        return SimpleNamespace(z=z, b=np.full_like(z, 8.0), Q=np.full_like(z, 2.0), u=np.full_like(z, 0.3), melt=np.full_like(z, 1e-5),
                               z_detach=zd, reached_top=False)

    monkeypatch.setattr(f1, "theory_solution", fake)
    zc = -(np.arange(50) + 0.5) * 10.0
    prof = {"z": zc, "z_detach": -300.0, "w_mean": np.full(50, 0.3), "melt": np.full(50, 1e-5),
            "w_face": np.full(50, 0.3), "w_max": np.full(50, 0.3), "width": np.full(50, 10.0),
            "Q": np.full(50, 3.0), "M": np.full(50, 0.9)}
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
        return SimpleNamespace(z=z, b=np.full_like(z, 8.0), Q=np.full_like(z, 2.0), u=np.full_like(z, u), melt=np.full_like(z, 1e-5),
                               z_detach=-300.0, reached_top=False)

    monkeypatch.setattr(f1, "theory_solution", fake)
    zc = -(np.arange(50) + 0.5) * 10.0
    prof = {"z": zc, "z_detach": -300.0, "w_mean": np.full(50, 0.3), "melt": np.full(50, 1e-5),
            "w_face": np.full(50, 0.3), "w_max": np.full(50, 0.3), "width": np.full(50, 10.0),
            "Q": np.full(50, 3.0), "M": np.full(50, 0.9)}
    out = f1.evaluate(prof, {"H": 500.0})
    assert out["modes"]["jenkins"]["u_ratio"] == pytest.approx(1.0)
    assert out["modes"]["mitgcm"]["u_ratio"] == pytest.approx(0.5)
    assert out["pass"]["C2p"] is False


def test_plume_extent_interpolated_edge_R22():
    from exp1_fjord_plumes.analysis.plume_diag import plume_extent

    x = np.array([2.5, 7.5, 12.5, 17.5])
    w = np.array([1.0, 1.0, 0.6, 0.0])
    width, Q, M = plume_extent(w, x, frac=0.1)
    edge = 12.5 + 5.0 * (0.6 - 0.1) / 0.6  # analytic crossing of 0.1 = 0.1 w_max
    assert width == pytest.approx(edge)
    # flat w = 1 on [0, 2.5], piecewise linear after; trapezoids (exact)
    q = 2.5 + 5.0 + 0.5 * (1.0 + 0.6) * 5.0 + 0.5 * (0.6 + 0.1) * (edge - 12.5)
    assert Q == pytest.approx(q)
    m = 2.5 + 5.0 + (1.0 + 0.6 + 0.36) / 3 * 5.0 + (0.36 + 0.06 + 0.01) / 3 * (edge - 12.5)
    assert M == pytest.approx(m)
    # whole row above threshold: edge is the last cell face (x_c[-1] + half cell is not
    # known to the helper, so it takes the domain edge passed in)
    width2, _, _ = plume_extent(np.array([1.0, 0.9, 0.8]), x[:3], frac=0.1, x_end=15.0)
    assert width2 == pytest.approx(15.0)


def test_post_hoc_and_declarations(monkeypatch):
    from types import SimpleNamespace

    from exp1_fjord_plumes.analysis import fig1_plume_validation as f1

    z = np.linspace(-500, 0, 501)
    calls = []

    def fake(params, mode, ambient=None):
        calls.append(ambient is not None)
        zd = -150.0 if ambient is not None else -100.0
        return SimpleNamespace(z=z, b=np.full_like(z, 8.0), Q=np.full_like(z, 2.0),
                               u=np.full_like(z, 0.3), melt=np.where(z < zd, 1e-5, np.nan),
                               z_detach=zd, reached_top=False)

    monkeypatch.setattr(f1, "theory_solution", fake)
    zc = -(np.arange(50) + 0.5) * 10.0
    prof = {"z": zc, "z_detach": -100.0, "w_mean": np.full(50, 0.3),
            "melt": np.where(zc < -100, 1e-5, 0.0),
            "w_face": np.full(50, 0.45), "w_max": np.full(50, 0.45), "width": np.full(50, 10.0),
            "Q": np.full(50, 3.0), "M": np.full(50, 0.9),
            "T_amb": np.full(50, 2.0), "S_amb": np.full(50, 34.0)}
    out = f1.evaluate(prof, {"H": 500.0})
    ph = out["post_hoc"]
    assert ph["w_face_peak_ratio"] == pytest.approx(1.5)
    assert ph["z_detach_theory_evolved"] == -150.0
    assert ph["z_detach_gap_over_H_evolved"] == pytest.approx(0.1)
    assert ph["melt_ratio_common_range"] == pytest.approx(1.0)
    assert ph["melt_ratio_whole_face"] == pytest.approx(1.0)
    assert "below its own z_detach" in f1.evaluate.__doc__
    assert out["melt_average_range"]


def _make_run3d(root: Path, nx=4, ny=6, nz=5, dy=10.0, dz=10.0, m0=2.0e-7, nan=False) -> Path:
    """Tiny 3-D run: uniform face melt m0 on all wet face cells (rows j = 0 and ny-1 are land)."""
    run = root / "run3d"
    run.mkdir()
    dt = 10.0
    (run / "run_status.json").write_text(json.dumps({"state": "OK", "last_iter": 720, "n_timesteps": 720, "wall_s": 1}))
    (run / "case.json").write_text(json.dumps({"params": {"deltaT": dt, "nx": nx, "ny": ny, "nz": nz}}))
    for name, arr in (("delX.bin", np.full(nx, 10.0)), ("delY.bin", np.full(ny, dy)), ("delZ.bin", np.full(nz, dz))):
        np.ascontiguousarray(arr, dtype=">f8").tofile(run / name)
    hfac = np.ones((nz, ny, nx))
    hfac[:, [0, -1], :] = 0.0
    hfac[-1, :, 0] = 0.5  # a partial bottom cell on the face
    np.ascontiguousarray(hfac, dtype=">f8").tofile(run / "hFacC.data")
    flux = np.zeros((2, nz, ny, nx))
    flux[0, :, 1:-1, 0] = -m0 * RHO_I_T
    flux[0, :, [0, -1], 0] = 0.0
    if nan:
        flux[:] = np.nan
    for it in (180, 360, 540, 720):  # 1800 s averages (timeInterval written like the model)
        stem = run / f"icefrontDiag.{it:010d}"
        np.ascontiguousarray(flux, dtype=">f8").tofile(f"{stem}.data")
        Path(f"{stem}.meta").write_text(
            f" nDims = [ 3 ];\n dimList = [\n {nx}, 1, {nx},\n {ny}, 1, {ny},\n {nz}, 1, {nz}\n ];\n"
            f" dataprec = [ 'float64' ];\n nrecords = [ 2 ];\n timeStepNumber = [ {it} ];\n"
            f" timeInterval = [ {it*dt-1800.0:.12E} {it*dt:.12E} ];\n nFlds = [ 2 ];\n"
            f" fldList = {{\n 'ICFfwFlx' 'ICFhtFlx'\n }};\n")
    return run


RHO_I_T = 917.0


def test_face_total_melt_uniform(tmp_path):
    from exp1_fjord_plumes.analysis.plume_diag import face_total_melt
    run = _make_run3d(tmp_path)
    # wet face area: 4 rows x (4 full levels + 0.5 partial) x dy dz = 4*4.5*100 = 1800 m2
    A = 4 * (4 + 0.5) * 10.0 * 10.0
    # 1800 s averages end at 0.5, 1, 1.5, 2 h; window (0.5, 2.0) h keeps those starting >= 0.5 h
    assert face_total_melt(run, (0.5, 2.0)) == pytest.approx(2.0e-7 * A, rel=1e-12)


def test_face_total_melt_nan_and_refuses(tmp_path):
    from exp1_fjord_plumes.analysis.plume_diag import face_total_melt
    run = _make_run3d(tmp_path, nan=True)
    assert np.isnan(face_total_melt(run, (0.5, 2.0)))
    with pytest.raises(ValueError):
        face_total_melt(run, (5.0, 6.0))  # no average inside the window
    (run / "run_status.json").write_text(json.dumps({"state": "FAILED"}))
    with pytest.raises(RuntimeError):
        face_total_melt(run, (0.5, 0.6))


def test_face_w_profile_max_over_y(tmp_path):
    from exp1_fjord_plumes.analysis.plume_diag import face_w_profile
    run = _make_run3d(tmp_path)
    nx, ny, nz = 4, 6, 5
    w = np.zeros((nz, ny, nx))
    w[2:, 3, 0] = 0.4  # plume in column j = 3 below level 2; another weaker plume at j = 4
    w[3:, 4, 0] = 0.1
    for it in (360, 720):
        stem = run / f"WVEL.{it:010d}"
        np.ascontiguousarray(w, dtype=">f8").tofile(f"{stem}.data")
        Path(f"{stem}.meta").write_text(
            f" nDims = [ 3 ];\n dimList = [\n {nx}, 1, {nx},\n {ny}, 1, {ny},\n {nz}, 1, {nz}\n ];\n"
            f" dataprec = [ 'float64' ];\n nrecords = [ 1 ];\n timeStepNumber = [ {it} ];\n")
    p = face_w_profile(run, (0.0, 2.0))
    assert p["w_ymax"] == pytest.approx([0, 0, 0.4, 0.4, 0.4])
    assert list(p["j_max"][2:]) == [3, 3, 3]
    assert p["t_snapshots_h"] == [1.0, 2.0]
    # detach_depth interpolates the 5 % crossing between top faces z = -20 m (w = 0.4) and -10 m (w = 0)
    assert p["z_detach"] == pytest.approx(-10.5)
