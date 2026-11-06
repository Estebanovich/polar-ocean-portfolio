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


def _meta3d(nx, ny, nz, nrec, it, extra=""):
    return (f" nDims = [ 3 ];\n dimList = [\n {nx}, 1, {nx},\n {ny}, 1, {ny},\n {nz}, 1, {nz}\n ];\n"
            f" dataprec = [ 'float64' ];\n nrecords = [ {nrec} ];\n timeStepNumber = [ {it} ];\n{extra}")


RHO_I_T = 917.0
NX3, NY3, NZ3, DT3 = 4, 6, 5, 10.0


def _make_run3d(root: Path, m0s=(2.0e-7,) * 4, period=1800.0, interval=True) -> Path:
    """Tiny 3-D run. Face melt is uniform on the wet face cells (rows j = 0, ny-1 are land; the
    bottom face cell is half-wet). One icefrontDiag average per entry of ``m0s``, the k-th ending at
    (k+1) * period (so 0.5, 1.0, 1.5, 2.0 h for period = 1800 s) with melt rate m0s[k]."""
    run = root / "run3d"
    run.mkdir()
    (run / "run_status.json").write_text(json.dumps({"state": "OK", "last_iter": 720, "n_timesteps": 720, "wall_s": 1}))
    (run / "case.json").write_text(json.dumps({"params": {"deltaT": DT3, "nx": NX3, "ny": NY3, "nz": NZ3}}))
    for name, arr in (("delX.bin", np.full(NX3, 10.0)), ("delY.bin", np.full(NY3, 10.0)),
                      ("delZ.bin", np.full(NZ3, 10.0))):
        np.ascontiguousarray(arr, dtype=">f8").tofile(run / name)
    hfac = np.ones((NZ3, NY3, NX3))
    hfac[:, [0, -1], :] = 0.0
    hfac[-1, :, 0] = 0.5
    np.ascontiguousarray(hfac, dtype=">f8").tofile(run / "hFacC.data")
    for k, m0 in enumerate(m0s):
        it = int(round((k + 1) * period / DT3))
        flux = np.zeros((2, NZ3, NY3, NX3))
        flux[0, :, 1:-1, 0] = -m0 * RHO_I_T
        stem = run / f"icefrontDiag.{it:010d}"
        np.ascontiguousarray(flux, dtype=">f8").tofile(f"{stem}.data")
        ti = f" timeInterval = [ {it*DT3-period:.12E} {it*DT3:.12E} ];\n" if interval else ""
        Path(f"{stem}.meta").write_text(_meta3d(NX3, NY3, NZ3, 2, it, ti + " nFlds = [ 2 ];\n fldList = {\n 'ICFfwFlx' 'ICFhtFlx'\n };\n"))
    return run


FACE_AREA3 = 4 * (4 + 0.5) * 10.0 * 10.0  # wet face: 4 rows x (4 full + 1 half level) x dy dz = 1800 m2


def test_face_total_melt_uniform(tmp_path):
    from exp1_fjord_plumes.analysis.plume_diag import face_total_melt
    run = _make_run3d(tmp_path)
    assert face_total_melt(run, (0.5, 2.0)) == pytest.approx(2.0e-7 * FACE_AREA3, rel=1e-12)


def test_face_total_melt_window_selection(tmp_path):
    """Different melt in each average (ending 0.5, 1, 1.5, 2 h; covering (end-0.5 h, end]).
    Window (1, 2] h must use exactly the averages ending at 1.5 and 2.0 h: not the one ending
    at 1.0 h (covers 0.5-1.0 h) and not the one ending at 0.5 h."""
    from exp1_fjord_plumes.analysis.plume_diag import face_total_melt
    m0s = (1.0e-7, 2.0e-7, 4.0e-7, 8.0e-7)
    run = _make_run3d(tmp_path, m0s=m0s)
    assert face_total_melt(run, (1.0, 2.0)) == pytest.approx(np.mean(m0s[2:]) * FACE_AREA3, rel=1e-12)
    assert face_total_melt(run, (0.5, 1.5)) == pytest.approx(np.mean(m0s[1:3]) * FACE_AREA3, rel=1e-12)
    assert face_total_melt(run, (0.0, 2.0)) == pytest.approx(np.mean(m0s) * FACE_AREA3, rel=1e-12)
    with pytest.raises(ValueError):
        face_total_melt(run, (2.0, 3.0))  # nothing averaged entirely inside (2, 3] h


def test_face_total_melt_reads_period_from_meta(tmp_path):
    """900 s averages (as in (b)): ends 0.25, 0.5, 0.75, 1.0 h; window (0.5, 1.0] -> the last two."""
    from exp1_fjord_plumes.analysis.plume_diag import face_total_melt
    m0s = (1.0e-7, 2.0e-7, 4.0e-7, 8.0e-7)
    run = _make_run3d(tmp_path, m0s=m0s, period=900.0)
    assert face_total_melt(run, (0.5, 1.0)) == pytest.approx(np.mean(m0s[2:]) * FACE_AREA3, rel=1e-12)


def test_face_total_melt_missing_timeinterval_raises(tmp_path):
    from exp1_fjord_plumes.analysis.plume_diag import face_total_melt
    run = _make_run3d(tmp_path, interval=False)
    with pytest.raises(ValueError, match="timeInterval"):
        face_total_melt(run, (0.5, 2.0))


def test_face_total_melt_nonfinite_and_refuses(tmp_path):
    from exp1_fjord_plumes.analysis.plume_diag import face_total_melt
    run = _make_run3d(tmp_path, m0s=(1e-7, 1e-7, float("nan"), 1e-7))
    with pytest.raises(ValueError, match="run3d"):  # names the run
        face_total_melt(run, (0.5, 2.0))
    assert face_total_melt(run, (0.5, 1.0)) == pytest.approx(1e-7 * FACE_AREA3)  # NaN file outside window
    (run / "run_status.json").write_text(json.dumps({"state": "FAILED"}))
    with pytest.raises(RuntimeError):
        face_total_melt(run, (0.5, 1.0))


def _write_field(run: Path, name: str, it: int, arr):
    stem = run / f"{name}.{it:010d}"
    np.ascontiguousarray(arr, dtype=">f8").tofile(f"{stem}.data")
    Path(f"{stem}.meta").write_text(_meta3d(NX3, NY3, NZ3, 1, it))


def test_face_w_profile_time_mean_then_max_over_y(tmp_path):
    """Two snapshots whose maxima are in different j: window-mean first (0.2 at j=2 and j=4 ->
    mean 0.1 each... plus a steady 0.15 at j=3) then max over y. max-then-mean would give 0.3."""
    from exp1_fjord_plumes.analysis.plume_diag import face_w_profile
    run = _make_run3d(tmp_path)
    w1 = np.zeros((NZ3, NY3, NX3)); w2 = np.zeros_like(w1)
    w1[2:, 2, 0] = 0.4; w2[2:, 4, 0] = 0.4      # moving plume: window mean 0.2 in each of j = 2, 4
    w1[2:, 3, 0] = 0.15; w2[2:, 3, 0] = 0.15    # steady plume: mean 0.15 at j = 3
    _write_field(run, "WVEL", 360, w1)
    _write_field(run, "WVEL", 720, w2)
    p = face_w_profile(run, (0.0, 2.0))
    assert p["w_ymax"] == pytest.approx([0, 0, 0.2, 0.2, 0.2])  # time-mean then max (reverse: 0.4)
    assert list(p["j_max"][2:]) == [2, 2, 2]
    assert p["t_snapshots_h"] == [1.0, 2.0]
    # detach_depth: 5 % crossing between top faces z = -20 m (0.2) and -10 m (0)
    assert p["z_detach"] == pytest.approx(-10.5)


def test_face_w_profile_nonfinite_raises(tmp_path):
    from exp1_fjord_plumes.analysis.plume_diag import face_w_profile
    run = _make_run3d(tmp_path)
    w = np.zeros((NZ3, NY3, NX3)); w[1, 2, 0] = np.nan
    _write_field(run, "WVEL", 360, w)
    with pytest.raises(ValueError, match="run3d"):
        face_w_profile(run, (0.0, 2.0))


def test_ambient_far_from_face(tmp_path):
    """Window mean over snapshots, then mean over the WET y cells of the column nearest x_m
    (land rows j = 0, ny-1 have hFacC = 0 and must not enter)."""
    from exp1_fjord_plumes.analysis.plume_diag import ambient_far_from_face
    run = _make_run3d(tmp_path)
    T1 = np.full((NZ3, NY3, NX3), 1.0); T2 = np.full_like(T1, 3.0)
    for T in (T1, T2):
        T[:, [0, -1], :] = 999.0  # land rows: must be ignored
    T1[0, 1:-1, 2] = 5.0  # level 0 of column i = 2, wet rows: window mean (5+3)/2 = 4
    S = np.full((NZ3, NY3, NX3), 34.0)
    for it, T in ((360, T1), (720, T2)):
        _write_field(run, "THETA", it, T)
        _write_field(run, "SALT", it, S)
    a = ambient_far_from_face(run, (0.0, 2.0), x_m=24.0)  # centres 5, 15, 25, 35 m -> i = 2
    assert a["i"] == 2 and a["x"] == pytest.approx(25.0)
    assert a["T"] == pytest.approx([4.0, 2.0, 2.0, 2.0, 2.0])
    assert a["S"] == pytest.approx(np.full(NZ3, 34.0))
    assert a["z"] == pytest.approx([-5, -15, -25, -35, -45])
