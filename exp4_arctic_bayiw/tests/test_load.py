"""Tests for exp4_arctic_bayiw.analysis.load on a synthetic two-stage MDS run directory."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from exp4_arctic_bayiw.analysis import load

NZ, NY, NX = 2, 3, 4
RHONIL, TALPHA, SBETA, TREF, SREF = 999.8, 2e-4, 7.4e-4, 16.582, 33.6


def _meta(path: Path, shape, it, flds, nrec=None):
    dims = ",\n".join(f"{n:5d}, 1, {n:5d}" for n in reversed(shape))
    fl = " ".join(f"'{f:<8s}'" for f in flds)
    nrec = len(flds) if nrec is None else nrec
    path.write_text(f" nDims = [ {len(shape)} ];\n dimList = [\n{dims}\n ];\n"
                    f" dataprec = [ 'float64' ];\n nrecords = [ {nrec} ];\n"
                    f" timeStepNumber = [ {it} ];\n nFlds = [ {len(flds)} ];\n"
                    f" fldList = {{\n {fl}\n }};\n")


def _write(d: Path, prefix: str, it: int, arr, flds):
    """MDS diagnostics file; arr = (nrec, [nz,] ny, nx)."""
    stem = d / f"{prefix}.{it:010d}"
    np.ascontiguousarray(arr, ">f8").tofile(f"{stem}.data")
    _meta(Path(f"{stem}.meta"), arr.shape[1:], it, flds)


DATA = """ &PARM01
 eosType='LINEAR',
 tAlpha=2.E-4,
 sBeta=7.4E-4,
 rhoNil=999.8,
 gravity=9.81,
 tRef= 2*16.582,
 sRef= 2*33.6,
 &
 &PARM03
 nIter0={n0},
 nTimeSteps=60,
 deltaT=30.,
{forc} &
"""
FORC = " periodicExternalForcing=.TRUE.,\n externForcingPeriod=600.,\n externForcingCycle=1800.,\n"


def _field(it, k=0):
    return np.arange(NZ * NY * NX, dtype=float).reshape(NZ, NY, NX) * 0.01 + it + k


@pytest.fixture
def run(tmp_path):
    s1, s2 = tmp_path, tmp_path / "stage2"
    s2.mkdir()
    for d, n0, forc in ((s1, 0, FORC + " meridWindFile='tauy.bin',\n"), (s2, 60, "")):
        (d / "data").write_text(DATA.format(n0=n0, forc=forc))
        (d / "run_status.json").write_text(json.dumps(
            {"state": "OK", "last_iter": n0 + 60, "n_timesteps": 60}))
    for name in ("XC", "YC", "DXG", "DYG", "DXF", "DYF", "RAC", "RAS", "RAW", "Depth"):
        np.full((NY, NX), 2.0).astype(">f8").tofile(s1 / f"{name}.data")
        _meta(s1 / f"{name}.meta", (NY, NX), 0, [name], 1)
    for name in ("hFacC", "hFacW", "hFacS"):
        np.ones((NZ, NY, NX)).astype(">f8").tofile(s1 / f"{name}.data")
        _meta(s1 / f"{name}.meta", (NZ, NY, NX), 0, [name], 1)
    for name, v in (("DRF", [10.0, 20.0]), ("RC", [-5.0, -20.0])):
        np.array(v).reshape(NZ, 1, 1).astype(">f8").tofile(s1 / f"{name}.data")
        _meta(s1 / f"{name}.meta", (NZ, 1, 1), 0, [name], 1)
    np.full((3, NY, NX), 0.1).astype(">f8").tofile(s1 / "tauy.bin")
    for d, its in ((s1, (0, 30)), (s2, (60, 90))):
        for it in its:
            _write(d, "state3d", it, np.stack([_field(it, 0), _field(it, 1), _field(it, 2)]),
                   ["UVEL", "VVEL", "THETA"])
            _write(d, "diag_rho", it, _field(it, 3)[None], ["RHOAnoma"])
            _write(d, "state2d", it, _field(it, 4)[0][None], ["ETAN"])
    # final pickup (iteration 120): U, V, T, S, 4 tendencies (3-D), EtaN, dEtaHdt, EtaH (2-D)
    T, S = _field(120, 2), 33.0 + _field(120, 5) * 0.01
    recs = [_field(120, 0), _field(120, 1), T, S] + [np.zeros((NZ, NY, NX))] * 4
    flat = np.concatenate([r.reshape(-1, NY, NX) for r in recs] + [_field(120, 4)[:1]] * 3)
    flat.astype(">f8").tofile(s2 / "pickup.ckptA.data")
    _meta(s2 / "pickup.ckptA.meta", (NY, NX), 120,
          ["Uvel", "Vvel", "Theta", "Salt", "GuNm1", "GvNm1", "GtNm1", "GsNm1", "EtaN",
           "dEtaHdt", "EtaH"], flat.shape[0])
    return tmp_path, T, S


def test_run_reads_stages_lazily_and_appends_the_final_pickup(run):
    d, T, S = run
    ds = load.run(d)
    np.testing.assert_allclose(ds["time"].values, [0, 900, 1800, 2700, 3600])
    np.testing.assert_allclose(ds["U"].isel(time=2).values, _field(60, 0))
    np.testing.assert_allclose(ds["V"].isel(time=1, k=1, j=slice(0, 2)).values, _field(30, 1)[1, :2])
    np.testing.assert_allclose(ds["ETA"].isel(time=3).values, _field(90, 4)[0])
    # rho = rhoConst + RHOAnoma (model/src/find_rho.F:101-108; rhoConst = rhoNil by default)
    np.testing.assert_allclose(ds["RHO"].isel(time=0).values, RHONIL + _field(0, 3))
    # last record from the pickup; rho from the linear EOS of Theta, Salt
    assert ds.attrs["final_record_source"].startswith("pickup")
    rho_eos = RHONIL * (1 - TALPHA * (T - TREF) + SBETA * (S - SREF))
    np.testing.assert_allclose(ds["RHO"].isel(time=4).values, rho_eos, rtol=1e-14)
    np.testing.assert_allclose(ds["drF"].values, [10.0, 20.0])
    assert ds.attrs["rho0"] == RHONIL and ds.attrs["g"] == 9.81


def test_run_refuses_a_failed_stage(run):
    d, _, _ = run
    (d / "stage2" / "run_status.json").write_text(json.dumps(
        {"state": "FAILED", "last_iter": 70, "n_timesteps": 60}))
    with pytest.raises(RuntimeError):
        load.run(d)


def test_forcing_stress_interpolates_like_mitgcm_and_stops_in_stage2(run):
    d, _, _ = run
    ds = load.run(d)
    taux, tauy = load.forcing_stress(ds, 900.0)
    np.testing.assert_allclose(tauy, 0.1)
    assert np.all(taux == 0)
    _, tauy2 = load.forcing_stress(ds, 2700.0)              # stage 2 has no wind file
    assert np.all(tauy2 == 0)


def test_forcing_stress_wraps_on_the_namelist_cycle(run):
    """Periodic record sequence of eesupp/src/get_periodic_interval.F:83, 117-127: nbRec =
    cycle / period (here 3 of the 4 file records), record k centred at (k + 1/2) period."""
    d, _, _ = run
    R = np.stack([np.full((NY, NX), v) for v in (1.0, 2.0, 3.0, 99.0)])
    R.astype(">f8").tofile(d / "tauy.bin")
    ds = load.run(d)
    _, ty = load.forcing_stress(ds, 300.0)                  # centre of record 0
    np.testing.assert_allclose(ty, 1.0)
    _, ty = load.forcing_stress(ds, 1500.0)                 # centre of record 2
    np.testing.assert_allclose(ty, 3.0)
    _, ty = load.forcing_stress(ds, 1650.0)                 # 1/4 of the way from record 2 to 0
    np.testing.assert_allclose(ty, 0.75 * 3.0 + 0.25 * 1.0)  # wraps to record 0, not to 3 (99)


def test_pickup_layout_checks(tmp_path):
    from pfcommon.mitgcm_io import _parse_meta
    nz = 2
    a = np.zeros((2 * nz + 1, NY, NX))
    a.astype(">f8").tofile(tmp_path / "p.data")
    _meta(tmp_path / "p.meta", (NY, NX), 1, ["Uvel", "Vvel", "EtaN"], a.shape[0])
    f = load._pickup_fields(tmp_path / "p.data", _parse_meta((tmp_path / "p.meta").read_text()), nz)
    assert f["Uvel"].shape == (nz, NY, NX) and f["EtaN"].shape == (NY, NX)
    _meta(tmp_path / "p.meta", (NY, NX), 1, ["Uvel", "Vvel", "EtaN", "EtaH"], a.shape[0])
    with pytest.raises(ValueError, match="records"):
        load._pickup_fields(tmp_path / "p.data", _parse_meta((tmp_path / "p.meta").read_text()), nz)
    _meta(tmp_path / "p.meta", (NY, NX), 1, ["Uvel", "Vvel", "Mystery"], a.shape[0])
    with pytest.raises(ValueError, match="unknown"):
        load._pickup_fields(tmp_path / "p.data", _parse_meta((tmp_path / "p.meta").read_text()), nz)
