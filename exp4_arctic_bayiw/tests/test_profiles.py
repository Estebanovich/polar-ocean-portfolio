"""Tests for exp4_arctic_bayiw.gendata.profiles and analysis.theory_numbers (Task 2).

Synthetic tests need no external data. The tests marked "real data" use the UNIS HD profiles,
IBCAO/NPI files in data_external/ and the read-only thesis inputs, and are skipped if missing.
"""
from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np
import pytest
import yaml

from exp4_arctic_bayiw.gendata import profiles, vgrid
from exp4_arctic_bayiw.theory import modes

PKG = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((PKG / "cases.yaml").read_text())["configs"]["fjord_vmf"]
EOS = profiles.EOS_EXP4
UNIS = PKG / "data_external" / "unis_hd" / "profiles"
THESIS_IN = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/input")
G = 9.81


def make_prof(z, T, S, lat=77.80, lon=15.5, time=dt.datetime(2010, 8, 15)):
    z = np.asarray(z, float)
    return {"lat": lat, "lon": lon, "time": time, "z": z, "p": -z, "T": np.asarray(T, float),
            "S": np.asarray(S, float), "file": "synthetic"}


def summer_like(zmax=100.0, dz=1.0, dS=0.0, dT=0.0):
    """Stable fjord-like profile: fresh warm surface layer over saltier cold water."""
    z = -np.arange(1.0, zmax + 0.5, dz)
    S = 34.4 - 3.0 * np.exp(z / 8.0) + dS          # z < 0: S increases with depth
    T = 1.0 + 4.0 * np.exp(z / 10.0) + dT
    return make_prof(z, T, S)


def n2_linear(zc, T, S):
    """N^2 = g (alpha dT/dz - beta dS/dz) at the interfaces between centres, z positive up (s^-2)."""
    dz = np.diff(zc)
    return G * (EOS["tAlpha"] * np.diff(T) / dz - EOS["sBeta"] * np.diff(S) / dz)


# ----------------------------------------------------------------------------- qc
def test_qc_accepts_clean_profile():
    p = profiles.qc(summer_like())
    assert p is not None
    assert p["qc"]["inversion_fraction"] < 0.01


def test_qc_rejects_salinity_out_of_range():
    p = summer_like()
    p["S"][20] = 40.0
    assert profiles.qc(p) is None


def test_qc_rejects_temperature_out_of_range():
    p = summer_like()
    p["T"][5] = 12.0
    assert profiles.qc(p) is None


def test_qc_rejects_large_inversions():
    # salinity decreasing with depth over the bottom half: a 50 % inverted column
    z = -np.arange(1.0, 101.0)
    S = np.where(z > -50, 33.0 + 0.02 * (-z), 34.0 - 0.02 * (-z - 50))
    T = np.full_like(z, 1.0)
    assert profiles.qc(make_prof(z, T, S)) is None


def test_qc_tolerates_small_noise():
    rng = np.random.default_rng(0)
    p = summer_like()
    p["S"] = p["S"] + rng.normal(0, 5e-4, p["S"].size)   # CTD-like salinity noise
    assert profiles.qc(p) is not None


def test_qc_drops_nan_levels_and_rejects_too_short():
    p = summer_like()
    p["S"][::2] = np.nan
    q = profiles.qc(p)
    assert q is not None and np.all(np.isfinite(q["S"]))
    short = make_prof([-1.0, -2.0, -3.0], [1, 1, 1], [34, 34.1, 34.2])
    assert profiles.qc(short) is None


# ----------------------------------------------------------------------------- select
def test_select_polygon_months_years():
    poly = [(15.0, 77.7), (16.0, 77.7), (16.0, 77.9), (15.0, 77.9)]
    a = make_prof([-1, -2], [1, 1], [34, 34], lat=77.8, lon=15.5, time=dt.datetime(2010, 8, 1))
    b = make_prof([-1, -2], [1, 1], [34, 34], lat=77.8, lon=14.5, time=dt.datetime(2010, 8, 1))  # outside
    c = make_prof([-1, -2], [1, 1], [34, 34], lat=77.8, lon=15.5, time=dt.datetime(2010, 3, 1))  # month
    d = make_prof([-1, -2], [1, 1], [34, 34], lat=77.8, lon=15.5, time=dt.datetime(1985, 8, 1))  # year
    out = profiles.select([a, b, c, d], poly, {8, 9}, (1990, 2019))
    assert out == [a]
    assert len(profiles.select([a, b, c, d], poly, {8, 9}, None)) == 2


# ----------------------------------------------------------------------------- representative
def zgrid():
    f = vgrid.faces(CFG["delR"])
    return -0.5 * (f[1:] + f[:-1])


def test_representative_median_and_extension():
    zc = zgrid()
    ps = [summer_like(zmax=80.0, dS=d) for d in (-0.1, 0.0, 0.1, 0.2, -0.2)]
    info = {}
    T, S = profiles.representative(ps, zc, info=info)
    ref = summer_like(zmax=80.0)
    k = (zc < -2) & (zc > -75)
    assert np.allclose(S[k], np.interp(-zc[k], -ref["z"], ref["S"]), atol=1e-6)
    # below the last datum (80 m) the last value is repeated
    deep = zc < -81
    assert deep.any() and np.allclose(S[deep], S[~deep][-1]) and np.allclose(T[deep], T[~deep][-1])
    assert info["n_profiles"] == 5 and info["n_levels_extended_below"] == int(deep.sum())


def test_representative_stable_and_finite():
    """Review Focus 5: N^2 >= 0 at every interface and no NaN, even from unstable inputs."""
    zc = zgrid()
    rng = np.random.default_rng(1)
    ps = []
    for i in range(9):
        p = summer_like(zmax=rng.uniform(30, 110), dS=rng.normal(0, 0.1))
        # a 10 m salinity inversion at a random depth and a warm (destabilising under the
        # linear EOS) intrusion: the median keeps part of them
        k = slice(15 + i, 25 + i)
        p["S"][k] -= 0.3
        p["T"][k] += 0.5
        ps.append(p)
    info = {}
    T, S = profiles.representative(ps, zc, info=info)
    assert np.all(np.isfinite(T)) and np.all(np.isfinite(S))
    assert np.all(n2_linear(zc, T, S) >= -1e-12)
    assert info["n_levels_adjusted"] >= 1


def test_representative_stable_input_untouched():
    zc = zgrid()
    info = {}
    profiles.representative([summer_like()] * 5, zc, info=info)
    assert info["n_levels_adjusted"] == 0


def test_stabilize_is_minimal_l2():
    """Pool-adjacent-violators: one inverted pair is replaced by its thickness-weighted mean."""
    zc = np.array([-1.0, -2.0, -3.0, -4.0])
    T = np.zeros(4)
    S = np.array([34.0, 34.2, 34.1, 34.3])
    S2, n = profiles.stabilize(zc, T, S, EOS)
    assert n == 2
    assert np.allclose(S2, [34.0, 34.15, 34.15, 34.3])


def test_representative_needs_profiles():
    with pytest.raises(ValueError):
        profiles.representative([], zgrid())


# ----------------------------------------------------------------------------- load
def test_load_unis_hd_synthetic(tmp_path):
    netCDF4 = pytest.importorskip("netCDF4")
    fn = tmp_path / "UNIS_HD_20100815T120000_77.8_15.5.nc"
    with netCDF4.Dataset(fn, "w") as d:
        d.createDimension("TIME", 1)
        d.createDimension("PRES", 4)
        t = d.createVariable("TIME", "f8", ("TIME",))
        t.units = "days since 1950-01-01 00:00:00"
        t[:] = [(dt.datetime(2010, 8, 15, 12) - dt.datetime(1950, 1, 1)).total_seconds() / 86400]
        pr = d.createVariable("PRES", "f8", ("PRES",))
        pr[:] = [1.0, 2.0, 3.0, 4.0]
        for name, val in (("LATITUDE", 77.8), ("LONGITUDE", 15.5)):
            v = d.createVariable(name, "f8", ("TIME",))
            v[:] = [val]
        te = d.createVariable("TEMP", "f8", ("TIME", "PRES"), fill_value=np.nan)
        te[:] = [[2.0, 1.9, np.nan, 1.7]]
        sa = d.createVariable("PSAL", "f8", ("TIME", "PRES"), fill_value=np.nan)
        sa[:] = [[33.0, 33.1, 33.2, 33.3]]
    out = profiles.load_unis_hd(tmp_path)
    assert len(out) == 1
    p = out[0]
    assert p["lat"] == pytest.approx(77.8) and p["lon"] == pytest.approx(15.5)
    assert p["time"].year == 2010 and p["time"].month == 8
    assert p["z"].size == 3 and np.all(p["z"] < 0) and np.all(np.diff(p["z"]) < 0)
    assert p["z"][0] == pytest.approx(-0.99, abs=0.01)          # 1 dbar ~ 0.99 m at 78 N
    assert np.allclose(p["T"], [2.0, 1.9, 1.7])


# ----------------------------------------------------------------------------- theory numbers
def test_modes_column_constant_N():
    """Linear-EOS column with dT/dz constant: c_n = N H / (n pi) (Gill 1982 Sec. 6.11 [verificar])."""
    from exp4_arctic_bayiw.analysis import theory_numbers as tn
    delR = np.full(100, 1.0)
    H = 100.0
    zc = -(np.arange(100) + 0.5)
    N2 = 1e-4
    T = 5.0 + N2 / (G * EOS["tAlpha"]) * zc        # alpha g dT/dz = N2
    S = np.full(100, 34.0)
    z, n2 = tn.modes_column(delR, T, S, H, EOS)
    assert z[0] == 0.0 and z[-1] == pytest.approx(-H)
    assert np.allclose(n2, N2)
    c, _ = modes.vertical_modes(z, n2, H, n_modes=3)
    assert np.allclose(c, np.sqrt(N2) * H / (np.arange(1, 4) * np.pi), rtol=2e-3)


def test_modes_column_truncates_at_H():
    from exp4_arctic_bayiw.analysis import theory_numbers as tn
    delR = np.array(CFG["delR"])
    zc = zgrid()
    T = 1.0 + 4.0 * np.exp(zc / 10.0)
    S = 34.4 - 3.0 * np.exp(zc / 8.0)
    z, n2 = tn.modes_column(delR, T, S, 119.0, EOS)
    assert z[-1] == pytest.approx(-119.0) and np.all(np.diff(z) < 0) and np.all(n2 > 0)


def test_mask_geometry_rectangle():
    """Rectangle of 20 x 5 cells (dx = 200 m): coast contour, staircase edges and water path."""
    from exp4_arctic_bayiw.analysis import theory_numbers as tn
    dx = 200.0
    m = np.zeros((9, 26), bool)
    m[2:7, 3:23] = True
    mouth = np.zeros_like(m)
    mouth[2:7, 2] = True          # open "mouth" column on the west side, wet but not basin
    wet = m | mouth
    g = tn.mask_geometry(m, wet, mouth, dx)
    L, W = 20 * dx, 5 * dx
    assert g["staircase_coast_m"] == pytest.approx(2 * L + W)
    assert g["contour_total_m"] == pytest.approx(2 * (L + W), rel=0.12)
    assert g["contour_coast_m"] == pytest.approx(2 * L + W, rel=0.12)
    assert g["along_axis_length_m"] == pytest.approx(L, rel=0.06)
    assert g["axis_angle_deg"] == pytest.approx(0.0, abs=1.0)   # east-west axis


def test_bands_non_overlapping_and_resolved():
    from exp4_arctic_bayiw.analysis import theory_numbers as tn
    f_cph = modes.coriolis(78.0) / (2 * np.pi) * 3600
    b = tn.spectral_bands(f_cph, record_h=240.0, dt_out_s=1200.0, merian_lowest_cph=0.3)
    names = [x["name"] for x in b["bands"]]
    assert names == ["low", "kelvin_subinertial", "near_inertial", "superinertial", "seiche"]
    edges = [(x["lo_cph"], x["hi_cph"]) for x in b["bands"]]
    for (lo1, hi1), (lo2, hi2) in zip(edges, edges[1:]):
        assert hi1 == pytest.approx(lo2) and lo1 < hi1
    assert all(x["n_bins"] >= 2 for x in b["bands"][1:])
    assert edges[-1][1] == pytest.approx(1 / (2 * 1200 / 3600))     # Nyquist
    with pytest.raises(ValueError):
        tn.spectral_bands(f_cph, 240.0, 1200.0, merian_lowest_cph=0.1)   # seiche band would hit NI


# ----------------------------------------------------------------------------- real data
@pytest.mark.skipif(not UNIS.is_dir(), reason="UNIS HD profiles not downloaded")
def test_real_unis_summer_winter():
    pc = CFG["profiles"]
    allp = profiles.load_unis_hd(UNIS)
    assert len(allp) == 1547
    for season in ("summer", "winter"):
        s = pc["seasons"][season]
        sel = profiles.select(allp, pc["polygon_lonlat"], set(s["months"]), tuple(s["years"]))
        ok = [q for q in (profiles.qc(p, **pc["qc"]) for p in sel) if q is not None]
        assert len(ok) >= pc["min_profiles"], season
        info = {}
        T, S = profiles.representative(ok, zgrid(), min_count=pc["min_count_per_level"], info=info)
        assert np.all(np.isfinite(T)) and np.all(np.isfinite(S))
        assert np.all(n2_linear(zgrid(), T, S) >= -1e-12)


@pytest.mark.skipif(not (PKG / "figures" / "theory.json").is_file(), reason="theory.json not written")
def test_theory_json_contents():
    import json
    d = json.loads((PKG / "figures" / "theory.json").read_text())
    for k in ("vmf_summer", "vmf_winter", "thesis_bay_thesis_f", "thesis_bay_78N_vmf_summer"):
        c = d["cases"][k]
        assert len(c["c_n"]) >= 3 and c["c_n"][0] > c["c_n"][1] > c["c_n"][2] > 0
        assert c["inertial_period_h"] > 0 and c["kelvin_circuit_period_h"][0] > 0
    assert d["cases"]["vmf_summer"]["inertial_period_h"] == pytest.approx(12.23, abs=0.01)
    assert d["cases"]["thesis_bay_thesis_f"]["f_s-1"] == pytest.approx(6.97e-5, rel=1e-3)


def test_parse_namelist_fortran_numbers():
    """Regression: '2.E-4' (Fortran style, as in the thesis data.stage1) must parse as 2e-4."""
    from exp4_arctic_bayiw.analysis import theory_numbers as tn
    txt = " &PARM01\n tAlpha=2.E-4,\n sBeta=7.4E-4,\n delR=1.0, 1.1,\n\t2.0,\n f0=6.97E-5,\n &\n"
    assert tn.parse_namelist_list(txt, "tAlpha") == [2e-4]
    assert tn.parse_namelist_list(txt, "delR") == [1.0, 1.1, 2.0]
    assert tn.parse_namelist_list(txt, "f0") == [6.97e-5]
