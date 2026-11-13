"""Tests for exp4_arctic_bayiw.analysis.validate (Task 6: V2-V4 against the PD0) and
analysis.ice_check (F4 ice-open-boundary distance, notes/pd1_seaice.md section 8).

Every evaluator is tested on synthetic outputs on both sides of its frozen threshold."""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from exp4_arctic_bayiw.analysis import budget, ice_check, load, validate as va

PKG = Path(__file__).resolve().parents[1]
PD0 = PKG / "notes" / "pd0_preregistro.md"
F = 1.4266e-4                     # s^-1, f at 78 N (PD0 section 2.2)
DT = 1200.0                       # s, record interval
T_FIN = 280 * 3600.0              # s


# ---- thresholds against the PD0 text ------------------------------------------------------------

SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")


def _num(pattern: str, text: str, n: int = 1):
    m = re.search(pattern, text, re.S)
    assert m, f"PD0 text not found: {pattern}"
    vals = [float(m.group(k).translate(SUP)) for k in range(1, n + 1)]
    return vals[0] if n == 1 else tuple(vals)


def test_thresholds_are_the_pd0_text():
    """A threshold cannot change in validate.THRESHOLDS without this test failing (plan Task 6
    Step 2): every number is parsed from the frozen PD0 sections 3-5."""
    t = PD0.read_text()
    T = va.THRESHOLDS
    m = re.search(r"máx \|Δη\| ≤ (\d+)×10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+) m", t)
    assert m and float(m.group(1)) * 10 ** int(m.group(2).translate(SUP)) == T["V1"]["exact_max_abs_deta_m"]
    assert _num(r"la DEP en el pico difiere ≤ (\d+) %", t) / 100 == T["V1"]["spectral_psd_max_rel_diff"]
    assert _num(r"\*\*cercana a la inercial\*\* \| \[(\d\.\d+) f, (\d\.\d+) f\)", t, 2) == \
        T["V2a"]["near_inertial_band_f"]
    assert _num(r"está a ≤ (\d+)Δf de f", t) == T["V2a"]["max_peak_offset_df"]
    assert _num(r"\[(\d\.\d+) c₁\(H.{1,2}\), (\d\.\d+) c₁\(H_máx\)\]", t, 2) == \
        (T["V2b"]["c1_factor_mean_depth"], T["V2b"]["c1_factor_max_depth"])
    assert _num(r"Con la EOS ajustada: \[(\d\.\d+), (\d\.\d+)\] m s⁻¹", t, 2) == T["V2b"]["range_m_s"]
    assert _num(r"a ≥ (\d+) estaciones costeras sucesivas", t) == T["V2b"]["min_stations"]
    assert _num(r"escala de decaimiento hacia fuera de la costa de (\d+)–(\d+) km", t, 2) == \
        T["V3"]["decay_scale_km"]
    assert _num(r"dentro de un factor (\d+) de (\d\.\d+) m s⁻¹, es decir (\d\.\d+)–(\d\.\d+) m s⁻¹",
                t, 4) == (T["V3"]["c_factor"], T["V3"]["c_ref_m_s"]) + T["V3"]["c_range_m_s"]
    assert "**sentido ciclónico:** costa a la derecha, hacia el este por la costa sur" in t
    assert T["V3"]["sense"] == "cyclonic: eastward along the south coast"
    assert _num(r"\*\*V4a:\*\*[^\n]*difiere ≤ (\d+) %", t) / 100 == T["V4a"]["max_rel_diff"]
    assert _num(r"\*\*V4b:\*\*[^\n]*difiere ≤ (\d\.\d+) en valor absoluto", t) == T["V4b"]["max_abs_diff"]
    assert _num(r"\*\*factor mínimo es (\d\.\d+)\*\*", t) == T["factor_min"]
    assert _num(r"\*\*Apoyada\*\* \(como la tesis\) si φ_NI > (\d\.\d+)", t) == T["H4a"]["phi_ni_support"]
    assert _num(r"\*\*Alternativa\*\* \(domina Kelvin\) si φ_K > (\d\.\d+)", t) == T["H4a"]["phi_k_alternative"]
    assert _num(r"la banda baja lleva más del (\d+) %", t) / 100 == T["H4a"]["low_band_max_frac"]
    assert _num(r"Ē_in\(F2\)/Ē_in\(F3\) ≤ (\d\.\d+)", t) == T["H4b"]["energy_ratio_max"]
    assert _num(r"W_NI\(F1\)/W_NI\(F3\) ≥ (\d\.\d+)", t) == T["H4c"]["ratio_min"]
    assert _num(r"F1/F3, es ≥ (\d\.\d+)", t) == T["H4c"]["ratio_min"]
    assert _num(r"η\(KPP\)/η\(cierre de la tesis\) ≥ (\d\.\d+)", t) == T["H4d"]["ratio_min"]


def test_every_pd0_amendment_keeps_the_thresholds():
    """A dated amendment may change methods, never thresholds silently: each one must state that
    no threshold changes (Enmiendas 1 and 2 do); a new amendment that does not fails here."""
    t = PD0.read_text().split("## Enmiendas", 1)[1]
    parts = re.split(r"\n### Enmienda ", t)[1:]
    assert len(parts) >= 2
    for p in parts:
        assert re.search(r"Ningún umbral ni tolerancia cambia|No cambia ninguno", p), p[:40]


def test_v2b_range_of_theory_json_rounds_to_the_pd0_range():
    th = va.load_theory()
    lo, hi = va.v2b_range(th)
    assert (round(lo, 2), round(hi, 2)) == va.THRESHOLDS["V2b"]["range_m_s"]


# ---- V2a: inertial peak at the shelf station --------------------------------------------------

def _t(n=None):
    n = n or int(T_FIN / DT) + 1
    return np.arange(n) * DT


@pytest.mark.parametrize("freq_factor, ok", [(1.0, True), (1.0 + 0.4 / (F / (2 * np.pi) * T_FIN), True),
                                              (1.10, False)])
def test_v2a_peak_within_one_df_of_f(freq_factor, ok):
    t = _t()
    w = F * freq_factor
    decay = np.exp(-t / (100 * 3600))
    u, v = decay * np.cos(w * t), -decay * np.sin(w * t)
    r = va.v2a_inertial_peak(u, v, t, F)
    assert r["pass"] is ok
    assert r["offset_df"] == pytest.approx(abs(r["f_peak_hz"] - F / (2 * np.pi)) / r["df_hz"])


def test_v2a_only_looks_inside_the_near_inertial_band():
    t = _t()
    u = 5 * np.cos(2 * np.pi * t / (24 * 3600)) + np.cos(1.15 * F * t)     # big diurnal outside
    r = va.v2a_inertial_peak(u, 0 * u, t, F)
    assert 0.85 <= r["f_peak_hz"] / (F / 2 / np.pi) < 1.2 and r["pass"] is False


# ---- V2b: Kelvin front speed --------------------------------------------------------------------

def _front(c, s):
    t = _t()
    t0 = 10 * 3600.0
    z = np.stack([np.where(t >= t0 + si / c, 1.0 - np.exp(-(t - t0 - si / c).clip(0) / 3600), 0.0)
                  for si in s], axis=1)
    return t, z


@pytest.mark.parametrize("c, ok", [(0.45, True), (0.30, True), (0.9, False), (0.2, False)])
def test_v2b_arrival_speed_against_the_range(c, ok):
    s = np.array([0.0, 5e3, 10e3, 15e3])
    t, z = _front(c, s)
    r = va.v2b_kelvin_speed(z, t, s, (0.28, 0.62))
    assert r["c_m_s"] == pytest.approx(c, rel=0.03)
    assert r["pass"] is ok


def test_v2b_needs_three_stations_with_a_front():
    s = np.array([0.0, 5e3, 10e3])
    t, z = _front(0.45, s)
    z[:, 2] = 0.0                                         # no front at the third station
    r = va.v2b_kelvin_speed(z, t, s, (0.28, 0.62))
    assert r["pass"] is False and r["n_stations"] == 2


def test_v2b_westward_front_fails():
    s = np.array([0.0, 5e3, 10e3])
    t, z = _front(0.45, s[::-1])                          # arrives first at the east station
    r = va.v2b_kelvin_speed(z, t, s, (0.28, 0.62))
    assert r["c_m_s"] < 0 and r["pass"] is False


def test_pycnocline_displacement_sign():
    z = np.array([-1.0, -3.0, -5.0, -7.0])
    rho_b = np.array([1020.0, 1021.0, 1024.0, 1024.5])          # max N2 between -3 and -5
    k = va.pycnocline_level(rho_b, z, 1025.0)
    assert k in (1, 2)
    up = rho_b[k] + 0.3                                          # denser water: isopycnals rise
    zeta = va.displacement(np.array([rho_b[k], up]), rho_b, z, k)
    assert zeta[0] == 0 and zeta[1] > 0


# ---- V3: M2 internal Kelvin wave -----------------------------------------------------------------

W_M2 = 2 * np.pi / (12.4206012 * 3600)


def _kelvin(c=0.42, L=3e3, n_st=12, ds=2e3):
    t = _t()
    s = np.arange(n_st) * ds
    d = np.arange(0, 4200, 200.0)
    k = W_M2 / c
    coast = np.cos(W_M2 * t[:, None] - k * s[None]) * 0.5
    trans = 0.5 * np.exp(-d[None] / L) * np.cos(W_M2 * t[:, None] - 0.3)
    return t, s, coast, d, trans


@pytest.mark.parametrize("c, L, ok", [(0.42, 3e3, (True, True, True)),
                                      (-0.42, 3e3, (True, False, True)),
                                      (1.2, 3e3, (True, True, False)),
                                      (0.42, 0.5e3, (False, True, True)),
                                      (0.42, 20e3, (False, True, True))])
def test_v3_three_conditions(c, L, ok):
    t, s, coast, d, trans = _kelvin(c=c, L=L)
    r = va.v3_kelvin(coast, s, trans, d, t, W_M2, ramp_s=12.42 * 3600)
    assert (r["presence"]["pass"], r["sense"]["pass"], r["speed"]["pass"]) == ok
    assert r["pass"] is all(ok)
    assert r["speed"]["c_m_s"] == pytest.approx(c, rel=0.02)
    if ok[0]:
        assert r["presence"]["decay_scale_km"] == pytest.approx(L / 1e3, rel=0.05)


def test_harmonic_fit_recovers_amplitude_and_lag():
    t = _t()
    x = 2.0 * np.cos(W_M2 * t - 1.0) + 0.3
    a, ph = va.harmonic(x, t, W_M2)
    assert a == pytest.approx(2.0, rel=1e-6) and ph == pytest.approx(1.0, abs=1e-6)


# ---- V4: band partition of the flux and resolution test ---------------------------------------

BANDS_HZ = {"low": (0.0, 2 / T_FIN), "kelvin_subinertial": (2 / T_FIN, 0.85 * F / 2 / np.pi),
            "near_inertial": (0.85 * F / 2 / np.pi, 1.2 * F / 2 / np.pi),
            "superinertial": (1.2 * F / 2 / np.pi, 0.1775 / 3600), "seiche": (0.1775 / 3600, 1 / (2 * DT))}


def test_band_cospectrum_partition_is_exact_and_separates_frequencies():
    t = _t(840)
    nz, nf = 2, 3
    w_ni, w_se = F, 2 * np.pi / (3 * 3600)
    up = (np.cos(w_ni * t) * 2 + np.cos(w_se * t))[:, None, None] * np.ones((1, nz, nf))
    pp = (np.cos(w_ni * t) * 3 + np.cos(w_se * t) * 1)[:, None, None] * np.ones((1, nz, nf))
    dz, wdt = np.full((nz, nf), 5.0), np.array([100.0, 200.0, -100.0])
    out = va.band_flux(up, pp, dz, wdt, DT, BANDS_HZ)
    total = float(((up * pp).mean(0) * dz * wdt).sum())
    assert sum(out.values()) == pytest.approx(total, rel=1e-10)
    assert out["near_inertial"] == pytest.approx(0.5 * 6 * 10 * 200, rel=0.05)
    assert out["seiche"] == pytest.approx(0.5 * 1 * 10 * 200, rel=0.05)


def _v4(total, fr):
    return {"F_bc_out_J": total, "band_out_J": {b: total * f for b, f in fr.items()}}


FR = {"kelvin_subinertial": 0.5, "near_inertial": 0.3, "superinertial": 0.15, "seiche": 0.05}


@pytest.mark.parametrize("ratio, shift, ok", [(1.2, 0.0, (True, True)), (0.76, 0.09, (True, True)),
                                              (1.3, 0.0, (False, True)), (1.0, 0.11, (True, False))])
def test_v4_resolution(ratio, shift, ok):
    fr2 = dict(FR, kelvin_subinertial=FR["kelvin_subinertial"] - shift,
               near_inertial=FR["near_inertial"] + shift)
    r = va.v4_resolution(_v4(1e12, FR), _v4(ratio * 1e12, fr2))
    assert (r["V4a"]["pass"], r["V4b"]["pass"]) == ok
    assert r["pass"] is all(ok)


def _channel_ds(nt=60):
    """A channel x in [0, 9) with a closed 'fjord' (columns 5..8) and an open face at i = 5;
    a mode-1 baroclinic wave through it (u' and p' in phase, energy into the box)."""
    nz, ny, nx, d = 3, 5, 10, 200.0
    drF = np.array([5.0, 10.0, 15.0])
    hC = np.zeros((nz, ny, nx)); hC[:, 1:4, 1:9] = 1.0
    wet = hC[0] > 0
    hW = np.zeros_like(hC); hS = np.zeros_like(hC)
    hW[:, :, 1:] = (wet[:, 1:] & wet[:, :-1])[None]
    hS[:, 1:, :] = (wet[1:] & wet[:-1])[None]
    t = np.arange(nt) * DT
    w = 2 * np.pi / (6 * 3600)
    mode = np.array([1.0, 0.2, -0.6])
    U = 0.05 * np.cos(w * t)[:, None, None, None] * mode[None, :, None, None] * hW[None]
    V = np.zeros_like(U)
    rho_b = np.array([1025.0, 1025.5, 1026.0])
    rho = rho_b[None, :, None, None] + 0.02 * np.cos(w * t)[:, None, None, None] \
        * np.array([-1.0, 0.0, 1.0])[None, :, None, None] * hC[None]
    eta = 0.01 * np.cos(w * t)[:, None, None] * wet[None]
    one = np.ones((ny, nx))
    ds = xr.Dataset(
        {"U": (("time", "k", "j", "i"), U), "V": (("time", "k", "j", "i"), V),
         "RHO": (("time", "k", "j", "i"), rho), "ETA": (("time", "j", "i"), eta),
         "hFacC": (("k", "j", "i"), hC), "hFacW": (("k", "j", "i"), hW),
         "hFacS": (("k", "j", "i"), hS), "drF": (("k",), drF),
         "RC": (("k",), -(np.cumsum(drF) - 0.5 * drF)),
         "rA": (("j", "i"), d * d * one), "rAs": (("j", "i"), d * d * one),
         "rAw": (("j", "i"), d * d * one), "dxG": (("j", "i"), d * one),
         "dyG": (("j", "i"), d * one)},
        coords={"time": t}, attrs={"rho0": 1026.0, "g": 9.81})
    mask = np.zeros((ny, nx), bool); mask[1:4, 5:9] = True
    groups = {"mouth": {"x": np.zeros((ny, nx), bool), "y": np.zeros((ny, nx), bool)}}
    groups["mouth"]["x"][1:4, 5] = True
    return ds, mask, groups, rho_b


def test_face_profiles_reproduce_the_kang_fringer_baroclinic_flux():
    ds, mask, groups, rho_b = _channel_ds()
    r = budget.kf_budget(ds, mask, rho_b=rho_b, face_groups=groups)
    fp = va.face_profiles(ds, mask, groups, rho_b=rho_b)
    Fbc = (fp["up"] * fp["pp"] * fp["dz"][None]).sum(1) @ fp["w"]
    np.testing.assert_allclose(Fbc, r["F_bc"], rtol=1e-10, atol=1e-6)
    out = va.flux_summary(fp)
    np.testing.assert_allclose(out["F_bc_out_J"], budget.f_bc_pressure_out(r)[-1], rtol=1e-10)
    assert out["per_section_out_J"]["mouth"] == pytest.approx(out["F_bc_out_J"])


# ---- stations ------------------------------------------------------------------------------------

def _fjord_grid():
    """East-west fjord (rows 10..19, columns 10..69) opening west onto a 'shelf' (columns 0..9)."""
    ny, nx, dx = 30, 80, 200.0
    wet = np.zeros((ny, nx), bool)
    wet[:, :10] = True
    wet[10:20, 10:70] = True
    inner = np.zeros_like(wet); inner[10:20, 20:70] = True
    xc = np.arange(nx) * dx; yc = np.arange(ny) * dx
    return wet, inner, xc, yc, dx


def test_south_coast_stations_are_ordered_eastward_and_on_the_south_coast():
    wet, inner, xc, yc, dx = _fjord_grid()
    deep = wet.copy()
    st = va.coastal_stations(wet, deep, inner, xc, yc, axis_deg=0.0, spacing_m=2000.0, side="south")
    js, is_ = zip(*st["cells"])
    assert all(j == 10 for j in js)                       # first wet row above the south wall
    assert list(is_) == sorted(is_) and len(is_) >= 5
    assert np.all(np.diff(st["s_m"]) > 0)
    np.testing.assert_allclose(np.diff(st["s_m"]), np.diff(np.array(is_)) * dx)


def test_shelf_station_is_away_from_the_open_boundaries_and_the_fjord():
    wet, inner, xc, yc, dx = _fjord_grid()
    wet2 = np.zeros((60, 80), bool); wet2[:, :40] = True; wet2[25:35, 40:80] = True
    inner2 = np.zeros_like(wet2); inner2[25:35, 50:80] = True
    j, i = va.shelf_station(wet2, inner2, ("W", "N", "S"), n_sponge=5)
    assert i >= 10 and 10 <= j <= 49 and not inner2[j, i]


# ---- F4: ice never within 10 cells of an open boundary ----------------------------------------

def _meta(path: Path, shape, it, flds, nrec):
    dims = ",\n".join(f"{n:5d}, 1, {n:5d}" for n in reversed(shape))
    fl = " ".join(f"'{f:<8s}'" for f in flds)
    path.write_text(f" nDims = [ {len(shape)} ];\n dimList = [\n{dims}\n ];\n"
                    f" dataprec = [ 'float64' ];\n nrecords = [ {nrec} ];\n"
                    f" timeStepNumber = [ {it} ];\n nFlds = [ {len(flds)} ];\n"
                    f" fldList = {{\n {fl}\n }};\n")


def _ice_run(root: Path, i_min_stage2: int, mode="free_drift"):
    """Two-stage ice run, 40 x 40 cells, ice west edge at i = 20 in stage 1 and i_min_stage2 in
    stage 2 (the W column is an open boundary)."""
    ny, nx = 40, 40
    (root / "input").mkdir(parents=True)
    (root / "input" / "case.json").write_text(json.dumps(
        {"params": {"case": "F4", "nx": nx, "ny": ny, "n_sponge": 10, "ob_sides": ["W", "N", "S"]},
         "actual": {"ice": {"mode": mode}}}))
    np.full((ny, nx), -100.0).astype(">f8").tofile(root / "bathy.bin")
    for d, its, imin in ((root, (0, 10), 20), (root / "stage2", (20, 30), i_min_stage2)):
        d.mkdir(exist_ok=True)
        (d / "run_status.json").write_text(json.dumps({"state": "OK", "last_iter": its[-1] + 10,
                                                         "n_timesteps": 20}))
        for it in its:
            a = np.zeros((4, ny, nx)); a[0, 18:22, imin:30] = 1.0
            a.astype(">f8").tofile(d / f"seaice.{it:010d}.data")
            _meta(d / f"seaice.{it:010d}.meta", (ny, nx), it, ["SIarea", "SIheff", "SIuice", "SIvice"], 4)
    (root / "run_status_final.json").write_text(json.dumps({"state": "OK", "last_iter": 40,
                                                              "expected_last_iter": 40}))
    return root


def test_ice_check_reads_every_stage_and_fails_near_the_boundary(tmp_path):
    ok = ice_check.ice_ob_distance(_ice_run(tmp_path / "a", 15))
    assert ok["passes"] and ok["min_cells"] == 15 and ok["records"] == 4
    bad = ice_check.ice_ob_distance(_ice_run(tmp_path / "b", 9))       # only in stage 2
    assert not bad["passes"] and bad["min_cells"] == 9 and bad["at_iter"] == 20


def test_ice_check_cli_exit_code_and_stamp(tmp_path):
    stamps = tmp_path / "ice_ob_check.json"
    good, badr = _ice_run(tmp_path / "a", 15), _ice_run(tmp_path / "b", 9)
    assert va.main(["ice-check", str(good), "--stamps", str(stamps)]) == 0
    assert va.main(["ice-check", str(badr), "--stamps", str(stamps)]) == 1
    rec = json.loads(stamps.read_text())
    assert rec[str(good.resolve())]["passes"] and not rec[str(badr.resolve())]["passes"]
    ice_check.require_ice_ob_check(good, stamps)
    with pytest.raises(RuntimeError, match="within"):
        ice_check.require_ice_ob_check(badr, stamps)


def test_ice_run_is_refused_without_a_current_passing_check(tmp_path):
    stamps = tmp_path / "ice_ob_check.json"
    good = _ice_run(tmp_path / "a", 15)
    with pytest.raises(RuntimeError, match="ice-check"):
        ice_check.require_ice_ob_check(good, stamps)
    ice_check.check_and_stamp(good, stamps)
    ice_check.require_ice_ob_check(good, stamps)
    # the run changed after the check (status rewritten): the stamp is stale
    (good / "run_status_final.json").write_text(json.dumps({"state": "OK", "last_iter": 40,
                                                              "expected_last_iter": 40, "x": 1}))
    with pytest.raises(RuntimeError, match="stale"):
        ice_check.require_ice_ob_check(good, stamps)


def test_load_run_enforces_the_ice_check(tmp_path, monkeypatch):
    run = _ice_run(tmp_path / "a", 15)
    monkeypatch.setattr(ice_check, "STAMPS", tmp_path / "none.json")
    with pytest.raises(RuntimeError, match="ice-check"):
        load.run(run)


def test_ice_stamp_goes_stale_when_the_seaice_files_change(tmp_path):
    stamps = tmp_path / "s.json"
    good = _ice_run(tmp_path / "a", 15)
    ice_check.check_and_stamp(good, stamps)
    f = good / "stage2" / "seaice.0000000030.data"
    f.write_bytes(f.read_bytes())                               # same bytes, new mtime
    import os
    os.utime(f, (1, 1))
    with pytest.raises(RuntimeError, match="stale"):
        ice_check.require_ice_ob_check(good, stamps)


def test_ice_check_fails_if_the_records_do_not_cover_the_run(tmp_path):
    r = _ice_run(tmp_path / "a", 15)
    for x in (r / "stage2").glob("seaice.0000000030.*"):       # last record 20 < 40 - 10
        x.unlink()
    res = ice_check.ice_ob_distance(r)
    assert not res["passes"] and "cover" in res["reason"]
    r2 = _ice_run(tmp_path / "b", 15)
    for x in (r2 / "stage2").glob("seaice.*"):                  # a stage without records
        x.unlink()
    res2 = ice_check.ice_ob_distance(r2)
    assert not res2["passes"] and "stage2" in res2["reason"]


def test_stage_dirs_ignore_renamed_stage_directories(tmp_path):
    for n in ("stage2", "stage2_old", "stage10", "stagex"):
        (tmp_path / n).mkdir()
    assert [p.name for p in ice_check.stage_dirs(tmp_path)][1:] == ["stage2", "stage10"]
    assert [p.name for p in load.stage_dirs(tmp_path)][1:] == ["stage2", "stage10"]


def test_ice_area_refuses_other_precisions(tmp_path):
    r = _ice_run(tmp_path / "a", 15)
    m = r / "seaice.0000000000.meta"
    m.write_text(m.read_text().replace("float64", "float32"))
    with pytest.raises(ValueError, match="float64"):
        ice_check.ice_ob_distance(r)


def test_open_water_run_needs_no_ice_check(tmp_path):
    assert ice_check.is_ice_run(tmp_path) is False


# ---- CLI on real run directories refuses runs that are not OK ----------------------------------

def test_evaluate_cli_refuses_a_failed_run(tmp_path):
    bad = tmp_path / "F1"
    bad.mkdir()
    (bad / "run_status.json").write_text(json.dumps({"state": "FAILED", "last_iter": 3, "n_timesteps": 10}))
    with pytest.raises(RuntimeError, match="FAILED"):
        va.main(["evaluate", "--f1", str(bad), "--out", str(tmp_path / "v.json")])


def test_evaluate_writes_validation_json(tmp_path, monkeypatch):
    """Wiring of the CLI: each criterion gets the runs it needs, results go to the JSON with
    the thresholds and a pending classification."""
    calls = []
    monkeypatch.setattr(va, "open_run", lambda p: calls.append(Path(p).name) or Path(p).name)
    monkeypatch.setattr(va, "v2_from_run", lambda ds, th: {"V2a": {"pass": True}, "V2b": {"pass": False}})
    monkeypatch.setattr(va, "v3_from_run", lambda ds: {"pass": True})
    monkeypatch.setattr(va, "v4_from_runs", lambda a, b: {"pass": True, "V4a": {}, "V4b": {}})
    out = tmp_path / "v.json"
    assert va.main(["evaluate", "--f1", str(tmp_path / "F1"), "--v3", str(tmp_path / "V3"),
                    "--v4", str(tmp_path / "V4"), "--out", str(out)]) == 0
    r = json.loads(out.read_text())
    assert calls == ["F1", "V3", "V4"]
    assert r["V2b"]["pass"] is False and r["V3"]["pass"] is True and r["V4"]["pass"] is True
    assert r["thresholds"] == json.loads(json.dumps(va.THRESHOLDS))
    assert r["summary"]["failed"] == ["V2b"] and r["summary"]["go_to_PD2"] is True
    assert r["V2b"]["classification"] == "pending"


# ---- Figure 1 on synthetic data -----------------------------------------------------------------

def test_fig1_renders_from_a_synthetic_validation_json(tmp_path):
    from exp4_arctic_bayiw.analysis import fig1_validation as f1
    t, s, coast, d, trans = _kelvin()
    v3 = va.v3_kelvin(coast, s, trans, d, t, W_M2, ramp_s=12.42 * 3600)
    tt = _t()
    u = np.exp(-tt / 3.6e5) * np.cos(F * tt)
    v2a = va.v2a_inertial_peak(u, -np.exp(-tt / 3.6e5) * np.sin(F * tt), tt, F)
    tf, zf = _front(0.45, np.array([0.0, 5e3, 10e3, 15e3]))
    v2b = va.v2b_kelvin_speed(zf, tf, np.array([0.0, 5e3, 10e3, 15e3]), (0.28, 0.62))
    v4 = va.v4_resolution(_v4(1e12, FR), _v4(1.1e12, FR))
    val = {"V2a": v2a, "V2b": v2b, "V3": v3, "V4": v4, "thresholds": va.THRESHOLDS}
    p = tmp_path / "validation.json"
    p.write_text(json.dumps(val, default=va.jsonable))
    files = f1.main(["--validation", str(p), "--out", str(tmp_path / "fig1_validation")])
    assert files == 0
    assert (tmp_path / "fig1_validation.png").stat().st_size > 10_000
    assert (tmp_path / "fig1_validation.pdf").exists()


# ---- the dataset paths (wiring of v2/v3/v4_from_run on a synthetic fjord) -----------------------

def _synthetic_fjord_ds(nt=240, kind="tide"):
    """Shelf (columns 0..29) and an east-west fjord (rows 25..36, columns 30..159) 60 m deep; M2
    Kelvin wave eastward along the south coast (kind 'tide') or an inertial shelf current and a
    Kelvin front (kind 'wind')."""
    ny, nx, nz, dx = 60, 160, 4, 200.0
    wet = np.zeros((ny, nx), bool); wet[:, :30] = True; wet[25:37, 30:160] = True
    inner = np.zeros_like(wet); inner[25:37, 40:159] = True
    drF = np.full(nz, 15.0)
    RC = -(np.cumsum(drF) - 0.5 * drF)
    hC = np.repeat(wet[None].astype(float), nz, 0)
    rho_b = np.array([1024.0, 1025.5, 1025.7, 1025.8])                     # max N2 at k = 0-1
    k = va.pycnocline_level(rho_b, RC, 1026.0)
    drdz = np.gradient(rho_b, RC)[k]
    t = np.arange(nt) * DT
    X = np.arange(nx) * dx; Y = np.arange(ny) * dx
    XX, YY = np.meshgrid(X, Y)
    d = (YY - 24.5 * dx).clip(0)                                  # distance from the south wall
    if kind == "tide":
        zeta = 0.5 * np.exp(-d / 3e3)[None] * np.cos(W_M2 * t[:, None, None] - W_M2 / 0.42 * XX[None])
    else:
        c = 0.45
        arr = 10 * 3600 + (XX - 30 * dx).clip(0) / c
        zeta = np.where(t[:, None, None] >= arr[None], 1.0, 0.0) * (d < 1000)[None]
    zeta = zeta * inner[None] * (t > 0)[:, None, None]
    RHO = rho_b[None, :, None, None] - zeta[:, None] * drdz * hC[None]
    U = np.zeros((nt, nz, ny, nx)); V = np.zeros_like(U)
    if kind == "wind":
        U[:, 0] = np.cos(F * t)[:, None, None] * wet
        V[:, 0] = -np.sin(F * t)[:, None, None] * wet
    one = np.ones((ny, nx))
    ds = xr.Dataset(
        {"U": (("time", "k", "j", "i"), U), "V": (("time", "k", "j", "i"), V),
         "RHO": (("time", "k", "j", "i"), RHO), "ETA": (("time", "j", "i"), np.zeros((nt, ny, nx))),
         "hFacC": (("k", "j", "i"), hC), "hFacW": (("k", "j", "i"), hC), "hFacS": (("k", "j", "i"), hC),
         "drF": (("k",), drF), "RC": (("k",), RC), "Depth": (("j", "i"), 60.0 * wet),
         "XC": (("j", "i"), XX), "YC": (("j", "i"), YY),
         "rA": (("j", "i"), dx * dx * one), "rAs": (("j", "i"), dx * dx * one),
         "rAw": (("j", "i"), dx * dx * one), "dxG": (("j", "i"), dx * one), "dyG": (("j", "i"), dx * one)},
        coords={"time": t}, attrs={"rho0": 1026.0, "g": 9.81, "run_dir": "synthetic"})
    groups = {"mouth": {"x": np.zeros((ny, nx), bool), "y": np.zeros((ny, nx), bool)}}
    groups["mouth"]["x"][25:37, 40] = True
    return ds, inner, groups


@pytest.fixture
def synth(monkeypatch):
    params = {"ob_sides": ["W", "N", "S"], "n_sponge": 5, "tide": {"period_s": 12.4206012 * 3600}}
    actual = {"tide": {"period_s": 12.4206012 * 3600, "ramp_s": 2 * 12.4206012 * 3600}}
    monkeypatch.setattr(va, "_case_params", lambda ds: params)
    monkeypatch.setattr(va, "_case_actual", lambda ds: actual)

    def make(kind):
        ds, inner, groups = _synthetic_fjord_ds(kind=kind)
        monkeypatch.setattr(budget, "fjord_box", lambda d, **kw: (inner, groups))
        return ds
    return make


def test_v3_from_run_on_a_synthetic_kelvin_wave(synth):
    r = va.v3_from_run(synth("tide"))
    assert r["sense"]["pass"] and r["speed"]["c_m_s"] == pytest.approx(0.42, rel=0.05)
    assert r["presence"]["decay_scale_km"] == pytest.approx(3.0, rel=0.1) and r["pass"]
    assert len(r["stations"]["cells"]) >= 5
    assert r["ramp_s"] == 2 * 12.4206012 * 3600                    # actual.tide.ramp_s (review 6)


def test_v2_from_run_on_a_synthetic_wind_response(synth):
    th = va.load_theory()
    r = va.v2_from_run(synth("wind"), th)
    assert r["V2a"]["pass"] and r["V2a"]["station_ji"][1] >= 10
    assert r["V2b"]["c_m_s"] == pytest.approx(0.45, rel=0.15) and r["V2b"]["n_stations"] >= 3


def test_flux_from_run_partitions_the_section_flux(synth):
    """Wind synthetic: baroclinic surface current and a density front at the mouth, so the flux
    is not zero; the band sum equals the record integral up to the trapezoid end terms."""
    ds = synth("wind")
    th = va.load_theory()
    mask, groups = budget.fjord_box(ds)
    fp = va.face_profiles(ds, mask, groups)
    Fq = ((fp["up"] * fp["pp"] * fp["dz"][None]).sum(1) * fp["w"][None]).sum(1)      # W, + in
    out = va.flux_from_run(ds, th)
    assert abs(out["F_bc_out_J"]) > 1.0                           # J: not a vacuous zero
    end = 0.5 * DT * (abs(Fq[0]) + abs(Fq[-1]))
    assert abs(sum(out["band_out_J"].values()) - out["F_bc_out_J"]) <= end + 1e-9 * abs(out["F_bc_out_J"])
    assert sum(out["band_out_J"].values()) == pytest.approx(-Fq.sum() * DT, rel=1e-10)
    assert set(out["per_section_out_J"]) == {"mouth"}


def test_flux_summary_requires_the_record_to_start_at_zero(synth):
    ds = synth("wind")
    mask, groups = budget.fjord_box(ds)
    fp = va.face_profiles(ds, mask, groups, time=slice(1, None))
    with pytest.raises(ValueError, match="t\\[0\\]"):
        va.flux_summary(fp)


def test_go_to_pd2_is_undecided_while_criteria_are_not_evaluated(tmp_path, monkeypatch):
    monkeypatch.setattr(va, "open_run", lambda p: Path(p).name)
    monkeypatch.setattr(va, "v2_from_run", lambda ds, th: {"V2a": {"pass": True}, "V2b": {"pass": True}})
    r = va.evaluate(f1=tmp_path / "F1")
    assert r["summary"]["not_evaluated"] == ["V3", "V4"] and r["summary"]["go_to_PD2"] is None


def test_v2b_flags_gaps_in_the_station_sequence():
    wet, inner, xc, yc, dx = _fjord_grid()
    deep = wet.copy()
    deep[:, 40:46] = False                                     # no deep water for one target
    st = va.coastal_stations(wet, deep, inner, xc, yc, axis_deg=0.0, spacing_m=1000.0, side="south")
    assert st["gaps"] and st["successive"] is False
    st2 = va.coastal_stations(wet, wet, inner, xc, yc, axis_deg=0.0, spacing_m=1000.0, side="south")
    assert st2["gaps"] == [] and st2["successive"] is True


def test_face_profiles_with_dry_bottom_levels_still_match_kang_fringer():
    """face_profiles reads only the wet levels of each strait rectangle (I/O); a dry bottom level
    must not change the flux."""
    ds, mask, groups, rho_b = _channel_ds()
    for v in ("hFacC", "hFacW", "hFacS"):
        ds[v].values[2] = 0.0
    ds["U"].values[:, 2] = 0.0
    r = budget.kf_budget(ds, mask, rho_b=rho_b, face_groups=groups)
    fp = va.face_profiles(ds, mask, groups, rho_b=rho_b)
    Fbc = (fp["up"] * fp["pp"] * fp["dz"][None]).sum(1) @ fp["w"]
    np.testing.assert_allclose(Fbc, r["F_bc"], rtol=1e-10, atol=1e-6)
    assert np.abs(Fbc).max() > 0


def test_coarse_box_maps_the_fine_box_and_its_sound_faces():
    """V4 (400 m): the inner-fjord box is the 200 m box mapped onto the coarse grid (cells with at
    least half of their 2 x 2 children in the fine box); every open face is assigned to the
    nearest sound, and a face far from every sound is an error."""
    nyf, nxf, dxf = 20, 40, 200.0
    wf = np.zeros((nyf, nxf), bool); wf[6:14, :38] = True               # channel closed at the east
    mf = np.zeros_like(wf); mf[6:14, 20:38] = True                      # fine box east of x = 20
    xc = (np.arange(nxf // 2) + 0.5) * 2 * dxf; yc = (np.arange(nyf // 2) + 0.5) * 2 * dxf
    wc = wf.reshape(nyf // 2, 2, nxf // 2, 2).any((1, 3))
    secs = {"S": ([4000.0, 1000.0], [4000.0, 3000.0])}
    m, g = va.coarse_box(mf, wc, xc, yc, secs)
    assert m.sum() == 4 * 9 and g["S"]["x"].sum() == 4 and g["S"]["y"].sum() == 0
    with pytest.raises(RuntimeError, match="far from"):
        va.coarse_box(mf, wc, xc, yc, {"S": ([0.0, 0.0], [0.0, 100.0])})


def test_evaluate_resumes_from_a_partial_json(tmp_path, monkeypatch):
    out = tmp_path / "v.json"
    out.write_text(json.dumps({"V2a": {"pass": True}, "V2b": {"pass": True}, "V3": {"pass": True}}))
    calls = []
    monkeypatch.setattr(va, "open_run", lambda p: Path(p).name)
    monkeypatch.setattr(va, "v2_from_run", lambda ds, th: calls.append("v2") or {})
    monkeypatch.setattr(va, "v3_from_run", lambda ds: calls.append("v3") or {})
    monkeypatch.setattr(va, "v4_from_runs", lambda a, b: {"pass": False, "V4a": {}, "V4b": {}})
    assert va.main(["evaluate", "--f1", "F1", "--v3", "V3", "--v4", "V4", "--out", str(out), "--resume"]) == 0
    r = json.loads(out.read_text())
    assert calls == [] and r["summary"]["failed"] == ["V4"] and r["summary"]["go_to_PD2"] is True
