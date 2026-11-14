"""Tests for exp4_arctic_bayiw.analysis.partition (Task 7, Step 1): band partition of the strait
fluxes, box energy without the low band, wind work in a band, station band energy, trapped
fraction and the run bands.

Every kernel is checked against budget.kf_budget (Kang & Fringer 2012, the thesis budget) on a
synthetic box whose fields hold known frequencies, so the band split and the time integrals can
be compared exactly."""
from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from exp4_arctic_bayiw.analysis import budget, partition as pa, validate as va

DT = 1200.0                          # s
NT = 120                             # 40 h
T_REC = NT * DT
H = 3600.0
# test bands [Hz]: low [0, 2 df), mid [2 df, 0.1 cph), high [0.1 cph, Nyquist]
BANDS = {"low": (0.0, 2.0 / T_REC), "mid": (2.0 / T_REC, 0.1 / H), "high": (0.1 / H, 0.5 / DT)}
F_LOW, F_MID, F_HIGH = 1.0 / T_REC, 0.05 / H, 0.25 / H           # bins 1, 2 and 10


def _synthetic(seed=0, nt=NT, wind=True):
    """Basin 4 x 10 x 12 (200 m cells), a dry corner and partial bottom cells; the box is rows
    3..7, columns 4..9, open on all four sides (x faces at i = 4, 10; y faces at j = 3, 8).
    U, V and rho' are sums of a constant, a low-band, a mid-band and a high-band cosine with
    random amplitudes and phases per point."""
    rng = np.random.default_rng(seed)
    nz, ny, nx, d = 4, 10, 12, 200.0
    drF = np.array([5.0, 10.0, 15.0, 20.0])
    hC = np.ones((nz, ny, nx))
    hC[:, 0:2, 0:3] = 0.0                                      # land corner
    hC[3, :, 6:] = 0.5                                         # partial bottom cells
    hC[3, 5, 7] = 0.0                                          # a dry bottom cell in the box
    wet3 = hC > 0
    hW = np.zeros_like(hC); hS = np.zeros_like(hC)
    hW[:, :, 1:] = np.minimum(hC[:, :, 1:], hC[:, :, :-1]) * (wet3[:, :, 1:] & wet3[:, :, :-1])
    hS[:, 1:, :] = np.minimum(hC[:, 1:, :], hC[:, :-1, :]) * (wet3[:, 1:, :] & wet3[:, :-1, :])
    t = np.arange(nt) * DT

    def field(scale, mask3):
        out = np.zeros((nt, nz, ny, nx))
        for f in (0.0, F_LOW, F_MID, F_HIGH):
            a = scale * rng.standard_normal((nz, ny, nx))
            ph = rng.uniform(0, 2 * np.pi, (nz, ny, nx))
            out += a[None] * np.cos(2 * np.pi * f * t[:, None, None, None] + ph[None])
        return out * mask3[None]

    U = field(0.05, hW > 0)
    V = field(0.05, hS > 0)
    rho_b = np.array([1025.0, 1025.4, 1025.9, 1026.5])
    RHO = rho_b[None, :, None, None] + field(0.02, wet3)
    ETA = 0.01 * rng.standard_normal((nt, ny, nx)) * wet3[0][None]
    one = np.ones((ny, nx))
    X, Y = np.meshgrid(np.arange(nx) * d, np.arange(ny) * d)
    ds = xr.Dataset(
        {"U": (("time", "k", "j", "i"), U), "V": (("time", "k", "j", "i"), V),
         "RHO": (("time", "k", "j", "i"), RHO), "ETA": (("time", "j", "i"), ETA),
         "hFacC": (("k", "j", "i"), hC), "hFacW": (("k", "j", "i"), hW),
         "hFacS": (("k", "j", "i"), hS), "drF": (("k",), drF),
         "RC": (("k",), -(np.cumsum(drF) - 0.5 * drF)),
         "rA": (("j", "i"), d * d * one), "rAs": (("j", "i"), d * d * one),
         "rAw": (("j", "i"), d * d * one), "dxG": (("j", "i"), d * one),
         "dyG": (("j", "i"), d * one), "XC": (("j", "i"), X), "YC": (("j", "i"), Y)},
        coords={"time": t}, attrs={"rho0": 1026.0, "g": 9.81})
    mask = np.zeros((ny, nx), bool)
    mask[3:8, 4:10] = True
    f = budget.box_faces(mask)
    groups = {n: {"x": np.zeros((ny, nx), bool), "y": np.zeros((ny, nx), bool)}
              for n in ("xfaces", "yfaces")}
    for j, i, _ in f["x"]:
        groups["xfaces"]["x"][j, i] = True
    for j, i, _ in f["y"]:
        groups["yfaces"]["y"][j, i] = True
    tx0 = 0.1 * rng.standard_normal((ny, nx))
    ty0 = 0.1 * rng.standard_normal((ny, nx))
    tau = (lambda tt: (tx0 * np.exp(-((tt - 10 * H) / (3 * H)) ** 2),
                       ty0 * np.exp(-((tt - 10 * H) / (3 * H)) ** 2))) if wind else None
    return ds, mask, groups, rho_b, tau


@pytest.fixture(scope="module")
def syn():
    ds, mask, groups, rho_b, tau = _synthetic()
    r = budget.kf_budget(ds, mask, tau=tau, rho_b=rho_b, face_groups=groups)
    return ds, mask, groups, rho_b, tau, r


# ---- strait fluxes ------------------------------------------------------------------------------

def test_face_series_reproduce_kf_budget_per_section(syn):
    """F_bc, A_bc, F_bt, A_bt of each section from the face series equal kf_budget's (+ in)."""
    ds, mask, groups, rho_b, tau, r = syn
    fs = pa.face_series(ds, mask, groups, rho_b=rho_b)
    ser = pa.section_series(fs)
    for g in groups:
        for k in ("F_bc", "A_bc", "F_bt", "A_bt"):
            np.testing.assert_allclose(ser[g][k], r[f"{k}_{g}"], rtol=1e-10,
                                       atol=1e-9 * np.abs(r[f"{k}_{g}"]).max())
        assert np.abs(ser[g]["A_bc"]).max() > 0 and np.abs(ser[g]["F_bc"]).max() > 0


def test_band_flux_of_pressure_work_and_advection_adds_up_to_the_dt_sum(syn):
    """The bands of F_bc, A_bc and F_bc + A_bc (+ out) add up to -dt sum over the records
    (Parseval, ruling: band sums reproduce the time integral, dt-sum equality as in Task 6)."""
    ds, mask, groups, rho_b, tau, r = syn
    fs = pa.face_series(ds, mask, groups, rho_b=rho_b)
    out = pa.band_flux_series(fs, BANDS)
    for g in list(groups) + ["total"]:
        sl = [g] if g != "total" else list(groups)
        for k in ("F_bc", "A_bc"):
            dts = -DT * sum(r[f"{k}_{s}"].sum() for s in sl)
            assert sum(out[g][f"{k}_out_J"].values()) == pytest.approx(dts, rel=1e-9)
        fa = {b: out[g]["F_bc_out_J"][b] + out[g]["A_bc_out_J"][b] for b in BANDS}
        for b in BANDS:
            assert out[g]["FA_bc_out_J"][b] == pytest.approx(fa[b], rel=1e-12, abs=1e-6)
        assert out[g]["dt_sum_out_J"]["FA_bc"] == pytest.approx(sum(fa.values()), rel=1e-9)


def test_band_flux_puts_a_pure_frequency_in_its_band():
    """u' and p' at one mid-band frequency: all the flux is in 'mid' (cospectrum)."""
    t = np.arange(NT) * DT
    up = np.cos(2 * np.pi * F_MID * t)[:, None, None] * np.ones((1, 2, 3))
    pp = 2 * np.cos(2 * np.pi * F_MID * t + 0.3)[:, None, None] * np.ones((1, 2, 3))
    fs = {"up": up, "pp": pp, "un": up, "e1": pp, "dz": np.ones((2, 3)), "w": np.ones(3),
          "group": np.array(["a", "a", "b"]), "t": t,
          "fbt": np.zeros((NT, 3)), "abt": np.zeros((NT, 3))}
    out = pa.band_flux_series(fs, BANDS)
    tot = out["total"]["F_bc_out_J"]
    assert abs(tot["mid"]) > 0 and abs(tot["low"]) < 1e-9 * abs(tot["mid"])
    assert abs(tot["high"]) < 1e-9 * abs(tot["mid"])
    assert tot["mid"] == pytest.approx(-np.cos(0.3) * 2 * 3 * NT * DT, rel=1e-9)   # <up pp> = cos/2 x 2


def test_band_flux_requires_t0_zero(syn):
    ds, mask, groups, rho_b, tau, r = syn
    fs = pa.face_series(ds, mask, groups, rho_b=rho_b, time=slice(1, None))
    with pytest.raises(ValueError, match="start at 0"):
        pa.band_flux_series(fs, BANDS)


# ---- spectral helpers ---------------------------------------------------------------------------

def test_band_power_is_parseval_and_matches_the_cospectrum():
    rng = np.random.default_rng(3)
    x = rng.standard_normal((NT, 5))
    p = pa.band_power(x, DT, BANDS)
    np.testing.assert_allclose(sum(p.values()), (x ** 2).mean(0), rtol=1e-12)
    for c in range(5):
        cs = va.band_flux(x[:, None, c:c + 1], x[:, None, c:c + 1], np.ones((1, 1)), np.ones(1),
                          DT, BANDS)
        for b in BANDS:
            assert p[b][c] == pytest.approx(cs[b], rel=1e-10, abs=1e-14)


def test_band_pass_and_high_pass():
    t = np.arange(NT) * DT
    x = 3.0 + np.cos(2 * np.pi * F_LOW * t) + 2 * np.cos(2 * np.pi * F_MID * t) \
        + 0.5 * np.sin(2 * np.pi * F_HIGH * t)
    np.testing.assert_allclose(pa.band_pass(x, DT, *BANDS["mid"]), 2 * np.cos(2 * np.pi * F_MID * t),
                               atol=1e-12)
    np.testing.assert_allclose(pa.high_pass(x, DT, BANDS["low"][1]),
                               2 * np.cos(2 * np.pi * F_MID * t) + 0.5 * np.sin(2 * np.pi * F_HIGH * t),
                               atol=1e-12)
    xx = np.stack([x, 2 * x], 1)
    np.testing.assert_allclose(pa.band_pass(xx, DT, *BANDS["mid"])[:, 1],
                               4 * np.cos(2 * np.pi * F_MID * t), atol=1e-12)


# ---- box energy (H4b(i)) ------------------------------------------------------------------------

def test_level_energy_all_bands_equals_the_kf_budget_time_mean(syn):
    """Sum over bands of the baroclinic KE and the exact APE of every record (no filter) equal the
    record means of kf_budget's ke_bc and ape."""
    ds, mask, groups, rho_b, tau, r = syn
    e = pa.level_energy(ds, mask, rho_b, BANDS, progress=False)
    assert e["totals"]["ke_bc_all_J"] == pytest.approx(r["ke_bc"].mean(), rel=1e-10)
    assert e["totals"]["ape_exact_all_J"] == pytest.approx(r["ape"].mean(), rel=1e-10)
    assert sum(e["totals"]["ke_bc_band_J"].values()) == pytest.approx(r["ke_bc"].mean(), rel=1e-10)


def test_level_energy_without_the_low_band_equals_kf_on_high_passed_fields(syn):
    """H4b(i) (ruling): u, v and rho' are high-passed at 2 df per point before forming KE_bc and
    the exact APE; the result equals kf_budget run on the high-passed fields."""
    ds, mask, groups, rho_b, tau, r = syn
    cut = BANDS["low"][1]
    hp = ds.copy()
    hp["U"] = (ds["U"].dims, pa.high_pass(ds["U"].values, DT, cut))
    hp["V"] = (ds["V"].dims, pa.high_pass(ds["V"].values, DT, cut))
    rp = ds["RHO"].values - rho_b[None, :, None, None]
    hp["RHO"] = (ds["RHO"].dims, rho_b[None, :, None, None] + pa.high_pass(rp, DT, cut))
    r_hp = budget.kf_budget(hp, mask, rho_b=rho_b)
    e = pa.level_energy(ds, mask, rho_b, BANDS, progress=False)
    assert e["totals"]["ke_bc_hp_J"] == pytest.approx(r_hp["ke_bc"].mean(), rel=1e-10)
    assert e["totals"]["ape_exact_hp_J"] == pytest.approx(r_hp["ape"].mean(), rel=1e-10)
    assert e["totals"]["E_in_hp_J"] == pytest.approx(r_hp["ke_bc"].mean() + r_hp["ape"].mean(), rel=1e-10)
    assert e["totals"]["ke_bc_hp_J"] == pytest.approx(
        sum(v for b, v in e["totals"]["ke_bc_band_J"].items() if b != "low"), rel=1e-10)
    assert e["totals"]["ke_bc_hp_J"] < e["totals"]["ke_bc_all_J"]


# ---- wind work (H4c(i)) -------------------------------------------------------------------------

def test_wind_work_all_frequencies_equals_kf_wind_power(syn):
    """W over every frequency = int (P_bt + P_bc) dt of kf_budget (same points and weights)."""
    ds, mask, groups, rho_b, tau, r = syn
    pts = pa.wind_points(ds, mask)
    u, v = pa.surface_series(ds, pts)
    W = pa.wind_work(ds["time"].values, tau, pts, u, v)
    assert W == pytest.approx(budget.cumulative(r["t"], r["P_bt"] + r["P_bc"])[-1], rel=1e-10)


def test_wind_work_in_a_band_uses_the_band_passed_surface_velocity(syn):
    ds, mask, groups, rho_b, tau, r = syn
    pts = pa.wind_points(ds, mask)
    u, v = pa.surface_series(ds, pts)
    W_mid = pa.wind_work(ds["time"].values, tau, pts, u, v, band=BANDS["mid"])
    ub = pa.band_pass(u, DT, *BANDS["mid"])
    vb = pa.band_pass(v, DT, *BANDS["mid"])
    W_ref = pa.wind_work(ds["time"].values, tau, pts, ub, vb)
    assert W_mid == pytest.approx(W_ref, rel=1e-12)
    W_parts = sum(pa.wind_work(ds["time"].values, tau, pts, u, v, band=b) for b in BANDS.values())
    W_all = pa.wind_work(ds["time"].values, tau, pts, u, v)
    assert W_parts == pytest.approx(W_all, rel=1e-9)


# ---- stations and trapped fraction ---------------------------------------------------------------

def test_station_band_energy_finds_the_band_of_a_pure_oscillation():
    t = np.arange(NT) * DT
    u = np.cos(2 * np.pi * F_HIGH * t)
    v = -np.sin(2 * np.pi * F_HIGH * t)
    s = pa.station_band_energy(u, v, DT, BANDS)
    assert s["fraction"]["high"] > 0.99
    assert sum(s["band_m2_s-2"].values()) == pytest.approx(s["total_m2_s-2"], rel=1e-12)


def test_trapped_fraction():
    e = np.array([[1.0, 0.0], [1.0, 2.0], [2.0, 2.0]])       # (cells, bands) J
    d = np.array([100.0, 900.0, 5000.0])                       # m from land
    out = pa.trapped_fraction(e, d, r1_m=500.0, names=("a", "b"))
    assert out["a"] == pytest.approx(0.5) and out["b"] == pytest.approx(0.5)
    out2 = pa.trapped_fraction(e, d, r1_m=50.0, names=("a", "b"))
    assert out2["a"] == pytest.approx(0.0)


# ---- bands per run ------------------------------------------------------------------------------

def test_run_bands():
    th = va.load_theory()
    f1 = pa.run_bands("F1", th)
    assert f1["bands_hz"] == va.bands_hz(th, 280) and f1["resolvable"]["near_inertial"]
    r2 = pa.run_bands("R2", th)
    edges = {b: (lo * 3600, hi * 3600) for b, (lo, hi) in r2["bands_hz"].items()}
    assert edges["seiche"][0] == pytest.approx(0.149, abs=5e-4)               # PD0 3, bay at 78 N
    assert edges["near_inertial"][0] == pytest.approx(0.0695, abs=1e-4)
    assert edges["seiche"][1] == pytest.approx(2.0)                          # 900 s records
    r1k = pa.run_bands("R1k", th)
    assert r1k["bands_hz"]["near_inertial"][0] * 3600 == pytest.approx(0.0339, abs=1e-4)
    r1 = pa.run_bands("R1", th)                                              # 120 h: PD0 3
    assert not r1["resolvable"]["near_inertial"]
    assert r1["bands_hz"]["low"][1] * 3600 == pytest.approx(2 / 120)
    assert pa.run_bands("R1k_0_120h", th)["bands_hz"] == r1["bands_hz"]


def test_al_band_is_the_gullmar_band():
    assert pa.AL_BAND_CPH == (0.01, 0.07)


# ---- driver pieces ------------------------------------------------------------------------------

def test_bay_box_has_one_mouth_group_and_refuses_a_dry_box():
    nz, ny, nx = 2, 40, 30
    wet = np.zeros((ny, nx), bool)
    wet[:30] = True                                           # shelf
    wet[30:38, 10:20] = True                                  # bay, mouth on row 30
    hC = np.repeat(wet[None].astype(float), nz, 0)
    ds = xr.Dataset({"hFacC": (("k", "j", "i"), hC)})
    mask, groups = pa.bay_box(ds, rows=(30, 38), cols=(10, 20))
    assert mask.sum() == 80 and list(groups) == ["mouth"]
    assert groups["mouth"]["y"][30, 10:20].all() and groups["mouth"]["y"].sum() == 10
    assert groups["mouth"]["x"].sum() == 0
    shelf = (np.arange(ny)[:, None] < 30) * np.ones((1, nx))   # no bay: the box is land
    ds2 = xr.Dataset({"hFacC": (("k", "j", "i"), np.repeat(shelf[None], nz, 0))})
    with pytest.raises(ValueError, match="no wet cell"):
        pa.bay_box(ds2, rows=(30, 38), cols=(10, 20))


def test_tau_from_diagnostics_returns_the_record_at_that_time():
    t = np.arange(4) * DT
    tx = np.arange(4)[:, None, None] * np.ones((1, 3, 2))
    ds = xr.Dataset({"TAUX": (("time", "j", "i"), tx), "TAUY": (("time", "j", "i"), -tx)},
                    coords={"time": t})
    tau = pa.tau_from_diagnostics(ds)
    a, b = tau(2 * DT)
    assert np.all(a == 2) and np.all(b == -2)
    with pytest.raises(KeyError):
        tau(1.5 * DT)


def test_budget_summary_and_window(syn):
    """budget_summary integrates kf_budget over the window; a window ending earlier is the same as
    the budget of the shorter record (R1k over 0-120 h for H4d)."""
    ds, mask, groups, rho_b, tau, r = syn
    s = pa.budget_summary(r, list(groups))
    c = budget.closure(r)
    assert s["closure_J"]["bc"]["D"] == pytest.approx(c["bc"]["D"][-1])
    assert s["F_bc_pressure_out_J"] == pytest.approx(budget.f_bc_pressure_out(r)[-1])
    assert s["totals_in_J"]["FA_bc"] == pytest.approx(
        budget.cumulative(r["t"], r["F_bc"] + r["A_bc"])[-1])
    assert s["P_J"] == pytest.approx(budget.cumulative(r["t"], r["P_bt"] + r["P_bc"])[-1])
    rw = pa.window(r, 20 * H)
    r2 = budget.kf_budget(ds, mask, tau=tau, rho_b=rho_b, face_groups=groups, time=slice(0, 61))
    s2 = pa.budget_summary(r2, list(groups))
    sw = pa.budget_summary(rw, list(groups))
    assert sw["F_bc_pressure_out_J"] == pytest.approx(s2["F_bc_pressure_out_J"], rel=1e-12)
    assert sw["window_h"] == [0.0, 20.0]
    with pytest.raises(ValueError):
        pa.window(r, 20.5 * H)                                 # not a record time


# ---- local store of the box (one LaCie pass per run) ---------------------------------------------

def test_local_store_gives_the_same_kernels_as_the_full_dataset(syn, tmp_path):
    """kf_budget, face_series and level_energy on the stored box (cells within 2 of the box, wet
    levels only, cropped grid) equal the same kernels on the full dataset, exactly."""
    ds, mask, groups, rho_b, tau, r = syn
    st = pa.build_store(ds, mask, tmp_path / "store", progress=False)
    dl, ml, gl, crop = pa.store_dataset(ds, st, mask, groups)
    taul = pa.crop_tau(tau, crop)
    rl = budget.kf_budget(dl, ml, tau=taul, rho_b=rho_b, face_groups=gl)
    for k in ("ke_bc", "ape", "F_bc", "A_bc", "F_bt", "A_bt", "C", "P_bt", "P_bc", "F_bc_xfaces",
              "A_bc_yfaces"):
        np.testing.assert_array_equal(rl[k], r[k])
    a = pa.section_series(pa.face_series(dl, ml, gl, rho_b=rho_b, progress=False))
    b = pa.section_series(pa.face_series(ds, mask, groups, rho_b=rho_b, progress=False))
    for g in groups:
        for k in ("F_bc", "A_bc", "F_bt", "A_bt"):
            np.testing.assert_array_equal(a[g][k], b[g][k])
    e1 = pa.level_energy(dl, ml, rho_b, BANDS, progress=False)["totals"]
    e2 = pa.level_energy(ds, mask, rho_b, BANDS, progress=False)["totals"]
    assert e1 == e2
    pa.remove_store(st)
    assert not (tmp_path / "store").exists()


def test_assemble_from_synthetic_caches(syn, tmp_path, monkeypatch):
    """partition.json / budget.json from caches written as compute() writes them."""
    import json
    ds, mask, groups, rho_b, tau, r = syn
    monkeypatch.setattr(pa, "run_bands", lambda name, th: {
        "bands_hz": BANDS, "resolvable": {b: True for b in BANDS}, "n_bins": {b: 3.0 for b in BANDS},
        "record_h": 40.0, "dt_out_s": DT, "f_cph": 0.08})
    d = {k: v for k, v in r.items() if isinstance(v, np.ndarray)}
    d.update(mask=mask, sections=np.array(sorted(groups)), tau_source=np.array("forcing files"))
    pa._save(pa.cache_file(tmp_path, "kf", "F1"), d)
    pa._save(pa.cache_file(tmp_path, "faces", "F1"), pa.face_series(ds, mask, groups, rho_b=rho_b, progress=False))
    t = ds["time"].values
    rng = np.random.default_rng(1)
    pa._save(pa.cache_file(tmp_path, "surface", "F1"), {
        "station_names": np.array(["shelf", "coast00", "coast01"]), "station_ji": np.array([[1, 1], [4, 5], [5, 6]]),
        "station_u": rng.standard_normal((t.size, 3)), "station_v": rng.standard_normal((t.size, 3)), "t": t,
        "wind_all_J": 5.0, "wind_band_names": np.array(list(BANDS)), "wind_band_J": np.array([1.0, 2.0, 2.0]),
        "tau_source": np.array("forcing files")})
    part, bud = pa.assemble(tmp_path, {}, runs=("F1", "F2"), out_dir=tmp_path)
    assert part["missing"] == ["F2"]
    p = part["runs"]["F1"]
    assert p["check_face_series_vs_kf_budget_max_rel"] < 1e-9
    tot = p["straits"]["total"]
    assert sum(tot["FA_bc_out_J"].values()) == pytest.approx(-DT * (r["F_bc"] + r["A_bc"]).sum(), rel=1e-9)
    assert p["wind_work_J"]["bands"]["mid"] == 2.0
    assert p["wind_work_J"]["check_all_vs_kf_P_J"] == pytest.approx(
        budget.cumulative(r["t"], r["P_bt"] + r["P_bc"])[-1])
    assert set(p["stations"]) >= {"shelf", "coast00", "coast_south_mean"}
    al = p["arneborg_liljebladh_band_post_hoc"]
    assert al["label"].startswith("POST HOC") and "total" in al
    b = bud["runs"]["F1"]
    assert b["F_bc_pressure_out_J"] == pytest.approx(budget.f_bc_pressure_out(r)[-1])
    assert len(b["series"]["t_h"]) == len(range(0, NT, 3))
    saved = json.loads((tmp_path / "partition.json").read_text())
    assert saved["runs"]["F1"]["straits"]["total"]["FA_bc_out_J"]["mid"] == pytest.approx(tot["FA_bc_out_J"]["mid"])


def test_cache_keys_are_frozen_at_import(tmp_path, monkeypatch):
    """A process keeps the keys of the code it imported: an edit of partition.py while it runs
    must not rename its caches (inspect.getsource reads the current file with the old line
    numbers; this renamed the first F1 caches)."""
    k = pa.cache_file(tmp_path, "faces", "F1").name
    monkeypatch.setitem(pa.KEYS, "faces", lambda: "000000000000")
    assert pa.cache_file(tmp_path, "faces", "F1").name == k


def test_cache_keys_include_the_run_driver(monkeypatch):
    """Review fix 8: the box, the background profile, the stress choice and the store are part of
    what makes a kf/faces/energy/surface cache, so their source keys it."""
    for f in (pa.compute, pa._box, pa.bay_box, pa.tau_from_diagnostics, pa.build_store,
              pa.store_dataset, pa._StoreArray, pa.crop_tau):
        assert f in pa.DRIVER
    before = {k: f() for k, f in pa.KEYS.items()}
    monkeypatch.setattr(pa, "DRIVER", (pa.window,))
    after = {k: f() for k, f in pa.KEYS.items()}
    assert all(before[k] != after[k] for k in before)
