"""Tests for exp2_isomip_cavity/gendata/make_case.py (Task 1; ruling R1 placeholder templates)."""
from __future__ import annotations

import json

import numpy as np
import pytest

from exp2_isomip_cavity.gendata import geometry as geo
from exp2_isomip_cavity.gendata import make_case as mc
from exp2_isomip_cavity.gendata import profiles as pr
from pfcommon.mitgcm_io import read_bin

CAL = ["cal2_gT0.015", "cal2_gT0.03", "cal2_gT0.06", "cal4_gT0.015", "cal4_gT0.03", "cal4_gT0.06",
       "cal2_gT0.0075", "cal4_gT0.0075"]   # 0.0075: PD2 sampling amendment (2026-11-09)


def _make(name, tmp_path, synth_geom, templates, **kw):
    out = tmp_path / name
    meta = mc.make_case(name, out, template_dir=templates, geom_path=synth_geom, **kw)
    return out, meta


# ---- registry ----------------------------------------------------------------------

def test_registry_has_all_cases():
    cases = mc.load_cases()["cases"]
    for n in ["ocean0_2km", "ocean1_2km", "ocean0_4km", *CAL]:
        assert n in cases
    assert cases["ocean0_2km"]["config"] == "ocean0_2km"
    assert all(cases[n]["config"] == "ocean0_4km" for n in cases if n.startswith(("ens_", "cal4")))


def test_ensemble_has_29_members_in_groups():
    spec = mc.load_cases()
    ens = mc.ensemble_cases(spec)
    assert len(ens) == 29 and len(set(ens)) == 29
    assert ens[0] == "ocean0_4km"
    cases = spec["cases"]
    groups = [cases[n]["ensemble_group"] for n in ens]
    assert groups.count("reference") == 1 and groups.count("factorial") == 15
    assert groups.count("geometry") == 10 and groups.count("extrapolation") == 3
    fact = {(cases[n]["Tbot"], cases[n]["z_c"]) for n in ens if cases[n]["ensemble_group"] == "factorial"}
    assert fact == {(T, z) for T in (-0.5, 0.0, 0.5, 1.0, 1.5) for z in (250.0, 400.0, 550.0)}
    ext = [cases[n] for n in ens if cases[n]["ensemble_group"] == "extrapolation"]
    assert sorted((c["Tbot"], c["z_c"], c["s"]) for c in ext) == [(2.0, 250.0, 1.0), (2.0, 400.0, 1.0), (2.0, 550.0, 1.0)]
    for n in ens[1:]:
        c = cases[n]
        assert n == f"ens_T{mc.fmt_num(c['Tbot'])}_zc{c['z_c']:.0f}_s{c['s']:.1f}"


def test_latin_hypercube_deterministic_and_stratified():
    a, b = mc.load_cases(), mc.load_cases()
    geo_a = [n for n in mc.ensemble_cases(a) if a["cases"][n]["ensemble_group"] == "geometry"]
    geo_b = [n for n in mc.ensemble_cases(b) if b["cases"][n]["ensemble_group"] == "geometry"]
    assert geo_a == geo_b
    for s in (0.8, 1.2):
        sub = [a["cases"][n] for n in geo_a if a["cases"][n]["s"] == s]
        assert len(sub) == 5
        T = np.array([c["Tbot"] for c in sub])
        z = np.array([c["z_c"] for c in sub])
        assert T.min() >= -0.5 and T.max() <= 1.5 and z.min() >= 250 and z.max() <= 550
        # one sample per fifth of each range (Latin hypercube), up to the name rounding
        assert sorted(np.minimum(((T + 0.5) / 2.0 * 5).astype(int), 4)) == [0, 1, 2, 3, 4]
        assert sorted(np.minimum(((z - 250) / 300.0 * 5).astype(int), 4)) == [0, 1, 2, 3, 4]
        assert all(c["lhs_seed"] == 20261108 for c in sub)


# ---- generated input ----------------------------------------------------------------

@pytest.mark.parametrize("name,nx,ny", [("cal2_gT0.03", 240, 41), ("cal4_gT0.03", 120, 21)])
def test_make_case_writes_inputs(tmp_path, synth_geom, templates, name, nx, ny):
    out, meta = _make(name, tmp_path, synth_geom, templates)
    for f in ("bathy.bin", "icetopo.bin", "shelficeLoadAnomaly.bin"):
        assert read_bin(out / f, (ny, nx)).shape == (ny, nx)
    for f in ("T.init.bin", "S.init.bin", "T.restore.bin", "S.restore.bin", "rbcs_mask.bin"):
        assert read_bin(out / f, (36, ny, nx)).shape == (36, ny, nx)
    data = (out / "data").read_text()
    assert f"delX = {nx}*{2000.0 * 240 / nx}" in data
    assert "nIter0 = 0," in data
    assert f"nTimeSteps = {meta['params']['nTimeSteps']}," in data
    assert "${" not in data
    shelf = (out / "data.shelfice").read_text()
    assert "SHELFICEheatTransCoeff = 0.03," in shelf
    assert f"SHELFICEsaltTransCoeff = {0.03 / 35}," in shelf
    on_disk = json.loads((out / "case.json").read_text())
    assert on_disk == json.loads(json.dumps(meta, default=float))


def test_case_json_records_actual_parameters(tmp_path, synth_geom, templates):
    out, meta = _make("cal2_gT0.03", tmp_path, synth_geom, templates)
    a = meta["actual"]
    assert a["s"] == 1.0 and a["gammaT"] == 0.03 and a["gammaS"] == pytest.approx(0.03 / 35)
    assert a["Tbot"] == 1.0 and a["z_c"] == 360.0 and a["d"] == 720.0
    assert a["bed_source"] == "GFZ Ocean1_input_geom_v1.01.nc bedrockTopography"
    wet = read_bin(out / "bathy.bin", (41, 240)) < 0
    assert a["n_wet_columns"] == int(wet.sum())
    assert a["n_wet_cells"] > a["n_wet_columns"]
    assert 100.0 < a["grounding_line_depth_m"] < 720.0
    assert len(a["geom_sha256"]) == 64
    assert a["deepest_wet_ice_base_depth_m"] == pytest.approx(-a["draft_min_m"])
    assert a["n_floating_columns"] > 0


def test_geometry_files_match_geometry_module(tmp_path, synth_geom, templates):
    out, meta = _make("cal4_gT0.03", tmp_path, synth_geom, templates)
    xc, yc = mc.grid_centres(4.0)
    bed, draft = geo.load_bed(str(synth_geom), xc, yc), geo.load_draft(str(synth_geom), xc, yc)
    fl = geo.load_floating(str(synth_geom), xc, yc)
    b2, d2, wet = geo.apply_scale_and_mask(bed, draft, 1.0, np.full(36, 20.0), xc_m=xc, floating=fl,
                                           hfac_min=meta["params"]["hFacMin"])
    assert meta["actual"]["geometry_snapped_to_hfac"] and meta["params"]["hFacMin"] == 0.2
    np.testing.assert_array_equal(read_bin(out / "bathy.bin", (21, 120))[:20], b2)
    np.testing.assert_array_equal(read_bin(out / "icetopo.bin", (21, 120))[:20], d2)
    assert not read_bin(out / "bathy.bin", (21, 120))[20].any()           # land row (R12)
    T, S = pr.ocean0(-(np.arange(36) + 0.5) * 20.0)
    pf = pr.load_anomaly(T, S, np.full(36, 20.0))
    np.testing.assert_allclose(read_bin(out / "shelficeLoadAnomaly.bin", (21, 120))[:20],
                               np.where(wet, pr.load_anomaly_at_draft(pf, np.full(36, 20.0), d2), 0.0))


def test_rbcs_mask_linear_ramp(tmp_path, synth_geom, templates):
    # AD16 Eq. 20: gamma = gamma0 max(0, (x - 790 km) / 10 km), at cell centres
    out, meta = _make("cal2_gT0.03", tmp_path, synth_geom, templates)
    m = read_bin(out / "rbcs_mask.bin", (36, 41, 240))
    xc, _ = mc.grid_centres(2.0)
    np.testing.assert_allclose(m[0, 0], np.clip((xc - 790e3) / 10e3, 0, None))
    assert np.all(m[:, :40] == m[:1, :1, :]) and not m[:, 40].any()     # land row (R12) is 0
    assert meta["params"]["tauRelax"] == pytest.approx(8640.0)   # 0.1 day: gamma0 = 10 day^-1


def test_ocean1_cold_initial_warm_restoring(tmp_path, synth_geom, templates):
    out, meta = _make("ocean1_2km", tmp_path, synth_geom, templates, overrides={"gammaT": 0.03})
    z = -(np.arange(36) + 0.5) * 20.0
    Ti, Si = read_bin(out / "T.init.bin", (36, 41, 240)), read_bin(out / "S.init.bin", (36, 41, 240))
    Tr, Sr = read_bin(out / "T.restore.bin", (36, 41, 240)), read_bin(out / "S.restore.bin", (36, 41, 240))
    np.testing.assert_allclose(Ti[:, 5, 5], pr.cold(z)[0]); np.testing.assert_allclose(Si[:, 5, 5], pr.cold(z)[1])
    np.testing.assert_allclose(Tr[:, 5, 5], pr.ocean0(z)[0]); np.testing.assert_allclose(Sr[:, 5, 5], pr.ocean0(z)[1])
    assert meta["params"]["endTime"] == pytest.approx(20 * 365 * 86400.0)
    assert meta["actual"]["init_profile"] == "cold" and meta["actual"]["restore_profile"] == "ocean0"


def test_ensemble_case_ramp_and_scale(tmp_path, synth_geom, templates):
    spec = mc.load_cases()
    name = next(n for n in mc.ensemble_cases(spec) if spec["cases"][n].get("s") == 1.2)
    out, meta = _make(name, tmp_path, synth_geom, templates, overrides={"gammaT": 0.03})
    c = spec["cases"][name]
    z = -(np.arange(36) + 0.5) * 20.0
    T = read_bin(out / "T.init.bin", (36, 21, 120))
    np.testing.assert_allclose(T[:, 3, 3], pr.ramp(z, c["Tbot"], c["z_c"], d=300.0)[0])
    assert meta["actual"]["s"] == 1.2 and meta["actual"]["Tbot"] == c["Tbot"]
    ref, _ = _make("cal4_gT0.03", tmp_path, synth_geom, templates)
    n12 = int((read_bin(out / "bathy.bin", (21, 120)) < 0).sum())
    n10 = int((read_bin(ref / "bathy.bin", (21, 120)) < 0).sum())
    assert n12 <= n10


def test_uncalibrated_gamma_refused(tmp_path, synth_geom, templates):
    with pytest.raises(ValueError, match="gammaT"):
        _make("ocean0_2km", tmp_path, synth_geom, templates)
    assert not (tmp_path / "ocean0_2km").exists()


def test_unknown_placeholder_raises(tmp_path, synth_geom, templates):
    (templates / "data.extra").write_text("x = ${no_such_parameter}\n")
    with pytest.raises(KeyError):
        _make("cal2_gT0.03", tmp_path, synth_geom, templates)


def test_missing_binary_reference_raises_before_writing(tmp_path, synth_geom, templates):
    (templates / "data.extra").write_text(" file = 'not_written.bin',\n")
    with pytest.raises(RuntimeError, match="not_written.bin"):
        _make("cal2_gT0.03", tmp_path, synth_geom, templates)
    assert not (tmp_path / "cal2_gT0.03").exists()


def test_nonintegral_step_count_rejected(tmp_path, synth_geom, templates):
    with pytest.raises(ValueError, match="deltaT"):
        _make("cal2_gT0.03", tmp_path, synth_geom, templates, overrides={"deltaT": 7.0})


def test_default_template_dir_is_config(tmp_path, synth_geom, monkeypatch):
    monkeypatch.setattr(mc, "PKG", tmp_path / "pkg")
    # both resolutions share config/ocean0/input_templates (cases.yaml template_config)
    with pytest.raises(FileNotFoundError, match="config/ocean0/input_templates"):
        mc.make_case("cal2_gT0.03", tmp_path / "x", geom_path=synth_geom)
