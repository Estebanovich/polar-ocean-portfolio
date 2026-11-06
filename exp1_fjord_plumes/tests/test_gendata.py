"""Tests for the Exp. 1 input generators (Task 3)."""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pytest

from pfcommon.mitgcm_io import read_bin

from exp1_fjord_plumes.gendata.grids import icefront_fields, stretched_x
from exp1_fjord_plumes.gendata.make_case import load_cases, make_case
from exp1_fjord_plumes.gendata.plume_source_sink import plume_source_sink
from exp1_fjord_plumes.gendata.profiles import continuous_sermilik, two_layer
from exp1_fjord_plumes.gendata.sources import discharge_addmass, place_sources
from exp1_fjord_plumes.theory.plume import LinearEOS

PKG_DIR = Path(__file__).resolve().parents[1]

Z = -(np.arange(60) + 0.5) * 10.0  # 600 m, dz = 10


def _n2(T, S, z):
    eos = LinearEOS()
    rho = eos.rho(T, S)
    return -eos.g / eos.rho0 * np.gradient(rho, z)


def test_profiles_statically_stable():
    for T, S in (two_layer(Z, 200.0, -1.5, 32.9, 4.0, 34.6, 10.0),
                 continuous_sermilik(Z)):
        assert np.all(_n2(T, S, Z) >= 0.0)


def test_two_layer_interface_depth():
    z_int, dz = 200.0, 10.0
    T, S = two_layer(Z, z_int, -1.5, 32.9, 4.0, 34.6, 5.0)
    n2 = _n2(T, S, Z)
    assert abs(-Z[np.argmax(n2)] - z_int) <= dz


def test_addmass_sums_to_discharge():
    shape = (5, 4, 6)
    hf = np.ones(shape)
    cells = [(4, 1, 0), (4, 2, 0)]
    a = discharge_addmass(shape, cells, 10.0, 1027.0, hFacC=hf, front_i=0)
    assert abs(a.sum() - 1027.0 * 10.0) < 1e-9
    assert a[4, 1, 0] == a[4, 2, 0]


def test_addmass_rejects_land_and_nonadjacent():
    shape = (5, 4, 6)
    hf = np.ones(shape)
    hf[4, 0, 0] = 0.0
    kw = dict(hFacC=hf, front_i=0)
    with pytest.raises(ValueError, match=re.escape("(4, 0, 0)")):
        discharge_addmass(shape, [(4, 0, 0)], 1.0, 1027.0, **kw)
    with pytest.raises(ValueError, match=re.escape("(4, 1, 3)")):
        discharge_addmass(shape, [(4, 1, 3)], 1.0, 1027.0, **kw)
    with pytest.raises(ValueError, match=re.escape("(2, 1, 0)")):  # above bottom
        discharge_addmass(shape, [(2, 1, 0)], 1.0, 1027.0, **kw)


def test_place_sources_reports_actual_separation():
    yc = (np.arange(120) + 0.5) * 10.0
    js, x0 = place_sources(yc, 25.0, 10.0)
    assert x0 in (20.0, 30.0)
    assert yc[js[1]] - yc[js[0]] == x0
    assert place_sources(yc, 40.0, 10.0)[1] == 40.0
    assert place_sources(yc, 0.0, 10.0)[0][0] == place_sources(yc, 0.0, 10.0)[0][1]


def test_stretched_x_geometric_then_capped():
    dx = stretched_x(5.0, 200.0, 1.03, 20000.0)
    assert dx[0] == 5.0 and dx.max() == 200.0
    assert np.allclose(dx[1] / dx[0], 1.03) and dx.sum() >= 20000.0


def test_icefront_length_definition():
    dx = stretched_x(5.0, 200.0, 1.03, 2000.0)
    R, length = icefront_fields(len(dx), 3, 50, dx, 10.0, 0, 500.0)
    assert length[0, 0] == 1.0 / dx[0]
    assert np.all(R[:, 0] == 500.0) and np.all(R[:, 1:] == 0.0)
    assert np.all(length[:, 1:] == 0.0)


def _templates(d):
    d.mkdir()
    (d / "data").write_text(
        "&PARM01\n tAlpha=${tAlpha}, deltaT=${deltaT},\n nTimeSteps=${nTimeSteps},\n"
        " hydrogThetaFile='T.init.bin',\n hydrogSaltFile='S.init.bin',\n&\n"
        "&PARM05\n delXfile='delX.bin', delYfile='delY.bin', delRfile='delZ.bin',\n"
        " bathyFile='topog.bin',\n&\n")
    (d / "data.icefront").write_text(
        "&ICEFRONT_PARM01\n ICEFRONTdepthFile='R_icefront.bin',\n"
        " ICEFRONTlengthFile='icefrontlength.bin',\n&\n")
    (d / "data.addmass").write_text(" addMassFile='addmass.bin',\n")


def test_make_case_writes_case_json(tmp_path):
    tdir = tmp_path / "tmpl"
    _templates(tdir)
    out = make_case("b_3d_twolayer_x0_40", tmp_path / "run", template_dir=tdir)
    meta = json.loads((out / "case.json").read_text())
    assert meta["x0_actual"] == 40.0
    assert abs(meta["sum_addmass_kg_s"] - 1027.0 * 20.0) < 1e-6
    for f in re.findall(r"'([^']+\.bin)'", "".join(p.read_text() for p in out.iterdir()
                                                   if p.name.startswith("data"))):
        assert (out / f).is_file()
    assert "deltaT=2.0" in (out / "data").read_text()  # PD1(b): dt = 2 s


def test_make_case_missing_file_raises(tmp_path):
    tdir = tmp_path / "tmpl"
    _templates(tdir)
    (tdir / "data.extra").write_text(" foo='nothere.bin',\n")
    with pytest.raises(RuntimeError, match="nothere.bin"):
        make_case("b_3d_twolayer_single", tmp_path / "run", template_dir=tdir)


def test_make_case_x0_aligned_to_grid(tmp_path):
    tdir = tmp_path / "tmpl"
    _templates(tdir)
    out = make_case("b_3d_sermilik_x0_80", tmp_path / "run", template_dir=tdir)
    assert json.loads((out / "case.json").read_text())["x0_actual"] == 80.0


def test_cases_yaml_lists_all_cases():
    names = set(load_cases()["cases"])
    exp = {f"a_line2d_{s}_q{q}" for s in ("twolayer", "sermilik") for q in ("0.01", "0.03", "0.1")}
    exp |= {f"b_3d_{s}_x0_{x}" for s in ("twolayer", "sermilik") for x in (0, 20, 40, 80, 160, 320)}
    exp |= {f"b_3d_{s}_single" for s in ("twolayer", "sermilik")}
    exp |= {f"c_fjord_{k}" for k in ("plume", "wind", "shelf", "all")}
    exp |= {"conv_dx10", "conv_dx5", "conv_dx2p5", "sens_freeslip", "sens_xu2012"}  # R20
    assert names == exp


@pytest.mark.parametrize("name", ["a_line2d_sermilik_q0.03", "b_3d_twolayer_x0_0"])
def test_make_case_other_configs(tmp_path, name):
    tdir = tmp_path / "tmpl"
    _templates(tdir)
    out = make_case(name, tmp_path / "run", template_dir=tdir)
    meta = json.loads((out / "case.json").read_text())
    assert meta["sum_addmass_kg_s"] > 0


def test_make_case_rejects_source_on_land(tmp_path):
    import yaml
    spec = load_cases()
    spec["cases"]["b_land"] = {**spec["cases"]["b_3d_twolayer_single"],
                               "land_columns": [[60, 0]]}  # (j, i) of the source
    cf = tmp_path / "cases.yaml"
    cf.write_text(yaml.safe_dump(spec))
    tdir = tmp_path / "tmpl"
    _templates(tdir)
    with pytest.raises(ValueError, match="land"):
        make_case("b_land", tmp_path / "run", template_dir=tdir, cases_file=cf)


def test_b_sermilik_pycnocline_at_lower_layer_top(tmp_path):
    tdir = tmp_path / "tmpl"
    _templates(tdir)
    out = make_case("b_3d_sermilik_single", tmp_path / "run", template_dir=tdir)
    p = json.loads((out / "case.json").read_text())["params"]
    T = read_bin(out / "T.init.bin", (p["nz"], p["ny"], p["nx"]))[:, 0, 0]
    S = read_bin(out / "S.init.bin", (p["nz"], p["ny"], p["nx"]))[:, 0, 0]
    zc = -(np.arange(p["nz"]) + 0.5) * p["dz"]
    n2 = _n2(T, S, zc)
    assert abs(-zc[np.argmax(n2)] - 100.0) <= p["dz"]


def test_make_case_writes_nothing_when_file_missing(tmp_path):
    tdir = tmp_path / "tmpl"
    _templates(tdir)
    (tdir / "data.extra").write_text(" foo='nothere.bin',\n")
    with pytest.raises(RuntimeError):
        make_case("b_3d_twolayer_single", tmp_path / "run", template_dir=tdir)
    assert not (tmp_path / "run").exists() or not any((tmp_path / "run").iterdir())


# ---- Task 6: config (b) walls, x0 = 0 and _single definitions -------------

def _b_case(tmp_path, name):
    tdir = tmp_path / "tmpl"
    _templates(tdir)
    out = make_case(name, tmp_path / "run", template_dir=tdir)
    meta = json.loads((out / "case.json").read_text())
    p = meta["params"]
    return out, meta, p, (p["nz"], p["ny"], p["nx"])


def test_b_cross_fjord_walls_and_closed_x_wrap(tmp_path):
    """Land rows j = 0 and j = Ny-1 (side walls) and land column i = Nx-1."""
    out, meta, p, shape = _b_case(tmp_path, "b_3d_twolayer_x0_40")
    depth = -read_bin(out / "topog.bin", (p["ny"], p["nx"]))
    assert np.all(depth[0, :] == 0) and np.all(depth[-1, :] == 0)
    assert np.all(depth[:, -1] == 0)
    assert np.all(depth[1:-1, :-1] == p["H"])


def test_b_land_rows_cols_reject_source_on_wall(tmp_path):
    import yaml
    spec = load_cases()
    spec["cases"]["b_wall"] = {**spec["cases"]["b_3d_twolayer_single"],
                               "land_j_rows": [60]}  # wall through the source
    cf = tmp_path / "cases.yaml"
    cf.write_text(yaml.safe_dump(spec))
    tdir = tmp_path / "tmpl"
    _templates(tdir)
    with pytest.raises(ValueError, match="land"):
        make_case("b_wall", tmp_path / "run", template_dir=tdir, cases_file=cf)


def test_b_x0_0_is_one_cell_of_2Q0(tmp_path):
    out, meta, p, shape = _b_case(tmp_path, "b_3d_twolayer_x0_0")
    add = read_bin(out / "addmass.bin", shape)
    assert np.count_nonzero(add) == 1
    assert abs(add.max() - 1027.0 * 2 * p["Q0"]) < 1e-9
    assert meta["x0_actual"] == 0.0 and meta["Q0_total"] == 2 * p["Q0"]


def test_b_single_is_one_source_of_Q0(tmp_path):
    out, meta, p, shape = _b_case(tmp_path, "b_3d_twolayer_single")
    add = read_bin(out / "addmass.bin", shape)
    assert np.count_nonzero(add) == 1
    assert abs(add.sum() - 1027.0 * p["Q0"]) < 1e-9
    assert meta["Q0_total"] == p["Q0"] and meta["x0_actual"] is None


def test_b_widest_separation_far_from_side_walls(tmp_path):
    """x0 = 320 m: both sources >= 100 m from the wet-side faces of the walls."""
    out, meta, p, shape = _b_case(tmp_path, "b_3d_twolayer_x0_320")
    assert meta["x0_actual"] == 320.0
    dy = p["dy"]
    y_wall_s, y_wall_n = 1 * dy, (p["ny"] - 1) * dy   # faces of rows 0 and Ny-1
    for (_, j, _) in meta["source_cells"]:
        yc = (j + 0.5) * dy
        assert yc - y_wall_s >= 100.0 and y_wall_n - yc >= 100.0


# ---------------------------------------------------------------- Task 7 (c)
# Plume as entrainment sinks + outflow sources (Cowton et al. 2015, Sec. 2.1.3)

ZC = -(np.arange(50) + 0.5) * 10.0   # config (c): H = 500 m, dz = 10 m


def _c_plume(Q0=100.0, uniform=False):
    if uniform:
        Ta, Sa = np.full(ZC.shape, 4.0), np.full(ZC.shape, 34.6)
    else:
        Ta, Sa = continuous_sermilik(ZC)
    return plume_source_sink(ZC, Ta, Sa, Q0, 10.0, return_solution=True)


def test_source_sink_sums_to_Q0():
    """C5: sum of the per-level volume fluxes = Q0 within 1 %."""
    for Q0 in (10.0, 100.0, 1000.0):
        q, T_out, S_out, _ = _c_plume(Q0)
        assert abs(q.sum() - Q0) <= 0.01 * Q0


def test_source_sink_signs_and_outflow_levels():
    """Cowton et al. (2015, Sec. 2.1.3), whole cells: sinks (< 0) in every
    cell below the one holding z_neutral, the whole outflow (> 0) in that
    cell, nothing above. Outflow T/S = plume T/S at its bottom face."""
    q, T_out, S_out, sol = _c_plume(100.0)
    assert np.isfinite(sol.z_neutral)
    top, bot = ZC + 5.0, ZC - 5.0
    k_n = int(np.flatnonzero((bot < sol.z_neutral) & (sol.z_neutral <= top))[0])
    assert q[k_n] > 0 and np.all(q[k_n + 1:] < 0) and np.all(q[:k_n] == 0)
    z_e = bot[k_n]
    assert np.isclose(T_out[k_n], np.interp(z_e, sol.z, sol.T), atol=1e-6)
    assert np.isclose(S_out[k_n], np.interp(z_e, sol.z, sol.S), atol=1e-6)
    others = np.arange(ZC.size) != k_n
    assert np.all(np.isnan(T_out[others])) and np.all(np.isnan(S_out[others]))
    # outflow is a mixture of discharge and ambient: inside their T, S range
    Ta, Sa = continuous_sermilik(ZC)
    assert 0.0 < S_out[k_n] < Sa.max() and min(0.0, Ta.min()) < T_out[k_n] < Ta.max()


def test_source_sink_plume_reaching_surface():
    """Uniform ambient: the plume reaches the surface and all outflow goes to
    the top cell; the sum is still Q0."""
    q, T_out, S_out, sol = _c_plume(100.0, uniform=True)
    assert sol.reached_top
    assert q[0] > 0 and np.all(q[1:] <= 0)
    assert abs(q.sum() - 100.0) <= 1.0


def test_source_sink_outflow_over_n_levels():
    """n_out_levels = 2: the outflow Q(z_lo) is split equally between the
    z_neutral cell and the cell above it (into the overshoot); sinks and the
    sum are unchanged."""
    Ta, Sa = continuous_sermilik(ZC)
    q1, T1, S1 = plume_source_sink(ZC, Ta, Sa, 100.0, 10.0)
    q2, T2, S2 = plume_source_sink(ZC, Ta, Sa, 100.0, 10.0, n_out_levels=2)
    k_n = int(np.flatnonzero(q1 > 0)[0])
    assert np.allclose(q2[k_n - 1:k_n + 1], q1[k_n] / 2)
    assert np.allclose(q2[k_n + 1:], q1[k_n + 1:]) and np.all(q2[:k_n - 1] == 0)
    assert abs(q2.sum() - 100.0) < 1e-9
    assert np.allclose(T2[k_n - 1:k_n + 1], T1[k_n]) and np.isnan(T2[k_n - 2])


# ------------------------------------------------- Task 7 Step 4: config (c)
def _active(text):
    """Namelist text without comment lines."""
    return "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("#"))


def _c_case(tmp_path, name):
    out = make_case(name, tmp_path / name)   # real c_fjord3d templates
    meta = json.loads((out / "case.json").read_text())
    p = meta["params"]
    return out, meta, p, (p["nz"], p["ny"], p["nx"])


def test_c_plume_net_addmass_is_rho0_Q0(tmp_path):
    """C5: net addMass = rho0 Q0 within 1 %, recorded in case.json."""
    out, meta, p, shape = _c_case(tmp_path, "c_fjord_plume")
    add = read_bin(out / "addmass.bin", shape)
    assert abs(add.sum() - p["rho0"] * p["Q0"]) <= 0.01 * p["rho0"] * p["Q0"]
    assert meta["plume_representation"] == "bpt_sink_source_rbcs"
    assert abs(meta["net_addmass_kg_s"] - add.sum()) < 1e-6
    assert abs(meta["net_addmass_over_rho0_Q0"] - 1.0) <= 0.01


def test_c_sinks_only_below_outflow_in_front_column(tmp_path):
    out, meta, p, shape = _c_case(tmp_path, "c_fjord_plume")
    add = read_bin(out / "addmass.bin", shape)
    k, j, i = np.nonzero(add)
    assert np.all(i == p["front_i"])
    assert set(j) == set(meta["plume_j"]) and len(meta["plume_j"]) == p["n_plume_cols"]
    kpos = np.nonzero(add > 0)[0]
    kneg = np.nonzero(add < 0)[0]
    assert kneg.min() > kpos.max()                      # sinks only below outflow
    assert sorted(set(kpos)) == meta["outflow_k"]
    assert len(meta["outflow_k"]) == p["n_out_levels"]


def test_c_rbcs_mask_only_in_outflow_cells(tmp_path):
    """mask = q tau / V in the outflow cells (RBCS tendency = q (T_out - T)/V,
    PD3 option iii), 0 elsewhere, <= 1."""
    out, meta, p, shape = _c_case(tmp_path, "c_fjord_plume")
    add = read_bin(out / "addmass.bin", shape)
    mask = read_bin(out / "rbcs_mask.bin", shape)
    assert np.array_equal(mask > 0, add > 0)
    V = p["dx"] * p["dy"] * p["dz"]
    q = add / p["rho0"]
    assert np.allclose(mask[add > 0], q[add > 0] * p["rbcs_tau"] / V)
    assert mask.max() <= 1.0


def test_c_rbcs_targets_are_plume_outflow(tmp_path):
    out, meta, p, shape = _c_case(tmp_path, "c_fjord_plume")
    mask = read_bin(out / "rbcs_mask.bin", shape)
    Tr = read_bin(out / "rbcs_T.bin", shape)
    Sr = read_bin(out / "rbcs_S.bin", shape)
    zc = -(np.arange(p["nz"]) + 0.5) * p["dz"]
    Ta, Sa = continuous_sermilik(zc)
    _, T_out, S_out = plume_source_sink(zc, Ta, Sa, p["Q0"], p["dz"],
                                        T0=p["T0_src"], S0=p["S0_src"],
                                        alpha=p["alpha"], n_out_levels=p["n_out_levels"])
    kk = np.nonzero(mask > 0)[0]
    assert np.allclose(Tr[mask > 0], T_out[kk]) and np.allclose(Sr[mask > 0], S_out[kk])
    assert np.isclose(meta["T_out"], np.nanmax(T_out)) and np.isclose(meta["S_out"], np.nanmax(S_out))
    data = _active((out / "data").read_text())
    assert "temp_addMass" not in data and "salt_addMass" not in data
    assert "salt_EvPrRn=1.234567E+05" in data


def test_c_wind_case_has_no_plume_and_two_stages(tmp_path):
    out, meta, p, shape = _c_case(tmp_path, "c_fjord_wind")
    assert not np.any(read_bin(out / "addmass.bin", shape))
    assert not np.any(read_bin(out / "rbcs_mask.bin", shape))
    nrec = int(round(p["externForcingCycle"] / p["externForcingPeriod"]))
    w = read_bin(out / "windx.bin", (nrec, p["ny"], p["nx"]))
    n_on = int(round(p["wind_hours"] * 3600 / p["externForcingPeriod"]))
    assert np.all(w[:n_on] == p["wind_tau"]) and np.all(w[n_on:] == 0)
    s1, s2 = (out / "data").read_text(), (out / "data.stage2").read_text()
    n1 = int(round(p["wind_stage1_s"] / p["deltaT"]))
    assert "zonalWindFile='windx.bin'" in s1 and "nTimeSteps=%d," % n1 in s1
    assert "zonalWindFile" not in _active(s2)
    assert "nIter0=%d," % n1 in s2
    assert "nTimeSteps=%d," % (int(round(p["endTime"] / p["deltaT"])) - n1) in s2


def test_c_plume_case_single_stage(tmp_path):
    out, meta, p, shape = _c_case(tmp_path, "c_fjord_plume")
    assert not (out / "data.stage2").exists() and not (out / "windx.bin").exists()
    assert "nTimeSteps=%d," % int(round(p["endTime"] / p["deltaT"])) in (out / "data").read_text()


def test_c_shelf_ob_series_heaves_pycnocline(tmp_path):
    """OB T/S records: Sermilik profile with the pycnocline displaced by
    A sin(2 pi t/P) at the record times (n - 1/2) dt_rec (records centred,
    eesupp/src/get_periodic_interval.F)."""
    out, meta, p, shape = _c_case(tmp_path, "c_fjord_shelf")
    nrec = int(round(p["externForcingCycle"] / p["externForcingPeriod"]))
    assert p["externForcingCycle"] == p["shelf_period_days"] * 86400
    T = read_bin(out / "OBEt.bin", (nrec, p["nz"], p["ny"]))
    zc = -(np.arange(p["nz"]) + 0.5) * p["dz"]
    for n in (0, nrec // 4, nrec // 2):
        t = (n + 0.5) * p["externForcingPeriod"]
        eta = p["shelf_amp"] * np.sin(2 * np.pi * t / p["externForcingCycle"])
        Te, _ = continuous_sermilik(zc, z_c=200.0 - eta)
        assert np.allclose(T[n, :, 5], Te)
    # plume-only case: constant (initial) OB profile
    out2, _, p2, _ = _c_case(tmp_path, "c_fjord_plume")
    T2 = read_bin(out2 / "OBEt.bin", (nrec, p2["nz"], p2["ny"]))
    assert np.allclose(T2, T2[0])


def test_c_grid_matches_size_h_and_5_ranks(tmp_path):
    out, meta, p, shape = _c_case(tmp_path, "c_fjord_all")
    txt = (PKG_DIR / "config" / "c_fjord3d" / "code" / "SIZE.h").read_text()
    v = {k: int(x) for k, x in re.findall(r"(sNx|sNy|nSx|nSy|nPx|nPy|Nr) *= *(\d+)", txt)}
    assert v["sNx"] * v["nSx"] * v["nPx"] == p["nx"]
    assert v["sNy"] * v["nSy"] * v["nPy"] == p["ny"]
    assert v["Nr"] == p["nz"] and v["nPx"] * v["nPy"] == 5
    depth = -read_bin(out / "topog.bin", (p["ny"], p["nx"]))
    assert np.all(depth[0] == 0) and np.all(depth[-1] == 0) and np.all(depth[:, -1] == 0)
    obcs = _active((out / "data.obcs").read_text())
    assert "OB_Ieast=%d*-2," % p["ny"] in obcs and "useOBCSbalance" not in obcs


def test_source_sink_rejects_outflow_levels_above_surface():
    """R19: n_out_levels that would reach above the surface raises instead of
    being clipped silently (plume reaching the surface -> outflow cell k = 0)."""
    Ta, Sa = np.full(ZC.shape, 4.0), np.full(ZC.shape, 34.6)
    with pytest.raises(ValueError, match="n_out_levels"):
        plume_source_sink(ZC, Ta, Sa, 100.0, 10.0, n_out_levels=2)


def test_c_outflow_one_level_six_columns_not_denser_than_ambient(tmp_path):
    """R18: lateral spreading only (n_plume_cols = 6, n_out_levels = 1); the
    outflow is not denser than the ambient water of its cell (no inversion,
    Cowton et al. 2015 terminate at neutral buoyancy)."""
    out, meta, p, shape = _c_case(tmp_path, "c_fjord_plume")
    assert p["n_plume_cols"] == 6 and p["n_out_levels"] == 1
    assert len(meta["plume_j"]) == 6 and len(meta["outflow_k"]) == 1
    k = meta["outflow_k"][0]
    eos = LinearEOS()
    zc = -(np.arange(p["nz"]) + 0.5) * p["dz"]
    Ta, Sa = continuous_sermilik(zc)
    assert eos.rho(meta["T_out"], meta["S_out"]) <= eos.rho(Ta[k], Sa[k])
    V = p["dx"] * p["dy"] * p["dz"]
    assert abs(meta["outflow_turnover_s"] - V / (meta["Q_outflow"] / 6)) < 1e-6


def test_c_case_json_records_melt_volume_from_ambient(tmp_path):
    """R19 (ii): sum(q) = Q0 is reached with the plume melt volume taken from
    the ambient sinks; case.json states it."""
    out, meta, p, shape = _c_case(tmp_path, "c_fjord_plume")
    assert 0.1 < meta["melt_volume_in_sinks_m3_s"] < 2.0
    assert "melt" in meta["sum_rule_note"]


def test_c_wind_run_not_longer_than_stage1_is_single_stage(tmp_path):
    """endTime <= wind_stage1_s (e.g. a short run): one stage with the wind file."""
    out = make_case("c_fjord_wind", tmp_path / "w", overrides={"endTime": 43200.0})
    d = (out / "data").read_text()
    assert not (out / "data.stage2").exists()
    assert "zonalWindFile='windx.bin'" in d and "nTimeSteps=720," in d


def test_c_wind_stage1_must_align_with_dt_and_diagnostics(tmp_path):
    """Stage boundary on a multiple of deltaT, the 3-h mean and the hourly snapshot."""
    for bad in (5400.0, 1830.0):          # not a multiple of 10800 / of 3600 and 60
        with pytest.raises(ValueError, match="multiple"):
            make_case("c_fjord_wind", tmp_path / f"w{bad}",
                      overrides={"endTime": 86400.0, "wind_stage1_s": bad})
    from exp1_fjord_plumes.gendata.make_case import C_DIAG_PERIODS
    t = (PKG_DIR / "config" / "c_fjord3d" / "input_templates" / "data.diagnostics").read_text()
    assert f"frequency(1) = {C_DIAG_PERIODS[0]:.0f}.," in t
    assert f"frequency(5) = -{C_DIAG_PERIODS[1]:.0f}.," in t


# ---- Task 6 fix: (b) eastern OB = Stevens (1990) with prescribed depth-mean u ----

def _b_real_case(tmp_path, name="b_3d_twolayer_x0_40"):
    """Generate a (b) case with the REAL config templates (data.obcs etc.)."""
    out = make_case(name, tmp_path / "run")
    meta = json.loads((out / "case.json").read_text())
    p = meta["params"]
    return out, meta, p, (p["nz"], p["ny"], p["nx"])


def _wet_ob_area_col(out, p):
    """hFacW * drF at the OB u-point (0-based column nx-2) from the written topography."""
    nz, ny, nx = p["nz"], p["ny"], p["nx"]
    depth = -read_bin(out / "topog.bin", (ny, nx))
    k = (np.arange(nz) * p["dz"])[:, None, None]
    hFacC = (k < depth[None]).astype(float)
    i = nx - 2  # Fortran i = Nx-1 (OB_Ieast=-2)
    return np.minimum(hFacC[:, :, i], hFacC[:, :, i - 1]) * p["dz"]  # [m], (nz, ny)


@pytest.mark.parametrize("name", ["b_3d_twolayer_x0_40", "b_3d_twolayer_single",
                                  "b_3d_sermilik_x0_0"])
def test_b_stevens_obeu_carries_Q0_total(tmp_path, name):
    out, meta, p, (nz, ny, nx) = _b_real_case(tmp_path, name)
    u = read_bin(out / "OBEu.bin", (nz, ny))
    A_col = _wet_ob_area_col(out, p)
    Q = float(np.sum(u * A_col) * p["dy"])
    assert abs(Q - meta["Q0_total"]) / meta["Q0_total"] < 1e-10
    assert np.all(u >= 0.0) and np.all(u[A_col == 0] == 0.0)
    assert abs(meta["ob_area_m2"] - A_col.sum() * p["dy"]) < 1e-6
    assert abs(meta["ob_u_m_per_s"] - meta["Q0_total"] / meta["ob_area_m2"]) < 1e-14
    assert np.all(u[A_col > 0] == meta["ob_u_m_per_s"])


def test_b_stevens_obev_zero_and_tracers_are_ambient(tmp_path):
    out, meta, p, (nz, ny, nx) = _b_real_case(tmp_path)
    assert np.all(read_bin(out / "OBEv.bin", (nz, ny)) == 0.0)
    T0 = read_bin(out / "T.init.bin", (nz, ny, nx))[:, 0, 0]
    S0 = read_bin(out / "S.init.bin", (nz, ny, nx))[:, 0, 0]
    for f, ref in (("OBEt.bin", T0), ("OBEs.bin", S0)):
        a = read_bin(out / f, (nz, ny))
        assert np.all(a == ref[:, None])


def test_b_stevens_namelist_and_cpp(tmp_path):
    out, meta, p, shape = _b_real_case(tmp_path)
    txt = (out / "data.obcs").read_text()
    live = "\n".join(l for l in txt.splitlines() if not l.lstrip().startswith("#"))
    assert re.search(r"useStevensEast\s*=\s*\.TRUE\.", live)
    assert not re.search(r"useOrlanskiEast\s*=\s*\.TRUE\.", live)
    assert "OBCS_PARM02" not in live
    assert re.search(r"useOBCSprescribe\s*=\s*\.TRUE\.", live)
    for key in ("OBEuFile='OBEu.bin'", "OBEvFile='OBEv.bin'", "OBEtFile='OBEt.bin'",
                "OBEsFile='OBEs.bin'", "TrelaxStevens=86400.", "SrelaxStevens=86400."):
        assert key in live.replace(" ", ""), key
    assert "OBCS_monSelect" in live
    # constant single-record files are read once at nIter0 (obcs_fields_load.F:347-363):
    # (b) must not turn on periodic forcing
    d = (out / "data").read_text()
    assert not re.search(r"periodicExternalForcing\s*=\s*\.TRUE\.", d)
    h = (PKG_DIR / "config" / "b_twoplume3d" / "code" / "OBCS_OPTIONS.h").read_text()
    assert re.search(r"^#define ALLOW_OBCS_STEVENS", h, re.M)
    assert re.search(r"^#define ALLOW_OBCS_PRESCRIBE", h, re.M)


# ---- (c) Stevens eastern OB (task-7 fix) ---------------------------------------
@pytest.mark.parametrize("name", ["c_fjord_all", "c_fjord_plume", "c_fjord_shelf"])
def test_c_stevens_obeu_carries_Q0_total_every_record(tmp_path, name):
    out, meta, p, (nz, ny, nx) = _c_case(tmp_path, name)
    nrec = int(round(p["externForcingCycle"] / p["externForcingPeriod"]))
    u = read_bin(out / "OBEu.bin", (nrec, nz, ny))
    # same record count as the time-varying OBEt/OBEs (periodic read, obcs_fields_load.F:102,159-173)
    sizes = {f: (out / f).stat().st_size for f in ("OBEu.bin", "OBEv.bin", "OBEt.bin", "OBEs.bin")}
    assert len(set(sizes.values())) == 1 and sizes["OBEu.bin"] == nrec * nz * ny * 8, sizes
    assert nrec > 1 and np.all(u == u[0])                    # identical records
    A_col = _wet_ob_area_col(out, p)
    Q = float(np.sum(u[0] * A_col) * p["dy"])
    if p["Q0"] > 0:
        assert abs(Q - meta["Q0_total"]) / meta["Q0_total"] < 1e-10
    else:
        assert Q == 0.0
    assert meta["Q0_total"] == p["Q0"]
    assert np.all(u >= 0.0) and np.all(u[:, A_col == 0] == 0.0)
    assert abs(meta["ob_area_m2"] - A_col.sum() * p["dy"]) < 1e-6
    assert meta["ob_column_0based"] == nx - 2
    assert abs(meta["ob_u_m_per_s"] - p["Q0"] / meta["ob_area_m2"]) < 1e-14
    assert np.all(u[:, A_col > 0] == meta["ob_u_m_per_s"])
    assert np.all(read_bin(out / "OBEv.bin", (nrec, nz, ny)) == 0.0)


def test_c_stevens_namelist_and_cpp(tmp_path):
    out, meta, p, shape = _c_case(tmp_path, "c_fjord_all")
    live = _active((out / "data.obcs").read_text())
    flat = live.replace(" ", "")
    assert re.search(r"useStevensEast\s*=\s*\.TRUE\.", live)
    assert not re.search(r"useOrlanskiEast\s*=\s*\.TRUE\.", live)
    assert "OBCS_PARM02" not in live and "cvelTimeScale" not in live
    assert "&OBCS_PARM04" in flat
    for key in ("OBEuFile='OBEu.bin'", "OBEvFile='OBEv.bin'", "OBEtFile='OBEt.bin'",
                "OBEsFile='OBEs.bin'", "TrelaxStevens=86400.", "SrelaxStevens=86400.",
                "OBCSsponge_Theta=.TRUE.", "OBCSsponge_Salt=.TRUE.", "Urelaxobcsinner=21600.",
                "Urelaxobcsbound=600."):
        assert key in flat, key
    assert "OBCS_monitorFreq" in live
    d = _active((out / "data").read_text())
    assert re.search(r"periodicExternalForcing\s*=\s*\.TRUE\.", d)   # needed for the T/S series
    h = (PKG_DIR / "config" / "c_fjord3d" / "code" / "OBCS_OPTIONS.h").read_text()
    assert re.search(r"^#define ALLOW_OBCS_STEVENS", h, re.M)
    assert re.search(r"^#define ALLOW_OBCS_PRESCRIBE", h, re.M)
