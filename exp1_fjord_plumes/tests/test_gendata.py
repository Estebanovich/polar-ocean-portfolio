"""Tests for the Exp. 1 input generators (Task 3)."""
from __future__ import annotations

import json
import re

import numpy as np
import pytest

from pfcommon.mitgcm_io import read_bin

from exp1_fjord_plumes.gendata.grids import icefront_fields, stretched_x
from exp1_fjord_plumes.gendata.make_case import load_cases, make_case
from exp1_fjord_plumes.gendata.profiles import continuous_sermilik, two_layer
from exp1_fjord_plumes.gendata.sources import discharge_addmass, place_sources
from exp1_fjord_plumes.theory.plume import LinearEOS

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
    assert "deltaT=1.0" in (out / "data").read_text()


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
    assert names == exp


@pytest.mark.parametrize("name", ["a_line2d_sermilik_q0.03", "b_3d_twolayer_x0_0"])
def test_make_case_other_configs(tmp_path, name):
    tdir = tmp_path / "tmpl"
    _templates(tdir)
    out = make_case(name, tmp_path / "run", template_dir=tdir)
    meta = json.loads((out / "case.json").read_text())
    assert meta["sum_addmass_kg_s"] > 0


def test_c_addmass_is_bottom_source_of_rho0_Q0(tmp_path):
    tdir = tmp_path / "tmpl"
    _templates(tdir)
    out = make_case("c_fjord_plume", tmp_path / "run", template_dir=tdir)
    meta = json.loads((out / "case.json").read_text())
    p = meta["params"]
    shape = (p["nz"], p["ny"], p["nx"])
    add = read_bin(out / "addmass.bin", shape)
    assert abs(add.sum() - 1027.0 * 100.0) < 1e-6
    k, j, i = np.nonzero(add)
    assert np.all(k == p["nz"] - 1) and np.all(i == p["front_i"])
    assert meta["plume_representation"] == "bottom_source_placeholder"
    assert "sink_cells" not in meta


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
from exp1_fjord_plumes.gendata.plume_source_sink import plume_source_sink  # noqa: E402

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
