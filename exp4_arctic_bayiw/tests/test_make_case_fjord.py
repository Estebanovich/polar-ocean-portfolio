"""Tests for the fjord_vmf cases of gendata.make_case (Task 4): F1-F5, F1c, F1g, V3, V4.

MITgcm is never run. The bathymetry is rebuilt from data_external/ (IBCAO, NPI); the tests
are skipped if it is not present. Large binaries (the 181-record wind files) are written
for one case only.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import yaml

from exp4_arctic_bayiw.gendata import make_case as mc
from pfcommon import campaign as cp

PKG = Path(__file__).resolve().parents[1]
DATA = PKG / "data_external"
SPEC = yaml.safe_load((PKG / "cases.yaml").read_text())
CFG = SPEC["configs"]["fjord_vmf"]
pytestmark = pytest.mark.skipif(not (DATA / CFG["bathy"]["ibcao_dir"]).exists(),
                                reason="IBCAO / NPI data not present in data_external/")
SMALL = ("bathy.bin", "T.init.bin", "S.init.bin", "HEFF.bin", "AREA.bin")
FJORD_CASES = ("F1", "F2", "F3", "F4", "F5", "F1c", "F1g", "V3", "V4")


@pytest.fixture(scope="module")
def cases(tmp_path_factory):
    root = tmp_path_factory.mktemp("fjord")
    out = {}
    for n in FJORD_CASES:
        wb = True if n == "F1" else SMALL + tuple(f for f in mc.TIDE_FILES)
        out[n] = (root / n, mc.make_case(n, root / n, write_bins=wb))
    return out


def _nml(d: Path, f: str = "data") -> dict:
    return mc.parse_namelist((d / f).read_text())


def test_all_fjord_cases_registered():
    for n in FJORD_CASES:
        assert SPEC["cases"][n]["config"].startswith("fjord_vmf")


def test_f1_pd0_mixing_eos_and_f_plane(cases):
    d, meta = cases["F1"]
    p1 = _nml(d)["PARM01"]
    assert float(p1["f0"]) == pytest.approx(1.4266e-4) and float(p1["beta"]) == 0.0
    assert float(p1["viscar"]) == 1e-4 and float(p1["diffkrt"]) == 1e-5 and float(p1["diffkrs"]) == 1e-5
    assert float(p1["viscc2smag"]) == 2.2 and float(p1["viscc4smag"]) == 1.0
    assert "viscah" not in p1 and float(p1["diffkht"]) == 0.0 and float(p1["diffkhs"]) == 0.0
    assert float(p1["bottomdragquadratic"]) == 2.5e-3
    assert p1["tempadvscheme"] == "77" and p1["saltadvscheme"] == "77"
    assert p1["eostype"] == "'LINEAR'"
    e = CFG["eos"]
    assert float(p1["talpha"]) == e["tAlpha"] and float(p1["sbeta"]) == e["sBeta"]
    assert float(p1["rhonil"]) == e["rhoNil"] and float(p1["rhoconst"]) == e["rhoNil"]
    assert p1["tref"] == f"40*{e['tRef']}" and p1["sref"] == f"40*{e['sRef']}"
    assert float(p1["hfacmin"]) == 0.2 and float(p1["hfacmindr"]) == 2.0
    assert "ivdc_kappa" not in p1 and p1["implicitviscosity"] == ".TRUE."
    pk = _nml(d, "data.pkg")["PACKAGES"]
    assert pk["usekpp"] == ".TRUE." and pk.get("useggl90", ".FALSE.") == ".FALSE."
    assert (d / "data.kpp").exists() and not (d / "data.ggl90").exists()
    assert meta["actual"]["profile"] == "vmf_summer"


def test_f1_grid_matches_bathymetry_and_size_h(cases):
    d, meta = cases["F1"]
    p4 = _nml(d)["PARM04"]
    delR = mc.parse_namelist_list((d / "data").read_text(), "delR")
    assert np.allclose(delR, CFG["delR"])
    assert p4["delx"] == "500*200.0" and p4["dely"] == "300*200.0"
    assert float(p4["xgorigin"]) == CFG["bbox_utm"][0] and float(p4["ygorigin"]) == CFG["bbox_utm"][2]
    b = np.fromfile(d / "bathy.bin", ">f8")
    ref = np.fromfile(DATA / "derived" / "bathy_fjord_vmf.bin", ">f8")
    assert np.array_equal(b, ref)                               # Task 1b bathymetry
    cp.check_size(d, PKG / "config" / "fjord_vmf" / "code" / "SIZE.h")
    cp.check_np(PKG / "config" / "fjord_vmf" / "code" / "SIZE.h", 5)


def test_f1_initial_state_is_the_summer_profile(cases):
    import json
    d, _ = cases["F1"]
    pj = json.loads((PKG / "figures" / "profiles_vmf.json").read_text())
    T = np.fromfile(d / "T.init.bin", ">f8").reshape(40, 300, 500)
    S = np.fromfile(d / "S.init.bin", ">f8").reshape(40, 300, 500)
    assert np.ptp(T, axis=(1, 2)).max() == 0.0 and np.ptp(S, axis=(1, 2)).max() == 0.0
    assert np.array_equal(T[:, 0, 0], pj["summer"]["T"]) and np.array_equal(S[:, 0, 0], pj["summer"]["S"])


def test_f1_wind_files_down_fjord(cases):
    d, meta = cases["F1"]
    n = 181
    tx = np.fromfile(d / "taux.bin", ">f8").reshape(n, 300, 500)
    ty = np.fromfile(d / "tauy.bin", ">f8").reshape(n, 300, 500)
    k = 7
    ang = np.radians(CFG["wind"]["angle_deg"])
    assert tx[k, 0, 0] == pytest.approx(0.1 * np.cos(ang)) and ty[k, 0, 0] == pytest.approx(0.1 * np.sin(ang))
    assert tx[k, 0, 0] < 0 and ty[k, 0, 0] < 0                  # towards WSW, down-fjord
    assert np.hypot(tx, ty).max() == pytest.approx(0.1)
    assert np.all(tx[-2:] == 0) and np.all(ty[-2:] == 0)
    p3 = _nml(d)["PARM03"]
    assert float(p3["externforcingperiod"]) == 1200.0
    assert float(p3["externforcingcycle"]) == n * 1200.0
    p5 = _nml(d)["PARM05"]
    assert p5["zonalwindfile"] == "'taux.bin'" and p5["meridwindfile"] == "'tauy.bin'"
    assert meta["actual"]["wind"]["max_tau_N_m-2"] == pytest.approx(0.1)


def test_free_stage_has_same_physics_without_wind(cases):
    d, _ = cases["F1"]
    a, b = _nml(d), _nml(d, "data.free")
    assert a["PARM01"] == b["PARM01"] and a["PARM04"] == b["PARM04"]
    assert "zonalwindfile" not in b["PARM05"] and "meridwindfile" not in b["PARM05"]
    assert b["PARM03"].get("periodicexternalforcing", ".FALSE.") == ".FALSE."
    st = cases["F1"][1]["actual"]["stages"]
    dt = float(a["PARM03"]["deltat"])
    assert st[0] == [0, round(60 * 3600 / dt)]
    assert float(a["PARM03"]["pchkptfreq"]) == 60 * 3600.0


def test_obcs_west_north_south(cases):
    d, _ = cases["F1"]
    o = _nml(d, "data.obcs")["OBCS_PARM01"]
    assert o["ob_iwest"] == "300*1" and o["ob_jsouth"] == "500*1" and o["ob_jnorth"] == "500*-1"
    assert "ob_ieast" not in o
    for s in ("west", "north", "south"):
        assert o[f"useorlanski{s}"] == ".TRUE."
    s3 = _nml(d, "data.obcs")["OBCS_PARM03"]
    assert int(s3["spongethickness"]) == CFG["n_sponge"]


def test_f5_cross_axis(cases):
    _, m1 = cases["F1"]
    _, m5 = cases["F5"]
    assert m5["actual"]["wind"]["angle_deg"] - m1["actual"]["wind"]["angle_deg"] == pytest.approx(90.0)
    assert m5["actual"]["profile"] == "vmf_summer"


def test_f1c_thesis_closure(cases):
    d, _ = cases["F1c"]
    p1 = _nml(d)["PARM01"]
    pk = _nml(d, "data.pkg")["PACKAGES"]
    assert pk["usekpp"] == ".FALSE." and pk.get("useggl90", ".FALSE.") == ".FALSE."
    assert float(p1["viscar"]) == 1e-5 and float(p1["diffkrt"]) == 1e-5
    assert float(p1["ivdc_kappa"]) == 1.0                       # thesis data.stage1
    assert not (d / "data.kpp").exists()
    a = _nml(cases["F1"][0])["PARM01"]
    diff = {k for k in set(a) | set(p1) if a.get(k) != p1.get(k)}
    assert diff == {"viscar", "ivdc_kappa"}                      # only the vertical closure


def test_f1g_ggl90(cases):
    d, _ = cases["F1g"]
    pk = _nml(d, "data.pkg")["PACKAGES"]
    assert pk["useggl90"] == ".TRUE." and pk["usekpp"] == ".FALSE."
    assert (d / "data.ggl90").exists() and not (d / "data.kpp").exists()
    assert _nml(d)["PARM01"] == _nml(cases["F1"][0])["PARM01"]


def test_v3_tide_only(cases):
    d, meta = cases["V3"]
    p5 = _nml(d)["PARM05"]
    assert "zonalwindfile" not in p5 and "meridwindfile" not in p5
    o = _nml(d, "data.obcs")["OBCS_PARM01"]
    assert o["useobcstides"] == ".TRUE."
    assert float(o["obcs_tidalperiod"].split(",")[0]) == pytest.approx(44714.16)
    tw = np.fromfile(d / "OBW_uTideSn.bin", ">f8")
    tn = np.fromfile(d / "OBN_vTideSn.bin", ">f8")
    ts = np.fromfile(d / "OBS_vTideSn.bin", ">f8")
    assert tw.size == 300 and tn.size == 500 and ts.size == 500
    u0 = meta["actual"]["tide"]["U0_m_s-1"]
    assert u0 > 0
    assert set(np.unique(tw)) <= {0.0, u0} and set(np.unique(ts)) <= {0.0, u0}
    assert set(np.unique(tn)) <= {0.0, -u0}                      # inward at the north edge
    t = meta["actual"]["tide"]
    # mass balance: wet OB face area x U0 = omega x eta0 x wet surface area
    assert t["U0_m_s-1"] * t["ob_face_area_m2"] == pytest.approx(
        2 * np.pi / 44714.16 * t["eta0_m"] * t["wet_area_m2"], rel=1e-12)


def test_v4_400m_grid(cases):
    d, meta = cases["V4"]
    p4 = _nml(d)["PARM04"]
    assert p4["delx"] == "250*400.0" and p4["dely"] == "150*400.0"
    assert np.fromfile(d / "bathy.bin", ">f8").size == 250 * 150
    cp.check_size(d, PKG / "config" / "fjord_vmf_400" / "code" / "SIZE.h")
    a = _nml(cases["F1"][0])["PARM01"]
    assert _nml(d)["PARM01"] == a


def test_fjord_configs_share_code_but_size_h():
    base = PKG / "config" / "fjord_vmf" / "code"
    for other in ("fjord_vmf_400", "fjord_vmf_100", "fjord_vmf_100s"):
        o = PKG / "config" / other / "code"
        names = sorted(p.name for p in base.iterdir())
        assert names == sorted(p.name for p in o.iterdir())
        for nme in names:
            if nme != "SIZE.h":
                assert (base / nme).read_text() == (o / nme).read_text(), (other, nme)


def test_namelists_name_only_written_inputs(cases):
    for n, (d, _) in cases.items():
        for f in d.iterdir():
            if f.suffix == ".bin" or f.name == "case.json":
                continue
            for line in f.read_text().splitlines():
                if line.strip().startswith(("#", "!")):
                    continue
                for ref in __import__("re").findall(r"'([^']+\.bin)'", line):
                    assert ref in mc.FJORD_BINS, (n, f.name, ref)


def test_run_campaign_dry_run_fjord(capsys):
    spec = importlib.util.spec_from_file_location("rc4", PKG / "scripts" / "run_campaign.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    assert m.main(["--config", "fjord_vmf", "--cases", "F*", "--dry-run"]) == 0
    assert capsys.readouterr().out.split() == ["F1", "F2", "F3", "F4", "F5", "F1c", "F1g"]


@pytest.mark.parametrize("name,cfg", [("F2_100", "fjord_vmf_100"), ("F2_100s", "fjord_vmf_100s")])
def test_100m_timing_cases(tmp_path, name, cfg):
    """R13 timing cases: the 100 m grid is the 200 m bathymetry split 2 x 2 (full box) or a cut
    of it (sub-box); SIZE.h matches; the ice mask is the same polygon."""
    meta = mc.make_case(name, tmp_path, overrides={"forced_h": 2.0, "t_end_h": 1.0},
                        write_bins=("bathy.bin",))
    cp.check_size(tmp_path, PKG / "config" / cfg / "code" / "SIZE.h")
    b = np.fromfile(tmp_path / "bathy.bin", ">f8").reshape(meta["params"]["ny"], meta["params"]["nx"])
    b200 = np.fromfile(DATA / "derived" / "bathy_fjord_vmf.bin", ">f8").reshape(300, 500)
    full = np.kron(b200, np.ones((2, 2)))
    if name == "F2_100":
        assert np.array_equal(b, full)
    else:
        x0, _, y0, _ = SPEC["configs"][cfg]["bbox_utm"]
        i0, j0 = (x0 - 455000) // 100, (y0 - 8594000) // 100
        assert np.array_equal(b, full[j0:j0 + 380, i0:i0 + 730])
    assert meta["actual"]["ice"]["mode"] == "landfast"
    assert meta["actual"]["stages"] == [[0, 120]]
