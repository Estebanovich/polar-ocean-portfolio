"""Tests for exp4_arctic_bayiw.gendata.make_case and scripts/run_campaign.py (Task 3).

The bay_thesis cases read the read-only thesis inputs (BayIW_Rectan/input); those tests are
skipped if the volume is not mounted. MITgcm is never run: run_campaign gets a fake runner.
"""
from __future__ import annotations

import importlib.util
import json
import os
import re
from pathlib import Path

import numpy as np
import pytest
import yaml

from exp4_arctic_bayiw.gendata import make_case as mc
from exp4_arctic_bayiw.gendata import profiles, vgrid

PKG = Path(__file__).resolve().parents[1]
SPEC = yaml.safe_load((PKG / "cases.yaml").read_text())
THESIS = Path(SPEC["configs"]["bay_thesis"]["thesis_dir"])
THESIS_IN = THESIS / "input"
TCOPY = PKG / "thesis_copy" / "BayIW_Rectan"
NOTES = PKG / "notes" / "pd1_puente.md"
real = pytest.mark.skipif(not (THESIS_IN / "bahia_01_expand_bat.bin").exists(),
                          reason="read-only thesis inputs (BayIW_Rectan/input) not mounted")
NAMELISTS = ("data", "data.free", "data.pkg", "data.obcs", "data.diagnostics", "eedata", "data.kpp")


@pytest.fixture(scope="module")
def cases(tmp_path_factory):
    """R1 with every binary; R1k, R2, R3 with the small binaries only (namelists, bathy)."""
    if not (THESIS_IN / "bahia_01_expand_bat.bin").exists():
        pytest.skip("thesis inputs not mounted")
    root = tmp_path_factory.mktemp("cases")
    out = {}
    out["R1"] = (root / "R1", mc.make_case("R1", root / "R1"))
    out["R2"] = (root / "R2", mc.make_case("R2", root / "R2", write_bins=("bathy.bin", "T.init.bin",
                                                                            "S.init.bin")))
    for n in ("R1k", "R3"):
        out[n] = (root / n, mc.make_case(n, root / n, write_bins=("bathy.bin",)))
    return out


def _nml(d: Path) -> dict:
    return {f: mc.parse_namelist((d / f).read_text()) for f in NAMELISTS if (d / f).exists()}


def _thesis_column(case_dir: Path, name: str, nz=50, ny=352, nx=560):
    a = np.fromfile(case_dir / name, ">f8").reshape(nz, ny, nx)
    assert np.ptp(a, axis=(1, 2)).max() == 0.0          # horizontally uniform
    return a[:, 0, 0]


# ---- namelist helpers (synthetic) -----------------------------------------------------------

TXT = """# comment
 &PARM01
 viscAr=1.E-5,
 tRef=1., 2.,
\t3.,

 sRef= 4.,
# f0=1.,
 f0=6.97E-5,
 &
 &PARM03
 nIter0=0,
 &
"""


def test_parse_namelist_multiline_values_and_comments():
    p = mc.parse_namelist(TXT)
    assert p["PARM01"] == {"viscar": "1.E-5", "tref": "1., 2., 3.", "sref": "4.", "f0": "6.97E-5"}
    assert p["PARM03"] == {"niter0": "0"}


def test_set_nml_key_replace_insert_remove():
    t = mc.set_nml_key(TXT, "PARM01", "tRef", "2*0.5")
    assert mc.parse_namelist(t)["PARM01"]["tref"] == "2*0.5"
    assert "# f0=1.," in t                                       # comments kept
    t = mc.set_nml_key(t, "PARM01", "beta", "0.")
    assert mc.parse_namelist(t)["PARM01"]["beta"] == "0."
    t = mc.set_nml_key(t, "PARM01", "viscAr", None)
    p = mc.parse_namelist(t)
    assert "viscar" not in p["PARM01"] and p["PARM01"]["sref"] == "4."
    assert p["PARM03"] == {"niter0": "0"}
    with pytest.raises(KeyError):
        mc.set_nml_key(TXT, "PARM09", "x", "1")


def test_thesis_wind_pulse_shape_synthetic():
    tau = mc.thesis_wind_patch(nx=40, ny=30, cx=20, cy=15, radius=6)
    assert tau.shape == (180, 30, 40)
    t = np.linspace(0, 60, 180)
    k = np.argmax(tau[:, 15, 20])
    assert t[k] <= 2.5 and np.all(tau[t > 2.5] == 0.0)               # half Gaussian
    assert tau[:, 15, 20].max() <= 0.10 and tau[:, 0, 0].max() == 0.0  # patch cut at the radius


# ---- R1 = thesis (Task 3 Step 1) --------------------------------------------------------------

def _thesis_f0():
    txt = (TCOPY / "run_expand" / "data.stage1").read_text()
    return mc.parse_namelist(txt)["PARM01"]["f0"]


@real
def test_r1_f0_is_the_thesis_data_stage1(cases):
    d, meta = cases["R1"]
    p = _nml(d)
    assert p["data"]["PARM01"]["f0"] == _thesis_f0() == "6.97E-5"
    assert p["data.free"]["PARM01"]["f0"] == "6.97E-5"
    assert "beta" not in p["data"]["PARM01"]                  # MITgcm default 1e-11, as the thesis
    a = meta["actual"]
    assert a["f0_s-1"] == 6.97e-5 and a["beta_m-1s-1"] == 1.0e-11
    # beta-plane: f = f0 + beta yC (model/src/ini_cori.F:77) is 5.36e-5 at the south edge
    assert a["f_domain_s-1"][0] == pytest.approx(5.3556e-5, rel=1e-4)


@real
def test_r1_bathymetry_byte_identical_to_thesis(cases):
    d, _ = cases["R1"]
    assert (d / "bathy.bin").read_bytes() == (THESIS_IN / "bahia_01_expand_bat.bin").read_bytes()
    assert (d / "delX.bin").read_bytes() == (THESIS_IN / "bahia_01_expand_dx.bin").read_bytes()
    assert (d / "delY.bin").read_bytes() == (THESIS_IN / "bahia_01_expand_dy.bin").read_bytes()
    assert (d / "T.init.bin").read_bytes() == (THESIS_IN / "feb_temp_50zlev_560x352.bin").read_bytes()
    assert (d / "S.init.bin").read_bytes() == (THESIS_IN / "feb_salt_50zlev_560x352.bin").read_bytes()


def _thesis_wind_function():
    os.environ.setdefault("MPLBACKEND", "Agg")
    spec = importlib.util.spec_from_file_location("thesis_wind_patch",
                                                  TCOPY / "input" / "wind_patch.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.make_meridional_wind_patch


@real
def test_r1_wind_equals_thesis_function(cases):
    d, meta = cases["R1"]
    w = SPEC["configs"]["bay_thesis"]["wind"]
    tau, _ = _thesis_wind_function()(nt=w["nt"], end_time=w["end_time_h"], max_tau=w["max_tau"],
                                     sigma_t=w["sigma_t_h"], center_t=w["center_t_h"], nx=560,
                                     ny=352, cx=w["cx"], cy=w["cy"], radius=w["radius"],
                                     write_bin=False)
    ref = tau.transpose(2, 1, 0)                              # (nt, ny, nx), as the thesis wrote it
    got = np.fromfile(d / "tauy.bin", ">f8").reshape(ref.shape)
    assert np.max(np.abs(got - ref)) < 1e-15
    assert meta["actual"]["wind"]["max_tau_N_m-2"] == pytest.approx(ref.max())


@real
def test_r1_wind_byte_identical_to_stored_thesis_file(cases):
    d, _ = cases["R1"]
    f = THESIS_IN / "merid_local_bay_560x352.bin"
    assert (d / "tauy.bin").read_bytes() == f.read_bytes()


# keys of R1 that differ from the thesis namelists: I/O only (MNC -> MDS, fixed file names,
# startTime/endTime -> nIter0/nTimeSteps with the same iterations)
R1_VS_THESIS = {
    ("data", "PARM03", "starttime"), ("data", "PARM03", "endtime"), ("data", "PARM03", "ntimesteps"),
    ("data", "PARM03", "dumpfreq"), ("data", "PARM04", "delyfile"), ("data", "PARM04", "delxfile"),
    ("data", "PARM05", "bathyfile"), ("data", "PARM05", "meridwindfile"),
    ("data", "PARM05", "hydrogthetafile"), ("data", "PARM05", "hydrogsaltfile"),
    ("data.free", "PARM03", "starttime"), ("data.free", "PARM03", "endtime"),
    ("data.free", "PARM03", "ntimesteps"), ("data.free", "PARM03", "dumpfreq"),
    ("data.free", "PARM04", "delyfile"), ("data.free", "PARM04", "delxfile"),
    ("data.free", "PARM05", "bathyfile"),
    ("data.pkg", "PACKAGES", "usemnc"),
}
R1_DIAG_ADDED = {"fields(1:3,2)", "filename(2)", "frequency(2)", "timephase(2)",
                 "fields(1,3)", "filename(3)", "frequency(3)", "timephase(3)"}


def _diff(a: dict, b: dict) -> dict:
    """{(file, group, key): (value in a, value in b)} over every key of either; '*' marks a
    group present in one only; None = absent."""
    out = {}
    for f in set(a) | set(b):
        ga, gb = a.get(f, {}), b.get(f, {})
        for g in set(ga) | set(gb):
            if (g in ga) != (g in gb):
                out[(f, g, "*")] = ("present" if g in ga else None, "present" if g in gb else None)
            ka, kb = ga.get(g, {}), gb.get(g, {})
            for k in set(ka) | set(kb):
                if ka.get(k) != kb.get(k):
                    out[(f, g, k)] = (ka.get(k), kb.get(k))
    return out


@real
def test_r1_namelists_differ_from_thesis_only_in_io(cases):
    d, _ = cases["R1"]
    r1 = _nml(d)
    th = {f: mc.parse_namelist((TCOPY / "run_expand" / s).read_text()) for f, s in
          (("data", "data.stage1"), ("data.free", "data.stage2"), ("data.pkg", "data.pkg"),
           ("data.obcs", "data.obcs"), ("data.diagnostics", "data.diagnostics"),
           ("eedata", "eedata"))}
    diff = _diff(th, r1)
    diag = {k for (f, g, k) in diff if f == "data.diagnostics"}
    assert diag == R1_DIAG_ADDED                                    # only new output streams
    assert {k for k in diff if k[0] != "data.diagnostics"} == R1_VS_THESIS
    assert r1["data"]["PARM03"]["ntimesteps"] == "7200" and r1["data"]["PARM03"]["niter0"] == "0"
    assert r1["data.free"]["PARM03"]["niter0"] == "7200"
    assert r1["data.free"]["PARM03"]["ntimesteps"] == "7200"
    assert "data.kpp" not in r1


@real
def test_r1_case_json_stage_plan_and_sources(cases):
    _, meta = cases["R1"]
    a = meta["actual"]
    assert a["stages"] == [[0, 7200], [7200, 7200]]
    assert a["profile"] == "thesis_realistic" and a["surface"] == "open water"
    assert a["sources"]["bathy.bin"]["sha256"] == \
        "b3201d2af900f3fdfba768835518b600e9eaee6acb8d97677395f3e02307f155"
    assert a["column_depth_m"] == pytest.approx(599.9) and a["bathy_max_depth_m"] == pytest.approx(1000.0)
    assert meta["params"]["nx"] == 560 and meta["params"]["ny"] == 352 and meta["params"]["nz"] == 50


# ---- R1k (Task 3 Step 7) -----------------------------------------------------------------------

def _notes_table(tag: str) -> dict:
    """Rows of the table between <!-- tag --> and <!-- /tag --> in notes/pd1_puente.md:
    {(file, group, key): (R1 value, other value)}; '—' = absent."""
    txt = NOTES.read_text()
    m = re.search(rf"<!-- {tag} -->(.*?)<!-- /{tag} -->", txt, re.S)
    assert m, f"table {tag} not found in {NOTES}"
    rows = {}
    for ln in m.group(1).splitlines():
        cells = [c.strip().strip("`") for c in ln.strip().strip("|").split("|")]
        if len(cells) < 5 or cells[0] in ("archivo", "") or set(cells[0]) <= {"-"}:
            continue
        val = [None if c == "—" else c for c in cells[3:5]]
        rows[(cells[0], cells[1], cells[2].lower())] = tuple(val)
    return rows


@real
def test_r1k_activates_kpp_and_differs_from_r1_only_in_listed_keys(cases):
    r1, r1k = _nml(cases["R1"][0]), _nml(cases["R1k"][0])
    assert r1k["data.pkg"]["PACKAGES"]["usekpp"] == ".TRUE."
    assert "data.kpp" in r1k
    listed = _notes_table("R1k-R1")
    assert _diff(r1, r1k) == listed


@real
def test_r1k_case_json(cases):
    a = cases["R1k"][1]["actual"]
    assert a["stages"] == [[0, 7200], [7200, 21600]]
    assert a["f0_s-1"] == 6.97e-5 and a["kpp"] is True and a["eos"]["tAlpha"] == 2.0e-4


# ---- R2 / R3 (bridge, generated only) ---------------------------------------------------------

def _summer_on_thesis_levels():
    pj = json.loads((PKG / "figures" / "profiles_vmf.json").read_text())
    zc = np.array(pj["z_c_m"])
    s = pj["summer"]
    delR = np.array(mc.parse_namelist_list((TCOPY / "run_expand" / "data.stage1").read_text(), "delR"))
    f = vgrid.faces(delR)
    tzc = -0.5 * (f[1:] + f[:-1])
    T = np.interp(-tzc, -zc, s["T"])
    S = np.interp(-tzc, -zc, s["S"])
    eos = SPEC["configs"]["fjord_vmf"]["eos"]
    S, _ = profiles.stabilize(tzc, T, S, eos)
    return T, S, pj["summer_on_thesis_levels"]


@real
def test_r2_f0_78n_and_summer_profile_on_thesis_levels(cases):
    d, meta = cases["R2"]
    p = _nml(d)
    for f in ("data", "data.free"):
        assert p[f]["PARM01"]["f0"] == "1.4266E-4" and p[f]["PARM01"]["beta"] == "0."
    T, S, ref = _summer_on_thesis_levels()
    Tm, Sm = _thesis_column(d, "T.init.bin"), _thesis_column(d, "S.init.bin")
    assert np.max(np.abs(Tm - T)) < 1e-12 and np.max(np.abs(Sm - S)) < 1e-12
    # and the same column as Task 2 (theory_numbers, figures/profiles_vmf.json)
    assert np.max(np.abs(Tm - np.array(ref["T"]))) < 1e-12
    assert np.max(np.abs(Sm - np.array(ref["S"]))) < 1e-12
    a = meta["actual"]
    assert a["f0_s-1"] == 1.4266e-4 and a["beta_m-1s-1"] == 0.0
    assert a["profile"] == "vmf_summer_thesis_levels"


@real
def test_r2_fitted_eos_kpp_and_pd0_mixing(cases):
    p = _nml(cases["R2"][0])["data"]["PARM01"]
    eos = SPEC["configs"]["fjord_vmf"]["eos"]
    assert float(p["talpha"]) == eos["tAlpha"] and float(p["sbeta"]) == eos["sBeta"]
    assert float(p["rhonil"]) == float(p["rhoconst"]) == eos["rhoNil"]
    assert p["tref"] == f"50*{eos['tRef']}" and p["sref"] == f"50*{eos['sRef']}"
    assert p["eostype"] == "'LINEAR'"
    assert p["viscar"] == "1.E-4" and p["diffkrt"] == p["diffkrs"] == "1.E-5"
    assert "ivdc_kappa" not in p and p["implicitviscosity"] == ".TRUE."
    assert "viscah" not in p and p["viscc2smag"] == "2.2" and p["viscc4smag"] == "1."
    assert p["diffkht"] == p["diffkhs"] == "0." and p["bottomdragquadratic"] == "2.5E-3"
    assert p["tempadvscheme"] == p["saltadvscheme"] == "77"
    assert _nml(cases["R2"][0])["data.pkg"]["PACKAGES"]["usekpp"] == ".TRUE."
    assert cases["R2"][1]["actual"]["stages"] == [[0, 7200], [7200, 26400]]


@real
def test_r3_is_r2_without_the_bay(cases):
    d3, _ = cases["R3"]
    assert (d3 / "bathy.bin").read_bytes() == (THESIS_IN / "nobahia_01_expand_bat.bin").read_bytes()
    assert _nml(d3) == _nml(cases["R2"][0])


# ---- errors ----------------------------------------------------------------------------------

@real
def test_unknown_case_and_edit_set(tmp_path):
    with pytest.raises(KeyError):
        mc.make_case("R9", tmp_path / "x")
    spec = yaml.safe_load((PKG / "cases.yaml").read_text())
    spec["cases"]["Rbad"] = {"config": "bay_thesis", "edits": ["nope"]}
    cf = tmp_path / "cases.yaml"
    cf.write_text(yaml.safe_dump(spec))
    with pytest.raises(KeyError, match="nope"):
        mc.make_case("Rbad", tmp_path / "y", cases_file=cf, write_bins=())


@real
def test_short_run_has_a_single_stage(tmp_path):
    meta = mc.make_case("R1", tmp_path / "s", overrides={"t_end_h": 2.0}, write_bins=())
    assert meta["actual"]["stages"] == [[0, 240]]
    assert not (tmp_path / "s" / "data.free").exists()
    assert mc.parse_namelist((tmp_path / "s" / "data").read_text())["PARM03"]["ntimesteps"] == "240"


# ---- run_campaign (fake runner) ----------------------------------------------------------------

def _load_rc():
    spec = importlib.util.spec_from_file_location("exp4_run_campaign", PKG / "scripts" / "run_campaign.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_run_campaign_dry_run_lists_cases(capsys):
    rc = _load_rc()
    assert rc.main(["--config", "bay_thesis_kpp", "--cases", "R*", "--dry-run"]) == 0
    assert capsys.readouterr().out.split() == ["R1b", "R1k", "R2", "R3"]


def test_run_campaign_np_must_match_size_h(tmp_path):
    rc = _load_rc()
    with pytest.raises(RuntimeError, match="nPx"):
        rc.run_campaign("bay_thesis", "R1", 5, runs_root=tmp_path, build_root=tmp_path,
                        sh=lambda c: 0, pgrep=lambda: "")


def test_run_campaign_two_stages_with_free_namelist(tmp_path, monkeypatch):
    """Generation stubbed (small files); the fake runner checks what each stage gets."""
    rc = _load_rc()

    def fake_make_case(name, out, **kw):
        out = Path(out)
        out.mkdir(parents=True, exist_ok=True)
        (out / "data").write_text(" &PARM03\n nIter0=0,\n nTimeSteps=4,\n deltaT=30.,\n"
                                  " pChkptFreq=120.,\n &\n")
        (out / "data.free").write_text(" &PARM03\n nIter0=4,\n nTimeSteps=8,\n deltaT=30.,\n"
                                       " chkptFreq=120.,\n &\n")
        meta = {"params": {"nx": 560, "ny": 352, "nz": 50, "t_end_h": 0.1, "deltaT": 30.0,
                           "forced_h": 60.0},
                "actual": {"stages": [[0, 4], [4, 8]]}}
        (out / "case.json").write_text(json.dumps(meta))
        return meta

    calls = []

    def fake_sh(cmd):
        cmd = [str(c) for c in cmd]
        if cmd[0].endswith("build_mitgcm.sh"):
            b = Path(cmd[2])
            b.mkdir(parents=True, exist_ok=True)
            (b / "mitgcmuv").write_bytes(b"exe")
            return 0
        _, build, inp, run, np_ = cmd
        inp, run = Path(inp), Path(run)
        d = (inp / "data").read_text()
        calls.append((run.name, d, sorted(p.name for p in inp.iterdir())))
        run.mkdir(parents=True, exist_ok=True)
        n0 = int(re.search(r"nIter0=(\d+)", d).group(1))
        n = int(re.search(r"nTimeSteps=(\d+)", d).group(1))
        if n0 == 0:
            (run / f"pickup.{n0 + n:010d}.data").write_bytes(b"p")
            (run / f"pickup.{n0 + n:010d}.meta").write_text("m")
            (run / f"pickup_orlanskiE.{n0 + n:010d}.001.001.data").write_bytes(b"o")
        (run / "run_status.json").write_text(json.dumps(
            {"state": "OK", "last_iter": n0 + n, "n_timesteps": n, "wall_s": 0}))
        return 0

    monkeypatch.setattr(rc, "make_case", fake_make_case)
    monkeypatch.setattr(rc.cp, "check_free_space", lambda *a, **k: None)
    res = rc.run_campaign("bay_thesis", "R1", 4, runs_root=tmp_path / "runs",
                          build_root=tmp_path / "b", sh=fake_sh, pgrep=lambda: "")
    assert res == {"R1": "OK"}
    assert [c[0] for c in calls] == ["R1", "stage2"]
    assert "pChkptFreq" in calls[0][1] and "chkptFreq=120." in calls[1][1]
    assert "nIter0=4" in calls[1][1] and "nTimeSteps=8" in calls[1][1]
    assert "pickup.0000000004.data" in calls[1][2]
    assert "pickup_orlanskiE.0000000004.001.001.data" in calls[1][2]      # OB pickups too
    fin = json.loads((tmp_path / "runs" / "R1" / "run_status_final.json").read_text())
    assert fin["state"] == "OK" and fin["last_iter"] == 12
    cj = json.loads((tmp_path / "runs" / "R1" / "input" / "case.json").read_text())
    assert len(cj["binary_sha256"]) == 64
