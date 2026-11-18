import json
import re

import numpy as np
import pytest

from exp6_southern_carbon.gendata import fields, make_case
from exp6_southern_carbon.gendata.grid import Grid

CO2 = {y: 285.0 + 0.5 * (y - 1850) for y in range(1850, 2101)}
G50 = Grid(20, 40, 50e3)


def _nml(text, key):
    m = re.search(rf"^\s*{key}\s*=\s*([^,\n]+)", text, re.M)
    return m.group(1).strip().strip("'") if m else None


def test_transient_branch(tmp_path):
    out = make_case.write_case("GT_w2060", tmp_path, {"co2": CO2, "iter_1850": 1000, "nIter0": 1000 + 210 * 15768})
    w = np.fromfile(out / "zonal_wind.bin", ">f4").reshape(12, 40, 20)   # periodicExternalForcing: 12 records
    for m in range(12):
        np.testing.assert_allclose(w[m], (1.2 * fields.taux(G50)).astype(">f4"))
    assert "gchem_secondsPerYear=31536000." in (out / "data.gchem").read_text().replace(" ", "")
    d = (out / "data.dic").read_text()
    assert (_nml(d, "dic_int1"), _nml(d, "dic_int2"), _nml(d, "dic_int3"), _nml(d, "dic_int4")) == ("2", "253", "1000", "15768")
    assert (out / "co2atmos.dat").read_text().count("\n") == 253
    data = (out / "data").read_text()
    assert int(_nml(data, "nIter0")) == 1000 + 210 * 15768 and int(_nml(data, "nTimeSteps")) == 25 * 15768
    assert float(_nml(data, "deltaT")) == 2000.0 and _nml(data, "deltaTtracer") is None
    c = json.loads((out / "case.json").read_text())
    assert c["branch_year"] == 2060 and c["wind_factor"] == 1.2


def test_control_fixed_1850(tmp_path):
    out = make_case.write_case("GC", tmp_path, {"co2": CO2, "iter_1850": 1000, "nIter0": 1000})
    d = (out / "data.dic").read_text()
    assert _nml(d, "dic_int1") == "1" and float(_nml(d, "dic_pCO2")) == pytest.approx(285e-6)
    assert not (out / "co2atmos.dat").exists()


def test_fixed_year_fine(tmp_path):
    out = make_case.write_case("FT_2060_c1", tmp_path, {"co2": CO2, "nIter0": 0})
    d = (out / "data.dic").read_text()
    assert float(_nml(d, "dic_pCO2")) == pytest.approx(CO2[2060] * 1e-6)
    assert "useGMRedi=.FALSE." in (out / "data.pkg").read_text().replace(" ", "")


def test_all_inputs_present_and_sized(tmp_path):
    out = make_case.write_case("G0", tmp_path, {"co2": CO2, "n_timesteps": 3942})
    sizes = {"bathy.bin": 800, "zonal_wind.bin": 12 * 800, "sst_relax.bin": 12 * 800, "t_relax_mask.bin": 49 * 800,
             "temperature.bin": 49 * 800, "par.bin": 12 * 800, "dic_wind.bin": 12 * 800, "fice.bin": 12 * 800,
             "silica.bin": 12 * 800}
    for i in range(1, 7):
        sizes[f"init_ptr{i}.bin"] = 49 * 800
    for i in range(1, 4):
        sizes[f"relax_ptr{i}.bin"] = 49 * 800
    for f, n in sizes.items():
        assert (out / f).stat().st_size == 4 * n, f
    for f in ("data", "data.pkg", "data.gmredi", "data.rbcs", "data.ptracers", "data.gchem", "data.dic",
              "data.diagnostics", "data.layers", "eedata"):
        assert (out / f).exists(), f
    assert "debugLevel=-1," in (out / "data").read_text().replace(" ", "")   # no per-step log lines
    wind = np.fromfile(out / "dic_wind.bin", ">f4")
    assert (wind == 10.0).all()
    po4, pre = (np.fromfile(out / f"init_ptr{i}.bin", ">f4") for i in (3, 6))
    np.testing.assert_array_equal(po4, pre)


def test_unknown_placeholder_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(make_case, "template_values", lambda p, g: {})
    with pytest.raises(KeyError):
        make_case.write_case("G0", tmp_path, {"co2": CO2, "n_timesteps": 10})


def test_g0_needs_length(tmp_path):
    with pytest.raises(ValueError):
        make_case.write_case("G0", tmp_path, {"co2": CO2})


def test_g0_accelerated_async(tmp_path):
    out = make_case.write_case("G0", tmp_path, {"co2": CO2, "n_timesteps": 3942})
    data = (out / "data").read_text()
    assert float(_nml(data, "deltaTMom")) == 2000.0
    assert float(_nml(data, "deltaTtracer")) == 8000.0 and float(_nml(data, "deltaTClock")) == 8000.0
    assert _nml(data, "deltaT") is None
    c = json.loads((out / "case.json").read_text())
    assert c["steps_per_year"] == 3942 and c["config"] == "g50p5"


def test_g0s_synchronous_after_g0(tmp_path):
    out = make_case.write_case("G0s", tmp_path, {"co2": CO2, "nIter0": 3942 * 1000})
    data = (out / "data").read_text()
    assert float(_nml(data, "deltaT")) == 2000.0 and int(_nml(data, "nTimeSteps")) == 48 * 15768
