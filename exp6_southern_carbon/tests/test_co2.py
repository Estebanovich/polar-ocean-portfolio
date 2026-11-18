from pathlib import Path

import numpy as np
import pytest

from exp6_southern_carbon.gendata import co2

I4 = 31536
START = 1000
DATA = Path(__file__).resolve().parents[2] / "literatura" / "data" / "cmip6_co2"


def _series():
    return {y: float(y - 1850) for y in range(1850, 2101)}


def _file_array(tmp_path):
    out = tmp_path / "co2atmos.dat"
    n = co2.write_co2atmos(_series(), out)
    vals = np.loadtxt(out)
    assert vals.size == n
    return vals * 1e6


def test_mid_year_value_is_that_year(tmp_path):
    arr = _file_array(tmp_path)
    i3 = co2.co2_offset_iter(START, I4)
    for y in (1850, 1851, 1900, 2060, 2100):
        it = START + int((y - 1850 + 0.5) * I4)
        assert co2.atm_pco2_at(it, arr, i3, I4) == pytest.approx(y - 1850, abs=1e-9)


def test_whole_run_stays_inside_array(tmp_path):
    arr = _file_array(tmp_path)
    i3 = co2.co2_offset_iter(START, I4)
    for it in (START, START + 1, START + 251 * I4 - 1, START + 251 * I4):
        co2.atm_pco2_at(it, arr, i3, I4)  # raises IndexError if it reads past the file


def test_file_pads_both_ends(tmp_path):
    out = tmp_path / "co2atmos.dat"
    n = co2.write_co2atmos({y: 280.0 + (y - 1850) for y in range(1850, 2101)}, out)
    v = np.loadtxt(out)
    assert n == 253 and v.size == 253
    assert v[0] == v[1] == pytest.approx(280e-6) and v[-1] == v[-2] == pytest.approx(530e-6)


def test_dic_namelist_values():
    assert co2.dic_namelist(START, I4) == {"dic_int1": 2, "dic_int2": 253, "dic_int3": START, "dic_int4": I4}


@pytest.mark.skipif(not DATA.exists() or not any(DATA.iterdir()), reason="CMIP6 CO2 file not downloaded (D1)")
def test_load_cmip6():
    c = co2.load_cmip6_co2(DATA)
    assert sorted(c) == list(range(1850, 2101))
    assert 284.0 <= c[1850] <= 286.0 and c[2100] > 800.0
