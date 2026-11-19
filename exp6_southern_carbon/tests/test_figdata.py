"""analysis.figdata (plan Task 10): pure helpers behind figures/{spinup,validation,figdata}.json."""
import numpy as np
import pytest

from exp6_southern_carbon.analysis import figdata as fd


def test_calendar_years_from_end_iterations():
    spy, i1850 = 31536, 1_000_000
    ends = [i1850 + k * spy for k in (1, 2, 51)]
    assert fd.calendar_years(ends, i1850, spy) == [1850, 1851, 1900]


def test_year_index_of_calendar_year():
    assert fd.year_index(1850) == 1 and fd.year_index(1900) == 51


def test_xmean_section_nan_where_row_is_dry():
    a = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4)
    h = np.ones_like(a)
    h[1, 2, :] = 0.0
    h[0, 0, :2] = 0.0
    s = fd.xmean_section(a, h)
    assert s.shape == (2, 3) and np.isnan(s[1, 2])
    assert s[0, 0] == pytest.approx(a[0, 0, 2:].mean())


def test_g200_by_row_matches_region_gradient():
    z = np.array([5.0, 100.0, 300.0])
    col = np.array([2.0, 2.05, 2.1])                    # mol m-3
    dic = np.broadcast_to(col[:, None, None], (3, 2, 2)).copy()
    dic[:, 1] += 0.01 * np.arange(3)[:, None]
    h = np.ones_like(dic)
    g = fd.g200_by_row(dic, z, h)
    want0 = (np.interp(200.0, z, col) - col[0]) * 1e6 / 1035.0 / 200.0
    assert g[0] == pytest.approx(want0) and g[1] > g[0]


def test_v6b_passes_only_if_every_year_below_one_percent():
    ok = {"surf": 100.0, "residual": 0.5}
    bad = {"surf": 100.0, "residual": -2.0}
    v = fd.v6b({"GC": {1900: ok, 2060: ok}, "GT": {1900: ok, 2060: ok}})
    assert v["pass"] and v["max_frac"] == pytest.approx(5e-3)
    v = fd.v6b({"GC": {1900: ok}, "GT": {1900: bad}})
    assert not v["pass"] and v["frac"]["GT"]["1900"] == pytest.approx(0.02)


def test_nan_to_none_is_json_ready():
    import json
    out = fd._nan_to_none({1900: np.array([1.0, np.nan]), "b": np.float64(2.0), "c": np.bool_(True)})
    assert json.loads(json.dumps(out)) == {"1900": [1.0, None], "b": 2.0, "c": True}
