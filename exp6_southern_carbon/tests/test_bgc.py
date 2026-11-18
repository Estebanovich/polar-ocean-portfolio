import numpy as np
import pytest

from exp6_southern_carbon.gendata import bgc
from exp6_southern_carbon.gendata.grid import Grid

TUT_REF = {  # tutorial_global_oce_biogeo/input/data.ptracers, first and last of 15 levels
    "DIC": (2.0282, 2.2758), "Alk": (2.3086, 2.3867), "PO4": (0.5438e-3, 2.2296e-3),
    "DOP": (0.0, 0.0), "O2": (0.2457, 0.2051)}


def test_ref_profiles_at_source_levels():
    zc = bgc.source_levels()
    assert zc.size == 15 and zc[0] == 25.0
    p = bgc.ref_profiles(zc)
    for k, (a, b) in TUT_REF.items():
        assert p[k][0] == pytest.approx(a, abs=1e-12) and p[k][-1] == pytest.approx(b, abs=1e-12)


def test_ref_profiles_on_model_grid():
    p = bgc.ref_profiles(Grid(20, 40, 50e3).z)
    assert set(p) == {"DIC", "Alk", "PO4", "DOP", "O2"}
    assert all(v.shape == (49,) for v in p.values())
    assert all((p[k] > 0).all() for k in ("DIC", "Alk", "PO4", "O2"))
    assert p["DIC"][0] == pytest.approx(2.0282)          # above the first source level: constant


def test_par_seasonal_cycle():
    par = bgc.par_monthly(lat_deg=-55.0, par_frac=0.4)
    assert par.shape == (12,) and (par >= 0).all()
    assert int(np.argmax(par)) in (11, 0) and int(np.argmin(par)) == 5
    assert 5.0 < par.min() < par.max() < 0.4 * 1361.0


def test_anth_sponge():
    z = np.array([0.0, 1000.0, 2000.0])
    pi = np.array([2.0, 2.2, 2.3])
    s = bgc.anth_sponge(pi, {y: 0.01 for y in range(1850, 2101)}, z, efold_m=1000.0)
    assert s.shape == (251, 3)
    assert s[0, 0] == pytest.approx(2.01) and s[0, 1] - 2.2 == pytest.approx(0.01 * np.exp(-1))


def test_silica_from_tutorial_south_of_50S():
    si = bgc.silica_surface()
    assert 0.0 < si < 0.2   # mol m-3
