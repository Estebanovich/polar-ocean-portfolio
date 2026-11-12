"""Tests for exp4_arctic_bayiw.gendata.ice (Task 4) and the ice switch of make_case.

Review Focus 3: the ice is only inside the inner-fjord polygon and on wet cells; outside it
the wind acts on open water (the shelf, Bellsund, the open boundaries and the sponge).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from exp4_arctic_bayiw.gendata import bathy, ice
from exp4_arctic_bayiw.gendata import make_case as mc

PKG = Path(__file__).resolve().parents[1]
DATA = PKG / "data_external"
SPEC = yaml.safe_load((PKG / "cases.yaml").read_text())
CFG = SPEC["configs"]["fjord_vmf"]
have_data = pytest.mark.skipif(not (DATA / CFG["bathy"]["ibcao_dir"]).exists(),
                               reason="IBCAO / NPI data not present in data_external/")


def test_mask_synthetic_polygon_and_wet():
    xc = np.arange(10) * 100.0 + 50.0
    yc = np.arange(8) * 100.0 + 50.0
    wet = np.ones((8, 10), bool)
    wet[3, 4] = False
    poly = [[300, 200], [800, 200], [800, 600], [300, 600]]
    m = ice.mask(xc, yc, wet, poly)
    assert m.shape == (8, 10) and m.dtype == bool
    XC, YC = np.meshgrid(xc, yc)
    inside = (XC > 300) & (XC < 800) & (YC > 200) & (YC < 600)
    assert np.array_equal(m, inside & wet)
    assert not m[3, 4]                                         # land inside the polygon


def test_fields_files_and_values():
    m = np.zeros((4, 5), bool)
    m[1:3, 2:4] = True
    f = ice.fields(m, h=1.0, A=1.0)
    assert set(f) == {"HEFF.bin", "AREA.bin"}
    for k in f:
        assert f[k].shape == (4, 5) and f[k].dtype == np.float64
    assert np.array_equal(f["HEFF.bin"], m.astype(float))
    assert np.array_equal(f["AREA.bin"], m.astype(float))
    f2 = ice.fields(m, h=0.5, A=0.8)
    assert f2["HEFF.bin"].max() == 0.5 and f2["AREA.bin"].max() == 0.8
    with pytest.raises(ValueError):
        ice.fields(m, A=1.2)                                   # concentration must be in [0, 1]


@pytest.fixture(scope="module")
def grid():
    if not (DATA / CFG["bathy"]["ibcao_dir"]).exists():
        pytest.skip("data_external not present")
    _, g = bathy.build(CFG, DATA)
    return g


@have_data
def test_real_mask_inner_fjord_only(grid):
    """The ice polygon (the profile polygon of cases.yaml, inner fjord east of Akseloya) holds
    ice only on wet cells of the inner fjord and of the two sound closure disks: no ice on the
    shelf / Bellsund (sea west of the sounds), none at an open boundary or in the sponge."""
    poly = CFG[CFG["ice"]["polygon"]]["polygon_utm"] if isinstance(CFG["ice"]["polygon"], str) \
        else CFG["ice"]["polygon"]
    m = ice.mask(grid["xc"], grid["yc"], grid["wet"], poly)
    inner, outer, disks = bathy.inner_fjord_mask(grid)
    both = np.logical_or.reduce(list(disks.values()))
    assert m.sum() > 0.99 * inner.sum()
    assert not np.any(m & ~grid["wet"])
    assert not np.any(m & outer)                               # Bellsund and the shelf
    assert not np.any(m & ~(inner | both))
    n = CFG["n_sponge"]
    band = np.zeros_like(m)
    band[:, :n + 1] = True                                     # W
    band[:n + 1, :] = True                                     # S
    band[-n - 1:, :] = True                                    # N
    assert not np.any(m & band)
    assert grid["xc"][np.nonzero(m.any(axis=0))[0].min()] > 488000.0


@pytest.fixture(scope="module")
def f2_f3(tmp_path_factory):
    if not (DATA / CFG["bathy"]["ibcao_dir"]).exists():
        pytest.skip("data_external not present")
    root = tmp_path_factory.mktemp("ice_cases")
    small = ("bathy.bin", "HEFF.bin", "AREA.bin")
    return {n: (root / n, mc.make_case(n, root / n, write_bins=small)) for n in ("F2", "F3", "F4")}


@have_data
def test_make_case_f2_landfast_ice_dynamics_off(f2_f3):
    d, meta = f2_f3["F2"]
    pkg = mc.parse_namelist((d / "data.pkg").read_text())["PACKAGES"]
    assert pkg["useseaice"] == ".TRUE."
    si = mc.parse_namelist((d / "data.seaice").read_text())["SEAICE_PARM01"]
    assert si["seaiceusedynamics"] == ".FALSE."
    assert si.get("seaiceusefreedrift", ".FALSE.") == ".FALSE."
    assert si["usepw79thermodynamics"] == ".FALSE."            # no thermodynamics
    assert float(si["seaice_waterdrag"]) == 5.5e-3              # PD0 section 8
    assert float(si["seaice_drag"]) == float(si["ocean_drag"])  # wind on ice = prescribed stress
    assert si["hefffile"] == "'HEFF.bin'" and si["areafile"] == "'AREA.bin'"
    a = np.fromfile(d / "AREA.bin", ">f8")
    h = np.fromfile(d / "HEFF.bin", ">f8")
    assert set(np.unique(a)) == {0.0, 1.0} and set(np.unique(h)) == {0.0, 1.0}
    assert meta["actual"]["ice"]["mode"] == "landfast"
    assert meta["actual"]["ice"]["cells"] == int(a.sum())
    assert meta["actual"]["profile"] == "vmf_winter"
    assert "pickup_seaice" in meta["params"]["pickup_stems"]


@have_data
def test_make_case_f3_has_no_ice(f2_f3):
    d, meta = f2_f3["F3"]
    pkg = mc.parse_namelist((d / "data.pkg").read_text())["PACKAGES"]
    assert pkg.get("useseaice", ".FALSE.") == ".FALSE."
    assert not (d / "data.seaice").exists()
    assert not (d / "HEFF.bin").exists() and not (d / "AREA.bin").exists()
    assert meta["actual"]["ice"] is None
    assert meta["actual"]["profile"] == "vmf_winter"


@have_data
def test_f2_f3_differ_only_by_the_ice(f2_f3):
    """pkg/seaice needs tracForcingOutAB = 1 (pkg/seaice/seaice_check.F:82-86); every fjord case
    sets it, so the data namelists of F2 and F3 are identical."""
    a = mc.parse_namelist((f2_f3["F2"][0] / "data").read_text())
    b = mc.parse_namelist((f2_f3["F3"][0] / "data").read_text())
    assert a == b
    assert a["PARM03"]["tracforcingoutab"] == "1"


@have_data
def test_make_case_f4_free_drift(f2_f3):
    d, meta = f2_f3["F4"]
    si = mc.parse_namelist((d / "data.seaice").read_text())["SEAICE_PARM01"]
    assert si["seaiceusedynamics"] == ".TRUE." and si["seaiceusefreedrift"] == ".TRUE."
    assert si["usepw79thermodynamics"] == ".FALSE."
    assert meta["actual"]["ice"]["mode"] == "free_drift"


@have_data
def test_f4_ice_equals_f2_ice(f2_f3):
    """Ruling R15: F4 (free drift) starts from the same ice as F2 (same polygon, h and A)."""
    d2, m2 = f2_f3["F2"]
    d4, m4 = f2_f3["F4"]
    for f in ("HEFF.bin", "AREA.bin"):
        assert (d2 / f).read_bytes() == (d4 / f).read_bytes()
    for k in ("cells", "area_km2", "h_m", "A", "polygon"):
        assert m2["actual"]["ice"][k] == m4["actual"]["ice"][k]
