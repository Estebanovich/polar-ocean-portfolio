from pathlib import Path

import yaml

CASES = Path(__file__).resolve().parents[1] / "cases.yaml"
KEYS = {"config", "pickup_from", "branch_year", "wind_factor", "pco2_mode", "sponge", "n_years"}


def _load():
    return yaml.safe_load(CASES.read_text())


def test_configs():
    c = _load()["configs"]
    assert (c["g50"]["nx"], c["g50"]["ny"], c["g50"]["nr"], c["g50"]["dx_m"], c["g50"]["deltaT"], c["g50"]["np"]) == (
        20, 40, 49, 50000.0, 1000.0, 1)
    assert (c["f10"]["nx"], c["f10"]["ny"], c["f10"]["nr"], c["f10"]["dx_m"], c["f10"]["np"]) == (100, 200, 49, 10000.0, 5)


def test_every_case_has_keys_and_valid_parent():
    cases = _load()["cases"]
    for name, p in cases.items():
        assert KEYS <= set(p), name
        assert p["pickup_from"] is None or p["pickup_from"] in cases, name
        assert p["pco2_mode"] in {"fixed1850", "transient", "fixedYear"}, name
        assert p["sponge"] in {"pi", "anth"}, name


def test_branches():
    cases = _load()["cases"]
    for name, p in cases.items():
        if "_w" in name:
            assert p["wind_factor"] == 1.2, name
        if p["branch_year"] is not None:
            assert p["branch_year"] in {1900, 1980, 2020, 2060}, name
        if name.startswith("FT_") or name.startswith("F0_"):
            assert p["config"] == "f10", name
    assert cases["GC"]["wind_factor"] == cases["GT"]["wind_factor"] == 1.0
    for y in (1900, 1980, 2020, 2060):
        assert cases[f"GT_w{y}"]["pickup_from"] == "GT" and cases[f"GC_w{y}"]["pickup_from"] == "GC"
    for y in (1900, 2060):
        for m in (1, 2, 3):
            assert cases[f"FT_{y}_c{m}"]["pickup_from"] == f"F0_{y}"
            assert cases[f"FT_{y}_w{m}"]["wind_factor"] == 1.2
