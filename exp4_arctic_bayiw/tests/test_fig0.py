"""Figure 0 (Task 8): the schematic is drawn from committed inputs only (satellite images under
figures/src/sat, the case registry, bathy_report.json and validation.json), never from model output."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
SRC = PKG / "figures" / "src"


def _mod(name):
    spec = importlib.util.spec_from_file_location(f"exp4_{name}", SRC / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_places_only_the_committed_satellite_images():
    f0 = _mod("fig0_schematic")
    cmds = f0.build()
    placed = [Path(c["params"]["path"]) for c in cmds if c["command"] == "file.place"]
    assert {p.name for p in placed} == {f0.S2, f0.MODIS}
    for p in placed:
        assert p.parent == SRC / "sat" and p.exists() and p.stat().st_size < 1_000_000


def test_no_model_output_is_read():
    src = (SRC / "fig0_schematic.py").read_text()
    for bad in ("load.run", "open_run", "EXP4_RUNS", "/Volumes/LaCie", "partition", "hypotheses.json"):
        assert bad not in src, bad


def test_facts_match_the_fetched_images():
    f0, fs = _mod("fig0_schematic"), _mod("fetch_sat")
    assert f0.BBOX_DOMAIN == fs.BBOX_DOMAIN and f0.SIZE_DOMAIN == fs.SIZE_DOMAIN
    assert f0.BBOX_SVALBARD == fs.BBOX_SVALBARD
    assert f0.S2 in fs.SCENES
    assert "2024-03-25" in f0.S2 and f0.FACTS["sentinel2"]["date_text"] == "25 March 2024"
    assert "2023-07-05" in f0.MODIS and f0.FACTS["modis"]["date_text"] == "5 July 2023"
    committed = json.loads((PKG / "figures" / "fig0_schematic.json").read_text())
    assert committed == f0.FACTS


def test_domain_inside_the_image():
    import yaml
    f0 = _mod("fig0_schematic")
    x0, x1, y0, y1 = yaml.safe_load((PKG / "cases.yaml").read_text())["configs"]["fjord_vmf"]["bbox_utm"]
    X0, Y0, X1, Y1 = f0.BBOX_DOMAIN
    assert X0 < x0 and x1 < X1 and Y0 < y0 and y1 < Y1
