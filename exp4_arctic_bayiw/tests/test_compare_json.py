"""scripts/compare_json.py (Task 8): the committed figures/*.json against a regenerated set."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]


def _cj():
    spec = importlib.util.spec_from_file_location("exp4_compare_json", PKG / "scripts" / "compare_json.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_diff_lists_paths_and_relative_size():
    cj = _cj()
    a = {"x": 1.0, "y": {"z": [1.0, 2.0], "s": "a"}, "only_a": 1}
    b = {"x": 1.0, "y": {"z": [1.0, 2.0 * (1 + 1e-9)], "s": "b"}, "only_b": 2}
    d = {p: r for p, r in cj.diff(a, b)}
    assert set(d) == {"y.z[1]", "y.s", "only_a", "only_b"}
    assert 0.9e-9 < d["y.z[1]"] < 1.1e-9
    assert d["y.s"] is None and d["only_a"] is None
    assert cj.diff(a, json.loads(json.dumps(a))) == []


def test_compare_dirs(tmp_path):
    cj = _cj()
    old, new = tmp_path / "old", tmp_path / "new"
    old.mkdir(), new.mkdir()
    (old / "a.json").write_text(json.dumps({"v": 1.0}))
    (new / "a.json").write_text(json.dumps({"v": 1.0}))
    (old / "b.json").write_text(json.dumps({"v": 1.0}))
    (new / "b.json").write_text(json.dumps({"v": 1.5}))
    (old / "c.json").write_text("{}")
    rep = cj.compare_dirs(old, new)
    assert rep["identical"] == ["a.json"]
    assert rep["differ"]["b.json"]["n"] == 1 and rep["differ"]["b.json"]["max_rel"] == 0.5
    assert rep["missing_in_new"] == ["c.json"]
    assert cj.main([str(old), str(new)]) == 1
    assert cj.main([str(old / "a.json").rsplit("/", 1)[0], str(old)]) == 0
