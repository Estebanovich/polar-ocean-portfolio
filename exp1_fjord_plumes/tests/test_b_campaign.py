"""Tests for the config (b) campaign preparation (Task 6, R29): animation
diagnostics in the data.diagnostics template and the eta-retention check."""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pytest

from exp1_fjord_plumes.analysis.b_checks import RETENTION_MAX, eta_retention

TEMPLATE = (Path(__file__).resolve().parents[1] / "config" / "b_twoplume3d"
            / "input_templates" / "data.diagnostics")

# The seven PD1-approved lists, exactly as before R29 (must stay unchanged)
OLD_LISTS = {
    1: ("'THETA   '", "'THETA'", "-1800."),
    2: ("'SALT    '", "'SALT'", "-1800."),
    3: ("'UVEL    '", "'UVEL'", "-1800."),
    4: ("'VVEL    '", "'VVEL'", "-1800."),
    5: ("'WVEL    '", "'WVEL'", "-1800."),
    6: ("'ETAN    '", "'ETAN'", "-1800."),
    7: ("'ICFfwFlx','ICFhtFlx'", "'icefrontDiag'", "900."),
}


def _lists(text: str) -> dict:
    out: dict = {}
    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            continue
        m = re.match(r"\s*(fields|fileName|frequency|timePhase|levels)\(([^)]*)\)\s*=\s*(.*?),?\s*$", line)
        if m:
            key, idx, val = m.groups()
            n = int(idx.split(",")[-1])
            out.setdefault(n, {})[key] = val.strip().rstrip(",")
    return out


def test_diagnostics_template_old_lists_unchanged():
    lists = _lists(TEMPLATE.read_text())
    for n, (fld, fname, freq) in OLD_LISTS.items():
        assert lists[n]["fields"] == fld
        assert lists[n]["fileName"] == fname
        assert lists[n]["frequency"] == freq
        assert "levels" not in lists[n]


def test_diagnostics_template_animation_lists_R29():
    lists = _lists(TEMPLATE.read_text())
    assert set(lists) == set(range(1, 10))
    w, t = lists[8], lists[9]
    assert w["fields"] == "'WVEL    '" and w["fileName"] == "'animW'"
    assert w["frequency"] == "-120." and "levels" not in w  # full 3-D (y-z face section)
    assert t["fields"] == "'THETA   '" and t["fileName"] == "'animT'"
    assert t["frequency"] == "-120."
    assert [float(v) for v in t["levels"].split(",")] == [8.0, 16.0]


def _write_eta(run: Path, it: int, eta: np.ndarray, name: str = "ETAN") -> None:
    ny, nx = eta.shape
    stem = run / f"{name}.{it:010d}"
    np.ascontiguousarray(eta, dtype=">f8").tofile(f"{stem}.data")
    Path(f"{stem}.meta").write_text(
        f" nDims = [ 2 ];\n dimList = [\n {nx}, 1, {nx},\n {ny}, 1, {ny}\n ];\n"
        f" dataprec = [ 'float64' ];\n nrecords = [ 1 ];\n timeStepNumber = [ {it} ];\n")


def _make_run(tmp_path: Path, retained: float, state: str = "OK") -> Path:
    """Box 4 x 5 cells of 10 m, land row j = 0; Q = 2 m3/s; eta rises at
    retained * Q / A (uniform over the wet cells, 0 on land)."""
    run = tmp_path / "run"
    run.mkdir()
    ny, nx, dt = 4, 5, 2.0
    (run / "run_status.json").write_text(json.dumps(
        {"state": state, "last_iter": 5400, "n_timesteps": 5400, "wall_s": 1}))
    (run / "case.json").write_text(json.dumps(
        {"params": {"deltaT": dt, "nx": nx, "ny": ny, "dx": 10.0, "dy": 10.0},
         "Q0_total": 2.0}))
    depth = np.full((ny, nx), 300.0)
    depth[0] = 0.0
    np.ascontiguousarray(-depth, dtype=">f8").tofile(run / "topog.bin")
    area = 3 * 5 * 100.0
    for it in (900, 1800, 2700, 3600, 4500):
        _write_eta(run, it, np.where(depth > 0, retained * 2.0 * it * dt / area, 0.0))
    _write_eta(run, 5400, np.where(depth > 0, retained * 2.0 * 5400 * dt / area, 0.0), "Eta")
    return run


def test_eta_retention_known_fraction(tmp_path):
    r = eta_retention(_make_run(tmp_path, 0.3))
    assert r["t_end_h"] == pytest.approx(3.0)
    assert r["eta_trapped_m"] == pytest.approx(2.0 * 10800 / 1500.0)
    assert r["retained_cum"] == pytest.approx(0.3)
    assert r["retained_last_hour"] == pytest.approx(0.3)
    assert r["flag"] is False


def test_eta_retention_flags_above_threshold(tmp_path):
    r = eta_retention(_make_run(tmp_path, RETENTION_MAX + 0.1))
    assert r["flag"] is True


def test_eta_retention_refuses_failed_run(tmp_path):
    with pytest.raises(RuntimeError):
        eta_retention(_make_run(tmp_path, 0.3, state="FAILED"))


def test_animation_frames_prefer_R29_output(tmp_path):
    from exp1_fjord_plumes.analysis.animate_b import frame_index

    for name, its in (("animW", (60, 120)), ("animT", (60, 120)), ("W", (0, 900)),
                      ("T", (0, 900))):
        for it in its:
            (tmp_path / f"{name}.{it:010d}.meta").write_text("")
    assert [i for i, _ in frame_index(tmp_path)] == [60, 120]
    for f in tmp_path.glob("anim*"):
        f.unlink()
    fr = frame_index(tmp_path)  # pre-R29 runs: state dumps / 3-D snapshots
    assert [i for i, _ in fr] == [0, 900] and fr[1][1]["W"] == ("W", 900)
