"""Tests for scripts/b_gate.py (R30): the first-case gate of the (b) campaign."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

PKG = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("b_gate", PKG / "scripts" / "b_gate.py")
bg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bg)

NY, NX, NZ, DT, IT = 4, 5, 3, 2.0, 900


def _mds(d: Path, name: str, it: int, arr: np.ndarray) -> None:
    arr = np.asarray(arr, dtype=float)
    stem = d / f"{name}.{it:010d}"
    np.ascontiguousarray(arr, dtype=">f8").tofile(f"{stem}.data")
    dims = "".join(f" {n}, 1, {n},\n" for n in reversed(arr.shape))
    Path(f"{stem}.meta").write_text(
        f" nDims = [ {arr.ndim} ];\n dimList = [\n{dims} ];\n dataprec = [ 'float64' ];\n"
        f" nrecords = [ 1 ];\n timeStepNumber = [ {it} ];\n")


def _common(d: Path, n_steps: int, state: str = "OK") -> None:
    d.mkdir()
    (d / "run_status.json").write_text(json.dumps(
        {"state": state, "last_iter": n_steps, "n_timesteps": n_steps, "wall_s": 1}))
    (d / "case.json").write_text(json.dumps(
        {"params": {"deltaT": DT, "nx": NX, "ny": NY, "nz": NZ, "dx": 10.0, "dy": 10.0},
         "Q0_total": 2.0}))
    depth = np.full((NY, NX), 30.0)
    depth[0] = 0.0
    np.ascontiguousarray(-depth, dtype=">f8").tofile(d / "topog.bin")


def _fields(seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    return {k: rng.normal(size=(NZ, NY, NX)) for k in ("ICF", "T", "S", "U", "V")} | {
        "Eta": rng.normal(size=(NY, NX))}


def _ref(tmp: Path, f: dict) -> Path:
    d = tmp / "ref"
    _common(d, IT)
    _mds(d, "icefrontDiag", IT, f["ICF"])
    for dump in ("T", "S", "U", "V", "Eta"):
        _mds(d, dump, IT, f[dump])
    return d


def _run(tmp: Path, f: dict, retained: float) -> Path:
    d = tmp / "run"
    _common(d, 5400)
    _mds(d, "icefrontDiag", IT, f["ICF"])
    for diag, dump in (("THETA", "T"), ("SALT", "S"), ("UVEL", "U"), ("VVEL", "V"),
                       ("ETAN", "Eta")):
        _mds(d, diag, IT, f[dump])
    area = 3 * 5 * 100.0
    wet = np.ones((NY, NX), bool)
    wet[0] = False
    for it in (900, 1800, 2700, 3600, 4500):
        _mds(d, "ETAN", it, np.where(wet, retained * 2.0 * it * DT / area, 0.0)
             if it != IT else f["Eta"])
    _mds(d, "Eta", 5400, np.where(wet, retained * 2.0 * 5400 * DT / area, 0.0))
    return d


def test_gate_passes(tmp_path):
    f = _fields()
    ok, reasons, det = bg.gate(_run(tmp_path, f, 0.3), _ref(tmp_path, f))
    assert ok, reasons
    assert all(v["identical"] for v in det["identity"].values())


@pytest.mark.parametrize("retained", [0.7, -0.7])
def test_gate_fails_on_retention(tmp_path, retained):
    f = _fields()
    ok, reasons, det = bg.gate(_run(tmp_path, f, retained), _ref(tmp_path, f))
    assert not ok and any("retained_last_hour" in r for r in reasons)


def test_gate_fails_and_reports_max_diff_when_not_identical(tmp_path):
    f = _fields()
    ref = _ref(tmp_path, f)
    g = dict(f)
    g["ICF"] = f["ICF"].copy()
    g["ICF"][1, 2, 0] += 1e-9
    ok, reasons, det = bg.gate(_run(tmp_path, g, 0.3), ref)
    assert not ok
    assert det["identity"]["icefrontDiag"]["identical"] is False
    assert det["identity"]["icefrontDiag"]["max_abs_diff"] == pytest.approx(1e-9, rel=1e-3)
    assert any("icefrontDiag" in r for r in reasons)


def test_gate_fails_when_reference_file_missing(tmp_path):
    f = _fields()
    ref = _ref(tmp_path, f)
    (ref / f"S.{IT:010d}.data").unlink()
    ok, reasons, _ = bg.gate(_run(tmp_path, f, 0.3), ref)
    assert not ok and any("missing" in r for r in reasons)


def test_gate_main_exit_codes(tmp_path, capsys):
    f = _fields()
    run, ref = _run(tmp_path, f, 0.3), _ref(tmp_path, f)
    assert bg.main([str(run), "--ref", str(ref)]) == 0
    assert "PASS" in capsys.readouterr().out
    (ref / f"U.{IT:010d}.data").unlink()
    assert bg.main([str(run), "--ref", str(ref)]) == 1
    assert "FAIL" in capsys.readouterr().out
