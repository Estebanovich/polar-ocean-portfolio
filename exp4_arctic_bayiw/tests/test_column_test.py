"""PD2 column test of the model's inertial frequency (scripts/column_test.py; post hoc)."""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import numpy as np

from pfcommon import campaign as cp

PKG = Path(__file__).resolve().parents[1]
NX = NY = 10


def _load():
    spec = importlib.util.spec_from_file_location("exp4_column_test", PKG / "scripts" / "column_test.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _fake_f1(d: Path, nx=12, ny=8, nz=40, nrec=181):
    """A small F1-like input directory with the F1 namelists (copied from the template rendering)."""
    d.mkdir(parents=True)
    src = PKG / "tests" / "data_f1_namelists"
    for f in src.iterdir():
        (d / f.name).write_text(f.read_text())
    t = (d / "data").read_text()
    t = re.sub(r"delX=\d+\*200\.0", f"delX={nx}*200.0", t)
    t = re.sub(r"delY=\d+\*200\.0", f"delY={ny}*200.0", t)
    (d / "data").write_text(t)
    T = np.broadcast_to(np.linspace(4, -1, nz)[:, None, None], (nz, ny, nx))
    S = np.broadcast_to(np.linspace(31, 34.5, nz)[:, None, None], (nz, ny, nx))
    np.ascontiguousarray(T, ">f8").tofile(d / "T.init.bin")
    np.ascontiguousarray(S, ">f8").tofile(d / "S.init.bin")
    pulse = 0.1 * np.exp(-0.5 * ((np.arange(nrec) + 0.5) * 1200 / 3600 - 2.5) ** 2)
    for name, comp in (("taux.bin", -0.977), ("tauy.bin", -0.215)):
        np.ascontiguousarray(np.broadcast_to((comp * pulse)[:, None, None], (nrec, ny, nx)), ">f8").tofile(d / name)
    b = -np.full((ny, nx), 150.0); b[0] = 0.0
    np.ascontiguousarray(b, ">f8").tofile(d / "bathy.bin")
    (d / "case.json").write_text(json.dumps({"params": {"nx": nx, "ny": ny, "nz": nz}}))


def test_column_inputs_keep_f1_physics_and_remove_boundaries(tmp_path):
    ct = _load()
    f1 = tmp_path / "F1in"
    _fake_f1(f1)
    out = tmp_path / "col"
    info = ct.make_column_inputs(f1, out, depth_m=153.0, ji=(4, 5))
    d = (out / "data").read_text()
    for key in ("staggerTimeStep=.TRUE.", "deltaT=120", "f0=0.00014266", "eosType='LINEAR'",
                "tAlpha=3.0456e-05", "viscAr=1.E-4", "externForcingCycle=217200"):
        assert key in d.replace(" ", ""), key
    assert f"delX={NX}*200.0" in d and f"delY={NY}*200.0" in d
    assert "useOBCS=.FALSE." in (out / "data.pkg").read_text()
    assert (out / "data.free").exists() and "zonalWindFile" not in (out / "data.free").read_text().split("&PARM05")[1]
    b = np.fromfile(out / "bathy.bin", ">f8").reshape(NY, NX)
    assert np.all(b == -153.0)                                    # no walls: doubly periodic
    T = np.fromfile(out / "T.init.bin", ">f8").reshape(40, NY, NX)
    assert np.allclose(T, np.linspace(4, -1, 40)[:, None, None])
    tx = np.fromfile(out / "taux.bin", ">f8").reshape(-1, NY, NX)
    assert tx.shape[0] == 181 and np.allclose(tx.max(), 0.0) and tx.min() < -0.09
    assert info["boundaries"].startswith("doubly periodic")
    cp.check_size(out, PKG / "config" / "fjord_vmf_column" / "code" / "SIZE.h")


def test_column_run_refuses_other_model_processes_unless_allowed(tmp_path, monkeypatch):
    ct = _load()
    monkeypatch.setattr(ct.cp, "_pgrep", lambda: "4242 model-process")
    import pytest
    with pytest.raises(RuntimeError, match="allow-concurrent"):
        ct.main(["run", "--input", str(tmp_path), "--run", str(tmp_path / "r")])
    called = []
    monkeypatch.setattr(ct, "run", lambda *a, **k: called.append(1) or "OK")
    assert ct.main(["run", "--input", str(tmp_path), "--run", str(tmp_path / "r"), "--allow-concurrent"]) == 0
    assert called == [1]


def test_column_analysis_reports_within_1df_not_pass():
    import inspect
    ct = _load()
    src = inspect.getsource(ct.analyse)
    assert "within_1df" in src and 'r["pass"]' not in src.replace('r.pop("pass")', "")
