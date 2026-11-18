import json
from pathlib import Path

import pytest

from exp6_southern_carbon.scripts import run_campaign as rc

SPY = 15768


def _run(root: Path, case: str, last: int, pickups=()):
    d = root / case
    d.mkdir(parents=True)
    (d / "run_status.json").write_text(json.dumps({"state": "OK", "last_iter": last, "n_timesteps": last}))
    for it in pickups:
        for stem in rc.PICKUP_STEMS:
            for ext in (".data", ".meta"):
                (d / f"{stem}.{it:010d}{ext}").write_text("x")


def test_branch_iter():
    assert rc.branch_iter(1900, 5_000_000) == 5_000_000 + 50 * SPY
    assert rc.branch_iter(1850, 7) == 7


def test_plan_branch_from_transient(tmp_path):
    i1850 = 4_000_000
    _run(tmp_path, "G0s", i1850)
    it = rc.branch_iter(2060, i1850)
    _run(tmp_path, "GT", i1850 + 251 * SPY, pickups=[it])
    (p,) = rc.plan(["GT_w2060"], tmp_path)
    assert p["nIter0"] == it and p["nTimeSteps"] == 25 * SPY and p["np"] == 1
    assert {f.name for f in p["pickups"]} == {f"{s}.{it:010d}{e}" for s in rc.PICKUP_STEMS for e in (".data", ".meta")}
    assert p["extra"]["iter_1850"] == i1850


def test_plan_parent_end(tmp_path):
    _run(tmp_path, "G0", 3942 * 1000, pickups=[3942 * 1000])
    _run(tmp_path, "G0s", 3942 * 1000 + 48 * SPY, pickups=[3942 * 1000 + 48 * SPY])
    (g0s,) = rc.plan(["G0s"], tmp_path)
    assert g0s["nIter0"] == 3942 * 1000 and g0s["np"] == 5
    (gc,) = rc.plan(["GC"], tmp_path)
    assert gc["nIter0"] == 3942 * 1000 + 48 * SPY


def test_plan_missing_pickup_names_file(tmp_path):
    _run(tmp_path, "G0s", 100)
    _run(tmp_path, "GT", 100 + 251 * SPY)
    with pytest.raises(FileNotFoundError, match=r"pickup\.%010d" % rc.branch_iter(1980, 100)):
        rc.plan(["GT_w1980"], tmp_path, wait_s=0)


def test_dry_run_writes_nothing(tmp_path, capsys):
    assert rc.main(["--cases", "GC_w19*", "--dry-run", "--runs-root", str(tmp_path)]) == 0
    assert list(tmp_path.iterdir()) == []
    assert "GC_w1900" in capsys.readouterr().out


def test_prepare_cold_start_writes_inputs(tmp_path):
    (item,) = rc.plan(["G0"], tmp_path, n_timesteps=3942)
    inp = rc.prepare(item, tmp_path)
    assert (inp / "case.json").exists() and (inp / "data").exists()
    assert "nTimeSteps=3942" in (inp / "data").read_text().replace(" ", "")


def test_prepare_size_mismatch_raises(tmp_path):
    (item,) = rc.plan(["G0"], tmp_path, n_timesteps=3942)
    item["config"] = "f10"                  # f10 SIZE.h is 100 x 200, the case is 20 x 40
    with pytest.raises(RuntimeError, match="Nx"):
        rc.prepare(item, tmp_path)


def test_rmtree_tolerates_vanishing_appledouble(tmp_path, monkeypatch):
    d = tmp_path / "input"
    d.mkdir()
    for n in ("bathy.bin", "._bathy.bin"):
        (d / n).write_text("x")
    real = __import__("os").unlink

    def unlink(path, *a, **k):              # exFAT: removing bathy.bin also removes ._bathy.bin
        real(path, *a, **k)
        if str(path) == "bathy.bin":
            real("._bathy.bin", *a, **k)
    monkeypatch.setattr(rc.os, "unlink", unlink)
    rc.rmtree_exfat(d)
    assert not d.exists()
