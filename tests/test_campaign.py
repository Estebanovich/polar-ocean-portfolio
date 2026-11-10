"""Tests for pfcommon/campaign.py: shared campaign helpers, concurrency pool and stage chaining.

MITgcm is never run: a fake runner stands in for scripts/run_mitgcm.sh.
"""
from __future__ import annotations

import importlib
import json
import struct
import threading
import time
from pathlib import Path

import pytest

from pfcommon import campaign as cp

DATA = (" &PARM03\n nIter0 = 0,\n nTimeSteps = 120,\n deltaT = 100.,\n"
        " pChkptFreq = 4000.,\n &\n")


# ---- paths ---------------------------------------------------------------------------

def test_paths_exp2_defaults_and_overrides(monkeypatch):
    import pfcommon.paths as paths
    for k in ("EXP1_RUNS", "EXP1_BUILDS", "EXP2_RUNS", "EXP2_BUILDS"):
        monkeypatch.delenv(k, raising=False)
    paths = importlib.reload(paths)
    assert paths.EXP1_RUNS == Path("/Volumes/LaCie/portfolio_runs/exp1")
    assert paths.EXP1_BUILDS == Path("/Volumes/ESTEBAN_MIT/portfolio_builds/exp1")
    assert paths.EXP2_RUNS == Path("/Volumes/LaCie/portfolio_runs/exp2")
    assert paths.EXP2_BUILDS == Path("/Volumes/ESTEBAN_MIT/portfolio_builds/exp2")
    monkeypatch.setenv("EXP2_RUNS", "/tmp/r2")
    monkeypatch.setenv("EXP2_BUILDS", "/tmp/b2")
    paths = importlib.reload(paths)
    assert paths.EXP2_RUNS == Path("/tmp/r2") and paths.EXP2_BUILDS == Path("/tmp/b2")
    monkeypatch.delenv("EXP2_RUNS")
    monkeypatch.delenv("EXP2_BUILDS")
    importlib.reload(paths)


# ---- copied helpers ------------------------------------------------------------------

CASES = {"a1": {"config": "A"}, "a2": {"config": "A"}, "b1": {"config": "B"}}


def test_select_cases():
    assert cp.select_cases(CASES, "A", "*") == ["a1", "a2"]
    with pytest.raises(ValueError, match="b1"):
        cp.select_cases(CASES, "A", "b*")
    with pytest.raises(ValueError, match="no case"):
        cp.select_cases(CASES, "A", "zz*")


def test_other_processes_by_name():
    assert cp.other_mitgcm_processes(pgrep=lambda: "") == []
    assert cp.other_mitgcm_processes(pgrep=lambda: "77 mpirun\n") == ["77 mpirun"]
    assert "-l" in cp.PGREP_ARGS and "-f" not in cp.PGREP_ARGS


def test_parse_sets_and_overrides():
    assert cp.parse_sets(["a=1", "b=2.5"]) == {"a": 1, "b": 2.5}
    with pytest.raises(ValueError):
        cp.parse_sets(["a"])
    spec = {"configs": {"A": {"deltaT": 100.0, "n": 3, "gammaT": None, "name": "x"}},
            "cases": {"a1": {"config": "A"}}}
    assert cp.check_overrides(["a1"], {"deltaT": 50, "n": 4, "gammaT": 0.03}, spec) == \
        {"deltaT": 50.0, "n": 4, "gammaT": 0.03}
    for bad in ({"nope": 1}, {"name": 1}, {"n": 1.5}):
        with pytest.raises(ValueError):
            cp.check_overrides(["a1"], bad, spec)


def test_check_size(tmp_path):
    (tmp_path / "case.json").write_text(json.dumps({"params": {"nx": 240, "ny": 40, "nz": 36}}))
    size = tmp_path / "SIZE.h"
    size.write_text("     &           sNx =  48,\n     &           sNy =  40,\n     &           nSx =   1,\n"
                    "     &           nSy =   1,\n     &           nPx =   5,\n     &           nPy =   1,\n"
                    "     &           Nr  =  36)\n")
    cp.check_size(tmp_path, size)
    size.write_text(size.read_text().replace("nPx =   5", "nPx =   4"))
    with pytest.raises(RuntimeError, match="Nx"):
        cp.check_size(tmp_path, size)


def test_finished_ok(tmp_path):
    assert not cp.finished_ok(tmp_path)
    (tmp_path / "run_status.json").write_text(json.dumps({"state": "OK", "last_iter": 10, "n_timesteps": 10}))
    assert cp.finished_ok(tmp_path)
    (tmp_path / "run_status_final.json").write_text(json.dumps(
        {"state": "OK", "last_iter": 10, "expected_last_iter": 30}))
    assert not cp.finished_ok(tmp_path)          # the final status wins when present


def test_build_up_to_date(tmp_path):
    cfg, b = tmp_path / "cfg", tmp_path / "build"
    (cfg / "code").mkdir(parents=True)
    b.mkdir()
    (cfg / "code" / "SIZE.h").write_text("x")
    assert not cp.build_up_to_date(cfg, b, mpi=False)
    time.sleep(0.01)
    (b / cp.EXE_NAME).write_text("bin")
    (b / "Makefile").write_text("DEFINES = -DALLOW_USE_MPI\n")
    assert cp.build_up_to_date(cfg, b, mpi=True)
    assert not cp.build_up_to_date(cfg, b, mpi=False)


def test_check_free_space(tmp_path):
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "f").write_bytes(b"0" * 1000)

    class U:
        def __init__(self, free):
            self.free = free
    cp.check_free_space(tmp_path, [tmp_path / "in"], disk_usage=lambda p: U(30 * 10**9))
    with pytest.raises(RuntimeError, match="free space"):
        cp.check_free_space(tmp_path, [tmp_path / "in"], disk_usage=lambda p: U(10 * 10**9))


# ---- concurrency -------------------------------------------------------------------

@pytest.mark.parametrize("np_each,mc", [(1, 1), (1, 3), (1, 5), (5, 1)])
def test_concurrency_allowed(np_each, mc):
    cp.check_concurrency(np_each, mc)


@pytest.mark.parametrize("np_each,mc", [(1, 0), (1, 6), (5, 2), (2, 1), (3, 1), (6, 1), (0, 1)])
def test_concurrency_rejected(np_each, mc):
    with pytest.raises(ValueError):
        cp.check_concurrency(np_each, mc)


@pytest.mark.parametrize("mc", [1, 2, 5])
def test_run_pool_never_exceeds_max_concurrent(mc):
    lock, live, peak = threading.Lock(), [0], [0]

    def run_one(case):
        with lock:
            live[0] += 1
            peak[0] = max(peak[0], live[0])
        time.sleep(0.03)
        with lock:
            live[0] -= 1
        return "OK"

    cases = [f"c{i}" for i in range(12)]
    res = cp.run_pool(cases, 1, mc, run_one)
    assert res == {c: "OK" for c in cases} and list(res) == cases
    assert peak[0] <= mc
    assert peak[0] == mc                     # the pool really runs cases side by side


def test_run_pool_reports_errors_and_continues():
    def run_one(case):
        if case == "bad":
            raise RuntimeError("boom")
        return "OK"
    res = cp.run_pool(["a", "bad", "c"], 1, 2, run_one)
    assert res["a"] == res["c"] == "OK" and res["bad"].startswith("ERROR")


def test_run_pool_rejects_bad_concurrency():
    with pytest.raises(ValueError):
        cp.run_pool(["a"], 5, 2, lambda c: "OK")


# ---- stage chaining (ruling R2) ----------------------------------------------------

def test_stage_plan_and_namelist_rewrite():
    assert cp.stage_plan(0, 120, 3) == [(0, 40), (40, 40), (80, 40)]
    assert cp.stage_plan(0, 120, 1) == [(0, 120)]
    with pytest.raises(ValueError):
        cp.stage_plan(0, 100, 3)
    t = cp.set_namelist_ints(DATA, {"nIter0": 40, "nTimeSteps": 40})
    assert cp.nml_int(t, "nIter0") == 40 and cp.nml_int(t, "nTimeSteps") == 40
    assert "deltaT = 100." in t
    with pytest.raises(ValueError):
        cp.set_namelist_ints(" &P\n deltaT = 1.,\n &\n", {"nIter0": 1})
    with pytest.raises(ValueError, match="endTime"):
        cp.set_namelist_ints(DATA + " endTime = 5.,\n", {"nIter0": 1})


class FakeModel:
    """Fake run_mitgcm.sh: a deterministic float64 'model' x <- x * 1.0001 + 1e-3 per step.

    It restarts from pickup.<nIter0>.data (in its input dir) when nIter0 > 0, writes
    pickup files every pChkptFreq and run_status.json like the real runner.
    """

    def __init__(self, fail_stage_iter0=None, extra_pickup=True):
        self.calls, self.fail_stage_iter0, self.extra = [], fail_stage_iter0, extra_pickup

    def __call__(self, cmd):
        _, build, inp, run, np_ = (str(c) for c in cmd)
        inp, run = Path(inp), Path(run)
        self.calls.append((inp, run))
        run.mkdir(parents=True, exist_ok=True)
        for f in inp.iterdir():
            if f.is_file():
                (run / f.name).write_bytes(f.read_bytes())
        d = (run / "data").read_text()
        n0, n = cp.nml_int(d, "nIter0"), cp.nml_int(d, "nTimeSteps")
        dt, ck = cp.nml_float(d, "deltaT"), cp.nml_float(d, "pChkptFreq")
        x = 1.0 if n0 == 0 else struct.unpack(">d", (run / f"pickup.{n0:010d}.data").read_bytes())[0]
        it = n0
        for _ in range(n):
            x = x * 1.0001 + 1e-3
            it += 1
            if (it * dt) % ck == 0:
                for stem in ("pickup",) + (("pickup_ptracers",) if self.extra else ()):
                    (run / f"{stem}.{it:010d}.data").write_bytes(struct.pack(">d", x))
                    (run / f"{stem}.{it:010d}.meta").write_text("meta")
        state = "FAILED" if n0 == self.fail_stage_iter0 else "OK"
        (run / "run_status.json").write_text(json.dumps(
            {"state": state, "last_iter": it, "n_timesteps": n, "wall_s": 0}))
        return 0


def _inp(tmp_path, data=DATA):
    inp = tmp_path / "input"
    inp.mkdir()
    (inp / "data").write_text(data)
    (inp / "field.bin").write_bytes(b"\0" * 16)
    return inp


def test_run_stages_chain_restarts_from_pickups(tmp_path):
    inp = _inp(tmp_path)
    fake = FakeModel()
    state = cp.run_stages(tmp_path / "run3", inp, tmp_path / "b", 1, 3, fake, "run_mitgcm.sh",
                          pickup_stems=("pickup", "pickup_ptracers"))
    assert state == "OK"
    assert [cp.nml_int((i / "data").read_text(), "nIter0") for i, _ in fake.calls] == [0, 40, 80]
    assert [r.name for _, r in fake.calls] == ["run3", "stage2", "stage3"]
    for k in (2, 3):
        si = tmp_path / "run3" / f"input_stage{k}"
        assert (si / f"pickup.{40 * (k - 1):010d}.data").exists()
        assert (si / f"pickup_ptracers.{40 * (k - 1):010d}.meta").exists()
        assert not any(p.is_symlink() for p in si.rglob("*"))      # copies, never links (exFAT)
    fin = json.loads((tmp_path / "run3" / "run_status_final.json").read_text())
    assert fin["state"] == "OK" and fin["last_iter"] == fin["expected_last_iter"] == 120
    assert set(fin["stages"]) == {"1", "2", "3"}
    assert (inp / "data").read_text() == DATA                      # generated input untouched
    assert cp.finished_ok(tmp_path / "run3")


def test_staged_run_bit_identical_to_single_run(tmp_path):
    inp = _inp(tmp_path)
    assert cp.run_stages(tmp_path / "one", inp, tmp_path / "b", 1, 1, FakeModel(), "run_mitgcm.sh") == "OK"
    assert cp.run_stages(tmp_path / "three", inp, tmp_path / "b", 1, 3, FakeModel(), "run_mitgcm.sh") == "OK"
    tag = f"{120:010d}"
    assert cp.compare_bitwise(tmp_path / "one", tmp_path / "three" / "stage3", f"pickup.{tag}.data") == []
    # the helper does see a one-bit difference
    f = tmp_path / "three" / "stage3" / f"pickup.{tag}.data"
    b = bytearray(f.read_bytes()); b[-1] ^= 1; f.write_bytes(bytes(b))
    assert cp.compare_bitwise(tmp_path / "one", tmp_path / "three" / "stage3", f"pickup.{tag}.data") \
        == [f"pickup.{tag}.data"]
    assert cp.compare_bitwise(tmp_path / "one", tmp_path / "nowhere", "pickup.*.data") != []


def test_run_stages_stops_on_failed_stage(tmp_path):
    inp = _inp(tmp_path)
    fake = FakeModel(fail_stage_iter0=40)
    assert cp.run_stages(tmp_path / "r", inp, tmp_path / "b", 1, 3, fake, "run_mitgcm.sh") == "FAILED"
    assert len(fake.calls) == 2
    fin = json.loads((tmp_path / "r" / "run_status_final.json").read_text())
    assert fin["state"] == "FAILED" and "stage 2" in fin["reason"]


def test_run_stages_missing_required_pickup(tmp_path):
    inp = _inp(tmp_path)
    fake = FakeModel(extra_pickup=False)
    st = cp.run_stages(tmp_path / "r", inp, tmp_path / "b", 1, 3, fake, "run_mitgcm.sh",
                       pickup_stems=("pickup", "pickup_ptracers"))
    assert st == "FAILED" and len(fake.calls) == 1
    fin = json.loads((tmp_path / "r" / "run_status_final.json").read_text())
    assert "pickup_ptracers" in fin["reason"]


def test_run_stages_rejects_stage_end_off_checkpoint(tmp_path):
    inp = _inp(tmp_path, DATA.replace("pChkptFreq = 4000.", "pChkptFreq = 3000."))
    with pytest.raises(ValueError, match="pChkptFreq"):
        cp.run_stages(tmp_path / "r", inp, tmp_path / "b", 1, 3, FakeModel(), "run_mitgcm.sh")


@pytest.mark.parametrize("n_stages", [1, 3])
def test_run_stages_stale_ok_status_with_failing_runner(tmp_path, n_stages):
    """The runner exits 2 without writing anything (e.g. LaCie not mounted) over a run dir
    that holds an earlier OK attempt: the case must be FAILED, never the stale OK."""
    inp = _inp(tmp_path)
    run = tmp_path / "r"
    run.mkdir()
    ok = {"state": "OK", "last_iter": 120 if n_stages == 1 else 40, "n_timesteps": 40, "wall_s": 0}
    (run / "run_status.json").write_text(json.dumps(ok))
    calls = []

    def dead_runner(cmd):
        calls.append(cmd)
        return 2

    assert cp.run_stages(run, inp, tmp_path / "b", 1, n_stages, dead_runner, "run_mitgcm.sh") == "FAILED"
    assert len(calls) == 1
    fin = json.loads((run / "run_status_final.json").read_text())
    assert fin["state"] == "FAILED" and "code 2" in fin["reason"]
    assert fin["stages"]["1"] is None                  # the stale status was removed, not read
    assert not cp.finished_ok(run)


def test_run_stages_nonzero_exit_fails_even_with_ok_status(tmp_path):
    inp = _inp(tmp_path)
    fake = FakeModel()
    st = cp.run_stages(tmp_path / "r", inp, tmp_path / "b", 1, 1, lambda c: fake(c) or 1, "run_mitgcm.sh")
    assert st == "FAILED"
    assert "code 1" in json.loads((tmp_path / "r" / "run_status_final.json").read_text())["reason"]
