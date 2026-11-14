"""scripts/restart_test.py run() end to end with a fake runner (no MITgcm is launched)."""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import pytest

from pfcommon import campaign as cp

PKG = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location("exp4_restart_test", PKG / "scripts" / "restart_test.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def fake_make_case(case, out, overrides=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    n = round(overrides["t_end_h"] * 3600 / 120)
    (out / "data").write_text(f" &PARM03\n nIter0=0,\n nTimeSteps={n},\n deltaT=120.,\n"
                              f" pChkptFreq=21600.,\n &\n")
    return {"params": {"case": case}}


def fake_runner(calls, perturb_stage2=False):
    """Stand-in for run_mitgcm.sh: writes run_status.json, the pickups at every 10th step after
    nIter0 and the diagnostics at every 10th step from nIter0, with content that depends only on
    the iteration (a perfect restart), or not (perturb_stage2)."""
    def sh(cmd):
        _, build, inp, rdir, np_ = cmd
        calls.append((Path(inp).name, Path(rdir).name, np_))
        text = (Path(inp) / "data").read_text()
        n0, n = cp.nml_int(text, "nIter0"), cp.nml_int(text, "nTimeSteps")
        assert cp.nml_float(text, "pChkptFreq") == 1200.0          # restart_test sets it
        rdir = Path(rdir)
        rdir.mkdir(parents=True, exist_ok=True)
        tag = b"x" if (perturb_stage2 and n0 > 0) else b""
        for it in range(n0 + 10, n0 + n + 1, 10):
            for stem in ("pickup", "pickup_seaice"):
                for ext in ("data", "meta"):
                    (rdir / f"{stem}.{it:010d}.{ext}").write_bytes(f"{stem}{it}".encode() + tag)
            (rdir / f"pickup_orlanskiW.{it:010d}.001.001.data").write_bytes(f"ob{it}".encode() + tag)
        for it in range(n0, n0 + n, 10):
            for stem in ("state3d", "surf", "seaice"):
                (rdir / f"{stem}.{it:010d}.data").write_bytes(f"{stem}{it}".encode() + tag)
        (rdir / "run_status.json").write_text(json.dumps(
            {"state": "OK", "last_iter": n0 + n, "n_timesteps": n}))
        return 0
    return sh


@pytest.mark.parametrize("perturb, identical", [(False, True), (True, False)])
def test_restart_run_end_to_end_with_a_fake_runner(tmp_path, perturb, identical):
    rt = _load()
    calls = []
    res = rt.run("F4", runs=tmp_path, builds=tmp_path / "builds", sh=fake_runner(calls, perturb),
                 pgrep=lambda: "", make_case=fake_make_case)
    assert res == {"single": "OK", "split": "OK"}
    # one 20-step run, then 10 + pickup + 10, all on 5 ranks
    assert [c[1:] for c in calls] == [("F4_restarttest_single", 5), ("F4_restarttest_split", 5),
                                      ("stage2", 5)]
    st2 = tmp_path / "F4_restarttest_split" / "input_stage2"
    assert (st2 / "pickup_seaice.0000000010.data").exists()          # pickups were copied
    assert "nIter0=10" in re.sub(r"\s", "", (st2 / "data").read_text())
    cmp = rt.compare("F4", runs=tmp_path)
    assert cmp["bit_identical"] is identical


def test_restart_run_refuses_while_a_model_runs(tmp_path):
    rt = _load()
    with pytest.raises(RuntimeError, match="another"):
        rt.run("F4", runs=tmp_path, builds=tmp_path, sh=lambda c: 0,
               pgrep=lambda: "123 model-process", make_case=fake_make_case)
