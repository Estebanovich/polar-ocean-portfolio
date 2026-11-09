"""Tests for exp2_isomip_cavity/scripts/run_campaign.py with a fake shell (no MITgcm run)."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("exp2_run_campaign", PKG / "scripts" / "run_campaign.py")
rc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rc)

SIZE = ("     &           sNx = {nx},\n     &           sNy = {ny},\n     &           nSx =   1,\n"
        "     &           nSy =   1,\n     &           nPx =   1,\n     &           nPy =   1,\n"
        "     &           Nr  =  36)\n")


class Shell:
    def __init__(self):
        self.calls = []

    def __call__(self, cmd):
        self.calls.append([str(c) for c in cmd])
        if Path(cmd[0]).name == "build_mitgcm.sh":           # fake build: an executable file
            from pfcommon.campaign import EXE_NAME
            Path(cmd[2]).mkdir(parents=True, exist_ok=True)
            (Path(cmd[2]) / EXE_NAME).write_bytes(b"fake binary " + str(cmd[1]).encode())
        if Path(cmd[0]).name == "run_mitgcm.sh":
            from pfcommon.campaign import nml_int
            inp, run = Path(cmd[2]), Path(cmd[3])
            run.mkdir(parents=True, exist_ok=True)
            d = (inp / "data").read_text()
            n0, n = nml_int(d, "nIter0"), nml_int(d, "nTimeSteps")
            for stem in ("pickup", "pickup_ptracers"):
                for ext in ("data", "meta"):
                    (run / f"{stem}.{n0 + n:010d}.{ext}").write_text("x")
            (run / "run_status.json").write_text(json.dumps(
                {"state": "OK", "last_iter": n0 + n, "n_timesteps": n, "wall_s": 0}))
        return 0


@pytest.fixture
def cfg_root(tmp_path):
    root = tmp_path / "config"
    for cfg, nx, ny in (("ocean0_2km", 240, 41), ("ocean0_4km", 120, 21)):
        (root / cfg / "code").mkdir(parents=True)
        (root / cfg / "code" / "SIZE.h").write_text(SIZE.format(nx=nx, ny=ny))
    return root


def _run(tmp_path, cfg_root, templates, synth_geom, config, pattern, **kw):
    sh = Shell()
    res = rc.run_campaign(config, pattern, kw.pop("np_", 1), max_concurrent=kw.pop("mc", 1),
                          runs_root=tmp_path / "runs", build_root=tmp_path / "builds",
                          config_root=cfg_root, template_dir=templates, geom_path=synth_geom,
                          sh=sh, pgrep=kw.pop("pgrep", lambda: ""), **kw)
    return res, sh


def test_refuses_when_model_running(tmp_path, cfg_root, templates, synth_geom):
    with pytest.raises(RuntimeError, match="refusing"):
        _run(tmp_path, cfg_root, templates, synth_geom, "ocean0_2km", "cal2_*",
             pgrep=lambda: "99 mpirun")
    assert not (tmp_path / "runs").exists()


def test_concurrency_rules_enforced(tmp_path, cfg_root, templates, synth_geom):
    for np_, mc in ((1, 6), (5, 2), (3, 1)):
        with pytest.raises(ValueError):
            _run(tmp_path, cfg_root, templates, synth_geom, "ocean0_2km", "cal2_*", np_=np_, mc=mc)


def test_campaign_builds_once_and_runs_all(tmp_path, cfg_root, templates, synth_geom):
    res, sh = _run(tmp_path, cfg_root, templates, synth_geom, "ocean0_4km", "cal4_*", mc=3,
                   overrides={"endTime": 4800.0, "deltaT": 1200.0})
    # four calibration cases since the PD2 sampling amendment (0.0075 added, 2026-11-09)
    assert res == {"cal4_gT0.015": "OK", "cal4_gT0.03": "OK", "cal4_gT0.06": "OK",
                   "cal4_gT0.0075": "OK"}
    kinds = [Path(c[0]).name for c in sh.calls]
    assert kinds.count("build_mitgcm.sh") == 1 and kinds[0] == "build_mitgcm.sh"
    assert "--mpi" not in sh.calls[0]
    assert kinds.count("run_mitgcm.sh") == 4
    for n in res:
        fin = json.loads((tmp_path / "runs" / n / "run_status_final.json").read_text())
        assert fin["state"] == "OK" and fin["last_iter"] == 4


def test_campaign_stages_and_skip_ok(tmp_path, cfg_root, templates, synth_geom):
    ov = {"endTime": 4800.0, "deltaT": 1200.0, "gammaT": 0.03}
    res, sh = _run(tmp_path, cfg_root, templates, synth_geom, "ocean0_2km", "ocean1_2km",
                   overrides=ov, n_stages=2)
    assert res == {"ocean1_2km": "OK"}
    run = tmp_path / "runs" / "ocean1_2km"
    assert (run / "stage2" / "run_status.json").exists()
    assert "nIter0 = 2," in (run / "input_stage2" / "data").read_text()
    res2, sh2 = _run(tmp_path, cfg_root, templates, synth_geom, "ocean0_2km", "ocean1_2km",
                     overrides=ov, n_stages=2, skip_ok=True)
    assert res2 == {"ocean1_2km": "OK"}
    assert not any(Path(c[0]).name == "run_mitgcm.sh" for c in sh2.calls)
    with pytest.raises(RuntimeError, match="other parameters"):
        _run(tmp_path, cfg_root, templates, synth_geom, "ocean0_2km", "ocean1_2km",
             overrides={**ov, "deltaT": 600.0}, n_stages=2, skip_ok=True)


def test_mpi_build_for_np5(tmp_path, cfg_root, templates, synth_geom):
    (cfg_root / "ocean0_4km" / "code" / "SIZE.h").write_text(
        SIZE.format(nx=24, ny=21).replace("nPx =   1", "nPx =   5"))
    res, sh = _run(tmp_path, cfg_root, templates, synth_geom, "ocean0_4km", "cal4_gT0.03",
                   np_=5, mc=1, overrides={"endTime": 4800.0, "deltaT": 1200.0})
    assert res == {"cal4_gT0.03": "OK"} and "--mpi" in sh.calls[0]


def test_dry_run_lists_cases(capsys):
    assert rc.main(["--config", "ocean0_4km", "--cases", "ens_*", "--dry-run"]) == 0
    out = capsys.readouterr().out.split()
    assert len(out) == 28 and all(n.startswith("ens_") for n in out)


def test_skip_ok_never_reuses_a_void_law_run(tmp_path, cfg_root, templates, synth_geom):
    """R21: an OK run whose case.json has no melt_law (made with the void pre-R21 u* law,
    Cd u^2 + u_tidal^2) is neither skipped nor overwritten: --skip-ok refuses to start."""
    ov = {"endTime": 4800.0, "deltaT": 1200.0, "gammaT": 0.03}
    res, _ = _run(tmp_path, cfg_root, templates, synth_geom, "ocean0_4km", "ocean0_4km",
                  overrides=ov)
    assert res == {"ocean0_4km": "OK"}
    cj = tmp_path / "runs" / "ocean0_4km" / "input" / "case.json"
    meta = json.loads(cj.read_text())
    assert meta["params"]["melt_law"] == "AD16_eq27_r21"        # new runs carry the tag
    del meta["params"]["melt_law"]                              # what a void run looks like
    cj.write_text(json.dumps(meta))
    with pytest.raises(RuntimeError, match="void pre-R21"):
        _run(tmp_path, cfg_root, templates, synth_geom, "ocean0_4km", "ocean0_4km",
             overrides=ov, skip_ok=True)
    assert "melt_law" not in json.loads(cj.read_text())["params"]   # left untouched


def test_production_run_suffix_defaults_to_r21(monkeypatch):
    seen = {}
    monkeypatch.setattr(rc, "run_campaign", lambda *a, **k: seen.update(k) or {})
    rc.main(["--config", "ocean0_4km", "--cases", "ocean0_4km"])
    assert rc.DEFAULT_RUN_SUFFIX == "_r21" and seen["run_suffix"] == "_r21"


def test_binary_sha256_recorded_at_run_time(tmp_path, cfg_root, templates, synth_geom):
    """Provenance: every run's input/case.json carries the sha256 of the executable of
    the build it ran with, written by the campaign runner before the run starts."""
    import hashlib
    from pfcommon.campaign import EXE_NAME
    res, sh = _run(tmp_path, cfg_root, templates, synth_geom, "ocean0_4km", "cal4_gT0.03",
                   overrides={"endTime": 4800.0, "deltaT": 1200.0})
    assert res == {"cal4_gT0.03": "OK"}
    exe = tmp_path / "builds" / "ocean0_4km" / EXE_NAME
    meta = json.loads((tmp_path / "runs" / "cal4_gT0.03" / "input" / "case.json").read_text())
    assert meta["binary_sha256"] == hashlib.sha256(exe.read_bytes()).hexdigest()
    assert meta["binary_path"] == str(exe) and meta["binary_recorded_at"]
    assert meta["params"]["melt_law"] == "AD16_eq27_r21"
