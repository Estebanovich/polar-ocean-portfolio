"""run_protocol.py --freeze-config (Task 7): frozen PD0 configuration, never overwritten."""
import hashlib
import json

import pytest

pytest.importorskip("sklearn")
pytest.importorskip("torch")

from exp3_melt_emulator.folds import drop_n1, make_folds  # noqa: E402
from exp3_melt_emulator.methods import REGISTRY  # noqa: E402
from exp3_melt_emulator.scripts import run_protocol as rp  # noqa: E402
from exp3_melt_emulator.tests.synth import synthetic_runs, write_npz  # noqa: E402

FAKE_GIT = {"commit": "0" * 40, "dirty_paths": []}


@pytest.fixture
def dataset(tmp_path):
    runs = synthetic_runs()
    path = tmp_path / "ds.npz"
    write_npz(runs, path)
    return runs, path, hashlib.sha256(path.read_bytes()).hexdigest()


def _freeze(tmp_path, dataset, out=None, noise=None):
    _, path, sha = dataset
    out = out or tmp_path / "figures" / "protocol_config.json"
    rp.freeze_config(out, dataset_path=path, dataset_sha=sha, noise_json=noise,
                     git_state=lambda: dict(FAKE_GIT))
    return out


def test_freeze_records_folds_methods_thresholds(tmp_path, dataset):
    runs, _, sha = dataset
    cfg = json.loads(_freeze(tmp_path, dataset).read_text())
    assert cfg["dataset"]["sha256"] == sha
    assert cfg["git"] == FAKE_GIT
    folds = make_folds(runs)
    assert [f["name"] for f in cfg["folds"]] == [f.name for f in folds]
    for rec, f in zip(cfg["folds"], folds):
        assert rec["kind"] == f.kind and rec["train"] == list(f.train) and rec["test"] == list(f.test)
    cv = [f for f in folds if f.kind != "final"]
    assert [r["name"] for r in cfg["folds_sensitivity_n1"]] == [f.name for f in cv]
    for rec, f in zip(cfg["folds_sensitivity_n1"], cv):
        assert rec["train"] == list(drop_n1(f, runs).train) and rec["test"] == list(f.test)
    assert set(cfg["methods"]) == set(REGISTRY)
    for mid, rec in cfg["methods"].items():
        assert rec == json.loads(json.dumps(REGISTRY[mid]().hparams()))   # constructor defaults
    assert cfg["methods"]["M6"]["max_epochs"] == 5000 and cfg["methods"]["M6"]["seeds"] == [0, 1, 2, 3, 4]
    assert cfg["methods"]["M7"]["tf_kind_folds"] == 5
    assert cfg["methods"]["M5"]["var_fraction"] == 0.99 and cfg["methods"]["M5"]["n_restarts"] == 5
    assert cfg["methods"]["M4"]["degrees"] == [1, 2]
    assert cfg["thresholds"] == json.loads(json.dumps(rp.PD0_THRESHOLDS))
    assert cfg["thresholds"]["H3b"]["max_ratio_M7_over_M6"] == 1.10
    assert cfg["thresholds"]["gl_frac"] == 0.10


def test_freeze_refuses_to_overwrite(tmp_path, dataset):
    out = _freeze(tmp_path, dataset)
    before = out.read_bytes()
    with pytest.raises(FileExistsError):
        _freeze(tmp_path, dataset, out=out)
    assert out.read_bytes() == before


def test_freeze_checks_dataset_sha(tmp_path, dataset):
    runs, path, _ = dataset
    with pytest.raises(ValueError, match="sha256"):
        rp.freeze_config(tmp_path / "c.json", dataset_path=path, dataset_sha="f" * 64,
                         git_state=lambda: dict(FAKE_GIT))
    assert not (tmp_path / "c.json").exists()


def test_freeze_records_noise_floor(tmp_path, dataset):
    noise = tmp_path / "noise_floor.json"
    noise.write_text(json.dumps({"sha256": "a" * 64, "pooled_median_rms": 0.16,
                                 "median_abs_integrated_rel_diff": 0.013}))
    cfg = json.loads(_freeze(tmp_path, dataset, noise=noise).read_text())
    nfl = cfg["noise_floor"]
    assert nfl["pooled_median_rms"] == 0.16 and nfl["npz_sha256"] == "a" * 64
    assert nfl["json_sha256"] == hashlib.sha256(noise.read_bytes()).hexdigest()


def test_cli_freeze_and_refuse(tmp_path, dataset, monkeypatch):
    _, path, sha = dataset
    monkeypatch.setattr(rp, "git_state", lambda: dict(FAKE_GIT))
    out = tmp_path / "p.json"
    args = ["--freeze-config", "--out", str(out), "--dataset", str(path), "--dataset-sha", sha,
            "--noise-json", str(tmp_path / "missing.json")]
    with pytest.raises(FileNotFoundError):
        rp.main(args)                                  # noise floor must exist
    args = args[:-2]
    assert rp.main(args + ["--noise-json", ""]) == 0
    assert json.loads(out.read_text())["git"] == FAKE_GIT
    with pytest.raises(FileExistsError):
        rp.main(args + ["--noise-json", ""])
    with pytest.raises(SystemExit):
        rp.main(["--out", str(tmp_path / "q.json")])   # no action requested


def test_git_state_of_this_repo():
    st = rp.git_state()
    assert len(st["commit"]) == 40 and all(c in "0123456789abcdef" for c in st["commit"])
    assert isinstance(st["dirty_paths"], list)
