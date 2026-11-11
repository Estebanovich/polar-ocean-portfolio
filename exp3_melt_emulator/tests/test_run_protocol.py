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
    assert cfg["methods"]["M7"]["tf_kind_folds"] is None          # PD0 decision 3: full LORO
    assert cfg["methods"]["M5"]["var_fraction"] == 0.99 and cfg["methods"]["M5"]["n_restarts"] == 5
    assert cfg["methods"]["M4"]["degrees"] == [1, 2]
    from exp3_melt_emulator.methods.nets import HYBRID_G_FEATURES, N0, N_VAL, NET_FEATURES
    m6, m7 = cfg["methods"]["M6"], cfg["methods"]["M7"]
    assert m6["input_features"] == list(NET_FEATURES)
    assert m7["input_features"] == list(HYBRID_G_FEATURES)
    assert m6["n_val"] == m7["n_val"] == N_VAL == 3
    assert "default_rng([split_seed, member_seed])" in m6["validation_rule"]
    assert m7["n0"] == N0 == 2.0 and "M2" in m7["K0_rule"] and "M3" in m7["K0_rule"]
    assert "stopgrad" in m7["variance"] and "R11" in m7["variance"]
    assert "softplus" in m6["variance"]
    assert "tie -> local" in m7["tf_kind_rule"] and "split_seed + 1" in m7["tf_kind_rule"]
    assert cfg["methods"]["M1"]["input_features"] == ["TF"]
    assert cfg["methods"]["M2"]["input_features"] == ["TF"]
    assert cfg["methods"]["M3"]["input_features"] == ["TF", "TF_mean"]
    for mid in ("M4", "M5"):
        assert cfg["methods"][mid]["input_features"] == ["Tbot", "z_c", "s"]
    assert cfg["methods"]["M5"]["white_bounds"] == [1e-8, 1e-2]   # PD0 decision 4
    assert cfg["methods"]["M5"]["white_init"] == 1e-3
    assert "8 h" in cfg["thresholds"]["stop_rule"]
    assert cfg["profile_test_methods"] == ["M0", "M1", "M2", "M3", "M6", "M7"]
    assert "drop_n1" in cfg["thresholds"]["H3b"]["cases"]
    assert cfg["thresholds"] == json.loads(json.dumps(rp.PD0_THRESHOLDS))
    assert cfg["thresholds"]["H3b"]["max_ratio_M7_over_M6"] == 1.10
    assert cfg["thresholds"]["gl_frac"] == 0.10


def test_freeze_refuses_to_overwrite(tmp_path, dataset):
    out = _freeze(tmp_path, dataset)
    before = out.read_bytes()
    with pytest.raises(FileExistsError):
        _freeze(tmp_path, dataset, out=out)
    assert out.read_bytes() == before


def test_freeze_refuses_dirty_tree_outside_notes(tmp_path, dataset):
    _, path, sha = dataset
    out = tmp_path / "c.json"
    dirty = {"commit": "0" * 40, "dirty_paths": ["exp3_melt_emulator/notes/pd0_preregistro.md",
                                                  "exp3_melt_emulator/methods/nets.py"]}
    with pytest.raises(RuntimeError, match="methods/nets.py"):
        rp.freeze_config(out, dataset_path=path, dataset_sha=sha, noise_json=None,
                         git_state=lambda: dict(dirty))
    assert not out.exists()
    notes_only = {"commit": "0" * 40, "dirty_paths": ["exp3_melt_emulator/notes/pd0_preregistro.md"]}
    rp.freeze_config(out, dataset_path=path, dataset_sha=sha, noise_json=None,
                     git_state=lambda: dict(notes_only))
    assert json.loads(out.read_text())["git"] == notes_only


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


@pytest.mark.skipif(not (rp.REPO / ".git").exists(), reason="not a git checkout")
def test_git_state_of_this_repo():
    st = rp.git_state()
    assert len(st["commit"]) == 40 and all(c in "0123456789abcdef" for c in st["commit"])
    assert isinstance(st["dirty_paths"], list)


# --- Task 8: --cv and --final ------------------------------------------------------------
import dataclasses  # noqa: E402

import numpy as np  # noqa: E402

from exp3_melt_emulator.metrics import integrated_rel_error, rmse_gl, rmse_map  # noqa: E402

PROTO_CASES = {   # i: (Tbot, z_c, group); index 6 is the split == "test" run (synth: last index)
    0: (-0.5, 250.0, "factorial"), 1: (0.0, 400.0, "factorial"), 2: (0.5, 550.0, "factorial"),
    3: (1.5, 400.0, "factorial"), 4: (1.0, 300.0, "geometry"), 5: (0.0, 400.0, "reference"),
    6: (2.0, 400.0, "factorial"),
}


def _proto_runs():
    """7 synthetic runs giving non-empty loro (5), zc_out (3), warm_out, profile_test and
    final folds; run 2 fails N1, run 5 is the reference, run 6 the split == "test" run."""
    runs = synthetic_runs(n=7)
    out = []
    for i, r in enumerate(runs):
        tb, zc, grp = PROTO_CASES[i]
        out.append(dataclasses.replace(r, Tbot=tb, z_c=zc, group=grp,
                                       split="reference" if grp == "reference" else r.split))
    return out


@pytest.fixture(scope="module")
def proto(tmp_path_factory):
    """Synthetic dataset + frozen config; EXP3_TEST_MODE=1 for the module (tiny methods)."""
    d = tmp_path_factory.mktemp("proto")
    runs = _proto_runs()
    ds = d / "ds.npz"
    write_npz(runs, ds)
    sha = hashlib.sha256(ds.read_bytes()).hexdigest()
    cfg = d / "protocol_config.json"
    rp.freeze_config(cfg, dataset_path=ds, dataset_sha=sha, noise_json=None,
                     git_state=lambda: dict(FAKE_GIT))
    csha = hashlib.sha256(cfg.read_bytes()).hexdigest()
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("EXP3_TEST_MODE", "1")
        yield {"dir": d, "runs": runs, "ds": ds, "cfg": cfg, "csha": csha}


def _cv(proto, root, **kw):
    kw.setdefault("workers", 2)
    return rp.run_cv(proto["cfg"], proto["csha"], units_dir=root / "units",
                     progress_log=root / "progress.log", out=root / "cv_results.json",
                     git_state=lambda: dict(FAKE_GIT), **kw)


@pytest.fixture(scope="module")
def cv_done(proto):
    root = proto["dir"] / "cv"
    status = _cv(proto, root)
    return root, status


def _expected_units(proto):
    cfg = json.loads(proto["cfg"].read_text())
    exp = []
    for variant, folds in (("main", cfg["folds"]), ("n1", cfg["folds_sensitivity_n1"])):
        for f in folds:
            if f["kind"] == "final":
                continue
            for mid in sorted(cfg["methods"]):
                if f["kind"] == "profile_test" and mid not in cfg["profile_test_methods"]:
                    continue
                exp.append((variant, f["name"], mid, tuple(f["test"])))
    return exp


def test_cv_results_schema(proto, cv_done):
    root, status = cv_done
    assert status["status"] == "complete"
    res = json.loads((root / "cv_results.json").read_text())
    assert res["config_sha256"] == proto["csha"]
    assert res["dataset_sha256"] == hashlib.sha256(proto["ds"].read_bytes()).hexdigest()
    assert res["test_overrides"] == rp.TEST_OVERRIDES
    exp = _expected_units(proto)
    want = {(v, f, m, c) for v, f, m, cases in exp for c in cases}
    got = [(e["variant"], e["fold"], e["method"], e["case"]) for e in res["entries"]]
    assert len(got) == len(set(got)) and set(got) == want          # one entry each
    assert not any(f == "final" for _, f, _, _ in got)
    test_cases = {r.case for r in proto["runs"] if r.split == "test"}
    assert not {c for *_, c in got} & test_cases                   # never the final test runs
    assert {u["unit"] for u in res["units"]} == {f"{v}__{f}__{m}" for v, f, m, _ in exp}
    assert {"M4", "M5"}.isdisjoint({m for _, f, m, _ in got if f == "profile_test"})
    runs = {r.case: r for r in proto["runs"]}
    for e in res["entries"]:
        mt = e["metrics"]
        assert set(mt) >= {"rmse_map", "integrated_rel_error", "rmse_gl", "integrated_pred",
                           "integrated_true"}
        assert all(np.isfinite(v) for v in mt.values())
        assert set(e["skill"]) == {"rmse_map", "rmse_gl", "abs_integrated_rel_error"}
        assert e["n1_pass"] == runs[e["case"]].n1_pass
        assert e["label"] == ("sensitivity" if e["variant"] == "n1" else "main")
        if e["method"] == "M0":
            assert all(v == 0.0 for v in e["skill"].values())
        if e["method"] in ("M5", "M6", "M7"):
            assert set(e["prob"]) == {"coverage_68", "coverage_95", "crps", "mean_sigma"}
            assert 0 <= e["prob"]["coverage_68"] <= e["prob"]["coverage_95"] <= 1
        else:
            assert e["prob"] is None
        if e["method"] in ("M6", "M7"):
            mb = e["members"]
            assert len(mb["integrated"]) == rp.TEST_OVERRIDES[e["method"]]["n_members"]
            assert mb["integrated_mixture"] == pytest.approx(np.mean(mb["integrated"]))
            assert mb["integrated_mixture"] == pytest.approx(mt["integrated_pred"])
            assert mb["spread_rel"] == pytest.approx(
                np.std(mb["integrated"]) / abs(mb["integrated_mixture"]))
            assert mb["mean_aleatoric"] >= 0 and mb["mean_epistemic"] >= 0
        else:
            assert e["members"] is None
    infos = {u["unit"]: u for u in res["units"]}
    m7 = infos["main__loro_synth_0__M7"]["info"]
    assert m7["tf_kind"] in ("local", "nonlocal") and "members" in m7
    assert {"K", "n", "epochs_run", "hit_max_epochs"} <= set(m7["members"][0])
    assert "gp" in infos["main__loro_synth_0__M5"]["info"]
    assert infos["main__loro_synth_0__M4"]["info"]["degree"] in (1, 2)
    u = infos["n1__loro_synth_0__M1"]
    assert "synth_2" not in u["train"] and u["wall_s"] > 0


def test_cv_metrics_recomputable(proto, cv_done):
    """Spot check: the M1 entry equals metrics recomputed from an in-process fit."""
    root, _ = cv_done
    res = json.loads((root / "cv_results.json").read_text())
    runs = {r.case: r for r in proto["runs"]}
    e = next(e for e in res["entries"] if (e["variant"], e["fold"], e["method"])
             == ("main", "zc_out_400", "M1") and e["case"] == "synth_3")
    u = next(u for u in res["units"] if u["unit"] == "main__zc_out_400__M1")
    m = REGISTRY["M1"]()
    m.fit([runs[c] for c in u["train"]])
    p = m.predict(runs["synth_3"])
    assert e["metrics"]["rmse_map"] == pytest.approx(rmse_map(p, runs["synth_3"]), rel=1e-12)
    assert e["metrics"]["rmse_gl"] == pytest.approx(rmse_gl(p, runs["synth_3"]), rel=1e-12)
    assert e["metrics"]["integrated_rel_error"] == pytest.approx(
        integrated_rel_error(p, runs["synth_3"]), rel=1e-12)


def test_cv_progress_log(cv_done):
    root, _ = cv_done
    log = (root / "progress.log").read_text().splitlines()
    done = [ln for ln in log if " DONE " in ln]
    assert len(done) == len(list((root / "units").glob("*.json")))
    assert all("wall" in ln for ln in done)
    assert any("ESTIMATE" in ln for ln in log)
    first = [ln for ln in done[:2]]
    assert all("__M7 " in ln for ln in first)                    # M7 units scheduled first


def test_resume_skips_completed_units(proto, cv_done):
    root, _ = cv_done
    units = sorted((root / "units").glob("*.json"))
    before = {p.name: p.read_bytes() for p in units}
    redo = root / "units" / "main__loro_synth_1__M2.json"
    redo.unlink()
    status = _cv(proto, root)
    assert status["status"] == "complete" and status["ran"] == ["main__loro_synth_1__M2"]
    after = {p.name: p.read_bytes() for p in (root / "units").glob("*.json")}
    assert set(after) == set(before)
    assert all(after[k] == v for k, v in before.items() if k != redo.name)
    status = _cv(proto, root)
    assert status["ran"] == []
    assert not list((root / "units").glob("*.tmp*"))


def test_resume_refuses_stale_unit(proto, tmp_path):
    units = tmp_path / "units"
    units.mkdir()
    (units / "main__loro_synth_1__M2.json").write_text(json.dumps(
        {"unit": "main__loro_synth_1__M2", "config_sha256": "f" * 64,
         "test_overrides": rp.TEST_OVERRIDES}))
    with pytest.raises(rp.ProtocolError, match="config"):
        _cv(proto, tmp_path)


def test_stop_rule_triggers(proto, tmp_path):
    calls = []

    def estimate(**kw):
        calls.append(kw)
        return {"estimate_h": 99.0}
    status = _cv(proto, tmp_path, estimate_fn=estimate)
    assert status["status"] == "stopped" and status["estimate"]["estimate_h"] == 99.0
    assert len(calls) == 1 and len(calls[0]["completed"]) >= 2      # after the first batch
    assert not (tmp_path / "cv_results.json").exists()
    n = len(list((tmp_path / "units").glob("*.json")))
    assert 2 <= n <= 3                                             # no new units after the check
    log = (tmp_path / "progress.log").read_text()
    assert "STOP" in log and "99.0" in log
    # below the limit: runs on (resume from the stopped state)
    status = _cv(proto, tmp_path, estimate_fn=lambda **kw: {"estimate_h": 1.0})
    assert status["status"] == "complete"


def test_default_estimate_lpt():
    pend = [{"method": "M7", "n_train": 24}] * 4 + [{"method": "M1", "n_train": 24}]
    done = [({"method": "M7", "n_train": 24}, 490.0)]                 # 10 s per inner fit
    est = rp.estimate_wall(completed=done, in_flight=[({"method": "M7", "n_train": 24}, 100.0)],
                           pending=pend, elapsed_s=500.0, workers=2)
    # in flight: 390 s left; pending M7 490 s each, LPT on 2 workers -> 390+490+490 = 1370
    assert est["estimate_h"] == pytest.approx((500.0 + 1370.0) / 3600, rel=1e-9)


def test_config_mismatch_refuses(proto, tmp_path):
    cfg = proto["cfg"]
    with pytest.raises(rp.ProtocolError, match="sha256"):
        rp.load_protocol(cfg, "0" * 64)
    # an edited config whose sha is "expected" still fails if hparams differ from the code
    rec = json.loads(cfg.read_text())
    rec["methods"]["M6"]["patience"] = 1
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(rec))
    with pytest.raises(rp.ProtocolError, match="M6"):
        rp.load_protocol(bad, hashlib.sha256(bad.read_bytes()).hexdigest())
    # folds that differ from make_folds
    rec = json.loads(cfg.read_text())
    rec["folds"][0]["train"] = rec["folds"][0]["train"][:-1]
    bad.write_text(json.dumps(rec))
    with pytest.raises(rp.ProtocolError, match="folds"):
        rp.load_protocol(bad, hashlib.sha256(bad.read_bytes()).hexdigest())
    # dataset sha
    rec = json.loads(cfg.read_text())
    rec["dataset"]["sha256"] = "e" * 64
    bad.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match="sha256"):
        rp.load_protocol(bad, hashlib.sha256(bad.read_bytes()).hexdigest())
    for fn in (rp.run_cv, rp.run_final):
        with pytest.raises(rp.ProtocolError, match="sha256"):
            fn(cfg, "0" * 64, git_state=lambda: dict(FAKE_GIT), **(
                {"units_dir": tmp_path / "u", "progress_log": tmp_path / "l",
                 "out": tmp_path / "o.json"} if fn is rp.run_cv else
                {"out": tmp_path / "f.json", "maps_out": tmp_path / "m.npz",
                 "marker": tmp_path / "started.json", "progress_log": tmp_path / "l"}))
    assert not (tmp_path / "f.json").exists() and not (tmp_path / "started.json").exists()


def test_cli_uses_frozen_sha(monkeypatch, tmp_path):
    assert rp.CONFIG_SHA256 == "56a9ffe5b57fd69da0462e63cb79cb1cbc7538e8444946e9bde5101f9a462486"
    seen = {}
    monkeypatch.setattr(rp, "run_cv", lambda cfg, sha, **kw: seen.update(sha=sha, **kw) or
                        {"status": "complete", "ran": []})
    assert rp.main(["--cv"]) == 0 and seen["sha"] == rp.CONFIG_SHA256
    assert seen["wall_limit_hours"] is None and seen["amendment"] is None   # frozen 8 h rule
    assert rp.main(["--cv", "--wall-limit-hours", "20", "--amendment", "E1"]) == 0
    assert seen["wall_limit_hours"] == 20.0 and seen["amendment"] == "E1"


def test_test_mode_refused_on_real_dataset(monkeypatch):
    from exp3_melt_emulator.data import DATASET_SHA256
    runs = _proto_runs()
    monkeypatch.setenv("EXP3_TEST_MODE", "1")
    with pytest.raises(rp.ProtocolError, match="EXP3_TEST_MODE"):
        rp.test_overrides({"dataset": {"sha256": DATASET_SHA256}}, runs)
    real_like = [dataclasses.replace(r, case="ens_" + r.case) for r in runs]
    with pytest.raises(rp.ProtocolError, match="EXP3_TEST_MODE"):
        rp.test_overrides({"dataset": {"sha256": "a" * 64}}, real_like)
    assert rp.test_overrides({"dataset": {"sha256": "a" * 64}}, runs) == rp.TEST_OVERRIDES
    monkeypatch.delenv("EXP3_TEST_MODE")
    assert rp.test_overrides({"dataset": {"sha256": "a" * 64}}, runs) == {}


def test_final_refuses_second_run(proto, tmp_path):
    kw = dict(out=tmp_path / "final_test.json", maps_out=tmp_path / "final_maps.npz",
              marker=tmp_path / "final_started.json", progress_log=tmp_path / "final.log",
              git_state=lambda: dict(FAKE_GIT), workers=2)
    rec = rp.run_final(proto["cfg"], proto["csha"], **kw)
    res = json.loads((tmp_path / "final_test.json").read_text())
    assert res == json.loads(json.dumps(rec))
    assert res["maps_sha256"] == hashlib.sha256((tmp_path / "final_maps.npz").read_bytes()).hexdigest()
    got = {(e["method"], e["case"]) for e in res["entries"]}
    assert got == {(m, "synth_6") for m in REGISTRY}
    assert all(e["fold"] == "final" for e in res["entries"])
    with np.load(tmp_path / "final_maps.npz", allow_pickle=False) as z:
        assert "M7__synth_6__mean" in z.files and "M7__synth_6__sigma" in z.files
        assert z["M7__synth_6__members_mu"].shape[0] == rp.TEST_OVERRIDES["M7"]["n_members"]
        assert "M2__synth_0__mean_in_sample" in z.files
    before = (tmp_path / "final_test.json").read_bytes()
    with pytest.raises(FileExistsError):
        rp.run_final(proto["cfg"], proto["csha"], **kw)
    assert (tmp_path / "final_test.json").read_bytes() == before
    # a started-but-unfinished run (marker present) also refuses
    (tmp_path / "final_test.json").unlink()
    with pytest.raises(rp.ProtocolError, match="started"):
        rp.run_final(proto["cfg"], proto["csha"], **kw)
    assert not (tmp_path / "final_test.json").exists()


# --- Enmienda 1 (PD0): wall-clock budget override, main-first scheduling -------------------
def test_cv_units_main_before_n1(proto):
    cfg = json.loads(proto["cfg"].read_text())
    units = rp.cv_units(cfg, {r.case: r for r in proto["runs"]})
    variants = [u["variant"] for u in units]
    assert variants == sorted(variants, key=lambda v: v != "main")    # all main, then all n1
    for v in ("main", "n1"):
        grp = [u for u in units if u["variant"] == v]
        meth = [u["method"] for u in grp]
        assert meth == sorted(meth, key=lambda m: m != "M7")           # M7 first within each
        for sub in ([u for u in grp if u["method"] == "M7"], [u for u in grp if u["method"] != "M7"]):
            costs = [rp.prior_seconds(u) for u in sub]
            assert costs == sorted(costs, reverse=True)                # longest first


def test_wall_limit_requires_amendment(proto, tmp_path):
    for kw in ({"wall_limit_hours": 20.0}, {"amendment": "text"}):
        with pytest.raises(rp.ProtocolError, match="together"):
            _cv(proto, tmp_path, **kw)
    assert not (tmp_path / "units").exists()
    with pytest.raises(rp.ProtocolError, match="amendment"):
        _cv(proto, tmp_path, wall_limit_hours=20.0, amendment="  ")
    with pytest.raises(SystemExit):
        rp.main(["--cv", "--wall-limit-hours", "20"])
    with pytest.raises(SystemExit):
        rp.main(["--cv", "--amendment", "x"])
    with pytest.raises(SystemExit):
        rp.main(["--final", "--wall-limit-hours", "20", "--amendment", "x"])


def test_wall_limit_amendment_recorded(proto, tmp_path):
    text = "Enmienda de prueba: presupuesto 20 h"
    # 9 h: above the frozen 8 h limit, so the frozen rule stops ...
    st = _cv(proto, tmp_path / "a", estimate_fn=lambda **kw: {"estimate_h": 9.0})
    assert st["status"] == "stopped"
    # ... and under the amendment it runs on; the earlier units are resumed
    st = _cv(proto, tmp_path / "a", estimate_fn=lambda **kw: {"estimate_h": 9.0},
             wall_limit_hours=20.0, amendment=text)
    assert st["status"] == "complete"
    root = tmp_path / "a"
    log = (root / "progress.log").read_text()
    assert "AMENDMENT" in log and text in log and "frozen 8 h" in log and "limit 20 h" in log
    res = json.loads((root / "cv_results.json").read_text())
    (am,) = res["amendments"]
    assert am["text"] == text and am["wall_limit_hours"] == 20.0
    assert am["frozen_wall_limit_hours"] == rp.STOP_HOURS == 8.0 and am["timestamp"]
    recs = [json.loads(p.read_text()) for p in (root / "units").glob("*.json")]
    under = [r for r in recs if r.get("amendment")]
    assert under and all(r["amendment"] == am for r in under)
    assert len(under) == len(st["ran"])                       # every unit written under it
    assert len(recs) - len(under) >= 2                        # pre-amendment units kept


def test_default_has_no_amendment(cv_done):
    root, _ = cv_done
    res = json.loads((root / "cv_results.json").read_text())
    assert res["amendments"] == [] and res["wall_limit_hours"] == 8.0
    assert all(json.loads(p.read_text())["amendment"] is None
               for p in (root / "units").glob("*.json"))
