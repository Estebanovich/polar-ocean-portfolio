"""Noise floor (ruling R4): pyocean-native, no torch / scikit-learn import."""
import dataclasses
import json
from pathlib import Path

import numpy as np
import pytest

from exp2_isomip_cavity.analysis import metrics as mt
from exp3_melt_emulator.analysis import noise_floor as nf
from exp3_melt_emulator.tests.synth import synthetic_runs

NY, NX = 5, 12


def _stub(shelf, melt_of_month, case="stub"):
    """RunData with 24 one-month records; record k carries SHIfwFlx such that the melt is
    ``melt_of_month(k + 1)`` (a constant or a map) everywhere."""
    recs = [mt.Record(iter=k + 1, t0=k * mt.MONTH_S, t1=(k + 1) * mt.MONTH_S,
                      stage_dir=Path("/nonexistent")) for k in range(24)]
    z2 = np.zeros((NY, NX))
    wet3 = np.zeros((2, NY, NX))
    wet3[:, :4, 3:] = 1.0
    draft = np.where(shelf, -200.0, 0.0)

    class Stub(mt.RunData):
        def field(self, kind, name, k):
            assert (kind, name) == ("shiDiag", "SHIfwFlx")
            m = np.broadcast_to(np.asarray(melt_of_month(k + 1), float), (NY, NX))
            return -m * mt.RHO_FW / mt.SECONDS_PER_YEAR        # inverse of melt_rate_m_per_a

    return Stub(run_dir=Path(f"/stub/{case}"), case={}, params={}, actual={}, stage_dirs=[],
                records=recs, nx=NX, ny=NY, nz=2, hfac=wet3, hfacW=wet3, hfacS=wet3,
                rA=np.full((NY, NX), 1.6e7), drF=np.ones(2), dyG=z2, dxG=z2, xc=z2, yc=z2,
                xg=z2, yg=z2, rc=np.array([-5.0, -15.0]), rf=np.array([0.0, -10.0, -20.0]),
                draft=draft, t_start=0.0, t_end=24 * mt.MONTH_S)


def test_noise_maps_windows(monkeypatch):
    shelf = synthetic_runs()[0].shelf_mask
    stub = _stub(shelf, lambda month: float(month))
    monkeypatch.setattr(nf.mt, "load_run", lambda d: stub)
    prev, last = nf.noise_maps("/stub/run")
    np.testing.assert_allclose(prev, np.mean(np.arange(13, 19)), rtol=1e-12)   # months 13-18
    np.testing.assert_allclose(last, np.mean(np.arange(19, 25)), rtol=1e-12)   # months 19-24


def test_noise_maps_rejects_wrong_record_count(monkeypatch):
    shelf = synthetic_runs()[0].shelf_mask
    stub = _stub(shelf, lambda month: 1.0)
    stub.records = stub.records[:23]
    monkeypatch.setattr(nf.mt, "load_run", lambda d: stub)
    with pytest.raises(ValueError):
        nf.noise_maps("/stub/run")


def test_area_rms_uses_mask_and_weights():
    mask = np.array([[True, True, False]])
    diff = np.array([[1.0, 3.0, 100.0]])
    area = np.array([[1.0, 3.0, 1.0]])
    assert nf.area_rms(diff, mask, area) == pytest.approx(np.sqrt((1 * 1 + 3 * 9) / 4))
    assert nf.area_rms(diff, mask, np.ones((1, 3))) == pytest.approx(np.sqrt(5.0))


def test_training_cases_excludes_reference_and_test():
    runs = synthetic_runs(n=6)
    runs[0] = dataclasses.replace(runs[0], group="reference", case="ocean0_4km")
    cases = nf.training_cases(runs)
    assert "ocean0_4km" not in cases
    assert all(r.split == "train" for r in runs if r.case in cases)
    assert len(cases) == 4                    # 6 - reference - 1 test run


def _stub_runs(runs, delta):
    """Stub RunData per run: months 13-18 = run melt + delta(run), months 19-24 = run melt."""
    def melt_of(r):
        base = np.where(r.shelf_mask, r.melt, 0.0)
        return lambda month: base + (delta(r) if month <= 18 else 0.0)
    return {r.case: _stub(r.shelf_mask, melt_of(r), r.case) for r in runs}


def test_build_rms_and_pooled_median(monkeypatch, tmp_path):
    runs = [dataclasses.replace(r, split="train") for r in synthetic_runs(n=5)]
    deltas = {r.case: 0.1 * (i + 1) for i, r in enumerate(runs)}     # uniform offsets
    stubs = _stub_runs(runs, lambda r: deltas[r.case])
    monkeypatch.setattr(nf.mt, "load_run", lambda d: stubs[Path(d).name.removesuffix("_r21")])
    res = nf.build(runs, tmp_path, "_r21")
    rms = {p["case"]: p["rms"] for p in res["summary"]["per_run"]}
    for c, d in deltas.items():
        assert rms[c] == pytest.approx(d, rel=1e-9)      # RMS of a uniform offset = |offset|
    assert res["summary"]["pooled_median_rms"] == pytest.approx(0.3, rel=1e-9)
    assert res["arrays"]["prev6"].shape == (5, NY, NX)


def test_build_refuses_dataset_mismatch(monkeypatch, tmp_path):
    runs = [dataclasses.replace(r, split="train") for r in synthetic_runs(n=3)]
    stubs = _stub_runs(runs, lambda r: 0.0)
    bad = runs[1].case
    shifted = _stub(runs[1].shelf_mask,
                    lambda month: np.where(runs[1].shelf_mask, runs[1].melt, 0.0) + 1.0, bad)
    stubs[bad] = shifted                                  # months 19-24 differ from the dataset
    monkeypatch.setattr(nf.mt, "load_run", lambda d: stubs[Path(d).name.removesuffix("_r21")])
    with pytest.raises(ValueError, match=bad):
        nf.build(runs, tmp_path, "_r21")


def test_write_outputs(monkeypatch, tmp_path):
    runs = [dataclasses.replace(r, split="train") for r in synthetic_runs(n=3)]
    stubs = _stub_runs(runs, lambda r: 0.2)
    monkeypatch.setattr(nf.mt, "load_run", lambda d: stubs[Path(d).name.removesuffix("_r21")])
    res = nf.build(runs, tmp_path, "_r21")
    npz, js = tmp_path / "out" / "noise_floor.npz", tmp_path / "noise_floor.json"
    s = nf.write_outputs(res, npz, js)
    rec = json.loads(js.read_text())
    import hashlib
    assert rec["sha256"] == s == hashlib.sha256(npz.read_bytes()).hexdigest()
    with np.load(npz, allow_pickle=False) as d:
        assert list(d["case"]) == [r.case for r in runs]
        np.testing.assert_allclose(d["rms"], 0.2, rtol=1e-9)


def test_module_imports_no_torch_or_sklearn():
    import subprocess
    import sys
    code = ("import sys; import exp3_melt_emulator.analysis.noise_floor; "
            "bad = [m for m in ('torch', 'sklearn') if m in sys.modules]; "
            "assert not bad, bad")
    root = Path(__file__).resolve().parents[2]
    subprocess.run([sys.executable, "-c", code], cwd=root, check=True)


def test_integrated_rel_diff(monkeypatch, tmp_path):
    runs = [dataclasses.replace(r, split="train") for r in synthetic_runs(n=3)]
    stubs = _stub_runs(runs, lambda r: 0.5)
    monkeypatch.setattr(nf.mt, "load_run", lambda d: stubs[Path(d).name.removesuffix("_r21")])
    res = nf.build(runs, tmp_path, "_r21")
    for p, r in zip(res["summary"]["per_run"], runs):
        mean = r.melt[r.shelf_mask].mean()
        assert p["integrated_rel_diff"] == pytest.approx(-0.5 / mean, rel=1e-9)


def test_build_refuses_nonpositive_mean_melt(monkeypatch, tmp_path):
    runs = [dataclasses.replace(r, split="train") for r in synthetic_runs(n=2)]
    r0 = runs[0]
    runs[0] = dataclasses.replace(r0, melt=np.where(r0.shelf_mask, -1.0, r0.melt))
    stubs = _stub_runs(runs, lambda r: 0.0)
    monkeypatch.setattr(nf.mt, "load_run", lambda d: stubs[Path(d).name.removesuffix("_r21")])
    with pytest.raises(ValueError, match="relative noise undefined"):
        nf.build(runs, tmp_path, "_r21")
