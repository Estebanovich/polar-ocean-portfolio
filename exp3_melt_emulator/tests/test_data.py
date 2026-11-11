import dataclasses
from pathlib import Path

import numpy as np
import pytest

from exp2_isomip_cavity.theory.melt import LAMBDA1, LAMBDA2, LAMBDA3, PA_TO_DBAR
from exp3_melt_emulator.data import (
    DATASET_PATH, FEATURES, cell_features, load_runs, pressure_dbar, union_mask,
)
from exp3_melt_emulator.tests.synth import synthetic_runs, write_npz


def test_sha_mismatch_raises(tmp_path):
    p = tmp_path / "copy.npz"
    write_npz(synthetic_runs(), p)
    with pytest.raises(ValueError):
        load_runs(p, sha256="0" * 64)


def test_load_roundtrip_with_own_sha(tmp_path):
    import hashlib
    runs = synthetic_runs()
    p = tmp_path / "copy.npz"
    write_npz(runs, p)
    got = load_runs(p, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    assert [r.case for r in got] == [r.case for r in runs]
    assert np.array_equal(got[1].shelf_mask, runs[1].shelf_mask)
    assert np.allclose(got[3].melt, runs[3].melt, equal_nan=True)


def test_cell_features_keys_and_mask(synth_runs):
    for r in synth_runs:
        f = cell_features(r)
        assert set(f) == set(FEATURES)
        for k, a in f.items():
            assert a.shape == r.shelf_mask.shape
            assert np.array_equal(np.isnan(a), ~r.shelf_mask), k


def test_TF_matches_liquidus(synth_runs):
    r = synth_runs[0]
    nz = len(r.z)
    r = dataclasses.replace(r, T_profile=np.full(nz, 1.0), S_profile=np.full(nz, 34.5),
                            draft=np.where(r.shelf_mask, -500.0, r.draft))
    p = pressure_dbar(-500.0)
    expected = 1.0 - (LAMBDA1 * 34.5 + LAMBDA2 + LAMBDA3 / PA_TO_DBAR * p)
    tf = cell_features(r)["TF"]
    assert np.allclose(tf[r.shelf_mask], expected, rtol=0, atol=1e-12)
    assert np.allclose(cell_features(r)["TF_mean"][r.shelf_mask], expected, rtol=0, atol=1e-12)


def test_T_S_interpolated_at_draft(synth_runs):
    r = synth_runs[0]
    f = cell_features(r)
    m = r.shelf_mask
    assert np.allclose(f["T"][m], np.interp(r.draft[m], r.z[::-1], r.T_profile[::-1]))
    assert np.allclose(f["S"][m], np.interp(r.draft[m], r.z[::-1], r.S_profile[::-1]))


def test_distances_ignore_wall_row(synth_runs):
    r = synth_runs[0]            # shelf columns 3-7, grounded columns 0-2, open from column 8
    f = cell_features(r)
    dx = 4000.0
    for c in range(3, 8):
        # the wall row (last row) is only dx away from row 3, but it is not grounded ice
        assert np.allclose(f["dist_gl"][:4, c], (c - 2) * dx)
    assert np.allclose(f["dist_front"][:4, 7], dx)
    assert np.allclose(f["dist_front"][:4, 3], 5 * dx)


def test_slope_and_scalars(synth_runs):
    r = synth_runs[0]
    f = cell_features(r)
    m = r.shelf_mask
    assert np.all(np.isfinite(f["slope"][m])) and np.all(f["slope"][m] > 0)
    # draft = -(100+60 col) - 5 row: |grad| = hypot(60/dx, 5/dy)
    assert np.allclose(f["slope"][m], np.hypot(60.0 / 4000.0, 5.0 / 4000.0))
    assert np.allclose(f["draft_minus_zc"][m], -r.draft[m] - r.z_c)
    assert np.all(f["s"][m] == r.s)


def test_synth_contents(synth_runs):
    assert len(synth_runs) == 6
    assert synth_runs[1].shelf_mask.sum() < synth_runs[0].shelf_mask.sum()
    assert union_mask(synth_runs).sum() == synth_runs[0].shelf_mask.sum()
    assert (synth_runs[3].melt[synth_runs[3].shelf_mask] < 0).any()
    assert [r.case for r in synth_runs if r.split == "test"] == ["synth_5"]


@pytest.mark.skipif(not Path(DATASET_PATH).exists(), reason="LaCie not mounted")
def test_real_dataset_shape():
    runs = load_runs()
    assert len(runs) == 29
    assert sorted(r.case for r in runs if r.split == "test") == [
        "ens_T2.0_zc250_s1.0", "ens_T2.0_zc400_s1.0", "ens_T2.0_zc550_s1.0"]
    counts = {s: {int(r.shelf_mask.sum()) for r in runs if r.s == s} for s in (0.8, 1.0, 1.2)}
    assert counts == {0.8: {726}, 1.0: {696}, 1.2: {620}}
    for r in runs:
        f = cell_features(r)
        assert np.all(np.isfinite(f["TF"][r.shelf_mask]))
