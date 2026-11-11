"""Pre-registered folds, checked against the real case list (no dataset needed)."""
from types import SimpleNamespace

import pytest

from exp3_melt_emulator.folds import Fold, drop_n1, make_folds

N1_FAIL = {"ens_T-0.5_zc250_s1.0", "ens_T-0.5_zc550_s1.0", "ens_T1.5_zc550_s1.0",
           "ens_T-0.42_zc402_s0.8", "ens_T-0.39_zc503_s1.2"}

GEOMETRY = [
    (0.75, 527, 0.8), (1.43, 295, 0.8), (-0.05, 328, 0.8), (-0.42, 402, 0.8),
    (0.34, 444, 0.8), (1.42, 481, 1.2), (0.02, 291, 1.2), (-0.39, 503, 1.2),
    (0.97, 386, 1.2), (0.66, 311, 1.2),
]


def _fixture():
    rows = [("ocean0_4km", "reference", "train", 1.0, 360.0, 1.0)]
    for T in (-0.5, 0.0, 0.5, 1.0, 1.5):
        for zc in (250, 400, 550):
            rows.append((f"ens_T{T}_zc{zc}_s1.0", "factorial", "train", T, float(zc), 1.0))
    for T, zc, s in GEOMETRY:
        rows.append((f"ens_T{T}_zc{zc}_s{s}", "geometry", "train", T, float(zc), s))
    for zc in (250, 400, 550):
        rows.append((f"ens_T2.0_zc{zc}_s1.0", "extrapolation", "test", 2.0, float(zc), 1.0))
    return [SimpleNamespace(case=c, group=g, split=sp, Tbot=T, z_c=zc, s=s,
                            n1_pass=c not in N1_FAIL)
            for c, g, sp, T, zc, s in rows]


@pytest.fixture
def runs():
    return _fixture()


@pytest.fixture
def folds(runs):
    return make_folds(runs)


def _by_name(folds):
    return {f.name: f for f in folds}


def test_fixture_size(runs):
    assert len(runs) == 29


def test_fold_counts(folds):
    counts = {}
    for f in folds:
        counts[f.kind] = counts.get(f.kind, 0) + 1
    assert counts == {"loro": 25, "zc_out": 3, "warm_out": 1, "profile_test": 1, "final": 1}


def test_reference_never_trains(folds):
    for f in folds:
        assert "ocean0_4km" not in f.train
        assert not any(c.startswith("ens_T2.0") for c in f.train)
    holders = [f for f in folds if "ocean0_4km" in f.test]
    assert [f.kind for f in holders] == ["profile_test"]


def test_loro_holds_out_each_training_run_once(folds):
    loro = [f for f in folds if f.kind == "loro"]
    held = [f.test[0] for f in loro]
    assert all(len(f.test) == 1 and len(f.train) == 24 for f in loro)
    assert len(set(held)) == 25
    for f in loro:
        assert f.test[0] not in f.train


def test_zc_out_partition(folds):
    by = _by_name(folds)
    tests = [set(by[n].test) for n in ("zc_out_250", "zc_out_400", "zc_out_550")]
    assert sum(len(t) for t in tests) == 25
    assert set().union(*tests) == set(by["profile_test"].train)
    for n in ("zc_out_250", "zc_out_400", "zc_out_550"):
        assert not set(by[n].train) & set(by[n].test)
        assert len(by[n].train) + len(by[n].test) == 25


def test_zc_out_550_excludes_527_and_503(folds):
    t = set(_by_name(folds)["zc_out_550"].test)
    assert "ens_T0.75_zc527_s0.8" in t
    assert "ens_T-0.39_zc503_s1.2" in t
    assert {f"ens_T{T}_zc550_s1.0" for T in (-0.5, 0.0, 0.5, 1.0, 1.5)} <= t
    assert "ens_T1.42_zc481_s1.2" in t  # 481 > 475


def test_zc_out_boundaries(folds):
    by = _by_name(folds)
    assert "ens_T-0.05_zc328_s0.8" not in by["zc_out_250"].test  # 328 >= 325
    assert "ens_T-0.05_zc328_s0.8" in by["zc_out_400"].test
    assert "ens_T0.66_zc311_s1.2" in by["zc_out_250"].test


def test_warm_out_test_set(folds):
    f = _by_name(folds)["warm_out"]
    assert set(f.test) == {
        "ens_T1.5_zc250_s1.0", "ens_T1.5_zc400_s1.0", "ens_T1.5_zc550_s1.0",
        "ens_T1.43_zc295_s0.8", "ens_T1.42_zc481_s1.2"}
    assert len(f.test) == 5
    train_T = {c: T for c, T in _case_T().items() if c in f.train}
    assert train_T and max(train_T.values()) <= 1.0


def _case_T():
    return {r.case: r.Tbot for r in _fixture()}


def test_profile_test_fold(folds):
    f = _by_name(folds)["profile_test"]
    assert len(f.train) == 25
    assert f.test == ("ocean0_4km",)


def test_final_test_is_T2(folds):
    f = _by_name(folds)["final"]
    assert set(f.test) == {f"ens_T2.0_zc{zc}_s1.0" for zc in (250, 400, 550)}
    assert len(f.train) == 25


def test_drop_n1(runs, folds):
    f = _by_name(folds)["final"]
    g = drop_n1(f, runs)
    assert isinstance(g, Fold)
    assert set(f.train) - set(g.train) == N1_FAIL
    assert g.test == f.test and g.name == f.name and g.kind == f.kind
    # a fold whose test already holds N1 runs keeps them in test
    h = drop_n1(_by_name(folds)["zc_out_550"], runs)
    assert not set(h.train) & N1_FAIL
    assert N1_FAIL & set(h.test)


def test_deterministic_order(runs):
    assert make_folds(runs) == make_folds(list(runs))
