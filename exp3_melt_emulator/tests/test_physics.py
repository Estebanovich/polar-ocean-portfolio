import dataclasses

import numpy as np
import pytest

from exp3_melt_emulator.data import cell_features
from exp3_melt_emulator.methods import REGISTRY, Method
from exp3_melt_emulator.methods.physics import QuadLocal, QuadNonLocal, LinearLocal
from exp3_melt_emulator.tests.synth import synthetic_runs


def _with_melt(runs, method):
    return [dataclasses.replace(r, melt=method.predict(r)) for r in runs]


@pytest.mark.parametrize("cls,K", [(QuadLocal, 3e-5), (LinearLocal, 0.7), (QuadNonLocal, 2e-4)])
def test_closed_form_recovers_K(cls, K):
    runs = synthetic_runs()
    gen = cls(); gen.coef = K
    train = _with_melt(runs, gen)
    m = cls(); m.fit(train)
    assert m.coef == pytest.approx(K, rel=1e-10)


def test_predict_nan_outside_mask():
    runs = synthetic_runs()
    for name in ("M1", "M2", "M3"):
        m = REGISTRY[name](); m.fit(runs)
        for r in runs:
            p = m.predict(r)
            assert p.shape == r.shelf_mask.shape
            assert np.isnan(p[~r.shelf_mask]).all()
            assert np.isfinite(p[r.shelf_mask]).all()


def test_registry_has_M1_M3():
    for k, name in (("M1", "M1"), ("M2", "M2"), ("M3", "M3")):
        m = REGISTRY[k]()
        assert isinstance(m, Method) and m.name == name


def test_predict_before_fit_and_empty_fit():
    r = synthetic_runs()[0]
    for k in ("M1", "M2", "M3"):
        with pytest.raises(RuntimeError):
            REGISTRY[k]().predict(r)
        with pytest.raises(ValueError):
            REGISTRY[k]().fit([])


def test_signed_refreezing():
    r = synthetic_runs()[0]
    cold = dataclasses.replace(r, T_profile=np.full_like(r.T_profile, -2.2))
    for cls in (LinearLocal, QuadLocal, QuadNonLocal):
        m = cls(); m.coef = 1.0
        p = m.predict(cold)[cold.shelf_mask]
        tf = cell_features(cold)["TF"][cold.shelf_mask]
        assert (tf < 0).all()
        assert (p < 0).all()  # signed: refreezing is not clipped (ruling R1)


def test_nonlocal_uses_mean():
    r = synthetic_runs()[0]
    sh = np.argwhere(r.shelf_mask)[0]
    i, j = sh
    # TF is interpolated from a shared profile, so double TF in one cell by patching
    # cell_features as seen by the methods.
    import exp3_melt_emulator.methods.physics as ph
    orig = ph.cell_features

    def patched(run):
        g = orig(run)
        TF = g["TF"].copy(); TF[i, j] *= 2
        g["TF"] = TF
        g["TF_mean"] = np.where(run.shelf_mask, np.nanmean(TF[run.shelf_mask]), np.nan)
        return g

    m3 = QuadNonLocal(); m3.coef = 1.0
    m2 = QuadLocal(); m2.coef = 1.0
    b3, b2 = m3.predict(r), m2.predict(r)
    ph.cell_features = patched
    try:
        a3, a2 = m3.predict(r), m2.predict(r)
    finally:
        ph.cell_features = orig
    other = r.shelf_mask.copy(); other[i, j] = False
    assert not np.allclose(a3[other], b3[other])
    np.testing.assert_array_equal(a2[other], b2[other])
