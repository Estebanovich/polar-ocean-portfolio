"""Tests for exp4_arctic_bayiw.analysis.posthoc_reference (Task 7 fix round 1, POST HOC): the KF
baroclinic pressure work with p' built from rho' relative to the evolving box-mean density
rho_box(z, t) instead of the fixed initial rho_b (method fixed in notes/pd3_hipotesis.md before
computing)."""
from __future__ import annotations

import numpy as np
import pytest

from exp4_arctic_bayiw.analysis import partition as pa, posthoc_reference as ph
from exp4_arctic_bayiw.tests.test_partition import BANDS, DT, _synthetic


def _drifted(ds, rho_b, amp=0.05):
    """The same dataset plus a horizontally uniform density drift delta(z, t) in every wet cell."""
    t = ds["time"].values
    nz = ds.sizes["k"]
    delta = amp * (t / t[-1])[:, None] * np.linspace(1.0, -0.5, nz)[None]          # (nt, nz)
    out = ds.copy()
    wet = ds["hFacC"].values > 0
    out["RHO"] = (ds["RHO"].dims, ds["RHO"].values + delta[:, :, None, None] * wet[None])
    return out


def test_box_mean_profile_is_the_volume_weighted_mean_of_the_wet_box_cells():
    ds, mask, groups, rho_b, tau = _synthetic()
    prof = ph.box_mean_profile(ds, mask, rho_b)
    R = ds["RHO"].values
    vol = ds["drF"].values[:, None, None] * ds["hFacC"].values * ds["rA"].values[None]
    w = vol * mask[None]
    ref = (R * w[None]).sum((2, 3)) / w.sum((1, 2))[None]
    np.testing.assert_allclose(prof, ref, rtol=1e-13)


def test_fixed_reference_reproduces_the_face_series():
    """With rho_ref(z, t) = rho_b for every record the moving pass is partition.face_series."""
    ds, mask, groups, rho_b, tau = _synthetic()
    nt = ds.sizes["time"]
    m = ph.face_pressure(ds, mask, groups, np.tile(rho_b, (nt, 1)), progress=False)
    f = pa.face_series(ds, mask, groups, rho_b=rho_b, progress=False)
    np.testing.assert_allclose(m["up"], f["up"], rtol=0, atol=1e-15)
    np.testing.assert_allclose(m["pp"], f["pp"], rtol=1e-13, atol=1e-9)


def test_moving_reference_removes_a_uniform_background_drift():
    """A horizontally uniform drift delta(z, t) changes the fixed-reference F_bc but not the
    moving-reference one."""
    ds, mask, groups, rho_b, tau = _synthetic()
    dd = _drifted(ds, rho_b)
    a = ph.moving_pressure_work(ds, mask, groups, rho_b, progress=False)
    b = ph.moving_pressure_work(dd, mask, groups, rho_b, progress=False)
    for g in groups:
        np.testing.assert_allclose(b["series_moving"][g], a["series_moving"][g], rtol=1e-9,
                                   atol=1e-9 * np.abs(a["series_moving"][g]).max())
    fa = pa.section_series(pa.face_series(ds, mask, groups, rho_b=rho_b, progress=False))
    fb = pa.section_series(pa.face_series(dd, mask, groups, rho_b=rho_b, progress=False))
    g = list(groups)[0]
    assert np.abs(fb[g]["F_bc"] - fa[g]["F_bc"]).max() > 1e-3 * np.abs(fa[g]["F_bc"]).max()


def test_summary_band_split_adds_up_and_eta_uses_the_budget_denominator():
    ds, mask, groups, rho_b, tau = _synthetic()
    r = ph.moving_pressure_work(ds, mask, groups, rho_b, progress=False)
    s = ph.summarize(r, BANDS, P_J=100.0)
    for g in list(groups) + ["total"]:
        assert sum(s[g]["band_out_J"].values()) == pytest.approx(s[g]["dt_sum_out_J"], rel=1e-9)
    assert s["eta"] == pytest.approx(s["total"]["out_J"] / 100.0)
    assert s["label"].startswith("POST HOC")
    res = sum(s["total"]["band_out_J"][b] for b in ("mid", "high"))
    assert s["eta_resolved"] == pytest.approx(res / 100.0)


def test_moving_pass_on_the_local_store_equals_the_full_dataset(tmp_path):
    ds, mask, groups, rho_b, tau = _synthetic()
    st = pa.build_store(ds, mask, tmp_path / "s", progress=False)
    dl, ml, gl, _ = pa.store_dataset(ds, st, mask, groups)
    a = ph.moving_pressure_work(dl, ml, gl, rho_b, progress=False)
    b = ph.moving_pressure_work(ds, mask, groups, rho_b, progress=False)
    for g in groups:
        np.testing.assert_array_equal(a["series_moving"][g], b["series_moving"][g])
