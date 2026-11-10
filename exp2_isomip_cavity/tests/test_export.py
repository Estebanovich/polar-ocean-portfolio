"""Tests for analysis/export_exp3.py: the Exp. 2 -> Exp. 3 ensemble dataset (Task 7).

Synthetic runs come from test_metrics.make_run (8 x 5 x 4 grid, shelf x < 640 km, draft
-40 m). One column is made dry (grounded ice) so the NaN rule can be checked. Each run
carries a uniform melt rate m under the ice shelf, so its exported melt map is m on the ice shelf, 0 on open
ocean and NaN on dry cells, exactly.
"""
from __future__ import annotations

import hashlib

import numpy as np
import pytest

from exp2_isomip_cavity.analysis import export_exp3 as ex
from exp2_isomip_cavity.gendata import profiles as pr
from test_metrics import NX, NY, NZ, default_geometry, make_run

DRY = (2, 1)                       # (j, i) of the dry (grounded) column


def _geometry():
    draft, _, hfac = default_geometry()
    hfac = hfac.copy()
    hfac[:, DRY[0], DRY[1]] = 0.0              # grounded: no wet cell in this column
    return draft, hfac


def _fw(m):
    """SHIfwFlx [kg m^-2 s^-1] giving melt m [m/a] (m = -fw / 1000 * 31 536 000) under the
    ice shelf and 0 elsewhere, as the model writes it."""
    draft, hfac = _geometry()
    shelf = (hfac > 0).any(0) & (draft < 0)
    return np.where(shelf, -m * 1000.0 / 31536000.0, 0.0)


def write_case(root, name, Tbot, z_c, s, m, *, group="factorial", status="OK", d=300.0,
               drift=0.0):
    """Run with melt m in months 19-24 and m (1 - drift) in months 13-18, so N1 = drift."""
    draft, hfac = _geometry()
    fw = lambda k: {"SHIfwFlx": _fw(m if k >= 18 else m * (1.0 - drift))}   # noqa: E731
    make_run(root / f"{name}_r21", nrec=24, status=status, hfac=hfac, draft=draft, shi=fw,
             case_extra={"params": {"case": name, "Tbot": Tbot, "z_c": z_c, "s": s, "d": d,
                                    "init_profile": "ramp", "restore_profile": "ramp",
                                    "ensemble_group": group},
                         "actual": {"Tbot": Tbot, "z_c": z_c, "s": s, "d": d}})
    return {"case": name, "group": group, "Tbot": Tbot, "z_c": z_c, "s": s}


@pytest.fixture
def ensemble(tmp_path):
    specs = [write_case(tmp_path, "ens_a", -0.5, 250.0, 1.0, m=4.0, drift=0.05),
             write_case(tmp_path, "ens_b", 1.0, 400.0, 0.8, m=9.0),
             write_case(tmp_path, "ens_fail", 0.5, 550.0, 1.0, m=2.0, status="FAILED"),
             {"case": "ens_missing", "group": "factorial", "Tbot": 0.0, "z_c": 400.0,
              "s": 1.0},
             write_case(tmp_path, "ens_T2.0", 2.0, 400.0, 1.0, m=20.0, group="extrapolation")]
    return tmp_path, specs


def test_runs_failing_n3_or_missing_are_excluded_and_listed(ensemble):
    root, specs = ensemble
    ds = ex.build_dataset(root, cases=specs)
    assert list(ds["case"]) == ["ens_a", "ens_b", "ens_T2.0"]
    assert list(ds["excluded_case"]) == ["ens_fail", "ens_missing"]
    assert ds["excluded_reason"][0].startswith("GAP_N3")
    assert ds["excluded_reason"][1] == "GAP_MISSING"


def test_shapes(ensemble):
    root, specs = ensemble
    ds = ex.build_dataset(root, cases=specs)
    n = 3
    assert ds["melt"].shape == (n, NY, NX)
    assert ds["wet_mask"].shape == ds["shelf_mask"].shape == ds["draft"].shape == (n, NY, NX)
    assert ds["z"].shape == (NZ,)
    assert ds["T_profile"].shape == ds["S_profile"].shape == (n, NZ)
    assert ds["x"].shape == (NX,) and ds["y"].shape == (NY,)
    for k in ("Tbot", "z_c", "s", "d", "N1", "N1_pass", "split", "group", "mean_melt_shelf"):
        assert ds[k].shape == (n,), k


def test_melt_is_nan_on_dry_cells_and_exact_elsewhere(ensemble):
    root, specs = ensemble
    ds = ex.build_dataset(root, cases=specs)
    j, i = DRY
    assert np.isnan(ds["melt"][:, j, i]).all()
    assert np.isnan(ds["draft"][:, j, i]).all()
    assert not ds["wet_mask"][:, j, i].any() and not ds["shelf_mask"][:, j, i].any()
    wet = ds["wet_mask"]
    assert np.isfinite(ds["melt"][wet]).all()
    assert np.isnan(ds["melt"][~wet]).all()
    for r, m in enumerate([4.0, 9.0, 20.0]):
        shelf, open_ = ds["shelf_mask"][r], wet[r] & ~ds["shelf_mask"][r]
        np.testing.assert_allclose(ds["melt"][r][shelf], m, rtol=1e-12)
        np.testing.assert_allclose(ds["melt"][r][open_], 0.0, atol=1e-12)
        assert ds["mean_melt_shelf"][r] == pytest.approx(m, rel=1e-12)


def test_split_n1_flags_and_parameters(ensemble):
    root, specs = ensemble
    ds = ex.build_dataset(root, cases=specs)
    assert list(ds["split"]) == ["train", "train", "test"]          # Tbot = 2.0 -> test
    assert ds["N1"][0] == pytest.approx(0.05, rel=1e-9)
    assert list(ds["N1_pass"]) == [False, True, True]               # pd0 N1 < 3 %
    np.testing.assert_allclose(ds["Tbot"], [-0.5, 1.0, 2.0])
    np.testing.assert_allclose(ds["z_c"], [250.0, 400.0, 400.0])
    np.testing.assert_allclose(ds["s"], [1.0, 0.8, 1.0])
    np.testing.assert_allclose(ds["d"], [300.0, 300.0, 300.0])


def test_profiles_are_the_initial_ramp_on_model_levels(ensemble):
    root, specs = ensemble
    ds = ex.build_dataset(root, cases=specs)
    T, S = pr.ramp(ds["z"], -0.5, 250.0, d=300.0)
    np.testing.assert_allclose(ds["T_profile"][0], T, rtol=0, atol=1e-14)
    np.testing.assert_allclose(ds["S_profile"][0], S, rtol=0, atol=1e-14)
    assert (ds["z"] < 0).all()


def test_reference_case_uses_the_ocean0_profile(tmp_path):
    draft, hfac = _geometry()
    make_run(tmp_path / "ocean0_4km_r21", nrec=24, hfac=hfac, draft=draft,
             shi=lambda k: {"SHIfwFlx": _fw(10.0)},
             case_extra={"params": {"case": "ocean0_4km", "init_profile": "ocean0"}})
    ds = ex.build_dataset(tmp_path, cases=[{"case": "ocean0_4km", "group": "reference",
                                             "Tbot": 1.0, "z_c": None, "s": 1.0}])
    T, S = pr.ocean0(ds["z"])
    np.testing.assert_allclose(ds["T_profile"][0], T, atol=1e-14)
    np.testing.assert_allclose(ds["S_profile"][0], S, atol=1e-14)
    assert ds["z_c"][0] == 360.0 and ds["d"][0] == 720.0 and ds["Tbot"][0] == 1.0


def test_write_is_byte_reproducible_and_loads_without_pickle(ensemble, tmp_path):
    root, specs = ensemble
    ds = ex.build_dataset(root, cases=specs)
    a, b = tmp_path / "a.npz", tmp_path / "b.npz"
    sha_a, sha_b = ex.write_npz(ds, a), ex.write_npz(ds, b)
    assert sha_a == sha_b == hashlib.sha256(a.read_bytes()).hexdigest()
    with np.load(a, allow_pickle=False) as z:
        assert set(z.files) == set(ds)
        np.testing.assert_array_equal(z["melt"], ds["melt"])
        assert list(z["case"]) == ["ens_a", "ens_b", "ens_T2.0"]


def test_all_runs_excluded_raises(tmp_path):
    with pytest.raises(ValueError, match="no run"):
        ex.build_dataset(tmp_path, cases=[{"case": "ens_missing", "group": "factorial",
                                           "Tbot": 0.0, "z_c": 400.0, "s": 1.0}])
