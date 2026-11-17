"""gendata.make_case.write_case (Exp. 5 Task 1)."""
import json

import numpy as np
import pytest

from exp5_adjoint_slope.gendata import make_case
from pfcommon.mitgcm_io import read_bin


def test_write_case_g0(tmp_path):
    out = make_case.write_case("G0", tmp_path / "G0")
    meta = json.loads((out / "case.json").read_text())
    nx, ny, nr = meta["nx"], meta["ny"], meta["nr"]
    # J line at y = 170 km, ~920 m isobath (notes/pd1b_config_dettling.md).
    assert (nx, ny, nr) == (40, 45, 35) and meta["j_qb"] == 17
    assert 900 < meta["y_qb_depth_m"] < 940
    sizes = {"bathy.bin": ny * nx, "T.init.bin": nr * ny * nx, "S.init.bin": nr * ny * nx,
             "tau_x.bin": ny * nx, "sst_relax.bin": ny * nx, "sss_relax.bin": ny * nx,
             "lambda_t.bin": ny * nx, "lambda_s.bin": ny * nx, "salt_flux.bin": ny * nx,
             "rbcs_mask.bin": nr * ny * nx, "rbcs_mask_uv.bin": nr * ny * nx,
             "rbcs_T.bin": nr * ny * nx, "rbcs_S.bin": nr * ny * nx,
             "gm_map.bin": ny * nx, "delZ.bin": nr}
    for name, n in sizes.items():
        assert (out / name).stat().st_size == n * 8, name
    b = read_bin(out / "bathy.bin", (ny, nx))
    assert b.min() < -2900 and b.max() > -600
    raw = np.fromfile(out / "delZ.bin", dtype=">f8")
    assert raw.sum() == np.float64(3000.0) or abs(raw.sum() - 3000) < 1e-9


def test_write_case_fine_j_qb(tmp_path):
    meta = json.loads((make_case.write_case("F0", tmp_path / "F0") / "case.json").read_text())
    assert (meta["nx"], meta["ny"], meta["j_qb"]) == (160, 180, 68)


def test_registry_numbers_are_numbers():
    # PyYAML (YAML 1.1) reads 90.0e3 as a string; every numeric key must load as int/float.
    for case in ("G0", "F0", "GA"):
        p = make_case.load_case(case)
        for k, v in p.items():
            if k not in ("config", "case", "input_dirs", "pickup_from", "init_from", "kind", "grdchk", "tau_pert"):
                assert isinstance(v, (int, float)), (case, k, v)


def test_tap_case_merges_inputs_and_writes_weights(tmp_path):
    out = make_case.write_case("GS", tmp_path / "GS")
    assert (out / "data.gmredi").exists() and (out / "data.ctrl").exists()
    assert "useGrdchk" in (out / "data.pkg").read_text()        # chan10_tap overrides chan10
    assert (out / "ones_64b.bin").stat().st_size == 35 * 45 * 40 * 8


def test_tau_pert_adds_interpolated_coarse_profile(tmp_path):
    # FP members (H5b): the coarse dtau_x(y) profile (JSON y_m / dtau) interpolated to the fine
    # rows and added to the base wind, uniform in x and constant in time.
    import numpy as np
    from pfcommon.mitgcm_io import read_bin
    f = tmp_path / "pat.json"
    f.write_text(json.dumps({"y_m": [5e3, 445e3], "dtau": [0.0, 0.01]}))
    base = make_case.write_case("FC1", tmp_path / "C")
    pert = make_case.write_case("FC1", tmp_path / "P", {"tau_pert": str(f)})
    tc = read_bin(base / "tau_x.bin", (180, 160))
    tp = read_bin(pert / "tau_x.bin", (180, 160))
    yc = (np.arange(180) + 0.5) * 2.5e3
    expect = np.interp(yc, [5e3, 445e3], [0.0, 0.01])
    assert np.allclose(tp - tc, expect[:, None])
    assert json.loads((pert / "case.json").read_text())["tau_pert"] == str(f)
