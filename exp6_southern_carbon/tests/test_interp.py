import json
from pathlib import Path

import numpy as np
import pytest

from exp6_southern_carbon.analysis import fine
from exp6_southern_carbon.gendata import interp
from exp6_southern_carbon.gendata.grid import Grid

GC, GF = Grid(20, 40, 50e3), Grid(100, 200, 10e3)


def _vol(g, wall_rows):
    v = np.ones((g.nr, g.ny, g.nx)) * (g.dx_m ** 2) * g.dr[:, None, None]
    v[:, :wall_rows, :] = 0.0
    ridge = (g.xc[None, :] < 250e3) & ~((g.yc[:, None] > 900e3) & (g.yc[:, None] < 1100e3))
    v[-5:][:, ridge] = 0.0                 # deep ridge with a wet notch, as in the tutorial
    return v


def test_tracers_to_fine_conserves_and_stays_positive():
    vc, vf = _vol(GC, 1), _vol(GF, 5)
    rng = np.random.default_rng(1)
    tr = np.stack([2.0 + 0.1 * rng.random((GC.nr, GC.ny, GC.nx)) for _ in range(6)])
    tr[2] = 1e-4 * rng.random((GC.nr, GC.ny, GC.nx))           # small, positive (PO4-like)
    tr[:, vc == 0] = 0.0                                         # MITgcm writes 0 on land
    out = interp.tracers_to_fine(tr, GC, GF, vc, vf)
    assert out.shape == (6, GF.nr, GF.ny, GF.nx)
    for i in range(6):
        inv_c, inv_f = (tr[i] * vc).sum(), (out[i] * vf).sum()
        assert abs(inv_f - inv_c) / inv_c < 0.005
        assert (out[i][vf > 0] >= 0).all()


def test_h6c():
    r = fine.h6c({1900: [0.0, 0.0, 0.0], 2060: [1.0, 1.1, 0.9]}, ds_coarse=1.5)
    assert r["dS"] == pytest.approx(1.0) and r["ratio"] == pytest.approx(1.0 / 1.5)
    assert r["detectable"] and r["sign_ok"] and r["pass"] and not r["indicative"]
    neg = fine.h6c({1900: [0.0, 0.0, 0.0], 2060: [1.0, 1.1, 0.9]}, ds_coarse=-1.5)
    assert not neg["sign_ok"] and not neg["pass"]
    two = fine.h6c({1900: [0.0, 0.1], 2060: [1.0, 1.1]}, ds_coarse=1.5)
    assert two["indicative"]
    noisy = fine.h6c({1900: [0.0, 2.0, -2.0], 2060: [0.5, 2.5, -1.5]}, ds_coarse=0.5)
    assert not noisy["detectable"] and not noisy["pass"]


def test_v6d():
    rising = list(np.geomspace(1.0, 10.0, 36))               # ~20 % per quarter at the end
    assert not fine.v6d(rising)["pass"]
    flat = [1.0] * 36
    assert fine.v6d(flat)["pass"]
    assert not fine.v6d([1.0] * 24)["pass"]                    # < 3 years minimum


def test_land_column_without_wet_cell_raises():
    vc, vf = _vol(GC, 1), _vol(GF, 5)
    vc[-1, :, 0] = 0.0                      # a whole y-column of land at one level and x
    tr = np.ones((1, GC.nr, GC.ny, GC.nx))
    with pytest.raises(ValueError, match="no wet cell"):
        interp.tracers_to_fine(tr, GC, GF, vc, vf)


def test_wet_volume_partial_cells():
    g, k = Grid(2, 3, 10e3), 30                              # dr ~ 120 m: hFacMin = 0.1 governs
    bathy = np.full((3, 2), -(g.rf[k] + 0.5 * g.dr[k]))      # half of level k wet
    bathy[0, :] = 0.0                                        # wall row: land
    bathy[1, 0] = -(g.rf[k] + 0.04 * g.dr[k])                # 4 % < hFacMin/2 -> dry
    bathy[1, 1] = -(g.rf[k] + 0.07 * g.dr[k])                # 7 %: rounded up to hFacMin
    v = interp.wet_volume(g, bathy, hfac_min=0.1, hfac_min_dr=5.0)
    area = g.dx_m ** 2
    assert (v[:, 0, :] == 0).all()
    assert v[k - 1, 2, 1] == pytest.approx(area * g.dr[k - 1])
    assert v[k, 2, 1] == pytest.approx(area * 0.5 * g.dr[k]) and v[k + 1, 2, 1] == 0.0
    assert v[k, 1, 0] == 0.0 and v[k, 1, 1] == pytest.approx(area * 0.1 * g.dr[k])
    thin = interp.wet_volume(g, np.full((3, 2), -(g.rf[2] + 0.5 * g.dr[2])))
    assert thin[2, 0, 0] == pytest.approx(area * 5.0)        # 7-m level: hFacMinDr = 5 m governs


def test_fine_init_conserves_and_perturbs_theta():
    vc, vf = _vol(GC, 1), _vol(GF, 5)
    rng = np.random.default_rng(2)
    th = np.where(vc > 0, 5.0 - GC.z[:, None, None] / 1000.0 + 0 * vc, 0.0)
    tr = np.stack([np.where(vc > 0, 2.0 + 0.1 * rng.random(vc.shape), 0.0) for _ in range(6)])
    out = interp.fine_init(th, tr, GC, GF, vc, vf, member=0)
    assert set(out) == {"theta", 1, 2, 3, 4, 5, 6}
    base = interp.fine_init(th, tr, GC, GF, vc, vf, member=None)["theta"]
    d = (out["theta"] - base)[vf > 0]
    assert np.abs(d).max() <= 1e-3 and np.abs(d).max() > 0
    assert (out["theta"] == base)[vf == 0].all()
    for i in range(1, 7):
        assert abs((out[i] * vf).sum() / (tr[i - 1] * vc).sum() - 1) < 0.005


def test_stabilize_wet_removes_inversions_conserving_heat():
    g = Grid(1, 2, 10e3)
    vol = np.ones((g.nr, 2, 1)) * g.dr[:, None, None]
    vol[10:, 1, :] = 0.0                                       # column 1: land below level 9
    th = np.repeat(np.linspace(5.0, 0.0, g.nr)[:, None, None], 2, axis=1)
    th[6, 0, 0] = 9.0                                          # warm blob at depth: unstable
    th[12, 1, 0] = 50.0                                        # land value: must not be touched
    out = interp.stabilize_wet(th, vol)
    wetpair = (vol[1:] > 0) & (vol[:-1] > 0)
    assert (np.diff(out, axis=0)[wetpair] <= 1e-12).all()
    assert (out * vol).sum() == pytest.approx((th * vol).sum())
    assert out[12, 1, 0] == 50.0 and (out[:, 1] == th[:, 1]).all()


def _mds(path_stem, arrs, names, itr):
    nr, ny, nx = arrs[0].shape
    np.stack(arrs).astype(">f8").tofile(f"{path_stem}.{itr:010d}.data")
    flds = " ".join(f"'{n:<8}'" for n in names)
    Path(f"{path_stem}.{itr:010d}.meta").write_text(
        f" nDims = [ 3 ];\n dimList = [ {nx}, 1, {nx}, {ny}, 1, {ny}, {nr}, 1, {nr} ];\n"
        f" dataprec = [ 'float64' ];\n nrecords = [ {len(arrs)} ];\n timeStepNumber = [ {itr} ];\n"
        f" nFlds = [ {len(arrs)} ];\n fldList = {{ {flds} }};\n")


def test_eke_series_reads_monthly_records(tmp_path):
    g = Grid(4, 3, 10e3)
    run = tmp_path / "F0_1900"
    (run / "input").mkdir(parents=True)
    (run / "input" / "case.json").write_text(json.dumps({"nr": g.nr, "ny": 3, "nx": 4, "pickup_from": "GT"}))
    hf = np.ones((g.nr, 3, 4)); hf[:, 0] = 0.0
    hf.astype(">f8").tofile(run / "hFacC.data")
    x = np.arange(4.0)[None, None, :] * np.ones((g.nr, 3, 1))
    for m, a in enumerate((1.0, 2.0), start=1):
        u, v = a * x, np.zeros_like(x)
        _mds(run / "eke_mon", [u, v, u ** 2 + 0.5, v ** 2], ["UVEL", "VVEL", "UVELSQ", "VVELSQ"], 100 * m)
    s = fine.eke_series(run)
    # 1/2 <UVELSQ - ubar^2> with ubar the zonal mean: 1/2 (a^2 var(x) + 0.5 + a^2 mean(x)^2 - a^2 mean(x)^2) on wet cells
    assert s == pytest.approx([0.5 * (1.0 * 1.25 + 0.5), 0.5 * (4.0 * 1.25 + 0.5)])


def test_f0_next_years():
    assert fine.f0_next_years(0, False) == 3
    assert fine.f0_next_years(24, False) == 1
    assert fine.f0_next_years(36, True) == 0
    assert fine.f0_next_years(36, False) == 1
    assert fine.f0_next_years(96, False) == 0
