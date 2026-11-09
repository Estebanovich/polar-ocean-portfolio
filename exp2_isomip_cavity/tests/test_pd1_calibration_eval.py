"""Tests for exp2_isomip_cavity/analysis/pd1_calibration_eval.py (R15, R17, resolution test)."""
from __future__ import annotations

import numpy as np
import pytest

from exp2_isomip_cavity.analysis import pd1_calibration_eval as ev
from exp2_isomip_cavity.theory import melt as ml


def test_freezing_point_is_ad16_liquidus():
    # AD16 Eq. 25: Tf = lam1 S + lam2 + lam3 p, p = rho0 g |z|
    tf = ev.freezing_point(34.5, -500.0)
    p = 1027.51 * 9.81 * 500.0
    assert tf == pytest.approx(ml.LAMBDA1 * 34.5 + ml.LAMBDA2 + ml.LAMBDA3 * p, rel=1e-14)


def test_cfl_maxima_parses_and_flags():
    log = "\n".join([
        "(PID.TID 0000.0001) %MON advcfl_uvel_max              =   1.0E-01",
        "(PID.TID 0000.0001) %MON advcfl_vvel_max              =   2.0E-01",
        "(PID.TID 0000.0001) %MON advcfl_wvel_max              =   3.0E-01",
        "(PID.TID 0000.0001) %MON advcfl_W_hf_max              =   3.1E-01",
        "(PID.TID 0000.0001) %MON advcfl_uvel_max              =   0.6E+00",
        "(PID.TID 0000.0001) %MON advcfl_vvel_max              =   2.0E-01",
        "(PID.TID 0000.0001) %MON advcfl_wvel_max              =   3.0E-01",
        "(PID.TID 0000.0001) %MON advcfl_W_hf_max              =   3.1E-01"])
    c = ev.cfl_maxima(log)
    assert c["advcfl_uvel_max"] == 0.6 and c["max"] == 0.6 and c["above_0.5"]
    assert c["n_monitor_lines"] == 2
    assert not ev.cfl_maxima(log.replace("0.6E+00", "0.4E+00"))["above_0.5"]
    with pytest.raises(ValueError):
        ev.cfl_maxima("no monitor")


def _box(nz=3, ny=5, nx=6):
    hfac = np.ones((nz, ny, nx))
    hfac[0, :, :2] = 0.0            # ice: column top is level 1 for i < 2
    hfac[1, :, :2] = 0.3            # thin top cells
    rc = -10.0 - 20.0 * np.arange(nz)
    S = np.full((nz, ny, nx), 34.5)
    T = ev.freezing_point(S, rc[:, None, None]) + 1.0   # all warm
    return T, S, hfac, rc


def test_supercooling_counts_and_location():
    T, S, hfac, rc = _box()
    tf = ev.freezing_point(S, rc[:, None, None])
    T[1, 2, 0] = tf[1, 2, 0] - 0.2      # interior row, thin top cell, lateral dry? i=0 edge
    T[2, 0, 3] = tf[2, 0, 3] - 0.05     # wall row (j = 0), interior column, normal cell
    T[2, 2, 3] = tf[2, 2, 3] - 0.01     # interior, no dry neighbour, normal cell
    st = ev.supercooling_stats(T, S, hfac, rc, wall_rows=(0, 4))
    assert st["count"] == 3
    assert st["min_dT"] == pytest.approx(-0.2, abs=1e-12)
    assert st["n_below_minus_0.1"] == 1
    assert st["wall_rows"]["count"] == 1 and st["interior_rows"]["count"] == 2
    assert st["thin_top_cells"]["count"] == 1 and st["other_cells"]["count"] == 2
    assert st["thin_top_cells"]["min_dT"] == pytest.approx(-0.2, abs=1e-12)
    # (1,2,0): domain edge i = 0 counts as a dry neighbour; (2,2,3) has none
    assert st["lateral_dry_neighbour"]["count"] == 1 and st["no_dry_neighbour"]["count"] == 1
    assert st["coldest"]["k"] == 1 and st["coldest"]["hFacC"] == 0.3


def test_supercooling_ignores_dry_cells():
    T, S, hfac, rc = _box()
    T[0, 1, 0] = -50.0                   # dry cell (hfac = 0)
    assert ev.supercooling_stats(T, S, hfac, rc, wall_rows=(0, 4))["count"] == 0


def test_coarsen_2x2_area_weighted_and_drops_land_row():
    f = np.arange(5 * 4, dtype=float).reshape(5, 4)     # 5 rows: 4 wet + 1 land row
    a = np.ones((5, 4))
    a[0, 0] = 3.0
    c = ev.coarsen_2x2(f, a, 2, 2)
    assert c.shape == (2, 2)
    assert c[0, 0] == pytest.approx((0 * 3 + 1 + 4 + 5) / 6.0)
    assert c[1, 1] == pytest.approx((10 + 11 + 14 + 15) / 4.0)
