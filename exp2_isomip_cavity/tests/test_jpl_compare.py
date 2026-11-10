"""Tests for analysis/jpl_compare.psi_stats (post hoc circulation extrema, hand values)."""
from __future__ import annotations

import numpy as np
import pytest

from exp2_isomip_cavity.analysis import jpl_compare as jc


def test_psi_stats_regions_and_extrema():
    xb = np.array([600e3, 635e3, 645e3, 700e3])           # x < 640: first two; front: 635, 645
    yb = np.array([0.0, 40e3, 60e3])                      # cavity rows: y > 50 km
    bt = np.array([[0.0, 0.0, 0.0, 9.0],
                   [-0.3, 0.1, 0.0, 9.0],
                   [0.2, -0.5, 0.7, 9.0]])
    moc = np.array([[0.0, 0.1, 5.0, 5.0],
                    [-0.25, 0.05, 5.0, 5.0]])
    s = jc.psi_stats(bt, moc, xb, yb, xb)
    assert s["psi_bt_max_abs_cavity_Sv"] == pytest.approx(0.5)      # row y = 60, x < 640
    assert s["psi_bt_max_abs_front_Sv"] == pytest.approx(0.7)       # 630 < x < 650, all y
    assert s["bt_strength_max_minus_min_Sv"] == pytest.approx(0.2 + 0.5)
    assert s["psi_moc_max_abs_Sv"] == pytest.approx(0.25)
    assert s["moc_strength_max_minus_min_Sv"] == pytest.approx(0.1 + 0.25)


def test_psi_stats_ignores_nan_land():
    xb, yb = np.array([600e3, 610e3]), np.array([60e3, 70e3])
    bt = np.array([[np.nan, 0.4], [0.1, np.nan]])
    moc = np.array([[np.nan, 0.3]])
    s = jc.psi_stats(bt, moc, xb, yb, xb)
    assert s["psi_bt_max_abs_cavity_Sv"] == pytest.approx(0.4)
    assert s["psi_moc_max_abs_Sv"] == pytest.approx(0.3)
    assert s["psi_bt_max_abs_front_Sv"] is None                     # no x in (630, 650) km


def test_psi_stats_full_domain_overturning_max():
    xb = np.array([600e3, 700e3, 766e3])
    yb = np.array([60e3])
    bt = np.zeros((1, 3))
    moc = np.array([[0.185, 0.0, -0.2224]])
    s = jc.psi_stats(bt, moc, xb, yb, xb)
    assert s["psi_moc_max_abs_Sv"] == pytest.approx(0.185)             # x < 640 km
    assert s["psi_moc_max_abs_full_domain_Sv"] == pytest.approx(0.2224)
