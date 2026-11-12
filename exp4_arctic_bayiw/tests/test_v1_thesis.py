"""Tests for exp4_arctic_bayiw.analysis.v1_thesis helpers (synthetic data, no model output)."""
from __future__ import annotations

import numpy as np

from exp4_arctic_bayiw.analysis import v1_thesis as v


def test_psd_hann_peak_and_resolution():
    dt = 900.0
    t = np.arange(480) * dt / 3600.0                       # 120 h
    x = np.sin(2 * np.pi * 0.25 * t)
    f, P, df = v.psd_hann(x, dt)
    assert df == 1.0 / 120.0
    k, fp, _ = v.band_peak(f, P, (0.18, 0.35))
    assert abs(fp - 0.25) <= df / 2


def test_mds_field_reads_one_record(tmp_path):
    a = np.arange(3 * v.NZ * v.NY * v.NX, dtype=float).reshape(3, v.NZ, v.NY, v.NX)
    a.astype(">f8").tofile(tmp_path / "state3d.0000000030.data")
    (tmp_path / "state3d.0000000030.meta").write_text(
        " nDims = [   3 ];\n nrecords = [          3 ];\n nFlds = [    3 ];\n")
    np.testing.assert_array_equal(v.mds_field(tmp_path, "state3d", 30, rec=2), a[2])
    np.testing.assert_array_equal(v.mds_field(tmp_path, "state3d", 30, rec=1, level=0), a[1, 0])
