"""analysis.adjsens: reading adxx files, tau_x sensitivity per record, lead times, scaling."""
import numpy as np
import pytest

from exp5_adjoint_slope.analysis import adjsens
from exp5_adjoint_slope.gendata.grid import Grid, vertical_spacing

G = Grid(40, 45, 10e3, vertical_spacing())


def test_lead_time_centre_to_centre():
    # Record m covers [(m-1)P, mP]; J is the mean over the last 30 days of the window. Lead time
    # = centre of the J window - centre of the record (PD0 amendment: GA = 360 d, records 30 d).
    assert adjsens.lead_time_days(1, 360.0) == pytest.approx(330.0)
    assert adjsens.lead_time_days(12, 360.0) == pytest.approx(0.0)
    assert adjsens.lead_time_days(6, 180.0) == pytest.approx(0.0)


def test_tau_x_sensitivity_sums_over_x():
    ad = np.full((12, G.ny, G.nx), 2.0)
    s = adjsens.tau_x_sensitivity(ad, G)
    assert s.shape == (12, G.ny) and np.allclose(s, 2.0 * G.nx)


def test_scaled():
    assert np.allclose(adjsens.scaled(np.array([1.0, -2.0]), 0.5), [0.5, -1.0])


def test_load_adxx_reads_records(tmp_path):
    from pfcommon.mitgcm_io import write_bin
    a = np.arange(2 * G.ny * G.nx, dtype=float).reshape(2, G.ny, G.nx)
    write_bin(tmp_path / "adxx_fu.0000000000.data", a)
    (tmp_path / "adxx_fu.0000000000.meta").write_text(
        f" nDims = [   2 ];\n dimList = [\n {G.nx}, 1, {G.nx},\n {G.ny}, 1, {G.ny}\n ];\n"
        " dataprec = [ 'float64' ];\n nrecords = [     2 ];\n timeStepNumber = [ 0 ];\n")
    assert np.array_equal(adjsens.load_adxx(tmp_path, "fu"), a)
