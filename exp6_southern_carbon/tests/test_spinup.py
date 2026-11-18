from pathlib import Path

import numpy as np
import pytest

from exp6_southern_carbon.analysis import spinup

T10 = Path("/Volumes/LaCie/portfolio_runs/exp6/T10")
DONE = (T10 / "run_status.json").exists() and '"OK"' in (T10 / "run_status.json").read_text()


def test_residual_moc_from_layer_transports():
    nl, ny, nx = 3, 4, 2
    vh = np.zeros((nl, ny, nx))            # v*h per layer (m2 s-1), layer 0 = coldest
    vh[2, 1:, :] = 1.0                     # warm water flows north from row 1 on
    vh[0, 1:, :] = -1.0                    # cold water returns south
    dxg = np.full((ny, nx), 1e6)
    psi = spinup.residual_moc(vh, dxg)     # Sv, psi[l, j] = northward transport warmer than bound l
    assert psi.shape == (nl, ny)
    assert psi[2, 1] == pytest.approx(2.0) and psi[1, 1] == pytest.approx(2.0) and psi[0, 1] == pytest.approx(0.0)
    assert spinup.moc_out_of_R(psi, j_edge=2) == pytest.approx(2.0)


def test_upwelling_fraction_in_R():
    psi_upper = np.array([0.0, 1.0, 3.0, 3.5, 3.6, 3.6])
    assert spinup.upwelling_fraction(psi_upper, j_edge=3) == pytest.approx(3.0 / 3.6)


@pytest.mark.skipif(not DONE, reason="T10 not finished")
def test_report_on_t10():
    r = spinup.report(T10, j_edge=20)
    assert set(r) >= {"v6a", "moc_R_sv", "upwelling_fraction_R", "mld_winter_R_m", "years"}
    assert 0.0 < r["mld_winter_R_m"] < 4000.0 and np.isfinite(r["moc_R_sv"])


def test_run_years_from_first_annual_record():
    spy = 7884
    staged = [spy * k for k in range(1, 201)]          # G0 after two 100-year stages (cold start at 0)
    assert list(spinup.run_years(staged, spy)) == list(range(1, 201))
    child = [5_000_000 + spy * k for k in range(1, 49)]  # a run started from a parent's pickup
    assert list(spinup.run_years(child, spy)) == list(range(1, 49))
