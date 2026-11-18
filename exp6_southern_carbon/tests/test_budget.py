from pathlib import Path

import numpy as np
import pytest

from exp6_southern_carbon.analysis import budget

T10 = Path("/Volumes/LaCie/portfolio_runs/exp6/T10")


def test_conv_moves_mass_between_cells_and_sums_to_surface_flux():
    nr, ny, nx = 2, 3, 4
    fx = np.zeros((nr, ny, nx)); fy = np.zeros((nr, ny, nx)); fr = np.zeros((nr, ny, nx))
    fx[0, 1, 2] = 5.0            # east-going flux through the west face of (0,1,2)
    fy[1, 2, 0] = 3.0            # north-going flux through the south face of (1,2,0)
    fr[0, 0, 0] = 2.0            # upward flux through the sea surface above (0,0,0)
    c = budget._conv(fx, fy, fr)
    assert c[0, 1, 2] == 5.0 and c[0, 1, 1] == -5.0
    assert c[1, 2, 0] == 3.0 and c[1, 1, 0] == -3.0
    assert c[0, 0, 0] == -2.0
    assert c.sum() == pytest.approx(-2.0)      # only the surface face leaves the domain


@pytest.mark.skipif(not (T10 / "run_status.json").exists() or '"OK"' not in (T10 / "run_status.json").read_text(),
                    reason="T10 not finished")
def test_t10_whole_domain_closes():
    b = budget.dic_budget(T10, year=10)
    assert abs(b["residual"]) < 0.01 * abs(b["surf"]), b
