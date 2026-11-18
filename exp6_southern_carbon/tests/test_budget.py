from pathlib import Path

import numpy as np
import pytest

from exp6_southern_carbon.analysis import budget

T10 = Path("/Volumes/LaCie/portfolio_runs/exp6/T10")
T2 = Path("/Volumes/LaCie/portfolio_runs/exp6/T2")      # 2 years with the DICGDIC diagnostic


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


@pytest.mark.skipif(not (T2 / "run_status.json").exists() or '"OK"' not in (T2 / "run_status.json").read_text(),
                    reason="T2 not finished")
def test_t2_whole_domain_closes():
    b = budget.dic_budget(T2, year=2)
    assert abs(b["residual"]) < 0.01 * abs(b["surf"]), b


@pytest.mark.skipif(not (T10 / "run_status.json").exists(), reason="T10 not run")
def test_pickup_fields_split_3d_and_2d():
    from exp6_southern_carbon.analysis import io
    it = io.year_iter(T10, 10)
    p = io.fields(T10, "pickup", it)
    assert p["Theta"].shape == (49, 40, 20) and p["EtaN"].shape == (40, 20)
    assert io.fields(T10, "pickup_ptracers", it)["pTr01"].shape == (49, 40, 20)
    assert io.fields(T10, "pickup_dic", it)["DIC_pH2d"].shape == (40, 20)
