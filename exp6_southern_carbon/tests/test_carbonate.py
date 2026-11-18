from pathlib import Path

import numpy as np
import pytest

from exp6_southern_carbon.analysis import carbonate, io

T10 = Path("/Volumes/LaCie/portfolio_runs/exp6/T10")
DONE = (T10 / "run_status.json").exists() and '"OK"' in (T10 / "run_status.json").read_text()


def test_warming_raises_pco2_about_4_percent_per_degree():
    # Takahashi et al. (1993): d ln pCO2 / dT ~ 0.0423 per degC (Jarnikova et al. 2025 use ~4 %)
    p1, _ = carbonate.pco2(2.10, 2.30, 2.0)
    p2, _ = carbonate.pco2(2.10, 2.30, 3.0)
    assert 0.035 < np.log(p2 / p1) < 0.048


def test_dic_at_pco2_inverts_pco2():
    d = carbonate.dic_at_pco2(400e-6, 2.30, 5.0)
    p, _ = carbonate.pco2(d, 2.30, 5.0)
    assert p == pytest.approx(400e-6, rel=1e-6)


@pytest.mark.skipif(not DONE, reason="T10 not finished")
def test_matches_model_pco2_in_t10():
    it = io.year_iter(T10, 10)
    ptr = io.fields(T10, "pickup_ptracers", it)
    pk = io.fields(T10, "pickup", it)
    dec = io.fields(T10, "surf_mon", it)["DICPCO2"]
    g = io.grid(T10)
    wet = g["hfac"][0] > 0
    sil = io.case(T10)["silica_mol_m3"]
    p, _ = carbonate.pco2(ptr["pTr01"][0][wet], ptr["pTr02"][0][wet], pk["Theta"][0][wet], 34.5,
                          ptr["pTr03"][0][wet], sil)
    assert abs(np.mean(p) - np.mean(dec[wet])) < 5e-6
