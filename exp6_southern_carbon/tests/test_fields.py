import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from exp6_southern_carbon.gendata import fields
from exp6_southern_carbon.gendata.grid import Grid

TUT = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification/tutorial_reentrant_channel/input")
G50 = Grid(20, 40, 50e3)
G10 = Grid(100, 200, 10e3)
FILES = {"bathy": "bathy.50km.bin", "taux": "zonal_wind.50km.bin", "sst": "SST_relax.50km.bin",
         "mask": "T_relax_mask.50km.bin", "temp": "temperature.50km.bin"}


def test_grid_is_tutorial_vertical_grid():
    assert G50.dr.size == 49 and G50.dr.sum() == pytest.approx(3982.26086, abs=1e-4)
    assert G50.yc[0] == 25e3 and G10.yc[-1] == 1995e3


@pytest.mark.parametrize("key", sorted(FILES))
def test_g50_matches_tutorial_bitwise(key):
    ours = fields.tutorial_fields(G50)[key].astype(">f4").tobytes()
    assert ours == (TUT / FILES[key]).read_bytes()


def test_fine5_port_matches_tutorial_script(tmp_path):
    shutil.copy(TUT / "gendata.5km.py", tmp_path)
    subprocess.run([sys.executable, "-I", "gendata.5km.py"], cwd=tmp_path, check=True,
                   env={"MPLBACKEND": "Agg", "PATH": "/usr/bin:/bin"})
    ours = fields.fine5_fields()
    for key, name in FILES.items():
        f = tmp_path / name.replace("50km", "5km")
        assert ours[key].astype(">f4").tobytes() == f.read_bytes(), key


def test_f10_is_2x2_mean_of_5km():
    f5, f10 = fields.fine5_fields(), fields.tutorial_fields(G10)
    b = f5["bathy"].reshape(200, 2, 100, 2).mean(axis=(1, 3))
    np.testing.assert_array_equal(f10["bathy"], b)
    assert f10["bathy"].shape == (200, 100) and (f10["bathy"][0] == 0).all()
    assert f10["mask"].shape == (49, 200, 100)


def test_wind_factor():
    np.testing.assert_array_equal(fields.taux(G50, 1.2), 1.2 * fields.taux(G50))


def test_seasonal_sst():
    m = fields.sst_relax_monthly(G50, amp=1.5, month_min=8)
    assert m.shape == (12, 40, 20)
    np.testing.assert_allclose(m.mean(0), fields.tutorial_fields(G50)["sst"], atol=1e-12)
    assert int(np.argmin(m[:, 10, 0])) == 7 and int(np.argmax(m[:, 10, 0])) == 1
    assert (m.max(0) - m.min(0)).max() == pytest.approx(3.0)
