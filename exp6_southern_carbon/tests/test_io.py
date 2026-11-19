"""io.grid: a run restarted from a pickup has no grid files (MITgcm writes them only on a cold
start), so the grid is read from the nearest ancestor along pickup_from."""
import json

import numpy as np
import pytest

from exp6_southern_carbon.analysis import io

NR, NY, NX = 3, 4, 2


def _case(d, name, parent):
    (d / name).mkdir()
    (d / name / "case.json").write_text(json.dumps({"nr": NR, "ny": NY, "nx": NX, "pickup_from": parent}))
    return d / name


def _grid_files(run, scale=1.0):
    for f, shape in (("hFacC", (NR, NY, NX)), ("RAC", (NY, NX)), ("DRF", (NR,)), ("RC", (NR,)), ("YC", (NY, NX)), ("DXG", (NY, NX))):
        (scale * np.arange(1, np.prod(shape) + 1)).astype(">f8").tofile(run / f"{f}.data")


def test_grid_of_restarted_run_comes_from_its_cold_start_ancestor(tmp_path):
    g0 = _case(tmp_path, "G0", None)
    _grid_files(g0)
    _case(tmp_path, "G0s", "G0")
    gc = _case(tmp_path, "GC", "G0s")
    g = io.grid(gc)
    assert g["hfac"].shape == (NR, NY, NX)
    np.testing.assert_array_equal(g["drf"], [1, 2, 3])
    assert g["dxg"].shape == (NY, NX)


def test_grid_prefers_the_runs_own_files(tmp_path):
    g0 = _case(tmp_path, "G0", None)
    _grid_files(g0)
    own = _case(tmp_path, "X", "G0")
    _grid_files(own, scale=2.0)
    np.testing.assert_array_equal(io.grid(own)["drf"], [2, 4, 6])


def test_grid_missing_everywhere_raises(tmp_path):
    _case(tmp_path, "G0", None)
    with pytest.raises(FileNotFoundError):
        io.grid(_case(tmp_path, "GC", "G0"))
