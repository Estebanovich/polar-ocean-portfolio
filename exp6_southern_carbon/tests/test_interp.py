import numpy as np
import pytest

from exp6_southern_carbon.analysis import fine
from exp6_southern_carbon.gendata import interp
from exp6_southern_carbon.gendata.grid import Grid

GC, GF = Grid(20, 40, 50e3), Grid(100, 200, 10e3)


def _vol(g, wall_rows):
    v = np.ones((g.nr, g.ny, g.nx)) * (g.dx_m ** 2) * g.dr[:, None, None]
    v[:, :wall_rows, :] = 0.0
    ridge = (g.xc[None, :] < 250e3) & ~((g.yc[:, None] > 900e3) & (g.yc[:, None] < 1100e3))
    v[-5:][:, ridge] = 0.0                 # deep ridge with a wet notch, as in the tutorial
    return v


def test_tracers_to_fine_conserves_and_stays_positive():
    vc, vf = _vol(GC, 1), _vol(GF, 5)
    rng = np.random.default_rng(1)
    tr = np.stack([2.0 + 0.1 * rng.random((GC.nr, GC.ny, GC.nx)) for _ in range(6)])
    tr[2] = 1e-4 * rng.random((GC.nr, GC.ny, GC.nx))           # small, positive (PO4-like)
    tr[:, vc == 0] = 0.0                                         # MITgcm writes 0 on land
    out = interp.tracers_to_fine(tr, GC, GF, vc, vf)
    assert out.shape == (6, GF.nr, GF.ny, GF.nx)
    for i in range(6):
        inv_c, inv_f = (tr[i] * vc).sum(), (out[i] * vf).sum()
        assert abs(inv_f - inv_c) / inv_c < 0.005
        assert (out[i][vf > 0] >= 0).all()


def test_h6c():
    r = fine.h6c({1900: [0.0, 0.0, 0.0], 2060: [1.0, 1.1, 0.9]}, ds_coarse=1.5)
    assert r["dS"] == pytest.approx(1.0) and r["ratio"] == pytest.approx(1.0 / 1.5)
    assert r["detectable"] and r["sign_ok"] and r["pass"] and not r["indicative"]
    neg = fine.h6c({1900: [0.0, 0.0, 0.0], 2060: [1.0, 1.1, 0.9]}, ds_coarse=-1.5)
    assert not neg["sign_ok"] and not neg["pass"]
    two = fine.h6c({1900: [0.0, 0.1], 2060: [1.0, 1.1]}, ds_coarse=1.5)
    assert two["indicative"]
    noisy = fine.h6c({1900: [0.0, 2.0, -2.0], 2060: [0.5, 2.5, -1.5]}, ds_coarse=0.5)
    assert not noisy["detectable"] and not noisy["pass"]


def test_v6d():
    rising = list(np.geomspace(1.0, 10.0, 36))               # ~20 % per quarter at the end
    assert not fine.v6d(rising)["pass"]
    flat = [1.0] * 36
    assert fine.v6d(flat)["pass"]
    assert not fine.v6d([1.0] * 24)["pass"]                    # < 3 years minimum


def test_land_column_without_wet_cell_raises():
    vc, vf = _vol(GC, 1), _vol(GF, 5)
    vc[-1, :, 0] = 0.0                      # a whole y-column of land at one level and x
    tr = np.ones((1, GC.nr, GC.ny, GC.nx))
    with pytest.raises(ValueError, match="no wet cell"):
        interp.tracers_to_fine(tr, GC, GF, vc, vf)
