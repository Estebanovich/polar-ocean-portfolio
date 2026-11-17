"""gendata.interp: coarse (G0) state -> fine grid, and ensemble-member perturbations (Task 4)."""
import numpy as np
import pytest

from exp5_adjoint_slope.gendata import fields, interp
from exp5_adjoint_slope.gendata.grid import Grid, vertical_spacing

GC = Grid(40, 45, 10e3, vertical_spacing())
GF = Grid(160, 180, 2.5e3, vertical_spacing())
ALPHA, BETA = 3.0e-5, 7.4e-4


def _coarse_state():
    T, S = fields.initial_ts(GC)
    T = T + 0.05 * np.sin(GC.xc / 4e4)[None, None, :]       # some x structure
    wet = GC.zc[:, None, None] < fields.depth(GC.yc)[None, :, None]
    T, S = np.where(wet, T, 0.0), np.where(wet, S, 0.0)   # MITgcm writes 0 in land
    return T, S


def _wet(g):
    return np.broadcast_to(g.zc[:, None, None] < fields.depth(g.yc)[None, :, None],
                           (g.nr, g.ny, g.nx))


def test_no_nan_and_no_land_values_in_wet_fine_cells():
    T, S = interp.coarse_to_fine(*_coarse_state(), GC, GF)
    w = _wet(GF)
    assert T.shape == (35, 180, 160)
    assert np.all(np.isfinite(T[w])) and np.all(np.isfinite(S[w]))
    assert S[w].min() > 30.0                              # no land zeros leaked in


def test_shelf_heat_content_within_one_percent():
    Tc, Sc = _coarse_state()
    Tf, _ = interp.coarse_to_fine(Tc, Sc, GC, GF)

    # Volume-mean heat relative to -1.9 degC over the region both grids cover the same way
    # (y < 70 km: whole coarse cells; the bathymetry is sampled at different centres near the
    # shelf break, so total volumes differ and totals are not comparable).
    def shelf_heat(T, g):
        w = _wet(g) & (g.yc < 70e3)[None, :, None]
        vol = g.dz_m[:, None, None] * g.dx_m ** 2 * w
        return float(((T + 1.9) * vol).sum() / vol.sum())
    assert shelf_heat(Tf, GF) == pytest.approx(shelf_heat(Tc, GC), rel=0.01)


def test_statically_stable():
    Tf, Sf = interp.coarse_to_fine(*_coarse_state(), GC, GF)
    rho = -ALPHA * Tf + BETA * Sf
    w = _wet(GF)
    d = np.diff(rho, axis=0)
    both = w[1:] & w[:-1]
    assert np.all(d[both] >= -1e-12)


def test_perturb_member():
    T = np.zeros((35, 180, 160))
    a = interp.perturb_member(T, k=1, amp=1e-3, seed=7)
    b = interp.perturb_member(T, k=1, amp=1e-3, seed=7)
    c = interp.perturb_member(T, k=2, amp=1e-3, seed=7)
    assert np.array_equal(a, b) and not np.array_equal(a, c)
    assert np.abs(a).max() <= 1e-3 and np.abs(a).max() > 0


def _write_pickup(d, T, S, it=788400):
    # MITgcm pickup layout (G0 stage30 meta): Uvel Vvel Theta Salt GuNm1 GvNm1 (nr records each),
    # then EtaN dEtaHdt EtaH (one record each), big-endian float64.
    from pfcommon.mitgcm_io import write_bin
    nr, ny, nx = T.shape
    z = np.zeros((nr, ny, nx))
    recs = np.concatenate([z, z + 0.01, T, S, z, z, np.zeros((3, ny, nx))])
    d.mkdir(parents=True, exist_ok=True)
    write_bin(d / f"pickup.{it:010d}.data", recs)
    (d / f"pickup.{it:010d}.meta").write_text(
        " nFlds = [    9 ];\n fldList = {\n 'Uvel    ' 'Vvel    ' 'Theta   ' 'Salt    ' 'GuNm1   '"
        " 'GvNm1   ' 'EtaN    ' 'dEtaHdt ' 'EtaH    '\n };\n")


def test_pickup_ts_reads_theta_and_salt(tmp_path):
    T, S = _coarse_state()
    _write_pickup(tmp_path, T, S)
    t, s = interp.pickup_ts(tmp_path / "pickup.0000788400.data", *T.shape)
    assert np.array_equal(t, T) and np.array_equal(s, S)


def test_write_case_init_from_coarse_pickup(tmp_path, monkeypatch):
    # F0 starts from the last G0 pickup interpolated to the fine grid; the sponge target stays
    # the initial profile (rbcs_T), not the G0 state.
    from pfcommon import paths
    from pfcommon.mitgcm_io import read_bin
    from exp5_adjoint_slope.gendata import make_case
    T, S = _coarse_state()
    _write_pickup(tmp_path / "runs" / "G0" / "stage30", T + 0.1, S)
    monkeypatch.setattr(paths, "EXP5_RUNS", tmp_path / "runs")
    out = make_case.write_case("F0", tmp_path / "in")
    tf = read_bin(out / "T.init.bin", (GF.nr, GF.ny, GF.nx))
    exp, _ = interp.coarse_to_fine(T + 0.1, S, GC, GF)
    assert np.allclose(tf, exp)
    t0, _ = fields.initial_ts(GF)
    assert np.allclose(read_bin(out / "rbcs_T.bin", (GF.nr, GF.ny, GF.nx)), t0)
    import json
    assert json.loads((out / "case.json").read_text())["init_pickup"].endswith("pickup.0000788400.data")


def _wet_walls(g):
    # The real bathymetry: land rows 0 and ny-1 (walls), so whole coarse columns are 0.
    b = fields.bathymetry(g)                                  # (ny, nx), negative depth, 0 land
    return g.zc[:, None, None] < -b[None, :, :]


def test_walls_leave_no_land_values_in_wet_fine_cells():
    # F0, 2026-11-16: INI_THETA found 2240 wet points with theta = 0 next to the walls.
    T, S = fields.initial_ts(GC)
    w = _wet_walls(GC)
    Tf, Sf = interp.coarse_to_fine(np.where(w, T, 0.0), np.where(w, S, 0.0), GC, GF)
    wf = _wet_walls(GF)
    assert Sf[wf].min() > 30.0 and np.all(Tf[wf] != 0.0)


def test_perturb_pickup_theta_in_place(tmp_path):
    # FC/FP members (spec Sec. 4): theta of the F0 pickup + noise <= 1e-3 degC in wet cells only;
    # salinity, velocities and land are untouched.
    T, S = _coarse_state()
    _write_pickup(tmp_path, T, S)
    f = tmp_path / "pickup.0000788400.data"
    before = np.fromfile(f, dtype=">f8").copy()
    interp.perturb_pickup_theta(f, *T.shape, k=2, amp=1e-3, seed=7)
    t, s = interp.pickup_ts(f, *T.shape)
    wet = S != 0
    d = t - T
    assert np.abs(d[wet]).max() <= 1e-3 and np.abs(d[wet]).max() > 0
    assert np.all(d[~wet] == 0) and np.array_equal(s, S)
    after = np.fromfile(f, dtype=">f8")
    n3 = T.size
    assert np.array_equal(after[:2 * n3], before[:2 * n3])          # Uvel, Vvel unchanged
    interp.perturb_pickup_theta(f, *T.shape, k=2, amp=0.0, seed=7)  # amp 0: identity
    t2, _ = interp.pickup_ts(f, *T.shape)
    assert np.array_equal(t2, t)
