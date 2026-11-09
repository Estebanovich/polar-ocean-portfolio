"""Tests for exp2_isomip_cavity/analysis/metrics.py: the ISOMIP+ metrics library (Task 3).

Every test builds a synthetic MITgcm run directory (MDS diagnostics, grid files, case.json,
run_status files) with known fields, so each metric has an analytic answer. Definitions:
notes/pd0_preregistro.md (binding); diagnostics and signs: notes/pd1_config.md.

Review Focus 3 (melt sign/units), 4 (calendar) and 5 (failed/NaN runs refused before any
metric) of the Exp. 2 context are covered by the tests marked "RF3", "RF4", "RF5".
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from exp2_isomip_cavity.analysis import metrics as mt
from exp2_isomip_cavity.theory import melt as ml

MONTH = 2628000.0
NX, NY, NZ = 8, 5, 4
DX = 10e3
DY = 20e3                  # y corners 0-100 km, so the cavity band y > 50 km exists
DZ = 20.0
X0 = 600e3                 # XG of the first column: columns 600-680 km, front at 640 km
SHI_FIELDS = ["SHIfwFlx", "SHIhtFlx", "SHIuStar", "SHIgammT", "SHIgammS", "SHIForcT",
              "SHIForcS", "SHI_TauX", "SHI_TauY", "ETAN"]
DYN_FIELDS = ["UVEL", "VVEL", "WVEL", "THETA", "SALT", "TRAC01", "UTRAC01"]


# ---- synthetic run writer ----------------------------------------------------------------

def _write_mds(stem: Path, fields: dict, names: list, itr: int, t0: float, t1: float):
    arrs = [np.asarray(fields.get(n, 0.0), dtype=float) for n in names]
    shape = np.broadcast_shapes(*[a.shape for a in arrs])
    data = np.stack([np.broadcast_to(a, shape) for a in arrs])
    data.astype(">f8").tofile(f"{stem}.{itr:010d}.data")
    dims = list(reversed(shape))
    dimlist = ",\n".join(f" {n:5d}, 1, {n:5d}" for n in dims)
    flds = " ".join(f"'{n:<8s}'" for n in names)
    Path(f"{stem}.{itr:010d}.meta").write_text(
        f" nDims = [ {len(dims)} ];\n dimList = [\n{dimlist}\n ];\n"
        f" dataprec = [ 'float64' ];\n nrecords = [ {len(names)} ];\n"
        f" timeStepNumber = [ {itr} ];\n timeInterval = [ {t0:.12E} {t1:.12E} ];\n"
        f" nFlds = [ {len(names)} ];\n fldList = {{\n {flds}\n }};\n")


def _write_state(stem: Path, arr, itr: int):
    arr = np.asarray(arr, dtype=float)
    arr.astype(">f8").tofile(f"{stem}.{itr:010d}.data")
    dims = list(reversed(arr.shape))
    dimlist = ",\n".join(f" {n:5d}, 1, {n:5d}" for n in dims)
    Path(f"{stem}.{itr:010d}.meta").write_text(
        f" nDims = [ {len(dims)} ];\n dimList = [\n{dimlist}\n ];\n"
        f" dataprec = [ 'float64' ];\n nrecords = [ 1 ];\n timeStepNumber = [ {itr} ];\n")


def default_geometry():
    """Shelf (draft -40 m) for x < 640 km, open ocean beyond; flat bed at -80 m."""
    xc = X0 + DX * (np.arange(NX) + 0.5)
    draft = np.where(xc < 640e3, -40.0, 0.0)[None, :].repeat(NY, 0)
    bathy = np.full((NY, NX), -80.0)
    rf = -DZ * np.arange(NZ + 1)
    hfac = np.zeros((NZ, NY, NX))
    for k in range(NZ):
        top, bot = rf[k], rf[k + 1]
        hfac[k] = np.clip((np.minimum(top, draft) - np.maximum(bot, bathy)) / DZ, 0, 1)
    return draft, bathy, hfac


def make_run(root: Path, nrec: int = 24, shi=None, dyn=None, *, stages: int = 1,
             hfac=None, draft=None, mask=None, c_start=None, c_end=None, eta_end=None,
             eta_start=None, status="OK", final=True, case_extra=None, t_end_state=None,
             s_end_state=None) -> Path:
    """Write a synthetic run of ``nrec`` monthly records (one time step per month).

    ``shi(k)`` / ``dyn(k)`` return dicts of fields for record k (missing fields are 0).
    With ``stages`` > 1 the records are split evenly over run_dir and run_dir/stage<k>,
    as pfcommon.campaign.run_stages lays them out."""
    root.mkdir(parents=True, exist_ok=True)
    d0, b0, h0 = default_geometry()
    draft = d0 if draft is None else draft
    hfac = h0 if hfac is None else hfac
    xg = X0 + DX * np.arange(NX)
    xc = xg + DX / 2
    yg = DY * np.arange(NY)
    yc = yg + DY / 2
    rf = -DZ * np.arange(NZ + 1)
    rc = 0.5 * (rf[:-1] + rf[1:])
    hW = np.minimum(hfac, np.roll(hfac, 1, axis=2))
    hS = np.minimum(hfac, np.roll(hfac, 1, axis=1))
    grid = {"hFacC": hfac, "hFacW": hW, "hFacS": hS, "RAC": np.full((NY, NX), DX * DY),
            "DRF": np.full(NZ, DZ), "DXG": np.full((NY, NX), DX), "DYG": np.full((NY, NX), DY),
            "XC": np.broadcast_to(xc, (NY, NX)), "YC": np.broadcast_to(yc[:, None], (NY, NX)),
            "XG": np.broadcast_to(xg, (NY, NX)), "YG": np.broadcast_to(yg[:, None], (NY, NX)),
            "RC": rc, "RF": rf}
    for n, a in grid.items():
        np.ascontiguousarray(a, dtype=">f8").tofile(root / f"{n}.data")
    np.ascontiguousarray(draft, dtype=">f8").tofile(root / "icetopo.bin")
    np.ascontiguousarray(b0, dtype=">f8").tofile(root / "bathy.bin")
    if mask is None:
        mask = np.zeros((NZ, NY, NX))
        mask[:, :, -1] = 1.0                     # sponge in the last column
    if mask is not False:
        np.ascontiguousarray(mask, dtype=">f8").tofile(root / "rbcs_mask.bin")
    params = {"nx": NX, "ny": NY, "nz": NZ, "dx": DX, "dy": DY, "dz": DZ, "deltaT": MONTH,
              "nIter0": 0, "nTimeSteps": nrec, "tauRelax": 8640.0, "Cd": 2.5e-3,
              "u_tidal": 0.01, "gammaT": 0.03, "rho_ref": 1027.51, "gravity": 9.81,
              "init_profile": "ocean0", "restore_profile": "ocean0", "case": "synthetic",
              "config": "synthetic"}
    actual = {"grounding_line_x_m": 999e3, "grounding_line_depth_m": 1.0,
              "deepest_wet_ice_base_depth_m": 572.1674970323206, "Tbot": 1.0, "z_c": 360.0,
              "d": 720.0, "s": 1.0}
    if case_extra:
        params.update(case_extra.get("params", {}))
        actual.update(case_extra.get("actual", {}))
    (root / "case.json").write_text(json.dumps({"params": params, "actual": actual}))
    per = nrec // stages
    stage_status = {}
    for s in range(stages):
        sdir = root if s == 0 else root / f"stage{s + 1}"
        sdir.mkdir(exist_ok=True)
        a, n = s * per, per
        for k in range(a, a + n):
            _write_mds(sdir / "shiDiag", (shi or (lambda k: {}))(k), SHI_FIELDS, k + 1,
                       k * MONTH, (k + 1) * MONTH)
            _write_mds(sdir / "dynDiag", (dyn or (lambda k: {}))(k), DYN_FIELDS, k + 1,
                       k * MONTH, (k + 1) * MONTH)
        st = {"state": status if s == stages - 1 else "OK", "last_iter": a + n,
              "n_timesteps": n, "wall_s": 1}
        (sdir / "run_status.json").write_text(json.dumps(st))
        (sdir / "STDOUT.0000").write_text("(PID.TID 0000.0001) %MON time_tsnumber = 1\n")
        stage_status[str(s + 1)] = st
    last = root if stages == 1 else root / f"stage{stages}"
    zero3 = np.zeros((NZ, NY, NX))
    _write_state(root / "PTRACER01", zero3 if c_start is None else c_start, 0)
    _write_state(root / "Eta", np.zeros((NY, NX)) if eta_start is None else eta_start, 0)
    _write_state(last / "PTRACER01", zero3 if c_end is None else c_end, nrec)
    _write_state(last / "Eta", np.zeros((NY, NX)) if eta_end is None else eta_end, nrec)
    _write_state(last / "T", zero3 if t_end_state is None else t_end_state, nrec)
    _write_state(last / "S", zero3 + 34.5 if s_end_state is None else s_end_state, nrec)
    if final:
        (root / "run_status_final.json").write_text(json.dumps(
            {"state": status, "last_iter": nrec, "expected_last_iter": nrec,
             "stages": stage_status}))
    return root


def shelf_mask():
    draft, _, hfac = default_geometry()
    return (hfac > 0).any(0) & (draft < 0)


# ---- RF3: melt sign and units --------------------------------------------------------------

def test_rf3_melt_rate_uniform_flux():
    """1e-3 kg m^-2 s^-1 of melting is 31.536 m/yr of freshwater (rho_fw = 1000)."""
    fw = np.full((3, 4), -1e-3)                 # SHIfwFlx < 0 = melting (pd1_config.md 1.3)
    np.testing.assert_allclose(mt.melt_rate_m_per_a(fw), 31.536, rtol=1e-14)
    np.testing.assert_allclose(mt.melt_rate_m_per_a(-fw), -31.536, rtol=1e-14)   # freezing < 0


def test_rf3_sign_matches_three_equation_law():
    """Warm water under the ice gives SHIfwFlx < 0 in the model law and melt > 0 here."""
    fw = ml.shifwflx(1.0, 34.5, 5.0e6, 0.01, 0.03)
    assert fw < 0
    assert mt.melt_rate_m_per_a(fw) > 0
    np.testing.assert_allclose(mt.melt_rate_m_per_a(fw), ml.melt_m_per_yr(fw), rtol=1e-15)


def test_mean_melt_area_weighted_and_masks(tmp_path):
    draft = default_geometry()[0].copy()
    draft[:, 0] = -400.0                        # one deep column (base < -300 m)
    _, _, hfac = default_geometry()
    fw = np.where(draft < 0, -1e-3, 0.0)
    fw[:, 0] = -2e-3
    rd = mt.load_run(make_run(tmp_path / "r", nrec=24, draft=draft, hfac=hfac,
                              shi=lambda k: {"SHIfwFlx": fw}))
    shelf = rd.shelf
    expect = (31.536 * 3 + 63.072) / 4          # 4 shelf columns per row, equal areas
    assert shelf.sum() == 4 * NY
    assert mt.mean_melt(rd, mt.LAST6) == pytest.approx(expect, rel=1e-13)
    assert mt.mean_melt(rd, mt.LAST6, mask="deep300") == pytest.approx(63.072, rel=1e-13)
    # total melt flux [kg/s]: sum of -SHIfwFlx * rA
    assert mt.total_melt_flux(rd, mt.LAST6) == pytest.approx(
        (3 * 1e-3 + 2e-3) * NY * DX * DY, rel=1e-13)


# ---- RF4: calendar ------------------------------------------------------------------------

def test_rf4_last_six_months_are_records_18_to_23(tmp_path):
    rd = mt.load_run(make_run(tmp_path / "r", nrec=24))
    assert list(mt.calendar_window(rd, "last6")) == [18, 19, 20, 21, 22, 23]
    assert list(mt.calendar_window(rd, "prev6")) == list(range(12, 18))
    assert list(mt.LAST6) == [18, 19, 20, 21, 22, 23]
    assert list(mt.YEAR20) == list(range(228, 240)) and list(mt.YEAR19) == list(range(216, 228))
    assert rd.records[18].t0 == pytest.approx(18 * MONTH)
    assert rd.records[23].t1 == pytest.approx(24 * MONTH)
    assert mt.MONTH_S == 2628000.0


def test_rf4_last6_picks_exactly_those_records(tmp_path):
    """Record k carries melt k+1 m/yr equivalent: last 6 = mean of 19..24."""
    rate = lambda k: -(k + 1) * 1000.0 / ml.SECONDS_PER_YEAR  # noqa: E731 - m = k+1 m/yr
    shelf = shelf_mask()
    rd = mt.load_run(make_run(tmp_path / "r", nrec=24,
                              shi=lambda k: {"SHIfwFlx": np.where(shelf, rate(k), 0.0)}))
    assert mt.mean_melt(rd, mt.LAST6) == pytest.approx(np.mean(np.arange(19, 25)), rel=1e-12)
    assert mt.mean_melt(rd, mt.calendar_window(rd, "prev6")) == pytest.approx(
        np.mean(np.arange(13, 19)), rel=1e-12)


def test_rf4_calendar_window_rejects_wrong_length(tmp_path):
    rd = mt.load_run(make_run(tmp_path / "r", nrec=12))
    with pytest.raises(ValueError, match="24"):
        mt.calendar_window(rd, "last6")
    with pytest.raises(ValueError, match="240"):
        mt.calendar_window(rd, "year20")


def test_records_must_be_contiguous(tmp_path):
    root = make_run(tmp_path / "r", nrec=24)
    for suf in ("data", "meta"):                 # a missing month must not be interpolated
        (root / f"shiDiag.{12:010d}.{suf}").unlink()
        (root / f"dynDiag.{12:010d}.{suf}").unlink()
    with pytest.raises(RuntimeError, match="gap|contiguous|cover"):
        mt.load_run(root)


# ---- N1 ----------------------------------------------------------------------------------

def test_n1_equilibrium_drift_two_year_form(tmp_path):
    shelf = shelf_mask()
    m_of = lambda k: 10.0 if k >= 18 else (9.0 if k >= 12 else 50.0)  # noqa: E731
    rd = mt.load_run(make_run(tmp_path / "r", nrec=24, shi=lambda k: {
        "SHIfwFlx": np.where(shelf, -m_of(k) / ml.SECONDS_PER_YEAR * 1000.0, 0.0)}))
    assert mt.equilibrium_drift(rd) == pytest.approx(0.1, rel=1e-12)


def test_n1_ocean1_form_uses_year19_and_year20(tmp_path):
    shelf = shelf_mask()
    m_of = lambda k: 8.0 if k >= 228 else (7.6 if k >= 216 else 1.0)  # noqa: E731
    rd = mt.load_run(make_run(tmp_path / "r", nrec=240, shi=lambda k: {
        "SHIfwFlx": np.where(shelf, -m_of(k) / ml.SECONDS_PER_YEAR * 1000.0, 0.0)}))
    assert mt.equilibrium_drift(rd) == pytest.approx(0.05, rel=1e-12)
    assert mt.equilibrium_drift(rd, form="ocean1") == pytest.approx(0.05, rel=1e-12)
    with pytest.raises(ValueError):
        mt.equilibrium_drift(rd, form="2yr")


# ---- streamfunctions ---------------------------------------------------------------------

def _all_wet():
    return np.ones((NZ, NY, NX))


def test_barotropic_streamfunction_uniform_channel_flow(tmp_path):
    """Uniform u0 over depth H in a channel: psi(y) = -u0 H y (AD16: U = -d psi/dy)."""
    u0 = 0.05
    rd = mt.load_run(make_run(tmp_path / "r", nrec=2, hfac=_all_wet(),
                              draft=np.zeros((NY, NX)), dyn=lambda k: {"UVEL": u0}))
    psi = mt.barotropic_streamfunction(rd, [0, 1])
    assert psi.shape == (NY + 1, NX)
    y = DY * np.arange(NY + 1)
    H = NZ * DZ
    np.testing.assert_allclose(psi, (-u0 * H * y / 1e6)[:, None] * np.ones(NX), rtol=1e-13,
                               atol=1e-18)


def test_overturning_streamfunction_two_layer(tmp_path):
    """u = +u0 in the top two levels, -u0 below: psi(x, z) = int_z^0 int u dy dz'."""
    u0 = 0.02
    U = np.zeros((NZ, NY, NX))
    U[:2], U[2:] = u0, -u0
    rd = mt.load_run(make_run(tmp_path / "r", nrec=1, hfac=_all_wet(),
                              draft=np.zeros((NY, NX)), dyn=lambda k: {"UVEL": U}))
    psi = mt.overturning_streamfunction(rd, [0])
    assert psi.shape == (NZ + 1, NX)
    W = NY * DY
    expect = np.array([0, 1, 2, 1, 0]) * u0 * DZ * W / 1e6
    np.testing.assert_allclose(psi, expect[:, None] * np.ones(NX), atol=1e-15)
    assert mt.overturning_max(rd, [0]) == pytest.approx(2 * u0 * DZ * W / 1e6, rel=1e-13)


def test_psi_regions_and_circulation_strength(tmp_path):
    """pd0 V2c: cavity y > 50 km, x < 640 km; V2d: 630 < x < 650 km; Y26 Fig. 13: max-min
    of psi for x < 640 km, one month."""
    U = np.zeros((NZ, NY, NX))
    U[:, 2, :] = -0.1                           # westward jet in row 2 (y 40-60 km)
    U[:, 3, 5] = -0.5                           # extra jet at x = 650 km (outside both)
    rd = mt.load_run(make_run(tmp_path / "r", nrec=2, hfac=_all_wet(),
                              draft=np.zeros((NY, NX)), dyn=lambda k: {"UVEL": U * (k + 1)}))
    psi = mt.barotropic_streamfunction(rd, [0])
    H = NZ * DZ
    jet = 0.1 * H * DY / 1e6                    # psi rises by this north of the jet
    # cavity: psi corners with y > 50 km (y = 60, 80, 100 km) and x < 640 km
    yg = DY * np.arange(NY + 1)
    xg = X0 + DX * np.arange(NX)
    sel = (yg[:, None] > 50e3) & (xg[None, :] < 640e3)
    assert mt.psi_bt_max(rd, [0], "cavity") == pytest.approx(np.abs(psi[sel]).max(), rel=1e-13)
    assert mt.psi_bt_max(rd, [0], "cavity") == pytest.approx(jet, rel=1e-13)
    front = (xg > 630e3) & (xg < 650e3)
    assert front.sum() == 1                     # only XG = 640 km
    assert mt.psi_bt_max(rd, [0], "front") == pytest.approx(np.abs(psi[:, front]).max(), rel=1e-13)
    with pytest.raises(ValueError):
        mt.psi_bt_max(rd, [0], "everywhere")
    cs0 = mt.circulation_strength(rd, 0, kind="bt")
    assert cs0 == pytest.approx(jet, rel=1e-12)     # the x = 650 km jet is excluded
    assert mt.circulation_strength(rd, 1, kind="bt") == pytest.approx(2 * jet, rel=1e-12)
    assert mt.circulation_strength(rd, 0, kind="moc") >= 0


def test_r13_psi_bt_max_is_absolute_for_reversed_jet(tmp_path):
    """Ruling R13: an eastward (positive-U) jet gives psi <= 0 everywhere, so the signed
    maximum is 0, but V2c/V2d must read max |psi| = the jet transport."""
    U = np.zeros((NZ, NY, NX))
    U[:, 2, :] = +0.1
    rd = mt.load_run(make_run(tmp_path / "r", nrec=1, hfac=_all_wet(),
                              draft=np.zeros((NY, NX)), dyn=lambda k: {"UVEL": U}))
    psi = mt.barotropic_streamfunction(rd, [0])
    jet = 0.1 * NZ * DZ * DY / 1e6
    assert psi.max() == 0.0 and psi.min() == pytest.approx(-jet, rel=1e-13)
    assert mt.psi_bt_max(rd, [0], "cavity") == pytest.approx(jet, rel=1e-13)
    assert mt.psi_bt_max(rd, [0], "front") == pytest.approx(jet, rel=1e-13)


def test_land_row_adds_nothing_to_psi_or_masks(tmp_path):
    """A dry northern row (the Task 2 walls fix) carries no transport: psi in it equals the
    value at its southern edge, and it holds no shelf or x_front cell."""
    hfac = _all_wet()
    hfac[:, -1, :] = 0.0
    U = np.zeros((NZ, NY, NX))
    U[:, 1, :] = -0.1
    U[:, -1, :] = 5.0                            # garbage velocity on land: masked by hFacW
    draft = default_geometry()[0]
    rd = mt.load_run(make_run(tmp_path / "r", nrec=1, hfac=hfac, draft=draft,
                              dyn=lambda k: {"UVEL": U}))
    psi = mt.barotropic_streamfunction(rd, [0])
    np.testing.assert_array_equal(psi[-1], psi[-2])
    assert not rd.shelf[-1].any()
    assert mt.x_front_index(rd)[-1] == -1


# ---- N2: tracer budget -------------------------------------------------------------------

def _budget_case(tmp_path, negative_c=False):
    """Records with different melt each month and a nonzero sponge sink; the end inventory
    is built to close the budget exactly, with r* from a nonzero final eta."""
    draft, _, hfac = default_geometry()
    shelf = (hfac > 0).any(0) & (draft < 0)
    nrec = 24
    tau = 8640.0
    mask = np.zeros((NZ, NY, NX))
    mask[:, :, -1] = 0.5
    rA = DX * DY
    vol = hfac * DZ * rA
    hcol = (hfac * DZ).sum(0)
    fw_of = lambda k: np.where(shelf, -1e-4 * (1 + 0.1 * k), 0.0)  # noqa: E731
    c_sponge = (lambda k: -1e-6 * (k + 1)) if negative_c else (lambda k: 1e-6 * (k + 1))
    eta_rec = lambda k: np.full((NY, NX), 0.01 * k)                # noqa: E731
    trac_of = lambda k: np.where(mask > 0, c_sponge(k), 0.0)       # noqa: E731
    src = sum((-fw_of(k) * rA).sum() / 1000.0 * MONTH for k in range(nrec))
    sink = 0.0
    for k in range(nrec):
        rs = (eta_rec(k) + hcol) / hcol
        sink += (mask / tau * trac_of(k) * vol * rs[None]).sum() * MONTH
    eta_end = np.full((NY, NX), 0.3)
    rs_end = (eta_end + hcol) / hcol
    # put the remaining meltwater uniformly in the wet cells of the first column
    c_end = np.zeros((NZ, NY, NX))
    tgt = (hfac > 0) & (np.arange(NX)[None, None, :] == 0)
    c_end[tgt] = (src - sink) / (vol * rs_end[None])[tgt].sum()
    root = make_run(tmp_path / "r", nrec=nrec, mask=mask, c_end=c_end, eta_end=eta_end,
                    shi=lambda k: {"SHIfwFlx": fw_of(k), "ETAN": eta_rec(k)},
                    dyn=lambda k: {"TRAC01": trac_of(k)})
    return root, src, sink, fw_of


def test_n2_budget_closes_with_nonzero_sink_over_all_records(tmp_path):
    root, src, sink, fw_of = _budget_case(tmp_path)
    rd = mt.load_run(root)
    b = mt.tracer_budget(rd)
    assert b["rbcs_sink_m3"] == pytest.approx(sink, rel=1e-13)
    assert abs(sink) > 1e-3 * src                       # the sink leg is really exercised
    assert b["source_m3"] == pytest.approx(src, rel=1e-13)
    assert mt.tracer_budget_residual(rd) < 1e-12
    # the Task 2 shortcut (last record x T_run) would not close this budget
    last_x_T = (-fw_of(23) * DX * DY).sum() / 1000.0 * 24 * MONTH
    assert abs(last_x_T - src) / src > 0.1


def test_n2_negative_tracer_is_not_clipped(tmp_path):
    """Refreezing makes C < 0 (volume-fraction budget); it must enter the sink unclipped."""
    root, src, sink, _ = _budget_case(tmp_path, negative_c=True)
    rd = mt.load_run(root)
    b = mt.tracer_budget(rd)
    assert sink < 0
    assert b["rbcs_sink_m3"] == pytest.approx(sink, rel=1e-13)
    assert mt.tracer_budget_residual(rd) < 1e-12


def test_n2_requires_both_legs(tmp_path):
    root = make_run(tmp_path / "r", nrec=2, mask=False)          # no rbcs_mask.bin
    rd = mt.load_run(root)
    with pytest.raises(RuntimeError, match="rbcs"):
        mt.tracer_budget(rd)
    root2 = make_run(tmp_path / "r2", nrec=2)                    # no melt at all: no source leg
    with pytest.raises(RuntimeError, match="source"):
        mt.tracer_budget(mt.load_run(root2))


# ---- H2: outflow depth -------------------------------------------------------------------

def test_h2_outflow_depth_of_delta_transport(tmp_path):
    """A delta of outgoing u C at level 2 gives z_out = RC[2]; inflow (u C < 0) is ignored."""
    UT = np.zeros((NZ, NY, NX))
    i_f = 4                                    # first open-ocean column (shelf max i = 3)
    UT[2, :, i_f] = UT[2, :, i_f + 1] = 1e-3   # outgoing at both faces of the x_front column
    UT[0, :, i_f] = UT[0, :, i_f + 1] = -5e-3  # incoming at the surface: excluded
    rd = mt.load_run(make_run(tmp_path / "r", nrec=24, dyn=lambda k: {"UTRAC01": UT}))
    z = mt.outflow_depth(rd, mt.LAST6)
    assert z == pytest.approx(-50.0, abs=1e-12)                  # RC[2] = -50 m
    assert mt.outflow_depth(rd, mt.LAST6, face="west") == pytest.approx(-50.0, abs=1e-12)
    np.testing.assert_array_equal(mt.x_front_index(rd), np.full(NY, i_f))
    assert mt.outflow_depth(rd, mt.LAST6, x_front_m=X0 + 4.5 * DX) == pytest.approx(-50.0)


def test_h2_weighting_between_two_levels(tmp_path):
    UT = np.zeros((NZ, NY, NX))
    UT[1, :, 4:6], UT[3, :, 4:6] = 1e-3, 3e-3
    rd = mt.load_run(make_run(tmp_path / "r", nrec=24, dyn=lambda k: {"UTRAC01": UT}))
    assert mt.outflow_depth(rd, mt.LAST6) == pytest.approx((-30 * 1 - 70 * 3) / 4, rel=1e-13)


# ---- RF5: failed / diverged runs are refused before any metric ----------------------------

def test_rf5_failed_run_raises(tmp_path):
    root = make_run(tmp_path / "r", nrec=2, status="FAILED")
    with pytest.raises(RuntimeError, match="FAILED"):
        mt.load_run(root)


def test_rf5_nan_in_log_raises(tmp_path):
    root = make_run(tmp_path / "r", nrec=2)
    (root / "STDOUT.0000").write_text("(PID.TID 0000.0001) %MON dynstat_theta_max = NaN\n")
    with pytest.raises(RuntimeError, match="non-finite"):
        mt.load_run(root)


def test_rf5_require_ok_is_called_first(tmp_path, monkeypatch):
    """No file of the run is read before require_ok has accepted it."""
    root = make_run(tmp_path / "r", nrec=2)
    calls = []

    def fake_require_ok(d):
        calls.append(Path(d))
        raise RuntimeError("refused by require_ok")

    monkeypatch.setattr(mt, "require_ok", fake_require_ok)
    (root / "case.json").write_text("not json")          # would fail if read first
    with pytest.raises(RuntimeError, match="refused by require_ok"):
        mt.load_run(root)
    assert calls == [root]


def test_rf5_nan_in_diagnostic_raises(tmp_path):
    shelf = shelf_mask()
    fw = np.where(shelf, -1e-4, 0.0)
    fw[0, 0] = np.nan
    rd = mt.load_run(make_run(tmp_path / "r", nrec=24, shi=lambda k: {"SHIfwFlx": fw}))
    with pytest.raises(RuntimeError, match="non-finite"):
        mt.mean_melt(rd, mt.LAST6)


def test_rf5_final_status_failed_raises(tmp_path):
    root = make_run(tmp_path / "r", nrec=4, stages=2)
    fin = json.loads((root / "run_status_final.json").read_text())
    fin["state"] = "FAILED"
    (root / "run_status_final.json").write_text(json.dumps(fin))
    with pytest.raises(RuntimeError, match="run_status_final"):
        mt.load_run(root)


def test_rf5_failed_later_stage_raises(tmp_path):
    root = make_run(tmp_path / "r", nrec=4, stages=2)
    st = json.loads((root / "stage2" / "run_status.json").read_text())
    st["state"] = "FAILED"
    (root / "stage2" / "run_status.json").write_text(json.dumps(st))
    with pytest.raises(RuntimeError):
        mt.load_run(root)


def test_staged_run_without_final_status_raises(tmp_path):
    root = make_run(tmp_path / "r", nrec=4, stages=2, final=False)
    with pytest.raises(RuntimeError, match="run_status_final"):
        mt.load_run(root)


# ---- staged runs: records concatenated in time order ---------------------------------------

def test_stage_concatenation_in_time_order(tmp_path):
    shelf = shelf_mask()
    rate = lambda k: np.where(shelf, -(k + 1) / ml.SECONDS_PER_YEAR * 1000.0, 0.0)  # noqa: E731
    root = make_run(tmp_path / "r", nrec=24, stages=3, shi=lambda k: {"SHIfwFlx": rate(k)})
    rd = mt.load_run(root)
    assert len(rd.records) == 24
    assert [r.iter for r in rd.records] == list(range(1, 25))
    assert rd.records[0].stage_dir == root and rd.records[8].stage_dir == root / "stage2"
    assert rd.records[23].stage_dir == root / "stage3"
    for k in (0, 7, 8, 15, 16, 23):
        assert mt.mean_melt(rd, [k]) == pytest.approx(k + 1, rel=1e-12)
    assert mt.mean_melt(rd, mt.LAST6) == pytest.approx(21.5, rel=1e-12)


# ---- u* check ----------------------------------------------------------------------------

def test_ustar_check_consistent_fields_pass(tmp_path):
    """Constant velocity: the monthly mean of u* equals u*(mean u) and the check passes."""
    shelf = shelf_mask()
    u = 0.1
    us = np.where(shelf, np.sqrt(2.5e-3 * (u ** 2 + 0.01 ** 2)), 0.0)    # AD16 Eq. 27 (R21)
    rd = mt.load_run(make_run(tmp_path / "r", nrec=2, hfac=_all_wet(),
                              shi=lambda k: {"SHIuStar": us, "SHIgammT": 0.03 * us},
                              dyn=lambda k: {"UVEL": u}))
    res = mt.ustar_check(rd, 1)
    assert res["pass"]
    assert res["floor"]["pass"] and res["gammT_identity"]["pass"]
    assert res["jensen_bound"]["min_ratio"] == pytest.approx(1.0, rel=1e-12)
    assert "not available" in res["equality"]


def test_ustar_check_detects_floor_violation(tmp_path):
    shelf = shelf_mask()
    us = np.where(shelf, 4.0e-4, 0.0)       # below sqrt(Cd) u_tidal = 5.0e-4 m/s (AD16 Eq. 27)
    rd = mt.load_run(make_run(tmp_path / "r", nrec=1,
                              shi=lambda k: {"SHIuStar": us, "SHIgammT": 0.03 * us}))
    res = mt.ustar_check(rd, 0)
    assert not res["floor"]["pass"] and not res["pass"]
    assert res["floor"]["floor_m_s"] == pytest.approx(5.0e-4, rel=1e-12)


def test_ustar_check_floor_is_sqrt_cd_times_u_tidal(tmp_path):
    """R21: the floor of u*^2 = Cd (|u|^2 + u_tidal^2) is sqrt(Cd) u_tidal = 5.0e-4 m/s, so
    u* = 0.005 m/s (below u_tidal, the old wrong floor) passes the floor test."""
    shelf = shelf_mask()
    us = np.where(shelf, 0.005, 0.0)
    rd = mt.load_run(make_run(tmp_path / "r", nrec=1,
                              shi=lambda k: {"SHIuStar": us, "SHIgammT": 0.03 * us}))
    res = mt.ustar_check(rd, 0)
    assert res["floor"]["pass"]


def test_ustar_check_detects_ustar_below_jensen_bound(tmp_path):
    shelf = shelf_mask()
    us = np.where(shelf, 0.0101, 0.0)                   # far below sqrt(Cd u^2 + ut^2)
    rd = mt.load_run(make_run(tmp_path / "r", nrec=1, hfac=_all_wet(),
                              shi=lambda k: {"SHIuStar": us, "SHIgammT": 0.03 * us},
                              dyn=lambda k: {"UVEL": 0.3}))
    res = mt.ustar_check(rd, 0)
    assert not res["jensen_bound"]["pass"] and not res["pass"]


# ---- H1 freezing point at the deepest wet ice base ----------------------------------------

def test_tf_gl_uses_deepest_wet_ice_base_not_wet_edge(tmp_path):
    rd = mt.load_run(make_run(tmp_path / "r", nrec=1))
    D = 572.1674970323206
    p = 1027.51 * 9.81 * D
    expect = -0.0573 * 34.7 + 0.0832 - 7.53e-8 * p          # AD16 Eq. 25, Table 4
    res = mt.freezing_point_gl(rd)
    assert res["T_f_GL_degC"] == pytest.approx(expect, rel=1e-14)
    assert res["depth_m"] == D


# ---- R10 (post hoc): Gade line -----------------------------------------------------------

def test_gade_line_cells_on_the_mixing_line_pass(tmp_path):
    """Cells built on the virtual-flux meltwater mixing line are within tolerance; a cell
    shifted off the line is not; C < 0 is used as is (no clipping)."""
    from exp2_isomip_cavity.gendata import make_case as mc
    rf = -DZ * np.arange(NZ + 1)
    rc = 0.5 * (rf[:-1] + rf[1:])
    Ta, Sa, *_ = mc._profile("ocean0", rc, {})
    Ta3 = np.broadcast_to(Ta[:, None, None], (NZ, NY, NX)).copy()
    Sa3 = np.broadcast_to(Sa[:, None, None], (NZ, NY, NX)).copy()
    C = np.zeros((NZ, NY, NX))
    C[:, 1, :] = 0.03                         # far outside the profile range: on the line
    C[:, 3, :] = -0.01                        # refreezing water: S above the profile
    r = ml.RHO_FW / ml.RHO_SW
    p = ml.RHO_SW * 9.81 * (-rc)
    Tf = ml.LAMBDA1 * Sa + ml.LAMBDA2 + ml.LAMBDA3 * p
    Teff = (Tf - ml.LATENT_HEAT / ml.CW)[:, None, None]
    T = Ta3 - r * C * (Ta3 - Teff)
    S = Sa3 * (1 - r * C)
    T[2, 1, 2] += 0.5                          # one cell off the line
    hfac = _all_wet()
    rd = mt.load_run(make_run(tmp_path / "r", nrec=1, hfac=hfac, draft=np.zeros((NY, NX)),
                              dyn=lambda k: {"THETA": T, "SALT": S, "TRAC01": C}))
    res = mt.gade_line_check(rd, 0)
    n_out = res["n_out_of_range"]
    assert n_out >= NZ * NX                    # the whole C = 0.03 row leaves the range
    assert res["n_within_tol"] == n_out - 1
    assert res["worst"]["T"]["residual"] == pytest.approx(0.5, abs=1e-12)
    assert res["post_hoc"] is True
    # the source-depth variant also accepts the on-line cells and rejects the shifted one
    inf = res["inferred_source_depth"]
    assert inf["n_within_tol"] == n_out - 1
    assert inf["worst_T"]["residual"] == pytest.approx(0.5, abs=1e-9)
    # the end-of-run state can be checked too
    res2 = mt.gade_line_check(mt.load_run(make_run(
        tmp_path / "r2", nrec=1, hfac=hfac, draft=np.zeros((NY, NX)), c_end=C,
        t_end_state=T, s_end_state=S)), None)
    assert res2["n_within_tol"] == res2["n_out_of_range"] - 1


def test_gade_inferred_source_accepts_water_lifted_from_depth(tmp_path):
    """Meltwater mixed into ambient water from 70 m but found at 10 m: the local-ambient
    form rejects it, the source-depth form (S -> z0) accepts it."""
    from exp2_isomip_cavity.gendata import make_case as mc
    rf = -DZ * np.arange(NZ + 1)
    rc = 0.5 * (rf[:-1] + rf[1:])
    Ta, Sa, *_ = mc._profile("ocean0", rc, {})
    r = ml.RHO_FW / ml.RHO_SW
    T = np.broadcast_to(Ta[:, None, None], (NZ, NY, NX)).copy()
    S = np.broadcast_to(Sa[:, None, None], (NZ, NY, NX)).copy()
    C = np.zeros((NZ, NY, NX))
    c = 0.03
    Tf0 = ml.LAMBDA1 * Sa[0] + ml.LAMBDA2 + ml.LAMBDA3 * ml.RHO_SW * 9.81 * (-rc[0])
    Teff0 = Tf0 - ml.LATENT_HEAT / ml.CW
    C[0, 2, 1] = c
    T[0, 2, 1] = (1 - r * c) * Ta[-1] + r * c * Teff0       # ambient from rc[-1] = -70 m
    S[0, 2, 1] = (1 - r * c) * Sa[-1]
    rd = mt.load_run(make_run(tmp_path / "r", nrec=1, hfac=_all_wet(), draft=np.zeros((NY, NX)),
                              dyn=lambda k: {"THETA": T, "SALT": S, "TRAC01": C}))
    res = mt.gade_line_check(rd, 0)
    assert res["n_out_of_range"] == 1
    assert res["model_virtual_flux"]["n_within_tol"] == 0
    assert res["inferred_source_depth"]["n_within_tol"] == 1
    assert abs(res["inferred_source_depth"]["worst_T"]["residual"]) < 1e-3


# ---- summary -----------------------------------------------------------------------------

def test_isomip_summary_writes_metrics_json(tmp_path):
    root, *_ = _budget_case(tmp_path)
    rd = mt.load_run(root)
    out = mt.isomip_summary(rd)
    f = root / "metrics.json"
    assert f.exists()
    js = json.loads(f.read_text())
    for key in ("N1", "N2", "mean_melt_shelf", "mean_melt_deep300", "total_melt_flux",
                "psi_bt_max_cavity", "psi_bt_max_front", "overturning_max",
                "circulation_strength", "z_out", "D_out", "T_f_GL", "window"):
        assert key in js, key
    assert js["window"]["records_0based"] == list(range(18, 24))
    assert js["mean_melt_shelf"]["units"] == "m a-1 (freshwater, positive = melting)"
    assert js["N2"]["value"] == pytest.approx(out["N2"]["value"])
    assert js["N2"]["value"] < 1e-12
