"""Tests for the config (c) analysis (Task 7 campaign prep): run parts and
two-stage status (c_runs), Fig. 3 exchange partition, animate_c. Synthetic
MDS output on a tiny grid; MITgcm is never run."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from exp1_fjord_plumes.analysis import c_runs
from exp1_fjord_plumes.analysis.exchange import exchange_flux

NX, NY, NZ, D, DZ, DT = 24, 5, 4, 200.0, 10.0, 60.0
I_SEC = 6
U0 = 0.05


def _write_mds(run: Path, name: str, it: int, recs: np.ndarray) -> None:
    """recs: (nrec, nz, ny, nx) or (nrec, ny, nx)."""
    recs = np.asarray(recs, dtype=float)
    nrec, *shape = recs.shape
    dims = list(reversed(shape))
    stem = run / f"{name}.{it:010d}"
    np.ascontiguousarray(recs, dtype=">f8").tofile(f"{stem}.data")
    dl = ",\n".join(f" {n}, 1, {n}" for n in dims)
    Path(f"{stem}.meta").write_text(
        f" nDims = [ {len(dims)} ];\n dimList = [\n{dl}\n ];\n dataprec = [ 'float64' ];\n"
        f" nrecords = [ {nrec} ];\n timeStepNumber = [ {it} ];\n")


def _params(case, **kw):
    p = {"case": case, "nx": NX, "ny": NY, "nz": NZ, "dx": D, "dy": D, "dz": DZ,
         "deltaT": DT, "f0": 1.37e-4, "H": NZ * DZ, "front_i": 0, "Q0": 100.0,
         "endTime": 7200.0}
    p.update(kw)
    return p


def _make_run(root: Path, case: str, iters, *, two_stage=False, shelf=False,
              t_mon=(-1.5, 4.0), s_mon=(32.9, 34.6), state="OK", u_scale=1.0) -> Path:
    """Run dir with mean3h (UVEL, THETA), animV (VVEL at 1 level), case.json,
    topog (land rows j = 0, NY-1 and column NX-1), T/S init, monitor log.
    u = +U0 in the upper half, -U0 in the lower half (two-layer exchange)."""
    run = root / case
    run.mkdir(parents=True)
    depth = np.full((NY, NX), NZ * DZ)
    depth[0], depth[-1], depth[:, -1] = 0.0, 0.0, 0.0
    np.ascontiguousarray(-depth, ">f8").tofile(run / "topog.bin")
    z = -(np.arange(NZ) + 0.5) * DZ
    T1 = np.linspace(-1.5, 4.0, NZ)
    S1 = np.linspace(32.9, 34.6, NZ)
    for nm, a in (("T.init.bin", T1), ("S.init.bin", S1)):
        np.ascontiguousarray(np.broadcast_to(a[:, None, None], (NZ, NY, NX)), ">f8").tofile(run / nm)
    ob = np.broadcast_to(T1[None, :, None], (3, NZ, NY)).copy()
    obs = np.broadcast_to(S1[None, :, None], (3, NZ, NY)).copy()
    if shelf:
        ob[1] -= 0.01  # heaved profile slightly outside the initial range
    np.ascontiguousarray(ob, ">f8").tofile(run / "OBEt.bin")
    np.ascontiguousarray(obs, ">f8").tofile(run / "OBEs.bin")
    meta = {"params": _params(case, externForcingCycle=3 * 3600.0, externForcingPeriod=3600.0),
            "T_out": 3.0, "S_out": 33.5}
    (run / "case.json").write_text(json.dumps(meta))
    u = np.where(np.arange(NZ)[:, None, None] < NZ // 2, U0, -U0) * u_scale
    u = np.broadcast_to(u, (NZ, NY, NX)).copy()
    T = np.broadcast_to(T1[:, None, None], (NZ, NY, NX)).copy()
    parts = [run] if not two_stage else [run, run / "stage2"]
    if two_stage:
        (run / "stage2").mkdir()
        (run / "data.stage2").write_text("&PARM03\n nIter0=%d,\n&\n" % iters[len(iters) // 2 - 1])
    for n, it in enumerate(iters):
        d = parts[0] if (not two_stage or n < len(iters) // 2) else parts[1]
        _write_mds(d, "mean3h", it, np.stack([u, T]))
        _write_mds(d, "animV", it, np.zeros((1, 1, NY, NX)))
    for d in parts:
        (d / "output.txt").write_text(
            f"(PID.TID 0000.0001) %MON dynstat_theta_max = {t_mon[1]:.13E}\n"
            f"(PID.TID 0000.0001) %MON dynstat_theta_min = {t_mon[0]:.13E}\n"
            f"(PID.TID 0000.0001) %MON dynstat_salt_max = {s_mon[1]:.13E}\n"
            f"(PID.TID 0000.0001) %MON dynstat_salt_min = {s_mon[0]:.13E}\n")
        (d / "run_status.json").write_text(json.dumps(
            {"state": state, "last_iter": iters[-1], "n_timesteps": iters[-1], "wall_s": 1}))
    if two_stage:
        (run / "run_status_final.json").write_text(json.dumps(
            {"state": state, "last_iter": iters[-1], "expected_last_iter": iters[-1]}))
    return run


ITERS = [180, 360, 540, 720]   # 3-h means with dt = 60 s


def _A_wet():
    return (NY - 2) * NZ * D * DZ


# ---- c_runs --------------------------------------------------------------------

def test_case_parts_two_stage_and_status(tmp_path):
    run = _make_run(tmp_path, "c_fjord_wind", ITERS, two_stage=True)
    assert c_runs.case_parts(run) == [run, run / "stage2"]
    fin = json.loads((run / "run_status_final.json").read_text())
    fin["state"] = "FAILED"
    (run / "run_status_final.json").write_text(json.dumps(fin))
    with pytest.raises(RuntimeError, match="FAILED"):
        c_runs.case_parts(run)


def test_case_parts_refuses_failed_single_stage(tmp_path):
    run = _make_run(tmp_path, "c_fjord_plume", ITERS, state="FAILED")
    with pytest.raises(RuntimeError):
        c_runs.case_parts(run)


def test_section_inside_sponge_refused(tmp_path):
    run = _make_run(tmp_path, "c_fjord_plume", ITERS)
    with pytest.raises(ValueError, match="sponge"):
        c_runs.section_area(run, NX - 2 - 10)
    dA = c_runs.section_area(run, I_SEC)
    assert dA.shape == (NZ, NY) and np.all(dA[:, 0] == 0) and dA.sum() == _A_wet()


def test_exchange_series_two_layer_and_two_stages(tmp_path):
    run = _make_run(tmp_path, "c_fjord_all", ITERS, two_stage=True)
    ts = c_runs.exchange_series(run, I_SEC, T_ref=-1.5)
    assert np.allclose(ts["t_s"], np.array(ITERS) * DT - 5400.0)
    assert np.allclose(ts["Q_ex"], 0.5 * U0 * _A_wet())
    # heat toward the front: rho0 cp sum(-u)(T - Tref) dA (upper layer warm? no: T rises with depth)
    T1 = np.linspace(-1.5, 4.0, NZ)
    u = np.where(np.arange(NZ) < NZ // 2, U0, -U0)
    H = 1027.0 * 3974.0 * np.sum(-u * (T1 + 1.5)) * (NY - 2) * D * DZ
    assert np.allclose(ts["heat_W"], H)


def test_c5p_range_includes_outflow_and_ob_series(tmp_path):
    ok = _make_run(tmp_path / "a", "c_fjord_plume", ITERS)
    assert c_runs.c5p_check(ok)["pass"]
    bad = _make_run(tmp_path / "b", "c_fjord_plume", ITERS, t_mon=(-1.6, 4.0))
    r = c_runs.c5p_check(bad)
    assert not r["pass"] and r["allowed"]["T"][0] == pytest.approx(-1.5)
    sh = _make_run(tmp_path / "c", "c_fjord_shelf", ITERS, shelf=True, t_mon=(-1.505, 4.0))
    r = c_runs.c5p_check(sh)
    assert r["pass"] and r["allowed"]["T"][0] == pytest.approx(-1.51)


def test_w_over_rd(tmp_path):
    run = _make_run(tmp_path, "c_fjord_plume", ITERS)
    r = c_runs.w_over_rd(run, h1=20.0)
    from exp1_fjord_plumes.theory.plume import LinearEOS
    e = LinearEOS()
    gp = 9.81 * (e.rho(4.0, 34.6) - e.rho(-1.5, 32.9)) / 1027.0
    c = np.sqrt(gp * 20.0 * 20.0 / 40.0)
    assert r["W_m"] == (NY - 2) * D
    assert r["Rd_m"] == pytest.approx(c / 1.37e-4, rel=1e-6)
    assert r["W_over_Rd"] == pytest.approx((NY - 2) * D * 1.37e-4 / c, rel=1e-6)


# ---- Fig. 3 --------------------------------------------------------------------

def test_fig3_writes_metrics_and_figure(tmp_path):
    from exp1_fjord_plumes.analysis import fig3_exchange_partition as f3
    root = tmp_path / "runs"
    for case, sc, two in (("c_fjord_plume", 1.0, False), ("c_fjord_wind", 2.0, True),
                          ("c_fjord_shelf", 3.0, False), ("c_fjord_all", 4.0, True)):
        _make_run(root, case + "_x", ITERS, two_stage=two, u_scale=sc)
    out = tmp_path / "fig"
    m = f3.main(root, out, suffix="_x", i_section=I_SEC)
    met = json.loads((out / "fig3_metrics.json").read_text())
    assert (out / "fig3_exchange_partition.png").exists()
    assert (out / "fig3_exchange_partition.pdf").exists()
    for case, sc in (("c_fjord_plume", 1), ("c_fjord_wind", 2), ("c_fjord_shelf", 3), ("c_fjord_all", 4)):
        assert met["cases"][case]["Q_ex_mean_m3_s"] == pytest.approx(0.5 * sc * U0 * _A_wet())
        assert met["cases"][case]["C5p"]["pass"]
    assert met["i_section"] == I_SEC and "W_over_Rd" in met and met == m


def test_fig3_skips_missing_cases(tmp_path):
    from exp1_fjord_plumes.analysis import fig3_exchange_partition as f3
    root = tmp_path / "runs"
    _make_run(root, "c_fjord_plume", ITERS)
    m = f3.main(root, tmp_path / "fig", i_section=I_SEC)
    assert list(m["cases"]) == ["c_fjord_plume"] and m["missing"] == [
        "c_fjord_wind", "c_fjord_shelf", "c_fjord_all"]


# ---- animate_c -----------------------------------------------------------------

def test_animate_c_frames_across_stages_and_gif_budget(tmp_path):
    from exp1_fjord_plumes.analysis import animate_c as ac
    run = _make_run(tmp_path, "c_fjord_all", ITERS, two_stage=True)
    fr = ac.frame_index(run)
    assert [it for it, _ in fr] == ITERS
    assert fr[0][1] == run and fr[-1][1] == run / "stage2"
    res = ac.make_gif(run, tmp_path / "a.gif", dpi=40, max_mb=8.0, i_section=I_SEC, k_plan=1)
    assert res["mb"] <= 8.0 and res["frames"] == len(ITERS) and (tmp_path / "a.gif").exists()


def test_animate_c_check_scales(tmp_path):
    from exp1_fjord_plumes.analysis import animate_c as ac
    run = _make_run(tmp_path, "c_fjord_plume", ITERS)
    r = ac.check_scales(run, i_section=I_SEC, k_plan=1)
    assert set(r) >= {"a_T_plan", "b_T_section", "c_Q_ex", "limits"}
    assert r["a_T_plan"]["max"] == 0.0           # T within the fixed scale
    assert r["c_Q_ex"]["max_over_lim"] < 1.0


def test_c_diagnostics_template_has_anim_v():
    t = (Path(__file__).resolve().parents[1] / "config" / "c_fjord3d" / "input_templates"
         / "data.diagnostics").read_text()
    assert "'animV'" in t and "levels(1,4) = 19." in t
