"""exp5_adjoint_slope/scripts/run_adjoint.py helpers (plan Task 3)."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("run_adjoint", PKG / "scripts" / "run_adjoint.py")
ra = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ra)


def test_subruns_expand_grdchk_points_and_eps():
    p = {"grdchk": [{"var": "xx_fu", "i": 20, "j": 10, "k": 1, "rec": 1, "eps": [1e-3, 1e-4]},
                    {"var": "xx_theta", "i": 20, "j": 12, "k": 5, "rec": 1, "eps": [1e-2]}]}
    subs = ra.subruns(p)
    assert [s["name"] for s in subs] == ["xx_fu_e0", "xx_fu_e1", "xx_theta_e0"]
    assert subs[1]["eps"] == 1e-4 and subs[2]["k"] == 5


def test_subruns_without_grdchk_is_single_run():
    assert ra.subruns({}) == [{"name": None}]


def test_last_fc():
    log = ("(PID.TID 0000.0001)  global fc =  0.100000000000000D+01\n"
           "(PID.TID 0000.0001)  global fc = -0.605596159366832D+12\n")
    assert ra.last_fc(log) == pytest.approx(-0.605596159366832e12)
    assert ra.last_fc("nothing") is None


def test_namelists_for_gs_like_case(tmp_path):
    import shutil
    from exp5_adjoint_slope.gendata.make_case import write_case
    inp = write_case("GS", tmp_path / "in")
    p = {"window_days": 10, "last_days": 5, "ctrl_period_days": 0, "pickup_from": "T2Y"}
    sub = {"name": "xx_fu_e0", "var": "xx_fu", "i": 20, "j": 10, "k": 1, "rec": 1, "eps": 1e-3}
    ra.edit_namelists(inp, p, sub)
    d = (inp / "data").read_text()
    assert " nTimeSteps=720," in d and " pickupSuff='G0endstate'," in d
    # periodic forcing always (test_forcing_reloaded_every_step); the old assertion that
    # externForcingCycle is absent encoded the accumulating-wind-control bug.
    assert " externForcingCycle=864000.0," in d and " periodicExternalForcing=.TRUE.," in d
    assert " lastinterval=432000.0," in (inp / "data.cost").read_text()
    g = (inp / "data.grdchk").read_text()
    assert " grdchkvarname='xx_fu'," in g and " grdchk_eps=0.001," in g and " jGloPos=10," in g
    # nbeg=0 makes grdchk locate the point from iGloPos/jGloPos/kGloPos (grdchk_main.F:227);
    # nbeg>=1 would walk the wet-point index and ignore the position.
    assert " nbeg=0," in g and " nend=0," in g
    assert " useGrdchk=.TRUE.," in (inp / "data.pkg").read_text()


def test_monthly_controls_set_cycle(tmp_path):
    from exp5_adjoint_slope.gendata.make_case import write_case
    inp = write_case("GS", tmp_path / "in")
    ra.edit_namelists(inp, {"window_days": 360, "last_days": 30, "ctrl_period_days": 30}, {"name": None})
    d = (inp / "data").read_text()
    assert " externForcingCycle=31104000.0," in d and " nTimeSteps=25920," in d
    c = (inp / "data.ctrl").read_text()
    assert " xx_gentim2d_period(1)=2592000.0," in c and " xx_gentim2d_period(2)=2592000.0," in c
    assert " useGrdchk=.FALSE.," in (inp / "data.pkg").read_text()


def test_window_must_be_multiple_of_period(tmp_path):
    from exp5_adjoint_slope.gendata.make_case import write_case
    inp = write_case("GS", tmp_path / "in")
    with pytest.raises(ValueError):
        ra.edit_namelists(inp, {"window_days": 365, "last_days": 30, "ctrl_period_days": 30}, {"name": None})


def test_snapshots_for_skew_anchor(tmp_path):
    # J = resolved + bolus + Redi (option A): the anchor rebuilds bolus and Redi fluxes from
    # GM_PsiY, GM_Kvy, GM_Kvz; GM_vbT and DFyE_TH are the model's own fluxes, a second check.
    from exp5_adjoint_slope.gendata.make_case import write_case
    inp = write_case("GJ", tmp_path / "in")
    p = {"window_days": 2, "last_days": 1, "ctrl_period_days": 0, "pickup_from": "G0",
         "snapshots": True}
    ra.edit_namelists(inp, p, {})
    d = (inp / "data.diagnostics").read_text()
    assert ("'VVEL    ','THETA   ','GM_PsiY ','GM_Kvy  ','GM_Kvz  ','GM_vbT  ','DFyE_TH '"
            in d)


def test_copy_inputs_skips_appledouble(tmp_path):
    # exFAT (LaCie) keeps macOS metadata in ._<name> files that vanish when <name> is rewritten:
    # GJ died on FileNotFoundError '._eedata' (2026-11-16).
    inp, rd = tmp_path / "input", tmp_path / "run"
    inp.mkdir(); rd.mkdir()
    (inp / "data").write_text("x")
    (inp / "._data").write_text("meta")
    ra.copy_inputs(inp, rd)
    assert sorted(f.name for f in rd.iterdir()) == ["data"]


def test_remove_tree_tolerates_vanished_appledouble(tmp_path, monkeypatch):
    # On exFAT, deleting <name> also deletes ._<name>; rmtree then fails on the listed ._<name>
    # (GJ rerun, FileNotFoundError '._data.gmredi'). Only a vanished file is tolerated.
    rd = tmp_path / "run"
    rd.mkdir()
    (rd / "data").write_text("x")
    (rd / "._data").write_text("meta")
    real_unlink = ra.os.unlink

    def unlink(path, *a, **k):                    # ._<name> is listed but already gone
        real_unlink(path, *a, **k)
        if ra.os.path.basename(ra.os.fspath(path)).startswith("._"):
            raise FileNotFoundError(path)

    monkeypatch.setattr(ra.os, "unlink", unlink)
    ra.remove_tree(rd)
    assert not rd.exists()


def test_subruns_repeated_control_get_distinct_names():
    # GSU checks xx_fu at three rows; names must not collide (run directories would be reused).
    p = {"grdchk": [{"var": "xx_fu", "i": 20, "j": 18, "k": 1, "rec": 1, "eps": [1e-5]},
                    {"var": "xx_fu", "i": 20, "j": 10, "k": 1, "rec": 1, "eps": [1e-3, 1e-4]}]}
    names = [s["name"] for s in ra.subruns(p)]
    assert names == ["xx_fu_i20j18k1_e0", "xx_fu_i20j10k1_e0", "xx_fu_i20j10k1_e1"]


def test_xx_from_installs_controls_and_keeps_them(tmp_path, monkeypatch):
    # V2/V3 forward runs with perturbed controls: pkg/ctrl overwrites xx_* with zeros unless
    # doInitXX = .FALSE. (pkg/ctrl/ctrl_readparms.F:171, ctrl_init_ctrlvar.F:140).
    from exp5_adjoint_slope.gendata.make_case import write_case
    inp = write_case("GS", tmp_path / "in")
    src = tmp_path / "pert"
    src.mkdir()
    (src / "xx_fu.0000000000.data").write_bytes(b"\0" * 8)
    (src / "xx_fu.0000000000.meta").write_text("m")
    p = {"window_days": 10, "last_days": 5, "ctrl_period_days": 0, "xx_from": str(src)}
    ra.edit_namelists(inp, p, {"name": None})
    ra.install_xx(inp, p)
    assert " doInitXX=.FALSE.," in (inp / "data.ctrl").read_text()
    assert (inp / "xx_fu.0000000000.data").exists() and (inp / "xx_fu.0000000000.meta").exists()


def test_without_xx_from_controls_start_at_zero(tmp_path):
    from exp5_adjoint_slope.gendata.make_case import write_case
    inp = write_case("GS", tmp_path / "in")
    p = {"window_days": 10, "last_days": 5, "ctrl_period_days": 0}
    ra.edit_namelists(inp, p, {"name": None})
    ra.install_xx(inp, p)
    assert "doInitXX" not in (inp / "data.ctrl").read_text()
    assert not list(inp.glob("xx_*"))


def test_nml_overrides_apply_after_case_edits(tmp_path):
    # Adjoint-instability diagnostics (2026-11-16): per-case namelist overrides,
    # {file: {group: {key: value}}}, applied last so they win over the window edits.
    from exp5_adjoint_slope.gendata.make_case import write_case
    inp = write_case("GS", tmp_path / "in")
    p = {"window_days": 10, "last_days": 5, "ctrl_period_days": 0,
         "nml": {"data": {"PARM01": {"viscC4leith": "0."}},
                 "data.gmredi": {"GM_PARM01": {"GM_taper_scheme": "'clipping'"}}}}
    ra.edit_namelists(inp, p, {"name": None})
    assert " viscC4leith=0.," in (inp / "data").read_text()
    assert " GM_taper_scheme='clipping'," in (inp / "data.gmredi").read_text()


def test_jobs_expand_case_globs_with_subruns():
    # One pool for many cases (forward ensembles): every case matching the glob, and every
    # grdchk subrun of each, with unique labels <case>/<subrun>.
    jobs = ra.jobs("GV2d_*")
    assert [j[0] for j in jobs] == [f"GV2d_{n}_{s}" for n in ("fu", "theta", "diffkr") for s in "pm"]
    assert all(j[2]["name"] is None for j in jobs)
    labels = [ra.label(j) for j in ra.jobs("GV1")]
    assert len(labels) == 14 and len(set(labels)) == 14 and labels[0].startswith("GV1/xx_fu_")
    import pytest
    with pytest.raises(KeyError):
        ra.jobs("NOPE*")


@pytest.mark.parametrize("per", [0, 30])
def test_forcing_reloaded_every_step(tmp_path, per):
    # pkg/ctrl adds xx_fu to fu every step (ctrl_map_forcing.F:119, forward_step.F:572), and fu
    # is reset from the file only with periodicExternalForcing (external_fields_load.F:268).
    # Without it the wind control accumulates (V3E blew up at day 7 with 0.01 N m-2).
    from exp5_adjoint_slope.gendata.make_case import write_case
    inp = write_case("GS", tmp_path / "in")
    ra.edit_namelists(inp, {"window_days": 90, "last_days": 30, "ctrl_period_days": per},
                      {"name": None})
    d = (inp / "data").read_text()
    assert " periodicExternalForcing=.TRUE.," in d
    assert " externForcingPeriod=7776000.0," in d and " externForcingCycle=7776000.0," in d
