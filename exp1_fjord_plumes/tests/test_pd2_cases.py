"""PD2 / ruling R20: convergence and sensitivity cases for config (a).

Grid refinement (dx0 = dz = 10, 5, 2.5 m) and sensitivity at the production
grid (free slip; "Xu-style closure (A4 grid-rescaled)", key sens_xu2012). Checks that make_case takes dz/dx0
from the case, scales the ice-face fields and the source, picks the shared
study templates, and that those templates reproduce production exactly when
given the production values.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from pfcommon.mitgcm_io import read_bin
from exp1_fjord_plumes.gendata.make_case import load_cases, make_case

PKG = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("run_campaign", PKG / "scripts" / "run_campaign.py")
rc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rc)

# name: (dx0, dz, deltaT)   deltaT = 0.5 s/m * min(dx0, dz) (pd2_convergence.md)
CONV = {"conv_dx10": (10.0, 10.0, 5.0),
        "conv_dx5": (5.0, 5.0, 2.5),
        "conv_dx2p5": (2.5, 2.5, 1.25)}
H, RHO0, Q0_LINE, DY = 500.0, 1027.0, 0.03, 10.0


def _data_lines(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")]


@pytest.mark.parametrize("name", list(CONV))
def test_conv_grid_icefront_and_source_scale(tmp_path, name):
    dx0, dz, dt = CONV[name]
    d = make_case(name, tmp_path / name)
    p = json.loads((d / "case.json").read_text())["params"]
    nx, ny, nz = p["nx"], p["ny"], p["nz"]
    assert nz == round(H / dz) and ny == 1
    assert np.allclose(read_bin(d / "delZ.bin", (nz,)), dz)
    dx = read_bin(d / "delX.bin", (nx,))
    assert dx[0] == dx0 and dx.max() == 200.0
    assert p["deltaT"] == dt and p["endTime"] == 86400.0
    assert p["nTimeSteps"] == round(86400.0 / dt)
    R = read_bin(d / "R_icefront.bin", (ny, nx))
    L = read_bin(d / "icefrontlength.bin", (ny, nx))
    assert R[0, 0] == H and np.count_nonzero(R) == 1
    assert L[0, 0] == pytest.approx(1.0 / dx0) and np.count_nonzero(L) == 1
    add = read_bin(d / "addmass.bin", (nz, ny, nx))
    assert add[nz - 1, 0, 0] == pytest.approx(RHO0 * Q0_LINE * DY)  # bottom cell at the face
    assert add.sum() == pytest.approx(add[nz - 1, 0, 0])
    topo = read_bin(d / "topog.bin", (ny, nx))
    assert topo[0, -1] == 0.0 and np.all(topo[0, :-1] == -H)  # land column closes the x-wrap
    rc.check_size(d, PKG / "config" / name / "code" / "SIZE.h")  # SIZE.h matches


@pytest.mark.parametrize("name", list(CONV) + ["sens_freeslip", "sens_xu2012"])
def test_study_cases_selectable_and_physics_unchanged(tmp_path, name):
    cfg = load_cases()["cases"][name]["config"]
    assert rc.select_cases(cfg, name) == [name]
    data = "\n".join(_data_lines(make_case(name, tmp_path / name) / "data"))
    for line in ("eosType='LINEAR',", "tAlpha=3.87e-05,", "rhoConst=1027.0,",
                 "HeatCapacity_Cp=3974.,", "nonHydrostatic=.TRUE.,", "selectAddFluid=1,",
                 "cg3dMaxIters=2000,", "cg3dTargetResidual=1.E-8,"):  # R21: cap 2000
        assert line in data
    # per-step cg3d/CFL validity check (pd2_convergence.md): monitor every step
    p = json.loads((tmp_path / name / "case.json").read_text())["params"]
    assert f"monitorFreq={p['monitorFreq']}," in data
    assert float(p["monitorFreq"]) == p["deltaT"]


def test_production_cases_not_selected_with_study_cases():
    assert len(rc.select_cases("a_line2d", "*")) == 6


def test_study_templates_reproduce_production(tmp_path):
    """With production values, a_line2d_study renders production's data, plus only
    the explicit zero biharmonic lines (MITgcm defaults are 0)."""
    prod = _data_lines(make_case("a_line2d_twolayer_q0.03", tmp_path / "prod") / "data")
    ov = {"template_config": "a_line2d_study", "viscAh": "0.25", "viscAr": "1.E-3",
          "diffKh": "0.25", "diffKr": "2.E-5", "viscA4": "0.", "diffK4": "0.",
          "no_slip_sides": ".TRUE.", "no_slip_bottom": ".TRUE.",
          "cg3dMaxIters": "500", "monitorFreq": "600."}
    study = _data_lines(make_case("a_line2d_twolayer_q0.03", tmp_path / "study", overrides=ov) / "data")
    extra = [ln for ln in study if ln not in prod]
    assert extra == ["viscA4=0.,", "diffK4T=0.,", "diffK4S=0.,"]
    assert [ln for ln in study if ln not in extra] == prod
    for f in ("data.pkg", "data.icefront", "data.obcs", "data.diagnostics", "eedata"):
        assert (tmp_path / "study" / f).read_text() == (tmp_path / "prod" / f).read_text()


def test_sens_freeslip_only_changes_wall_conditions(tmp_path):
    """Against production rendered for 24 h: only the two wall lines change, apart
    from the study-wide R21 cap / per-step monitor and the explicit zero biharmonic."""
    prod = _data_lines(make_case("a_line2d_twolayer_q0.03", tmp_path / "p",
                                 overrides={"endTime": 86400.0}) / "data")
    fs = _data_lines(make_case("sens_freeslip", tmp_path / "s") / "data")
    only_prod = [ln for ln in prod if ln not in fs]
    only_fs = [ln for ln in fs if ln not in prod]
    assert sorted(only_prod) == sorted(["no_slip_sides=.TRUE.,", "no_slip_bottom=.TRUE.,",
                                        "cg3dMaxIters=500,", "monitorFreq=600.,"])
    assert sorted(only_fs) == sorted(["no_slip_sides=.FALSE.,", "no_slip_bottom=.FALSE.,",
                                      "cg3dMaxIters=2000,", "monitorFreq=3.0,",
                                      "viscA4=0.,", "diffK4T=0.,", "diffK4S=0.,"])
    p = json.loads((tmp_path / "s" / "case.json").read_text())["params"]
    assert (p["dx0"], p["dz"], p["deltaT"], p["endTime"]) == (5.0, 10.0, 3.0, 86400.0)


def test_sens_xu2012_mixing(tmp_path):
    d = _data_lines(make_case("sens_xu2012", tmp_path / "x") / "data")
    # "Xu-style closure (A4 grid-rescaled)". Xu et al. 2012, Ann. Glaciol. 53(60), p. 231:
    # nu_v = 0.1 m2/s, kappa = nu (Pr = 1); horizontal biharmonic
    # 300 m4/s at 20 m, scaled to dx0 = 5 m with the same grid damping time:
    # 300 (5/20)^4 = 1.171875 m4/s (literal 300 exceeds the stability limit
    # L^4/(32 dt) = 6.5 m4/s at dx = 5 m, dt = 3 s; mom_calc_visc.F:48-49).
    for line in ("viscAh=0.,", "diffKhT=0.,", "diffKhS=0.,", "viscAr=0.1,", "diffKrT=0.1,",
                 "diffKrS=0.1,", "viscA4=1.171875,", "diffK4T=1.171875,", "diffK4S=1.171875,",
                 "no_slip_sides=.TRUE.,", "no_slip_bottom=.TRUE.,"):
        assert line in d
    assert 1.171875 < 5.0 ** 4 / (32 * 3.0)


def test_case_json_records_git_commit(tmp_path):
    meta = json.loads((make_case("conv_dx10", tmp_path / "g") / "case.json").read_text())
    assert len(meta["git_commit"]) == 40 and int(meta["git_commit"], 16) >= 0
    assert isinstance(meta["git_dirty"], bool)
