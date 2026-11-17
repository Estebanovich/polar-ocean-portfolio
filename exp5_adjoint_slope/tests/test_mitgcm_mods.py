"""Namelists and SIZE.h of the Exp. 5 MITgcm configurations (plan Tasks 2-4)."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from pfcommon.campaign import _size_param, nml_float

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "config"


def _txt(cfg, name):
    return (CFG / cfg / "input" / name).read_text()


def _live(text):
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("#"))


def _has(text, key, val):
    return re.search(rf"^\s*{key}\s*=\s*{re.escape(val)}\s*,?\s*$", _live(text), re.M | re.I)


def test_chan10_gmredi():
    t = _txt("chan10", "data.gmredi")
    # Option A (notes/v4_revision_posthoc.md sec. 4): advective form, floor 1e-6.
    assert _has(t, "GM_AdvForm", ".TRUE.") and _has(t, "GM_Small_Number", "1.E-6")
    # Dettling et al. (2024): kappa_GM = 130/(1 + 800 s_topo) (Eq. 24, Fig. 5) through the bolus
    # map; uniform kappa_Redi = 15 (p. 1682); Gerdes et al. (1991) taper (p. 1678).
    assert _has(t, "GM_background_K", "130.") and _has(t, "GM_isopycK", "15.")
    assert _has(t, "GM_bol2dFile", "'gm_map.bin'") and "GM_iso2dFile" not in _live(t)
    assert _has(t, "GM_taper_scheme", "'gkw91'")


def test_chan10_monthly_diags_have_j_terms():
    # J from monthly means: Eulerian VVELTH/VVEL, bolus GM_vbT, Redi DFyE_TH (diffKhT = 0).
    t = _live(_txt("chan10", "data.diagnostics"))
    for f in ("'VVELTH  '", "'VVEL    '", "'GM_vbT  '", "'DFyE_TH '"):
        assert f in t


DETTLING_PHYSICS = [("viscAh", "12."), ("viscA4Grid", "0.1"), ("viscC4leith", "1."),
                    ("viscC4leithD", "1."), ("viscAr", "3.E-4"), ("bottomDragLinear", "1.E-3"),
                    ("diffKrT", "5.E-6"), ("diffKrS", "5.E-6"),
                    ("saltFluxFile", "'salt_flux.bin'"), ("lambdaThetaFile", "'lambda_t.bin'"),
                    ("lambdaSaltFile", "'lambda_s.bin'")]


@pytest.mark.parametrize("cfg", ["chan10", "chan2p5"])
def test_dettling_physics(cfg):
    # Dettling et al. (2024) Table 1; notes/pd1b_config_dettling.md.
    t = _txt(cfg, "data")
    assert nml_float(t, "f0") == pytest.approx(-1.31e-4)
    for k, v in DETTLING_PHYSICS:
        assert _has(t, k, v), (cfg, k)
    live = _live(t)
    # tau*ClimRelax must stay: it switches the relaxation on (test_surface_relaxation_switched_on).
    for gone in ("bottomDragQuadratic", "viscC2Leith"):
        assert not re.search(rf"^\s*{gone}\s*=", live, re.M), (cfg, gone)


def test_chan10_data_core():
    t = _txt("chan10", "data")
    assert _has(t, "eosType", "'LINEAR'") and _has(t, "useSingleCpuIO", ".TRUE.")
    assert _has(t, "readBinaryPrec", "64") and _has(t, "writeBinaryPrec", "64")
    assert _has(t, "tempAdvScheme", "30") and _has(t, "saltAdvScheme", "30")
    assert _has(t, "delX", "40*10.E3") and _has(t, "delY", "45*10.E3")


@pytest.mark.parametrize("cfg", ["chan10", "chan2p5"])
def test_rbcs_sponge_dettling(cfg):
    # Table 1: T and S restored in 54 days, velocity to zero in 27 days, 50 km sponge.
    assert _has(_txt(cfg, "data.pkg"), "useRBCS", ".TRUE.")
    t = _txt(cfg, "data.rbcs")
    assert _has(t, "useRBCtemp", ".TRUE.") and _has(t, "useRBCsalt", ".TRUE.")
    assert _has(t, "useRBCuVel", ".TRUE.") and _has(t, "useRBCvVel", ".TRUE.")
    assert _has(t, "tauRelaxT", "4665600.") and _has(t, "tauRelaxS", "4665600.")
    assert _has(t, "tauRelaxU", "2332800.") and _has(t, "tauRelaxV", "2332800.")
    assert "relaxMaskUFile='rbcs_mask_uv.bin'" in t and "relaxMaskVFile='rbcs_mask_uv.bin'" in t


@pytest.mark.parametrize("cfg,nx,ny", [("chan10", 40, 45), ("chan2p5", 160, 180)])
def test_size_h(cfg, nx, ny):
    t = (CFG / cfg / "code" / "SIZE.h").read_text()
    p = {k: _size_param(t, k) for k in ("sNx", "sNy", "nSx", "nSy", "nPx", "nPy", "Nr")}
    assert p["sNx"] * p["nSx"] * p["nPx"] == nx and p["sNy"] * p["nSy"] * p["nPy"] == ny
    assert p["Nr"] == 35 and p["nPx"] * p["nPy"] <= 5


def test_chan2p5_no_gm():
    assert "gmredi" not in (CFG / "chan2p5" / "code" / "packages.conf").read_text()
    t = _txt("chan2p5", "data")
    assert _has(t, "delX", "160*2.5E3") and _has(t, "delY", "180*2.5E3")
    assert "useGMRedi" not in _txt("chan2p5", "data.pkg")


def test_tap_cost_row_matches_registry():
    import yaml
    from exp5_adjoint_slope.gendata import fields
    from exp5_adjoint_slope.gendata.make_case import grid_of, load_case
    h = (CFG / "chan10_tap" / "code_tap" / "EXP5_COST.h").read_text()
    jq = int(re.search(r"jqb_glob\s*=\s*(\d+)", h).group(1))
    p = load_case("GS")
    assert jq == fields.j_qb(grid_of(p), p["y_qb_m"]) + 1
    assert re.search(r"thetaF\s*=\s*-1\.9", h) and p["theta_f"] == -1.9


def test_fixed_form_fortran_fits_72_columns():
    # Fixed-form Fortran ignores columns > 72 after cpp: _RL expands to Real*8, so a
    # declaration that fit in the source truncated 'psiBelow' to 'psiBe' and broke the build.
    for f in (CFG / "chan10_tap" / "code_tap").glob("*.F"):
        for n, ln in enumerate(_fortran_lines(f.read_text()), 1):
            assert len(ln) <= 72, f"{f.name}:{n} has {len(ln)} columns after cpp"


def _fortran_lines(text):
    out = []
    for ln in text.splitlines():
        if ln[:1] in ("C", "c", "!", "#") or not ln.strip():
            continue
        out.append(ln.replace("_RL", "Real*8").replace("_RS", "Real*8").replace(" _d ", "D"))
    return out


@pytest.mark.parametrize("cfg", ["chan10", "chan2p5"])
def test_surface_relaxation_switched_on(cfg):
    # MITgcm applies surface relaxation only if tau*ClimRelax > 0 (model/src/set_parms.F:250-253,
    # external_forcing_surf.F:173); lambda*File then overrides the rate. Without tau, G0 ran 30
    # years with no shelf heat sink (2026-11-16).
    from pfcommon.campaign import nml_float
    data = (PKG / "config" / cfg / "input" / "data").read_text()
    assert "lambdaThetaFile" in data and "lambdaSaltFile" in data
    assert nml_float(data, "tauThetaClimRelax") > 0
    assert nml_float(data, "tauSaltClimRelax") > 0


def test_the_main_loop_copy_only_adds_binomial_checkpointing():
    # GM30 (2026-11-16): without checkpointing the Tapenade adjoint stores the whole trajectory
    # (8.06 GB for 30 days; GA would need ~97 GB). The copy in code_tap differs from the MITgcm
    # source only by the $AD BINOMIAL-CKP directive in front of the time loop
    # (tools/TAP_support/ADFirstAidKit/README.md: "Triggered by Tapenade directive $AD BINOMIAL-CKP").
    src = Path("/Volumes/ESTEBAN_MIT/MITgcm/model/src/the_main_loop.F")
    if not src.exists():
        pytest.skip("MITgcm source not mounted")
    ours = (CFG / "chan10_tap" / "code_tap" / "the_main_loop.F").read_text().splitlines()
    orig = src.read_text().splitlines()
    added = [l for l in ours if l not in orig]
    assert [l.split()[0] for l in added if not l.startswith("C ")] == ["C$AD"]
    d = [l for l in added if l.startswith("C$AD")][0]
    assert d.split()[1] == "BINOMIAL-CKP" and d.split()[2] == "nTimeSteps+1"
    i = ours.index(d)
    assert ours[i + 1].strip() == "DO iloop = 1, nTimeSteps"
    assert all(len(l) <= 72 for l in added)
    assert [l for l in orig if l not in ours] == []
