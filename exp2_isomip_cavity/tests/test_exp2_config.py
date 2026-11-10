"""Tests for the Exp. 2 MITgcm configuration (Task 2): code/ layout, SIZE.h per resolution,
the rendered namelists and the np / SIZE.h check of the campaign."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from pfcommon import campaign as cp
from exp2_isomip_cavity.gendata import make_case as mc

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "config"
SHARED = CFG / "ocean0" / "code"
LAYOUT = {  # config: (sNx, sNy, nPx, Nx, Ny)
    "ocean0_2km": (48, 41, 5, 240, 41),
    "ocean0_4km": (120, 21, 1, 120, 21),
    "ocean0_4km_mpi": (24, 21, 5, 120, 21),
}


def _par(text, name):
    return cp._size_param(text, name)


@pytest.mark.parametrize("cfg", sorted(LAYOUT))
def test_size_h_per_resolution_and_shared_code(cfg):
    sNx, sNy, nPx, Nx, Ny = LAYOUT[cfg]
    code = CFG / cfg / "code"
    t = (code / "SIZE.h").read_text()
    assert (_par(t, "sNx"), _par(t, "sNy"), _par(t, "nPx"), _par(t, "nPy")) == (sNx, sNy, nPx, 1)
    assert _par(t, "sNx") * _par(t, "nSx") * nPx == Nx and _par(t, "sNy") * _par(t, "nSy") == Ny
    assert _par(t, "Nr") == 36
    # every other file is a link to the one shared code/ (only SIZE.h differs)
    others = sorted(p.name for p in code.iterdir() if p.name != "SIZE.h")
    assert others == sorted(p.name for p in SHARED.iterdir())
    for n in others:
        assert (code / n).is_symlink() and (code / n).resolve() == (SHARED / n).resolve()


def test_grid_matches_cases_yaml():
    spec = mc.load_cases()
    for cfg in ("ocean0_2km", "ocean0_4km"):
        xc, yc = mc.grid_centres(spec["configs"][cfg]["dx_km"])
        _, _, _, Nx, Ny = LAYOUT[cfg]
        assert (xc.size, yc.size + 1) == (Nx, Ny)     # + land row (R12)


def test_fortran_changes_are_marked():
    """Each copied .F says what it is a copy of, and every change is tagged."""
    th = (SHARED / "shelfice_thermodynamics.F").read_text()
    assert "copy of pkg/shelfice/shelfice_thermodynamics.F" in th
    assert len(re.findall(r"^#ifdef SHI_ISOMIP_PLUS_MELT", th, re.MULTILINE)) == 3
    assert len(re.findall(r"^#ifndef SHI_ISOMIP_PLUS_MELT", th, re.MULTILINE)) == 1
    assert th.count("ISOMIP+ (Exp. 2)") >= 5
    pa = (SHARED / "ptracers_apply_forcing.F").read_text()
    assert "copy of pkg/ptracers/ptracers_apply_forcing.F" in pa and "SHELFICE_FORCING_MWTR" in pa
    assert "#define SHI_ISOMIP_PLUS_MELT" in (SHARED / "SHELFICE_OPTIONS.h").read_text()
    assert "#define NONLIN_FRSURF" in (SHARED / "CPP_OPTIONS.h").read_text()
    assert (SHARED / "packages.conf").read_text().split("\n")[2:7] == [
        "gfd", "shelfice", "ptracers", "rbcs", "diagnostics"]


@pytest.fixture
def rendered(tmp_path, synth_geom):
    out = tmp_path / "case"
    meta = mc.make_case("cal2_gT0.03", out, geom_path=synth_geom)   # real templates
    return out, meta


def _nml(text, key):
    m = re.search(rf"^\s*{re.escape(key)}\s*=\s*([^,\n]+)", text, re.MULTILINE | re.IGNORECASE)
    assert m, f"{key} not set"
    return m.group(1).strip()


def test_real_templates_render(rendered):
    out, meta = rendered
    for f in ("data", "data.pkg", "data.shelfice", "data.rbcs", "data.ptracers",
              "data.diagnostics", "eedata"):
        assert (out / f).is_file()
        assert "${" not in (out / f).read_text()
    d = (out / "data").read_text()
    assert int(_nml(d, "nIter0")) == 0
    assert int(_nml(d, "nTimeSteps")) == meta["params"]["nTimeSteps"]
    assert not re.search(r"^\s*(endTime|startTime)\s*=", d, re.MULTILINE | re.IGNORECASE)
    assert float(_nml(d, "hFacMin")) == 0.2 and _nml(d, "useMin4hFacEdges") == ".TRUE."
    # namelist group: useMin4hFacEdges is PARM04 (model/src/ini_parms.F:302); a wrong group
    # stops the model at start-up (smoke run, first attempt)
    code = "\n".join(l for l in d.splitlines() if not l.lstrip().startswith("#"))
    groups = dict(re.findall(r"&(PARM0\d)(.*?)\n &\n", code, re.DOTALL))
    assert "useMin4hFacEdges" in groups["PARM04"] and "useMin4hFacEdges" not in groups["PARM01"]
    # no mixing of z- and r-named inputs (ini_parms.F: 'Cannot mix z, p and r', 2nd attempt)
    assert "delR" in groups["PARM04"] and not re.search(r"\b(delZ|viscAz|diffKz[TS]|hFacMinDz)\b", code)
    for k in ("hFacMin", "hFacInf", "select_rStar", "rhoConstFresh", "useSingleCpuIO"):
        assert k in groups["PARM01"]
    assert _nml(d, "useRealFreshWaterFlux") == ".FALSE."
    assert _nml(d, "tempAdvScheme") == "77" and _nml(d, "saltAdvScheme") == "77"
    assert float(_nml(d, "rhoConst")) == 1027.51 and float(_nml(d, "rhoConstFresh")) == 1000.0
    assert float(_nml(d, "HeatCapacity_Cp")) == 3974.0
    assert _nml(d, "eosType") == "'LINEAR'"
    assert _nml(d, "nonlinFreeSurf") == "4" and _nml(d, "select_rStar") == "2"
    # the pickup checkpoint divides the 2-year run (stage contract of run_stages)
    T = meta["params"]["endTime"]
    assert T % float(_nml(d, "pChkptFreq")) == 0


def test_shelfice_namelist_selects_isomip_plus_law(rendered):
    out, meta = rendered
    s = (out / "data.shelfice").read_text()
    for k, v in (("SHELFICEuseGammaFrict", ".TRUE."), ("SHELFICEconserve", ".TRUE."),
                 ("SHELFICEboundaryLayer", ".TRUE."), ("useISOMIPTD", ".FALSE."),
                 ("SHELFICEadvDiffHeatFlux", ".FALSE."), ("SHI_withBL_realFWflux", ".FALSE."),
                 ("SHELFICEMassStepping", ".FALSE.")):
        assert _nml(s, k) == v, k
    assert float(_nml(s, "SHELFICEkappa")) == 0.0
    assert float(_nml(s, "SHELFICEheatTransCoeff")) == 0.03
    assert float(_nml(s, "SHELFICEsaltTransCoeff")) == pytest.approx(0.03 / 35, rel=1e-15)
    assert float(_nml(s, "shiCdrag")) == 2.5e-3 and float(_nml(s, "SHELFICEDragQuadratic")) == 2.5e-3
    assert float(_nml(s, "SHELFICElatentHeat")) == 3.34e5


def test_rbcs_and_tracer_namelists(rendered):
    out, _ = rendered
    r = (out / "data.rbcs").read_text()
    for k in ("useRBCtemp", "useRBCsalt", "useRBCptracers", "useRBCpTrNum(1)"):
        assert _nml(r, k) == ".TRUE."
    for i in (1, 2, 3):
        assert _nml(r, f"relaxMaskFile({i})") == "'rbcs_mask.bin'"
    assert float(_nml(r, "tauRelaxPTR(1)")) == 8640.0
    code = "\n".join(l for l in r.splitlines() if not l.lstrip().startswith("#"))
    assert "relaxPtracerFile" not in code         # target 0 (rbcs_init_varia.F:72)
    p = (out / "data.ptracers").read_text()
    assert _nml(p, "PTRACERS_numInUse") == "1" and _nml(p, "PTRACERS_Iter0") == "0"
    assert _nml(p, "PTRACERS_advScheme(1)") == "77"
    assert "PTRACERS_EvPrRn" not in p and "PTRACERS_initialFile" not in p


def test_diagnostics_monthly_streams(rendered):
    """R4: monthly 3-D UVEL, VVEL, THETA, SALT, TRAC01 (+ UTRAC01, WVEL) and the SHI 2-D
    fields for the whole run; months of 2628000 s."""
    out, _ = rendered
    t = (out / "data.diagnostics").read_text()
    streams = {}
    for m in re.finditer(r"(?<![\w])fields\(1:\d+,(\d+)\)\s*=\s*((?:'[^']*'\s*,\s*)+)", t):
        streams[int(m.group(1))] = [f.strip() for f in re.findall(r"'([^']*)'", m.group(2))]
    freq = {int(n): float(v) for n, v in re.findall(r"frequency\((\d+)\)\s*=\s*([-0-9.eE+]+)", t)}
    assert freq[1] == freq[2] == 2628000.0
    assert {"SHIfwFlx", "SHIhtFlx", "SHIuStar", "SHIForcT", "ETAN"} <= set(streams[1])
    assert {"UVEL", "VVEL", "WVEL", "THETA", "SALT", "TRAC01", "UTRAC01"} <= set(streams[2])
    assert freq[3] < 0 and {"SHIfwFlx", "SHIuStar", "SHIgammT"} <= set(streams[3])
    assert freq[6] == freq[3] and {"THETA", "SALT", "TRAC01"} <= set(streams[6])
    # a list takes the levels of its first field: never mix 2-D and 3-D fields (smoke
    # run 1: the mixed snapshot list kept only level 1 of THETA / SALT)
    two_d = {"SHIfwFlx", "SHIhtFlx", "SHIuStar", "SHIgammT", "SHIgammS", "SHIForcT",
             "SHIForcS", "SHI_TauX", "SHI_TauY", "ETAN"}
    for n, flds in streams.items():
        assert set(flds) <= two_d or not (set(flds) & two_d), f"list {n} mixes 2-D and 3-D"
    # snapshots at multiples of |snapFreq| (default phase would be -freq/2)
    for n in (3, 6):
        assert re.search(rf"timePhase\({n}\)\s*=\s*0\.", t)
    assert "levels(1:3,5)" in t
    assert _nml(t, "dumpAtLast") == ".TRUE."
    n_fields = sum(len(v) for v in streams.values())
    assert n_fields <= 50                         # DIAGNOSTICS_SIZE.h numperlist per list


def test_campaign_np_must_match_size_h(tmp_path):
    for cfg, (_, _, nPx, _, _) in LAYOUT.items():
        cp.check_np(CFG / cfg / "code" / "SIZE.h", nPx)
        with pytest.raises(RuntimeError, match="nPx"):
            cp.check_np(CFG / cfg / "code" / "SIZE.h", 6 - nPx if nPx == 5 else 5)


def test_load_anomaly_at_snapped_draft(rendered):
    """The load anomaly is evaluated at the snapped draft = the model's R_shelfIce."""
    from pfcommon.mitgcm_io import read_bin
    from exp2_isomip_cavity.gendata import geometry as geo, profiles as pr
    out, meta = rendered
    ny, nx = meta["params"]["ny"], meta["params"]["nx"]
    d = read_bin(out / "icetopo.bin", (ny, nx))
    b = read_bin(out / "bathy.bin", (ny, nx))
    h, r_low, ro_surf = geo.mitgcm_hfac(b, d, np.full(36, 20.0), 0.2)
    np.testing.assert_array_equal(ro_surf[d < 0], d[d < 0])
    T, S = pr.ocean0(-(np.arange(36) + 0.5) * 20.0)
    pf = pr.load_anomaly(T, S, np.full(36, 20.0))
    wet = h.sum(axis=0) > 0
    np.testing.assert_allclose(read_bin(out / "shelficeLoadAnomaly.bin", (ny, nx)),
                               np.where(wet, pr.load_anomaly_at_draft(pf, np.full(36, 20.0), ro_surf), 0.0))


def test_smoke_check_boundary_layer_mean():
    """Losch (2008) average of shelfice_thermodynamics.F:341-389: top cell K (hFacC 0.3)
    plus 0.7 dz of cell K+1."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("pd1_smoke_check", PKG / "analysis" / "pd1_smoke_check.py")
    sc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sc)
    h = np.zeros((4, 1, 2))
    h[1:, 0, 0] = [0.3, 1.0, 1.0]           # column 0: ice base in level 1
    h[2:, 0, 1] = [0.5, 0.2]                # column 1: K+1 thinner than the gap
    f = np.zeros((4, 1, 2))
    f[:, 0, 0] = [9.0, 1.0, 2.0, 5.0]
    f[:, 0, 1] = [9.0, 9.0, 1.0, -3.0]
    ktop = np.array([[1, 2]])
    out = sc.boundary_layer_mean(f, h, ktop, 20.0)
    assert out[0, 0] == pytest.approx((1.0 * 6.0 + 2.0 * 14.0) / 20.0, rel=1e-15)
    assert out[0, 1] == pytest.approx((1.0 * 10.0 - 3.0 * 4.0) / 14.0, rel=1e-15)
    assert sc.boundary_layer_mean(f, h, ktop, 20.0, clip0=True)[0, 1] == pytest.approx(10.0 / 14.0)


# ---- R12: closed lateral boundaries (AD16 sect. 3.1.4: walls on all lateral sides) -------
# MITgcm's Cartesian domain is periodic in x and y (eesupp/src/ini_procs.F:123-125); faces use
# hFacW/S = min of the two adjacent hFacC (useMin4hFacEdges, model/src/ini_masks_etc.F:251-262),
# with the periodic halo (ini_masks_etc.F:402). Checked on generated input, not model output.

def _model_hfac(out, meta):
    from pfcommon.mitgcm_io import read_bin
    from exp2_isomip_cavity.gendata import geometry as geo
    ny, nx, nz = meta["params"]["ny"], meta["params"]["nx"], meta["params"]["nz"]
    b = read_bin(out / "bathy.bin", (ny, nx))
    d = read_bin(out / "icetopo.bin", (ny, nx))
    h = geo.mitgcm_hfac(b, d, np.full(nz, meta["params"]["dz"]), meta["params"]["hFacMin"])[0]
    return (h,) + geo.mitgcm_hfac_faces(h)


@pytest.mark.parametrize("case", ["cal2_gT0.03", "cal4_gT0.03"])
def test_lateral_boundaries_are_closed(tmp_path, synth_geom, case):
    from exp2_testdata import REAL_GEOM
    for geom in [synth_geom] + ([REAL_GEOM] if REAL_GEOM.exists() else []):
        out = tmp_path / f"{case}_{Path(geom).stem}"
        meta = mc.make_case(case, out, geom_path=geom)
        h, hW, hS = _model_hfac(out, meta)
        ny = meta["params"]["ny"]
        assert not hS[:, 0, :].any(), "y = 0 face (periodic link to the last row) is open"
        assert not hS[:, ny - 1, :].any(), "y = 80 km face is open"
        assert not hW[:, :, 0].any(), "x-periodic link face (x = 320 km) is open"
        # no wet cell adjacent across a periodic seam
        assert not ((h[:, 0, :] > 0) & (h[:, -1, :] > 0)).any()
        assert not ((h[:, -1, :] > 0)).any()                  # land row y = 80-82 km
        assert not ((h[:, :, 0] > 0) & (h[:, :, -1] > 0)).any()
        assert not (h[:, :, 0] > 0).any()                     # x = 320-322 km column dry
        # the ISOMIP+ wet domain is unchanged: rows 0..ny-2 on the same centres
        assert meta["params"]["ny_wet_domain"] == ny - 1
        assert int((h.sum(axis=0) > 0).sum()) == meta["actual"]["n_wet_columns"]
        d = (out / "data").read_text()
        assert _nml(d, "ygOrigin") == "0." and f"delY={ny}*" in d.replace(" ", "")


def test_mitgcm_hfac_faces_periodic():
    from exp2_isomip_cavity.gendata import geometry as geo
    h = np.zeros((1, 3, 4))
    h[0, 0, :] = 1.0                       # row 0 wet; row 2 wet only at i = 3
    h[0, 2, 3] = 0.5
    hW, hS = geo.mitgcm_hfac_faces(h)
    assert hS[0, 0, 3] == 0.5 and hS[0, 0, 0] == 0.0          # wraps to row 2
    assert hW[0, 0, 0] == 1.0                                  # wraps to column 3
    assert hW[0, 2, 3] == 0.0 and hS[0, 1, 0] == 0.0


def test_wet_first_column_is_refused():
    from exp2_isomip_cavity.gendata import geometry as geo
    bed = np.full((2, 3), -300.0)
    draft = np.zeros((2, 3))
    with pytest.raises(ValueError, match="x-periodic"):
        geo.close_lateral_boundaries(bed, draft, np.ones((2, 3), bool))
    draft[:, 0] = -300.0
    b, d, w = geo.close_lateral_boundaries(np.where([[0, 1, 1]] * 2, bed, 0.0), np.where([[0, 1, 1]] * 2, 0.0, 0.0),
                                           np.array([[False, True, True]] * 2))
    assert b.shape == (3, 3) and not w[-1].any() and (b[-1] == 0).all() and (d[-1] == 0).all()


def test_r18_check_ini_temp_off_without_freezing(rendered):
    """R18: T = 0 exactly (ramp, Tbot = 0.0) must not stop the model. checkIniTemp only
    stops the run (model/src/ini_theta.F:110-119); it alters theta only together with
    allowFreezing (ini_theta.F:130), which must stay at its .FALSE. default."""
    out, _ = rendered
    d = (out / "data").read_text()
    code = "\n".join(l for l in d.splitlines() if not l.lstrip().startswith("#"))
    parm05 = re.search(r"&PARM05(.*?)\n &", code, re.DOTALL).group(1)
    assert re.search(r"^\s*checkIniTemp\s*=\s*\.FALSE\.", parm05, re.MULTILINE)
    assert not re.search(r"allowFreezing", code, re.IGNORECASE)
