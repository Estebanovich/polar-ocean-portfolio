"""Fixtures for Exp. 2 tests: synthetic 1 km ISOMIP+ geometry and placeholder templates."""
from __future__ import annotations

import pytest

from exp2_testdata import synthetic_fields, write_geom


@pytest.fixture(scope="session")
def synth_geom(tmp_path_factory):
    d = tmp_path_factory.mktemp("geom")
    return write_geom(d / "synthetic_geom.nc", synthetic_fields())


@pytest.fixture
def templates(tmp_path):
    """Placeholder templates (ruling R1): the real namelists are written in Task 2."""
    d = tmp_path / "templates"
    d.mkdir()
    (d / "data").write_text(
        " &PARM01\n tRef = ${nz}*${tRef},\n rhoNil = ${rho_ref},\n f0 = ${f0},\n &\n"
        " &PARM03\n nIter0 = ${nIter0},\n nTimeSteps = ${nTimeSteps},\n deltaT = ${deltaT},\n &\n"
        " &PARM05\n bathyFile = 'bathy.bin',\n hydrogThetaFile = 'T.init.bin',\n"
        " hydrogSaltFile = 'S.init.bin',\n delX = ${nx}*${dx},\n delY = ${ny}*${dy},\n"
        " delZ = ${nz}*${dz},\n &\n")
    (d / "data.shelfice").write_text(
        " &SHELFICE_PARM01\n SHELFICEtopoFile = 'icetopo.bin',\n"
        " SHELFICEloadAnomalyFile = 'shelficeLoadAnomaly.bin',\n"
        " SHELFICEheatTransCoeff = ${gammaT},\n SHELFICEsaltTransCoeff = ${gammaS},\n &\n")
    (d / "data.rbcs").write_text(
        " &RBCS_PARM01\n relaxMaskFile(1) = 'rbcs_mask.bin',\n relaxTFile = 'T.restore.bin',\n"
        " relaxSFile = 'S.restore.bin',\n tauRelaxT = ${tauRelax},\n &\n")
    return d
