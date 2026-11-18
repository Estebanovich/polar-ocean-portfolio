import re
from pathlib import Path

import pytest
import yaml

from pfcommon import campaign

PKG = Path(__file__).resolve().parents[1]
SRC = Path("/Volumes/ESTEBAN_MIT/MITgcm/pkg/gchem/gchem_forcing_sep.F")
CFGS = ("g50", "g50p5", "f10")


def _code(cfg, name):
    return (PKG / "config" / cfg / "code" / name).read_text()


@pytest.mark.parametrize("cfg", CFGS)
def test_gchem_forcing_sep_single_exp6_block(cfg):
    ours = _code(cfg, "gchem_forcing_sep.F")
    blocks = re.findall(r"C EXP6 begin\n.*?C EXP6 end\n", ours, re.S)
    assert len(blocks) == 1
    assert ours.replace(blocks[0], "") == SRC.read_text()
    assert "CALL EXP6_PREFORMED( bi, bj, myThid )" in blocks[0]
    sub = _code(cfg, "exp6_preformed.F")
    assert "pTracer(i,j,1,bi,bj,6) = pTracer(i,j,1,bi,bj,3)" in sub


@pytest.mark.parametrize("cfg", CFGS)
def test_dic_options(cfg):
    o = _code(cfg, "DIC_OPTIONS.h")
    assert "#undef ALLOW_FE" in o and "#define READ_PAR" in o and "#define ALLOW_FE" not in o


@pytest.mark.parametrize("cfg", CFGS)
def test_ptracers_size(cfg):
    assert re.search(r"PTRACERS_num\s*=\s*6\b", _code(cfg, "PTRACERS_SIZE.h"))


@pytest.mark.parametrize("cfg", CFGS)
def test_layers_size_matches_bounds(cfg):
    n = int(re.search(r"Nlayers\s*=\s*(\d+)", _code(cfg, "LAYERS_SIZE.h")).group(1))
    t = (PKG / "config" / "input_templates" / "data.layers").read_text()
    assert n == int(re.search(r"layers_bounds\(1:(\d+),1\)", t).group(1)) - 1


@pytest.mark.parametrize("cfg", CFGS)
def test_size_matches_registry(cfg):
    c = yaml.safe_load((PKG / "cases.yaml").read_text())["configs"][cfg]
    s = _code(cfg, "SIZE.h")
    get = lambda k: int(re.search(rf"{k}\s*=\s*(\d+)", s).group(1))
    assert get("sNx") * get("nSx") * get("nPx") == c["nx"]
    assert get("sNy") * get("nSy") * get("nPy") == c["ny"]
    assert get("Nr") == c["nr"] and get("nPx") * get("nPy") == c["np"]
    pk = _code(cfg, "packages.conf").split()
    assert {"gfd", "gmredi", "rbcs", "layers", "diagnostics", "ptracers", "gchem", "dic"} <= set(pk)


DIC = Path("/Volumes/ESTEBAN_MIT/MITgcm/pkg/dic")


@pytest.mark.parametrize("cfg", CFGS)
@pytest.mark.parametrize("name,needle", [
    ("dic_biotic_forcing.F", "CALL DIAGNOSTICS_FILL(GDIC   ,'DICGDIC ',0,Nr,2,bi,bj,myThid)"),
    ("dic_diagnostics_init.F", "diagName  = 'DICGPO4 '")])
def test_dic_copies_single_exp6_block(cfg, name, needle):
    ours = _code(cfg, name)
    blocks = re.findall(r"C EXP6 begin\n.*?C EXP6 end\n", ours, re.S)
    assert len(blocks) == 1 and needle in blocks[0]
    assert ours.replace(blocks[0], "") == (DIC / name).read_text()
