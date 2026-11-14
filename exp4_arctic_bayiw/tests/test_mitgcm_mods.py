"""The modified MITgcm source copies of config/fjord_vmf/code (PD1, 2026-11-12) differ from the
originals only in the documented blocks (Task 4 PD1 finalization).

Each copy starts with a header (``C--- Exp. 4 header`` ... ``C--- Exp. 4 header end``) naming
the original and its sha256; every change sits between ``C--- Exp. 4 BEGIN`` and
``C--- Exp. 4 END`` lines. Removing the header and those blocks must give the original byte
for byte. The originals are read from the read-only MITgcm tree (skipped if not mounted).
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pytest
import yaml

PKG = Path(__file__).resolve().parents[1]
CODE = PKG / "config" / "fjord_vmf" / "code"
MITGCM = Path("/Volumes/ESTEBAN_MIT/MITgcm")
COPIES = {"obcs_check.F": "pkg/obcs/obcs_check.F",
          "seaice_dynsolver.F": "pkg/seaice/seaice_dynsolver.F",
          "obcs_add_tides.F": "pkg/obcs/obcs_add_tides.F"}
mounted = pytest.mark.skipif(not (MITGCM / "pkg").exists(), reason="MITgcm source not mounted")


def strip_mods(text: str) -> tuple[str, str, list[str]]:
    """(text without header and Exp. 4 blocks, header, list of block bodies)."""
    head, sep, rest = text.partition("C--- Exp. 4 header end\n")
    assert sep, "no header"
    blocks = re.findall(r"C--- Exp\. 4 BEGIN[^\n]*\n(.*?)C--- Exp\. 4 END\n", rest, re.S)
    clean = re.sub(r"C--- Exp\. 4 BEGIN[^\n]*\n.*?C--- Exp\. 4 END\n", "", rest, flags=re.S)
    return clean, head, blocks


@mounted
@pytest.mark.parametrize("name,src", sorted(COPIES.items()))
def test_copy_equals_original_outside_blocks(name, src):
    text = (CODE / name).read_text()
    clean, head, blocks = strip_mods(text)
    orig = (MITGCM / src).read_bytes()
    assert clean.encode() == orig
    assert src in head
    assert hashlib.sha256(orig).hexdigest() in head
    assert blocks


def test_blocks_are_the_documented_changes():
    _, _, b = strip_mods((CODE / "obcs_check.F").read_text())
    assert b[0].strip().startswith("#ifdef ALLOW_SEAICE") and "SEAICE_OPTIONS.h" in b[0]
    assert [x.strip() for x in b[1:]] == ["#ifndef DISABLE_SEAICE_OBCS",
                                         "#endif /* DISABLE_SEAICE_OBCS */",
                                         "#ifndef DISABLE_SEAICE_OBCS",
                                         "#endif /* DISABLE_SEAICE_OBCS */"]
    _, _, b = strip_mods((CODE / "seaice_dynsolver.F").read_text())
    assert len(b) == 1 and "CALL SEAICE_OCEANDRAG_COEFFS(" in b[0]
    assert ".NOT.SEAICEuseDYNAMICS .OR. SEAICEuseFREEDRIFT" in b[0]
    _, _, b = strip_mods((CODE / "obcs_add_tides.F").read_text())
    assert len(b) == 2 and "rampFac" in b[0] and "COS( PI*myTime/EXP4_TIDE_RAMP_S )" in b[1]


def test_options_enable_the_mods():
    so = (CODE / "SEAICE_OPTIONS.h").read_text()
    assert re.search(r"^#define DISABLE_SEAICE_OBCS", so, re.M)
    oo = (CODE / "OBCS_OPTIONS.h").read_text()
    m = re.search(r"^#define EXP4_TIDE_RAMP_S ([0-9.]+) _d 0", oo, re.M)
    spec = yaml.safe_load((PKG / "cases.yaml").read_text())["configs"]["fjord_vmf"]["tide"]
    assert m and float(m.group(1)) == spec["ramp_s"] == spec["period_s"]
