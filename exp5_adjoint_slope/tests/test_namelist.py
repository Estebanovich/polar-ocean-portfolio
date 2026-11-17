"""gendata.namelist.set_param (plan Task 3): set or add one key in a namelist group."""
import pytest

from exp5_adjoint_slope.gendata.namelist import set_param

TXT = """# comment deltaT=1.
 &PARM03
 nIter0=0,
 deltaT=1200.0,
# pickupSuff='old',
 &
 &PARM04
 delX=40*10.E3,
 &
"""


def test_replace_existing_key():
    out = set_param(TXT, "PARM03", "deltaT", "600.0")
    assert " deltaT=600.0,\n" in out and "1200.0" not in out
    assert out.count("deltaT") == 2          # the comment line is untouched


def test_add_missing_key_inside_group():
    out = set_param(TXT, "PARM03", "pickupSuff", "'G0end'")
    i_grp, i_end = out.index("&PARM03"), out.index(" &\n &PARM04")
    i_key = out.index(" pickupSuff='G0end',")
    assert i_grp < i_key < i_end
    assert "# pickupSuff='old'," in out


def test_unknown_group_raises():
    with pytest.raises(ValueError):
        set_param(TXT, "PARM09", "x", "1")


def test_duplicate_key_raises():
    with pytest.raises(ValueError):
        set_param(TXT.replace(" nIter0=0,", " nIter0=0,\n nIter0=1,"), "PARM03", "nIter0", "2")
