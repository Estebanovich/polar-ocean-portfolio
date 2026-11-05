from pathlib import Path

import numpy as np

from pfcommon.mitgcm_io import read_bin, read_mds, write_bin


def test_bin_roundtrip_big_endian(tmp_path):
    rng = np.random.default_rng(0)
    a = rng.standard_normal((3, 4, 5))
    p = tmp_path / "a.bin"
    write_bin(p, a)
    assert np.array_equal(read_bin(p, a.shape), a)
    raw = p.read_bytes()
    assert raw[:8] == np.float64(a[0, 0, 0]).byteswap().tobytes()
    assert len(raw) == a.size * 8


def test_read_mds_parses_meta(tmp_path):
    arr = np.arange(2 * 3 * 4 * 5, dtype=float).reshape(2, 3, 4, 5)
    prefix = tmp_path / "Eta"
    stem = tmp_path / "Eta.0000000020"
    Path(f"{stem}.data").write_bytes(arr.astype(">f8").tobytes())
    Path(f"{stem}.meta").write_text(
        " nDims = [   3 ];\n"
        " dimList = [\n"
        "    5,    1,    5,\n"
        "    4,    1,    4,\n"
        "    3,    1,    3\n"
        " ];\n"
        " dataprec = [ 'float64' ];\n"
        " nrecords = [     2 ];\n"
        " timeStepNumber = [ 20 ];\n"
        " fldList = {\n"
        " 'Theta   ' 'Salt    '\n"
        " };\n"
    )
    out, meta = read_mds(prefix, itr=20)
    assert out.shape == (2, 3, 4, 5)
    assert np.array_equal(out, arr)
    assert meta["fldList"] == ["Theta", "Salt"]
    out2, _ = read_mds(prefix)  # latest iteration
    assert np.array_equal(out2, arr)
