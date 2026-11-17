"""exp5_adjoint_slope/scripts/run_campaign.py --dry-run planning (plan Task 2)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
SCRIPT = PKG / "scripts" / "run_campaign.py"


def _dry(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args, "--dry-run"], capture_output=True,
                          text=True)


def test_g0_stages_sum_to_ntimesteps():
    r = _dry("--cases", "G0", "--np", "5")
    assert r.returncode == 0, r.stderr
    assert r.stdout.split() == ["G0"]
    import yaml
    reg = yaml.safe_load((PKG / "cases.yaml").read_text())
    n, k = reg["cases"]["G0"]["n_timesteps"], reg["cases"]["G0"]["n_stages"]
    assert f"stages={k}" in r.stderr and f"nTimeSteps={n}" in r.stderr
    assert n % k == 0


def test_rejects_np_over_5():
    r = _dry("--cases", "G0", "--np", "6")
    assert r.returncode != 0 and "np" in r.stderr


def test_rejects_np_not_matching_size_h():
    r = _dry("--cases", "G0", "--np", "1")
    assert r.returncode != 0 and "nPx" in r.stderr


def test_prepare_member_perturbs_copied_pickup(tmp_path):
    # FC1..3: the F1b pickup is copied, then theta gets member-k noise (<= 1e-3 degC, wet only).
    import importlib.util
    import numpy as np
    from exp5_adjoint_slope.gendata import interp, make_case
    spec = importlib.util.spec_from_file_location("run_campaign", SCRIPT)
    rc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rc)
    g = make_case.grid_of(make_case.load_case("F0"))
    T = np.full((g.nr, g.ny, g.nx), -1.0)
    S = np.full((g.nr, g.ny, g.nx), 34.5)
    p = make_case.load_case("FC1")
    assert p["pickup_from"] == "F1b" and int(p["pickup_iter"]) == 328500
    pk = f"pickup.{int(p['pickup_iter']):010d}"
    src = tmp_path / "runs" / p["pickup_from"]
    src.mkdir(parents=True)
    z = np.zeros_like(T)
    np.concatenate([z, z, T, S, z, z, np.zeros((3, g.ny, g.nx))]).astype(">f8").tofile(
        src / f"{pk}.data")
    (src / f"{pk}.meta").write_text(
        " fldList = {\n 'Uvel    ' 'Vvel    ' 'Theta   ' 'Salt    ' 'GuNm1   ' 'GvNm1   '"
        " 'EtaN    ' 'dEtaHdt ' 'EtaH    '\n };\n")
    inp = rc.prepare("FC1", p, tmp_path / "FC1", tmp_path / "runs")
    t, s = interp.pickup_ts(inp / f"{pk}.data", g.nr, g.ny, g.nx)
    d = t - T
    assert 0 < np.abs(d).max() <= 1e-3 and np.array_equal(s, S)
    t0, _ = interp.pickup_ts(src / f"{pk}.data", g.nr, g.ny, g.nx)
    assert np.array_equal(t0, T)                                     # source untouched
