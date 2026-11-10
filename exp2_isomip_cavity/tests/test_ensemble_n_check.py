"""Tests for analysis/ensemble_n_check.py (pd0 N1-N3 summary of the production queue)."""
from __future__ import annotations

import json

from exp2_isomip_cavity.analysis import ensemble_n_check as enc
from test_metrics import make_run


def test_failed_run_is_a_gap_and_ok_run_passes_n3(tmp_path, capsys):
    make_run(tmp_path / "ens_a_r21", nrec=24)
    make_run(tmp_path / "ens_b_r21", nrec=24, status="FAILED")
    res = enc.check_ensemble(tmp_path, "_r21", names=["ens_a", "ens_b"])
    assert res["gaps_N3"] == ["ens_b"]
    assert res["runs"]["ens_a"]["N3"]["pass"] and not res["runs"]["ens_a"]["gap"]
    assert set(res["runs"]["ens_a"]) >= {"N1", "N2"}
    assert "N1" not in res["runs"]["ens_b"]                    # nothing read from a gap


def test_missing_run_is_a_gap_and_cli_never_fails(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr(enc, "ensemble_cases", lambda spec: ["ens_missing"])
    out = tmp_path / "s.json"
    assert enc.main(["--runs-root", str(tmp_path), "--out", str(out)]) == 0
    res = json.loads(out.read_text())
    assert res["gaps_N3"] == ["ens_missing"]
    assert "SUMMARY: 1 runs; N3 gaps 1" in capsys.readouterr().out
