"""Tests for exp4_arctic_bayiw.analysis.pd1_timing (Task 4): parsing of the MITgcm monitor."""
from __future__ import annotations

import json

from exp4_arctic_bayiw.analysis import pd1_timing as pt

LOG = """(PID.TID 0000.0001) %MON time_tsnumber                =                     0
(PID.TID 0000.0001) %MON advcfl_uvel_max              =   0.0000000000000E+00
(PID.TID 0000.0001) %MON advcfl_W_hf_max              =   0.0000000000000E+00
(PID.TID 0000.0001) %MON dynstat_wvel_max             =  -0.0000000000000E+00
(PID.TID 0000.0001) %MON time_tsnumber                =                    20
(PID.TID 0000.0001) %MON advcfl_uvel_max              =   1.5000000000000E-02
(PID.TID 0000.0001) %MON advcfl_W_hf_max              =   2.5000000000000E-01
(PID.TID 0000.0001) %MON dynstat_wvel_max             =   1.0000000000000E-03
(PID.TID 0000.0001) %MON time_tsnumber                =                    40
(PID.TID 0000.0001) %MON advcfl_uvel_max              =   2.0000000000000E-02
(PID.TID 0000.0001) %MON advcfl_W_hf_max              =   7.0000000000000E-01
(PID.TID 0000.0001) %MON dynstat_wvel_max             =   8.0000000000000E-03
(PID.TID 0000.0001)   Seconds in section "MAIN LOOP           [THE_MAIN_LOOP]":
(PID.TID 0000.0001)           User time:   141.4
(PID.TID 0000.0001)         System time:  0.93
(PID.TID 0000.0001)     Wall clock time:   181.9104
"""


def test_monitor_series_and_main_loop():
    m = pt.monitor_series(LOG)
    assert m["advcfl_W_hf_max"] == [(0, 0.0), (20, 0.25), (40, 0.7)]
    assert m["dynstat_wvel_max"][-1] == (40, 8e-3)
    assert pt.main_loop_wall(LOG) == 181.9104


def test_stage_info_cfl_and_noise_onset(tmp_path):
    (tmp_path / "STDOUT.0000").write_text(LOG)
    (tmp_path / "run_status.json").write_text(json.dumps(
        {"state": "OK", "last_iter": 40, "n_timesteps": 40, "wall_s": 20}))
    s = pt.stage_info(tmp_path, dt=60.0)
    assert s["cfl_max"] == 0.7 and s["cfl_max_iter"] == 40 and s["cfl_max_key"] == "advcfl_W_hf_max"
    assert s["noise_onset_iter"] == 40 and s["noise_onset_h"] == 40 * 60 / 3600
    assert s["s_per_step"] == 0.5
