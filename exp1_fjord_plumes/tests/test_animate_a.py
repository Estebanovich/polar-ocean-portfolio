"""animate_a: frame selection (snapshots with both THETA and WVEL, t > 0)."""
from exp1_fjord_plumes.analysis.animate_a import frames


def test_frames_need_theta_and_wvel_and_skip_t0(tmp_path):
    for name, its in (("THETA", (0, 1200, 2400, 3600)), ("WVEL", (0, 1200, 3600))):
        for it in its:
            (tmp_path / f"{name}.{it:010d}.meta").write_text("")
    assert frames(tmp_path) == [1200, 3600]
