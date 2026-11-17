"""analysis.pd2_eval: assembly of V1d points (best eps per point, adjoint vs FD and vs TLM)."""
import pytest

from exp5_adjoint_slope.analysis import pd2_eval as pe


def _g(grad, fd, rel):
    return {"grad": grad, "fd_grad": fd, "rel_err": rel}


def test_v1d_points_best_eps_and_tlm():
    adj = {"xx_fu_i26j17k1_e0": _g(1.0, 1.02, 0.02), "xx_fu_i26j17k1_e1": _g(1.0, 1.001, 0.001),
           "xx_theta_i1j17k18_e0": _g(2.0, 2.0, 1e-6)}
    tlm = {"xx_fu_i26j17k1_e0": _g(1.0000001, 1.001, 0.001), "xx_theta_i1j17k18_e0": _g(2.1, 2.0, 0.05)}
    pts = pe.v1d_points(adj, tlm)
    assert [p["point"] for p in pts] == ["xx_fu_i26j17k1", "xx_theta_i1j17k18"]
    assert pts[0]["rel_err_fd"] == pytest.approx(0.001) and pts[0]["eps_run"] == "e1"
    assert pts[0]["rel_err_tlm"] == pytest.approx(1e-7, rel=1e-3)
    assert pts[1]["rel_err"] == pytest.approx(abs(2.0 - 2.1) / 2.1)      # worst of FD and TLM
    assert pe.v1d_points(adj, {})[1]["rel_err_tlm"] is None


def test_drop_periodic_end_requires_a_zero_last_record():
    import numpy as np
    S = np.ones((2, 4, 5))
    S[:, -1] = 0.0
    assert pe.drop_periodic_end(S).shape == (2, 3, 5)
    S[0, -1, 2] = 1e-3
    with pytest.raises(ValueError, match="last record"):
        pe.drop_periodic_end(S)
