"""analysis.figures (plan Task 7): every figure renders from synthetic JSON without error."""
import numpy as np
import pytest

from exp5_adjoint_slope.analysis import figures as fg

Y = ((np.arange(45) + 0.5) * 10e3).tolist()
Z = (-np.linspace(5, 2900, 35)).tolist()


def _ms():
    th = np.tile(np.linspace(-1.9, 0.8, 35)[:, None], (1, 45))
    return {"y_m": Y, "z_m": Z, "depth_m": np.linspace(500, 3000, 45).tolist(),
            "tau_x": (-0.075 * np.sin(np.linspace(0, np.pi, 45))).tolist(), "theta": th.tolist(),
            "salt": (th * 0 + 34.5).tolist(), "sigma_anom": (th * 0.1).tolist(), "y_qb_m": 170e3,
            "j_qb": 17, "J": {"J_TW": -2e-4, "J_resolved_TW": -2e-3, "J_bolus_TW": 1.7e-3,
                              "J_redi_TW": 2e-4, "shelf_theta_degC": -1.9, "shelf_salt": 34.4}}


def _sens():
    r = np.random.default_rng(0)
    return {"leads_days": [60.0, 30.0, 0.0], "y_m": Y,
            "dJdtaux_rows_mean": r.normal(size=(3, 45)).tolist(),
            "dJdtaux_rows_se": np.full((3, 45), 0.1).tolist(), "fp_pattern": np.zeros(45).tolist()}


def _val():
    pts = [{"point": f"xx_fu_i{i}j17k1", "grad_adj": 1e-3, "fd_grad": 1.0001e-3, "grad_tlm": 1e-3,
            "rel_err": 1e-4} for i in range(10)]
    w = {k: {"dj_uniform": -1e-3, "dj_dipole": 2e-3, "pass": True, "paired_diff": 1e-3,
             "paired_se": 1e-4, "significant": True} for k in ("37.5", "18.75", "75")}
    return {"A": {"V1d": {"pass": True, "value": 1e-4, "points": pts}},
            "B": {"V3e": {"pred_members": [1e-3] * 10, "fwd_members": [1.1e-3] * 10},
                  "H5a": {"widths": w, "robust": True, "direct_37.5km_years26-30": {
                      "i": {"fwd": [-9e-4, 1e-4]}, "ii": {"fwd": [1.7e-3, 1e-4]}}},
                  "H5c": {"leads_days": [0, 30, 60], "sens_mean": [0.01, 0.03, 0.05],
                          "sens_se": [0.002] * 3, "verdict": "not testable beyond 60 days"}}}


def _diag():
    return {"adjoint_growth": {"30": {"max_abs_dJdtheta0": 8e-3}, "90": {"max_abs_dJdtheta0": 2e-3},
                               "180": {"max_abs_dJdtheta0": 2e11}, "360": {"max_abs_dJdtheta0": 9e43}},
            "ensemble_J_TW": {f"E{k}": np.linspace(0, 1e-3 * k, 12).tolist() for k in range(5)},
            "ensemble_std_TW": np.linspace(1e-8, 1e-3, 12).tolist()}


def _xc():
    return {"J_control_TW": {"FC1": 4e-3, "FC2": 3e-3, "FC3": 1e-3},
            "J_pert_TW": {"FP1": 3.7e-3, "FP2": 4.1e-3, "FP3": 9.6e-3},
            "H5b": {"dj": 3.1e-3, "sigma_delta": 2.1e-3, "dj_pred": 1.4e-3, "pass": False},
            "dj_pred_se_TW": 1.6e-4, "V5_eddy_fraction": {"FC1": 1.6, "FC2": 2.1, "FC3": 1.8}}


@pytest.mark.parametrize("name", ["fig0_schematic", "fig1_mean_state", "fig2_sensitivity",
                                  "fig3_validation", "fig4_hypotheses", "fig5_crosscheck"])
def test_each_figure_renders(tmp_path, name):
    data = {"mean_state": _ms(), "sensitivity": _sens(), "validation": _val(), "diagnostics": _diag(),
            "crosscheck": _xc()}
    out = fg.FIGS[name](data, tmp_path / name)
    assert all(p.exists() and p.stat().st_size > 1000 for p in out)
