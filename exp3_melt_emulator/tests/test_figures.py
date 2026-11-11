"""Figures 1, 2, 4 and the physics diagnostics (Figure 3) of Exp. 3 on synthetic inputs."""
import dataclasses
import json

import numpy as np
import pytest

pytest.importorskip("cmocean", reason="pfcommon.style needs cmocean (absent in anemoi): "
                    "run this module under pyocean")

from exp3_melt_emulator.analysis import figures as fg  # noqa: E402
from exp3_melt_emulator.analysis import physics_diag as pdg  # noqa: E402
from exp3_melt_emulator.tests.synth import synthetic_runs  # noqa: E402
from exp3_melt_emulator.tests.test_hypotheses import METHODS, _results  # noqa: E402

NOISE = {"pooled_median_rms": 0.157, "median_abs_integrated_rel_diff": 0.013}


@pytest.fixture
def no_npz(monkeypatch):
    """Any attempt to read an npz fails the test."""
    def boom(*a, **k):
        raise AssertionError("figure code read an npz")
    monkeypatch.setattr(np, "load", boom)


def _final_maps(run, methods=METHODS):
    rng = np.random.default_rng(1)
    return {f"{m}__{run.case}__mean": np.where(run.shelf_mask, run.melt, np.nan)
            * (1 + 0.1 * rng.standard_normal(run.melt.shape)) for m in methods}


def test_figures_from_json_only(tmp_path, no_npz):
    cv, fin = _results({"M7": (0.1, 0.3, -0.2)},
                       h3c_points={("M7", "zc_out_0"): (0.02, 0.4)})
    j1 = fg.fig1_method_errors(cv, fin, NOISE, tmp_path)
    assert (tmp_path / "fig1_method_errors.png").exists()
    assert (tmp_path / "fig1_method_errors.pdf").exists()
    saved = json.loads((tmp_path / "fig1_method_errors.json").read_text())
    assert saved["panels"]["final"]["M7"]["median_abs_ire"] == pytest.approx(0.2)
    assert j1["noise_floor"]["rmse_map"] == pytest.approx(0.157)
    j4 = fg.fig4_spread_error(cv, fin, tmp_path, expected=None)
    assert (tmp_path / "fig4_spread_error.png").exists()
    assert [v["case"] for v in j4["M7"]["violations"]] == ["zc_out_0"]

    # Figure 2 takes the maps explicitly (no npz read inside)
    runs = synthetic_runs()
    run = dataclasses.replace(runs[-1], case="final_1", z_c=400.0)
    maps = _final_maps(run)
    j2 = fg.fig2_maps(fin, maps, run, tmp_path)
    assert (tmp_path / "fig2_maps_T2.png").exists()
    assert j2["best_physics"] in ("M1", "M2", "M3")
    assert j2["colour_scale"]["vmin"] <= np.nanmin(run.melt[run.shelf_mask])


def test_fig2_best_physics_by_rmse(tmp_path):
    cv, fin = _results({})
    for e in fin["entries"]:
        if e["case"] == "final_1":
            e["metrics"]["rmse_map"] = {"M1": 3.0, "M2": 2.0, "M3": 5.0}.get(e["method"], 1.0)
    run = dataclasses.replace(synthetic_runs()[-1], case="final_1", z_c=400.0)
    assert fg.fig2_maps(fin, _final_maps(run), run, tmp_path)["best_physics"] == "M2"


def test_powerlaw_fit_recovers_exponent():
    T = np.linspace(-0.5, 2.0, 6)
    m = 3.0 * (T + 1.9) ** 2.5
    fit = pdg.powerlaw_fit(T + 1.9, m)
    assert fit["n"] == pytest.approx(2.5, abs=1e-6)
    assert fit["k"] == pytest.approx(3.0, rel=1e-6)
    assert fit["R2_log"] == pytest.approx(1.0, abs=1e-12)


def test_powerlaw_fit_excludes_nonpositive():
    fit = pdg.powerlaw_fit([-0.1, 1.0, 2.0, 3.0], [1.0, 2.0, 4.0, 8.0], labels="abcd")
    assert fit["excluded"] == ["a"] and fit["n_runs"] == 3


def _physics_inputs():
    """Synthetic factorial: z_c 250/400/550 x Tbot 0.0/1.0/2.0 (2.0 = test), with
    cv loro entries and final entries carrying the integrated melt of each run."""
    base = synthetic_runs(n=2)[0]
    runs, cv, fin, maps = [], [], [], {}
    for zc in (250.0, 400.0, 550.0):
        for T in (0.0, 1.0, 2.0):
            case = f"ens_T{T}_zc{zc:.0f}_s1.0"
            melt = np.where(base.shelf_mask, (T + 2.0) ** 2 * (1 + 0.001 * zc), np.nan)
            T_prof = -1.9 + (T + 1.9) * np.clip((-base.z - zc + 100) / 200, 0, 1)
            r = dataclasses.replace(base, case=case, Tbot=T, z_c=zc, melt=melt,
                                    T_profile=T_prof,
                                    split="test" if T == 2.0 else "train")
            runs.append(r)
            tot = float(melt[base.shelf_mask].sum())
            for m in METHODS:
                pred = np.where(base.shelf_mask, melt * (1.1 if m == "M2" else 0.9), np.nan)
                key = f"{m}__{case}__" + ("mean" if T == 2.0 else "mean_in_sample")
                maps[key] = pred
                e = {"case": case, "method": m, "Tbot": T, "z_c": zc, "s": 1.0,
                     "metrics": {"integrated_pred": float(pred[base.shelf_mask].sum()),
                                 "integrated_true": tot}}
                if T == 2.0:
                    fin.append({**e, "fold": "final"})
                else:
                    cv.append({**e, "kind": "loro", "fold": "loro", "variant": "main"})
    return runs, {"entries": cv}, {"entries": fin}, maps


def test_physics_diag(tmp_path, no_npz):
    runs, cv, fin, maps = _physics_inputs()
    res = pdg.make(runs, cv, fin, maps, tmp_path)
    assert (tmp_path / "fig3_physics.png").exists()
    saved = json.loads((tmp_path / "fig3_physics.json").read_text())
    assert "post hoc" in saved["b_power_law"]["label"].lower()
    a = res["a_integrated_vs_Tbot"]["400"]
    assert a["Tbot"] == [0.0, 1.0, 2.0]
    assert a["M2"] == pytest.approx([1.1 * x for x in a["MITgcm"]])
    lam = res["c_lambert"]["400"]
    # melt = c (T + 2)^2: dM/dT over (0, 1) = 5c, over (1, 2) = 7c, d2M/dT2 = 2c
    c = lam["MITgcm"]["dMdT_ref_p1"] / 5.0
    assert lam["MITgcm"]["dMdT_p1_p2"] == pytest.approx(7 * c)
    assert lam["MITgcm"]["d2MdT2"] == pytest.approx(2 * c)
    assert lam["M2"]["dMdT_p1_p2"] == pytest.approx(1.1 * 7 * c)
    # uniform maps: deep amplification is 1
    assert lam["MITgcm"]["DA_ref_p1"] == pytest.approx(1.0)
    # labelling (fix round 1): in-sample caveat and small-denominator flags in the JSON
    assert "in-sample" in saved["c_lambert"]["caveat"].lower()
    assert lam["MITgcm"]["DA_ref_p1_small_denominator"] == (abs(5 * c) < pdg.DA_DENOM_MIN)
    assert "poorly_constrained" in saved["b_power_law"]["fits"]["TF_shelf_mean"]["pooled_all_runs"]


def test_powerlaw_fit_flags_narrow_forcing():
    assert pdg.powerlaw_fit([1.0, 1.1, 1.2], [1.0, 2.0, 3.0])["poorly_constrained"]
    assert not pdg.powerlaw_fit([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])["poorly_constrained"]


def test_physics_diag_checks_integrated_true(tmp_path):
    runs, cv, fin, maps = _physics_inputs()
    cv["entries"][0]["metrics"]["integrated_true"] *= 1.01
    with pytest.raises(ValueError, match="integrated_true"):
        pdg.make(runs, cv, fin, maps, tmp_path)
