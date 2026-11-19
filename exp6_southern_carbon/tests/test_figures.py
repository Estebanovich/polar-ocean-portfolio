"""analysis.figures (plan Task 10): every figure renders from synthetic JSON; main() skips figures
whose source JSON is missing (coarse README first, fine run pending)."""
import json

import numpy as np
import pytest

from exp6_southern_carbon.analysis import figures as fg
from exp6_southern_carbon.analysis import hypotheses

Y = (1900, 1980, 2020, 2060)
NY, NR = 40, 49


def _figdata():
    y = ((np.arange(NY) + 0.5) * 50).tolist()
    z = np.linspace(5, 3900, NR)
    dic = 2100 + 150 * (1 - np.exp(-z / 500))[:, None] * np.ones(NY)
    years = list(range(1850, 2101))
    ser = {"GC": {"years": years, "flux": [0.45] * len(years)},
           "GT": {"years": years, "flux": np.linspace(0.45, 3.0, len(years)).tolist()}}
    for b in Y:
        for p in ("GC", "GT"):
            yy = list(range(b, b + 25))
            ser[f"{p}_w{b}"] = {"years": yy, "flux": (np.array(ser[p]["flux"][b - 1850:b - 1825]) - 0.1).tolist()}
    return {"geometry": {"y_km": y, "depth_m": np.r_[0, np.full(NY - 1, 4000.0)].tolist(),
                         "tau_x": (0.2 * np.sin(np.linspace(0, np.pi, NY))).tolist(),
                         "sst_mean": np.linspace(-1.85, 9.85, NY).tolist(), "sst_amp": 1.5,
                         "sponge_y_km": [1900.0, 2000.0], "R_edge_km": 1000.0},
            "sections": {"z_m": z.tolist(), "y_km": y, "dic": {"1900": dic.tolist(), "2060": (dic - 20).tolist()},
                         "g200_y": {"1900": np.full(NY, 0.3).tolist(), "2060": np.full(NY, 0.2).tolist()},
                         "g200_R": {"1900": 0.3, "2060": 0.2}},
            "series": ser}


def _hyp():
    s_gt = {1900: -0.10, 1980: -0.08, 2020: -0.05, 2060: -0.02}
    s_gc = {y: -0.10 for y in Y}
    terms = {y: {"mld": 1e10, "moc": 3e10, "bio": -1e10, "sum_circ_bio": 3e10, "temp": -2e9,
                 "conductance": 5e4, "parent": {"g200": 0.3 - 1e-4 * (y - 1900)}, "branch": {}} for y in Y}
    return hypotheses.assemble(s_gt, s_gc, {y: 3e10 for y in Y}, {y: 2.8e10 for y in Y}, terms)


def _xc():
    return {"S_fine_members": {"1900": [-0.11, -0.09, -0.12], "2060": [-0.03, -0.01, -0.02]},
            "H6c": {"dS": 0.087, "sigma": 0.01, "detectable": True, "sign_ok": True, "ratio": 1.1,
                    "pass": True, "indicative": False, "n": 3}}


def _data():
    return {"figdata": _figdata(), "hypotheses": json.loads(json.dumps(_hyp())), "crosscheck": _xc()}


@pytest.mark.parametrize("name", list(fg.FIGS))
def test_each_figure_renders(tmp_path, name):
    out = fg.FIGS[name](_data(), tmp_path / name)
    assert all(p.exists() and p.stat().st_size > 1000 for p in out)


def test_main_skips_figures_without_their_json(tmp_path, monkeypatch):
    d = _data()
    for n in ("figdata", "hypotheses"):
        (tmp_path / f"{n}.json").write_text(json.dumps(d[n]))
    monkeypatch.setattr(fg, "FIG", tmp_path)
    made = fg.main()
    assert "fig5_fine" not in made and {"fig0_schematic", "fig3_central"} <= set(made)
    assert (tmp_path / "fig3_central.png").exists() and not (tmp_path / "fig5_fine.png").exists()
