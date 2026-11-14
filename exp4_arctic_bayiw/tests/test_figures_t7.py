"""Tests for exp4_arctic_bayiw.analysis.figures (Task 7, Step 3): Figures 2-4 render from
synthetic partition/budget JSONs and synthetic maps, and write their numbers to a JSON of the
same stem."""
from __future__ import annotations

import json

import numpy as np

from exp4_arctic_bayiw.analysis import figures as fg

B5 = ("low", "kelvin_subinertial", "near_inertial", "superinertial", "seiche")


def _bands(v):
    return dict(zip(B5, v))


def _spec(peak):
    f = np.linspace(0, 1.5, 200)
    return {"freq_cph": f.tolist(), "psd_m2_s-2_cph-1": (1e-3 + np.exp(-((f - peak) / 0.01) ** 2)).tolist()}


def _run(peak, st_names, straits=True):
    st = {n: {"fraction": _bands([0.1, 0.3, 0.4, 0.1, 0.1]), "band_m2_s-2": _bands([1, 3, 4, 1, 1])}
          for n in st_names}
    r = {"bands_cph": {b: [0.01 * k, 0.01 * (k + 1)] for k, b in enumerate(B5)},
         "stations": st, "spectra": {n: _spec(peak) for n in st_names}}
    if straits:
        r["straits"] = {"total": {"FA_bc_out_J": _bands([-1e9, 2e8, -1e8, 5e7, 1e7]),
                                  "F_bc_out_J": _bands([-1e9, 1e8, -1e8, 5e7, 1e7])}}
    return r


def _partition():
    runs = {"F1": _run(0.09, ["shelf", "coast_south_mean"]), "F3": _run(0.085, ["shelf", "coast_south_mean"]),
            "F2": _run(0.085, ["shelf", "coast_south_mean"]), "F1c": _run(0.09, ["shelf"]),
            "R1": _run(0.04, ["Seiche", "Omega"]), "R1k": _run(0.04, ["Seiche", "Omega"]),
            "R2": _run(0.085, ["Seiche", "Omega"]), "R3": _run(0.085, ["Omega"], straits=False)}
    runs["R3"]["straits"] = None
    return {"runs": runs}


def _budget():
    t = list(np.arange(0, 280, 1.0))
    def b(scale):
        return {"series": {"t_h": t,
                           "FA_bc_in_W": {"Akselsundet": list(scale * np.sin(np.array(t) / 10)),
                                          "Mariasundet": list(0.3 * scale * np.cos(np.array(t) / 10))},
                           "F_bc_in_W": {"Akselsundet": list(scale * np.sin(np.array(t) / 12)),
                                         "Mariasundet": list(0.2 * scale * np.cos(np.array(t) / 12))},
                           "bc_cumulative_J": {k: list(scale * 1e3 * np.array(t) * (q + 1))
                                               for q, k in enumerate(("P", "C", "F", "dE", "D"))}},
                "closure_J": {"bc": {"P": 1.0, "C": 1.0, "F": 1.0, "dE": 1.0, "D": 1.0}}}
    return {"runs": {"F1": b(1e6), "F2": b(1e5), "F3": b(5e5), "F1c": b(8e5)}}


def _maps():
    jj, ii = np.meshgrid(np.arange(140, 160), np.arange(200, 260), indexing="ij")
    jj, ii = jj.ravel(), ii.ravel()
    rng = np.random.default_rng(0)
    m = {"band_names": list(B5), "x0_m": 455000.0, "y0_m": 8594000.0, "dx_m": 200.0}
    for r, s in (("F2", 1.0), ("F3", 10.0)):
        m[r] = {"jj": jj, "ii": ii, "ke_bc_band": s * rng.random((jj.size, 5)) + 1e-3}
    return m


def test_fig2_renders(tmp_path):
    part = _partition()
    part["runs"]["F1"]["straits"]["total"]["FA_bc_out_J"] = {}      # Enmienda 3: F_bc is plotted
    rec = fg.fig2_spectra_partition(part, tmp_path)
    assert (tmp_path / "fig2_spectra_partition.png").exists() and (tmp_path / "fig2_spectra_partition.pdf").exists()
    saved = json.loads((tmp_path / "fig2_spectra_partition.json").read_text())
    assert saved == json.loads(json.dumps(rec)) and "F1" in saved["flux_share"]


def test_fig3_renders(tmp_path):
    bud = _budget()
    for r in bud["runs"].values():
        r["series"].pop("FA_bc_in_W")                                # Enmienda 3: F_bc is plotted
    rec = fg.fig3_strait_flux_budget(bud, tmp_path)
    assert (tmp_path / "fig3_strait_flux_budget.png").exists()
    assert set(rec["final_cumulative_bc_J"]) >= {"F1", "F2"}


def test_fig4_renders(tmp_path):
    rec = fg.fig4_ice_maps(_maps(), tmp_path)
    assert (tmp_path / "fig4_ice_maps.png").exists()
    assert rec["box_total_ke_bc_J"]["F3"]["near_inertial"] > rec["box_total_ke_bc_J"]["F2"]["near_inertial"]


def test_fig2_takes_f_from_theory_json(tmp_path):
    from exp4_arctic_bayiw.analysis import validate as va
    th = va.load_theory()
    rec = fg.fig2_spectra_partition(_partition(), tmp_path)
    assert rec["f_cph"]["fjord"] == th["bands"]["fjord_78N_280h"]["f_cph"]
    assert rec["f_cph"]["thesis"] == th["bands"]["bay_thesis_f_240h"]["f_cph"]
    assert not hasattr(fg, "F_CPH")
