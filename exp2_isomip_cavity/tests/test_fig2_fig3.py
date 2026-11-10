"""Tests for fig2_melt_sensitivity.py, fig3_meltwater_export.py and animate_tracer.py.

The figure functions take the ensemble table (ensemble.collect) and, for Fig. 3 and the GIF,
synthetic runs from test_metrics.make_run with known tracer fields, so the plotted data have
analytic values. The figures themselves are checked for existence; what each panel plots
(and with which role: core, test, reference, gap) is returned and checked.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from exp2_isomip_cavity.analysis import animate_tracer as at
from exp2_isomip_cavity.analysis import ensemble as en
from exp2_isomip_cavity.analysis import fig2_melt_sensitivity as f2
from exp2_isomip_cavity.analysis import fig3_meltwater_export as f3
from exp2_isomip_cavity.analysis import metrics as mt
from test_metrics import DZ, NX, NY, NZ, default_geometry, make_run

TF = -2.34
CORE_T = list(en.CORE_TBOT)


def synth_table(gap=None):
    """A table as ensemble.collect returns it: core (m = 2 (T - TF)^1.6), 3 test runs,
    2 geometry runs and the Ocean0 reference."""
    rows = []

    def add(case, group, T, z, s, m, D, status="OK"):
        ok = status == "OK"
        r = {c: None for c in en.COLUMNS}
        r.update(case=case, group=group, Tbot=T, z_c=z, s=s, split=en.split_of(T),
                 status=status, N3_pass=ok, notes=[])
        if ok:
            r.update(mean_melt_shelf=m, mean_melt_deep300=1.5 * m, total_flux=1e4 * m,
                     z_out=-D, D_out=D, N1=0.001, N1_pass=True, N2=1e-6, N2_pass=True,
                     T_f_GL=TF)
        rows.append(r)

    add("ocean0_4km", "reference", 1.0, None, 1.0, 10.0, 300.0)
    for z in en.CORE_ZC:
        for T in CORE_T:
            name = f"ens_T{T}_zc{z:.0f}_s1.0"
            st = "GAP_N3: failed" if name == gap else "OK"
            add(name, "factorial", T, z, 1.0, 2.0 * (T - TF) ** 1.6 * (z / 400.0),
                500.0 - 50 * T - 0.1 * z, status=st)
        add(f"ens_T2.0_zc{z:.0f}_s1.0", "extrapolation", 2.0, z, 1.0,
            2.0 * (2.0 - TF) ** 1.6, 380.0)
    add("ens_T0.34_zc444_s0.8", "geometry", 0.34, 444.0, 0.8, 8.0, 350.0)
    add("ens_T0.97_zc386_s1.2", "geometry", 0.97, 386.0, 1.2, 14.0, 320.0)
    return {k: [r[k] for r in rows] for k in en.COLUMNS}


# ---- Fig. 2 ---------------------------------------------------------------------------------

def test_fig2_writes_png_pdf_json_and_roles(tmp_path):
    tab = synth_table(gap="ens_T0.5_zc250_s1.0")
    res = f2.run(tab, tmp_path,
                 resolution={"shelf_mean_melt": {"rel_diff_4_vs_2": -0.0005},
                             "psi_bt_max_cavity": {"rel_diff_4_vs_2": 1.2458},
                             "melt_map_correlation": {"r": 0.9146}})
    for ext in ("png", "pdf"):
        assert (tmp_path / f"fig2_melt_sensitivity.{ext}").stat().st_size > 0
    js = json.loads((tmp_path / "fig2_h1.json").read_text())
    assert js["H1"]["melt_key"] == "mean_melt_shelf"
    assert js["H1"]["per_z_c"]["250"]["verdict"] == "NOT_EVALUABLE"
    assert js["H1"]["per_z_c"]["400"]["verdict"] == "PASS"
    assert js["H1"]["per_z_c"]["400"]["power_law"]["n"] == pytest.approx(1.6, abs=1e-10)
    cap = js["caption"]
    assert "4 km" in cap and "resolution test" in cap
    # review fix 2: the melt agreement is (largely) the calibration target; cite psi_bt, r
    assert "calibration target" in cap and "+125 %" in cap and "r = 0.91" in cap
    # review fix 1: both halves of H1 and the verdict are visible in panel (a)
    leg = res["plotted"]["a"]["legend"]
    assert any("z_c = 400 m" in l and "n = 1.60" in l and "H1 PASS" in l for l in leg)
    assert any("z_c = 250 m" in l and "H1 not evaluable" in l for l in leg)
    assert "power law" in cap and "[1.3, 2.0]" in cap
    assert js["panel_c_label"] == "descriptive"
    p = res["plotted"]
    # (a): 14 core points (one gap), 3 test points marked as such, the Ocean0 point
    assert sorted(p["a"]["core"]) == sorted(
        c for c, g, st in zip(tab["case"], tab["group"], tab["status"])
        if g == "factorial" and st == "OK")
    assert "ens_T0.5_zc250_s1.0" not in p["a"]["core"]
    assert sorted(p["a"]["test"]) == [f"ens_T2.0_zc{z:.0f}_s1.0" for z in (250, 400, 550)]
    assert p["a"]["reference"] == ["ocean0_4km"]
    assert p["a"]["fit_curves"] == ["400", "550"]          # no fit drawn for the gap z_c
    # (b): heat map over Tbot (core + test column) x z_c, gap cell left empty (NaN)
    grid = np.array(p["b"]["grid"], dtype=float)
    assert grid.shape == (3, 6)
    assert np.isnan(grid[0, 2])                             # z_c 250, Tbot 0.5
    assert np.isfinite(grid[1]).all()
    assert p["b"]["test_column_Tbot"] == 2.0
    # (c): s = 0.8, 1.0 and 1.2 shown, test runs excluded
    assert sorted(p["c"]) == ["0.8", "1.0", "1.2"]
    assert all("T2.0" not in c for c in p["c"]["1.0"])


def test_fig2_default_is_shelf_and_deep_is_secondary(tmp_path):
    """R22: the H1 verdict uses mean_melt_shelf; the deep300 fit is written as secondary."""
    p = en.write_table(synth_table(), tmp_path / "t.json")
    assert f2.main(["--table", str(p), "--out-dir", str(tmp_path)]) == 0
    js = json.loads((tmp_path / "fig2_h1.json").read_text())
    assert js["H1"]["melt_key"] == "mean_melt_shelf" and js["H1"]["role"] == "criterion"
    sec = js["H1_secondary_deep300"]
    assert sec["melt_key"] == "mean_melt_deep300" and sec["role"].startswith("secondary")
    assert all(v["verdict"] == "PASS" for v in sec["per_z_c"].values())
    assert "R22" in js["melt_key_note"]


def test_fig2_legend_shows_n_and_fail(tmp_path):
    """A z_c whose quadratic fits but whose n is outside [1.3, 2.0] is shown as H1 FAIL."""
    tab = synth_table()
    for i, c in enumerate(tab["case"]):
        if tab["group"][i] == "factorial" and tab["z_c"][i] == 550.0:
            tab["mean_melt_shelf"][i] = 0.1 * (tab["Tbot"][i] - TF) ** 2.7
    res = f2.run(tab, tmp_path)
    leg = res["plotted"]["a"]["legend"]
    assert any("z_c = 550 m" in l and "n = 2.70" in l and "H1 FAIL (n > 2.0)" in l
               for l in leg)


# ---- Fig. 3 ---------------------------------------------------------------------------------

def _tracer_run(tmp_path, name="r"):
    """C = 0.01 (k + 1) + 0.001 j in every cell (time-constant); UTRAC01 > 0 at level 2."""
    k, j = np.meshgrid(np.arange(NZ), np.arange(NY), indexing="ij")
    C = (0.01 * (k + 1) + 0.001 * j)[:, :, None] * np.ones((NZ, NY, NX))
    ut = np.zeros((NZ, NY, NX)); ut[2] = 1.0
    return make_run(tmp_path / name, nrec=24, shi=lambda r: {"SHIfwFlx": -1e-5},
                    dyn=lambda r: {"TRAC01": C, "UTRAC01": ut}), C


def test_front_section_reads_x_front_column_and_masks_dry(tmp_path):
    run, C = _tracer_run(tmp_path)
    rd = mt.load_run(run)
    sec = f3.front_section(rd, mt.LAST6)
    assert list(sec["i_front"]) == [4] * NY                 # shelf x < 640 km: columns 0-3
    np.testing.assert_allclose(sec["C"].filled(np.nan), C[:, :, 4], rtol=1e-14)
    _, _, hfac = default_geometry()
    assert (sec["C"].mask == (hfac[:, :, 4] == 0)).all()


def test_column_integral_analytic(tmp_path):
    run, C = _tracer_run(tmp_path)
    rd = mt.load_run(run)
    col = f3.column_tracer(rd, mt.LAST6)
    _, _, hfac = default_geometry()
    expect = (C * hfac * DZ).sum(axis=0)
    np.testing.assert_allclose(col.filled(np.nan), expect, rtol=1e-14)
    # shelf column: draft -40, bed -80 -> levels 2, 3 wet: (0.03 + 0.04 + 2 * 0.001 j) * 20
    assert col[0, 0] == pytest.approx((0.03 + 0.04) * DZ)


def test_fig3_writes_outputs_and_h2(tmp_path):
    run, _ = _tracer_run(tmp_path)
    rd = mt.load_run(run)
    tab = synth_table(gap="ens_T1.0_zc550_s1.0")
    res = f3.run(tab, tmp_path, sections={"cold": ("ens_T-0.5_zc400_s1.0", rd),
                                          "warm": ("ens_T1.5_zc400_s1.0", rd)},
                 map_run=("ens_T1.5_zc400_s1.0", rd))
    for ext in ("png", "pdf"):
        assert (tmp_path / f"fig3_meltwater_export.{ext}").stat().st_size > 0
    js = json.loads((tmp_path / "fig3_h2.json").read_text())
    # synthetic D_out = 500 - 50 T - 0.1 z: strictly decreasing in Tbot -> rho = -1
    assert js["H2"]["per_z_c"]["400"]["rho"] == pytest.approx(-1.0)
    assert js["H2"]["per_z_c"]["400"]["verdict"] == "PASS"
    assert js["H2"]["per_z_c"]["550"]["verdict"] == "NOT_EVALUABLE"
    p = res["plotted"]
    assert sorted(p["b"]["test"]) == [f"ens_T2.0_zc{z:.0f}_s1.0" for z in (250, 400, 550)]
    assert "ens_T1.0_zc550_s1.0" not in p["b"]["core"]
    assert p["a"] == ["ens_T-0.5_zc400_s1.0", "ens_T1.5_zc400_s1.0"]
    assert "psi" not in json.dumps(js).lower()              # psi_bt not used (4 km caveat)


# ---- GIF ------------------------------------------------------------------------------------

def test_gif_within_budget_and_stride_doubles(tmp_path):
    k = np.arange(NZ)[:, None, None]
    make_run(tmp_path / "o", nrec=24, shi=lambda r: {"SHIfwFlx": -1e-5},
             dyn=lambda r: {"TRAC01": 1e-3 * (r + 1) * (k + 1) * np.ones((NZ, NY, NX))})
    out = tmp_path / "o.gif"
    res = at.make_gif(tmp_path / "o", out, dpi=40, max_mb=5.0)
    assert out.exists() and res["mb"] <= 5.0 and not res["over_budget"]
    assert res["frames"] == 24 and res["stride"] == 1
    # fixed colour scale = max over all frames of the column integral (record 24, wet cells)
    _, _, hfac = default_geometry()
    assert res["vmax_column"] == pytest.approx(
        (24e-3 * (k + 1) * hfac * DZ).sum(axis=0).max(), rel=1e-12)
    small = at.make_gif(tmp_path / "o", tmp_path / "s.gif", dpi=40, max_mb=res["mb"] * 0.6)
    assert small["stride"] >= 2 and small["frames"] < 24
