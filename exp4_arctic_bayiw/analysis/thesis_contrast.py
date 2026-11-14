"""Task 5, Step 5: contrast of the analysis library with the MSc-thesis scripts on R1.

Two contrasts, both on the realistic bay run (R1 is bit-identical to the thesis run
``BayIW_Rectan/run_expand``, figures/v1_thesis.json):

* **Nash flux** (Nash et al. 2005): ``flux.section_flux`` on the R1 MDS output against the
  thesis script ``thesis_copy/Output_IWBay/energy_flux_IWBay.py`` run on the thesis MNC output,
  section j = 323, i = 274..285 (ruling R16), per stage (forced 0-60 h, free 60-120 h).
  Output: ``figures/flux_vs_thesis.json``.
* **Kang & Fringer (2012) budget** of the bay box (j 324-347, i 274-285, 0-120 h):
  ``budget.kf_budget`` on R1 against ``thesis_copy/Output_IWBay/balance_kang_fringer.py``
  (``calcular``) run on the thesis MNC output, term by term; the residual D as in the thesis
  ``figura_energia_balance_kf.terminos``. Also the explicit Laplacian dissipation
  (``disipacion_viscosa_bahia.py``). Output: ``figures/budget_vs_thesis.json``.

How the thesis scripts are run (stated in the JSON): their own functions are imported from the
verbatim copies in ``thesis_copy/Output_IWBay`` and called on the thesis MNC files; nothing is
re-implemented on the reference side. The pyocean environment has no ``MITgcmutils``, which
``MIT_DatProcesW.py`` imports at module level but none of the called functions use; it is
replaced by an empty stub module. ``figura_energia_balance_kf.py`` imports the thesis plotting
module ``estructura_vertical`` only for the npz folder; it is stubbed with that folder. The
thesis ``main()`` functions are not called (they would write into the read-only thesis tree);
the reference results are cached as npz in ``EXP4_RUNS/analysis_cache`` (or ``--cache``). Every cache
file is keyed on its kind, the run (or thesis scenario) and a hash of the source files that make it
(``CACHE_SOURCES``), so a different run or changed code misses the cache; ``--recompute`` forces
the library series to be recomputed.

Sign conventions: project + = into the bay (and into the fjord). energy_flux_IWBay.py uses
+ = into the bay (+y, "Positivo = flujo neto entrando a la bahía"); balance_kang_fringer.py uses
+ = out of the bay ("+ = sale de la bahía"), so its F and A are compared with a minus sign.

Usage::

    thesis_contrast.py reference-nash | reference-kf | reference-diss   # thesis side (cache)
    thesis_contrast.py contrast [--run R1]                              # both JSON files
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import sys
import types
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parents[1]
REPO = PKG.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pfcommon.paths import EXP4_RUNS  # noqa: E402

THESIS_SCRIPTS = PKG / "thesis_copy" / "Output_IWBay"
THESIS_RUN = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand")
CACHE = EXP4_RUNS / "analysis_cache"
SCEN = "realistic"                          # the thesis scenario that R1 reproduces
J_SEC, I_SEC = 323, slice(274, 286)         # ruling R16 (energy_flux_IWBay.J_BOCA, X_BOCA)
J_BAY, I_BAY = (324, 348), (274, 286)       # balance_kang_fringer.J_BAHIA, I_BAHIA
T_STAGE2 = 216000.0                         # s; thesis ITER_ETAPA2 = 7200 x 30 s
TOL_REL = 0.01                              # Task 5 Step 5: < 1 %
ANALYSIS = PKG / "analysis"
_TS = THESIS_SCRIPTS
# Sources whose content keys each cache file (a change in any of them misses the cache).
CACHE_SOURCES = {
    "thesis_nash": [_TS / "energy_flux_IWBay.py"],
    "thesis_kf": [_TS / f for f in ("balance_kang_fringer.py", "conversion_bt_2D.py",
                                    "fondo_kang_fringer.py", "energia_viento.py",
                                    "ape_dos_capas.py", "conversion_baroclinica_2D.py",
                                    "MIT_DatProcesW.py", "energy_flux_IWBay.py")],
    "thesis_diss": [_TS / f for f in ("disipacion_viscosa_bahia.py", "balance_kang_fringer.py",
                                      "energy_flux_IWBay.py", "MIT_DatProcesW.py")],
    "library_kf": [ANALYSIS / "budget.py", ANALYSIS / "load.py"],
}


def source_hash(paths) -> str:
    """First 12 hex digits of the sha256 of the concatenated file contents."""
    h = hashlib.sha256()
    for p in paths:
        h.update(Path(p).read_bytes())
    return h.hexdigest()[:12]


def cache_file(cache: Path, kind: str, run_name: str, sources) -> Path:
    """Cache file keyed on the kind, the run (or thesis scenario) and the source hash."""
    return Path(cache) / f"{kind}_{run_name}_{source_hash(sources)}.npz"


def cached(f: Path, compute, recompute: bool = False) -> dict:
    """np.load(f) if it exists and not ``recompute``; else compute() -> dict, saved to f."""
    f = Path(f)
    if f.exists() and not recompute:
        return dict(np.load(f))
    r = compute()
    f.parent.mkdir(parents=True, exist_ok=True)
    np.savez(f, **r)
    return dict(r)


def _ref_file(cache, kind):
    return cache_file(cache, kind, f"{SCEN}_bay", CACHE_SOURCES[kind])


# ---- thesis side -------------------------------------------------------------------------------

def thesis_modules():
    """Import the verbatim thesis copies (stubbing MITgcmutils / estructura_vertical if absent)."""
    sys.dont_write_bytecode = True
    p = str(THESIS_SCRIPTS)
    if p not in sys.path:
        sys.path.insert(0, p)
    try:
        importlib.import_module("MITgcmutils")
    except ImportError:
        sys.modules["MITgcmutils"] = types.ModuleType("MITgcmutils")
    if "estructura_vertical" not in sys.modules:
        try:
            importlib.import_module("estructura_vertical")
        except ImportError:
            ev = types.ModuleType("estructura_vertical")
            ev.BASE = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay")
            sys.modules["estructura_vertical"] = ev
    names = ("energy_flux_IWBay", "balance_kang_fringer", "disipacion_viscosa_bahia",
             "figura_energia_balance_kf")
    return {n: importlib.import_module(n) for n in names}


def reference_nash(cache: Path = CACHE) -> Path:
    """energy_flux_IWBay.py on the thesis run, realistic scenario, the steps of its main()."""
    m = thesis_modules()["energy_flux_IWBay"]
    cfg = m.SCENARIOS[SCEN]
    run_dir = cfg["base"] / "run_expand"
    eos = {k: cfg[k] for k in ("tRef", "sRef", "tAlpha", "sBeta", "rhoNil")}
    groups = m.discover_stage_groups(run_dir)
    todos = [g for grp in groups for g in grp]
    dx_sec = m._leer_dx_seccion(run_dir)
    nz_sec, _, H_sec = m.niveles_humedos_seccion(run_dir)
    rho_ref, u_mean, v_mean, drF = m.calcular_medias_run(run_dir, todos, nz_sec=nz_sec, **eos)
    out = dict(rho_ref=rho_ref, u_mean_sec=u_mean, v_mean_sec=v_mean, drF=drF, dx_sec=dx_sec,
               H_sec=H_sec, j_boca=m.J_BOCA, x_boca=[m.X_BOCA.start, m.X_BOCA.stop])
    for n, grp in enumerate(groups, start=1):
        r = m.procesar_etapa(run_dir, grp, rho_ref, u_mean, v_mean, drF, **eos)
        for k in ("t_seconds", "Fy_sec", "Fy_bt_sec"):
            out[f"{k}_e{n}"] = r[k]
    f = _ref_file(cache, "thesis_nash")
    cached(f, lambda: out, recompute=True)
    return f


def reference_kf(cache: Path = CACHE) -> Path:
    """balance_kang_fringer.calcular('realistic') on the thesis run (time series per record)."""
    m = thesis_modules()["balance_kang_fringer"]
    res, meta = m.calcular(SCEN)
    f = _ref_file(cache, "thesis_kf")
    cached(f, lambda: dict(res, meta_json=json.dumps(meta)), recompute=True)
    return f


def reference_dissipation(cache: Path = CACHE) -> Path:
    """disipacion_viscosa_bahia.calcular('realistic') on the thesis run."""
    m = thesis_modules()["disipacion_viscosa_bahia"]
    res = m.calcular(SCEN)
    f = _ref_file(cache, "thesis_diss")
    cached(f, lambda: res, recompute=True)
    return f


# ---- library side and contrast -------------------------------------------------------------------

def _rel(a, b):
    """|a - b| / |b| (b the thesis value); 0 if both are 0."""
    a, b = float(a), float(b)
    return 0.0 if a == b else abs(a - b) / abs(b) if b != 0 else float("inf")


def _series_rel(a, b):
    """max_t |a - b| / max_t |b|."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b))) if np.max(np.abs(b)) > 0 else 0.0


def contrast_nash(ds, ref: dict, stored: Path | None = None) -> dict:
    from exp4_arctic_bayiw.analysis import flux
    sec = flux.section_from_dataset(ds, "y", J_SEC, I_SEC)       # +y = into the bay (both)
    t = sec["t"]
    r = flux.section_flux(sec["u_n"], sec["rho"], sec["eta"], sec["dz"], sec["dl"], t,
                          rho0=ds.attrs["rho0"], g=ds.attrs["g"],
                          stages=[(0.0, T_STAGE2), (T_STAGE2, t[-1])])
    out = {"section": {"j": J_SEC, "i": [I_SEC.start, I_SEC.stop - 1]},
           "n_records": int(t.size), "records_thesis": int(len(np.unique(np.concatenate(
               [ref["t_seconds_e1"], ref["t_seconds_e2"]])))),
           "dx_max_rel_diff": float(np.max(np.abs(sec["dl"] - ref["dx_sec"]) / ref["dx_sec"])),
           "rho_ref_max_abs_diff_kg_m-3": float(np.max(np.abs(r["rho_ref"] - ref["rho_ref"]))),
           "u_mean_max_abs_diff_m_s-1": float(np.max(np.abs(r["u_mean"] - ref["v_mean_sec"]))),
           "stages": {}}
    worst = 0.0
    for n, s in enumerate(r["stages"], start=1):
        tr = ref[f"t_seconds_e{n}"]
        if not np.array_equal(s["t"], tr):
            raise RuntimeError(f"stage {n}: record times differ from the thesis")
        Fbc_th, Fbt_th = ref[f"Fy_sec_e{n}"], ref[f"Fy_bt_sec_e{n}"]
        dx = ref["dx_sec"]
        e = {"t_h": [float(tr[0] / 3600), float(tr[-1] / 3600)], "n_records": int(tr.size)}
        for k, mine, th in (("F_bc", s["F_bc"], Fbc_th), ("F_bt", s["F_bt"], Fbt_th)):
            tot_m, tot_t = mine @ sec["dl"], th @ dx
            e[k] = {"total_mean_W": [float(tot_m.mean()), float(tot_t.mean())],
                    "total_mean_rel_diff": _rel(tot_m.mean(), tot_t.mean()),
                    "column_mean_max_rel_diff": float(np.max(np.abs(mine.mean(0) - th.mean(0)))
                                                      / np.max(np.abs(th.mean(0)))),
                    "series_total_max_rel_diff": _series_rel(tot_m, tot_t)}
            worst = max(worst, e[k]["total_mean_rel_diff"], e[k]["series_total_max_rel_diff"])
        out["stages"][f"stage{n}"] = e
    if stored is not None and stored.exists():
        d = np.load(stored)
        out["stored_thesis_npz"] = {"file": str(stored), "stage_total_mean_rel_diff": {
            f"stage{n}": {k: _rel((r["stages"][n - 1][k] @ sec["dl"]).mean(),
                                  (d[f"{kk}_e{n}"] @ d["dx_sec"]).mean())
                          for k, kk in (("F_bc", "Fy_sec"), ("F_bt", "Fy_bt_sec"))}
            for n in (1, 2)}}
    out["max_rel_diff"] = worst
    out["pass_lt_1pct"] = bool(worst < TOL_REL)
    return out


def bay_mask(ds, rows=J_BAY, cols=I_BAY):
    ny, nx = ds["hFacC"].shape[1:]
    m = np.zeros((ny, nx), bool)
    m[rows[0]:rows[1], cols[0]:cols[1]] = True
    return m


TERMS = ("ke_bt", "pe_bt", "ke_bc", "ape", "ape_linear", "F_bt", "F_bc", "A_bt", "A_bc", "C",
         "P_bt", "P_bc")
FLUX_TERMS = ("F_bt", "F_bc", "A_bt", "A_bc")                 # thesis: + out of the bay
THESIS_KEY = {"ape_linear": "ape_lineal"}


def contrast_kf(ds, ref: dict, figs, stored: Path | None = None, cache: Path = CACHE,
                recompute: bool = False) -> dict:
    from exp4_arctic_bayiw.analysis import budget, load
    tau = (lambda t: load.forcing_stress(ds, t))
    def compute():
        r = budget.kf_budget(ds, bay_mask(ds), tau=tau, conversion_masks={
            "escalon": bay_mask(ds, rows=(J_BAY[0] - 1, J_BAY[0]))})
        return {k: v for k, v in r.items() if isinstance(v, np.ndarray)}

    f = cache_file(cache, "library_kf", Path(ds.attrs["run_dir"]).name, CACHE_SOURCES["library_kf"])
    r = cached(f, compute, recompute)        # the library series (reading R1 is slow)
    esc = {"t": r["t"], "C": r["C_escalon"]}
    if not np.allclose(r["t"] / 3600.0, ref["t_horas"], atol=1e-9, rtol=0):
        raise RuntimeError("record times differ from the thesis")
    th = {k: (-ref[THESIS_KEY.get(k, k)] if k in FLUX_TERMS else ref[THESIS_KEY.get(k, k)])
          for k in TERMS}                                     # thesis -> + into the bay
    out = {"box": {"j": [J_BAY[0], J_BAY[1] - 1], "i": [I_BAY[0], I_BAY[1] - 1]},
           "n_records": int(r["t"].size), "terms": {}}
    worst = 0.0
    c_me = budget.closure(r)
    c_th = figs.terminos({k: v for k, v in ref.items()}, ape_campo="ape")
    for k in TERMS:
        e = {"series_max_rel_diff": _series_rel(r[k], th[k])}
        if k in ("ke_bt", "pe_bt", "ke_bc", "ape", "ape_linear"):
            e["value_120h_J"] = [float(r[k][-1]), float(th[k][-1])]
            e["rel_diff_120h"] = _rel(r[k][-1], th[k][-1])
        else:
            a, b = budget.cumulative(r["t"], r[k]), budget.cumulative(ref["t_horas"] * 3600, th[k])
            e["integral_0_120h_J"] = [float(a[-1]), float(b[-1])]
            e["rel_diff_integral"] = _rel(a[-1], b[-1])
        out["terms"][k] = e
        if k != "ape_linear":
            worst = max(worst, e["series_max_rel_diff"])
    for comp in ("bt", "bc"):
        for k in ("dE", "F", "D"):
            mine = c_me[comp][k][-1]
            # thesis F is + out; D is the same quantity
            thv = -c_th[comp][k][-1] if k == "F" else c_th[comp][k][-1]
            out["terms"][f"{k}_{comp}_cumulative"] = {"value_120h_J": [float(mine), float(thv)],
                                                      "rel_diff": _rel(mine, thv)}
            worst = max(worst, _rel(mine, thv))
    a, b = budget.cumulative(esc["t"], esc["C"]), budget.cumulative(ref["t_horas"] * 3600,
                                                                    ref["C_escalon"])
    out["terms"]["C_escalon"] = {"series_max_rel_diff": _series_rel(esc["C"], ref["C_escalon"]),
                                 "integral_0_120h_J": [float(a[-1]), float(b[-1])],
                                 "rel_diff_integral": _rel(a[-1], b[-1])}
    out["rho_b_max_abs_diff_kg_m-3"] = float(np.max(np.abs(r["rho_b"] - np.array(
        json.loads(str(ref["meta_json"]))["rho_b"]))))
    if stored is not None and stored.exists():
        d = np.load(stored)
        out["stored_thesis_npz"] = {"file": str(stored), "series_max_rel_diff_vs_rerun": {
            k: _series_rel(d[THESIS_KEY.get(k, k)], ref[THESIS_KEY.get(k, k)]) for k in TERMS}}
    out["max_rel_diff"] = worst
    out["pass_lt_1pct"] = bool(worst < TOL_REL)
    return out, r, c_me


def contrast_dissipation(ds, ref: dict, c_me: dict) -> dict:
    """Explicit Laplacian dissipation of the bay (descriptive, R17) against the thesis
    disipacion_viscosa_bahia.py; and its share of the residual D."""
    from exp4_arctic_bayiw.analysis import budget
    th = thesis_modules()["disipacion_viscosa_bahia"]
    r = budget.explicit_dissipation(ds, J_BAY, I_BAY, th.NU_H, th.NU_V)
    out = {"nu_h_m2_s-1": th.NU_H, "nu_v_m2_s-1": th.NU_V, "terms": {}}
    worst = 0.0
    for k_me, k_th in (("bt", "bt"), ("bc_h", "bc_h"), ("bc_v", "bc_v"), ("sup", "sup")):
        a = budget.cumulative(r["t"], r[k_me])[-1]
        b = budget.cumulative(ref["t"], ref[k_th])[-1]
        e = {"series_max_rel_diff": _series_rel(r[k_me], ref[k_th]),
             "integral_0_120h_J": [float(a), float(b)], "rel_diff_integral": _rel(a, b)}
        out["terms"][k_me] = e
        worst = max(worst, e["series_max_rel_diff"])
    eps_bt = budget.cumulative(r["t"], r["bt"])[-1]
    eps_bc = budget.cumulative(r["t"], r["bc_h"] + r["bc_v"])[-1]
    out["share_of_residual"] = {"eps_bt_over_D_bt": float(eps_bt / c_me["bt"]["D"][-1]),
                                "eps_bc_over_D_bc": float(eps_bc / c_me["bc"]["D"][-1])}
    out["max_rel_diff"] = worst
    return out


def _load_ref(cache, kind):
    f = _ref_file(cache, kind)
    if not f.exists():
        raise FileNotFoundError(f"{f}: run `thesis_contrast.py reference-{kind[7:]}` first")
    return dict(np.load(f))


def contrast(run: Path, cache: Path = CACHE, recompute: bool = False) -> dict:
    from exp4_arctic_bayiw.analysis import load
    mods = thesis_modules()
    ds = load.run(run)
    base = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay")
    common = {"run": str(run), "thesis_run": str(THESIS_RUN), "final_record": ds.attrs[
        "final_record_source"], "how_the_thesis_ran": (
        "thesis functions imported from the verbatim copies in thesis_copy/Output_IWBay and "
        "called on the thesis MNC output (MITgcmutils stubbed: not installed in pyocean and not "
        "used by the called functions); main() not called; reference cached in " + str(cache)),
        "tolerance_rel": TOL_REL}
    nash = dict(common, script="thesis_copy/Output_IWBay/energy_flux_IWBay.py",
                sign="+ into the bay in both (thesis +y)",
                **contrast_nash(ds, _load_ref(cache, "thesis_nash"),
                                base / "npz_flux" / "realistic_bay_flujo.npz"))
    (PKG / "figures" / "flux_vs_thesis.json").write_text(json.dumps(nash, indent=1))
    kf_ref = _load_ref(cache, "thesis_kf")
    kfres, r, c = contrast_kf(ds, kf_ref, mods["figura_energia_balance_kf"],
                              base / "npz_energy" / "realistic_balance_kf_bahia.npz",
                              cache=cache, recompute=recompute)
    kf = dict(common, script="thesis_copy/Output_IWBay/balance_kang_fringer.py",
              sign="project F, A + into the bay = - thesis (+ out of the bay); D identical",
              **kfres)
    dref = _ref_file(cache, "thesis_diss")
    if dref.exists():
        kf["explicit_dissipation_descriptive"] = contrast_dissipation(ds, dict(np.load(dref)), c)
    (PKG / "figures" / "budget_vs_thesis.json").write_text(json.dumps(kf, indent=1))
    return {"summary": {"nash_max_rel_diff": nash["max_rel_diff"], "nash_pass": nash["pass_lt_1pct"],
                        "kf_max_rel_diff": kf["max_rel_diff"], "kf_pass": kf["pass_lt_1pct"]}}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("what", choices=("reference-nash", "reference-kf", "reference-diss", "contrast"))
    ap.add_argument("--run", type=Path, default=EXP4_RUNS / "R1")
    ap.add_argument("--cache", type=Path, default=CACHE)
    ap.add_argument("--recompute", action="store_true",
                    help="ignore the cached library series of this run and recompute them")
    a = ap.parse_args(argv)
    if a.what == "reference-nash":
        print(reference_nash(a.cache))
    elif a.what == "reference-kf":
        print(reference_kf(a.cache))
    elif a.what == "reference-diss":
        print(reference_dissipation(a.cache))
    else:
        res = contrast(a.run, a.cache, a.recompute)
        print(json.dumps(res["summary"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
