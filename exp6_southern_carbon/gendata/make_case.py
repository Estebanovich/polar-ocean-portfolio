"""Write the MITgcm inputs of one Exp. 6 case (plan Task 2).

write_case(case, out_dir, extra) merges cases.yaml (config + case), writes the binaries (>f4,
read with readBinaryPrec = 32) under fixed names, renders config/input_templates/* with
string.Template (an unknown placeholder raises KeyError) and writes case.json with the actual
parameters. `extra` carries what the registry cannot know:
  co2         year -> ppm (CMIP6; default: gendata.co2.load_cmip6_co2 of literatura/data/cmip6_co2)
  iter_1850   first iteration of GC/GT (end of G0); needed by pco2_mode = transient
  nIter0      start iteration (default 0)
  n_timesteps overrides n_years * steps per year (required when n_years is null: G0, F0, TF)
  pco2_ppm    fixed pCO2 for the smoke/timing runs T10/TF only
  init        {tracer index 1..6 or 'theta': 3-D array} initial fields (fine runs, Task 9)
  cant_sup    year -> C_ant,sup (mol m-3) for sponge = anth (GTa)"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from string import Template

import numpy as np
import yaml

from . import bgc, co2 as co2mod, fields
from .grid import Grid

PKG = Path(__file__).resolve().parents[1]
TEMPLATES = PKG / "config" / "input_templates"
YEAR_S = 31536000.0
WIND_GAS = 10.0          # m s-1, uniform and fixed gas-transfer wind (spec Sec. 3.3)
S_UNIFORM = 34.5


def load_case(case: str) -> dict:
    reg = yaml.safe_load((PKG / "cases.yaml").read_text())
    over = dict(reg["cases"][case])
    p = dict(reg["configs"][over["config"]])
    p.update(over)
    p["case"] = case
    return p


def grid_of(p: dict) -> Grid:
    return Grid(p["nx"], p["ny"], p["dx_m"])


def clock_dt(p: dict) -> float:
    """Model-clock step: deltaT, or deltaT_tracer (= deltaTClock) in the accelerated spin-up."""
    return float(p.get("deltaT_tracer") or p["deltaT"])


def steps_per_year(p: dict) -> int:
    n = YEAR_S / clock_dt(p)
    assert n == int(n), "deltaT must divide the 365-day year"
    return int(n)


def _w(path: Path, a) -> None:
    np.asarray(a, dtype=">f4").tofile(path)


def _co2(p: dict, extra: dict) -> dict[int, float] | None:
    if "co2" in extra:
        return extra["co2"]
    if "pco2_ppm" in extra:
        return None
    return co2mod.load_cmip6_co2(PKG.parent / "literatura" / "data" / "cmip6_co2")


def dic_co2_block(p: dict, extra: dict, out: Path) -> str:
    c = _co2(p, extra)
    mode = p["pco2_mode"]
    if c is None:
        ppm = float(extra["pco2_ppm"])
    elif mode == "transient":
        nml = co2mod.dic_namelist(int(extra["iter_1850"]), steps_per_year(p))
        co2mod.write_co2atmos(c, out / "co2atmos.dat")
        return "\n".join(f"  {k}={v}," for k, v in nml.items())
    elif mode == "fixed1850":
        ppm = c[1850]
    elif mode == "fixedYear":
        ppm = c[int(p["branch_year"])]
    else:
        raise ValueError(mode)
    return f"  dic_int1=1,\n  dic_pCO2={ppm * 1e-6:.10e},"


def template_values(p: dict, g: Grid) -> dict:
    fine = p["config"] == "f10"
    visc = (" viscC2Leith = 1.,\n useFullLeith=.TRUE.,\n viscAhGridMax = 0.5," if fine else " viscAh=2000.,")
    extra_diag = (" fields(1:4,5) = 'UVEL    ','VVEL    ','UVELSQ  ','VVELSQ  ',\n frequency(5) = 2628000.,\n"
                  " filename(5) = 'eke_mon',") if fine else (
                  " fields(1:2,5) = 'GM_PsiX ','GM_PsiY ',\n frequency(5) = 31536000.,\n filename(5) = 'gm_ann',")
    if p.get("deltaT_tracer"):   # Bryan (1984) acceleration: tracers 4x faster than momentum
        dt = (f" deltaTMom={float(p['deltaT']):.1f},\n deltaTtracer={clock_dt(p):.1f},\n"
              f" deltaTClock={clock_dt(p):.1f},\n deltaTFreeSurf={float(p['deltaT']):.1f},")
    else:
        dt = f" deltaT={float(p['deltaT']):.1f},"
    return {"visc_block": visc, "dt_block": dt, "nx": g.nx, "ny": g.ny, "dx": f"{g.dx_m:.1f}",
            "delR": ",\n      ".join(", ".join(f"{d:.8f}" for d in g.dr[i:i + 6]) for i in range(0, g.nr, 6)),
            "useGMRedi": ".FALSE." if fine else ".TRUE.", "ptr_gm": ".FALSE." if fine else ".TRUE.",
            "extra_diags": extra_diag}


def write_case(case: str, out_dir: Path, extra: dict | None = None) -> Path:
    extra = dict(extra or {})
    p = load_case(case)
    g = grid_of(p)
    spy = steps_per_year(p)
    if "n_timesteps" in extra:
        nts = int(extra["n_timesteps"])
    elif p["n_years"] is not None:
        nts = int(p["n_years"]) * spy
    else:
        raise ValueError(f"{case}: n_years is null; pass extra['n_timesteps']")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    f = fields.tutorial_fields(g)
    _w(out / "bathy.bin", f["bathy"])
    _w(out / "zonal_wind.bin", float(p["wind_factor"]) * f["taux"][None] * np.ones((12, 1, 1)))  # periodic: 12 records
    _w(out / "sst_relax.bin", fields.sst_relax_monthly(g, amp=float(extra.get("sst_amp", 1.5))))
    _w(out / "t_relax_mask.bin", f["mask"])
    _w(out / "temperature.bin", f["temp"])
    ones = np.ones((12, g.ny, g.nx))
    _w(out / "par.bin", bgc.par_monthly()[:, None, None] * ones)
    _w(out / "dic_wind.bin", WIND_GAS * ones)
    _w(out / "fice.bin", 0.0 * ones)
    _w(out / "silica.bin", bgc.silica_surface() * ones)
    prof = bgc.ref_profiles(g.z)
    order = ["DIC", "Alk", "PO4", "DOP", "O2", "PO4"]
    init = extra.get("init", {})
    for i, name in enumerate(order, start=1):
        a = init.get(i, prof[name][:, None, None] * np.ones((g.nr, g.ny, g.nx)))
        _w(out / f"init_ptr{i}.bin", a)
    rbcs_period = ""
    for i, name in enumerate(order[:3], start=1):
        base = prof[name][:, None, None] * np.ones((g.nr, g.ny, g.nx))
        if i == 1 and p["sponge"] == "anth":
            rec = bgc.anth_sponge(prof["DIC"], extra["cant_sup"], g.z)          # (251, nr)
            _w(out / "relax_ptr1.bin", rec[:, :, None, None] * np.ones((1, 1, g.ny, g.nx)))
            rbcs_period = " rbcsForcingPeriod=31536000.,\n rbcsForcingCycle=0.,"
        else:
            _w(out / f"relax_ptr{i}.bin", base)
    vals = template_values(p, g)
    vals.update({"nIter0": int(extra.get("nIter0", 0)), "nTimeSteps": nts,
                 "dic_co2": dic_co2_block(p, extra, out), "rbcs_period": rbcs_period})
    for t in sorted(TEMPLATES.iterdir()):
        (out / t.name).write_text(Template(t.read_text()).substitute(vals))
    meta = dict(p, n_timesteps=nts, nIter0=int(extra.get("nIter0", 0)), steps_per_year=spy,
                iter_1850=extra.get("iter_1850"), sst_amp=float(extra.get("sst_amp", 1.5)),
                wind_gas_ms=WIND_GAS, salinity=S_UNIFORM, silica_mol_m3=bgc.silica_surface())
    (out / "case.json").write_text(json.dumps(meta, indent=1) + "\n")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("out_dir")
    ap.add_argument("--n-timesteps", type=int)
    ap.add_argument("--pco2-ppm", type=float)
    a = ap.parse_args(argv)
    extra = {}
    if a.n_timesteps:
        extra["n_timesteps"] = a.n_timesteps
    if a.pco2_ppm:
        extra["pco2_ppm"] = a.pco2_ppm
    print(write_case(a.case, Path(a.out_dir), extra))


if __name__ == "__main__":
    main()
