"""Atmospheric CO2 forcing for pkg/dic (dic_int1 = 2, file co2atmos.dat).

pkg/dic reads dic_int2 values (atm, one per line; pkg/dic/dic_init_fixed.F:66-74) and, at
iteration myIter, interpolates linearly between entries (pkg/dic/dic_atmos.F:71-85):
    ntim  = int((myIter - dic_int3) / dic_int4) + 1          (integer division)
    aWght = 0.5 + (myIter - dic_int3)/dic_int4 - (ntim - 1)
    if aWght > 1: ntim = ntim + 1, aWght = aWght - 1
    pCO2  = co2atmos(ntim) (1 - aWght) + co2atmos(ntim + 1) aWght
so entry k carries full weight at myIter = dic_int3 + (k - 1.5) dic_int4. With dic_int3 = the
first iteration of the run (1 January 1850) and one padding entry at each end, entry k = 2 is
1850 at mid-1850, and the run never reads outside the file.
Source of the concentrations: CMIP6 global-mean annual CO2 (Meinshausen et al. 2017, historical;
Meinshausen et al. 2020, SSP3-7.0), file downloaded by Esteban (literatura/data/cmip6_co2)."""
from __future__ import annotations

from pathlib import Path

import numpy as np

YEARS = range(1850, 2101)


def write_co2atmos(co2: dict[int, float], out: Path) -> int:
    """co2: year -> ppm for 1850..2100. Writes atm values with the first and last repeated."""
    v = [co2[y] * 1e-6 for y in YEARS]
    v = [v[0]] + v + [v[-1]]
    Path(out).write_text("".join(f"{x:.10e}\n" for x in v))
    return len(v)


def co2_offset_iter(start_iter: int, steps_per_year: int) -> int:
    """dic_int3 such that the value of year Y applies at mid-Y (module docstring)."""
    return int(start_iter)


def dic_namelist(start_iter: int, steps_per_year: int) -> dict[str, int]:
    return {"dic_int1": 2, "dic_int2": len(YEARS) + 2,
            "dic_int3": co2_offset_iter(start_iter, steps_per_year), "dic_int4": int(steps_per_year)}


def atm_pco2_at(it: int, co2atm, i3: int, i4: int) -> float:
    """Replica of pkg/dic/dic_atmos.F:71-85 (1-based co2atmos -> 0-based array)."""
    a = np.asarray(co2atm, float)
    d = it - i3
    ntim = int(d / i4) + 1 if d >= 0 else -int(-d / i4) + 1  # Fortran INT truncates toward 0
    aw = 0.5 + d / i4 - (ntim - 1)
    if aw > 1.0:
        ntim += 1
        aw -= 1.0
    if ntim < 1 or ntim + 1 > a.size:
        raise IndexError(f"co2atmos({ntim}+1) outside 1..{a.size}")
    return float(a[ntim - 1] * (1.0 - aw) + a[ntim] * aw)


def load_cmip6_co2(path, first: int = 1850, last: int = 2100) -> dict[int, float]:
    """Annual global-mean CO2 (ppm) for first..last from the input4MIPs files in ``path``.

    Reads the yearly ``gr1-GMNHSH`` files of UoM-CMIP-1-2-0 (historical) and
    UoM-AIM-ssp370-1-2-1 (SSP3-7.0); sector 0 is the global mean, units 1e-6 (ppm).
    The year comes from the decoded time midpoint; the historical file wins on overlap."""
    import cftime
    import xarray as xr

    path = Path(path)
    out: dict[int, float] = {}
    for pat in ("*UoM-AIM-ssp370-1-2-1_gr1-GMNHSH_*.nc", "*UoM-CMIP-1-2-0_gr1-GMNHSH_*.nc"):
        files = sorted(f for f in path.glob(pat) if not f.name.startswith("._"))
        if len(files) != 1:
            raise FileNotFoundError(f"expected one file {pat} in {path}, found {len(files)}")
        with xr.open_dataset(files[0], decode_times=False) as d:
            t = d["time"]
            cal = t.attrs.get("calendar", "standard")
            y0 = " 0-1-1" in t.attrs["units"]           # historical file: "days since 0-1-1"
            years = [x.year for x in cftime.num2date(t.values, t.attrs["units"],
                                                     "proleptic_gregorian" if y0 else cal,
                                                     has_year_zero=True if y0 else None)]
            v = d["mole_fraction_of_carbon_dioxide_in_air"].isel(sector=0).values.astype(float)
        if d["mole_fraction_of_carbon_dioxide_in_air"].attrs.get("units") not in ("1.e-6", "1e-6", "ppm"):
            raise ValueError(f"{files[0].name}: unexpected units")
        out.update({int(y): float(x) for y, x in zip(years, v)})
    missing = [y for y in range(first, last + 1) if y not in out]
    if missing:
        raise ValueError(f"CMIP6 CO2 missing years {missing[:5]}...")
    return {y: out[y] for y in range(first, last + 1)}
