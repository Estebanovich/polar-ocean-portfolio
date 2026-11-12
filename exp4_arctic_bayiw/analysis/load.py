"""Lazy reader of an Exp. 4 MITgcm run (MDS diagnostics) as an xarray.Dataset.

``run(run_dir)`` checks every stage with ``pfcommon.runctl.require_ok`` (the run directory and
its ``stage<N>`` subdirectories, as written by scripts/run_campaign.py) and returns

* grid (eager): drF, drC, RC (k); hFacC, hFacW, hFacS (k, j, i); rA, rAs, rAw, rAz, dxG, dyG,
  dxF, dyF, dxC, dyC, dxV, dyU, XC, YC, Depth (j, i). MDS face indexing: U[..., j, i] and hFacW
  are on the west face of cell (j, i), V and hFacS on its south face.
* time series (lazy, memory-mapped; only the slices asked for are read): U, V [m s^-1],
  RHO [kg m^-3] = rhoConst + RHOAnoma (model/src/find_rho.F:101-108, RHOAnoma = rho - rhoConst,
  pkg/diagnostics/diagnostics_main_init.F:500-501; rhoConst = rhoNil if unset,
  model/src/ini_parms.F:483), ETA [m], and T, TAUX, TAUY [N m^-2] if present (oceTAUX at U
  points, oceTAUY at V points). time [s] = iteration x deltaT.
* attrs: rho0 (= rhoConst, the Boussinesq reference), g, eos (eosType, tAlpha, sBeta, tRef,
  sRef, rhoNil), deltaT, stages (dir, nIter0, nTimeSteps, wind files), final_record_source.

Final record. A diagnostics snapshot labelled iteration n is filled at the start of step n + 1
(DO_STATEVARS_DIAGS, model/src/forward_step.F:509-510), so the state after the last step is
never written as a diagnostic. If the last stage holds a pickup at its last iteration, one record
spacing after the last diagnostic, that state is appended (U, V, T, ETA from Uvel, Vvel, Theta,
EtaN; RHO from the linear EOS of Theta and Salt). This gives the same 0-120 h records as the
thesis MNC state files for R1 (481). Disable with ``final_from_pickup=False``. Variables the
pickup does not hold (TAUX, TAUY) disable the append.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from xarray.backends import BackendArray
from xarray.core import indexing

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pfcommon.mitgcm_io import _parse_meta  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402

GRID2D = {"XC": "XC", "YC": "YC", "RAC": "rA", "RAS": "rAs", "RAW": "rAw", "RAZ": "rAz",
          "DXG": "dxG", "DYG": "dyG", "DXF": "dxF", "DYF": "dyF", "DXC": "dxC", "DYC": "dyC",
          "DXV": "dxV", "DYU": "dyU", "Depth": "Depth"}
GRID3D = {"hFacC": "hFacC", "hFacW": "hFacW", "hFacS": "hFacS"}
GRID1D = {"DRF": "drF", "DRC": "drC", "RC": "RC"}
DIAG = {"UVEL": "U", "VVEL": "V", "THETA": "T", "RHOAnoma": "RHO", "ETAN": "ETA",
        "oceTAUX": "TAUX", "oceTAUY": "TAUY"}
PICKUP = {"U": "Uvel", "V": "Vvel", "T": "Theta", "ETA": "EtaN"}


def _grid(stem: Path) -> np.ndarray:
    """A grid file without iteration number (``DRF.meta``/``DRF.data``), dims as in its meta."""
    meta = _parse_meta(Path(f"{stem}.meta").read_text())
    return np.fromfile(Path(f"{stem}.data"), ">f8").reshape(meta["shape"]).astype(np.float64)


def stage_dirs(run_dir: Path) -> list[Path]:
    run_dir = Path(run_dir)
    return [run_dir] + sorted((p for p in run_dir.glob("stage[0-9]*") if p.is_dir()),
                              key=lambda p: int(p.name[5:]))


def _nml(d: Path) -> dict:
    from exp4_arctic_bayiw.gendata.make_case import parse_namelist
    return parse_namelist((Path(d) / "data").read_text())


def _num(v: str | None, default=None):
    if v is None:
        return default
    return float(v.replace("D", "E").replace("d", "e").rstrip(","))


def _list(v: str | None, nz: int):
    """MITgcm namelist list with repeat counts (``40*0.4676``, ``16.582, 16.582``)."""
    if v is None:
        return None
    out = []
    for tok in re.findall(r"(?:(\d+)\*)?([-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?)", v):
        out += [float(tok[1].replace("D", "E").replace("d", "e"))] * (int(tok[0]) if tok[0] else 1)
    if len(out) == 1:
        out *= nz
    return np.array(out[:nz])


def _str(v: str | None):
    return None if v is None else v.strip().strip("'\"")


class _Records(BackendArray):
    """(time, [k,] j, i) view of one field over many MDS files, read on indexing."""

    def __init__(self, items, shape3):
        self.items = items                    # [(data path, nrec, rec, shape) | ("array", a)]
        self.shape = (len(items),) + tuple(shape3)
        self.dtype = np.dtype("float64")

    def __getitem__(self, key):
        return indexing.explicit_indexing_adapter(key, self.shape, indexing.IndexingSupport.OUTER,
                                                  self._get)

    def _record(self, n):
        it = self.items[n]
        if it[0] == "array":
            return it[1]
        path, nrec, rec, shape, offset = it
        m = np.memmap(path, dtype=">f8", mode="r", shape=(nrec,) + shape)
        a = m[rec]
        return a + offset if offset else a

    def _get(self, key):
        tkey, rest = key[0], key[1:]
        tidx = np.arange(self.shape[0])[tkey]
        scalar_t = np.ndim(tidx) == 0
        out = []
        for n in np.atleast_1d(tidx):
            a = self._record(int(n))
            basic = tuple(k if isinstance(k, (slice, int, np.integer)) else slice(None) for k in rest)
            a = a[basic]
            ax = 0
            for k in rest:                    # outer (array) keys after the basic ones
                if isinstance(k, (int, np.integer)):
                    continue
                if not isinstance(k, slice):
                    a = np.take(a, k, axis=ax)
                ax += 1
            out.append(np.array(a, dtype=np.float64))
        return out[0] if scalar_t else np.stack(out)


def _diag_files(dirs):
    """{field: {iteration: (data path, nrec, rec, shape)}} of every MDS diagnostics file."""
    out: dict = {}
    for d in dirs:
        for m in sorted(d.glob("*.meta")):
            parts = m.name.split(".")
            if len(parts) != 3 or not re.fullmatch(r"\d{10}", parts[1]):
                continue
            meta = _parse_meta(m.read_text())
            flds = meta["fldList"]
            if not flds or meta["timeStepNumber"] is None:
                continue
            shape = meta["shape"]
            nper = 1
            for k, f in enumerate(flds):
                if f not in DIAG:
                    continue
                out.setdefault(DIAG[f], {}).setdefault(
                    int(parts[1]), (m.with_suffix(".data"), meta["nrecords"], k * nper, shape))
    return out


def _pickup(d: Path, it: int):
    for m in sorted(d.glob("pickup*.meta")):
        if not re.fullmatch(r"pickup\.(ckpt[AB]|\d{10})\.meta", m.name):
            continue
        meta = _parse_meta(m.read_text())
        if meta["timeStepNumber"] == it:
            return m.with_suffix(".data"), meta
    return None


PICKUP_3D = ("Uvel", "Vvel", "Theta", "Salt", "GuNm1", "GvNm1", "GtNm1", "GsNm1")
PICKUP_2D = ("EtaN", "dEtaHdt", "EtaH")


def _pickup_fields(path, meta, nz):
    """{field: array} of a MITgcm pickup: the 3-D fields of PICKUP_3D take nz records each, the
    2-D fields of PICKUP_2D one. ValueError for an unknown field or if the records consumed do not
    equal the meta's nrecords."""
    ny, nx = meta["shape"][-2:]
    unknown = [f for f in meta["fldList"] if f not in PICKUP_3D + PICKUP_2D]
    if unknown:
        raise ValueError(f"{path}: unknown pickup field layout {unknown}")
    n_rec = sum(nz if f in PICKUP_3D else 1 for f in meta["fldList"])
    if n_rec != meta["nrecords"]:
        raise ValueError(f"{path}: fields need {n_rec} records, the meta says {meta['nrecords']}")
    a = np.fromfile(path, ">f8").reshape(meta["nrecords"], ny, nx)
    out, r = {}, 0
    for f in meta["fldList"]:
        n = nz if f in PICKUP_3D else 1
        out[f] = a[r:r + n].astype(float) if n > 1 else a[r].astype(float)
        r += n
    return out


def eos_linear(T, S, eos: dict):
    """rho = rhoNil (1 - tAlpha (T - tRef) + sBeta (S - sRef)) [kg m^-3]
    (model/src/find_rho.F:101-108, LINEAR; tRef, sRef per level)."""
    tr = np.asarray(eos["tRef"])[:, None, None]
    sr = np.asarray(eos["sRef"])[:, None, None]
    return eos["rhoNil"] * (1.0 - eos["tAlpha"] * (T - tr) + eos["sBeta"] * (S - sr))


def run(run_dir: Path, final_from_pickup: bool = True) -> xr.Dataset:
    """The run as an xarray.Dataset (grid eager, time series lazy). See the module docstring."""
    dirs = stage_dirs(run_dir)
    for d in dirs:
        require_ok(d)
    g0 = dirs[0]
    nml = _nml(g0)
    p1 = nml.get("PARM01", {})
    data_vars, coords = {}, {}
    nz = _grid(g0 / "DRF").size
    for f, name in GRID1D.items():
        if (g0 / f"{f}.meta").exists():
            a = np.asarray(_grid(g0 / f)).reshape(-1)
            data_vars[name] = (("k" if a.size == nz else "kp1",), a)
    for f, name in GRID2D.items():
        if (g0 / f"{f}.meta").exists():
            a = _grid(g0 / f)
            data_vars[name] = (("j", "i"), np.asarray(a).reshape(a.shape[-2:]))
    for f, name in GRID3D.items():
        a = _grid(g0 / f)
        data_vars[name] = (("k", "j", "i"), np.asarray(a).reshape((nz,) + a.shape[-2:]))
    ny, nx = data_vars["hFacC"][1].shape[1:]
    rho_nil = _num(p1.get("rhonil"), 999.8)
    rho_const = _num(p1.get("rhoconst"), rho_nil)                 # ini_parms.F:483
    eos = {"eosType": _str(p1.get("eostype")) or "LINEAR", "tAlpha": _num(p1.get("talpha"), 2e-4),
           "sBeta": _num(p1.get("sbeta"), 7.4e-4), "rhoNil": rho_nil,
           "tRef": _list(p1.get("tref"), nz), "sRef": _list(p1.get("sref"), nz)}
    p3 = nml.get("PARM03", {})
    dt = _num(p3.get("deltat"), None) or _num(p3.get("deltatclock"), None) or \
        _num(p3.get("deltatmom"), None)
    stages = []
    for d in dirs:
        n = _nml(d)
        q3 = n.get("PARM03", {})
        q5 = {k: v for g in n.values() for k, v in g.items()}       # file names: any group
        stages.append({"dir": str(d), "nIter0": int(_num(q3.get("niter0"), 0)),
                       "nTimeSteps": int(_num(q3.get("ntimesteps"), 0)),
                       "periodic": (q3.get("periodicexternalforcing", ".FALSE.").upper()
                                    .strip(".") == "TRUE"),
                       "period": _num(q3.get("externforcingperiod"), 0.0),
                       "cycle": _num(q3.get("externforcingcycle"), 0.0),
                       "zonalWindFile": _str(q5.get("zonalwindfile")),
                       "meridWindFile": _str(q5.get("meridwindfile"))})

    files = _diag_files(dirs)
    need = [v for v in ("U", "V", "RHO", "ETA") if v in files]
    its = sorted(set.intersection(*(set(files[v]) for v in need))) if need else []
    if not its:
        raise RuntimeError(f"{run_dir}: no common diagnostics records of {need}")
    present = [v for v in files if set(its) <= set(files[v])]
    items = {v: [] for v in present}
    for it in its:
        for v in present:
            path, nrec, rec, shape = files[v][it]
            off = rho_const if v == "RHO" else 0.0
            items[v].append((path, nrec, rec, tuple(shape), off))
    src = "diagnostics"
    last = stages[-1]
    it_end = last["nIter0"] + last["nTimeSteps"]
    step = its[-1] - its[-2] if len(its) > 1 else None
    if final_from_pickup and step and it_end == its[-1] + step and \
            not ({"TAUX", "TAUY"} & set(present)) and eos["eosType"] == "LINEAR":
        pk = _pickup(Path(last["dir"]), it_end)
        if pk is not None:
            f = _pickup_fields(*pk, nz)
            vals = {v: f[PICKUP[v]] for v in PICKUP if v in present}
            if "RHO" in present:
                vals["RHO"] = eos_linear(f["Theta"], f["Salt"], eos)
            if set(vals) == set(present):
                for v in present:
                    items[v].append(("array", vals[v]))
                its = its + [it_end]
                src = f"pickup {pk[0].name} (iteration {it_end}) for the last record"
    for v in present:
        shape3 = items[v][0][3] if items[v][0][0] != "array" else items[v][0][1].shape
        dims = ("time", "k", "j", "i") if len(shape3) == 3 else ("time", "j", "i")
        data_vars[v] = xr.Variable(dims, indexing.LazilyIndexedArray(_Records(items[v], shape3)))
    coords["time"] = np.array(its, float) * dt
    coords["iter"] = ("time", np.array(its, int))
    ds = xr.Dataset(data_vars, coords=coords)
    ds.attrs.update(run_dir=str(run_dir), rho0=rho_const, g=_num(p1.get("gravity"), 9.81),
                    eos=eos, deltaT=dt, stages=stages, final_record_source=src, shape=(nz, ny, nx))
    return ds


def forcing_stress(ds: xr.Dataset, t: float):
    """(taux at U points, tauy at V points) [N m^-2] that the model applies at time t from the
    stage's forcing files: record k centred at (k + 1/2) period, linear interpolation, wrapped
    on nbRec = externForcingCycle / externForcingPeriod from the namelist (not on the file's record
    count; eesupp/src/get_periodic_interval.F:83, 117-127;
    model/src/external_fields_load.F:268-277). The stage is the one with
    nIter0 deltaT <= t < (nIter0 + nTimeSteps) deltaT (the last stage holds its end); a stage
    without wind files gives 0 (the thesis stage 2, balance_kang_fringer.tau_modelo)."""
    nz, ny, nx = ds.attrs["shape"]
    dt = ds.attrs["deltaT"]
    st = ds.attrs["stages"]
    s = st[-1]
    for c in st:
        if c["nIter0"] * dt <= t < (c["nIter0"] + c["nTimeSteps"]) * dt:
            s = c
            break
    out = []
    for key in ("zonalWindFile", "meridWindFile"):
        f = s[key]
        if not f:
            out.append(np.zeros((ny, nx)))
            continue
        if not s["periodic"]:
            # external_fields_load.F:63 only interpolates with periodicExternalForcing
            raise NotImplementedError(f"{key} without periodicExternalForcing is not handled")
        p = Path(s["dir"]) / f
        R = np.fromfile(p, ">f8").reshape(-1, ny, nx)
        per, cyc = s["period"], s["cycle"]
        if cyc > 0:
            # periodic sequence, eesupp/src/get_periodic_interval.F:83 (nbRec) and 117-127
            nb = int(round(cyc / per))
            if abs(nb * per - cyc) > 1e-9 * cyc:
                raise ValueError(f"externForcingCycle {cyc} is not a multiple of the period {per}")
            if R.shape[0] < nb:
                raise ValueError(f"{p}: {R.shape[0]} records < cycle/period = {nb}")
            loc = t - 0.5 * per + cyc * (2 - round(t / cyc))
            tm = loc % cyc
            k1 = int(tm // per)
            w2 = (tm - per * k1) / per
            out.append((1 - w2) * R[k1] + w2 * R[(k1 + 1) % nb])
        else:
            # non-periodic sequence, get_periodic_interval.F:98-107
            loc = t - 0.5 * per
            k1 = int(np.floor(loc / per))
            w2 = (loc - per * k1) / per
            out.append((1 - w2) * R[k1] + w2 * R[k1 + 1])
    return out[0], out[1]
