#!/opt/anaconda3/envs/pyocean/bin/python
r"""ISOMIP+ metrics for Exp. 2 (MITgcm ice-shelf cavity).

The binding definitions are those of ``notes/pd0_preregistro.md`` (pd0); diagnostic names,
units and signs come from ``notes/pd1_config.md`` (pd1). AD16 = Asay-Davis et al. (2016),
Y26 = Yung et al. (2026, final version; "L" = line of the local .md).

Usage::

    metrics.py <run_dir>                 # isomip_summary -> <run_dir>/metrics.json
    metrics.py <run_dir> --gade [--record K]   # post hoc Gade-line check (R10)
    metrics.py <run_dir> --ustar [--record K]  # u* consistency check

Conventions
-----------
* Calendar (pd0 "Definiciones comunes"): 365-day year, month = 2 628 000 s. Records are the
  monthly-mean diagnostics files in time order, 0-based. "Last 6 months" of a 2-year run =
  records 18-23 (``LAST6``); "year 20" of Ocean1 = records 228-239 (``YEAR20``).
* Melt (pd0; pd1 sect. 1.3): ``SHIfwFlx`` [kg m^-2 s^-1] is negative for melting;
  m [m a^-1 of freshwater, positive = melting] = -SHIfwFlx / rho_fw * 31 536 000,
  rho_fw = 1000 kg m^-3 (AD16 Table 4).
* Shelf (pd0, AD16 sect. 3.1.2): wet columns whose ice base (``icetopo.bin``, after calving
  and hFac snapping) is < 0. Deep shelf: base < -300 m (AD16 Eq. 37). Means are area-weighted.
* Time means over several records are weighted by the record length (timeInterval of the
  .meta file); for calendar windows all weights are equal.
* Tracer TRAC01 = meltwater volume fraction C (pd1 sect. 4). It is never clipped: C < 0
  under refreezing is physical for a volume-fraction budget.

Every run is accepted by ``pfcommon.runctl.require_ok`` before anything of it is read
(Review Focus 5); staged runs also need ``run_status_final.json`` (pfcommon.campaign).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pfcommon.mitgcm_io import _parse_meta, read_mds  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402
from exp2_isomip_cavity.theory import melt as ml  # noqa: E402

# ---- constants -----------------------------------------------------------------------------
MONTH_S = 2628000.0                      # pd0: 365-day calendar, month = 2 628 000 s
SECONDS_PER_YEAR = ml.SECONDS_PER_YEAR   # 31 536 000 s
RHO_FW = ml.RHO_FW                       # 1000 kg m^-3, AD16 Table 4
SV = 1.0e6                               # m^3 s^-1 per Sverdrup
LAST6 = range(18, 24)                    # pd0: months 19-24 of a 2-year run (0-based 18-23)
PREV6 = range(12, 18)                    # pd0 N1: months 13-18
YEAR20 = range(228, 240)                 # pd0: Ocean1 year 20 = records 229-240
YEAR19 = range(216, 228)                 # pd0 N1 (Ocean1 form): year 19
N_RECORDS_2YR = 24
N_RECORDS_OCEAN1 = 240
X_CAVITY_MAX_M = 640e3                   # pd0 V2b/V2c, Y26 Fig. 13: cavity x < 640 km
Y_CAVITY_MIN_M = 50e3                    # pd0 V2c, Y26 L1118: y > 50 km
X_FRONT_BAND_M = (630e3, 650e3)          # pd0 V2d, Y26 L1155: 630 < x < 650 km
DEEP_DRAFT_M = -300.0                    # pd0, AD16 Eq. 37: base < -300 m
S_BOT_TF = 34.7                          # pd0 H1: Sbot = 34.7 for T_f,GL
TOL_RANGE = 0.1                          # pd1 S3: profile range +- 0.1

_SHI = "shiDiag"
_DYN = "dynDiag"
_TIME_INTERVAL = re.compile(r"timeInterval\s*=\s*\[(.*?)\]\s*;", re.S)


# ---- run loading ----------------------------------------------------------------------------

@dataclass
class Record:
    """One monthly-mean record: diagnostics files ``<kind>.<iter:010d>`` in ``stage_dir``."""
    iter: int
    t0: float
    t1: float
    stage_dir: Path

    @property
    def period(self) -> float:
        return self.t1 - self.t0


@dataclass
class RunData:
    run_dir: Path
    case: dict
    params: dict
    actual: dict
    stage_dirs: list
    records: list
    nx: int
    ny: int
    nz: int
    hfac: np.ndarray            # hFacC (nz, ny, nx)
    hfacW: np.ndarray
    hfacS: np.ndarray
    rA: np.ndarray              # RAC (ny, nx) [m^2]
    drF: np.ndarray             # (nz,) [m]
    dyG: np.ndarray             # (ny, nx) [m], length of the west (U) face
    dxG: np.ndarray             # (ny, nx) [m], length of the south (V) face
    xc: np.ndarray
    yc: np.ndarray
    xg: np.ndarray
    yg: np.ndarray
    rc: np.ndarray              # (nz,) cell-centre z [m, negative down]
    rf: np.ndarray              # (nz+1,)
    draft: np.ndarray           # ice base (ny, nx) [m], icetopo.bin
    t_start: float
    t_end: float
    _meta_cache: dict = field(default_factory=dict, repr=False)

    @property
    def wet3(self) -> np.ndarray:
        return self.hfac > 0

    @property
    def wet(self) -> np.ndarray:
        return self.wet3.any(axis=0)

    @property
    def shelf(self) -> np.ndarray:
        """pd0: wet columns with ice base < 0 after calving (AD16 sect. 3.1.2)."""
        return self.wet & (self.draft < 0)

    @property
    def ktop(self) -> np.ndarray:
        """Index of the top wet cell of each column (-1 where dry)."""
        return np.where(self.wet, np.argmax(self.wet3, axis=0), -1)

    @property
    def n_records(self) -> int:
        return len(self.records)

    @property
    def vol(self) -> np.ndarray:
        """Static cell volume hFacC drF rA [m^3]."""
        return self.hfac * self.drF[:, None, None] * self.rA[None]

    @property
    def hcol(self) -> np.ndarray:
        """Static water-column thickness sum(hFacC drF) [m]."""
        return (self.hfac * self.drF[:, None, None]).sum(axis=0)

    def rstar(self, eta: np.ndarray) -> np.ndarray:
        """rStarFacC = (eta + H) / H, H = static column thickness
        (model/src/calc_r_star.F:103-105); 1 in dry columns."""
        h = self.hcol
        return np.where(h > 0, (eta + h) / np.where(h > 0, h, 1.0), 1.0)

    def field(self, kind: str, name: str, k: int) -> np.ndarray:
        """Field ``name`` of record ``k`` of diagnostics list ``kind`` (shiDiag | dynDiag).
        Raises if it holds a non-finite value (a diverged record must never be averaged)."""
        rec = self.records[k]
        stem = rec.stage_dir / f"{kind}.{rec.iter:010d}"
        key = str(stem)
        if key not in self._meta_cache:
            self._meta_cache[key] = _parse_meta(Path(f"{stem}.meta").read_text())
        meta = self._meta_cache[key]
        names = meta["fldList"]
        if name not in names:
            raise KeyError(f"{stem}: field {name} not in {names}")
        shape = tuple(meta["shape"])
        n = int(np.prod(shape))
        mm = np.memmap(f"{stem}.data", dtype=">f8", mode="r", offset=names.index(name) * n * 8,
                       shape=(n,))
        arr = np.array(mm, dtype=float).reshape(shape)
        if not np.isfinite(arr).all():
            raise RuntimeError(f"{stem}: non-finite values in {name}; refusing to use it")
        return arr

    def mean_field(self, kind: str, name: str, records) -> np.ndarray:
        """Record-length-weighted time mean over ``records``."""
        recs = _as_records(self, records)
        w = np.array([self.records[k].period for k in recs])
        acc = None
        for k, wk in zip(recs, w):
            f = self.field(kind, name, k) * wk
            acc = f if acc is None else acc + f
        return acc / w.sum()


def _as_records(rd: RunData, records) -> list:
    if isinstance(records, (int, np.integer)):
        records = [int(records)]
    recs = [int(k) for k in records]
    if not recs:
        raise ValueError("empty record list")
    for k in recs:
        if not 0 <= k < rd.n_records:
            raise IndexError(f"record {k} outside 0..{rd.n_records - 1}")
    return recs


def _stage_dirs(run_dir: Path) -> list:
    """Stage directories in time order (pfcommon.campaign.run_stages layout: stage 1 in
    run_dir, stage k > 1 in run_dir/stage<k>). For staged runs run_status_final.json must be
    OK and complete, and every stage must pass require_ok."""
    fin = run_dir / "run_status_final.json"
    if not fin.exists():
        if any(run_dir.glob("stage[0-9]*")):
            raise RuntimeError(f"{run_dir}: staged run without run_status_final.json; refusing it")
        return [run_dir]
    try:
        st = json.loads(fin.read_text())
    except json.JSONDecodeError as e:
        raise RuntimeError(f"{run_dir}: corrupt run_status_final.json ({e})") from e
    if st.get("state") != "OK" or st.get("last_iter") != st.get("expected_last_iter"):
        raise RuntimeError(f"{run_dir}: run_status_final.json is {st.get('state')} "
                           f"(last_iter {st.get('last_iter')}, expected "
                           f"{st.get('expected_last_iter')}); refusing it")
    keys = sorted(int(k) for k in st.get("stages", {}))
    if keys != list(range(1, len(keys) + 1)) or not keys:
        raise RuntimeError(f"{run_dir}: run_status_final.json stages {keys} are not 1..n")
    dirs = [run_dir if k == 1 else run_dir / f"stage{k}" for k in keys]
    for k, d in zip(keys, dirs):
        s = st["stages"][str(k)]
        if not s or s.get("state") != "OK":
            raise RuntimeError(f"{run_dir}: stage {k} is not OK in run_status_final.json")
        if k > 1:
            require_ok(d)
    return dirs


def _time_interval(meta_file: Path) -> tuple[float, float]:
    m = _TIME_INTERVAL.search(meta_file.read_text())
    if not m:
        raise RuntimeError(f"{meta_file}: no timeInterval (not a time-averaged diagnostic)")
    t0, t1 = (float(v.replace("D", "E").replace("d", "e")) for v in m.group(1).split())
    return t0, t1


def _records(run_dir: Path, stage_dirs: list, t_start: float, t_end: float) -> list:
    recs = []
    for sd in stage_dirs:
        its = {}
        for kind in (_SHI, _DYN):
            its[kind] = sorted(int(f.name.split(".")[1]) for f in sd.glob(f"{kind}.[0-9]*.meta"))
        if its[_SHI] != its[_DYN]:
            raise RuntimeError(f"{sd}: {_SHI} and {_DYN} records differ "
                               f"({len(its[_SHI])} vs {len(its[_DYN])})")
        for it in its[_SHI]:
            t0, t1 = _time_interval(sd / f"{_SHI}.{it:010d}.meta")
            if _time_interval(sd / f"{_DYN}.{it:010d}.meta") != (t0, t1):
                raise RuntimeError(f"{sd}: {_SHI} and {_DYN} intervals differ at {it}")
            recs.append(Record(it, t0, t1, sd))
    recs.sort(key=lambda r: r.t1)
    if not recs:
        raise RuntimeError(f"{run_dir}: no monthly diagnostics records")
    tol = 1e-6
    # Review Focus 5: a gap must be reported, never interpolated
    if abs(recs[0].t0 - t_start) > tol or abs(recs[-1].t1 - t_end) > tol:
        raise RuntimeError(f"{run_dir}: records cover {recs[0].t0:g}-{recs[-1].t1:g} s, "
                           f"run is {t_start:g}-{t_end:g} s")
    for a, b in zip(recs[:-1], recs[1:]):
        if abs(a.t1 - b.t0) > tol:
            raise RuntimeError(f"{run_dir}: gap between records ending {a.t1:g} s and "
                               f"starting {b.t0:g} s (records are not contiguous)")
    return recs


def _grid(rd: Path, name: str, shape) -> np.ndarray:
    """Model grid file <name>.data (big-endian float64, useSingleCpuIO)."""
    return np.fromfile(rd / f"{name}.data", dtype=">f8").reshape(shape).astype(float)


def load_run(run_dir) -> RunData:
    """Load a finished run: require_ok first (Review Focus 5), then the stage layout,
    case.json, the model grid and the list of monthly records (concatenated over stages
    in time order). Fields are read lazily, record by record."""
    run_dir = Path(run_dir)
    require_ok(run_dir)                                   # before anything else is read
    stage_dirs = _stage_dirs(run_dir)
    case = json.loads((run_dir / "case.json").read_text())
    p, actual = case["params"], case.get("actual", {})
    if p.get("melt_law") != ml.MELT_LAW:                  # R21: void u* law
        raise RuntimeError(f"{run_dir}: melt_law {p.get('melt_law')!r} != {ml.MELT_LAW!r}: run "
                           "made with the void pre-R21 u* law (Cd u^2 + u_tidal^2), refused")
    nx, ny, nz = int(p["nx"]), int(p["ny"]), int(p["nz"])
    dt = float(p["deltaT"])
    n0, nts = int(p["nIter0"]), int(p["nTimeSteps"])
    fin = run_dir / "run_status_final.json"
    if fin.exists():
        exp = json.loads(fin.read_text()).get("expected_last_iter")
        if exp != n0 + nts:
            raise RuntimeError(f"{run_dir}: expected_last_iter {exp} != nIter0 + nTimeSteps "
                               f"{n0 + nts} of case.json")
    t_start, t_end = n0 * dt, (n0 + nts) * dt
    recs = _records(run_dir, stage_dirs, t_start, t_end)
    g3, g2 = (nz, ny, nx), (ny, nx)
    return RunData(
        run_dir=run_dir, case=case, params=p, actual=actual, stage_dirs=stage_dirs,
        records=recs, nx=nx, ny=ny, nz=nz,
        hfac=_grid(run_dir, "hFacC", g3), hfacW=_grid(run_dir, "hFacW", g3),
        hfacS=_grid(run_dir, "hFacS", g3), rA=_grid(run_dir, "RAC", g2),
        drF=_grid(run_dir, "DRF", (nz,)), dyG=_grid(run_dir, "DYG", g2),
        dxG=_grid(run_dir, "DXG", g2), xc=_grid(run_dir, "XC", g2), yc=_grid(run_dir, "YC", g2),
        xg=_grid(run_dir, "XG", g2), yg=_grid(run_dir, "YG", g2),
        rc=_grid(run_dir, "RC", (nz,)), rf=_grid(run_dir, "RF", (nz + 1,)),
        draft=np.fromfile(run_dir / "icetopo.bin", dtype=">f8").reshape(g2).astype(float),
        t_start=t_start, t_end=t_end)


# ---- calendar -------------------------------------------------------------------------------

_WINDOWS = {"last6": (LAST6, N_RECORDS_2YR), "prev6": (PREV6, N_RECORDS_2YR),
            "year20": (YEAR20, N_RECORDS_OCEAN1), "year19": (YEAR19, N_RECORDS_OCEAN1)}


def calendar_window(rd: RunData, name: str) -> range:
    """Pre-registered record window (pd0 "Definiciones comunes", N1): 'last6' / 'prev6' of a
    24-record (2-year) run, 'year20' / 'year19' of a 240-record Ocean1 run. Asserts the run
    has exactly that many records, each one month long, starting at the run start."""
    if name not in _WINDOWS:
        raise ValueError(f"unknown window {name!r}; one of {sorted(_WINDOWS)}")
    win, n = _WINDOWS[name]
    if rd.n_records != n:
        raise ValueError(f"window {name!r} needs a run of exactly {n} monthly records, "
                         f"this one has {rd.n_records}")
    for k, r in enumerate(rd.records):
        if abs(r.period - MONTH_S) > 1e-6 or abs(r.t0 - (rd.t_start + k * MONTH_S)) > 1e-6:
            raise ValueError(f"record {k} spans {r.t0:g}-{r.t1:g} s, not month {k + 1} "
                             f"of a {MONTH_S:g} s calendar")
    return win


# ---- melt -----------------------------------------------------------------------------------

def melt_rate_m_per_a(fw_flux) -> np.ndarray:
    """Melt rate [m a^-1 of freshwater, positive = melting] from SHIfwFlx [kg m^-2 s^-1,
    negative = melting] (pd1 sect. 1.3; AD16 Eq. 26 gives SHIfwFlx = -rho_fw m):
    m = -SHIfwFlx / rho_fw * 31 536 000."""
    return -np.asarray(fw_flux, dtype=float) / RHO_FW * SECONDS_PER_YEAR


def _mask(rd: RunData, mask: str) -> np.ndarray:
    if mask == "shelf":
        return rd.shelf
    if mask == "deep300":
        return rd.shelf & (rd.draft < DEEP_DRAFT_M)      # AD16 Eq. 37
    raise ValueError(f"mask must be 'shelf' or 'deep300', not {mask!r}")


def mean_melt(rd: RunData, records, mask: str = "shelf") -> float:
    """Area-weighted mean melt [m a^-1] over ``mask`` ('shelf' | 'deep300', pd0) of the
    time-mean SHIfwFlx over ``records``."""
    sel = _mask(rd, mask)
    if not sel.any():
        raise ValueError(f"mask {mask!r} is empty")
    m = melt_rate_m_per_a(rd.mean_field(_SHI, "SHIfwFlx", records))
    return float((m * rd.rA)[sel].sum() / rd.rA[sel].sum())


def total_melt_flux(rd: RunData, records) -> float:
    """Total meltwater mass flux [kg s^-1, positive = melting]: sum of -SHIfwFlx rA over the
    shelf (SHIfwFlx is zero elsewhere, pd1 sect. 1.3), time mean over ``records``."""
    fw = rd.mean_field(_SHI, "SHIfwFlx", records)
    return float((-fw * rd.rA)[rd.shelf].sum())


def mean_melt_per_record(rd: RunData, mask: str = "shelf") -> list:
    return [mean_melt(rd, [k], mask) for k in range(rd.n_records)]


# ---- streamfunctions ------------------------------------------------------------------------

def _u_transport(rd: RunData, records) -> np.ndarray:
    """Zonal volume transport per cell face, UVEL hFacW drF dyG [m^3 s^-1] at U points
    (west faces), from the time-mean UVEL. Static hFacW: the r* stretching (relative size
    eta/H, below 1e-3 here) is not applied because UVEL is not a thickness-weighted mean."""
    u = rd.mean_field(_DYN, "UVEL", records)
    return u * rd.hfacW * rd.drF[:, None, None] * rd.dyG[None]


def barotropic_streamfunction(rd: RunData, records) -> np.ndarray:
    """psi(x, y) = -int_0^y U dy' [Sv], U = vertically integrated zonal transport
    (pd0; AD16 sect. 3.3 output list: U = -d psi / dy). Shape (ny + 1, nx): psi at the
    cell corners (XG[i], y_j), y_j = j dy, psi(y = 0) = 0."""
    U = _u_transport(rd, records).sum(axis=0)             # (ny, nx), m^3 s^-1
    psi = np.zeros((rd.ny + 1, rd.nx))
    psi[1:] = -np.cumsum(U, axis=0)
    return psi / SV


def overturning_streamfunction(rd: RunData, records) -> np.ndarray:
    """psi(x, z) = int_z^0 int u dy dz' [Sv], so that the y-integrated u = -d psi / dz
    (AD16 sect. 3.3 output list); integrated from the top, psi = 0 at z = 0. Shape
    (nz + 1, nx): psi at (XG[i], RF[k])."""
    Uy = _u_transport(rd, records).sum(axis=1)            # (nz, nx), m^3 s^-1
    psi = np.zeros((rd.nz + 1, rd.nx))
    psi[1:] = np.cumsum(Uy, axis=0)
    return psi / SV


def _corner_y(rd: RunData) -> np.ndarray:
    """y of the psi corners: YG of each row plus the northern edge of the last row."""
    return np.concatenate([rd.yg[:, 0], [rd.yg[-1, 0] + rd.dyG[-1, 0]]])


def psi_bt_max(rd: RunData, records, region: str) -> float:
    """max |psi_bt| [Sv] of the time-mean barotropic streamfunction over a pd0 region:
    'cavity' = V2c, y > 50 km and x < 640 km (Y26 L1118); 'front' = V2d, 630 < x < 650 km
    (Y26 L1155). Ruling R13: the absolute value, so that the result does not depend on the
    gyre's sense of rotation or on the sign convention. Regions are selected by corner
    coordinates (XG, y); land rows carry no transport (hFacW = 0), so they add no extremum."""
    psi = barotropic_streamfunction(rd, records)
    x = rd.xg[0][None, :]
    y = _corner_y(rd)[:, None]
    if region == "cavity":
        sel = (y > Y_CAVITY_MIN_M) & (x < X_CAVITY_MAX_M)
    elif region == "front":
        sel = (x > X_FRONT_BAND_M[0]) & (x < X_FRONT_BAND_M[1]) & np.ones_like(y, bool)
    else:
        raise ValueError(f"region must be 'cavity' or 'front', not {region!r}")
    if not sel.any():
        raise ValueError(f"region {region!r} selects no psi points")
    return float(np.abs(psi[sel]).max())


def overturning_max(rd: RunData, records) -> float:
    """pd0 V2b: max |psi(x, z)| [Sv] in the cavity, x < 640 km (Y26 L1100-1102, Fig. 8)."""
    psi = overturning_streamfunction(rd, records)
    return float(np.abs(psi[:, rd.xg[0] < X_CAVITY_MAX_M]).max())


def circulation_strength(rd: RunData, record: int, kind: str = "bt") -> float:
    """Circulation strength of one month [Sv] = max - min of psi for x < 640 km
    (pd0; Y26 Fig. 13 caption, L1365-1369). kind 'bt' = barotropic, 'moc' = overturning."""
    if not isinstance(record, (int, np.integer)):
        raise TypeError("circulation_strength is defined for a single monthly record")
    if kind == "bt":
        psi = barotropic_streamfunction(rd, [record])
    elif kind == "moc":
        psi = overturning_streamfunction(rd, [record])
    else:
        raise ValueError(f"kind must be 'bt' or 'moc', not {kind!r}")
    sub = psi[:, rd.xg[0] < X_CAVITY_MAX_M]
    return float(sub.max() - sub.min())


# ---- N1, N2 ---------------------------------------------------------------------------------

def equilibrium_drift(rd: RunData, form: str | None = None) -> float:
    """pd0 N1: |<m>(late) - <m>(early)| / <m>(late), <m> = whole-shelf mean melt.
    form '2yr': late = months 19-24, early = 13-18 (24 records); 'ocean1': late = year 20,
    early = year 19 (240 records). ``None`` picks the form from the number of records."""
    if form is None:
        form = {N_RECORDS_2YR: "2yr", N_RECORDS_OCEAN1: "ocean1"}.get(rd.n_records)
        if form is None:
            raise ValueError(f"N1 needs 24 or 240 monthly records, run has {rd.n_records}")
    if form == "2yr":
        late, early = calendar_window(rd, "last6"), calendar_window(rd, "prev6")
    elif form == "ocean1":
        late, early = calendar_window(rd, "year20"), calendar_window(rd, "year19")
    else:
        raise ValueError(f"form must be '2yr' or 'ocean1', not {form!r}")
    a, b = mean_melt(rd, late), mean_melt(rd, early)
    return abs(a - b) / a


def _read_state(stage_dir: Path, name: str, itr: int) -> np.ndarray:
    arr = read_mds(stage_dir / name, itr)[0][0]
    if not np.isfinite(arr).all():
        raise RuntimeError(f"{stage_dir}/{name}.{itr:010d}: non-finite values")
    return arr


def _inventory(rd: RunData, stage_dir: Path, itr: int) -> float:
    C = _read_state(stage_dir, "PTRACER01", itr)          # never clipped (C < 0 is physical)
    eta = _read_state(stage_dir, "Eta", itr)
    return float((C * rd.vol * rd.rstar(eta)[None]).sum())


def tracer_budget(rd: RunData) -> dict:
    """pd0 N2 legs over the whole run [m^3 of meltwater]:

    * inventory I = sum C hFacC drF rA rStarFacC (pd1 sect. 6 S2), at nIter0 and at the end;
    * source = sum_k period_k * sum(-<SHIfwFlx>_k rA) / rho_fw (pd1 sect. 4: the tracer
      source is -SHIfwFlx / rho_fw per unit ice-base area, config/ocean0/code/
      shelfice_forcing_mwtr.F);
    * rbcs sink = sum_k period_k * sum(mask / tau <C>_k hFacC drF rA rStarFacC_k), with
      rStarFacC_k from the record-mean ETAN (the restoring term is linear in C; the product
      of record means of C and rStarFacC is the only approximation, of order the eta
      variance over H^2).

    Integrated over ALL records (deferred Task 2 finding i). Both legs must be present."""
    mask_f = rd.run_dir / "rbcs_mask.bin"
    tau = rd.params.get("tauRelax")
    if not mask_f.exists() or not tau or float(tau) <= 0:
        raise RuntimeError(f"{rd.run_dir}: rbcs sink leg missing (rbcs_mask.bin / tauRelax)")
    mask = np.fromfile(mask_f, dtype=">f8").reshape(rd.nz, rd.ny, rd.nx).astype(float)
    tau = float(tau)
    src = sink = 0.0
    for k, r in enumerate(rd.records):
        fw = rd.field(_SHI, "SHIfwFlx", k)
        C = rd.field(_DYN, "TRAC01", k)                   # unclipped
        eta = rd.field(_SHI, "ETAN", k)
        src += float((-fw * rd.rA).sum()) / RHO_FW * r.period
        sink += float((mask / tau * C * rd.vol * rd.rstar(eta)[None]).sum()) * r.period
    if not src > 0:
        raise RuntimeError(f"{rd.run_dir}: tracer source leg is {src:g} m^3; N2 is undefined")
    n0 = int(rd.params["nIter0"])
    n_end = n0 + int(rd.params["nTimeSteps"])
    inv0 = _inventory(rd, rd.stage_dirs[0], n0)
    inv1 = _inventory(rd, rd.stage_dirs[-1], n_end)
    d_inv = inv1 - inv0
    res = abs(d_inv - (src - sink)) / src
    return {"value": res, "pass": bool(res < 0.01), "inventory_start_m3": inv0,
            "inventory_end_m3": inv1, "delta_inventory_m3": d_inv, "source_m3": src,
            "rbcs_sink_m3": sink, "n_records": rd.n_records,
            "definition": "pd0 N2: |dInv - (int source - int rbcs sink)| / int source < 1 %"}


def tracer_budget_residual(rd: RunData) -> float:
    """pd0 N2 relative residual (see tracer_budget)."""
    return tracer_budget(rd)["value"]


# ---- H1, H2 ---------------------------------------------------------------------------------

def freezing_point_gl(rd: RunData, S: float = S_BOT_TF) -> dict:
    """pd0 H1: T_f,GL = lam1 S + lam2 + lam3 p (AD16 Eq. 25, Table 4) with S = Sbot = 34.7
    and p = rho_ref g D [Pa] at D = the deepest wet ice base after calving
    (case.json actual.deepest_wet_ice_base_depth_m; pd0 "la base mas profunda tras el
    vaciado"). Not case.json's grounding_line_* (those are the first wet ice column, the wet
    edge: deferred Task 1 finding)."""
    D = rd.actual.get("deepest_wet_ice_base_depth_m")
    if D is None:
        raise RuntimeError(f"{rd.run_dir}: case.json has no deepest_wet_ice_base_depth_m")
    p = float(rd.params["rho_ref"]) * float(rd.params["gravity"]) * float(D)
    tf = ml.LAMBDA1 * S + ml.LAMBDA2 + ml.LAMBDA3 * p
    return {"T_f_GL_degC": float(tf), "depth_m": D, "pressure_Pa": p, "S": S,
            "source": "AD16 Eq. 25 (Table 4); p = rho_ref g D, D = deepest wet ice base"}


def x_front_index(rd: RunData) -> np.ndarray:
    """pd0 H2: per row j, the open-ocean column just past the front = (max i of the shelf)
    + 1; -1 in rows without shelf."""
    out = np.full(rd.ny, -1)
    for j in range(rd.ny):
        idx = np.nonzero(rd.shelf[j])[0]
        if idx.size:
            out[j] = idx.max() + 1
    return out


def outflow_depth(rd: RunData, records, x_front_m: float | None = None,
                  face: str = "centre") -> float:
    """pd0 H2: z_out [m, negative] = transport-weighted mean depth of the outgoing tracer
    transport at x_front, counting only u C > 0, over ``records`` (last 6 months).

    UTRAC01 = u C hFacW at U points (pd1 sect. 8), so the meltwater transport per face is
    UTRAC01 drF dyG [m^3 s^-1]. The time mean is taken first, then u C > 0 is selected.
    x_front: pd0's per-row column (x_front_index) or, if ``x_front_m`` is given, the column
    whose XC is nearest. face 'centre' interpolates to the centre of that column (mean of
    its west and east faces, the C-grid interpolation to x = x_front); 'west' uses its west
    face (the ice-front face). Depth weight: RC of the level."""
    ut = rd.mean_field(_DYN, "UTRAC01", records) * rd.drF[:, None, None] * rd.dyG[None]
    if x_front_m is None:
        cols = x_front_index(rd)
    else:
        cols = np.full(rd.ny, int(np.argmin(np.abs(rd.xc[0] - x_front_m))))
    num = den = 0.0
    for j, i in enumerate(cols):
        if i < 0:
            continue
        if face == "centre":
            f = 0.5 * (ut[:, j, i] + ut[:, j, (i + 1) % rd.nx])
        elif face == "west":
            f = ut[:, j, i]
        else:
            raise ValueError(f"face must be 'centre' or 'west', not {face!r}")
        pos = np.where(f > 0, f, 0.0)
        num += float((pos * rd.rc).sum())
        den += float(pos.sum())
    if den <= 0:
        raise RuntimeError("no outgoing tracer transport at x_front")
    return num / den


# ---- u* check (deferred Task 2 finding iii) -------------------------------------------------

def boundary_layer_velocity_sq(rd: RunData, U: np.ndarray, V: np.ndarray) -> np.ndarray:
    """|u_BL|^2 at C points as config/ocean0/code/shelfice_thermodynamics.F computes it with
    SHELFICEboundaryLayer (which forces SHELFICE_oldCalcUStar, shelfice_readparms.F:208):
    uLoc(K) = recip_hFacC (u_i hFacW_i + u_i+1 hFacW_i+1) / 2 (lines 335-347), averaged over
    the top cell plus drKp1 = min(drF(K)(1 - hFacC(K)), drF(K+1) hFacC(K+1)) of the next
    (lines 373-418). Periodic neighbours (np.roll), as the model's exchanges."""
    def cpoint(F, h, axis):
        Fp = np.roll(F * h, -1, axis=axis)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(rd.hfac > 0, 1.0 / rd.hfac, 0.0)
        return r * 0.5 * (F * h + Fp)

    uc = cpoint(U, rd.hfacW, 2)
    vc = cpoint(V, rd.hfacS, 1)
    out = np.zeros((rd.ny, rd.nx))
    kt = rd.ktop
    for j, i in zip(*np.nonzero(rd.shelf)):
        k = kt[j, i]
        u, v = uc[k, j, i], vc[k, j, i]
        if k < rd.nz - 1:
            h0 = rd.drF[k] * rd.hfac[k, j, i]
            dr = max(min(rd.drF[k] * (1 - rd.hfac[k, j, i]), rd.drF[k + 1] * rd.hfac[k + 1, j, i]), 0.0)
            u = (u * h0 + dr * uc[k + 1, j, i]) / (h0 + dr)
            v = (v * h0 + dr * vc[k + 1, j, i]) / (h0 + dr)
        out[j, i] = u * u + v * v
    return out


def ustar_check(rd: RunData, record: int, rtol: float = 1e-3) -> dict:
    """Consistency of the friction velocity with AD16 Eq. 27 (PDF p. 15, ruling R21),
    u*^2 = Cd (|u_BL|^2 + u_tidal^2), on the shelf, for one monthly record:

    * floor: <SHIuStar> >= sqrt(Cd) u_tidal (exact: every instantaneous u* >= sqrt(Cd) u_tidal);
    * gammT identity: <SHIgammT> = Gamma_T <SHIuStar> (exact, linear in u*);
    * Jensen bound: u* is a convex function of the (linear) boundary-layer velocity, so
      <u*> >= sqrt(Cd (|BL(<u>)|^2 + u_tidal^2)) with <u> the monthly-mean UVEL/VVEL; checked
      to ``rtol`` (static hFac; the model's r*-weighted hFac differ by O(eta/H)).

    The equality itself cannot be checked: the diagnostics have no boundary-layer U/V
    (AD16's uBoundaryLayer) and no velocity snapshot coincides with a SHIuStar snapshot."""
    p = rd.params
    cd, ut, gT = float(p["Cd"]), float(p["u_tidal"]), float(p["gammaT"])
    sh = rd.shelf
    us = rd.field(_SHI, "SHIuStar", record)
    gm = rd.field(_SHI, "SHIgammT", record)
    floor = np.sqrt(cd) * ut                                        # AD16 Eq. 27, |u| = 0
    floor_ok = bool(us[sh].min() >= floor * (1 - 1e-12))
    den = np.where(us[sh] > 0, gT * us[sh], 1.0)
    rel_g = np.where(us[sh] > 0, np.abs(gm[sh] - gT * us[sh]) / den, np.inf)
    g_ok = bool(rel_g.max() < 1e-10)
    vsq = boundary_layer_velocity_sq(rd, rd.field(_DYN, "UVEL", record),
                                     rd.field(_DYN, "VVEL", record))
    bound = np.sqrt(cd * (vsq + ut * ut))
    ratio = us[sh] / bound[sh]
    j_ok = bool(ratio.min() >= 1 - rtol)
    return {"pass": floor_ok and g_ok and j_ok, "record": int(record),
            "n_shelf_columns": int(sh.sum()),
            "floor": {"pass": floor_ok, "min_ustar_m_s": float(us[sh].min()), "u_tidal": ut,
                      "floor_m_s": float(floor)},
            "gammT_identity": {"pass": g_ok, "max_rel_err": float(rel_g.max())},
            "jensen_bound": {"pass": j_ok, "rtol": rtol, "min_ratio": float(ratio.min()),
                             "median_ratio": float(np.median(ratio)),
                             "frac_below_1": float((ratio < 1).mean())},
            "equality": "not available: no boundary-layer U/V diagnostic and no velocity "
                        "snapshot simultaneous with SHIuStar (no model output added in Task 3)"}


# ---- R10 (post hoc): Gade line ---------------------------------------------------------------

def gade_line_check(rd: RunData, record: int | None, tol_T: float = 0.05,
                    tol_S: float = 0.05, ambient: str = "init") -> dict:
    """POST HOC (ruling R10; not a pre-registered criterion, S3 stays FAIL). For wet cells
    whose T or S lies outside the initial-profile range +- 0.1 (pd1 S3), test whether
    (T, S) is a mixture of the local ambient profile and meltwater, given C = TRAC01.

    Derivation for this model's virtual fluxes (verified against AD16 Eq. 24-30 and
    config/ocean0/code/shelfice_thermodynamics.F, conserve branch lines 700-724, and
    shelfice_forcing_mwtr.F):

    * salt: SHIForcS = (gS - q/rho0)(Sb - Sw), q = SHIfwFlx; with Eq. 26,
      gS (Sb - Sw) = q Sb / rho0, so SHIForcS = q Sw / rho0 = -(rho_fw/rho0) m Sw
      (= AD16 Eq. 30);
    * heat: SHIForcT = (gT - q/rho0)(Tb - Tw); with Eq. 24 (kappa_i = 0),
      gT (Tb - Tw) = q L / (rho0 cw), so SHIForcT = -(rho_fw/rho0) m (Tw - T_eff),
      T_eff = Tb - L/cw (= AD16 Eq. 29; Gade 1979 effective meltwater temperature with
      kappa_i = 0 [verificar: Gade 1979 not in the local library]);
    * tracer: source m [m s^-1] = -q / rho_fw, spread with the same weights.

    So per unit C the cell gains -(rho_fw/rho0) Sw of salt and -(rho_fw/rho0)(Tw - T_eff)
    of temperature; to first order in C:
        S_mix = (1 - r C) S_amb,  T_mix = (1 - r C) T_amb + r C T_eff,  r = rho_fw / rho0,
    with T_eff from the AD16 liquidus at the cell depth (p = rho0 g |z|) and S_amb.
    The controller's form (r = 1, the volume-flux mixing line) is reported alongside.
    ``record`` None uses the end-of-run state (T, S, PTRACER01); otherwise the monthly
    means (THETA, SALT, TRAC01), for which the linear mixing relation also holds."""
    from exp2_isomip_cavity.gendata import make_case as mc
    p = rd.params
    if record is None:
        n_end = int(p["nIter0"]) + int(p["nTimeSteps"])
        sd = rd.stage_dirs[-1]
        T, S, C = (_read_state(sd, n, n_end) for n in ("T", "S", "PTRACER01"))
        src = f"end-of-run state, iteration {n_end}"
    else:
        T, S, C = (rd.field(_DYN, n, record) for n in ("THETA", "SALT", "TRAC01"))
        src = f"monthly-mean record {record}"
    prof = p["init_profile"] if ambient == "init" else p["restore_profile"]
    zz = np.linspace(0.0, float(rd.rf[-1]), 7201)         # 0 to the bottom, as pd1 S3
    Tp, Sp, *_ = mc._profile(p["init_profile"], zz, p)
    lo_T, hi_T = Tp.min() - TOL_RANGE, Tp.max() + TOL_RANGE
    lo_S, hi_S = Sp.min() - TOL_RANGE, Sp.max() + TOL_RANGE
    Ta, Sa, *_ = mc._profile(prof, rd.rc, p)
    Ta3 = np.broadcast_to(np.asarray(Ta)[:, None, None], T.shape)
    Sa3 = np.broadcast_to(np.asarray(Sa)[:, None, None], T.shape)
    rho0 = float(p["rho_ref"])
    pres = rho0 * float(p["gravity"]) * (-rd.rc)
    Tf = ml.LAMBDA1 * np.asarray(Sa) + ml.LAMBDA2 + ml.LAMBDA3 * pres      # AD16 Eq. 25
    Teff3 = np.broadcast_to((Tf - ml.LATENT_HEAT / ml.CW)[:, None, None], T.shape)
    wet = rd.wet3
    out = wet & ((T < lo_T) | (T > hi_T) | (S < lo_S) | (S > hi_S))
    res = {"post_hoc": True, "label": "POST HOC (R10): Gade-line diagnostic, not a criterion",
           "state": src, "ambient_profile": prof, "T_bounds": [lo_T, hi_T],
           "S_bounds": [lo_S, hi_S], "tol_T": tol_T, "tol_S": tol_S,
           "n_wet": int(wet.sum()), "n_out_of_range": int(out.sum())}
    idx = np.nonzero(out)
    zc = rd.rc[idx[0]]
    for label, r in (("model_virtual_flux", RHO_FW / rho0), ("controller_r1", 1.0)):
        Tm = (1 - r * C) * Ta3 + r * C * Teff3
        Sm = (1 - r * C) * Sa3
        dT, dS = (T - Tm)[out], (S - Sm)[out]
        ok = (np.abs(dT) <= tol_T) & (np.abs(dS) <= tol_S)
        blk = {"r": r, "n_within_tol": int(ok.sum()),
               "fraction_within_tol": float(ok.mean()) if ok.size else None}
        for nm, d in (("T", dT), ("S", dS)):
            if d.size:
                w = int(np.argmax(np.abs(d)))
                blk[f"worst_{nm}"] = {
                    "residual": float(d[w]), "k": int(idx[0][w]), "j": int(idx[1][w]),
                    "i": int(idx[2][w]), "z_m": float(zc[w]),
                    "x_km": float(rd.xc[idx[1][w], idx[2][w]] / 1e3),
                    "y_km": float(rd.yc[idx[1][w], idx[2][w]] / 1e3),
                    "C": float(C[out][w]), "T": float(T[out][w]), "S": float(S[out][w])}
                blk[f"abs_residual_{nm}_p50_p95_p99"] = [float(np.percentile(np.abs(d), q))
                                                         for q in (50, 95, 99)]
        res[label] = blk
    # Extension (post hoc, beyond R10's local-ambient form): the ambient water need not come
    # from the cell's own depth (meltwater-laden water rises along the ice base). Invert the
    # salt relation for the ambient salinity, S_amb = S / (1 - r C), find the profile depth
    # z0 with S_a(z0) = S_amb, and predict T from T_a(z0). For the Ocean0 and ramp profiles
    # T_a and S_a are both linear in the same r(z), so the ambient T-S set is a straight
    # segment and any mixture of ambient waters from several depths lies on it: this test
    # is then exact to first order in r C (T_eff still taken at the cell depth; its
    # variation, < 0.5 degC, enters multiplied by C).
    r = RHO_FW / rho0
    order = np.argsort(Sp)
    Sq, Tq = Sp[order], Tp[order]
    s_amb = (S / (1 - r * C))[out]
    inside = (s_amb >= Sq[0] - tol_S) & (s_amb <= Sq[-1] + tol_S)
    t_amb = np.interp(np.clip(s_amb, Sq[0], Sq[-1]), Sq, Tq)
    t_mix = (1 - r * C[out]) * t_amb + r * C[out] * Teff3[out]
    dT = T[out] - t_mix
    ok = inside & (np.abs(dT) <= tol_T)
    blk = {"r": r, "n_within_tol": int(ok.sum()),
           "fraction_within_tol": float(ok.mean()) if ok.size else None,
           "n_S_amb_outside_profile": int((~inside).sum()),
           "note": "extension: ambient source depth inferred from S (not R10's local form)"}
    if dT.size:
        w = int(np.argmax(np.abs(dT)))
        blk["worst_T"] = {"residual": float(dT[w]), "k": int(idx[0][w]), "j": int(idx[1][w]),
                          "i": int(idx[2][w]), "z_m": float(zc[w]), "C": float(C[out][w]),
                          "S_amb": float(s_amb[w])}
        blk["abs_residual_T_p50_p95_p99"] = [float(np.percentile(np.abs(dT), q))
                                             for q in (50, 95, 99)]
    res["inferred_source_depth"] = blk
    m = res["model_virtual_flux"]
    res["n_within_tol"] = m["n_within_tol"]
    res["fraction_within_tol"] = m["fraction_within_tol"]
    res["worst"] = {"T": m.get("worst_T"), "S": m.get("worst_S")}
    if res["n_out_of_range"]:
        cc = C[out]
        res["C_out_of_range_min_median_max"] = [float(cc.min()), float(np.median(cc)),
                                                float(cc.max())]
    return res


# ---- summary --------------------------------------------------------------------------------

UNITS = {"melt": "m a-1 (freshwater, positive = melting)", "flux": "kg s-1 (positive = melting)",
         "psi": "Sv", "z_out": "m (negative down)", "D_out": "m (positive depth)",
         "N1": "1 (relative)", "N2": "1 (relative)", "T_f": "degC"}


def _q(value, units, source, **kw):
    return {"value": value, "units": units, "source": source, **kw}


def isomip_summary(rd: RunData, write: bool = True) -> dict:
    """Every pre-registered quantity computable for this run (pd0), with units and the record
    indices used, written to <run_dir>/metrics.json. The window is 'last6' for a 24-record
    run, 'year20' for a 240-record Ocean1 run; any other length uses all records and says
    that it is not a pre-registered window."""
    n = rd.n_records
    if n == N_RECORDS_2YR:
        wname, win = "last6", list(calendar_window(rd, "last6"))
    elif n == N_RECORDS_OCEAN1:
        wname, win = "year20", list(calendar_window(rd, "year20"))
    else:
        wname, win = "all", list(range(n))
    out = {"run_dir": str(rd.run_dir), "case": rd.params.get("case"),
           "config": rd.params.get("config"),
           "case_actual": {k: rd.actual.get(k) for k in ("Tbot", "z_c", "d", "s", "gammaT")},
           "n_records": n, "record_iters": [r.iter for r in rd.records],
           "window": {"name": wname, "records_0based": win,
                      "records_1based": [k + 1 for k in win],
                      "pre_registered": wname != "all",
                      "note": {"last6": "pd0: last 6 months of a 2-year run",
                               "year20": "pd0: Ocean1 year 20",
                               "all": "not a pre-registered window (run is neither 24 nor "
                                      "240 monthly records)"}[wname]},
           "units": UNITS}
    out["mean_melt_shelf"] = _q(mean_melt(rd, win, "shelf"), UNITS["melt"],
                                "pd0 Melt/Plataforma; V2a; H1 <m>")
    deep = rd.shelf & (rd.draft < DEEP_DRAFT_M)
    out["mean_melt_deep300"] = _q(mean_melt(rd, win, "deep300") if deep.any() else None,
                                  UNITS["melt"], "pd0 K1; AD16 Eq. 37")
    out["total_melt_flux"] = _q(total_melt_flux(rd, win), UNITS["flux"], "sum -SHIfwFlx rA")
    out["psi_bt_max_cavity"] = _q(psi_bt_max(rd, win, "cavity"), "Sv",
                                  "pd0 V2c (y > 50 km, x < 640 km; Y26 L1118)")
    out["psi_bt_max_front"] = _q(psi_bt_max(rd, win, "front"), "Sv",
                                 "pd0 V2d (630 < x < 650 km; Y26 L1155)")
    out["overturning_max"] = _q(overturning_max(rd, win), "Sv",
                                "pd0 V2b: max |psi(x,z)|, x < 640 km (Y26 L1100-1102)")
    melt_m = mean_melt_per_record(rd)
    cs_bt = [circulation_strength(rd, k, "bt") for k in range(n)]
    cs_moc = [circulation_strength(rd, k, "moc") for k in range(n)]
    out["circulation_strength"] = {
        "units": "Sv", "source": "pd0; Y26 Fig. 13 (L1365-1369): max - min psi, x < 640 km, "
                                 "per month",
        "window_mean": {"melt_shelf": float(np.mean([melt_m[k] for k in win])),
                        "bt": float(np.mean([cs_bt[k] for k in win])),
                        "moc": float(np.mean([cs_moc[k] for k in win]))},
        "monthly": {"melt_shelf": melt_m, "bt": cs_bt, "moc": cs_moc}}
    try:
        z = outflow_depth(rd, win)
        zw = outflow_depth(rd, win, face="west")
        out["z_out"] = _q(z, UNITS["z_out"], "pd0 H2 (x_front centre, u C > 0)",
                          west_face_informative=zw)
        out["D_out"] = _q(-z, UNITS["D_out"], "pd0 H2: D_out = -z_out",
                          west_face_informative=-zw)
    except RuntimeError as e:
        out["z_out"] = _q(None, UNITS["z_out"], "pd0 H2", reason=str(e))
        out["D_out"] = _q(None, UNITS["D_out"], "pd0 H2", reason=str(e))
    try:
        out["N1"] = _q(equilibrium_drift(rd), UNITS["N1"], "pd0 N1",
                       form="2yr" if n == N_RECORDS_2YR else "ocean1")
        out["N1"]["pass"] = bool(out["N1"]["value"] < 0.03)
    except ValueError as e:
        out["N1"] = _q(None, UNITS["N1"], "pd0 N1", reason=str(e))
    b = tracer_budget(rd)
    out["N2"] = {**b, "units": UNITS["N2"], "source": "pd0 N2"}
    out["T_f_GL"] = {**freezing_point_gl(rd), "units": UNITS["T_f"]}
    out["ustar_check_last_record (informative)"] = ustar_check(rd, n - 1)
    if write:
        (rd.run_dir / "metrics.json").write_text(json.dumps(out, indent=2, default=float))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="ISOMIP+ metrics (pd0) for one Exp. 2 run")
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--gade", action="store_true", help="post hoc Gade-line check (R10)")
    ap.add_argument("--ustar", action="store_true", help="u* consistency check")
    ap.add_argument("--record", type=int, default=None,
                    help="record (0-based) for --gade/--ustar; --gade default: end state")
    ap.add_argument("--json", type=Path, default=None, help="also write the result here")
    a = ap.parse_args(argv)
    rd = load_run(a.run_dir)
    if a.gade:
        res = gade_line_check(rd, a.record)
    elif a.ustar:
        res = ustar_check(rd, rd.n_records - 1 if a.record is None else a.record)
    else:
        res = isomip_summary(rd)
    txt = json.dumps(res, indent=2, default=float)
    print(txt)
    if a.json:
        a.json.write_text(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
