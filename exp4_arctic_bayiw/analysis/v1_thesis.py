"""V1 (PD0 section 4): does R1 reproduce the stored thesis run BayIW_Rectan/run_expand?

Compares the R1 run (MDS diagnostics snapshots every 900 s: state2d = ETAN, state3d = UVEL,
VVEL, THETA) with the thesis MNC state files (Eta, U, V, Temp every 900 s, 8 tiles of
140 x 176) at the same iterations, and the model pickups byte for byte (pickup.0000007200 at
the end of stage 1; pickup.ckptA at iteration 14400, end of stage 2).

- V1-exacto (frozen): max |d eta| <= 1e-10 m over every common record.
- V1-espectral (frozen; evaluated in any case, descriptive if V1-exacto passes): eta at the mouth
  station, Hann-windowed periodogram of 0-120 h (df = 1/120 cph): (i) the bay-seiche peak in the
  same frequency bin; (ii) PSD at the peak within 10 %.
- f0 (Task 3 Step 6): f = f0 + beta yC (model/src/ini_cori.F:77) with beta = 1e-11 (MITgcm
  default, model/src/set_defaults.F:112; the thesis data does not set beta); the near-inertial
  peak of the surface kinetic energy at stations on the shelf is compared with the local f.

Usage::

    v1_thesis.py --run /Volumes/LaCie/portfolio_runs/exp4/R1 [--out figures/v1_thesis.json]
                 [--every-3d 12]

``--every-3d k`` compares U, V, T at every k-th common record (all records of eta are always
compared). The run must pass pfcommon.runctl.require_ok (each stage).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parents[1]
REPO = PKG.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pfcommon.mitgcm_io import read_mds  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402

THESIS_RUN = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand")
NX, NY, NZ, SNX, SNY, NPX = 560, 352, 50, 140, 176, 4
TOL_ETA = 1.0e-10                      # m, PD0 section 4, V1-exacto
DF_TOL_BINS = 0                        # same bin, PD0 V1-espectral (i)
PSD_TOL = 0.10                         # PD0 V1-espectral (ii)
SEICHE_BAND_CPH = (0.18, 0.35)         # thesis picos_seiche_inercial.BANDA_SEICHE
NI_BAND_CPH = (0.02, 0.07)             # thesis picos_seiche_inercial.BANDA_NI
# thesis stations (j, i), picos_seiche_inercial.ESTACIONES; the coast is j = 323, so Seiche2
# (j = 325, i = 280) is the station just inside the bay mouth, on the bay axis.
STATIONS = {"Seiche": (345, 280), "Seiche2": (325, 280), "Omega": (315, 280),
            "A": (315, 230), "B": (315, 330), "C": (265, 280)}
MOUTH_STATION = "Seiche2"
# extra shelf stations along the bay axis i = 280 for the f(y) test (post hoc, Step 6)
F_STATIONS = {"S300": (300, 280), "S200": (200, 280), "S120": (120, 280), "S60": (60, 280)}


# ---- thesis MNC -------------------------------------------------------------------------------

class ThesisMNC:
    """Index of the thesis state.<iter>.t<tile>.nc files: (var, iteration) -> global field."""

    def __init__(self, run_dir: Path = THESIS_RUN):
        import netCDF4  # only needed for the thesis MNC files
        self.files: dict[int, list[tuple[np.ndarray, Path]]] = {}
        for f in sorted(Path(run_dir).glob("mnc_*/state.*.t*.nc")):
            tile = int(re.search(r"\.t(\d+)\.nc$", f.name).group(1))
            self.files.setdefault(tile, []).append((None, f))
        if sorted(self.files) != list(range(1, 9)):
            raise RuntimeError(f"thesis run: expected tiles 1..8, found {sorted(self.files)}")
        for tile, lst in self.files.items():
            out = []
            for _, f in lst:
                with netCDF4.Dataset(f) as d:
                    out.append((np.asarray(d["iter"][:]), f))
            self.files[tile] = out

    def iterations(self) -> np.ndarray:
        its = [set(np.concatenate([a for a, _ in lst]).tolist()) for lst in self.files.values()]
        return np.array(sorted(set.intersection(*its)))

    def _where(self, tile: int, it: int):
        for a, f in self.files[tile]:
            k = np.flatnonzero(a == it)
            if k.size:
                return f, int(k[0])
        raise KeyError(f"tile {tile}: iteration {it} not stored")

    def field(self, var: str, it: int, level: int | slice | None = None) -> np.ndarray:
        """Global (ny, nx) or (nz, ny, nx) field; U/V on the MITgcm west/south faces."""
        import netCDF4
        out = None
        for tile in range(1, 9):
            f, k = self._where(tile, it)
            ti, tj = (tile - 1) % NPX, (tile - 1) // NPX
            with netCDF4.Dataset(f) as d:
                v = d[var]
                a = np.asarray(v[k] if level is None else v[k, level])
            if var == "U":
                a = a[..., :SNX]
            elif var == "V":
                a = a[..., :SNY, :]
            if out is None:
                out = np.empty(a.shape[:-2] + (NY, NX))
            out[..., tj * SNY:(tj + 1) * SNY, ti * SNX:(ti + 1) * SNX] = a
        return out


# ---- R1 MDS -----------------------------------------------------------------------------------

def stage_dirs(run: Path) -> list[Path]:
    run = Path(run)
    return [run] + sorted((p for p in run.glob("stage[0-9]*") if p.is_dir()),
                          key=lambda p: int(p.name[5:]))


def mds_iterations(run: Path, prefix: str) -> dict[int, Path]:
    out = {}
    for d in stage_dirs(run):
        for m in d.glob(f"{prefix}.[0-9]*.meta"):
            out.setdefault(int(m.name.split(".")[1]), d)    # the first stage that wrote it
    return dict(sorted(out.items()))


def mds_field(d: Path, prefix: str, it: int, rec: int = 0, level: int | None = None) -> np.ndarray:
    """One record of an MDS file (memory map: only the requested record/level is read)."""
    meta_txt = (d / f"{prefix}.{it:010d}.meta").read_text()
    nrec = int(re.search(r"nrecords\s*=\s*\[\s*(\d+)", meta_txt).group(1))    # one per field
    nd = int(re.search(r"nDims\s*=\s*\[\s*(\d+)", meta_txt).group(1))
    shape = (NZ, NY, NX) if nd == 3 else (NY, NX)
    m = np.memmap(d / f"{prefix}.{it:010d}.data", dtype=">f8", mode="r", shape=(nrec,) + shape)
    a = m[rec] if level is None else m[rec, level]
    return np.array(a, dtype=np.float64)


# ---- spectra ----------------------------------------------------------------------------------

def psd_hann(x: np.ndarray, dt_s: float):
    """One-sided periodogram with a Hann window (single window), units^2 / cph; f in cph."""
    x = np.asarray(x, float) - np.mean(x)
    n = x.size
    w = np.hanning(n)
    X = np.fft.rfft(x * w)
    fs = 3600.0 / dt_s                                  # samples per hour
    P = (np.abs(X) ** 2) / (fs * np.sum(w ** 2))
    P[1:-1 if n % 2 == 0 else None] *= 2
    f = np.fft.rfftfreq(n, d=1.0 / fs)                  # cph
    return f, P, 1.0 / (n * dt_s / 3600.0)


def band_peak(f, P, band):
    m = (f >= band[0]) & (f <= band[1])
    k = np.flatnonzero(m)[np.argmax(P[m])]
    return int(k), float(f[k]), float(P[k])


# ---- main -------------------------------------------------------------------------------------

def compare(run: Path, every_3d: int = 12, thesis_run: Path = THESIS_RUN) -> dict:
    for d in stage_dirs(run):
        require_ok(d)
    th = ThesisMNC(thesis_run)
    ti = th.iterations()
    eta_its = mds_iterations(run, "state2d")
    s3_its = mds_iterations(run, "state3d")
    common = [int(i) for i in ti if int(i) in eta_its]
    if not common:
        raise RuntimeError("no common iterations between R1 and the thesis run")
    out: dict = {"run": str(run), "thesis_run": str(thesis_run), "n_common_eta_records": len(common),
                 "iterations": [common[0], common[-1]], "tolerance_eta_m": TOL_ETA}

    # eta: every common record, full domain; plus station series for the spectra
    deta, series_r1, series_th = [], {}, {}
    sts = {**STATIONS, **F_STATIONS}
    for it in common:
        a = mds_field(eta_its[it], "state2d", it)
        b = th.field("Eta", it)
        deta.append(float(np.max(np.abs(a - b))))
        for s, (j, i) in sts.items():
            series_r1.setdefault(s, []).append(a[j, i])
            series_th.setdefault(s, []).append(b[j, i])
    deta = np.array(deta)
    out["eta_max_abs_diff_m"] = float(deta.max())
    out["eta_max_abs_diff_at_iter"] = int(common[int(np.argmax(deta))])
    out["eta_n_records_bitwise_equal"] = int(np.sum(deta == 0.0))
    out["V1_exacto_pass"] = bool(deta.max() <= TOL_ETA)

    # U, V, T on a subset of the common records (all levels, full domain)
    sub = [it for it in common[::every_3d] if it in s3_its] + \
          ([common[-1]] if common[-1] in s3_its and common[-1] not in common[::every_3d] else [])
    d3 = {"U": [], "V": [], "T": []}
    for it in sub:
        for rec, (key, var) in enumerate((("U", "U"), ("V", "V"), ("T", "Temp"))):
            a = mds_field(s3_its[it], "state3d", it, rec=rec)
            b = th.field(var, it)
            d3[key].append(float(np.max(np.abs(a - b))))
    out["uvt_records_compared"] = sub
    out["max_abs_diff"] = {"U_m_s-1": max(d3["U"]), "V_m_s-1": max(d3["V"]), "T_degC": max(d3["T"])}

    # pickups, byte for byte
    pk = {}
    for name, ours in (("pickup.0000007200.data", run / "pickup.0000007200.data"),
                       ("pickup.ckptA.data", run / "stage2" / "pickup.ckptA.data")):
        pk[name] = pickup_check(ours, Path(thesis_run) / name)
    out["pickups"] = pk

    # V1-espectral at the mouth station and spectra at all stations
    dt = 900.0
    spec = {}
    for s in sts:
        f, P1, df = psd_hann(np.array(series_r1[s]), dt)
        _, P2, _ = psd_hann(np.array(series_th[s]), dt)
        k1, f1, p1 = band_peak(f, P1, SEICHE_BAND_CPH)
        k2, f2, p2 = band_peak(f, P2, SEICHE_BAND_CPH)
        spec[s] = {"seiche_peak_cph": [f1, f2], "seiche_peak_bin": [k1, k2],
                   "psd_at_peak_m2_cph-1": [p1, p2],
                   "psd_rel_diff": float(abs(p1 - p2) / p2) if p2 > 0 else None}
    out["df_cph"] = float(df)
    out["spectra_eta"] = spec
    m = spec[MOUTH_STATION]
    out["V1_espectral"] = {"station": MOUTH_STATION, "jj_ii": STATIONS[MOUTH_STATION],
                           "same_bin": m["seiche_peak_bin"][0] == m["seiche_peak_bin"][1],
                           "psd_rel_diff": m["psd_rel_diff"],
                           "pass": bool(m["seiche_peak_bin"][0] == m["seiche_peak_bin"][1]
                                        and m["psd_rel_diff"] is not None and m["psd_rel_diff"] <= PSD_TOL)}
    return out


def pickup_check(ours: Path, theirs: Path) -> dict:
    """Byte-for-byte comparison of two pickups. An absent file is recorded as status "missing"
    (bitwise_equal None), never skipped."""
    if not (ours.exists() and theirs.exists()):
        return {"bitwise_equal": None, "status": "missing",
                "missing": [str(p) for p in (ours, theirs) if not p.exists()]}
    r = {"bitwise_equal": ours.read_bytes() == theirs.read_bytes(), "status": "compared",
         "iter_ours": _meta_iter(ours), "iter_thesis": _meta_iter(theirs)}
    if not r["bitwise_equal"]:
        x = np.fromfile(ours, ">f8")
        y = np.fromfile(theirs, ">f8")
        r["max_abs_diff"] = float(np.max(np.abs(x - y)))
    return r


def _meta_iter(data_file: Path):
    m = Path(str(data_file)[: -len(".data")] + ".meta")
    if not m.exists():
        return None
    r = re.search(r"timeStepNumber\s*=\s*\[\s*(\d+)", m.read_text())
    return int(r.group(1)) if r else None


def inertial_test(run: Path, f0: float, beta: float, yg0: float, dely: np.ndarray) -> dict:
    """Near-inertial peak of the surface KE spectrum (u, v of level 1, 0-120 h) at stations
    along i = 280, against the local f = f0 + beta yC (ini_cori.F:77) and against f0."""
    its = mds_iterations(run, "state3d")
    keys = sorted(its)
    yc = yg0 + np.cumsum(dely) - 0.5 * dely
    sts = {**STATIONS, **F_STATIONS}
    u = {s: [] for s in sts}
    v = {s: [] for s in sts}
    for it in keys:
        a = mds_field(its[it], "state3d", it, rec=0, level=0)
        b = mds_field(its[it], "state3d", it, rec=1, level=0)
        for s, (j, i) in sts.items():
            u[s].append(0.5 * (a[j, i] + a[j, i + 1]))
            v[s].append(0.5 * (b[j, i] + b[j + 1, i]))
    res = {}
    for s, (j, i) in sts.items():
        f, Pu, df = psd_hann(np.array(u[s]), 900.0)
        _, Pv, _ = psd_hann(np.array(v[s]), 900.0)
        k, fp, _ = band_peak(f, Pu + Pv, NI_BAND_CPH)
        f_loc = (f0 + beta * yc[j]) * 3600 / (2 * np.pi)
        res[s] = {"j": j, "i": i, "y_km": float(yc[j] / 1e3), "ke_peak_cph": fp,
                  "f_local_cph": float(f_loc), "f0_cph": float(f0 * 3600 / (2 * np.pi)),
                  "peak_minus_flocal_bins": float((fp - f_loc) / df),
                  "peak_minus_f0_bins": float((fp - f0 * 3600 / (2 * np.pi)) / df)}
    return {"df_cph": float(df), "stations": res}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=PKG / "figures" / "v1_thesis.json")
    ap.add_argument("--every-3d", type=int, default=12)
    a = ap.parse_args(argv)
    res = compare(a.run, a.every_3d)
    cj = json.loads((a.run / "input" / "case.json").read_text())
    from exp4_arctic_bayiw.gendata.make_case import parse_namelist_list
    txt = (a.run / "input" / "data").read_text()
    dely = np.fromfile(a.run / "input" / "delY.bin", ">f8")
    res["inertial"] = inertial_test(a.run, cj["actual"]["f0_s-1"], cj["actual"]["beta_m-1s-1"],
                                    parse_namelist_list(txt, "ygOrigin")[0], dely)
    a.out.write_text(json.dumps(res, indent=2))
    print(json.dumps({k: res[k] for k in ("eta_max_abs_diff_m", "V1_exacto_pass", "max_abs_diff",
                                          "pickups", "V1_espectral")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
