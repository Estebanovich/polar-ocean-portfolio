"""Plume diagnostics from MITgcm output of the 2D line-plume configuration (a).

Conventions
-----------
* Arrays are in MITgcm vertical order: index k = 0 is the surface cell. ``z``
  is the cell-centre elevation [m], negative downward.
* The ice face is the west side of column i = 0 (it carries R_icefront); the
  "face column" is therefore i = 0, whose cells are the 5 m cells touching the ice.
* Melt sign: pkg/icefront's ICFfwFlx [kg m-2 s-1] is NEGATIVE for melting
  (Task 2, verified in the MITgcm source), so ``m_ice = -ICFfwFlx / rho_i``
  [m s-1 of ice], positive = melting. This is the only place the conversion
  is done (`melt_rate_m_per_s`).
* Final-state fields: the diagnostics snapshots are written one step late, so
  there is no snapshot at t = endTime. State dumps (T, S, U, W at
  dumpInitAndLast) are used as additional instantaneous snapshots when their
  time lies in the window; the iteration-0 state is never used.
* Time windows are in hours. Snapshot (instantaneous) files are used when
  t0 <= t <= t1; hourly-average files (ICFfwFlx) labelled t cover (t-1 h, t] and
  are used when t0 <= t-1 h and t <= t1.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from pfcommon.mitgcm_io import read_bin, read_mds
from pfcommon.runctl import require_ok

RHO_I = 917.0  # kg m-3, icefront_readparms.F:86
WIDTH_FRAC = 0.10  # plume edge: w >= 10 % of the level's w_max
DETACH_FRAC = 0.05  # detachment: face-column w falls below 5 % of its maximum
AVG_PERIOD_S = 3600.0  # frequency of ICFfwFlx averages in data.diagnostics
_DUMP = {"THETA": "T", "SALT": "S", "UVEL": "U", "WVEL": "W"}
_EPS = 1e-6


def melt_rate_m_per_s(fw_flux_kg_m2_s: np.ndarray, rho_i: float = RHO_I) -> np.ndarray:
    """ICFfwFlx [kg m-2 s-1, < 0 = melt] -> ice melt rate [m s-1], positive = melting."""
    return -np.asarray(fw_flux_kg_m2_s, dtype=float) / rho_i


def detach_depth(z: np.ndarray, w: np.ndarray, frac: float = DETACH_FRAC) -> float:
    """Elevation where the rising plume ends.

    Starting at the level of maximum w and moving upward (towards the surface),
    the first level where w < frac * w_max; the crossing is linearly
    interpolated. ``z`` surface-first (decreasing), same layout as ``w``.
    Returns the surface elevation z[0] if the plume never decays, NaN if w_max <= 0.
    """
    z = np.asarray(z, dtype=float)
    w = np.asarray(w, dtype=float)
    kmax = int(np.nanargmax(w))
    wmax = w[kmax]
    if not wmax > 0:
        return float("nan")
    thr = frac * wmax
    for k in range(kmax - 1, -1, -1):
        if w[k] < thr:
            f = (w[k + 1] - thr) / (w[k + 1] - w[k])  # 0 at k+1, 1 at k
            return float(z[k + 1] + f * (z[k] - z[k + 1]))
    return float(z[0])


def _iterations(run_dir: Path, name: str) -> list[int]:
    out = []
    for f in run_dir.glob(f"{name}.[0-9]*.meta"):
        m = re.fullmatch(rf"{re.escape(name)}\.(\d+)\.meta", f.name)
        if m:
            out.append(int(m.group(1)))
    return sorted(out)


def _snapshots(run_dir, field, dt, t0, t1):
    """{iteration: 2D array (nz, nx)} of instantaneous `field` with t in [t0, t1] hours."""
    its = {}
    for name in (field, _DUMP[field]):
        for it in _iterations(run_dir, name):
            t = it * dt / 3600.0
            if it > 0 and t0 - _EPS <= t <= t1 + _EPS and it not in its:
                its[it] = read_mds(run_dir / name, it)[0][0]
    return its


def plume_profile(run_dir: Path, t_window_h: tuple[float, float]) -> dict:
    """Plume profiles in the face column, averaged over a time window.

    Raises RuntimeError if the run is not OK (Review Focus 3) and ValueError if
    the window contains no output.

    Returned keys (arrays surface-first, length nz unless noted):
      z, zf (nz+1 faces), x, xf        grid [m]
      w_face   w in the face column, interpolated to cell centres [m s-1]
      w_max    max over x of w at each level
      width    plume width [m]: from the ice face to the outer edge of the
               contiguous run of cells around argmax_x(w) with w >= 10 % of
               w_max(z); NaN at/above z_detach
      w_mean   mean of w over [0, width] ; Q = w_mean * width [m2 s-1]
      T, S, u_face   face-column THETA, SALT, UVEL (UVEL is on the west face of
               cell i, i.e. u = 0 at a wall: ~0 in column 0 by construction)
      melt     face melt rate [m s-1 of ice], positive = melting
      z_detach plume detachment elevation from the face-column w at its native
               (top-face) levels (see `detach_depth`)
      t_snapshots_h, t_means_h   hours used
    """
    run_dir = Path(run_dir)
    require_ok(run_dir)
    t0, t1 = map(float, t_window_h)
    params = json.loads((run_dir / "case.json").read_text())["params"]
    dt, nx, nz = float(params["deltaT"]), int(params["nx"]), int(params["nz"])
    dx = read_bin(run_dir / "delX.bin", (nx,))
    dz = read_bin(run_dir / "delZ.bin", (nz,))
    xf = np.concatenate([[0.0], np.cumsum(dx)])
    zf = -np.concatenate([[0.0], np.cumsum(dz)])
    x, z = 0.5 * (xf[1:] + xf[:-1]), 0.5 * (zf[1:] + zf[:-1])

    fields, times = {}, {}
    for f in _DUMP:
        snaps = _snapshots(run_dir, f, dt, t0, t1)
        if not snaps:
            raise ValueError(f"{run_dir}: no {f} snapshot in window {t_window_h} h")
        fields[f] = np.mean([snaps[i] for i in sorted(snaps)], axis=0)
        times[f] = [round(i * dt / 3600.0, 6) for i in sorted(snaps)]

    fw = []
    t_means = []
    for it in _iterations(run_dir, "icefrontDiag"):
        t = it * dt / 3600.0
        if t0 - _EPS <= t - AVG_PERIOD_S / 3600.0 and t <= t1 + _EPS:
            a, meta = read_mds(run_dir / "icefrontDiag", it)
            fw.append(a[meta["fldList"].index("ICFfwFlx")])
            t_means.append(round(t, 6))
    if not fw:
        raise ValueError(f"{run_dir}: no ICFfwFlx average inside window {t_window_h} h")
    melt = melt_rate_m_per_s(np.mean(fw, axis=0))[:, 0]

    W = fields["WVEL"]  # on the top face of each cell; bottom face is the sea floor
    Wc = 0.5 * (W + np.vstack([W[1:], np.zeros((1, nx))]))
    w_face = Wc[:, 0]
    z_detach = detach_depth(zf[:-1], W[:, 0])  # native w levels (top faces): no cell-centre smoothing

    w_max = Wc.max(axis=1)
    width = np.full(nz, np.nan)
    w_mean = np.full(nz, np.nan)
    for k in range(nz):
        if not (z[k] < z_detach) or w_max[k] <= 0:
            continue
        thr = WIDTH_FRAC * w_max[k]
        i = int(np.argmax(Wc[k]))
        while i + 1 < nx and Wc[k, i + 1] >= thr:
            i += 1
        width[k] = xf[i + 1]
        w_mean[k] = float(np.sum(Wc[k, : i + 1] * dx[: i + 1]) / width[k])

    return {
        "z": z, "zf": zf, "x": x, "xf": xf,
        "w_face": w_face, "w_max": w_max, "width": width, "w_mean": w_mean,
        "Q": w_mean * width,
        "T": fields["THETA"][:, 0], "S": fields["SALT"][:, 0], "u_face": fields["UVEL"][:, 0],
        "melt": melt, "z_detach": z_detach,
        "t_snapshots_h": times["WVEL"], "t_means_h": t_means,
    }
