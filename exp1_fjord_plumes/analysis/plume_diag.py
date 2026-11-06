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
AMBIENT_X = 1000.0  # m, column whose window-mean T/S is the "evolved ambient" (post-hoc)
WIDTH_METHOD = ("10 % of w_max(z); outer edge by linear interpolation between cell centres "
                "(R22); w piecewise linear, constant between the face and the first centre")
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


def plume_extent(w: np.ndarray, x: np.ndarray, frac: float = WIDTH_FRAC,
                 x_end: float | None = None) -> tuple[float, float, float]:
    """Plume width, volume flux and momentum flux of one level (R22).

    ``w`` at cell centres ``x`` (x = 0 is the ice face). The plume is the
    contiguous run of cells around argmax(w) with w >= frac * w_max; its outer
    edge is located by LINEAR INTERPOLATION between the centres of the last cell
    inside and the first cell outside (the same 10 % definition as before, but
    not quantised to stretched-grid cell faces). w is piecewise linear between
    centres and constant (= w[0]) between the face and the first centre.
    Returns (width, Q = int_0^edge w dx, M = int_0^edge w^2 dx); each segment is
    integrated exactly. If the run reaches the last cell, the edge is ``x_end``
    (default x[-1]) with w constant there.
    """
    w = np.asarray(w, dtype=float)
    x = np.asarray(x, dtype=float)
    thr = frac * w.max()
    i = int(np.argmax(w))
    while i + 1 < w.size and w[i + 1] >= thr:
        i += 1
    xs = np.concatenate([[0.0], x[: i + 1]])
    ws = np.concatenate([[w[0]], w[: i + 1]])
    if i + 1 < w.size:
        edge = x[i] + (w[i] - thr) / (w[i] - w[i + 1]) * (x[i + 1] - x[i])
        xs, ws = np.append(xs, edge), np.append(ws, thr)
    else:
        edge = x[-1] if x_end is None else float(x_end)
        xs, ws = np.append(xs, edge), np.append(ws, w[i])
    d = np.diff(xs)
    a, b = ws[:-1], ws[1:]
    return float(edge), float(np.sum(0.5 * (a + b) * d)), float(np.sum((a * a + a * b + b * b) / 3.0 * d))


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
               w_max(z), edge interpolated linearly between cell centres
               (`plume_extent`, R22); NaN at/above z_detach
      w_mean   mean of w over [0, width] ; Q = w_mean * width [m2 s-1];
               M = int w^2 dx over the same width [m3 s-2] (M/Q = momentum-weighted
               mean velocity, post-hoc metric)
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
    mom = np.full(nz, np.nan)
    for k in range(nz):
        if not (z[k] < z_detach) or w_max[k] <= 0:
            continue
        width[k], Q, mom[k] = plume_extent(Wc[k], x, WIDTH_FRAC, x_end=xf[-1])
        w_mean[k] = Q / width[k]

    ic = int(np.argmin(np.abs(x - AMBIENT_X)))  # ambient column used for the evolved-ambient diagnostic
    return {
        "z": z, "zf": zf, "x": x, "xf": xf,
        "w_face": w_face, "w_max": w_max, "width": width, "w_mean": w_mean,
        "Q": w_mean * width, "M": mom,
        "x_amb": float(x[ic]), "T_amb": fields["THETA"][:, ic], "S_amb": fields["SALT"][:, ic],
        "width_method": WIDTH_METHOD,
        "T": fields["THETA"][:, 0], "S": fields["SALT"][:, 0], "u_face": fields["UVEL"][:, 0],
        "melt": melt, "z_detach": z_detach,
        "t_snapshots_h": times["WVEL"], "t_means_h": t_means,
    }


def _avg_period_s(run_dir: Path, it: int) -> float:
    """Averaging period [s] of a diagnostics file from ``timeInterval = [t_start t_end]``
    in its .meta (pkg/diagnostics writes the interval; (b) averages ICFfwFlx over 900 s,
    (a) over 3600 s). Raises ValueError if the meta has no timeInterval: guessing a period
    would silently change which averages enter a window."""
    text = (run_dir / f"icefrontDiag.{it:010d}.meta").read_text()
    m = re.search(r"timeInterval\s*=\s*\[\s*([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*\]", text)
    if not m:
        raise ValueError(f"{run_dir}: icefrontDiag.{it:010d}.meta has no timeInterval")
    return float(m.group(2)) - float(m.group(1))


def _require_finite(run_dir: Path, what: str, arr) -> None:
    """Raise ValueError naming the run if ``arr`` (a loaded diagnostic) has a non-finite value."""
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{Path(run_dir).name}: non-finite values in {what} inside the window "
                         "(diverged run?); no analysis is done on it")


def _grid3d(run_dir: Path, params: dict):
    """(dy[ny], dz[nz], hFacC(nz, ny, nx)) of a 3-D run; hFacC.data is written by the model
    (MDS, big-endian float64, array shape (nz, ny, nx))."""
    nx, ny, nz = int(params["nx"]), int(params["ny"]), int(params["nz"])
    dy = read_bin(run_dir / "delY.bin", (ny,))
    dz = read_bin(run_dir / "delZ.bin", (nz,))
    hfac = read_bin(run_dir / "hFacC.data", (nz, ny, nx))
    return dy, dz, hfac


def face_melt_field(run_dir: Path, t_window_h: tuple[float, float]) -> tuple[np.ndarray, list[float]]:
    """Window-mean ice melt rate on the face column i = 0 of a 3-D run.

    Returns (m_ice (nz, ny) [m s-1 of ice, positive = melting], t_means_h). Uses the
    ICFfwFlx averages whose whole averaging interval lies in the window
    (t_end - period >= t0 and t_end <= t1; period read from the .meta). ICFfwFlx is
    filled per unit ice-face area (pkg/icefront/icefront_thermodynamics.F:188, before the
    icefrontlength scaling at :156-:175), so m_ice = -ICFfwFlx / rho_i.
    """
    run_dir = Path(run_dir)
    require_ok(run_dir)
    t0, t1 = map(float, t_window_h)
    params = json.loads((run_dir / "case.json").read_text())["params"]
    dt = float(params["deltaT"])
    fw, t_means = [], []
    for it in _iterations(run_dir, "icefrontDiag"):
        t = it * dt / 3600.0
        if t0 - _EPS <= t - _avg_period_s(run_dir, it) / 3600.0 and t <= t1 + _EPS:
            a, meta = read_mds(run_dir / "icefrontDiag", it)
            fw.append(a[meta["fldList"].index("ICFfwFlx")][..., 0])  # face column i = 0
            t_means.append(round(t, 6))
    if not fw:
        raise ValueError(f"{run_dir}: no ICFfwFlx average inside window {t_window_h} h")
    _require_finite(run_dir, "ICFfwFlx", fw)
    return melt_rate_m_per_s(np.mean(fw, axis=0)), t_means


def face_total_melt(run_dir: Path, t_window_h: tuple[float, float]) -> float:
    """Total ice loss on the face, M = sum_{k,j} m_ice dy dz hFacC over i = 0 [m3 s-1 of ice],
    window-mean (pre-registered in notes/pd1_b_twoplume3d.md; the face area of cell (k, j) is
    dy_j * dz_k * hFacC, as icefrontlength = 1/dx scales the per-area flux to the cell,
    pkg/icefront/icefront_thermodynamics.F:156-175). Refuses runs that are not OK."""
    run_dir = Path(run_dir)
    m_ice, _ = face_melt_field(run_dir, t_window_h)
    params = json.loads((run_dir / "case.json").read_text())["params"]
    dy, dz, hfac = _grid3d(run_dir, params)
    area = dz[:, None] * dy[None, :] * hfac[:, :, 0]  # m2 of ice face per cell
    return float(np.sum(m_ice * area))


def face_w_profile(run_dir: Path, t_window_h: tuple[float, float]) -> dict:
    """Face-column vertical velocity of a 3-D run, for the C4 detachment depth.

    ``w_ymax(k)`` = max over y of the window-mean (over the WVEL snapshots / W dumps in the
    window, `_snapshots`) of w at i = 0, at the native (top-face) levels, so that
    ``detach_depth(zf[:-1], w_ymax)`` is the same estimator as in (a). Keys: z (centres),
    zf (faces), w_ymax, j_max (y index of the maximum at each level), z_detach, t_snapshots_h.
    """
    run_dir = Path(run_dir)
    require_ok(run_dir)
    t0, t1 = map(float, t_window_h)
    params = json.loads((run_dir / "case.json").read_text())["params"]
    dt, nz = float(params["deltaT"]), int(params["nz"])
    dz = read_bin(run_dir / "delZ.bin", (nz,))
    zf = -np.concatenate([[0.0], np.cumsum(dz)])
    snaps = _snapshots(run_dir, "WVEL", dt, t0, t1)
    if not snaps:
        raise ValueError(f"{run_dir}: no WVEL snapshot in window {t_window_h} h")
    _require_finite(run_dir, "WVEL", [snaps[i][..., 0] for i in snaps])
    Wm = np.mean([snaps[i][..., 0] for i in sorted(snaps)], axis=0)  # (nz, ny) at i = 0
    w_ymax = Wm.max(axis=1)
    return {"z": 0.5 * (zf[1:] + zf[:-1]), "zf": zf, "w_ymax": w_ymax, "j_max": Wm.argmax(axis=1),
            "z_detach": detach_depth(zf[:-1], w_ymax),
            "t_snapshots_h": [round(i * dt / 3600.0, 6) for i in sorted(snaps)]}


def ambient_far_from_face(run_dir: Path, t_window_h: tuple[float, float], x_m: float) -> dict:
    """Window-mean THETA/SALT of the column nearest x = ``x_m`` from the ice face, averaged over
    the wet y cells (hFacC > 0): the "evolved ambient far from the face" of the C4 theory
    (same snapshots as the other window means). Keys z (centres, surface first), T, S, x, i."""
    run_dir = Path(run_dir)
    require_ok(run_dir)
    t0, t1 = map(float, t_window_h)
    params = json.loads((run_dir / "case.json").read_text())["params"]
    dt, nx, nz = float(params["deltaT"]), int(params["nx"]), int(params["nz"])
    dx = read_bin(run_dir / "delX.bin", (nx,))
    dz = read_bin(run_dir / "delZ.bin", (nz,))
    _, _, hfac = _grid3d(run_dir, params)
    xc = np.cumsum(dx) - 0.5 * dx
    i = int(np.argmin(np.abs(xc - x_m)))
    wet = hfac[:, :, i] > 0  # (nz, ny)
    out = {"z": -(np.cumsum(dz) - 0.5 * dz), "x": float(xc[i]), "i": i}
    for f, key in (("THETA", "T"), ("SALT", "S")):
        snaps = _snapshots(run_dir, f, dt, t0, t1)
        if not snaps:
            raise ValueError(f"{run_dir}: no {f} snapshot in window {t_window_h} h")
        _require_finite(run_dir, f, [snaps[k][..., i] for k in snaps])
        col = np.mean([snaps[k][..., i] for k in sorted(snaps)], axis=0)  # (nz, ny)
        out[key] = np.array([col[k][wet[k]].mean() if wet[k].any() else np.nan for k in range(nz)])
    return out
