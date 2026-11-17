"""Coarse (chan10, G0) hydrography -> fine grid (chan2p5) for the F0 spin-up, and the small
initial-condition perturbations that separate ensemble members (spec Sec. 4; plan Task 4).

Velocities are not interpolated: the fine model starts at rest with the G0 hydrography."""
from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .grid import Grid

ALPHA, BETA = 3.0e-5, 7.4e-4          # linear EOS of config/*/input/data


def _fill_land(F: np.ndarray, wet: np.ndarray) -> np.ndarray:
    """Replace land values (MITgcm writes 0 in land) so that interpolation never mixes them in:
    first from the nearest wet row in y on the same level (this covers the wall rows, which
    are land from the surface down), then from the value above."""
    F, filled = F.copy(), wet.copy()
    for k in range(F.shape[0]):
        while not filled[k].all() and filled[k].any():
            for shift in (1, -1):
                src = np.roll(filled[k], shift, axis=0)
                src[0 if shift == 1 else -1] = False          # no wrap in y (closed channel)
                take = src & ~filled[k]
                F[k][take] = np.roll(F[k], shift, axis=0)[take]
                filled[k] |= take
    for k in range(1, F.shape[0]):
        F[k] = np.where(filled[k], F[k], F[k - 1])
    return F


def _regrid(F: np.ndarray, gc: Grid, gf: Grid) -> np.ndarray:
    """Bilinear in (y, x) on each level; periodic in x, constant beyond the outer y centres."""
    xs = np.concatenate([[gc.xc[-1] - gc.nx * gc.dx_m], gc.xc, [gc.xc[0] + gc.nx * gc.dx_m]])
    yq = np.clip(gf.yc, gc.yc[0], gc.yc[-1])
    Y, X = np.meshgrid(yq, gf.xc, indexing="ij")
    out = np.empty((F.shape[0], gf.ny, gf.nx))
    for k in range(F.shape[0]):
        Fk = np.concatenate([F[k][:, -1:], F[k], F[k][:, :1]], axis=1)
        out[k] = RegularGridInterpolator((gc.yc, xs), Fk)((Y, X))
    return out


def _stabilize(T: np.ndarray, S: np.ndarray, max_pass: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """Convective adjustment: mix (equal-thickness weighted) any statically unstable pair."""
    T, S = T.copy(), S.copy()
    for _ in range(max_pass):
        rho = -ALPHA * T + BETA * S
        bad = np.diff(rho, axis=0) < 0
        if not bad.any():
            break
        for k in np.unique(np.nonzero(bad)[0]):
            m = bad[k]
            tm, sm = 0.5 * (T[k][m] + T[k + 1][m]), 0.5 * (S[k][m] + S[k + 1][m])
            T[k][m], T[k + 1][m], S[k][m], S[k + 1][m] = tm, tm, sm, sm
    return T, S


def coarse_to_fine(theta_c, salt_c, grid_c: Grid, grid_f: Grid) -> tuple[np.ndarray, np.ndarray]:
    """Fine-grid theta and salinity (nr, ny_f, nx_f) from the coarse state (nr, ny_c, nx_c)."""
    wet = np.asarray(salt_c) != 0.0
    T = _regrid(_fill_land(np.asarray(theta_c, float), wet), grid_c, grid_f)
    S = _regrid(_fill_land(np.asarray(salt_c, float), wet), grid_c, grid_f)
    return _stabilize(T, S)


def perturb_member(theta: np.ndarray, k: int, amp: float = 1e-3, seed: int = 0) -> np.ndarray:
    """theta plus uniform noise in [-amp, amp] (degC), reproducible per (seed, member k)."""
    rng = np.random.default_rng([seed, k])
    return np.asarray(theta, float) + rng.uniform(-amp, amp, size=np.shape(theta))


_PICKUP_2D = {"EtaN", "dEtaHdt", "EtaH"}


def pickup_ts(data_path, nr: int, ny: int, nx: int) -> tuple[np.ndarray, np.ndarray]:
    """Theta and Salt (nr, ny, nx) from a MITgcm pickup (pickup.<it>.data + .meta). The fields
    are stacked in the order of fldList, nr records for 3-D fields and one for EtaN, dEtaHdt,
    EtaH (model/src/write_pickup.F)."""
    import re
    from pathlib import Path
    data_path = Path(data_path)
    meta = data_path.with_suffix(".meta").read_text()
    names = re.findall(r"'(\w+)\s*'", meta[meta.index("fldList"):])
    raw = np.fromfile(data_path, dtype=">f8")
    out, off = {}, 0
    for n in names:
        size = (1 if n in _PICKUP_2D else nr) * ny * nx
        out[n] = raw[off:off + size]
        off += size
    if off != raw.size:
        raise ValueError(f"{data_path}: {raw.size} values, fldList accounts for {off}")
    return out["Theta"].reshape(nr, ny, nx), out["Salt"].reshape(nr, ny, nx)


def last_pickup(run_dir):
    """Newest pickup.<it>.data under run_dir (input copies excluded)."""
    from pathlib import Path
    hits = [h for h in Path(run_dir).rglob("pickup.[0-9]*.data")
            if not h.parent.name.startswith("input")]
    if not hits:
        raise FileNotFoundError(f"no pickup under {run_dir}")
    return max(hits, key=lambda f: int(f.name.split(".")[1]))


def perturb_pickup_theta(data_path, nr: int, ny: int, nx: int, k: int, amp: float = 1e-3,
                         seed: int = 0) -> None:
    """Add perturb_member noise to Theta of a pickup, in place, in wet cells only (Salt != 0),
    for ensemble member k."""
    import re
    from pathlib import Path
    data_path = Path(data_path)
    meta = data_path.with_suffix(".meta").read_text()
    names = re.findall(r"'(\w+)\s*'", meta[meta.index("fldList"):])
    raw = np.fromfile(data_path, dtype=">f8")
    off, where = 0, {}
    for n in names:
        size = (1 if n in _PICKUP_2D else nr) * ny * nx
        where[n] = slice(off, off + size)
        off += size
    T = raw[where["Theta"]].reshape(nr, ny, nx)
    wet = raw[where["Salt"]].reshape(nr, ny, nx) != 0
    Tp = np.where(wet, perturb_member(T, k, amp, seed), T)
    raw[where["Theta"]] = Tp.ravel()
    raw.astype(">f8").tofile(data_path)
