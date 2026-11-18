"""Biogeochemical inputs (spec Sec. 3.3-3.4).

Reference profiles: PTRACERS_ref of verification/tutorial_global_oce_biogeo/input/data.ptracers
on the 15 levels of its data (delZ, line 67-69), interpolated linearly in depth, constant beyond
the end levels. Silica: mean of tutorial_global_oce_biogeo/input/sillev1.bin (12 monthly
records, 128 x 64, 2.8125 deg, from 90 S) over wet points south of 50 S. PAR: daily-mean
top-of-atmosphere insolation (solar constant 1361 W m-2; declination of Cooper 1969,
delta = 23.45 deg sin(2 pi (284 + d)/365); sunset hour angle h0 = arccos(-tan phi tan delta))
times a PAR fraction, at mid-month days of a 365-day year."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

BIOGEO = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification/tutorial_global_oce_biogeo/input")
NAMES = ("DIC", "Alk", "PO4", "DOP", "O2")
S0 = 1361.0


def _floats(block: str) -> list[float]:
    out = []
    for tok in re.findall(r"(\d+)\*([-+0-9.Ee]+)|([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)", block):
        if tok[0]:
            out += [float(tok[1])] * int(tok[0])
        else:
            out.append(float(tok[2]))
    return out


def source_levels() -> np.ndarray:
    txt = (BIOGEO / "data").read_text()
    m = re.search(r"delZ\s*=\s*(.*?)\n\s*[a-zA-Z]", txt, re.S)
    dz = np.array(_floats(m.group(1))[:15])
    rf = np.insert(np.cumsum(dz), 0, 0.0)
    return 0.5 * (rf[1:] + rf[:-1])


def _refs() -> dict[str, np.ndarray]:
    txt = (BIOGEO / "data.ptracers").read_text()
    out = {}
    for i, name in enumerate(NAMES, start=1):
        m = re.search(rf"PTRACERS_ref\(1:15,{i}\)\s*=\s*(.*?)(?:\n\s*[#A-Za-z&])", txt, re.S)
        v = _floats(m.group(1))
        assert len(v) == 15, (name, len(v))
        out[name] = np.array(v)
    return out


def ref_profiles(z_m) -> dict[str, np.ndarray]:
    zc, refs = source_levels(), _refs()
    return {k: np.interp(np.asarray(z_m, float), zc, v) for k, v in refs.items()}


def par_monthly(lat_deg: float = -55.0, par_frac: float = 0.4) -> np.ndarray:
    phi = np.deg2rad(lat_deg)
    mid = np.cumsum([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) - np.array(
        [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) / 2.0
    delta = np.deg2rad(23.45) * np.sin(2 * np.pi * (284 + mid) / 365.0)
    h0 = np.arccos(np.clip(-np.tan(phi) * np.tan(delta), -1.0, 1.0))
    q = S0 / np.pi * (h0 * np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.sin(h0))
    return par_frac * q


def anth_sponge(dic_pi_z, cant_sup_by_year: dict[int, float], z, efold_m: float = 1000.0) -> np.ndarray:
    """(251, nr): DIC_pi(z) + C_ant,sup(Y) exp(-z/efold_m), Y = 1850..2100 (spec Sec. 3.4)."""
    p = np.exp(-np.asarray(z, float) / efold_m)
    return np.array([np.asarray(dic_pi_z, float) + cant_sup_by_year[y] * p for y in range(1850, 2101)])


def silica_surface() -> float:
    a = np.fromfile(BIOGEO / "sillev1.bin", ">f4").reshape(12, 64, 128)
    lat = -90.0 + 2.8125 * (np.arange(64) + 0.5)
    s = a[:, lat <= -50.0, :]
    return float(s[s > 0].mean())
