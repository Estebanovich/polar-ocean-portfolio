"""POST HOC check (Task 7 fix round 1): the Kang & Fringer baroclinic pressure work through the
straits with p' built from rho' relative to an EVOLVING reference, the box-mean density
rho_box(z, t), instead of the fixed initial rho_b of PD0 Enmienda 1.

Why: p~ = rho0 g eta + g int rho' dz (Nash et al. 2005 eq. 3 [Pa]) uses the same fixed rho_b
whose drift inflates the exact APE (notes/pd3_hipotesis.md), so a drift of the background
profile adds int p'_drift(z) Q'(z) dz to F_bc wherever there is a mean baroclinic exchange Q'.
Method (fixed in notes/pd3_hipotesis.md, commit 2af0db1, BEFORE computing):

1. rho_box(z_k, t) = volume-weighted (dz rA) mean of rho over the wet cells of the box at level k
   and record t; the fixed rho_b where the box has no wet cell at level k;
2. rho'_m = rho - rho_box(z, t) in the wet cells;
3. p~_m = rho0 g eta + g sum rho'_m dz (budget.pressure_tilde), p'_m = p~_m minus its depth mean
   on the face (KF12 eqs. 15-16); u' as always;
4. F_bc,m = sum u_n' p'_m dz [W] per face, summed per strait; record integral (trapezoid, + out),
   band split (validate.band_flux, the same bands) and eta_m = F_bc,m,out / (P_bt + P_bc) with
   the denominator of budget.json;
5. F1 and F1c over 0-279.67 h, R1 and R1k over 0-120 h; the box-mean profile at the first and the
   last record of F1, F1c, F2, F3, R1, R1k and R2.

Descriptive only: no frozen verdict changes. Output: figures/posthoc_moving_reference.json.

Usage: posthoc_reference.py [--out figures/posthoc_moving_reference.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp4_arctic_bayiw.analysis import budget, validate as va  # noqa: E402

PKG = Path(__file__).resolve().parents[1]
LABEL = ("POST HOC check (Task 7 fix round 1; method fixed in notes/pd3_hipotesis.md before "
         "computing): pressure work with an evolving box-mean reference rho_box(z, t); descriptive, "
         "never a verdict")


def _progress(what, r, n, on=True, every=120):
    if on and (r % every == 0 or r == n - 1):
        print(f"[posthoc] {what}: {r + 1}/{n}", file=sys.stderr, flush=True)


def box_mean_profile(ds, mask, rho_b, time=slice(None), progress=True) -> np.ndarray:
    """rho_box(z_k, t) [kg m^-3] (nt, nz): volume-weighted (drF hFacC rA) mean of RHO over the wet
    cells of ``mask`` at each level and record; rho_b where the box has no wet cell."""
    m = np.asarray(mask, bool)
    jj, ii = np.nonzero(m)
    j0, j1, i0, i1 = jj.min(), jj.max() + 1, ii.min(), ii.max() + 1
    hC = ds["hFacC"].values
    nz = hC.shape[0]
    vol = ds["drF"].values[:, None] * hC[:, jj, ii] * ds["rA"].values[jj, ii][None]     # (nz, nc)
    W = vol.sum(1)
    nk = int(np.flatnonzero(W > 0).max()) + 1
    tidx = np.arange(ds.sizes["time"])[time]
    out = np.tile(np.asarray(rho_b, float), (tidx.size, 1))
    for q, n in enumerate(tidx):
        R = ds["RHO"].isel(time=int(n), k=slice(0, nk), j=slice(j0, j1), i=slice(i0, i1)).values
        v = R[:, jj - j0, ii - i0]
        ok = W[:nk] > 0
        out[q, :nk][ok] = (v * vol[:nk]).sum(1)[ok] / W[:nk][ok]
        _progress("box mean profile", q, tidx.size, progress)
    return out


def face_pressure(ds, mask, face_groups: dict, rho_ref, time=slice(None), progress=True) -> dict:
    """u_n' and p'_m on the open faces of ``face_groups`` with rho' = rho - rho_ref[t] (rho_ref
    (nt, nz)), as budget.kf_budget forms the baroclinic pressure work otherwise. One rectangle
    holding every face (+1/+2 halo) is read per record. Returns up, pp (nt, nz, nf), dz, w, group,
    t (as partition.face_series)."""
    rho0, g = float(ds.attrs["rho0"]), float(ds.attrs.get("g", budget.G))
    drF = ds["drF"].values
    hC, hW, hS = ds["hFacC"].values, ds["hFacW"].values, ds["hFacS"].values
    dxG, dyG = ds["dxG"].values, ds["dyG"].values
    ny, nx = hC.shape[1:]
    faces = budget.box_faces(np.asarray(mask, bool), wet=None)
    sel = []
    for kind in ("x", "y"):
        for j, i, s in faces[kind]:
            hit = [n for n, gr in face_groups.items() if gr[kind][j, i]]
            if hit:
                sel.append((kind, j, i, s, hit[0]))
    if not sel:
        raise ValueError("no face of the box is in face_groups")
    tidx = np.arange(ds.sizes["time"])[time]
    rho_ref = np.asarray(rho_ref, float)
    if rho_ref.shape[0] != tidx.size:
        raise ValueError("rho_ref must have one profile per record")
    nt, nz, nf = tidx.size, drF.size, len(sel)
    up, pp = np.zeros((nt, nz, nf)), np.zeros((nt, nz, nf))
    dz, w = np.zeros((nz, nf)), np.zeros(nf)
    for q, (kind, j, i, s, _) in enumerate(sel):
        dz[:, q] = drF * (hS[:, j, i] if kind == "y" else hW[:, j, i])
        w[q] = s * (dxG[j, i] if kind == "y" else dyG[j, i])
    js = [x[1] for x in sel]
    is_ = [x[2] for x in sel]
    j0, j1 = max(min(js) - 1, 0), min(max(js) + 2, ny)
    i0, i1 = max(min(is_) - 1, 0), min(max(is_) + 2, nx)
    nk = int(np.flatnonzero((hC[:, j0:j1, i0:i1] > 0).any(axis=(1, 2))).max()) + 1
    dzC = drF[:nk, None, None] * hC[:nk, j0:j1, i0:i1]
    wetC = dzC > 0
    mW, mS = hW[:nk, j0:j1, i0:i1] > 0, hS[:nk, j0:j1, i0:i1] > 0
    for r, n in enumerate(tidx):
        win = dict(time=int(n), k=slice(0, nk), j=slice(j0, j1), i=slice(i0, i1))
        U = ds["U"].isel(**win).values * mW
        V = ds["V"].isel(**win).values * mS
        R = ds["RHO"].isel(**win).values
        E = ds["ETA"].isel(time=int(n), j=slice(j0, j1), i=slice(i0, i1)).values
        rp = np.where(wetC, R - rho_ref[r, :nk, None, None], 0.0)
        pt = budget.pressure_tilde(rp, E, dzC, rho0, g)                  # Pa, Nash 2005 eq. 3
        for q, (kind, j, i, _, _) in enumerate(sel):
            a, b = j - j0, i - i0
            if kind == "y":
                un, p = V[:, a, b], 0.5 * (pt[:, a - 1, b] + pt[:, a, b])
            else:
                un, p = U[:, a, b], 0.5 * (pt[:, a, b - 1] + pt[:, a, b])
            dzq = dz[:nk, q]
            H = dzq.sum()
            Hs = H if H > 0 else 1.0
            up[r, :nk, q] = un - (un * dzq).sum() / Hs                    # KF12: u'
            pp[r, :nk, q] = p - (p * dzq).sum() / Hs                      # KF12: p'
        _progress("face pressure", r, nt, progress)
    return {"up": up, "pp": pp, "dz": dz, "w": w, "group": np.array([x[4] for x in sel]),
            "t": ds["time"].values[tidx].astype(float)}


def moving_pressure_work(ds, mask, face_groups: dict, rho_b, time=slice(None), progress=True) -> dict:
    """The moving-reference pass: rho_box(z, t), then the face pressure work. series_moving
    {group: F_bc,m (nt,) [W], + INTO the box}."""
    ref = box_mean_profile(ds, mask, rho_b, time=time, progress=progress)
    fp = face_pressure(ds, mask, face_groups, ref, time=time, progress=progress)
    Fq = (fp["up"] * fp["pp"] * fp["dz"][None]).sum(1) * fp["w"][None]
    ser = {str(n): Fq[:, fp["group"] == n].sum(1) for n in np.unique(fp["group"])}
    return {"fp": fp, "rho_ref": ref, "series_moving": ser, "t": fp["t"]}


def summarize(r: dict, bands: dict, P_J: float) -> dict:
    """Per strait and total: record integral out [J] (trapezoid), dt-sum and band split out [J]
    (validate.band_flux), and eta = total out / P_J, eta_resolved = (bands other than the one
    starting at 0) / P_J."""
    fp, t = r["fp"], r["t"]
    if t[0] != 0.0:
        raise ValueError("the record must start at t = 0")
    dt, n = float(t[1] - t[0]), t.size
    low = min(bands, key=lambda b: bands[b][0])
    out = {"label": LABEL}
    names = list(r["series_moving"])
    for g in names + ["total"]:
        s = np.ones(fp["group"].size, bool) if g == "total" else fp["group"] == g
        F = sum(r["series_moving"][x] for x in names) if g == "total" else r["series_moving"][g]
        b = va.band_flux(fp["up"][:, :, s], fp["pp"][:, :, s], fp["dz"][:, s], fp["w"][s], dt, bands)
        out[g] = {"out_J": float(-budget.cumulative(t, F)[-1]), "dt_sum_out_J": float(-dt * F.sum()),
                  "band_out_J": {k: -v * n * dt for k, v in b.items()}}
    out["P_J"] = P_J
    out["eta"] = out["total"]["out_J"] / P_J
    out["eta_resolved"] = sum(v for k, v in out["total"]["band_out_J"].items() if k != low) / P_J
    out["low_band_out_J"] = out["total"]["band_out_J"][low]
    return out


def main(argv=None) -> int:
    from exp4_arctic_bayiw.analysis import load, partition as pa
    from pfcommon.paths import EXP4_RUNS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=PKG / "figures" / "posthoc_moving_reference.json")
    ap.add_argument("--runs-root", type=Path, default=EXP4_RUNS)
    ap.add_argument("--store-root", type=Path, default=None, help="local disk for the box store")
    a = ap.parse_args(argv)
    th = va.load_theory()
    bud = json.loads((PKG / "figures" / "budget.json").read_text())["runs"]
    res = json.loads(a.out.read_text()) if a.out.exists() else {}
    res["label"] = LABEL
    res.setdefault("runs", {})
    res.setdefault("profile_drift", {})
    jobs = (("F1", "F1", None), ("F1c", "F1c", None), ("R1", "R1", None), ("R1k_0_120h", "R1k", 481))
    cache = {}
    for name, run, n in jobs:
        if name in res["runs"]:
            continue
        ds = load.run(a.runs_root / run)
        mask, groups = pa._box(run, ds, a.runs_root)
        cache[run] = (ds, mask)
        rho_b = budget.background_profile(ds["RHO"].isel(time=0).values, ds["hFacC"].values > 0)
        tsl = slice(0, n) if n else slice(None)
        # one LaCie pass into the local box store (partition.build_store; identical results, tested)
        import tempfile
        st = pa.build_store(ds, mask, Path(a.store_root or tempfile.gettempdir()) / "exp4_posthoc" / run)
        dl, ml, gl, _ = pa.store_dataset(ds, st, mask, groups)
        r = moving_pressure_work(dl, ml, gl, rho_b, time=tsl)
        pa.remove_store(st)
        if n:
            assert abs(r["t"][-1] - 120 * 3600.0) < 1e-6
        s = summarize(r, pa.run_bands(name, th)["bands_hz"], bud[name]["P_J"])
        s["window_h"] = [0.0, float(r["t"][-1] / 3600.0)]
        s["n_records"] = int(r["t"].size)
        s["fixed_reference"] = {"F_bc_pressure_out_J": bud[name]["F_bc_pressure_out_J"],
                                "eta": bud[name]["F_bc_pressure_out_J"] / bud[name]["P_J"]}
        res["runs"][name] = s
        a.out.write_text(json.dumps(res, indent=1, default=va.jsonable) + "\n")
    for run in ("F1", "F1c", "F2", "F3", "R1", "R1k", "R2"):
        if run in res["profile_drift"]:
            continue
        ds, mask = cache.get(run, (None, None))
        if ds is None:
            ds = load.run(a.runs_root / run)
            mask, _ = pa._box(run, ds, a.runs_root)
        rho_b = budget.background_profile(ds["RHO"].isel(time=0).values, ds["hFacC"].values > 0)
        last = ds.sizes["time"] - 1
        prof = box_mean_profile(ds, mask, rho_b, time=[0, last], progress=False)
        nk = int(np.flatnonzero((ds["hFacC"].values[:, np.asarray(mask, bool)] > 0).any(1)).max()) + 1
        res["profile_drift"][run] = {
            "t_h": [0.0, float(ds["time"].values[last] / 3600.0)], "z_m": ds["RC"].values[:nk].tolist(),
            "rho_box_first_kg_m-3": prof[0, :nk].tolist(), "rho_box_last_kg_m-3": prof[1, :nk].tolist(),
            "max_abs_change_kg_m-3": float(np.abs(prof[1, :nk] - prof[0, :nk]).max())}
        a.out.write_text(json.dumps(res, indent=1, default=va.jsonable) + "\n")
    print(json.dumps({k: {kk: res["runs"][k][kk] for kk in ("eta", "eta_resolved", "low_band_out_J")}
                      for k in res["runs"]}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
