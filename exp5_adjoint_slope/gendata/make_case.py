"""Write the MITgcm inputs of one Exp. 5 case (plan Task 1).

`write_case(case, out_dir)` merges cases.yaml (configs + case overrides), writes the binaries
under fixed names and copies config/<cfg>/input/* when that directory exists. case.json holds the
actual parameters, including the J-line index j_qb (0-based; Fortran row j_qb + 1)."""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import yaml

from pfcommon.mitgcm_io import write_bin

from . import fields
from .grid import Grid, vertical_spacing

PKG = Path(__file__).resolve().parents[1]


def load_case(case: str) -> dict:
    reg = yaml.safe_load((PKG / "cases.yaml").read_text())
    over = dict(reg["cases"][case])
    cfg = over.pop("config")
    p = dict(reg["configs"][cfg])
    p.update(over)
    p["config"], p["case"] = cfg, case
    return p


def grid_of(p: dict) -> Grid:
    return Grid(int(p["nx"]), int(p["ny"]), float(p["dx_m"]),
                vertical_spacing(int(p["nr"]), float(p["dz_top_m"]), float(p["depth_total_m"])))


def write_case(case: str, out_dir: Path, overrides: dict | None = None) -> Path:
    p = load_case(case)
    p.update(overrides or {})
    g = grid_of(p)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name in p.get("input_dirs", [p["config"]]):      # later directories override earlier
        src = PKG / "config" / name / "input"
        if src.is_dir():
            for f in src.iterdir():
                if f.is_file():
                    shutil.copy2(f, out / f.name)
    T, S = fields.initial_ts(g)
    T0, S0 = T, S                                       # sponge target: initial profiles
    init_pickup = None
    if p.get("init_from"):                              # start from a coarse spin-up (F0)
        from pfcommon import paths
        from . import interp
        src = load_case(p["init_from"])
        gc = grid_of(src)
        init_pickup = interp.last_pickup(paths.EXP5_RUNS / p["init_from"])
        Tc, Sc = interp.pickup_ts(init_pickup, gc.nr, gc.ny, gc.nx)
        T, S = interp.coarse_to_fine(Tc, Sc, gc, g)
    sst, sss, lam_t, lam_s = fields.surface_relaxation(g)
    tau = np.repeat(fields.wind_tau_x(g.yc)[:, None], g.nx, axis=1)
    if p.get("tau_pert"):                               # FP members (H5b): coarse dtau_x(y) profile
        pat = json.loads(Path(p["tau_pert"]).read_text())
        tau = tau + np.interp(g.yc, pat["y_m"], pat["dtau"])[:, None]
    gm = np.repeat(fields.gm_factor(g.yc)[:, None], g.nx, axis=1)
    sponge = fields.sponge_mask(g)
    files = {"bathy.bin": fields.bathymetry(g), "T.init.bin": T, "S.init.bin": S,
             "tau_x.bin": tau, "sst_relax.bin": sst, "sss_relax.bin": sss,
             "lambda_t.bin": lam_t, "lambda_s.bin": lam_s,
             "salt_flux.bin": fields.polynya_salt_flux(g),
             "rbcs_mask.bin": sponge, "rbcs_mask_uv.bin": sponge, "rbcs_T.bin": T0, "rbcs_S.bin": S0,
             "gm_map.bin": gm, "delZ.bin": g.dz_m}
    if p["config"].endswith("_tap"):                      # control weights (data.ctrl)
        files["ones_64b.bin"] = np.ones((g.nr, g.ny, g.nx))
    for name, arr in files.items():
        write_bin(out / name, arr)
    jq = fields.j_qb(g, float(p["y_qb_m"]))
    meta = {**p, "params": {"nx": g.nx, "ny": g.ny, "nz": g.nr}, "j_qb": jq, "y_qb_depth_m": float(fields.depth(np.array([g.yg[jq]]))[0]),
            "init_pickup": str(init_pickup) if init_pickup else None}
    (out / "case.json").write_text(json.dumps(meta, indent=1, sort_keys=True))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("case")
    ap.add_argument("out_dir", type=Path)
    a = ap.parse_args(argv)
    print(write_case(a.case, a.out_dir))


if __name__ == "__main__":
    main()
