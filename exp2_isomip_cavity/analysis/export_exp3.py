#!/opt/anaconda3/envs/pyocean/bin/python
r"""Export the Exp. 2 ensemble as the training / test dataset of Exp. 3 (melt emulator).

Usage::

    export_exp3.py [--runs-root /Volumes/LaCie/portfolio_runs/exp2] [--run-suffix _r21]
                   [--out /Volumes/LaCie/portfolio_runs/exp2/exp3_dataset.npz]
                   [--summary exp2_isomip_cavity/figures/exp3_dataset.json]

Members: the ensemble of cases.yaml (``ensemble.registry_cases``: the reference ocean0_4km
and the 28 ``ens_*`` runs, 4 km, Gamma_T,cal,4km). The post hoc ``hol_*`` runs (R24) are not
ensemble members and are not exported.

* A run that fails N3 (``metrics.load_run`` refuses it: not OK, incomplete, NaN, void melt
  law) or that does not exist is **excluded**, never imputed; it is listed in
  ``excluded_case`` / ``excluded_reason`` (pd0 N3, R22).
* Runs failing N1 stay in the dataset, flagged by ``N1`` / ``N1_pass`` (pd0 N1, R22).
* ``split`` = 'test' for Tbot = 2.0 (the extrapolation group), 'train' otherwise
  (``ensemble.split_of``).

Arrays (n = exported runs, model grid ny x nx incl. the land row, nz levels):

    case, group, split            (n,) str
    Tbot, z_c, d, s               (n,) float   profile / geometry parameters (degC, m, m, -)
    N1, N1_pass                   (n,) float, bool
    mean_melt_shelf               (n,) float   whole-shelf mean melt, m a^-1 (H1 <m>, R22)
    z                             (nz,) float  cell-centre heights RC, m (negative down)
    T_profile, S_profile          (n, nz)      initial (= restoring) profile on z (degC, PSU)
    x, y                          (nx,), (ny,) cell centres, m
    melt                          (n, ny, nx)  last-6-month (records 19-24) mean melt rate,
                                               m a^-1 freshwater, + = melting; NaN on dry cells
    draft                         (n, ny, nx)  ice base (icetopo.bin), m; NaN on dry cells
    wet_mask, shelf_mask          (n, ny, nx)  bool; shelf = wet & draft < 0 (pd0)
    excluded_case, excluded_reason (k,) str

The .npz is written byte-reproducibly (fixed zip timestamps, no pickle), so its sha256 is a
fingerprint of the content; the summary JSON records it.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import zipfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp2_isomip_cavity.analysis import ensemble as en  # noqa: E402
from exp2_isomip_cavity.analysis import metrics as mt  # noqa: E402

DEFAULT_OUT = en.DEFAULT_RUNS / "exp3_dataset.npz"
DEFAULT_SUMMARY = REPO / "exp2_isomip_cavity" / "figures" / "exp3_dataset.json"
ZIP_DATE = (1980, 1, 1, 0, 0, 0)          # fixed member timestamp: byte-reproducible archive


def _profile(rd: mt.RunData) -> tuple:
    """(T, S, Tbot, z_c, d) of the run's initial profile on the model levels (case.json)."""
    from exp2_isomip_cavity.gendata.make_case import _profile as profile_of
    p = rd.params
    return profile_of(p["init_profile"], rd.rc, p)


def _member(rd: mt.RunData) -> dict:
    win = list(mt.calendar_window(rd, "last6"))                 # pd0: records 19-24
    wet, shelf = rd.wet, rd.shelf
    melt = mt.melt_rate_m_per_a(rd.mean_field("shiDiag", "SHIfwFlx", win))
    T, S, Tbot, z_c, d = _profile(rd)
    n1 = float(mt.equilibrium_drift(rd, "2yr"))
    return {"melt": np.where(wet, melt, np.nan), "draft": np.where(wet, rd.draft, np.nan),
            "wet_mask": wet, "shelf_mask": shelf, "T_profile": T, "S_profile": S,
            "Tbot": Tbot, "z_c": z_c, "d": d, "s": float(rd.params.get("s", 1.0)),
            "N1": n1, "N1_pass": bool(n1 < en.N1_MAX),
            "mean_melt_shelf": float(mt.mean_melt(rd, win, "shelf")),
            "z": rd.rc, "x": rd.xc[0], "y": rd.yc[:, 0]}


def build_dataset(runs_root, suffix: str = "_r21", cases: list[dict] | None = None) -> dict:
    """Dataset dict of numpy arrays (module doc). Raises if no run can be exported or if the
    exported runs do not share one grid."""
    cases = en.registry_cases() if cases is None else cases
    rows, excluded = [], []
    for spec in cases:
        run_dir = Path(runs_root) / f"{spec['case']}{suffix}"
        if not run_dir.exists():
            excluded.append((spec["case"], "GAP_MISSING"))
            continue
        try:
            rd = mt.load_run(run_dir)                           # N3 gate, before any read
        except Exception as e:                                  # noqa: BLE001 - N3 gap
            excluded.append((spec["case"], f"GAP_N3: {type(e).__name__}: {e}"[:300]))
            continue
        row = _member(rd)
        row.update(case=spec["case"], group=str(spec["group"]), split=en.split_of(row["Tbot"]))
        rows.append(row)
    if not rows:
        raise ValueError(f"no run of {len(cases)} could be exported from {runs_root}")
    for key in ("z", "x", "y"):
        if any(not np.array_equal(r[key], rows[0][key]) for r in rows):
            raise ValueError(f"exported runs do not share one grid ({key} differs)")
    ds = {"case": np.array([r["case"] for r in rows], dtype=str),
          "group": np.array([r["group"] for r in rows], dtype=str),
          "split": np.array([r["split"] for r in rows], dtype=str)}
    for key in ("Tbot", "z_c", "d", "s", "N1", "mean_melt_shelf"):
        ds[key] = np.array([r[key] for r in rows], dtype=float)
    ds["N1_pass"] = np.array([r["N1_pass"] for r in rows], dtype=bool)
    for key in ("z", "x", "y"):
        ds[key] = np.asarray(rows[0][key], dtype=float)
    for key in ("T_profile", "S_profile", "melt", "draft"):
        ds[key] = np.stack([np.asarray(r[key], dtype=float) for r in rows])
    for key in ("wet_mask", "shelf_mask"):
        ds[key] = np.stack([np.asarray(r[key], dtype=bool) for r in rows])
    ds["excluded_case"] = np.array([c for c, _ in excluded], dtype=str)
    ds["excluded_reason"] = np.array([r for _, r in excluded], dtype=str)
    return ds


def write_npz(ds: dict, path) -> str:
    """Write ``ds`` as an uncompressed .npz (np.load-compatible, no pickle) with fixed member
    timestamps, so identical content gives identical bytes. Returns the file's sha256."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as zf:
        for key in sorted(ds):
            buf = io.BytesIO()
            np.lib.format.write_array(buf, np.asarray(ds[key]), allow_pickle=False)
            info = zipfile.ZipInfo(f"{key}.npy", date_time=ZIP_DATE)
            info.external_attr = 0o644 << 16
            zf.writestr(info, buf.getvalue())
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summary(ds: dict, path: Path, sha: str) -> dict:
    """Small JSON record of the export (tracked in the repo; the .npz is not)."""
    split = list(ds["split"])
    return {"file": str(path), "sha256": sha, "bytes": Path(path).stat().st_size,
            "window": "last 6 months, records 19-24 (pd0)",
            "n_runs": len(split), "n_train": split.count("train"), "n_test": split.count("test"),
            "test_cases": [c for c, s in zip(ds["case"], split) if s == "test"],
            "n_N1_fail": int((~ds["N1_pass"]).sum()),
            "N1_fail_cases": [str(c) for c, p in zip(ds["case"], ds["N1_pass"]) if not p],
            "excluded": [{"case": str(c), "reason": str(r)}
                         for c, r in zip(ds["excluded_case"], ds["excluded_reason"])],
            "grid": {"nx": int(ds["x"].size), "ny": int(ds["y"].size), "nz": int(ds["z"].size)},
            "arrays": {k: {"dtype": str(np.asarray(v).dtype), "shape": list(np.shape(v))}
                       for k, v in sorted(ds.items())}}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs-root", type=Path, default=en.DEFAULT_RUNS)
    ap.add_argument("--run-suffix", default="_r21")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    a = ap.parse_args(argv)
    ds = build_dataset(a.runs_root, a.run_suffix)
    sha = write_npz(ds, a.out)
    rec = summary(ds, a.out, sha)
    a.summary.parent.mkdir(parents=True, exist_ok=True)
    a.summary.write_text(json.dumps(rec, indent=1) + "\n")
    print(f"wrote {a.out} ({rec['n_runs']} runs: {rec['n_train']} train, {rec['n_test']} test; "
          f"{len(rec['excluded'])} excluded; {rec['n_N1_fail']} N1-flagged) sha256 {sha}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
