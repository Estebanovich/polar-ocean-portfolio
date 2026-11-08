"""D1 (post-hoc, pre-registered in notes/pd1_b_twoplume3d.md): plume volume flux Q_m(z) near the face
vs half-cone theory Q_t(z), same ambient as C4 (window-mean far field), window 2-3 h."""
import sys, glob, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pfcommon.paths import EXP1_RUNS
from pfcommon.mitgcm_io import read_mds
from pfcommon.runctl import require_ok
from exp1_fjord_plumes.analysis.plume_diag import ambient_far_from_face
from exp1_fjord_plumes.analysis.fig2_two_plume_melt import AMBIENT_X_M
from exp1_fjord_plumes.theory.plume import point_plume
from exp1_fjord_plumes.theory.melt import MeltParams

R = EXP1_RUNS
out = {}
for s in ("twolayer", "sermilik"):
    run = R / f"b_3d_{s}_single"
    require_ok(run)
    meta = json.loads((run / "case.json").read_text()); p = meta["params"]
    H, dx, dy, dz = p["H"], p["dx"], p["dy"], p["dz"]
    its = sorted(int(f.split(".")[-2]) for f in glob.glob(str(run / "WVEL.*.meta")))
    its = [i for i in its if 7200 <= i * p["deltaT"] <= 10800]  # snapshots in the 2-3 h window
    W = np.mean([read_mds(run / "WVEL", i)[0][0] for i in its], axis=0)  # (nz, ny, nx), w at cell top
    assert np.isfinite(W).all()
    hf = np.fromfile(run / "hFacC.data", ">f8").reshape(W.shape)
    wpos = np.clip(W[:, :, :10], 0, None) * (hf[:, :, :10] > 0)
    Qm = wpos.sum(axis=(1, 2)) * dx * dy               # m3/s through the top face of each level, i < 10
    zf = -np.arange(W.shape[0]) * dz                   # elevation of the top face of level k (k=0 surface)
    amb = ambient_far_from_face(run, (2, 3), AMBIENT_X_M)
    zg = np.linspace(-H, 0.0, int(H) + 1)
    zc, Tm, Sm = (np.asarray(a, float)[::-1] for a in (amb["z"], amb["T"], amb["S"]))
    sol = point_plume(zg, np.interp(zg, zc, Tm), np.interp(zg, zc, Sm), p["Q0"], p["T0_src"], p["S0_src"],
                      MeltParams(ice_depth=H), "mitgcm", alpha=p.get("alpha", 0.1), geometry="half_cone")
    Qt = np.interp(zf, sol.z, np.nan_to_num(sol.Q, nan=np.nan), left=np.nan, right=np.nan)
    rows = []
    for k in range(W.shape[0] - 1, 0, -1):  # bottom to top, skip the surface face (w = 0)
        rows.append((float(zf[k]), float(Qm[k]), float(Qt[k]) if np.isfinite(Qt[k]) else None))
    lower = [r for r in rows if r[0] <= -H / 2 and r[2]]
    ratio_lower = float(np.median([r[1] / r[2] for r in lower])) if lower else None
    out[s] = {"z_detach_theory": float(sol.z_detach), "rows_z_Qm_Qt": rows, "median_Qm_over_Qt_lower_half": ratio_lower,
              "n_snapshots": len(its)}
    print(f"{s}: theory z_detach {sol.z_detach:.1f} m; median Qm/Qt (z <= -H/2) = {ratio_lower:.3f}; n_snap {len(its)}")
    for z, qm, qt in rows[::3]:
        print(f"   z={z:7.1f}  Qm={qm:8.2f}  Qt={'%8.2f' % qt if qt else '   (none)'}")
(Path(__file__).resolve().parents[1] / "figures" / "d1_volume_flux.json").write_text(json.dumps(out, indent=1))
