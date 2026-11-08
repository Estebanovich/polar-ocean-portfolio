"""Plume-excess metric with clean no-discharge background (pre-registered, pd1_b_twoplume3d.md)."""
import sys, glob, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pfcommon.paths import EXP1_RUNS
from pfcommon.mitgcm_io import read_mds
from pfcommon.runctl import require_ok
from exp1_fjord_plumes.analysis.plume_diag import face_total_melt
R = EXP1_RUNS; out = {}
for s in ("twolayer", "sermilik"):
    bg = R / f"b_bg_{s}_nodischarge"; require_ok(bg)
    Mb, Mb1 = face_total_melt(bg, (2, 3)), face_total_melt(bg, (1, 2))
    eps = abs(Mb - Mb1)
    Ms = face_total_melt(R / f"b_3d_{s}_single", (2, 3)); Es = Ms - Mb
    res = {"M_bg": Mb, "M_bg_1_2h": Mb1, "eps": eps, "M_single": Ms, "E_single": Es,
           "Mbg_over_Msingle": Mb / Ms, "resolvable": Es >= 10 * eps, "R_E": {}}
    for x in (0, 20, 40, 80, 160, 320):
        run = R / f"b_3d_{s}_x0_{x}"; require_ok(run)
        c = json.loads((run / "case.json").read_text())
        E = face_total_melt(run, (2, 3)) - Mb
        res["R_E"][x] = {"x0_actual": c.get("x0_actual"), "x0c": c.get("x0c"), "E": E, "R_E": E / (2 * Es)}
    # face-wide upwelling without a plume: window-mean w at i=0 (post-hoc check from D1)
    p = c["params"]; its = sorted(int(f.split(".")[-2]) for f in glob.glob(str(bg / "WVEL.*.meta")))
    its = [i for i in its if 7200 <= i * p["deltaT"] <= 10800]
    W = np.mean([read_mds(bg / "WVEL", i)[0][0] for i in its], axis=0)
    hf = np.fromfile(bg / "hFacC.data", ">f8").reshape(W.shape)
    w0 = np.where(hf[:, :, 0] > 0, W[:, :, 0], np.nan)
    res["bg_w_face_median_by_depth"] = {f"{-k*p['dz']:.0f}": float(np.nanmedian(w0[k])) for k in range(0, W.shape[0], 3)}
    out[s] = res
    print(f"== {s}: M_bg {Mb:.4f} (1-2h {Mb1:.4f}) eps {eps:.2e}; M_single {Ms:.4f}; E_single {Es:.4f}; "
          f"10eps {10*eps:.2e} -> resolvable {Es >= 10*eps}; Mbg/Ms {Mb/Ms:.3f}")
    for x, r in res["R_E"].items():
        print(f"   x0={x:4d} actual {r['x0_actual']} x0c {r['x0c']}  E {r['E']:.4f}  R_E {r['R_E']:.3f}")
    print("   bg median w at face by z:", {k: f"{v:+.4f}" for k, v in res["bg_w_face_median_by_depth"].items()})
(Path(__file__).resolve().parents[1] / "figures" / "fig2_excess_clean.json").write_text(json.dumps(out, indent=1, default=float))
