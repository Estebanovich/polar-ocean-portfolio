"""Plume-excess ratio R_E at dx = 5 m (pre-registered subset, notes/pd1_b_twoplume3d.md, 2026-11-06).

Same metric as analysis/excess_clean_bg.py at 10 m, on the twolayer 5 m runs:
  R_E(x0) = E(x0) / (2 E_single),  E = M - M_bg (face_total_melt, window 2-3 h);
  eps = |M_bg(2-3 h) - M_bg(1-2 h)|; resolvable only if E_single >= 10 eps.
Pre-registered reading: if |R_E(5 m) - R_E(10 m)| > 0.1 for any x0, the plume interaction depends on
resolution and Fig. 2 is presented at 5 m; otherwise "the plumes hardly interact" is robust in this range.
C4 of each two-plume run is reported as context only (z_detach of the y-max of w at the face; theory with
the far-field window-mean ambient, Q0 per plume and, for the merged x0 = 0 source, also 2 Q0).
Output: figures/fig2_excess_dx5.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from exp1_fjord_plumes.analysis.fig2_two_plume_melt import AMBIENT_X_M, theory_zdetach  # noqa: E402
from exp1_fjord_plumes.analysis.plume_diag import (  # noqa: E402
    ambient_far_from_face, face_total_melt, face_w_profile)
from pfcommon.runctl import require_ok  # noqa: E402

R = Path("/Volumes/LaCie/portfolio_runs/exp1")
REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "exp1_fjord_plumes" / "figures" / "fig2_excess_dx5.json"
W, W1 = (2, 3), (1, 2)
X0S = (0, 20, 80)
DIFF_TOL = 0.1


def c4_context(run: Path) -> dict:
    p = json.loads((run / "case.json").read_text())["params"]
    zm = face_w_profile(run, W)["z_detach"]
    amb = ambient_far_from_face(run, W, AMBIENT_X_M)
    a3 = (amb["z"], amb["T"], amb["S"])
    zt = theory_zdetach(p, a3)
    zt2 = theory_zdetach({**p, "Q0": 2 * p["Q0"]}, a3)
    return {"z_detach_model": zm, "z_detach_theory_Q0": zt, "diff_over_H_Q0": (zm - zt) / p["H"],
            "z_detach_theory_2Q0": zt2, "diff_over_H_2Q0": (zm - zt2) / p["H"],
            "label": "context only (pre-registered, no threshold)"}


def main() -> dict:
    ten = json.loads((REPO / "exp1_fjord_plumes" / "figures" / "fig2_excess_clean.json").read_text())["twolayer"]
    bg = R / "b_dx5_bg_twolayer_nodischarge"
    single = R / "b_dx5_twolayer_single"
    require_ok(bg)
    require_ok(single)
    Mb, Mb1 = face_total_melt(bg, W), face_total_melt(bg, W1)
    eps = abs(Mb - Mb1)
    Ms = face_total_melt(single, W)
    Es = Ms - Mb
    res = {"source": "analysis/excess_dx5.py; notes/pd1_b_twoplume3d.md (5 m subset, pre-registered)",
           "M_bg": Mb, "M_bg_1_2h": Mb1, "eps": eps, "M_single": Ms, "E_single": Es,
           "Mbg_over_Msingle": Mb / Ms, "resolvable": bool(Es >= 10 * eps), "R_E": {}}
    for x in X0S:
        run = R / f"b_dx5_twolayer_x0_{x}"
        require_ok(run)
        c = json.loads((run / "case.json").read_text())
        E = face_total_melt(run, W) - Mb
        r5 = E / (2 * Es)
        r10 = ten["R_E"][str(x)]["R_E"]
        res["R_E"][str(x)] = {"x0_actual": c.get("x0_actual"), "E": E, "R_E": r5, "R_E_10m": r10,
                              "diff_5m_minus_10m": r5 - r10, "exceeds_0p1": bool(abs(r5 - r10) > DIFF_TOL),
                              "c4": c4_context(run)}
    res["resolution_dependent"] = any(v["exceeds_0p1"] for v in res["R_E"].values())
    res["reading"] = ("interaction depends on resolution: Fig. 2 at 5 m" if res["resolution_dependent"]
                      else "plumes hardly interact: robust to resolution in this range")
    OUT.write_text(json.dumps(res, indent=1, default=float))
    print(f"M_bg {Mb:.4f} (1-2 h {Mb1:.4f}) eps {eps:.2e}; M_single {Ms:.4f}; E_single {Es:.4f}; "
          f"10 eps {10 * eps:.2e} -> resolvable {res['resolvable']}; Mbg/Ms {Mb / Ms:.3f}")
    for x, v in res["R_E"].items():
        c = v["c4"]
        print(f"x0={x:>3}: E {v['E']:.4f}  R_E 5 m {v['R_E']:.3f}  10 m {v['R_E_10m']:.3f}  "
              f"diff {v['diff_5m_minus_10m']:+.3f}  >0.1 {v['exceeds_0p1']}  | z_det model {c['z_detach_model']:.1f} "
              f"theory(Q0) {c['z_detach_theory_Q0']:.1f} theory(2Q0) {c['z_detach_theory_2Q0']:.1f}")
    print("reading:", res["reading"])
    return res


if __name__ == "__main__":
    main()
