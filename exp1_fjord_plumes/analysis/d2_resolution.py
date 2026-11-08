"""D2: C4 of the single (b) plume at dx = 10 m and 5 m (pre-registered, notes/pd1_b_twoplume3d.md, D2).

Same code (fig2_two_plume_melt.c4_for) and window (2-3 h) at both resolutions.
Output: figures/d2_resolution.json.  Usage: d2_resolution.py [--runs-root DIR]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pfcommon.paths import EXP1_RUNS  # noqa: E402
from exp1_fjord_plumes.analysis.fig2_two_plume_melt import WIN_MAIN, c4_for  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402

RUNS = {"dx10": "b_3d_twolayer_single", "dx5": "b_dx5_twolayer_single"}
OUT = Path(__file__).resolve().parents[1] / "figures" / "d2_resolution.json"


def main(root: Path) -> dict:
    res = {"source": "c4_for(run, (2, 3)); notes/pd1_b_twoplume3d.md, D2"}
    for k, name in RUNS.items():
        run = root / name
        require_ok(run)
        c = c4_for(run, WIN_MAIN)
        res[k] = {"run": name, "z_detach_model": c["z_detach_model"], "z_detach_theory": c["z_detach_theory"],
                  "diff_over_H": c["diff_over_H"]}
    OUT.write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs-root", type=Path, default=EXP1_RUNS)
    main(ap.parse_args().runs_root)
