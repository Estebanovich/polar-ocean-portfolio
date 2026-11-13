"""POST HOC diagnostics of the V4 failure (Task 6 part B, written after V4a and V4b failed).

Not a criterion; it backs the (a)-(d) classification in notes/pd2_validacion.md. For F1 (200 m,
budget.fjord_box) and V4 (400 m, validate.coarse_box), per sound: the outward Kang & Fringer
baroclinic pressure work F_bc(t) [W] (validate.face_profiles, + out), its record integral (net)
and the integral of its absolute value (gross) [J], the net-to-gross ratio, the band split of the
net [J] (validate.band_flux) and the number of open faces and of wet levels on them. A net flux
that is a small residual of a large gross flux makes a ratio of nets ill-conditioned.

Usage: pd2_v4_diagnostics.py --f1 RUN --v4 RUN [--out figures/pd2_v4_diagnostics.json]
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


def section_series(fp: dict) -> dict:
    """{section: F_bc(t) [W], + out of the box} from validate.face_profiles output."""
    Fq = -(fp["up"] * fp["pp"] * fp["dz"][None]).sum(1) * fp["w"][None]        # (nt, nf), + out
    return {str(n): Fq[:, fp["group"] == n].sum(1) for n in np.unique(fp["group"])}


def summarise(fp: dict, bands: dict) -> dict:
    t = fp["t"]
    dt = float(t[1] - t[0])
    out = {}
    for n, F in section_series(fp).items():
        sel = fp["group"] == n
        net = float(budget.cumulative(t, F)[-1])
        gross = float(budget.cumulative(t, np.abs(F))[-1])
        b = va.band_flux(fp["up"][:, :, sel], fp["pp"][:, :, sel], fp["dz"][:, sel],
                         fp["w"][sel], dt, bands)
        out[n] = {"net_out_J": net, "gross_J": gross, "net_over_gross": net / gross if gross else None,
                  "band_out_J": {k: -v * t.size * dt for k, v in b.items()},
                  "open_faces": int(sel.sum()), "wet_levels_max": int((fp["dz"][:, sel] > 0).sum(0).max()),
                  "mean_out_W": float(F.mean()), "std_W": float(F.std())}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--f1", type=Path, required=True)
    ap.add_argument("--v4", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=PKG / "figures" / "pd2_v4_diagnostics.json")
    a = ap.parse_args(argv)
    th = va.load_theory()
    bands = va.bands_hz(th)
    d1, d4 = va.open_run(a.f1), va.open_run(a.v4)
    m1, g1 = budget.fjord_box(d1)
    m4, g4 = va.coarse_box(m1, d4["hFacC"].values[0] > 0, d4["XC"].values[0], d4["YC"].values[:, 0],
                           va._sections())
    res = {"post_hoc": "diagnostics written after V4 failed; not criteria"}
    for name, ds, m, g in (("F1", d1, m1, g1), ("V4", d4, m4, g4)):
        res[name] = summarise(va.face_profiles(ds, m, g), bands)
        a.out.write_text(json.dumps(res, indent=1, default=va.jsonable))
    print(json.dumps(res, indent=1, default=va.jsonable))
    return 0


if __name__ == "__main__":
    sys.exit(main())
