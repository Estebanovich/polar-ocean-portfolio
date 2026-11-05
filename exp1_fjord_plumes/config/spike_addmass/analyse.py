"""PD3 spike: read the two runs and print the evidence (sink/source T, S,
volume conservation, flow divergence in the sink cell).
Usage: python analyse.py /Volumes/LaCie/portfolio_runs/exp1/spike_addmass
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from pfcommon.mitgcm_io import read_mds  # noqa: E402

from gendata import DX, DZ, Q, SINK, SRC, T_OUT, S_OUT  # noqa: E402

V = DX * DX * DZ


def series(run, name, cell):
    out = []
    for it in range(0, 301, 30):
        a = read_mds(run / name, it)[0][0]
        out.append(a[cell])
    return np.array(out)


def main(root):
    root = Path(root)
    res = {}
    for v in ("tset", "local"):
        run = root / v
        res[v] = {f"{n}_{c}": series(run, n, cell) for n in ("T", "S")
                  for c, cell in (("sink", SINK), ("src", SRC))}
        eta = np.array([read_mds(run / "Eta", it)[0][0].sum() * DX * DX
                        for it in range(0, 301, 30)])
        k, j, i = SINK
        U = read_mds(run / "U", 300)[0][0]
        Vv = read_mds(run / "V", 300)[0][0]
        W = read_mds(run / "W", 300)[0][0]
        # W(k) is at the top face of cell k (positive up); W(k+1) the bottom
        div = ((U[k, j, i + 1] - U[k, j, i]) * DX * DZ
               + (Vv[k, j + 1, i] - Vv[k, j, i]) * DX * DZ
               + (W[k, j, i] - W[k + 1, j, i]) * DX * DX)
        res[v]["eta_vol"] = eta
        res[v]["outflux_sink"] = div
    t = np.arange(0, 301, 30) * 60.0
    np.set_printoptions(precision=6, suppress=False, linewidth=120)
    for v, r in res.items():
        print(f"== {v}")
        print(" t [h]       ", t / 3600)
        for key in ("T_sink", "S_sink", "T_src", "S_src"):
            print(f" {key:12s}", r[key])
        print(" sum(eta)*A [m^3] min/max", r["eta_vol"].min(), r["eta_vol"].max())
        print(f" net outflow through sink-cell faces at t_end [m^3/s]: {r['outflux_sink']:.6f}"
              f" (expected -Q = {-Q})")
    # removal at temp_addMass, no advection: dT/dt = -(Q/V)(T_out - T)
    lam = Q / V
    T_pred = T_OUT + (2.0 - T_OUT) * np.exp(lam * t)
    S_pred = S_OUT + (34.0 - S_OUT) * np.exp(lam * t)
    print("== prediction if the sink removes water at temp/salt_addMass (no inflow mixing)")
    print(" T_sink_pred ", T_pred)
    print(" S_sink_pred ", S_pred)
    print("== max |local - initial| in sink cell: T %.3e, S %.3e" % (
        np.abs(res["local"]["T_sink"] - 2).max(), np.abs(res["local"]["S_sink"] - 34).max()))
    print("== max |T_src(local) - T_src(tset)| = %.3e, S: %.3e" % (
        np.abs(res["local"]["T_src"] - res["tset"]["T_src"]).max(),
        np.abs(res["local"]["S_src"] - res["tset"]["S_src"]).max()))


if __name__ == "__main__":
    main(sys.argv[1])
