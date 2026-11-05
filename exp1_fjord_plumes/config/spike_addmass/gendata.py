"""PD3 spike: inputs for the addMass sink/source test box (10 x 10 x 10).

Two variants written to <out>/tset and <out>/local:
  tset  : temp_addMass = 10, salt_addMass = 30 (one scalar for all cells),
          no RBCS.
  local : temp_addMass unset and salt_EvPrRn = UNSET_RL so that
          salt_addMass stays unset (ini_parms.F:1579-1580; default
          salt_EvPrRn = 0, set_defaults.F:265, would make it 0 psu) ->
          the addMass T/S term is skipped (apply_forcing.F:505, :875) and
          added/removed water has the local T, S; the outflow T/S is then
          imposed with RBCS in the source cell only.
Uniform T = 2 degC, S = 34 psu, so advection and (zero) diffusion give no
tracer tendency: any change of T, S in the sink cell is the addMass term.
Usage: python gendata.py <out_dir>
"""
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NX = NY = NZ = 10
DX, DZ, RHO0 = 200.0, 10.0, 1027.0
Q = 2.0                      # m^3 s^-1 removed at SINK, added at SRC
SRC = (2, 4, 7)              # (k, j, i), 0-based, k = 0 at the surface
SINK = (7, 4, 2)
T_OUT, S_OUT = 10.0, 30.0


def wbin(path, a):
    np.ascontiguousarray(a, dtype=">f8").tofile(path)


def main(out):
    out = Path(out)
    V = DX * DX * DZ
    for variant in ("tset", "local"):
        d = out / variant
        d.mkdir(parents=True, exist_ok=True)
        wbin(d / "topog.bin", -NZ * DZ * np.ones((NY, NX)))
        wbin(d / "T.init.bin", 2.0 * np.ones((NZ, NY, NX)))
        wbin(d / "S.init.bin", 34.0 * np.ones((NZ, NY, NX)))
        am = np.zeros((NZ, NY, NX))
        am[SRC], am[SINK] = RHO0 * Q, -RHO0 * Q   # kg s^-1 per cell
        wbin(d / "addmass.bin", am)
        parm = (HERE / "input" / "data.common").read_text()
        pkg = " &PACKAGES\n &\n"
        if variant == "tset":
            parm += f" temp_addMass={T_OUT},\n salt_addMass={S_OUT},\n"
        else:
            parm += " salt_EvPrRn=1.234567E+05,\n"
            pkg = " &PACKAGES\n useRBCS=.TRUE.,\n &\n"
            tau = 2.0e5
            mask = np.zeros((NZ, NY, NX))
            mask[SRC] = Q * tau / V          # = 1.0
            wbin(d / "rbcs_mask.bin", mask)
            wbin(d / "rbcs_T.bin", T_OUT * np.ones((NZ, NY, NX)))
            wbin(d / "rbcs_S.bin", S_OUT * np.ones((NZ, NY, NX)))
            shutil.copy(HERE / "input" / "data.rbcs", d / "data.rbcs")
        parm += (HERE / "input" / "data.tail").read_text()
        (d / "data").write_text(f"# PD3 spike, variant {variant}\n" + parm)
        (d / "data.pkg").write_text(pkg)
        shutil.copy(HERE / "input" / "eedata", d / "eedata")


if __name__ == "__main__":
    main(sys.argv[1])
