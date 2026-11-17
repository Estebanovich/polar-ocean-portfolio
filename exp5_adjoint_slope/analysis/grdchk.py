"""Parse MITgcm gradient-check output (pkg/grdchk, grdchk_print.F).

Each control point prints two ``grad-res`` lines: (proc, #, i, j, k, bi, bj, iobc, fc ref,
fc + eps, fc - eps) and (proc, #, ..., adj or tlm grad, fd grad, 1 - fd/grad). The header
line names ``adj grad`` or ``tlm grad``."""
from __future__ import annotations

import re

_NUM = r"[-+]?\d+\.\d+E[-+]\d+"
# Integer columns are fixed width (grdchk_print.F) and can run together ("124260"), so a row is
# any grad-res line ending in three reals; the integers are read from the first line of a pair.
_ROW = re.compile(rf"^\s*grad-res\s+([\d\s]+?)\s+({_NUM})\s+({_NUM})\s+({_NUM})\s*$", re.M)


def parse(text: str) -> list[dict]:
    """List of points: ind, i, j, k, fc_ref, fc_plus, fc_minus, grad, fd_grad, rel_err
    (= |1 - fd/grad|), kind ('adj' or 'tlm')."""
    kind = "tlm" if "tlm grad" in text else "adj"
    rows = _ROW.findall(text)
    pts = []
    for a, b in zip(rows[0::2], rows[1::2]):
        ints = [int(x) for x in a[0].split()]          # proc, #, i, j, k, bi, bj, iobc
        pts.append({"ind": ints[1], "i": ints[2], "j": ints[3], "k": ints[4],
                    "fc_ref": float(a[1]), "fc_plus": float(a[2]), "fc_minus": float(a[3]),
                    "grad": float(b[1]), "fd_grad": float(b[2]), "rel_err": abs(float(b[3])),
                    "kind": kind})
    return pts
