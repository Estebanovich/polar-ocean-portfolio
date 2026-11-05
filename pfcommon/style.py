"""Common matplotlib style: 9 pt sans, 7.0 x 4.3 in figures, colour-blind safe lines."""
from __future__ import annotations

from pathlib import Path

import cmocean
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGSIZE = (7.0, 4.3)
# Okabe & Ito (2008) colour-blind-safe palette
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#F0E442", "#000000"]
CMAPS = {"temperature": cmocean.cm.thermal, "salinity": cmocean.cm.haline, "diverging": cmocean.cm.balance}


def apply() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "figure.figsize": FIGSIZE,
        "axes.prop_cycle": matplotlib.cycler(color=OKABE_ITO),
        "image.cmap": "cmo.thermal",
        "savefig.dpi": 200,
    })


def save(fig, stem: Path) -> list[Path]:
    """Write `<stem>.png` (200 dpi) and `<stem>.pdf`; return both paths."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    out = [stem.with_suffix(".png"), stem.with_suffix(".pdf")]
    fig.savefig(out[0], dpi=200, bbox_inches="tight")
    fig.savefig(out[1], bbox_inches="tight")
    return out
